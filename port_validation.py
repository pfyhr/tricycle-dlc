"""Is the JavaScript port faithful? Run the SAME car on the SAME line in both plants.

Takes the baked Knutstorp tables (geometry, racing line, steer feedforward) and the
Elise's reference speed straight out of the built live simulator (outputs/index.html),
writes them as the OpenModelica track table, and laps Tricycle.Examples.TrackLap with
the Elise parameter set. The JS plant laps the identical inputs via node. Both closed
loops then face the same line, the same v_ref and the same car - what differs is the
implementation (and the documented deltas: the JS dropped tyre relaxation, and its
driver gained the telemetry-calibrated late-hard braking the Modelica TrackDriver
never had, so braking points differ slightly by design).

Writes outputs/svg/port_validation.svg and prints the agreement metrics quoted in the
README. Requires OpenModelica (omc) and node on PATH.
"""
import json, os, re, subprocess, sys, tempfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'tracks'))
from speed_profile import write_track_table
from cars import build_config, override_string

TRACK, CAR = 'knutstorp', 'elise'
MOD = os.path.join(HERE, 'modelica')
OMC = 'omc'

# ---- 1. pull the baked tables + the JS lap out of the built simulator --------------
html = open(os.path.join(HERE, 'outputs', 'index.html')).read()
a = html.index('const D=')
b = html.index('// ---- world chase camera')
block = html[a:b]

node_src = block + f"""
T = D.tracks.{TRACK};
const P = Object.assign({{m0: D.presets.{CAR}.m}}, D.presets.{CAR});
const vRef = profile(P);
const R = simulate(P);
const o = R.o, out = {{v: [], del: [], n: [], s: []}};
for (let i = 0; i < o.s.length; i += 5) {{     // thin to ~2700 samples
  const s = o.s[i], xc = interp(T.x, s), yc = interp(T.y, s), pc = interp(T.psiU, s);
  out.s.push(s); out.v.push(o.v[i]); out.del.push(o.del[i]);
  out.n.push((o.y[i] - yc)*Math.cos(pc) - (o.x[i] - xc)*Math.sin(pc));
}}
console.log(JSON.stringify({{track: {{s: T.s, kap_c: T.kap_c, nRef: T.nRef, psiRef: T.psiRef,
  kapLine: T.kapLine, dsSeg: T.dsSeg, deltaFF: T.deltaFF, L: T.L, n: T.n}},
  vRef: vRef, lap: R.lap, js: out}}));
"""
with tempfile.NamedTemporaryFile('w', suffix='.js', delete=False) as f:
    f.write(node_src)
    tmp = f.name
data = json.loads(subprocess.run(['node', tmp], capture_output=True, text=True,
                                 check=True).stdout)
os.unlink(tmp)
trk, js, vRef = data['track'], data['js'], np.array(data['vRef'])
N = trk['n']
s = np.array(trk['s'][:N]);      kap = np.array(trk['kap_c'][:N])
nRef = np.array(trk['nRef'][:N]); psiRef = np.array(trk['psiRef'][:N])
kapL = np.array(trk['kapLine'][:N]); dsSeg = np.array(trk['dsSeg'][:N])
dFF = np.array(trk['deltaFF'][:N]); L = trk['L']
v = vRef[:N]
axFF = (np.roll(v, -1)**2 - v**2)/(2*dsSeg)   # feedforward accel between knots

# ---- 2. lap the OpenModelica plant on the identical table --------------------------
cfg = build_config(CAR, 'base')
# NOTE the drivers are compared as they ship: the Modelica TrackDriver runs its original
# gains (raising them to the JS's KMUL-multiplied level reproduces the straight-line
# steering limit cycle the JS driver later fixed with its lookahead cap + gain knee -
# a bug the Modelica side still has at high gain, which is itself evidence the two
# implementations share their dynamics).
os.makedirs(os.path.join(MOD, 'build'), exist_ok=True)
LTRK = write_track_table(os.path.join(MOD, 'build', 'track.txt'),
                         s, kap, v, axFF, nRef=nRef, psiRef=psiRef,
                         kappaLine=kapL, deltaFF=dFF)
override = override_string(cfg, LTRK, v[0])
mos = ('loadModel(Modelica); loadFile("Tricycle.mo");\n'
       'simulate(Tricycle.Examples.TrackLap, stopTime=300, numberOfIntervals=12000, '
       'outputFormat="csv", variableFilter="time|s|n|vKmh|deltaDeg", '
       f'fileNamePrefix="build/portval", simflags="{override}"); getErrorString();\n')
open('/tmp/portval.mos', 'w').write(mos)
subprocess.run([OMC, '/tmp/portval.mos'], cwd=MOD, check=True,
               capture_output=True, text=True)
d = np.genfromtxt(os.path.join(MOD, 'build', 'portval_res.csv'),
                  delimiter=',', names=True)
m = d['s'] <= L
tMo = d['time'][m][-1]

# ---- 3. resample both laps onto a common s grid and compare ------------------------
sg = np.linspace(0, L - 5, 900)
def at(sq, ss, ff):
    return np.interp(sq, ss, ff)
vJ, vM = at(sg, js['s'], js['v']), at(sg, d['s'][m], d['vKmh'][m])
dJ = at(sg, js['s'], np.degrees(js['del']))
dM = at(sg, d['s'][m], d['deltaDeg'][m])
nJ, nM = at(sg, js['s'], js['n']), at(sg, d['s'][m], d['n'][m])
rms = lambda e: float(np.sqrt(np.mean(e**2)))
print(f'lap time: Modelica {tMo:.1f} s vs JS {data["lap"]:.1f} s '
      f'(Δ {abs(tMo - data["lap"]):.1f} s)')
print(f'v(s)  rms {rms(vJ - vM):.1f} km/h   delta(s) rms {rms(dJ - dM):.2f} deg   '
      f'n(s) rms {rms(nJ - nM):.2f} m')

# ---- 4. the figure ------------------------------------------------------------------
fig, ax = plt.subplots(3, 1, figsize=(9.5, 7.2), sharex=True)
fig.patch.set_facecolor('white')
for a_ in ax:
    a_.set_facecolor('white'); a_.grid(alpha=0.25)
ax[0].plot(sg, vM, lw=1.6, color='#2b6cb0', label='OpenModelica TrackTricycle')
ax[0].plot(sg, vJ, lw=1.2, color='#dd6b20', ls='--', label='JavaScript port (live sim)')
ax[0].set_ylabel('v  [km/h]'); ax[0].legend(loc='lower right', fontsize=9)
ax[1].plot(sg, dM, lw=1.6, color='#2b6cb0')
ax[1].plot(sg, dJ, lw=1.2, color='#dd6b20', ls='--')
ax[1].set_ylabel('road-wheel δ  [°]')
ax[2].plot(sg, nM, lw=1.6, color='#2b6cb0')
ax[2].plot(sg, nJ, lw=1.2, color='#dd6b20', ls='--')
ax[2].set_ylabel('lateral offset n  [m]'); ax[2].set_xlabel('distance s  [m]')
fig.suptitle(f'Same car (Elise), same racing line + v_ref, both plants — Knutstorp\n'
             f'Modelica {tMo:.1f} s  vs  JS {data["lap"]:.1f} s;  '
             f'v rms {rms(vJ - vM):.1f} km/h,  δ rms {rms(dJ - dM):.2f}°,  '
             f'n rms {rms(nJ - nM):.2f} m', fontsize=10)
fig.tight_layout()
out = os.path.join(HERE, 'outputs', 'svg', 'port_validation.svg')
fig.savefig(out, bbox_inches='tight')
print('wrote', out)
