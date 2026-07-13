"""Both plants drive the SAME ISO 3888-1 double lane change - animated, one as a ghost.

Builds a straight 'track' whose racing line is the ISO 3888-1 path (sections
15/30/25/25/15/15 m, 3.5 m offset, cosine transitions) at a constant 80 km/h
reference, runs the OpenModelica TrackTricycle and the JavaScript plant on the
identical table with the identical (unified) driver and Elise parameters, and renders

  outputs/gif/port_dlc.gif       top-down animation, JS car as a translucent ghost
                                 over the Modelica car
  outputs/svg/port_dlc_traj.svg  CG trajectory y(x) of both plants + the reference

Requires OpenModelica (omc) and node on PATH.
"""
import json, os, subprocess, sys, tempfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import animation, transforms

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'tracks'))
from speed_profile import write_track_table
from racing_line import offset_geometry
from cars import build_config, override_string

MOD = os.path.join(HERE, 'modelica')
V0 = 80/3.6                      # ISO 3888-1 recommended entry speed
DS, L = 0.5, 250.0               # straight test track
RUNIN = 40.0                     # maneuver starts here

# ---- 1. the ISO 3888-1 path as a racing line on a straight centerline --------------
s = np.arange(0, L, DS)
x, y, psi, kap = s.copy(), np.zeros_like(s), np.zeros_like(s), np.zeros_like(s)
sec = np.cumsum([0, 15, 30, 25, 25, 15, 15]) + RUNIN     # section boundaries
nRef = np.zeros_like(s)
m2 = (s >= sec[1]) & (s < sec[2])                        # transition out
nRef[m2] = 3.5*(1 - np.cos(np.pi*(s[m2] - sec[1])/30))/2
nRef[(s >= sec[2]) & (s < sec[3])] = 3.5                 # offset lane
m4 = (s >= sec[3]) & (s < sec[4])                        # transition back
nRef[m4] = 3.5*(1 + np.cos(np.pi*(s[m4] - sec[3])/25))/2
psiRef, kapLine, dsSeg = offset_geometry(x, y, psi, DS, nRef)
vRef = np.full_like(s, V0)
axFF = np.zeros_like(s)
dFF = np.zeros_like(s)

# ---- 2. OpenModelica lap ------------------------------------------------------------
cfg = build_config('elise', 'base')
os.makedirs(os.path.join(MOD, 'build'), exist_ok=True)
LTRK = write_track_table(os.path.join(MOD, 'build', 'track.txt'),
                         s, kap, vRef, axFF, nRef=nRef, psiRef=psiRef,
                         kappaLine=kapLine, deltaFF=dFF)
override = override_string(cfg, LTRK, V0) + ',trike.n0=0,trike.dpsi0=0'
mos = ('loadModel(Modelica); loadFile("Tricycle.mo");\n'
       'simulate(Tricycle.Examples.TrackLap, stopTime=11.5, method="rungekutta", '
       'numberOfIntervals=1920, outputFormat="csv", '
       'variableFilter="time|s|n|vKmh|deltaDeg|dpsiDeg", '
       f'fileNamePrefix="build/portdlc", simflags="{override}"); getErrorString();\n')
open('/tmp/portdlc.mos', 'w').write(mos)
res = os.path.join(MOD, 'build', 'portdlc_res.csv')
if os.path.exists(res):
    os.unlink(res)
r = subprocess.run(['omc', '/tmp/portdlc.mos'], cwd=MOD, check=True,
                   capture_output=True, text=True)
if not os.path.exists(res):
    sys.exit('OMC produced no result:\n' + r.stdout[-3000:])
dm = np.genfromtxt(res, delimiter=',', names=True)

# ---- 3. JavaScript lap on the identical table ---------------------------------------
html = open(os.path.join(HERE, 'outputs', 'index.html')).read()
a, b = html.index('const D='), html.index('// ---- world chase camera')
T = dict(s=list(np.round(np.append(s, L), 3)), L=L, ds=DS, n=len(s), width=20.0,
         kap_c=list(np.append(kap, kap[0])), nRef=list(np.round(np.append(nRef, 0), 4)),
         psiRef=list(np.round(np.append(psiRef, 0), 5)),
         kapLine=list(np.round(np.append(kapLine, 0), 6)),
         dsSeg=list(np.round(np.append(dsSeg, DS), 4)),
         deltaFF=list(np.append(dFF, 0)),
         x=list(np.round(np.append(x, L), 2)), y=list(np.append(y, 0.0)),
         psiU=list(np.append(psi, 0.0)))
node_src = html[a:b] + f"""
T = {json.dumps(T)};
const P = Object.assign({{m0: D.presets.elise.m}}, D.presets.elise);
const vRef = T.s.map(_ => {V0});
const dt = 0.006;
let yv = [0, 0, 0, {V0}, 0, 0, 0, 0, 0];
const out = {{t: [], s: [], n: [], dp: [], del: [], v: []}};
let t = 0;
while (yv[0] < T.L - 2 && t < 11.5) {{
  out.t.push(t); out.s.push(yv[0]); out.n.push(yv[1]); out.dp.push(yv[2]);
  out.del.push(yv[6]); out.v.push(yv[3]*3.6);
  const k1 = deriv(yv, P, vRef), y2 = yv.map((v, i) => v + 0.5*dt*k1.dy[i]),
        k2 = deriv(y2, P, vRef), y3 = yv.map((v, i) => v + 0.5*dt*k2.dy[i]),
        k3 = deriv(y3, P, vRef), y4 = yv.map((v, i) => v + dt*k3.dy[i]),
        k4 = deriv(y4, P, vRef);
  yv = yv.map((v, i) => v + dt/6*(k1.dy[i] + 2*k2.dy[i] + 2*k3.dy[i] + k4.dy[i]));
  t += dt;
}}
console.log(JSON.stringify(out));
"""
with tempfile.NamedTemporaryFile('w', suffix='.js', delete=False) as f:
    f.write(node_src)
    tmp = f.name
js = json.loads(subprocess.run(['node', tmp], capture_output=True, text=True,
                               check=True).stdout)
os.unlink(tmp)

# common time grid over the maneuver
mask = (dm['s'] >= 5) & (dm['s'] <= 230)
tg = np.arange(0.3, min(dm['time'][mask][-1], js['t'][-1]) - 0.1, 1/24)
def series(tq, tt, ff):
    return np.interp(tq, tt, ff)
sM, nM = series(tg, dm['time'], dm['s']), series(tg, dm['time'], dm['n'])
hM = np.radians(series(tg, dm['time'], dm['dpsiDeg']))
sJ, nJ = series(tg, js['t'], js['s']), series(tg, js['t'], js['n'])
hJ = series(tg, js['t'], js['dp'])
print(f'DLC done: Modelica exit n = {nM[-1]:+.3f} m, JS exit n = {nJ[-1]:+.3f} m; '
      f'peak |n error| between plants '
      f'{np.max(np.abs(series(tg, dm["time"], dm["n"]) - np.interp(sM, sJ, nJ))):.3f} m')

# ---- 4. CG trajectory figure --------------------------------------------------------
# Same presentation as the original DLC study (dlc_maneuver.py): distance axis
# compressed SQX:1 for readability, ISO 3888-1 cone gates drawn as shaded lanes with
# solid boundaries. Gate widths per ISO for a BW = 1.8 m vehicle.
SQX = 2
BW = 1.8
GATES = [(sec[0], sec[1], 1.1*BW + 0.25, 0.0),   # section 1: entry lane
         (sec[2], sec[3], 1.2*BW + 0.25, 3.5),   # section 3: offset lane
         (sec[4], sec[6], 1.3*BW + 0.25, 0.0)]   # sections 5+6: exit lane
fig, ax = plt.subplots(figsize=(9.5, 3.0))
fig.patch.set_facecolor('white'); ax.set_facecolor('white'); ax.grid(alpha=0.3)
for x0, x1, w, c in GATES:
    ax.fill_between([x0/SQX, x1/SQX], c - w/2, c + w/2, color='tab:blue',
                    alpha=0.10, lw=0)
    ax.plot([x0/SQX, x1/SQX], [c - w/2]*2, 'k-', lw=1.6)
    ax.plot([x0/SQX, x1/SQX], [c + w/2]*2, 'k-', lw=1.6)
ax.plot(s/SQX, nRef, color='#9aa0a6', lw=1.0, ls=':', label='ISO 3888-1 reference path')
ax.plot(sM/SQX, nM, color='#2b6cb0', lw=1.8, label='OpenModelica TrackTricycle')
ax.plot(sJ/SQX, nJ, color='#dd6b20', lw=1.2, ls='--', label='JavaScript port')
ticks = np.arange(RUNIN, sec[-1] + 1, 20)
ax.set_xticks(ticks/SQX); ax.set_xticklabels(f'{v:.0f}' for v in ticks)
ax.set_xlim((RUNIN - 8)/SQX, (sec[-1] + 8)/SQX); ax.set_ylim(-2.2, 5.2)
ax.set_aspect('equal')
ax.set_xlabel(f'distance [m]  (axis compressed {SQX}:1)')
ax.set_ylabel('lateral position [m]')
ax.legend(loc='upper left', fontsize=8)
ax.set_title('ISO 3888-1 double lane change at 80 km/h — CG path vs cone gates, '
             'Elise, both plants', fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(HERE, 'outputs', 'svg', 'port_dlc_traj.svg'),
            bbox_inches='tight')
print('wrote outputs/svg/port_dlc_traj.svg')

# ---- 5. the animation ---------------------------------------------------------------
CARL, CARW = 3.9, 1.7
figA, axA = plt.subplots(figsize=(10, 2.5))
figA.patch.set_facecolor('white'); axA.set_facecolor('#3b3e44')
axA.set_xlim(RUNIN - 12, sec[-1] + 12); axA.set_ylim(-3.2, 6.9)
axA.set_aspect('equal'); axA.set_yticks([])
axA.set_xlabel('x  [m]', fontsize=8); axA.tick_params(labelsize=7)
axA.plot(s, nRef, color='#8f9498', lw=0.9, ls=':')
for x0, x1, w, c in GATES:                     # the ISO cone gates
    axA.fill_between([x0, x1], c - w/2, c + w/2, color='#5a7e9e', alpha=0.30, lw=0)
    axA.plot([x0, x1], [c - w/2]*2, color='#e8e6df', lw=1.3)
    axA.plot([x0, x1], [c + w/2]*2, color='#e8e6df', lw=1.3)
axA.axhline(-2.6, color='#707377', lw=0.8)
axA.axhline(6.3, color='#707377', lw=0.8)
carM = plt.Rectangle((0, 0), CARL, CARW, fc='#f2c40e', ec='#111111', lw=0.8)
carJ = plt.Rectangle((0, 0), CARL, CARW, fc='#e6e7e9', ec='#111111', lw=0.8,
                     alpha=0.45)
axA.add_patch(carM); axA.add_patch(carJ)
lbl = axA.text(0.015, 0.93, '', transform=axA.transAxes, fontsize=8,
               color='#e8e6df', va='top')
axA.text(0.015, 0.10, 'solid: OpenModelica   ghost: JavaScript port',
         transform=axA.transAxes, fontsize=7, color='#c9c7c0')

def place(car, xc, yc, h):
    tr = (transforms.Affine2D().translate(-CARL/2, -CARW/2).rotate(h)
          .translate(xc, yc) + axA.transData)
    car.set_transform(tr)

def frame(i):
    place(carM, sM[i], nM[i], hM[i])
    place(carJ, sJ[i], nJ[i], hJ[i])
    lbl.set_text(f't = {tg[i]:4.1f} s   ISO 3888-1 @ 80 km/h — '
                 f'Elise, both plants, one driver law')
    return carM, carJ, lbl

ani = animation.FuncAnimation(figA, frame, frames=len(tg), blit=True)
os.makedirs(os.path.join(HERE, 'outputs', 'gif'), exist_ok=True)
out = os.path.join(HERE, 'outputs', 'gif', 'port_dlc.gif')
ani.save(out, writer=animation.PillowWriter(fps=24), dpi=90)
print('wrote', out)
