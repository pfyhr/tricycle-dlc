"""Minimum-manageable-time lap of a planar OSM circuit with the power- and
grip-limited tricycle: generate the track + speed-profile table for the vehicle setup,
run Tricycle.Examples.TrackLap, and produce the lap figures (track map colored
by speed, speed trace vs reference, g-g diagram, tie-rod forces) plus
outputs/<track>_lap_summary.csv. The result CSV is kept in
modelica/build/<track>_lap_res.csv for track_render.py (chase-camera viewer).

Run:  python3 track_lap.py [--track=nordschleife|knutstorp|anderstorp|gelleras|kinnekulle]
"""
import argparse, os, shutil, subprocess, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'tracks'))
from speed_profile import load_centerline, speed_profile, write_track_table
from fetch_track import TRACKS

ap = argparse.ArgumentParser()
ap.add_argument('--track', default='nordschleife', choices=list(TRACKS))
TRACK = ap.parse_args().track
CFG = TRACKS[TRACK]
PFX, DISPLAY = CFG['prefix'], CFG['display']

OMC = shutil.which('omc') or '/Users/pontus/opt/openmodelica/bin/omc'
os.environ.setdefault('OPENMODELICAHOME',
                      os.path.dirname(os.path.dirname(os.path.realpath(OMC))))
MOD = os.path.join(HERE, 'modelica')
SVG = os.path.join(HERE, 'outputs/svg')
os.makedirs(os.path.join(MOD, 'build'), exist_ok=True)
os.makedirs(SVG, exist_ok=True)

# ---- vehicle setup (must mirror Tricycle.Track.TrackTricycle defaults) ------------
SETUP = dict(m=1650.0, Pmax=150e3, CdA=0.72, Crr=0.012, rho=1.20,
             mu=0.95, ayFrac=0.90, hcg=0.55, a=1.20, b=1.60)
LAP_VARS = ('time|s|n|vKmh|vRefKmh|ayG|axG|deltaDeg|dpsiDeg|yawRateDegS|betaDeg|'
            'FtieL|FtieR|FxR|Pdrive|FzFL|FzFR')

# ---- 1. track table (geometry + minimum-time speed profile for this setup) --------
s, x, y, psi, kap = load_centerline(os.path.join(HERE, f'tracks/{TRACK}.csv'))
vRef, axFF = speed_profile(s, kap, **SETUP)
LTRK = write_track_table(os.path.join(MOD, 'build/track.txt'), s, kap, vRef, axFF)
ds = s[1] - s[0]
tIdeal = np.sum(ds/vRef)
print(f'{DISPLAY}: L = {LTRK:.0f} m; quasi-steady ideal lap '
      f'{int(tIdeal//60)}:{tIdeal % 60:04.1f} (driver will be a bit slower)')

# ---- 2. simulate -------------------------------------------------------------------
mos = (f'loadModel(Modelica); loadFile("Tricycle.mo");\n'
       f'simulate(Tricycle.Examples.TrackLap, stopTime=900, '
       f'numberOfIntervals=18000, outputFormat="csv", variableFilter="{LAP_VARS}", '
       f'fileNamePrefix="build/{PFX}_lap", '
       f'simflags="-override sLap={LTRK:.1f},u0={vRef[0]:.2f}"); getErrorString();\n')
open('/tmp/trike_nslap.mos', 'w').write(mos)
r = subprocess.run([OMC, '/tmp/trike_nslap.mos'], cwd=MOD, check=True,
                   capture_output=True, text=True)
RES = os.path.join(MOD, f'build/{PFX}_lap_res.csv')
if not os.path.exists(RES):
    sys.exit('simulation produced no result file:\n' + r.stdout[-2000:])
d = np.genfromtxt(RES, delimiter=',', names=True)
t = d['time']
if d['s'][-1] < LTRK - 10:
    sys.exit(f"lap NOT completed: s ended at {d['s'][-1]:.0f} of {LTRK:.0f} m")
tLap = t[-1]
print(f"lap time {int(tLap//60)}:{tLap % 60:04.1f}  "
      f"(vmax {d['vKmh'].max():.0f} km/h, vmin {d['vKmh'][t > 5].min():.0f} km/h, "
      f"|n|max {np.abs(d['n']).max():.2f} m, "
      f"peak ay {np.abs(d['ayG']).max():.2f} g, peak brake {-d['axG'].min():.2f} g, "
      f"peak tie-rod {max(np.abs(d['FtieL']).max(), np.abs(d['FtieR']).max()):.0f} N)")

# centerline interpolators for mapping (s,n) -> global x,y
sq = np.concatenate([s, [LTRK]])
per = lambda f: np.concatenate([f, [f[0]]])
xC, yC, psiC = (lambda ss: np.interp(ss, sq, per(x)),
                lambda ss: np.interp(ss, sq, per(y)),
                None)
psiU = per(psi)  # already unwrapped and periodic-consistent (ends 2*pi apart)
sMod = np.mod(d['s'], LTRK)
xCar = np.interp(sMod, sq, per(x)) - d['n']*np.sin(np.interp(sMod, sq, psiU))
yCar = np.interp(sMod, sq, per(y)) + d['n']*np.cos(np.interp(sMod, sq, psiU))

# ---- 3. summary CSV ----------------------------------------------------------------
with open(os.path.join(HERE, f'outputs/{PFX}_lap_summary.csv'), 'w') as f:
    f.write('quantity,value,unit\n')
    f.write(f'track,{DISPLAY},-\n')
    f.write(f'lap_time,{tLap:.1f},s\n')
    f.write(f'track_length,{LTRK:.1f},m\n')
    f.write(f'ideal_quasi_steady_lap,{tIdeal:.1f},s\n')
    f.write(f'v_max,{d["vKmh"].max():.1f},km/h\n')
    f.write(f'v_min,{d["vKmh"][t > 5].min():.1f},km/h\n')
    f.write(f'abs_n_max,{np.abs(d["n"]).max():.2f},m\n')
    f.write(f'ay_peak,{np.abs(d["ayG"]).max():.2f},g\n')
    f.write(f'brake_peak,{-d["axG"].min():.2f},g\n')
    f.write(f'tie_rod_peak,{max(np.abs(d["FtieL"]).max(), np.abs(d["FtieR"]).max()):.0f},N\n')
    f.write(f'drive_power_peak,{d["Pdrive"].max()/1e3:.0f},kW\n')

# ---- 4. map colored by speed -------------------------------------------------------
fig, a = plt.subplots(figsize=(9, 7.2))
pts = np.column_stack([xCar, yCar]).reshape(-1, 1, 2)
segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
lc = LineCollection(segs, cmap='viridis', array=d['vKmh'][:-1], lw=2.2)
a.add_collection(lc)
a.plot(xCar[0], yCar[0], 'o', color='tab:red', ms=6, zorder=5)
a.annotate('start/finish', (xCar[0], yCar[0]), textcoords='offset points',
           xytext=(8, 6), fontsize=8)
cb = fig.colorbar(lc, ax=a, shrink=0.8); cb.set_label('speed [km/h]')
a.set_xlim(xCar.min() - 200, xCar.max() + 200)
a.set_ylim(yCar.min() - 200, yCar.max() + 200)
a.set_aspect('equal'); a.grid(alpha=.3)
a.set_xlabel('x [m]'); a.set_ylabel('y [m]')
a.set_title(f'{DISPLAY} (planar, OSM centerline) — lap {int(tLap//60)}:{tLap % 60:04.1f} '
            f'at {SETUP["Pmax"]/1e3:.0f} kW / {SETUP["m"]:.0f} kg')
plt.tight_layout(); plt.savefig(f'{SVG}/{PFX}_map.svg', facecolor='white'); plt.close()

# ---- 5. speed trace vs reference + lateral offset ---------------------------------
fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True,
                       gridspec_kw={'height_ratios': [2.4, 1]})
ax[0].plot(d['s']/1000, d['vRefKmh'], 'k--', lw=0.9, label='minimum-time reference (previewed)')
ax[0].plot(d['s']/1000, d['vKmh'], color='tab:blue', lw=1.2, label='vehicle')
ax[0].set_ylabel('speed [km/h]'); ax[0].legend(loc='lower right', fontsize=8)
ax[0].grid(alpha=.3)
ax[0].set_title('speed tracking and lateral offset over the lap')
ax[1].plot(d['s']/1000, d['n'], color='tab:red', lw=1)
ax[1].axhline(0, color='k', lw=.6)
ax[1].set_ylabel('offset n [m]'); ax[1].set_xlabel('distance s [km]')
ax[1].grid(alpha=.3)
plt.tight_layout(); plt.savefig(f'{SVG}/{PFX}_speed.svg', facecolor='white'); plt.close()

# ---- 6. g-g diagram ----------------------------------------------------------------
fig, a = plt.subplots(figsize=(5.6, 5.6))
a.scatter(d['ayG'], d['axG'], s=2, c=d['vKmh'], cmap='viridis', alpha=0.5)
th = np.linspace(0, 2*np.pi, 200)
a.plot(0.855*np.cos(th), 0.855*np.sin(th), 'k--', lw=0.8,
       label='profile envelope 0.9·μ·g')
a.set_xlabel('lateral acceleration [g]'); a.set_ylabel('longitudinal acceleration [g]')
a.set_aspect('equal'); a.grid(alpha=.3); a.legend(fontsize=8, loc='upper right')
a.set_title('g-g diagram over the lap (color: speed)')
plt.tight_layout(); plt.savefig(f'{SVG}/{PFX}_gg.svg', facecolor='white'); plt.close()

# ---- 7. tie-rod forces over the lap ------------------------------------------------
fig, a = plt.subplots(figsize=(10, 4))
a.plot(d['s']/1000, d['FtieL']/1000, color='tab:green', lw=0.9, label='tie rod left')
a.plot(d['s']/1000, d['FtieR']/1000, color='tab:red', lw=0.9, alpha=0.8, label='tie rod right')
a.set_xlabel('distance s [km]'); a.set_ylabel('tie-rod force [kN]')
a.grid(alpha=.3); a.legend(fontsize=8)
a.set_title(f'steering-link loads over one {DISPLAY} lap')
plt.tight_layout(); plt.savefig(f'{SVG}/{PFX}_tierod.svg', facecolor='white'); plt.close()

print(f'wrote outputs/{PFX}_lap_summary.csv and 4 SVGs to outputs/svg/')
print(f'renderer: python3 track_render.py --track={TRACK}')
