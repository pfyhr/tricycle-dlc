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
from racing_line import min_curvature_line
from fetch_track import TRACKS
from cars import CARS, override_string, ocp_params
from telemetry import load_trace

ap = argparse.ArgumentParser()
ap.add_argument('--track', default='nordschleife', choices=list(TRACKS))
ap.add_argument('--car', default='tourer', choices=list(CARS),
                help='vehicle parameter set (tourer = original 1650 kg/150 kW, '
                     'elise = track-day Lotus Elise tuned to real telemetry)')
ap.add_argument('--line', default='optimal', choices=['optimal', 'center', 'ocp'],
                help='minimum-curvature racing line (optimal), centerline following '
                     '(center), or provably min-time optimal control (ocp; needs casadi)')
ap.add_argument('--width', type=float, default=None,
                help='track width [m] override (default: per-track value)')
args = ap.parse_args()
TRACK, LINE, CAR = args.track, args.line, args.car
CFG = TRACKS[TRACK]
PFX, DISPLAY = CFG['prefix'], CFG['display']
STEM = f'{CAR}_{PFX}'          # car-tagged output stem, e.g. elise_knutstorp
CAR_DISPLAY = CARS[CAR]['display']
CAR_HALF = 0.9      # BW/2, must mirror the vehicle footprint
EDGE_MARGIN = 0.2   # keep the tyres just inside the white line
DRIVER_MARGIN = 0.4 # min-curve line: leave room for the driver's tracking overshoot
                    # (the OCP line pulls in via its speed-dependent corridor instead)

OMC = shutil.which('omc') or '/Users/pontus/opt/openmodelica/bin/omc'
os.environ.setdefault('OPENMODELICAHOME',
                      os.path.dirname(os.path.dirname(os.path.realpath(OMC))))
MOD = os.path.join(HERE, 'modelica')
SVG = os.path.join(HERE, 'outputs/svg')
os.makedirs(os.path.join(MOD, 'build'), exist_ok=True)
os.makedirs(SVG, exist_ok=True)

# ---- vehicle setup (the car's profile params; the same car is pushed to the plant
# via -override below, so the simulated tyre/mass/power match the line we build) ------
SETUP = dict(CARS[CAR]['profile'])
LAP_VARS = ('time|s|n|vKmh|vRefKmh|ayG|axG|deltaDeg|dpsiDeg|yawRateDegS|betaDeg|'
            'FtieL|FtieR|FxR|Pdrive|FzFL|FzFR')

# ---- 1. track table (geometry + racing line + minimum-time speed profile) ---------
s, x, y, psi, kap = load_centerline(os.path.join(HERE, f'tracks/{TRACK}.csv'))
ds = s[1] - s[0]
width = args.width if args.width is not None else CFG.get('width', 10.0)
wMax = max(0.0, width/2 - CAR_HALF - EDGE_MARGIN)

tOcp = None
deltaFF = None
if LINE in ('optimal', 'ocp') and wMax > 0.3:
    # minimum-curvature line first (both the 'optimal' answer and the OCP warm start)
    from racing_line import offset_geometry, _gauss_periodic, apply_driver_margin
    nRef, psiRef, kapLine, dsSeg = min_curvature_line(x, y, psi, kap, ds, wMax)
    if LINE == 'ocp':
        from opt_lap import solve_min_time_dyn
        vMC, _ = speed_profile(s, kapLine, ds_seg=dsSeg, **SETUP)
        # solve on a coarser grid (~10 m spacing, capped ~700 nodes for the long tracks)
        # then interpolate the optimal offset back to the full centerline grid
        stride = max(int(round(10.0/ds)), int(np.ceil(len(s)/700)))
        sc, kc = s[::stride], kap[::stride]
        # SPEED-DEPENDENT corridor: pull the optimal line in from the edge where the car is
        # fast, since it runs wider on exit the faster it goes (a uniform margin can't - tight
        # enough for fast corners over-conservatises slow ones). Margin grows ~v^2 (grip-limited
        # run-wide), from ~0.5 m in slow corners to ~2.1 m at top speed.
        vN = vMC[::stride]
        wOcp = np.clip(wMax - 0.5 - 1.6*(vN/vN.max())**2, 0.3, wMax)
        print(f'{DISPLAY}: solving 3-DOF minimum-time OCP ({len(sc)} nodes'
              f'{f" (every {stride}th of {len(s)})" if stride > 1 else ""}, '
              f'corridor +/-{wOcp.min():.1f}-{wOcp.max():.1f} m)...')
        # grip_frac 0.93 matches the plant's realised peak; w_reg gentles the steer rate
        # so the optimal transitions stay within the preview driver's tracking bandwidth
        nOpt_c, ocp = solve_min_time_dyn(sc, kc, wOcp, p=ocp_params(CAR),
                                         grip_frac=CARS[CAR]['grip_frac'], w_reg=4e-2,
                                         v_init=vMC[::stride], n_init=nRef[::stride],
                                         dpsi_init=psiRef[::stride])
        tOcp = ocp['T']
        print(f'  IPOPT {"converged" if ocp["converged"] else "did NOT converge"}: '
              f'T_opt = {int(tOcp//60)}:{tOcp % 60:04.1f} '
              f'(min-time for the 3-DOF yaw+sideslip vehicle)')
        Lc = sc[-1] + (sc[1] - sc[0])
        nRef = np.interp(s, np.append(sc, Lc), np.append(nOpt_c, nOpt_c[0]))
        # light final smoothing, just to de-ripple the interpolation (the yaw dynamics
        # already make the line smooth - no path-curvature penalty needed)
        nRef = np.clip(_gauss_periodic(nRef, ds, 6.0), -wMax, wMax)
        psiRef, kapLine, dsSeg = offset_geometry(x, y, psi, ds, nRef)
        # dynamic steer feedforward: the OCP's optimal steer minus the kinematic part
        # (the driver adds its own (Lwb+Kus*vx^2)*kapLine), i.e. the sideslip-steer the
        # geometric feedforward misses - this is what lets the real car actually hold
        # the aggressive optimal line (Kapania-Gerdes feedforward+feedback).
        dOcp = np.interp(s, np.append(sc, Lc), np.append(ocp['delta'], ocp['delta'][0]))
        uOcp = np.interp(s, np.append(sc, Lc), np.append(ocp['u'], ocp['u'][0]))
        drv = CARS[CAR]['driver']
        deltaFF = dOcp - (drv['Lwb'] + drv['Kus']*uOcp**2)*kapLine
        # de-ripple and bound the feedforward: keep the useful low-frequency sideslip
        # steer, drop the spikes that would just saturate the road-wheel clamp
        deltaFF = np.clip(_gauss_periodic(deltaFF, ds, 6.0), -0.08, 0.08)
        tag = f'OCP min-time line, corridor +/-{wMax:.1f} m'
    else:
        # inset the geometric line where the car is fast, so the preview driver's
        # tracking overshoot stays on asphalt (see racing_line.apply_driver_margin)
        vMC, _ = speed_profile(s, kapLine, ds_seg=dsSeg, **SETUP)
        nRef, psiRef, kapLine, dsSeg = apply_driver_margin(
            x, y, psi, ds, nRef, wMax, vMC, margin=DRIVER_MARGIN)
        tag = f'min-curvature line, corridor +/-{wMax - DRIVER_MARGIN - 1.4:.1f}-{wMax:.1f} m'
else:
    nRef = psiRef = None
    kapLine, dsSeg = kap, None
    tag = 'centerline'
vRef, axFF = speed_profile(s, kapLine, ds_seg=dsSeg, **SETUP)
LTRK = write_track_table(os.path.join(MOD, 'build/track.txt'), s, kap, vRef, axFF,
                         nRef=nRef, psiRef=psiRef, kappaLine=kapLine, deltaFF=deltaFF)
# lap distance/time along the driven path (line if optimal, else centerline)
pathLen = float(np.sum(dsSeg)) if dsSeg is not None else LTRK
tIdeal = np.sum((dsSeg if dsSeg is not None else np.full(len(s), ds))/vRef)
print(f'{DISPLAY} ({tag}): centerline L = {LTRK:.0f} m, path {pathLen:.0f} m; '
      f'quasi-steady ideal {int(tIdeal//60)}:{tIdeal % 60:04.1f} '
      f'(driver will be a bit slower)')

# ---- 2. simulate -------------------------------------------------------------------
override = override_string(CAR, LTRK, vRef[0])
mos = (f'loadModel(Modelica); loadFile("Tricycle.mo");\n'
       f'simulate(Tricycle.Examples.TrackLap, stopTime=900, '
       f'numberOfIntervals=18000, outputFormat="csv", variableFilter="{LAP_VARS}", '
       f'fileNamePrefix="build/{STEM}_lap", '
       f'simflags="{override}"); getErrorString();\n')
open('/tmp/trike_nslap.mos', 'w').write(mos)
r = subprocess.run([OMC, '/tmp/trike_nslap.mos'], cwd=MOD, check=True,
                   capture_output=True, text=True)
RES = os.path.join(MOD, f'build/{STEM}_lap_res.csv')
if not os.path.exists(RES):
    sys.exit('simulation produced no result file:\n' + r.stdout[-2000:])
d = np.genfromtxt(RES, delimiter=',', names=True)
t = d['time']
if d['s'][-1] < LTRK - 10:
    sys.exit(f"lap NOT completed: s ended at {d['s'][-1]:.0f} of {LTRK:.0f} m")
tLap = t[-1]
# racing-line offset target at the car's current station, and how well it is tracked
nRefT = np.zeros_like(d['s']) if nRef is None else \
    np.interp(np.mod(d['s'], LTRK), np.append(s, LTRK),
              np.append(nRef, nRef[0]))
trackRms = np.sqrt(((d['n'] - nRefT)[t > 5]**2).mean())
print(f"lap time {int(tLap//60)}:{tLap % 60:04.1f}  "
      f"(vmax {d['vKmh'].max():.0f} km/h, vmin {d['vKmh'][t > 5].min():.0f} km/h, "
      f"line-tracking rms {trackRms:.2f} m, |n|max {np.abs(d['n']).max():.2f} m, "
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
with open(os.path.join(HERE, f'outputs/{STEM}_lap_summary.csv'), 'w') as f:
    f.write('quantity,value,unit\n')
    f.write(f'track,{DISPLAY},-\n')
    f.write(f'line,{LINE},-\n')
    if tOcp is not None:
        f.write(f'ocp_min_time_lower_bound,{tOcp:.1f},s\n')
    f.write(f'lap_time,{tLap:.1f},s\n')
    f.write(f'track_length,{LTRK:.1f},m\n')
    f.write(f'path_length,{pathLen:.1f},m\n')
    f.write(f'ideal_quasi_steady_lap,{tIdeal:.1f},s\n')
    f.write(f'v_max,{d["vKmh"].max():.1f},km/h\n')
    f.write(f'v_min,{d["vKmh"][t > 5].min():.1f},km/h\n')
    f.write(f'abs_n_max,{np.abs(d["n"]).max():.2f},m\n')
    f.write(f'line_tracking_rms,{trackRms:.2f},m\n')
    f.write(f'ay_peak,{np.abs(d["ayG"]).max():.2f},g\n')
    f.write(f'brake_peak,{-d["axG"].min():.2f},g\n')
    f.write(f'tie_rod_peak,{max(np.abs(d["FtieL"]).max(), np.abs(d["FtieR"]).max()):.0f},N\n')
    f.write(f'drive_power_peak,{d["Pdrive"].max()/1e3:.0f},kW\n')

# ---- 4. map colored by speed, with track edges and centerline ----------------------
fig, a = plt.subplots(figsize=(9, 7.2))
nx, ny = -np.sin(psi), np.cos(psi)
hw = width/2
for sgn in (+1, -1):                        # track edges (asphalt corridor)
    a.plot(x + sgn*hw*nx, y + sgn*hw*ny, color='0.7', lw=0.6)
a.plot(np.append(x, x[0]), np.append(y, y[0]), color='0.55', lw=0.7, ls=(0, (6, 6)),
       label='centerline')
pts = np.column_stack([xCar, yCar]).reshape(-1, 1, 2)
segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
LINE_LABEL = {'optimal': 'min-curvature line', 'ocp': 'min-time (OCP) line',
              'center': 'centerline follow'}[LINE]
lc = LineCollection(segs, cmap='viridis', array=d['vKmh'][:-1], lw=2.2,
                    label=LINE_LABEL)
a.add_collection(lc)
a.plot(xCar[0], yCar[0], 'o', color='tab:red', ms=6, zorder=5)
a.annotate('start/finish', (xCar[0], yCar[0]), textcoords='offset points',
           xytext=(8, 6), fontsize=8)
cb = fig.colorbar(lc, ax=a, shrink=0.8); cb.set_label('speed [km/h]')
a.set_xlim(xCar.min() - 200, xCar.max() + 200)
a.set_ylim(yCar.min() - 200, yCar.max() + 200)
a.set_aspect('equal'); a.grid(alpha=.3); a.legend(loc='best', fontsize=8)
a.set_xlabel('x [m]'); a.set_ylabel('y [m]')
tag_ocp = f' (OCP bound {int(tOcp//60)}:{tOcp % 60:04.1f})' if tOcp is not None else ''
a.set_title(f'{DISPLAY} — {CAR_DISPLAY} — {LINE_LABEL}\nlap {int(tLap//60)}:{tLap % 60:04.1f}'
            f'{tag_ocp} at {SETUP["Pmax"]/1e3:.0f} kW / {SETUP["m"]:.0f} kg')
plt.tight_layout(); plt.savefig(f'{SVG}/{STEM}_map.svg', facecolor='white'); plt.close()

# ---- 5. speed trace vs reference + lateral offset ---------------------------------
fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True,
                       gridspec_kw={'height_ratios': [2.4, 1]})
ax[0].plot(d['s']/1000, d['vRefKmh'], 'k--', lw=0.9, label='minimum-time reference (previewed)')
ax[0].plot(d['s']/1000, d['vKmh'], color='tab:blue', lw=1.2, label='vehicle')
tel = load_trace(CAR, TRACK)          # real logger trace, if we have one for this car+track
if tel is not None:
    frac, vReal = tel
    ax[0].plot(frac*LTRK/1000, vReal, color='tab:orange', lw=1.1, ls=(0, (5, 2)),
               alpha=0.9, label='real telemetry (fastest logged lap)')
ax[0].set_ylabel('speed [km/h]'); ax[0].legend(loc='lower right', fontsize=8)
ax[0].grid(alpha=.3)
ax[0].set_title('speed tracking and lateral offset over the lap')
if nRef is not None:
    ax[1].plot(d['s']/1000, nRefT, 'k--', lw=0.9, label='racing-line target')
    ax[1].fill_between(s/1000, -wMax, wMax, color='0.85', lw=0, zorder=0,
                       label='corridor')
ax[1].plot(d['s']/1000, d['n'], color='tab:red', lw=1, label='vehicle')
ax[1].axhline(0, color='k', lw=.6)
ax[1].set_ylabel('offset n [m]'); ax[1].set_xlabel('distance s [km]')
ax[1].grid(alpha=.3); ax[1].legend(loc='upper right', fontsize=7, ncol=3)
plt.tight_layout(); plt.savefig(f'{SVG}/{STEM}_speed.svg', facecolor='white'); plt.close()

# ---- 6. g-g diagram ----------------------------------------------------------------
fig, a = plt.subplots(figsize=(5.6, 5.6))
a.scatter(d['ayG'], d['axG'], s=2, c=d['vKmh'], cmap='viridis', alpha=0.5)
th = np.linspace(0, 2*np.pi, 200)
a.plot(0.855*np.cos(th), 0.855*np.sin(th), 'k--', lw=0.8,
       label='profile envelope 0.9·μ·g')
a.set_xlabel('lateral acceleration [g]'); a.set_ylabel('longitudinal acceleration [g]')
a.set_aspect('equal'); a.grid(alpha=.3); a.legend(fontsize=8, loc='upper right')
a.set_title('g-g diagram over the lap (color: speed)')
plt.tight_layout(); plt.savefig(f'{SVG}/{STEM}_gg.svg', facecolor='white'); plt.close()

# ---- 7. tie-rod forces over the lap ------------------------------------------------
fig, a = plt.subplots(figsize=(10, 4))
a.plot(d['s']/1000, d['FtieL']/1000, color='tab:green', lw=0.9, label='tie rod left')
a.plot(d['s']/1000, d['FtieR']/1000, color='tab:red', lw=0.9, alpha=0.8, label='tie rod right')
a.set_xlabel('distance s [km]'); a.set_ylabel('tie-rod force [kN]')
a.grid(alpha=.3); a.legend(fontsize=8)
a.set_title(f'steering-link loads over one {DISPLAY} lap')
plt.tight_layout(); plt.savefig(f'{SVG}/{STEM}_tierod.svg', facecolor='white'); plt.close()

print(f'wrote outputs/{STEM}_lap_summary.csv and 4 SVGs to outputs/svg/')
print(f'renderer: python3 track_render.py --track={TRACK} --car={CAR}')
