"""ISO 3888-1 double lane change with the tricycle vehicle and a manual (unassisted)
rack-and-pinion steering: run Tricycle.Examples.DoubleLaneChange plus a StepSteer speed
sweep, and produce the link-force plots (tie-rod forces, kingpin-moment trail
decomposition), the steering-feel plot (handwheel angle/torque), model-validation plots
(brush tire curves, understeer gradient) and outputs/dlc_summary.csv."""
import numpy as np, subprocess, shutil, os, csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OMC = shutil.which('omc') or '/Users/pontus/opt/openmodelica/bin/omc'
# pin OPENMODELICAHOME to the omc install so generated code always finds its headers
# (openmodelica.h / omc_simulation_settings.h), even via symlinks or relocated installs
os.environ.setdefault('OPENMODELICAHOME',
                      os.path.dirname(os.path.dirname(os.path.realpath(OMC))))
HERE = os.path.dirname(os.path.abspath(__file__)); MOD = os.path.join(HERE, 'modelica')
SVG = os.path.join(HERE, 'outputs/svg')
GIF = os.path.join(HERE, 'outputs/gif')

SPEEDS_KMH = [40, 60, 80, 100, 120]
DLC_VARS = ('time|X|Y|psiDeg|ayG|yawRateDegS|deltaLdeg|deltaRdeg|dCmdDeg|'
            'FtieL|FtieR|MkpL|MkpR|MkpMechL|MkpMechR|MkpPneuL|MkpPneuR|rackForce|'
            'FyFL|FyFR|FyR|FzFL|FzFR|MzFL|MzFR|alphaFLdeg|alphaFRdeg|alphaRdeg|'
            'betaDeg|hwaDeg|hwTorque|rackDisp')
STEP_VARS = ('time|yawRateDegS|ayG|deltaLdeg|betaDeg|FtieL|FtieR|FzFL|FzFR|'
             'alphaFLdeg|FyFL|MzFL')

# ---- vehicle / tire parameters (must mirror Tricycle defaults) -------------------
M, IZZ, A, B = 1650.0, 2700.0, 1.20, 1.60
L = A + B; G = 9.80665
TF, HCG, XIF = 1.55, 0.55, 0.6
TMECH, LARM = 0.025, 0.11
IS = 20.0                                       # overall steering ratio
TIRE_F = dict(c1=7.0e4, c2=4000.0, mu=0.95, ap0=0.06, FzNom=4000.0)
TIRE_R = dict(c1=1.4e5, c2=8000.0, mu=0.95, ap0=0.085, FzNom=8000.0)
BW = 1.8   # vehicle width for the ISO gate drawing


def brush(alpha, Fz, d):
    """NumPy twin of Tricycle.brushForces (same regularizations)."""
    Fzc = np.maximum(Fz, 50.0)
    Ca = d['c1']*np.sin(2*np.arctan(Fzc/d['c2']))
    ap = d['ap0']*np.sqrt(Fzc/d['FzNom'])
    sy = np.tan(alpha)
    sAbs = np.sqrt(sy**2 + 1e-12)
    sgn = sy/sAbs
    th = Ca/(3*d['mu']*Fzc)
    uu = 0.5*(th*sAbs + 1 - np.sqrt((th*sAbs - 1)**2 + 1e-4))
    lam = 1 - uu
    Fy = d['mu']*Fzc*(1 - lam**3)*sgn
    Mz = -d['mu']*Fzc*ap*uu*lam**3*sgn
    return Fy, Mz


def yref_iso(x):
    """ISO 3888-1 reference centerline (mirror of Tricycle.yRefIso3888)."""
    x = np.asarray(x, float)
    y = np.zeros_like(x)
    m = (x > 15) & (x <= 45); y[m] = 1.75*(1 - np.cos(np.pi*(x[m] - 15)/30))
    y[(x > 45) & (x <= 70)] = 3.5
    m = (x > 70) & (x <= 95); y[m] = 1.75*(1 + np.cos(np.pi*(x[m] - 70)/25))
    return y


GATES = [(0, 15, 1.1*BW + 0.25, 0.0),     # section 1: entry lane
         (45, 70, 1.2*BW + 0.25, 3.5),    # section 3: offset lane
         (95, 125, 1.3*BW + 0.25, 0.0)]   # sections 5+6: exit lane (ISO total 125 m)

# ---- simulate --------------------------------------------------------------------
lines = ['loadModel(Modelica); loadFile("Tricycle.mo");']
runs = [('dlc', 'dlc', '/tmp/trikedlc.csv', '', 8, 8000,
         'Tricycle.Examples.DoubleLaneChange', DLC_VARS)]
lines.append(
    f'simulate(Tricycle.Examples.DoubleLaneChange, stopTime=8, numberOfIntervals=8000, '
    f'outputFormat="csv", variableFilter="{DLC_VARS}", fileNamePrefix="build/dlc", '
    f'simflags="-r=/tmp/trikedlc.csv");')
for kmh in SPEEDS_KMH:
    u = kmh/3.6
    rf = f'/tmp/trikestep_{kmh}.csv'
    lines.append(
        f'simulate(Tricycle.Examples.StepSteer, stopTime=6, numberOfIntervals=3000, '
        f'outputFormat="csv", variableFilter="{STEP_VARS}", fileNamePrefix="build/step", '
        f'simflags="-override=u={u:.4f} -r={rf}");')
    runs.append(('step', kmh, rf, f'-override=u={u:.4f}', 6, 3000,
                 'Tricycle.Examples.StepSteer', STEP_VARS))
open('/tmp/trike_all.mos', 'w').write('\n'.join(lines))
subprocess.run([OMC, '/tmp/trike_all.mos'], cwd=MOD, check=True, capture_output=True)

# retry any result file that didn't get written in the batch run
for kind, key, rf, ov, stop, n, model, vf in runs:
    if not os.path.exists(rf):
        m = ('loadModel(Modelica); loadFile("Tricycle.mo");\n'
             f'simulate({model}, stopTime={stop}, numberOfIntervals={n}, '
             f'outputFormat="csv", variableFilter="{vf}", fileNamePrefix="build/retry", '
             f'simflags="{ov} -r={rf}");')
        open('/tmp/trike_retry.mos', 'w').write(m)
        subprocess.run([OMC, '/tmp/trike_retry.mos'], cwd=MOD, check=True, capture_output=True)

d = np.genfromtxt('/tmp/trikedlc.csv', delimiter=',', names=True)
S = {k: np.genfromtxt(f'/tmp/trikestep_{k}.csv', delimiter=',', names=True) for k in SPEEDS_KMH}
t = d['time']

# ---- gate check: full vehicle footprint (yawed corners), not just the CG ---------
CORNERS_BODY = np.array([[A + 0.9, -BW/2], [A + 0.9, BW/2],
                         [-(B + 0.9), BW/2], [-(B + 0.9), -BW/2]])


def gate_margins(d):
    """Worst corner clearance to each gate boundary [m]; positive = inside."""
    X, Y, psi = d['X'], d['Y'], np.radians(d['psiDeg'])
    out = []
    for x0, x1, w, c in GATES:
        lo, hi = c - w/2, c + w/2
        over = -9e9
        for cx, cy in CORNERS_BODY:
            xc = X + cx*np.cos(psi) - cy*np.sin(psi)
            yc = Y + cx*np.sin(psi) + cy*np.cos(psi)
            z = (xc >= x0) & (xc <= x1)
            over = max(over, np.maximum(yc[z] - hi, lo - yc[z]).max())
        out.append(-over)
    return out


marg = gate_margins(d)
print(f"ISO 3888-1 gates {'PASS' if min(marg) > 0 else 'FAIL'} "
      f"(footprint margins {'/'.join(f'{1000*m:+.0f}' for m in marg)} mm), "
      f"peak ay {np.abs(d['ayG']).max():.2f} g, "
      f"peak tie-rod {max(np.abs(d['FtieL']).max(), np.abs(d['FtieR']).max()):.0f} N, "
      f"peak handwheel torque {np.abs(d['hwTorque']).max():.1f} N.m")

# ---- 1. path + gates + vehicle outlines ------------------------------------------
# The x axis is compressed SQX:1 for readability. A naive aspect change would shear a
# rotated rectangle, so the outline is drawn as a TRUE rectangle whose yaw is mapped
# into the compressed geometry (psi' = atan(SQX*tan(psi)) keeps it tangent to the
# drawn path); the labels state the true yaw.
SQX = 2
CAR_FRONT, CAR_REAR, CAR_HALFW = A + 0.9, B + 0.9, BW/2   # body frame, CG at origin


def car_poly(x, y, psi):
    """Undistorted vehicle rectangle at compressed CG position, display-mapped yaw."""
    dpsi = np.arctan(SQX*np.tan(psi))
    pts = np.array([[CAR_FRONT, -CAR_HALFW], [CAR_FRONT, CAR_HALFW],
                    [-CAR_REAR, CAR_HALFW], [-CAR_REAR, -CAR_HALFW]])
    c, s = np.cos(dpsi), np.sin(dpsi)
    return pts @ np.array([[c, s], [-s, c]]) + [x/SQX, y]


def draw_course(a):
    xs = np.linspace(0, 125, 500)
    a.plot(xs/SQX, yref_iso(xs), 'k:', lw=1, label='reference centerline')
    for x0, x1, w, c in GATES:
        a.fill_between([x0/SQX, x1/SQX], c - w/2, c + w/2, color='tab:blue', alpha=0.10, lw=0)
        a.plot([x0/SQX, x1/SQX], [c - w/2]*2, 'k-', lw=1.6)
        a.plot([x0/SQX, x1/SQX], [c + w/2]*2, 'k-', lw=1.6)
    ticks = np.arange(0, 121, 20)
    a.set_xticks(ticks/SQX); a.set_xticklabels(ticks)
    a.set_xlabel(f'distance [m]  (axis compressed {SQX}:1)')
    a.set_ylabel('lateral position [m]')
    a.set_aspect('equal'); a.grid(alpha=.3)
    a.set_xlim(-3/SQX, 128/SQX); a.set_ylim(-2.6, 6.2)


FIGSZ_PATH = (16, 3.2)
fig, a = plt.subplots(figsize=FIGSZ_PATH)
draw_course(a)
a.plot(d['X']/SQX, d['Y'], color='tab:red', lw=1.6, label='vehicle CG')
for xst in (20, 70, 100):
    i = np.argmin(np.abs(d['X'] - xst))
    psi = np.radians(d['psiDeg'][i])
    a.add_patch(plt.Polygon(car_poly(d['X'][i], d['Y'][i], psi), closed=True,
                            fc='tab:red', alpha=0.15, ec='tab:red', lw=1.4, zorder=5))
    a.annotate(f"$\\psi$ = {d['psiDeg'][i]:+.1f}°", (d['X'][i]/SQX, d['Y'][i]),
               textcoords='offset points', xytext=(0, -26), ha='center', fontsize=8)
a.set_title('ISO 3888-1 double lane change at 80 km/h — CG path vs cone gates, '
            'vehicle outline at 20 / 70 / 100 m (labels: true yaw)')
a.legend(fontsize=8, loc='upper right')
plt.tight_layout(); plt.savefig(f'{SVG}/dlc_path.svg', facecolor="white"); plt.close()

# ---- 1b. animated trajectory (GIF) ------------------------------------------------
from matplotlib import animation

figA, aA = plt.subplots(figsize=FIGSZ_PATH)
draw_course(aA)
trail, = aA.plot([], [], color='tab:red', lw=1.4)
car = plt.Polygon(car_poly(0, 0, 0), closed=True, fc='tab:red', alpha=0.25,
                  ec='tab:red', lw=1.4, zorder=5)
aA.add_patch(car)
info = aA.text(0.99, 0.93, '', transform=aA.transAxes, ha='right', va='top', fontsize=8)
aA.set_title('ISO 3888-1 double lane change at 80 km/h — manual steering', fontsize=9)
figA.tight_layout()

mask = t <= 6.2
tA, XA, YA = t[mask], d['X'][mask], d['Y'][mask]
psiA, ayA = np.radians(d['psiDeg'][mask]), d['ayG'][mask]
# one frame per 40 ms of simulated time (CSV carries extra event rows)
frames = np.searchsorted(tA, np.arange(0, tA[-1], 0.04))


def _update(i):
    trail.set_data(XA[:i+1]/SQX, YA[:i+1])
    car.set_xy(car_poly(XA[i], YA[i], psiA[i]))
    info.set_text(f't = {tA[i]:.2f} s   $\\psi$ = {np.degrees(psiA[i]):+.1f}°   '
                  f'$a_y$ = {ayA[i]:+.2f} g')
    return trail, car, info


ani = animation.FuncAnimation(figA, _update, frames=frames, blit=True)
os.makedirs(GIF, exist_ok=True)
ani.save(f'{GIF}/dlc_anim.gif', writer=animation.PillowWriter(fps=25), dpi=80)
plt.close(figA)

# ---- 2. vehicle states ------------------------------------------------------------
fig, ax = plt.subplots(3, 1, figsize=(9, 6.4), sharex=True)
ax[0].plot(t, d['ayG'], color='tab:blue', lw=1.3)
ax[0].set_ylabel('lateral accel [g]')
ax[1].plot(t, d['yawRateDegS'], color='tab:purple', lw=1.3)
ax[1].set_ylabel('yaw rate [deg/s]')
ax[2].plot(t, d['dCmdDeg'], 'k--', lw=1, label='driver command')
ax[2].plot(t, d['deltaLdeg'], color='tab:green', lw=1.2, label='road wheel (left)')
ax[2].plot(t, d['alphaFLdeg'], color='tab:orange', lw=1, alpha=0.8, label='slip angle FL')
ax[2].set_ylabel('[deg]'); ax[2].set_xlabel('time [s]'); ax[2].legend(fontsize=8)
for a_ in ax: a_.grid(alpha=.3); a_.axhline(0, color='k', lw=.4)
ax[0].set_title('Double lane change — vehicle response')
plt.tight_layout(); plt.savefig(f'{SVG}/dlc_states.svg', facecolor="white"); plt.close()

# ---- 3. tie-rod forces (headline) --------------------------------------------------
fig, a = plt.subplots(figsize=(9, 4.2))
a.plot(t, d['FtieL']/1000, color='tab:green', lw=1.5, label='tie rod left')
a.plot(t, d['FtieR']/1000, color='tab:red', lw=1.5, label='tie rod right')
a.plot(t, d['rackForce']/1000, color='k', lw=1, alpha=0.5, ls='--', label='rack (sum + wheel inertia)')
for sig, name in [(d['FtieL'], 'left'), (d['FtieR'], 'right')]:
    i = np.argmax(np.abs(sig))
    a.annotate(f"{sig[i]/1000:+.2f} kN", (t[i], sig[i]/1000),
               textcoords='offset points', xytext=(8, 6), fontsize=8)
a.set_xlabel('time [s]'); a.set_ylabel('axial force [kN]')
a.set_title('Tie-rod / steering-arm axial forces — ISO 3888-1 DLC at 80 km/h')
a.grid(alpha=.3); a.axhline(0, color='k', lw=.4); a.legend(fontsize=8)
plt.tight_layout(); plt.savefig(f'{SVG}/dlc_tierod.svg', facecolor="white"); plt.close()

# ---- 4. kingpin axis: road-wheel angle, moment decomposition, trail ----------------
fig, ax = plt.subplots(4, 1, figsize=(9, 10.5), sharex=True)
ax[0].plot(t, d['dCmdDeg'], 'k--', lw=1, label='driver command')
ax[0].plot(t, d['deltaLdeg'], color='tab:green', lw=1.3, label='road wheel left')
ax[0].plot(t, d['deltaRdeg'], color='tab:red', lw=1.3, alpha=0.6, label='road wheel right')
ax[0].set_ylabel('road-wheel angle [deg]'); ax[0].legend(fontsize=8)
for a_, side in zip(ax[1:3], 'LR'):
    a_.plot(t, d[f'Mkp{side}'], color='k', lw=1.5, label='total kingpin moment')
    a_.plot(t, d[f'MkpMech{side}'], color='tab:blue', lw=1.2,
            label=f'mechanical trail (Fy x {TMECH*1000:.0f} mm)')
    a_.plot(t, d[f'MkpPneu{side}'], color='tab:orange', lw=1.2,
            label='pneumatic trail (-Mz, collapses near the grip limit)')
    a_.set_ylabel(f'{"left" if side == "L" else "right"} wheel [N.m]')
ax[1].legend(fontsize=8)
tp0mm = 1000*TIRE_F['ap0']/3
for side, cc in [('L', 'tab:green'), ('R', 'tab:red')]:
    fy = d[f'FyF{side}']
    ok = np.abs(fy) > 500          # trail = -Mz/Fy is ill-conditioned near force reversals
    tp = np.where(ok, -d[f'MzF{side}']/np.where(ok, fy, 1)*1000, np.nan)
    ax[3].plot(t, tp, color=cc, lw=1.3, alpha=0.8 if side == 'R' else 1.0,
               label=f'pneumatic trail {side}')
ax[3].axhline(tp0mm, color='k', ls=':', lw=0.9)
ax[3].text(0.03, tp0mm + 0.8, f'zero-slip trail $a_p/3$ = {tp0mm:.0f} mm at nominal load',
           fontsize=7.5, color='k')
ax[3].set_ylabel('pneumatic trail [mm]'); ax[3].set_ylim(-2, 30)
ax[3].legend(fontsize=8, loc='lower right'); ax[3].set_xlabel('time [s]')
for a_ in ax: a_.grid(alpha=.3); a_.axhline(0, color='k', lw=.4)
ax[0].set_title('Kingpin (steering-axis) view: road-wheel angle, moment split, trail collapse')
plt.tight_layout(); plt.savefig(f'{SVG}/dlc_kingpin_decomp.svg', facecolor="white"); plt.close()

# ---- 5. steering feel: handwheel angle and torque ----------------------------------
fig, ax = plt.subplots(2, 1, figsize=(9, 5.6), sharex=True)
ax[0].plot(t, d['hwaDeg'], color='tab:blue', lw=1.3, label='handwheel angle')
ax[0].plot(t, d['dCmdDeg']*IS, 'k--', lw=0.9, label='driver reference (i_S x command)')
ax[0].set_ylabel('handwheel angle [deg]'); ax[0].legend(fontsize=8)
ax[1].plot(t, d['hwTorque'], color='tab:red', lw=1.4)
i = np.argmax(np.abs(d['hwTorque']))
ax[1].annotate(f"{d['hwTorque'][i]:+.1f} N.m", (t[i], d['hwTorque'][i]),
               textcoords='offset points', xytext=(8, 6), fontsize=8)
ax[1].set_ylabel('handwheel torque [N.m]'); ax[1].set_xlabel('time [s]')
for a_ in ax: a_.grid(alpha=.3); a_.axhline(0, color='k', lw=.4)
ax[0].set_title('Manual (unassisted) steering: handwheel angle and steering-feel torque')
plt.tight_layout(); plt.savefig(f'{SVG}/dlc_handwheel.svg', facecolor="white"); plt.close()

# ---- 6. brush tire curves + Modelica cross-check -----------------------------------
fig, ax = plt.subplots(1, 2, figsize=(9, 3.8))
al = np.radians(np.linspace(-12, 12, 400))
for Fz, cc in [(2000, 'tab:blue'), (4000, 'tab:orange'), (6000, 'tab:green')]:
    Fy, Mz = brush(al, Fz, TIRE_F)
    ax[0].plot(np.degrees(al), Fy/1000, color=cc, lw=1.3, label=f'Fz = {Fz/1000:.0f} kN')
    ax[1].plot(np.degrees(al), Mz, color=cc, lw=1.3)
ss = S[80]; i = ss['time'] > 5.0
a_ss = np.radians(ss['alphaFLdeg'][i].mean())
fz_ss, fy_ss, mz_ss = ss['FzFL'][i].mean(), ss['FyFL'][i].mean(), ss['MzFL'][i].mean()
fy_np, mz_np = brush(a_ss, fz_ss, TIRE_F)
errF, errM = abs(fy_np/fy_ss - 1), abs(mz_np/mz_ss - 1)
assert errF < 1e-3 and errM < 1e-3, f'NumPy/Modelica mismatch: {errF:.2e}, {errM:.2e}'
print(f'brush cross-check (StepSteer 80 km/h steady state): dFy {errF*100:.4f} %, dMz {errM*100:.4f} %')
ax[0].plot(np.degrees(a_ss), fy_ss/1000, 'k*', ms=10, label='Modelica steady state')
ax[1].plot(np.degrees(a_ss), mz_ss, 'k*', ms=10)
ax[0].set_xlabel('slip angle [deg]'); ax[0].set_ylabel('lateral force Fy [kN]')
ax[1].set_xlabel('slip angle [deg]'); ax[1].set_ylabel('aligning moment Mz [N.m]')
ax[0].legend(fontsize=8)
for a_ in ax: a_.grid(alpha=.3); a_.axhline(0, color='k', lw=.4)
fig.suptitle('Brush tire model (front): saturation and aligning-moment collapse', fontsize=10)
plt.tight_layout(); plt.savefig(f'{SVG}/tire_brush.svg', facecolor="white"); plt.close()

# ---- 7. understeer gradient vs analytic bicycle (with aligning moments) ------------
Ca_f = lambda Fz: TIRE_F['c1']*np.sin(2*np.arctan(Fz/TIRE_F['c2']))
Ca_r = lambda Fz: TIRE_R['c1']*np.sin(2*np.arctan(Fz/TIRE_R['c2']))
FzF0 = M*G*B/(2*L); FzR0 = M*G*A/L
Cf, Cr = 2*Ca_f(FzF0), Ca_r(FzR0)
tpf = TIRE_F['ap0']*np.sqrt(FzF0/TIRE_F['FzNom'])/3
tpr = TIRE_R['ap0']*np.sqrt(FzR0/TIRE_R['FzNom'])/3


def bicycle_gain(u):
    """Linear steady-state yaw gain r/delta including aligning-moment trails."""
    Amat = np.array([
        [-(Cf + Cr)/u, -Cf*A/u + Cr*B/u - M*u],
        [(-(A - tpf)*Cf + (B + tpr)*Cr)/u, -((A - tpf)*Cf*A + (B + tpr)*Cr*B)/u]])
    rhs = np.array([-Cf, -(A - tpf)*Cf])
    return np.linalg.solve(Amat, rhs)[1]


fig, a = plt.subplots(figsize=(7, 4))
uu = np.linspace(8, 36, 100)
a.plot(uu*3.6, [bicycle_gain(x) for x in uu], 'k-', lw=1.2,
       label='linear bicycle + aligning moments (analytic)')
gains = []
for kmh in SPEEDS_KMH:
    ss = S[kmh]; i = ss['time'] > 5.0
    gains.append(np.radians(ss['yawRateDegS'][i].mean())/np.radians(ss['deltaLdeg'][i].mean()))
a.plot(SPEEDS_KMH, gains, 'o', color='tab:blue', ms=6, label='tricycle model (0.5 deg step)')
Kus_degg = np.degrees(M/L*(B/Cf - A/Cr)*G)
a.set_xlabel('speed [km/h]'); a.set_ylabel('yaw-rate gain r/delta [1/s]')
a.set_title(f'Steady-state handling check (force-only understeer gradient {Kus_degg:.2f} deg/g)')
a.grid(alpha=.3); a.legend(fontsize=8)
plt.tight_layout(); plt.savefig(f'{SVG}/dlc_understeer.svg', facecolor="white"); plt.close()

# ---- summary CSV -------------------------------------------------------------------
with open(os.path.join(HERE, 'outputs/dlc_summary.csv'), 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['peak_FtieL_N', 'peak_FtieR_N', 'peak_MkpL_Nm', 'peak_MkpR_Nm',
                'peak_rack_N', 'peak_ay_g', 'peak_hwTorque_Nm', 'peak_hwa_deg',
                'gate_margin_entry_mm', 'gate_margin_offset_mm', 'gate_margin_exit_mm'])
    w.writerow([round(np.abs(d[k]).max(), 1) for k in
                ('FtieL', 'FtieR', 'MkpL', 'MkpR', 'rackForce')] +
               [round(np.abs(d['ayG']).max(), 3),
                round(np.abs(d['hwTorque']).max(), 2),
                round(np.abs(d['hwaDeg']).max(), 1)] +
               [round(1000*m, 0) for m in marg])
print('wrote outputs/dlc_summary.csv and 7 SVGs + 1 GIF to outputs/')
