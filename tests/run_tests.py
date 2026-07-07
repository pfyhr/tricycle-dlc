"""Test suite for the Tricycle package: checkModel on every class, simulate every
example, and assert the physics against analytic references. Exits nonzero on failure.

Run locally or in CI:  python3 tests/run_tests.py
Requires: omc (OpenModelica) on PATH with the Modelica Standard Library installed,
Python 3 with numpy.
"""
import numpy as np, subprocess, shutil, os, sys, tempfile

OMC = shutil.which('omc') or '/Users/pontus/opt/openmodelica/bin/omc'
# pin OPENMODELICAHOME to the omc install so generated code always finds its headers
# (openmodelica.h / omc_simulation_settings.h), even via symlinks or relocated installs
os.environ.setdefault('OPENMODELICAHOME',
                      os.path.dirname(os.path.dirname(os.path.realpath(OMC))))
HERE = os.path.dirname(os.path.abspath(__file__))
MOD = os.path.join(os.path.dirname(HERE), 'modelica')
TMP = tempfile.mkdtemp(prefix='tricycle_tests_')
os.makedirs(os.path.join(os.path.dirname(HERE), 'modelica', 'build'), exist_ok=True)

CLASSES = ['Tricycle.PlanarTricycle', 'Tricycle.ManualSteering', 'Tricycle.Iso3888Path',
           'Tricycle.Examples.StepSteer', 'Tricycle.Examples.DoubleLaneChange',
           'Tricycle.Examples.OpenLoopDLC', 'Tricycle.Track.TrackTricycle',
           'Tricycle.Track.TrackDriver', 'Tricycle.Examples.NordschleifeLap']

# ---- parameters mirrored from Tricycle.mo defaults --------------------------------
M, A, B = 1650.0, 1.20, 1.60
L = A + B; G = 9.80665
BW = 1.8
TIRE_F = dict(c1=7.0e4, c2=4000.0, mu=0.95, ap0=0.06, FzNom=4000.0)
TIRE_R = dict(c1=1.4e5, c2=8000.0, mu=0.95, ap0=0.085, FzNom=8000.0)
GATES = [(0, 15, 1.1*BW + 0.25, 0.0), (45, 70, 1.2*BW + 0.25, 3.5),
         (95, 125, 1.3*BW + 0.25, 0.0)]
CORNERS = np.array([[A + 0.9, -BW/2], [A + 0.9, BW/2],
                    [-(B + 0.9), BW/2], [-(B + 0.9), -BW/2]])

failures = []


def check(name, cond, detail=''):
    status = 'ok  ' if cond else 'FAIL'
    print(f'[{status}] {name}' + (f'  ({detail})' if detail else ''))
    if not cond:
        failures.append(name)


def run_mos(script):
    path = os.path.join(TMP, 'script.mos')
    open(path, 'w').write(script)
    r = subprocess.run([OMC, path], cwd=MOD, capture_output=True, text=True)
    return r.stdout + r.stderr


def brush(alpha, Fz, d):
    """NumPy twin of Tricycle.brushForces."""
    Fzc = np.maximum(Fz, 50.0)
    Ca = d['c1']*np.sin(2*np.arctan(Fzc/d['c2']))
    ap = d['ap0']*np.sqrt(Fzc/d['FzNom'])
    sy = np.tan(alpha)
    sAbs = np.sqrt(sy**2 + 1e-12)
    sgn = sy/sAbs
    th = Ca/(3*d['mu']*Fzc)
    uu = 0.5*(th*sAbs + 1 - np.sqrt((th*sAbs - 1)**2 + 1e-4))
    lam = 1 - uu
    return d['mu']*Fzc*(1 - lam**3)*sgn, -d['mu']*Fzc*ap*uu*lam**3*sgn


def bicycle_gain(u):
    """Linear steady-state yaw gain r/delta including aligning-moment trails."""
    Ca_f = lambda Fz: TIRE_F['c1']*np.sin(2*np.arctan(Fz/TIRE_F['c2']))
    Ca_r = lambda Fz: TIRE_R['c1']*np.sin(2*np.arctan(Fz/TIRE_R['c2']))
    FzF0 = M*G*B/(2*L); FzR0 = M*G*A/L
    Cf, Cr = 2*Ca_f(FzF0), Ca_r(FzR0)
    tpf = TIRE_F['ap0']*np.sqrt(FzF0/TIRE_F['FzNom'])/3
    tpr = TIRE_R['ap0']*np.sqrt(FzR0/TIRE_R['FzNom'])/3
    Amat = np.array([
        [-(Cf + Cr)/u, -Cf*A/u + Cr*B/u - M*u],
        [(-(A - tpf)*Cf + (B + tpr)*Cr)/u, -((A - tpf)*Cf*A + (B + tpr)*Cr*B)/u]])
    rhs = np.array([-Cf, -(A - tpf)*Cf])
    return np.linalg.solve(Amat, rhs)[1]


def gate_margins(X, Y, psiDeg):
    psi = np.radians(psiDeg)
    out = []
    for x0, x1, w, c in GATES:
        lo, hi = c - w/2, c + w/2
        over = -9e9
        for cx, cy in CORNERS:
            xc = X + cx*np.cos(psi) - cy*np.sin(psi)
            yc = Y + cx*np.sin(psi) + cy*np.cos(psi)
            z = (xc >= x0) & (xc <= x1)
            if z.any():
                over = max(over, np.maximum(yc[z] - hi, lo - yc[z]).max())
        out.append(-over)
    return out


# ---- 1. checkModel on every class -------------------------------------------------
script = 'loadModel(Modelica); getErrorString();\nloadFile("Tricycle.mo"); getErrorString();\n'
script += '\n'.join(f'checkModel({c}); getErrorString();' for c in CLASSES)
out = run_mos(script)
# omc prints one line per statement: true / "" / true / "" ... for the two loads
check('package loads', out.splitlines()[:4].count('true') == 2, 'loadModel + loadFile')
for c in CLASSES:
    check(f'checkModel {c}', f'Check of {c} completed successfully' in out)

# ---- 2. StepSteer: understeer vs analytic bicycle + brush cross-check --------------
rf = os.path.join(TMP, 'step.csv')
out = run_mos(
    'loadModel(Modelica); loadFile("Tricycle.mo");\n'
    f'simulate(Tricycle.Examples.StepSteer, stopTime=6, numberOfIntervals=3000, '
    f'outputFormat="csv", variableFilter="time|yawRateDegS|deltaLdeg|alphaFLdeg|FzFL|FyFL|MzFL", '
    f'fileNamePrefix="build/step", simflags="-r={rf}"); getErrorString();')
check('StepSteer simulates', os.path.exists(rf), rf if not os.path.exists(rf) else '')
if os.path.exists(rf):
    s = np.genfromtxt(rf, delimiter=',', names=True)
    i = s['time'] > 5.0
    gain = np.radians(s['yawRateDegS'][i].mean())/np.radians(s['deltaLdeg'][i].mean())
    ref = bicycle_gain(80/3.6)
    check('StepSteer yaw gain vs analytic bicycle+Mz', abs(gain/ref - 1) < 0.03,
          f'sim {gain:.3f} vs analytic {ref:.3f} 1/s')
    a_ss = np.radians(s['alphaFLdeg'][i].mean())
    fy_np, mz_np = brush(a_ss, s['FzFL'][i].mean(), TIRE_F)
    check('brush tire NumPy twin matches Modelica',
          abs(fy_np/s['FyFL'][i].mean() - 1) < 1e-3 and abs(mz_np/s['MzFL'][i].mean() - 1) < 1e-3,
          f'dFy {abs(fy_np/s["FyFL"][i].mean()-1):.2e}')

# ---- 3. DoubleLaneChange: ISO gates (footprint), stability, sane steering effort ---
rf = os.path.join(TMP, 'dlc.csv')
out = run_mos(
    'loadModel(Modelica); loadFile("Tricycle.mo");\n'
    f'simulate(Tricycle.Examples.DoubleLaneChange, stopTime=8, numberOfIntervals=8000, '
    f'outputFormat="csv", variableFilter="time|X|Y|psiDeg|ayG|hwTorque|FtieL|FtieR", '
    f'fileNamePrefix="build/dlc", simflags="-r={rf}"); getErrorString();')
check('DoubleLaneChange simulates', os.path.exists(rf))
if os.path.exists(rf):
    d = np.genfromtxt(rf, delimiter=',', names=True)
    marg = gate_margins(d['X'], d['Y'], d['psiDeg'])
    check('DLC passes ISO 3888-1 gates (footprint)', min(marg) > 0,
          '/'.join(f'{1000*m:+.0f} mm' for m in marg))
    check('DLC settles (no spin-out)', abs(d['Y'][-1]) < 0.05, f'Y end {d["Y"][-1]:+.3f} m')
    ay = np.abs(d['ayG']).max()
    check('DLC peak lateral acceleration in severe-maneuver band', 0.5 < ay < 0.95, f'{ay:.2f} g')
    hwt = np.abs(d['hwTorque']).max()
    check('DLC handwheel torque plausible for unassisted rack', 4 < hwt < 20, f'{hwt:.1f} N.m')
    ftie = max(np.abs(d['FtieL']).max(), np.abs(d['FtieR']).max())
    check('DLC tie-rod force in expected range', 500 < ftie < 4000, f'{ftie:.0f} N')

# ---- 4. OpenLoopDLC: runs, bounded, lane-change-like ------------------------------
rf = os.path.join(TMP, 'ol.csv')
out = run_mos(
    'loadModel(Modelica); loadFile("Tricycle.mo");\n'
    f'simulate(Tricycle.Examples.OpenLoopDLC, stopTime=5, numberOfIntervals=5000, '
    f'outputFormat="csv", variableFilter="time|Y|ayG|deltaLdeg", '
    f'fileNamePrefix="build/ol", simflags="-r={rf}"); getErrorString();')
check('OpenLoopDLC simulates', os.path.exists(rf))
if os.path.exists(rf):
    o = np.genfromtxt(rf, delimiter=',', names=True)
    check('OpenLoopDLC bounded (no spin-out)',
          np.abs(o['Y']).max() < 8 and np.abs(o['ayG']).max() < 1.0,
          f'Ypk {np.abs(o["Y"]).max():.2f} m, ay {np.abs(o["ayG"]).max():.2f} g')
    check('OpenLoopDLC steer tracks command (arm filter attenuation < 15%)',
          2.5 < np.abs(o['deltaLdeg']).max() <= 3.0,
          f'{np.abs(o["deltaLdeg"]).max():.2f} of 3 deg')

# ---- 5. NordschleifeLap: full lap on the planar OSM track --------------------------
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'tracks'))
from speed_profile import load_centerline, speed_profile, write_track_table

sT, xT, yT, psiT, kapT = load_centerline(
    os.path.join(os.path.dirname(HERE), 'tracks', 'nordschleife.csv'))
vRef, axFF = speed_profile(sT, kapT)
LTRK = write_track_table(os.path.join(MOD, 'build', 'ns_track.txt'),
                         sT, kapT, vRef, axFF)
check('Nordschleife track data sane', 20000 < LTRK < 22000 and
      abs(np.abs(kapT).max() - 1/19.5) < 0.02, f'L {LTRK:.0f} m')
tIdeal = np.sum((sT[1] - sT[0])/vRef)

rf = os.path.join(TMP, 'nslap.csv')
out = run_mos(
    'loadModel(Modelica); loadFile("Tricycle.mo");\n'
    f'simulate(Tricycle.Examples.NordschleifeLap, stopTime=900, numberOfIntervals=9000, '
    f'outputFormat="csv", variableFilter="time|s|n|vKmh|vRefKmh|ayG|axG|FtieL|FtieR", '
    f'fileNamePrefix="build/nslap_test", simflags="-r={rf}"); getErrorString();')
check('NordschleifeLap simulates', os.path.exists(rf))
if os.path.exists(rf):
    nl = np.genfromtxt(rf, delimiter=',', names=True)
    tl, m = nl['time'][-1], nl['time'] > 15
    check('lap completes', nl['s'][-1] >= LTRK - 10,
          f"s {nl['s'][-1]:.0f} of {LTRK:.0f} m")
    check('lap time near the quasi-steady ideal (driver overhead < 5%)',
          tIdeal < tl < 1.05*tIdeal, f'{tl:.1f} s vs ideal {tIdeal:.1f} s')
    check('stays on the road (|n| < 2.5 m)', np.abs(nl['n']).max() < 2.5,
          f"|n|max {np.abs(nl['n']).max():.2f} m")
    check('speed tracks the minimum-time profile', np.sqrt(
          ((nl['vKmh'] - nl['vRefKmh'])[m]**2).mean()) < 15,
          f"rms {np.sqrt(((nl['vKmh'] - nl['vRefKmh'])[m]**2).mean()):.1f} km/h")
    check('lap uses the grip (peak |ay| in 0.8-1.0 g)',
          0.8 < np.abs(nl['ayG']).max() < 1.0, f"{np.abs(nl['ayG']).max():.2f} g")
    check('lap tie-rod force in expected range', 500 < max(
          np.abs(nl['FtieL']).max(), np.abs(nl['FtieR']).max()) < 4000,
          f"{max(np.abs(nl['FtieL']).max(), np.abs(nl['FtieR']).max()):.0f} N")

# ---- summary -----------------------------------------------------------------------
print()
if failures:
    print(f'{len(failures)} test(s) FAILED: ' + ', '.join(failures))
    sys.exit(1)
print('all tests passed')
