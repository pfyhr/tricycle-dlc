"""Build the live in-browser simulator (proof of concept): outputs/webgui.html.

Ports the planar plant (Track.TrackTricycle, reduced: steer + roll/pitch load-transfer lags
kept, tyre relaxation dropped), the lookahead driver (Track.TrackDriver), and the
quasi-steady speed profile to JavaScript, so a full lap runs in the browser in a few ms.
The racing line (min-curvature geometry) is fixed and car-independent; sliders change the
CAR (mass, power, grip µ, downforce) plus a dry/wet toggle, and the lap re-solves live.
Real telemetry is overlaid on the speed trace. One track (Knutstorp) for the POC.

The UI adopts the Modelica chase-cam player's look (track_render.py): a full-viewport
heading-up follow camera with floating HUD panels. The JS driver is tuned past the Modelica
reference for clean on-line tracking (honest grip budget + mass-aware cornering limit +
tighter lookahead feedback — see GRIP_DERATE/MASS_EXP/KMUL in webgui_template.html); it no
longer matches Modelica lap times to the tenth by design (it drives a tidier, saner line).
"""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'tracks'))
import numpy as np
from speed_profile import load_centerline, speed_profile
from racing_line import min_curvature_line, offset_geometry, _gauss_periodic
from opt_lap import solve_min_time_dyn
from cars import build_config, CARS
from telemetry import load_trace

TRACK = 'knutstorp'
s, x, y, psi, kap = load_centerline(os.path.join(HERE, f'tracks/{TRACK}.csv'))
ds = float(s[1] - s[0]); LTRK = float(s[-1] + ds)
wMax = 10.0/2 - 0.9 - 0.2
psiU = np.unwrap(psi)

# Fixed, car-independent racing line: the 3-DOF minimum-TIME line (CasADi/IPOPT), not the
# geometric min-curvature line. Min-time positions wide on straights to straighten the next
# corner (e.g. sits LEFT onto Knutstorp's main straight, where min-curvature would centre it)
# and picks true late apexes. The JS driver follows it as the car is varied; deltaFF gives it
# the OCP's sideslip steer feedforward so it can hold the aggressive line. Solved once for the
# Elise; the line geometry is close enough across the (fixed-line) car range.
cfg0 = build_config('elise', 'base')
nMC, psiMC, kMC, dsMC = min_curvature_line(x, y, psi, kap, ds, wMax)     # OCP warm start
vMC, _ = speed_profile(s, kMC, ds_seg=dsMC, **cfg0['profile'])
stride = max(int(round(10.0/ds)), int(np.ceil(len(s)/700)))
sc, kc, vN = s[::stride], kap[::stride], vMC[::stride]
wOcp = np.clip(wMax - 0.5 - 1.6*(vN/vN.max())**2, 0.3, wMax)   # inset where fast (run-wide ~v^2)
print(f'solving webgui racing line: 3-DOF min-time OCP ({len(sc)} nodes)...')
nOpt_c, ocp = solve_min_time_dyn(sc, kc, wOcp, p=cfg0['ocp'], grip_frac=cfg0['grip_frac'],
                                 w_reg=4e-2, v_init=vN, n_init=nMC[::stride],
                                 dpsi_init=psiMC[::stride])
print(f'  IPOPT {"converged" if ocp["converged"] else "did NOT converge"}, T={ocp["T"]:.1f}s')
Lc = sc[-1] + (sc[1] - sc[0])
nRef = np.clip(_gauss_periodic(np.interp(s, np.append(sc, Lc), np.append(nOpt_c, nOpt_c[0])),
                               ds, 6.0), -wMax, wMax)
psiRef, kapLine, dsSeg = offset_geometry(x, y, psi, ds, nRef)
# dynamic (sideslip) steer feedforward = OCP steer minus the kinematic part the driver adds
dOcp = np.interp(s, np.append(sc, Lc), np.append(ocp['delta'], ocp['delta'][0]))
uOcp = np.interp(s, np.append(sc, Lc), np.append(ocp['u'], ocp['u'][0]))
d0 = cfg0['driver']
deltaFF = np.clip(_gauss_periodic(dOcp - (d0['Lwb'] + d0['Kus']*uOcp**2)*kapLine, ds, 6.0),
                  -0.08, 0.08)

wrap = lambda a: np.append(a, a[0]).tolist()          # periodic: append s=L point
sW = np.append(s, LTRK).tolist()
# psiU is the UNWRAPPED heading (winds ~2pi over the lap); at s=L continue it smoothly
# (extrapolate) rather than snapping back to psiU[0], else interp swings 2pi across the seam
psiUW = np.append(psiU, 2*psiU[-1] - psiU[-2]).tolist()
TRK = dict(L=LTRK, ds=ds, width=10.0, n=len(s),
           s=sW, x=wrap(x), y=wrap(y), psiU=psiUW, kap_c=wrap(kap),
           nRef=wrap(nRef), psiRef=wrap(psiRef), kapLine=wrap(kapLine), dsSeg=wrap(dsSeg),
           deltaFF=wrap(deltaFF))


def flat(car):
    c = build_config(car, 'base'); t = c['trike']; d = c['driver']
    return dict(name=CARS[car]['display'], m=t['m'], Izz=t['Izz'], a=t['a'], b=t['b'],
                tf=t['tf'], hcg=t['hcg'], xiF=t['xiF'], kBf=t['kBf'], CdA=t['CdA'],
                Crr=t['Crr'], rho=t['rho'], muF=t['muF'], muR=t['muR'], mu=t['muF'],
                c1F=t['c1F'], c1R=t['c1R'], c2F=t['c2F'], c2R=t['c2R'],
                FzNomF=t['FzNomF'], FzNomR=t['FzNomR'], ap0F=t['ap0F'], ap0R=t['ap0R'],
                Pmax=c['Pmax'], ayFrac=c['profile']['ayFrac'], ClA=0.0, aeroBal=0.42,
                Lwb=d['Lwb'], Kus=d['Kus'], KLA=d['KLA'], Kr=d['Kr'])


PRESETS = {c: flat(c) for c in ('elise', 'miata', 'tourer')}
# telemetry per car (fraction, kmh) + real lap seconds from the CSV header
TEL = {}
for car in ('elise', 'miata'):
    tr = load_trace(car, TRACK)
    if tr is not None:
        hdr = open(os.path.join(HERE, 'tracks', 'telemetry', f'{car}_{TRACK}.csv')).readline()
        import re
        m = re.search(r'lap\s+([0-9.]+)s', hdr)
        TEL[car] = dict(frac=tr[0].tolist(), v=tr[1].tolist(),
                        lap=float(m.group(1)) if m else None)

DATA = dict(track=TRK, presets=PRESETS, tel=TEL, trackName='Ring Knutstorp')
html = open(os.path.join(HERE, 'webgui_template.html')).read()
html = html.replace('__DATA__', json.dumps(DATA, separators=(',', ':')))
out = os.path.join(HERE, 'outputs', 'webgui.html')
open(out, 'w').write(html)
print(f'wrote {out} ({len(html)//1024} kB)')
