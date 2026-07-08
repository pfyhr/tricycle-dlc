"""Build the live in-browser simulator (proof of concept): outputs/webgui.html.

Ports the planar plant (Track.TrackTricycle, reduced: steer + roll/pitch load-transfer lags
kept, tyre relaxation dropped), the lookahead driver (Track.TrackDriver), and the
quasi-steady speed profile to JavaScript, so a full lap runs in the browser in a few ms.
The racing line (min-curvature geometry) is fixed and car-independent; sliders change the
CAR (mass, power, grip, downforce, brake bias, balance) and the lap re-solves live. Real
telemetry is overlaid on the speed trace. One track (Knutstorp) for the POC.
"""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'tracks'))
import numpy as np
from speed_profile import load_centerline, speed_profile
from racing_line import min_curvature_line, apply_driver_margin
from cars import build_config, CARS
from telemetry import load_trace

TRACK = 'knutstorp'
s, x, y, psi, kap = load_centerline(os.path.join(HERE, f'tracks/{TRACK}.csv'))
ds = float(s[1] - s[0]); LTRK = float(s[-1] + ds)
wMax = 10.0/2 - 0.9 - 0.2
psiU = np.unwrap(psi)

# fixed, car-independent racing line (min-curvature, inset for the driver overshoot using
# the Elise-base speed as a reference); the JS driver follows it as the car is varied
cfg0 = build_config('elise', 'base')
nRef, psiRef, kapLine, dsSeg = min_curvature_line(x, y, psi, kap, ds, wMax)
vMC, _ = speed_profile(s, kapLine, ds_seg=dsSeg, **cfg0['profile'])
nRef, psiRef, kapLine, dsSeg = apply_driver_margin(x, y, psi, ds, nRef, wMax, vMC, margin=0.4)

wrap = lambda a: np.append(a, a[0]).tolist()          # periodic: append s=L point
sW = np.append(s, LTRK).tolist()
TRK = dict(L=LTRK, ds=ds, width=10.0, n=len(s),
           s=sW, x=wrap(x), y=wrap(y), psiU=wrap(psiU), kap_c=wrap(kap),
           nRef=wrap(nRef), psiRef=wrap(psiRef), kapLine=wrap(kapLine), dsSeg=wrap(dsSeg))


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
