"""Build the live in-browser simulator: outputs/index.html (the site landing page).

Ports the planar plant (Track.TrackTricycle, reduced: steer + roll/pitch load-transfer lags
kept, tyre relaxation dropped), the lookahead driver (Track.TrackDriver), and the
quasi-steady speed profile to JavaScript, so a full lap runs in the browser in a few ms.

For each sim track the fixed, car-independent racing line is the 3-DOF minimum-TIME line
(CasADi/IPOPT, solved once for the Elise), not the geometric min-curvature line: min-time
positions wide on straights to straighten the next corner (e.g. sits LEFT onto Knutstorp's
main straight) and picks true late apexes. Its sideslip steer feedforward (deltaFF) is baked
too so the JS driver holds the aggressive line. Sliders change the CAR (mass, power, grip mu,
downforce) plus a dry/wet toggle; a dropdown switches track. The lap re-solves live. Real
telemetry (where logged) is overlaid on the speed trace.

The UI adopts the Modelica chase-cam player's look (track_render.py): a full-viewport
heading-up follow camera with floating HUD panels. The JS driver reflexes are tuned to a
capable club driver (see GRIP_DERATE/MASS_EXP/KMUL + C.TpFF/tauSteer in webgui_template.html).
"""
import json, os, sys, re, pickle
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'tracks'))
import numpy as np
from speed_profile import load_centerline, speed_profile
from racing_line import min_curvature_line, offset_geometry, _gauss_periodic, apply_driver_margin
from opt_lap import solve_min_time_dyn
from cars import build_config, CARS
from fetch_track import TRACKS
from telemetry import load_trace

# cache the (slow) per-track OCP solves so line post-processing can be re-tuned cheaply;
# invalidated by bumping CACHE_VER (delete .webgui_ocp_cache.pkl to force a full re-solve)
CACHE_VER = 1
CACHE = os.path.join(HERE, '.webgui_ocp_cache.pkl')
_cache = {}
if os.path.exists(CACHE):
    with open(CACHE, 'rb') as f:
        _cache = pickle.load(f)

# sim tracks (Knutstorp first = default); short tracks then the Nordschleife
TLIST = ['knutstorp', 'anderstorp', 'gelleras', 'kinnekulle', 'nordschleife']
CAR_HALF, EDGE = 0.9, 0.05
cfg0 = build_config('elise', 'base')          # racing line solved for the Elise
d0 = cfg0['driver']


def rnd(a, n):
    return [round(float(v), n) for v in a]


def build_track(key):
    """Load a centreline, solve the min-time racing line, return the JS track dict."""
    s, x, y, psi, kap = load_centerline(os.path.join(HERE, f'tracks/{key}.csv'))
    ds = float(s[1] - s[0]); L = float(s[-1] + ds)
    width = TRACKS[key].get('width', 10.0)
    wMax = width/2 - CAR_HALF - EDGE
    psiU = np.unwrap(psi)
    # min-curvature warm start, then the 3-DOF minimum-time OCP line
    nMC, psiMC, kMC, dsMC = min_curvature_line(x, y, psi, kap, ds, wMax)
    vMC, _ = speed_profile(s, kMC, ds_seg=dsMC, **cfg0['profile'])
    stride = max(int(round(10.0/ds)), int(np.ceil(len(s)/700)))
    sc, kc, vN = s[::stride], kap[::stride], vMC[::stride]
    # corridor narrows with LATERAL DEMAND, not speed: a flat-out straight (kappa~0) keeps
    # full width - entry positioning before a turn-in point should hug the edge - while
    # fast kinks (where the driver's overshoot lives) stay pulled in
    ayb = 0.9*cfg0['profile']['ayFrac']*cfg0['profile']['mu']*9.80665
    demN = (np.clip(np.abs(kc)*vN**2/ayb, 0.0, 1.0)     # loaded...
            * np.clip((vN/38.0)**2, 0.0, 1.0))          # ...AND fast (slow corners keep width)
    wOcp = np.clip(wMax - 0.2 - 0.8*demN, 0.3, wMax)
    o = _cache.get(key)                             # OCP inputs fixed per track; rm cache to re-solve
    if o is not None and 'n' in o:
        nOpt_c, dOcp_c, uOcp_c = o['n'], o['d'], o['u']
        print(f'  {key}: OCP line (cached, T={o["T"]:.1f}s)', flush=True)
    else:
        print(f'  {key}: OCP min-time line ({len(sc)} nodes)...', flush=True)
        nOpt_c, ocp = solve_min_time_dyn(sc, kc, wOcp, p=cfg0['ocp'], grip_frac=cfg0['grip_frac'],
                                         w_reg=4e-2, v_init=vN, n_init=nMC[::stride],
                                         dpsi_init=psiMC[::stride])
        print(f'    IPOPT {"converged" if ocp["converged"] else "did NOT converge"}, '
              f'T={ocp["T"]:.1f}s', flush=True)
        nOpt_c, dOcp_c, uOcp_c = np.array(nOpt_c), np.array(ocp['delta']), np.array(ocp['u'])
        _cache[key] = dict(n=nOpt_c, d=dOcp_c, u=uOcp_c, T=ocp['T'])
        with open(CACHE, 'wb') as f:
            pickle.dump(_cache, f)
    Lc = sc[-1] + (sc[1] - sc[0])
    nRef = np.clip(_gauss_periodic(np.interp(s, np.append(sc, Lc), np.append(nOpt_c, nOpt_c[0])),
                                   ds, 6.0), -wMax, wMax)
    psiRef, kapLine, dsSeg = offset_geometry(x, y, psi, ds, nRef)
    # light driver margin: the JS lookahead driver overshoots the aggressive OCP line by up to
    # ~0.6 m on the grippier cars; cap the line a touch (more where fast) so the tyres stay on
    vLine, _ = speed_profile(s, kapLine, ds_seg=dsSeg, **cfg0['profile'])
    # per-track driver margin: the short Swedish tracks are wide and slow enough that the
    # driver's apex overshoot lands on the kerb; the Nordschleife is narrow with 200 km/h
    # corners (where the steering-gain knee limits tracking) and needs more line in hand
    mg, mk = (0.4, 3.4) if key == 'nordschleife' else (0.2, 2.6)
    nRef, psiRef, kapLine, dsSeg = apply_driver_margin(x, y, psi, ds, nRef, wMax, vLine,
                                                       margin=mg, k=mk, ay_budget=ayb)
    dOcp = np.interp(s, np.append(sc, Lc), np.append(dOcp_c, dOcp_c[0]))
    uOcp = np.interp(s, np.append(sc, Lc), np.append(uOcp_c, uOcp_c[0]))
    deltaFF = np.clip(_gauss_periodic(dOcp - (d0['Lwb'] + d0['Kus']*uOcp**2)*kapLine, ds, 6.0),
                      -0.08, 0.08)
    wrap = lambda a: np.append(a, a[0])
    return dict(
        L=round(L, 2), ds=round(ds, 3), width=width, n=len(s),
        s=rnd(np.append(s, L), 2), x=rnd(wrap(x), 2), y=rnd(wrap(y), 2),
        # psiU winds ~2pi/lap: continue it (extrapolate) at s=L, not snap to psiU[0]
        psiU=rnd(np.append(psiU, 2*psiU[-1] - psiU[-2]), 4), kap_c=rnd(wrap(kap), 5),
        nRef=rnd(wrap(nRef), 3), psiRef=rnd(wrap(psiRef), 4), kapLine=rnd(wrap(kapLine), 5),
        dsSeg=rnd(wrap(dsSeg), 3), deltaFF=rnd(wrap(deltaFF), 4))


def flat(car):
    c = build_config(car, 'base'); t = c['trike']; d = c['driver']
    return dict(name=CARS[car]['display'], m=t['m'], Izz=t['Izz'], a=t['a'], b=t['b'],
                tf=t['tf'], hcg=t['hcg'], xiF=t['xiF'], kBf=t['kBf'], CdA=t['CdA'],
                Crr=t['Crr'], rho=t['rho'], muF=t['muF'], muR=t['muR'], mu=t['muF'],
                c1F=t['c1F'], c1R=t['c1R'], c2F=t['c2F'], c2R=t['c2R'],
                FzNomF=t['FzNomF'], FzNomR=t['FzNomR'], ap0F=t['ap0F'], ap0R=t['ap0R'],
                Pmax=c['Pmax'], ayFrac=c['profile']['ayFrac'],
                ClA=t.get('ClA', 0.0), aeroBal=t.get('aeroBal', 0.42),
                Lwb=d['Lwb'], Kus=d['Kus'], KLA=d['KLA'], Kr=d['Kr'])


print('solving racing lines:')
TRK = {t: build_track(t) for t in TLIST}
NAMES = {t: TRACKS[t]['display'] for t in TLIST}
PRESETS = {c: flat(c) for c in ('elise', 'miata', 'tourer', 'm140', 'clubman')}

# real telemetry per (track, car): fraction, kmh, real lap seconds from the CSV header
TEL = {}
for t in TLIST:
    for car in ('elise', 'miata'):
        tr = load_trace(car, t)
        if tr is None:
            continue
        hdr = open(os.path.join(HERE, 'tracks', 'telemetry', f'{car}_{t}.csv')).readline()
        m = re.search(r'lap\s+([0-9.]+)s', hdr)
        TEL.setdefault(t, {})[car] = dict(frac=rnd(tr[0], 4), v=rnd(tr[1], 1),
                                          lap=float(m.group(1)) if m else None)

DATA = dict(tracks=TRK, order=TLIST, names=NAMES, default='knutstorp',
            presets=PRESETS, tel=TEL)
html = open(os.path.join(HERE, 'webgui_template.html')).read()
html = html.replace('__DATA__', json.dumps(DATA, separators=(',', ':')))
out = os.path.join(HERE, 'outputs', 'index.html')
open(out, 'w').write(html)
print(f'wrote {out} ({len(html)//1024} kB, {len(TLIST)} tracks)')
