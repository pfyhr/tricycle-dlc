"""Fetch a race-track centerline from OpenStreetMap and write a planar,
arc-length-parameterized track file tracks/<key>.csv with columns
s [m], x [m], y [m], psi [rad, unwrapped], kappa [1/m].

Tracks are declared in TRACKS below - either an OSM route relation (ordered member
ways, e.g. the Nordschleife) or a bounding-box query whose raceway ways are assembled
into the unique closed circuit loop (pit lanes, karting tracks, and closed infield
loops excluded). The polyline is projected onto a local tangent plane (planar:
elevation is dropped by design), resampled, low-pass filtered with a periodic
Gaussian, and rotated so s = 0 sits mid-straight (except when the raw start is kept).

Data (c) OpenStreetMap contributors, ODbL - see sources/SOURCES.md.

Usage:  python3 tracks/fetch_track.py --track=knutstorp   (or --track=all)
Raw Overpass responses are cached in tracks/.<key>_osm.json, so re-processing
needs no network.
"""
import argparse, json, os, sys, urllib.parse, urllib.request

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DS = 2.0          # internal resample step [m]
DS_OUT = 5.0      # output step [m]
SIGMA = 6.0       # Gaussian smoothing of the centerline [m]
R_EARTH = 6371008.8

TRACKS = {
    'nordschleife': dict(
        display='Nürburgring Nordschleife', prefix='ns',
        query='relation(38566);out geom;', relation=True,
        lmin=20000, lmax=22000, rotate_start=False),  # committed CSV predates rotation
    'knutstorp': dict(
        display='Ring Knutstorp', prefix='knutstorp',
        query='way(229220546);out geom;',
        lmin=1900, lmax=2200),
    'anderstorp': dict(
        display='Anderstorp Raceway', prefix='anderstorp',
        query='way["highway"="raceway"](57.23,13.55,57.30,13.65);out geom;',
        seed='Flygrakan', lmin=3900, lmax=4150),
    'gelleras': dict(
        display='Gelleråsen Arena', prefix='gelleras',
        query='way["highway"="raceway"](59.33,14.45,59.40,14.55);out geom;',
        lmin=2250, lmax=2450),
    'kinnekulle': dict(
        display='Kinnekulle Ring', prefix='kinnekulle',
        query='way["highway"="raceway"](58.52,13.33,58.60,13.45);out geom;',
        seed='Kinnekulle Ring', lmin=1950, lmax=2200),
}


def fetch_raw(key, query):
    cache = os.path.join(HERE, f'.{key}_osm.json')
    if os.path.exists(cache):
        return json.load(open(cache))
    req = urllib.request.Request(
        'https://overpass-api.de/api/interpreter',
        data=('data=' + urllib.parse.quote(f'[out:json][timeout:120];{query}')).encode(),
        headers={'User-Agent': 'tricycle-dlc-tracksim/0.1'})
    with urllib.request.urlopen(req, timeout=150) as r:
        data = json.load(r)
    json.dump(data, open(cache, 'w'))
    return data


def stitch_relation(data):
    """Order/orient a route relation's member ways into one closed node loop."""
    rel = next(e for e in data['elements'] if e['type'] == 'relation')
    segs = [[(p['lat'], p['lon']) for p in m['geometry']]
            for m in rel['members'] if m['type'] == 'way' and 'geometry' in m]
    loop = segs.pop(0)
    while segs:
        for i, sgm in enumerate(segs):
            if sgm[0] == loop[-1]:
                loop += sgm[1:]; break
            if sgm[-1] == loop[-1]:
                loop += sgm[::-1][1:]; break
        else:
            raise RuntimeError(f'no way continues from {loop[-1]}; {len(segs)} left')
        segs.pop(i)
    if loop[0] != loop[-1]:
        raise RuntimeError('stitched polyline is not closed')
    return loop[:-1], []


def assemble_loop(data, seed_name, lmin, lmax):
    """Find the unique closed circuit loop among a bbox's raceway ways."""
    def rough_len0(pts, closed=True):
        lat0 = np.radians(np.mean([p[0] for p in pts]))
        x = R_EARTH*np.cos(lat0)*np.radians([p[1] for p in pts])
        y = R_EARTH*np.radians([p[0] for p in pts])
        if closed:
            x, y = np.append(x, x[0]), np.append(y, y[0])
        return np.hypot(np.diff(x), np.diff(y)).sum()

    ways, oneway = {}, {}
    for w in data['elements']:
        t = w.get('tags', {})
        if 'geometry' not in w or t.get('sport') == 'karting' or t.get('name') == 'Pit Lane':
            continue
        g = [(p['lat'], p['lon']) for p in w['geometry']]
        if g[0] == g[-1]:
            if lmin < rough_len0(g[:-1]) < lmax:  # a single closed way IS the circuit
                return g[:-1], [w['id']]
            continue  # other closed ways are infield side loops, not circuit parts
        ways[w['id']] = dict(g=g, name=t.get('name'))
        oneway[w['id']] = t.get('oneway') == 'yes'
    seed = ([i for i, w in ways.items() if seed_name and w['name'] == seed_name] or
            [max(ways, key=lambda i: rough_len0(ways[i]['g'], closed=False))])[0]
    ends = {i: (w['g'][0], w['g'][-1]) for i, w in ways.items()}
    found = []

    def dfs(chain, node, used):
        if len(chain) > 30:
            return
        if node == ends[seed][0]:
            found.append(list(chain)); return
        for j, (a, b) in ends.items():
            if j in used:
                continue
            if a == node:
                dfs(chain + [(j, False)], b, used | {j})
            elif b == node:
                dfs(chain + [(j, True)], a, used | {j})
    dfs([(seed, False)], ends[seed][1], {seed})
    if not found:
        raise RuntimeError('no closed loop found')
    loops = []
    for chain in found:
        pts = []
        for j, rev in chain:
            g = ways[j]['g'][::-1] if rev else ways[j]['g']
            pts += g[:-1]
        loops.append((chain, pts))
    # keep loops in the expected length band; they must agree
    def rough_len(pts):
        lat0 = np.radians(np.mean([p[0] for p in pts]))
        x = R_EARTH*np.cos(lat0)*np.radians([p[1] for p in pts])
        y = R_EARTH*np.radians([p[0] for p in pts])
        return np.hypot(np.diff(np.append(x, x[0])), np.diff(np.append(y, y[0]))).sum()
    loops = [(c, p) for c, p in loops if lmin < rough_len(p) < lmax]
    if not loops:
        raise RuntimeError('no loop in the expected length band')
    chain, pts = loops[0]
    # orient by the mapped one-way (racing) direction where tagged
    fwd = sum(1 for j, rev in chain if oneway[j] and not rev)
    bwd = sum(1 for j, rev in chain if oneway[j] and rev)
    if bwd > fwd:
        pts = pts[::-1]
    return pts, [j for j, _ in chain]


def gauss_periodic(f, ds, sigma):
    n = int(np.ceil(4*sigma/ds))
    k = np.exp(-0.5*(np.arange(-n, n + 1)*ds/sigma)**2)
    k /= k.sum()
    return np.convolve(np.concatenate([f[-n:], f, f[:n]]), k, mode='same')[n:-n]


def process(key, cfg):
    data = fetch_raw(key, cfg['query'])
    if cfg.get('relation'):
        latlon, way_ids = stitch_relation(data)
    else:
        latlon, way_ids = assemble_loop(data, cfg.get('seed'), cfg['lmin'], cfg['lmax'])
    latlon = np.array(latlon)
    print(f"{key}: {len(latlon)} nodes" +
          (f" from ways {way_ids}" if way_ids else ' (relation)'))

    lat0, lon0 = np.radians(latlon.mean(axis=0))
    lat, lon = np.radians(latlon).T
    x = R_EARTH*np.cos(lat0)*(lon - lon0)
    y = R_EARTH*(lat - lat0)

    dx, dy = np.diff(np.append(x, x[0])), np.diff(np.append(y, y[0]))
    sRaw = np.concatenate([[0], np.cumsum(np.hypot(dx, dy))])
    sU = np.arange(0, sRaw[-1], DS)
    xU = np.interp(sU, sRaw, np.append(x, x[0]))
    yU = np.interp(sU, sRaw, np.append(y, y[0]))

    xS = gauss_periodic(xU, DS, SIGMA)
    yS = gauss_periodic(yU, DS, SIGMA)
    xp = (np.roll(xS, -1) - np.roll(xS, 1))/(2*DS)
    yp = (np.roll(yS, -1) - np.roll(yS, 1))/(2*DS)
    xpp = (np.roll(xS, -1) - 2*xS + np.roll(xS, 1))/DS**2
    ypp = (np.roll(yS, -1) - 2*yS + np.roll(yS, 1))/DS**2
    kap = (xp*ypp - yp*xpp)/np.maximum((xp**2 + yp**2)**1.5, 1e-12)

    # start the lap at the ENTRY of the longest straight (corner exit, like a real
    # start/finish line): a flying lap then begins accelerating on settled states
    # instead of braking into the first corner with cold ones
    if cfg.get('rotate_start', True):
        w = int(120/DS)
        score = np.convolve(np.abs(np.concatenate([kap, kap[:w]])),
                            np.ones(w)/w, 'valid')[:len(kap)]
        thr = 1.3*score.min() + 1e-5
        iMin = int(np.argmin(score))
        i0 = iMin
        while score[(i0 - 1) % len(kap)] < thr:   # walk back to the straight's entry
            i0 = (i0 - 1) % len(kap)
            if i0 == iMin:
                break
        i0 = (i0 + int(30/DS)) % len(kap)         # small margin past the corner exit
        xS, yS, kap = np.roll(xS, -i0), np.roll(yS, -i0), np.roll(kap, -i0)

    xp = (np.roll(xS, -1) - np.roll(xS, 1))/(2*DS)
    yp = (np.roll(yS, -1) - np.roll(yS, 1))/(2*DS)
    psi = np.unwrap(np.arctan2(yp, xp))

    dsS = np.hypot(np.roll(xS, -1) - xS, np.roll(yS, -1) - yS)
    sS = np.concatenate([[0], np.cumsum(dsS)])[:-1]
    L = sS[-1] + dsS[-1]
    sOut = np.arange(0, L, DS_OUT)
    cols = [np.interp(sOut, sS, f) for f in (xS, yS, psi, kap)]

    out = os.path.join(HERE, f'{key}.csv')
    hdr = (f'# {cfg["display"]} planar centerline, L = {L:.1f} m\n'
           f'# from OSM ({cfg["query"].split(";")[0]}); (c) OpenStreetMap contributors, ODbL 1.0\n'
           f'# smoothed with a periodic Gaussian, sigma = {SIGMA:.0f} m; elevation dropped\n'
           's_m,x_m,y_m,psi_rad,kappa_1pm')
    np.savetxt(out, np.column_stack([sOut] + cols), delimiter=',',
               header=hdr, comments='', fmt='%.6f')
    turn = (psi[-1] - psi[0])/(2*np.pi)
    print(f'  L = {L:.1f} m, min radius = {1/np.abs(kap).max():.1f} m, '
          f'net turn = {turn:+.2f} rev ({"CW" if turn < 0 else "CCW"}), wrote {out}')
    if not cfg['lmin'] < L < cfg['lmax']:
        sys.exit(f'unexpected length {L:.0f} m for {key}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--track', default='all', help='track key or "all": '
                    + ', '.join(TRACKS))
    args = ap.parse_args()
    keys = list(TRACKS) if args.track == 'all' else [args.track]
    for k in keys:
        if k not in TRACKS:
            sys.exit(f'unknown track {k!r}; known: {", ".join(TRACKS)}')
        process(k, TRACKS[k])
