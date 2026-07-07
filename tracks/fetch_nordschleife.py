"""Fetch the Nurburgring Nordschleife centerline from OpenStreetMap and write a
planar, arc-length-parameterized track file: tracks/nordschleife.csv with columns
s [m], x [m], y [m], psi [rad, unwrapped], kappa [1/m].

The OSM route relation 38566 ("Nurburgring Nordschleife", route=road) holds the
ordered raceway ways of the ~20.8 km loop. Nodes are stitched into one closed
polyline, projected onto a local tangent plane (planar: elevation is dropped by
design), resampled to uniform arc length, and low-pass filtered with a periodic
Gaussian before computing heading and curvature (raw polyline curvature is
noise).

Data (c) OpenStreetMap contributors, ODbL - see sources/SOURCES.md.

Usage:  python3 tracks/fetch_nordschleife.py            # writes tracks/nordschleife.csv
The raw Overpass response is cached in tracks/.nordschleife_osm.json so the
script only needs network on first run.
"""
import json, os, sys, urllib.request

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, '.nordschleife_osm.json')
OUT = os.path.join(HERE, 'nordschleife.csv')

RELATION = 38566  # OSM route relation "Nurburgring Nordschleife"
DS = 2.0          # internal resample step [m]
DS_OUT = 5.0      # output step [m]
SIGMA = 6.0       # Gaussian smoothing of the centerline [m]
R_EARTH = 6371008.8


def fetch_raw():
    if os.path.exists(RAW):
        return json.load(open(RAW))
    query = f'[out:json][timeout:120];relation({RELATION});out geom;'
    req = urllib.request.Request(
        'https://overpass-api.de/api/interpreter',
        data=('data=' + urllib.parse.quote(query)).encode(),
        headers={'User-Agent': 'tricycle-dlc-tracksim/0.1'})
    with urllib.request.urlopen(req, timeout=120) as r:
        data = json.load(r)
    json.dump(data, open(RAW, 'w'))
    return data


def stitch(ways):
    """Order/orient the member ways into one closed node loop."""
    segs = [[(p['lat'], p['lon']) for p in w['geometry']] for w in ways]
    loop = segs.pop(0)
    while segs:
        for i, s in enumerate(segs):
            if s[0] == loop[-1]:
                loop += s[1:]; break
            if s[-1] == loop[-1]:
                loop += s[::-1][1:]; break
        else:
            raise RuntimeError(f'no way continues from {loop[-1]}; {len(segs)} left')
        segs.pop(i)
    if loop[0] != loop[-1]:
        raise RuntimeError('stitched polyline is not closed')
    return np.array(loop[:-1])  # drop duplicate closing node


def gauss_periodic(f, ds, sigma):
    n = int(np.ceil(4*sigma/ds))
    k = np.exp(-0.5*(np.arange(-n, n + 1)*ds/sigma)**2)
    k /= k.sum()
    return np.convolve(np.concatenate([f[-n:], f, f[:n]]), k, mode='same')[n:-n]


def main():
    data = fetch_raw()
    rel = next(e for e in data['elements'] if e['type'] == 'relation')
    ways = [m for m in rel['members'] if m['type'] == 'way' and 'geometry' in m]
    latlon = stitch(ways)
    print(f'{len(ways)} ways, {len(latlon)} nodes')

    # local tangent plane at the centroid; x east, y north (planar by design)
    lat0, lon0 = np.radians(latlon.mean(axis=0))
    lat, lon = np.radians(latlon).T
    x = R_EARTH*np.cos(lat0)*(lon - lon0)
    y = R_EARTH*(lat - lat0)

    # ensure the lap runs in OSM way direction starting at the relation's first node,
    # resample raw polyline to uniform ds
    dx, dy = np.diff(np.append(x, x[0])), np.diff(np.append(y, y[0]))
    sRaw = np.concatenate([[0], np.cumsum(np.hypot(dx, dy))])
    Lraw = sRaw[-1]
    sU = np.arange(0, Lraw, DS)
    xU = np.interp(sU, sRaw, np.append(x, x[0]))
    yU = np.interp(sU, sRaw, np.append(y, y[0]))

    # periodic Gaussian smoothing, then heading/curvature by central differences
    xS = gauss_periodic(xU, DS, SIGMA)
    yS = gauss_periodic(yU, DS, SIGMA)
    xp = (np.roll(xS, -1) - np.roll(xS, 1))/(2*DS)
    yp = (np.roll(yS, -1) - np.roll(yS, 1))/(2*DS)
    xpp = (np.roll(xS, -1) - 2*xS + np.roll(xS, 1))/DS**2
    ypp = (np.roll(yS, -1) - 2*yS + np.roll(yS, 1))/DS**2
    kap = (xp*ypp - yp*xpp)/np.maximum((xp**2 + yp**2)**1.5, 1e-12)
    psi = np.unwrap(np.arctan2(yp, xp))

    # re-parameterize by the smoothed arc length and thin to the output step
    dsS = np.hypot(np.roll(xS, -1) - xS, np.roll(yS, -1) - yS)
    sS = np.concatenate([[0], np.cumsum(dsS)])[:-1]
    L = sS[-1] + dsS[-1]
    sOut = np.arange(0, L, DS_OUT)
    cols = [np.interp(sOut, sS, f) for f in (xS, yS, psi, kap)]

    hdr = (f'# Nurburgring Nordschleife planar centerline, L = {L:.1f} m\n'
           f'# from OSM relation {RELATION} (route=road); (c) OpenStreetMap contributors, ODbL 1.0\n'
           f'# smoothed with a periodic Gaussian, sigma = {SIGMA:.0f} m; elevation dropped\n'
           's_m,x_m,y_m,psi_rad,kappa_1pm')
    np.savetxt(OUT, np.column_stack([sOut] + cols), delimiter=',',
               header=hdr, comments='', fmt='%.6f')
    rmin = 1/np.abs(kap).max()
    print(f'L = {L:.1f} m, min radius = {rmin:.1f} m, wrote {OUT} ({len(sOut)} rows)')
    if not 20000 < L < 22000:
        sys.exit('unexpected track length - check the OSM relation members')


if __name__ == '__main__':
    main()
