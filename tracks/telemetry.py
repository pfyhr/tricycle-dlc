"""Real-telemetry speed traces for validation overlays.

The raw VBOX logs live outside the repo; this module extracts one clean flying lap from
each and stores a small speed-vs-lap-fraction CSV under tracks/telemetry/ so the overlay
is self-contained and reproducible without the logs. `load_trace(car, track)` returns the
committed trace for the speed-figure overlay in track_lap.py.

Speed comes straight from the logger; distance is the integral of speed over time from the
start/finish crossing (no GPS projection needed), then normalized to lap fraction so it maps
onto the simulated centerline s-axis (0..L). Corner locations align to within the small
difference between the real driven line and the centerline.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STORE = os.path.join(HERE, 'telemetry')
LOGDIR = os.path.expanduser('~/Downloads/RaceLogs')

# (car, track) -> raw VBOX log used to build the committed trace
SOURCES = {
    ('elise', 'knutstorp'): 'session_knutstorp_20180922_1031.vbo',
    ('miata', 'knutstorp'): 'oskar_knutstorp_lap19_20180512_1039.vbo',
}
# first-block VBOX column indices: velocity(kmh), long, lat, time(hhmmss.ss)
IV, ILON, ILAT, ITIME = 3, 2, 1, 15


def _read(path):
    gate, rows = None, []
    with open(path, errors='replace') as f:
        indata = False
        for ln in f:
            if ln.startswith('Start') and gate is None:
                gate = tuple(float(x) for x in ln.split()[1:5])
            if ln.startswith('[data]'):
                indata = True
                continue
            if indata:
                p = ln.split()
                if len(p) >= 16:
                    try:
                        rows.append((float(p[ITIME]), float(p[ILON]),
                                     float(p[ILAT]), float(p[IV])))
                    except ValueError:
                        pass
    a = np.array(rows)
    T = a[:, 0]
    tsec = np.floor(T/10000)*3600 + np.floor((T % 10000)/100)*60 + (T % 100)
    return gate, tsec, a[:, 1], a[:, 2], a[:, 3]     # gate, t[s], lon, lat, v[kmh]


def _fastest_lap(gate, t, lon, lat, v):
    """Return (t, v) for the fastest flying lap; whole file if no gate laps resolve."""
    gm = ((gate[0] + gate[2])/2, (gate[1] + gate[3])/2)
    d = np.hypot(lon - gm[0], lat - gm[1])
    below = d < 0.02                                  # ~0.02 min ~ 37 m capture radius
    marks, i = [], 0
    while i < len(d):
        if below[i]:
            j = i
            while j < len(d) and below[j]:
                j += 1
            marks.append(i + int(np.argmin(d[i:j])))
            i = j
        else:
            i += 1
    laps = [(marks[k], marks[k+1]) for k in range(len(marks)-1)
            if 55 < t[marks[k+1]] - t[marks[k]] < 85]
    if not laps:
        m = v > 30
        return t[m], v[m]
    i0, i1 = min(laps, key=lambda ab: t[ab[1]] - t[ab[0]])
    return t[i0:i1+1], v[i0:i1+1]


def build():
    """Extract every SOURCES trace and write tracks/telemetry/<car>_<track>.csv."""
    os.makedirs(STORE, exist_ok=True)
    for (car, track), fn in SOURCES.items():
        gate, t, lon, lat, v = _read(os.path.join(LOGDIR, fn))
        tl, vl = _fastest_lap(gate, t, lon, lat, v)
        tl = tl - tl[0]
        dist = np.concatenate([[0], np.cumsum(0.5*(vl[1:] + vl[:-1])/3.6 * np.diff(tl))])
        frac = dist/dist[-1]
        out = os.path.join(STORE, f'{car}_{track}.csv')
        hdr = f'lap_fraction,speed_kmh   # {fn}, lap {tl[-1]:.1f}s, {dist[-1]:.0f} m'
        np.savetxt(out, np.column_stack([frac, vl]), fmt='%.5f', delimiter=',', header=hdr)
        print(f'{car:6s} {track:10s} lap {int(tl[-1]//60)}:{tl[-1] % 60:04.1f}  '
              f'{dist[-1]:.0f} m, vmax {vl.max():.0f} km/h -> {os.path.relpath(out, HERE)}')


def load_trace(car, track):
    """Committed real speed trace: (lap_fraction[0..1], speed_kmh) or None if absent."""
    p = os.path.join(STORE, f'{car}_{track}.csv')
    if not os.path.exists(p):
        return None
    d = np.loadtxt(p, delimiter=',')
    return d[:, 0], d[:, 1]


if __name__ == '__main__':
    build()
