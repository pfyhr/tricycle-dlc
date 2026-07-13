"""RaceChrono VBO session parser: every lap from every log, with GPS + IMU channels.

telemetry.py extracts ONE speed trace per (car, track) for the committed overlay CSVs.
This module is the raw-data layer under the model-tuning work: it parses whole sessions
(23- and 25-column RaceChrono VBO variants, 10 Hz GPS with logger long/lat accel),
splits them into laps at the session's start/finish gate, and returns per-lap traces
parameterized by driven distance so they can be compared - and fitted - against the
simulator's speed profile on the same 0..1 lap-fraction axis.

Raw logs live outside the repo (~/Downloads/RaceLogs); nothing here is imported by the
build. Use catalog() for a session/lap-time survey and laps_of() for the trace data.
"""
import os
import glob
import numpy as np

LOGDIR = os.path.expanduser('~/Downloads/RaceLogs')


def read_vbo(path):
    """Parse one RaceChrono VBO -> dict of channel arrays + the start/finish gate.

    Columns (first block, both variants): satellites, lat, lon, v_kmh, heading,
    height, long_g, lat_g, update_rate, lean, combined_g, fix, ... ; second block
    repeats satellites then time(hhmmss.ss). Uses lat/lon/v/accel from block 1 and
    time from block 2 (index 15), same layout telemetry.py relies on.
    """
    import re
    gate, rows, cols = None, [], None
    tpat = re.compile(r'^\d{5,6}\.\d{1,3}$')                 # hhmmss.cc wall-clock field
    with open(path, errors='replace') as f:
        indata = False
        for ln in f:
            if ln.startswith('Start') and gate is None:
                p = ln.split()
                gate = tuple(float(x) for x in p[1:5])       # lon1 lat1 lon2 lat2
            if ln.startswith('[data]'):
                indata = True
                continue
            if indata:
                p = ln.split()
                if len(p) < 9:
                    continue
                if cols is None:                              # layouts differ between logger
                    it = next((k for k in range(1, len(p)) if tpat.match(p[k])), None)
                    if it is None:
                        continue
                    # standard VBOX order (RaceChrono v7+): sat time lat lon v hdg alt ax ay
                    # old two-block export: sat lat lon v hdg alt ax ay ... sat time ...
                    cols = ((it, it+1, it+2, it+3, it+6, it+7) if it == 1 else
                            (it, 1, 2, 3, 6, 7))              # (t, lat, lon, v, ax, ay)
                try:
                    rows.append((float(p[cols[0]]), float(p[cols[2]]), float(p[cols[1]]),
                                 float(p[cols[3]]), float(p[cols[4]]), float(p[cols[5]])))
                except (ValueError, IndexError):
                    pass
    if not rows:
        return None
    a = np.array(rows)
    T = a[:, 0]
    t = np.floor(T/10000)*3600 + np.floor((T % 10000)/100)*60 + (T % 100)
    t = np.unwrap(t, period=86400)                            # sessions may cross midnight
    return dict(gate=gate, t=t, lon=a[:, 1], lat=a[:, 2], v=a[:, 3],
                ax=a[:, 4], ay=a[:, 5])


def split_laps(d, tmin=30, tmax=900):
    """Gate crossings -> list of (i0, i1) sample index ranges, one per flying lap."""
    g = d['gate']
    gm = ((g[0] + g[2])/2, (g[1] + g[3])/2)
    r = np.hypot(d['lon'] - gm[0], d['lat'] - gm[1])
    below = r < 0.02                                          # ~37 m capture radius
    marks, i = [], 0
    while i < len(r):
        if below[i]:
            j = i
            while j < len(r) and below[j]:
                j += 1
            marks.append(i + int(np.argmin(r[i:j])))
            i = j
        else:
            i += 1
    t = d['t']
    return [(marks[k], marks[k+1]) for k in range(len(marks)-1)
            if tmin < t[marks[k+1]] - t[marks[k]] < tmax]


def lap_trace(d, i0, i1):
    """One lap -> dict(t0-based t, frac 0..1 by driven distance, v kmh, ax/ay g,
    lat/lon [NMEA minutes], length m)."""
    t = d['t'][i0:i1+1] - d['t'][i0]
    v = d['v'][i0:i1+1]
    dist = np.concatenate([[0], np.cumsum(0.5*(v[1:] + v[:-1])/3.6 * np.diff(t))])
    return dict(t=t, frac=dist/max(dist[-1], 1e-9), dist=dist, v=v,
                ax=d['ax'][i0:i1+1], ay=d['ay'][i0:i1+1],
                lat=d['lat'][i0:i1+1], lon=d['lon'][i0:i1+1],
                time=float(t[-1]), length=float(dist[-1]))


def laps_of(path, tmin=30, tmax=900, lref=None, ltol=0.05):
    """All flying laps of one session, fastest first. With lref, keep only laps whose
    driven length is within ltol of it (rejects gate double-captures = half laps)."""
    d = read_vbo(path)
    if d is None or d['gate'] is None:
        return []
    laps = [lap_trace(d, i0, i1) for i0, i1 in split_laps(d, tmin, tmax)]
    if lref is None and laps:                                 # default: the session's own
        lref = float(np.median([L['length'] for L in laps]))  # median lap length
    laps = [L for L in laps if abs(L['length'] - lref) <= ltol*lref]
    return sorted(laps, key=lambda L: L['time'])


def catalog(pattern='*.vbo'):
    """Survey every session: lap count, best/median time, lap length, peak |ay|."""
    out = []
    for p in sorted(glob.glob(os.path.join(LOGDIR, pattern))):
        laps = laps_of(p)
        if not laps:
            out.append((os.path.basename(p), 0, None, None, None, None))
            continue
        times = [L['time'] for L in laps]
        best = laps[0]
        out.append((os.path.basename(p), len(laps), min(times),
                    float(np.median(times)), best['length'],
                    float(np.abs(best['ay']).max())))
    return out


if __name__ == '__main__':
    fmt = lambda s: '   -  ' if s is None else f'{int(s//60)}:{s % 60:04.1f}'
    print(f'{"session":52s} laps  best   median  len_m  pk_ay')
    for name, n, best, med, ln, ay in catalog():
        print(f'{name:52s} {n:4d}  {fmt(best)}  {fmt(med)}  '
              f'{"" if ln is None else f"{ln:5.0f}"}  {"" if ay is None else f"{ay:.2f}"}')
