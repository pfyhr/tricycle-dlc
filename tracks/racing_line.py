"""Minimum-curvature racing line for a closed planar track.

Given the centerline (x, y, heading, curvature on a uniform arc-length grid) and the
usable half-width, solve for the lateral offset profile w(s) that minimizes the
path curvature the car has to take, subject to staying within the track:

    minimize   integral (kappa_c + w'')^2 ds     (linearized offset-path curvature)
    over w(s)  periodic,  |w(s)| <= w_max

This is the classic minimum-curvature line (Braghin et al. 2008; Heilmeier et al.,
Veh. Syst. Dyn. 2020) - a convex bounded least-squares problem. It is not the true
minimum-TIME line (that also trades a little curvature for a shorter path and depends
on the power/grip budget), but for a fixed corridor it is very close and needs no
nonlinear optimal-control machinery.

The linearized objective picks the offsets; the returned curvature, heading offset and
per-segment length are then recomputed from the EXACT geometry of the resulting offset
path, so the speed profile and the driver feedforward see the real line, not the
approximation.

Sign convention matches the model: n / w positive to the LEFT of the centerline
(ISO 8855 y-left), offset point = (x,y) + w*(-sin psi, cos psi).
"""
import numpy as np
from scipy import sparse
from scipy.optimize import lsq_linear


def _gauss_periodic(f, ds, sigma):
    n = int(np.ceil(4*sigma/ds))
    k = np.exp(-0.5*(np.arange(-n, n + 1)*ds/sigma)**2)
    k /= k.sum()
    return np.convolve(np.concatenate([f[-n:], f, f[:n]]), k, mode='same')[n:-n]


def min_curvature_line(x, y, psi, kappa_c, ds, w_max, alpha=5e-6, smooth=3.0):
    """Return (n_ref, psi_ref, kappa_line, ds_seg) on the centerline grid.

    Minimum-curvature offset as the BOX-CONSTRAINED least-squares problem

        minimize  || D2 w + kappa_c ||^2 + alpha || w ||^2   s.t.  |w| <= w_max

    where (D2 w)[i] ~ w''(s_i), so D2 w + kappa_c is the (linearized) curvature of
    the offset line. Solving with the corridor as a hard box constraint lets the line
    ride the inside edge through an apex and swing wide on entry/exit - a real racing
    line. This replaces an earlier formulation that solved unconstrained and then
    scaled the WHOLE line down to make its single widest point fit the corridor, which
    squashed every other corner toward the centre (a timid, centre-hugging line). The
    `alpha` ridge pins the otherwise free constant lateral offset (D2's null space) and
    keeps the system well posed; because the line's useful large offsets are LOW-curvature
    (cheap under the D2 objective), alpha also trades how eagerly the line rides the
    corridor edge (near 0 it hugs the edge bang-bang; ~5e-6 relaxes it toward the more
    selective apexing of the true minimum-TIME line - keep it small or the line goes
    timid). `smooth` gently rounds the junctions where the line meets the corridor edge
    so the steer feedforward stays continuous.

    n_ref     lateral offset of the racing line [m], |n_ref| <= w_max
    psi_ref   heading of the line relative to the centerline [rad]
    kappa_line curvature of the line [1/m] (exact, from the offset geometry)
    ds_seg    arc length of each line segment i->i+1 [m] (for the speed profile)
    """
    N = len(x)
    # periodic second-difference operator D2, so (D2 w)[i] ~ w''(s_i)
    off = np.ones(N)
    D2 = sparse.diags([off, -2*off, off], [-1, 0, 1], (N, N), format='lil')
    D2[0, N - 1] = 1.0
    D2[N - 1, 0] = 1.0
    D2 = (D2.tocsr())/ds**2

    # stack the curvature objective with a sqrt(alpha) centering ridge, box-bounded solve
    A = sparse.vstack([D2, np.sqrt(alpha)*sparse.eye(N, format='csr')]).tocsr()
    b = np.concatenate([-kappa_c, np.zeros(N)])
    w = lsq_linear(A, b, bounds=(-w_max, w_max), max_iter=300, tol=1e-10).x
    if smooth:
        w = np.clip(_gauss_periodic(w, ds, smooth), -w_max, w_max)
    dpsi_ref, kappa_line, ds_seg = offset_geometry(x, y, psi, ds, w)
    return w, dpsi_ref, kappa_line, ds_seg


def apply_driver_margin(x, y, psi, ds, w, w_max, vRef, margin=0.4, k=1.4,
                        ay_budget=None, smooth=4.0, asym=True, bite=0.0, kHi=0.0):
    """Pull a racing line in from the corridor edge where the driver needs slack.

    The preview driver overshoots the line under LATERAL LOAD - through corners,
    worst when fast - but tracks a straight essentially exactly, so margin taken on
    straights (e.g. the flat-out run to a turn-in point) is pure width wasted. With
    ay_budget the inset scales with the line's own lateral demand |kappa|*v^2 /
    ay_budget, smoothed ~15 m so it eases in around turn-in; without it, the legacy
    speed-only (v/vmax)^2 scaling. Returns the reshaped (w, dpsi_ref, kappa_line,
    ds_seg) from the exact offset geometry of the inset line."""
    if ay_budget is not None:
        # overshoot needs speed AND load: every grip-limited corner rides the full
        # lateral budget, so the load term alone saturates in slow hairpins too -
        # gate it by speed (full above ~137 km/h) so slow corners keep their width
        _, kap0, _ = offset_geometry(x, y, psi, ds, w)
        dem = (np.clip(np.abs(kap0)*vRef*vRef/ay_budget, 0.0, 1.0)
               * np.clip((vRef/38.0)**2, 0.0, 1.0))
        dem = _gauss_periodic(dem, ds, 15.0)
        # hiV: how far past the driver's full-authority speed this point is (the steering
        # gain knee sits at ~150 km/h). With kHi the margin gets a high-speed surcharge and
        # the apex freedom / bite fade out - lets a long track (Nordschleife) run the full
        # aggressive treatment through its slow-medium corners while its 180+ km/h sections
        # keep the protection that validation showed they need.
        hiV = _gauss_periodic(np.clip((vRef - 38.0)/10.0, 0.0, 1.0)**2, ds, 20.0)
        wCap = np.clip(w_max - margin - (k + kHi*hiV)*dem, 0.3, w_max)
        # DIRECTION-AWARE: overshoot only ever pushes the car OUTWARD - it cannot
        # overshoot into an apex (it arrives wide, never deep). Margin therefore
        # belongs on the outside of a corner only; the apex side keeps full width.
        if asym:
            aX = _gauss_periodic(np.clip(kap0/0.004, -1.0, 1.0), ds, 8.0)
            fade = 1.0 - hiV if kHi else np.ones_like(hiV)   # apex freedom fades where fast
            wIn = wCap + (w_max - margin - wCap)*fade
            hi = wCap + (wIn - wCap)*np.maximum(aX, 0.0)     # +n bound: apex side of a LEFT
            lo = -(wCap + (wIn - wCap)*np.maximum(-aX, 0.0)) # -n bound: apex side of a RIGHT
            w = np.clip(w, lo, hi)
            if bite:
                # aim past the apex to hit it: the saturated driver converges onto a deep
                # apex ~0.5 m short no matter the gains (measured; feedback is tanh-limited
                # at the grip limit) - so bias the REFERENCE deeper at the dip bottoms by
                # the known undershoot, and the driven car lands on the kerb
                wB = w_max - margin
                apx = (np.maximum(aX, 0.0)*np.clip((w - 0.6*wB)/(0.4*wB), 0, 1)
                       - np.maximum(-aX, 0.0)*np.clip((-w - 0.6*wB)/(0.4*wB), 0, 1))
                w = w + bite*apx*fade
            if kHi:   # speed-blended smoothing: sharp apex dips survive at slow corners,
                      # fast flowing sections stay calm (heavy sigma = the weave guard)
                w = (1.0 - hiV)*_gauss_periodic(w, ds, 4.0) + hiV*_gauss_periodic(w, ds, smooth)
            else:
                w = _gauss_periodic(w, ds, smooth)
        else:
            w = _gauss_periodic(np.clip(w, -wCap, wCap), ds, smooth)
    else:
        wCap = np.clip(w_max - margin - k*(vRef/vRef.max())**2, 0.3, w_max)
        w = _gauss_periodic(np.clip(w, -wCap, wCap), ds, smooth)
    dpsi_ref, kappa_line, ds_seg = offset_geometry(x, y, psi, ds, w)
    return w, dpsi_ref, kappa_line, ds_seg


def offset_geometry(x, y, psi, ds, w):
    """Exact heading (vs centerline), curvature and per-segment length of the path
    offset laterally by w(s) from the centerline. Shared by the minimum-curvature
    line and the optimal-control line so both feed the driver the same way."""
    nx, ny = -np.sin(psi), np.cos(psi)
    xl = x + w*nx
    yl = y + w*ny
    xp = (np.roll(xl, -1) - np.roll(xl, 1))/(2*ds)
    yp = (np.roll(yl, -1) - np.roll(yl, 1))/(2*ds)
    xpp = (np.roll(xl, -1) - 2*xl + np.roll(xl, 1))/ds**2
    ypp = (np.roll(yl, -1) - 2*yl + np.roll(yl, 1))/ds**2
    kappa_line = (xp*ypp - yp*xpp)/np.maximum((xp**2 + yp**2)**1.5, 1e-12)
    psi_line = np.unwrap(np.arctan2(yp, xp))
    dpsi_ref = np.arctan2(np.sin(psi_line - psi), np.cos(psi_line - psi))
    ds_seg = np.hypot(np.roll(xl, -1) - xl, np.roll(yl, -1) - yl)
    return dpsi_ref, kappa_line, ds_seg
