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
from scipy.sparse.linalg import spsolve


def _gauss_periodic(f, ds, sigma):
    n = int(np.ceil(4*sigma/ds))
    k = np.exp(-0.5*(np.arange(-n, n + 1)*ds/sigma)**2)
    k /= k.sum()
    return np.convolve(np.concatenate([f[-n:], f, f[:n]]), k, mode='same')[n:-n]


def min_curvature_line(x, y, psi, kappa_c, ds, w_max, alpha=3e-6, smooth=10.0):
    """Return (n_ref, psi_ref, kappa_line, ds_seg) on the centerline grid.

    Minimum-curvature offset via the regularized normal equations
        (D2^T D2 + alpha I) w = -D2^T kappa_c ,
    a periodic pentadiagonal system solved directly (fast, and smooth by
    construction - no bang-bang against the corridor). The offset is then scaled
    uniformly to fit within +/- w_max (scaling preserves smoothness), so the line
    stays inside the track without introducing curvature kinks. `alpha` trades
    curvature reduction against how far the raw solution reaches (a larger value
    keeps it gentler); the scaling makes the exact value non-critical.

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

    H = (D2.T @ D2 + alpha*sparse.eye(N, format='csr')).tocsc()
    w = spsolve(H, -(D2.T @ kappa_c))
    w = _gauss_periodic(w, ds, smooth) if smooth else w
    peak = np.abs(w).max()
    if peak > w_max:
        w *= w_max/peak

    # exact geometry of the offset path (derivatives wrt centerline arc length)
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
    return w, dpsi_ref, kappa_line, ds_seg
