"""Quasi-steady minimum-time speed profile for a closed track, and the OpenModelica
table file the Tricycle.Track models read.

Classic three-step lap-time-simulation profile (e.g. Brayshaw & Harrison 2005):
  1. corner-speed limit  v = sqrt(ayCap/|kappa|), capped at the drag-limited top speed
  2. forward pass:  acceleration limited by engine power, driven-axle traction, and
     the friction-ellipse remainder of the previewed lateral acceleration
  3. backward pass: braking limited the same way (drag helps)
Passes run around the loop until converged (periodic track). The result is what a
driver can manage on the centerline - not a formal minimum-time optimum.
"""
import numpy as np

G = 9.80665


def load_centerline(path):
    """Read tracks/nordschleife.csv (or same format): s,x,y,psi,kappa."""
    d = np.genfromtxt(path, delimiter=',', skip_header=3, names=True)
    return (d['s_m'], d['x_m'], d['y_m'], d['psi_rad'], d['kappa_1pm'])


def top_speed(m, Pmax, CdA, Crr, rho):
    v = 50.0
    for _ in range(60):
        f = 0.5*rho*CdA*v**3 + Crr*m*G*v - Pmax
        v -= f/(1.5*rho*CdA*v**2 + Crr*m*G)
    return v


def speed_profile(s, kappa, m=1650.0, Pmax=150e3, CdA=0.72, Crr=0.012, rho=1.20,
                  mu=0.95, ayFrac=0.90, hcg=0.55, a=1.20, b=1.60, nPass=4):
    """Return (vRef, axFF) on the same s grid (closed track, uniform spacing)."""
    N = len(s)
    ds = s[1] - s[0]
    L = a + b
    ayCap = ayFrac*mu*G
    vTop = top_speed(m, Pmax, CdA, Crr, rho)
    # steady RWD traction limit including longitudinal load transfer
    axTrac = mu*G*(a/L)/(1 - mu*hcg/L)
    v = np.minimum(vTop, np.sqrt(ayCap/np.maximum(np.abs(kappa), 1e-5)))

    def ellipse(vj, kj):
        return np.sqrt(max(0.0, 1 - min((vj**2*abs(kj)/ayCap)**2, 1.0)))

    for _ in range(nPass):
        for i in range(N):            # forward: power/traction-limited acceleration
            j, jn = i, (i + 1) % N
            axP = (Pmax/max(v[j], 5) - 0.5*rho*CdA*v[j]**2 - Crr*m*G)/m
            ax = min(axP, axTrac)*ellipse(v[j], kappa[j])
            v[jn] = min(v[jn], np.sqrt(max(v[j]**2 + 2*ds*ax, 1.0)))
        for i in range(N - 1, -1, -1):  # backward: grip-limited braking (drag helps)
            j, jn = i, (i + 1) % N
            axB = ayFrac*mu*G*ellipse(v[jn], kappa[jn]) \
                + (0.5*rho*CdA*v[jn]**2 + Crr*m*G)/m
            v[j] = min(v[j], np.sqrt(v[jn]**2 + 2*ds*axB))
    axFF = (np.roll(v, -1)**2 - np.roll(v, 1)**2)/(4*ds)
    return v, axFF


def write_track_table(path, s, kappa, vRef, axFF, extend=1500.0):
    """Write the OpenModelica text table 'track' (s, kappa, vRef, axFF), with the
    first `extend` meters repeated beyond s = L so preview lookups never leave the
    table."""
    ds = s[1] - s[0]
    Ltrk = s[-1] + ds
    nExt = int(extend/ds) + 1
    sAll = np.concatenate([s, Ltrk + s[:nExt]])
    wrap = lambda f: np.concatenate([f, f[:nExt]])
    tab = np.column_stack([sAll, wrap(kappa), wrap(vRef), wrap(axFF)])
    with open(path, 'w') as f:
        f.write('#1\n')
        f.write(f'double track({tab.shape[0]},4)  '
                f'# s [m], kappa [1/m], vRef [m/s], axFF [m/s2]\n')
        np.savetxt(f, tab, fmt='%.6f')
    return Ltrk
