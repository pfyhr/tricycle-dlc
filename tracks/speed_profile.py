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
                  mu=0.95, ayFrac=0.90, hcg=0.55, a=1.20, b=1.60, nPass=4,
                  ds_seg=None):
    """Return (vRef, axFF) on the s grid (closed track). `s` and `kappa` are indexed
    on the centerline grid; `ds_seg[i]` is the arc length of the driven path over the
    i->i+1 step (pass the racing line's per-segment length here, else the uniform
    centerline ds is used)."""
    N = len(s)
    dsc = s[1] - s[0]
    dseg = np.full(N, dsc) if ds_seg is None else np.asarray(ds_seg)
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
            v[jn] = min(v[jn], np.sqrt(max(v[j]**2 + 2*dseg[j]*ax, 1.0)))
        for i in range(N - 1, -1, -1):  # backward: grip-limited braking (drag helps)
            j, jn = i, (i + 1) % N
            axB = ayFrac*mu*G*ellipse(v[jn], kappa[jn]) \
                + (0.5*rho*CdA*v[jn]**2 + Crr*m*G)/m
            v[j] = min(v[j], np.sqrt(v[jn]**2 + 2*dseg[j]*axB))
    # acceleration feedforward along the driven path (central difference over 2 segments)
    dsFF = dseg + np.roll(dseg, 1)
    axFF = (np.roll(v, -1)**2 - np.roll(v, 1)**2)/(2*dsFF)
    return v, axFF


def write_track_table(path, s, kappa, vRef, axFF, nRef=None, psiRef=None,
                      kappaLine=None, deltaFF=None, extend=1500.0):
    """Write the OpenModelica text table 'track' with 8 columns
    (s, kappa_c, vRef, axFF, n_ref, psi_ref, kappa_line, delta_ff), repeating the first
    `extend` meters beyond s = L so preview lookups never leave the table.

    n_ref / psi_ref / kappa_line describe the racing line the driver follows (offset,
    heading vs centerline, line curvature). Omit them for centerline following:
    n_ref = psi_ref = 0 and kappa_line = kappa_c, which the driver reduces to exactly.
    delta_ff is an optional dynamic steer feedforward (the OCP line's sideslip-steer
    beyond the kinematic term); 0 for the geometric lines (the driver's own kinematic
    feedforward then carries it)."""
    ds = s[1] - s[0]
    Ltrk = s[-1] + ds
    nExt = int(extend/ds) + 1
    z = np.zeros_like(kappa)
    nRef = z if nRef is None else nRef
    psiRef = z if psiRef is None else psiRef
    kappaLine = kappa if kappaLine is None else kappaLine
    deltaFF = z if deltaFF is None else deltaFF
    sAll = np.concatenate([s, Ltrk + s[:nExt]])
    wrap = lambda f: np.concatenate([f, f[:nExt]])
    tab = np.column_stack([sAll, wrap(kappa), wrap(vRef), wrap(axFF),
                           wrap(nRef), wrap(psiRef), wrap(kappaLine), wrap(deltaFF)])
    with open(path, 'w') as f:
        f.write('#1\n')
        f.write(f'double track({tab.shape[0]},8)  '
                f'# s, kappa_c, vRef, axFF, n_ref, psi_ref, kappa_line, delta_ff\n')
        np.savetxt(f, tab, fmt='%.6f')
    return Ltrk
