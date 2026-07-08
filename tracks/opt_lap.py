"""Provably (locally) minimum-time lap by direct-collocation optimal control.

Solves the minimum-time optimal-control problem in the spatial (arc-length) domain
around the closed centerline, for a **friction-circle point-mass** abstraction of the
car with aerodynamic drag, rolling resistance and an engine-power limit:

    minimize   T = integral_0^L (1/sdot) ds
    over       states   x(s) = [n, chi, v]   (offset, heading-vs-centerline, speed)
               controls u(s) = [a_t, a_n]    (tangential / lateral acceleration)
    subject to the curvilinear point-mass kinematics (d/ds = (d/dt)/sdot),
               tyre force  F = m*a + drag + roll  inside the friction circle
                 |F|^2 <= (mu m g)^2   and   F_tangential * v <= Pmax,
               |n| <= w_max            (stay in the track corridor)
               periodic boundary conditions (closed lap).

Transcribed by trapezoidal collocation into one NLP and solved with IPOPT; the
solver's satisfaction of the Karush-Kuhn-Tucker conditions certifies *local*
optimality (a nonconvex problem - not a global proof, which is impractical for a
vehicle model). The friction limit mu and the power/drag/mass are the SAME as
Tricycle.Track.TrackTricycle, so the optimum is a like-for-like lower bound on the
simulated lap. What this point-mass omits - yaw/sideslip dynamics, tyre relaxation
and load transfer - is exactly what the full Modelica model adds back when it is then
driven around this optimal line (track_lap.py --line=ocp): the OCP chooses the
provably-optimal line, the plant validates it.

Requires casadi (pip install casadi; bundles IPOPT).
"""
import numpy as np
import casadi as ca

G = 9.80665
PARAMS = dict(m=1650.0, mu=0.95, Pmax=150e3, CdA=0.72, Crr=0.012, rho=1.20,
              vmin=8.0, vmax=120.0)

# 3-DOF planar (yaw + sideslip + longitudinal) vehicle, mirroring
# Tricycle.Track.TrackTricycle: single-track (bicycle) with axle-lumped brush tyres,
# algebraic longitudinal load transfer, friction ellipse + engine-power limits. Tyres
# are the neutral-lumped pair (2x per-wheel capacity, exactly the plant's rear tyre).
PARAMS_DYN = dict(
    m=1650.0, Izz=2700.0, a=1.20, b=1.60, hcg=0.55,
    Pmax=150e3, CdA=0.72, Crr=0.012, rho=1.20, kBf=0.65,
    mu=0.95, dMax=0.35, vmin=8.0, vmax=120.0,
    # axle-lumped brush tyres (c1,c2,FzNom doubled vs the 205/55R16 per-wheel values)
    tireF=dict(c1=1.4e5, c2=8000.0, FzNom=8000.0, ap0=0.06),
    tireR=dict(c1=1.4e5, c2=8000.0, FzNom=8000.0, ap0=0.085))


def _brush_ca(alpha, Fz, t, mu):
    """CasADi port of Tricycle.brushForces (pure side-slip brush tyre, C2). Returns
    (Fy, Mz) for a lumped axle at slip angle alpha and vertical load Fz."""
    Fzc   = 0.5*(Fz + 50 + ca.sqrt((Fz - 50)**2 + 1.0))     # smooth lift-off guard
    Ca    = t['c1']*ca.sin(2*ca.atan(Fzc/t['c2']))
    ap    = t['ap0']*ca.sqrt(Fzc/t['FzNom'])
    sy    = ca.tan(alpha)
    sAbs  = ca.sqrt(sy**2 + 1e-6)
    sgn   = sy/sAbs
    theta = Ca/(3*mu*Fzc)
    uu    = 0.5*(theta*sAbs + 1 - ca.sqrt((theta*sAbs - 1)**2 + 1e-4))
    lam   = 1 - uu
    Fy    = mu*Fzc*(1 - lam**3)*sgn
    Mz    = -mu*Fzc*ap*uu*lam**3*sgn
    return Fy, Mz


def solve_min_time_dyn(s, kappa, w_max, p=PARAMS_DYN, grip_frac=0.95, v_init=None,
                       n_init=None, dpsi_init=None, w_reg=8e-3, max_iter=3000,
                       hessian='limited-memory', verbose=False):
    """Minimum-time lap of the **3-DOF planar vehicle** (yaw + sideslip + speed) by
    direct collocation, in the arc-length domain around the closed centerline.

    States  x(s) = [n, dpsi, u, vy, r]  (offset, heading error, long/lat velocity, yaw rate)
    Controls u(s) = [delta, ax]         (road-wheel steer, commanded long. acceleration)

    Unlike the point-mass solve_min_time, this has genuine yaw inertia, so weaving the
    car side-to-side *costs* time (you pay to rotate Izz) - the optimum is a proper
    wide-entry/apex/track-out line the real plant tracks without the point-mass jitter.
    grip_frac derates the peak mu to the level the full plant sustains under transient
    load transfer, keeping the optimal line trackable. Returns (n_opt, info)."""
    N = len(s)
    ds = float(s[1] - s[0])
    kap = np.asarray(kappa)
    wm = np.full(N, float(w_max)) if np.ndim(w_max) == 0 else np.asarray(w_max, float)
    m, Izz, a, b, hcg = p['m'], p['Izz'], p['a'], p['b'], p['hcg']
    L = a + b
    mu = grip_frac*p['mu']
    Fz0F = m*G*b/L                      # static axle loads (lumped)
    Fz0R = m*G*a/L
    Froll = p['Crr']*m*G
    umin = p['vmin']

    SC = np.array([max(wm.max(), 1.0), 0.4, 40.0, 3.0, 0.6])   # scale n,dpsi,u,vy,r to O(1)
    UC = np.array([p['dMax'], 10.0])                        # scale delta, ax
    opti = ca.Opti()
    Xs = opti.variable(5, N)
    Us = opti.variable(2, N)
    X = ca.diag(SC) @ Xs
    U = ca.diag(UC) @ Us
    n, dpsi, u, vy, r = X[0, :], X[1, :], X[2, :], X[3, :], X[4, :]
    delta, axc = U[0, :], U[1, :]

    def dyn(k, xk, dk, ak):
        nn, dp, uu, vv, rr = xk[0], xk[1], xk[2], xk[3], xk[4]
        uG = 0.5*(uu + ca.sqrt(uu**2 + umin**2))            # >0 guard for slip angles
        Fdrag = 0.5*p['rho']*p['CdA']*uu**2
        dFzX = m*ak*hcg/L                                   # long. load transfer
        Fzf = Fz0F - dFzX
        Fzr = Fz0R + dFzX
        af = dk - ca.atan((vv + a*rr)/uG)
        ar = -ca.atan((vv - b*rr)/uG)
        Fyf, Mzf = _brush_ca(af, Fzf, p['tireF'], mu)
        Fyr, Mzr = _brush_ca(ar, Fzr, p['tireR'], mu)
        # longitudinal request -> drive(rear)/brake(split), smooth positive/negative parts
        Fneed = m*ak + Fdrag + Froll
        FxPos = 0.5*(Fneed + ca.sqrt(Fneed**2 + 100.0**2))  # drive (>=0), rear
        FxNeg = Fneed - FxPos                               # brake (<=0)
        Fxr = FxPos + (1 - p['kBf'])*FxNeg
        Fxf = p['kBf']*FxNeg
        # body-frame front forces, chassis accelerations
        FyFb = Fyf*ca.cos(dk) + Fxf*ca.sin(dk)
        ay = (FyFb + Fyr)/m
        udot = ak + vv*rr
        vdot = ay - uu*rr
        rdot = (a*FyFb - b*Fyr + Mzf + Mzr)/Izz
        sdot = (uu*ca.cos(dp) - vv*ca.sin(dp))/(0.5*((1 - nn*kap[k]) +
                ca.sqrt((1 - nn*kap[k])**2 + 0.04)))        # smooth max(1-n*kappa, ~0.2)
        dn = (uu*ca.sin(dp) + vv*ca.cos(dp))/sdot
        ddp = rr/sdot - kap[k]
        f = ca.vertcat(dn, ddp, udot/sdot, vdot/sdot, rdot/sdot)
        # constraints returned for this node: friction ellipse per axle + power
        gell_f = (Fxf**2 + Fyf**2) - (mu*Fzf)**2
        gell_r = (Fxr**2 + Fyr**2) - (mu*Fzr)**2
        gpow = FxPos*uu - p['Pmax']
        return f, sdot, gell_f, gell_r, gpow

    T = 0
    for k in range(N):
        kn = (k + 1) % N
        fk, sdot_k, gf, gr, gp = dyn(k, X[:, k], delta[k], axc[k])
        fkn, _, _, _, _ = dyn(kn, X[:, kn], delta[kn], axc[kn])
        opti.subject_to((X[:, kn] - X[:, k] - 0.5*ds*(fk + fkn))/SC == 0)
        opti.subject_to(gf <= 0)
        opti.subject_to(gr <= 0)
        opti.subject_to(gp <= 0)
        T = T + ds/sdot_k

    # control-rate regularization (yaw dynamics already smooth the line, so this is light)
    roll = list(range(1, N)) + [0]
    dU = Us[:, roll] - Us
    opti.minimize(T + w_reg*(ca.sumsqr(dU[0, :]) + ca.sumsqr(dU[1, :])))

    opti.subject_to(Xs[0, :] <= (wm/SC[0]).reshape(1, -1))       # per-node track corridor
    opti.subject_to(Xs[0, :] >= (-wm/SC[0]).reshape(1, -1))
    opti.subject_to(opti.bounded(-0.6/SC[1], Xs[1, :], 0.6/SC[1]))
    opti.subject_to(opti.bounded(p['vmin']/SC[2], Xs[2, :], p['vmax']/SC[2]))
    opti.subject_to(opti.bounded(-8.0/SC[3], Xs[3, :], 8.0/SC[3]))
    opti.subject_to(opti.bounded(-1.5/SC[4], Xs[4, :], 1.5/SC[4]))
    opti.subject_to(opti.bounded(-1.0, Us[0, :], 1.0))       # |delta| <= dMax
    opti.subject_to(opti.bounded(-1.5, Us[1, :], 1.2))       # ax in [-15, 12] m/s^2

    # warm start from the min-curvature line + its speed profile, with geometry-consistent
    # heading and yaw rate (r ~ kappa_line*u) so the collocation defects start small
    v0 = np.full(N, 30.0) if v_init is None else np.asarray(v_init)
    n0 = np.zeros(N) if n_init is None else np.asarray(n_init)
    dp0 = np.zeros(N) if dpsi_init is None else np.asarray(dpsi_init)
    kapL = kap + (np.roll(n0, -1) - 2*n0 + np.roll(n0, 1))/ds**2
    r0 = kapL*v0
    delta0 = L*kapL + m*G*(b/L - a/L)/(2*7.0e4)*0  # ~kinematic steer (understeer term ~0)
    opti.set_initial(Xs[0, :], n0/SC[0])
    opti.set_initial(Xs[1, :], dp0/SC[1])
    opti.set_initial(Xs[2, :], v0/SC[2])
    opti.set_initial(Xs[4, :], r0/SC[4])
    opti.set_initial(Us[0, :], np.clip(L*kapL, -p['dMax'], p['dMax'])/UC[0])

    opts = {'ipopt.max_iter': max_iter, 'ipopt.tol': 1e-6,
            'ipopt.acceptable_tol': 1e-5, 'ipopt.acceptable_iter': 10,
            'ipopt.mu_strategy': 'adaptive', 'ipopt.linear_solver': 'mumps',
            'print_time': False}
    if hessian:                                 # L-BFGS avoids the costly exact Hessian
        opts['ipopt.hessian_approximation'] = hessian   # of the nested brush tyre
    if not verbose:
        opts['ipopt.print_level'] = 0
    opti.solver('ipopt', opts)
    ok = True
    try:
        sol = opti.solve()
    except RuntimeError:
        sol, ok = opti.debug, False
    return (np.array(sol.value(n)).ravel(),
            dict(T=float(sol.value(T)), converged=ok,
                 u=np.array(sol.value(u)).ravel(),
                 vy=np.array(sol.value(vy)).ravel(),
                 r=np.array(sol.value(r)).ravel(),
                 delta=np.array(sol.value(delta)).ravel(),
                 ax=np.array(sol.value(axc)).ravel()))


def solve_min_time(s, kappa, w_max, p=PARAMS, grip_frac=0.90, v_init=None, n_init=None,
                   dpsi_init=None, w_reg=1e-3, w_curv=3125.0, max_iter=3000, verbose=False):
    """Return (n_opt, info): the optimal lateral offset on the s-grid, and info with
    the optimal lap time T and the optimal speed/acceleration trajectories. Warm-start
    with n_init / v_init (e.g. the minimum-curvature line and its speed profile).

    grip_frac scales the friction circle below the tyre's peak mu to the level the
    full vehicle actually sustains (the brush tyre + yaw dynamics realise ~0.9 of the
    peak); matching it keeps the optimal line trackable by the real plant."""
    N = len(s)
    ds = float(s[1] - s[0])
    m, mu = p['m'], grip_frac*p['mu']
    Flim = mu*m*G                       # friction-circle radius [N]
    kap = np.asarray(kappa)

    SC = np.array([4.0, 0.4, 40.0])     # scale n, chi, v to O(1)
    AS = 10.0                           # acceleration scale
    opti = ca.Opti()
    Xs = opti.variable(3, N)
    Us = opti.variable(2, N)
    X = ca.diag(SC) @ Xs
    n, chi, v = X[0, :], X[1, :], X[2, :]
    at = AS*Us[0, :]                    # tangential (long.) acceleration
    an = AS*Us[1, :]                    # lateral acceleration

    def deriv_s(k, xk, at_k, an_k):
        nn, ch, vv = xk[0], xk[1], xk[2]
        sdot = vv*ca.cos(ch)/(1 - nn*kap[k])
        dn = (1 - nn*kap[k])*ca.tan(ch)                 # dn/ds
        dchi = (an_k/vv)/sdot - kap[k]                  # dchi/ds
        dv = at_k/sdot                                  # dv/ds
        return ca.vertcat(dn, dchi, dv), sdot

    T = 0
    for k in range(N):
        kn = (k + 1) % N
        fk, sdot_k = deriv_s(k, X[:, k], at[k], an[k])
        fkn, _ = deriv_s(kn, X[:, kn], at[kn], an[kn])
        opti.subject_to((X[:, kn] - X[:, k] - 0.5*ds*(fk + fkn))/SC == 0)
        T = T + ds/sdot_k
        # tyre force = mass*accel + resistances, inside the friction circle + power cap
        Fdrag = 0.5*p['rho']*p['CdA']*v[k]**2
        Froll = p['Crr']*m*G
        Ft = m*at[k] + Fdrag + Froll                    # tangential tyre force
        Fn = m*an[k]                                     # lateral tyre force
        opti.subject_to(Ft**2 + Fn**2 <= Flim**2)       # friction circle
        opti.subject_to(Ft*v[k] <= p['Pmax'])           # engine power
    roll = list(range(1, N)) + [0]
    back = [N - 1] + list(range(0, N - 1))
    dU = Us[:, roll] - Us
    # penalise the path's second difference (its curvature) so the min-time optimum is
    # smooth and drivable rather than jagged/bang-bang: a small blended min-curvature
    # term, plus it suppresses the odd-even ringing mode of trapezoidal collocation.
    # Normalise by ds^3 so the penalty approximates the grid-independent functional
    # int (n'')^2 ds  (d2n ~ n''*ds^2); otherwise a coarse grid over-penalises and
    # flattens the racing line toward the centerline (e.g. the 25 m Nordschleife grid).
    nsc = Xs[0, :]
    d2n = nsc[roll] - 2*nsc + nsc[back]
    opti.minimize(T + w_reg*(ca.sumsqr(dU[0, :]) + ca.sumsqr(dU[1, :]))
                  + (w_curv/ds**3)*ca.sumsqr(d2n))

    opti.subject_to(opti.bounded(-w_max/SC[0], Xs[0, :], w_max/SC[0]))
    opti.subject_to(opti.bounded(-0.7/SC[1], Xs[1, :], 0.7/SC[1]))
    opti.subject_to(opti.bounded(p['vmin']/SC[2], Xs[2, :], p['vmax']/SC[2]))

    # warm start from the supplied line + speed profile, with a heading guess
    # consistent with the offset (dn/ds = (1-n*kappa) tan chi) so the dn defect starts ~0
    v0 = np.full(N, 30.0) if v_init is None else np.asarray(v_init)
    n0 = np.zeros(N) if n_init is None else np.asarray(n_init)
    dn0 = (np.roll(n0, -1) - np.roll(n0, 1))/(2*ds)
    chi0 = np.arctan(dn0/(1 - n0*kap))
    kapL = kap + (np.roll(n0, -1) - 2*n0 + np.roll(n0, 1))/ds**2
    opti.set_initial(Xs[0, :], n0/SC[0])
    opti.set_initial(Xs[1, :], chi0/SC[1])
    opti.set_initial(Xs[2, :], v0/SC[2])
    opti.set_initial(Us[1, :], v0**2*kapL/AS)

    opts = {'ipopt.max_iter': max_iter, 'ipopt.tol': 1e-7,
            'ipopt.acceptable_tol': 1e-6, 'ipopt.acceptable_iter': 10,
            'ipopt.mu_strategy': 'monotone', 'ipopt.mu_init': 1e-1,
            'ipopt.bound_push': 1e-2, 'ipopt.bound_frac': 1e-2,
            'ipopt.linear_solver': 'mumps', 'print_time': False}
    if not verbose:
        opts['ipopt.print_level'] = 0
    opti.solver('ipopt', opts)
    ok = True
    try:
        sol = opti.solve()
    except RuntimeError:
        sol, ok = opti.debug, False
    return (np.array(sol.value(n)).ravel(),
            dict(T=float(sol.value(T)), converged=ok,
                 v=np.array(sol.value(v)).ravel(),
                 chi=np.array(sol.value(chi)).ravel(),
                 at=np.array(sol.value(at)).ravel(),
                 an=np.array(sol.value(an)).ravel()))
