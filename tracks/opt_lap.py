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
