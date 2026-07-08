"""Vehicle parameter sets ("cars") for the track-lap pipeline.

One car drives two consumers that must agree:
  * `profile` - the quasi-steady speed profile and the 3-DOF min-time OCP
    (tracks/speed_profile.py, tracks/opt_lap.py): mass, power, drag, grip, geometry.
  * `trike`   - the OpenModelica plant `Track.TrackTricycle` (component `trike` inside
    Examples.TrackLap). These become `-override trike.<name>=...` at simulate time, so
    the simulated car matches the line/speed the profile was built for.

`Pmax` is a TrackLap-level parameter (overridden as `Pmax=`, not `trike.Pmax`).

The `tourer` values are exactly the Tricycle.mo defaults, so selecting it reproduces the
original 1650 kg / 150 kW results. `elise` is tuned to a track-day Lotus Elise from real
VBOX telemetry at Knutstorp (peak ~1.45 g lateral, ~1.1 g braking, power-limited
acceleration; best real lap ~1:11): 862 kg, 118 hp, mid-engine ~38% front, µ ≈ 1.45.
"""
G = 9.80665


def _car(display, m, Pmax, CdA, Crr, rho, mu, ayFrac, hcg, a, b,
         Izz, tf, xiF, kBf, c1F, c1R, c2F, c2R, FzNomF, FzNomR, ap0F, ap0R,
         Kus, KLA=0.10, Kr=0.6, grip_frac=0.93):
    # steady understeer gradient from the bicycle model, Kus = Wf/Caf - Wr/Car
    # [rad/(m/s^2)], used by the driver's steer feedforward delta = (L + Kus*v^2)*kappa.
    # grip_frac derates nominal mu to what the full plant sustains transiently (the OCP
    # target); tyre load sensitivity means the axle realizes < mu*Fz at the limit.
    return dict(
        display=display,
        Pmax=Pmax,                       # TrackLap-level
        grip_frac=grip_frac,
        profile=dict(m=m, Pmax=Pmax, CdA=CdA, Crr=Crr, rho=rho, mu=mu,
                     ayFrac=ayFrac, hcg=hcg, a=a, b=b),
        trike=dict(m=m, Izz=Izz, a=a, b=b, tf=tf, hcg=hcg, xiF=xiF, kBf=kBf,
                   CdA=CdA, Crr=Crr, rho=rho, muF=mu, muR=mu,
                   c1F=c1F, c1R=c1R, c2F=c2F, c2R=c2R,
                   FzNomF=FzNomF, FzNomR=FzNomR, ap0F=ap0F, ap0R=ap0R),
        driver=dict(Lwb=a + b, Kus=Kus, KLA=KLA, Kr=Kr),
    )


def ocp_params(car):
    """Build the 3-DOF OCP parameter dict (tracks/opt_lap.solve_min_time_dyn) from a car.
    The OCP tyres are per-AXLE: front axle = 2x the Modelica per-wheel front, rear axle =
    the lumped rear. For `tourer` this reproduces opt_lap.PARAMS_DYN exactly."""
    c = CARS[car]; t = c['trike']
    return dict(
        m=t['m'], Izz=t['Izz'], a=t['a'], b=t['b'], hcg=t['hcg'],
        Pmax=c['Pmax'], CdA=t['CdA'], Crr=t['Crr'], rho=t['rho'], kBf=t['kBf'],
        mu=t['muF'], dMax=0.35, vmin=8.0, vmax=120.0,
        tireF=dict(c1=2*t['c1F'], c2=2*t['c2F'], FzNom=2*t['FzNomF'], ap0=t['ap0F']),
        tireR=dict(c1=t['c1R'], c2=t['c2R'], FzNom=t['FzNomR'], ap0=t['ap0R']))


CARS = {
    'tourer': _car(
        'Sport Tourer (1650 kg, 150 kW)',
        m=1650, Pmax=150e3, CdA=0.72, Crr=0.012, rho=1.20, mu=0.95, ayFrac=0.90,
        hcg=0.55, a=1.20, b=1.60,
        Izz=2700, tf=1.55, xiF=0.60, kBf=0.65,
        c1F=7.0e4, c1R=1.4e5, c2F=4000, c2R=8000,
        FzNomF=4000, FzNomR=8000, ap0F=0.06, ap0R=0.085,
        Kus=1.63e-3),

    'elise': _car(
        'Lotus Elise (862 kg, 118 hp)',
        m=862, Pmax=88e3, CdA=0.68, Crr=0.012, rho=1.20, mu=1.60, ayFrac=0.90,
        hcg=0.46, a=1.43, b=0.87,          # L=2.30 m, mid-engine ~38% front
        Izz=1070, tf=1.46, xiF=0.50, kBf=0.62,
        # track tyres: high µ; nominal loads scaled to the Elise's axle loads
        # (front per wheel ~1600 N, lumped rear ~5250 N), stiffness kept ~proportional
        c1F=2.8e4, c1R=9.2e4, c2F=1600, c2R=5250,
        FzNomF=1600, FzNomR=5250, ap0F=0.055, ap0R=0.075,
        Kus=1.0e-3,           # near-neutral mid-engine, mild understeer for stability
        KLA=0.30,             # higher lookahead gain: the grippier car wants tighter feedback
        grip_frac=0.88),      # plant realizes ~1.40 g of the 1.60 nominal mu
}


def override_string(car, sLap, u0):
    """Build the OMC simflags -override list for a car (plus run-specific sLap/u0)."""
    c = CARS[car]
    parts = [f'sLap={sLap:.1f}', f'u0={u0:.2f}', f'Pmax={c["Pmax"]:.0f}']
    parts += [f'trike.{k}={v:g}' for k, v in c['trike'].items()]
    parts += [f'driver.{k}={v:g}' for k, v in c['driver'].items()]
    return '-override ' + ','.join(parts)
