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
         Kus, KLA=0.10, Kr=0.6, grip_frac=0.93, ClA=0.0, aeroBal=0.42):
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
                   FzNomF=FzNomF, FzNomR=FzNomR, ap0F=ap0F, ap0R=ap0R,
                   ClA=ClA, aeroBal=aeroBal),
        driver=dict(Lwb=a + b, Kus=Kus, KLA=KLA, Kr=Kr, mu=mu),
    )


def ocp_params(c):
    """Build the 3-DOF OCP parameter dict (tracks/opt_lap.solve_min_time_dyn) from a config
    (build_config output). The OCP tyres are per-AXLE: front axle = 2x the Modelica per-wheel
    front, rear axle = the lumped rear. For `tourer`/base this reproduces PARAMS_DYN exactly.
    ClA/aeroBal (downforce) pass through so the friction ellipse grows with speed."""
    t = c['trike']
    return dict(
        m=t['m'], Izz=t['Izz'], a=t['a'], b=t['b'], hcg=t['hcg'],
        Pmax=c['Pmax'], CdA=t['CdA'], Crr=t['Crr'], rho=t['rho'], kBf=t['kBf'],
        mu=t['muF'], dMax=0.35, vmin=8.0, vmax=120.0,
        ClA=t.get('ClA', 0.0), aeroBal=t.get('aeroBal', 0.45),
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
        'Lotus Elise (862 kg, 88 kW)',
        # mu/ayFrac calibrated against the owner's RaceLogs (racelogs.py): best real
        # Knutstorp lap 1:09.5 over 4 consistent sessions, sustained |ay| p99 = 1.62 g,
        # vmax 163 km/h, braking p90 0.9 g. R-compound at light Elise wheel loads ->
        # mu 1.80; ayFrac 0.96 puts the profile budget at what the plant realizes (1.5 g).
        m=862, Pmax=88e3, CdA=0.68, Crr=0.012, rho=1.20, mu=1.80, ayFrac=0.96,
        hcg=0.46, a=1.43, b=0.87,          # L=2.30 m, mid-engine ~38% front
        Izz=1070, tf=1.46, xiF=0.50, kBf=0.62,
        # track tyres: high µ; nominal loads scaled to the Elise's axle loads
        # (front per wheel ~1600 N, lumped rear ~5250 N), stiffness kept ~proportional
        c1F=2.8e4, c1R=9.2e4, c2F=1600, c2R=5250,
        FzNomF=1600, FzNomR=5250, ap0F=0.055, ap0R=0.075,
        Kus=1.0e-3,           # near-neutral mid-engine, mild understeer for stability
        KLA=0.30,             # higher lookahead gain: the grippier car wants tighter feedback
        grip_frac=0.88),      # plant realizes ~1.5 g of the 1.80 nominal mu

    'miata': _car(
        'Mazda MX-5 "Oskar" (980 kg, 104 kW)',
        # calibrated against Oskar's real Knutstorp lap 19 (racelogs.py): 1:08.0,
        # sustained |ay| p99 = 1.67 g, vmax 167 km/h - the fastest real lap on file.
        m=980, Pmax=104e3, CdA=0.66, Crr=0.012, rho=1.20, mu=1.82, ayFrac=0.97,
        hcg=0.48, a=1.15, b=1.15,          # L=2.30 m, front-engine RWD ~50/50
        Izz=1300, tf=1.41, xiF=0.50, kBf=0.62,
        # track tyres: realized ~1.55 g lateral; loads scaled to ~2400 N front / 4800 N rear
        c1F=4.2e4, c1R=8.4e4, c2F=2400, c2R=4800,
        FzNomF=2400, FzNomR=4800, ap0F=0.055, ap0R=0.075,
        Kus=1.2e-3, KLA=0.30, grip_frac=0.88),

    # BMW M140i: 340 hp straight-six, ~1530 kg with driver, sport street tyres.
    # Spec-based; mu 1.30 matches the car's logged |ay| p99 = 1.28 g (Falkenberg VBO).
    # NOTE the 2024 GRC Kinnekulle lap (0:57.6, racelogs.py) was NOT used to fit:
    # its speed trace shows a corner at 163 km/h where our centerline has a 103 km/h
    # bend - the track was rebuilt since our layout data; see falkenberg for a fit
    # target once that track is added.
    'm140': _car(
        'BMW M140i (1530 kg, 250 kW)',
        m=1530, Pmax=250e3, CdA=0.66, Crr=0.011, rho=1.20, mu=1.30, ayFrac=0.95,
        hcg=0.52, a=1.29, b=1.40,           # L=2.69 m, front-engine RWD ~52/48
        Izz=2600, tf=1.57, xiF=0.55, kBf=0.64,
        c1F=8.0e4, c1R=1.5e5, c2F=3700, c2R=7300,
        FzNomF=3700, FzNomR=7300, ap0F=0.06, ap0R=0.08,
        Kus=1.5e-3, KLA=0.26, grip_frac=0.90),

    # Scandinavian Clubman sports-prototype: 1800 Ford Zetec ~155 hp, 580 kg race weight
    # (real Swedish class figures), slicks and a big rear wing -> one of the fastest
    # sportscar classes in Scandinavia. ~200 W/kg, twice the Elise, plus downforce.
    'clubman': _car(
        'Clubman Racer (580 kg, 116 kW)',
        m=580, Pmax=116e3, CdA=0.85, Crr=0.012, rho=1.20, mu=1.75, ayFrac=0.92,
        hcg=0.30, a=1.20, b=1.20,           # L=2.40 m, ~50/50, very low CG
        Izz=780, tf=1.55, xiF=0.48, kBf=0.60,
        # slicks; nominal loads scaled to the light axle loads (~1420 N front, ~2840 N rear)
        c1F=2.6e4, c1R=5.2e4, c2F=1420, c2R=2840,
        FzNomF=1420, FzNomR=2840, ap0F=0.05, ap0R=0.07,
        Kus=0.8e-3,            # near-neutral, twitchy
        KLA=0.32, grip_frac=0.90,
        ClA=0.5, aeroBal=0.40),   # wing + flat floor; light downforce the lookahead driver holds at speed
}


import copy
import math

# Aerodynamic cost of downforce: classical lift-induced drag, CdA_i = ClA^2/(pi*e*b^2)
# per lifting surface, with the span b capped at the car's overall width - race wings
# are width-limited, which is why they run aspect ratios near 2 (Katz 2006, Annu. Rev.
# Fluid Mech. 38:27-63). e ~ 0.9 with endplates (endplates raise the effective span;
# Hoerner 1985 via Katz 2006). The package is modeled as TWO such surfaces - front and
# rear, carrying aeroBal and (1-aeroBal) of the lift - so the induced drag is the sum
# of two per-wing terms. This makes aero balance a real trade-off: concentrating the
# lift on one span-limited wing costs more drag than splitting it. At the Clubman's
# ClA = 0.5 m^2 / 40% front this costs ~0.015 m^2 (package L/D ~ 33, the efficient
# wing-plus-floor end); a slider-maxed 5 m^2 at 45% front costs ~1.5 m^2 (package
# L/D ~ 3.4 - the whole-car efficiency of a modern F1 package).
AERO_E = 0.9
AERO_SPAN_PAD = 0.20      # overall width ~ front track + tyre width


def induced_drag(ClA, tf, aeroBal=0.45):
    """Lift-induced drag area [m^2] of a front+rear downforce package totaling ClA
    on a car of front track tf, with aeroBal of the lift on the front surface."""
    b = tf + AERO_SPAN_PAD
    qF, qR = aeroBal*ClA, (1.0 - aeroBal)*ClA
    return (qF**2 + qR**2)/(math.pi*AERO_E*b**2)


# Setup options layered on any car. Each is a single-axis variation from baseline: pure
# parameter changes that the Modelica plant already supports. (Downforce lives only in the
# live JS simulator - build_webgui.py - since it adds speed-dependent grip the Modelica
# plant doesn't model yet. ClA/aeroBal are still threaded through so the OCP could use it.)
SETUPS = {
    'base':       dict(label='Baseline'),
    'ballast_lo': dict(label='Ballast −50 kg', dm=-50.0),
    'ballast_hi': dict(label='Ballast +50 kg', dm=+50.0),
    'oversteer':  dict(label='Loose (oversteer)', dxiF=-0.12),
    'understeer': dict(label='Tight (understeer)', dxiF=+0.12),
}


def build_config(car, setup='base'):
    """Return a deep copy of a car with a setup applied, plus its derived OCP params.
    Keys: display, setup, setup_label, Pmax, grip_frac, profile, trike, driver, ocp."""
    c = copy.deepcopy(CARS[car])
    s = SETUPS[setup]
    muS = s.get('mu_scale', 1.0)
    dm, dxiF, dCdA = s.get('dm', 0.0), s.get('dxiF', 0.0), s.get('dCdA', 0.0)
    ClA, aeroBal = s.get('ClA', 0.0), s.get('aeroBal', 0.45)
    p, t = c['profile'], c['trike']
    p['mu'] *= muS;  p['m'] += dm;  p['CdA'] += dCdA
    t['muF'] *= muS; t['muR'] *= muS; t['m'] += dm; t['CdA'] += dCdA
    t['xiF'] = min(0.9, max(0.1, t['xiF'] + dxiF))
    if ClA:                       # downforce: profile + plant + OCP all see the aero map
        p['ClA'], p['aeroBal'] = ClA, aeroBal
        t['ClA'], t['aeroBal'] = ClA, aeroBal
    # Downforce is not free: add the lift-induced drag of the final ClA. The presets
    # exported to the web sim carry the BASE CdA (c['CdA_base']) and the browser
    # applies the identical formula live as the slider moves (webgui liveParams) -
    # keep the two in sync.
    c['CdA_base'] = t['CdA']
    dCdAi = induced_drag(t.get('ClA', 0.0), t['tf'], t.get('aeroBal', 0.45))
    p['CdA'] += dCdAi
    t['CdA'] += dCdAi
    c['ocp'] = ocp_params(c)
    c['setup'], c['setup_label'] = setup, s['label']
    return c


def override_string(c, sLap, u0):
    """Build the OMC simflags -override list from a config (plus run-specific sLap/u0)."""
    parts = [f'sLap={sLap:.1f}', f'u0={u0:.2f}', f'Pmax={c["Pmax"]:.0f}']
    parts += [f'trike.{k}={v:g}' for k, v in c['trike'].items()]
    parts += [f'driver.{k}={v:g}' for k, v in c['driver'].items()]
    return '-override ' + ','.join(parts)
