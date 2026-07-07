# tricycle-dlc

![tests](https://github.com/pfyhr/tricycle-dlc/actions/workflows/ci.yml/badge.svg)

Minimal, defensible planar vehicle model for studying tire forces on the suspension
links and steering arm during an **ISO 3888-1 double lane change**, with a traditional
unassisted rack-and-pinion** steering.

## Model

`modelica/Tricycle.mo` (single-file Modelica package, simulated with OpenModelica):

- **`PlanarTricycle`** — planar three-wheel ("tricycle") vehicle: individual front
  wheels (own slip angles, quasi-static lateral load transfer with a roll-mode lag),
  lumped rear wheel — the architecture used for front-axle force estimation in
  WO 2025/113783 (Marzbanrad & Jonasson, Volvo). Chassis: 2 DOF (lateral velocity, yaw
  rate) at constant speed + path kinematics. Kingpin moment per side =
  Fy·(mechanical trail) − Mz; tie-rod force = M_kp/L_arm. `toeL`/`toeR` inputs are
  hooks for an active-toe actuator (wired to 0 here).
- **`TireData` / `brushForces`** — Pacejka brush tire (Tire and Vehicle Dynamics,
  Ch. 3): analytic Fy(α, Fz) and Mz(α, Fz) with a pneumatic trail that starts at
  a_p/3 ≈ 20 mm and collapses to zero at the grip limit; degressive load sensitivity;
  first-order relaxation lag. Smooth and event-free.
- **`ManualSteering`** — handwheel/column inertia + ideal pinion
  (r_p = L_arm/i_S, i_S = 20 typical for unassisted steering) + rack mass. No assist:
  the full kingpin reaction reflects to the handwheel, τ_HW = F_rack·r_p.
- **`Iso3888Path`** — ISO 3888-1 reference centerline (sections 15/30/25/25/15/15 m,
  3.5 m offset, 125 m total) + MacAdam-class single-point preview driver with yaw-rate
  damping. The driver turns the handwheel through a 2 Hz arm filter.
- **Examples** — `StepSteer` (understeer validation), `DoubleLaneChange` (headline,
  closed-loop at the ISO-recommended 80 km/h), `OpenLoopDLC` (prescribed one-period
  sine, repeatable sweeps).

![Double lane change animation](outputs/gif/dlc_anim.gif)

*Closed-loop ISO 3888-1 double lane change at 80 km/h (real-time playback): vehicle
outline at true yaw with the running time, heading and lateral-acceleration readout.
The distance axis is compressed 2:1; gate compliance is checked on the true footprint.*

## Headline results (defaults: D-segment sedan, 80 km/h)

| Quantity | Value |
|---|---|
| ISO 3888-1 gates (full-footprint check) | PASS, margins +91/+68/+151 mm |
| Peak lateral acceleration | 0.76 g |
| Peak tie-rod force | ≈ 1.4 kN per side (left/right split by load transfer) |
| Peak kingpin moment | ≈ 150 N·m per wheel |
| Peak handwheel angle / torque | ≈ 110° / **10 N·m** (unassisted!) |
| Understeer gradient (tires-only) | 0.96 deg/g — low edge of the production band, as expected with compliance/roll steer omitted |

Note the driver tuning (Tp = 0.55 s, Kdrv = 0.22, Kr = 0.25) is deliberately slow and
well damped: the 2 Hz arm filter adds phase lag that the preview must compensate, and
tighter tunings destabilize the driver–steering loop — a genuine interaction, not a
numerical artifact.

## Track sim: minimum-time laps of planar circuits

The `Tricycle.Track` sub-package re-expresses the tricycle in **track (Frenet)
coordinates** (s, n, Δψ) and adds a **longitudinal degree of freedom**: rear-wheel
drive limited by engine power (P_max/u) and by the rear friction-ellipse remainder
√((μF_z)² − F_y²), brakes split front/rear under the same per-axle ellipse limit
("ideal TC/ABS"), aero drag, rolling resistance, and pitch-lagged longitudinal load
transfer. Centerlines come from OpenStreetMap (elevation dropped — planar by design;
© OpenStreetMap contributors, ODbL), smoothed and tabulated as κ(s) in
`tracks/<key>.csv`.

The driver drives *his* fastest lap for the given setup, not a formal optimum: a
quasi-steady minimum-time speed profile v_ref(s) is precomputed for the exact vehicle
parameters (corner-speed limit → power/traction-limited forward pass →
braking-limited backward pass, `tracks/speed_profile.py`), and the two-channel
`TrackDriver` tracks it — curvature feedforward + gain-scheduled offset feedback +
yaw-rate-error damping for steering, and a preview-consistent constant-acceleration
law for throttle/brake. He stays on the centerline instead of using the track width,
so on tight hairpins he runs a metre or two wide — realistic for a line that ignores
the racing groove.

Lap times for the default setup (150 kW / 1650 kg / μ = 0.95):

| Track (`--track=`) | Length | Lap | Quasi-steady ideal | v_max |
|---|--:|--:|--:|--:|
| `nordschleife` — Nürburgring Nordschleife | 20.72 km | **11:02.3** | 10:53.7 | 229 km/h |
| `anderstorp` — Anderstorp Raceway | 4.01 km | **2:20.7** | 2:18.5 | 175 km/h |
| `gelleras` — Gelleråsen Arena (Karlskoga) | 2.33 km | **1:40.4** | 1:38.3 | 156 km/h |
| `knutstorp` — Ring Knutstorp | 2.06 km | **1:33.6** | 1:31.4 | 159 km/h |
| `kinnekulle` — Kinnekulle Ring | 2.06 km | **1:19.2** | 1:17.9 | 153 km/h |

The lap is always within ~2 % of the quasi-steady ideal — the driver's overhead for
steering the real line rather than sitting exactly on the profile.

![Nordschleife map colored by speed](outputs/svg/ns_map.svg)

```
python3 tracks/fetch_track.py --track=all    # (re)build centerlines from OSM - needs network
python3 track_lap.py    --track=knutstorp    # speed profile + lap sim + figures + summary CSV
python3 track_render.py --track=knutstorp    # chase-camera HTML viewer (outputs/<key>_chase.html):
                                             # GTA-style follow cam, minimap, speed/yaw/accel HUD
```

Adding a track is one entry in the `TRACKS` registry in `tracks/fetch_track.py`
(an OSM route relation, or a bounding box whose raceway ways are auto-assembled into
the closed circuit loop). Vehicle setup is sweepable: `track_lap.py` mirrors the
`Tricycle.Track.TrackTricycle` defaults, and per-run overrides pass straight through
to `simulate(..., simflags="-override Pmax=...")`.

## Run

```
python3 dlc_maneuver.py      # simulations, figures (outputs/svg), animation (outputs/gif), summary CSV
python3 render_diagram.py    # component diagram from the OpenModelica instance
python3 tests/run_tests.py   # test suite (also runs in GitHub Actions on every push)
```

Requires OpenModelica (`omc` on PATH), Python 3 with numpy + matplotlib.
Everything sweepable (speed, driver, steering ratio, tire and chassis parameters) is a
top-level parameter of the examples, overridable per run via `-override=...`.

## Validation

- Steady-state yaw-rate gain matches the analytic linear bicycle *including aligning
  moments* to < 1 % across 40–120 km/h:

  ![Yaw-rate gain vs speed](outputs/svg/dlc_understeer.svg)

- A NumPy twin of the brush tire is asserted against the Modelica steady state to
  < 10⁻³ on every run of `dlc_maneuver.py`.
- Gate compliance is checked on the full yawed vehicle footprint (all four corners),
  matching the ISO "no cones displaced" intent.

## Sources

Parameter provenance and reference anchors in `sources/SOURCES.md`
(Pacejka, Heydinger SAE 1999-01-1336, ISO 3888-1:1999 incl. the recommended
(80 ± 3) km/h entry speed, Milliken, MacAdam, Reimpell, a measured DLC field test
[Bîndac et al. 2022], WO 2025/113783).

## Documented omissions

Constant speed; parallel steer (no Ackermann); no KPI/caster jacking; no scrub×Fx;
lumped rear tire (slightly understeer-optimistic; front link loads unaffected); no roll
DOF (roll-mode lag on load transfer instead); no steering column compliance or rack
friction; tire parameters are textbook-typical, not fitted to a specific tire.
