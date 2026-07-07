within ;
package Tricycle
  "Planar three-wheel (tricycle) vehicle with manual rack-and-pinion steering, for suspension-link force studies (ISO 3888-1 double lane change)"
  extends Modelica.Icons.Package;


  record TireData
    "Brush-model tire parameters (Pacejka, Tire and Vehicle Dynamics, Ch. 3)"
    extends Modelica.Icons.Record;
    Real c1(unit="N/rad") = 7.0e4
      "Cornering-stiffness scale: Calpha = c1*sin(2*atan(Fz/c2))";
    Modelica.Units.SI.Force c2 = 4000 "Load at which Calpha peaks";
    Real mu = 0.95 "Tire-road friction coefficient (dry asphalt)";
    Modelica.Units.SI.Length ap0 = 0.06
      "Contact-patch half length at FzNom (pneumatic trail at zero slip = ap0/3)";
    Modelica.Units.SI.Force FzNom = 4000 "Nominal load for patch-length scaling";
    Modelica.Units.SI.Length sigmaRel = 0.5
      "Relaxation length (first-order transient slip-angle lag)";
    annotation (Documentation(info="<html><p>Parameters of the single-friction
<b>brush</b> tire (Pacejka, <i>Tire and Vehicle Dynamics</i>, 3rd ed., Ch. 3). The
cornering stiffness uses the standard degressive load law
C<sub>&alpha;</sub> = c1&middot;sin(2&middot;atan(Fz/c2)); the contact patch half length
grows as &radic;Fz so the low-slip pneumatic trail is ap0/3 (&asymp;20 mm) at nominal
load. Defaults are 205/55R16-class passenger-car values.</p></html>"));
  end TireData;

  function brushForces
    "Steady-state brush tire: lateral force and aligning moment vs slip angle and load"
    extends Modelica.Icons.Function;
    input Modelica.Units.SI.Angle alpha "Slip angle";
    input Modelica.Units.SI.Force Fz "Vertical load";
    input TireData d "Tire parameters";
    output Modelica.Units.SI.Force Fy
      "Lateral force (sign of alpha: alpha>0 -> Fy>0, leftward)";
    output Modelica.Units.SI.Torque Mz
      "Aligning moment (opposes alpha: alpha>0 -> Mz<0)";
  protected
    constant Real epsS = 1e-6 "Slip regularization [rad]";
    constant Real epsC = 1e-4 "Smooth-clamp regularization";
    Real Fzc, Ca, ap, sy, sAbs, sgn, theta, uu, lam;
  algorithm
    Fzc   := max(Fz, 50);  // lift-off guard
    Ca    := d.c1*sin(2*atan(Fzc/d.c2));
    ap    := d.ap0*sqrt(Fzc/d.FzNom);
    sy    := tan(alpha);
    sAbs  := sqrt(sy^2 + epsS^2);
    sgn   := sy/sAbs;
    theta := Ca/(3*d.mu*Fzc);
    // smooth min(theta*|tan(alpha)|, 1): full sliding at tan(alpha) = 3*mu*Fz/Calpha
    uu    := 0.5*(theta*sAbs + 1 - sqrt((theta*sAbs - 1)^2 + epsC));
    lam   := 1 - uu;
    Fy    := d.mu*Fzc*(1 - lam^3)*sgn;
    Mz    := -d.mu*Fzc*ap*uu*lam^3*sgn;
    annotation (smoothOrder=2, Documentation(info="<html><p>Pacejka brush model,
pure side slip (eqs. 3.9&ndash;3.12): with &theta; = C<sub>&alpha;</sub>/(3&mu;Fz) and
u = min(&theta;|tan&alpha;|, 1),</p>
<blockquote><code>|Fy| = &mu;Fz(1-(1-u)&sup3;), &nbsp;
|Mz| = &mu;Fz&middot;a<sub>p</sub>&middot;u(1-u)&sup3;</code></blockquote>
<p>so the pneumatic trail -Mz/Fy starts at a<sub>p</sub>/3 and collapses to zero at
full sliding &mdash; the mechanism that shapes the tie-rod load mid-maneuver. The
min and sign are smoothly regularized, so the function is C&sup2; and event-free.</p></html>"));
  end brushForces;

  function yRefIso3888
    "ISO 3888-1 reference centerline: half-cosine lane-offset blends"
    extends Modelica.Icons.Function;
    input Modelica.Units.SI.Position x "Distance along the course";
    input Modelica.Units.SI.Length offsetY = 3.5 "Lateral offset between lanes";
    output Modelica.Units.SI.Position y "Reference lateral position";
    output Real dydx "Reference path slope dy/dx";
  protected
    constant Real pi = Modelica.Constants.pi;
    // ISO 3888-1:1999 Table 1: 15 m lane / 30 m transition / 25 m offset lane
    // (3.5 m lane offset) / 25 m transition / 2 x 15 m exit lane, total 125 m
    constant Modelica.Units.SI.Length x1 = 15, x2 = 45, x3 = 70, x4 = 95;
  algorithm
    if x <= x1 then
      y := 0;
      dydx := 0;
    elseif x <= x2 then
      y := offsetY/2*(1 - cos(pi*(x - x1)/(x2 - x1)));
      dydx := offsetY/2*pi/(x2 - x1)*sin(pi*(x - x1)/(x2 - x1));
    elseif x <= x3 then
      y := offsetY;
      dydx := 0;
    elseif x <= x4 then
      y := offsetY/2*(1 + cos(pi*(x - x3)/(x4 - x3)));
      dydx := -offsetY/2*pi/(x4 - x3)*sin(pi*(x - x3)/(x4 - x3));
    else
      y := 0;
      dydx := 0;
    end if;
    annotation (smoothOrder=1);
  end yRefIso3888;

  block Iso3888Path
    "ISO 3888-1 path reference + single-point preview driver -> road-wheel steer command"
    parameter Modelica.Units.SI.Velocity u = 80/3.6 "Forward speed";
    parameter Modelica.Units.SI.Time Tp = 0.45 "Preview time (MacAdam-class single point)";
    parameter Real Kdrv(unit="rad/m") = 0.3 "Road-wheel steer per meter of preview error";
    parameter Real Kr(unit="rad.s/rad") = 0.08
      "Yaw-rate damping (vestibular feedback term; zero steady-state offset on straights)";
    parameter Modelica.Units.SI.Angle dMax = 0.35 "Smooth steer-command clamp (~20 deg)";
    parameter Modelica.Units.SI.Length offsetY = 3.5 "Lane offset";
    parameter Modelica.Units.SI.Position xStart = 0 "Vehicle X at course entry";
    Modelica.Blocks.Interfaces.RealInput X "Vehicle global x [m]"
      annotation (Placement(transformation(extent={{-120,50},{-100,70}})));
    Modelica.Blocks.Interfaces.RealInput Y "Vehicle global y [m]"
      annotation (Placement(transformation(extent={{-120,10},{-100,30}})));
    Modelica.Blocks.Interfaces.RealInput psi "Vehicle heading [rad]"
      annotation (Placement(transformation(extent={{-120,-30},{-100,-10}})));
    Modelica.Blocks.Interfaces.RealInput vy "Vehicle lateral velocity [m/s]"
      annotation (Placement(transformation(extent={{-120,-70},{-100,-50}})));
    Modelica.Blocks.Interfaces.RealInput r "Vehicle yaw rate [rad/s]"
      annotation (Placement(transformation(extent={{-120,-110},{-100,-90}})));
    Modelica.Blocks.Interfaces.RealOutput dCmd "Road-wheel steer command [rad]"
      annotation (Placement(transformation(extent={{100,-10},{120,10}})));
    Modelica.Units.SI.Position yRef "Reference lateral position at the preview point";
    Real dyRef "Reference path slope at the preview point (for logging)";
    Modelica.Units.SI.Position yPrev "Predicted lateral position at the preview point";
    Modelica.Units.SI.Position ePrev "Preview position error";
  equation
    (yRef, dyRef) = yRefIso3888(X - xStart + u*Tp, offsetY);
    // consistent single-point preview: vehicle predicted Tp ahead, path read Tp ahead
    yPrev = Y + Tp*(vy*cos(psi) + u*sin(psi));
    ePrev = yRef - yPrev;
    dCmd  = dMax*tanh((Kdrv*ePrev - Kr*r)/dMax);
    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Rectangle(extent={{-100,100},{100,-100}}, lineColor={0,0,0},
          fillColor={245,245,235}, fillPattern=FillPattern.Solid),
        Line(points={{-90,-40},{-40,-40},{-10,40},{40,40},{70,-40},{90,-40}},
          color={0,0,255}, smooth=Smooth.Bezier, thickness=0.5),
        Rectangle(extent={{-90,-24},{-60,-30}}, lineColor={128,128,128}),
        Rectangle(extent={{-24,56},{24,50}}, lineColor={128,128,128}),
        Rectangle(extent={{60,-24},{90,-30}}, lineColor={128,128,128}),
        Text(extent={{-96,-60},{96,-90}}, textColor={0,0,0}, textString="ISO 3888-1"),
        Text(extent={{-150,140},{150,110}}, textColor={0,0,255}, textString="%name")}),
      Documentation(info="<html><p>Single-point preview driver (MacAdam-class):
the lateral position the vehicle will have after the preview time <code>Tp</code> is
compared with the ISO 3888-1 reference centerline evaluated at the previewed station,
and the error is turned into a road-wheel steer command with gain <code>Kdrv</code>,
smoothly clamped at <code>dMax</code>.</p></html>"));
  end Iso3888Path;

  model PlanarTricycle
    "Planar 2-DOF (sideslip + yaw) three-wheel vehicle: individual front wheels, lumped rear; brush tires; kingpin trail -> tie-rod loads"
    // Sign conventions (ISO 8855): x forward, y left, z up.
    // delta > 0, yaw rate r > 0, ay > 0  =  LEFT turn; load transfers to the RIGHT (outer) wheel.
    // Positive rack travel s steers LEFT: delta = s/Larm (parallel steer).
    parameter Modelica.Units.SI.Mass m = 1650
      "Vehicle mass (D-segment sedan, Heydinger et al. SAE 1999-01-1336)";
    parameter Modelica.Units.SI.Inertia Izz = 2700 "Yaw moment of inertia";
    parameter Modelica.Units.SI.Length a = 1.20 "CG to front axle";
    parameter Modelica.Units.SI.Length b = 1.60 "CG to rear axle";
    parameter Modelica.Units.SI.Length tf = 1.55 "Front track width";
    parameter Modelica.Units.SI.Length hcg = 0.55 "CG height above ground";
    parameter Real xiF(min=0, max=1) = 0.6
      "Front share of the lateral load transfer (roll-stiffness fraction)";
    parameter Modelica.Units.SI.Time tauRoll = 0.15
      "Roll-mode lag applied to the load transfer";
    parameter Modelica.Units.SI.Velocity u = 80/3.6 "Constant forward speed";
    parameter Modelica.Units.SI.Velocity uMin = 1
      "Low-speed guard for slip-angle denominators";
    parameter Modelica.Units.SI.Length tMech = 0.025
      "Mechanical (caster) trail: 5-6 deg caster x ~0.29 m rolling radius";
    parameter Modelica.Units.SI.Length Larm = 0.11
      "Steering-arm length: rack travel per steer radian";
    parameter Modelica.Units.SI.Inertia Jkp = 0.79
      "Road-wheel inertia about the kingpin, per wheel (Yin et al. 2024)";
    parameter TireData tireF "Front tire, per wheel";
    parameter TireData tireR(c1=1.4e5, c2=8000, FzNom=8000, ap0=0.085)
      "Lumped rear-axle tire (2x per-wheel capacity)";

    Modelica.Mechanics.Translational.Interfaces.Flange_a rack
      "Steering-rack / tie-rod connection"
      annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
    Modelica.Blocks.Interfaces.RealInput toeL(unit="rad")
      "Left-wheel toe offset (active-toe actuator hook)"
      annotation (Placement(transformation(extent={{-120,50},{-100,70}})));
    Modelica.Blocks.Interfaces.RealInput toeR(unit="rad")
      "Right-wheel toe offset (active-toe actuator hook)"
      annotation (Placement(transformation(extent={{-120,-70},{-100,-50}})));

    // chassis states
    Modelica.Units.SI.Velocity vy(start=0, fixed=true) "Lateral velocity at CG";
    Modelica.Units.SI.AngularVelocity r(start=0, fixed=true) "Yaw rate";
    Modelica.Units.SI.Angle psi(start=0, fixed=true) "Heading";
    Modelica.Units.SI.Position X(start=0, fixed=true) "Global x position";
    Modelica.Units.SI.Position Y(start=0, fixed=true) "Global y position";
    Modelica.Units.SI.Force dFz(start=0, fixed=true)
      "Front lateral load transfer (roll-lagged)";
    // tire states: relaxation-lagged slip angles (one per tire; Fy and Mz are the
    // steady-state brush outputs at the lagged angle, so they stay mutually consistent)
    Modelica.Units.SI.Angle aLagFL(start=0, fixed=true) "Front-left lagged slip angle";
    Modelica.Units.SI.Angle aLagFR(start=0, fixed=true) "Front-right lagged slip angle";
    Modelica.Units.SI.Angle aLagR(start=0, fixed=true) "Rear lagged slip angle";
    Modelica.Units.SI.Force FyFL, FyFR, FyR;
    Modelica.Units.SI.Torque MzFL, MzFR, MzR;
    // kinematics and loads
    Modelica.Units.SI.Angle dL "Left road-wheel steer angle";
    Modelica.Units.SI.Angle dR "Right road-wheel steer angle";
    Modelica.Units.SI.Angle aFL "Front-left slip angle";
    Modelica.Units.SI.Angle aFR "Front-right slip angle";
    Modelica.Units.SI.Angle aR "Rear slip angle";
    Modelica.Units.SI.Force FzFL, FzFR, FzR;
    Modelica.Units.SI.Acceleration ay "Lateral acceleration (= der(vy) + u*r)";
    Modelica.Units.SI.Force FyF "Front-axle lateral force in the body frame (both wheels)";
    // kingpin / tie-rod
    Modelica.Units.SI.Torque MkpL "Left kingpin moment resisting steer";
    Modelica.Units.SI.Torque MkpR "Right kingpin moment resisting steer";
    Modelica.Units.SI.Torque MkpMechL, MkpMechR "Mechanical-trail contribution";
    Modelica.Units.SI.Torque MkpPneuL, MkpPneuR "Pneumatic-trail contribution";
    Modelica.Units.SI.Force FtieL "Left tie-rod axial force";
    Modelica.Units.SI.Force FtieR "Right tie-rod axial force";
    Modelica.Units.SI.Position s;
    Modelica.Units.SI.Velocity v;
    Modelica.Units.SI.Acceleration accRack;
  protected
    constant Modelica.Units.SI.Acceleration g = Modelica.Constants.g_n;
    final parameter Modelica.Units.SI.Length L = a + b;
    final parameter Modelica.Units.SI.Force Fz0F = m*g*b/(2*L);
    final parameter Modelica.Units.SI.Force Fz0R = m*g*a/L;
    // u is constant, so the low-speed guard and relaxation time constants are parameters
    final parameter Modelica.Units.SI.Velocity uGuard = max(u, uMin);
    final parameter Modelica.Units.SI.Time tauF = tireF.sigmaRel/uGuard;
    final parameter Modelica.Units.SI.Time tauR = tireR.sigmaRel/uGuard;
  equation
    s = rack.s;
    v = der(s);
    accRack = der(v);
    dL = s/Larm + toeL;
    dR = s/Larm + toeR;
    // per-wheel slip angles (individual front wheels: the tricycle content)
    aFL = dL - atan((vy + a*r)/noEvent(max(u - r*tf/2, uMin)));
    aFR = dR - atan((vy + a*r)/noEvent(max(u + r*tf/2, uMin)));
    aR  =    - atan((vy - b*r)/uGuard);
    // quasi-static front load transfer, lagged with the roll mode (breaks the Fz<->Fy loop)
    tauRoll*der(dFz) + dFz = xiF*m*ay*hcg/tf;
    FzFL = Fz0F - dFz;
    FzFR = Fz0F + dFz;
    FzR  = Fz0R;
    // brush tires with first-order relaxation (sigma/u lag on the slip angle)
    tauF*der(aLagFL) + aLagFL = aFL;
    tauF*der(aLagFR) + aLagFR = aFR;
    tauR*der(aLagR)  + aLagR  = aR;
    (FyFL, MzFL) = brushForces(aLagFL, FzFL, tireF);
    (FyFR, MzFR) = brushForces(aLagFR, FzFR, tireF);
    (FyR,  MzR)  = brushForces(aLagR,  FzR,  tireR);
    // planar chassis: lateral + yaw at constant forward speed, plus path kinematics
    FyF = FyFL*cos(dL) + FyFR*cos(dR);
    m*ay = FyF + FyR;
    der(vy) = ay - u*r;
    Izz*der(r) = a*FyF - b*FyR
               + (tf/2)*(FyFL*sin(dL) - FyFR*sin(dR))
               + MzFL + MzFR + MzR;
    der(psi) = r;
    der(X) = u*cos(psi) - vy*sin(psi);
    der(Y) = u*sin(psi) + vy*cos(psi);
    // kingpin moments -> tie-rod forces -> rack reaction
    // (Mz < 0 for alpha > 0, so -Mz is the positive restoring pneumatic part)
    MkpMechL = FyFL*tMech;
    MkpPneuL = -MzFL;
    MkpL = MkpMechL + MkpPneuL;
    MkpMechR = FyFR*tMech;
    MkpPneuR = -MzFR;
    MkpR = MkpMechR + MkpPneuR;
    FtieL = MkpL/Larm;
    FtieR = MkpR/Larm;
    rack.f = FtieL + FtieR + (2*Jkp/Larm^2)*accRack;
    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Rectangle(extent={{-60,80},{60,-80}}, lineColor={0,0,0},
          fillColor={235,235,245}, fillPattern=FillPattern.Solid, radius=20),
        Rectangle(extent={{-58,72},{-30,44}}, lineColor={0,0,0},
          fillColor={64,64,64}, fillPattern=FillPattern.Solid, radius=6),
        Rectangle(extent={{30,72},{58,44}}, lineColor={0,0,0},
          fillColor={64,64,64}, fillPattern=FillPattern.Solid, radius=6),
        Rectangle(extent={{-14,-76},{14,-48}}, lineColor={0,0,0},
          fillColor={64,64,64}, fillPattern=FillPattern.Solid, radius=6),
        Text(extent={{-150,120},{150,90}}, textColor={0,0,255}, textString="%name")}),
      Documentation(info="<html>
<p>Planar three-wheel (&quot;tricycle&quot;) vehicle for studying tire forces on the
suspension links and steering arm: <b>individual front-left/right wheels</b> (own slip
angles, own vertical loads via quasi-static lateral load transfer) and a <b>lumped rear
wheel</b> &mdash; the same architecture used for the front-axle force estimation in
WO&nbsp;2025/113783 (Marzbanrad &amp; Jonasson). Chassis: 2 DOF (lateral velocity,
yaw rate) at constant forward speed, plus path kinematics. Tires: Pacejka brush model
with degressive-load cornering stiffness and relaxation-length lag.</p>
<p>The kingpin moment per side is
M<sub>kp</sub> = F<sub>y</sub>&middot;t<sub>mech</sub> &minus; M<sub>z</sub>
(mechanical + pneumatic trail); the tie-rod force is M<sub>kp</sub>/L<sub>arm</sub>,
and both sides sum into the rack flange together with the reflected road-wheel inertia.
<code>toeL</code>/<code>toeR</code> add per-wheel steer offsets &mdash; the hooks for an
active toe-control actuator (wire to 0 when unused).</p>
<p><b>Documented omissions</b> (negligible for a constant-speed double lane change at
road-wheel angles of a few degrees): Ackermann split (&lt;0.1&deg;), Fz-jacking from
KPI/caster (&lt;2 N&middot;m), scrub-radius&times;Fx (no drive/brake force), rear tire
load sensitivity (lumped axle, slightly understeer-optimistic; front link loads
unaffected), roll DOF (load transfer carries a first-order roll-mode lag instead).</p>
<p>Sources: Pacejka <i>Tire and Vehicle Dynamics</i> Ch. 3 &amp; 7; Heydinger et al.
SAE 1999-01-1336 (mass/inertia); Milliken &amp; Milliken RCVD (load transfer, roll
stiffness fraction); Reimpell (caster trail); Yin et al. 2024 (kingpin wheel inertia).</p>
</html>"));
  end PlanarTricycle;

  model ManualSteering
    "Unassisted rack-and-pinion steering: handwheel inertia + ideal pinion + rack mass"
    parameter Real iS = 20 "Overall steering ratio, handwheel angle / road-wheel angle";
    parameter Modelica.Units.SI.Length Larm = 0.11
      "Steering-arm length (must match the vehicle): rack travel per road-wheel radian";
    parameter Modelica.Units.SI.Inertia Jhw = 0.035
      "Handwheel + column inertia";
    parameter Modelica.Units.SI.Mass mRack = 3.0
      "Rack + tie-rod translating mass";
    final parameter Modelica.Units.SI.Length rPinion = Larm/iS
      "Effective pinion radius: rack travel per handwheel radian";

    Modelica.Mechanics.Rotational.Interfaces.Flange_a handwheel
      "Handwheel connection (driver side)"
      annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
    Modelica.Mechanics.Translational.Interfaces.Flange_b rack
      "Rack / tie-rod connection (vehicle side)"
      annotation (Placement(transformation(extent={{90,-10},{110,10}})));

    Modelica.Mechanics.Rotational.Components.Inertia hw(J=Jhw)
      annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
    Modelica.Mechanics.Rotational.Components.IdealGearR2T pinion(ratio=1/rPinion)
      "Rack and pinion: handwheel angle <-> rack travel"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
    Modelica.Mechanics.Translational.Components.Mass rackMass(m=mRack)
      annotation (Placement(transformation(extent={{40,-10},{60,10}})));
  equation
    connect(handwheel, hw.flange_a);
    connect(hw.flange_b, pinion.flangeR);
    connect(pinion.flangeT, rackMass.flange_a);
    connect(rackMass.flange_b, rack);
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Ellipse(extent={{-90,30},{-30,-30}}, lineColor={0,0,0}, lineThickness=0.5),
        Ellipse(extent={{-66,6},{-54,-6}}, lineColor={0,0,0}, fillColor={64,64,64},
          fillPattern=FillPattern.Solid),
        Rectangle(extent={{-30,6},{90,-6}}, lineColor={64,64,64},
          fillPattern=FillPattern.HorizontalCylinder, fillColor={192,192,192}),
        Text(extent={{-150,80},{150,50}}, textColor={0,0,255}, textString="%name"),
        Text(extent={{-150,-40},{150,-70}}, textColor={0,0,0}, textString="i=%iS")}),
      Documentation(info="<html><p>Traditional <b>unassisted</b> rack-and-pinion
steering: a handwheel/column inertia drives the rack through an ideal (lossless,
backlash-free) pinion of effective radius r<sub>p</sub> = L<sub>arm</sub>/i<sub>S</sub>,
so the kinematics are &delta; = &phi;<sub>HW</sub>/i<sub>S</sub>. With no assistance,
the full kingpin reaction reflects to the handwheel:
&tau;<sub>HW</sub> = F<sub>rack</sub>&middot;r<sub>p</sub> &mdash; the steering-feel
torque is a primary output. The default ratio i<sub>S</sub> = 20 is typical for manual
(unassisted) passenger-car steering (rack travel &asymp; 35 mm per handwheel rev);
column compliance and friction are omitted (document if they matter for your use).</p>
</html>"));
  end ManualSteering;

  package Examples "Vehicle validation and the ISO 3888-1 double lane change"
    extends Modelica.Icons.ExamplesPackage;

    model StepSteer
      "Validation: ideal-position steer step into the tricycle vehicle (no steering hardware) for understeer-gradient checks"
      extends Modelica.Icons.Example;
      parameter Modelica.Units.SI.Velocity u = 80/3.6 "Forward speed";
      parameter Real deltaStepDeg = 0.5 "Road-wheel steer step [deg]";
      parameter Modelica.Units.SI.Length Larm = 0.11 "Steering-arm length";

      output Real yawRateDegS = trike.r*180/pi "Yaw rate [deg/s]";
      output Real ayG = trike.ay/Modelica.Constants.g_n "Lateral acceleration [g]";
      output Real deltaLdeg = trike.dL*180/pi "Left road-wheel angle [deg]";
      output Real betaDeg = atan(trike.vy/u)*180/pi "Body sideslip [deg]";
      output Modelica.Units.SI.Force FtieL = trike.FtieL "Left tie-rod force [N]";
      output Modelica.Units.SI.Force FtieR = trike.FtieR "Right tie-rod force [N]";
      output Modelica.Units.SI.Force FzFL = trike.FzFL "Front-left vertical load [N]";
      output Modelica.Units.SI.Force FzFR = trike.FzFR "Front-right vertical load [N]";
      output Real alphaFLdeg = trike.aFL*180/pi "Front-left slip angle [deg]";
      output Modelica.Units.SI.Force FyFL = trike.FyFL "Front-left lateral force [N]";
      output Modelica.Units.SI.Torque MzFL = trike.MzFL "Front-left aligning moment [N.m]";
    protected
      constant Real pi = Modelica.Constants.pi;
    public
      Tricycle.PlanarTricycle trike(u=u, Larm=Larm)
        annotation (Placement(transformation(extent={{20,-10},{40,10}})));
      Modelica.Mechanics.Translational.Sources.Position pos(exact=false, f_crit=5)
        "Ideal (filtered) rack-position source - no steering hardware"
        annotation (Placement(transformation(extent={{-20,-10},{0,10}})));
      Modelica.Blocks.Sources.Step sRef(
        height=deltaStepDeg*Modelica.Constants.pi/180*Larm, startTime=1)
        annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
      Modelica.Blocks.Sources.Constant toe0L(k=0)
        annotation (Placement(transformation(extent={{-20,30},{0,50}})));
      Modelica.Blocks.Sources.Constant toe0R(k=0)
        annotation (Placement(transformation(extent={{-20,-50},{0,-30}})));
    equation
      connect(sRef.y, pos.s_ref);
      connect(pos.flange, trike.rack);
      connect(toe0L.y, trike.toeL);
      connect(toe0R.y, trike.toeR);
      annotation (experiment(StopTime=6, Interval=0.002, Tolerance=1e-6),
        Documentation(info="<html><p>Steer-step response of the
<a href=\"modelica://Tricycle.PlanarTricycle\">PlanarTricycle</a> with an ideal
position source, so vehicle dynamics can be validated in isolation: sweeping
<code>u</code> and reading the steady-state yaw-rate gain checks the understeer
gradient against the analytic bicycle model.</p></html>"));
    end StepSteer;

    model DoubleLaneChange
      "ISO 3888-1 double lane change: preview driver turns the handwheel of the manual rack, steering the tricycle vehicle"
      extends Modelica.Icons.Example;
      parameter Modelica.Units.SI.Velocity u = 80/3.6 "Vehicle speed";
      parameter Modelica.Units.SI.Time Tp = 0.55
        "Driver preview time (compensates the arm lag; shorter values destabilize the loop)";
      parameter Real Kdrv(unit="rad/m") = 0.22 "Driver preview gain";
      parameter Real Kr(unit="rad.s/rad") = 0.25 "Driver yaw-rate damping";
      parameter Real iS = 20 "Overall steering ratio";
      parameter Modelica.Units.SI.Length Larm = 0.11 "Steering-arm length";
      parameter Modelica.Units.SI.Frequency fArm = 2
        "Driver arm/neuromuscular bandwidth (position-tracking filter)";

      output Modelica.Units.SI.Force FtieL = trike.FtieL "Left tie-rod force [N]";
      output Modelica.Units.SI.Force FtieR = trike.FtieR "Right tie-rod force [N]";
      output Modelica.Units.SI.Torque MkpL = trike.MkpL "Left kingpin moment [N.m]";
      output Modelica.Units.SI.Torque MkpR = trike.MkpR "Right kingpin moment [N.m]";
      output Modelica.Units.SI.Torque MkpMechL = trike.MkpMechL "Mechanical-trail part, left";
      output Modelica.Units.SI.Torque MkpMechR = trike.MkpMechR "Mechanical-trail part, right";
      output Modelica.Units.SI.Torque MkpPneuL = trike.MkpPneuL "Pneumatic-trail part, left";
      output Modelica.Units.SI.Torque MkpPneuR = trike.MkpPneuR "Pneumatic-trail part, right";
      output Modelica.Units.SI.Force rackForce = trike.rack.f "Rack axial reaction [N]";
      output Modelica.Units.SI.Force FyFL = trike.FyFL "Front-left lateral force [N]";
      output Modelica.Units.SI.Force FyFR = trike.FyFR "Front-right lateral force [N]";
      output Modelica.Units.SI.Force FyR = trike.FyR "Rear (lumped) lateral force [N]";
      output Modelica.Units.SI.Force FzFL = trike.FzFL "Front-left vertical load [N]";
      output Modelica.Units.SI.Force FzFR = trike.FzFR "Front-right vertical load [N]";
      output Modelica.Units.SI.Torque MzFL = trike.MzFL "Front-left aligning moment [N.m]";
      output Modelica.Units.SI.Torque MzFR = trike.MzFR "Front-right aligning moment [N.m]";
      output Real alphaFLdeg = trike.aFL*180/pi "Front-left slip angle [deg]";
      output Real alphaFRdeg = trike.aFR*180/pi "Front-right slip angle [deg]";
      output Real alphaRdeg = trike.aR*180/pi "Rear slip angle [deg]";
      output Real yawRateDegS = trike.r*180/pi "Yaw rate [deg/s]";
      output Real ayG = trike.ay/Modelica.Constants.g_n "Lateral acceleration [g]";
      output Real betaDeg = atan(trike.vy/u)*180/pi "Body sideslip angle [deg]";
      output Modelica.Units.SI.Position X = trike.X "Global x [m]";
      output Modelica.Units.SI.Position Y = trike.Y "Global y [m]";
      output Real psiDeg = trike.psi*180/pi "Heading [deg]";
      output Real deltaLdeg = trike.dL*180/pi "Left road-wheel angle [deg]";
      output Real deltaRdeg = trike.dR*180/pi "Right road-wheel angle [deg]";
      output Real dCmdDeg = driver.dCmd*180/pi "Driver steer command [deg]";
      output Real hwaDeg = steering.hw.phi*180/pi "Handwheel angle [deg]";
      output Modelica.Units.SI.Torque hwTorque = hwTorqueSensor.tau
        "Handwheel (steering-feel) torque [N.m]";
      output Modelica.Units.SI.Position rackDisp = steering.rackMass.s "Rack displacement [m]";
    protected
      constant Real pi = Modelica.Constants.pi;
    public
      Tricycle.Iso3888Path driver(u=u, Tp=Tp, Kdrv=Kdrv, Kr=Kr)
        annotation (Placement(transformation(extent={{-120,40},{-100,60}})));
      Modelica.Blocks.Math.Gain refGain(k=iS)
        "Road-wheel steer command [rad] -> handwheel angle reference [rad]"
        annotation (Placement(transformation(extent={{-86,44},{-74,56}})));
      Modelica.Mechanics.Rotational.Sources.Position arm(exact=false, f_crit=fArm)
        "Driver arm: filtered handwheel-position tracking"
        annotation (Placement(transformation(extent={{-60,40},{-40,60}})));
      Modelica.Mechanics.Rotational.Sensors.TorqueSensor hwTorqueSensor
        annotation (Placement(transformation(extent={{-30,40},{-10,60}})));
      Tricycle.ManualSteering steering(iS=iS, Larm=Larm)
        annotation (Placement(transformation(extent={{0,40},{20,60}})));
      Tricycle.PlanarTricycle trike(u=u, Larm=Larm)
        annotation (Placement(transformation(extent={{40,40},{60,60}})));
      Modelica.Blocks.Sources.RealExpression Xe(y=trike.X)
        annotation (Placement(transformation(extent={{-156,54},{-136,68}})));
      Modelica.Blocks.Sources.RealExpression Ye(y=trike.Y)
        annotation (Placement(transformation(extent={{-156,42},{-136,56}})));
      Modelica.Blocks.Sources.RealExpression psie(y=trike.psi)
        annotation (Placement(transformation(extent={{-156,30},{-136,44}})));
      Modelica.Blocks.Sources.RealExpression vye(y=trike.vy)
        annotation (Placement(transformation(extent={{-156,18},{-136,32}})));
      Modelica.Blocks.Sources.RealExpression re(y=trike.r)
        annotation (Placement(transformation(extent={{-156,6},{-136,20}})));
      Modelica.Blocks.Sources.Constant toe0L(k=0)
        annotation (Placement(transformation(extent={{20,70},{40,90}})));
      Modelica.Blocks.Sources.Constant toe0R(k=0)
        annotation (Placement(transformation(extent={{20,10},{40,30}})));
    equation
      connect(Xe.y, driver.X);
      connect(Ye.y, driver.Y);
      connect(psie.y, driver.psi);
      connect(vye.y, driver.vy);
      connect(re.y, driver.r);
      connect(driver.dCmd, refGain.u);
      connect(refGain.y, arm.phi_ref);
      connect(arm.flange, hwTorqueSensor.flange_a);
      connect(hwTorqueSensor.flange_b, steering.handwheel);
      connect(steering.rack, trike.rack);
      connect(toe0L.y, trike.toeL);
      connect(toe0R.y, trike.toeR);
      annotation (experiment(StopTime=7, Interval=0.001, Tolerance=1e-6),
        Documentation(info="<html><p>ISO 3888-1 double lane change at constant speed:
the preview driver tracks the reference centerline and commands a road-wheel angle,
converted to a handwheel reference (ratio i<sub>S</sub>) and applied by a filtered
position source representing the driver's arm. The <b>unassisted</b> manual rack
reflects the full kingpin reaction back to the handwheel, so the steering-feel torque
<code>hwTorque</code> is a headline output alongside the tie-rod forces.</p></html>"));
    end DoubleLaneChange;

    model OpenLoopDLC
      "Open-loop lane-change-shaped steer (one-period sine) through the manual rack - repeatable sweeps"
      extends Modelica.Icons.Example;
      parameter Modelica.Units.SI.Velocity u = 80/3.6 "Vehicle speed";
      parameter Real iS = 20 "Overall steering ratio";
      parameter Modelica.Units.SI.Length Larm = 0.11 "Steering-arm length";
      parameter Modelica.Units.SI.Frequency fArm = 2 "Driver arm bandwidth";
      parameter Real ampDeg = 3 "Road-wheel steer amplitude [deg]";
      parameter Modelica.Units.SI.Time period = 2 "Steer period [s]";
      parameter Modelica.Units.SI.Time t0 = 1 "Steer start time [s]";

      output Modelica.Units.SI.Force FtieL = trike.FtieL "Left tie-rod force [N]";
      output Modelica.Units.SI.Force FtieR = trike.FtieR "Right tie-rod force [N]";
      output Modelica.Units.SI.Torque MkpL = trike.MkpL "Left kingpin moment [N.m]";
      output Modelica.Units.SI.Torque MkpR = trike.MkpR "Right kingpin moment [N.m]";
      output Modelica.Units.SI.Force rackForce = trike.rack.f "Rack axial reaction [N]";
      output Real ayG = trike.ay/Modelica.Constants.g_n "Lateral acceleration [g]";
      output Real yawRateDegS = trike.r*180/pi "Yaw rate [deg/s]";
      output Modelica.Units.SI.Position Y = trike.Y "Global y [m]";
      output Real deltaLdeg = trike.dL*180/pi "Left road-wheel angle [deg]";
      output Modelica.Units.SI.Torque hwTorque = hwTorqueSensor.tau "Handwheel torque [N.m]";
    protected
      constant Real pi = Modelica.Constants.pi;
    public
      Modelica.Blocks.Sources.RealExpression steerCmd(
        y=if time > t0 and time < t0 + period then
            iS*ampDeg*Modelica.Constants.pi/180*sin(2*Modelica.Constants.pi*(time - t0)/period)
          else 0.0) "One-period-sine handwheel angle command [rad]"
        annotation (Placement(transformation(extent={{-120,42},{-100,58}})));
      Modelica.Mechanics.Rotational.Sources.Position arm(exact=false, f_crit=fArm)
        annotation (Placement(transformation(extent={{-60,40},{-40,60}})));
      Modelica.Mechanics.Rotational.Sensors.TorqueSensor hwTorqueSensor
        annotation (Placement(transformation(extent={{-30,40},{-10,60}})));
      Tricycle.ManualSteering steering(iS=iS, Larm=Larm)
        annotation (Placement(transformation(extent={{0,40},{20,60}})));
      Tricycle.PlanarTricycle trike(u=u, Larm=Larm)
        annotation (Placement(transformation(extent={{40,40},{60,60}})));
      Modelica.Blocks.Sources.Constant toe0L(k=0)
        annotation (Placement(transformation(extent={{20,70},{40,90}})));
      Modelica.Blocks.Sources.Constant toe0R(k=0)
        annotation (Placement(transformation(extent={{20,10},{40,30}})));
    equation
      connect(steerCmd.y, arm.phi_ref);
      connect(arm.flange, hwTorqueSensor.flange_a);
      connect(hwTorqueSensor.flange_b, steering.handwheel);
      connect(steering.rack, trike.rack);
      connect(toe0L.y, trike.toeL);
      connect(toe0R.y, trike.toeR);
      annotation (experiment(StopTime=5, Interval=0.001, Tolerance=1e-6),
        Documentation(info="<html><p>Same steering and vehicle as
<a href=\"modelica://Tricycle.Examples.DoubleLaneChange\">DoubleLaneChange</a> but with
a prescribed one-period-sine handwheel command instead of the closed-loop driver:
perfectly repeatable for amplitude/frequency sweeps.</p></html>"));
    end OpenLoopDLC;

    model NordschleifeLap
      "Minimum-manageable-time lap of the (planar) Nurburgring Nordschleife: power- and grip-limited tricycle, preview driver with speed control"
      extends Modelica.Icons.Example;
      parameter String fileName = "build/ns_track.txt"
        "Track table (generate with track_lap.py: s, kappa, vRef, axFF)";
      parameter Modelica.Units.SI.Length sLap = 20718.5 "Lap length (terminate there)";
      parameter Modelica.Units.SI.Velocity u0 = 30 "Rolling-start speed";
      parameter Modelica.Units.SI.Power Pmax = 150e3 "Peak drive power";

      output Modelica.Units.SI.Position s = trike.s "Distance along centerline [m]";
      output Modelica.Units.SI.Length n = trike.n "Lateral offset from centerline [m]";
      output Real vKmh = trike.u*3.6 "Speed [km/h]";
      output Real vRefKmh = driver.vRefPrev*3.6 "Previewed speed reference [km/h]";
      output Real ayG = trike.ay/Modelica.Constants.g_n "Lateral acceleration [g]";
      output Real axG = trike.ax/Modelica.Constants.g_n "Longitudinal acceleration [g]";
      output Real deltaDeg = trike.delta*180/Modelica.Constants.pi "Road-wheel angle [deg]";
      output Real dpsiDeg = trike.dpsi*180/Modelica.Constants.pi "Heading error [deg]";
      output Real yawRateDegS = trike.r*180/Modelica.Constants.pi "Yaw rate [deg/s]";
      output Real betaDeg = atan(trike.vy/max(trike.u, 1))*180/Modelica.Constants.pi
        "Body sideslip [deg]";
      output Modelica.Units.SI.Force FtieL = trike.FtieL "Left tie-rod force [N]";
      output Modelica.Units.SI.Force FtieR = trike.FtieR "Right tie-rod force [N]";
      output Modelica.Units.SI.Force FxR = trike.FxR "Rear drive/brake force [N]";
      output Modelica.Units.SI.Power Pdrive = trike.FxR*trike.u "Drive power [W]";
      output Modelica.Units.SI.Force FzFL = trike.FzFL "Front-left vertical load [N]";
      output Modelica.Units.SI.Force FzFR = trike.FzFR "Front-right vertical load [N]";

      Tricycle.Track.TrackTricycle trike(fileName=fileName, u0=u0, Pmax=Pmax)
        annotation (Placement(transformation(extent={{20,-10},{40,10}})));
      Tricycle.Track.TrackDriver driver(fileName=fileName)
        annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));
      Modelica.Blocks.Sources.RealExpression se(y=trike.s)
        annotation (Placement(transformation(extent={{-80,22},{-60,36}})));
      Modelica.Blocks.Sources.RealExpression ne(y=trike.n)
        annotation (Placement(transformation(extent={{-80,8},{-60,22}})));
      Modelica.Blocks.Sources.RealExpression dpsie(y=trike.dpsi)
        annotation (Placement(transformation(extent={{-80,-6},{-60,8}})));
      Modelica.Blocks.Sources.RealExpression ue(y=trike.u)
        annotation (Placement(transformation(extent={{-80,-20},{-60,-6}})));
      Modelica.Blocks.Sources.RealExpression vye(y=trike.vy)
        annotation (Placement(transformation(extent={{-80,-34},{-60,-20}})));
      Modelica.Blocks.Sources.RealExpression re(y=trike.r)
        annotation (Placement(transformation(extent={{-80,-48},{-60,-34}})));
    equation
      connect(se.y, driver.s);
      connect(ne.y, driver.n);
      connect(dpsie.y, driver.dpsi);
      connect(ue.y, driver.vx);
      connect(vye.y, driver.vy);
      connect(re.y, driver.r);
      connect(driver.dCmd, trike.dCmd);
      connect(driver.axCmd, trike.axCmd);
      when trike.s >= sLap then
        terminate("lap complete");
      end when;
      annotation (experiment(StopTime=900, Interval=0.05, Tolerance=1e-6),
        Documentation(info="<html><p>One flying lap of the planar Nordschleife
centerline (OSM data, ~20.7 km). The driver tracks the quasi-steady minimum-time
speed profile computed for exactly this vehicle setup (mass, power, grip, drag), so
the lap time is the minimum <i>he</i> can manage on the centerline. The simulation
terminates when s reaches <code>sLap</code>; the lap time is the final simulation
time plus the small correction for the rolling start. Run via
<code>track_lap.py</code>, which generates the track/speed-profile table first.</p></html>"));
    end NordschleifeLap;
  end Examples;

  package Track
    "Curvilinear (Frenet) track following: power- and traction-limited tricycle driven around a closed circuit at the driver's minimum manageable time"
    extends Modelica.Icons.Package;

    function smoothMin "C2 smooth minimum, 0.5*(x+y-sqrt((x-y)^2+eps^2))"
      extends Modelica.Icons.Function;
      input Real x;
      input Real y;
      input Real eps = 1 "Blend width (units of x,y)";
      output Real z;
    algorithm
      z := 0.5*(x + y - sqrt((x - y)^2 + eps^2));
      annotation (smoothOrder=2);
    end smoothMin;

    function smoothMax "C2 smooth maximum, 0.5*(x+y+sqrt((x-y)^2+eps^2))"
      extends Modelica.Icons.Function;
      input Real x;
      input Real y;
      input Real eps = 1 "Blend width (units of x,y)";
      output Real z;
    algorithm
      z := 0.5*(x + y + sqrt((x - y)^2 + eps^2));
      annotation (smoothOrder=2);
    end smoothMax;

    model TrackTricycle
      "Tricycle in track coordinates (s,n,dpsi) with a longitudinal DOF: RWD, power- and friction-ellipse-limited drive/brake, brush tires, tie-rod loads"
      // Track coordinates: s along the centerline, n leftward lateral offset,
      // dpsi = vehicle heading - centerline heading. kappa > 0 = left curve.
      parameter String fileName = "build/ns_track.txt"
        "Track table file: columns s [m], kappa [1/m], vRef [m/s], axFF [m/s2]";
      parameter Modelica.Units.SI.Velocity u0 = 30 "Initial (rolling-start) speed";
      // chassis (identical defaults to PlanarTricycle)
      parameter Modelica.Units.SI.Mass m = 1650 "Vehicle mass";
      parameter Modelica.Units.SI.Inertia Izz = 2700 "Yaw moment of inertia";
      parameter Modelica.Units.SI.Length a = 1.20 "CG to front axle";
      parameter Modelica.Units.SI.Length b = 1.60 "CG to rear axle";
      parameter Modelica.Units.SI.Length tf = 1.55 "Front track width";
      parameter Modelica.Units.SI.Length hcg = 0.55 "CG height above ground";
      parameter Real xiF(min=0, max=1) = 0.6 "Front share of lateral load transfer";
      parameter Modelica.Units.SI.Time tauRoll = 0.15 "Roll-mode lag on lateral transfer";
      parameter Modelica.Units.SI.Time tauPitch = 0.20 "Pitch-mode lag on longitudinal transfer";
      parameter Modelica.Units.SI.Velocity uMin = 1 "Low-speed guard";
      parameter Modelica.Units.SI.Length tMech = 0.025 "Mechanical (caster) trail";
      parameter Modelica.Units.SI.Length Larm = 0.11 "Steering-arm length";
      parameter Modelica.Units.SI.Time tauSteer = 0.12
        "Steer actuation lag (driver arm + linkage), road-wheel level";
      // powertrain / resistances
      parameter Modelica.Units.SI.Power Pmax = 150e3 "Peak drive power at the rear wheel";
      parameter Real kBf(min=0, max=1) = 0.65 "Front share of the brake force";
      parameter Real CdA(unit="m2") = 0.72 "Drag area Cd*A";
      parameter Real Crr = 0.012 "Rolling-resistance coefficient";
      parameter Modelica.Units.SI.Density rho = 1.20 "Air density";
      parameter TireData tireF "Front tire, per wheel";
      parameter TireData tireR(c1=1.4e5, c2=8000, FzNom=8000, ap0=0.085)
        "Lumped rear-axle tire";

      Modelica.Blocks.Interfaces.RealInput dCmd(unit="rad") "Road-wheel steer command"
        annotation (Placement(transformation(extent={{-120,30},{-100,50}})));
      Modelica.Blocks.Interfaces.RealInput axCmd(unit="m/s2")
        "Longitudinal acceleration request (drive > 0, brake < 0)"
        annotation (Placement(transformation(extent={{-120,-50},{-100,-30}})));

      // track-coordinate states
      Modelica.Units.SI.Position s(start=0, fixed=true) "Distance along centerline";
      Modelica.Units.SI.Length n(start=0, fixed=true) "Lateral offset, left of centerline";
      Modelica.Units.SI.Angle dpsi(start=0, fixed=true) "Heading error to centerline";
      Modelica.Units.SI.Velocity u(start=u0, fixed=true) "Forward speed (state)";
      Modelica.Units.SI.Velocity vy(start=0, fixed=true) "Lateral velocity at CG";
      Modelica.Units.SI.AngularVelocity r(start=0, fixed=true) "Yaw rate";
      Modelica.Units.SI.Angle delta(start=0, fixed=true) "Road-wheel steer angle";
      Modelica.Units.SI.Force dFz(start=0, fixed=true) "Front lateral load transfer";
      Modelica.Units.SI.Force dFzX(start=0, fixed=true) "Longitudinal load transfer";
      Modelica.Units.SI.Angle aLagFL(start=0, fixed=true) "Front-left lagged slip angle";
      Modelica.Units.SI.Angle aLagFR(start=0, fixed=true) "Front-right lagged slip angle";
      Modelica.Units.SI.Angle aLagR(start=0, fixed=true) "Rear lagged slip angle";
      // kinematics, tires, loads
      Real kappa(unit="1/m") "Centerline curvature at s";
      Modelica.Units.SI.Velocity sdot "Centerline progress rate";
      Modelica.Units.SI.Angle dL, dR;
      Modelica.Units.SI.Angle aFL, aFR, aR;
      Modelica.Units.SI.Force FzFL, FzFR, FzR;
      Modelica.Units.SI.Force FyFL, FyFR, FyR;
      Modelica.Units.SI.Torque MzFL, MzFR, MzR;
      Modelica.Units.SI.Acceleration ax "Longitudinal acceleration (= der(u) - vy*r)";
      Modelica.Units.SI.Acceleration ay "Lateral acceleration (= der(vy) + u*r)";
      Modelica.Units.SI.Force FyF "Front-axle lateral force in the body frame";
      // longitudinal forces
      Modelica.Units.SI.Force Fdrag, Freq;
      Modelica.Units.SI.Force FxR "Rear drive/brake force (applied)";
      Modelica.Units.SI.Force FxFL, FxFR "Front brake forces (applied, <= 0)";
      Modelica.Units.SI.Force FxRMax, FxFLMax, FxFRMax "Friction-ellipse remainders";
      Modelica.Units.SI.Force FxbFL, FxbFR, FybFL, FybFR "Front wheel forces, body frame";
      // kingpin / tie-rod (same decomposition as PlanarTricycle)
      Modelica.Units.SI.Torque MkpL, MkpR;
      Modelica.Units.SI.Force FtieL, FtieR;
    protected
      constant Modelica.Units.SI.Acceleration g = Modelica.Constants.g_n;
      final parameter Modelica.Units.SI.Length L = a + b;
      final parameter Modelica.Units.SI.Force Fz0F = m*g*b/(2*L);
      final parameter Modelica.Units.SI.Force Fz0R = m*g*a/L;
      final parameter Modelica.Units.SI.Force Froll = Crr*m*g;
      Modelica.Units.SI.Velocity uG "Guarded speed";
      Modelica.Units.SI.Force FreqPos, FreqNeg "Drive / brake split of the request";
      Modelica.Blocks.Tables.CombiTable1Ds kappaT(
        tableOnFile=true, tableName="track", fileName=fileName, columns={2},
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative)
        "Centerline curvature lookup";
    equation
      uG = noEvent(max(u, uMin));
      // track (Frenet) kinematics
      kappaT.u = s;
      kappa = kappaT.y[1];
      sdot = (u*cos(dpsi) - vy*sin(dpsi))/noEvent(max(1 - n*kappa, 0.2));
      der(s) = sdot;
      der(n) = u*sin(dpsi) + vy*cos(dpsi);
      der(dpsi) = r - kappa*sdot;
      // steering: first-order actuation lag, parallel steer
      tauSteer*der(delta) + delta = dCmd;
      dL = delta;
      dR = delta;
      aFL = dL - atan((vy + a*r)/noEvent(max(u - r*tf/2, uMin)));
      aFR = dR - atan((vy + a*r)/noEvent(max(u + r*tf/2, uMin)));
      aR  =    - atan((vy - b*r)/uG);
      // quasi-static load transfer, roll- and pitch-lagged
      tauRoll*der(dFz) + dFz = xiF*m*ay*hcg/tf;
      tauPitch*der(dFzX) + dFzX = m*ax*hcg/L;
      FzFL = Fz0F - dFz - dFzX/2;
      FzFR = Fz0F + dFz - dFzX/2;
      FzR  = Fz0R + dFzX;
      // brush tires at relaxation-lagged slip angles
      (tireF.sigmaRel/uG)*der(aLagFL) + aLagFL = aFL;
      (tireF.sigmaRel/uG)*der(aLagFR) + aLagFR = aFR;
      (tireR.sigmaRel/uG)*der(aLagR)  + aLagR  = aR;
      (FyFL, MzFL) = brushForces(aLagFL, FzFL, tireF);
      (FyFR, MzFR) = brushForces(aLagFR, FzFR, tireF);
      (FyR,  MzR)  = brushForces(aLagR,  FzR,  tireR);
      // longitudinal request -> drive (rear) / brake (split), limited by engine power
      // and by the friction-ellipse remainder per axle ("ideal TC/ABS"; no wheel-spin
      // model). Lateral force is not reduced by Fx - documented omission.
      Fdrag = 0.5*rho*CdA*u^2;
      Freq = m*axCmd + Fdrag + Froll;
      FreqPos = 0.5*(Freq + sqrt(Freq^2 + 100^2));
      FreqNeg = Freq - FreqPos;
      FxRMax  = tireR.mu*FzR*sqrt(smoothMax(1 - (FyR/(tireR.mu*FzR))^2, 0, 0.02));
      FxFLMax = tireF.mu*FzFL*sqrt(smoothMax(1 - (FyFL/(tireF.mu*FzFL))^2, 0, 0.02));
      FxFRMax = tireF.mu*FzFR*sqrt(smoothMax(1 - (FyFR/(tireF.mu*FzFR))^2, 0, 0.02));
      FxR  = smoothMin(smoothMin(FreqPos, Pmax/uG, 100), FxRMax, 100)
           + smoothMax((1 - kBf)*FreqNeg, -FxRMax, 100);
      FxFL = smoothMax(0.5*kBf*FreqNeg, -FxFLMax, 100);
      FxFR = smoothMax(0.5*kBf*FreqNeg, -FxFRMax, 100);
      // chassis: longitudinal + lateral + yaw in the body frame
      FxbFL = FxFL*cos(dL) - FyFL*sin(dL);
      FxbFR = FxFR*cos(dR) - FyFR*sin(dR);
      FybFL = FyFL*cos(dL) + FxFL*sin(dL);
      FybFR = FyFR*cos(dR) + FxFR*sin(dR);
      FyF = FybFL + FybFR;
      m*ax = FxbFL + FxbFR + FxR - Fdrag - Froll;
      m*ay = FyF + FyR;
      der(u)  = ax + vy*r;
      der(vy) = ay - u*r;
      Izz*der(r) = a*FyF - b*FyR - (tf/2)*(FxbFL - FxbFR)
                 + MzFL + MzFR + MzR;
      // kingpin moments -> tie-rod forces (scrub-radius x Fx omitted, as before)
      MkpL = FyFL*tMech - MzFL;
      MkpR = FyFR*tMech - MzFR;
      FtieL = MkpL/Larm;
      FtieR = MkpR/Larm;
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(extent={{-60,80},{60,-80}}, lineColor={0,0,0},
            fillColor={235,245,235}, fillPattern=FillPattern.Solid, radius=20),
          Rectangle(extent={{-58,72},{-30,44}}, lineColor={0,0,0},
            fillColor={64,64,64}, fillPattern=FillPattern.Solid, radius=6),
          Rectangle(extent={{30,72},{58,44}}, lineColor={0,0,0},
            fillColor={64,64,64}, fillPattern=FillPattern.Solid, radius=6),
          Rectangle(extent={{-14,-76},{14,-48}}, lineColor={0,0,0},
            fillColor={64,64,64}, fillPattern=FillPattern.Solid, radius=6),
          Text(extent={{-150,120},{150,90}}, textColor={0,0,255}, textString="%name")}),
        Documentation(info="<html>
<p>The planar tricycle re-expressed in <b>track (Frenet) coordinates</b>
(s,&nbsp;n,&nbsp;&Delta;&psi;) around a closed centerline given as a table
&kappa;(s), with a <b>longitudinal degree of freedom</b>: rear-wheel drive limited by
engine power P<sub>max</sub>/u and by the rear friction-ellipse remainder
&radic;((&mu;F<sub>z</sub>)&sup2;&nbsp;&minus;&nbsp;F<sub>y</sub>&sup2;), brakes split
front/rear and limited the same way (&quot;ideal TC/ABS&quot;), plus aerodynamic drag
and rolling resistance. Lateral tire behavior, load-transfer lags, and the
kingpin/tie-rod decomposition are carried over from
<a href=\"modelica://Tricycle.PlanarTricycle\">PlanarTricycle</a>; longitudinal load
transfer (pitch-lagged) is added because braking/traction depend on it.</p>
<p><b>Documented omissions</b>: no combined-slip reduction of F<sub>y</sub> by
F<sub>x</sub> (grip-optimistic at corner exit/entry), no wheel-spin/lock dynamics,
no elevation or banking (planar by design), steering hardware replaced by a
first-order road-wheel lag.</p></html>"));
    end TrackTricycle;

    block TrackDriver
      "Preview driver for a closed track: curvature feedforward + offset feedback steering, and speed control tracking a precomputed minimum-time profile vRef(s)"
      parameter String fileName = "build/ns_track.txt"
        "Track table file: columns s, kappa, vRef, axFF";
      parameter Modelica.Units.SI.Length Lwb = 2.80 "Wheelbase (steer feedforward)";
      parameter Real Kus(unit="rad.s2/m") = 1.63e-3
        "Understeer gradient for the feedforward (analytic value of the tricycle)";
      parameter Modelica.Units.SI.Time TpS = 0.45 "Steering preview time";
      parameter Modelica.Units.SI.Time TpV = 1.0 "Speed preview time";
      parameter Modelica.Units.SI.Time TpN = 0.8 "Lateral-offset prediction time";
      parameter Real Kn(unit="rad/m") = 0.25
        "Steer per meter of predicted offset at low speed";
      parameter Modelica.Units.SI.Velocity vKn = 30
        "Gain-scheduling speed: KnEff = Kn/(1+(vx/vKn)^2), the path response gain grows ~vx^2";
      parameter Real Kpsi = 0.5 "Steer per rad of heading error";
      parameter Real Kr(unit="rad.s/rad") = 0.3
        "Damping on the yaw-rate error r - kappa*vx (weave damping that does not fight steady cornering)";
      parameter Modelica.Units.SI.Angle dMax = 0.35 "Smooth steer clamp";
      parameter Modelica.Units.SI.Velocity uPrevMin = 5 "Preview distance floor";

      Modelica.Blocks.Interfaces.RealInput s(unit="m") "Distance along centerline"
        annotation (Placement(transformation(extent={{-120,70},{-100,90}})));
      Modelica.Blocks.Interfaces.RealInput n(unit="m") "Lateral offset"
        annotation (Placement(transformation(extent={{-120,30},{-100,50}})));
      Modelica.Blocks.Interfaces.RealInput dpsi(unit="rad") "Heading error"
        annotation (Placement(transformation(extent={{-120,-10},{-100,10}})));
      Modelica.Blocks.Interfaces.RealInput vx(unit="m/s") "Forward speed"
        annotation (Placement(transformation(extent={{-120,-50},{-100,-30}})));
      Modelica.Blocks.Interfaces.RealInput vy(unit="m/s") "Lateral velocity"
        annotation (Placement(transformation(extent={{-120,-90},{-100,-70}})));
      Modelica.Blocks.Interfaces.RealInput r(unit="rad/s") "Yaw rate"
        annotation (Placement(transformation(extent={{-120,-120},{-100,-100}})));
      Modelica.Blocks.Interfaces.RealOutput dCmd(unit="rad") "Road-wheel steer command"
        annotation (Placement(transformation(extent={{100,30},{120,50}})));
      Modelica.Blocks.Interfaces.RealOutput axCmd(unit="m/s2") "Acceleration request"
        annotation (Placement(transformation(extent={{100,-50},{120,-30}})));
      Real kapPrev(unit="1/m") "Previewed curvature";
      Real vRefPrev "Previewed speed reference [m/s] (unit-free: shares a table output vector)";
      Real axFFPrev "Previewed acceleration feedforward [m/s2] (logged; axCmd uses the preview-consistent law)";
      Modelica.Units.SI.Length nPrev "Predicted lateral offset";
      Modelica.Units.SI.Length dPrev "Speed preview distance";
    protected
      Modelica.Blocks.Tables.CombiTable1Ds steerT(
        tableOnFile=true, tableName="track", fileName=fileName, columns={2},
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative);
      Modelica.Blocks.Tables.CombiTable1Ds speedT(
        tableOnFile=true, tableName="track", fileName=fileName, columns={3,4},
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative);
      Modelica.Blocks.Tables.CombiTable1Ds nowT(
        tableOnFile=true, tableName="track", fileName=fileName, columns={2},
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative)
        "Curvature at the current station (yaw-rate reference)";
    equation
      dPrev = noEvent(max(vx, uPrevMin))*TpV;
      steerT.u = s + noEvent(max(vx, uPrevMin))*TpS;
      speedT.u = s + dPrev;
      kapPrev = steerT.y[1];
      vRefPrev = speedT.y[1];
      axFFPrev = speedT.y[2];
      nowT.u = s;
      nPrev = n + TpN*(vx*sin(dpsi) + vy*cos(dpsi));
      dCmd = dMax*tanh(((Lwb + Kus*vx^2)*kapPrev
                        - (Kn/(1 + (vx/vKn)^2))*nPrev - Kpsi*dpsi
                        - Kr*(r - nowT.y[1]*vx))/dMax);
      // constant-acceleration law to meet the previewed reference: self-correcting
      // (equals the profile's own feedforward when exactly on the profile)
      axCmd = (vRefPrev^2 - vx^2)/(2*dPrev);
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(extent={{-100,100},{100,-100}}, lineColor={0,0,0},
            fillColor={245,235,235}, fillPattern=FillPattern.Solid),
          Ellipse(extent={{-40,40},{40,-40}}, lineColor={0,0,0}, lineThickness=0.5),
          Line(points={{0,0},{22,32}}, color={0,0,0}, thickness=0.5),
          Text(extent={{-96,-60},{96,-90}}, textColor={0,0,0}, textString="min-time driver"),
          Text(extent={{-150,140},{150,110}}, textColor={0,0,255}, textString="%name")}),
        Documentation(info="<html><p>Two-channel driver: <b>steering</b> is
kinematic-plus-understeer curvature feedforward at a previewed station plus feedback
on the predicted lateral offset and heading error; <b>speed</b> tracks a
quasi-steady minimum-time profile v<sub>ref</sub>(s) (computed offline by
forward/backward passes under the same power, grip, and drag limits as the vehicle)
with its acceleration feedforward and a proportional speed error term. The lap time
this driver achieves is &quot;the minimum he can manage&quot; for the given setup
&mdash; not a formal optimum: he stays on the centerline instead of using the full
track width.</p></html>"));
    end TrackDriver;
  end Track;

  annotation (
    version="0.1.0",
    Documentation(info="<html>
<p>Minimal, defensible planar vehicle model for studying tire forces on the suspension
links and steering arm: a <b>three-wheel (tricycle)</b> chassis (individual front
wheels, lumped rear &mdash; the architecture of WO&nbsp;2025/113783), Pacejka brush
tires with collapsing pneumatic trail, kingpin trail decomposition, and a traditional
<b>unassisted rack-and-pinion</b> steering driven by a preview driver through the
ISO&nbsp;3888-1 double lane change. See <code>sources/SOURCES.md</code> for the
provenance of every parameter.</p>
</html>"));
end Tricycle;
