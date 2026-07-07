"""Chase-camera lap viewer: package the NordschleifeLap simulation results into a
self-contained HTML canvas player (outputs/ns_chase.html) - GTA-style follow camera
at road level zoom, a north-up minimap, and a telemetry HUD (speed, yaw rate,
lateral/longitudinal acceleration, g-g dot, steer and pedal bars).

Run AFTER track_lap.py (which leaves modelica/build/nslap_res.csv behind):
    python3 track_lap.py && python3 track_render.py
"""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'tracks'))
from speed_profile import load_centerline

RES = os.path.join(HERE, 'modelica/build/nslap_res.csv')
OUT = os.path.join(HERE, 'outputs/ns_chase.html')
FPS = 10          # embedded telemetry rate; the player interpolates
PMAX = 150e3      # must mirror the run (throttle bar = Pdrive/Pmax)

if not os.path.exists(RES):
    sys.exit('no simulation result - run track_lap.py first')

s, xc, yc, psic, kap = load_centerline(os.path.join(HERE, 'tracks/nordschleife.csv'))
ds = s[1] - s[0]
LTRK = s[-1] + ds
d = np.genfromtxt(RES, delimiter=',', names=True)

# uniform playback timeline
t = np.arange(0, d['time'][-1], 1.0/FPS)
def at(name): return np.interp(t, d['time'], d[name])

# map (s, n) -> global pose; psi tables are unwrapped so interpolation is safe
sq = np.concatenate([s, [LTRK]])
per = lambda f, closer=None: np.concatenate([f, [f[0] if closer is None else closer]])
sMod = np.mod(at('s'), LTRK)
psiTrk = np.interp(sMod, sq, per(psic, psic[0] + (psic[-1] - psic[0])*len(s)/(len(s)-1)))
# (use the actual wrap: last->first differs by the lap's -2pi; append first+total turn)
psiTrk = np.interp(sMod, sq, np.concatenate([psic, [psic[0] - 2*np.pi]]))
n = at('n')
xCar = np.interp(sMod, sq, per(xc)) - n*np.sin(psiTrk)
yCar = np.interp(sMod, sq, per(yc)) + n*np.cos(psiTrk)
psiCar = psiTrk + np.radians(at('dpsiDeg'))

r2 = lambda a: [round(float(v), 2) for v in a]
r4 = lambda a: [round(float(v), 4) for v in a]
data = {
    'fps': FPS,
    'L': round(float(LTRK), 1),
    'lapTime': round(float(d['time'][-1]), 2),
    'track': {'x': r2(xc[::1]), 'y': r2(yc[::1]), 'k': r4(kap[::1]), 'ds': ds},
    'car': {
        'x': r2(xCar), 'y': r2(yCar), 'psi': r4(psiCar), 's': r2(at('s')),
        'v': r2(at('vKmh')), 'vRef': r2(at('vRefKmh')), 'yaw': r2(at('yawRateDegS')),
        'ax': [round(float(v), 3) for v in at('axG')],
        'ay': [round(float(v), 3) for v in at('ayG')],
        'st': r2(at('deltaDeg')),
        'thr': [round(float(v), 3) for v in np.clip(at('Pdrive')/PMAX, 0, 1)],
        'brk': [round(float(v), 3) for v in np.clip(-at('axG')/0.9, 0, 1)],
    },
}

HTML = r"""<title>Nordschleife lap — chase cam</title>
<style>
:root{
  --ground:#151a16; --asphalt:#3b3e44; --edge:#b9b4a6; --kerb:#a8473c;
  --panel:rgba(12,14,16,.78); --panel-line:rgba(255,255,255,.08);
  --ink:#e8e6df; --ink-dim:#9a9891; --accent:#e8c15a;
  --good:#7dbf7f; --bad:#d06355;
  --mono:ui-monospace,'SF Mono',Menlo,Consolas,monospace;
  --sans:system-ui,-apple-system,'Segoe UI',sans-serif;
}
html,body{margin:0;height:100%;overflow:hidden;background:#0c0e10;color:var(--ink);font-family:var(--sans)}
canvas#world{position:fixed;inset:0;width:100vw;height:100vh;display:block}
.hud{position:fixed;padding:10px 14px;background:var(--panel);border:1px solid var(--panel-line);
     border-radius:10px;backdrop-filter:blur(6px);user-select:none}
#timer{top:14px;left:14px;font-family:var(--mono);font-variant-numeric:tabular-nums}
#timer .big{font-size:22px;letter-spacing:.5px}
#timer .sub{font-size:11px;color:var(--ink-dim);margin-top:2px}
#mini{top:14px;right:14px;padding:8px}
#mini canvas{display:block}
#tele{bottom:64px;left:14px;display:flex;gap:18px;align-items:flex-end}
#speed{font-family:var(--mono);font-variant-numeric:tabular-nums;line-height:1}
#speed .v{font-size:52px;font-weight:600}
#speed .u{font-size:11px;color:var(--ink-dim);letter-spacing:1.5px;text-transform:uppercase}
.rows{display:grid;grid-template-columns:auto auto;gap:2px 10px;
      font-family:var(--mono);font-size:13px;font-variant-numeric:tabular-nums}
.rows .l{color:var(--ink-dim);font-size:10px;letter-spacing:1px;text-transform:uppercase;align-self:center}
#gg{position:relative}
#gg canvas{display:block}
#gg .cap{font-size:9px;color:var(--ink-dim);letter-spacing:1px;text-transform:uppercase;text-align:center;margin-top:3px}
#bars{display:flex;gap:8px;align-items:flex-end;height:86px}
.bar{width:14px;height:70px;background:rgba(255,255,255,.07);border-radius:4px;position:relative;overflow:hidden}
.bar i{position:absolute;bottom:0;left:0;right:0;border-radius:4px}
.bar.t i{background:var(--good)} .bar.b i{background:var(--bad)}
.bar .bl{position:absolute;top:100%;left:50%;transform:translateX(-50%);margin-top:4px;
         font-size:9px;color:var(--ink-dim);letter-spacing:1px}
#steerwrap{display:flex;flex-direction:column;gap:4px;align-items:center}
#steer{width:120px;height:10px;background:rgba(255,255,255,.07);border-radius:5px;position:relative}
#steer i{position:absolute;top:1px;bottom:1px;width:8px;border-radius:4px;background:var(--accent);left:calc(50% - 4px)}
#steer:after{content:'';position:absolute;left:50%;top:-2px;bottom:-2px;width:1px;background:rgba(255,255,255,.25)}
#steerwrap .bl{font-size:9px;color:var(--ink-dim);letter-spacing:1px;text-transform:uppercase}
#ctrl{bottom:14px;left:14px;right:14px;display:flex;gap:12px;align-items:center;padding:8px 14px}
#ctrl button{background:none;border:1px solid var(--panel-line);color:var(--ink);border-radius:6px;
  font-family:var(--mono);font-size:13px;padding:4px 12px;cursor:pointer}
#ctrl button:hover{border-color:var(--accent);color:var(--accent)}
#ctrl button:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
#scrub{flex:1;accent-color:var(--accent)}
#ctrl .spd{color:var(--ink-dim);font-size:11px}
#ctrl select{background:#1a1d21;color:var(--ink);border:1px solid var(--panel-line);border-radius:6px;
  font-family:var(--mono);padding:3px 6px}
#zooml{font-size:10px;color:var(--ink-dim);letter-spacing:1px}
</style>
<canvas id="world"></canvas>
<div class="hud" id="timer"><div class="big" id="tval">0:00.0</div>
  <div class="sub"><span id="sval">0.00</span> km of <span id="ltot"></span> km — Nürburgring Nordschleife (planar sim)</div></div>
<div class="hud" id="mini"><canvas id="minimap" width="230" height="230"></canvas></div>
<div class="hud" id="tele">
  <div id="speed"><div class="v" id="vval">0</div><div class="u">km/h</div></div>
  <div class="rows">
    <span class="l">yaw rate</span><span id="yawval">+0.0 °/s</span>
    <span class="l">lat acc</span><span id="ayval">+0.00 g</span>
    <span class="l">long acc</span><span id="axval">+0.00 g</span>
    <span class="l">v ref</span><span id="vrefval">0 km/h</span>
  </div>
  <div id="gg"><canvas id="ggc" width="92" height="92"></canvas><div class="cap">g-g</div></div>
  <div id="bars">
    <div class="bar t"><i id="thr"></i><span class="bl">THR</span></div>
    <div class="bar b"><i id="brk"></i><span class="bl">BRK</span></div>
  </div>
  <div id="steerwrap"><div id="steer"><i id="sti"></i></div><span class="bl">steer</span></div>
</div>
<div class="hud" id="ctrl">
  <button id="play" aria-label="play/pause">▶</button>
  <input id="scrub" type="range" min="0" max="1000" value="0">
  <span class="spd">speed <select id="rate"><option>0.5</option><option selected>1</option><option>2</option><option>4</option><option>8</option></select>×</span>
  <span class="spd" id="zooml">zoom <input id="zoom" type="range" min="50" max="260" value="120" style="width:90px;vertical-align:middle;accent-color:var(--accent)"></span>
</div>
<script>
const D = __DATA__;
const TX = D.track.x, TY = D.track.y, TK = D.track.k, N = TX.length;
const C = D.car, FPS = D.fps, T = D.lapTime;
const cv = document.getElementById('world'), ctx = cv.getContext('2d');
let W, H, dpr = Math.min(devicePixelRatio || 1, 2);
function resize(){ W = innerWidth; H = innerHeight; cv.width = W*dpr; cv.height = H*dpr; }
addEventListener('resize', resize); resize();

// world geometry as Path2D (world meters); kerb segments where |kappa| is high
const road = new Path2D();
road.moveTo(TX[0], TY[0]);
for (let i = 1; i < N; i++) road.lineTo(TX[i], TY[i]);
road.closePath();
const kerbs = new Path2D();
for (let i = 0; i < N; i++){
  if (Math.abs(TK[i]) > 0.008){
    const j = (i+1) % N;
    kerbs.moveTo(TX[i], TY[i]); kerbs.lineTo(TX[j], TY[j]);
  }
}
const ROADW = 7.0;

// playback state
let playing = !matchMedia('(prefers-reduced-motion: reduce)').matches;
let tSim = 0, rate = 1, viewH = 120, last = performance.now();
const lerp = (a, f, i) => a[i] + (a[Math.min(i+1, a.length-1)] - a[i])*f;
function sample(tq){
  const fi = Math.min(Math.max(tq*FPS, 0), C.x.length - 1.001);
  const i = Math.floor(fi), f = fi - i;
  return { x: lerp(C.x, f, i), y: lerp(C.y, f, i), psi: lerp(C.psi, f, i),
           s: lerp(C.s, f, i), v: lerp(C.v, f, i), vRef: lerp(C.vRef, f, i),
           yaw: lerp(C.yaw, f, i), ax: lerp(C.ax, f, i), ay: lerp(C.ay, f, i),
           st: lerp(C.st, f, i), thr: lerp(C.thr, f, i), brk: lerp(C.brk, f, i) };
}

function drawWorld(p){
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.fillStyle = getComputedStyle(document.documentElement).getPropertyValue('--ground');
  ctx.fillRect(0, 0, W, H);
  const k = H/viewH;                       // px per meter
  const cx = W/2, cy = H*0.62;
  // world -> screen: rotate so car heading points up (canvas y down)
  const cs = Math.cos(p.psi), sn = Math.sin(p.psi);
  // screen = k * [ right; -fwd ] * (world - car) + (cx, cy)
  ctx.setTransform(k*sn*dpr, -k*cs*dpr, -k*cs*dpr, -k*sn*dpr,
                   dpr*(cx - k*(sn*p.x - cs*p.y)), dpr*(cy + k*(cs*p.x + sn*p.y)));
  ctx.lineJoin = 'round'; ctx.lineCap = 'butt';
  const css = getComputedStyle(document.documentElement);
  ctx.strokeStyle = css.getPropertyValue('--kerb'); ctx.lineWidth = ROADW + 1.9;
  ctx.stroke(kerbs);
  ctx.strokeStyle = css.getPropertyValue('--edge'); ctx.lineWidth = ROADW + 0.7;
  ctx.stroke(road);
  ctx.strokeStyle = css.getPropertyValue('--asphalt'); ctx.lineWidth = ROADW;
  ctx.stroke(road);
  ctx.strokeStyle = 'rgba(255,255,255,.28)'; ctx.lineWidth = 0.14;
  ctx.setLineDash([3, 6]); ctx.stroke(road); ctx.setLineDash([]);
  // start/finish
  ctx.save();
  ctx.translate(TX[0], TY[0]);
  const p0 = Math.atan2(TY[1]-TY[0], TX[1]-TX[0]);
  ctx.rotate(p0);
  ctx.fillStyle = 'rgba(255,255,255,.75)';
  for (let i = -3; i < 3; i++)
    for (let j = 0; j < 2; j++)
      if ((i + j) % 2 === 0) ctx.fillRect(j*0.7, i*1.16, 0.7, 1.16);
  ctx.restore();
  // car (world frame)
  ctx.save();
  ctx.translate(p.x, p.y); ctx.rotate(p.psi);
  ctx.fillStyle = 'rgba(0,0,0,.35)';
  ctx.beginPath(); ctx.roundRect(-2.3, -1.02, 4.6, 2.04, 0.5); ctx.fill();
  ctx.fillStyle = '#d8dade';
  ctx.beginPath(); ctx.roundRect(-2.2, -0.9, 4.4, 1.8, 0.5); ctx.fill();
  ctx.fillStyle = '#23262b';
  ctx.beginPath(); ctx.roundRect(0.15, -0.72, 1.15, 1.44, 0.25); ctx.fill(); // glasshouse
  ctx.fillStyle = css.getPropertyValue('--accent');
  ctx.fillRect(1.9, -0.62, 0.35, 1.24);                                     // nose band
  ctx.restore();
}

// minimap (static path prerendered)
const mm = document.getElementById('minimap'), mctx = mm.getContext('2d');
let mmBase, mmap;
{
  const xmin = Math.min(...TX), xmax = Math.max(...TX);
  const ymin = Math.min(...TY), ymax = Math.max(...TY);
  const pad = 10, sz = 230;
  const kk = (sz - 2*pad)/Math.max(xmax - xmin, ymax - ymin);
  mmap = (x, y) => [pad + (x - xmin)*kk + (sz - 2*pad - (xmax - xmin)*kk)/2,
                    sz - pad - (y - ymin)*kk - (sz - 2*pad - (ymax - ymin)*kk)/2];
  mmBase = document.createElement('canvas');
  mmBase.width = mmBase.height = sz;
  const b = mmBase.getContext('2d');
  b.strokeStyle = 'rgba(232,230,223,.55)'; b.lineWidth = 1.6; b.beginPath();
  const [x0, y0] = mmap(TX[0], TY[0]); b.moveTo(x0, y0);
  for (let i = 1; i < N; i++){ const [xx, yy] = mmap(TX[i], TY[i]); b.lineTo(xx, yy); }
  b.closePath(); b.stroke();
}
function drawMini(p){
  mctx.clearRect(0, 0, 230, 230);
  mctx.drawImage(mmBase, 0, 0);
  const [xx, yy] = mmap(p.x, p.y);
  mctx.fillStyle = getComputedStyle(document.documentElement).getPropertyValue('--accent');
  mctx.beginPath(); mctx.arc(xx, yy, 4, 0, 7); mctx.fill();
  mctx.strokeStyle = 'rgba(0,0,0,.5)'; mctx.lineWidth = 1;
  mctx.beginPath(); mctx.arc(xx, yy, 4, 0, 7); mctx.stroke();
}

// g-g dot with short trail
const gg = document.getElementById('ggc').getContext('2d');
const trail = [];
function drawGG(p){
  gg.clearRect(0, 0, 92, 92);
  gg.strokeStyle = 'rgba(255,255,255,.18)';
  gg.beginPath(); gg.arc(46, 46, 40, 0, 7); gg.stroke();
  gg.beginPath(); gg.arc(46, 46, 20, 0, 7); gg.stroke();
  gg.beginPath(); gg.moveTo(6, 46); gg.lineTo(86, 46); gg.moveTo(46, 6); gg.lineTo(46, 86); gg.stroke();
  trail.push([p.ay, p.ax]); if (trail.length > 25) trail.shift();
  for (let i = 0; i < trail.length; i++){
    const a = i/trail.length;
    gg.fillStyle = `rgba(232,193,90,${0.12 + 0.5*a})`;
    gg.beginPath();
    gg.arc(46 + trail[i][0]*40, 46 - trail[i][1]*40, i === trail.length-1 ? 3.5 : 2, 0, 7);
    gg.fill();
  }
}

const $ = id => document.getElementById(id);
const fmt = (v, dg, u) => (v >= 0 ? '+' : '−') + Math.abs(v).toFixed(dg) + ' ' + u;
function drawHUD(p, tq){
  $('vval').textContent = Math.round(p.v);
  $('vrefval').textContent = Math.round(p.vRef) + ' km/h';
  $('yawval').textContent = fmt(p.yaw, 1, '°/s');
  $('ayval').textContent = fmt(p.ay, 2, 'g');
  $('axval').textContent = fmt(p.ax, 2, 'g');
  $('thr').style.height = (p.thr*100).toFixed(0) + '%';
  $('brk').style.height = (Math.min(p.brk, 1)*100).toFixed(0) + '%';
  $('sti').style.left = `calc(50% - 4px + ${(-p.st/20*56).toFixed(1)}px)`;
  const mn = Math.floor(tq/60), sc = tq - 60*mn;
  $('tval').textContent = `${mn}:${sc < 10 ? '0' : ''}${sc.toFixed(1)}`;
  $('sval').textContent = (p.s/1000).toFixed(2);
  if (!scrubbing) $('scrub').value = Math.round(1000*tq/T);
  drawGG(p); drawMini(p);
}

function frame(now){
  const dt = Math.min((now - last)/1000, 0.1); last = now;
  if (playing){ tSim += dt*rate; if (tSim >= T){ tSim = 0; } }
  const p = sample(tSim);
  drawWorld(p); drawHUD(p, tSim);
  requestAnimationFrame(frame);
}
$('ltot').textContent = (D.L/1000).toFixed(2);
let scrubbing = false;
$('play').onclick = () => { playing = !playing; $('play').textContent = playing ? '❚❚' : '▶'; };
$('rate').onchange = e => rate = parseFloat(e.target.value);
$('zoom').oninput = e => viewH = 310 - parseFloat(e.target.value);
$('scrub').addEventListener('pointerdown', () => scrubbing = true);
$('scrub').addEventListener('pointerup', () => scrubbing = false);
$('scrub').addEventListener('input', e => { tSim = e.target.value/1000*T; });
addEventListener('keydown', e => {
  if (e.code === 'Space'){ e.preventDefault(); $('play').click(); }
  if (e.code === 'ArrowRight') tSim = Math.min(tSim + 5, T - 0.01);
  if (e.code === 'ArrowLeft') tSim = Math.max(tSim - 5, 0);
});
$('play').textContent = playing ? '❚❚' : '▶';
requestAnimationFrame(t0 => { last = t0; requestAnimationFrame(frame); });
</script>
"""

html = HTML.replace('__DATA__', json.dumps(data, separators=(',', ':')))
os.makedirs(os.path.dirname(OUT), exist_ok=True)
open(OUT, 'w').write(html)
kb = os.path.getsize(OUT)/1024
print(f'wrote {OUT} ({kb:.0f} KB, {len(t)} telemetry frames at {FPS} Hz)')
