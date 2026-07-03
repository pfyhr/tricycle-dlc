"""Render a clean diagram of Tricycle.Examples.DoubleLaneChange using the *real*
OpenModelica component icons (from the model-instance JSON) but placed on a tidy grid
with orthogonal wiring — authentic icons, controlled layout.
Produces outputs/svg/dlc_diagram.svg."""
import json, subprocess, shutil, re, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Ellipse, FancyBboxPatch, Arc

MODEL = 'Tricycle.Examples.DoubleLaneChange'
OMC = shutil.which('omc') or '/Users/pontus/opt/openmodelica/bin/omc'
# pin OPENMODELICAHOME to the omc install so generated code always finds its headers
# (openmodelica.h / omc_simulation_settings.h), even via symlinks or relocated installs
os.environ.setdefault('OPENMODELICAHOME',
                      os.path.dirname(os.path.dirname(os.path.realpath(OMC))))
HERE = os.path.dirname(os.path.abspath(__file__)); MODELICA = os.path.join(HERE, 'modelica')

mos = ('loadModel(Modelica); loadFile("Tricycle.mo");\n'
       f'writeFile("/tmp/mi.json", getModelInstance({MODEL}, prettyPrint=true));\n')
open('/tmp/_render.mos', 'w').write(mos)
subprocess.run([OMC, '/tmp/_render.mos'], cwd=MODELICA, check=True, capture_output=True)
d = json.loads(open('/tmp/mi.json').read())

# ---- tidy manual layout: (col, row); forward chain on row 0, feedback on row -1 -----
SP, RSP, W = 2.25, 3.0, 2.0             # column spacing, row spacing, box size
LABEL = {'driver': 'ISO 3888-1\npreview driver', 'refGain': 'ratio $i_S$',
         'arm': 'driver arm\n(2 Hz)', 'hwTorqueSensor': 'torque\nsensor',
         'steering': 'manual rack\n($i_S$=20, unassisted)', 'trike': 'tricycle\nvehicle'}
GRID = {'driver': (0, 0), 'refGain': (1, 0), 'arm': (2, 0), 'hwTorqueSensor': (3, 0),
        'steering': (4, 0), 'trike': (5, 0)}
POS = {k: (c*SP, r*RSP) for k, (c, r) in GRID.items()}


def col(c):
    if not isinstance(c, list) or len(c) != 3 or c[0] < 0:
        return None
    return tuple(max(0, min(255, v))/255 for v in c)


def name_of(e):
    n = e.get('name'); return n.get('name') if isinstance(n, dict) else n


def collect_graphics(t):
    """Compose a class icon: inherited (extends) graphics first, then the class's own
    on top — OM/OMEdit don't merge inheritance into Icon.graphics, so we do it here."""
    if not isinstance(t, dict):
        return []
    out = []
    for e in t.get('elements', []):
        if isinstance(e, dict) and e.get('$kind') == 'extends':
            out += collect_graphics(e.get('baseClass', {}))
    icon = t.get('annotation', {}).get('Icon')
    if isinstance(icon, dict):
        out += [g for g in icon.get('graphics', []) if isinstance(g, dict)]
    return out


fig, ax = plt.subplots(figsize=(10.5, 4.2))


def mapper(cx, cy, icon_ext):
    (ix1, iy1), (ix2, iy2) = icon_ext
    def f(p, o=(0, 0)):
        px, py = p[0]+o[0], p[1]+o[1]
        return (cx + (px-(ix1+ix2)/2)/(ix2-ix1)*W,
                cy + (py-(iy1+iy2)/2)/(iy2-iy1)*W)
    return f


def draw_prim(g, f):
    k = g.get('name'); el = g.get('elements', [])
    if not el or el[0] is False:
        return
    o = el[1] if len(el) > 1 and isinstance(el[1], list) else (0, 0)
    if k in ('Rectangle', 'Ellipse'):
        line, fill = col(el[3]), col(el[4])
        fp = el[6]['name'] if isinstance(el[6], dict) else 'None'
        ext = el[9] if k == 'Rectangle' else el[8]   # Rectangle has borderPattern at 8, Ellipse's extent is at 8
        (x1, y1), (x2, y2) = f(ext[0], o), f(ext[1], o)
        fc = fill if (fill and 'None' not in fp) else 'none'
        ec = line if line else 'none'
        if k == 'Rectangle':
            ax.add_patch(Polygon([(x1, y1), (x2, y1), (x2, y2), (x1, y2)], closed=True,
                                 fc=fc, ec=ec, lw=0.7, zorder=3))
        else:
            cx0, cy0 = (x1+x2)/2, (y1+y2)/2
            sa = el[9] if len(el) > 9 and not isinstance(el[9], (list, dict)) else 0
            ea = el[10] if len(el) > 10 and not isinstance(el[10], (list, dict)) else 360
            rot = el[2] if isinstance(el[2], (int, float)) else 0
            if abs(ea - sa) < 359.5 and not (sa == 0 and ea == 360):   # partial -> open arc (e.g. torque arrows)
                ax.add_patch(Arc((cx0, cy0), abs(x2-x1), abs(y2-y1), angle=rot,
                                 theta1=sa, theta2=ea, edgecolor=ec if ec != 'none' else 'k',
                                 lw=0.9, zorder=3))
            else:
                ax.add_patch(Ellipse((cx0, cy0), abs(x2-x1), abs(y2-y1), angle=rot,
                                     fc=fc, ec=ec, lw=0.7, zorder=3))
    elif k in ('Line', 'Polygon'):
        pts = el[3] if k == 'Line' else el[8]
        P = [f(p, o) for p in pts]
        if k == 'Polygon':
            fill = col(el[4]); fp = el[6]['name'] if isinstance(el[6], dict) else 'None'
            ax.add_patch(Polygon(P, closed=True, fc=fill if (fill and 'None' not in fp) else 'none',
                                 ec=col(el[3]) or 'k', lw=0.7, zorder=3))
        else:
            ax.plot([p[0] for p in P], [p[1] for p in P], color=col(el[4]) or 'k',
                    lw=0.9, solid_capstyle='round', zorder=3)
    elif k == 'Text':
        if '%name' in (el[9] or ''):
            return                                 # skip OM name labels; we draw our own
        s = re.sub(r'%\w+', '', el[9] or '').strip()
        if not s or s.endswith('='):
            return
        (x1, y1), (x2, y2) = f(el[8][0], o), f(el[8][1], o)
        ax.text((x1+x2)/2, (y1+y2)/2, s, ha='center', va='center', fontsize=6,
                color=col(el[11]) or (0, 0, 0), zorder=4)


types = {name_of(e): e.get('type') for e in d['elements'] if e.get('$kind') == 'component'}
for name, (cx, cy) in POS.items():
    t = types.get(name)
    if isinstance(t, dict):
        icon = t.get('annotation', {}).get('Icon', {}) or {}
        ext = (icon.get('coordinateSystem', {}) or {}).get('extent', [[-100, -100], [100, 100]])
        f = mapper(cx, cy, ext)
        for g in collect_graphics(t):
            try:
                draw_prim(g, f)
            except Exception:
                pass
    ax.text(cx, cy+W/2+0.32, LABEL.get(name, name), ha='center', va='bottom',
            fontsize=7.5, color='tab:blue', zorder=5)

# ---- orthogonal wiring from the real connection list --------------------------------
def cname(ref):
    return ref['parts'][0]['name']

def edge(name, side):
    cx, cy = POS[name]
    return {'r': (cx+W/2, cy), 'l': (cx-W/2, cy), 't': (cx, cy+W/2), 'b': (cx, cy-W/2)}[side]

import math

def kind_of(color):
    r, g, bl = color
    if bl > 0.35 and bl >= r and bl >= g:
        return 'sig'                 # signal (Real) connector
    if g > 0.3 and g >= r and g >= bl:
        return 'trans'               # translational flange
    return 'rot'                     # rotational flange

def connector(p, color, ang=0.0):
    """Draw the proper Modelica connector glyph: signal=triangle, rotational=circle,
    translational=square (the shape carries the domain meaning)."""
    k = kind_of(color)
    if k == 'sig':                                  # causal signal -> filled triangle
        base = [(-0.07, -0.10), (0.11, 0.0), (-0.07, 0.10)]
        ca, sa = math.cos(ang), math.sin(ang)
        pts = [(p[0]+x*ca-y*sa, p[1]+x*sa+y*ca) for x, y in base]
        ax.add_patch(Polygon(pts, closed=True, fc=(0.0, 0.0, 0.55), ec=(0.0, 0.0, 0.55), zorder=7))
    elif k == 'rot':                                # rotational flange -> filled grey disc
        ax.add_patch(Ellipse(p, 0.16, 0.16, fc=(0.5, 0.5, 0.5), ec='k', lw=0.5, zorder=7))
    else:                                           # translational flange -> filled square
        ax.add_patch(plt.Rectangle((p[0]-0.08, p[1]-0.08), 0.16, 0.16,
                                   fc=(0.45, 0.6, 0.4), ec='k', lw=0.5, zorder=7))

for c in d.get('connections', []):
    a, b = cname(c['lhs']), cname(c['rhs'])
    if a not in POS or b not in POS:
        continue
    color = col(c.get('annotation', {}).get('Line', {}).get('color')) or (0.1, 0.2, 0.5)
    (ax_, ay_), (bx_, by_) = POS[a], POS[b]
    if ay_ == by_:                                  # same row: straight horizontal (flow ->)
        lo, hi = (a, b) if ax_ < bx_ else (b, a)
        p1, p2 = edge(lo, 'r'), edge(hi, 'l')
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color, lw=1.6, zorder=2)
        connector(p1, color, 0.0); connector(p2, color, 0.0)
    elif 'angleSensor' in (a, b):                   # feedback loop
        other = b if a == 'angleSensor' else a
        if other == 'rotor':                        # rotor.flange -> sensor (vertical, rotational)
            p1, p2 = edge('rotor', 'b'), edge('angleSensor', 't')
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color, lw=1.6, zorder=2)
            connector(p1, color); connector(p2, color)
        else:                                       # sensor.phi -> controller.u_m (down, left, up)
            s = edge('angleSensor', 'b'); cpt = edge('controller', 'b')
            ylow = s[1] - 0.9
            ax.plot([s[0], s[0], cpt[0], cpt[0]], [s[1], ylow, ylow, cpt[1]],
                    color=color, lw=1.6, zorder=2)
            connector(s, color, -math.pi/2)         # signal leaves sensor downward
            connector(cpt, color, math.pi/2)        # signal enters controller upward


# ---- hand-drawn extra: vehicle-state feedback to the driver (outermost loop) --------
SIG = (0.0, 0.0, 0.55)
tb = edge('trike', 'b'); db = edge('driver', 'b')
ylow2 = edge('trike', 'b')[1] - 1.3
ax.plot([tb[0], tb[0], db[0], db[0]], [tb[1], ylow2, ylow2, db[1]], color=SIG, lw=1.6, zorder=2)
connector(tb, SIG, -math.pi/2); connector(db, SIG, math.pi/2)
ax.text((tb[0]+db[0])/2, ylow2-0.15, 'vehicle states  $X,\\ Y,\\ \\psi,\\ v_y,\\ r$',
        ha='center', va='top', fontsize=7.5, color=SIG)
ax.set_aspect('equal'); ax.axis('off')
xs = [p[0] for p in POS.values()]; ys = [p[1] for p in POS.values()]
ax.set_xlim(min(xs)-W, max(xs)+W); ax.set_ylim(min(ys)-W-1.6, max(ys)+W+1.2)
ax.set_title(f'{MODEL}  —  OpenModelica components', fontweight='bold', fontsize=11)
plt.tight_layout()
os.makedirs(os.path.join(HERE, 'outputs/svg'), exist_ok=True)
plt.savefig(os.path.join(HERE, 'outputs/svg/dlc_diagram.svg'), facecolor="white")
plt.savefig('/tmp/om_diagram.png', dpi=140)
print('rendered outputs/svg/dlc_diagram.svg (clean layout, real OM icons)')
