"""Print the README's car-parameter appendix as markdown, straight from cars.py.

Regenerate after any preset change and paste over the table in the README appendix:
    python3 car_table.py
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tracks'))
from cars import CARS

ORDER = ['elise', 'miata', 'm140', 'clubman']
ROWS = [
    ('m',      'kg',      'mass (incl. driver)',                 lambda c: f"{c['trike']['m']:.0f}"),
    ('I_zz',   'kg·m²',   'yaw inertia',                         lambda c: f"{c['trike']['Izz']:.0f}"),
    ('a / b',  'm',       'CG to front / rear axle',             lambda c: f"{c['trike']['a']:.2f} / {c['trike']['b']:.2f}"),
    ('L',      'm',       'wheelbase',                           lambda c: f"{c['trike']['a']+c['trike']['b']:.2f}"),
    ('front weight', '%', 'static front axle share = b/L',       lambda c: f"{100*c['trike']['b']/(c['trike']['a']+c['trike']['b']):.0f}"),
    ('t_f',    'm',       'front track width',                   lambda c: f"{c['trike']['tf']:.2f}"),
    ('h_cg',   'm',       'CG height',                           lambda c: f"{c['trike']['hcg']:.2f}"),
    ('ξ_F',    '–',       'front share of the roll couple',      lambda c: f"{c['trike']['xiF']:.2f}"),
    ('k_Bf',   '–',       'front share of brake force',          lambda c: f"{c['trike']['kBf']:.2f}"),
    ('P_max',  'kW',      'peak drive power (rear wheel)',       lambda c: f"{c['Pmax']/1e3:.0f}"),
    ('C_dA',   'm²',      'drag area',                           lambda c: f"{c['trike']['CdA']:.2f}"),
    ('C_lA',   'm²',      'lift (downforce) area',               lambda c: f"{c['trike'].get('ClA',0.0):.1f}"),
    ('aeroBal','–',       'front share of downforce',            lambda c: f"{c['trike'].get('aeroBal',0.42):.2f}" if c['trike'].get('ClA',0) else '–'),
    ('C_rr',   '–',       'rolling resistance',                  lambda c: f"{c['trike']['Crr']:.3f}"),
    ('µ',      '–',       'tyre friction coefficient',           lambda c: f"{c['trike']['muF']:.2f}"),
    ('ayFrac', '–',       'driver share of the lateral budget',  lambda c: f"{c['profile']['ayFrac']:.2f}"),
    ('c1F / c1R', 'N/rad','tyre stiffness scale, front wheel / lumped rear', lambda c: f"{c['trike']['c1F']:.0f} / {c['trike']['c1R']:.0f}"),
    ('c2F / c2R', 'N',    'stiffness load scale (C_α = c1·sin(2·atan(F_z/c2)))', lambda c: f"{c['trike']['c2F']:.0f} / {c['trike']['c2R']:.0f}"),
    ('F_z,nom F / R', 'N','nominal wheel loads (trail + contact length scale)', lambda c: f"{c['trike']['FzNomF']:.0f} / {c['trike']['FzNomR']:.0f}"),
    ('a_p0 F / R', 'm',   'contact half-length at F_z,nom',      lambda c: f"{c['trike']['ap0F']:.3f} / {c['trike']['ap0R']:.3f}"),
    ('K_us',   'rad·s²/m','understeer gradient (steer FF)',      lambda c: f"{c['driver']['Kus']*1e3:.2f}e-3"),
    ('K_LA',   'rad/m',   'lookahead feedback gain (× KMUL 1.7)',lambda c: f"{c['driver']['KLA']:.2f}"),
    ('K_r',    'rad·s/rad','yaw damping gain',                   lambda c: f"{c['driver']['Kr']:.2f}"),
]

names = {k: CARS[k]['display'].split(' (')[0] for k in ORDER}
hdr = '| parameter | unit | meaning | ' + ' | '.join(names[k] for k in ORDER) + ' |'
sep = '|---|---|---|' + '--:|'*len(ORDER)
print(hdr)
print(sep)
for sym, unit, meaning, fn in ROWS:
    cells = ' | '.join(fn(CARS[k]) for k in ORDER)
    print(f'| {sym} | {unit} | {meaning} | {cells} |')
