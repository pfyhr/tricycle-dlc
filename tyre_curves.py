"""Elise brush-tyre characteristics in Pacejka's Chapter-3 conventions.

Evaluates the exact force/moment forms the plants use - lambda = 1 - theta*|sigma|
(3.8), Fy = mu*Fz*(1 - lambda^3)*sgn(alpha) (3.11), Mz = -mu*Fz*a*lambda^3*(1-lambda)
*sgn(alpha) (3.12), pneumatic trail t0 = a/3 at vanishing slip (3.14) - with the
Magic-Formula load function C_alpha = c1*sin(2*atan(Fz/c2)) for the cornering
stiffness, at the calibrated Elise's front/rear parameters and loads.

Writes outputs/svg/elise_tyres.svg for the README's tyre section.
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'tracks'))
from cars import CARS

t = CARS['elise']['trike']
MU = t['muF']

def brush(al, Fz, c1, c2, FzNom, ap0, mu):
    """Exact port of the plants' brushForces (Pacejka ch. 3 forms)."""
    Fz = np.maximum(Fz, 50.0)
    Ca = c1*np.sin(2*np.arctan(Fz/c2))
    ap = ap0*np.sqrt(Fz/FzNom)
    sy = np.tan(al)
    sA = np.sqrt(sy**2 + 1e-12)
    sg = sy/sA
    th = Ca/(3*mu*Fz)
    u = 0.5*(th*sA + 1 - np.sqrt((th*sA - 1)**2 + 1e-4))   # smooth clip of theta*sigma
    lam = 1 - u
    Fy = mu*Fz*(1 - lam**3)*sg
    Mz = -mu*Fz*ap*u*lam**3*sg
    return Fy, Mz, ap

alp = np.radians(np.linspace(0, 14, 400))
FRONT = dict(c1=t['c1F'], c2=t['c2F'], FzNom=t['FzNomF'], ap0=t['ap0F'], mu=MU)
REAR = dict(c1=t['c1R'], c2=t['c2R'], FzNom=t['FzNomR'], ap0=t['ap0R'], mu=MU)

fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.2))
fig.patch.set_facecolor('white')
for a_ in ax:
    a_.set_facecolor('white'); a_.grid(alpha=0.25)

# --- Fy(alpha): front tyre at loads spanning the lateral load transfer -------------
loads = [0.5, 1.0, 1.6]
blues = ['#9dc3e6', '#2b6cb0', '#173f66']
for f, c in zip(loads, blues):
    Fz = f*t['FzNomF']
    Fy, _, _ = brush(alp, Fz, **FRONT)
    ax[0].plot(np.degrees(alp), Fy/1e3, color=c, lw=1.6,
               label=f'front, $F_z$ = {f:.1f}·$F_{{z,nom}}$')
    ax[0].axhline(MU*Fz/1e3, color=c, lw=0.6, ls=':')
FyR, _, _ = brush(alp, t['FzNomR'], **REAR)
ax[0].plot(np.degrees(alp), FyR/1e3, color='#dd6b20', lw=1.6, ls='--',
           label='rear (lumped), $F_{z,nom}$')
ax[0].axhline(MU*t['FzNomR']/1e3, color='#dd6b20', lw=0.6, ls=':')
ax[0].set_xlabel('slip angle α  [°]'); ax[0].set_ylabel('$F_y$  [kN]')
ax[0].set_title('side force (dotted: µ$F_z$)', fontsize=9)
ax[0].legend(fontsize=7, loc='center right')

# --- Mz(alpha): the aligning-moment dip that "lightens" the steering ----------------
for f, c in zip(loads, blues):
    Fz = f*t['FzNomF']
    _, Mz, _ = brush(alp, Fz, **FRONT)
    ax[1].plot(np.degrees(alp), Mz, color=c, lw=1.6)
_, MzR, _ = brush(alp, t['FzNomR'], **REAR)
ax[1].plot(np.degrees(alp), MzR, color='#dd6b20', lw=1.6, ls='--')
ax[1].set_xlabel('slip angle α  [°]'); ax[1].set_ylabel('$M_z$  [N·m]')
ax[1].set_title('aligning moment ($M_z$ = −t·$F_y$)', fontsize=9)

# --- pneumatic trail t(alpha)/a: a/3 at zero slip, collapses at the limit -----------
for f, c in zip(loads, blues):
    Fz = f*t['FzNomF']
    Fy, Mz, ap = brush(alp, Fz, **FRONT)
    tr = np.where(np.abs(Fy) > 1, -Mz/np.maximum(np.abs(Fy), 1), ap/3)/ap
    ax[2].plot(np.degrees(alp), tr, color=c, lw=1.6)
ax[2].axhline(1/3, color='#999999', lw=0.7, ls=':')
ax[2].text(9.5, 0.345, 't$_0$ = a/3', fontsize=8, color='#666666')
ax[2].set_xlabel('slip angle α  [°]'); ax[2].set_ylabel('pneumatic trail  t/a')
ax[2].set_ylim(-0.02, 0.4)
ax[2].set_title('trail collapse → steering goes light', fontsize=9)

fig.suptitle('Calibrated Elise brush tyres — exact Pacejka ch. 3 forms '
             '(λ = 1−θ|σ|, F$_y$ = µF$_z$(1−λ³), M$_z$ = −µF$_z$·a·λ³(1−λ)), '
             'MF load function for C$_α$', fontsize=9)
fig.tight_layout()
out = os.path.join(HERE, 'outputs', 'svg', 'elise_tyres.svg')
fig.savefig(out, bbox_inches='tight')
print('wrote', out)
