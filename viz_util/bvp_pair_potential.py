"""
Plot: compute the pair potential using the linear superposition approx.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import matplotlib
import numdifftools as nd
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from viz_util.bvp_solver import (
    Parameters, solve_auto, d2fb_dphi2, dfb_dphi, fb
)

# ── Parameters ───────────────────────────────────────────────────────────────
p = Parameters(
    nb=1.0,
    chi_b=2.5,      # moderate Flory-Huggins bulk interaction (above spinodal threshold)
    omega_b=-0.2,   # surface attraction to bulk species
    kappa=1,
    phi_inf=0.15,    # far-field volume fraction (below spinodal)
    L=3.0,
    eps = 1e-6,
    R=5
)

# ── Free energy ───────────────────────────────────────────────────────────────

def F(phi, dphi, p):
    dr = p.L / len(phi)

    return (np.sum(dr * (fb(phi, p) - phi * dfb_dphi(phi, p) / 2) )
            + 2 * np.pi * p.R**2 * p.omega_b * phi[0])

# ── Thierry's formula ─────────────────────────────────────────────────────────

def phi_linear(z, p):
    alpha = d2fb_dphi2(p.phi_inf, p)

    h0 = p.omega_b
    lam = np.sqrt(alpha/p.kappa)

    phi = (p.phi_inf - h0 * p.R * np.exp(- lam * z) 
            / ((z + p.R) * p.kappa * (lam + 1 / p.R)))

    dphi = (phi - p.phi_inf) * (- lam - 1 / (z + p.R))

    return phi, dphi, lam

# ── Solve ─────────────────────────────────────────────────────────────────────
print("Solving BVP …")
sol = solve_auto(p, verbose=True)

if sol is None or not sol.success:
    print("Solver failed – try adjust parameters.")
    sys.exit(1)

nsteps = 1000
z = np.linspace(0, p.L, nsteps)
distances = np.linspace(0, p.L, nsteps)
phi, dphi = sol.sol(z)
phi_lin, dphi_lin, lam = phi_linear(z, p)
print(dfb_dphi(p.phi_inf, p))

# ── Compute potential ─────────────────────────────────────────────────────────

Fvals = np.zeros(nsteps)
Fvals_lin = np.zeros(nsteps)

Uvals = np.zeros(nsteps)
Uvals_lin = np.zeros(nsteps)

for i, d in enumerate(distances):
    # z = np.linspace(0, d, nsteps)
    u = np.abs(d - z)

    phi2, dphi2 = sol.sol(u)
    phi_lin2, dphi_lin2, lam = phi_linear(u, p)
    dphi_lin2 *= -np.sign(d-z)

    Uvals_lin[i] = (F(phi_lin + phi_lin2 - p.phi_inf, dphi_lin + dphi_lin2, p) 
                - F(phi_lin, dphi_lin, p) - F(phi_lin2, dphi_lin2, p))

    Uvals[i] = (F(phi + phi2 - p.phi_inf, dphi + dphi2, p) 
                - F(phi, dphi, p) - F(phi2, dphi2, p))
    
    Fvals[i] = F(phi, dphi, p)
    Fvals_lin[i] = F(phi_lin, dphi_lin, p)

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(2, 1, figsize=(5, 7))

minmax = np.where(np.diff(np.sign(np.diff(Fvals))))[0]

dFvals = np.diff(Fvals) * nsteps / p.L
d2Fvals = np.diff(dFvals) * nsteps / p.L

# ax.axvline(x=distances[minmax[0]], linestyle="dashed", color="black", lw=1)
# ax.axvline(x=distances[minmax[1]], linestyle="dashed", color="black", lw=1)

ax[0].plot(distances, Fvals)
ax[0].plot(distances, Fvals_lin)
# ax[0].axhline(y=- 2 * np.pi * p.R**2 * p.omega_b * phi[0], 
#               linestyle="dashed", color="black")

ax[1].plot(z, phi_lin, color="black", label="linearized")
ax[1].plot(z, phi, color="red", label="full")
ax[1].legend()

for d in [5]:
    u = np.abs(d - z)
    # z = np.linspace(0, d, int(d/p.L)*nsteps)
    # z = np.concatenate(z, np.zeros((1-int(d/p.L))*nsteps))
    # phi_lin2, dphi_lin2, lam = phi_linear(u, p)

    # ax[1].plot(z, phi_lin2, color="C{}".format(d))
    # ax[1].plot(z, phi_lin + phi_lin2 - p.phi_inf, color="C{}".format(d))

# ax.set_yscale('log')
# ax.set_xlim(10, None)

fig.tight_layout()
out = "media/figures/20260512/bvp_pair_potential.png"
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"\nSaved figure → {out}")
