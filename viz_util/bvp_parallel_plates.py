"""
Plot: compute the pair potential using the linear superposition approx.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from viz_util.bvp_solver import (
    Parameters, solve_auto, d2fb_dphi2, dfb_dphi, fb, compute_spinodal,
    compute_binodal, psi_linear, F, compute_mu_inf
)

# ── Parameters ───────────────────────────────────────────────────────────────
p = Parameters(
    nb=1.0,
    chi_b=2.5,      # moderate Flory-Huggins bulk interaction (above spinodal threshold)
    omega_b=-1,   # surface attraction to bulk species
    kappa=1,
    phi_inf=0.25,    # far-field volume fraction (below spinodal)
    L=5.0,
    eps = 1e-6,
    R=np.inf,
    conserve_volume=True,
    parallel_plates = True
)


# ── Solve ─────────────────────────────────────────────────────────────────────
print("Solving BVP …")
sol = solve_auto(p, verbose=True)

if sol is None or not sol.success:
    print("Solver failed – try adjust parameters.")
    sys.exit(1)

nsteps = 1000
z = np.linspace(0, p.L, nsteps)
distances = np.linspace(0, p.L, nsteps)

# full solution, left sphere
phi, dphi, iphi = sol.sol(z)

print(p.is_spherical)

# ── Compute potential ─────────────────────────────────────────────────────────

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1,1, figsize=(5, 3.5))

# Left: concentration profile
ax = axes
ax.plot(z, phi, color="#2563EB", lw=2, label=r"$\phi(z)$")
# ax.set_yscale('log')

ax.axhline(p.phi_inf, color="gray", ls="--", lw=1, label=r"$\phi_\infty$")
phi_spin_low, phi_spin_high = compute_spinodal(p)
if phi_spin_low is not None:
    ax.axhspan(phi_spin_low, phi_spin_high, alpha=0.10, color="red",
               label=f"spinodal ({phi_spin_low:.2f}, {phi_spin_high:.2f})")
phi_bin_low, phi_bin_high = compute_binodal(p)
if phi_bin_low is not None:
    ax.axhline(phi_bin_low, color="purple", ls="-.", lw=1.2,
               label=f"binodal ({phi_bin_low:.2f}, {phi_bin_high:.2f})")
    ax.axhline(phi_bin_high, color="purple", ls="-.", lw=1.2)
ax.axvline(0, color="k", lw=0.8, ls=":")
ax.scatter([0], np.array([phi[0]]), color="#2563EB", zorder=5, s=50)
ax.set_xlabel(r"Distance from surface $z$", fontsize=12)
ax.set_ylabel(r"Volume fraction $\phi$", fontsize=12)
ax.set_title("Bulk concentration profile", fontsize=13)
ax.set_xlim(0, p.L)
# ax.set_ylim(0.0, 0.2)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)


fig.tight_layout()
out = "media/figures/20260515/bvp_parallel_plates.png"
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"\nSaved figure → {out}")
