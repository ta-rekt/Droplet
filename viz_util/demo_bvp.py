"""
Demo: run the BVP solver and plot the bulk concentration profile φ(z)
alongside the equilibrium membrane composition φ_m.

Output: media/figures/20260508/bvp_concentration_profile.png
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from viz_util.bvp_solver import (
    Parameters, solve_auto, compute_spinodal, compute_binodal,
    dfb_dphi, d2fb_dphi2, compute_mu_inf
)

# ── Parameters ───────────────────────────────────────────────────────────────
p = Parameters(
    nb=1.0,
    chi_b=2.5,      # moderate Flory-Huggins bulk interaction (above spinodal threshold)
    omega_b=-0.3,   # surface attraction to bulk species
    kappa=1,
    phi_inf=0.25,    # far-field volume fraction (below spinodal)
    L=10.0,
    R=50,
    conserve_volume=True
)

# ── Thierry's formula ─────────────────────────────────────────────────────────

def phi_profile(z, p):

    alpha = d2fb_dphi2(p.phi_inf, p)

    h0 = p.omega_b
    lam = np.sqrt(alpha/p.kappa)

    return (p.phi_inf - h0 * p.R * np.exp(- lam * z) 
            / ((z + p.R) * p.kappa * (lam + 1 / p.R)))

# ── Solve ─────────────────────────────────────────────────────────────────────
print("Solving BVP …")
sol = solve_auto(p, verbose=True)

if sol is None or not sol.success:
    print("Solver failed – try adjust parameters.")
    sys.exit(1)

z = np.linspace(0, p.L, 500)
if p.conserve_volume:
    phi, dphi, iphi = sol.sol(z)
else:
    phi, dphi = sol.sol(z)

print(f"\nSolution summary:")
print(f"  φ(0)  = {phi[0]:.4f}   (surface bulk composition)")
print(f"  φ(L)  = {phi[-1]:.4f}  (should ≈ φ_inf = {p.phi_inf})")
print("mu_infinity = {}".format(-dfb_dphi(p.phi_inf, p)))

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

# Left: concentration profile
ax = axes[0]
ax.plot(z, phi, color="#2563EB", lw=2, label=r"$\phi(z)$")
ax.plot(z, phi_profile(z, p), color="orange", linestyle="dashed", label="analytical")
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

# Right: gradient profile
ax2 = axes[1]
ax2.plot(z, dphi, color="#059669", lw=2)
ax2.axhline(0, color="gray", ls="--", lw=1)
ax2.set_xlabel(r"Distance from surface $z$", fontsize=12)
ax2.set_ylabel(r"$\partial_z\phi$", fontsize=12)
ax2.set_title("Composition gradient", fontsize=13)
ax2.set_xlim(0, p.L)
ax2.grid(True, alpha=0.3)

# Parameter annotation
param_text = (
    rf"$\chi_b={p.chi_b}$,"
    "\n"
    rf"$\omega_b={p.omega_b}$, $\mu_\infty={compute_mu_inf(p)}$"
    "\n"
    rf"$\kappa={p.kappa}$, $\phi_\infty={p.phi_inf}$, $L={p.L}$"
)
fig.text(0.5, -0.04, param_text, ha="center", fontsize=9, color="0.4")

fig.tight_layout()
out = "media/figures/20260509/bvp_concentration_profile_no_mem.png"
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"\nSaved figure → {out}")
