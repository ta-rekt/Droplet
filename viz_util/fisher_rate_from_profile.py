"""
Polymer prewetting as a concentration sensor 

A. linear approximation (analytical)
    1. plotting the equilibrium distance distribution between poly beads
    2. integrating over it for \psi_probe
    3. computing the fisher rate from \psi_probe

B. nonlinear profile from solver
    4. get free energy profile from numerical solution of \psi(z), Derjaguin
    5. repeat 1-3
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
import matplotlib
from dataclasses import dataclass
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from viz_util.bvp_solver import (
    Parameters, solve_auto, d2fb_dphi2, dfb_dphi, fb, compute_spinodal,
    compute_binodal, psi_linear, F, compute_mu_inf
)

@dataclass
class Probe:
    kb: float = 1.0
    ku: float = 1.0

# ── Parameters ────────────────────────────────────────────────────────────
p = Parameters(
    nb=1.0,
    chi_b=2.5,     
    omega_b=-0.4,   
    kappa=1,
    phi_inf=0.1,   
    L=20.0,
    eps=1e-6,
    R=2.0,
    conserve_volume=True,
    parallel_plates = False
)

# polymer params
s = 100   # inter-bead distance
a = 1   # Kuhn length

# ── Distance distribution ─────────────────────────────────────────────────
nsteps = 1000
r = np.linspace(0, p.L, nsteps)     # distance from center
d = np.linspace(0, p.L-p.R, nsteps) # distance from surface

def compute_P_eq(p: Parameters, d: np.ndarray) -> np.ndarray:
    """Compute the equilibrium distribution of bead separation on a polymer"""

    sol = solve_auto(p, verbose=False)
    
    # profile
    phi = sol.sol(d)[0]
    psi_lin = psi_linear(d, p)[0]

    # potential
    U_lin = 4 * np.pi * p.R**2 * p.omega_b * psi_lin

    # equilibrium distribution
    P_eq_unnorm = (d+p.R)**2 * np.exp(-3*(d+p.R)**2/(2*s*a**2) - U_lin)
    P_eq = P_eq_unnorm / np.sum(P_eq_unnorm)

    # concentration at the probe as a function of distance
    phi_boundary = phi[-1]
    Phi_d = phi_boundary + psi_lin[0] + psi_lin 

    return P_eq, Phi_d, U_lin


def compute_P_eq_full(p: Parameters, d: np.ndarray) -> np.ndarray:
    """Compute the equilibrium distribution of bead separation on a polymer
    using the nonlinear profile from the BVP solver"""

    sol = solve_auto(p, verbose=False)
    
    # profile
    phi = sol.sol(d)[0]
    psi_lin = psi_linear(d, p)[0]

    # potential
    U_lin = 4 * np.pi * p.R**2 * p.omega_b * psi_lin

    # equilibrium distribution
    P_eq_unnorm = (d+p.R)**2 * np.exp(-3*(d+p.R)**2/(2*s*a**2) - U_lin)
    P_eq = P_eq_unnorm / np.sum(P_eq_unnorm)

    # concentration at the probe as a function of distance
    phi_boundary = phi[-1]
    Phi_d = phi_boundary + psi_lin[0] + psi_lin 

    return P_eq, Phi_d, U_lin


# ── Mean probe concentration varying total phi volume ─────────────────────
n_params = 100
phi_probe = np.zeros(n_params)
phi_tf_vals = np.linspace(0.003, 0.15, n_params)
dr = (p.L-p.R) / nsteps

for i, _phi in enumerate(phi_tf_vals):
    p.phi_inf = _phi
    P_eq, Phi_d, U_lin = compute_P_eq(p, d)
    phi_probe[i] = np.sum(P_eq * Phi_d)


# ── Fisher rate ───────────────────────────────────────────────────────────
def fisher_rate(p2: Probe, 
                phi_probe: np.ndarray,
                phi_tf: np.ndarray) -> np.ndarray:
    
    # fisher rate without microcondensation, 2 state probe
    fisher_base = phi_tf * p2.kb * p2.ku / (p2.kb + p2.ku)

    # fisher rate accounting for probe-enhancer microcondensation
    dphi_probe = np.gradient(phi_probe, phi_tf)
    inv_taub = phi_probe * p2.kb * p2.ku / (p2.kb + p2.ku)
    fisher_probe = phi_tf**2 * inv_taub * (dphi_probe / phi_probe)**2

    return fisher_base, fisher_probe

probe = Probe(
    ku = 1.0, 
    kb = 1.0,
)
fisher_base, fisher_probe = fisher_rate(probe, phi_probe, phi_tf_vals)

# ── Plot ──────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(3,1, figsize=(5, 8.0))

p.phi_inf = 0.1
P_eq, Phi_d, U_lin = compute_P_eq(p, d)

# Left: concentration profile
axl = axes[0]
axl.set_title("Linear")
axl.set_ylabel(r"$U(d)$", color="black")
axl.plot(d, U_lin, color="black", lw=2, label=r"$U(d)$")
axl.set_xlabel("d (distance between surfaces)", fontsize=12)
axl.set_ylim(-1.10*np.max(np.abs(U_lin)), 1.10*np.max(np.abs(U_lin)))

axr = axl.twinx()
axr.set_ylabel(r"$P_{eq}(d)$", color="#2563EB")
axr.plot(d, P_eq, color="#2563EB", lw=2, label=r"$\psi(d)$")
axr.set_ylim(-1.10*np.max(P_eq), 1.10*np.max(P_eq))

axl = axes[1]
axl.set_ylabel(r"$\phi_{probe}$", color="black")
axl.plot(phi_tf_vals, phi_probe, color="C2", lw=2)
axl.plot(phi_tf_vals, phi_tf_vals, color="C2", lw=1.5, linestyle="dashed")
axl.set_xlabel(r"$\phi_{TF}$ (total volume fraction)", fontsize=12)

axl = axes[2]
axl.set_ylabel(r"Fisher information rate", color="black")
axl.plot(phi_tf_vals, fisher_base, color="red", lw=1.5, label="control",linestyle="dashed")
axl.plot(phi_tf_vals, fisher_probe, color="red", lw=2, label="with interaction")
axl.set_xlabel(r"$\phi_{TF}$ (total volume fraction)", fontsize=12)
axl.set_yscale('log')
axl.legend()

# ── Save ──────────────────────────────────────────────────────────────────
fig.tight_layout()
out = "media/figures/20260516/fisher_rate_from_profile.png"
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"\nSaved figure → {out}")
