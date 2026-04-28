import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp

# Lattice sizes / molecular sizes (paper uses nb, nm)  (Eq. 15a, 15b)  :contentReference[oaicite:6]{index=6}
nb = 1.0   # bulk molecule volume in units of solvent molecular volume
nm = 1.0   # membrane molecule area in units of solvent molecular area

# Bulk and membrane interaction parameters (Eq. 15a, 15b) :contentReference[oaicite:7]{index=7}
chi_b = 3
chi_m = 2.5
omega_m = -0.3

# Coupling free energy parameters (Eq. 15c) :contentReference[oaicite:8]{index=8}
omega_b = -0.3
chi_bb = -0.5
chi_bm = -1.0

# Gradient penalty in bulk (paper's kappa in Eq. 3 / Eq. 9) :contentReference[oaicite:9]{index=9}
kappa = 1.0

# Far-field bulk composition (reservoir)
phi_inf = 0.06  # must be in (0,1) for FH logs

# Domain truncation: solve on z in [0, L]
L = 10.0

# ============================================================
# Parameters (dimensionless units: set kBT=1, solvent nu=1, nu~ = 1)
# ============================================================

# Numerical stability for logs
EPS = 1e-10

# ============================================================
# Free energies and derivatives (Eq. 15a–15c) :contentReference[oaicite:10]{index=10}
# ============================================================

def clip01(x):
    return np.clip(x, EPS, 1.0 - EPS)

def dfb_dphi(phi):
    """
    d/dphi of bulk free-energy density fb(phi):
    fb = (1/nb)[phi ln phi + (1-phi) ln(1-phi)] + chi_b phi(1-phi)
    (overall factor kBT/nu set to 1)
    """
    phi = clip01(phi)
    return (1.0/nb) * (np.log(phi) - np.log(1.0 - phi)) + chi_b * (1.0 - 2.0*phi)

def dfm_dphim(phim):
    """
    d/dphim of membrane free-energy density fm(phim):
    fm = (1/nm)[phim ln phim + (1-phim) ln(1-phim)] + chi_m phim(1-phim) + omega_m phim
    (overall factor kBT/nu~ set to 1)
    """
    phim = clip01(phim)
    return (1.0/nm) * (np.log(phim) - np.log(1.0 - phim)) + chi_m * (1.0 - 2.0*phim) + omega_m

def J(phi, phim):
    """
    Coupling free energy J(phi, phim) at the surface (Eq. 15c) :contentReference[oaicite:11]{index=11}
    J = omega_b*phi + chi_bb*phi^2 + chi_bm*phi*phim
    """
    return omega_b*phi + chi_bb*phi**2 + chi_bm*phi*phim

def dJ_dphi(phi, phim):
    return omega_b + 2.0*chi_bb*phi + chi_bm*phim

def dJ_dphim(phi, phim):
    return chi_bm*phi

# ============================================================
# Chemical potential of the reservoir μ∞
# Paper: μ = νb δF/δφ, and Eq. 9a uses μ∞/νb. :contentReference[oaicite:12]{index=12}
# With our fb in "per-volume" units, μ∞ = νb * dfb/dphi(phi_inf)
# and νb = nb (since solvent nu=1 and molecule volume = nb*nu).
# ============================================================

nu_b = nb
nu_m = nm
mu_inf = nu_b * dfb_dphi(phi_inf)

# ============================================================
# BVP formulation
# y[0]=phi, y[1]=dphi/dz; parameter p[0]=phim (constant)
# Bulk ODE from Eq. 9a: kappa*phi'' = dfb/dphi - mu_inf/nu_b :contentReference[oaicite:13]{index=13}
# ============================================================

def ode(z, y, p):
    phi = y[0]
    dphi = y[1]
    phim = p[0]

    # second derivative from rearranged Eq. 9a
    ddphi = (dfb_dphi(phi) - (mu_inf / nu_b)) / kappa
    return np.vstack((dphi, ddphi))

# Boundary conditions:
# 1) Eq. 9c: -kappa*phi'(0) + dJ/dphi(phi(0),phim)=0 :contentReference[oaicite:14]{index=14}
# 2) Eq. 9b at z=L: phi(L)=phi_inf :contentReference[oaicite:15]{index=15}
# 3) Eq. 9d: dfm/dphim + dJ/dphim - mu_inf/nu_m = 0 :contentReference[oaicite:16]{index=16}
def bc(ya, yb, p):
    phim = p[0]
    phi0 = ya[0]
    dphi0 = ya[1]
    return np.array([
        -kappa*dphi0 + dJ_dphi(phi0, phim),             # Eq. 9c
        yb[0] - phi_inf,                                 # Eq. 9b at z=L
        dfm_dphim(phim) + dJ_dphim(phi0, phim) - (mu_inf/nu_m)  # Eq. 9d
    ])

# ============================================================
# Initial guess (important for convergence)
# ============================================================

z = np.linspace(0.0, L, 400)

# Start with a mild enrichment near the surface
phi0_guess = min(0.95, phi_inf + 0.1)
ell = 5.0
phi_guess = phi_inf + (phi0_guess - phi_inf) * np.exp(-z/ell)
dphi_guess = -(phi0_guess - phi_inf)/ell * np.exp(-z/ell)

y_guess = np.vstack((phi_guess, dphi_guess))

# Initial guess for phim (must be in (0,1))
p_guess = np.array([0.2])

# ============================================================
# Solve
# ============================================================

sol = solve_bvp(ode, bc, z, y_guess, p=p_guess, tol=1e-6, max_nodes=50000)

print("status:", sol.status, sol.message)
print("phi_m =", sol.p[0])

# Evaluate solution
z_fine = np.linspace(0.0, L, 800)
phi_fine, dphi_fine = sol.sol(z_fine)

# ============================================================
# Plot
# ============================================================

plt.figure()
plt.plot(z_fine, phi_fine, label=r'$\phi(z)$')
plt.axhline(phi_inf, linestyle='--', label=r'$\phi_\infty$')
plt.xlabel("z")
plt.ylabel("phi")
plt.legend()
plt.tight_layout()
plt.ylim((0,1))
plt.show()
