"""
2D axisymmetric BVP solver for two spheres in a Flory-Huggins bulk.

Two identical spheres of radius R sit on the z-axis at z = 0 and z = 2R + d,
so the surface-to-surface separation is d.  The equilibrium volume fraction
phi(rho, z) satisfies

    kappa * Laplacian(phi) = dfb/dphi - mu_inf / nb

on the cylindrical (rho, z) domain with axial symmetry, where the Laplacian
is the cylindrical (rho, z) one:

    Laplacian phi = d2 phi/drho2 + (1/rho) d phi/drho + d2 phi/dz2

Boundary conditions:
  - Robin (adsorption) on each sphere:  -kappa * d phi/dn + omega_b = 0
    with n the outward normal pointing away from the sphere centre.
  - Dirichlet phi = phi_inf on the outer cylinder rho = rho_max and the z
    end-caps z = -z_buf, z = 2R + d + z_buf.  This is the "absorbing"
    truncation of the far-field reservoir.
  - dphi/drho = 0 on the symmetry axis rho = 0.

Sphere interiors are masked out of the FD grid; their surfaces are handled
with a ghost-cell immersed-boundary method (linear extrapolation from the
mirror-image fluid point using bilinear interpolation, with the slope at
the surface pinned by the Robin BC).

Nonlinear solve: scipy.optimize.newton_krylov, warm-started from the linear-
response superposition of two single-sphere psi_linear profiles.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from dataclasses import dataclass, replace
from typing import Optional

import numpy as np
from scipy.optimize import newton_krylov, NoConvergence

from viz_util.bvp_solver import (
    Parameters, dfb_dphi, d2fb_dphi2, fb, compute_mu_inf, clip01, psi_linear,
)


# ------------------------------------------------------------------
# Geometry and node classification
# ------------------------------------------------------------------

@dataclass
class Geometry:
    p: Parameters
    d: float
    Nrho: int
    Nz: int
    rho_grid: np.ndarray
    z_grid: np.ndarray
    drho: float
    dz: float
    z_c1: float
    z_c2: float
    inside_mask: np.ndarray
    fluid_mask: np.ndarray
    dir_mask: np.ndarray
    unk_mask: np.ndarray
    ghost_ij: np.ndarray
    ghost_i0j0: np.ndarray
    ghost_wts: np.ndarray
    ghost_2delta: np.ndarray


def _decay_length(p: Parameters) -> float:
    d2 = float(d2fb_dphi2(np.asarray([p.phi_inf]), p)[0])
    if d2 > 0.05:
        return float(np.sqrt(p.kappa / d2))
    return float(np.sqrt(p.kappa * p.nb / 0.05))


def build_geometry(p: Parameters, d: float,
                   rho_max: Optional[float] = None,
                   z_buf:   Optional[float] = None,
                   Nrho: int = 80, Nz: int = 200) -> Geometry:
    R = p.R
    ell = _decay_length(p)
    if rho_max is None:
        rho_max = R + 4.0 * ell
    if z_buf is None:
        z_buf = 4.0 * ell

    z_c1 = 0.0
    z_c2 = 2.0 * R + d

    rho_grid = np.linspace(0.0, rho_max, Nrho)
    # z-box bounds include a buffer beyond each sphere's pole (not just past
    # each centre), so the spheres are fully enclosed with z_buf of clear
    # reservoir on either side.
    z_grid   = np.linspace(z_c1 - R - z_buf, z_c2 + R + z_buf, Nz)
    drho = float(rho_grid[1] - rho_grid[0])
    dz   = float(z_grid[1]   - z_grid[0])

    Rho, Z = np.meshgrid(rho_grid, z_grid, indexing="ij")
    r1_sq = Rho ** 2 + (Z - z_c1) ** 2
    r2_sq = Rho ** 2 + (Z - z_c2) ** 2
    inside_mask = (r1_sq < R ** 2) | (r2_sq < R ** 2)
    fluid_mask  = ~inside_mask

    dir_mask = np.zeros_like(fluid_mask)
    dir_mask[-1, :] = True
    dir_mask[:,  0] = True
    dir_mask[:, -1] = True
    dir_mask &= fluid_mask
    unk_mask = fluid_mask & ~dir_mask

    # Ghost cells: inside-sphere nodes with at least one 4-connected fluid
    # neighbour.  Each gets a bilinear-interp prescription that fills its
    # value from the mirror-image fluid point.
    ghost_records = []
    nbr_offsets = ((1, 0), (-1, 0), (0, 1), (0, -1))
    ii_inside, jj_inside = np.where(inside_mask)
    for i, j in zip(ii_inside, jj_inside):
        has_fluid = False
        for di, dj in nbr_offsets:
            ii, jj = i + di, j + dj
            if 0 <= ii < Nrho and 0 <= jj < Nz and fluid_mask[ii, jj]:
                has_fluid = True
                break
        if not has_fluid:
            continue
        rho_G = rho_grid[i]
        z_G   = z_grid[j]
        d1 = np.hypot(rho_G, z_G - z_c1)
        d2 = np.hypot(rho_G, z_G - z_c2)
        if d1 <= d2:
            zc, Gc = z_c1, d1
        else:
            zc, Gc = z_c2, d2
        if Gc < 1e-12:
            continue
        u_rho = rho_G / Gc
        u_z   = (z_G - zc) / Gc
        delta = R - Gc
        # Image point: start at the mirror across the surface (distance
        # 2*delta from G).  If the natural image's bilinear stencil falls on
        # any inside-sphere cell, push the image further out along the normal
        # until the stencil is clean.  We then use the forward-FD relation
        # phi(G) = phi(I) - d_GI * omega_b/kappa, where d_GI = |G - I|.
        step = max(drho, dz)
        d_GI = 2.0 * delta
        rho_I = rho_G + d_GI * u_rho
        z_I   = z_G   + d_GI * u_z
        chosen_i0 = chosen_j0 = 0
        chosen_w = np.zeros(4)
        for _ in range(40):
            rho_I_c = float(np.clip(rho_I, rho_grid[0], rho_grid[-1] - 1e-12 * drho))
            z_I_c   = float(np.clip(z_I,   z_grid[0],   z_grid[-1]   - 1e-12 * dz))
            ti = (rho_I_c - rho_grid[0]) / drho
            tj = (z_I_c   - z_grid[0])   / dz
            i0 = int(np.clip(int(np.floor(ti)), 0, Nrho - 2))
            j0 = int(np.clip(int(np.floor(tj)), 0, Nz   - 2))
            t_rho = ti - i0
            t_z   = tj - j0
            stencil_clean = not (inside_mask[i0, j0] or inside_mask[i0 + 1, j0]
                                 or inside_mask[i0, j0 + 1]
                                 or inside_mask[i0 + 1, j0 + 1])
            chosen_i0, chosen_j0 = i0, j0
            chosen_w = np.array([
                (1.0 - t_rho) * (1.0 - t_z),
                t_rho         * (1.0 - t_z),
                (1.0 - t_rho) * t_z,
                t_rho         * t_z,
            ])
            if stencil_clean:
                break
            d_GI += step
            rho_I = rho_G + d_GI * u_rho
            z_I   = z_G   + d_GI * u_z
        ghost_records.append((i, j, chosen_i0, chosen_j0, chosen_w, d_GI))

    if ghost_records:
        ghost_ij     = np.array([(r[0], r[1]) for r in ghost_records], dtype=int)
        ghost_i0j0   = np.array([(r[2], r[3]) for r in ghost_records], dtype=int)
        ghost_wts    = np.array([r[4] for r in ghost_records], dtype=float)
        ghost_2delta = np.array([r[5] for r in ghost_records], dtype=float)
    else:
        ghost_ij     = np.zeros((0, 2), dtype=int)
        ghost_i0j0   = np.zeros((0, 2), dtype=int)
        ghost_wts    = np.zeros((0, 4), dtype=float)
        ghost_2delta = np.zeros((0,), dtype=float)

    return Geometry(
        p=p, d=d, Nrho=Nrho, Nz=Nz,
        rho_grid=rho_grid, z_grid=z_grid,
        drho=drho, dz=dz, z_c1=z_c1, z_c2=z_c2,
        inside_mask=inside_mask, fluid_mask=fluid_mask,
        dir_mask=dir_mask, unk_mask=unk_mask,
        ghost_ij=ghost_ij, ghost_i0j0=ghost_i0j0,
        ghost_wts=ghost_wts, ghost_2delta=ghost_2delta,
    )


# ------------------------------------------------------------------
# Scatter (assemble the full phi field including Dirichlet + ghost fills)
# ------------------------------------------------------------------

def _scatter(x_unk: np.ndarray, geom: Geometry) -> np.ndarray:
    p = geom.p
    phi = np.full((geom.Nrho, geom.Nz), p.phi_inf)
    phi[geom.unk_mask] = x_unk
    if geom.ghost_ij.shape[0] > 0:
        i0 = geom.ghost_i0j0[:, 0]
        j0 = geom.ghost_i0j0[:, 1]
        gi = geom.ghost_ij[:, 0]
        gj = geom.ghost_ij[:, 1]
        # By construction the bilinear stencil of each ghost's image lies in
        # the fluid (or on the Dirichlet boundary), so a single pass suffices.
        phi_I = (geom.ghost_wts[:, 0] * phi[i0,     j0]
               + geom.ghost_wts[:, 1] * phi[i0 + 1, j0]
               + geom.ghost_wts[:, 2] * phi[i0,     j0 + 1]
               + geom.ghost_wts[:, 3] * phi[i0 + 1, j0 + 1])
        phi[gi, gj] = phi_I - geom.ghost_2delta * p.omega_b / p.kappa
    return phi


# ------------------------------------------------------------------
# Cylindrical Laplacian on the structured grid
# ------------------------------------------------------------------

def _laplacian(phi: np.ndarray, geom: Geometry) -> np.ndarray:
    drho, dz = geom.drho, geom.dz
    lap = np.zeros_like(phi)
    rho_int = geom.rho_grid[1:-1, None]
    lap[1:-1, 1:-1] = (
        (phi[2:, 1:-1] - 2.0 * phi[1:-1, 1:-1] + phi[:-2, 1:-1]) / drho ** 2
        + (phi[2:, 1:-1] - phi[:-2, 1:-1]) / (2.0 * drho * rho_int)
        + (phi[1:-1, 2:] - 2.0 * phi[1:-1, 1:-1] + phi[1:-1, :-2]) / dz ** 2
    )
    # Axis: dphi/drho = 0 by symmetry; L'Hopital gives (1/rho)dphi/drho ->
    # d2phi/drho2, so the axial Laplacian is 2 d2phi/drho2 + d2phi/dz2.
    # With the symmetric ghost phi[-1, j] = phi[1, j], 2 d2phi/drho2 reduces
    # to 4 (phi[1, j] - phi[0, j]) / drho2.
    lap[0, 1:-1] = (
        4.0 * (phi[1, 1:-1] - phi[0, 1:-1]) / drho ** 2
        + (phi[0, 2:] - 2.0 * phi[0, 1:-1] + phi[0, :-2]) / dz ** 2
    )
    return lap


# ------------------------------------------------------------------
# Residual + solve
# ------------------------------------------------------------------

def _make_residual(geom: Geometry, mu_inf: float):
    p = geom.p
    unk = geom.unk_mask

    def residual(x_unk: np.ndarray) -> np.ndarray:
        phi = _scatter(x_unk, geom)
        lap = _laplacian(phi, geom)
        driving = dfb_dphi(phi, p) - mu_inf / p.nb
        return (p.kappa * lap - driving)[unk]

    return residual


def _initial_guess(geom: Geometry) -> np.ndarray:
    p = geom.p
    Rho, Z = np.meshgrid(geom.rho_grid, geom.z_grid, indexing="ij")
    s1 = np.sqrt(Rho ** 2 + (Z - geom.z_c1) ** 2) - p.R
    s2 = np.sqrt(Rho ** 2 + (Z - geom.z_c2) ** 2) - p.R
    s1c = np.clip(s1, 0.0, None)
    s2c = np.clip(s2, 0.0, None)
    psi1, _, _ = psi_linear(s1c, p)
    psi2, _, _ = psi_linear(s2c, p)
    phi0 = p.phi_inf + psi1 + psi2
    phi0 = clip01(phi0, p.eps)
    return phi0[geom.unk_mask]


def solve_two_spheres(p: Parameters, d: float,
                      rho_max: Optional[float] = None,
                      z_buf:   Optional[float] = None,
                      Nrho: int = 80, Nz: int = 200,
                      f_tol: float = 1e-6,
                      maxiter: int = 200,
                      verbose: bool = False,
                      omega_continuation: bool = True) -> dict:
    """
    Solve the two-sphere axisymmetric BVP.

    Returns dict with keys:
      'phi'      : (Nrho, Nz) array including ghost-filled inside cells
                   (and phi_inf in non-ghost inside cells)
      'geom'     : Geometry struct (grids, masks, etc.)
      'p', 'd'   : echo back inputs
      'mu_inf'   : reservoir chemical potential
      'success'  : True if newton_krylov converged
      'x'        : flat solution vector (unknowns only)
    """
    geom   = build_geometry(p, d, rho_max=rho_max, z_buf=z_buf,
                            Nrho=Nrho, Nz=Nz)
    mu_inf = compute_mu_inf(p)

    x0 = _initial_guess(geom)
    success = True
    try:
        x_sol = newton_krylov(_make_residual(geom, mu_inf), x0,
                              f_tol=f_tol, maxiter=maxiter, verbose=verbose)
    except NoConvergence as e:
        x_sol = np.asarray(e.args[0])
        success = False

    # If the cold start failed and continuation is on, ramp omega_b from a
    # small value up to the target, reusing the previous solution each step.
    if not success and omega_continuation and p.omega_b != 0.0:
        if verbose:
            print("Cold start failed; trying omega_b continuation ...")
        n_steps = 6
        x_cur = _initial_guess(build_geometry(replace(p, omega_b=0.0), d,
                                              rho_max=rho_max, z_buf=z_buf,
                                              Nrho=Nrho, Nz=Nz))
        for k in range(1, n_steps + 1):
            frac = k / n_steps
            p_k = replace(p, omega_b=frac * p.omega_b)
            geom_k = build_geometry(p_k, d, rho_max=rho_max, z_buf=z_buf,
                                    Nrho=Nrho, Nz=Nz)
            try:
                x_cur = newton_krylov(
                    _make_residual(geom_k, compute_mu_inf(p_k)),
                    x_cur, f_tol=f_tol, maxiter=maxiter, verbose=False,
                )
            except NoConvergence as e:
                x_cur = np.asarray(e.args[0])
                if verbose:
                    print(f"  continuation step {k}/{n_steps} (frac={frac:.2f}) "
                          f"did not converge")
                break
            if verbose:
                print(f"  continuation step {k}/{n_steps} (frac={frac:.2f}) ok")
        x_sol = x_cur
        geom = geom_k
        mu_inf = compute_mu_inf(p)
        # Final residual check
        res = _make_residual(geom, mu_inf)(x_sol)
        success = float(np.max(np.abs(res))) < 10.0 * f_tol

    phi = _scatter(x_sol, geom)
    return dict(phi=phi, geom=geom, p=p, d=d, mu_inf=mu_inf,
                success=success, x=x_sol)


# ------------------------------------------------------------------
# Bilinear sampling utility (used by the surface integral)
# ------------------------------------------------------------------

def _interp_phi(phi: np.ndarray, geom: Geometry,
                rho_q: np.ndarray, z_q: np.ndarray) -> np.ndarray:
    rho_g, z_g = geom.rho_grid, geom.z_grid
    drho, dz = geom.drho, geom.dz
    rho_c = np.clip(rho_q, rho_g[0], rho_g[-1] - 1e-12 * drho)
    z_c   = np.clip(z_q,   z_g[0],   z_g[-1]   - 1e-12 * dz)
    ti = (rho_c - rho_g[0]) / drho
    tj = (z_c   - z_g[0])   / dz
    i0 = np.clip(np.floor(ti).astype(int), 0, geom.Nrho - 2)
    j0 = np.clip(np.floor(tj).astype(int), 0, geom.Nz   - 2)
    t_rho = ti - i0
    t_z   = tj - j0
    p00 = phi[i0,     j0]
    p10 = phi[i0 + 1, j0]
    p01 = phi[i0,     j0 + 1]
    p11 = phi[i0 + 1, j0 + 1]
    return ((1.0 - t_rho) * (1.0 - t_z) * p00
          + t_rho         * (1.0 - t_z) * p10
          + (1.0 - t_rho) * t_z         * p01
          + t_rho         * t_z         * p11)


# ------------------------------------------------------------------
# Excess free energy
# ------------------------------------------------------------------

def excess_free_energy_two_spheres(result: dict) -> float:
    """
    Grand-potential excess over the bulk reservoir at phi_inf:

      Delta Omega
        = int_V [ df(phi) + (kappa/2) |grad phi|^2 ] 2 pi rho drho dz
        + omega_b * sum_{both spheres} int_{sphere surface} (phi - phi_inf) dA

    where df(phi) = fb(phi) - fb(phi_inf) - mu_ref/nb (phi - phi_inf) and the
    surface integral uses spherical-polar quadrature dA = R^2 sin theta
    dtheta dphi (closed form in phi gives 2 pi R^2).  Mirrors the open-mode
    branch of excess_free_energy in bvp_solver.py.
    """
    geom = result["geom"]
    phi  = result["phi"]
    p    = geom.p
    phi_inf = p.phi_inf

    f_inf  = float(fb(np.asarray([phi_inf]), p)[0])
    mu_ref = p.nb * float(dfb_dphi(np.asarray([phi_inf]), p)[0])

    # Use the ghost-filled phi for the gradient: this gives smooth FD across
    # the sphere surface (the ghost values encode the BC slope).  Mask the
    # integrand to zero on non-fluid cells AFTER computing it.
    df = (fb(phi, p) - f_inf - mu_ref / p.nb * (phi - phi_inf))
    dphi_drho = np.gradient(phi, geom.drho, axis=0)
    dphi_dz   = np.gradient(phi, geom.dz,   axis=1)
    grad_sq = dphi_drho ** 2 + dphi_dz ** 2
    integrand = (df + 0.5 * p.kappa * grad_sq)
    integrand = np.where(geom.fluid_mask, integrand, 0.0)

    weight = 2.0 * np.pi * geom.rho_grid[:, None]
    bulk = float(np.trapezoid(np.trapezoid(weight * integrand,
                                           geom.z_grid, axis=1),
                              geom.rho_grid, axis=0))

    # Surface integral.  Parametrise each sphere by polar angle theta in [0, pi];
    # symmetric in the azimuthal angle (factor 2 pi).
    R = p.R
    theta = np.linspace(0.0, np.pi, 400)
    sin_t = np.sin(theta)
    cos_t = np.cos(theta)
    rho_s = R * sin_t

    surf_total = 0.0
    for zc in (geom.z_c1, geom.z_c2):
        z_s = zc + R * cos_t
        phi_s = _interp_phi(phi, geom, rho_s, z_s)
        surf_total += float(np.trapezoid((phi_s - phi_inf) * sin_t, theta))
    surf = p.omega_b * 2.0 * np.pi * R ** 2 * surf_total
    return bulk + surf


# ------------------------------------------------------------------
# __main__: verification tests + U(d) plot
# ------------------------------------------------------------------

def _run_main():
    import datetime
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from viz_util.fisher_rate_from_profile import (
        solve_isolated_sphere, total_U, L_ISO,
    )

    # Parameters at the spinodal threshold chi_b = 2 with a noticeably
    # attractive surface.  Strong enough to produce a clear pair potential
    # but not so strong (chi_b > 2 with large omega_b) that the cold-start
    # initial guess fails to converge to the wetting-layer branch.
    p = Parameters(
        nb=1.0, chi_b=2.0, omega_b=-0.5, kappa=1.0,
        phi_inf=0.1, L=20.0, R=5.0, eps=1e-6,
        conserve_volume=False, parallel_plates=False,
    )

    # ------------------------------------------------------------------
    # Test 1: far-separation matches isolated single-sphere on axis
    # ------------------------------------------------------------------
    print("\n[1] Far-separation comparison (d = 4 R) ...")
    d_far = 4.0 * p.R
    res = solve_two_spheres(p, d_far, Nrho=60, Nz=200, verbose=False)
    print(f"  converged: {res['success']}")
    z_axis   = res["geom"].z_grid
    phi_axis = res["phi"][0, :]
    psi_i, _, s_max = solve_isolated_sphere(p)
    # Region above sphere 1 (z in (R, R + d_far/2))
    region = (z_axis > p.R + 1e-3) & (z_axis < p.R + d_far / 2.0)
    s_ref  = np.clip(z_axis[region] - p.R, 0.0, s_max)
    phi_ref = p.phi_inf + psi_i(s_ref)
    err = float(np.max(np.abs(phi_axis[region] - phi_ref)))
    print(f"  max |two-sphere - single-sphere| on axis, near sphere 1: {err:.3e}")

    # ------------------------------------------------------------------
    # Test 2: linear-response superposition (small omega_b)
    # ------------------------------------------------------------------
    print("\n[2] Linear-response limit (omega_b = -0.05, chi_b = 0) ...")
    p_lin = replace(p, chi_b=0.0, omega_b=-0.05)
    d_lin = 2.0 * p_lin.R
    res_lin = solve_two_spheres(p_lin, d_lin, Nrho=60, Nz=200, verbose=False)
    print(f"  converged: {res_lin['success']}")
    geom = res_lin["geom"]
    Rho, Z = np.meshgrid(geom.rho_grid, geom.z_grid, indexing="ij")
    s1 = np.clip(np.sqrt(Rho ** 2 + (Z - geom.z_c1) ** 2) - p_lin.R, 0.0, None)
    s2 = np.clip(np.sqrt(Rho ** 2 + (Z - geom.z_c2) ** 2) - p_lin.R, 0.0, None)
    psi1 = psi_linear(s1, p_lin)[0]
    psi2 = psi_linear(s2, p_lin)[0]
    phi_lin_pred = p_lin.phi_inf + psi1 + psi2
    diff = np.abs(res_lin["phi"] - phi_lin_pred)
    diff = np.where(geom.fluid_mask, diff, 0.0)
    err_lin = float(np.max(diff))
    print(f"  max |numeric - linear superposition|: {err_lin:.3e}  "
          f"(expected O(omega_b^2) = {p_lin.omega_b ** 2:.3e})")

    # ------------------------------------------------------------------
    # Test 3: U(d) sweep vs Derjaguin / additive
    # ------------------------------------------------------------------
    # Reference F at a large separation: U(d) = F(d) - F(d_ref).  This
    # cancels the systematic discretisation error of the 2D FD method and
    # gives an interaction energy that vanishes at large d by construction.
    print("\n[3] U(d) sweep vs Derjaguin / additive references ...")
    d_grid = np.linspace(0.3 * p.R, 2.5 * p.R, 8)
    Nrho = 100
    dz_target = 0.10
    def _Nz(dd):
        return int((2 * p.R + dd + 8.0) / dz_target) + 1

    d_ref = 8.0 * p.R
    res_ref = solve_two_spheres(p, d_ref, Nrho=Nrho, Nz=_Nz(d_ref),
                                z_buf=4.0, rho_max=10.0, verbose=False)
    F_ref = excess_free_energy_two_spheres(res_ref)
    print(f"  reference F(d={d_ref}) = {F_ref:.4e}")
    U_full = np.full(d_grid.shape, np.nan)
    for k, dd in enumerate(d_grid):
        res = solve_two_spheres(p, float(dd), Nrho=Nrho, Nz=_Nz(dd),
                                z_buf=4.0, rho_max=10.0, verbose=False)
        if not res["success"]:
            print(f"  d={dd:.3f}: NOT converged, skipping")
            continue
        F = excess_free_energy_two_spheres(res)
        U_full[k] = F - F_ref
        print(f"  d={dd:.3f}: F={F:.4e}  U={U_full[k]:+.4e}")

    print("  Computing Derjaguin + additive reference (may take a moment) ...")
    try:
        U_blend, U_D, U_A, _ = total_U(d_grid, p)
    except Exception as exc:
        print(f"  reference failed ({exc}); plotting only U_full")
        U_blend = U_D = U_A = np.full_like(d_grid, np.nan)

    today = datetime.date.today().strftime("%Y%m%d")
    out_dir = f"media/figures/{today}"
    os.makedirs(out_dir, exist_ok=True)
    out = f"{out_dir}/bvp_two_spheres_U_d.png"

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    ax = axes[0]
    ax.plot(d_grid, U_full, "o-",  color="C0",       label="full 2D BVP",
            lw=2.0)
    ax.plot(d_grid, U_D,    "s--", color="#9333EA",  label="Derjaguin",
            alpha=0.7)
    ax.plot(d_grid, U_A,    "d:",  color="#10B981",  label="additive",
            alpha=0.7)
    ax.plot(d_grid, U_blend, "-",  color="#2563EB",  label="blended (current)",
            alpha=0.8)
    ax.set_xlabel("d (surface-to-surface)")
    ax.set_ylabel("U(d)")
    ax.set_title(f"chi_b={p.chi_b}, omega_b={p.omega_b}, phi_inf={p.phi_inf}")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # Right panel: midplane phi(rho) at d = R
    ax = axes[1]
    res_mid = solve_two_spheres(p, p.R, Nrho=Nrho, Nz=_Nz(p.R),
                                z_buf=4.0, rho_max=10.0, verbose=False)
    geom = res_mid["geom"]
    j_mid = np.argmin(np.abs(geom.z_grid - (geom.z_c1 + geom.z_c2) / 2.0))
    ax.plot(geom.rho_grid, res_mid["phi"][:, j_mid], "o-", color="C2",
            label="phi(rho, z=midplane) at d=R")
    ax.axhline(p.phi_inf, ls="--", color="gray", alpha=0.6, label="phi_inf")
    ax.set_xlabel("rho")
    ax.set_ylabel("phi")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"\nSaved -> {out}")


if __name__ == "__main__":
    _run_main()
