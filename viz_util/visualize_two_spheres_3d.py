"""
3D visualization of the two-sphere phi(rho, z) profile produced by
spec/fisher_two_spheres.jl.

The spec saves the full 2D axisymmetric phi field per converged d under
res["profile_phi"][k] with companion res["profile_rho"][k],
res["profile_z"][k], res["profile_d"][k].  This script:

  1. Reads one (alpha, beta, d-index) tuple from the JSON output.
  2. Rotates the 2D phi(rho, z) field around the z-axis to recover the full
     3D axisymmetric volume.
  3. Renders the two spheres as semi-transparent meshes plus a translucent
     point cloud sampled with density proportional to (phi - phi_inf).  The
     visual effect is "denser cloud where phi is higher" -- wetting layers
     around each sphere show up as bright shells in the otherwise dilute box.

Usage:
    python viz_util/visualize_two_spheres_3d.py path/to/output_cpuK.json \\
        --alpha 1 --beta 1 --d_index 0 [--out figure.png]

NaN values in the JSON (JSON.jl encodes NaN as null -> Python None) are
converted in-place to np.nan when reading scalar arrays.
"""

import argparse
import json
import os

import matplotlib.pyplot as plt
import numpy as np


def _to_float(x):
    """Convert None -> NaN, otherwise float."""
    return float("nan") if x is None else float(x)


def _to_array(seq):
    return np.array([_to_float(v) for v in seq], dtype=float)


def load_result(path, alpha_key, beta_key):
    with open(path) as f:
        data = json.load(f)
    return data[str(alpha_key)][str(beta_key)]


def bilinear_interp(phi, rho_grid, z_grid, rho_q, z_q, fill):
    """
    Vectorised bilinear interpolation of phi(rho, z) at (rho_q, z_q).  Points
    outside the grid get `fill`.  phi shape is (Nrho, Nz).
    """
    Nrho, Nz = phi.shape
    drho = rho_grid[1] - rho_grid[0]
    dz   = z_grid[1]   - z_grid[0]
    out_of = ((rho_q < rho_grid[0]) | (rho_q > rho_grid[-1]) |
              (z_q   < z_grid[0])   | (z_q   > z_grid[-1]))
    ti = (rho_q - rho_grid[0]) / drho
    tj = (z_q   - z_grid[0])   / dz
    i0 = np.clip(np.floor(ti).astype(int), 0, Nrho - 2)
    j0 = np.clip(np.floor(tj).astype(int), 0, Nz - 2)
    t_rho = ti - i0
    t_z   = tj - j0
    p00 = phi[i0,     j0]
    p10 = phi[i0 + 1, j0]
    p01 = phi[i0,     j0 + 1]
    p11 = phi[i0 + 1, j0 + 1]
    out = ((1 - t_rho) * (1 - t_z) * p00 +
           t_rho       * (1 - t_z) * p10 +
           (1 - t_rho) * t_z       * p01 +
           t_rho       * t_z       * p11)
    out[out_of] = fill
    return out


def sample_density_cloud(phi, rho_grid, z_grid, R, z_c1, z_c2, phi_inf,
                         n_samples, seed):
    """
    Rejection sampling: accept random Cartesian points (x, y, z) with
    probability proportional to max(0, phi(rho, z) - phi_inf).  Returns
    (xyz, phi_at_xyz) arrays of length n_samples.  Points inside either
    sphere are rejected.
    """
    rng = np.random.default_rng(seed)
    rho_max = rho_grid[-1]
    z_min, z_max = z_grid[0], z_grid[-1]
    weight_max = max(phi.max() - phi_inf, 1e-12)

    xs = np.empty(0); ys = np.empty(0); zs = np.empty(0); phs = np.empty(0)
    while xs.size < n_samples:
        n_try = max(2000, 4 * (n_samples - xs.size))
        # Uniform in a cylinder (rho, theta, z): rho^2 = uniform * rho_max^2
        rr = rho_max * np.sqrt(rng.uniform(0.0, 1.0, n_try))
        th = rng.uniform(0.0, 2.0 * np.pi, n_try)
        zz = rng.uniform(z_min, z_max, n_try)
        phi_at = bilinear_interp(phi, rho_grid, z_grid, rr, zz, phi_inf)
        d1sq = rr ** 2 + (zz - z_c1) ** 2
        d2sq = rr ** 2 + (zz - z_c2) ** 2
        outside = (d1sq > R ** 2) & (d2sq > R ** 2)
        w = np.clip(phi_at - phi_inf, 0.0, None)
        accept = (rng.uniform(0.0, weight_max, n_try) < w) & outside
        xs  = np.concatenate([xs,  rr[accept] * np.cos(th[accept])])
        ys  = np.concatenate([ys,  rr[accept] * np.sin(th[accept])])
        zs  = np.concatenate([zs,  zz[accept]])
        phs = np.concatenate([phs, phi_at[accept]])
    return (np.column_stack([xs[:n_samples], ys[:n_samples], zs[:n_samples]]),
            phs[:n_samples])


def render(result, d_index, out_path, n_samples, seed):
    R       = float(result["R"])
    phi_inf = float(result["phi_inf"])
    omega_b = float(result["omega_b"])

    profile_d   = _to_array(result["profile_d"])
    if d_index < 0 or d_index >= len(profile_d):
        raise IndexError(f"d_index {d_index} out of range "
                         f"(have {len(profile_d)} converged d-points)")
    d   = float(profile_d[d_index])
    rho = _to_array(result["profile_rho"][d_index])
    z   = _to_array(result["profile_z"][d_index])
    phi = np.array(result["profile_phi"][d_index], dtype=float)  # (Nrho, Nz)

    z_c1 = 0.0
    z_c2 = 2.0 * R + d

    pts, phi_pts = sample_density_cloud(
        phi, rho, z, R, z_c1, z_c2, phi_inf,
        n_samples=n_samples, seed=seed,
    )

    fig = plt.figure(figsize=(7.5, 11.0))
    ax  = fig.add_subplot(111, projection="3d")

    # Two spheres -- semi-opaque, light grey
    u  = np.linspace(0.0, 2.0 * np.pi, 48)
    v  = np.linspace(0.0, np.pi, 36)
    xs = R * np.outer(np.cos(u), np.sin(v))
    ys = R * np.outer(np.sin(u), np.sin(v))
    zs_c = np.outer(np.ones_like(u), np.cos(v))
    for zc in (z_c1, z_c2):
        ax.plot_surface(xs, ys, zc + R * zs_c,
                        color="#d0d0d0", alpha=0.55, linewidth=0,
                        edgecolor="none", shade=True, antialiased=True)

    # Density cloud -- denser where phi is higher; color also encodes phi
    sc = ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2],
                    c=phi_pts, cmap="magma",
                    s=2, alpha=0.18, edgecolor="none",
                    depthshade=False,
                    vmin=phi_inf, vmax=max(phi.max(), phi_inf + 1e-9))
    cb = plt.colorbar(sc, ax=ax, shrink=0.55, pad=0.02)
    cb.set_label(r"$\phi$")

    # Equal aspect in (x, y); z is naturally longer than rho_max in this setup
    rho_max = rho[-1]
    z_span  = z[-1] - z[0]
    ax.set_box_aspect((1.0, 1.0, z_span / (2.0 * rho_max)))
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title(
        f"two-sphere wet field: phi_inf={phi_inf}, omega_b={omega_b}, "
        f"R={R}, d={d:.3f}"
        f"\nphi range = [{phi.min():.4f}, {phi.max():.4f}]"
    )

    # Trim panes
    ax.xaxis.pane.set_alpha(0.0)
    ax.yaxis.pane.set_alpha(0.0)
    ax.zaxis.pane.set_alpha(0.0)

    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"saved -> {out_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("json_path",
                        help="Path to output_cpuK.json from spec/fisher_two_spheres.jl")
    parser.add_argument("--alpha", type=int, default=1,
                        help="Top-level alpha key (1-indexed)")
    parser.add_argument("--beta", type=int, default=1,
                        help="Inner beta key (1-indexed)")
    parser.add_argument("--d_index", type=int, default=0,
                        help="Index into res['profile_d'] (0-indexed; "
                             "only converged d-points are stored)")
    parser.add_argument("--n_samples", type=int, default=80000,
                        help="Number of accepted density-cloud points")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", default=None,
                        help="Output PNG path (default: alongside JSON)")
    args = parser.parse_args()

    result = load_result(args.json_path, args.alpha, args.beta)
    if args.out is None:
        base = os.path.splitext(os.path.basename(args.json_path))[0]
        outdir = os.path.dirname(os.path.abspath(args.json_path))
        args.out = os.path.join(outdir,
            f"{base}_3d_a{args.alpha}_b{args.beta}_d{args.d_index}.png")
    render(result, args.d_index, args.out,
           n_samples=args.n_samples, seed=args.seed)


if __name__ == "__main__":
    main()
