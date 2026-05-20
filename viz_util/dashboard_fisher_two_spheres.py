"""
Load output_cpu*.json files from a fisher_two_spheres sweep, assemble
phi_probe(phi_inf) per (omega_b, R) group, compute the Fisher information
rate, and emit both a JSON dump and a 3-panel figure:

  - top:    U_full(d) for a few representative phi_inf
  - middle: phi_probe vs phi_inf (with phi_inf as control reference)
  - bottom: Fisher rate vs phi_inf (log y), control + with-interaction

Output:
  media/figures/<date>/fisher_two_spheres.png
  media/figures/<date>/fisher_two_spheres.json
"""

import sys, os, glob, json
from datetime import date

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ── Config ───────────────────────────────────────────────────────────────────
EXP_NAME = "fisher_two_spheres_test"
OUT_DATE = date.today().strftime("%Y%m%d")
DATA_DIR = f"data/{OUT_DATE}/{EXP_NAME}"
OUT_PNG  = f"media/figures/{OUT_DATE}/{EXP_NAME}.pdf"
OUT_JSON = f"media/figures/{OUT_DATE}/{EXP_NAME}.json"

# Probe rate constants for the Fisher rate (match fisher_rate_from_profile)
KB = 1.0
KU = 1.0


# ── Load ─────────────────────────────────────────────────────────────────────
files = sorted(glob.glob(os.path.join(DATA_DIR, "output_cpu*.json")))
if not files:
    sys.exit(f"no output_cpu*.json under {DATA_DIR}")

records = []
for fp in files:
    with open(fp) as f:
        blob = json.load(f)
    for _, by_b in blob.items():
        for _, by_c in by_b.items():
            for _, rec in by_c.items():
                records.append(rec)

print(f"loaded {len(records)} records from {len(files)} files")


# ── Fisher rate ──────────────────────────────────────────────────────────────
def fisher_rate(phi_probe, phi_tf, kb=KB, ku=KU):
    inv_tau      = phi_probe * kb * ku / (kb + ku)
    fisher_base  = phi_tf * kb * ku / (kb + ku)
    with np.errstate(divide="ignore", invalid="ignore"):
        dphi_probe   = np.gradient(phi_probe, phi_tf)
        fisher_probe = phi_tf ** 2 * inv_tau * (dphi_probe / phi_probe) ** 2
    return fisher_base, fisher_probe


# ── Group by (omega_b, R), assemble phi_inf series ───────────────────────────
groups = {}
for r in records:
    key = (round(r["omega_b"], 12), round(r["R"], 12))
    groups.setdefault(key, []).append(r)

assembled = {}
for (omega_b, R), recs in groups.items():
    recs.sort(key=lambda x: x["phi_inf"])
    phi_inf  = np.array([r["phi_inf"]   for r in recs], dtype=float)
    phi_probe = np.array([r["phi_probe"] for r in recs], dtype=float)
    d_grid_list = [np.asarray(r["d_grid"], dtype=float) for r in recs]
    d_grid = d_grid_list[0]
    if not all(arr.shape == d_grid.shape and np.allclose(arr, d_grid)
               for arr in d_grid_list):
        print(f"  warn: d_grid varies within (omega_b={omega_b}, R={R}); "
              "using first record's grid")
    U_full_grid   = np.array([r["U_full"]   for r in recs], dtype=float)
    Phi_surf_grid = np.array([r["Phi_surf"] for r in recs], dtype=float)

    valid = np.isfinite(phi_probe) & (phi_probe > 0)
    fisher_base  = np.full_like(phi_inf, np.nan)
    fisher_probe = np.full_like(phi_inf, np.nan)
    if int(valid.sum()) >= 2:
        fb_v, fp_v = fisher_rate(phi_probe[valid], phi_inf[valid])
        fisher_base[valid]  = fb_v
        fisher_probe[valid] = fp_v

    assembled[(omega_b, R)] = dict(
        phi_inf=phi_inf,
        phi_probe=phi_probe,
        d_grid=d_grid,
        U_full_grid=U_full_grid,
        Phi_surf_grid=Phi_surf_grid,
        fisher_base=fisher_base,
        fisher_probe=fisher_probe,
    )
    print(f"(omega_b={omega_b:g}, R={R:g}): {len(recs)} phi_inf points, "
          f"{int(valid.sum())} valid for Fisher rate")


# ── Dump aggregated JSON ─────────────────────────────────────────────────────
os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)

def _jsonify(a):
    if isinstance(a, np.ndarray):
        return a.tolist()
    return a

any_rec = records[0]
out_blob = dict(
    blurb=any_rec.get("blurb", ""),
    chi_b=any_rec.get("chi_b"),
    s=any_rec.get("s"),
    a=any_rec.get("a"),
    kb=KB,
    ku=KU,
    groups=[
        dict(
            omega_b=omega_b,
            R=R,
            phi_inf=_jsonify(g["phi_inf"]),
            phi_probe=_jsonify(g["phi_probe"]),
            d_grid=_jsonify(g["d_grid"]),
            U_full_grid=_jsonify(g["U_full_grid"]),
            Phi_surf_grid=_jsonify(g["Phi_surf_grid"]),
            fisher_base=_jsonify(g["fisher_base"]),
            fisher_probe=_jsonify(g["fisher_probe"]),
        )
        for (omega_b, R), g in assembled.items()
    ],
)
with open(OUT_JSON, "w") as f:
    json.dump(out_blob, f, indent=2)
print(f"wrote -> {OUT_JSON}")


# ── Plot ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(3, 1, figsize=(5.4, 9.0))

# We plot the FIRST (omega_b, R) group on the figure.  If the sweep includes
# multiple operating points, additional groups go to additional colour cycles
# but share the axes.
(omega_b, R) = next(iter(assembled))
g = assembled[(omega_b, R)]
phi_inf  = g["phi_inf"]
phi_probe = g["phi_probe"]
d_grid   = g["d_grid"]
U_full_grid   = g["U_full_grid"]
fisher_base   = g["fisher_base"]
fisher_probe  = g["fisher_probe"]
chi_b = any_rec.get("chi_b")

# Top: U_full(d) for a few phi_inf values
ax = axes[0]
n_show = min(5, len(phi_inf))
sel = np.linspace(0, len(phi_inf) - 1, n_show).astype(int)
cmap = plt.get_cmap("viridis")
for k, idx in enumerate(sel):
    if not np.all(np.isnan(U_full_grid[idx])):
        ax.plot(d_grid, U_full_grid[idx], "o-",
                color=cmap(k / max(1, n_show - 1)),
                lw=1.5, ms=3,
                label=fr"$\phi_\infty$={phi_inf[idx]:.3f}")
ax.axhline(0.0, color="gray", lw=0.6, ls="--")
ax.set_xlabel(r"$d$ (surface-to-surface)")
ax.set_ylabel(r"$U_{\rm full}(d)$")
ax.set_title(fr"Full 2D pair potential ($\chi_b$={chi_b}, "
             fr"$\omega_b$={omega_b:g}, $R$={R:g})")
ax.legend(fontsize=8, loc="best")
ax.grid(True, alpha=0.3)

# Middle: phi_probe vs phi_inf
ax = axes[1]
ax.plot(phi_inf, phi_probe, "o-", color="C2", lw=2, label=r"$\phi_{\rm probe}$")
ax.plot(phi_inf, phi_inf, color="gray", lw=1.0, ls="--",
        label=r"$\phi_\infty$ (no interaction)")
ax.set_xlabel(r"$\phi_\infty$ (TF total volume fraction)")
ax.set_ylabel(r"$\phi_{\rm probe}$")
ax.legend(fontsize=9)
ax.set_yscale('log')
ax.grid(True, alpha=0.3)

# Bottom: Fisher rate vs phi_inf
ax = axes[2]
ax.plot(phi_inf, fisher_base,  color="red", lw=1.0, ls="--", label="control")
ax.plot(phi_inf, fisher_probe, color="red", lw=2.0,           label="with interaction")
ax.set_xlabel(r"$\phi_\infty$")
ax.set_ylabel("Fisher information rate")
ax.set_yscale("log")
ax.legend(fontsize=9)
ax.grid(True, which="both", alpha=0.3)

fig.tight_layout()
os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
print(f"saved -> {OUT_PNG}")
