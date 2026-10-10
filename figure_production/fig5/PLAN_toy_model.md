# Fig. 5 — switch back to the June toy model (plan)

Status: **done 2026-10-09** (Claude Code; see `../README.md` Fig. 5 and `droplet/CLAUDE.md`). Written 2026-10-09 by the Cowork session, to be picked up from a
Claude Code session in a terminal. Read `../../../CLAUDE.md` (the `droplet/`
notes, one level above the repo) first: house style, traps, standing
instructions. Do **not** edit the repo's own `Droplet/CLAUDE.md` (shared with
Tarek) and never write `fig5/fig5.pdf` yourself — the user runs the notebook.

## Why

The current `fig5.ipynb` rebuilds the old Fig. 4 of Tarek's draft
(`Tarek_Droplets/figs/main-fig4.pdf`) but feeds it the **s = 110 physical
sweep** (BVP surface potential, lattice units). In that model the Rouse
confinement is almost flat over the plotted range, so panels C, D and E come out
flat and look nothing like the original. The user chose to reproduce the
original model instead.

The original panels come from a **toy, dimensionless model**:

- Panel D is `viz-rouse.ipynb` cell 37 (same cell in the 28 July commit
  `ef5a572` "figs 1-4 v1 done", which also committed its output
  `media/figures/20260727/rouse_distance_distribution_lowc_hic.pdf`).
  It reads `data/20260615/rouse_mfpt_attractive` (cell 21), sorted by alpha,
  and plots indices **0 and 9** (low, high c), `bins=100, range=(R, L)`,
  `density=True`, against `exp(-U_tot/kBT)` normalised by `trapz` on
  `linspace(R, L, 500)`.
- That data folder no longer exists, but the spec that produced it is in git:
  `git show 9f78af2:spec/rouse_mfpt_attractive.jl` (Tarek, 2026-06-15,
  "about to launch a big sim to check rouse exit rate").

### The June model (from that spec; Python twin in `viz-rouse.ipynb` cell 23)

| | value |
|---|---|
| surface | nb = 1, chi_b = 2.5, omega_b = h0 = -0.5, kappa = 1, eps = 1e-6, **R = 2.5**, L = 10 (Parameters default) |
| concentrations | `alphas = range(0.003, stop=0.15, length=16)` = phi_inf; low = index 1 (0.003), high = index 10 (0.0912) (1-based) |
| lambda | `sqrt(d2fb_dphi2(phi_inf, p) / kappa)` |
| U_ads(d) | `-(4π R³ h0² exp(-λ(d-R))) / (κ (λ + 1/R) d)` |
| U_tot(d) | `U_ads(d) - 2 log d + 3 d² / (2 s b²)`, s = beadSep = N-1 = 99, b = 1 |
| Rouse | N = 100 beads, kBT = 1, friction = 1, b = 1, tracked beads 1 and 100, `r_min = R`, `dt = 0.001*friction/(2kBT/b²) = 5e-4`, nsteps = 5e7, `saveMode = :distance_only`, `sample_skip = 500`, `seed = keya + 42` |

Checked: the Rouse integrator's physics is unchanged since June
(`git diff 9f78af2 HEAD -- src/rouse.jl` only adds sojourn tracking, an O(N)
tridiagonal drift and `:sojourn_only`; stiffness is 3kBT/b² both then and now).
The current `RousePolymer` and Julia `Parameters` accept the June keywords
(`a`, `b` sojourn bounds are optional). `sim/sweep.jl` still only needs
`slist, checkpoint_path, alphas, betas, nparams, blurb, target!`.

Caption of the old figure (Tarek_Droplets `rouse.tex`, label `rouse-setup`) for
the landmarks: ridge d* per concentration (dashed lines in C),
**d^o = 1.2 max_{c in [c_min, c_max]} d*(c)**, **d^x = 2R**.

## Steps

### 1. Rouse spec — `fig5/rouse_toy.jl` (new)

Start from `git show 9f78af2:spec/rouse_mfpt_attractive.jl`. Change only:

- `alphas = collect(range(0.003, stop=0.15, length=16))[[1, 10]]` (2 runs).
  Optionally 4 seeds per concentration for smoother histograms
  (`[[1,1,1,1,10,10,10,10]]`); the original used one run each.
- seed: reproduce the June seeds — June keya 1 and 10 → seeds 43 and 52. With
  the 2-run list, use e.g. `seed = (alpha == alphas_all[1] ? 43 : 52) + 1000*(rep-1)`;
  keep it simple and document it. (Bitwise reproduction is not expected —
  different Julia/RNG versions — only statistical agreement.)
- Also store `res["ads_potential"]`, `res["pot_grid"] = dists`, and
  `res["x_ridge"]` is not needed (computed in Python).
- Header comment: provenance (commit 9f78af2, the data folder it reproduces).

Delete the s110 files this replaces once the toy pipeline works:
`fig5/param.jl`, `fig5/rouse_mfpt_attractive.jl` (both untracked). Deletion in
the user's folders needs his OK.

### 2. Run script — rewrite `fig5/run_fig5_rouse.sh`

Keep its structure (repo-root check, julia check, `Pkg.instantiate()`, ncpus
choice, final verification) but:

- drop the s110 input-copy step (the toy model needs no input files);
- `./runJL sim/sweep.jl figure_production/fig5/rouse_toy.jl $NCPUS data/<today>/fig5_rouse_toy`
  (ncpus must divide the number of runs: 1 or 2, or 1/2/4/8 with 8 runs);
- memory ≈ nsteps × 8 B = 0.4 GB per process (`:distance_only` allocates
  `zeros(nsteps)`); runtime unknown for N = 100 at 5e7 steps — time one run
  first (the s110 runs did 1e8 steps in ~12.6 min).

Julia must be **1.11.1** (Manifest). Run from the repo root.

### 3. `fig5/prepare_fig5_data.py` — rewrite for the toy model

Everything except the histograms is analytic, so compute it in Python with
`viz_util.bvp_solver.Parameters` and `d2fb_dphi2` (import with the repo root on
`sys.path`; the script already resolves the repo). Write to
`figure_production/processed_data/` (`%.17g`):

- `fig5_potentials.csv`: `d` = `linspace(R, L, 500)`, `U_ads_k{i}` and
  `U_tot_k{i}` for all 16 concentrations (small, and E/A need several).
- `fig5_scalars.json`: model parameters, the 16 phi_inf, indices used
  (low 0 / high 9, 0-based; panel A's three, e.g. 0, 7, 15 — choose so the
  curves are distinct), ridge d*(c) for each c (first sign change of
  F = -dU/dd past the contact well — same rule as cell 23), d^o = 1.2 max d*,
  d^x = R in code coordinates (labelled 2R in the figure, see below).
- `fig5_hist.csv` only when `data/<date>/fig5_rouse_toy/output_cpu*.json`
  exist: 100 bins on (R, L), density as in cell 37 (normalised over in-range
  samples), plus counts; burn-in: cell 37 used none — keep none, or 5 % and
  say so. Check each run's `total_potential` against the Python U_tot
  (same formula, so agreement to ~1e-10, not the 1e-3 tolerance of the s110
  version) and its `alpha` against the two phi_inf.

### 4. `fig5/fig5.ipynb` — adapt (edit the notebook directly; never generate it)

Keep: all artwork (`art_*.svg`, `GROUPS`, stamping), the geometry constants,
house style, `schematic_xlabel`, the overlap-free text placements. Change:

- Data cell: read the new CSVs/JSON; drop `NM_PER_SITE` and nm units
  everywhere. D's x-axis becomes dimensionless `Distance $d$`, R to L
  (2.5–10, ticks 4, 6, 8, 10 as in the original).
- A (plot): three U_ads curves, increasing c; arrow "Increasing $c$".
- C: U_tot low/high on [R, L] with dashed ridge lines at d*; y-limits from
  the data. **The inset** (Rouse + surface attraction sketch) currently sits in
  the lower right of panel C-left, which is free only because the s110 curves
  are flat. With the toy curves (barrier near contact, broad minimum mid-range,
  rise at large d) it will collide: move it to where the original had it —
  above the curves, top centre of C-left — and re-check `GROUPS["C"]`,
  `GROUP_SCALE["C"]` and the "Rouse"/"Surface attraction" label positions.
- D: as cell 37 (log y, histogram + theory, ridge dashed, "Closed"/"Open"
  above the high-c frame as now, legend in high-c panel).
- E: U_tot at high c; d^x, d*, d^o dashed lines and labels; re-tune the
  exit / re-entry arc endpoints to the new curve (original: exit from the
  contact well over the barrier to d^o; re-entry from beyond d^o back to d^x,
  arching above).

Open question to raise with the user, not to decide silently: code
coordinates have contact at d = R (D's axis starts at 2.5 = R), while the
schematic panels and the caption say d^x = 2R. The caption itself carries
"{\TM redefine d}". Keep the original's labelling unless told otherwise.

### 5. Checks

- Compare D (and C) with the original: `git show ef5a572:media/figures/20260727/rouse_distance_distribution_lowc_hic.pdf > /tmp/orig_D.pdf`
  and `Tarek_Droplets/figs/main-fig4.pdf` (low-c PDF: spike at contact, dip near
  2.8, broad maximum near 8–9; high-c: dominant contact peak, plateau ~0.05).
- Fonts/strokes with the pdfplumber snippet in `../README.md` (9 / 8 / 12 pt,
  ≤ 0.5 pt strokes). Text overlaps: render, then check every visible text
  artist's window extent pairwise.
- Render with real LaTeX to a scratch path only, never into the repo.

### 6. Docs

- `../README.md`, Fig. 5 section: replace the s110 description (inputs, k = 0/9,
  nm axis, runtimes, the data folder) with the toy model and the new commands.
- `droplet/CLAUDE.md` (outside the repo), "Fig. 5 specifics": same update;
  remove the s110 notes.
- Delete this plan file when done (or mark it done).
