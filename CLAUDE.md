# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Droplet is a Julia package for lattice Monte Carlo simulations of polymers and biomolecular condensates (transcription factors, RNAP, DNA polymers). The Julia core (`src/`) defines species types and Metropolis dynamics; experiments are written as small spec files that get executed in parallel by a sweep driver and post-processed by Python dashboards loaded from `viz.ipynb`.

## Running Simulations

The standard workflow is **sweep driver + spec file + output dir**, parallelized across CPUs by `runJL`:

```bash
./runJL sim/sweep.jl spec/<experiment>.jl <ncpus> data/YYYYMMDD/<experiment>_<info>
```

`runJL` builds a pipe of N `julia --project=@. sim/sweep.jl <specfile> <cpukey> <ncpus> <outdir>` calls, one per CPU. `sim/sweep.jl` `include`s the spec file (path is resolved as `../<specfile>` relative to where Julia is invoked from), then iterates `target!(res, slist, keya, alpha, beta, blurb)` over its slice of `alphas` × `betas`, writing `output_cpu<k>.json` per CPU.

Spec files therefore must define at module-top-level: `slist`, `alphas`, `betas`, `nparams`, `blurb`, `checkpoint_path`, and a mutating `target!(res, args...)`. See `context/spec_template_guide.md` for the full contract — the `args` unpacking is `slist, keya, alpha, beta, blurb` (always include `beta` even if unused).

## Tests

```bash
julia --project=@. -e 'using Pkg; Pkg.test()'
# or run a single file:
julia --project=@. test/test_hamiltonians.jl
```

`test/runtests.jl` is the orchestrator. **Do not run the full suite during iteration** — it is large and expensive. Iterate on individual `test/test_*.jl` files and only run the full suite at the end of an implementation phase to check global consistency (per `context/style_guide.md`).

## Architecture

`src/Droplet.jl` is a thin module that includes, in order: `types.jl`, `interactions.jl`, `hamiltonian.jl`, `metropolis.jl`, `move.jl`, `port.jl`, `init.jl`, `utils.jl`, `stats.jl`, `ctmc.jl`, `generate.jl`, `get_residence_time.jl`. The order matters — `types.jl` defines the structs (`Lattice`, `FastBulk`/`FastBulkMonomer`, `Bulk` (legacy), `Polymer`) that everything else dispatches on.

Key conventions baked into the API:
- Species carry their own state (`s.pos`, `s.ham`, `s.deltaH`, `s.intSpecies`, `s.intIdx`, `s.intMat`). `move!(s)` returns `Bool` (accepted), mutates `s.pos` and `s.deltaH`. `ham!(s)` recomputes `s.ham`. `intPot!(s, v)` and `intPot!(s1, s2, v)` register self/pair interactions.
- The Monte Carlo loop in `target!` is wall-clock-bounded (`while t_current <= t_max`, with `t_current = (time_ns() - t_start)/6e10` minutes), samples every `step` accepted moves, and only accumulates spatial density after `t_burn_in`.
- Polymers need `initMonomers!(poly, monomer_list)` before `init!(l, slist)`. `FastBulk` is the default; `Bulk` is legacy.
- Output dicts must include `species_name`, `nframes`, `step`, `boxlength`, `runtime_min`; common extras are `binding`, `transcription`, `*_density_lattice`, `hamiltonian`, `acceptance_rate`, `*_final_coordinates`.

Three-way (coupling) interactions and time-varying Gillespie dynamics are separate subsystems — see `test_couplingPot.jl` / `test_*_coupling.jl` and `test_generate_tv.jl` respectively.

## Authoring Experiments

Per `context/CRITICAL_PRINCIPLES.md` and `context/style_guide.md`:

- **Copy, don't write from scratch.** Start from the closest existing `spec/*.jl` (catalog in `context/experiment_types.md`) and modify minimally. Same for Python dashboards (`viz_util/dashboard_*.py`).
- **Start small.** First runs are always `t_max = 1` minute, `ncpus = 1`, regardless of what the source spec used. Scale up only after the proposal/spec/dashboard validate end-to-end.
- **Naming/locations:** proposals → `prompts/YYYYMMDD/proposal_<name>.md`; specs → `spec/<name>.jl`; dashboards → `viz_util/dashboard_<name>.py`; outputs → `data/YYYYMMDD/<name>_<boxlength>_<info>/`.
- **Frugal code.** Match existing style; avoid helper functions, abstractions, or refactors beyond what the task requires. Suggest higher-level reuse instead of duplicating logic.

The `data/` and `jobs/` directories are gitignored.

## Paper Figures

Figures for the manuscript are built in `figure_production/` — one notebook per
figure, reading `figure_production/processed_data/` and writing
`figure_production/figures/<fig>/`.

**Before producing or changing any figure, read `figure_production/FIGURE_STYLE.md`
and follow it.** It is the source of truth for figure typography and geometry, and
it exists because these mistakes are easy to make and invisible until the figure is
on the page:

- all text at the caption size (9 pt *as printed*), tick labels 8 pt, panel letters Arial 12 pt;
- build at the width the figure is included at, so LaTeX's scale factor is 1 —
  a column-width PDF dropped into a `figure*` prints at nearly twice the intended size;
- never `bbox_inches="tight"` on a print-sized figure; it re-crops the page and
  LaTeX silently rescales it;
- verify the output PDF's page size and embedded font sizes rather than trusting
  the rcParams.

## Context Files

When the task is generating proposals, specs, or visualizations, read the corresponding doc in `context/` first — it is the source of truth and more current than this file:
- `context/CRITICAL_PRINCIPLES.md` — proposal format, scaling rules, common mistakes
- `context/spec_template_guide.md` — full spec-file syntax and `target!` contract
- `context/visualization_guide.md` — dashboard structure and the spec→dashboard→`viz.ipynb` pipeline
- `context/experiment_types.md` — catalog of existing experiments to copy from
- `figure_production/FIGURE_STYLE.md` — **required reading before producing or
  modifying any paper figure** (type sizes, scale factors, panel letters)
- `context/domain_knowledge.md` — physics/biology background
- `context/heterogeneous_polymer_implementation.md` — heteropolymer subsystem details
- VERY IMPORTANT: please don't use skills without asking first. do not load a skill
autonomously, even if it seems relevant in the context of the prompt