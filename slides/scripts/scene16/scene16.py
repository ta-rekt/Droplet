"""
Scene 16 -- "This plot shows the maximum possible Fisher rate for each given
concentration. You can clearly see the two regimes: diffusion-limited, and
activation-time limited."

Figure 2a of the SPRT paper, re-drawn from its source (ConstrainedPromoters,
figures.ipynb, "Figure 2 / A / Plot") in the slides' font: the maximum Fisher
rate over four-state circuits ("full" decoding), the diffusion limit 4aD c
and the activation limit 0.16 kmax. Left out: the "acc" curve, its diffusion
limit (8aD c/3) and the vertical line at the crossover.

On top, the Fisher rate (scene 15's expression) of five circuits optimal at
concentrations c0 from the crossover up, with the paper's parameters and
optimal rates: kb = min(4aD, 0.503 kmax / c0), kon = ku = koff = kmax. Each
peaks at its own c0 on the maximum, at 0.808 = 0.1616 kmax (checked when the
script runs). Below the crossover, the optimal circuit is a single one,
kb = 4aD: the maximum there is the rising side of the first curve.

Data: fig2a_full.json, a copy of ConstrainedPromoters/data/fig2/a1/
ir_full_vs_ep_pareto_pid1.json (commit dea9418); one entry per c0: the eight
log10 rates of the optimal circuit, its drift V at dc/c = 0.1, and c0.

Output: slides/figures/scene16/scene16.{pdf,svg}, a full 128 x 96 mm slide.

Usage (from slides/):
    python3 scripts/scene16/scene16.py
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

import numpy as np
from scipy.optimize import minimize_scalar

from common.drawing import LABEL_PT, SLIDE_UNITS, TEXT_COLOR, new_slide, save
from scene15.scene15 import fisher_off

# ---------------------------------------------------------------------------
# The paper's parameters and data
# ---------------------------------------------------------------------------

FOUR_AD = 0.0888       # 4 a D: a = 3e-3 um, D = 7.4 um^2/s
K_MAX = 5.0            # maximum activation and unbinding rates (1/s)
X_STAR = 0.503427      # optimal c kb, in units of kmax (SI)
ACTIVATION_LIMIT = 1 / (1.0 + 2.5 ** 2 / 25)  # 0.8 = 0.16 kmax, as plotted in the paper
DRIFT_TO_FISHER = 2 / 0.1 ** 2  # V = 1/2 (dc/c)^2 Fdot, at dc/c = 0.1

DATA = Path(__file__).with_name("fig2a_full.json")


def paper_curve():
    """(c0, maximum Fisher rate) as plotted in the paper."""
    rows = json.loads(DATA.read_text())
    rows = [rows[str(i)] for i in range(1, len(rows) + 1)]
    return (np.array([r[9] for r in rows]),
            DRIFT_TO_FISHER * np.array([r[8] for r in rows]))


def fisher_rate(c, kb):
    """Fisher rate at c of the circuit binding at kb, other rates at kmax and
    counter-clockwise rates zero."""
    x = c * kb
    return fisher_off(x, K_MAX) / (1 / x + 3 / K_MAX)


def optimal_kb(c0):
    return min(FOUR_AD, X_STAR * K_MAX / c0)


# Five circuits, optimal from the crossover c0 = 0.503 kmax / 4aD up.
CROSSOVER = X_STAR * K_MAX / FOUR_AD
C0S = CROSSOVER * 2.5 ** np.arange(5)


def check_peaks(c, maximum, tolerance=0.002):
    """Each curve must peak at its c0, on the paper's maximum."""
    for c0 in C0S:
        kb = optimal_kb(c0)
        peak = minimize_scalar(lambda lc: -fisher_rate(np.exp(lc), kb),
                               bounds=(np.log(c0) - 3, np.log(c0) + 3), method="bounded")
        c_peak, value = np.exp(peak.x), -peak.fun
        on_curve = np.exp(np.interp(np.log(c_peak), np.log(c), np.log(maximum)))
        print(f"c0 = {c0:7.1f}: peak at c = {c_peak:7.1f}, Fdot = {value:.4f}, "
              f"paper's maximum there {on_curve:.4f} ({value / on_curve - 1:+.2%})")
        assert abs(c_peak / c0 - 1) < tolerance and abs(value / on_curve - 1) < tolerance


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

AXES = (3.2, 2.9, 15.6, 9.6)  # left, bottom, width, height (slide units)
X_LIMITS, Y_LIMITS = (8e-2, 1.2e3), (6e-3, 16)  # as in the paper
TICK_PT = 10
Y_LABEL = r"$\max_{\rm rates} \dot{\mathcal F}(c)$"  # the best circuit at each c
CURVE_COLOR = "#2F62C8"
LIMIT_STYLE = dict(color="black", alpha=0.5, lw=1.2)
DIFFUSION_LABEL_C = 0.25   # where the diffusion limit is labelled...
LABEL_LIFT = 2.0           # ...this factor above the line
ACTIVATION_LABEL_C = 0.1


def paper_axes(fig, box, x_limits, y_limits, y_label, tick_pt=TICK_PT):
    """Log-log axes at `box` (left, bottom, width, height in slide units), in
    the style of the paper's figures but the slides' font, with the y label
    horizontal on top of the axis."""
    left, bottom, width, height = box
    ax = fig.add_axes([left / SLIDE_UNITS[0], bottom / SLIDE_UNITS[1],
                       width / SLIDE_UNITS[0], height / SLIDE_UNITS[1]])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    # Ticks at the powers of ten in range.
    for limits, set_ticks in ((x_limits, ax.set_xticks), (y_limits, ax.set_yticks)):
        set_ticks(10.0 ** np.arange(np.ceil(np.log10(limits[0])),
                                    np.floor(np.log10(limits[1])) + 1))
    ax.tick_params(labelsize=tick_pt, colors=TEXT_COLOR)
    ax.set_xlabel("Concentration", fontsize=LABEL_PT, color=TEXT_COLOR)
    ax.text(0, 1.03, y_label, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=LABEL_PT, color=TEXT_COLOR)
    return ax


def label_along(ax, c, y, text, slope, size=LABEL_PT, below=False):
    """Text starting at (c, y), above (or `below`) a line of the given log-log
    slope, rotated along it."""
    (x0, y0), (x1, y1) = ax.transData.transform([(c, y), (10 * c, y * 10 ** slope)])
    angle = np.degrees(np.arctan2(y1 - y0, x1 - x0))
    ax.text(c, y, text, rotation=angle, rotation_mode="anchor", ha="left",
            va="top" if below else "bottom", fontsize=size, color=TEXT_COLOR)


if __name__ == "__main__":
    c, maximum = paper_curve()
    check_peaks(c, maximum)

    fig, _ = new_slide()
    ax = paper_axes(fig, AXES, X_LIMITS, Y_LIMITS, Y_LABEL)
    ax.plot(c, maximum, color="black", lw=1.8, zorder=4)
    ax.plot(c, FOUR_AD * c, **LIMIT_STYLE)
    ax.axhline(ACTIVATION_LIMIT, linestyle="dotted", **LIMIT_STYLE)
    cs = np.geomspace(*X_LIMITS, 400)
    for c0 in C0S:
        ax.plot(cs, [fisher_rate(ci, optimal_kb(c0)) for ci in cs], color=CURVE_COLOR, lw=1.3,
                zorder=3)
    label_along(ax, DIFFUSION_LABEL_C, LABEL_LIFT * FOUR_AD * DIFFUSION_LABEL_C,
                "Diffusion limit", slope=1)
    label_along(ax, ACTIVATION_LABEL_C, 1.15 * ACTIVATION_LIMIT, "Activation limit", slope=0)
    save(fig, "scene16", "scene16")
