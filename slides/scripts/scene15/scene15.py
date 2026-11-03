"""
Scene 15 -- "Up till now, we have used the simplest possible model of
transcription initiation [...] We can capture this in the following model
with 4 states instead of 2. You can notice the delay in activating
transcription, governed by this rate kon. [...] the Fisher rate: it decreases
at high c. [...]"

The cartoon of scene 11 (same molecules, same places); the four-state model
(TF unbound / bound x transcription OFF / ON), going round clockwise: binding,
activation (kon), unbinding, deactivation; and, from 15_5 on, its Fisher rate
against c (log scale), with the current c marked.

    scene15_1  nothing bound, no transcription; top-left state highlighted
    scene15_2  a TF binds (highlighted, "Binding"), no transcription yet; top right
    scene15_3  bound and transcribing; bottom right
    scene15_4  transcribing, the TF gone; bottom left
    scene15_5  the Fisher rate appears; c low, below its maximum (diffusion-limited):
               waiting for a TF, top left (nothing bound, no transcription)
    scene15_6  c high, past the maximum (activation-limited), more TFs: waiting
               for activation, top right (TF bound, no transcription)
    scene15_7  the same, with the Fisher rate at twice kon overlaid

The Fisher rate is Fdot(c) = F(c) / tau: the Fisher information of an OFF
time (TF binding at rate x = c kb, then activation at rate a = kon), per mean
cycle time tau, with the counter-clockwise rates all zero (SI of the SPRT
paper, eq. fisher_off):

    F(c) = 1 - alpha^2 + (2a/x) alpha^4 zeta(3, beta),
    alpha = x / (x - a),  beta = alpha if x > a, 1 - alpha if x < a,
    tau = 1/x + 1/kon + 1/ku + 1/koff.

Outputs: slides/figures/scene15/scene15_<k>.{pdf,svg}, full 128 x 96 mm slides.

Usage (from slides/):
    python3 scripts/scene15/scene15.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

import numpy as np
from matplotlib.patches import Circle
from scipy.special import zeta

from common.drawing import (AXIS_COLOR, LABEL_PT, draw_label, draw_plot_axes,
                            new_slide, save)
from common.telegraph import (K_BIND, K_UNBIND, STATE_LABELS, STATE_LABEL_X,
                              draw_state_graph)
from scene11.scene11 import draw_cartoon, molecules

# ---------------------------------------------------------------------------
# Fisher rate of the four-state model
# ---------------------------------------------------------------------------

K_ON = 1.0   # activation, while bound (a)
K_OFF = 1.0  # deactivation, once unbound


SINGULAR_WINDOW = 1e-3  # relative distance from x = a below which we interpolate


def fisher_off(x, a):
    """Fisher information on ln c of one OFF time: binding at rate x, then
    activation at rate a."""
    def closed_form(x):
        alpha = x / (x - a)
        beta = alpha if x > a else 1 - alpha
        return 1 - alpha ** 2 + 2 * a / x * alpha ** 4 * zeta(3, beta)

    if abs(x - a) < SINGULAR_WINDOW * a:
        # x = a is a removable singularity, near which the closed form loses
        # precision: interpolate across it.
        low, high = a * (1 - SINGULAR_WINDOW), a * (1 + SINGULAR_WINDOW)
        return np.interp(x, [low, high], [closed_form(low), closed_form(high)])
    return closed_form(x)


def fisher_rate(c, k_on=K_ON):
    x = c * K_BIND
    tau = 1 / x + 1 / k_on + 1 / K_UNBIND + 1 / K_OFF
    return fisher_off(x, k_on) / tau


# ---------------------------------------------------------------------------
# Slides
# ---------------------------------------------------------------------------

C_LOW, C_HIGH = 0.12, 4.0      # either side of the maximum (c = 0.5)
TF_CONCENTRATION = {C_LOW: 1.0, C_HIGH: 2.0}  # TFs drawn, as in scene 11 (relative to 1x)

# Cartoon (bound TF, transcription), highlighted state, Fisher plot: the
# concentration marked (None: no plot), and the kon of the extra curve.
SLIDES = {
    "scene15_1": dict(bound=False, transcribing=False, state="TL"),
    "scene15_2": dict(bound=True, transcribing=False, state="TR", binding=True),
    "scene15_3": dict(bound=True, transcribing=True, state="BR"),
    "scene15_4": dict(bound=False, transcribing=True, state="BL"),
    "scene15_5": dict(bound=False, transcribing=False, state="TL", c=C_LOW),
    "scene15_6": dict(bound=True, transcribing=False, state="TR", c=C_HIGH),
    "scene15_7": dict(bound=True, transcribing=False, state="TR", c=C_HIGH,
                      other_k_on=2 * K_ON),
}

# ---------------------------------------------------------------------------
# Layout (slide units)
# ---------------------------------------------------------------------------

BINDING = dict(radius=0.62, color="#FFD54A", alpha=0.45)  # as in scene 6
BINDING_LABEL_OFFSET = (0.6, 1.85)  # from the bound TF; above the DNA arc

# Four-state graph: the two-state graph's rows, two states each.
STATES = {"TL": (12.3, 13.22), "TR": (15.0, 13.22), "BL": (12.3, 10.10), "BR": (15.0, 10.10)}
STATE_LABELS_DX = 0.95  # state labels, moved right of the two-state graph's
# Both ways between neighbours (clockwise outside). Rate labels: binding and
# unbinding (rows), activation and deactivation of the bound promoter (right
# column).
CYCLE = ["TL", "TR", "BR", "BL"]
ARROWS = [arrow for a, b in zip(CYCLE, CYCLE[1:] + CYCLE[:1]) for arrow in ((a, b), (b, a))]
RATE_LABELS = {("TL", "TR"): r"$c\,k_b$", ("TR", "TL"): r"$k_u$",
               ("BL", "BR"): r"$c\,k_b$", ("BR", "BL"): r"$k_u$",
               ("TR", "BR"): r"$k_{\rm on}$", ("BR", "TR"): r"$k_{\rm off}$"}

# Fisher rate against c, on a log scale from C_RANGE[0] (at the vertical axis)
# to C_RANGE[1] (at C_RIGHT).
FISHER_ORIGIN, FISHER_TIPS = (12.0, 3.6), (19.2, 8.5)
C_RANGE = (0.03, 30.0)
C_RIGHT = 18.7
FISHER_UNIT = 16.0  # slide units per unit of Fdot (ku = 1)
OTHER_CURVE_COLOR = "#2F62C8"
# In 15_7, the curves are labelled: kon inside its hump, the other above its
# right flank, clear of the marked point.
CURVE_LABELS = [dict(c=0.5, rate=0.06, ha="center"), dict(c=2.2, rate=None, ha="left")]


def draw_four_state_graph(ax, highlight=None):
    draw_state_graph(ax, STATES, ARROWS, highlight, labels=RATE_LABELS)
    for labels in STATE_LABELS.values():
        for text, y in labels:
            draw_label(ax, (STATE_LABEL_X + STATE_LABELS_DX, y), text)


def draw_fisher_plot(ax, c, other_k_on=None):
    draw_plot_axes(ax, FISHER_ORIGIN, FISHER_TIPS, y_label=r"$\dot{\mathcal F}(c)$",
                   x_label="$c$")
    x0, y0 = FISHER_ORIGIN

    def at(c, rate):
        fraction = np.log(c / C_RANGE[0]) / np.log(C_RANGE[1] / C_RANGE[0])
        return x0 + fraction * (C_RIGHT - x0), y0 + rate * FISHER_UNIT

    cs = np.geomspace(*C_RANGE, 300)
    curves = [(K_ON, "#222222", r"$k_{\rm on}$")]
    if other_k_on:
        curves.append((other_k_on, OTHER_CURVE_COLOR, rf"${other_k_on / K_ON:g}\,k_{{\rm on}}$"))
    for (k_on, color, text), where in zip(curves, CURVE_LABELS):
        ax.plot(*at(cs, np.array([fisher_rate(ci, k_on) for ci in cs])), color=color, lw=1.6,
                zorder=4)
        if other_k_on:  # two curves: say which is which
            rate = where["rate"]
            x, y = at(where["c"], fisher_rate(where["c"], k_on) + 0.01 if rate is None else rate)
            ax.text(x, y, text, ha=where["ha"], va="bottom", fontsize=LABEL_PT, color=color)

    point = at(c, fisher_rate(c))
    ax.plot([point[0]] * 2, [y0, point[1]], color=AXIS_COLOR, lw=1.0, dashes=(3, 2), zorder=3)
    ax.add_patch(Circle(point, 0.12, fc="#222222", ec="none", zorder=5))


def draw_slide(name, mol, bound, transcribing, state=None, binding=False, c=None,
               other_k_on=None):
    fig, ax = new_slide()
    tfs = TF_CONCENTRATION.get(c, 1.0)
    draw_cartoon(ax, mol, tfs, transcribing=transcribing, bound_tfs=mol["bound"] if bound else [])
    if binding:
        tf = mol["bound"][0]
        ax.add_patch(Circle(tf, BINDING["radius"], fc=BINDING["color"], ec="none",
                            alpha=BINDING["alpha"], zorder=0.8))
        draw_label(ax, tf + np.array(BINDING_LABEL_OFFSET), "Binding")
    draw_four_state_graph(ax, state)
    if c is not None:
        draw_fisher_plot(ax, c, other_k_on)
    save(fig, "scene15", name)


if __name__ == "__main__":
    mol = molecules()
    for name, spec in SLIDES.items():
        draw_slide(name, mol, **spec)
