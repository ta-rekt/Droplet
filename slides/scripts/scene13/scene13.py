"""
Scene 13 -- "Let's plot this ratio in real time, while transcription is
happening, and vary concentration. [...] a larger concentration difference is
easier to tell apart. Larger absolute concentrations are also easier to
distinguish. [...]"

A transcription time series simulated at concentration c_low = c (top left);
below it, the log-likelihood ratio R(t) of that very trace, low vs high
concentration c_high = c + dc, with its trend V t; the two-state graph (top
right); the slope V against c (bottom right), with the current c marked and
the value it saturates at, 1/2 (dc/c)^2 ku, dashed; and the concentration-
difference slider, at c_high / c_low = 1 + dc/c.

    scene13_1  c = 1, dc/c = 1 (a 2x difference)
    scene13_2  c = 1, dc/c = 0.5 (1.5x): same trace, shallower trend
    scene13_3  low c (dc/c = 1): fewer bindings, shallower trend
    scene13_4  high c (dc/c = 1): steeper, close to saturation

V is plotted as 1/2 (dc/c)^2 c kb ku / (c kb + ku), the slope for small dc. All
slides use the same random seed: the first one whose log R(t) follows V t on
every slide (least-squares slope within SLOPE_TOLERANCE of V, and RMS distance
from V t within DISTANCE_TOLERANCE of its end value).

Outputs: slides/figures/scene13/scene13_<k>.{pdf,svg}, full 128 x 96 mm slides.

Usage (from slides/):
    python3 scripts/scene13/scene13.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

import numpy as np
from matplotlib.patches import Circle

from common.drawing import LABEL_PT, TEXT_COLOR, draw_arrow, draw_label, new_slide, save
from common.telegraph import (AXIS_COLOR, K_UNBIND, TimeAxes, draw_slider, draw_time_series,
                              draw_two_state_graph, first_seed, likelihood, log_likelihood_ratio,
                              slope, telegraph, use_trace_glyph)

use_trace_glyph()

# ---------------------------------------------------------------------------
# Slides and their traces
# ---------------------------------------------------------------------------

DURATION = 48.0  # time shown, as in 11_5-6
SLIDES = {
    "scene13_1": dict(c=1.0, dc_over_c=1.0),
    "scene13_2": dict(c=1.0, dc_over_c=0.5),
    "scene13_3": dict(c=0.3, dc_over_c=1.0),
    "scene13_4": dict(c=3.0, dc_over_c=1.0),
}
SLOPE_TOLERANCE = 0.2
DISTANCE_TOLERANCE = 0.16


def ratio(seed, c, dc_over_c):
    """The trace at c, and its log R(t) (low c vs high c + dc)."""
    trace = telegraph(c, DURATION, seed)
    return trace, log_likelihood_ratio(trace, c, c * (1 + dc_over_c))


def follows_trend(seed, c, dc_over_c):
    t, log_r = ratio(seed, c, dc_over_c)[1]
    grid = np.linspace(0, DURATION, 1000)
    on_grid = np.interp(grid, t, log_r)
    v = slope(c, dc_over_c)
    fitted = on_grid @ grid / (grid @ grid)  # least squares, through the origin
    distance = np.sqrt(np.mean((on_grid - v * grid) ** 2))
    return (abs(fitted / v - 1) < SLOPE_TOLERANCE
            and distance < DISTANCE_TOLERANCE * v * DURATION)


# ---------------------------------------------------------------------------
# Layout (slide units)
# ---------------------------------------------------------------------------

ACCENT = "#D9532B"  # the slope V: the trend line and its point on V(c)

# Left: the time series, and log R(t) below it, on the same time axis.
TIME_SERIES_AXES = TimeAxes(origin=(4.3, 10.7), tips=(11.2, 13.07), span=(4.37, 10.8))
TRACE_PT = 0.6
RATIO_AXES = TimeAxes(origin=(4.3, 2.6), tips=(11.2, 8.6), span=(4.37, 10.8))
LOG_R_UNIT = 0.3      # slide units per unit of log R (natural log)
LOG_R_BOTTOM = -2.0   # where the vertical axis starts
LOG_R_TICKS = [(0, "$1$"), (4, "$10^4$"), (8, "$10^8$")]  # powers of 10 of R
RATIO_PT = 1.0

# Right: the two-state graph (moved from its scene-11 place) and V(c) below it.
GRAPH_OFFSET = (0.9, 0.0)
SLOPE_ORIGIN = (13.2, 2.6)
SLOPE_TIPS = (19.0, 7.9)
C_MAX = 4.0           # c at the right end of the curve...
C_RIGHT = 18.5        # ...drawn at this x
V_UNIT = 8.0          # slide units per unit of V

# The slider, left of the plots, between them.
SLIDER_OFFSET = (-0.17, -3.0)


# ---------------------------------------------------------------------------
# Panels
# ---------------------------------------------------------------------------

def draw_ratio(ax, t, log_r, v):
    """log R(t) on a log scale of R, with the trend V t."""
    (x0, y0), (x_tip, y_tip) = RATIO_AXES.origin, RATIO_AXES.tips
    draw_arrow(ax, (x0, y0), (x_tip, y0), AXIS_COLOR, z=3, lw=1.3, head=9)
    draw_arrow(ax, (x0, y0 + LOG_R_BOTTOM * LOG_R_UNIT), (x0, y_tip), AXIS_COLOR, z=3, lw=1.3,
               head=9)
    draw_label(ax, (x_tip - 0.15, y0 - 0.5), "time")
    for power, text in LOG_R_TICKS:
        y = y0 + power * np.log(10) * LOG_R_UNIT
        ax.plot([x0 - 0.1, x0], [y, y], color=AXIS_COLOR, lw=1.3, zorder=3)
        draw_label(ax, (x0 - 0.2, y), text, ha="right")
    draw_label(ax, (x0 + 0.27, y_tip), rf"$\dfrac{{{likelihood('low')}}}{{{likelihood('high')}}}$",
               ha="left")

    ax.plot(RATIO_AXES.x(t, DURATION), y0 + log_r * LOG_R_UNIT, color="#222222", lw=RATIO_PT,
            solid_joinstyle="miter", zorder=4)
    end = (RATIO_AXES.span[1], y0 + v * DURATION * LOG_R_UNIT)
    ax.plot([RATIO_AXES.span[0], end[0]], [y0, end[1]], color=ACCENT, lw=1.8,
            solid_capstyle="round", zorder=5)
    ax.text(end[0] + 0.15, end[1], "$Vt$", ha="left", va="center", fontsize=LABEL_PT + 2,
            color=ACCENT)


def draw_slope(ax, c, dc_over_c):
    """V against c, saturating at 1/2 (dc/c)^2 ku (dashed), with c marked."""
    (x0, y0), (x_tip, y_tip) = SLOPE_ORIGIN, SLOPE_TIPS
    draw_arrow(ax, (x0, y0), (x_tip, y0), AXIS_COLOR, z=3, lw=1.3, head=9)
    draw_arrow(ax, (x0, y0), (x0, y_tip), AXIS_COLOR, z=3, lw=1.3, head=9)
    draw_label(ax, (x0 + 0.27, y_tip), r"slope $V$", ha="left")
    draw_label(ax, (x_tip + 0.25, y0), "$c$", ha="left")

    def at(c, v):
        return x0 + c / C_MAX * (C_RIGHT - x0), y0 + v * V_UNIT

    cs = np.linspace(0, C_MAX, 200)
    ax.plot(*at(cs, slope(cs, dc_over_c)), color="#222222", lw=1.6, zorder=4)
    saturation = at(C_MAX, 0.5 * dc_over_c ** 2 * K_UNBIND)[1]
    ax.plot([x0, C_RIGHT], [saturation] * 2, color=AXIS_COLOR, lw=1.1, dashes=(4, 3), zorder=3)

    point = at(c, slope(c, dc_over_c))
    ax.plot([point[0]] * 2, [y0, point[1]], color=ACCENT, lw=1.1, dashes=(3, 2), zorder=3)
    ax.add_patch(Circle(point, 0.13, fc=ACCENT, ec="none", zorder=5))


def draw_slide(name, c, dc_over_c, seed):
    fig, ax = new_slide()
    trace, (t, log_r) = ratio(seed, c, dc_over_c)
    draw_time_series(ax, trace, DURATION, TIME_SERIES_AXES, lw=TRACE_PT, time_label=False)
    draw_ratio(ax, t, log_r, slope(c, dc_over_c))
    draw_two_state_graph(ax, offset=GRAPH_OFFSET)
    draw_slope(ax, c, dc_over_c)
    draw_slider(ax, 1 + dc_over_c, offset=SLIDER_OFFSET)
    save(fig, "scene13", name)


if __name__ == "__main__":
    seed = first_seed(lambda s: all(follows_trend(s, **spec) for spec in SLIDES.values()), 1)
    print(f"seed {seed}")
    for name, spec in SLIDES.items():
        draw_slide(name, seed=seed, **spec)
