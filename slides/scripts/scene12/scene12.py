"""
Scene 12 -- "Let's think about what you were actually doing when you attempted
this problem. You looked at the transcription trace, and asked yourself: how
likely is it that this trace came from the low concentration versus the high
concentration? Which is the same as asking whether this ratio is larger than 1"

The storyboard's equations, typeset by LaTeX in the slides' font, with the
sketched little traces as a pictogram (\trace, see common/telegraph.py):

    scene12_1  P(trace | c_low) >? P(trace | c_high), and the ratio of the two > 1
    scene12_2  the same, naming the ratio R(t) =, and the drift-diffusion
               equation log R(t) = V t + sqrt(2D) xi

Outputs: slides/figures/scene12/scene12_<k>.{pdf,svg}, full 128 x 96 mm slides.

Usage (from slides/):
    python3 scripts/scene12/scene12.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

from common.drawing import TEXT_COLOR, new_slide, save
from common.telegraph import likelihood, use_trace_glyph

use_trace_glyph()

EQUATION_PT = 18
CENTER_X = 10.0
# Baselines of the three lines (slide units).
QUESTION_Y = 11.6
RATIO_Y = 8.0
DRIFT_Y = 4.4

LOW, HIGH = likelihood("low"), likelihood("high")
QUESTION = rf"${LOW} \overset{{\textstyle ?}}{{>}} {HIGH}$"
RATIO = rf"\dfrac{{{LOW}}}{{{HIGH}}} > 1"
DRIFT = r"$\log R(t) = V t + \sqrt{2D}\,\xi$"


def equation(ax, xy, tex, ha="center"):
    return ax.text(*xy, tex, ha=ha, va="baseline", fontsize=EQUATION_PT, color=TEXT_COLOR)


def draw_ratio(ax, named):
    """The ratio > 1, with "R(t) =" in front if `named`. The line is aligned on
    its right end, so the fraction stays put when "R(t) =" appears; it is
    centred with it."""
    named_line = equation(ax, (CENTER_X, RATIO_Y), rf"$R(t) = {RATIO}$")
    right = named_line.get_window_extent().transformed(ax.transData.inverted()).x1
    named_line.remove()
    equation(ax, (right, RATIO_Y), rf"$R(t) = {RATIO}$" if named else rf"${RATIO}$", ha="right")


def draw_slide(name, drift):
    fig, ax = new_slide()
    equation(ax, (CENTER_X, QUESTION_Y), QUESTION)
    draw_ratio(ax, named=drift)
    if drift:
        equation(ax, (CENTER_X, DRIFT_Y), DRIFT)
    save(fig, "scene12", name)


if __name__ == "__main__":
    draw_slide("scene12_1", drift=False)
    draw_slide("scene12_2", drift=True)
