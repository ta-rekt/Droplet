"""
Scene 23 -- "multi-scale temporal information in transcription activity"
(storyboard note).

Panel F of the paper's Fig. 1, its two transcription time series put on top of
each other as published: the polymer's probe (top), switching with the polymer
between closed (red) and open (green) states, and the free receptor (bottom).
Scene 22, in miniature, stays at the top of the frame.

The plots are cropped from the published figure (slides/figures/scene23/
fig1.pdf) and placed by pdflatex, so they stay vector and keep the figure's own
lettering, drawn at the slides' label size. The panel letter F is left out. The
miniature is scene 22's output, cropped to its ink.

Output: slides/figures/scene23/scene23.pdf, a full 128 x 96 mm slide.

Usage (from slides/):
    python3 scripts/scene22/scene22.py   # the miniature, if not yet rendered
    python3 scripts/scene23/scene23.py
Requires the `gs` and `pdflatex` command-line tools.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

from common.drawing import FIGURES_DIR, LABEL_PT, SLIDE_UNITS, UNIT_MM
from scene22.scene22 import NODE, PT_MM, ink_box, typeset

FIG1 = FIGURES_DIR / "scene23" / "fig1.pdf"

# Scene 22 in miniature, centred at the top of the frame: top edge and height
# of its ink (slide units).
MINIATURE = dict(top=14.6, height=4.0)

# Panel F's plots share one geometry, the free receptor's 243 pt right of the
# polymer's, so one viewport of fig1.pdf (x0, y0, x1, y1 in points, y up),
# shifted, crops each: from the y-label to the right spine, and from the line
# of labels above the axes to the T's under them. The panel letter, above the
# polymer's y-label, is outside it.
POLYMER_BOX = (4.0, 4.6, 230.8, 73.2)
RECEPTOR_SHIFT = 243.0

# Fig. 1's lettering is 9 pt; this scale draws it at the slides' label size.
FIG_SCALE = LABEL_PT / 9
PLOTS_TOP = 10.2   # top edge of the polymer's plot (slide units)
PLOTS_GAP = 0.3    # between the two plots


def miniature():
    """TikZ node: scene 22, cropped to its ink, where MINIATURE places it."""
    pdf = FIGURES_DIR / "scene22" / "scene22.pdf"
    box = ink_box(pdf)
    scale = MINIATURE["height"] * UNIT_MM / ((box[3] - box[1]) * PT_MM)
    width = (box[2] - box[0]) * scale * PT_MM / UNIT_MM
    return NODE % ((SLIDE_UNITS[0] - width) / 2 * UNIT_MM, MINIATURE["top"] * UNIT_MM,
                   scale, *box, pdf)


if __name__ == "__main__":
    x0, y0, x1, y1 = POLYMER_BOX
    width = (x1 - x0) * FIG_SCALE * PT_MM / UNIT_MM
    height = (y1 - y0) * FIG_SCALE * PT_MM / UNIT_MM
    left = (SLIDE_UNITS[0] - width) / 2
    nodes = [miniature()]
    for i, shift in enumerate((0, RECEPTOR_SHIFT)):
        top = PLOTS_TOP - i * (height + PLOTS_GAP)
        nodes.append(NODE % (left * UNIT_MM, top * UNIT_MM, FIG_SCALE,
                             x0 + shift, y0, x1 + shift, y1, FIG1))
    typeset(nodes, "scene23")
