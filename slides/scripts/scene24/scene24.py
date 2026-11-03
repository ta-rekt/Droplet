"""
Scene 24 -- Panel E of the paper's Fig. 1: one stretch of the polymer's
trajectory, its smoothed transcription activity (black) against the
probe-enhancer distance (orange); transcription rises while the two are in
contact. Scene 22's miniature stays at the top of the frame, as in scene 23.

The panel is cropped from the published figure (slides/figures/scene24/
fig1.pdf) and placed by pdflatex, as in scene 23, keeping its own lettering.
The panel letter E is left out.

Output: slides/figures/scene24/scene24.pdf, a full 128 x 96 mm slide.

Usage (from slides/):
    python3 scripts/scene22/scene22.py   # the miniature, if not yet rendered
    python3 scripts/scene24/scene24.py
Requires the `gs` and `pdflatex` command-line tools.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

from common.drawing import FIGURES_DIR, SLIDE_UNITS, UNIT_MM
from scene22.scene22 import NODE, PT_MM, typeset
from scene23.scene23 import PLOTS_TOP, miniature

FIG1 = FIGURES_DIR / "scene24" / "fig1.pdf"

# Panel E, a viewport of fig1.pdf (x0, y0, x1, y1 in points, y up): from one
# y-label to the other, and from the top tick labels to the x-label.
PANEL_BOX = (295.3, 79.9, 476.6, 192.0)
# The panel letter E stands left of the y-label, its lowest stroke level with
# the top of the axes, and the end of that stroke falls inside the viewport: it
# is painted white (box in points, as PANEL_BOX; nothing else is drawn there).
LETTER_MASK = (295.3, 189.3, 296.5, 192.0)
MASK = r"\fill[white] (%.3f, %.3f) rectangle (%.3f, %.3f);"

# Centred, its top edge where scene 23's plots start (slide units).
PANEL = dict(top=PLOTS_TOP, width=15.3)


if __name__ == "__main__":
    x0, y0, x1, y1 = PANEL_BOX
    scale = PANEL["width"] * UNIT_MM / ((x1 - x0) * PT_MM)
    left, top = (SLIDE_UNITS[0] - PANEL["width"]) / 2 * UNIT_MM, PANEL["top"] * UNIT_MM

    def to_mm(x, y):
        """Slide position (mm) of a point of fig1.pdf (points)."""
        return left + (x - x0) * scale * PT_MM, top - (y1 - y) * scale * PT_MM

    nodes = [miniature(),
             NODE % (left, top, scale, *PANEL_BOX, FIG1),
             MASK % (*to_mm(*LETTER_MASK[:2]), *to_mm(*LETTER_MASK[2:]))]
    typeset(nodes, "scene24")
