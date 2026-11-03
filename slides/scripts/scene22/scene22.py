"""
Scene 22 -- "We wondered what a gene regulation mechanism that accounts for
this biological reality would look like. I started by building intuition with
lattice simulations of polymers in 3D space interacting with dilute molecules,
representing regulatory proteins. [...]"

Panels B and C of the paper's Fig. 1, side by side: the lattice polymer (green)
with its promoter (orange) and enhancer (red), among ligands (blue) that are
dilute (B, c = 0.003), or condensed on the regulatory regions, holding promoter
and enhancer together (C, c = 0.017).

They are the published plots as extracted by figure_production/fig1/
extract_art.py (art_b.pdf, art_c.pdf), which carry no lettering: the axis and
tick labels that Fig. 1 sets on top of them are left out. Each is cropped to its
ink (Ghostscript's bounding box) and placed by pdflatex, so the plots stay
vector.

Output: slides/figures/scene22/scene22.pdf, a full 128 x 96 mm slide.

Usage (from slides/):
    python3 scripts/scene22/scene22.py
Requires the `gs` and `pdflatex` command-line tools.
"""

import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

from common.drawing import FIGURES_DIR, SLIDE_MM, UNIT_MM

FIG1_DIR = Path(__file__).resolve().parents[3] / "figure_production" / "fig1"
PT_MM = 25.4 / 72

# Top-left corner of each panel's ink (slide units, as in the storyboard). Both
# are drawn at one scale, the one that makes B PANEL_WIDTH wide.
PANELS = {"art_b.pdf": (1.1, 13.4), "art_c.pdf": (10.9, 13.4)}
PANEL_WIDTH = 8.6

TEX = r"""\documentclass[border=0pt]{standalone}
\usepackage{graphicx}
\usepackage{tikz}
\begin{document}
\begin{tikzpicture}[x=1mm, y=1mm]
\useasboundingbox (0, 0) rectangle (%(width)s, %(height)s);
%(nodes)s
\end{tikzpicture}
\end{document}
"""
NODE = (r"\node[anchor=north west, inner sep=0] at (%.3f, %.3f) "
        r"{\includegraphics[scale=%.5f, viewport=%.3f %.3f %.3f %.3f, clip]{%s}};")


def ink_box(pdf):
    """(x0, y0, x1, y1), in points, of everything drawn on the page."""
    out = subprocess.run(["gs", "-q", "-dBATCH", "-dNOPAUSE", "-sDEVICE=bbox", str(pdf)],
                         capture_output=True, text=True, check=True).stderr
    return tuple(map(float, re.search(r"%%HiResBoundingBox: (.*)", out).group(1).split()))


def typeset(nodes, scene):
    """Compile TikZ `nodes` on a full slide into figures/<scene>/<scene>.pdf."""
    out = FIGURES_DIR / scene / f"{scene}.pdf"
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        tex = Path(tmp) / f"{scene}.tex"
        tex.write_text(TEX % dict(width=SLIDE_MM[0], height=SLIDE_MM[1], nodes="\n".join(nodes)))
        subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", tex.name],
                       cwd=tmp, check=True, capture_output=True)
        shutil.copy(tex.with_suffix(".pdf"), out)
    print(f"wrote {out}")


if __name__ == "__main__":
    boxes = {name: ink_box(FIG1_DIR / name) for name in PANELS}
    x0, _, x1, _ = boxes["art_b.pdf"]
    scale = PANEL_WIDTH * UNIT_MM / ((x1 - x0) * PT_MM)
    nodes = [NODE % (x * UNIT_MM, y * UNIT_MM, scale, *boxes[name], FIG1_DIR / name)
             for name, (x, y) in PANELS.items()]
    typeset(nodes, "scene22")
