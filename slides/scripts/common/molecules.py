"""
Molecules for the gene-expression scenes: DNA (chromatin), RNA, RNA polymerase,
proteins.

DNA and RNA are `Polymer`s: chains of backbone nodes. Chunk i is the short
stretch of chain between nodes i and i+1, and every chunk has a kind ("coil"
for compacted chromatin, "linker" for open DNA, "gene", "rna", ...).

Everything drawn along a chain -- chromatin coils, the gene tube, RNA beads, a
polymerase riding on the DNA -- is computed from the node positions and placed
by node index (the chain's material coordinate), never by absolute position.
To animate thermal motion later, move the nodes (e.g. a small random-walk or
Rouse update per frame) and redraw: every detail follows its chunk.
"""

import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgb
from matplotlib.patches import Circle, PathPatch, Polygon
from matplotlib.path import Path
from scipy.interpolate import splev, splprep

from .drawing import (EDGE_PT, draw_silhouette, placer, ribbon, smooth_closed, smooth_open,
                      unit_normals)

# ---------------------------------------------------------------------------
# Palettes and sizes (drawing units unless stated)
# ---------------------------------------------------------------------------

BLUE = dict(fill="#2F62C8", edge="#1D3F85")  # same blue as the rounded cell's genes
GREEN = dict(fill="#4AA851", edge="#2B6E31")  # transcription factor, its gene and mRNA
ORANGE = dict(fill="#F39A2B", edge="#B5651D")  # promoter
RNAP_COLORS = dict(fill="#D957DB", edge="#8E1F90")
MOLECULE_EDGE_PT = 0.7 * EDGE_PT

DNA_BOND = 0.11            # node spacing along DNA
DNA_LINE = dict(back="#A5A5A5", front="#2B2B2B")
DNA_LINE_PT = (0.5, 1.1)   # line width (points) at the back / front of a coil loop
COIL_NODES_PER_LOOP = 2    # compacted chromatin: one loop every 2 nodes...
COIL_RADIUS = 0.17         # ...of this radius
TUBE_HALF_WIDTH = 0.075
TUBE_COLORS = {"gene": BLUE, "promoter": ORANGE}  # chunk kinds drawn as coloured tubes

RNA_BOND = 0.15            # node spacing along RNA = distance between beads
RNA_BEAD_RADIUS = 0.052

PROTEIN_HALF_WIDTH = 0.08
PROTEIN_RADIUS = 0.4       # outer radius of a protein drawn at scale 1


# ---------------------------------------------------------------------------
# Polymer chains
# ---------------------------------------------------------------------------

def resample(points, spacing):
    """Points evenly spaced (as close to `spacing` as fits) along a smooth curve
    through `points`, including both ends."""
    curve, s = smooth_open(points, n=400)
    length = np.linalg.norm(np.diff(curve, axis=0), axis=1).sum()
    ss = np.linspace(0, 1, max(1, round(length / spacing)) + 1)
    return np.column_stack([np.interp(ss, s, curve[:, 0]), np.interp(ss, s, curve[:, 1])])


def runs(mask):
    """(start, stop) index pairs of the contiguous True stretches of a boolean array."""
    edges = np.diff(np.r_[0, np.asarray(mask, int), 0])
    return list(zip(np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)))


def smoothstep(x):
    return x * x * (3 - 2 * x)


def circle_points(center, radius, n=24):
    """Points around a circle, e.g. to keep growing chains clear of a protein."""
    a = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.asarray(center) + radius * np.column_stack([np.cos(a), np.sin(a)])


class Polymer:
    """A chain of backbone nodes; chunk i joins nodes i and i+1 and has kind kinds[i]."""

    def __init__(self, nodes, kinds):
        self.nodes = np.asarray(nodes, float)        # (N, 2): the state to move when animating
        self.kinds = np.asarray(kinds, dtype=object)  # (N - 1,); editable, e.g. to mark a promoter

    @classmethod
    def from_trace(cls, pieces, bond):
        """Chain from traced pieces [(kind, points), ...] laid end to end.

        Each piece is resampled to nodes `bond` apart; consecutive pieces share
        their joining node.
        """
        nodes, kinds = [], []
        for kind, points in pieces:
            pts = resample(points, bond)
            nodes.extend(pts if not nodes else pts[1:])
            kinds.extend([kind] * (len(pts) - 1))
        return cls(nodes, kinds)

    def spline(self):
        u = np.arange(len(self.nodes), dtype=float)
        return splprep(self.nodes.T, u=u, s=0, k=min(3, len(self.nodes) - 1))[0]

    def curve(self, per_chunk=16):
        """Smooth curve through the nodes: (points, u), u = fractional node index."""
        u = np.linspace(0, len(self.nodes) - 1, (len(self.nodes) - 1) * per_chunk + 1)
        return np.column_stack(splev(u, self.spline())), u

    def at(self, u):
        """Position and unit tangent at fractional node index u."""
        tck = self.spline()
        d = np.array(splev(u, tck, der=1))
        return np.array(splev(u, tck)), d / np.linalg.norm(d)

    def kind_of(self, u):
        """Kind of the chunk containing each fractional node index in u."""
        return self.kinds[np.minimum(np.asarray(u).astype(int), len(self.kinds) - 1)]

    def span(self, kind):
        """(first, last) node index of the stretch of chunks of this kind."""
        idx = np.flatnonzero(self.kinds == kind)
        return float(idx[0]), float(idx[-1] + 1)

    def nearest(self, point):
        """Fractional node index of the point of the chain closest to `point`."""
        pts, u = self.curve()
        return u[np.argmin(np.linalg.norm(pts - np.asarray(point), axis=1))]


# ---------------------------------------------------------------------------
# DNA
# ---------------------------------------------------------------------------

def draw_dna(ax, dna, tube_colors=TUBE_COLORS, z=0, line_colors=DNA_LINE, line_pt=DNA_LINE_PT):
    """Chromatin: coiled over "coil" chunks, a plain line over open DNA, and a
    coloured tube over chunks whose kind is in `tube_colors` (e.g. the gene).

    Coils are loops around the backbone whose phase follows the node index, so
    they stay attached to the same chunks when the nodes move. The front half of
    each loop is drawn darker, thicker and on top, for a 3D coil look; its
    colours and widths are `line_colors` and `line_pt` (lighter colours dim the
    chromatin; thinner lines keep the coils open on a zoomed-out drawing).
    """
    center, u = dna.curve()
    n = unit_normals(center)
    t = np.column_stack([n[:, 1], -n[:, 0]])
    kinds = dna.kind_of(u)

    # Loop radius tapers to 0 over one loop at each end of a coiled stretch, so
    # coils grow smoothly out of the open DNA and out of the chain's free ends.
    radius = np.zeros(len(u))
    for a, b in runs(kinds == "coil"):
        from_end = np.minimum(u[a:b] - u[a], u[b - 1] - u[a:b]) / COIL_NODES_PER_LOOP
        radius[a:b] = COIL_RADIUS * smoothstep(np.clip(from_end, 0, 1))
    phase = 2 * np.pi * u / COIL_NODES_PER_LOOP
    path = center + radius[:, None] * (-np.sin(phase)[:, None] * t + np.cos(phase)[:, None] * n)
    depth = 1 - radius / COIL_RADIUS * 0.5 * (1 + np.sin(phase))  # 0 = back, 1 = front

    segments = np.stack([path[:-1], path[1:]], axis=1)
    seg_depth = 0.5 * (depth[:-1] + depth[1:])
    order = np.argsort(seg_depth, kind="stable")  # back first, front on top
    back, front = np.array(to_rgb(line_colors["back"])), np.array(to_rgb(line_colors["front"]))
    ax.add_collection(LineCollection(
        segments[order], colors=back + (front - back) * seg_depth[order, None],
        linewidths=line_pt[0] + (line_pt[1] - line_pt[0]) * seg_depth[order],
        capstyle="round", zorder=z))

    for kind, colors in tube_colors.items():
        for a, b in runs(kinds == kind):
            piece = center[a:b + 1]
            draw_silhouette(ax, [ribbon(piece, np.full(len(piece), TUBE_HALF_WIDTH))],
                            colors["fill"], colors["edge"], z + 1, edge_pt=MOLECULE_EDGE_PT)


# ---------------------------------------------------------------------------
# RNA
# ---------------------------------------------------------------------------

def draw_rna(ax, rna, color=BLUE["fill"], z=0):
    """RNA as a string of beads, one bead per node."""
    for p in rna.nodes:
        ax.add_patch(Circle(p, RNA_BEAD_RADIUS, fc=color, ec="none", zorder=z))


def grow_rna(start, direction, n_beads, rng, avoid=(), clearance=0.25, wiggle=0.45,
             pull=0.12, tries=40):
    """A transcript as a random, self-avoiding worm-like chain of `n_beads` beads,
    growing from `start` (the polymerase) roughly along `direction`.

    Each step turns by a random angle (`wiggle`) with a weak `pull` back to the
    initial heading; steps that land within `clearance` of the `avoid` points
    (e.g. the DNA) or on the chain itself are redrawn, with wider and wider
    turns, so a blocked chain bends around the obstacle.
    """
    avoid = np.asarray(avoid, float).reshape(-1, 2)
    goal = heading = np.arctan2(direction[1], direction[0])
    nodes = [np.asarray(start, float)]
    for i in range(n_beads - 1):
        for attempt in range(tries):
            spread = wiggle * (1 + attempt / 4)
            turn = heading + spread * rng.standard_normal() + pull * (goal - heading)
            p = nodes[-1] + RNA_BOND * np.array([np.cos(turn), np.sin(turn)])
            clear = (i < 3 or len(avoid) == 0
                     or np.linalg.norm(avoid - p, axis=1).min() > clearance)
            self_clear = (len(nodes) < 3
                          or np.linalg.norm(np.array(nodes[:-2]) - p, axis=1).min() > 1.5 * RNA_BOND)
            if clear and self_clear:
                break
        heading = turn
        nodes.append(p)
    return Polymer(nodes, ["rna"] * (len(nodes) - 1))


# ---------------------------------------------------------------------------
# Proteins
# ---------------------------------------------------------------------------

# RNA polymerase: a lumpy globule around its centre, in a frame where the DNA
# under it runs along +x; it turns with the DNA as it moves along it.
RNAP_OUTLINE = [(0.37, 0.0), (0.26, 0.19), (0.09, 0.29), (-0.11, 0.33), (-0.27, 0.19),
                (-0.32, 0.0), (-0.29, -0.21), (-0.1, -0.31), (0.09, -0.28), (0.28, -0.2)]


def draw_rnap_at(ax, center, angle=0.0, colors=RNAP_COLORS, scale=1.0, alpha=1.0, z=0):
    """RNA polymerase at `center`, turned by `angle` degrees (e.g. floating free)."""
    blob = placer(center, scale, angle)(smooth_closed(np.array(RNAP_OUTLINE)))
    if alpha < 1:  # one translucent patch (the silhouette trick would show its layers)
        ax.add_patch(Polygon(blob, closed=True, fc=colors["fill"], ec=colors["edge"],
                             lw=MOLECULE_EDGE_PT, alpha=alpha, zorder=z))
    else:
        draw_silhouette(ax, [blob], colors["fill"], colors["edge"], z, edge_pt=MOLECULE_EDGE_PT)


def draw_rnap(ax, dna, u, colors=RNAP_COLORS, z=0):
    """RNA polymerase sitting on the DNA at fractional node index u."""
    p, t = dna.at(u)
    draw_rnap_at(ax, p, np.degrees(np.arctan2(t[1], t[0])), colors=colors, z=z)
    return p


# A folded protein, drawn as in the storyboard: a thick closed loop.
PROTEIN_TRACE = [(-0.26, 0.17), (-0.07, 0.29), (0.15, 0.29), (0.24, 0.09), (0.25, -0.14),
                 (0.12, -0.28), (-0.13, -0.29), (-0.3, -0.17)]


def draw_protein(ax, center, colors=BLUE, scale=1.0, angle=0.0, alpha=1.0, z=0):
    loop = placer(center, scale, angle)(smooth_closed(np.array(PROTEIN_TRACE)))[:-1]
    n = unit_normals(np.vstack([loop[-1:], loop, loop[:1]]))[1:-1]  # periodic normals
    outer, inner = loop + n * PROTEIN_HALF_WIDTH * scale, loop - n * PROTEIN_HALF_WIDTH * scale
    # Outer and inner edges with opposite orientation: the inner one cuts a hole.
    verts = np.vstack([outer, outer[:1], inner[::-1], inner[-1:]])
    codes = ([Path.MOVETO] + [Path.LINETO] * (len(outer) - 1) + [Path.CLOSEPOLY]) * 2
    ax.add_patch(PathPatch(Path(verts, codes), fc=colors["fill"], ec=colors["edge"],
                           lw=MOLECULE_EDGE_PT, joinstyle="round", alpha=alpha, zorder=z))


# ---------------------------------------------------------------------------
# Proteins scattered in 3D
# ---------------------------------------------------------------------------
# Each protein gets a depth d in [0, 1] (0 = far back, 1 = front): far proteins
# are drawn smaller and fainter, and behind the DNA.

def depth_scale(d):
    return 0.6 + 0.5 * d


def depth_alpha(d):
    return 0.3 + 0.7 * d


def scatter_proteins(n, region, rng, avoid=(), existing=(), gap=0.1, tries=20000,
                     radius=PROTEIN_RADIUS):
    """`n` proteins at random positions, depths and orientations.

    `region` is (x0, x1, y0, y1); `avoid` is a list of (points, radius) pairs
    that proteins must stay clear of (DNA, polymerases, transcripts...).
    `existing` proteins (an earlier scatter) are kept, and new ones added to them.
    `radius` is the size of one protein at depth scale 1.
    Returns a list of dicts (center, depth, angle), ordered from back to front.
    """
    x0, x1, y0, y1 = region
    placed = list(existing)
    for _ in range(tries):
        if len(placed) == n:
            break
        d = rng.uniform()
        r = radius * depth_scale(d)
        c = np.array([rng.uniform(x0 + r, x1 - r), rng.uniform(y0 + r, y1 - r)])
        if all(np.linalg.norm(np.asarray(pts).reshape(-1, 2) - c, axis=1).min() > rad + r + gap
               for pts, rad in avoid) and \
           all(np.linalg.norm(q["center"] - c) > r + radius * depth_scale(q["depth"]) + gap
               for q in placed):
            placed.append(dict(center=c, depth=d, angle=rng.uniform(0, 360)))
    return sorted(placed, key=lambda q: q["depth"])


def draw_scattered_proteins(ax, proteins, colors=GREEN, z_back=0, z_front=10):
    """Draw scatter_proteins() output: the far half behind z_back, the near half at z_front."""
    for q in proteins:
        draw_protein(ax, q["center"], colors, scale=depth_scale(q["depth"]), angle=q["angle"],
                     alpha=depth_alpha(q["depth"]), z=z_back if q["depth"] < 0.5 else z_front)
