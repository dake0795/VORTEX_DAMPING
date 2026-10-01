r"""Figure style of the vortex-damping working notes (Dan, 1 Oct 2026).

The style is the long dipole paper's (ep_turbulence_paper/scripts/lib/paper_style.py, imported here, not copied),
with one departure that Dan asked for: axis labels, tick labels, legends and panel tags are ALL the body size of the
notes (10 pt in jpp.cls), where the long paper's module uses 9 pt labels and 8 pt ticks and legends.

Rules (also in STYLE_GUIDE.md, section 7a):
  * every figure is made at its final size, TEXTWIDTH_IN wide, and included with a bare \includegraphics{...};
  * at most TWO panels per row;
  * all text 10 pt (FS); Computer Modern through usetex;
  * panel tags "(a)", "(b)" inside the axes, top left, through tag();
  * in-panel legends frameless on white; a legend that does not fit goes below the figure, boxed (legend_below);
  * colours from PAL; lines 1.1 pt; no titles; units in the axis labels;
  * one directory per figure, figures/<name>/figure.py: `--rebuild` reads the run data and writes cache.npz
    (small, committed), the default reads cache.npz alone and writes figures/<name>.pdf.
"""
import os
import sys

LONG_PAPER_LIB = "/rds/project/rds-aSo1XX0UOlw/ir-kenn3/ep_turbulence_paper/scripts/lib"
sys.path.insert(0, LONG_PAPER_LIB)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import paper_style as _ps
from paper_style import TEXTWIDTH_IN, HALFWIDTH_IN, PANEL_H_FULL, PANEL_H_HALF, style_axes, legend_below  # noqa: F401

FS = 10                                   # body size of jpp.cls
FIG_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))     # <repo>/figures
RB = "/rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925"

# one colour per run / quantity, used the same way in every figure
PAL = {"orig": "k", "a03": "#1f77b4", "a003": "#d62728",          # original run, eps = 0.3, eps = 0.03 restarts
       "E": "k", "Z": "#d62728", "W": "#1f77b4", "zon": "0.55",
       "plus": "#d62728", "minus": "#1f77b4", "theory": "#2ca02c", "aux": "#ff7f0e"}


def apply_style():
    _ps.apply_style()
    plt.rcParams.update({"font.size": FS, "axes.titlesize": FS, "axes.labelsize": FS, "legend.fontsize": FS,
                         "xtick.labelsize": FS, "ytick.labelsize": FS, "legend.frameon": False,
                         "legend.handlelength": 1.6, "legend.borderaxespad": 0.3, "legend.labelspacing": 0.3})


def two_panel(height=2.55, left=0.62, right=0.08, bottom=0.50, top=0.10, wspace=0.72, sharey=False):
    """One row of two panels at TEXTWIDTH_IN; margins and gap in inches.  Returns fig, (ax_a, ax_b)."""
    fig = plt.figure(figsize=(TEXTWIDTH_IN, height))
    w = (TEXTWIDTH_IN - left - right - wspace) / 2.0
    h = height - bottom - top
    axa = fig.add_axes([left / TEXTWIDTH_IN, bottom / height, w / TEXTWIDTH_IN, h / height])
    axb = fig.add_axes([(left + w + wspace) / TEXTWIDTH_IN, bottom / height, w / TEXTWIDTH_IN, h / height],
                       sharey=axa if sharey else None)
    return fig, (axa, axb)


def grid_2xN(nrows, height=None, left=0.62, right=0.08, bottom=0.50, top=0.10, wspace=0.72, hspace=0.55, panel_h=1.75):
    """nrows rows of two panels at TEXTWIDTH_IN.  Returns fig, axes[nrows][2]."""
    if height is None:
        height = top + bottom + nrows * panel_h + (nrows - 1) * hspace
    panel_h = (height - top - bottom - (nrows - 1) * hspace) / nrows
    fig = plt.figure(figsize=(TEXTWIDTH_IN, height))
    w = (TEXTWIDTH_IN - left - right - wspace) / 2.0
    axes = []
    for r in range(nrows):
        y0 = bottom + (nrows - 1 - r) * (panel_h + hspace)
        axes.append([fig.add_axes([(left + c * (w + wspace)) / TEXTWIDTH_IN, y0 / height, w / TEXTWIDTH_IN,
                                   panel_h / height]) for c in range(2)])
    return fig, axes


def tag(ax, s, x=0.04, y=0.95, **kw):
    """Panel tag '(a)' inside the axes, top left, on an opaque white patch."""
    opts = dict(transform=ax.transAxes, ha="left", va="top", fontsize=FS, zorder=8,
                bbox=dict(facecolor="white", edgecolor="none", pad=1.2))
    opts.update(kw)
    return ax.text(x, y, s, **opts)


def save(fig, name):
    """figures/<name>.pdf at the exact figure size (no tight bbox) and a PNG beside it for quick viewing."""
    path = os.path.join(FIG_DIR, name + ".pdf")
    fig.savefig(path)
    fig.savefig(os.path.join(FIG_DIR, name + ".png"), dpi=200)
    print("saved", path)
    return path
