"""Thesis figure style — one place that decides what figures look like.

`apply_thesis_style()` sets matplotlib rcParams. `save_figure()` writes
both `<name>.pdf` (vector, for inclusion in the thesis) and `<name>.png`
(raster, 300 DPI, for quick previews) under `figures/out/`.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

# Okabe & Ito colorblind-safe palette (8 distinguishable colors).
OKABE_ITO: tuple[str, ...] = (
    "#000000",  # black
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
)

# Thesis page widths (millimetres). Single column ~ 90 mm, double ~ 180 mm.
SINGLE_COL_MM: float = 90.0
DOUBLE_COL_MM: float = 180.0


def _mm_to_in(mm: float) -> float:
    return mm / 25.4


def apply_thesis_style(width_mm: float = SINGLE_COL_MM, aspect: float = 3 / 4) -> None:
    """Apply the thesis-wide matplotlib style.

    Args:
        width_mm: Figure width in millimetres. 90 mm = single column, 180 mm = double.
        aspect: height / width ratio.
    """

    width_in = _mm_to_in(width_mm)
    height_in = width_in * aspect

    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 10,
            "legend.fontsize": 9,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "figure.figsize": (width_in, height_in),
            "figure.dpi": 100,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "axes.prop_cycle": mpl.cycler(color=list(OKABE_ITO)),
            # LaTeX rendering off — Windows LaTeX toolchains are fragile; revisit later.
            "text.usetex": False,
            # TrueType (Type 42) so PDF text stays editable in Illustrator / Inkscape.
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save_figure(fig: "plt.Figure", name: str, out_dir: Path | str | None = None) -> tuple[Path, Path]:
    """Write `<name>.pdf` and `<name>.png` to `out_dir` (default `figures/out/`)."""

    out_dir = Path(out_dir) if out_dir is not None else Path("figures") / "out"
    out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / f"{name}.pdf"
    png_path = out_dir / f"{name}.png"
    fig.savefig(pdf_path)
    fig.savefig(png_path)
    return pdf_path, png_path
