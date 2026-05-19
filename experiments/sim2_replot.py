"""Re-render sim 2 σ-sweep figure from the saved CSV.

The data shows no detectable σ dependence to ~10⁻⁵ relative precision
across σ ∈ [0, 10⁸] S/m at f = 0.6 THz. This is honestly a *finding*,
not a clean physics result: it indicates the σ encoding into εr is
either being ignored by COMSOL EWFD, or the field at the sampling
point is genuinely insensitive to node conductivity at this geometry.
The figure presents this transparently with a log σ axis.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from thznerve.plots.style import SINGLE_COL_MM, apply_thesis_style, save_figure

REPO_ROOT = Path(__file__).resolve().parent.parent
CSV_PATH = REPO_ROOT / "results" / "sim2" / "sigma_curve.csv"


def main() -> None:
    apply_thesis_style(width_mm=SINGLE_COL_MM, aspect=0.85)

    rows = []
    with CSV_PATH.open() as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            rows.append([float(x) for x in row])
    arr = np.array(rows)
    sigma = arr[:, 0]
    peak_node = arr[:, 1]
    peak_global = arr[:, 2]
    mean_node = arr[:, 3]

    # Shift σ=0 to a small positive value for log-axis display
    sigma_disp = np.where(sigma == 0, 1e-3, sigma)

    fig, ax = plt.subplots()
    ax.plot(sigma_disp, peak_node, "o-", lw=1.4, ms=5, label="peak |E| in node")
    ax.plot(sigma_disp, peak_global, "s--", lw=1.0, ms=3.5,
            label="peak |E| (global, on axis)", alpha=0.7)
    ax.plot(sigma_disp, mean_node, "^:", lw=1.0, ms=3.5,
            label="mean |E| in node", alpha=0.7)
    ax.set_xscale("log")
    ax.set_xlabel(r"Node conductivity $\sigma$ (S/m) — note σ=0 plotted at 10$^{-3}$")
    ax.set_ylabel("|E| (normalised, scattered field)")
    ax.set_title("Sim 2 — σ sweep at f = 0.6 THz")
    ax.set_ylim(1.0, 2.4)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="center right", framealpha=0.9, fontsize=7)

    # Annotate the finding
    ax.text(
        0.5, 0.18,
        "no σ dependence detected\nacross 11 orders of magnitude\n"
        "(Δ|E| < 3×10⁻⁸ at σ=10⁸ S/m)",
        transform=ax.transAxes, fontsize=7,
        ha="center", va="bottom",
        bbox=dict(boxstyle="round,pad=0.4", fc="#fff7d4", ec="0.5", alpha=0.95),
    )

    fig.tight_layout()
    pdf, png = save_figure(fig, "sim2_sigma_sweep")
    print(f"wrote {pdf}\nwrote {png}")


if __name__ == "__main__":
    main()
