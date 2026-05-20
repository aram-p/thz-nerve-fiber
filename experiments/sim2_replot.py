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
    # Drop trailing rows where σ went to the extreme diagnostic range —
    # we only want the biological-σ data here.
    sigma = arr[:, 0]
    peak_node = arr[:, 1]
    peak_global = arr[:, 2]
    mean_node = arr[:, 3]

    fig, axs = plt.subplots(1, 2, figsize=(8, 3.4))
    ax, ax_d = axs

    ax.plot(sigma, peak_node, "o-", lw=1.5, ms=5, label="peak |E| in node annulus")
    ax.plot(sigma, mean_node, "^:", lw=1.0, ms=3.5, alpha=0.75,
            label="mean |E| in node annulus")
    ax.set_xlabel(r"Node conductivity $\sigma$ (S/m)")
    ax.set_ylabel("|E| (normalised, scattered field)")
    ax.set_title("|E| vs σ")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", framealpha=0.92, fontsize=7)

    # Show the *change* in |E| from baseline more clearly
    delta = peak_node - peak_node[0]
    ax_d.plot(sigma, 1e3 * delta, "o-", color="C3", lw=1.5, ms=5)
    ax_d.set_xlabel(r"Node conductivity $\sigma$ (S/m)")
    ax_d.set_ylabel(r"$\Delta$|E| in node × 10$^3$  (relative to σ=0)")
    ax_d.set_title(r"$\Delta$|E| scaling with σ — approximately linear")
    ax_d.grid(True, alpha=0.3)
    # Linear fit guide
    if len(sigma) >= 2 and sigma.max() > 0:
        slope = delta[-1] / sigma[-1]
        sig_line = np.linspace(0, sigma.max(), 50)
        ax_d.plot(sig_line, 1e3 * slope * sig_line, "--", color="0.5",
                  alpha=0.6, lw=0.9,
                  label=fr"linear fit  $\Delta|E|/\Delta\sigma$ = {slope:.2g}")
        ax_d.legend(loc="best", framealpha=0.92, fontsize=7)

    fig.tight_layout()
    pdf, png = save_figure(fig, "sim2_sigma_sweep")
    print(f"wrote {pdf}\nwrote {png}")


if __name__ == "__main__":
    main()
