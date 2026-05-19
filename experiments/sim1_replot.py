"""Re-render sim 1's frequency-sweep figure from the saved CSV.

Loads results/sim1/spectrum.csv produced by sim1_freq_sweep.py and
re-plots with a smoothed trend line and an annotated peak. Lets the
plotting iteration happen without re-running the 5-minute COMSOL sweep.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from thznerve.plots.style import SINGLE_COL_MM, apply_thesis_style, save_figure

REPO_ROOT = Path(__file__).resolve().parent.parent
CSV_PATH = REPO_ROOT / "results" / "sim1" / "spectrum.csv"


def _smooth(y: np.ndarray, window: int = 3) -> np.ndarray:
    """Centered moving average."""
    k = np.ones(window) / window
    pad = window // 2
    return np.convolve(np.pad(y, pad, mode="edge"), k, mode="valid")


def main() -> None:
    from thznerve.plots.style import DOUBLE_COL_MM
    apply_thesis_style(width_mm=DOUBLE_COL_MM, aspect=0.45)

    rows = []
    with CSV_PATH.open() as f:
        r = csv.reader(f)
        next(r)  # header
        for row in r:
            rows.append([float(x) for x in row])
    arr = np.array(rows)
    f_thz = arr[:, 0]
    peak_node = arr[:, 1]
    peak_global = arr[:, 2]
    mean_node = arr[:, 3]

    smooth_node = _smooth(peak_node, window=3)
    peak_idx = int(np.argmax(peak_node))
    f_peak = f_thz[peak_idx]
    e_peak = peak_node[peak_idx]

    fig, ax = plt.subplots()
    ax.plot(f_thz, peak_global, "s--", lw=1.0, ms=2.5, alpha=0.5,
            label="peak |E| (global)")
    ax.plot(f_thz, mean_node, "^:", lw=1.0, ms=2.5, alpha=0.65,
            label="mean |E| in node")
    ax.plot(f_thz, peak_node, "o", lw=0, ms=4.0, label="peak |E| in node")
    ax.plot(f_thz, smooth_node, "-", lw=1.6, alpha=0.85,
            label="peak |E| (3-pt smooth)")

    ax.axvline(0.6, color="grey", ls=":", lw=0.8, alpha=0.7)
    ax.axvline(2.0, color="grey", ls=":", lw=0.8, alpha=0.7)
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel("|E| (normalised, scattered field)")
    ax.set_title("Sim 1 — Frequency sweep (σ_node = 0, baseline geometry)")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0.85, 3.05)
    ax.text(0.62, 2.98, "0.6 THz (exp.)", fontsize=7, color="0.35")
    ax.text(2.02, 2.98, "2 THz (exp.)", fontsize=7, color="0.35")
    ax.annotate(
        f"computed peak: f = {f_peak:.3f} THz\n|E| = {e_peak:.2f}",
        xy=(f_peak, e_peak), xytext=(0.85, 2.85),
        fontsize=7,
        arrowprops=dict(arrowstyle="->", color="0.25", lw=0.7),
    )
    ax.legend(loc="lower right", framealpha=0.92, fontsize=7, ncols=2)
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim1_freq_sweep")
    print(f"wrote {pdf}\nwrote {png}\npeak at f={f_peak:.3f} THz, |E|={e_peak:.3f}")


if __name__ == "__main__":
    main()
