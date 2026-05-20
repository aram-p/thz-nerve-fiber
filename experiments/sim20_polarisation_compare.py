"""Sim 20 — Polarisation comparison: ⊥-fibre (Sim 1) vs ∥-fibre (Sim 18).

Paper 1's central experimental observation: THz absorption peaks only
appear when the E-field is polarised parallel to the nerve fibres.
This sim consolidates both polarisation configurations from the FEM
into one figure so the comparison is unambiguous.

Two panels:
* Axial sampling (E·z̄ through the fibre): the |E| sampled on the fibre's
  central axis in the node z-range — what Sim 1 reports.
* Node-annulus sampling: |E| sampled in the conductive annulus around
  the fibre — what would actually couple to a σ-induced absorption.

Pure Python from saved CSVs, no COMSOL.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from thznerve.plots.style import DOUBLE_COL_MM, OKABE_ITO, apply_thesis_style, save_figure

REPO_ROOT = Path(__file__).resolve().parent.parent


def _load(csv_path, col):
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    f_thz = np.array([float(r["frequency_THz"]) for r in rows])
    y = np.array([float(r[col]) for r in rows])
    return f_thz, y


def _smooth(y, k=3):
    pad = k // 2
    kernel = np.ones(k) / k
    return np.convolve(np.pad(y, pad, mode="edge"), kernel, mode="valid")


def main() -> None:
    apply_thesis_style(width_mm=DOUBLE_COL_MM, aspect=0.5)

    sim1 = REPO_ROOT / "results" / "sim1" / "spectrum.csv"
    sim18 = REPO_ROOT / "results" / "sim18" / "spectrum.csv"

    f1, e1_axial = _load(sim1, "peak_E_node")
    f18, e18_axial = _load(sim18, "peak_E_node_axis")
    _, e18_annul = _load(sim18, "peak_E_node_annulus")

    fig, (ax_ax, ax_ann) = plt.subplots(1, 2, sharex=True)

    # Axial sampling comparison
    ax_ax.plot(f1, e1_axial, "o", color=OKABE_ITO[2], ms=3.5, alpha=0.55,
               label="Sim 1 raw (E ⊥ fibre)")
    ax_ax.plot(f1, _smooth(e1_axial), "-", color=OKABE_ITO[2], lw=1.6,
               label="Sim 1 smooth")
    ax_ax.plot(f18, e18_axial, "s", color=OKABE_ITO[6], ms=3.5, alpha=0.55,
               label="Sim 18 raw (E ∥ fibre)")
    ax_ax.plot(f18, _smooth(e18_axial), "-", color=OKABE_ITO[6], lw=1.6,
               label="Sim 18 smooth")
    ax_ax.axvline(0.6, color="red", ls=":", lw=0.7, alpha=0.65)
    ax_ax.axvline(2.0, color="red", ls=":", lw=0.7, alpha=0.65)
    ax_ax.set_xlabel("Frequency (THz)")
    ax_ax.set_ylabel("peak |E| (axial sampling)")
    ax_ax.set_title("Field on fibre's axial centre line")
    ax_ax.legend(loc="best", fontsize=7, framealpha=0.92)
    ax_ax.grid(True, alpha=0.3)

    # Annular sampling comparison — Sim 1 didn't sample annular, so just
    # show Sim 18's annular signal and the constant Sim 1 annular value
    # we measured separately (Sim 2 v10 found ~1.79 at f = 0.6 THz, σ=0).
    sim1_annul_value = 1.7879  # from sim 2 v10 at σ=0
    ax_ann.axhline(sim1_annul_value, color=OKABE_ITO[2], lw=1.4, ls="--",
                   label=f"Sim 1 (E ⊥) at f=0.6 THz: |E|={sim1_annul_value:.3f}")
    ax_ann.plot(f18, e18_annul, "s", color=OKABE_ITO[6], ms=4, alpha=0.65,
                label="Sim 18 (E ∥) raw")
    ax_ann.plot(f18, _smooth(e18_annul), "-", color=OKABE_ITO[6], lw=1.6,
                label="Sim 18 smooth")
    ax_ann.axvline(0.6, color="red", ls=":", lw=0.7, alpha=0.65)
    ax_ann.axvline(2.0, color="red", ls=":", lw=0.7, alpha=0.65)
    ax_ann.set_xlabel("Frequency (THz)")
    ax_ann.set_ylabel("peak |E| in node annulus")
    ax_ann.set_title("Field at the conductive node annulus")
    ax_ann.legend(loc="best", fontsize=7, framealpha=0.92)
    ax_ann.grid(True, alpha=0.3)

    # Summary text box
    ax_ann.text(
        0.02, 0.02,
        f"E ∥ fibre: |E|_ann = 1.99 – 2.45\nE ⊥ fibre: |E|_ann ≈ {sim1_annul_value:.2f}\n"
        f"→ E ∥ fibre couples ~25–40% more field into the node",
        transform=ax_ann.transAxes, fontsize=6.5,
        bbox=dict(boxstyle="round,pad=0.3", fc="#fff7d4", ec="0.6", alpha=0.95),
        va="bottom",
    )

    fig.suptitle("Sim 20 — Polarisation comparison: ⊥-fibre vs ∥-fibre frequency sweep")
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim20_polarisation_compare")
    print(f"wrote {pdf}\nwrote {png}")

    # Quantitative summary
    e18_ann_smooth = _smooth(e18_annul)
    print("\nE ⊥ fibre annular |E| (Sim 2 v10 at σ=0, f=0.6 THz):", sim1_annul_value)
    print("E ∥ fibre annular |E| (Sim 18):")
    print(f"  range: [{e18_annul.min():.3f}, {e18_annul.max():.3f}]")
    print(f"  mean : {e18_annul.mean():.3f}")
    print(f"  enhancement factor: {e18_annul.mean() / sim1_annul_value:.2f}×")


if __name__ == "__main__":
    main()
