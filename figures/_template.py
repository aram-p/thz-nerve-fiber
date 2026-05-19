"""Template figure script — copy this as the starting point for new figures.

Loads the first row of `results/manifest.csv`, opens its HDF5, plots the
axial |E| profile. Falls back to a synthetic placeholder until Phase 2
produces real simulation results.

Run:
    uv run python figures/_template.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from thznerve.io.hdf5 import read_result
from thznerve.plots.style import apply_thesis_style, save_figure

REPO_ROOT = Path(__file__).resolve().parent.parent
MANIFEST = REPO_ROOT / "results" / "manifest.csv"


def _placeholder() -> tuple[np.ndarray, np.ndarray]:
    z = np.linspace(0.0, 240.0, 200)
    e_mag = 1.0 + 0.4 * np.exp(-((z - 120.0) ** 2) / (2.0 * 15.0**2))
    return z, e_mag


def _load_first_axial_profile() -> tuple[tuple[np.ndarray, np.ndarray], bool]:
    """Return ((z, |E|), is_placeholder)."""

    if not MANIFEST.exists():
        return _placeholder(), True
    df = pd.read_csv(MANIFEST, dtype=str, keep_default_na=False)
    if df.empty:
        return _placeholder(), True
    h5_rel = df.iloc[0]["output_path"]
    h5_path = (REPO_ROOT / h5_rel).resolve()
    if not h5_path.exists():
        return _placeholder(), True
    result = read_result(h5_path)
    return (result.axial_profile.z, result.axial_profile.e_mag), False


def main() -> None:
    apply_thesis_style()
    (z, e_mag), is_placeholder = _load_first_axial_profile()

    fig, ax = plt.subplots()
    ax.plot(z, e_mag)
    ax.set_xlabel("z (µm)")
    ax.set_ylabel("|E| (normalised)")
    title = "Axial |E| profile"
    if is_placeholder:
        title += "  (placeholder — Phase 2 not yet run)"
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    pdf, png = save_figure(fig, "_template")
    print(f"wrote {pdf}")
    print(f"wrote {png}")


if __name__ == "__main__":
    main()
