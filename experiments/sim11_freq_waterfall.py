"""Sim 11 — 3D waterfall of |E|(z, f) from the sim 1 frequency sweep.

Loads all 26 HDF5 result files produced by sim 1, reads the axial
profile from each, and renders a 3-D surface with z along one axis,
frequency along the other, and |E| as the height (and colour). Also
shows a top-down heatmap projection in the second panel.

This is the *single most informative figure* from the COMSOL pipeline —
it converts a noisy 1-D spectrum (sim 1) into a 2-D pattern where
resonance features show up as ridges that span a range of z. Reveals
the geometric pattern of the resonant field that a per-frequency peak
sweep would average over.

Pure Python (no COMSOL), runs in < 5 s.
"""

from __future__ import annotations

from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import LightSource

from thznerve.plots.style import DOUBLE_COL_MM, apply_thesis_style, save_figure

REPO_ROOT = Path(__file__).resolve().parent.parent
SIM1_DIR = REPO_ROOT / "results" / "sim1"


def main() -> None:
    apply_thesis_style(width_mm=DOUBLE_COL_MM, aspect=0.5)

    h5_files = sorted(SIM1_DIR.glob("freq_*.h5"))
    if not h5_files:
        raise SystemExit(
            "no sim1 HDF5 files found — run experiments/sim1_freq_sweep.py first"
        )

    freqs = []
    profiles = []
    for p in h5_files:
        with h5py.File(p, "r") as f:
            freqs.append(float(f["scalars"].attrs["frequency_THz"]))
            z = f["axial_profile/z"][...]
            e = f["axial_profile/E_mag"][...]
            profiles.append(e)
    freqs = np.array(freqs)
    z_um = z  # all profiles share the same z grid
    Egrid = np.array(profiles)  # shape (n_freq, n_z)

    print(f"loaded {len(freqs)} freq points × {len(z_um)} z samples")
    print(f"|E| range: [{Egrid.min():.3f}, {Egrid.max():.3f}]")

    F, Z = np.meshgrid(freqs, z_um, indexing="ij")

    fig = plt.figure()
    ax_3d = fig.add_subplot(1, 2, 1, projection="3d")
    ax_2d = fig.add_subplot(1, 2, 2)

    # 3-D surface with shading
    ls = LightSource(azdeg=200, altdeg=30)
    norm = (Egrid - Egrid.min()) / (Egrid.max() - Egrid.min() + 1e-12)
    rgb = ls.shade(norm, cmap=cm.viridis, blend_mode="soft", vert_exag=1.0)
    surf = ax_3d.plot_surface(Z, F, Egrid, facecolors=rgb, rcount=80, ccount=80,
                               linewidth=0, antialiased=True)
    ax_3d.set_xlabel("z (µm)")
    ax_3d.set_ylabel("Frequency (THz)")
    ax_3d.set_zlabel("|E|")
    ax_3d.set_title("|E|(z, f) — full axial profile vs frequency")
    ax_3d.view_init(elev=28, azim=-62)
    ax_3d.set_box_aspect((1.4, 1.1, 0.7))

    # Mark experimental peaks as semi-transparent vertical planes
    for f_peak in (0.6, 2.0):
        z_line = np.linspace(z_um.min(), z_um.max(), 50)
        f_line = np.full_like(z_line, f_peak)
        e_line = np.full_like(z_line, Egrid.max())
        ax_3d.plot(z_line, f_line, e_line, color="red", lw=0.8, alpha=0.55)

    # 2-D heatmap projection
    im = ax_2d.pcolormesh(z_um, freqs, Egrid, cmap="viridis",
                          shading="auto")
    # Mark node region
    z_node_lo, z_node_hi = 100, 140
    ax_2d.axvline(z_node_lo, color="white", lw=0.6, alpha=0.5)
    ax_2d.axvline(z_node_hi, color="white", lw=0.6, alpha=0.5)
    ax_2d.text((z_node_lo + z_node_hi) / 2, freqs.min() + 0.02,
               "node", color="white", ha="center", fontsize=7, alpha=0.85)
    ax_2d.axhline(0.6, color="red", lw=0.7, alpha=0.65, ls="--")
    ax_2d.axhline(2.0, color="red", lw=0.7, alpha=0.65, ls="--")
    ax_2d.set_xlabel("z (µm)")
    ax_2d.set_ylabel("Frequency (THz)")
    ax_2d.set_title("Top-down view: |E|(z, f) heatmap")
    fig.colorbar(im, ax=ax_2d, shrink=0.85, pad=0.02).set_label("|E|")

    fig.suptitle("Sim 11 — Frequency-sweep waterfall of axial |E| profiles")
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim11_freq_waterfall")
    print(f"wrote {pdf}\nwrote {png}")


if __name__ == "__main__":
    main()
