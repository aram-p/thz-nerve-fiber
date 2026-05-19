"""Sim 7 — Wavelength vs geometry-scale diagram.

Plots vacuum wavelength λ₀ = c/f and the wavelength in water
λ_w = c/(n_w(f) f) over 0.1–3 THz, overlaid with the geometric scales
of the modelled fibre (axon, myelin, node, internode, total length,
unit cell). The figure makes immediately visible:

  * The wavelength in water at THz is comparable to the node length and
    internode length — that's the regime where dipole / half-wave
    resonances on the fibre geometry are expected.
  * The vacuum wavelength is ~2× larger than the water wavelength
    because Re(n_w) ≈ 2 across this band.

Pure Python, no COMSOL. Runs in < 1 s.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from thznerve.model.materials import debye_water_epsilon
from thznerve.plots.style import SINGLE_COL_MM, OKABE_ITO, apply_thesis_style, save_figure

C_0: float = 2.99792458e8  # m/s


def main() -> None:
    apply_thesis_style(width_mm=SINGLE_COL_MM, aspect=0.85)

    f_thz = np.linspace(0.1, 3.0, 600)
    f_hz = f_thz * 1e12

    lam_vac_um = (C_0 / f_hz) * 1e6  # µm
    eps_w = debye_water_epsilon(f_hz)
    n_w = np.sqrt(eps_w).real  # real part of refractive index
    lam_water_um = lam_vac_um / n_w

    geom_scales = [
        ("axon diameter", 10.0),
        ("myelin diameter", 14.0),
        ("node length", 40.0),
        ("internode length", 100.0),
        ("unit-cell side", 40.0),   # 2 * external_half_width_um
        ("fibre total length", 240.0),
    ]

    fig, ax = plt.subplots()

    ax.plot(f_thz, lam_vac_um, lw=1.6, color=OKABE_ITO[0],
            label=r"$\lambda_0 = c/f$ (vacuum)")
    ax.plot(f_thz, lam_water_um, lw=1.6, color=OKABE_ITO[2],
            label=r"$\lambda_w = c/(n_w f)$ (water)")
    ax.plot(f_thz, lam_water_um / 2, ":", lw=1.0, color=OKABE_ITO[2], alpha=0.6,
            label=r"$\lambda_w/2$ (half-wave)")

    # Geometry scale lines
    palette = [OKABE_ITO[1], OKABE_ITO[3], OKABE_ITO[5], OKABE_ITO[6], OKABE_ITO[7], "0.4"]
    for (name, val), c in zip(geom_scales, palette):
        ax.axhline(val, lw=0.7, ls="--", color=c, alpha=0.75)
        ax.text(2.95, val * 1.04, name, fontsize=6.5, ha="right", color=c)

    # Experimental peaks
    for f_peak in (0.6, 2.0):
        ax.axvline(f_peak, lw=0.7, color="grey", alpha=0.6)
    ax.text(0.61, 1.2, "0.6 THz", color="0.3", fontsize=7, rotation=90, va="bottom")
    ax.text(2.01, 1.2, "2 THz",  color="0.3", fontsize=7, rotation=90, va="bottom")

    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel("Length (µm)")
    ax.set_yscale("log")
    ax.set_xlim(0.1, 3.0)
    ax.set_ylim(1, 4000)
    ax.set_title("Sim 7 — Wavelength scales vs fibre geometry")
    ax.legend(loc="upper right", framealpha=0.95)
    ax.grid(True, which="both", alpha=0.25)

    fig.tight_layout()
    pdf, png = save_figure(fig, "sim7_wavelength_scales")
    print(f"wrote {pdf}\nwrote {png}")


if __name__ == "__main__":
    main()
