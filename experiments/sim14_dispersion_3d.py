"""Sim 14 — 3D complex-ε dispersion curves.

Plots the three nerve-fibre materials as curves in (Re ε, Im ε, freq)
space: water (double-Debye), myelin (constant 4.5 − 0.5 i), and node at
five σ values from 0 to 1000 S/m. The frequency axis is vertical; each
material traces out a path through the complex-ε plane as f varies.

Why this is interesting to a tutor
----------------------------------
A 2-panel (Re, Im) plot like sim 4 makes you flip your eyes back and
forth to relate frequency to position; a single 3-D plot puts both
ε-components and frequency in one view. The geometry of the
double-Debye relaxation (a curved trajectory in ε-space) becomes
visible at a glance, as does the way the σ term lifts the node curve
along Im ε at low frequency only.

Pure Python, no COMSOL.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from thznerve.model.materials import (
    MYELIN_EPS_IMAG, MYELIN_EPS_REAL,
    debye_water_epsilon, node_epsilon,
)
from thznerve.plots.style import DOUBLE_COL_MM, OKABE_ITO, apply_thesis_style, save_figure


def main() -> None:
    apply_thesis_style(width_mm=DOUBLE_COL_MM, aspect=0.62)

    f_thz = np.linspace(0.05, 2.5, 250)
    f_hz = f_thz * 1e12

    eps_water = debye_water_epsilon(f_hz)
    eps_myelin_const = MYELIN_EPS_REAL + 1j * MYELIN_EPS_IMAG

    sigma_values = [0.0, 1.0, 10.0, 100.0, 1000.0]
    sigma_colors = [OKABE_ITO[0], OKABE_ITO[3], OKABE_ITO[5], OKABE_ITO[6], OKABE_ITO[7]]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # Water curve
    ax.plot(eps_water.real, -eps_water.imag, f_thz,
            color=OKABE_ITO[2], lw=2.0, label="water (double Debye)")
    # Myelin — vertical line at the constant ε
    ax.plot([MYELIN_EPS_REAL] * len(f_thz),
            [-MYELIN_EPS_IMAG] * len(f_thz),
            f_thz,
            color=OKABE_ITO[1], lw=1.6, ls="--", label="myelin (constant)")

    # Node curves at different σ values
    for sigma, c in zip(sigma_values, sigma_colors):
        eps_n = node_epsilon(f_hz, sigma_S_per_m=sigma)
        ax.plot(eps_n.real, -eps_n.imag, f_thz,
                color=c, lw=1.4, alpha=0.85,
                label=fr"node, $\sigma = {sigma:g}$ S/m")

    # Vertical drop-lines at experimentally observed frequencies
    for f_peak, label in [(0.6, "0.6 THz exp."), (2.0, "2 THz exp.")]:
        ax.plot([0, 80], [0, 0], [f_peak, f_peak], color="0.7", lw=0.6, alpha=0.6)
        ax.text(80, 0.1, f_peak, label, fontsize=6.5, color="0.4")

    # Highlight ε at f=0.6 THz for each material
    f_target = 0.6e12
    eps_w_06 = complex(debye_water_epsilon(f_target))
    ax.scatter([eps_w_06.real], [-eps_w_06.imag], [0.6],
               color=OKABE_ITO[2], s=40, zorder=10, edgecolors="white", linewidths=0.6)
    ax.scatter([MYELIN_EPS_REAL], [-MYELIN_EPS_IMAG], [0.6],
               color=OKABE_ITO[1], s=40, zorder=10, edgecolors="white", linewidths=0.6)

    ax.set_xlabel(r"Re $\varepsilon_r$")
    ax.set_ylabel(r"$-$Im $\varepsilon_r$ (loss)")
    ax.set_zlabel("Frequency (THz)")
    ax.set_title("Sim 14 — Complex permittivity ε(f) of the three materials, 3-D view")
    ax.legend(loc="upper left", fontsize=6.5, framealpha=0.92)
    ax.view_init(elev=22, azim=-58)

    fig.tight_layout()
    pdf, png = save_figure(fig, "sim14_dispersion_3d")
    print(f"wrote {pdf}\nwrote {png}")


if __name__ == "__main__":
    main()
