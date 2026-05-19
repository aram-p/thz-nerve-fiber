"""Sim 4 — Material dispersion in the THz range.

Plots Re(ε) and Im(ε) of the three materials used in the model:

  * Water (double-Debye, Ellison/Jepsen parameters)
  * Myelin (constant 4.5 − 0.5 i at THz)
  * Node of Ranvier = water + i σ / (ω ε₀), at three σ values (0, 1, 10 S/m)

Why this matters for the tutor presentation
-------------------------------------------
The choice of ε(f) drives every other result. Showing the actual
dispersion curves makes the modelling choices explicit:

  * Water has a real part dropping from ~5 to ~3.5 across 0.1–2 THz,
    and a non-monotonic imaginary part with a broad loss peak — these
    are the well-known double-Debye relaxations of liquid water.
  * Myelin is intentionally constant — the (small) loss term is the
    only frequency dependence we assert at THz.
  * The node's σ-term shows up as a 1/f imaginary contribution on
    top of water, which is large at low frequency and dies out by 2 THz.

No COMSOL needed; this runs in < 1 s.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from thznerve.model.materials import (
    debye_water_epsilon,
    myelin_epsilon,
    node_epsilon,
)
from thznerve.plots.style import DOUBLE_COL_MM, apply_thesis_style, save_figure


def main() -> None:
    apply_thesis_style(width_mm=DOUBLE_COL_MM, aspect=0.45)

    f_thz = np.linspace(0.05, 2.5, 400)
    f_hz = f_thz * 1e12

    eps_water = debye_water_epsilon(f_hz)
    eps_myelin = np.full_like(f_hz, myelin_epsilon(0.6e12), dtype=complex)
    eps_node_0 = node_epsilon(f_hz, sigma_S_per_m=0.0)
    eps_node_1 = node_epsilon(f_hz, sigma_S_per_m=1.0)
    eps_node_10 = node_epsilon(f_hz, sigma_S_per_m=10.0)

    fig, (ax_re, ax_im) = plt.subplots(1, 2)

    ax_re.plot(f_thz, eps_water.real, label="water (double Debye)", lw=1.6)
    ax_re.plot(f_thz, eps_myelin.real, "--", label="myelin", lw=1.4)
    ax_re.plot(f_thz, eps_node_0.real, ":", label="node, σ=0", lw=1.4)
    ax_re.plot(f_thz, eps_node_1.real, "-.", label="node, σ=1 S/m", lw=1.2)
    ax_re.plot(f_thz, eps_node_10.real, "-.", label="node, σ=10 S/m", lw=1.0)
    ax_re.set_xlabel("Frequency (THz)")
    ax_re.set_ylabel(r"Re$\{\varepsilon_r(f)\}$")
    ax_re.set_title("Real part")
    ax_re.legend(loc="upper right", framealpha=0.95)
    ax_re.grid(True, alpha=0.3)

    ax_im.plot(f_thz, -eps_water.imag, label="water", lw=1.6)
    ax_im.plot(f_thz, -eps_myelin.imag, "--", label="myelin", lw=1.4)
    ax_im.plot(f_thz, -eps_node_0.imag, ":", label="node, σ=0", lw=1.4)
    ax_im.plot(f_thz, -eps_node_1.imag, "-.", label="node, σ=1 S/m", lw=1.2)
    ax_im.plot(f_thz, -eps_node_10.imag, "-.", label="node, σ=10 S/m", lw=1.0)
    ax_im.set_xlabel("Frequency (THz)")
    ax_im.set_ylabel(r"$-$Im$\{\varepsilon_r(f)\}$ (loss)")
    ax_im.set_title("Imaginary part (loss)")
    ax_im.set_yscale("log")
    ax_im.legend(loc="upper right", framealpha=0.95)
    ax_im.grid(True, which="both", alpha=0.3)

    fig.suptitle("Sim 4 — Material dispersion of nerve-fibre model")
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim4_dispersion")
    print(f"wrote {pdf}\nwrote {png}")


if __name__ == "__main__":
    main()
