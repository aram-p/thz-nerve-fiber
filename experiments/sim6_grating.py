"""Sim 6 — Diffraction-grating / wire-array analytical model.

Per Hovhannisyan & Makaryan 2024 paper 3, the spinal-cord sample under
applied voltage behaves as a periodic array of conductive nerve fibres
(a "diffraction grating of wires"). This sim implements the textbook
1D-array-of-wires equivalent-sheet model (Tretyakov 2003, "Analytical
Modelling in Applied Electromagnetics", ch. 4) and plots transmittance,
reflectance, and absorbance vs frequency as a function of node
conductivity σ.

Physics
-------
Array of parallel wires (parallel to E-field), period a, wire radius
r_wire, wire conductivity σ. Immersed in water (ε = water Debye), single
sheet. The sheet impedance for a wire grid with E ∥ wires is

    Z_sheet(ω) = R_grid + j X_L(ω),
    X_L(ω)     = (ω μ₀ a / 2π) · ln(a / (2π r_wire))   (Kontorovich/Tretyakov)
    R_grid     = a / (σ · π r_wire²)                    (DC limit, per square)

Normal-incidence transmittance through a single shunt sheet in a
medium of intrinsic impedance η:

    t = 2 Z_sheet / (2 Z_sheet + η)
    r = -η         / (2 Z_sheet + η)
    A = 1 - |t|² - |r|²

We use η_water = η₀ / √ε_water(f), so the dispersion of water is folded
in via the frequency-dependent permittivity from `materials.py`.

What the figure shows
---------------------
* Transmittance dips and absorbance peaks at frequencies where Z_sheet
  becomes comparable to η/2 — these are σ-tunable resonance-like
  features.
* At σ = 0 the sheet is a pure-reactive inductor and absorbance is
  zero; finite σ turns absorbance on.
* As σ grows, absorbance increases up to a critical value, then
  decreases as the wires become near-perfect reflectors.

This is an analytical first-order model; it does NOT include
node-segment finiteness or fibre-radius-dependent capacitive coupling.
For quantitative agreement with the experimental 0.6 / 2 THz peaks,
the full 3D FEM is needed (sim 1-3). This sim's role is to show the
mechanism by which a periodic fibre array produces frequency-dependent
THz absorption.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from thznerve.model.materials import debye_water_epsilon
from thznerve.plots.style import DOUBLE_COL_MM, apply_thesis_style, save_figure

MU_0: float = 4e-7 * np.pi          # H/m
EPS_0: float = 8.8541878128e-12     # F/m
C_0: float = 1 / np.sqrt(MU_0 * EPS_0)
ETA_0: float = np.sqrt(MU_0 / EPS_0)  # ~376.7 Ω


def wire_grid_sheet_impedance(
    f_hz: np.ndarray, *, period_um: float, wire_radius_um: float, sigma_S_per_m: float
) -> np.ndarray:
    """Z_sheet(ω) = R_grid + j X_L(ω) for a 1D wire array, E ∥ wires."""

    omega = 2 * np.pi * f_hz
    a = period_um * 1e-6
    r = wire_radius_um * 1e-6

    X_L = (omega * MU_0 * a / (2 * np.pi)) * np.log(a / (2 * np.pi * r))
    if sigma_S_per_m > 0:
        # DC-limit per-square resistance of the wire array.
        R_grid = a / (sigma_S_per_m * np.pi * r**2)
    else:
        R_grid = np.inf * np.ones_like(f_hz)  # effectively no current → no abs
    Z = R_grid + 1j * X_L
    return Z


def single_sheet_TR(Z_sheet: np.ndarray, eta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Voltage transmission/reflection of a shunt sheet between two half-spaces of
    intrinsic impedance ``eta`` (assumed equal on both sides).
    """

    denom = 2 * Z_sheet + eta
    t = 2 * Z_sheet / denom
    r = -eta / denom
    return t, r


def main() -> None:
    apply_thesis_style(width_mm=DOUBLE_COL_MM, aspect=0.5)

    f_thz = np.linspace(0.05, 2.5, 600)
    f_hz = f_thz * 1e12

    # Water-immersed sheet — intrinsic impedance of water at each frequency.
    eps_water = debye_water_epsilon(f_hz)
    eta_w = ETA_0 / np.sqrt(eps_water)  # complex impedance of water

    # The wire-array reactance scales like ωln(a/r) and the resistance like
    # a/(σ π r²). To see a clear resonance-like absorbance peak we need
    # σ values where the resistive part matches the reactance somewhere in
    # the band — that happens at much higher conductivities (the node-channel
    # interpretation in paper 3 imagines highly conductive segments) and
    # closer wire packing.
    period_um = 50.0        # closer fibre packing (matches dense white matter)
    wire_radius_um = 2.5    # effective node-channel radius (smaller than axon)

    sigmas = [10.0, 1e2, 1e3, 1e4, 1e5]  # S/m, log-spaced — node-channel regime

    fig, (ax_T, ax_A) = plt.subplots(1, 2, sharex=True)

    for sigma in sigmas:
        Z = wire_grid_sheet_impedance(
            f_hz, period_um=period_um, wire_radius_um=wire_radius_um,
            sigma_S_per_m=sigma,
        )
        t, r = single_sheet_TR(Z, eta_w)
        T = np.abs(t) ** 2
        R = np.abs(r) ** 2
        A = 1.0 - T - R
        A = np.clip(A, 0.0, 1.0)

        label = fr"$\sigma = {sigma:g}$ S/m"
        ax_T.plot(f_thz, T, label=label, lw=1.4)
        ax_A.plot(f_thz, A, label=label, lw=1.4)

    for ax in (ax_T, ax_A):
        ax.axvline(0.6, color="grey", ls=":", lw=0.8, alpha=0.7)
        ax.axvline(2.0, color="grey", ls=":", lw=0.8, alpha=0.7)
        ax.set_xlabel("Frequency (THz)")
        ax.grid(True, alpha=0.3)

    ax_T.set_ylabel("Transmittance |T|²")
    ax_T.set_title("Transmittance through wire-grid sheet (water host)")
    ax_T.legend(loc="lower right", framealpha=0.95)

    ax_A.set_ylabel("Absorbance A = 1 − |T|² − |R|²")
    ax_A.set_title("Absorbance — σ-tuned via node conductivity")
    ax_A.legend(loc="upper right", framealpha=0.95)

    # Annotate experimentally observed peaks
    ax_A.annotate(
        "0.6 THz (exp.)", xy=(0.6, 0.92), xytext=(0.7, 0.94),
        fontsize=7, color="grey",
    )
    ax_A.annotate(
        "2 THz (exp.)", xy=(2.0, 0.92), xytext=(1.65, 0.94),
        fontsize=7, color="grey",
    )

    fig.suptitle(
        f"Sim 6 — Wire-array model: period {period_um:g} µm, "
        f"wire r = {wire_radius_um:g} µm, water host"
    )
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim6_grating")
    print(f"wrote {pdf}\nwrote {png}")


if __name__ == "__main__":
    main()
