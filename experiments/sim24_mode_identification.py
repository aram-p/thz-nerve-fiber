"""Sim 24 — Mode identification: the two resonance peaks are half-wave
standing modes on different fibre segments.

The two peaks resolved by Sim 22 are:

  * Low window:  f₀ = 0.619 ± 0.008 THz, Q = 15.0  (annular sampling)
  * High window: f₀ = 1.938 ± 0.010 THz, Q = 30.5  (axial sampling)

This sim asks: do these correspond to half-wave standing modes on the
fibre geometry? Specifically:

  * λ_water(f) / 2 ≈ internode + half-node length  → predicts the low peak
  * λ_water(f) / 2 ≈ node length                   → predicts the high peak

Both predictions use the double-Debye water refractive index n(f) at the
relevant f (the wave inside the sheath is in a water-like medium).

Outputs two 1-D plots:
  * Predicted vs observed peak f as a function of segment length
  * λ_water(f) and λ_water/2 overlaid with the geometric segment lengths

Pure Python, no COMSOL.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from thznerve.model.materials import debye_water_epsilon
from thznerve.plots.style import DOUBLE_COL_MM, OKABE_ITO, apply_thesis_style, save_figure

REPO_ROOT = Path(__file__).resolve().parent.parent
C_0 = 2.99792458e8  # m/s

# Observed peaks from Sim 22 (refined mesh, dense sampling).
PEAKS_OBSERVED = [
    {"label": "low",  "f_THz": 0.619, "Q": 15.0, "sampling": "annular"},
    {"label": "high", "f_THz": 1.938, "Q": 30.5, "sampling": "axial"},
]

# Geometric segments (baseline geometry, all in µm).
NODE_LENGTH_UM = 40.0
INTERNODE_LENGTH_UM = 100.0
NODE_PLUS_HALF_INTERNODE_UM = NODE_LENGTH_UM + INTERNODE_LENGTH_UM  # = 140
TOTAL_FIBRE_LENGTH_UM = 2 * INTERNODE_LENGTH_UM + NODE_LENGTH_UM    # = 240


def n_water(f_hz):
    """Real part of the water refractive index from double-Debye permittivity."""
    return np.sqrt(debye_water_epsilon(f_hz)).real


def lambda_water_um(f_hz):
    """Wavelength in water at frequency f, in µm."""
    return (C_0 / (n_water(f_hz) * f_hz)) * 1e6


def predicted_half_wave_freq_THz(segment_length_um, n_guess=2.0):
    """Iteratively find f such that λ_water(f) / 2 = segment_length_um.

    n_water depends on f (Debye dispersion), so we iterate fixed-point.
    """
    f_hz = C_0 / (n_guess * 2 * segment_length_um * 1e-6)
    for _ in range(50):
        n = n_water(f_hz)
        f_new = C_0 / (n * 2 * segment_length_um * 1e-6)
        if abs(f_new - f_hz) / f_hz < 1e-8:
            break
        f_hz = f_new
    return f_hz / 1e12


def main() -> None:
    apply_thesis_style(width_mm=DOUBLE_COL_MM, aspect=0.55)

    # ---- Predicted half-wave frequencies for each segment ----
    predictions = {
        "node only ({:g} µm)".format(NODE_LENGTH_UM):
            predicted_half_wave_freq_THz(NODE_LENGTH_UM),
        "internode + half-node ({:g} µm)".format(0.5 * NODE_LENGTH_UM + INTERNODE_LENGTH_UM):
            predicted_half_wave_freq_THz(0.5 * NODE_LENGTH_UM + INTERNODE_LENGTH_UM),
        "node + 2 half-internode ({:g} µm)".format(NODE_LENGTH_UM + INTERNODE_LENGTH_UM):
            predicted_half_wave_freq_THz(NODE_LENGTH_UM + INTERNODE_LENGTH_UM),
        "full fibre ({:g} µm)".format(TOTAL_FIBRE_LENGTH_UM):
            predicted_half_wave_freq_THz(TOTAL_FIBRE_LENGTH_UM),
    }

    print("=== Half-wave standing-mode predictions ===")
    for label, f_pred in predictions.items():
        print(f"  λ/2 = {label:35s} → f = {f_pred:.3f} THz")
    print()
    print("=== Observed peaks ===")
    for p in PEAKS_OBSERVED:
        print(f"  {p['label']} ({p['sampling']:>7}): f = {p['f_THz']:.3f} THz, Q = {p['Q']:.1f}")

    # ---- Figure ----
    fig, (ax_wl, ax_match) = plt.subplots(1, 2)

    # Left: λ_water(f) and λ_water/2 vs frequency, with horizontal lines at
    # the relevant segment lengths.
    f_grid_THz = np.linspace(0.1, 3.0, 400)
    f_grid_Hz = f_grid_THz * 1e12
    lam_w = lambda_water_um(f_grid_Hz)

    ax_wl.plot(f_grid_THz, lam_w, color=OKABE_ITO[2], lw=1.6,
               label=r"$\lambda_w(f)$ — wavelength in water")
    ax_wl.plot(f_grid_THz, lam_w / 2, "--", color=OKABE_ITO[2], lw=1.3,
               alpha=0.7, label=r"$\lambda_w(f) / 2$ — half-wavelength")

    segments = [
        ("node (40 µm)", NODE_LENGTH_UM, OKABE_ITO[6]),
        ("internode + half-node (120 µm)", 0.5 * NODE_LENGTH_UM + INTERNODE_LENGTH_UM,
         OKABE_ITO[5]),
        ("node + 2×half-internode (140 µm)", NODE_LENGTH_UM + INTERNODE_LENGTH_UM,
         OKABE_ITO[1]),
    ]
    for label, L, c in segments:
        ax_wl.axhline(L, color=c, ls=":", lw=1.1, alpha=0.85)
        ax_wl.text(2.85, L * 1.05, label, fontsize=6.5, color=c, ha="right")

    # Mark the observed peaks
    for p in PEAKS_OBSERVED:
        ax_wl.axvline(p["f_THz"], color="red", ls=":", lw=0.8, alpha=0.7)
        ax_wl.text(p["f_THz"] + 0.02, ax_wl.get_ylim()[1] * 0.05,
                   f"observed {p['f_THz']:.3f} THz",
                   rotation=90, fontsize=6.5, color="0.3")

    ax_wl.set_xlabel("Frequency (THz)")
    ax_wl.set_ylabel("Length (µm)")
    ax_wl.set_yscale("log")
    ax_wl.set_xlim(0.1, 3.0)
    ax_wl.set_ylim(20, 4000)
    ax_wl.set_title(r"Wavelength in water vs fibre segment scales")
    ax_wl.legend(loc="upper right", fontsize=6.5, framealpha=0.92)
    ax_wl.grid(True, which="both", alpha=0.25)

    # Right: predicted vs observed peak frequency
    f_peak_low = PEAKS_OBSERVED[0]["f_THz"]
    f_peak_high = PEAKS_OBSERVED[1]["f_THz"]
    seg_lengths = np.array([s[1] for s in segments])
    seg_labels = [s[0] for s in segments]

    pred_freqs = [predictions[k] for k in [
        "node only (40 µm)",
        "internode + half-node (120 µm)",
        "node + 2 half-internode (140 µm)",
    ]]

    ax_match.scatter(seg_lengths, pred_freqs, s=90, color=OKABE_ITO[2],
                     marker="o", label="predicted (λ_w / 2 = L)", zorder=5,
                     edgecolor="white", linewidths=1.2)
    for L, p, lbl in zip(seg_lengths, pred_freqs, seg_labels):
        ax_match.text(L * 1.04, p, lbl.split(" (")[0], fontsize=6.5)

    # Observed peaks
    ax_match.axhline(f_peak_low, color=OKABE_ITO[6], ls="--", lw=1.0, alpha=0.7,
                     label=f"observed low: {f_peak_low:.3f} THz")
    ax_match.axhline(f_peak_high, color=OKABE_ITO[1], ls="--", lw=1.0, alpha=0.7,
                     label=f"observed high: {f_peak_high:.3f} THz")
    ax_match.scatter([NODE_LENGTH_UM], [f_peak_high], s=140, marker="*",
                     color=OKABE_ITO[1], edgecolor="white", linewidths=1.4,
                     label="observed @ node length", zorder=6)
    ax_match.scatter([NODE_LENGTH_UM + INTERNODE_LENGTH_UM], [f_peak_low], s=140,
                     marker="*", color=OKABE_ITO[6], edgecolor="white", linewidths=1.4,
                     label="observed @ node + 2 half-internode", zorder=6)

    ax_match.set_xlabel("Segment length L (µm)")
    ax_match.set_ylabel("Frequency (THz)")
    ax_match.set_title("Predicted half-wave f vs observed peaks")
    ax_match.set_xscale("log")
    ax_match.set_yscale("log")
    ax_match.set_xlim(20, 400)
    ax_match.legend(loc="lower left", fontsize=6.5, framealpha=0.92)
    ax_match.grid(True, which="both", alpha=0.25)

    fig.suptitle("Sim 24 — The two peaks are half-wave standing modes on different fibre segments")
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim24_mode_identification")
    print(f"\nwrote {pdf}\nwrote {png}")

    # Quantitative agreement table
    print()
    print("=== Match quality ===")
    print(f"  Low peak observed: {f_peak_low:.3f} THz")
    pred_for_low = predictions["node + 2 half-internode (140 µm)"]
    err_low = (pred_for_low - f_peak_low) / f_peak_low * 100
    print(f"  Predicted from L = 140 µm (node + 2 half-internodes): {pred_for_low:.3f} THz  "
          f"({err_low:+.1f}%)")

    print(f"  High peak observed: {f_peak_high:.3f} THz")
    pred_for_high = predictions["node only (40 µm)"]
    err_high = (pred_for_high - f_peak_high) / f_peak_high * 100
    print(f"  Predicted from L = 40 µm (node only):                  {pred_for_high:.3f} THz  "
          f"({err_high:+.1f}%)")


if __name__ == "__main__":
    main()
