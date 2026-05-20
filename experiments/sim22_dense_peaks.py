"""Sim 22 — Dense frequency sweeps with refined mesh near the two peaks.

Sim 1 (26 points, baseline mesh) found a peak around 0.6 THz but
nothing decisive at 2 THz. Two candidate explanations: (a) sampling
density too low to resolve a narrow peak there, (b) mesh too coarse at
2 THz where λ_water/2 ≈ 75 µm vs our baseline 30 µm mesh.

This sim does dense sweeps in two windows with **finer mesh**:
  * 0.55 – 0.70 THz, 10 points
  * 1.85 – 2.15 THz, 10 points
  * MeshParams(max_element_size_um=15, node_element_size_um=3, refinement_factor=1.0)

Outputs both spectra + a Lorentzian fit on each window. The expectation
is that the 0.6 THz peak shrinks slightly (mesh-refined value ~2.48
per Sim 16) and the 2 THz peak — if real — finally resolves.
"""

from __future__ import annotations

import csv
import time
from pathlib import Path

import matplotlib.pyplot as plt
import mph
import numpy as np
from scipy.optimize import curve_fit

from thznerve.io.provenance import get_git_sha
from thznerve.model.geometry import GeometryParams, build_geometry, total_length_um
from thznerve.model.materials import MaterialParams, apply_materials
from thznerve.model.mesh import MeshParams, build_mesh
from thznerve.model.study import (
    _evaluate_at_points, setup_physics, setup_study, solve_study,
)
from thznerve.plots.style import DOUBLE_COL_MM, OKABE_ITO, apply_thesis_style, save_figure

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "results" / "sim22"

GEOM = GeometryParams(
    axon_radius_um=5, myelin_radius_um=7,
    node_length_um=40, internode_length_um=100, external_half_width_um=20,
)
MAT = MaterialParams(node_sigma_S_per_m=0.0)
# Refined mesh: 15 µm max, 3 µm at node — captures λ_water/2 = 75 µm at 2 THz.
MESH = MeshParams(max_element_size_um=15.0, node_element_size_um=3.0)

WINDOW_LOW = np.linspace(0.55, 0.70, 10)   # near 0.6 THz peak
WINDOW_HIGH = np.linspace(1.85, 2.15, 10)  # near 2 THz peak
FREQ_THZ = np.concatenate([WINDOW_LOW, WINDOW_HIGH])


def lorentzian(f, A, f0, gamma, b, m):
    return A * (gamma / 2) ** 2 / ((f - f0) ** 2 + (gamma / 2) ** 2) + b + m * f


def fit_peak(f, y, *, f0_guess: float, window: float = 0.1):
    mask = (f >= f0_guess - window) & (f <= f0_guess + window)
    f_fit, y_fit = f[mask], y[mask]
    if len(f_fit) < 4:
        return None, None
    A0 = float(np.max(y_fit) - np.min(y_fit))
    p0 = [max(A0, 0.05), f0_guess, 0.07, float(np.min(y_fit)), 0.0]
    bounds = (
        [1e-4,        f0_guess - 0.15, 0.01, -10.0, -50.0],
        [10 * (A0 + 0.5), f0_guess + 0.15, window, 10.0, 50.0],
    )
    try:
        popt, pcov = curve_fit(lorentzian, f_fit, y_fit, p0=p0, bounds=bounds, maxfev=10000)
        perr = np.sqrt(np.diag(pcov))
    except Exception:
        return None, None
    return popt, perr


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    apply_thesis_style(width_mm=DOUBLE_COL_MM, aspect=0.5)
    t0 = time.monotonic()

    client = mph.start()
    model = client.create("sim22_dense")
    build_geometry(model, GEOM)
    apply_materials(model, MAT, freq_hz=FREQ_THZ[0] * 1e12)
    setup_physics(model, GEOM)
    n_elem = build_mesh(model, MESH)
    setup_study(model, FREQ_THZ[0] * 1e12)
    print(f"setup: {n_elem} elements, {time.monotonic() - t0:.1f}s")

    git_sha = get_git_sha()
    rows = []

    L = total_length_um(GEOM)
    z_node_lo = GEOM.internode_length_um
    z_node_hi = z_node_lo + GEOM.node_length_um
    r_node = 0.5 * (GEOM.axon_radius_um + GEOM.myelin_radius_um)

    # Annulus sampling points
    z_grid = np.linspace(z_node_lo + 1, z_node_hi - 1, 24)
    angles = np.linspace(0, 2 * np.pi, 6, endpoint=False)
    ann_pts = np.array([(r_node * np.cos(a), r_node * np.sin(a), z)
                        for z in z_grid for a in angles])
    z_axis = np.linspace(0.5, L - 0.5, 200)
    axial_pts = np.column_stack([np.zeros_like(z_axis), np.zeros_like(z_axis), z_axis])

    for i, f_thz in enumerate(FREQ_THZ, start=1):
        f_hz = float(f_thz) * 1e12
        # Update study step
        java = model.java
        java.study(str(java.study().tags()[0])).feature("freq").set("plist", str(f_hz))
        ts = time.monotonic()
        java.study(str(java.study().tags()[0])).run()
        e_axis = _evaluate_at_points(model, "ewfd.normE", axial_pts)
        e_ann = _evaluate_at_points(model, "ewfd.normE", ann_pts)
        dt = time.monotonic() - ts

        node_mask = (z_axis >= z_node_lo) & (z_axis <= z_node_hi)
        peak_axis = float(np.max(e_axis[node_mask])) if node_mask.any() else float("nan")
        peak_ann = float(np.max(e_ann))
        mean_ann = float(np.mean(e_ann))

        rows.append((f_thz, peak_axis, peak_ann, mean_ann, dt))
        print(f"  [{i:>2}/{len(FREQ_THZ)}] f={f_thz:.4f}  "
              f"|E|_ax={peak_axis:5.3f}  |E|_ann={peak_ann:5.3f}  ({dt:4.1f}s)")

    client.clear()

    csv_path = OUT_DIR / "spectrum.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["frequency_THz", "peak_E_node_axis", "peak_E_node_annulus",
                    "mean_E_node_annulus", "solve_s"])
        w.writerows(rows)

    arr = np.array(rows)
    f_all = arr[:, 0]
    p_ax = arr[:, 1]
    p_ann = arr[:, 2]

    # Split into the two windows
    low_mask = f_all < 1.0
    high_mask = f_all > 1.5

    f_low, ax_low, ann_low = f_all[low_mask], p_ax[low_mask], p_ann[low_mask]
    f_hi, ax_hi, ann_hi = f_all[high_mask], p_ax[high_mask], p_ann[high_mask]

    # Lorentzian fits
    print()
    print("== Low window (0.55-0.70 THz) ==")
    popt_low, perr_low = fit_peak(f_low, ax_low, f0_guess=0.62, window=0.08)
    if popt_low is not None:
        Q = popt_low[1] / abs(popt_low[2])
        print(f"  axial: f0 = {popt_low[1]:.4f} ± {perr_low[1]:.4f} THz, "
              f"γ = {abs(popt_low[2]):.4f} THz, Q = {Q:.2f}, A = {popt_low[0]:.3f}")
    popt_low_ann, perr_low_ann = fit_peak(f_low, ann_low, f0_guess=0.62, window=0.08)
    if popt_low_ann is not None:
        Q = popt_low_ann[1] / abs(popt_low_ann[2])
        print(f"  annular: f0 = {popt_low_ann[1]:.4f} ± {perr_low_ann[1]:.4f} THz, "
              f"γ = {abs(popt_low_ann[2]):.4f} THz, Q = {Q:.2f}, A = {popt_low_ann[0]:.3f}")

    print("== High window (1.85-2.15 THz) ==")
    popt_hi, perr_hi = fit_peak(f_hi, ax_hi, f0_guess=2.0, window=0.15)
    if popt_hi is not None:
        Q = popt_hi[1] / abs(popt_hi[2])
        print(f"  axial: f0 = {popt_hi[1]:.4f} ± {perr_hi[1]:.4f} THz, "
              f"γ = {abs(popt_hi[2]):.4f} THz, Q = {Q:.2f}, A = {popt_hi[0]:.3f}")
    popt_hi_ann, perr_hi_ann = fit_peak(f_hi, ann_hi, f0_guess=2.0, window=0.15)
    if popt_hi_ann is not None:
        Q = popt_hi_ann[1] / abs(popt_hi_ann[2])
        print(f"  annular: f0 = {popt_hi_ann[1]:.4f} ± {perr_hi_ann[1]:.4f} THz, "
              f"γ = {abs(popt_hi_ann[2]):.4f} THz, Q = {Q:.2f}, A = {popt_hi_ann[0]:.3f}")

    fig, (ax_l, ax_h) = plt.subplots(1, 2, figsize=(9, 3.6))

    def _plot(ax, f_w, ax_w, ann_w, popt_ax, popt_an, title, f0_marker):
        ax.plot(f_w, ax_w, "o-", color=OKABE_ITO[2], lw=1.4, ms=4, label="axial")
        ax.plot(f_w, ann_w, "s-", color=OKABE_ITO[6], lw=1.4, ms=4, label="annular")
        if popt_ax is not None:
            f_g = np.linspace(f_w.min(), f_w.max(), 200)
            ax.plot(f_g, lorentzian(f_g, *popt_ax), "--", color=OKABE_ITO[2],
                    lw=1.2, alpha=0.7,
                    label=f"axial fit: f₀={popt_ax[1]:.3f} THz, Q={popt_ax[1]/abs(popt_ax[2]):.1f}")
        if popt_an is not None:
            f_g = np.linspace(f_w.min(), f_w.max(), 200)
            ax.plot(f_g, lorentzian(f_g, *popt_an), "--", color=OKABE_ITO[6],
                    lw=1.2, alpha=0.7,
                    label=f"annular fit: f₀={popt_an[1]:.3f} THz, Q={popt_an[1]/abs(popt_an[2]):.1f}")
        ax.axvline(f0_marker, color="red", ls=":", lw=0.7, alpha=0.6)
        ax.set_xlabel("Frequency (THz)")
        ax.set_ylabel("|E| (normalised)")
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=6.5, framealpha=0.92)

    _plot(ax_l, f_low, ax_low, ann_low, popt_low, popt_low_ann,
          "Low window: 0.55–0.70 THz (refined mesh)", 0.6)
    _plot(ax_h, f_hi, ax_hi, ann_hi, popt_hi, popt_hi_ann,
          "High window: 1.85–2.15 THz (refined mesh)", 2.0)

    fig.suptitle("Sim 22 — dense frequency sweeps + Lorentzian fits in two peak windows")
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim22_dense_peaks")
    print(f"wrote {pdf}\nwrote {png}\nwrote {csv_path}")
    print(f"TOTAL: {time.monotonic() - t0:.1f}s")


if __name__ == "__main__":
    main()
