"""Sim 26 — Mesh convergence at both peak frequencies.

Sim 16 ran mesh-convergence only at f = 0.632 THz. Sim 22 found the
two real peaks at f = 0.619 and 1.938 THz with a refined (15 µm)
mesh. This sim verifies those peak |E| values survive *further* mesh
refinement, so the reported Q values aren't mesh artefacts.

Three mesh refinements × two frequencies = six solves, ~3-5 min total.
"""

from __future__ import annotations

import csv
import time
from pathlib import Path

import matplotlib.pyplot as plt
import mph
import numpy as np

from thznerve.model.geometry import GeometryParams, build_geometry, total_length_um
from thznerve.model.materials import MaterialParams, apply_materials
from thznerve.model.mesh import MeshParams, build_mesh
from thznerve.model.study import (
    _evaluate_at_points, setup_physics, setup_study, solve_study,
)
from thznerve.plots.style import DOUBLE_COL_MM, OKABE_ITO, apply_thesis_style, save_figure

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "results" / "sim26"

GEOM = GeometryParams(
    axon_radius_um=5, myelin_radius_um=7,
    node_length_um=40, internode_length_um=100, external_half_width_um=20,
)
MAT = MaterialParams(node_sigma_S_per_m=0.0)

FREQS_THZ = [0.619, 1.938]
MESHES = [
    ("coarse",   MeshParams(max_element_size_um=30.0, node_element_size_um=8.0,  refinement_factor=1.0)),
    ("baseline", MeshParams(max_element_size_um=15.0, node_element_size_um=3.0,  refinement_factor=1.0)),
    ("fine",     MeshParams(max_element_size_um=10.0, node_element_size_um=2.0,  refinement_factor=1.0)),
]


def sample_field(model):
    L = total_length_um(GEOM)
    z_node_lo = GEOM.internode_length_um
    z_node_hi = z_node_lo + GEOM.node_length_um
    r_node = 0.5 * (GEOM.axon_radius_um + GEOM.myelin_radius_um)

    z_grid = np.linspace(z_node_lo + 1, z_node_hi - 1, 24)
    angles = np.linspace(0, 2 * np.pi, 6, endpoint=False)
    ann_pts = np.array([(r_node * np.cos(a), r_node * np.sin(a), z)
                        for z in z_grid for a in angles])
    z_ax = np.linspace(0.5, L - 0.5, 200)
    ax_pts = np.column_stack([np.zeros_like(z_ax), np.zeros_like(z_ax), z_ax])

    e_ax = _evaluate_at_points(model, "ewfd.normE", ax_pts)
    e_ann = _evaluate_at_points(model, "ewfd.normE", ann_pts)
    mask = (z_ax >= z_node_lo) & (z_ax <= z_node_hi)
    return {
        "peak_E_node_axis": float(np.max(e_ax[mask])) if mask.any() else float("nan"),
        "peak_E_node_annulus": float(np.max(e_ann)),
        "mean_E_node_annulus": float(np.mean(e_ann)),
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    apply_thesis_style(width_mm=DOUBLE_COL_MM, aspect=0.5)
    t0 = time.monotonic()

    client = mph.start()
    rows = []
    for f_thz in FREQS_THZ:
        f_hz = float(f_thz) * 1e12
        for label, mesh_params in MESHES:
            model = client.create(f"sim26_f{f_thz:.3f}_{label}")
            build_geometry(model, GEOM)
            apply_materials(model, MAT, freq_hz=f_hz)
            setup_physics(model, GEOM)
            n_elem = build_mesh(model, mesh_params)
            setup_study(model, f_hz)
            ts = time.monotonic()
            solve_study(model)
            data = sample_field(model)
            dt = time.monotonic() - ts
            client.remove(model)
            row = (f_thz, label, n_elem, data["peak_E_node_axis"],
                   data["peak_E_node_annulus"], data["mean_E_node_annulus"], dt)
            rows.append(row)
            print(f"  f={f_thz:.3f}  {label:>8}  {n_elem:>6} elem  "
                  f"|E|_ax={data['peak_E_node_axis']:6.3f}  "
                  f"|E|_ann={data['peak_E_node_annulus']:6.3f}  "
                  f"({dt:5.1f}s)")
    client.clear()

    csv_path = OUT_DIR / "convergence.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["frequency_THz", "mesh", "n_elem",
                    "peak_E_node_axis", "peak_E_node_annulus",
                    "mean_E_node_annulus", "solve_s"])
        w.writerows(rows)

    arr = np.array(rows, dtype=object)

    fig, axes = plt.subplots(1, 2)
    for ax, f_thz in zip(axes, FREQS_THZ):
        mask = arr[:, 0] == f_thz
        sub = arr[mask]
        n_elem = sub[:, 2].astype(float)
        peak_ax = sub[:, 3].astype(float)
        peak_ann = sub[:, 4].astype(float)
        mesh_labels = sub[:, 1]
        ax.plot(n_elem, peak_ax, "o-", color=OKABE_ITO[2], lw=1.5, ms=6,
                label="peak |E| axial")
        ax.plot(n_elem, peak_ann, "s-", color=OKABE_ITO[6], lw=1.5, ms=6,
                label="peak |E| annular")
        for x, y, lbl in zip(n_elem, peak_ax, mesh_labels):
            ax.annotate(str(lbl), xy=(x, y), xytext=(4, 4), textcoords="offset points",
                        fontsize=6.5, color=OKABE_ITO[2])
        ax.set_xlabel("Number of elements")
        ax.set_ylabel("Peak |E| in node")
        ax.set_xscale("log")
        ax.set_title(f"f = {f_thz:.3f} THz")
        ax.legend(loc="best", fontsize=7, framealpha=0.92)
        ax.grid(True, which="both", alpha=0.25)

    fig.suptitle("Sim 26 — Mesh convergence at both peak frequencies")
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim26_convergence_both_peaks")
    print(f"wrote {pdf}\nwrote {png}\nwrote {csv_path}")
    print(f"TOTAL: {time.monotonic() - t0:.1f}s")


if __name__ == "__main__":
    main()
