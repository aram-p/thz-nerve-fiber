"""Sim 15 — 2-D parameter sweep: peak |E| as a surface over (f, L_node).

The headline scientific result. Runs a grid of solves spanning four
node lengths × six frequencies = 24 simulations. Renders peak |E| in
the node as a 3-D surface over the (frequency, node-length) plane.
Reveals whether the 0.6 THz peak is a robust feature of the
single-fibre geometry or a sampling artefact specific to one node
length.

Expect ~6–8 minutes of compute.
"""

from __future__ import annotations

import csv
import time
from pathlib import Path

import matplotlib.pyplot as plt
import mph
import numpy as np

from thznerve.io.hdf5 import write_result
from thznerve.io.provenance import get_git_sha
from thznerve.model.geometry import GeometryParams, build_geometry, total_length_um
from thznerve.model.materials import MaterialParams, apply_materials
from thznerve.model.mesh import MeshParams, build_mesh
from thznerve.model.study import (
    _evaluate_at_points, extract_axial_profile,
    setup_physics, setup_study, solve_study,
)
from thznerve.plots.style import DOUBLE_COL_MM, apply_thesis_style, save_figure

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "results" / "sim15"

# Compact sweep: enough to map the surface, small enough to finish overnight.
NODE_LENGTHS_UM = [20.0, 40.0, 60.0, 100.0]
FREQ_THZ = np.array([0.20, 0.45, 0.63, 0.85, 1.40, 2.00])
MAT = MaterialParams(node_sigma_S_per_m=0.0)


def _geom(node_L: float) -> GeometryParams:
    return GeometryParams(
        axon_radius_um=5, myelin_radius_um=7,
        node_length_um=node_L, internode_length_um=100, external_half_width_um=20,
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    apply_thesis_style(width_mm=DOUBLE_COL_MM, aspect=0.55)
    t0 = time.monotonic()

    client = mph.start()
    git_sha = get_git_sha()

    peak_grid = np.zeros((len(NODE_LENGTHS_UM), len(FREQ_THZ)))
    mean_grid = np.zeros_like(peak_grid)

    for i, node_L in enumerate(NODE_LENGTHS_UM):
        geom = _geom(node_L)
        # Build model once per node-length, sweep frequency on the same model.
        model = client.create(f"sim15_L{node_L:g}")
        build_geometry(model, geom)
        apply_materials(model, MAT)
        setup_physics(model, geom)
        build_mesh(model, MeshParams())
        setup_study(model, FREQ_THZ[0] * 1e12)

        z_node_lo = geom.internode_length_um
        z_node_hi = z_node_lo + geom.node_length_um
        r_node = 0.5 * (geom.axon_radius_um + geom.myelin_radius_um)
        z_grid_node = np.linspace(z_node_lo + 1, z_node_hi - 1, 24)
        annulus_pts = np.column_stack([
            np.full_like(z_grid_node, r_node),
            np.zeros_like(z_grid_node),
            z_grid_node,
        ])

        for j, f_thz in enumerate(FREQ_THZ):
            f_hz = float(f_thz) * 1e12
            model.java.study(str(model.java.study().tags()[0])).feature("freq").set("plist", str(f_hz))
            ts = time.monotonic()
            model.java.study(str(model.java.study().tags()[0])).run()
            e_annulus = _evaluate_at_points(model, "ewfd.normE", annulus_pts)
            dt = time.monotonic() - ts

            peak_grid[i, j] = float(np.max(e_annulus))
            mean_grid[i, j] = float(np.mean(e_annulus))
            print(f"  L={node_L:5.1f} um  f={f_thz:5.2f} THz  |E|peak={peak_grid[i,j]:5.3f}  ({dt:4.1f}s)")
        client.remove(model)

    client.clear()

    # CSV summary
    csv_path = OUT_DIR / "freq_nodelen.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["node_length_um", "frequency_THz", "peak_E_node", "mean_E_node"])
        for i, node_L in enumerate(NODE_LENGTHS_UM):
            for j, f_thz in enumerate(FREQ_THZ):
                w.writerow([node_L, f_thz, peak_grid[i, j], mean_grid[i, j]])

    # 3-D surface
    F_grid, L_grid = np.meshgrid(FREQ_THZ, NODE_LENGTHS_UM, indexing="xy")

    fig = plt.figure()
    ax_3d = fig.add_subplot(1, 2, 1, projection="3d")
    ax_2d = fig.add_subplot(1, 2, 2)

    surf = ax_3d.plot_surface(
        F_grid, L_grid, peak_grid,
        cmap="viridis", edgecolor="0.15", linewidth=0.4, antialiased=True, alpha=0.92,
    )
    ax_3d.set_xlabel("Frequency (THz)")
    ax_3d.set_ylabel("Node length (µm)")
    ax_3d.set_zlabel("peak |E| in node")
    ax_3d.set_title("peak |E| over (f, $L_{node}$)")
    ax_3d.view_init(elev=24, azim=-62)
    ax_3d.set_box_aspect((1.4, 1.0, 0.75))

    # Mark experimentally observed peaks as vertical sticks at min L
    z_top = peak_grid.max() * 1.05
    for f_peak in (0.6, 2.0):
        ax_3d.plot([f_peak, f_peak], [NODE_LENGTHS_UM[0]] * 2, [0, z_top],
                   color="red", lw=0.9, alpha=0.7)

    # 2-D heatmap
    im = ax_2d.pcolormesh(FREQ_THZ, NODE_LENGTHS_UM, peak_grid,
                          cmap="viridis", shading="nearest")
    ax_2d.set_xlabel("Frequency (THz)")
    ax_2d.set_ylabel("Node length (µm)")
    ax_2d.set_title("peak |E| heatmap")
    ax_2d.axvline(0.6, color="red", lw=0.7, alpha=0.7, ls="--")
    ax_2d.axvline(2.0, color="red", lw=0.7, alpha=0.7, ls="--")
    fig.colorbar(im, ax=ax_2d, shrink=0.85, pad=0.02).set_label("peak |E|")

    fig.suptitle("Sim 15 — peak |E|-at-node over the (f, $L_{node}$) plane",
                 y=0.99)
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim15_freq_nodelen_surface")
    print(f"wrote {pdf}\nwrote {png}\nwrote {csv_path}")
    print(f"TOTAL: {time.monotonic() - t0:.1f}s")


if __name__ == "__main__":
    main()
