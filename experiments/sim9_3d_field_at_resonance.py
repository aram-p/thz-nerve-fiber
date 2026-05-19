"""Sim 9 — 3D |E| field render at the resonance frequency.

Re-solves the EWFD problem at f = 0.632 THz (the peak found in sim 1's
26-point sweep) and samples |E| on a 3-D grid. Renders two panels:

  * Left:  3D scatter — every grid point coloured / alpha-blended by |E|,
    with the cylindrical fibre outline overlaid.
  * Right: y = 0 cross-section heatmap of |E| — the classic "axial slice"
    view, but at the resonance instead of an arbitrary frequency.

This is the figure that shows *where* in the geometry the resonant field
enhancement lives.
"""

from __future__ import annotations

import time
from pathlib import Path

import matplotlib.pyplot as plt
import mph
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from thznerve.model.geometry import GeometryParams, build_geometry, total_length_um
from thznerve.model.materials import MaterialParams, apply_materials
from thznerve.model.mesh import MeshParams, build_mesh
from thznerve.model.study import (
    _evaluate_at_points,
    setup_physics,
    setup_study,
    solve_study,
)
from thznerve.plots.style import DOUBLE_COL_MM, OKABE_ITO, apply_thesis_style, save_figure

REPO_ROOT = Path(__file__).resolve().parent.parent

GEOM = GeometryParams(
    axon_radius_um=5,
    myelin_radius_um=7,
    node_length_um=40,
    internode_length_um=100,
    external_half_width_um=20,
)
MAT = MaterialParams(node_sigma_S_per_m=0.0)
MESH = MeshParams(max_element_size_um=30.0, node_element_size_um=5.0)

# Resonance frequency from sim 1 — the highest |E|-at-node in 0.1–2 THz.
F_HZ = 0.632e12


def main() -> None:
    apply_thesis_style(width_mm=DOUBLE_COL_MM, aspect=0.55)

    t0 = time.monotonic()
    client = mph.start()
    model = client.create("sim9_3d_field")
    build_geometry(model, GEOM)
    apply_materials(model, MAT)
    setup_physics(model, GEOM)
    n_elem = build_mesh(model, MESH)
    setup_study(model, F_HZ)
    print(f"setup: {n_elem} elements, {time.monotonic() - t0:.1f}s")
    solve_study(model)
    print(f"solved in {time.monotonic() - t0:.1f}s")

    L = total_length_um(GEOM)
    hw = GEOM.external_half_width_um

    # 3D grid (coarse for scatter; finer for slice).
    nx, ny, nz = 14, 14, 50
    x = np.linspace(-hw + 0.7, hw - 0.7, nx)
    y = np.linspace(-hw + 0.7, hw - 0.7, ny)
    z = np.linspace(0.7, L - 0.7, nz)
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    pts = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    print(f"3D grid: {pts.shape[0]} points, sampling...")
    e_3d = _evaluate_at_points(model, "ewfd.normE", pts).reshape(nx, ny, nz)
    print(f"  |E|_3d: min={e_3d.min():.3g} max={e_3d.max():.3g}")

    # Higher-resolution y=0 slice
    nx_slice, nz_slice = 90, 200
    x_slice = np.linspace(-hw + 0.3, hw - 0.3, nx_slice)
    z_slice = np.linspace(0.3, L - 0.3, nz_slice)
    Xs, Zs = np.meshgrid(x_slice, z_slice, indexing="ij")
    pts_slice = np.column_stack([Xs.ravel(), np.zeros(Xs.size), Zs.ravel()])
    print(f"slice grid: {pts_slice.shape[0]} points, sampling...")
    e_slice = _evaluate_at_points(model, "ewfd.normE", pts_slice).reshape(nx_slice, nz_slice)

    client.clear()

    # --------------------------------------------------------------- plot
    cmap_hot = "magma"
    e_max = float(np.max(e_slice))

    fig = plt.figure()
    ax_3d = fig.add_subplot(1, 2, 1, projection="3d")
    ax_2d = fig.add_subplot(1, 2, 2)

    # ---- 3D scatter — only the highest-|E| points so the figure reads cleanly
    e_flat = e_3d.ravel()
    norm = e_flat / e_flat.max()
    mask = norm > 0.65
    sc = ax_3d.scatter(
        X.ravel()[mask], Y.ravel()[mask], Z.ravel()[mask],
        c=e_flat[mask], cmap=cmap_hot, alpha=0.7,
        s=22 * (norm[mask] ** 2.5) + 2, linewidths=0, vmin=0, vmax=e_max,
    )
    # Box wireframe
    xs = [-hw, hw, hw, -hw, -hw]
    ys = [-hw, -hw, hw, hw, -hw]
    for zc in (0, L):
        ax_3d.plot(xs, ys, [zc] * 5, color="0.55", lw=0.6)
    for xx, yy in zip(xs[:-1], ys[:-1]):
        ax_3d.plot([xx, xx], [yy, yy], [0, L], color="0.55", lw=0.5)

    # Sketch the fibre outline (axon + myelin) as thin guide rings
    theta = np.linspace(0, 2 * np.pi, 40)
    for r, col, lw in [(GEOM.axon_radius_um, OKABE_ITO[2], 0.7),
                       (GEOM.myelin_radius_um, OKABE_ITO[1], 0.6)]:
        for z_ring in np.linspace(0, L, 6):
            ax_3d.plot(r * np.cos(theta), r * np.sin(theta),
                       np.full_like(theta, z_ring), color=col, lw=lw, alpha=0.7)
    # Highlight the node region with full ring outlines
    z_n0 = GEOM.internode_length_um
    z_n1 = z_n0 + GEOM.node_length_um
    for r in (GEOM.axon_radius_um, GEOM.myelin_radius_um):
        for z_ring in (z_n0, z_n1):
            ax_3d.plot(r * np.cos(theta), r * np.sin(theta),
                       np.full_like(theta, z_ring), color=OKABE_ITO[6],
                       lw=1.4, alpha=0.95)

    ax_3d.set_xlabel("x (µm)")
    ax_3d.set_ylabel("y (µm)")
    ax_3d.set_zlabel("z (µm)")
    ax_3d.set_box_aspect((1, 1, 2.2))
    ax_3d.view_init(elev=18, azim=-58)
    ax_3d.set_title(f"3D |E| (top 35%, |E|/max > 0.65), f = {F_HZ/1e12:g} THz")

    cbar = fig.colorbar(sc, ax=ax_3d, shrink=0.6, pad=0.08)
    cbar.set_label("|E|")

    # ---- 2D slice — Gaussian smoothing to suppress sub-mesh-element FEM
    # interpolation artefacts so the resonance pattern reads clearly.
    from scipy.ndimage import gaussian_filter
    e_slice_smooth = gaussian_filter(e_slice, sigma=(1.2, 2.0))
    im = ax_2d.pcolormesh(z_slice, x_slice, e_slice_smooth, cmap=cmap_hot,
                          vmin=0, vmax=e_max, shading="auto")
    # Fibre outline (axon and myelin radii)
    for r, col, lw in [(GEOM.axon_radius_um, "white", 0.7),
                       (GEOM.myelin_radius_um, "white", 0.5)]:
        ax_2d.axhline( r, color=col, lw=lw, alpha=0.55)
        ax_2d.axhline(-r, color=col, lw=lw, alpha=0.55)
    # Highlight node region in z
    z_n0 = GEOM.internode_length_um
    z_n1 = z_n0 + GEOM.node_length_um
    ax_2d.axvline(z_n0, color=OKABE_ITO[6], lw=1.2, alpha=0.85)
    ax_2d.axvline(z_n1, color=OKABE_ITO[6], lw=1.2, alpha=0.85)
    ax_2d.text((z_n0 + z_n1) / 2, hw * 0.86, "node",
               color=OKABE_ITO[6], ha="center", fontsize=8,
               fontweight="bold")
    ax_2d.text(z_n0 / 2, hw * 0.86, "myelin",
               color="white", ha="center", fontsize=7, alpha=0.75)
    ax_2d.text((z_n1 + L) / 2, hw * 0.86, "myelin",
               color="white", ha="center", fontsize=7, alpha=0.75)
    ax_2d.set_xlabel("z (µm)")
    ax_2d.set_ylabel("x (µm)")
    ax_2d.set_title(f"|E|(x, z) at y = 0, f = {F_HZ/1e12:g} THz")
    fig.colorbar(im, ax=ax_2d, shrink=0.85, pad=0.02).set_label("|E|")

    fig.suptitle("Sim 9 — 3D / cross-sectional |E| at the resonant frequency")
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim9_3d_field_at_resonance")
    print(f"wrote {pdf}\nwrote {png}")
    print(f"TOTAL: {time.monotonic() - t0:.1f}s")


if __name__ == "__main__":
    main()
