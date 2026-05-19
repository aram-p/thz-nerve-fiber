"""Sim 10 — 3D nested isosurfaces of |E| at the resonant frequency.

Replaces sim 9's scatter rendering with proper iso-surfaces (computed
via marching cubes) so the spatial structure of the resonant scattered
field reads cleanly. Three nested transparent shells at |E|/max =
0.55, 0.70, 0.85 reveal the "skin" of the high-field region.
"""

from __future__ import annotations

import time
from pathlib import Path

import matplotlib.pyplot as plt
import mph
import numpy as np
from matplotlib import cm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage import measure

from thznerve.model.geometry import GeometryParams, build_geometry, total_length_um
from thznerve.model.materials import MaterialParams, apply_materials
from thznerve.model.mesh import MeshParams, build_mesh
from thznerve.model.study import (
    _evaluate_at_points, setup_physics, setup_study, solve_study,
)
from thznerve.plots.style import DOUBLE_COL_MM, OKABE_ITO, apply_thesis_style, save_figure

REPO_ROOT = Path(__file__).resolve().parent.parent

GEOM = GeometryParams(
    axon_radius_um=5, myelin_radius_um=7,
    node_length_um=40, internode_length_um=100, external_half_width_um=20,
)
MAT = MaterialParams(node_sigma_S_per_m=0.0)
MESH = MeshParams()

F_HZ = 0.632e12

# Sampling grid — finer than sim 9 since isosurfaces benefit from it.
NX, NY, NZ = 28, 28, 70


def _ring(ax, r, z, n=48, **kw):
    th = np.linspace(0, 2 * np.pi, n)
    ax.plot(r * np.cos(th), r * np.sin(th), np.full_like(th, z), **kw)


def main() -> None:
    apply_thesis_style(width_mm=DOUBLE_COL_MM, aspect=0.5)

    t0 = time.monotonic()
    client = mph.start()
    model = client.create("sim10_iso")
    build_geometry(model, GEOM)
    apply_materials(model, MAT)
    setup_physics(model, GEOM)
    build_mesh(model, MESH)
    setup_study(model, F_HZ)
    print(f"setup: {time.monotonic() - t0:.1f}s")
    solve_study(model)
    print(f"solved: {time.monotonic() - t0:.1f}s")

    L = total_length_um(GEOM)
    hw = GEOM.external_half_width_um
    x = np.linspace(-hw + 0.6, hw - 0.6, NX)
    y = np.linspace(-hw + 0.6, hw - 0.6, NY)
    z = np.linspace(0.6, L - 0.6, NZ)
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    pts = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    print(f"sampling {pts.shape[0]} points...")
    e_3d = _evaluate_at_points(model, "ewfd.normE", pts).reshape(NX, NY, NZ)
    print(f"|E| range: {e_3d.min():.3g} – {e_3d.max():.3g}")
    client.clear()

    e_max = float(e_3d.max())
    dx = float(x[1] - x[0])
    dy = float(y[1] - y[0])
    dz = float(z[1] - z[0])

    fig = plt.figure()
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    ax_2d = fig.add_subplot(1, 2, 2)

    levels = [0.55, 0.70, 0.85]
    colors = [OKABE_ITO[6], OKABE_ITO[5], OKABE_ITO[1]]
    alphas = [0.18, 0.32, 0.55]

    for level, col, alpha in zip(levels, colors, alphas):
        iso = level * e_max
        try:
            verts, faces, _normals, _vals = measure.marching_cubes(
                e_3d, level=iso, spacing=(dx, dy, dz),
            )
        except (ValueError, RuntimeError):
            continue
        # Shift verts to actual coordinates (grid starts at (x[0], y[0], z[0]))
        verts += np.array([x[0], y[0], z[0]])
        mesh = Poly3DCollection(
            verts[faces],
            facecolor=col, alpha=alpha, edgecolor="none",
        )
        ax.add_collection3d(mesh)

    # Box outline
    xs = [-hw, hw, hw, -hw, -hw]
    ys = [-hw, -hw, hw, hw, -hw]
    for zc in (0, L):
        ax.plot(xs, ys, [zc] * 5, color="0.55", lw=0.6)
    for xx, yy in zip(xs[:-1], ys[:-1]):
        ax.plot([xx, xx], [yy, yy], [0, L], color="0.55", lw=0.5)

    # Fibre outline: a few rings of the axon and myelin
    for z_ring in np.linspace(0, L, 7):
        _ring(ax, GEOM.axon_radius_um, z_ring,
              color=OKABE_ITO[2], lw=0.6, alpha=0.55)
        _ring(ax, GEOM.myelin_radius_um, z_ring,
              color=OKABE_ITO[3], lw=0.5, alpha=0.45)
    # Node-region full rings
    z_n0 = GEOM.internode_length_um
    z_n1 = z_n0 + GEOM.node_length_um
    for r in (GEOM.axon_radius_um, GEOM.myelin_radius_um):
        for z_ring in (z_n0, z_n1):
            _ring(ax, r, z_ring, color=OKABE_ITO[6], lw=1.4, alpha=0.95)

    # Legend proxies
    for level, col, alpha in zip(levels, colors, alphas):
        ax.plot([], [], color=col, lw=8, alpha=alpha,
                label=fr"|E|/max = {level}")
    ax.plot([], [], color=OKABE_ITO[6], lw=2, label="node boundary")
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")
    ax.set_zlabel("z (µm)")
    ax.set_box_aspect((1, 1, 2.2))
    ax.view_init(elev=20, azim=-58)
    ax.set_title(f"Nested |E| isosurfaces at f = {F_HZ/1e12:g} THz")
    ax.legend(loc="upper left", fontsize=6.5, framealpha=0.92)

    # Companion: integrate |E| over (x, y) for each z to get an axial-energy
    # profile that captures the full cross-section (not just on-axis).
    from scipy.integrate import trapezoid
    e_integrated_xy = trapezoid(trapezoid(e_3d ** 2, x, axis=0), y, axis=0)
    ax_2d.plot(z, e_integrated_xy, color=OKABE_ITO[6], lw=1.6,
               label=r"$\iint |E|^2 \, dx\,dy$")
    ax_2d.axvspan(z_n0, z_n1, color=OKABE_ITO[6], alpha=0.12)
    ax_2d.set_xlabel("z (µm)")
    ax_2d.set_ylabel(r"$\int\!\!\int |E|^2\, dx\, dy$  (µm$^2$)")
    ax_2d.set_title("Cross-section-integrated energy density vs z")
    ax_2d.grid(True, alpha=0.3)
    ax_2d.legend(loc="best", framealpha=0.92)

    fig.suptitle("Sim 10 — Nested |E| isosurfaces + axial-energy profile")
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim10_field_isosurfaces")
    print(f"wrote {pdf}\nwrote {png}")
    print(f"TOTAL: {time.monotonic() - t0:.1f}s")


if __name__ == "__main__":
    main()
