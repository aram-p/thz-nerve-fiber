"""Sim 16 — 3D mesh convergence at the resonance frequency.

A simulation result is only as trustworthy as its mesh. Sim 1's 0.6 THz
peak might be a real physics feature or it might be sub-resolution
mesh noise. This sim re-solves at f = 0.632 THz at three mesh
refinements — coarse, baseline, fine — and renders the |E| field for
each so the convergence (or lack of it) is visually obvious.

Refinement factor scales 1 / (mesh element size), so:
  factor 0.5 → larger elements, ~13 k elements (coarse)
  factor 1.0 → baseline,        ~26 k elements (sim 1 / sim 9 setting)
  factor 1.8 → smaller elements, ~80–100 k elements (fine)

3 fresh models, ~2–4 minutes total.
"""

from __future__ import annotations

import time
from pathlib import Path

import matplotlib.pyplot as plt
import mph
import numpy as np
from scipy.ndimage import gaussian_filter
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
F_HZ = 0.632e12

REFINEMENTS = [
    ("coarse",   MeshParams(max_element_size_um=40.0, node_element_size_um=10.0, refinement_factor=0.6)),
    ("baseline", MeshParams(max_element_size_um=30.0, node_element_size_um=5.0,  refinement_factor=1.0)),
    ("fine",     MeshParams(max_element_size_um=18.0, node_element_size_um=3.0,  refinement_factor=1.5)),
]

NX, NY, NZ = 16, 16, 50
NX_SLICE, NZ_SLICE = 80, 180


def _ring(ax, r, z, n=40, **kw):
    th = np.linspace(0, 2 * np.pi, n)
    ax.plot(r * np.cos(th), r * np.sin(th), np.full_like(th, z), **kw)


def _sample(model, hw, L):
    x = np.linspace(-hw + 0.7, hw - 0.7, NX)
    y = np.linspace(-hw + 0.7, hw - 0.7, NY)
    z = np.linspace(0.7, L - 0.7, NZ)
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    pts = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    e_3d = _evaluate_at_points(model, "ewfd.normE", pts).reshape(NX, NY, NZ)

    xs = np.linspace(-hw + 0.3, hw - 0.3, NX_SLICE)
    zs = np.linspace(0.3, L - 0.3, NZ_SLICE)
    Xs, Zs = np.meshgrid(xs, zs, indexing="ij")
    pts_slice = np.column_stack([Xs.ravel(), np.zeros(Xs.size), Zs.ravel()])
    e_slice = _evaluate_at_points(model, "ewfd.normE", pts_slice).reshape(NX_SLICE, NZ_SLICE)
    return x, y, z, e_3d, xs, zs, e_slice


def main() -> None:
    apply_thesis_style(width_mm=DOUBLE_COL_MM, aspect=0.65)
    t0 = time.monotonic()
    client = mph.start()

    L = total_length_um(GEOM)
    hw = GEOM.external_half_width_um

    results = []
    for label, mesh_params in REFINEMENTS:
        model = client.create(f"sim16_{label}")
        build_geometry(model, GEOM)
        apply_materials(model, MAT)
        setup_physics(model, GEOM)
        n_elem = build_mesh(model, mesh_params)
        setup_study(model, F_HZ)
        print(f"[{label}] {n_elem} elements, solving...")
        ts = time.monotonic()
        solve_study(model)
        x, y, z, e_3d, xs, zs, e_slice = _sample(model, hw, L)
        dt = time.monotonic() - ts
        # Sample on axis through node to get a 1D profile too
        z_node = np.linspace(0.7, L - 0.7, 200)
        axis_pts = np.column_stack([np.zeros_like(z_node), np.zeros_like(z_node), z_node])
        e_axis = _evaluate_at_points(model, "ewfd.normE", axis_pts)
        results.append({
            "label": label, "n_elem": n_elem, "dt": dt,
            "x": x, "y": y, "z": z, "e_3d": e_3d,
            "xs": xs, "zs": zs, "e_slice": e_slice,
            "z_axis": z_node, "e_axis": e_axis,
            "peak_axis": float(np.max(e_axis)),
            "peak_3d": float(np.max(e_3d)),
        })
        print(f"  [{label}] peak |E|_axis = {results[-1]['peak_axis']:.4f}, "
              f"peak |E|_3d = {results[-1]['peak_3d']:.4f}, "
              f"solved in {dt:.1f}s")
        client.remove(model)

    client.clear()

    # Use the fine result for color scale
    e_max_slice = max(r["e_slice"].max() for r in results)

    fig = plt.figure()
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 0.9])
    ax_2d = [fig.add_subplot(gs[0, k]) for k in range(3)]
    ax_3d = [fig.add_subplot(gs[1, k], projection="3d") for k in range(3)]

    for k, res in enumerate(results):
        label = res["label"]
        # 2D slice — Gaussian smoothing for visibility
        e_smooth = gaussian_filter(res["e_slice"], sigma=(1.0, 1.6))
        im = ax_2d[k].pcolormesh(res["zs"], res["xs"], e_smooth, cmap="magma",
                                  vmin=0, vmax=e_max_slice, shading="auto")
        z_n0 = GEOM.internode_length_um
        z_n1 = z_n0 + GEOM.node_length_um
        ax_2d[k].axvline(z_n0, color=OKABE_ITO[6], lw=0.9, alpha=0.85)
        ax_2d[k].axvline(z_n1, color=OKABE_ITO[6], lw=0.9, alpha=0.85)
        for r in (GEOM.axon_radius_um, GEOM.myelin_radius_um):
            ax_2d[k].axhline( r, color="white", lw=0.5, alpha=0.55)
            ax_2d[k].axhline(-r, color="white", lw=0.5, alpha=0.55)
        ax_2d[k].set_title(
            f"{label}\n{res['n_elem']:>6,} elements  ·  peak |E|_axis = {res['peak_axis']:.3f}",
            fontsize=8,
        )
        ax_2d[k].set_xlabel("z (µm)", fontsize=7)
        if k == 0:
            ax_2d[k].set_ylabel("x (µm)", fontsize=7)
        ax_2d[k].tick_params(labelsize=6)

        # 3D scatter (top-third |E|) so the resolution change reads at a glance
        e_3d = res["e_3d"]
        e_max = e_3d.max()
        x_g, y_g, z_g = np.meshgrid(res["x"], res["y"], res["z"], indexing="ij")
        e_flat = e_3d.ravel()
        norm = e_flat / e_max
        mask = norm > 0.7
        ax_3d[k].scatter(
            x_g.ravel()[mask], y_g.ravel()[mask], z_g.ravel()[mask],
            c=e_flat[mask], cmap="magma", alpha=0.75,
            s=18 * (norm[mask] ** 2.5) + 2, linewidths=0, vmin=0, vmax=e_max_slice,
        )
        # Highlight node rings
        for r in (GEOM.axon_radius_um, GEOM.myelin_radius_um):
            for z_ring in (z_n0, z_n1):
                _ring(ax_3d[k], r, z_ring, color=OKABE_ITO[6], lw=1.1, alpha=0.9)
        ax_3d[k].set_xlabel("x", fontsize=6, labelpad=-2)
        ax_3d[k].set_ylabel("y", fontsize=6, labelpad=-2)
        ax_3d[k].set_zlabel("z", fontsize=6, labelpad=-2)
        ax_3d[k].tick_params(labelsize=5)
        ax_3d[k].view_init(elev=22, azim=-58)
        ax_3d[k].set_box_aspect((1, 1, 2.0))
        ax_3d[k].set_title(f"3-D top-30% |E| ({label})", fontsize=7)

    fig.suptitle(
        f"Sim 16 — Mesh convergence at f = {F_HZ/1e12:g} THz "
        f"(coarse / baseline / fine, σ = 0)",
        y=0.99,
    )
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim16_mesh_convergence_3d")
    print(f"wrote {pdf}\nwrote {png}")
    print(f"TOTAL: {time.monotonic() - t0:.1f}s")


if __name__ == "__main__":
    main()
