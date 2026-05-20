"""Sim 25 — |E| field structure at both resonance peaks.

Solves at f = 0.619 THz (low peak, annular mode, half-wave on
~120 µm segment) and f = 1.938 THz (high peak, axial mode, half-wave
on ~40 µm node). Renders |E| 3D scatter + 2D x-z cross-section
side-by-side so the spatial structure of each mode is visible.

Expected visual confirmation (per Sim 24):
* 0.619 THz: one half-wave spanning the internode-node segment;
  bright lobes near z=50 and z=140 (sheath / node-junction).
* 1.938 THz: tighter half-wave confined to the node (z=100-140);
  bright lobe at the node centre.
"""

from __future__ import annotations

import time
from pathlib import Path

import matplotlib.pyplot as plt
import mph
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.ndimage import gaussian_filter
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
MESH = MeshParams(max_element_size_um=15.0, node_element_size_um=3.0)  # refined

PEAKS = [
    {"label": "low",  "f_THz": 0.619, "predicted_segment_um": 120.0,
     "mode_name": "λ/2 on internode + half-node (120 µm)"},
    {"label": "high", "f_THz": 1.938, "predicted_segment_um": 40.0,
     "mode_name": "λ/2 on node only (40 µm)"},
]

NX, NY, NZ = 16, 16, 60
NX_SLICE, NZ_SLICE = 90, 220


def _ring(ax, r, z, n=40, **kw):
    th = np.linspace(0, 2 * np.pi, n)
    ax.plot(r * np.cos(th), r * np.sin(th), np.full_like(th, z), **kw)


def _sample(model, hw, L):
    x = np.linspace(-hw + 0.6, hw - 0.6, NX)
    y = np.linspace(-hw + 0.6, hw - 0.6, NY)
    z = np.linspace(0.6, L - 0.6, NZ)
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    pts = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    e_3d = _evaluate_at_points(model, "ewfd.normE", pts).reshape(NX, NY, NZ)
    xs = np.linspace(-hw + 0.3, hw - 0.3, NX_SLICE)
    zs = np.linspace(0.3, L - 0.3, NZ_SLICE)
    Xs, Zs = np.meshgrid(xs, zs, indexing="ij")
    pts_s = np.column_stack([Xs.ravel(), np.zeros(Xs.size), Zs.ravel()])
    e_slice = _evaluate_at_points(model, "ewfd.normE", pts_s).reshape(NX_SLICE, NZ_SLICE)
    return x, y, z, e_3d, xs, zs, e_slice


def _render_iso(ax, x, y, z, e_3d, *, e_max, hw, L):
    dx, dy, dz = float(x[1] - x[0]), float(y[1] - y[0]), float(z[1] - z[0])
    levels = [0.55, 0.72, 0.88]
    colors = [OKABE_ITO[6], OKABE_ITO[5], OKABE_ITO[1]]
    alphas = [0.18, 0.30, 0.50]
    for level, col, alpha in zip(levels, colors, alphas):
        try:
            verts, faces, _, _ = measure.marching_cubes(e_3d, level=level * e_max,
                                                        spacing=(dx, dy, dz))
        except (ValueError, RuntimeError):
            continue
        verts += np.array([x[0], y[0], z[0]])
        ax.add_collection3d(Poly3DCollection(
            verts[faces], facecolor=col, alpha=alpha, edgecolor="none"))
    xs = [-hw, hw, hw, -hw, -hw]
    ys = [-hw, -hw, hw, hw, -hw]
    for zc in (0, L):
        ax.plot(xs, ys, [zc] * 5, color="0.55", lw=0.5)
    z_n0 = GEOM.internode_length_um
    z_n1 = z_n0 + GEOM.node_length_um
    for r in (GEOM.axon_radius_um, GEOM.myelin_radius_um):
        for z_ring in (z_n0, z_n1):
            _ring(ax, r, z_ring, color=OKABE_ITO[6], lw=1.3, alpha=0.9)
    ax.set_xlabel("x", fontsize=6, labelpad=-2)
    ax.set_ylabel("y", fontsize=6, labelpad=-2)
    ax.set_zlabel("z", fontsize=6, labelpad=-2)
    ax.tick_params(labelsize=5)
    ax.set_box_aspect((1, 1, 2.0))
    ax.view_init(elev=20, azim=-58)


def main() -> None:
    apply_thesis_style(width_mm=DOUBLE_COL_MM, aspect=0.85)
    t0 = time.monotonic()
    client = mph.start()

    L = total_length_um(GEOM)
    hw = GEOM.external_half_width_um
    results = {}

    for peak in PEAKS:
        label = peak["label"]
        f_hz = peak["f_THz"] * 1e12
        model = client.create(f"sim25_{label}")
        build_geometry(model, GEOM)
        apply_materials(model, MAT, freq_hz=f_hz)
        setup_physics(model, GEOM)
        n_elem = build_mesh(model, MESH)
        setup_study(model, f_hz)
        ts = time.monotonic()
        solve_study(model)
        print(f"[{label}] f = {f_hz/1e12:.3f} THz, solved in {time.monotonic() - ts:.1f}s "
              f"({n_elem} elements), sampling...")
        results[label] = _sample(model, hw, L)
        client.remove(model)

    client.clear()

    e_max_3d = max(results[p["label"]][3].max() for p in PEAKS)
    e_max_slice = max(results[p["label"]][6].max() for p in PEAKS)

    fig = plt.figure()
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.3], height_ratios=[1, 1])

    for row, peak in enumerate(PEAKS):
        label = peak["label"]
        f_thz = peak["f_THz"]
        x, y, z, e_3d, xs, zs, e_slice = results[label]

        ax_3d = fig.add_subplot(gs[row, 0], projection="3d")
        ax_2d = fig.add_subplot(gs[row, 1])
        _render_iso(ax_3d, x, y, z, e_3d, e_max=e_max_3d, hw=hw, L=L)
        ax_3d.set_title(
            f"{peak['label']} peak — f = {f_thz:.3f} THz\n{peak['mode_name']}",
            fontsize=7.5,
        )

        e_smooth = gaussian_filter(e_slice, sigma=(1.0, 1.6))
        im = ax_2d.pcolormesh(zs, xs, e_smooth, cmap="magma",
                              vmin=0, vmax=e_max_slice, shading="auto")
        for r in (GEOM.axon_radius_um, GEOM.myelin_radius_um):
            ax_2d.axhline( r, color="white", lw=0.6, alpha=0.55)
            ax_2d.axhline(-r, color="white", lw=0.6, alpha=0.55)
        z_n0 = GEOM.internode_length_um
        z_n1 = z_n0 + GEOM.node_length_um
        ax_2d.axvline(z_n0, color=OKABE_ITO[6], lw=1.0, alpha=0.85)
        ax_2d.axvline(z_n1, color=OKABE_ITO[6], lw=1.0, alpha=0.85)
        # Predicted standing-mode wavelength as horizontal lines
        seg = peak["predicted_segment_um"]
        # mark the segment along z
        if label == "high":
            z_seg_lo, z_seg_hi = z_n0, z_n1  # node only
        else:
            z_seg_lo = z_n0 - 0.5 * GEOM.node_length_um  # internode + half-node = 100 + 20 = 120
            z_seg_hi = z_seg_lo + 120
        ax_2d.annotate(
            f"λ/2 = {seg:g} µm", xy=((z_seg_lo + z_seg_hi) / 2, hw * 0.85),
            color="white", ha="center", fontsize=7,
            bbox=dict(boxstyle="round,pad=0.2", fc="black", alpha=0.45),
        )
        ax_2d.set_xlabel("z (µm)", fontsize=8)
        ax_2d.set_ylabel("x (µm)", fontsize=8)
        ax_2d.set_title(f"y = 0 cross-section, {peak['label']} peak",
                        fontsize=8)
        fig.colorbar(im, ax=ax_2d, shrink=0.85, pad=0.02).set_label("|E|")

    fig.suptitle("Sim 25 — |E| field structure at the two peaks (visual confirmation of mode identification)",
                 y=0.995)
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim25_field_at_both_peaks")
    print(f"wrote {pdf}\nwrote {png}")
    print(f"TOTAL: {time.monotonic() - t0:.1f}s")


if __name__ == "__main__":
    main()
