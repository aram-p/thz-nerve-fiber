"""Sim 13 — On-resonance vs off-resonance 3D field comparison.

Solves the EWFD problem at the resonant frequency 0.632 THz (the peak
from sim 1) and at an off-resonance frequency 0.328 THz (a local
minimum in sim 1's peak |E| trace). Renders both fields side-by-side
as 3D isosurfaces + x-z slice heatmaps. The visual difference is the
*operational definition* of "what makes this frequency resonant".
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
MESH = MeshParams()

F_ON  = 0.632e12   # resonance (sim 1 peak)
F_OFF = 0.328e12   # off-resonance (sim 1 local minimum)

NX, NY, NZ = 20, 20, 60
NX_SLICE, NZ_SLICE = 70, 160


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
    # 2D slice (finer)
    xs = np.linspace(-hw + 0.3, hw - 0.3, NX_SLICE)
    zs = np.linspace(0.3, L - 0.3, NZ_SLICE)
    Xs, Zs = np.meshgrid(xs, zs, indexing="ij")
    pts_slice = np.column_stack([Xs.ravel(), np.zeros(Xs.size), Zs.ravel()])
    e_slice = _evaluate_at_points(model, "ewfd.normE", pts_slice).reshape(NX_SLICE, NZ_SLICE)
    return x, y, z, e_3d, xs, zs, e_slice


def _render_iso(ax, x, y, z, e_3d, *, e_max, hw, L):
    dx, dy, dz = float(x[1] - x[0]), float(y[1] - y[0]), float(z[1] - z[0])
    levels = [0.6, 0.78, 0.92]
    colors = [OKABE_ITO[6], OKABE_ITO[5], OKABE_ITO[1]]
    alphas = [0.16, 0.28, 0.48]
    for level, col, alpha in zip(levels, colors, alphas):
        iso = level * e_max
        try:
            verts, faces, _, _ = measure.marching_cubes(e_3d, level=iso,
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
    for xx, yy in zip(xs[:-1], ys[:-1]):
        ax.plot([xx, xx], [yy, yy], [0, L], color="0.55", lw=0.5)
    z_n0 = GEOM.internode_length_um
    z_n1 = z_n0 + GEOM.node_length_um
    for r in (GEOM.axon_radius_um, GEOM.myelin_radius_um):
        for z_ring in (z_n0, z_n1):
            _ring(ax, r, z_ring, color=OKABE_ITO[6], lw=1.3, alpha=0.9)
    ax.set_box_aspect((1, 1, 2.0))
    ax.view_init(elev=20, azim=-58)


def main() -> None:
    apply_thesis_style(width_mm=DOUBLE_COL_MM, aspect=0.85)
    t0 = time.monotonic()
    client = mph.start()

    L = total_length_um(GEOM)
    hw = GEOM.external_half_width_um

    results = {}
    for label, f_hz in [("off", F_OFF), ("on", F_ON)]:
        model = client.create(f"sim13_{label}")
        build_geometry(model, GEOM)
        apply_materials(model, MAT)
        setup_physics(model, GEOM)
        build_mesh(model, MESH)
        setup_study(model, f_hz)
        print(f"[{label}] f = {f_hz/1e12:g} THz, solving...")
        ts = time.monotonic()
        solve_study(model)
        print(f"[{label}] solved in {time.monotonic() - ts:.1f}s, sampling...")
        results[label] = _sample(model, hw, L)
        client.remove(model)

    client.clear()

    # Shared color scale based on the on-resonance maximum
    e_max = max(results["on"][3].max(), results["off"][3].max())
    e_max_slice = max(results["on"][6].max(), results["off"][6].max())

    fig = plt.figure()
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.3], height_ratios=[1, 1])
    ax_off_3d = fig.add_subplot(gs[0, 0], projection="3d")
    ax_off_2d = fig.add_subplot(gs[0, 1])
    ax_on_3d  = fig.add_subplot(gs[1, 0], projection="3d")
    ax_on_2d  = fig.add_subplot(gs[1, 1])

    for label, ax_3d, ax_2d, freq_hz in [
        ("off", ax_off_3d, ax_off_2d, F_OFF),
        ("on",  ax_on_3d,  ax_on_2d,  F_ON),
    ]:
        x, y, z, e_3d, xs, zs, e_slice = results[label]
        _render_iso(ax_3d, x, y, z, e_3d, e_max=e_max, hw=hw, L=L)
        ax_3d.set_title(
            f"{'OFF' if label == 'off' else 'ON'}  f = {freq_hz/1e12:.3f} THz  "
            f"|E|_max = {e_3d.max():.2f}",
            fontsize=8,
        )
        # 2D slice
        e_slice_sm = gaussian_filter(e_slice, sigma=(1.0, 1.6))
        im = ax_2d.pcolormesh(zs, xs, e_slice_sm, cmap="magma",
                              vmin=0, vmax=e_max_slice, shading="auto")
        for r in (GEOM.axon_radius_um, GEOM.myelin_radius_um):
            ax_2d.axhline( r, color="white", lw=0.6, alpha=0.55)
            ax_2d.axhline(-r, color="white", lw=0.6, alpha=0.55)
        z_n0 = GEOM.internode_length_um
        z_n1 = z_n0 + GEOM.node_length_um
        ax_2d.axvline(z_n0, color=OKABE_ITO[6], lw=1.0, alpha=0.85)
        ax_2d.axvline(z_n1, color=OKABE_ITO[6], lw=1.0, alpha=0.85)
        ax_2d.set_xlabel("z (µm)")
        ax_2d.set_ylabel("x (µm)")
        ax_2d.set_title(
            f"y = 0 cross-section, {'OFF' if label == 'off' else 'ON'}  "
            f"({freq_hz/1e12:.3f} THz)", fontsize=8,
        )
        fig.colorbar(im, ax=ax_2d, shrink=0.85, pad=0.02).set_label("|E|")

    fig.suptitle(
        f"Sim 13 — Resonance contrast: off (f={F_OFF/1e12:g} THz) vs on (f={F_ON/1e12:g} THz)",
        y=0.995,
    )
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim13_resonance_contrast_3d")
    print(f"wrote {pdf}\nwrote {png}")
    print(f"TOTAL: {time.monotonic() - t0:.1f}s")


if __name__ == "__main__":
    main()
