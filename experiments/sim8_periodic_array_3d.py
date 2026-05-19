"""Sim 8 — 3D rendering of the periodic fibre array (diffraction grating).

The single-fibre unit cell sim (sim 5) doesn't communicate the
"diffraction grating of conductive fibres" framing from Hovhannisyan
& Makaryan paper 3. This sim tiles the unit cell into a 3 x 3 (x, y)
periodic array so the array geometry — fibre spacing, node alignment,
incident wave from below — is visible at a glance.

Pure matplotlib 3D, no COMSOL.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from thznerve.model.geometry import GeometryParams, total_length_um
from thznerve.plots.style import DOUBLE_COL_MM, OKABE_ITO, apply_thesis_style, save_figure


def _cyl_surface(radius, z0, z1, *, xc=0.0, yc=0.0, n_theta=36):
    theta = np.linspace(0, 2 * np.pi, n_theta)
    z = np.linspace(z0, z1, 2)
    th, zz = np.meshgrid(theta, z)
    return xc + radius * np.cos(th), yc + radius * np.sin(th), zz


def _cyl_endcap(radius, z, *, xc=0.0, yc=0.0, n_theta=36):
    theta = np.linspace(0, 2 * np.pi, n_theta)
    return xc + radius * np.cos(theta), yc + radius * np.sin(theta), np.full_like(theta, z)


def _add_cyl(ax, *, xc, yc, r, z0, z1, color, alpha):
    X, Y, Z = _cyl_surface(r, z0, z1, xc=xc, yc=yc)
    ax.plot_surface(X, Y, Z, color=color, alpha=alpha, linewidth=0)
    for z_cap in (z0, z1):
        xs, ys, zs = _cyl_endcap(r, z_cap, xc=xc, yc=yc)
        verts = [list(zip(xs, ys, zs))]
        ax.add_collection3d(
            Poly3DCollection(verts, facecolors=color, alpha=alpha,
                             edgecolor="0.3", linewidths=0.25)
        )


def _add_unit_cell(ax, *, p, xc, yc, faded=False):
    L = total_length_um(p)
    z_n0 = p.internode_length_um
    z_n1 = z_n0 + p.node_length_um
    a = 0.65 if not faded else 0.18
    # myelin (proximal)
    _add_cyl(ax, xc=xc, yc=yc, r=p.myelin_radius_um, z0=0,   z1=z_n0,
             color=OKABE_ITO[1], alpha=0.30 * (a / 0.65))
    # node
    _add_cyl(ax, xc=xc, yc=yc, r=p.myelin_radius_um, z0=z_n0, z1=z_n1,
             color=OKABE_ITO[6], alpha=0.55 * (a / 0.65))
    # myelin (distal)
    _add_cyl(ax, xc=xc, yc=yc, r=p.myelin_radius_um, z0=z_n1, z1=L,
             color=OKABE_ITO[1], alpha=0.30 * (a / 0.65))
    # axon
    _add_cyl(ax, xc=xc, yc=yc, r=p.axon_radius_um, z0=0, z1=L,
             color=OKABE_ITO[2], alpha=0.7 * (a / 0.65))


def _draw_cell_outline(ax, *, xc, yc, hw, L, color="0.45", lw=0.55):
    xs = [xc - hw, xc + hw, xc + hw, xc - hw, xc - hw]
    ys = [yc - hw, yc - hw, yc + hw, yc + hw, yc - hw]
    for z in (0, L):
        ax.plot(xs, ys, [z] * 5, color=color, lw=lw)
    for x, y in zip(xs[:-1], ys[:-1]):
        ax.plot([x, x], [y, y], [0, L], color=color, lw=lw)


def main() -> None:
    apply_thesis_style(width_mm=DOUBLE_COL_MM, aspect=0.55)

    p = GeometryParams(
        axon_radius_um=5,
        myelin_radius_um=7,
        node_length_um=40,
        internode_length_um=100,
        external_half_width_um=20,
    )
    L = total_length_um(p)
    hw = p.external_half_width_um
    cell_pitch = 2 * hw  # 40 µm between unit-cell centres

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # 3 × 3 grid centred on origin. The central cell is the "highlighted"
    # one; surrounding cells are faded to convey periodicity.
    for ix in (-1, 0, 1):
        for iy in (-1, 0, 1):
            xc = ix * cell_pitch
            yc = iy * cell_pitch
            faded = not (ix == 0 and iy == 0)
            _add_unit_cell(ax, p=p, xc=xc, yc=yc, faded=faded)
            _draw_cell_outline(ax, xc=xc, yc=yc, hw=hw, L=L,
                               color=("0.25" if not faded else "0.65"),
                               lw=0.7 if not faded else 0.4)

    # Plane wave: k along z (vertical arrow under the array), E along x.
    extent = 1.5 * cell_pitch
    arrow_z0 = -55.0
    arrow_z1 = -15.0
    ax.quiver(extent, -extent, arrow_z0,
              0, 0, arrow_z1 - arrow_z0,
              color="0.1", lw=2.0, arrow_length_ratio=0.2)
    ax.text(extent + 1, -extent + 2, arrow_z1 + 4,
            r"$\mathbf{k}$", fontsize=10, color="0.1")
    ax.quiver(extent, -extent, arrow_z0,
              20, 0, 0,
              color=OKABE_ITO[6], lw=2.0, arrow_length_ratio=0.3)
    ax.text(extent + 22, -extent, arrow_z0 + 2,
            r"$\mathbf{E}_0\,\hat{x}$", fontsize=10, color=OKABE_ITO[6])

    # Legend proxies
    ax.plot([], [], color=OKABE_ITO[2], lw=8, alpha=0.7, label="Axon (water)")
    ax.plot([], [], color=OKABE_ITO[1], lw=8, alpha=0.5, label="Myelin sheath")
    ax.plot([], [], color=OKABE_ITO[6], lw=8, alpha=0.7, label="Node of Ranvier")
    ax.plot([], [], color="0.45", lw=1.5, label="Unit-cell box")

    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")
    ax.set_zlabel("z (µm)")
    ax.set_title(
        "Sim 8 — Periodic 3×3 fibre array (diffraction-grating unit cells)"
    )
    ax.set_box_aspect((1, 1, 1.6))
    ax.view_init(elev=18, azim=-58)
    ax.legend(loc="upper left", framealpha=0.92, fontsize=7)

    fig.tight_layout()
    pdf, png = save_figure(fig, "sim8_periodic_array_3d")
    print(f"wrote {pdf}\nwrote {png}")


if __name__ == "__main__":
    main()
