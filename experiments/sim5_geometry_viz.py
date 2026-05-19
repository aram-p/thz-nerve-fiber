"""Sim 5 — 3D visualisation of the unit-cell geometry.

Render the rectangular periodic unit cell with the five labelled
domains (axon, myelin proximal, node of Ranvier, myelin distal, external
water) and the incident-wave geometry (k along z, E along x). This is
the "what is the simulation actually modelling" figure for the tutor.

Pure matplotlib, no COMSOL — runs in < 2 s.

If COMSOL field data is available later, a follow-up sim renders the
|E| heatmap on top of this geometry; for now this is the static
"setup" figure.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from thznerve.model.geometry import GeometryParams, total_length_um
from thznerve.plots.style import DOUBLE_COL_MM, OKABE_ITO, apply_thesis_style, save_figure


def _cylinder_surface(radius: float, z0: float, z1: float, n_theta: int = 48):
    theta = np.linspace(0, 2 * np.pi, n_theta)
    z = np.linspace(z0, z1, 2)
    th, zz = np.meshgrid(theta, z)
    return radius * np.cos(th), radius * np.sin(th), zz


def _cylinder_endcaps(radius: float, z: float, n_theta: int = 48):
    theta = np.linspace(0, 2 * np.pi, n_theta)
    return radius * np.cos(theta), radius * np.sin(theta), np.full_like(theta, z)


def _draw_cylinder(ax, *, r, z0, z1, color, alpha=0.45, label=None, edgecolor="k"):
    X, Y, Z = _cylinder_surface(r, z0, z1)
    ax.plot_surface(X, Y, Z, color=color, alpha=alpha, linewidth=0, antialiased=True)
    for z_cap in (z0, z1):
        xs, ys, zs = _cylinder_endcaps(r, z_cap)
        verts = [list(zip(xs, ys, zs))]
        ax.add_collection3d(
            Poly3DCollection(verts, facecolors=color, alpha=alpha, edgecolor=edgecolor, linewidths=0.4)
        )
    if label is not None:
        ax.plot([], [], [], color=color, lw=8, alpha=0.6, label=label)


def _draw_box_wireframe(ax, *, hw, L, color="k", lw=0.7):
    xs = [-hw, hw, hw, -hw, -hw]
    ys = [-hw, -hw, hw, hw, -hw]
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
    z_node0 = p.internode_length_um
    z_node1 = z_node0 + p.node_length_um

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # Axon — solid through full length
    _draw_cylinder(ax, r=p.axon_radius_um, z0=0, z1=L,
                   color=OKABE_ITO[2], alpha=0.6, label="Axon")

    # Myelin proximal & distal (annular — drawn as outer cylinder)
    _draw_cylinder(ax, r=p.myelin_radius_um, z0=0, z1=z_node0,
                   color=OKABE_ITO[1], alpha=0.35, label="Myelin")
    _draw_cylinder(ax, r=p.myelin_radius_um, z0=z_node1, z1=L,
                   color=OKABE_ITO[1], alpha=0.35)

    # Node of Ranvier
    _draw_cylinder(ax, r=p.myelin_radius_um, z0=z_node0, z1=z_node1,
                   color=OKABE_ITO[6], alpha=0.55, label="Node of Ranvier")

    # External box outline
    _draw_box_wireframe(ax, hw=p.external_half_width_um, L=L, color="0.4")
    ax.plot([], [], [], color="0.4", lw=1.5, label="External water (box)")

    # Incident wave arrow (k along z), E vector at start (along x)
    arrow_z0 = -25.0
    arrow_z1 = -5.0
    ax.quiver(
        p.external_half_width_um * 0.85, p.external_half_width_um * 0.85, arrow_z0,
        0, 0, arrow_z1 - arrow_z0,
        color="0.1", lw=2.0, arrow_length_ratio=0.25,
    )
    ax.text(
        p.external_half_width_um * 0.9, p.external_half_width_um * 0.9, arrow_z0 - 5,
        "k", fontsize=9, color="0.1",
    )
    ax.quiver(
        p.external_half_width_um * 0.85, p.external_half_width_um * 0.85, arrow_z0,
        15, 0, 0,
        color=OKABE_ITO[6], lw=2.0, arrow_length_ratio=0.3,
    )
    ax.text(
        p.external_half_width_um * 0.85 + 16, p.external_half_width_um * 0.85, arrow_z0,
        r"$E_0\,\hat{x}$", fontsize=9, color=OKABE_ITO[6],
    )

    ax.set_xlabel("x (µm)", labelpad=-2)
    ax.set_ylabel("y (µm)", labelpad=-2)
    ax.set_zlabel("z (µm)", labelpad=-2)
    fig.suptitle(
        "Sim 5 — Single-fibre periodic unit cell\n"
        f"axon r = {p.axon_radius_um:g} µm, myelin r = {p.myelin_radius_um:g} µm, "
        f"node L = {p.node_length_um:g} µm, total L = {L:g} µm, "
        f"box ± {p.external_half_width_um:g} µm",
        y=0.98, fontsize=9,
    )
    ax.set_box_aspect((1, 1, 1.8))
    ax.view_init(elev=18, azim=-58)
    ax.legend(loc="upper left", framealpha=0.92, fontsize=7,
              bbox_to_anchor=(-0.05, 0.95))

    fig.tight_layout()
    pdf, png = save_figure(fig, "sim5_geometry_viz")
    print(f"wrote {pdf}\nwrote {png}")


if __name__ == "__main__":
    main()
