"""Sim 12 — 3D stack of axial |E| profiles per node length.

Reads the sim 3 HDF5 results (one per node length) and stacks them
along a node-length axis as a 3D surface. Shows how the axial structure
of the resonant scattered field reorganises as the node grows from
10 µm to 100 µm.

A secondary plot extracts the peak-|E|-in-node trend (as in sim 3's
left panel) and overlays it on the 3D structure for context.

Pure Python from saved HDF5, no COMSOL.
"""

from __future__ import annotations

from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

from thznerve.plots.style import DOUBLE_COL_MM, OKABE_ITO, apply_thesis_style, save_figure

REPO_ROOT = Path(__file__).resolve().parent.parent
SIM3_DIR = REPO_ROOT / "results" / "sim3"


def main() -> None:
    apply_thesis_style(width_mm=DOUBLE_COL_MM, aspect=0.5)

    h5_files = sorted(SIM3_DIR.glob("nodeL_*.h5"))
    if not h5_files:
        raise SystemExit("no sim3 HDF5 files — run experiments/sim3_node_length.py")

    rows = []
    for p in h5_files:
        with h5py.File(p, "r") as f:
            node_L = float(f["scalars"].attrs["node_length_um"])
            peak = float(f["scalars"].attrs["peak_E_node"])
            z = f["axial_profile/z"][...]
            e = f["axial_profile/E_mag"][...]
            rows.append((node_L, z, e, peak))

    rows.sort(key=lambda r: r[0])
    node_Ls = np.array([r[0] for r in rows])
    peaks   = np.array([r[3] for r in rows])

    # All sims have the same n_z but different z extents (different total_L).
    # Normalise to z/total_L for a comparable axis.
    n_z = len(rows[0][1])
    z_norm = np.linspace(0, 1, n_z)
    Egrid = np.array([r[2] for r in rows])  # (n_L, n_z)

    fig = plt.figure()
    ax_3d = fig.add_subplot(1, 2, 1, projection="3d")
    ax_2d = fig.add_subplot(1, 2, 2)

    # 3-D ribbon — one polyline per node-length on its own y-slot
    for (node_L, _, e_profile, _), zorder in zip(rows, range(len(rows), 0, -1)):
        ax_3d.plot(z_norm, [node_L] * n_z, e_profile,
                   lw=1.4, alpha=0.85, label=fr"$L_{{node}} = {node_L:g}$ µm")
        # Fill under the curve for ribbon effect
        verts = [(z, node_L, 0) for z in z_norm]
        verts += [(z, node_L, e) for z, e in zip(z_norm[::-1], e_profile[::-1])]

    # Mark the node region (in normalised z, the node is centred at 0.5
    # with width node_L / total_L — varies per profile; we shade the
    # centre 0.45–0.55 band as a guide)
    for node_L, _, _, _ in rows:
        for z_marker in (0.45, 0.55):
            ax_3d.plot([z_marker, z_marker], [node_L, node_L], [0, 3.5],
                       color="0.7", lw=0.5, alpha=0.4)

    ax_3d.set_xlabel("z / total_L (normalised)")
    ax_3d.set_ylabel("node length (µm)")
    ax_3d.set_zlabel("|E|")
    ax_3d.set_title("Axial |E| profile per node length")
    ax_3d.legend(loc="upper left", fontsize=6.5, framealpha=0.92)
    ax_3d.view_init(elev=22, azim=-72)
    ax_3d.set_box_aspect((1.2, 1.0, 0.85))

    # 2-D heatmap projection
    im = ax_2d.pcolormesh(z_norm, node_Ls, Egrid, cmap="magma", shading="auto")
    ax_2d.axvspan(0.45, 0.55, color="white", alpha=0.08)
    ax_2d.axvline(0.5, color="white", lw=0.6, alpha=0.5, ls=":")
    ax_2d.text(0.5, node_Ls.min() + 4, "node centre",
               color="white", ha="center", fontsize=7, alpha=0.85)
    ax_2d.set_xlabel("z / total_L")
    ax_2d.set_ylabel("node length (µm)")
    ax_2d.set_title("|E|(z, $L_{node}$) heatmap")
    fig.colorbar(im, ax=ax_2d, shrink=0.85, pad=0.02).set_label("|E|")

    # Inset: peak vs node_L
    ax_inset = ax_2d.inset_axes([0.62, 0.65, 0.34, 0.32])
    ax_inset.plot(node_Ls, peaks, "o-", color=OKABE_ITO[6], lw=1.4, ms=4)
    ax_inset.set_xlabel("$L_{node}$ µm", fontsize=6.5)
    ax_inset.set_ylabel("peak |E|", fontsize=6.5)
    ax_inset.tick_params(labelsize=5.5)
    ax_inset.grid(True, alpha=0.3)
    ax_inset.set_facecolor("#ffffffcc")

    fig.suptitle("Sim 12 — Node-length sweep: axial structure in 3-D")
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim12_node_length_3d")
    print(f"wrote {pdf}\nwrote {png}")


if __name__ == "__main__":
    main()
