"""Sim 3 — Node-of-Ranvier length sensitivity at f = 0.6 THz, σ = 0.

Real biology has node lengths ~1 µm; the modelled value (40 µm) is
order-of-magnitude larger because finer features would explode the
mesh. The question this sim asks: how strongly does the |E|-at-node
response depend on node length? If the field-at-node is roughly flat
in node_L, the absolute number isn't critical and the wire-array /
diffraction-grating interpretation is robust. If it's a sharp
function of node_L, mesh-and-geometry refinement becomes essential
before any quantitative claim.
"""

from __future__ import annotations

import csv
import time
from pathlib import Path

import matplotlib.pyplot as plt
import mph
import numpy as np

from thznerve.io.hdf5 import write_result
from thznerve.io.provenance import config_hash, get_git_sha
from thznerve.model.geometry import GeometryParams, build_geometry, total_length_um
from thznerve.model.materials import MaterialParams, apply_materials
from thznerve.model.mesh import MeshParams, build_mesh
from thznerve.model.study import (
    extract_axial_profile,
    setup_physics,
    setup_study,
    solve_study,
)
from thznerve.plots.style import SINGLE_COL_MM, apply_thesis_style, save_figure

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "results" / "sim3"

FREQ_HZ = 0.6e12
NODE_LENGTHS_UM: list[float] = [10.0, 20.0, 40.0, 60.0, 100.0]
MAT = MaterialParams(node_sigma_S_per_m=0.0)


def _geom(node_L: float) -> GeometryParams:
    return GeometryParams(
        axon_radius_um=5,
        myelin_radius_um=7,
        node_length_um=node_L,
        internode_length_um=100,
        external_half_width_um=20,
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    apply_thesis_style(width_mm=SINGLE_COL_MM, aspect=0.75)

    t0 = time.monotonic()
    client = mph.start()

    git_sha = get_git_sha()
    rows = []
    profiles = []

    for i, node_L in enumerate(NODE_LENGTHS_UM, start=1):
        geom = _geom(node_L)
        model = client.create(f"sim3_nodeL_{node_L:g}")
        build_geometry(model, geom)
        apply_materials(model, MAT)
        setup_physics(model, geom)
        n_elem = build_mesh(model, MeshParams(max_element_size_um=30.0,
                                              node_element_size_um=min(5.0, node_L / 4)))
        setup_study(model, FREQ_HZ)
        t_solve = time.monotonic()
        solve_study(model)
        z, e = extract_axial_profile(model, geom, n_points=200)
        dt = time.monotonic() - t_solve

        L = total_length_um(geom)
        z_node_lo = geom.internode_length_um
        z_node_hi = z_node_lo + geom.node_length_um
        node_mask = (z >= z_node_lo) & (z <= z_node_hi)
        peak_node = float(np.max(e[node_mask])) if node_mask.any() else float("nan")
        peak_global = float(np.max(e))
        mean_node = float(np.mean(e[node_mask])) if node_mask.any() else float("nan")

        h5_path = OUT_DIR / f"nodeL_{i:02d}_{node_L:g}um.h5"
        write_result(
            h5_path,
            scalars={
                "node_length_um": node_L, "peak_E_node": peak_node,
                "peak_E_global": peak_global, "mean_E_node": mean_node,
                "n_elem": n_elem, "solve_seconds": dt,
            },
            axial_profile=(z, e),
            axial_slice=(np.array([0.0]), z, e[None, :]),
            metadata={
                "git_sha": git_sha,
                "sim": "sim3_node_length",
                "params": {**geom.model_dump(), **MAT.model_dump(),
                           "freq_THz": FREQ_HZ / 1e12},
                "comsol_version": "6.3",
            },
        )
        rows.append((node_L, peak_node, peak_global, mean_node, dt, n_elem))
        profiles.append((node_L, z, e, z_node_lo, z_node_hi))
        print(f"  [{i}/{len(NODE_LENGTHS_UM)}] node_L = {node_L:5.1f} µm   "
              f"elem = {n_elem:6d}   |E|_node = {peak_node:5.3f}   "
              f"({dt:4.1f}s)")
        # Clean up the model so client doesn't accumulate (we're creating fresh per L)
        client.remove(model)

    client.clear()

    csv_path = OUT_DIR / "node_length.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["node_length_um", "peak_E_node", "peak_E_global", "mean_E_node",
                    "solve_s", "n_elem"])
        w.writerows(rows)

    arr = np.array(rows)
    nL = arr[:, 0]
    peak_node_arr = arr[:, 1]
    mean_node_arr = arr[:, 3]

    fig, (ax_scan, ax_prof) = plt.subplots(1, 2, figsize=(8, 3))

    ax_scan.plot(nL, peak_node_arr, "o-", lw=1.4, ms=4, label="peak |E| in node")
    ax_scan.plot(nL, mean_node_arr, "^:", lw=1.0, ms=3, label="mean |E| in node", alpha=0.75)
    ax_scan.set_xlabel("Node length (µm)")
    ax_scan.set_ylabel("|E| (normalised)")
    ax_scan.set_title(f"|E| vs node length at f = {FREQ_HZ / 1e12:.2f} THz")
    ax_scan.legend(loc="best", framealpha=0.9)
    ax_scan.grid(True, alpha=0.3)

    for node_L, z, e, z_lo, z_hi in profiles:
        ax_prof.plot(z, e, lw=1.2, label=f"node L = {node_L:g} µm")
    # Shade the typical node region for the 40 µm reference geometry
    ax_prof.axvspan(100, 140, color="grey", alpha=0.12)
    ax_prof.set_xlabel("z (µm)")
    ax_prof.set_ylabel("|E| (normalised)")
    ax_prof.set_title("Axial profile per node length")
    ax_prof.legend(loc="best", framealpha=0.9, fontsize=7)
    ax_prof.grid(True, alpha=0.3)

    fig.suptitle("Sim 3 — Node-length sensitivity")
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim3_node_length")
    print(f"wrote {pdf}\nwrote {png}\nwrote {csv_path}")
    print(f"TOTAL: {time.monotonic() - t0:.1f}s")


if __name__ == "__main__":
    main()
