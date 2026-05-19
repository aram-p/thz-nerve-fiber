"""Sim 2 — Node-conductivity (σ) sweep at f = 0.6 THz.

For each σ in a logarithmically-spaced set, solve the EWFD problem and
record peak |E| at the node and in the external region. This directly
parameterises Hovhannisyan & Makaryan paper 3's "open ion channels →
nonlinear conductivity" mechanism: above some threshold voltage, σ at
the node grows from ~0 to finite values; the question is how the
THz field responds to that conductivity change at the experimentally
observed resonance frequency.
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
OUT_DIR = REPO_ROOT / "results" / "sim2"

GEOM = GeometryParams(
    axon_radius_um=5,
    myelin_radius_um=7,
    node_length_um=40,
    internode_length_um=100,
    external_half_width_um=20,
)
MESH = MeshParams(max_element_size_um=30.0, node_element_size_um=5.0)

FREQ_HZ: float = 0.6e12
SIGMA_VALUES: list[float] = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    apply_thesis_style(width_mm=SINGLE_COL_MM, aspect=0.75)

    t0 = time.monotonic()
    client = mph.start()
    model = client.create("sim2_sigma_sweep")

    build_geometry(model, GEOM)
    setup_physics(model, GEOM)
    n_elem = build_mesh(model, MESH)
    setup_study(model, FREQ_HZ)
    print(f"setup: {n_elem} elements, {time.monotonic() - t0:.1f}s")

    git_sha = get_git_sha()
    rows = []
    L = total_length_um(GEOM)
    z_node_lo = GEOM.internode_length_um
    z_node_hi = z_node_lo + GEOM.node_length_um

    for i, sigma in enumerate(SIGMA_VALUES, start=1):
        mat = MaterialParams(node_sigma_S_per_m=sigma)
        apply_materials(model, mat)
        t_solve = time.monotonic()
        # Need to rerun study because materials changed via the global parameter.
        model.java.study(str(model.java.study().tags()[0])).run()
        z, e = extract_axial_profile(model, GEOM, n_points=200)
        dt = time.monotonic() - t_solve

        node_mask = (z >= z_node_lo) & (z <= z_node_hi)
        peak_node = float(np.max(e[node_mask])) if node_mask.any() else float("nan")
        peak_global = float(np.max(e))
        mean_node = float(np.mean(e[node_mask])) if node_mask.any() else float("nan")

        h5_path = OUT_DIR / f"sigma_{i:02d}_{sigma:g}.h5"
        write_result(
            h5_path,
            scalars={
                "sigma_S_per_m": sigma, "peak_E_node": peak_node,
                "peak_E_global": peak_global, "mean_E_node": mean_node,
                "solve_seconds": dt, "frequency_THz": FREQ_HZ / 1e12,
            },
            axial_profile=(z, e),
            axial_slice=(np.array([0.0]), z, e[None, :]),
            metadata={
                "git_sha": git_sha,
                "sim": "sim2_sigma_sweep",
                "params": {**GEOM.model_dump(), **mat.model_dump(),
                           "freq_THz": FREQ_HZ / 1e12},
                "comsol_version": "6.3",
            },
        )
        rows.append((sigma, peak_node, peak_global, mean_node, dt))
        print(f"  [{i}/{len(SIGMA_VALUES)}] σ = {sigma:5.2f} S/m   "
              f"|E|_node = {peak_node:5.3f}   |E|_global = {peak_global:5.3f}   "
              f"({dt:4.1f}s)")

    client.clear()

    csv_path = OUT_DIR / "sigma_curve.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sigma_S_per_m", "peak_E_node", "peak_E_global", "mean_E_node", "solve_s"])
        w.writerows(rows)

    arr = np.array(rows)
    sig = arr[:, 0]
    peak_node_arr = arr[:, 1]
    peak_global_arr = arr[:, 2]
    mean_node_arr = arr[:, 3]

    fig, ax = plt.subplots()
    # Add small floor to avoid log(0) issues if we ever plot logarithmic
    ax.plot(sig, peak_node_arr, "o-", lw=1.4, ms=4, label="peak |E| in node")
    ax.plot(sig, peak_global_arr, "s--", lw=1.0, ms=3, label="peak |E| (global)", alpha=0.7)
    ax.plot(sig, mean_node_arr, "^:", lw=1.0, ms=3, label="mean |E| in node", alpha=0.7)
    ax.set_xlabel(r"Node conductivity $\sigma$ (S/m)")
    ax.set_ylabel("|E| (normalised, scattered field)")
    ax.set_title(f"Sim 2 — σ sweep at f = {FREQ_HZ / 1e12:.2f} THz")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", framealpha=0.9)
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim2_sigma_sweep")
    print(f"wrote {pdf}\nwrote {png}\nwrote {csv_path}")
    print(f"TOTAL: {time.monotonic() - t0:.1f}s")


if __name__ == "__main__":
    main()
