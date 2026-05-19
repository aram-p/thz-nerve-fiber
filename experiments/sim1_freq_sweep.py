"""Sim 1 — Frequency sweep at baseline geometry, σ_node = 0.

26 frequencies linearly spaced 0.1 → 2.0 THz. At each frequency:

    * Set up EWFD scattered-field problem in COMSOL (periodic-cell unit).
    * Solve.
    * Sample |E| along the z-axis (x=y=0), 200 points.
    * Record peak |E| in the node window (z ∈ [internode_L, internode_L+node_L])
      and the global peak.

Output: a CSV of (frequency, peak_E_node, peak_E_global) under
results/sim1/, an HDF5 with the full axial profiles, and a thesis-styled
spectrum figure.

Why this matters
----------------
This is the central experimental observation we're trying to reproduce:
Hovhannisyan & Makaryan saw resonant THz absorption near 0.6 THz and
2 THz in spinal-cord samples (paper 1). If our single-fibre unit-cell
model shows |E|-at-the-node enhancement at corresponding frequencies,
that supports the "diffraction-grating-of-conductive-fibres"
interpretation. If not, it tells us the resonance mechanism is somewhere
else (geometry orientation, multi-fibre coupling, voltage-driven
nonlinearity).
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
OUT_DIR = REPO_ROOT / "results" / "sim1"

GEOM = GeometryParams(
    axon_radius_um=5,
    myelin_radius_um=7,
    node_length_um=40,
    internode_length_um=100,
    external_half_width_um=20,
)
MAT = MaterialParams(node_sigma_S_per_m=0.0)
MESH = MeshParams(max_element_size_um=30.0, node_element_size_um=5.0)

FREQ_THZ: np.ndarray = np.linspace(0.1, 2.0, 26)


def _solve_one(model, freq_hz: float) -> tuple[np.ndarray, np.ndarray]:
    """Update study frequency, solve, extract axial |E| profile. Returns (z, |E|)."""

    java = model.java
    study_tag = str(java.study().tags()[0])
    java.study(study_tag).feature("freq").set("plist", str(freq_hz))
    java.study(study_tag).run()
    return extract_axial_profile(model, GEOM, n_points=200)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    apply_thesis_style(width_mm=SINGLE_COL_MM, aspect=0.7)

    t0 = time.monotonic()
    client = mph.start()
    model = client.create("sim1_freq_sweep")

    build_geometry(model, GEOM)
    apply_materials(model, MAT)
    setup_physics(model, GEOM)
    n_elem = build_mesh(model, MESH)
    setup_study(model, FREQ_THZ[0] * 1e12)
    print(f"setup: {n_elem} elements, {time.monotonic() - t0:.1f}s")

    git_sha = get_git_sha()
    cfg_hash = config_hash({
        "sim": "sim1_freq_sweep",
        "geom": GEOM.model_dump(),
        "mat": MAT.model_dump(),
        "mesh": MESH.model_dump(),
    })

    rows = []
    L = total_length_um(GEOM)
    z_node_lo = GEOM.internode_length_um
    z_node_hi = z_node_lo + GEOM.node_length_um

    for i, f_thz in enumerate(FREQ_THZ, start=1):
        f_hz = float(f_thz) * 1e12
        t_solve = time.monotonic()
        z, e = _solve_one(model, f_hz)
        dt = time.monotonic() - t_solve

        node_mask = (z >= z_node_lo) & (z <= z_node_hi)
        peak_node = float(np.max(e[node_mask])) if node_mask.any() else float("nan")
        peak_global = float(np.max(e))
        mean_node = float(np.mean(e[node_mask])) if node_mask.any() else float("nan")

        h5_path = OUT_DIR / f"freq_{i:03d}_{f_thz:.3f}THz.h5"
        write_result(
            h5_path,
            scalars={
                "frequency_THz": f_thz, "peak_E_node": peak_node,
                "peak_E_global": peak_global, "mean_E_node": mean_node,
                "solve_seconds": dt,
            },
            axial_profile=(z, e),
            axial_slice=(np.array([0.0]), z, e[None, :]),
            metadata={
                "git_sha": git_sha, "config_hash": cfg_hash,
                "sim": "sim1_freq_sweep",
                "params": {**GEOM.model_dump(), **MAT.model_dump(), "freq_THz": f_thz},
                "comsol_version": "6.3",
            },
        )

        rows.append((f_thz, peak_node, peak_global, mean_node, dt))
        print(f"  [{i:>2}/{len(FREQ_THZ)}] f = {f_thz:5.3f} THz  "
              f"|E|_peak_node = {peak_node:6.3f}  |E|_peak = {peak_global:6.3f}  "
              f"({dt:4.1f}s)")

    client.clear()

    # Write CSV summary
    csv_path = OUT_DIR / "spectrum.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["frequency_THz", "peak_E_node", "peak_E_global", "mean_E_node", "solve_s"])
        w.writerows(rows)

    # Plot
    arr = np.array(rows)
    f_thz_arr = arr[:, 0]
    peak_node = arr[:, 1]
    peak_global = arr[:, 2]
    mean_node = arr[:, 3]

    fig, ax = plt.subplots()
    ax.plot(f_thz_arr, peak_node, "o-", lw=1.3, ms=3.5, label="peak |E| in node")
    ax.plot(f_thz_arr, peak_global, "s--", lw=1.0, ms=2.5, label="peak |E| (global)", alpha=0.65)
    ax.plot(f_thz_arr, mean_node, "^:", lw=1.0, ms=2.5, label="mean |E| in node", alpha=0.65)
    ax.axvline(0.6, color="grey", ls=":", lw=0.7, alpha=0.7)
    ax.axvline(2.0, color="grey", ls=":", lw=0.7, alpha=0.7)
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel("|E| (normalised, scattered field)")
    ax.set_title(f"Sim 1 — Frequency sweep (σ_node = {MAT.node_sigma_S_per_m:g})")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", framealpha=0.9)
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim1_freq_sweep")
    print(f"wrote {pdf}\nwrote {png}\nwrote {csv_path}")
    print(f"TOTAL: {time.monotonic() - t0:.1f}s")


if __name__ == "__main__":
    main()
