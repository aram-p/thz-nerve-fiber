"""End-to-end smoke test: build geometry + materials + mesh + study, solve at 0.6 THz.

Run:
    uv run python scripts/study_smoke.py
"""

from __future__ import annotations

import time
from pathlib import Path

import mph
import numpy as np

from thznerve.model.geometry import GeometryParams, build_geometry
from thznerve.model.materials import MaterialParams, apply_materials
from thznerve.model.mesh import MeshParams, build_mesh
from thznerve.model.study import (
    extract_axial_profile,
    extract_axial_slice,
    setup_physics,
    setup_study,
    solve_study,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
ARTIFACT = REPO_ROOT / "results" / "study_smoke.mph"

GEOM = GeometryParams(
    axon_radius_um=5,
    myelin_radius_um=7,
    node_length_um=40,
    internode_length_um=100,
    external_half_width_um=20,
)
MAT = MaterialParams(node_sigma_S_per_m=0.0)
MESH = MeshParams(max_element_size_um=30.0, node_element_size_um=5.0)
FREQ_HZ = 0.6e12


def main() -> None:
    t0 = time.monotonic()

    client = mph.start()
    model = client.create("study_smoke")

    print("[1/5] geometry...")
    n_dom = build_geometry(model, GEOM)
    print(f"      {n_dom} domains")
    assert n_dom == 5

    print("[2/5] materials...")
    apply_materials(model, MAT)

    print("[3/5] physics + mesh...")
    setup_physics(model, GEOM)
    n_elem = build_mesh(model, MESH)
    print(f"      {n_elem} elements")

    print(f"[4/5] solve at {FREQ_HZ:g} Hz...")
    setup_study(model, FREQ_HZ)
    t_solve = time.monotonic()
    solve_study(model)
    print(f"      solved in {time.monotonic() - t_solve:.1f}s")

    print("[5/5] extract axial profile + slice...")
    z, e_axis = extract_axial_profile(model, GEOM, n_points=200)
    x, z_slice, e_slice = extract_axial_slice(model, GEOM, nx=60, nz=120)
    print(f"      |E| axis: min={e_axis.min():.3g}  max={e_axis.max():.3g}  mean={e_axis.mean():.3g}")
    print(f"      |E| slice shape: {e_slice.shape}")
    assert np.isfinite(e_axis).all(), "non-finite values in axial profile"
    assert e_axis.max() > 0, "axial profile is all zero — solve produced no field"

    ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
    model.save(ARTIFACT)
    print(f"saved artifact: {ARTIFACT}")
    print(f"TOTAL: {time.monotonic() - t0:.1f}s")
    client.clear()


if __name__ == "__main__":
    main()
