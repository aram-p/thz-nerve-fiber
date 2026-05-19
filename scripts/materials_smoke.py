"""Smoke test for `thznerve.model.materials.apply_materials`.

Builds the geometry, applies the three materials by named selection,
and saves the .mph for visual inspection.

Run:
    uv run python scripts/materials_smoke.py
"""

from __future__ import annotations

from pathlib import Path

import mph

from thznerve.model.geometry import GeometryParams, build_geometry
from thznerve.model.materials import MaterialParams, apply_materials

REPO_ROOT = Path(__file__).resolve().parent.parent
ARTIFACT = REPO_ROOT / "results" / "materials_smoke.mph"


def main() -> None:
    client = mph.start()
    model = client.create("materials_smoke")

    geom_params = GeometryParams(
        axon_radius_um=5,
        myelin_radius_um=7,
        node_length_um=40,
        internode_length_um=100,
        external_radius_um=20,
    )
    n_domains = build_geometry(model, geom_params)
    print(f"geometry: {n_domains} domains")
    assert n_domains == 5

    mat_params = MaterialParams(node_sigma_S_per_m=1.0)
    assignments = apply_materials(model, mat_params)
    for tag, label in assignments.items():
        print(f"  {tag:35s} -> {label}")

    ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
    model.save(ARTIFACT)
    print(f"saved artifact: {ARTIFACT}")
    print("OK: materials smoke test complete.")

    client.clear()


if __name__ == "__main__":
    main()
