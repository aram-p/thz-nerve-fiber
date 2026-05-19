"""Smoke test for `thznerve.model.geometry.build_geometry`.

Builds the geometry for several parameter sets and asserts the resulting
domain count is 5. Saves the first model as a .mph artifact for visual
inspection.

Run:
    uv run python scripts/geometry_smoke.py
"""

from __future__ import annotations

from pathlib import Path

import mph

from thznerve.model.geometry import GeometryParams, build_geometry

PARAM_SETS: list[GeometryParams] = [
    GeometryParams(
        axon_radius_um=5,
        myelin_radius_um=7,
        node_length_um=40,
        internode_length_um=100,
        external_radius_um=20,
    ),
    GeometryParams(
        axon_radius_um=3,
        myelin_radius_um=5,
        node_length_um=30,
        internode_length_um=80,
        external_radius_um=15,
    ),
    GeometryParams(
        axon_radius_um=8,
        myelin_radius_um=10,
        node_length_um=60,
        internode_length_um=120,
        external_radius_um=30,
    ),
]

REPO_ROOT = Path(__file__).resolve().parent.parent
ARTIFACT_PATH = REPO_ROOT / "results" / "geometry_smoke.mph"


def main() -> None:
    client = mph.start()
    expected_domains = 5

    for i, params in enumerate(PARAM_SETS, start=1):
        model = client.create(f"geom_smoke_{i}")
        n = build_geometry(model, params)
        print(
            f"[set {i}] axon_r={params.axon_radius_um} myelin_r={params.myelin_radius_um} "
            f"node_L={params.node_length_um} internode_L={params.internode_length_um} "
            f"ext_r={params.external_radius_um}  -->  {n} domains"
        )
        assert n == expected_domains, (
            f"set {i}: expected {expected_domains} domains, got {n}"
        )
        if i == 1:
            ARTIFACT_PATH.parent.mkdir(parents=True, exist_ok=True)
            model.save(ARTIFACT_PATH)
            print(f"  saved artifact: {ARTIFACT_PATH}")

    print("OK: geometry smoke test passed for all parameter sets.")
    client.clear()


if __name__ == "__main__":
    main()
