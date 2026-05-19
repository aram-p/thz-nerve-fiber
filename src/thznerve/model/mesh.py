"""Free-tetrahedral mesh for the nerve fiber model.

A single global FreeTet feature with a coarse default size (~30 µm) and a
finer Size sub-feature on the node selection (~5 µm). The mesh refinement
factor is exposed so a future mesh-convergence study (issue #14 follow-up)
can sweep it via the YAML config.
"""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, ConfigDict, Field


class MeshParams(BaseModel):
    model_config = ConfigDict(extra="forbid", protected_namespaces=())

    max_element_size_um: float = Field(default=30.0, gt=0)
    node_element_size_um: float = Field(default=5.0, gt=0)
    refinement_factor: float = Field(default=1.0, gt=0)


def _clear_mesh_features(mesh_java: Any) -> None:
    for tag in [str(t) for t in mesh_java.feature().tags()]:
        if tag == "size":  # auto-managed global Size
            continue
        mesh_java.feature().remove(tag)


def build_mesh(model: Any, params: MeshParams) -> int:
    """Build a free tetrahedral mesh with refinement at the node. Returns element count."""

    from thznerve.model.geometry import selection_tag

    java = model.java
    comp_tag = str(java.component().tags()[0])
    comp = java.component(comp_tag)
    geom_tag = str(comp.geom().tags()[0])

    if comp.mesh().tags().length == 0:
        comp.mesh().create("mesh1", geom_tag)
    mesh_tag = str(comp.mesh().tags()[0])
    mesh = comp.mesh(mesh_tag)
    _clear_mesh_features(mesh)

    coarse = params.max_element_size_um / params.refinement_factor
    fine = params.node_element_size_um / params.refinement_factor

    # Global size — sets the max element size everywhere.
    mesh.feature("size").set("custom", "on")
    mesh.feature("size").set("hmax", f"{coarse}[um]")
    mesh.feature("size").set("hmin", f"{fine / 5}[um]")

    # FreeTet over all domains.
    ftet = mesh.feature().create("ftet1", "FreeTet")
    ftet.label("Free Tetrahedral")

    # Size sub-feature: finer mesh in the node domain.
    size_node = ftet.feature().create("size_node", "Size")
    size_node.label("Node refinement")
    size_node.selection().geom(geom_tag, 3)
    size_node.selection().named(selection_tag("node"))
    size_node.set("custom", "on")
    size_node.set("hmax", f"{fine}[um]")
    size_node.set("hmaxactive", True)

    mesh.run()
    return int(mesh.getNumElem())
