"""3D cylindrical nerve-fiber geometry builder.

Produces exactly five domains after Form Union:

    label             approx. region
    ---------------   ----------------------------------------------------
    axon              r ∈ [0, axon_r],          z ∈ [0, total_L]
    myelin_proximal   r ∈ [axon_r, myelin_r],   z ∈ [0, internode_L]
    node              r ∈ [axon_r, myelin_r],   z ∈ [internode_L, internode_L+node_L]
    myelin_distal     r ∈ [axon_r, myelin_r],   z ∈ [internode_L+node_L, total_L]
    external          r ∈ [myelin_r, ext_r],    z ∈ [0, total_L]

with ``total_L = 2*internode_L + node_L``.

Construction strategy: the four annular regions (3 sheath segments + external)
are built as ``Difference(outer_cyl, inner_cyl)`` so their top/bottom faces
only exist within their own radial band — that keeps the inner axon cylinder
from being chopped into z-segments by sheath face cuts.
"""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, ConfigDict, Field

DOMAIN_LABELS: tuple[str, ...] = (
    "axon",
    "myelin_proximal",
    "node",
    "myelin_distal",
    "external",
)


class GeometryParams(BaseModel):
    model_config = ConfigDict(extra="forbid", protected_namespaces=())

    axon_radius_um: float = Field(gt=0)
    myelin_radius_um: float = Field(gt=0)
    node_length_um: float = Field(gt=0)
    internode_length_um: float = Field(gt=0)
    external_radius_um: float = Field(gt=0)


def total_length_um(params: GeometryParams) -> float:
    return 2 * params.internode_length_um + params.node_length_um


def domain_centers(params: GeometryParams) -> dict[str, tuple[float, float, float]]:
    """One (x, y, z) point in µm inside each labeled domain.

    Phase 2.2 (materials) uses these to look up COMSOL domain indices per label.
    """

    r_axon = params.axon_radius_um
    r_myelin = params.myelin_radius_um
    r_ext = params.external_radius_um
    L_inter = params.internode_length_um
    L_node = params.node_length_um
    L_total = total_length_um(params)

    r_sheath_mid = (r_axon + r_myelin) / 2
    r_ext_mid = (r_myelin + r_ext) / 2

    return {
        "axon": (0.0, 0.0, L_total / 2),
        "myelin_proximal": (r_sheath_mid, 0.0, L_inter / 2),
        "node": (r_sheath_mid, 0.0, L_inter + L_node / 2),
        "myelin_distal": (r_sheath_mid, 0.0, L_inter + L_node + L_inter / 2),
        "external": (r_ext_mid, 0.0, L_total / 2),
    }


def _ensure_component(model: Any) -> Any:
    java = model.java
    if java.component().tags().length == 0:
        java.component().create("comp1", True)
    comp_tag = str(java.component().tags()[0])
    return java.component(comp_tag)


def _ensure_geometry(component_java: Any) -> Any:
    if component_java.geom().tags().length == 0:
        component_java.geom().create("geom1", 3)
    geom_tag = str(component_java.geom().tags()[0])
    geom = component_java.geom(geom_tag)
    geom.lengthUnit("um")
    return geom


def _clear_geometry_features(geom_java: Any) -> None:
    """Remove every feature except the auto-managed 'fin' (Form Union) finalize step."""

    feats = geom_java.feature()
    for tag in [str(t) for t in feats.tags()]:
        if tag == "fin":
            continue
        feats.remove(tag)


def _add_cyl(geom_java: Any, tag: str, *, radius: float, height: float, z_pos: float) -> None:
    cyl = geom_java.feature().create(tag, "Cylinder")
    cyl.set("r", str(radius))
    cyl.set("h", str(height))
    cyl.setIndex("pos", "0", 0)
    cyl.setIndex("pos", "0", 1)
    cyl.setIndex("pos", str(z_pos), 2)


def _add_difference(
    geom_java: Any,
    tag: str,
    label: str,
    *,
    inputs: list[str],
    subtract: list[str],
) -> None:
    diff = geom_java.feature().create(tag, "Difference")
    # `input` / `input2` on Boolean features are Selections — set via the
    # selection accessor rather than `set(...)`.
    diff.selection("input").set(inputs)
    diff.selection("input2").set(subtract)
    diff.label(label)


def _add_ball_selection(
    geom_java: Any,
    tag: str,
    label: str,
    *,
    x: float,
    y: float,
    z: float,
    radius: float = 0.5,
) -> None:
    """Add a BallSelection picking the single domain that contains (x, y, z).

    A small ball at a known interior point intersects exactly one domain.
    The selection's tag at the model level becomes ``<geom>_<tag>`` (e.g.
    ``geom1_sel_axon``) — materials & boundary conditions reference it by that.
    """

    sel = geom_java.feature().create(tag, "BallSelection")
    sel.set("entitydim", "3")
    sel.set("posx", str(x))
    sel.set("posy", str(y))
    sel.set("posz", str(z))
    sel.set("r", str(radius))
    sel.set("condition", "intersects")
    sel.label(label)


def selection_tag(label: str, geom_tag: str = "geom1") -> str:
    """Model-level tag of the geometry BallSelection for the given domain label."""

    return f"{geom_tag}_sel_{label}"


def build_geometry(model: Any, params: GeometryParams) -> int:
    """Build the five-domain nerve geometry. Returns the resulting domain count."""

    comp = _ensure_component(model)
    geom = _ensure_geometry(comp)
    _clear_geometry_features(geom)

    L_inter = params.internode_length_um
    L_node = params.node_length_um
    L_total = total_length_um(params)
    r_axon = params.axon_radius_um
    r_myelin = params.myelin_radius_um
    r_ext = params.external_radius_um

    # ---- Inner solid cylinder (axon, full length, single domain) ----
    _add_cyl(geom, "cyl_axon", radius=r_axon, height=L_total, z_pos=0.0)
    geom.feature("cyl_axon").label("Axon")

    # ---- Outer + inner cylinders for the four annular regions ----
    # Proximal sheath segment
    _add_cyl(geom, "cyl_myp_out", radius=r_myelin, height=L_inter, z_pos=0.0)
    _add_cyl(geom, "cyl_myp_in", radius=r_axon, height=L_inter, z_pos=0.0)
    _add_difference(geom, "myp", "Myelin proximal",
                    inputs=["cyl_myp_out"], subtract=["cyl_myp_in"])

    # Node segment
    _add_cyl(geom, "cyl_nd_out", radius=r_myelin, height=L_node, z_pos=L_inter)
    _add_cyl(geom, "cyl_nd_in", radius=r_axon, height=L_node, z_pos=L_inter)
    _add_difference(geom, "nd", "Node",
                    inputs=["cyl_nd_out"], subtract=["cyl_nd_in"])

    # Distal sheath segment
    z_dist = L_inter + L_node
    _add_cyl(geom, "cyl_myd_out", radius=r_myelin, height=L_inter, z_pos=z_dist)
    _add_cyl(geom, "cyl_myd_in", radius=r_axon, height=L_inter, z_pos=z_dist)
    _add_difference(geom, "myd", "Myelin distal",
                    inputs=["cyl_myd_out"], subtract=["cyl_myd_in"])

    # External annulus (full length, r=myelin_r .. ext_r)
    _add_cyl(geom, "cyl_ext_out", radius=r_ext, height=L_total, z_pos=0.0)
    _add_cyl(geom, "cyl_ext_in", radius=r_myelin, height=L_total, z_pos=0.0)
    _add_difference(geom, "ext", "External medium",
                    inputs=["cyl_ext_out"], subtract=["cyl_ext_in"])

    # Named ball selections — one per labeled domain — for downstream
    # material and boundary-condition assignment.
    centers = domain_centers(params)
    for label, (x, y, z) in centers.items():
        _add_ball_selection(
            geom, f"sel_{label}", f"Selection: {label}",
            x=x, y=y, z=z, radius=0.5,
        )

    geom.run()
    return int(geom.getNDomains())
