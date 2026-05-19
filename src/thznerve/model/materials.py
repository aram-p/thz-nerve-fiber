"""Frequency-dependent materials for the THz nerve-fiber model.

Three materials per the thesis (issue #9):

* **Double-Debye water** — applied to *axon* and *external* domains.
  ε(f) = 3.17 + 73.5 / (1 + i·2π·f·9.36 ps) + 1.73 / (1 + i·2π·f·0.30 ps)

* **Myelin** — constant at THz, ε = 4.5 − 0.5 i. Applied to both *myelin_proximal*
  and *myelin_distal*.

* **Node of Ranvier** — water plus a node-conductivity term.
  ε_node(f, σ) = ε_water(f) + i · σ / (2π · f · ε₀).
  Applied to *node*.

Each material is exposed two ways:

* As a Python function for unit-testing and figure scripts.
* As a COMSOL-side analytic expression string for `apply_materials()` to bind
  into the material node's relative-permittivity property. COMSOL substitutes
  its built-in `freq` (set by the EWFD physics) and `epsilon0_const`.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from pydantic import BaseModel, ConfigDict, Field

# --- Physical constants -----------------------------------------------------

EPS_0: float = 8.8541878128e-12  # vacuum permittivity, F/m

# Double-Debye water parameters (Ellison 2007 / Jepsen et al. fits at THz).
DEBYE_EPS_INF: float = 3.17
DEBYE_DELTA_EPS_1: float = 73.5
DEBYE_DELTA_EPS_2: float = 1.73
DEBYE_TAU_1_S: float = 9.36e-12  # 9.36 ps
DEBYE_TAU_2_S: float = 0.30e-12  # 0.30 ps

# Myelin (constant at THz).
MYELIN_EPS_REAL: float = 4.5
MYELIN_EPS_IMAG: float = -0.5

# --- Material parameters ----------------------------------------------------


class MaterialParams(BaseModel):
    """Material parameters that vary per simulation (only σ at the node, currently)."""

    model_config = ConfigDict(extra="forbid", protected_namespaces=())

    node_sigma_S_per_m: float = Field(ge=0)


# --- Python implementations (testable) --------------------------------------


def debye_water_epsilon(freq_hz: float | np.ndarray) -> complex | np.ndarray:
    """Complex relative permittivity of water (double-Debye fit)."""

    omega = 2 * np.pi * np.asarray(freq_hz, dtype=float)
    eps = (
        DEBYE_EPS_INF
        + DEBYE_DELTA_EPS_1 / (1 + 1j * omega * DEBYE_TAU_1_S)
        + DEBYE_DELTA_EPS_2 / (1 + 1j * omega * DEBYE_TAU_2_S)
    )
    return eps


def myelin_epsilon(freq_hz: float | np.ndarray) -> complex:
    """Complex relative permittivity of myelin (constant at THz)."""

    return MYELIN_EPS_REAL + 1j * MYELIN_EPS_IMAG


def node_epsilon(freq_hz: float | np.ndarray, sigma_S_per_m: float) -> complex | np.ndarray:
    """Complex relative permittivity at the node — water plus conductivity term."""

    omega = 2 * np.pi * np.asarray(freq_hz, dtype=float)
    return debye_water_epsilon(freq_hz) + 1j * sigma_S_per_m / (omega * EPS_0)


# --- COMSOL analytic expressions --------------------------------------------

# Water — references COMSOL's built-in `freq` (Hz).
WATER_EPS_EXPR: str = (
    "3.17"
    " + 73.5/(1 + i*2*pi*freq*9.36e-12)"
    " + 1.73/(1 + i*2*pi*freq*0.30e-12)"
)

# Myelin — constant complex.
MYELIN_EPS_EXPR: str = "4.5 - 0.5*i"

# Node — water plus σ/(ωε₀) term. ``epsilon0_const`` is a COMSOL built-in.
# ``node_sigma_S_per_m`` was originally a global parameter; we now inject the
# numeric value directly via ``node_eps_expr(sigma)`` so the material
# expression text changes between σ iterations — defeating COMSOL's
# stationary-solver caching that previously returned identical solutions for
# every σ.
NODE_EPS_EXPR_TEMPLATE: str = (
    f"({WATER_EPS_EXPR})"
    " + i*{sigma}/(2*pi*freq*epsilon0_const)"
)
# Legacy symbol kept for any external references (uses 0 for σ).
NODE_EPS_EXPR: str = NODE_EPS_EXPR_TEMPLATE.format(sigma="0")


def node_eps_expr(sigma_S_per_m: float) -> str:
    """COMSOL analytic expression for the node ε with σ baked in."""

    return NODE_EPS_EXPR_TEMPLATE.format(sigma=f"{sigma_S_per_m:g}")


# --- COMSOL binding ---------------------------------------------------------


def _create_material(
    comp_java: Any,
    tag: str,
    label: str,
    eps_expr: str,
    *,
    sigma_expr: str = "0",
) -> Any:
    """Create a Common material with isotropic ``relpermittivity = eps_expr`` and
    optional electrical conductivity ``sigma_expr`` (S/m).

    EWFD treats ``electricconductivity`` and ``relpermittivity`` as
    independent constitutive parameters — encoding σ as the imaginary part
    of ε is ignored unless ``electricconductivity`` is also set. So we keep
    them separate.
    """

    materials = comp_java.material()
    existing = [str(t) for t in materials.tags()]
    if tag in existing:
        materials.remove(tag)
    mat = materials.create(tag, "Common")
    mat.label(label)

    basic = mat.propertyGroup("def")
    basic.set(
        "relpermittivity",
        [
            eps_expr, "0", "0",
            "0", eps_expr, "0",
            "0", "0", eps_expr,
        ],
    )
    basic.set(
        "electricconductivity",
        [
            sigma_expr, "0", "0",
            "0", sigma_expr, "0",
            "0", "0", sigma_expr,
        ],
    )
    basic.set("relpermeability", ["1", "0", "0", "0", "1", "0", "0", "0", "1"])
    return mat


def apply_materials(model: Any, params: MaterialParams) -> dict[str, str]:
    """Apply water / myelin / node materials to the geometry's labeled selections.

    Returns a dict mapping material tag → list of domain-label assignments
    (for logging / debugging by the caller).
    """

    from thznerve.model.geometry import selection_tag

    java = model.java
    comp_tag = str(java.component().tags()[0])
    comp = java.component(comp_tag)

    # Water — covers axon and external.
    mat_water = _create_material(comp, "mat_water", "Water (double Debye)", WATER_EPS_EXPR)
    mat_water.selection().named(selection_tag("axon"))
    mat_water_ext = _create_material(comp, "mat_water_ext", "Water (external)", WATER_EPS_EXPR)
    mat_water_ext.selection().named(selection_tag("external"))

    # Myelin — covers proximal & distal sheath segments.
    mat_myp = _create_material(comp, "mat_myelin_proximal", "Myelin (proximal)", MYELIN_EPS_EXPR)
    mat_myp.selection().named(selection_tag("myelin_proximal"))
    mat_myd = _create_material(comp, "mat_myelin_distal", "Myelin (distal)", MYELIN_EPS_EXPR)
    mat_myd.selection().named(selection_tag("myelin_distal"))

    # Node — water permittivity, with σ on the electricconductivity property
    # (NOT folded into Im(ε), because EWFD reads the two properties separately).
    mat_node = _create_material(
        comp, "mat_node", "Node of Ranvier",
        WATER_EPS_EXPR,
        sigma_expr=f"{params.node_sigma_S_per_m:g}",
    )
    mat_node.selection().named(selection_tag("node"))

    return {
        "mat_water": "axon",
        "mat_water_ext": "external",
        "mat_myelin_proximal": "myelin_proximal",
        "mat_myelin_distal": "myelin_distal",
        "mat_node": "node",
    }
