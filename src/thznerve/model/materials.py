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
#
# IMPORTANT: COMSOL EWFD silently drops the imaginary part of analytic
# expressions like ``73.5/(1 + i*2*pi*freq*9.36e-12)`` when used for
# `relpermittivity`. A probe verified that *literal* complex values
# ("3.99 - 2.94*i") produce different fields than purely-real expressions,
# but the analytic Debye form is treated as if its Im part were zero.
#
# So: for any sim that depends on Im(ε), pass ``freq_hz=...`` to
# ``apply_materials`` and the water εr is baked as a literal complex value
# computed at that frequency. Same trick as the σ-at-node fix.

# Analytic water εr — kept for completeness, but DON'T use this directly
# for εr assignment in EWFD. Use `water_eps_expr(freq_hz)` instead.
_WATER_EPS_EXPR_ANALYTIC: str = (
    "3.17"
    " + 73.5/(1 + i*2*pi*freq*9.36e-12)"
    " + 1.73/(1 + i*2*pi*freq*0.30e-12)"
)

# Backwards-compat alias — points at the analytic form for callers that don't
# yet know about the literal-at-frequency requirement. New code should use
# `water_eps_expr(freq_hz)`.
WATER_EPS_EXPR: str = _WATER_EPS_EXPR_ANALYTIC

# Myelin — constant complex. This one is already literal so works as-is.
MYELIN_EPS_EXPR: str = "4.5 - 0.5*i"


def water_eps_expr(freq_hz: float | None = None) -> str:
    """COMSOL water εr expression.

    If ``freq_hz`` is given, returns a literal complex value computed from
    the double-Debye formula at that frequency — this is the only encoding
    that actually delivers a non-zero Im(ε) into EWFD's solver.

    Without ``freq_hz``, returns the analytic expression. EWFD will treat
    its Im part as zero (lossless water) — use only when you know loss
    won't matter (e.g. checking geometry).
    """

    if freq_hz is None:
        return _WATER_EPS_EXPR_ANALYTIC
    eps = complex(debye_water_epsilon(float(freq_hz)))
    # COMSOL accepts "3.99 + (-2.94)*i" style; the sign goes on the literal.
    return f"({eps.real:.8g}) + ({eps.imag:.8g})*i"


def node_eps_expr(sigma_S_per_m: float, freq_hz: float | None = None) -> str:
    """COMSOL node εr = water + i·σ/(ω·ε₀), baked as a literal at ``freq_hz``.

    Same encoding requirement as ``water_eps_expr``: COMSOL only honours Im(ε)
    when the value is a literal complex number, so we precompute the whole
    thing numerically.
    """

    if freq_hz is None:
        # Analytic fallback (Im will be silently dropped). Useful only for
        # σ = 0 cases or geometry-only checks.
        return _WATER_EPS_EXPR_ANALYTIC
    eps_water = complex(debye_water_epsilon(float(freq_hz)))
    omega = 2 * np.pi * float(freq_hz)
    # In e^{-iωt} convention, loss has negative Im(ε); add the σ contribution.
    eps_node = eps_water + (-sigma_S_per_m / (omega * EPS_0)) * 1j
    return f"({eps_node.real:.8g}) + ({eps_node.imag:.8g})*i"


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


def apply_materials(
    model: Any,
    params: MaterialParams,
    *,
    freq_hz: float | None = None,
) -> dict[str, str]:
    """Apply water / myelin / node materials to the geometry's labeled selections.

    If ``freq_hz`` is supplied (and σ != 0), the σ contribution is baked into
    the node εr as a literal complex term at that frequency — the only
    encoding that actually changes the EWFD field (see materials.py docstring
    and the σ-encoding section of SESSION.md).
    """

    from thznerve.model.geometry import selection_tag

    java = model.java
    comp_tag = str(java.component().tags()[0])
    comp = java.component(comp_tag)

    # Water — covers axon and external. Bake the literal complex εr at the
    # simulation frequency so EWFD sees the imaginary (loss) part.
    water_expr = water_eps_expr(freq_hz=freq_hz)
    mat_water = _create_material(comp, "mat_water", "Water (double Debye)", water_expr)
    mat_water.selection().named(selection_tag("axon"))
    mat_water_ext = _create_material(comp, "mat_water_ext", "Water (external)", water_expr)
    mat_water_ext.selection().named(selection_tag("external"))

    # Myelin — covers proximal & distal sheath segments.
    mat_myp = _create_material(comp, "mat_myelin_proximal", "Myelin (proximal)", MYELIN_EPS_EXPR)
    mat_myp.selection().named(selection_tag("myelin_proximal"))
    mat_myd = _create_material(comp, "mat_myelin_distal", "Myelin (distal)", MYELIN_EPS_EXPR)
    mat_myd.selection().named(selection_tag("myelin_distal"))

    # Node — σ baked into εr as a literal complex term at `freq_hz`.
    mat_node = _create_material(
        comp, "mat_node", "Node of Ranvier",
        node_eps_expr(params.node_sigma_S_per_m, freq_hz=freq_hz),
    )
    mat_node.selection().named(selection_tag("node"))

    return {
        "mat_water": "axon",
        "mat_water_ext": "external",
        "mat_myelin_proximal": "myelin_proximal",
        "mat_myelin_distal": "myelin_distal",
        "mat_node": "node",
    }
