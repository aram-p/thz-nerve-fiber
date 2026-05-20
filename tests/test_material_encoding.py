"""Regression tests for the COMSOL material-expression encoders.

These guard against the class of bug where an εr expression evaluates
to a different value (or has its imaginary part silently dropped) by
COMSOL than the Python implementation computes. They don't run COMSOL —
they verify the *encoding* the COMSOL side will receive.

The σ encoding bug (analytic σ/(ω·ε₀) expression silently evaluating
to zero in EWFD) wasted hours of overnight session compute. A simple
unit test on the encoder would have caught it instantly. Specifically:

* ``water_eps_expr(freq_hz)`` must return a literal complex value at
  that frequency, not an analytic expression that COMSOL may strip.
* ``node_eps_expr(σ, freq_hz)`` must produce a different string for
  different σ values, with the imaginary part baked numerically.
* The literal-complex format must use a syntactic form that COMSOL
  parses (the legacy ``WATER_EPS_EXPR`` and analytic σ expressions
  did not).
"""

from __future__ import annotations

import re

import numpy as np
import pytest

from thznerve.model.materials import (
    EPS_0,
    _WATER_EPS_EXPR_ANALYTIC,
    debye_water_epsilon,
    node_eps_expr,
    water_eps_expr,
)


_NUM_PATTERN = r"[+-]?\d+\.?\d*(?:[eE][+-]?\d+)?"


def _parse_literal_complex(expr: str) -> complex:
    """Parse expressions like ``(3.94) + (-2.83)*i`` into a Python complex.

    Tolerates the few literal-complex syntaxes the encoders produce:
    ``A + B*i``, ``A - B*i``, ``A + (-B)*i``, ``(A) + (B)*i``, and
    plain real ``A``.
    """

    s = expr.replace(" ", "").replace("(", "").replace(")", "")
    # Normalise ``A+-B*i`` → ``A-B*i`` and ``A-+B*i`` → ``A-B*i`` (rare but
    # produced by the format string ``"{a:.8g} + {b:.8g}*i"`` when b is -ve).
    s = s.replace("+-", "-").replace("-+", "-")

    # A ± B*i (with optional *)
    m = re.fullmatch(rf"({_NUM_PATTERN})([+-]{_NUM_PATTERN})\*?i", s)
    if m:
        return complex(float(m.group(1)), float(m.group(2)))
    # Just A (real)
    m = re.fullmatch(rf"({_NUM_PATTERN})", s)
    if m:
        return complex(float(m.group(1)), 0.0)
    raise ValueError(f"could not parse {expr!r}")


@pytest.mark.parametrize("f_hz", [0.1e12, 0.6e12, 1.0e12, 2.0e12])
def test_water_eps_expr_matches_python_implementation(f_hz: float):
    """Literal-complex water εr matches debye_water_epsilon to 6 digits."""

    expr = water_eps_expr(freq_hz=f_hz)
    parsed = _parse_literal_complex(expr)
    expected = complex(debye_water_epsilon(f_hz))
    assert parsed.real == pytest.approx(expected.real, rel=1e-6)
    assert parsed.imag == pytest.approx(expected.imag, rel=1e-6)


def test_water_eps_expr_is_literal_at_frequency():
    """No ``freq`` variable reference — must be a literal at the given f."""

    expr = water_eps_expr(freq_hz=0.6e12)
    assert "freq" not in expr, (
        f"water_eps_expr({0.6e12:g}) returned analytic form referencing `freq`: {expr!r}"
    )


def test_water_eps_expr_no_freq_returns_analytic():
    """Backwards-compat: omitting freq_hz returns the analytic form."""

    expr = water_eps_expr(freq_hz=None)
    assert expr == _WATER_EPS_EXPR_ANALYTIC


@pytest.mark.parametrize("sigma", [0.0, 0.1, 1.0, 10.0, 100.0, 1000.0])
def test_node_eps_expr_changes_with_sigma(sigma: float):
    """Different σ values must produce *different* expression strings."""

    f_hz = 0.6e12
    expr_0 = node_eps_expr(0.0, freq_hz=f_hz)
    expr_sigma = node_eps_expr(sigma, freq_hz=f_hz)
    if sigma == 0:
        assert expr_sigma == expr_0
    else:
        assert expr_sigma != expr_0, (
            f"node_eps_expr at σ=0 vs σ={sigma} returned the same string: {expr_sigma!r}"
        )


def test_node_eps_expr_im_part_scales_linearly_with_sigma():
    """At a fixed frequency, the Im part of node εr is linear in σ
    (the only frequency-dependent thing is water's natural loss)."""

    f_hz = 0.6e12
    sigmas = np.array([0.0, 1.0, 10.0, 100.0])
    ims = []
    for s in sigmas:
        eps = _parse_literal_complex(node_eps_expr(float(s), freq_hz=f_hz))
        ims.append(eps.imag)
    ims = np.array(ims)
    expected_slope = -1.0 / (2 * np.pi * f_hz * EPS_0)
    # ε_node(σ) = ε_water + σ · (slope)·i  → ε_node.imag = ε_water.imag + σ·slope
    derived = (ims[1:] - ims[0]) / sigmas[1:]
    for d in derived:
        assert d == pytest.approx(expected_slope, rel=1e-3), (
            f"node Im scaling: got {d:.6g}, expected {expected_slope:.6g}"
        )


def test_water_and_node_match_at_sigma_zero():
    """At σ=0, node εr == water εr (they share the same Debye expression)."""

    f_hz = 0.6e12
    w = _parse_literal_complex(water_eps_expr(freq_hz=f_hz))
    n = _parse_literal_complex(node_eps_expr(0.0, freq_hz=f_hz))
    assert n.real == pytest.approx(w.real)
    assert n.imag == pytest.approx(w.imag)


def test_water_im_part_is_negative_throughout_thz():
    """Water at THz is lossy (Im(ε) < 0 in -iωt convention)."""

    for f_hz in [0.1e12, 0.6e12, 1.0e12, 2.0e12, 3.0e12]:
        eps = _parse_literal_complex(water_eps_expr(freq_hz=f_hz))
        assert eps.imag < 0, f"Im(ε_water) at {f_hz/1e12:g} THz = {eps.imag:g} (should be < 0)"


def test_node_loss_exceeds_water_at_finite_sigma():
    """At σ > 0, the node should be MORE lossy than water (more negative Im)."""

    f_hz = 0.6e12
    w = _parse_literal_complex(water_eps_expr(freq_hz=f_hz))
    for sigma in [0.1, 1.0, 10.0]:
        n = _parse_literal_complex(node_eps_expr(sigma, freq_hz=f_hz))
        assert n.imag < w.imag, (
            f"node Im at σ={sigma} = {n.imag:g}, expected < water Im = {w.imag:g}"
        )
