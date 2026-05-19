"""Unit tests for the Python ε(f) implementations.

The COMSOL-side analytic expressions in `materials.py` use the same
constants and same algebraic form as these functions, so matching here
is sufficient validation of the formulas themselves.
"""

from __future__ import annotations

import numpy as np
import pytest

from thznerve.model.materials import (
    EPS_0,
    MYELIN_EPS_IMAG,
    MYELIN_EPS_REAL,
    debye_water_epsilon,
    myelin_epsilon,
    node_epsilon,
)

# Hand-computed reference values at f = 0.6 THz.
# Derived from:
#   ω·τ1 = 2π·0.6e12·9.36e-12 = 35.286
#   term1 = 73.5/(1 + i·35.286) = 0.05897 - 2.0815i
#   ω·τ2 = 2π·0.6e12·0.30e-12 = 1.1310
#   term2 = 1.73/(1 + i·1.131)  = 0.7589  - 0.8584i
#   total = 3.17 + term1 + term2 ≈ 3.988 - 2.940i
_F_TEST_HZ = 0.6e12
_EPS_WATER_EXPECTED = 3.9882 - 2.9399j


def test_debye_water_at_0p6_thz_within_1pct():
    eps = complex(debye_water_epsilon(_F_TEST_HZ))
    # 1 % tolerance on |ε|
    rel = abs(eps - _EPS_WATER_EXPECTED) / abs(_EPS_WATER_EXPECTED)
    assert rel < 0.01, f"|Δε|/|ε| = {rel:.4f} (got {eps}, expected {_EPS_WATER_EXPECTED})"


def test_debye_water_vectorised():
    freqs = np.array([0.1e12, 0.6e12, 2.0e12])
    eps = debye_water_epsilon(freqs)
    assert eps.shape == freqs.shape
    # Real part of water ε(f) decreases monotonically through THz.
    assert eps.real[0] > eps.real[1] > eps.real[2]
    # Imaginary part is negative (lossy).
    assert (eps.imag < 0).all()


def test_myelin_is_constant():
    assert myelin_epsilon(0.6e12) == MYELIN_EPS_REAL + 1j * MYELIN_EPS_IMAG
    assert myelin_epsilon(2.0e12) == MYELIN_EPS_REAL + 1j * MYELIN_EPS_IMAG


def test_node_reduces_to_water_when_sigma_zero():
    eps_w = debye_water_epsilon(_F_TEST_HZ)
    eps_n = node_epsilon(_F_TEST_HZ, sigma_S_per_m=0.0)
    assert eps_n == pytest.approx(eps_w)


def test_node_sigma_imaginary_offset():
    sigma = 1.0
    eps_w = complex(debye_water_epsilon(_F_TEST_HZ))
    eps_n = complex(node_epsilon(_F_TEST_HZ, sigma_S_per_m=sigma))
    expected_offset = 1j * sigma / (2 * np.pi * _F_TEST_HZ * EPS_0)
    assert (eps_n - eps_w) == pytest.approx(expected_offset)
