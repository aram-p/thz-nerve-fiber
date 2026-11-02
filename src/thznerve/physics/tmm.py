"""Normal-incidence transfer-matrix method for planar multilayers (exp(+jwt))."""

from __future__ import annotations

import numpy as np

from thznerve.physics.materials import C0, refractive_index


def slab_rta(f_hz: float, n_in: complex, layers: list[tuple[complex, float]], n_out: complex):
    """Reflectance, transmittance, absorptance of a stack at normal incidence.

    layers: list of (relative permittivity, thickness in m). n_in must be real
    (lossless incidence medium) for R+T+A=1 to be meaningful.
    Returns (R, T, A, r, t) with t the field transmission coefficient.
    """
    k0 = 2 * np.pi * f_hz / C0
    ns = [complex(n_in)] + [complex(refractive_index(e)) for e, _ in layers] + [complex(n_out)]
    M = np.eye(2, dtype=complex)
    for i, (_, d) in enumerate(layers, start=1):
        n = ns[i]
        delta = k0 * n * d
        # characteristic matrix of a layer (Born & Wolf), exp(+jwt): use cos, j sin
        Mi = np.array([[np.cos(delta), 1j * np.sin(delta) / n],
                       [1j * n * np.sin(delta), np.cos(delta)]], dtype=complex)
        M = M @ Mi
    n1, n2 = ns[0], ns[-1]
    B, C = M @ np.array([1.0, n2])
    r = (n1 * B - C) / (n1 * B + C)
    t = 2 * n1 / (n1 * B + C)
    R = abs(r) ** 2
    T = (n2.real / n1.real) * abs(t) ** 2
    return R, T, 1 - R - T, r, t
