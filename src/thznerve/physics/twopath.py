"""Two-path (partial air path) model of a THz-TDS transmission measurement.

Hypothesis for the reported voltage-induced "resonances": a fraction of the beam
crosses the cuvette through air (gas bubble, crack, dried/shrunk region) instead
of tissue. The detected field is
    T_on(f) = (1 - p) T_tissue(f) + p T_air(f),   T_air = 1 (same as empty cuvette)
Since the tissue path is delayed by dt = (n - 1) d / c relative to air and is
strongly attenuated, the two paths cancel at f_k = (2k - 1) / (2 dt) when
p ~ |T_tissue(f_k)|. The standard retrieval n = 1 + c phi / (w d) then shows a
2 pi phase slip at the notch: n dips below 1 on one side and spikes on the other.
No fibre physics enters.
"""

from __future__ import annotations

import numpy as np

from thznerve.physics.materials import C0


def tissue_transmission(f_hz, n_real, alpha_np_per_m, d_m):
    """Bulk transmission of a slab relative to an equal air gap (no reflections)."""
    k0 = 2 * np.pi * np.asarray(f_hz) / C0
    n = np.asarray(n_real)
    t_fresnel = 4 * n / (1 + n) ** 2
    return t_fresnel * np.exp(-1j * (n - 1) * k0 * d_m) * np.exp(-0.5 * alpha_np_per_m * d_m)


def two_path(T_tissue, p):
    return (1 - p) * T_tissue + p


def retrieve_n(f_hz, T, d_m):
    """Textbook retrieval n = 1 + c * unwrap(-arg T) / (w d), unwrapped from low f."""
    phi = np.unwrap(-np.angle(T))
    w = 2 * np.pi * np.asarray(f_hz)
    return 1 + C0 * phi / (w * d_m)


def notch_frequencies(n_real, d_m, kmax=3):
    dt = (n_real - 1) * d_m / C0
    return np.array([(2 * k - 1) / (2 * dt) for k in range(1, kmax + 1)])
