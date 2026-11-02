"""THz dielectric models (exp(+jwt) convention: eps = eps' - j eps'').

All models are plain functions of frequency in Hz returning complex relative
permittivity. Parameters and sources are kept next to each model; the
``scale_*`` helpers implement the +/-20 % robustness checks.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

EPS0 = 8.8541878128e-12
MU0 = 1.25663706212e-6
C0 = 299792458.0
ETA0 = MU0 * C0


@dataclass(frozen=True)
class DoubleDebye:
    """eps = eps_inf + (eps_s - eps_1)/(1 + j w tau1) + (eps_1 - eps_inf)/(1 + j w tau2) + sigma/(j w eps0)."""

    eps_s: float
    eps_1: float
    eps_inf: float
    tau1_s: float
    tau2_s: float
    sigma_dc: float = 0.0
    source: str = ""

    def __call__(self, f_hz):
        w = 2 * np.pi * np.asarray(f_hz, dtype=float)
        eps = (
            self.eps_inf
            + (self.eps_s - self.eps_1) / (1 + 1j * w * self.tau1_s)
            + (self.eps_1 - self.eps_inf) / (1 + 1j * w * self.tau2_s)
        )
        if self.sigma_dc:
            eps = eps + self.sigma_dc / (1j * w * EPS0)
        return eps


# Kindt & Schmuttenmaer, J. Phys. Chem. 100, 10373 (1996), water at 25 C, THz-TDS 0.1-2 THz.
WATER_KS1996 = DoubleDebye(78.36, 5.16, 3.49, 8.24e-12, 0.18e-12, source="Kindt & Schmuttenmaer 1996")
# Legacy model used by Sims 1-26 (eps_s = 78.4, eps_1 = 4.90): kept as an alternative.
WATER_LEGACY = DoubleDebye(78.40, 4.90, 3.17, 9.36e-12, 0.30e-12, source="legacy (Sims 1-26)")


@dataclass(frozen=True)
class ConstantEps:
    eps: complex
    source: str = ""

    def __call__(self, f_hz):
        return np.full(np.shape(f_hz), self.eps, dtype=complex) if np.ndim(f_hz) else complex(self.eps)


def with_sigma(model, sigma_s_per_m: float):
    """Add an extra (frequency-independent) conductivity to any model."""

    def eps(f_hz):
        w = 2 * np.pi * np.asarray(f_hz, dtype=float)
        return model(f_hz) + sigma_s_per_m / (1j * w * EPS0)

    return eps


def scaled(model, re_scale: float = 1.0, im_scale: float = 1.0):
    """Scale eps' and eps'' independently (robustness checks)."""

    def eps(f_hz):
        e = model(f_hz)
        return re_scale * e.real + 1j * im_scale * e.imag

    return eps


def refractive_index(eps):
    """n = n' - j kappa with n' > 0 (principal sqrt in exp(+jwt) convention)."""
    n = np.sqrt(np.asarray(eps, dtype=complex))
    return np.where(n.real < 0, -n, n)


def alpha_per_cm(eps, f_hz):
    """Power absorption coefficient 2 k0 kappa in 1/cm."""
    n = refractive_index(eps)
    return 2 * (2 * np.pi * np.asarray(f_hz) / C0) * (-n.imag) / 100.0


def sigma_eff(eps, f_hz):
    """Effective THz conductivity w eps0 eps'' (S/m)."""
    return 2 * np.pi * np.asarray(f_hz) * EPS0 * (-np.imag(eps))
