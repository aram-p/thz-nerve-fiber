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


    def comsol(self) -> str:
        """The same model as a COMSOL expression in `freq` (exp(+jwt))."""
        s = (f"({self.eps_inf})+({self.eps_s - self.eps_1})/(1+i*2*pi*freq*{self.tau1_s}[s])"
             f"+({self.eps_1 - self.eps_inf})/(1+i*2*pi*freq*{self.tau2_s}[s])")
        if self.sigma_dc:
            s += f"-i*({self.sigma_dc}[S/m])/(2*pi*freq*epsilon0_const)"
        return f"({s})"


# Water. Values from docs/literature.md section 2a.
# Kindt & Schmuttenmaer, J. Phys. Chem. 100, 10373 (1996), 25 C, THz-TDS ~0.1-2 THz.
WATER_KS1996 = DoubleDebye(78.36, 4.93, 3.48, 8.24e-12, 0.180e-12, source="Kindt & Schmuttenmaer 1996, 25 C")
# Jepsen, Moller, Merbold, Opt. Express 15, 14717 (2007), room temperature.
WATER_JEPSEN2007 = DoubleDebye(78.36, 5.16, 3.49, 7.89e-12, 0.181e-12, source="Jepsen et al. 2007")
# Liebe, Hufford, Manabe, Int. J. IR MM Waves 12, 659 (1991) evaluated at 37 C (valid < 1 THz).
WATER_LIEBE37 = DoubleDebye(74.28, 4.98, 3.77, 6.28e-12, 0.158e-12, source="Liebe et al. 1991 at 37 C")
# Legacy model used by Sims 1-26: kept only to reproduce them.
WATER_LEGACY = DoubleDebye(78.40, 4.90, 3.17, 9.36e-12, 0.30e-12, source="legacy (Sims 1-26)")


@dataclass(frozen=True)
class ConstantEps:
    eps: complex
    source: str = ""

    def __call__(self, f_hz):
        return np.full(np.shape(f_hz), self.eps, dtype=complex) if np.ndim(f_hz) else complex(self.eps)

    def comsol(self) -> str:
        return f"({self.eps.real:.10g}+({self.eps.imag:.10g})*i)"


@dataclass(frozen=True)
class Looyenga:
    """Looyenga mixture: eps^(1/3) = phi eps_a^(1/3) + (1 - phi) eps_b^(1/3).

    Used for myelin (water fraction phi ~0.4 in lipid/protein) and axoplasm
    (phi ~0.85 water in protein), docs/literature.md sections 2b, 2d.
    """

    a: object
    b: object
    phi: float
    source: str = ""

    def __call__(self, f_hz):
        ea = np.asarray(self.a(f_hz), dtype=complex)
        eb = np.asarray(self.b(f_hz), dtype=complex)
        return (self.phi * ea ** (1 / 3) + (1 - self.phi) * eb ** (1 / 3)) ** 3

    def comsol(self) -> str:
        return f"(({self.phi})*({self.a.comsol()})^(1/3)+({1 - self.phi})*({self.b.comsol()})^(1/3))^3"


@dataclass(frozen=True)
class Scaled:
    """eps' and eps'' scaled independently (+/-20 % robustness checks)."""

    m: object
    re: float = 1.0
    im: float = 1.0

    def __call__(self, f_hz):
        e = np.asarray(self.m(f_hz), dtype=complex)
        return self.re * e.real + 1j * self.im * e.imag

    def comsol(self) -> str:
        e = self.m.comsol()
        return f"({self.re}*real({e})+i*{self.im}*imag({e}))"


@dataclass(frozen=True)
class PlusSigma:
    """Extra frequency-independent conductivity sigma (S/m): eps - j sigma/(w eps0)."""

    m: object
    sigma: float

    def __call__(self, f_hz):
        w = 2 * np.pi * np.asarray(f_hz, dtype=float)
        return self.m(f_hz) + self.sigma / (1j * w * EPS0)

    def comsol(self) -> str:
        return f"({self.m.comsol()}-i*({self.sigma}[S/m])/(2*pi*freq*epsilon0_const))"


# Dry lipid/protein of myelin and protein of axoplasm (docs/literature.md 2b, 2d).
DRY_LIPID = ConstantEps(2.3 - 0.05j, source="dry lipid/protein, n ~ 1.52 (Berry 2003 adipose, lipid TDS)")
PROTEIN = ConstantEps(2.5 - 0.1j, source="dry protein, eps ~ 2.5")


def myelin(water=WATER_KS1996, phi: float = 0.4):
    """Myelin: ~40 % water in lipid/protein (Basic Neurochemistry), Looyenga mix."""
    return Looyenga(water, DRY_LIPID, phi, source=f"Looyenga water {phi:.2f} + lipid")


def axoplasm(water=WATER_KS1996, phi: float = 0.85):
    """Axoplasm: ~80-85 % water plus protein, Looyenga mix."""
    return Looyenga(water, PROTEIN, phi, source=f"Looyenga water {phi:.2f} + protein")


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
