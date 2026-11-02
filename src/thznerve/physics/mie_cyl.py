"""Multilayer infinite circular cylinder, normal incidence, exact series (exp(+jwt)).

Layers: core (index 0) ... host (index N). ``radii`` are the N interface radii
(ascending, metres). The host may be lossy: the incident plane wave has
amplitude E0 = 1 V/m *at the cylinder axis*, and cross-sections are normalised
by the local plane-wave intensity I0 = |E0|^2 Re(n_host)/(2 eta0).

pol = "TM": E parallel to the axis (Ez).   pol = "TE": E perpendicular (Hz).
In layer l the axial field is A_l J_m(k_l r) + B_l H2_m(k_l r).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import h2vp, hankel2, jv, jvp

from thznerve.physics.materials import C0, EPS0, ETA0, MU0, refractive_index


@dataclass
class CylSolution:
    f_hz: float
    pol: str
    radii: np.ndarray
    eps: np.ndarray
    ms: np.ndarray
    coef: np.ndarray  # (n_m, n_layers, 2): A_l, B_l per order

    @property
    def omega(self):
        return 2 * np.pi * self.f_hz

    @property
    def k(self):
        return self.omega / C0 * refractive_index(self.eps)

    def layer_of(self, r):
        return int(np.searchsorted(self.radii, r, side="right"))

    def axial_and_deriv(self, r, l=None):
        """Per-order axial field F_m(r) and dF_m/dr in layer l (default: layer containing r)."""
        if l is None:
            l = self.layer_of(r)
        kl = self.k[l]
        A, B = self.coef[:, l, 0], self.coef[:, l, 1]
        x = kl * r
        F = A * jv(self.ms, x) + (B * hankel2(self.ms, x) if np.any(B) else 0)
        dF = kl * (A * jvp(self.ms, x) + (B * h2vp(self.ms, x) if np.any(B) else 0))
        return F, dF

    def inward_power(self, r, l=None, per_order=False):
        """Net power per unit length flowing inward through the circle of radius r (W/m)."""
        F, dF = self.axial_and_deriv(r, l)
        w = self.omega
        if self.pol == "TM":
            Hphi = dF / (1j * w * MU0)
            p = np.pi * r * np.real(F * np.conj(Hphi))
        else:
            lay = self.layer_of(r) if l is None else l
            Ephi = -dF / (1j * w * EPS0 * self.eps[lay])
            p = -np.pi * r * np.real(Ephi * np.conj(F))
        return p if per_order else float(np.sum(p))

    def absorbed_per_layer(self, per_order=False):
        """Absorbed power per unit length in each finite layer (W/m), core first.

        per_order=True returns shape (n_layers-1, n_orders) aligned with self.ms.
        """
        out = []
        prev = 0.0
        for i, r in enumerate(self.radii):
            p_in = self.inward_power(r, l=i, per_order=per_order)
            out.append(p_in - prev)
            prev = self.inward_power(r, l=i + 1, per_order=per_order)
        return np.array(out)

    @property
    def I0(self):
        nh = refractive_index(self.eps[-1])
        if self.pol == "TM":
            return 0.5 * nh.real / ETA0  # |E0| = 1
        return 0.5 * nh.real / ETA0  # TE also normalised to |E0| = 1

    def c_abs_per_layer(self):
        """Absorption cross-section per unit length for each layer (m)."""
        return self.absorbed_per_layer() / self.I0

    def field_E2(self, r, phi):
        """|E|^2 at (r, phi) (for checks / plots)."""
        l = self.layer_of(r)
        F, dF = self.axial_and_deriv(r, l)
        e = np.exp(1j * self.ms * phi)
        if self.pol == "TM":
            return abs(np.sum(F * e)) ** 2
        w = self.omega
        eps = self.eps[l]
        Er = np.sum((1j * self.ms / r) * F * e) / (1j * w * EPS0 * eps)
        Ephi = np.sum(-dF * e) / (1j * w * EPS0 * eps)
        return abs(Er) ** 2 + abs(Ephi) ** 2


def solve(f_hz: float, radii, eps, pol: str = "TM", m_max: int | None = None) -> CylSolution:
    radii = np.asarray(radii, dtype=float)
    eps = np.asarray(eps, dtype=complex)
    assert len(eps) == len(radii) + 1
    k = 2 * np.pi * f_hz / C0 * refractive_index(eps)
    if m_max is None:
        x = abs(k[-1]) * radii[-1]
        m_max = int(max(4, x + 4 * x ** (1 / 3) + 4))
    ms = np.arange(-m_max, m_max + 1)
    nl = len(eps)
    coef = np.zeros((len(ms), nl, 2), dtype=complex)
    nh = refractive_index(eps[-1])
    # incident amplitude of the axial field per order (|E0| = 1 at the axis)
    a_inc0 = 1.0 if pol == "TM" else nh / ETA0  # Hz = E0 n / eta0
    for im, m in enumerate(ms):
        p = (lambda l: 1.0) if pol == "TM" else (lambda l: 1.0 / eps[l])

        def W(l, r):
            x = k[l] * r
            return np.array([[jv(m, x), hankel2(m, x)],
                             [p(l) * k[l] * jvp(m, x), p(l) * k[l] * h2vp(m, x)]], dtype=complex)

        vecs = [np.array([1.0 + 0j, 0.0 + 0j])]
        for i, r in enumerate(radii):
            vecs.append(np.linalg.solve(W(i + 1, r), W(i, r) @ vecs[-1]))
        a_inc = a_inc0 * (-1j) ** m
        scale = a_inc / vecs[-1][0]
        for l in range(nl):
            coef[im, l] = vecs[l] * scale
    return CylSolution(f_hz, pol, radii, eps, ms, coef)
