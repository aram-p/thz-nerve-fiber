"""Closed-form law for THz absorption at a node of Ranvier (quasi-static, non-resonant).

    P_node = 1/2 w eps0 eps''_w |E0|^2 V_node F,   V_node = pi (b^2 - a^2) Ln

F interpolates between two exact limits in x = Ln / t (t = b - a, myelin thickness):
  * E parallel to the fibre. Thin gap (x -> 0): normal D is continuous across the gap faces and
    the myelin carries the outside field, so F0 = |eps_my / eps_w|^2. Long gap (x -> inf): the gap
    is bulk water in a tangential field, F = 1.
  * E perpendicular. Thin gap: the gap takes the tangential field of the myelin shell of an
    infinite coated cylinder, F0 = S (exact 2D series, `shell_factor`). Long gap: F = 1.
    F = F0 + (1 - F0) x^p / (x^p + c),   c = c0 |eps_w / eps_my|^q
The transition length depends on the water/myelin contrast: E parallel needs a longer gap when
water is more polarisable than myelin (q > 0), E perpendicular a shorter one (q < 0).
The three constants (c0, p, q) per polarisation were fitted once to the Laplace solver
(thznerve.physics.quasistatic, Richardson-extrapolated) over x = 0.2-8, g = 0.6-0.8,
0.1-3 THz (experiments/prod/fit_nodelaw.py): max error 2.1 % (parallel), 2.0 % (perpendicular),
rms 0.7 % / 0.5 %. They are not fitted to full-wave results or to any experiment.
The Laplace problem has no length scale, so F depends on x, g and the permittivities only.
No length in the law compares with the wavelength, and the internode length does not enter:
there is no resonance condition.
"""

from __future__ import annotations

import numpy as np

from thznerve.physics.materials import EPS0

C_PAR, P_PAR, Q_PAR = 1.066, 0.934, 0.742
C_PERP, P_PERP, Q_PERP = 0.999, 0.826, -0.784


def shell_factor(eps_ax, eps_my, eps_w, g):
    """<|E|^2>/|E0|^2 in the myelin shell of an infinite coated cylinder, E perpendicular (2D)."""
    a, b = g, 1.0
    M = np.array([[a, -a, -1 / a, 0], [eps_ax, -eps_my, eps_my / a ** 2, 0], [0, b, 1 / b, -1 / b],
                  [0, eps_my, -eps_my / b ** 2, eps_w / b ** 2]], complex)
    _, B, C, _ = np.linalg.solve(M, np.array([0, 0, -b, -eps_w], complex))
    r = np.linspace(a, b, 401)
    return float(np.trapezoid((abs(B) ** 2 + abs(C) ** 2 / r ** 4) * r, r) / np.trapezoid(r, r))


def node_factor_closed(a, b, Ln, eps_ax, eps_my, eps_w, pol="par"):
    x = Ln / (b - a)
    rho = abs(eps_w / eps_my)
    if pol == "par":
        F0, c, p = abs(eps_my / eps_w) ** 2, C_PAR * rho ** Q_PAR, P_PAR
    else:
        F0, c, p = shell_factor(eps_ax, eps_my, eps_w, a / b), C_PERP * rho ** Q_PERP, P_PERP
    return F0 + (1 - F0) * x ** p / (x ** p + c)


def node_power(f_hz, a_um, b_um, Ln_um, eps_ax, eps_my, eps_w, E0=1.0, pol="par"):
    """Absorbed power per node (W) for incident field amplitude E0 (V/m) at the fibre axis."""
    V = np.pi * (b_um ** 2 - a_um ** 2) * Ln_um * 1e-18
    F = node_factor_closed(a_um, b_um, Ln_um, eps_ax, eps_my, eps_w, pol)
    return 0.5 * 2 * np.pi * f_hz * EPS0 * (-np.imag(eps_w)) * abs(E0) ** 2 * V * F
