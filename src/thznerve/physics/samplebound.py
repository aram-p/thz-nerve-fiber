"""How much must the nodes change to move a tissue slab's THz transmission by a given dB?

A slab of thickness d contains nodes with volume fraction f_n. A change d_eps'' of the node
permittivity adds absorbed power per volume 1/2 w eps0 d_eps'' F f_n |E|^2, where F =
<|E|^2>_node / |E|^2 is the node field factor from the full-wave fibre model. Divided by the
intensity 1/2 n |E|^2 / eta0 this is an extra attenuation coefficient

    d_alpha = k0 d_eps'' F f_n / n        (power, 1/m)

and the transmission change is 10 log10(e) d_alpha d dB. Inverting gives the node change a
reported dip needs. First order in f_n (dilute nodes); local-field corrections in dense
tissue change F by O(1), not by orders of magnitude.
"""

from __future__ import annotations

import numpy as np

from thznerve.physics.materials import C0, EPS0

DB_PER_NEPER_POWER = 10 / np.log(10)   # 4.343 dB per unit of alpha*d


def node_volume_fraction(a_um, b_um, Ln_um, Lint_um, fibre_area_fraction=0.5):
    """Node volume / tissue volume for fibres of outer radius b packed at the given area fraction."""
    v_node = np.pi * (b_um ** 2 - a_um ** 2) * Ln_um
    v_cell = (Lint_um + Ln_um) * np.pi * b_um ** 2 / fibre_area_fraction
    return v_node / v_cell


def delta_db(f_hz, d_eps_im, F, f_n, n, d_m):
    """Transmission change (dB, positive = more loss) from a node eps'' change d_eps_im."""
    k0 = 2 * np.pi * np.asarray(f_hz) / C0
    return DB_PER_NEPER_POWER * k0 * d_eps_im * F * f_n / n * d_m


def required_delta_sigma(f_hz, db, F, f_n, n, d_m):
    """Node conductivity change (S/m) needed for a transmission change of `db` dB."""
    k0 = 2 * np.pi * np.asarray(f_hz) / C0
    d_eps_im = db * n / (DB_PER_NEPER_POWER * k0 * F * f_n * d_m)
    return 2 * np.pi * np.asarray(f_hz) * EPS0 * d_eps_im
