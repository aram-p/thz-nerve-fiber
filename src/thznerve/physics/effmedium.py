"""Effective permittivity of tissue made of aligned myelinated fibres in water (quasi-static).

Fibres: axon (eps_ax) inside myelin (eps_my), g = a/b, packed at area fraction p in
extracellular water (eps_w). The fibres are much thinner than the THz wavelength in water
(kb < 0.7 up to 3 THz for b <= 5.5 um), so mixing rules apply.

E parallel to the fibres: fields are tangential to every interface, so eps is the area average
(Wiener upper bound, exact for parallel cylinders).
E perpendicular: 2D Maxwell Garnett with the coated-cylinder equivalent permittivity
(exact for dilute coated cylinders, a good approximation up to p ~ 0.6).
Nodes occupy a volume fraction ~1e-4 to 1e-3 and are ignored here (see samplebound).
"""

from __future__ import annotations

import numpy as np


def coated_cylinder_eps(eps_core, eps_shell, g):
    """2D quasi-static equivalent permittivity of a core-shell cylinder (core radius / outer = g)."""
    q = g ** 2
    num = (eps_core + eps_shell) + q * (eps_core - eps_shell)
    den = (eps_core + eps_shell) - q * (eps_core - eps_shell)
    return eps_shell * num / den


def tissue_eps(eps_ax, eps_my, eps_w, g=0.7, p=0.6):
    """(eps_parallel, eps_perpendicular) of aligned fibres at area fraction p."""
    q = g ** 2
    e_par = p * (q * eps_ax + (1 - q) * eps_my) + (1 - p) * eps_w
    ec = coated_cylinder_eps(eps_ax, eps_my, g)
    e_perp = eps_w * ((ec + eps_w) + p * (ec - eps_w)) / ((ec + eps_w) - p * (ec - eps_w))
    return e_par, e_perp
