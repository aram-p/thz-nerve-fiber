import numpy as np

from thznerve.physics.effmedium import coated_cylinder_eps, tissue_eps


def test_limits():
    e = 4 - 2j
    assert abs(coated_cylinder_eps(e, e, 0.7) - e) < 1e-12          # homogeneous
    assert abs(coated_cylinder_eps(3, 2, 1.0) - 3) < 1e-12           # all core
    assert abs(coated_cylinder_eps(3, 2, 0.0) - 2) < 1e-12           # all shell
    par, perp = tissue_eps(e, e, e)
    assert abs(par - e) < 1e-12 and abs(perp - e) < 1e-12
    par, perp = tissue_eps(2.0, 2.0, 4.0, p=0.0)
    assert abs(par - 4) < 1e-12 and abs(perp - 4) < 1e-12


def test_wiener_order():
    # for real permittivities the perpendicular value lies below the parallel one (Wiener bounds)
    par, perp = tissue_eps(3.9, 2.0, 4.1)
    assert perp.real < par.real
