import numpy as np

from thznerve.physics.materials import WATER_KS1996, axoplasm, myelin
from thznerve.physics.nodelaw import node_factor_closed
from thznerve.physics.quasistatic import node_factor


def test_closed_law_vs_laplace_solver():
    for f in (0.3e12, 2e12):                         # not in the fit set (0.1, 1, 3 THz)
        e = [complex(axoplasm()(f)), complex(myelin()(f)), complex(WATER_KS1996(f))]
        for a, b, Ln in ((3.5, 5.0, 1.0), (2.0, 3.0, 1.5), (5.0, 7.0, 0.8)):
            for pol, field, tol in (("par", "axial", 0.08), ("perp", "transverse", 0.05)):
                F = node_factor(a, b, Ln, e[0], e[1], e[2], e[2], field=field)
                assert abs(node_factor_closed(a, b, Ln, *e, pol=pol) / F - 1) < tol


def test_limits():
    e = [4 - 2j, 3 - 0.7j, 4.1 - 2.1j]
    assert abs(node_factor_closed(3.5, 5, 1e4, *e) - 1) < 1e-3
    assert abs(node_factor_closed(3.5, 5, 1e-6, *e) - abs(e[1] / e[2]) ** 2) < 1e-3
