import numpy as np

from thznerve.physics.materials import WATER_KS1996, axoplasm, myelin
from thznerve.physics.nodelaw import shell_factor
from thznerve.physics.packing import layered_cylinder_factor, periodic_shell_factor_extrap


def _eps(f):
    return complex(axoplasm()(f)), complex(myelin()(f)), complex(WATER_KS1996(f))


def test_layered_matches_shell_factor():
    for f in (0.1e12, 1e12, 3e12):
        ea, em, ew = _eps(f)
        assert abs(layered_cylinder_factor([3.5, 5.0], [ea, em], ew, 1) / shell_factor(ea, em, ew, 0.7) - 1) < 1e-6


def test_layered_null_and_core():
    # all layers = host: field uniform, factor 1 in every layer
    assert abs(layered_cylinder_factor([1.0, 2.0, 3.0], [4 - 1j] * 3, 4 - 1j, 1) - 1) < 1e-9
    # bare cylinder core: |2 eh / (ec + eh)|^2
    ec, eh = 2.3 - 0.05j, 4.1 - 2.1j
    assert abs(layered_cylinder_factor([1.0], [ec], eh, 0) - abs(2 * eh / (ec + eh)) ** 2) < 1e-6


def test_periodic_dilute_limit():
    ea, em, ew = _eps(1e12)
    v = periodic_shell_factor_extrap(3.5, 5.0, ea, em, ew, 0.02, N=200)
    assert abs(v / shell_factor(ea, em, ew, 0.7) - 1) < 0.03
