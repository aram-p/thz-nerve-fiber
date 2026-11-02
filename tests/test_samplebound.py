import numpy as np

from thznerve.physics.materials import C0, EPS0
from thznerve.physics.samplebound import delta_db, node_volume_fraction, required_delta_sigma


def test_roundtrip_and_bulk_limit():
    f, F, fn, n, d = 0.6e12, 0.6, 2.5e-4, 2.25, 100e-6
    s = required_delta_sigma(f, 9.0, F, fn, n, d)
    de = s / (2 * np.pi * f * EPS0)
    assert abs(delta_db(f, de, F, fn, n, d) - 9.0) < 1e-9
    # f_n = F = 1: a uniform eps'' change of a medium with index n gives alpha = k0 eps''/n
    k0 = 2 * np.pi * f / C0
    assert abs(delta_db(f, 0.1, 1.0, 1.0, n, d) - 4.343 * k0 * 0.1 / n * d) < 1e-3


def test_node_fraction():
    # node pi (25 - 12.25) um^3, cell 1001 um * pi 25 / 0.5 um^2
    assert abs(node_volume_fraction(3.5, 5.0, 1.0, 1000.0) - 12.75 / 50050) < 1e-12
