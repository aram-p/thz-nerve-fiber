import numpy as np

from thznerve.physics.scattering import foldy_eps, q_ext_q_sca


def test_mie_lossless_sphere_matches_reference():
    q_ext, q_sca = q_ext_q_sca(1.5, 1.0)
    assert abs(q_ext - 0.21510) < 1e-4 and abs(q_ext - q_sca) < 1e-12      # Bohren-Huffman, m = 1.5, x = 1


def test_absorbing_sphere_extinguishes_more_than_it_scatters():
    q_ext, q_sca = q_ext_q_sca(1.5 + 0.1j, 2.0)
    assert q_ext > q_sca > 0


def test_foldy_small_spheres_is_maxwell_garnett_to_first_order():
    eh, es, f = 4.0 - 1.0j, 1.0, 1e-3
    beta = (es - eh) / (es + 2 * eh)
    mg_first = eh * (1 + 3 * f * beta)
    got = foldy_eps(0.3e12, eh, es, 0.5, f)
    assert abs(got - mg_first) < 0.02 * abs(mg_first - eh)


def test_foldy_zero_fraction_is_host():
    assert foldy_eps(1e12, 4.0 - 1.0j, 1.0, 50, 0.0) == 4.0 - 1.0j
