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


def test_absorbing_host_uses_sphere_wavenumber_inside():
    # a sphere of the host's own material changes nothing, at any size and host loss
    from thznerve.physics.scattering import mie_ab
    nh = 2.1 + 0.6j
    _, a, b = mie_ab(1.0 + 0j, (2 * np.pi * 1e12 / 3e8) * nh * 40e-6)
    assert np.max(np.abs(a)) < 1e-12 and np.max(np.abs(b)) < 1e-12


def test_absorbing_host_lossless_limit_continuous():
    # as the host loss goes to zero the complex-x result joins the real-x one
    a = foldy_eps(0.6e12, 4.4 - 1e-9j, 1.0, 30, 0.1)
    b = foldy_eps(0.6e12, 4.4 + 0j, 1.0, 30, 0.1)
    assert abs(a - b) < 1e-6
