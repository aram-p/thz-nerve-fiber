import numpy as np

from thznerve.waveform_fit import TwoPulse, f_test, n_eff


def _pulse(tt, t0, w=0.15):
    x = (tt - t0) / w
    return -x * np.exp(-x ** 2)


def test_fit_recovers_shift_and_second_pulse():
    tt = np.arange(0, 8, 0.005)
    B, K = _pulse(tt, 4.0), _pulse(tt, 3.4)
    tp = TwoPulse(tt, K)
    y = 0.7 * tp.shift(B, 0.05) + 0.15 * tp.shift(K, 0.05 - 0.03)
    p, misfit, _ = tp.fit(y, B)
    assert np.allclose(p, [0.7, 0.15, 0.05, -0.03], atol=2e-3)
    assert misfit < 1e-8


def test_apply_transfer_pure_delay_equals_shift():
    tt = np.arange(0, 8, 0.005)
    B = _pulse(tt, 4.0)
    tp = TwoPulse(tt, B)
    d = 0.05                                             # ps
    moved = tp.apply_transfer(B, np.exp(-2j * np.pi * tp.fr * d))
    assert np.max(np.abs(moved - tp.shift(B, d))) < 2e-3 * np.max(np.abs(B))


def test_n_eff_white_and_f_test():
    r = np.random.default_rng(0).normal(size=4000)
    ne, _ = n_eff(r)
    assert 0.7 * 4000 < ne < 1.3 * 4000
    F, p = f_test(1.0, 1.0, 20.0)
    assert F == 0.0 and p == 1.0


def test_a_max_bound_holds():
    tt = np.arange(0, 8, 0.005)
    B, K = _pulse(tt, 4.0), _pulse(tt, 3.4)
    tp = TwoPulse(tt, K)
    y = 1.5 * B + 0.1 * K
    p, _, _ = tp.fit(y, B, s=0.0, a_max=1.0)
    assert p[0] <= 1.0 + 1e-9


def test_node_transfer_reproduces_solved_points():
    from thznerve.waveform_fit import node_transfer
    st = {"f_thz": [0.6, 2.0], "loss_db_dipole": [0.5, 2.0], "delay_rad": [0.05, 0.3]}
    r = node_transfer(st, np.array([0.0, 0.6, 2.0]))
    assert abs(r[0] - 1) < 1e-12
    assert abs(-20 * np.log10(abs(r[1])) - 0.5) < 1e-9 and abs(-np.angle(r[1]) - 0.05) < 1e-9
    assert abs(-20 * np.log10(abs(r[2])) - 2.0) < 1e-9 and abs(-np.angle(r[2]) - 0.3) < 1e-9
