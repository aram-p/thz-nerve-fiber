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


def test_second_pulse_arrival_bound_holds():
    # the second pulse placed at the tissue delay (0.6 ps late) fits best unbounded; with t2_max it must stay near 0
    tt = np.arange(0, 8, 0.005)
    B, K = _pulse(tt, 4.0), _pulse(tt, 3.4)
    tp = TwoPulse(tt, K)
    y = 0.8 * B + 0.1 * tp.shift(K, 0.6) + 0.05 * K
    for s in (None, 0.0):
        p, _, _ = tp.fit(y, B, s=s, p0=[0.8, 0.1, 0.6] if s == 0.0 else [0.8, 0.1, 0.0, 0.6], t2_max=0.2)
        assert abs((p[2] + p[3]) if s is None else p[2]) <= 0.2 + 1e-3
    p, _, _ = tp.fit(y, B, s=0.0, p0=[0.8, 0.1, 0.6])
    assert abs(p[2] - 0.6) < 0.05


def test_multistart_keeps_fixed_s_with_loss():
    from thznerve.waveform_fit import multistart
    tt = np.arange(0, 8, 0.005)
    B, K = _pulse(tt, 4.0), _pulse(tt, 3.4)
    tp = TwoPulse(tt, K)
    y = 0.7 * tp.shift(tp.apply_transfer(B, tp.loss_transfer(4.0, "minphase")), 0.05) + 0.15 * tp.shift(K, 0.05)
    fit = multistart(lambda y_, x_, **kw: tp.fit(y_, x_, a_max=1.0, **kw))
    q, _, _ = fit(y, B, s=0.0, loss="minphase")
    assert len(q) == 4                                   # (a, b, tau, kappa): s stays fixed
    q, m, _ = fit(y, B, loss="minphase")
    assert len(q) == 5 and m < 1e-6


def test_node_transfer_reproduces_solved_points():
    from thznerve.waveform_fit import node_transfer
    st = {"f_thz": [0.6, 2.0], "loss_db_dipole": [0.5, 2.0], "delay_rad": [0.05, 0.3]}
    r = node_transfer(st, np.array([0.0, 0.6, 2.0]), delay_factor="loss")
    assert abs(r[0] - 1) < 1e-12
    assert abs(-20 * np.log10(abs(r[1])) - 0.5) < 1e-9 and abs(-np.angle(r[1]) - 0.05) < 1e-9
    assert abs(-20 * np.log10(abs(r[2])) - 2.0) < 1e-9 and abs(-np.angle(r[2]) - 0.3) < 1e-9
    r = node_transfer(st, np.array([0.6, 2.0]))                     # delay with the 2.2 path factor at 2 THz
    assert abs(-np.angle(r[0]) - 0.05) < 1e-9 and abs(-np.angle(r[1]) - 0.3 * 2.2 / 3.1) < 1e-9


def test_loss_term_recovered_and_minphase_is_causal():
    tt = np.arange(0, 8, 0.005)
    B, K = _pulse(tt, 4.0), _pulse(tt, 3.4)
    tp = TwoPulse(tt, K)
    y = 0.7 * tp.shift(tp.apply_transfer(B, tp.loss_transfer(6.0)), 0.04) + 0.15 * tp.shift(K, 0.0)
    p, misfit, _ = tp.fit(y, B, loss="zero")
    assert abs(p[2] - 0.04) < 2e-3 and abs(p[4] - 6.0) < 0.3
    h = tp.loss_transfer(6.0, "minphase")
    assert np.allclose(np.abs(h), tp.loss_transfer(6.0), rtol=1e-6)
    imp = np.fft.irfft(h, tp.np_)
    assert np.sum(imp[tp.np_ // 2:] ** 2) < 1e-6 * np.sum(imp ** 2)      # no response before t = 0


def test_loss_cutoff_no_gain_and_ab_bound():
    tt = np.arange(0, 8, 0.005)
    B, K = _pulse(tt, 4.0), _pulse(tt, 3.4)
    tp = TwoPulse(tt, K, f_loss=1.5)
    g = tp.loss_transfer(6.0)
    assert np.all(g <= 1 + 1e-12) and np.allclose(g[tp.fr <= 0.2], 1.0)
    assert np.allclose(g[tp.fr >= 1.5], 10 ** (-6.0 * 1.3 / 20))
    y = 0.9 * B + 0.3 * K
    p, _, _ = tp.fit(y, B, s=0.0, a_max=1.0, ab_max=1.0)
    assert p[0] + p[1] <= 1.0 + 1e-4
    F, pv = f_test(2.0, 1.0, 22.0, n_par=6, n_fixed=2)
    assert abs(F - 0.5 * 16) < 1e-12 and 0 < pv < 0.01


def test_mixed_loss_is_floor_times_causal_increment():
    tt = np.arange(0, 8, 0.005)
    tp = TwoPulse(tt, _pulse(tt, 3.4), kappa0=2.0)
    m = tp.loss_transfer(3.0, "mixed")
    assert np.allclose(np.abs(m), tp.loss_transfer(5.0), rtol=1e-6)
    assert np.allclose(np.angle(m), np.angle(tp.loss_transfer(3.0, "minphase")), atol=1e-9)


def test_supply_test_passes_true_cause_and_fails_none():
    from thznerve.waveform_fit import supply_test
    tt = np.arange(0, 8, 0.005)
    B, K = _pulse(tt, 4.0), _pulse(tt, 3.4)
    tp = TwoPulse(tt, K)
    tr = lambda v: np.exp(-v * np.clip(tp.fr, 0, 4) ** 2)          # a steep loss, zero phase
    y = 0.7 * tp.apply_transfer(B, tr(0.3)) + 0.15 * tp.shift(K, -0.03)
    y = y + np.random.default_rng(1).normal(scale=1e-3 * np.abs(y).max(), size=len(y))
    fit = lambda yy, x, **kw: tp.fit(yy, x, a_max=1.0, **kw)
    pr, s_ref, res = fit(y, B, loss="minphase")
    ne, _ = n_eff(res)
    q0, _, _ = fit(y, B, s=0.0)
    good = supply_test(tp, fit, y, B, tr, np.linspace(0, 0.6, 7), q0, pr, s_ref, ne)
    assert good["verdict"] == "pass" and abs(good["best"] - 0.3) < 1e-9
    flat = supply_test(tp, fit, y, B, lambda v: np.full(len(tp.fr), np.exp(-v)), np.linspace(0, 0.6, 7),
                       q0, pr, s_ref, ne)
    assert flat["verdict"] == "fail"


def test_null_calibrated_test_accepts_null_and_rejects_change():
    from thznerve.waveform_fit import null_calibrated_test
    tt = np.arange(0, 6, 0.01)
    B, K = _pulse(tt, 3.0), _pulse(tt, 2.5)
    tp = TwoPulse(tt, K, f_loss=2.0)
    fit = lambda y, x, **kw: tp.fit(y, x, a_max=1.0, **kw)
    model = lambda q, x, s: tp.model(q, x, s=s)
    rng = np.random.default_rng(1)
    noise = np.convolve(rng.normal(size=len(tt)), np.ones(5) / 5, "same") * 0.01
    y0 = 0.8 * B + 0.1 * tp.shift(K, -0.02) + noise                       # no change of the tissue pulse
    d0 = null_calibrated_test(fit, y0, B, [B], 0.0, 2, nb=40, rng=rng, model=model)
    assert d0["p_cal"] > 0.05
    lossy = tp.apply_transfer(B, tp.loss_transfer(8.0, "minphase"))
    y1 = 0.8 * lossy + 0.1 * tp.shift(K, -0.02) + noise                    # a causal loss and its delay
    d1 = null_calibrated_test(fit, y1, B, [B], 0.0, 2, nb=40, rng=rng, model=model)
    assert d1["p_cal"] < 0.05 and d1["p_nominal"] < d0["p_nominal"]


def test_multistart_never_worse_than_one_start():
    from thznerve.waveform_fit import multistart
    tt = np.arange(0, 6, 0.01)
    B, K = _pulse(tt, 3.0), _pulse(tt, 2.5)
    tp = TwoPulse(tt, K, f_loss=2.0)
    one = lambda y, x, **kw: tp.fit(y, x, a_max=1.0, **kw)
    fit = multistart(one)
    y = 0.8 * tp.apply_transfer(B, tp.loss_transfer(6.0, "minphase")) + 0.1 * tp.shift(K, -0.02)
    for kw in (dict(s=0.0), dict(), dict(loss="minphase")):
        assert fit(y, B, **kw)[1] <= one(y, B, **kw)[1] + 1e-12
    q, m, _ = fit(y, B, loss="minphase")
    assert len(q) == 5 and m < 1e-3


def test_null_calibrated_test_stops_early_when_clear():
    from thznerve.waveform_fit import null_calibrated_test
    tt = np.arange(0, 6, 0.01)
    B, K = _pulse(tt, 3.0), _pulse(tt, 2.5)
    tp = TwoPulse(tt, K, f_loss=2.0)
    fit = lambda y, x, **kw: tp.fit(y, x, a_max=1.0, **kw)
    model = lambda q, x, s: tp.model(q, x, s=s)
    rng = np.random.default_rng(2)
    noise = np.convolve(rng.normal(size=len(tt)), np.ones(5) / 5, "same") * 0.01
    y1 = 0.8 * tp.apply_transfer(B, tp.loss_transfer(8.0, "minphase")) + 0.1 * tp.shift(K, -0.02) + noise
    d = null_calibrated_test(fit, y1, B, [B], 0.0, 2, nb=200, rng=rng, model=model, batch=20)
    assert d["nb"] < 200 and d["p_cal"] < 0.05


def test_f_test_without_dof_has_no_power():
    from thznerve.waveform_fit import f_test
    assert f_test(2.0, 1.0, 5.0, 6, 1)[1] == 1.0 and f_test(2.0, 1.0, 50.0, 6, 1)[1] < 0.01


def test_phase_surrogate_keeps_spectrum():
    from thznerve.waveform_fit import phase_surrogate
    rng = np.random.default_rng(3)
    r = np.convolve(rng.normal(size=600), np.ones(9) / 9, "same")
    s = phase_surrogate(rng, r)
    assert np.allclose(np.abs(np.fft.rfft(s)), np.abs(np.fft.rfft(r)), rtol=1e-9, atol=1e-12)
    assert not np.allclose(s, r) and abs(n_eff(s)[0] / n_eff(r)[0] - 1) < 0.35


def test_size_fitted_test_passes_true_cause_and_fails_wrong_one():
    from thznerve.waveform_fit import multistart, size_fitted_test
    tt = np.arange(0, 6, 0.01)
    B, K = _pulse(tt, 3.0), _pulse(tt, 2.5)
    tp = TwoPulse(tt, K, f_loss=2.0)
    one = lambda y, x, **kw: tp.fit(y, x, a_max=1.0, **kw)
    model = lambda q, x, s: tp.model(q, x, s=s)
    rng = np.random.default_rng(4)
    noise = np.convolve(rng.normal(size=len(tt)), np.ones(5) / 5, "same") * 0.01
    lossy = tp.apply_transfer(B, tp.loss_transfer(8.0, "minphase"))
    y = 0.8 * lossy + 0.1 * tp.shift(K, -0.02) + noise
    true = [tp.apply_transfer(B, tp.loss_transfer(k, "minphase")) for k in (0.0, 4.0, 8.0, 12.0)]
    flat = [B * np.exp(-v) for v in (0.0, 0.2, 0.4)]                 # a frequency-flat loss cannot give the change
    d1 = size_fitted_test(multistart(one), one, y, true, 0.0, 2, model, rng=rng, nb_max=100, batch=25)
    d2 = size_fitted_test(multistart(one), one, y, flat, 0.0, 2, model, rng=rng, nb_max=100, batch=25)
    assert d1["verdict"] != "fail" and d1["size_index_null"] == 2
    assert d2["verdict"] == "fail"


def test_effect_noise_small_when_restricted_model_true():
    from thznerve.waveform_fit import effect_noise, effect_size
    rng = np.random.default_rng(1)
    t = np.linspace(0, 10, 400)
    x1, x2 = np.sin(t), np.exp(-(t - 5) ** 2)
    lsq = lambda X: (lambda y: (None, float(np.sum((y - X @ np.linalg.lstsq(X, y, rcond=None)[0]) ** 2))))
    fit_r, fit_f = lsq(x1[:, None]), lsq(np.c_[x1, x2])
    res = 0.05 * rng.standard_normal(len(t))
    d0 = effect_noise(fit_r, fit_f, 2 * x1, res, 40, rng)
    y = 2 * x1 + 0.5 * x2 + res
    d1 = effect_size(fit_r(y)[1], fit_f(y)[1], y)
    assert np.quantile(d0, 0.95) < 0.3 * d1
    assert effect_size(1.0, 2.0, y) == 0.0


def test_noise_scale_recovers_absorbed_noise():
    from thznerve.waveform_fit import noise_scale
    rng = np.random.default_rng(2)
    n, k = 60, 20
    X = rng.standard_normal((n, k))
    fit_f = lambda y: (None, float(np.sum((y - X @ np.linalg.lstsq(X, y, rcond=None)[0]) ** 2)))
    res = rng.standard_normal(n)
    s = noise_scale(fit_f, X[:, 0], res, 200, rng)
    assert abs(s - np.sqrt(n / (n - k))) < 0.05 * np.sqrt(n / (n - k))


def test_threshold_verdict_false_rate_near_level():
    from thznerve.waveform_fit import threshold_verdict
    rng = np.random.default_rng(3)
    draws = rng.standard_normal((3, 2000)) ** 2
    v = threshold_verdict([0.1, 0.1, 0.1], draws, 0.0, 1 - 0.05 / 3, nboot=20, rng=rng)
    assert not v["excluded"] and v["flip"] == 0.0
    assert 0.03 < v["false_rate"] < 0.07
    v = threshold_verdict([0.1, 50.0, 0.1], draws, 1.0, 1 - 0.05 / 3, nboot=20, rng=rng)
    assert v["excluded"] and v["largest"] > 1
