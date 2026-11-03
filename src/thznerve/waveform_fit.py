"""Two-pulse fits of the traced published waveforms (Sims 47-51, Sec. 5.6).

    biased(t) = a B(t - s) + b K(t - s - tau)

B the unbiased (tissue) pulse traced from the same figure as the biased one, K the empty-cuvette pulse
(H1 Fig. 3a). A candidate change of the slice is applied to B through its transfer ratio on the FFT grid
(`apply_transfer`), and the fit with s = 0 then tells how much of what a free s buys the change recovers.
The residual is correlated over ~0.3 ps, so significance uses N_eff (`n_eff`) less the fitted parameters.
"""
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import f as fdist


def load_scans(root):
    """Traced pulses on one time grid: (t_ps, K, {scan: (biased, unbiased)})."""
    root = Path(root)
    h1 = np.load(root / "results" / "sim42" / "digitised.npz")
    h3 = np.load(root / "results" / "sim39" / "digitised.npz")
    p2 = np.load(root / "results" / "sim43" / "paper2_2min.npz")
    tt = h1["t_ps"]
    h3i = {k: np.interp(tt, h3["t_ps"], h3["wave_" + k]) for k in ("blue", "red")}
    scans = {"30 s (H1)": (h1["30s"], h1["none"]), "2 min (H2)": (p2["biased"], p2["none"]),
             "4 min (H1)": (h1["4min"], h1["none"]), "5 min (H3)": (h3i["red"], h3i["blue"])}
    return tt, h1["empty"], scans


class TwoPulse:
    def __init__(self, tt, K, pad=4, f_loss=4.0, kappa0=0.0):
        self.tt, self.K = tt, K
        self.kappa0 = kappa0                                  # dB/THz: zero-phase floor of the 'mixed' loss
        self.F_LOSS = f_loss                                  # THz: the loss term is held constant above
        self.dt = tt[1] - tt[0]
        self.n = len(tt)
        self.np_ = pad * self.n
        self.fr = np.fft.rfftfreq(self.np_, self.dt)          # THz when t is in ps

    def shift(self, y, tau):
        return np.interp(self.tt - tau, self.tt, y, left=0.0, right=0.0)

    F_LO = 0.2                                                # THz: the loss term is 1 here (band edge), so a <= 1 keeps the
                                                              # tissue pulse at or below the unbiased one over the band

    def loss_transfer(self, kappa, kind="zero"):
        """Smooth extra loss of the tissue pulse, kappa dB per THz of field amplitude above F_LO (1 below
        F_LO, held above F_LOSS).
        kind 'zero': zero phase (a smoothing, e.g. of the traced curve); 'minphase': the causal (minimum-phase)
        transfer with the same magnitude, from the folded cepstrum, which adds its own delay; 'mixed': a zero-phase
        floor of kappa0 times a causal increment of kappa."""
        ramp = (np.clip(self.fr, self.F_LO, self.F_LOSS) - self.F_LO) * np.log(10) / 20
        lnmag = -kappa * ramp
        if kind == "zero":
            return np.exp(lnmag)
        if kind == "mixed":
            return np.exp(-self.kappa0 * ramp) * self.loss_transfer(kappa, "minphase")
        return np.exp(kappa * self._minphase_log())

    def _minphase_log(self):
        """Complex log of the causal transfer for kappa = 1 (the folded cepstrum is linear in ln|T|, so the
        transfer for any kappa is exp(kappa * this)); computed once."""
        if getattr(self, "_mp_log", None) is None:
            lnmag = -(np.clip(self.fr, self.F_LO, self.F_LOSS) - self.F_LO) * np.log(10) / 20
            full = np.r_[lnmag, lnmag[-2:0:-1]]
            c = np.fft.ifft(full).real
            m = len(full) // 2
            fold = np.r_[c[0], 2 * c[1:m], c[m], np.zeros(len(full) - m - 1)]
            self._mp_log = np.fft.fft(fold)[: len(self.fr)]
        return self._mp_log

    def model(self, p, x0, s=None, loss=None):
        """p = (a, b, [s], tau, [kappa]): s omitted when fixed, kappa present when loss is 'zero' or 'minphase'."""
        p = list(p)
        kap = p.pop() if loss else 0.0
        a, b, s_, tau = p if s is None else (p[0], p[1], s, p[2])
        xt = self.apply_transfer(x0, self.loss_transfer(kap, loss)) if loss else x0
        return a * self.shift(xt, s_) + b * self.shift(self.K, s_ + tau)

    def fit(self, y, x0, s=None, p0=None, a_max=np.inf, loss=None, ab_max=None):
        """Least squares over (a, b, s, tau), or (a, b, tau) with s fixed, plus kappa >= 0 when loss is given;
        a <= a_max (a passive slice cannot raise the tissue pulse above what its own transfer gives: a_max = 1);
        a + b <= ab_max if given (the second pulse takes its area from the tissue path), as a stiff penalty.
        Returns (params, misfit, residual)."""
        if p0 is not None:
            starts = [list(p0)]
        elif s is None:
            starts = [[0.6, b0, s0, 0.0] for b0 in (0.05, 0.2) for s0 in (-0.1, 0.0, 0.1)]
        else:
            starts = [[0.6, b0, t0] for b0 in (0.05, 0.2) for t0 in (-0.1, 0.0, 0.1)]
        if loss and p0 is None:
            starts = [st + [k0] for st in starts for k0 in (1.0, 6.0)]
        w = 1e3 * np.sqrt(len(y))
        if ab_max is None:
            r = lambda p: self.model(p, x0, s, loss) - y
        else:
            r = lambda p: np.r_[self.model(p, x0, s, loss) - y, w * max(p[0] + p[1] - ab_max, 0.0)]
        n = len(starts[0])
        ub = np.full(n, np.inf); ub[0] = a_max
        lb = np.full(n, -np.inf)
        if loss:
            lb[-1], ub[-1] = 0.0, 60.0
        starts = [np.clip(st, lb + 1e-9, np.minimum(ub, np.r_[a_max * 0.999, np.full(n - 1, np.inf)])) for st in starts]
        best = min((least_squares(r, st, bounds=(lb, ub)) for st in starts), key=lambda q: q.cost)
        res = best.fun[: len(y)]
        return best.x, float(np.sum(res ** 2)), res

    def apply_transfer(self, x0, ratio):
        """x0 with every frequency multiplied by ratio (array on self.fr, numpy sign: a delay makes arg negative).
        The spectrum of the last x0 is kept (fits apply many transfers to one x0)."""
        c = getattr(self, "_spec_cache", None)
        if c is None or c[0] is not x0:
            c = self._spec_cache = (x0, np.fft.rfft(x0, self.np_))
        return np.fft.irfft(c[1] * ratio, self.np_)[: self.n]


def n_eff(res):
    """Independent points in a correlated residual: N / (1 + 2 sum rho_k), rho summed to its first zero.
    Returns (N_eff, correlation length in samples)."""
    r = res - res.mean()
    ac = np.correlate(r, r, "full")[len(r) - 1:] / np.dot(r, r)
    k0 = int(np.argmax(ac <= 0)) if np.any(ac <= 0) else len(ac)
    return len(r) / (1 + 2 * ac[1:k0].sum()), k0


def f_test(misfit, misfit_free, ne, n_par=4, n_fixed=1):
    """F statistic and p of a nested fit with n_fixed parameters fewer against the free fit with n_par
    parameters (F(n_fixed, N_eff - n_par)). With N_eff <= n_par the test has no power: p = 1."""
    dof = ne - n_par
    if dof <= 0:
        return float("nan"), 1.0
    F = (misfit - misfit_free) / n_fixed / (misfit_free / dof)
    return F, float(fdist.sf(max(F, 0.0), n_fixed, dof))


def slice_factor(f_thz):
    """Every-extreme factor of Sec. 4.6 at f (THz): x6 x4 x4/3 and the slice factor 2.2 up to 1 THz, 3.1 from
    2 THz (Sim 45), linear between."""
    return 6 * 4 * 4 / 3 * np.interp(f_thz, [1.0, 2.0], [2.2, 3.1])


def node_transfer(state, fr_thz, lam=1.0, use="loss_db_dipole", delay_factor="path"):
    """exp(change of ln T) of the slice for a node state from results/sim30/sim36_signedband.json (values with
    every extreme at the solved frequencies): the factor is taken off at each solved frequency, the raw values
    interpolated linearly in f with 0 at f = 0 and held above the last, then multiplied by the factor at f and
    by lam (1/70.4: internode 1000 um, packing 0.5, 100 um slice). The loss keeps the larger slice factor
    (2.2 to 1 THz, 3.1 from 2 THz; loss-based, Sec. 4.6). The delay is scaled with the delay-based path factor
    2.2 at every frequency (delay_factor 'path', the 220 um optical thickness of Sec. 5.5) or, with 'loss',
    with the same factor as the loss."""
    f = np.array(state["f_thz"]); fac = slice_factor(f)
    raw_l = -np.array(state[use]) * np.log(10) / 20 / fac
    raw_d = np.array(state["delay_rad"]) / fac
    f = np.r_[0.0, f]; raw_l = np.r_[0.0, raw_l]; raw_d = np.r_[0.0, raw_d]
    fd = slice_factor(np.minimum(fr_thz, 1.0)) if delay_factor == "path" else slice_factor(fr_thz)
    d = lam * (slice_factor(fr_thz) * np.interp(fr_thz, f, raw_l) - 1j * fd * np.interp(fr_thz, f, raw_d))
    return np.exp(d)


def supply_test(tp, fit, y, xf, tr, grid, q0, pr, s_ref, ne, p_pass=0.05, n_size=1):
    """Can a candidate change supply a causal loss increment and a shift by itself? (Sims 53-54)

    Reference: fit(y, xf, loss="minphase") gave parameters pr and misfit s_ref (a, b, s, tau, kappa free;
    N_eff ne); fit(y, xf, s=0.0) gave q0. For each size v in grid the candidate tr(v) (a transfer on tp.fr)
    is applied to xf and fitted with kappa = s = 0 (candidate) and with both free (full model, six
    parameters with v). The candidate supplies the change if (kappa, s) = 0 is not rejected, F(2, N_eff - 6),
    with N_eff from the reference and from the candidate's own residual ('edge' if they disagree); sizes
    pass where (v, kappa, s) is not rejected, F(3, N_eff - 6). A candidate with n_size sizes (grid then lists
    combinations) is charged for them: the full model has 5 + n_size parameters."""
    S_c, S_f, R_c, prev_c, prev_f = [], [], [], list(q0), list(pr)
    for v in grid:
        xt = tp.apply_transfer(xf, tr(v))
        qc, mc, rc = min((fit(y, xt, s=0.0, p0=prev_c), fit(y, xt, s=0.0, p0=q0)), key=lambda t: t[1])
        qf, mf, _ = min((fit(y, xt, loss="minphase", p0=prev_f), fit(y, xt, loss="minphase", p0=pr)),
                        key=lambda t: t[1])
        prev_c, prev_f = list(qc), list(qf)
        S_c.append(mc); S_f.append(mf); R_c.append(rc)
    S_full = min(min(S_f), s_ref)
    ib = int(np.argmin(S_c))
    ne_own, _ = n_eff(R_c[ib])
    npar = 5 + n_size
    p1 = f_test(S_c[ib], S_full, ne, npar, 2)[1]
    p2 = f_test(S_c[ib], S_full, ne_own, npar, 2)[1]
    verdict = "pass" if min(p1, p2) >= p_pass else ("fail" if max(p1, p2) < p_pass else "edge")
    ok = [v for v, m in zip(grid, S_c) if f_test(m, S_full, ne, npar, 2 + n_size)[1] >= p_pass]
    rng = [float(v) for v in ok] if n_size == 1 else ok
    return dict(verdict=verdict, p=p1, p_own_neff=p2, best=grid[ib] if n_size > 1 else float(grid[ib]), misfit=float(S_c[ib]),
                misfit_full=float(S_full), pass_range=([min(rng), max(rng)] if rng else None) if n_size == 1 else rng)


def block_resample(rng, res, blk):
    """Circular block bootstrap of a residual: blocks of blk samples at random starts, cut to len(res)."""
    n = len(res)
    starts = rng.integers(0, n, size=n // max(blk, 1) + 1)
    return np.concatenate([np.take(res, np.arange(s0, s0 + blk), mode="wrap") for s0 in starts])[:n]


def null_calibrated_test(fit, y, x_null, x_alts, s_null, n_fixed, n_par=6, nb=100, rng=None, blk=None,
                         model=None, batch=None, level=0.05):
    """A nested F test calibrated on its own null by a block bootstrap.

    Null: fit(y, x_null, s=s_null) with no loss term (s_null 0.0 fixes the shift, None leaves it free).
    Alternative: the best of fit(y, x, loss='minphase') over x in x_alts (s free). The statistic is
    f_test(null misfit, alternative misfit, N_eff of the best alternative's residual, n_par, n_fixed), computed
    the same way on the data and on nb synthetic scans: the null model plus the best alternative's residual
    resampled in circular blocks of blk samples (default: its correlation length). model(params, x, s) gives
    the null model. Returns dict(p_nominal, p_cal, thr, rej_nominal, nb): p_cal = (1 + #{p_b <= p}) / (1 + nb).
    With batch set, draws come in batches of that size and stop early once the 99 % Clopper-Pearson interval
    of #{p_b <= p} / n lies wholly above or below level; nb is then the number of draws made.
    """
    rng = rng if rng is not None else np.random.default_rng(0)

    def stat(yy, p0n=None, p0a=None):
        qn, mn, _ = fit(yy, x_null, s=s_null, p0=p0n) if p0n is not None else fit(yy, x_null, s=s_null)
        alts = [fit(yy, x, loss="minphase", p0=p0a) if p0a is not None else fit(yy, x, loss="minphase") for x in x_alts]
        qa, ma, ra = min(alts, key=lambda t: t[1])
        ne, corr = n_eff(ra)
        return f_test(mn, ma, ne, n_par, n_fixed)[1], qn, qa, ra, corr, ne

    p_obs, qn, qa, ra, corr, ne_obs = stat(y)
    base = model(qn, x_null, s_null)
    blk = blk or max(corr, 2)
    from scipy.stats import beta
    pb, neb = [], []
    while len(pb) < nb:
        for _ in range(min(batch or nb, nb - len(pb))):
            st = stat(base + block_resample(rng, ra, blk), qn, qa)
            pb.append(st[0]); neb.append(st[5])
        k, n = int(np.sum(np.array(pb) <= p_obs)), len(pb)
        lo = beta.ppf(0.005, k, n - k + 1) if k > 0 else 0.0
        hi = beta.ppf(0.995, k + 1, n - k) if k < n else 1.0
        if batch and (hi < level or lo > level):
            break
    pb = np.array(pb)
    return dict(p_nominal=float(p_obs), p_cal=float((1 + np.sum(pb <= p_obs)) / (1 + len(pb))),
                thr=float(np.quantile(pb, 0.05)), rej_nominal=float(np.mean(pb < 0.05)), nb=len(pb),
                n_nan=int(np.sum(~np.isfinite(pb))), ne_obs=float(ne_obs), ne_draws=[float(np.min(neb)), float(np.median(neb))])


def multistart(one):
    """Wrap a fit(y, x, s=, loss=, p0=) so each call returns the best of several starts: the default, the
    given p0, and the nested model's solution (for a loss fit, the no-loss fit plus kappa = 1; for a no-loss
    fit, the loss fit without kappa, with s folded into tau when s is fixed)."""
    def fit(y, x, s=None, loss=None, p0=None):
        tries = [one(y, x, s=s, loss=loss)]
        if p0 is not None:
            tries.append(one(y, x, s=s, loss=loss, p0=p0))
        if loss:
            q, _, _ = one(y, x)
            tries.append(one(y, x, loss=loss, p0=list(q) + [1.0]))
        else:
            q, _, _ = one(y, x, loss="minphase")
            tries.append(one(y, x, s=s, p0=list(q[:4]) if s is None else [q[0], q[1], q[2] + q[3]]))
        return min(tries, key=lambda t: t[1])
    return fit


def phase_surrogate(rng, res):
    """A surrogate of res with the same amplitude spectrum (so the same circular autocorrelation and N_eff)
    and random Fourier phases (DC and Nyquist kept)."""
    R = np.fft.rfft(res)
    S = np.abs(R) * np.exp(1j * rng.uniform(0, 2 * np.pi, len(R)))
    S[0] = R[0]
    if len(res) % 2 == 0:
        S[-1] = R[-1]
    return np.fft.irfft(S, len(res))


def clopper_pearson(k, n, conf=0.95):
    from scipy.stats import beta
    a = (1 - conf) / 2
    lo = beta.ppf(a, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(1 - a, k + 1, n - k) if k < n else 1.0
    return float(lo), float(hi)


def size_fitted_test(fit_multi, fit_one, y, xs, s_null, n_fixed, model, n_par=6, rng=None, nb_max=400, batch=50,
                     level=0.05, conf=0.95):
    """Nested F test of 'a candidate of some size, no loss term' against 'the candidate of some size with the
    causal loss increment', the size fitted over the grid xs (xs[0] = no candidate) inside the statistic,
    calibrated on its own null by phase-randomised surrogates of the alternative's residual.

    Statistic Q(y, W): null misfit = min over g of fit_one(y, xs[g], s=s_null, p0=W0[g]); alternative misfit =
    min over g of the better of fit_one(y, xs[g], loss='minphase') from W1[g] and from the null solution at g
    plus kappa = 1; p = f_test(null, alternative, N_eff of the best alternative's residual on the data); the
    synthetic scans use the data's N_eff too, so p is a fixed monotone function of the misfit ratio and the
    surrogates alone supply its null distribution (their own N_eff is returned as a diagnostic). W0, W1 come from
    multi-start fits (fit_multi) on the data; Q is computed the same way on the data and on every synthetic
    scan (the null model at its best size plus a surrogate). Draws come in batches until the conf interval of
    the fraction of p_b <= p excludes level (verdict 'pass' if above, 'fail' if below) or nb_max is reached
    ('undecided'). Returns dict(p_nominal, p_cal, ci, nb, verdict, size_index_null, size_index_alt, ne_obs,
    ne_draws_median)."""
    rng = rng if rng is not None else np.random.default_rng(0)
    W0 = [fit_multi(y, x, s=s_null)[0] for x in xs]
    W1 = [fit_multi(y, x, loss="minphase")[0] for x in xs]

    def nested(q0):
        q = list(q0) if s_null is None else [q0[0], q0[1], s_null, q0[2]]
        return q + [1.0]

    def Q(yy, ne_use=None):
        nul = [fit_one(yy, x, s=s_null, p0=W0[g]) for g, x in enumerate(xs)]
        alt = [min((fit_one(yy, x, loss="minphase", p0=W1[g]), fit_one(yy, x, loss="minphase", p0=nested(nul[g][0]))),
                   key=lambda t: t[1]) for g, x in enumerate(xs)]
        g0 = int(np.argmin([t[1] for t in nul])); g1 = int(np.argmin([t[1] for t in alt]))
        ne, _ = n_eff(alt[g1][2])
        return f_test(nul[g0][1], alt[g1][1], ne_use or ne, n_par, n_fixed)[1], g0, g1, nul[g0][0], alt[g1][2], ne

    p_obs, g0, g1, q0, ra, ne_obs = Q(y)
    base = model(q0, xs[g0], s_null)
    pb, neb = [], []
    verdict = "undecided"
    while len(pb) < nb_max:
        for _ in range(min(batch, nb_max - len(pb))):
            r = Q(base + phase_surrogate(rng, ra), ne_obs)
            pb.append(r[0]); neb.append(r[5])
        k, n = int(np.sum(np.array(pb) <= p_obs)), len(pb)
        lo, hi = clopper_pearson(k, n, conf)
        if lo > level:
            verdict = "pass"; break
        if hi < level:
            verdict = "fail"; break
    k, n = int(np.sum(np.array(pb) <= p_obs)), len(pb)
    return dict(p_nominal=float(p_obs), p_cal=k / n, ci=clopper_pearson(k, n, conf), nb=n, verdict=verdict,
                size_index_null=g0, size_index_alt=g1, ne_obs=float(ne_obs), ne_draws_median=float(np.median(neb)))


def effect_size(S_r, S_f, y):
    """What a restricted fit leaves unexplained against a fuller fit, in % of the trace rms:
    sqrt(max(S_r - S_f, 0) / N) / rms(y). For linear models it is the rms difference of the two fitted curves."""
    return 100 * np.sqrt(max(S_r - S_f, 0.0) / len(y)) / np.sqrt(np.mean(y ** 2))


def effect_noise(fit_r, fit_f, y_model, res, nb, rng):
    """Effect sizes that fitting noise alone gives when the restricted model is true: y_model (the restricted
    model's best fit) plus a phase-randomised copy of res (same autocorrelation), refitted by both fits.
    fit_r and fit_f take a trace and return (params, misfit, ...). Returns an array of nb effect sizes."""
    out = []
    for _ in range(nb):
        yy = y_model + phase_surrogate(rng, res)
        out.append(effect_size(fit_r(yy)[1], fit_f(yy)[1], yy))
    return np.array(out)


def noise_scale(fit_f, y_model, res, n, rng):
    """Factor by which the full fit's residual understates the noise: inject phase-randomised copies of
    res into the full fit's curve y_model, refit, and return 1 / mean(sqrt(S_after / |noise|^2))."""
    r = [np.sqrt(fit_f(y_model + e)[1] / np.sum(e ** 2)) for e in (phase_surrogate(rng, res) for _ in range(n))]
    return float(1 / np.mean(r))


def threshold_verdict(D, draws, tracing, level, nboot=200, rng=None):
    """Any-of-scans rule. D: effect size per scan; draws: noise draws per scan (scans x nb); tracing: tracing D.
    Threshold per scan = hypot(quantile(draws, level), tracing). Returns dict with the ratio per scan, the
    largest, excluded (largest > 1), the false-exclusion rate counted on the draws, and the fraction of
    bootstrap resamples of the draws that flip the verdict."""
    draws = np.atleast_2d(np.asarray(draws, float)); D = np.asarray(D, float)
    th = np.hypot(np.quantile(draws, level, axis=1), tracing)
    ratio = D / th; exc = bool(ratio.max() > 1)
    false = float(np.mean(np.any(draws / th[:, None] > 1, axis=0)))
    rng = rng or np.random.default_rng(0); nb = draws.shape[1]; flips = 0
    for _ in range(nboot):
        idx = rng.integers(0, nb, nb)
        tb = np.hypot(np.quantile(draws[:, idx], level, axis=1), tracing)
        flips += bool((D / tb).max() > 1) != exc
    return dict(ratio=ratio.tolist(), largest=float(ratio.max()), excluded=exc, false_rate=false, flip=flips / nboot if nboot else float("nan"))
