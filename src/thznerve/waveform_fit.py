"""Two-pulse fits of the traced published waveforms (Sims 47-49, Sec. 5.6).

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
    def __init__(self, tt, K, pad=4):
        self.tt, self.K = tt, K
        self.dt = tt[1] - tt[0]
        self.n = len(tt)
        self.np_ = pad * self.n
        self.fr = np.fft.rfftfreq(self.np_, self.dt)          # THz when t is in ps

    def shift(self, y, tau):
        return np.interp(self.tt - tau, self.tt, y, left=0.0, right=0.0)

    def model(self, p, x0, s=None):
        a, b, s_, tau = p if s is None else (p[0], p[1], s, p[2])
        return a * self.shift(x0, s_) + b * self.shift(self.K, s_ + tau)

    def fit(self, y, x0, s=None, p0=None, a_max=np.inf):
        """Least squares over (a, b, s, tau), or (a, b, tau) with s fixed; a <= a_max (a passive slice cannot
        raise the tissue pulse above what its own transfer gives: a_max = 1). Returns (params, misfit, residual)."""
        if p0 is not None:
            starts = [p0]
        elif s is None:
            starts = [[0.6, b0, s0, 0.0] for b0 in (0.05, 0.2) for s0 in (-0.1, 0.0, 0.1)]
        else:
            starts = [[0.6, b0, t0] for b0 in (0.05, 0.2) for t0 in (-0.1, 0.0, 0.1)]
        r = lambda p: self.model(p, x0, s) - y
        n = len(starts[0])
        ub = np.full(n, np.inf); ub[0] = a_max
        starts = [np.minimum(st, np.r_[a_max * 0.999, np.full(n - 1, np.inf)]) for st in starts]
        best = min((least_squares(r, st, bounds=(-np.inf, ub)) for st in starts), key=lambda q: q.cost)
        return best.x, float(np.sum(best.fun ** 2)), best.fun

    def apply_transfer(self, x0, ratio):
        """x0 with every frequency multiplied by ratio (array on self.fr, numpy sign: a delay makes arg negative)."""
        return np.fft.irfft(np.fft.rfft(x0, self.np_) * ratio, self.np_)[: self.n]


def n_eff(res):
    """Independent points in a correlated residual: N / (1 + 2 sum rho_k), rho summed to its first zero.
    Returns (N_eff, correlation length in samples)."""
    r = res - res.mean()
    ac = np.correlate(r, r, "full")[len(r) - 1:] / np.dot(r, r)
    k0 = int(np.argmax(ac <= 0)) if np.any(ac <= 0) else len(ac)
    return len(r) / (1 + 2 * ac[1:k0].sum()), k0


def f_test(misfit, misfit_free, ne, n_par=4):
    """F statistic and p of a fit with one parameter fewer against the free fit (F(1, N_eff - n_par))."""
    dof = ne - n_par
    F = (misfit - misfit_free) / (misfit_free / dof)
    return F, float(fdist.sf(max(F, 0.0), 1, dof))


def slice_factor(f_thz):
    """Every-extreme factor of Sec. 4.6 at f (THz): x6 x4 x4/3 and the slice factor 2.2 up to 1 THz, 3.1 from
    2 THz (Sim 45), linear between."""
    return 6 * 4 * 4 / 3 * np.interp(f_thz, [1.0, 2.0], [2.2, 3.1])


def node_transfer(state, fr_thz, lam=1.0, use="loss_db_dipole"):
    """exp(change of ln T) of the slice for a node state from results/sim30/sim36_signedband.json (values with
    every extreme at the solved frequencies): the factor is taken off at each solved frequency, the raw values
    interpolated linearly in f with 0 at f = 0 and held above the last, then multiplied by the factor at f and
    by lam (1/70.4: internode 1000 um, packing 0.5, 100 um slice)."""
    f = np.array(state["f_thz"]); fac = slice_factor(f)
    raw = (-np.array(state[use]) * np.log(10) / 20 - 1j * np.array(state["delay_rad"])) / fac
    f = np.r_[0.0, f]; raw = np.r_[0.0, raw]
    d = lam * slice_factor(fr_thz) * (np.interp(fr_thz, f, raw.real) + 1j * np.interp(fr_thz, f, raw.imag))
    return np.exp(d)
