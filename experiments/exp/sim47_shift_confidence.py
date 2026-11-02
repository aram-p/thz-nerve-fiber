"""Sim 47 - Review round 18: a confidence interval for the shift s of the tissue pulse.

Model (Sec. 5.6): biased(t) = a * B(t - s) + b * K(t - s - tau), B the unbiased pulse traced from the same
figure as the biased one, K the empty-cuvette pulse (paper 1 Fig. 3a). s is measured between curves on the
same axes, so axis calibration cancels; its error comes from the curve tracing and from what the model
leaves out, both of which are in the residual. The residual is correlated over a fraction of a picosecond,
so two methods are used:

(1) profile: for s on a grid the other parameters are refitted; with the noise variance estimated from the
    best fit and N_eff = N / (1 + 2 sum rho_k) independent points (rho summed to its first zero),
    F(s) = (S(s) - S_min) / (S_min / (N_eff - 4)), compared with F(1, N_eff - 4) (the 4 fitted parameters
    taken off the degrees of freedom; review round 19).
(2) moving-block bootstrap: residual blocks of the correlation length are resampled (circularly) and added to
    the best-fit model, which is refitted; percentiles of s.

Node states give at most 0.063 rad of extra delay at 0.6 THz with every extreme (Sec. 5.7), a shift of
0.0168 ps; that value and s = 0 are tested per scan. No joint test: the scans share the empty-cuvette
pulse and come from one session, so they are not independent. Sim 48 applies each node state's whole-band
change inside the fit, which supersedes the 0.6 THz shift test.
Writes results/sim47/sim47.txt and results/sim47/sim47.json.
"""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import f as fdist, norm

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim47"; OUT.mkdir(parents=True, exist_ok=True)
lines = []


def say(s=""):
    print(s, flush=True); lines.append(s)


H1 = np.load(ROOT / "results" / "sim42" / "digitised.npz")
H3 = np.load(ROOT / "results" / "sim39" / "digitised.npz")
P2 = np.load(ROOT / "results" / "sim43" / "paper2_2min.npz")
tt = H1["t_ps"]; dt = tt[1] - tt[0]
h3 = {k: np.interp(tt, H3["t_ps"], H3["wave_" + k]) for k in ("blue", "red")}
K = H1["empty"]
S_NODE = 0.0634 / (2 * np.pi * 0.6e12) * 1e12        # ps: largest node-state delay at 0.6 THz, every extreme
SCANS = {"30 s (H1)": (H1["30s"], H1["none"]), "2 min (H2)": (P2["biased"], P2["none"]),
         "4 min (H1)": (H1["4min"], H1["none"]), "5 min (H3)": (h3["red"], h3["blue"])}
rng = np.random.default_rng(47)


def shift(y, tau):
    return np.interp(tt - tau, tt, y, left=0.0, right=0.0)


def model(p, x0, s=None):
    a, b, s_, tau = p if s is None else (p[0], p[1], s, p[2])
    return a * shift(x0, s_) + b * shift(K, s_ + tau)


def fit(y, x0, s=None, p0=None):
    if s is None:
        starts = [p0] if p0 is not None else [[0.6, b0, s0, 0.0] for b0 in (0.05, 0.2) for s0 in (-0.1, 0.0, 0.1)]
        r = lambda p: model(p, x0) - y
    else:
        starts = [p0] if p0 is not None else [[0.6, b0, t0] for b0 in (0.05, 0.2) for t0 in (-0.1, 0.0, 0.1)]
        r = lambda p: model(p, x0, s) - y
    best = min((least_squares(r, st) for st in starts), key=lambda q: q.cost)
    return best.x, float(np.sum(best.fun ** 2)), best.fun


def n_eff(res):
    r = res - res.mean()
    ac = np.correlate(r, r, "full")[len(r) - 1:] / np.dot(r, r)
    k0 = int(np.argmax(ac <= 0)) if np.any(ac <= 0) else len(ac)
    return len(r) / (1 + 2 * ac[1:k0].sum()), k0 * dt


say("Sim 47: confidence interval for the shift s of the tissue pulse (tau free)")
say(f"  node-state limit: 0.063 rad at 0.6 THz = s {S_NODE:.4f} ps")
out = {}
for name, (y, x0) in SCANS.items():
    p, smin, res = fit(y, x0)
    ne, lc = n_eff(res)
    grid = np.arange(p[2] - 0.08, p[2] + 0.0801, 0.002)
    prof = []
    q0 = [p[0], p[1], p[3]]
    for s in grid:
        q, ss, _ = fit(y, x0, s=s, p0=q0)
        prof.append(ss)
    prof = np.array(prof)
    dof = ne - 4
    Fs = lambda ss: (ss - smin) / (smin / dof)
    d = Fs(prof)
    ci = lambda lev: (grid[d <= fdist.ppf(lev, 1, dof)].min(), grid[d <= fdist.ppf(lev, 1, dof)].max())
    c68, c95 = ci(0.6827), ci(0.95)
    lo1 = grid[d <= fdist.ppf(0.90, 1, dof)].min()     # one-sided 95 % lower bound from the profile
    dnode = Fs(fit(y, x0, s=S_NODE, p0=q0)[1])
    dzero = Fs(fit(y, x0, s=0.0, p0=q0)[1])
    sig = lambda F: norm.isf(fdist.sf(max(F, 0), 1, dof) / 2)
    # moving-block bootstrap
    L = max(2, int(round(lc / dt)))
    base = model(p, x0)
    boot = []
    for _ in range(400):
        starts = rng.integers(0, len(res), size=int(np.ceil(len(res) / L)))
        rr = np.concatenate([np.take(res, np.arange(st, st + L), mode="wrap") for st in starts])[: len(res)]
        boot.append(fit(base + rr, x0, p0=p)[0][2])
    boot = np.array(boot)
    one_sided = lambda arr, lev: float(np.percentile(arr, 100 * (1 - lev)))
    say(f"  {name}: s = {p[2]:+.4f} ps, tau = {p[3]:+.4f} ps, residual {np.sqrt(smin / np.sum(y ** 2)):.0%}; "
        f"correlation length {lc:.2f} ps, N = {len(y)}, N_eff = {ne:.1f}")
    say(f"      profile: 68 % {c68[0]:+.3f} to {c68[1]:+.3f} ps, 95 % {c95[0]:+.3f} to {c95[1]:+.3f} ps; "
        f"one-sided 95 % lower bound {lo1:+.4f} ps; F (1, {dof:.1f} dof) at s = {S_NODE:.4f} (node limit) {dnode:.1f} "
        f"({sig(dnode):.1f} sigma two-sided), at s = 0 {dzero:.1f} ({sig(dzero):.1f} sigma)")
    say(f"      block bootstrap (block {L * dt:.2f} ps, 400 resamples): s {np.mean(boot):+.4f} +- {np.std(boot):.4f} ps; "
        f"one-sided 95 % lower bound {one_sided(boot, 0.95):+.4f} ps, 97.5 % {one_sided(boot, 0.975):+.4f} ps; "
        f"fraction of resamples with s <= node limit {np.mean(boot <= S_NODE):.3f}")
    out[name] = dict(s=float(p[2]), tau=float(p[3]), n_eff=float(ne), corr_len_ps=float(lc), ci68=list(map(float, c68)),
                     ci95=list(map(float, c95)), prof_lo95=float(lo1), F_node=float(dnode), F_zero=float(dzero),
                     sigma_node=float(sig(dnode)), sigma_zero=float(sig(dzero)), dof=float(dof),
                     boot_mean=float(np.mean(boot)), boot_sd=float(np.std(boot)),
                     boot_lo95=one_sided(boot, 0.95), boot_lo975=one_sided(boot, 0.975))
(OUT / "sim47.json").write_text(json.dumps(out, indent=1))
(OUT / "sim47.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
