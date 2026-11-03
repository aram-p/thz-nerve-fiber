"""Sim 48 - Review round 19: the waveform result as a phase lag at 0.6 THz, and node states inside the fit.

(A) Phase lag of the tissue pulse at 0.6 THz. Two model forms:
      shift:          biased(t) = a B(t - s) + b K(t - s - tau)
      shift + phase:  biased(t) = a B(t - s) + c H[B](t - s) + b K(t - s - tau)
    H the Hilbert transform; a B + c H[B] is B with every frequency delayed by the same phase atan2(c, a), so
    the lag at 0.6 THz is 2 pi 0.6 THz s + atan2(c, a). Three fit windows: the whole trace, t0 - 0.6 to t0 + 0.8 ps
    (t0 the peak of B), and the whole trace without t0 +- 0.4 ps. Spread from a moving-block bootstrap of the
    residual (block = its correlation length).
(B) Node states inside the fit. Each node state's computed change of ln T (loss and delay at every solved
    frequency, every extreme, Sec. 5.7; linear in f between solved points, 0 at f = 0, held above the last) is
    applied to B: B_node = irfft(rfft(B) exp(dlnT)); then biased = a B_node + b K(t - tau) with a, b, tau free.
    Compared with s = 0 (no node change) and s free: the fraction of the misfit reduction that s buys and the
    node state recovers, and an F test of the node fit against s free with N_eff - 4 degrees of freedom.
Writes results/sim48/sim48.txt and results/sim48/sim48.json.
"""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import hilbert
from scipy.stats import f as fdist, t as tdist

from thznerve.waveform_fit import node_transfer

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim48"; OUT.mkdir(parents=True, exist_ok=True)
lines = []


def say(s=""):
    print(s, flush=True); lines.append(s)


H1 = np.load(ROOT / "results" / "sim42" / "digitised.npz")
H3 = np.load(ROOT / "results" / "sim39" / "digitised.npz")
P2 = np.load(ROOT / "results" / "sim43" / "paper2_2min.npz")
tt = H1["t_ps"]; dt = tt[1] - tt[0]; N = len(tt)
h3 = {k: np.interp(tt, H3["t_ps"], H3["wave_" + k]) for k in ("blue", "red")}
K = H1["empty"]
F0 = 0.6                                                  # THz
SCANS = {"30 s (H1)": (H1["30s"], H1["none"]), "2 min (H2)": (P2["biased"], P2["none"]),
         "4 min (H1)": (H1["4min"], H1["none"]), "5 min (H3)": (h3["red"], h3["blue"])}
NODES = json.loads((ROOT / "results" / "sim30" / "sim36_signedband.json").read_text())["states"]
rng = np.random.default_rng(48)
NB = 200


def shift(y, tau):
    return np.interp(tt - tau, tt, y, left=0.0, right=0.0)


def model(p, x0, hx0, form, s=None):
    if form == "shift":
        a, b, s_, tau = p if s is None else (p[0], p[1], s, p[2]); c = 0.0
    else:
        a, b, s_, tau, c = p if s is None else (p[0], p[1], s, p[2], p[3])
    return a * shift(x0, s_) + c * shift(hx0, s_) + b * shift(K, s_ + tau)


def fit(y, x0, hx0, form, m, s=None, p0=None):
    n = 4 if form == "shift" else 5
    if s is not None:
        n -= 1
    if p0 is not None:
        starts = [p0]
    else:
        base = [[0.6, b0, s0, 0.0] for b0 in (0.05, 0.2) for s0 in (-0.1, 0.0, 0.1)]
        if s is not None:
            base = [[q[0], q[1], q[3] + t0] for q in base[:2] for t0 in (-0.1, 0.0, 0.1)]
        starts = [q + [0.0] for q in base] if form != "shift" else base
    r = lambda p: (model(p, x0, hx0, form, s) - y)[m]
    best = min((least_squares(r, st) for st in starts), key=lambda q: q.cost)
    return best.x, float(np.sum(best.fun ** 2)), best.fun


def lag(p, form):
    return 2 * np.pi * F0 * p[2] + (0.0 if form == "shift" else np.arctan2(p[4], p[0]))


def n_eff(res):
    r = res - res.mean()
    ac = np.correlate(r, r, "full")[len(r) - 1:] / np.dot(r, r)
    k0 = int(np.argmax(ac <= 0)) if np.any(ac <= 0) else len(ac)
    return len(r) / (1 + 2 * ac[1:k0].sum()), k0


def windows(x0):
    t0 = tt[np.argmax(np.abs(x0))]
    return {"whole trace": np.ones(N, bool), "t0-0.6 to t0+0.8 ps": (tt >= t0 - 0.6) & (tt <= t0 + 0.8),
            "without t0+-0.4 ps": np.abs(tt - t0) > 0.4}


NP = 4 * N
fr = np.fft.rfftfreq(NP, dt)                              # THz (dt in ps)


def apply_node(x0, st, use="loss_db_dipole", lam=1.0):
    return np.fft.irfft(np.fft.rfft(x0, NP) * node_transfer(st, fr, lam, use, delay_factor="loss"), NP)[:N]


out = {}
say("Sim 48 (A): phase lag of the tissue pulse at 0.6 THz, per scan, two model forms and three windows "
    f"(moving-block bootstrap, {NB} resamples)")
for name, (y, x0) in SCANS.items():
    hx0 = np.imag(hilbert(x0))
    out[name] = {"lag": {}}
    for form in ("shift", "shift + phase"):
        for wn, m in windows(x0).items():
            p, smin, res = fit(y, x0, hx0, form, m)
            ne, k0 = n_eff(res)
            L = max(2, k0)
            full_res = np.zeros(N); full_res[m] = res
            base = model(p, x0, hx0, form)
            idx = np.flatnonzero(m)
            boot = []
            for _ in range(NB):
                st = rng.integers(0, len(res), size=int(np.ceil(len(res) / L)))
                rr = np.concatenate([np.take(res, np.arange(q, q + L), mode="wrap") for q in st])[: len(res)]
                yb = base.copy(); yb[idx] = base[idx] + rr
                boot.append(lag(fit(yb, x0, hx0, form, m, p0=p)[0], form))
            boot = np.array(boot)
            sd = float(np.std(boot)); g = float(lag(p, form))
            dof = max(ne - len(p), 1.0)
            lo95 = g - tdist.ppf(0.95, dof) * sd
            out[name]["lag"][f"{form} | {wn}"] = dict(lag_rad=g, sd_rad=sd, lo95_rad=lo95, s_ps=float(p[2]), n_eff=float(ne),
                                                      residual=float(np.sqrt(smin / np.sum(y[m] ** 2))))
            say(f"  {name:11s} {form:13s} {wn:20s}: lag {g:+.3f} +- {sd:.3f} rad (one-sided 95 % lower bound "
                f"{lo95:+.3f}, t with N_eff - {len(p)} = {dof:.1f} dof), s {p[2]:+.4f} ps, residual "
                f"{np.sqrt(smin / np.sum(y[m] ** 2)):.0%}")
    lags = [v["lag_rad"] for v in out[name]["lag"].values()]
    los = [v["lo95_rad"] for v in out[name]["lag"].values()]
    out[name]["lag_range"] = [min(lags), max(lags)]; out[name]["lo95_min"] = min(los)
    say(f"  {name}: lag {min(lags):.3f} to {max(lags):.3f} rad over forms and windows; smallest 95 % lower bound {min(los):+.3f}")

nmax = max(NODES, key=lambda s: s["delay_rad"][s["f_thz"].index(F0)])
say(f"\n  largest node-state delay at 0.6 THz, every extreme: {nmax['state']} {nmax['delay_rad'][nmax['f_thz'].index(F0)]:.4f} rad")

say("\nSim 48 (B): each node state's change of ln T over the band applied to the unbiased pulse inside the fit "
    f"({len(NODES)} states, scaled to a fraction lam of every extreme, 1/70.4 being nominal geometry; loss of the dipole "
    "observable or the largest of three). Per family: the best state, the share of the s-free misfit reduction it "
    "recovers, and the F test p of its fit against s free")
FAM = ("sigma", "epsnode", "epsfine", "epsdense", "joint", "aqueous")
LAMS = (1 / 70.4, 1 / 5.33, 0.25, 0.5, 0.75, 1.0)    # 1/70.4: internode 1000 um, packing 0.5, 100 um; 1/5.33: 200 um, 0.6, 220 um


def aqueous(state):
    """States an aqueous gap can hold at THz: extra conductivity up to saturated NaCl (25 S/m), eps' >= 1, or
    negative eps' with eps'' >= 3 |eps'| (ion collisions make eps'' >> |eps'|, Sec. 4.6)."""
    fam, k = state.split(":")
    if fam == "sigma":
        return float(k) <= 25.0
    if fam == "joint":
        return False
    e1, e2 = (float(x) for x in k.split("_")[:2])
    return e1 >= 1.0 or e2 >= 3 * abs(e1)
for name, (y, x0) in SCANS.items():
    hx0 = np.imag(hilbert(x0)); m = np.ones(N, bool)
    p, s_free, res = fit(y, x0, hx0, "shift", m)
    ne, _ = n_eff(res)
    _, s_zero, _ = fit(y, x0, hx0, "shift", m, s=0.0)
    gain = s_zero - s_free
    Fp = lambda sn: ((sn - s_free) / (s_free / (ne - 4)), fdist.sf(max((sn - s_free) / (s_free / (ne - 4)), 0), 1, ne - 4))
    fz, pz = Fp(s_zero)
    res_ = {"n_eff": float(ne), "gain_frac": gain / s_zero, "F_zero": fz, "p_zero": pz, "by_lam": {}}
    say(f"  {name}: s free removes {gain / s_zero:.1%} of the s = 0 misfit (F = {fz:.1f}, p = {pz:.4f}, 1 and {ne - 4:.1f} dof)")
    if name.startswith("30 s"):
        out[name]["node_fit"] = res_; continue
    for lam in LAMS:
        fam_best = {}
        for use in ("loss_db_dipole", "loss_db"):
            for st in NODES:
                _, sn, _ = fit(y, apply_node(x0, st, use, lam), hx0, "shift", m, s=0.0)
                for fm in (st["state"].split(":")[0],) + (("aqueous",) if aqueous(st["state"]) else ()):
                    if fm not in fam_best or sn < fam_best[fm][1]:
                        fam_best[fm] = (st["state"], sn, use)
        row = {}
        for fm in FAM:
            st_, sn, use = fam_best[fm]
            F, pv = Fp(sn)
            row[fm] = dict(state=st_, loss_kind=use, recovered=(s_zero - sn) / gain, F=F, p=pv)
        res_["by_lam"][f"{lam:.4g}"] = row
        say(f"    x{lam:.3g} of every extreme: " + "; ".join(
            f"{fm} {row[fm]['state'].split(':', 1)[1] if fm != 'aqueous' else row[fm]['state']} {row[fm]['recovered']:.0%} (p {row[fm]['p']:.2f})" for fm in FAM))
    out[name]["node_fit"] = res_

(OUT / "sim48.json").write_text(json.dumps(out, indent=1))
(OUT / "sim48.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
