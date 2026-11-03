"""Sim 57 - Review round 24: the bootstrap of Sim 55b was built on the full fit (shift and loss present), so its
'no change' scans were not null and its thresholds were too strict. Here each test is calibrated on its own null.

(A) For every scan (reference: 30 s floor 2.07 dB/THz zero phase, causal increment held above 2 THz, a <= 1):
synthetic scans are the NULL fit of that scan plus the residual of its full fit resampled in circular blocks of
one and two correlation lengths. Shift test: null = s = 0 with the loss term; no-change test: null = s = 0, no
loss term. The rejection rate at nominal 0.05 and the 5th percentile of p (the calibrated threshold) are
reported per scan, and the observed p is converted to a calibrated p (fraction of null p at or below it).
(B) Candidate verdicts of Sims 53-54 (six references each) under the per-scan calibrated 2-dof threshold.
(C) Effect size against an empirical null: one traced unbiased pulse fitted with another (H1, H2, H3 are three
tracings of no-voltage pulses from three figures), same model, no floor; the loss slope, the lag at 0.6 THz
and the rms of the change term (full model minus the no-change model, % of the trace rms) are compared with
those of the biased scans.
Writes results/sim57/sim57.txt and results/sim57/sim57.json.
"""
import json
import sys
import time
from itertools import permutations
from pathlib import Path

import numpy as np

from thznerve.runlog import RunLog
from thznerve.waveform_fit import TwoPulse, f_test, load_scans, n_eff

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim57"; OUT.mkdir(parents=True, exist_ok=True)
NB = int(sys.argv[1]) if len(sys.argv) > 1 else 200
F0 = 0.6
lines = []


def say(s=""):
    print(s, flush=True); lines.append(s)


def blocks(rng, res, blk):
    n = len(res)
    starts = rng.integers(0, n, size=n // blk + 1)
    return np.concatenate([np.take(res, np.arange(s0, s0 + blk), mode="wrap") for s0 in starts])[:n]


tt, K, SCANS = load_scans(ROOT)
LATE = ("2 min (H2)", "4 min (H1)", "5 min (H3)")
tp = TwoPulse(tt, K, f_loss=2.0)
i06 = int(np.argmin(np.abs(tp.fr - F0)))
FLOOR = tp.loss_transfer(2.07, "zero")
fit = lambda y, x, **kw: tp.fit(y, x, a_max=1.0, loss=kw.pop("loss", None), **kw)
out = {"n_boot": NB, "A": {}, "B": {}, "C": {}}

with RunLog("sim57", "round 24: per-scan bootstrap on the null fit of each test; candidates re-scored; tracing-to-tracing null") as log:
    t0 = time.time()
    for scan, (y, x0) in SCANS.items():
        xf = tp.apply_transfer(x0, FLOOR)
        pr, m, res = fit(y, xf, loss="minphase")
        ne, corr = n_eff(res)
        q1, m1, _ = fit(y, xf, s=0.0, loss="minphase")
        q0, m0, _ = fit(y, xf, s=0.0)
        p1_obs, p2_obs = f_test(m1, m, ne, 5, 1)[1], f_test(m0, m, ne, 5, 2)[1]
        base1 = tp.model(q1, xf, s=0.0, loss="minphase")
        base0 = tp.model(q0, xf, s=0.0)
        rng = np.random.default_rng(57)
        d = dict(p_shift=p1_obs, p_none=p2_obs, n_eff=float(ne), corr=int(corr))
        p1s, p2s = [], []
        for mult in (1, 2):
            blk = max(mult * corr, 2)
            p1, p2 = [], []
            for _ in range(NB):
                yb = base1 + blocks(rng, res, blk)
                _, mb, rb = fit(yb, xf, loss="minphase", p0=pr)
                _, mb1, _ = fit(yb, xf, s=0.0, loss="minphase")
                p1.append(f_test(mb1, mb, n_eff(rb)[0], 5, 1)[1])
                yb = base0 + blocks(rng, res, blk)
                _, mb, rb = fit(yb, xf, loss="minphase", p0=pr)
                _, mb0, _ = fit(yb, xf, s=0.0)
                p2.append(f_test(mb0, mb, n_eff(rb)[0], 5, 2)[1])
            p1, p2 = np.array(p1), np.array(p2)
            d[f"block x{mult}"] = dict(block=int(blk), rej1=float(np.mean(p1 < 0.05)), rej2=float(np.mean(p2 < 0.05)),
                                       thr1=float(np.quantile(p1, 0.05)), thr2=float(np.quantile(p2, 0.05)))
            p1s += p1.tolist(); p2s += p2.tolist()
        p1s, p2s = np.array(p1s), np.array(p2s)
        d["thr1"] = min(d["block x1"]["thr1"], d["block x2"]["thr1"])
        d["thr2"] = min(d["block x1"]["thr2"], d["block x2"]["thr2"])
        d["p_shift_cal"] = float(np.mean(p1s <= p1_obs)); d["p_none_cal"] = float(np.mean(p2s <= p2_obs))
        out["A"][scan] = d
        say(f"{scan}: corr {corr} samples; nominal 5 % rejects shift {d['block x1']['rej1']:.2f}/{d['block x2']['rej1']:.2f}, "
            f"no change {d['block x1']['rej2']:.2f}/{d['block x2']['rej2']:.2f} (blocks x1/x2); thresholds shift p < {d['thr1']:.3f}, "
            f"2-dof p < {d['thr2']:.3f}; observed p shift {p1_obs:.2g} (calibrated {d['p_shift_cal']:.3f}), "
            f"no change {p2_obs:.2g} (calibrated {d['p_none_cal']:.3f}, floor {1 / len(p2s):.4f})  [{time.time() - t0:.0f} s]")

    say()
    summ = {}
    for sim, f in (("53", ROOT / "results" / "sim53" / "sim53.json"), ("54", ROOT / "results" / "sim54" / "sim54.json")):
        dd = json.loads(f.read_text())
        vs = [v for v in dd if v.startswith("increment")]
        names = sorted({c for v in vs for c in dd[v]["2 min (H2)"]["cands"] if not c.startswith("node ")})
        has_nodes = any(c.startswith("node ") for c in dd[vs[0]]["2 min (H2)"]["cands"])
        for c in names + (["node states (best of the three screened)"] if has_nodes else []):
            row, rown = [], []
            for sc in LATE:
                if c.startswith("node states"):
                    ps = [max(x["p"] for n, x in dd[v][sc]["cands"].items() if n.startswith("node ")) for v in vs]
                else:
                    ps = [dd[v][sc]["cands"][c]["p"] for v in vs]
                row.append(int(sum(q >= out["A"][sc]["thr2"] for q in ps))); rown.append(int(sum(q >= 0.05 for q in ps)))
            summ.setdefault(sim, {})[c] = dict(calibrated=row, nominal=rown)
            say(f"  Sim {sim} {c:42s} references passing at 2/4/5 min (of 6): calibrated {row}, nominal 0.05 {rown}")
    out["B"] = summ

    say()
    tp0 = TwoPulse(tt, K, f_loss=2.0)
    ref = {"H1": SCANS["30 s (H1)"][1], "H2": SCANS["2 min (H2)"][1], "H3": SCANS["5 min (H3)"][1]}
    pairs = [(f"{a} from {b}", ref[a], ref[b]) for a, b in permutations(ref, 2)]
    pairs += [(sc, SCANS[sc][0], tp.apply_transfer(SCANS[sc][1], FLOOR)) for sc in SCANS]
    for name, y, x in pairs:
        pr, m, res = fit(y, x, loss="minphase")
        q0, m0, _ = fit(y, x, s=0.0)
        inc = tp.loss_transfer(pr[4], "minphase")
        lag = float(2 * np.pi * F0 * pr[2] - np.angle(inc[i06]))
        dchg = tp.model(pr, x, loss="minphase") - tp.model(q0, x, s=0.0)
        rms = lambda v: float(np.sqrt(np.mean(v ** 2)))
        out["C"][name] = dict(kappa=float(pr[4]), lag=lag, s=float(pr[2]), change_rms=100 * rms(dchg) / rms(y),
                              resid_rms=100 * rms(res) / rms(y))
        say(f"{name:22s} kappa {pr[4]:5.2f} dB/THz, lag {lag:+.3f} rad, change term {out['C'][name]['change_rms']:5.1f} %, "
            f"residual {out['C'][name]['resid_rms']:4.1f} % of the trace rms")
    log.outcome = "; ".join(f"{sc[:5]} thr2 {out['A'][sc]['thr2']:.3f} p_none_cal {out['A'][sc]['p_none_cal']:.3f}" for sc in SCANS)
(OUT / "sim57.json").write_text(json.dumps(out, indent=1, default=float))
(OUT / "sim57.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
