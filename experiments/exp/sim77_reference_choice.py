"""Sim 77 - Review round 36: which empty-cuvette pulse the H1-H3 data prefer, and H4 sample 2 without a second pulse.

(1) H1-H3 scans fitted with H1's and with H4's empty-cuvette pulse as K, no floor, a <= 1: residual (% of the trace
    rms) of the full fit (causal loss term, s free) and of the no-change fit with s = 0, and the second pulse's
    arrival in the full fit (s + tau, ps; negative = early).
(2) H4 sample 2 refitted (Sim 73b's fits) with no second pulse (K = 0), base tracing and the 8 tracing shifts, and
    the no-change fit with a free.
Writes results/sim77/sim77.txt.
"""
import itertools
from pathlib import Path

import numpy as np

from thznerve.runlog import RunLog
from thznerve.waveform_fit import TwoPulse, load_scans, multistart

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim77"


def load(name, file):
    import importlib.util
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / "exp" / file)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


S76 = load("s76", "sim76_bounded_second_pulse.py")
S73, A = S76.S73, S76.A


def no_k(W):
    return {**W, "empty": np.zeros_like(W["empty"])}


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    lines = []
    say = lambda s="": (print(s, flush=True), lines.append(s))
    with RunLog("sim77", "round 36: H1-H3 residuals with H1's vs H4's empty pulse; H4 sample 2 with no second pulse") as log:
        tt, K1, SC = load_scans(ROOT)
        h4 = np.load(ROOT / "results" / "sim73" / "digitised.npz")
        say("(1) H1-H3 with each empty-cuvette pulse as K (no floor, a <= 1): residual % of the full fit / of no change with s = 0; second pulse in the full fit (ps)")
        res = {}
        for KK, nm in ((K1, "H1"), (h4["empty"], "H4")):
            tp = TwoPulse(tt, KK, f_loss=2.0)
            fit = multistart(lambda y, x, **kw: tp.fit(y, x, a_max=1.0, **kw))
            for sc in ("30 s (H1)", "2 min (H2)", "4 min (H1)", "5 min (H3)"):
                y, x0 = SC[sc]
                qf, Sf, _ = fit(y, x0, loss="minphase"); _, S0, _ = fit(y, x0, s=0.0)
                n = np.sum(y ** 2); res[(nm, sc)] = (np.sqrt(Sf / n), np.sqrt(S0 / n))
                say(f"  K = {nm} empty pulse, {sc:10s}: full {100 * np.sqrt(Sf / n):5.1f} %, no change {100 * np.sqrt(S0 / n):5.1f} %, second pulse at {qf[2] + qf[3]:+.3f}")
        say()
        say("(2) H4 sample 2 with no second pulse (K = 0): residual % (range over the 8 tracing shifts); best size; parameters")
        W, valid = A.traces()
        base = S73.run(no_k(W), valid)["s2"]; fa = S76.free_a(no_k(W), valid, None)
        sp, fsp = [], []
        for sh in itertools.product((-1, 1), repeat=3):
            d = {"empty": sh[0], "s1_none": sh[1], "s2_none": sh[1], "s1_6min": sh[2], "s2_6min": sh[2]}
            Wd, vd = A.traces(d)
            sp.append(S73.run(no_k(Wd), vd)["s2"]); fsp.append(S76.free_a(no_k(Wd), vd, None))
        for k in ("uniform | alone", "uniform | drift", "uneven | alone", "uneven | drift", "none | alone", "none | drift"):
            v = base[k]
            say(f"  {k:16s} residual {100 * v['resid']:5.1f} % ({100 * min(x[k]['resid'] for x in sp):.1f}-{100 * max(x[k]['resid'] for x in sp):.1f})  size {v['size']}  a {v['q'][0]:+.3f}" +
                (f", s {v['q'][2]:+.3f}" if len(v['q']) == 4 else ""))
        for mode in ("alone", "drift"):
            v = fa[mode]
            say(f"  no change, a free, {mode:5s}: residual {100 * v['resid']:5.1f} % ({100 * min(x[mode]['resid'] for x in fsp):.1f}-{100 * max(x[mode]['resid'] for x in fsp):.1f})  a {v['q'][0]:+.3f}" +
                (f", s {v['q'][2]:+.3f}" if mode == "drift" else ""))
        log.outcome = (f"full fit H1 K {min(100 * res[('H1', s)][0] for s in SC):.1f}-{max(100 * res[('H1', s)][0] for s in SC):.1f} %, H4 K "
                       f"{min(100 * res[('H4', s)][0] for s in SC):.1f}-{max(100 * res[('H4', s)][0] for s in SC):.1f} %; s2 K=0 thinner drift "
                       f"{100 * base['uniform | drift']['resid']:.1f} %, none drift {100 * base['none | drift']['resid']:.1f} %, a free drift {100 * fa['drift']['resid']:.1f} %")
    (OUT / "sim77.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
