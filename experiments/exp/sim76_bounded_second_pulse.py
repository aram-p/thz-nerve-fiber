"""Sim 76 - Review round 35: two checks of the second pulse.

(1) H4 refitted (Sim 73b's fits) with the second pulse held within 0.2 ps of the empty-cuvette delay
    (|s + tau| <= 0.2 ps; TwoPulse.fit t2_max), the base tracing and the 8 tracing shifts. Without the bound some
    no-change fits of sample 2 put the empty pulse at the tissue delay. For sample 2 also no change with a free
    (a not capped at 1): a stronger, earlier tissue pulse, the reading that needs no thinner path.
(2) The lead of the second pulse over the empty-cuvette pulse in the fits with the causal loss term (full fit,
    a <= 1, s free and s = 0) of the H1-H3 scans, with H1's and with H4's empty-cuvette pulse as K.
Writes results/sim76/sim76.{txt,json}.
"""
import itertools
import json
from pathlib import Path

import numpy as np

from thznerve.runlog import RunLog
from thznerve.waveform_fit import TwoPulse, effect_size, load_scans, multistart

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim76"
T2 = 0.2


def load(name, file):
    import importlib.util
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / "exp" / file)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


S73 = load("s73b", "sim73_h4_fit.py")
A = S73.A


def free_a(W, valid, t2_max):
    """Sample 2, no change, a free: (q, residual) alone and with a drift."""
    tp = TwoPulse(A.TT, W["empty"], f_loss=2.0)
    fit = multistart(lambda y, x, **kw: tp.fit(y, x, t2_max=t2_max, **kw))
    y, x0 = np.where(valid, W["s2_6min"], 0.0), np.where(valid, W["s2_none"], 0.0)
    out = {}
    for mode, sv in (("alone", 0.0), ("drift", None)):
        q, S, _ = fit(y, x0, s=sv)
        out[mode] = dict(q=[float(v) for v in q], resid=float(np.sqrt(S / np.sum(y ** 2))))
    return out


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    lines, js = [], {}
    say = lambda s="": (print(s, flush=True), lines.append(s))
    with RunLog("sim76", "round 35: H4 refitted with the second pulse within 0.2 ps of the empty-cuvette delay; second-pulse lead with the loss term, H1's and H4's empty pulse") as log:
        W, valid = A.traces()
        base = S73.run(W, valid, t2_max=T2); fa = free_a(W, valid, T2)
        spread, fsp = [], []
        for sh in itertools.product((-1, 1), repeat=3):
            d = {"empty": sh[0], "s1_none": sh[1], "s2_none": sh[1], "s1_6min": sh[2], "s2_6min": sh[2]}
            Wd, vd = A.traces(d)
            spread.append(S73.run(Wd, vd, t2_max=T2)); fsp.append(free_a(Wd, vd, T2))
        js["h4"] = {"base": base, "spread": spread, "free_a": fa, "free_a_spread": fsp}
        say(f"(1) H4, second pulse within {T2} ps of the empty-cuvette delay (residual = rms misfit / rms trace; D against the full fit, %)")
        for s in ("s1", "s2"):
            say(f"  sample {s[1]}, 6 min at 100 V")
            for k, v in base[s].items():
                sp = [x[s][k] for x in spread]
                say(f"    {k:16s} residual {100 * v['resid']:5.1f} % ({100 * min(x['resid'] for x in sp):.1f}-{100 * max(x['resid'] for x in sp):.1f} over tracing shifts)" +
                    (f"  D {v['D']:5.1f} % ({min(x['D'] for x in sp):.1f}-{max(x['D'] for x in sp):.1f})  size {v['size']}" if "D" in v else "") + f"  {S73.fmt(v['q'])}")
        for mode in ("alone", "drift"):
            v = fa[mode]; sp = [x[mode] for x in fsp]
            say(f"  sample 2, no change, a free, {mode:5s}: residual {100 * v['resid']:5.1f} % ({100 * min(x['resid'] for x in sp):.1f}-{100 * max(x['resid'] for x in sp):.1f})  {S73.fmt(v['q'])}")
        say()
        say("(2) second pulse relative to the empty-cuvette pulse in the full fit (causal loss term, a <= 1), ps; negative = early")
        tt, K1, SC = load_scans(ROOT)
        h4 = np.load(ROOT / "results" / "sim73" / "digitised.npz")
        js["lead"] = {}
        for KK, nm in ((K1, "H1 empty pulse"), (h4["empty"], "H4 empty pulse")):
            tp = TwoPulse(tt, KK, f_loss=2.0)
            fit = multistart(lambda y, x, **kw: tp.fit(y, x, a_max=1.0, **kw))
            for sc in ("30 s (H1)", "2 min (H2)", "4 min (H1)", "5 min (H3)"):
                y, x0 = SC[sc]
                qa = fit(y, x0, s=0.0, loss="minphase")[0]; qd = fit(y, x0, loss="minphase")[0]
                js["lead"][f"{nm} | {sc}"] = dict(alone=float(qa[2]), drift=float(qd[2] + qd[3]), b=[float(qa[1]), float(qd[1])])
                say(f"  {nm}, {sc}: s = 0: b {qa[1]:.3f}, at {qa[2]:+.3f}; s free: b {qd[1]:.3f}, s {qd[2]:+.3f}, at {qd[2] + qd[3]:+.3f}")
        L = js["lead"]
        log.outcome = (f"s2 no change drift {100 * base['s2']['none | drift']['resid']:.1f} %, uniform drift {100 * base['s2']['uniform | drift']['resid']:.1f} %, "
                       f"a free drift {100 * fa['drift']['resid']:.1f} %; loss-fit lead H1 K " +
                       f"{min(v['drift'] for k, v in L.items() if k.startswith('H1')):+.3f}..{max(v['drift'] for k, v in L.items() if k.startswith('H1')):+.3f}, H4 K " +
                       f"{min(v['drift'] for k, v in L.items() if k.startswith('H4')):+.3f}..{max(v['drift'] for k, v in L.items() if k.startswith('H4')):+.3f}")
    (OUT / "sim76.json").write_text(json.dumps(js, default=float))
    (OUT / "sim76.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
