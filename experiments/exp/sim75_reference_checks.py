"""Sim 75 - Review round 34: is H4's sample 1 a new scan, and is the lead of the second pulse real?

(1) Every pair of traced scans (H1-H4) compared after the best time shift (+-0.1 ps) and scale: rms difference
    over the window where H4 is traced, in % of the first trace's rms.
(2) The two-pulse fit of the H1-H3 biased scans (no loss term, a <= 1, alone and with a drift) with H1's empty-
    cuvette pulse and with H4's: where the second pulse sits relative to the empty-cuvette pulse (tau alone,
    s + tau with a drift; negative = early).
Writes results/sim75/sim75.txt.
"""
from pathlib import Path

import numpy as np

from thznerve.runlog import RunLog
from thznerve.waveform_fit import TwoPulse, load_scans, multistart

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim75"


def best_match(x, y, tt, v):
    best = None
    for sh in np.arange(-0.1, 0.1001, 0.005):
        ys = np.interp(tt - sh, tt, y, left=0, right=0)
        k = np.dot(x[v], ys[v]) / np.dot(ys[v], ys[v])
        r = np.sqrt(np.mean((x[v] - k * ys[v]) ** 2) / np.mean(x[v] ** 2))
        if best is None or r < best[0]:
            best = (r, k, sh)
    return best


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    lines = []
    say = lambda s="": (print(s, flush=True), lines.append(s))
    with RunLog("sim75", "round 34: H4 sample 1 against H1-H3 scans; second-pulse lead with H1's and H4's empty pulses") as log:
        tt, K1, SC = load_scans(ROOT)
        h4 = np.load(ROOT / "results" / "sim73" / "digitised.npz"); v = h4["valid"]
        tr = {"H1 empty": K1, "H4 empty": h4["empty"], "H1-H3 unbiased (H1)": SC["4 min (H1)"][1],
              "H3 unbiased": SC["5 min (H3)"][1], "H2 unbiased": SC["2 min (H2)"][1], "H4 sample 1 unbiased": h4["s1_none"],
              "H1 30 s": SC["30 s (H1)"][0], "H2 2 min": SC["2 min (H2)"][0], "H1 4 min": SC["4 min (H1)"][0],
              "H3 5 min": SC["5 min (H3)"][0], "H4 sample 1, '6 min'": h4["s1_6min"]}
        pairs = [("H1 empty", "H4 empty"), ("H1-H3 unbiased (H1)", "H3 unbiased"), ("H1-H3 unbiased (H1)", "H2 unbiased"),
                 ("H1-H3 unbiased (H1)", "H4 sample 1 unbiased"), ("H3 unbiased", "H4 sample 1 unbiased"),
                 ("H1 30 s", "H2 2 min"), ("H2 2 min", "H1 4 min"), ("H1 4 min", "H3 5 min"),
                 ("H1 4 min", "H4 sample 1, '6 min'"), ("H3 5 min", "H4 sample 1, '6 min'")]
        say("(1) scans compared after the best shift and scale (window where H4 is traced)")
        res = {}
        for a, b in pairs:
            r, k, sh = best_match(tr[a], tr[b], tt, v); res[(a, b)] = r
            say(f"  {a:22s} vs {b:22s}: {100 * r:5.1f} % (scale {k:.3f}, shift {sh:+.3f} ps)")
        say("")
        say("(2) second pulse relative to the empty-cuvette pulse (ps; negative = early), no loss term, a <= 1")
        for KK, nm in ((K1, "H1 empty pulse"), (h4["empty"], "H4 empty pulse")):
            tp = TwoPulse(tt, KK, f_loss=2.0)
            fit = multistart(lambda y, x, **kw: tp.fit(y, x, a_max=1.0, **kw))
            for sc in ("30 s (H1)", "2 min (H2)", "4 min (H1)", "5 min (H3)"):
                y, x0 = SC[sc]
                qa = fit(y, x0, s=0.0)[0]; qd = fit(y, x0)[0]
                say(f"  {nm}, {sc}: alone b {qa[1]:.3f}, at {qa[2]:+.3f}; with a drift b {qd[1]:.3f}, s {qd[2]:+.3f}, second pulse at {qd[2] + qd[3]:+.3f}")
        log.outcome = f"H4 '6 min' vs H3 5 min {100 * res[('H3 5 min', 'H4 sample 1, ' + chr(39) + '6 min' + chr(39))]:.1f} %; other biased pairs " \
                      f"{100 * min(res[p] for p in pairs[5:8]):.0f}-{100 * max(res[p] for p in pairs[5:8]):.0f} %"
    (OUT / "sim75.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
