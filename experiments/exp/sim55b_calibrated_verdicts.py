"""Sim 55b - Calibrated thresholds for the waveform F tests, and the verdicts of Sims 53-54 under them.

Sim 55 (B) found that on block-bootstrap 'no change' versions of the 30 s scan the nominal 5 % tests reject
far too often (shift test 31 %, no-change test 12 %): the misfit is systematic, not noise, and N_eff from
the residual autocorrelation does not cover it. Here the bootstrap p values are stored for blocks of one and
two correlation lengths (300 scans each), for the 1-dof shift test (F(1, N_eff - 5)) and the 2-dof
no-change test (F(2, N_eff - 5)). The calibrated threshold is the 5th percentile of the bootstrap p: a test
then rejects at a true 5 % rate. The 2-dof threshold is applied to the candidate tests of Sims 53 and 54
(F(2, N_eff - 6), reference N_eff), the 1-dof threshold to the shift tests of Sim 55 (A), and the observed
no-change p values are converted to a calibrated p (fraction of bootstrap p at or below them).
Writes results/sim55/sim55b.txt and results/sim55/sim55b.json.
"""
import json
from pathlib import Path

import numpy as np

from thznerve.runlog import RunLog
from thznerve.waveform_fit import TwoPulse, f_test, load_scans, n_eff

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim55"
lines = []


def say(s=""):
    print(s, flush=True); lines.append(s)


tt, K, SCANS = load_scans(ROOT)
LATE = ("2 min (H2)", "4 min (H1)", "5 min (H3)")
out = {}
with RunLog("sim55b", "round 23: calibrated thresholds of the waveform F tests (block bootstrap of the 30 s scan); Sims 53-54 reclassified") as log:
    tp = TwoPulse(tt, K, f_loss=2.0)
    floor = tp.loss_transfer(2.07, "zero")
    fit = lambda y, x, **kw: tp.fit(y, x, a_max=1.0, loss=kw.pop("loss", None), **kw)
    y, x0 = SCANS["30 s (H1)"]
    xf = tp.apply_transfer(x0, floor)
    pr, m, res = fit(y, xf, loss="minphase")
    base = tp.model(pr, xf, loss="minphase")
    _, corr = n_eff(res)
    n = len(res); rng = np.random.default_rng(56)
    for mult in (1, 2):
        blk = max(mult * corr, 2)
        p1, p2 = [], []
        for _ in range(300):
            starts = rng.integers(0, n, size=n // blk + 1)
            yb = base + np.concatenate([np.take(res, np.arange(s0, s0 + blk), mode="wrap") for s0 in starts])[:n]
            _, mb, resb = fit(yb, xf, loss="minphase", p0=pr)
            neb, _ = n_eff(resb)
            _, mb0, _ = fit(yb, xf, s=0.0, loss="minphase")
            _, mbn, _ = fit(yb, xf, s=0.0)
            p1.append(f_test(mb0, mb, neb, 5, 1)[1]); p2.append(f_test(mbn, mb, neb, 5, 2)[1])
        p1, p2 = np.sort(p1), np.sort(p2)
        out[f"block x{mult}"] = dict(block=int(blk), thr1=float(np.quantile(p1, 0.05)), thr2=float(np.quantile(p2, 0.05)),
                                     rej1=float(np.mean(p1 < 0.05)), rej2=float(np.mean(p2 < 0.05)), p1=p1.tolist(), p2=p2.tolist())
        say(f"block {blk} samples ({mult} correlation length{'s' if mult > 1 else ''}): nominal 5 % rejects {np.mean(p1 < 0.05):.2f} "
            f"(shift) and {np.mean(p2 < 0.05):.2f} (no change); calibrated thresholds p < {out[f'block x{mult}']['thr1']:.4f} "
            f"(shift), p < {out[f'block x{mult}']['thr2']:.4f} (2-dof)")
    (OUT / "sim55b.json").write_text(json.dumps(out, indent=1))
    thr1 = min(v["thr1"] for v in out.values()); thr2 = min(v["thr2"] for v in out.values())
    p2all = np.sort(np.r_[out["block x1"]["p2"], out["block x2"]["p2"]])
    say(f"using the stricter: shift p < {thr1:.4f}, 2-dof p < {thr2:.4f}")

    s55 = json.loads((OUT / "sim55.json").read_text())["A"]
    refs = [v for v in s55 if not v.startswith("summary")]
    pn = [s55[v][sc]["p_none"] for v in refs for sc in LATE]
    cal = [float(np.mean(p2all <= q)) for q in pn]
    say(f"no change at all, 36 references x 3 scans: nominal p <= {max(pn):.4f}; calibrated p <= {max(cal):.4f} "
        f"(bootstrap floor {1 / len(p2all):.4f})")
    ok30 = [v for v in refs if s55[v]["30 s (H1)"]["p_s0"] >= thr1]
    say(f"references passing the 30 s control at the calibrated threshold: {len(ok30)} of {len(refs)}; s = 0 rejected at "
        f"2-5 min in {sum(s55[v][sc]['p_s0'] < thr1 for v in ok30 for sc in LATE)} of {3 * len(ok30)}")
    out["no_change_cal_p_max"] = max(cal)

    summ = {}
    for sim, f in (("53", ROOT / "results" / "sim53" / "sim53.json"), ("54", ROOT / "results" / "sim54" / "sim54.json")):
        d = json.loads(f.read_text())
        vs = [v for v in d if v.startswith("increment")]
        names = sorted({c for v in vs for c in d[v]["2 min (H2)"]["cands"] if not c.startswith("node ")})
        has_nodes = any(c.startswith("node ") for c in d[vs[0]]["2 min (H2)"]["cands"])
        for c in names + (["node states (best of the three screened)"] if has_nodes else []):
            row = []
            for sc in LATE:
                if c.startswith("node states"):
                    ps = [max(x["p"] for n, x in d[v][sc]["cands"].items() if n.startswith("node ")) for v in vs]
                else:
                    ps = [d[v][sc]["cands"][c]["p"] for v in vs]
                row.append(sum(q >= thr2 for q in ps))
            summ.setdefault(sim, {})[c] = row
            say(f"  Sim {sim} {c:42s} references passing at 2/4/5 min (of 6): {row}")
    out["calibrated_counts"] = summ
    log.outcome = f"thresholds shift {thr1:.4f}, 2-dof {thr2:.4f}"
(OUT / "sim55b.json").write_text(json.dumps(out, indent=1))
(OUT / "sim55b.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
