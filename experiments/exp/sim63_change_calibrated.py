"""Sim 63 - Review round 27 follow-up: is there a change of the tissue pulse at all, with the calibration of Sim 61?

Sim 61 found that multi-start fits on the synthetic scans make the nominal F tests reject their null 10-25 % of
the time at 0.05, far more than the single-start bootstrap of Sim 57 showed. The no-change and shift-alone tests
of Sim 57 carry the paper's claim that the biased waveforms show a change of the tissue pulse, so they are
recalibrated here the same way (thznerve.waveform_fit.null_calibrated_test with multistart fits), on the eight
references of Sim 59 and all four scans:
  no change   null: the reference pulse, s = 0, no loss; alternative: loss increment and s free; F(2).
  shift alone null: the reference pulse, s free, no loss; alternative: loss increment and s free; F(1).
Up to 60 synthetic scans each, in batches of 20, stopping once the 99 % interval of the fraction lies wholly
above or below 0.05. A change is shown where p_cal < 0.05. Four worker processes; resumes from a partial
results/sim63/sim63.json. Writes results/sim63/sim63.txt and results/sim63/sim63.json.
"""
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from thznerve.runlog import RunLog
from thznerve.waveform_fit import TwoPulse, load_scans, multistart, null_calibrated_test

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim63"
NB = int(sys.argv[1]) if len(sys.argv) > 1 else 60
tt, K, SCANS = load_scans(ROOT)
REFS = [f"floor {k:g}, " + ab for k in (1.0, 1.5, 2.07, 2.5) for ab in ("a <= 1", "a + b <= 1")]
tp = TwoPulse(tt, K, f_loss=2.0)


def run(case):
    vn, sc, test = case
    k0 = float(vn.split()[1].rstrip(",")); ab = 1.0 if "a + b" in vn else None
    fit = multistart(lambda y, x, **kw: tp.fit(y, x, a_max=1.0, ab_max=ab, **kw))
    y, x0 = SCANS[sc]
    xf = tp.apply_transfer(x0, tp.loss_transfer(k0, "zero"))
    d = null_calibrated_test(fit, y, xf, [xf], 0.0 if test == "no change" else None, 2 if test == "no change" else 1,
                             nb=NB, batch=20, rng=np.random.default_rng(63), model=lambda q, x, s: tp.model(q, x, s=s))
    return case, d


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    part = OUT / "sim63.json"
    res = json.loads(part.read_text()).get("cases", {}) if part.exists() else {}
    with RunLog("sim63", f"round 27: no-change and shift-alone tests recalibrated with multistart fits, up to {NB} scans, 8 references x 4 scans") as log:
        cs = [(vn, sc, t) for t in ("no change", "shift alone") for sc in SCANS for vn in REFS
              if " | ".join((vn, sc, t)) not in res]
        print(f"{len(res)} cases done, {len(cs)} to run", flush=True)
        with Pool(4) as pool:
            for (vn, sc, t), d in pool.imap_unordered(run, cs):
                res[f"{vn} | {sc} | {t}"] = d
                print(f"{vn:22s} {sc[:5]} {t:11s}: nominal p {d['p_nominal']:.4f}, calibrated p {d['p_cal']:.3f} "
                      f"({d['nb']} scans; null rejects {d['rej_nominal']:.2f} at nominal 0.05)", flush=True)
                part.write_text(json.dumps({"cases": res}, indent=1, default=float))
        lines = [f"{k}: nominal p {d['p_nominal']:.4f}, calibrated p {d['p_cal']:.3f} ({d['nb']} scans), null rejects "
                 f"{d['rej_nominal']:.2f} at nominal 0.05" for k, d in sorted(res.items())] + [""]
        summ = {}
        for t in ("no change", "shift alone"):
            for sc in SCANS:
                ds = [res[f"{vn} | {sc} | {t}"] for vn in REFS]
                summ[f"{t} | {sc}"] = dict(rejected=sum(d["p_cal"] < 0.05 for d in ds), of=len(ds),
                                           p_cal=[min(d["p_cal"] for d in ds), max(d["p_cal"] for d in ds)],
                                           rej_nominal=[min(d["rej_nominal"] for d in ds), max(d["rej_nominal"] for d in ds)])
                s = summ[f"{t} | {sc}"]
                lines.append(f"{t:11s} {sc:11s}: rejected in {s['rejected']}/{s['of']} references, calibrated p "
                             f"{s['p_cal'][0]:.3f}-{s['p_cal'][1]:.3f}; null rejects {s['rej_nominal'][0]:.2f}-"
                             f"{s['rej_nominal'][1]:.2f} at nominal 0.05")
        print("\n".join(lines[-8:]), flush=True)
        log.outcome = "; ".join(f"{k.split(' | ')[0][:8]} {k.split(' | ')[1][:5]} {v['rejected']}/8" for k, v in summ.items())
    part.write_text(json.dumps({"cases": res, "summary": summ}, indent=1, default=float))
    (OUT / "sim63.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
