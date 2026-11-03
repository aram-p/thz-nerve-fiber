"""Sim 62 - Review round 27, item on the second-pulse lead: where the second pulse arrives in the fits that pass.

For every case of Sim 61 (eight references, floors 1.0-2.5 dB/THz, a <= 1 or a + b <= 1) whose candidate passes
the calibrated test (p_cal >= 0.05), the candidate is fitted at its size (alone: s = 0; with a drift: s free,
no loss term; multi-start) and the lead of the second pulse over the empty-cuvette pulse, -(s + tau) ps, is
recorded; so is the lead in each reference fit (loss increment, s free). Reported per candidate and scan as a
range with the references it comes from, against the 0.015 ps to which the no-voltage minima of the three
reports agree. Writes results/sim62/sim62.txt and results/sim62/sim62.json.
"""
import importlib.util
import json
from pathlib import Path

from thznerve.runlog import RunLog
from thznerve.waveform_fit import multistart

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim62"
spec = importlib.util.spec_from_file_location("s61", ROOT / "experiments" / "exp" / "sim61_calibrate_table2.py")
s61 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s61)
tp, SCANS, S59, transfer = s61.tp, s61.SCANS, s61.S59, s61.transfer
R61 = json.loads((ROOT / "results" / "sim61" / "sim61.json").read_text())["cases"]

OUT.mkdir(parents=True, exist_ok=True)
lines, out = [], {"reference": {}, "passing": {}}


def say(s=""):
    print(s, flush=True); lines.append(s)


def fit_for(vn):
    ab = 1.0 if "a + b" in vn else None
    return multistart(lambda y, x, **kw: tp.fit(y, x, a_max=1.0, ab_max=ab, **kw))


with RunLog("sim62", "round 27: second-pulse lead in the fits that pass the calibrated test (Sim 61), 8 references") as log:
    for vn in (v for v in S59 if v.startswith("floor")):
        fit = fit_for(vn); k0 = float(vn.split()[1].rstrip(","))
        for sc in s61.LATE:
            y, x0 = SCANS[sc]
            xf = tp.apply_transfer(x0, tp.loss_transfer(k0, "zero"))
            q, _, _ = fit(y, xf, loss="minphase")
            out["reference"][f"{vn} | {sc}"] = float(-(q[2] + q[3]))
    for key, d in R61.items():
        if d["p_cal"] < 0.05:
            continue
        vn, sc, c, mode = key.split(" | ")
        k0 = float(vn.split()[1].rstrip(","))
        y, x0 = SCANS[sc]
        xf = tp.apply_transfer(x0, tp.loss_transfer(k0, "zero"))
        xt = tp.apply_transfer(xf, transfer(c, d["size"]))
        q, _, _ = fit_for(vn)(y, xt, s=0.0) if mode == "alone" else fit_for(vn)(y, xt)
        out["passing"][key] = dict(lead_ps=float(-q[2] if mode == "alone" else -(q[2] + q[3])), p_cal=d["p_cal"],
                                   size=d["size"])
    say("Reference fits (loss increment, s free), lead in ps: " + "; ".join(
        f"{sc[:5]} {min(v for k, v in out['reference'].items() if sc in k):+.3f} to "
        f"{max(v for k, v in out['reference'].items() if sc in k):+.3f}" for sc in s61.LATE))
    groups = {}
    for key, d in out["passing"].items():
        vn, sc, c, mode = key.split(" | ")
        g = "node states" if c.startswith("node ") else c
        groups.setdefault((g, mode, sc), []).append((d["lead_ps"], vn))
    for (g, mode, sc), xs in sorted(groups.items()):
        say(f"{g[:28]:28s} {mode:5s} {sc[:5]}: {len(xs)} passing fits, lead {min(x for x, _ in xs):+.3f} to "
            f"{max(x for x, _ in xs):+.3f} ps ({'; '.join(v for _, v in xs)})")
    out["groups"] = {" | ".join(k): dict(n=len(xs), lo=min(x for x, _ in xs), hi=max(x for x, _ in xs),
                                         refs=[v for _, v in xs]) for k, xs in groups.items()}
    leads = [d["lead_ps"] for d in out["passing"].values()]
    log.outcome = (f"{len(leads)} passing fits; lead {min(leads):+.3f} to {max(leads):+.3f} ps" if leads else "no passing fits")
(OUT / "sim62.json").write_text(json.dumps(out, indent=1, default=float))
(OUT / "sim62.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
