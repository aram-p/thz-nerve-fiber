"""Sim 67 - Do two causes together reproduce the change of the tissue pulse?

Sim 66 judged each candidate alone by its effect size D (rms of what it leaves unexplained, in % of the trace rms,
against the fit with the causal loss increment and a free shift) and found none within the 1.63 / 1.18 % that
separates traced unbiased pulses. The three closest were a thicker tissue path, gas in a slice that swells
(r 30 um) and node states (wide range). Here each pair of them is applied together (product of the two
transfers) on a grid of both sizes, alone (s = 0) and with a drift (s free), for the eight references of Sim 59
and the scans 2, 4 and 5 min; node states are the one Sim 66 found best for that reference, scan and mode.
Resumes from results/sim67/sim67_partial.json. Writes results/sim67/sim67.txt and results/sim67/sim67.json.
"""
import json
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from thznerve.runlog import RunLog
from thznerve.waveform_fit import multistart

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim67"
LATE = ("2 min (H2)", "4 min (H1)", "5 min (H3)")
GAS = "gas, slice swells, r 30 um"
THICK = "thicker path (fraction)"


def s66():
    import importlib.util
    spec = importlib.util.spec_from_file_location("s66", ROOT / "experiments" / "exp" / "sim66_effect_size.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


S = s66()
R66 = json.loads((ROOT / "results" / "sim66" / "sim66.json").read_text())
GR = {THICK: np.linspace(0, 0.3, 11), GAS: np.linspace(0, 0.4, 11), "node": np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6.0])}
PAIRS = [(THICK, GAS), (THICK, "node"), (GAS, "node")]


def run(job):
    vn, sc = job
    tp, M = S.tp, S.M
    k0 = float(vn.split()[1].rstrip(",")); ab = 1.0 if "a + b" in vn else None
    fit = multistart(lambda y, x, **kw: tp.fit(y, x, a_max=1.0, ab_max=ab, **kw))
    y, x0 = S.SCANS[sc]
    xf = tp.apply_transfer(x0, tp.loss_transfer(k0, "zero"))
    Sf = fit(y, xf, loss="minphase")[1]
    case = R66["cases"][f"{vn} | {sc}"]
    row = {}
    for mode in ("alone", "drift"):
        node = case["node | wide"][mode]["name"]
        tr = {c: [M.transfer(node if c == "node" else c, v) if v else None for v in GR[c]] for c in GR}
        for c1, c2 in PAIRS:
            best = None
            for i, v1 in enumerate(GR[c1]):
                for j, v2 in enumerate(GR[c2]):
                    t = np.ones(len(tp.fr), complex)
                    for tt_ in (tr[c1][i], tr[c2][j]):
                        if tt_ is not None:
                            t = t * tt_
                    d = S.D(fit(y, tp.apply_transfer(xf, t), s=0.0 if mode == "alone" else None)[1], Sf, y)
                    if best is None or d < best[0]:
                        best = (d, float(v1), float(v2))
            row[f"{c1[:12]} + {c2[:12]} | {mode}"] = dict(D=best[0], sizes=best[1:], node=node)
    return job, row


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    T = R66["yardstick"]
    with RunLog("sim67", "pairs of the three closest candidates (thicker path, gas r 30, node wide) by effect size, 8 references, 2-5 min") as log:
        part = OUT / "sim67_partial.json"
        res = json.loads(part.read_text()) if part.exists() else {}
        jobs = [(vn, sc) for vn in S.REFS8 for sc in LATE if f"{vn} | {sc}" not in res]
        print(f"{len(res)} done, {len(jobs)} to run", flush=True)
        with Pool(4) as pool:
            for (vn, sc), row in pool.imap_unordered(run, jobs):
                res[f"{vn} | {sc}"] = row
                part.write_text(json.dumps(res, indent=1, default=float))
                print(f"{vn} | {sc}: " + "; ".join(f"{k} {d['D']:.2f}@{d['sizes'][0]:g},{d['sizes'][1]:g}" for k, d in row.items()), flush=True)
        lines = [f"{k}: " + "; ".join(f"{p} D {d['D']:.2f} at {d['sizes'][0]:g}, {d['sizes'][1]:g}" for p, d in row.items())
                 for k, row in sorted(res.items())] + [""]
        summ = {}
        for sc in LATE:
            for p in res[f"{S.REFS8[0]} | {sc}"]:
                ds = [res[f"{vn} | {sc}"][p]["D"] for vn in S.REFS8]
                mode = p.split(" | ")[1]
                n = sum(d <= T[mode] for d in ds)
                summ[f"{sc} | {p}"] = dict(reproduces=n, of=8, D=[min(ds), max(ds)])
                lines.append(f"{sc:11s} {p:40s}: reproduces in {n}/8 (D {min(ds):.2f}-{max(ds):.2f} %; yardstick {T[mode]:.2f} %)")
        print("\n".join(lines[-len(summ):]), flush=True)
        log.outcome = "best pair D " + ", ".join(f"{sc[:5]} {min(v['D'][0] for k, v in summ.items() if k.startswith(sc)):.2f} %" for sc in LATE) + \
            f"; reproduces in {sum(v['reproduces'] for v in summ.values())} of {8 * len(summ)}"
    (OUT / "sim67.json").write_text(json.dumps({"cases": res, "summary": summ, "yardstick": T}, indent=1, default=float))
    (OUT / "sim67.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
