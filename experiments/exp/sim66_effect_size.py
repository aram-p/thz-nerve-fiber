"""Sim 66 - The candidates judged by effect size against the differences between traced unbiased pulses.

Sim 65 showed that significance tests on the traced waveforms are not usable: a null that keeps the residual's
spectrum rejects 'no change' between the unbiased pulses of different reports in 10 of 12 comparisons. Every model
leaves a systematic misfit, so a large enough number of points finds any difference significant. Here the
criterion is the size of what a fit leaves unexplained:
  D = sqrt(max(S_restricted - S_full, 0) / N) / rms(y), in % of the trace rms,
where S is the summed squared misfit, S_full the fit with the causal loss increment and a free shift (no
candidate), and S_restricted the fit under test: no change (s = 0, no loss), a shift alone (s free, no loss), or
a candidate at its best size on a grid without the loss term, s = 0 (alone) or s free (with a drift). For linear
models D is the rms of the difference between the two fitted waveforms.
Yardstick: the same D for the six ordered pairs of traced unbiased pulses (H1, H2, H3; each fitted with another,
a <= 1 and a free), which have no tissue change but include tracing error and any difference between the three
reports' measurements. A restricted fit 'reproduces the change' if its D is no larger than the largest D of
those pairs in the same mode. References: the eight of Sim 59 (floors 1.0-2.5 dB/THz, increment held above 2 THz,
a <= 1 or a + b <= 1); scans 30 s, 2, 4, 5 min; candidates and grids as in Sim 64 (plausible and wide ranges,
node states: the best of the three screened in Sim 59 for that reference and scan). Multi-start fits.
Writes results/sim66/sim66.txt and results/sim66/sim66.json.
"""
import json
from itertools import permutations
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from thznerve.runlog import RunLog
from thznerve.waveform_fit import TwoPulse, load_scans, multistart

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim66"
tt, K, SCANS = load_scans(ROOT)
tp = TwoPulse(tt, K, f_loss=2.0)
S59 = json.loads((ROOT / "results" / "sim59" / "sim59.json").read_text())
REFS8 = [f"floor {k:g}, " + ab for k in (1.0, 1.5, 2.07, 2.5) for ab in ("a <= 1", "a + b <= 1")]


def s64():
    import importlib.util
    spec = importlib.util.spec_from_file_location("s64", ROOT / "experiments" / "exp" / "sim64_size_fitted_tests.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


M = s64()


def D(S_r, S_f, y):
    return 100 * np.sqrt(max(S_r - S_f, 0.0) / len(y)) / np.sqrt(np.mean(y ** 2))


def tracing_pairs():
    ref = {"H1": SCANS["30 s (H1)"][1], "H2": SCANS["2 min (H2)"][1], "H3": SCANS["5 min (H3)"][1]}
    out = {}
    for a, b in permutations(ref, 2):
        for amax in (1.0, np.inf):
            fit = multistart(lambda y, x, **kw: tp.fit(y, x, a_max=amax, **kw))
            y, x = ref[a], ref[b]
            Sf = fit(y, x, loss="minphase")[1]
            out[f"{a} from {b}, a <= {amax:g}"] = dict(alone=D(fit(y, x, s=0.0)[1], Sf, y), drift=D(fit(y, x)[1], Sf, y))
    return out


def run(job):
    vn, sc = job
    k0 = float(vn.split()[1].rstrip(",")); ab = 1.0 if "a + b" in vn else None
    fit = multistart(lambda y, x, **kw: tp.fit(y, x, a_max=1.0, ab_max=ab, **kw))
    y, x0 = SCANS[sc]
    xf = tp.apply_transfer(x0, tp.loss_transfer(k0, "zero"))
    Sf = fit(y, xf, loss="minphase")[1]
    row = {"none": dict(alone=D(fit(y, xf, s=0.0)[1], Sf, y), drift=D(fit(y, xf)[1], Sf, y))}
    if sc.startswith("30 s"):
        return job, row
    nodes = sorted((n for n in S59[vn][sc] if n.startswith("node ")), key=lambda n: -S59[vn][sc][n]["p_alone"])
    for c, (g, w) in M.GRIDS.items():
        for rn, grid in (("plausible", g), ("wide", w)):
            if grid is None:
                continue
            names = nodes if c == "node" else [c]
            best = {}
            for name in names:
                for v in grid:
                    xt = xf if v == 0 else tp.apply_transfer(xf, M.transfer(name, v))
                    for mode in ("alone", "drift"):
                        d = D(fit(y, xt, s=0.0 if mode == "alone" else None)[1], Sf, y)
                        if mode not in best or d < best[mode][0]:
                            best[mode] = (d, float(v), name)
            row[f"{c} | {rn}"] = {m: dict(D=b[0], size=b[1], name=b[2]) for m, b in best.items()}
    return job, row


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    with RunLog("sim66", "effect size of each restricted fit vs the traced-unbiased-pulse yardstick, 8 references, 4 scans, plausible and wide sizes") as log:
        tr = tracing_pairs()
        T = {m: max(d[m] for d in tr.values()) for m in ("alone", "drift")}
        lines = ["Yardstick, traced unbiased pulses (D in % of trace rms, no change / shift alone):"]
        lines += [f"  {k:22s} {d['alone']:5.2f} / {d['drift']:5.2f}" for k, d in tr.items()]
        lines.append(f"  largest: {T['alone']:.2f} (no change, s = 0) and {T['drift']:.2f} (s free)")
        print("\n".join(lines), flush=True)
        jobs = [(vn, sc) for vn in REFS8 for sc in SCANS]
        res = {}
        with Pool(4) as pool:
            for (vn, sc), row in pool.imap_unordered(run, jobs):
                res[f"{vn} | {sc}"] = row
                print(f"{vn} | {sc}: none {row['none']['alone']:.2f}/{row['none']['drift']:.2f}", flush=True)
        lines.append("")
        for k in sorted(res):
            row = res[k]
            lines.append(f"{k}: no change {row['none']['alone']:.2f}, shift alone {row['none']['drift']:.2f}" + "".join(
                f"; {c[:14]} {d['alone']['D']:.2f}@{d['alone']['size']:g} / {d['drift']['D']:.2f}@{d['drift']['size']:g}"
                for c, d in row.items() if c != "none"))
        lines.append("")
        summ = {}
        for sc in SCANS:
            rows = [res[f"{vn} | {sc}"] for vn in REFS8]
            for c in rows[0]:
                for mode in ("alone", "drift"):
                    ds = [r[c][mode] if c == "none" else r[c][mode]["D"] for r in rows]
                    n = sum(d <= T[mode] for d in ds)
                    summ[f"{c} | {mode} | {sc}"] = dict(reproduces=n, of=len(ds), D=[min(ds), max(ds)])
                    lines.append(f"{sc:11s} {c:42s} {mode:5s}: reproduces in {n}/{len(ds)} (D {min(ds):.2f}-{max(ds):.2f} %; "
                                 f"yardstick {T[mode]:.2f} %)")
        print("\n".join(lines[-len(summ):]), flush=True)
        log.outcome = f"yardstick {T['alone']:.2f}/{T['drift']:.2f} %; " + "; ".join(
            f"{k.split(' | ')[0][:10]} {k.split(' | ')[-1][:5]} {v['reproduces']}" for k, v in summ.items() if "| alone |" in k and "plausible" in k)[:200]
    (OUT / "sim66.json").write_text(json.dumps({"tracings": tr, "yardstick": T, "cases": res, "summary": summ}, indent=1, default=float))
    (OUT / "sim66.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
