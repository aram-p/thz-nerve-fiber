"""Sim 69 - Review round 30: the Sim 68 threshold with the noise level corrected, more draws, and the paper's rule
applied in the script.

Round 30 found: (1) the Sim 68 script judged at the 95th percentile of 30 draws while the paper used the largest
draw; (2) the full fit (5 parameters, N_eff 13-20) leaves a residual smaller than the noise it fits, so the draws
had ~12-15 % too little noise; (3) the largest of 30 draws is itself noisy. Here, for every restricted fit of
Sim 68 (same best sizes, node states and references):
  noise    the full fit's residual, scaled up so that the full fit leaves as much as the data's residual when
           noise of that size is injected into its own curve (10 calibration draws per case; the scale is
           compared with sqrt(N_eff / (N_eff - 5))), phase-randomised, added to the restricted fit's curve; both
           models refitted; NB = 100 draws, all stored (thznerve.waveform_fit.effect_noise).
  rule     per scan the threshold is the 1 - 0.05/3 quantile of the noise D combined in quadrature with the tracing
           D of Sim 68 (pairs built like each reference); a candidate is excluded in a reference if D exceeds the
           threshold at one of 2, 4 and 5 min, so a true cause is excluded by chance in ~5 % of references. The
           rate is also estimated directly from the draws (any of the three scans' draw i over its threshold).
  margin   the largest D / threshold over the three scans; and how often the verdict flips when the 100 draws
           are resampled (200 bootstrap resamples): above 10 % the verdict is called borderline.
Writes results/sim69/sim69.txt and results/sim69/sim69.json; resumes from results/sim69/sim69_partial.json.
"""
import json
import zlib
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from thznerve.runlog import RunLog
from thznerve.waveform_fit import effect_noise, effect_size, multistart, n_eff, phase_surrogate

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim69"
NB = 100
LEVEL = 1 - 0.05 / 3
LATE = ("2 min (H2)", "4 min (H1)", "5 min (H3)")


def s64():
    import importlib.util
    spec = importlib.util.spec_from_file_location("s64", ROOT / "experiments" / "exp" / "sim64_size_fitted_tests.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


M = s64()
tp, SCANS = M.tp, M.SCANS
S68 = json.loads((ROOT / "results" / "sim68" / "sim68.json").read_text())
R68, T68 = S68["cases"], S68["tracing"]


def run(job):
    vn, sc = job
    k0 = float(vn.split()[1].rstrip(",")); ab = 1.0 if "a + b" in vn else None
    fit = multistart(lambda y, x, **kw: tp.fit(y, x, a_max=1.0, ab_max=ab, **kw))
    y, x0 = SCANS[sc]
    xf = tp.apply_transfer(x0, tp.loss_transfer(k0, "zero"))
    qf, Sf, rf = fit(y, xf, loss="minphase")
    ne = n_eff(rf)[0]
    fit_f = lambda yy: fit(yy, xf, loss="minphase")
    rc = np.random.default_rng(zlib.crc32(f"69scale|{vn}|{sc}".encode())); ratios = []
    for _ in range(10):                                   # how much of injected noise the full fit leaves
        e = phase_surrogate(rc, rf); ratios.append(np.sqrt(fit_f(y + rf + e)[1] / np.sum(e ** 2)))
    scale = float(1 / np.mean(ratios))
    row = {"full": dict(ne=float(ne), scale=scale, scale_formula=float(np.sqrt(ne / (ne - 5))))}
    done = {}
    for k, d in R68[f"{vn} | {sc}"].items():
        if k == "full":
            continue
        name, mode = k.split(" | ")
        s = 0.0 if mode == "alone" else None
        if d["size"] == 0:
            xt = xf; key = ("none", mode)
        elif name == "node":
            xt = tp.apply_transfer(xf, M.transfer("node " + d["state"], d["size"])); key = (d["state"], d["size"], mode)
        else:
            xt = tp.apply_transfer(xf, M.transfer(name, d["size"])); key = (name, d["size"], mode)
        if key in done:
            row[k] = dict(done[key], state=d["state"]); continue
        _, Sr, rr = fit(y, xt, s=s)
        rng = np.random.default_rng(zlib.crc32(f"69|{vn}|{sc}|{k}".encode()))
        dn = effect_noise(lambda yy: fit(yy, xt, s=s), fit_f, y + rr, scale * rf, NB, rng)
        row[k] = done[key] = dict(D=effect_size(Sr, Sf, y), size=d["size"], state=d["state"], draws=[float(v) for v in dn])
    return job, row


def thr(draws, t):
    return float(np.hypot(np.quantile(draws, LEVEL), t))


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    part = OUT / "sim69_partial.json"
    res = json.loads(part.read_text()) if part.exists() else {}
    jobs = [(vn, sc) for sc in LATE + ("30 s (H1)",) for vn in M.REFS8 if f"{vn} | {sc}" not in res]
    with RunLog("sim69", f"round 30: noise scale calibrated on the full fit, {NB} draws, per-scan level {LEVEL:.4f}, any-of-three rule, bootstrap of the verdict") as log:
        print(f"{len(res)} done, {len(jobs)} to run", flush=True)
        with Pool(4) as pool:
            for job, row in pool.imap_unordered(run, jobs):
                res[" | ".join(job)] = row
                part.write_text(json.dumps(res, default=float))
                print(" | ".join(job), "done", flush=True)
        rng = np.random.default_rng(69)
        lines = [f"Noise scale (calibrated): {min(r['full']['scale'] for r in res.values()):.3f}-"
                 f"{max(r['full']['scale'] for r in res.values()):.3f} (sqrt(N_eff/(N_eff-5)): "
                 f"{min(r['full']['scale_formula'] for r in res.values()):.3f}-{max(r['full']['scale_formula'] for r in res.values()):.3f}); "
                 f"per-scan level {LEVEL:.4f}; {NB} draws", ""]
        summ = {}
        keys = [k for k in res[f"{M.REFS8[0]} | {LATE[0]}"] if k != "full"]
        for k in keys:
            mode = k.split(" | ")[1]
            per = {}
            for vn in M.REFS8:
                rows = [res[f"{vn} | {sc}"][k] for sc in LATE]
                th = [thr(r["draws"], T68[vn][mode]) for r in rows]
                ratio = [r["D"] / t for r, t in zip(rows, th)]
                exc = max(ratio) > 1
                dr = np.array([r["draws"] for r in rows])                    # 3 x NB
                false = float(np.mean(np.any(dr / np.array(th)[:, None] > 1, axis=0)))
                flips = 0
                for _ in range(200):
                    idx = rng.integers(0, NB, NB)
                    rb = max(r["D"] / thr(np.array(r["draws"])[idx], T68[vn][mode]) for r in rows)
                    flips += (rb > 1) != exc
                per[vn] = dict(ratio=ratio, largest=max(ratio), excluded=exc, false_rate=false, flip=flips / 200,
                               sizes=[r["size"] for r in rows], states=[r["state"] for r in rows])
            n_ex = sum(p["excluded"] for p in per.values()); n_bl = sum(p["flip"] > 0.1 for p in per.values())
            L = [p["largest"] for p in per.values()]; F = [p["false_rate"] for p in per.values()]
            summ[k] = dict(excluded=n_ex, borderline=n_bl, largest=[min(L), max(L)], false_rate=[min(F), max(F)], per=per)
            lines.append(f"{k:40s} excluded in {n_ex}/8 ({n_bl} borderline); largest D/threshold {min(L):.2f}-{max(L):.2f}; "
                         f"false exclusion from draws {min(F):.2f}-{max(F):.2f}")
            for vn, p in per.items():
                lines.append(f"    {vn:22s} " + " / ".join(f"{x:.2f}" for x in p["ratio"]) + f"  flip {p['flip']:.2f}  sizes {p['sizes']}")
        lines.append("")
        for k in ("none | alone", "none | drift"):
            rr = [(vn, res[f"{vn} | 30 s (H1)"][k]) for vn in M.REFS8]
            ratio = [r["D"] / thr(r["draws"], T68[vn][k.split(' | ')[1]]) for vn, r in rr]
            lines.append(f"30 s {k}: D/threshold " + ", ".join(f"{x:.2f}" for x in ratio))
        print("\n".join(lines), flush=True)
        log.outcome = "; ".join(f"{k.split(' | ')[0][:10]} {k.split(' | ')[1][:5]} {v['excluded']}/8 ({v['borderline']} bl)"
                                for k, v in summ.items())[:250]
    (OUT / "sim69.json").write_text(json.dumps({"cases": res, "summary": summ, "level": LEVEL}, default=float))
    (OUT / "sim69.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
