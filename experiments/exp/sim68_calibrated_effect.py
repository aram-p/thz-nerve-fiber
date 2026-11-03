"""Sim 68 - Review round 29: the effect size D judged against what fitting noise gives in a biased trace.

Round 29 found that the Sim 66 yardstick (traced unbiased pairs, N_eff 100-230) is on the wrong scale for
the biased traces (N_eff 13-20, full fits leaving 6-15 %): a true cause would often leave D above it. Here:
  noise     for every restricted fit (no change; each candidate at its Sim 66 best size, or for node states the
            best of all 80 states over lambda 0.5-6), synthetic traces = that fit's curve + a phase-randomised
            copy of the full fit's residual (same autocorrelation), refitted by the restricted fit (size held:
            this favours the candidate) and by the full fit (causal loss increment, s free): NB draws of D
            (thznerve.waveform_fit.effect_noise).
  tracing   the six ordered pairs of traced unbiased pulses, built like each reference (floor applied, a <= 1,
            a + b <= 1 where the reference has it): largest D per reference and mode.
  verdict   the candidate is excluded in a case if D > sqrt(q95(noise D)^2 + tracing D^2), i.e. beyond what
            noise and tracing differences together give.
Also recorded: absolute residual (% of trace rms) of the full fit, of each restricted fit and of the tracing
pairs, and the lead of the second pulse over the empty-cuvette pulse in each best fit.
References: the eight of Sim 59; scans 30 s (no change only), 2, 4, 5 min. Resumes from
results/sim68/sim68_partial.json. Writes results/sim68/sim68.txt and results/sim68/sim68.json.
"""
import json
import zlib
from itertools import permutations
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from thznerve.runlog import RunLog
from thznerve.waveform_fit import effect_noise, effect_size, multistart

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim68"
NB = 30
LAMS = [0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6]


def s64():
    import importlib.util
    spec = importlib.util.spec_from_file_location("s64", ROOT / "experiments" / "exp" / "sim64_size_fitted_tests.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


M = s64()
tp, SCANS = M.tp, M.SCANS
R66 = json.loads((ROOT / "results" / "sim66" / "sim66.json").read_text())["cases"]
CANDS = {"thicker path (fraction)": "plausible", "gas, slice swells, r 30 um": "wide",
         "ions: all tissue water (S/m)": "wide", "axoplasm conductivity (S/m)": "wide", M.HX: "plausible"}


def setup(vn):
    k0 = float(vn.split()[1].rstrip(",")); ab = 1.0 if "a + b" in vn else None
    one = lambda y, x, **kw: tp.fit(y, x, a_max=1.0, ab_max=ab, **kw)
    return k0, one, multistart(one)


def rms_pct(S, y):
    return 100 * np.sqrt(S / len(y)) / np.sqrt(np.mean(y ** 2))


def lead(q, mode):
    return float(-q[2] if mode == "alone" else -(q[2] + q[3]))


def yardstick(vn):
    k0, _, fit = setup(vn)
    ref = {"H1": SCANS["30 s (H1)"][1], "H2": SCANS["2 min (H2)"][1], "H3": SCANS["5 min (H3)"][1]}
    out = {}
    for a, b in permutations(ref, 2):
        y = ref[a]; x = tp.apply_transfer(ref[b], tp.loss_transfer(k0, "zero"))
        Sf = fit(y, x, loss="minphase")[1]
        out[f"{a} from {b}"] = dict(alone=effect_size(fit(y, x, s=0.0)[1], Sf, y), drift=effect_size(fit(y, x)[1], Sf, y),
                                    resid_full=rms_pct(Sf, y))
    return out


def run(job):
    if job[0] == "yard":
        return job, yardstick(job[1])
    vn, sc = job
    k0, one, fit = setup(vn)
    y, x0 = SCANS[sc]
    xf = tp.apply_transfer(x0, tp.loss_transfer(k0, "zero"))
    qf, Sf, rf = fit(y, xf, loss="minphase")
    fit_f = lambda yy: fit(yy, xf, loss="minphase")
    row = {"full": dict(resid=rms_pct(Sf, y), lead=float(-(qf[2] + qf[3])))}
    models = [("none", None, 0.0)]
    if not sc.startswith("30 s"):
        models += [(c, c, None) for c in CANDS] + [("node", "node", None)]
    for name, c, _ in models:
        for mode in ("alone", "drift"):
            s = 0.0 if mode == "alone" else None
            if c is None:
                xt, size, tag = xf, 0.0, "-"
            elif c == "node":
                scr = sorted((one(y, tp.apply_transfer(xf, M.transfer("node " + st, lam)), s=s)[1], st, lam)
                             for st in M.NODES for lam in LAMS)[:3]
                best = min((fit(y, tp.apply_transfer(xf, M.transfer("node " + st, lam)), s=s)[1], st, lam) for _, st, lam in scr)
                _, tag, size = best
                xt = tp.apply_transfer(xf, M.transfer("node " + tag, size))
            else:
                size = R66[f"{vn} | {sc}"][f"{c} | {CANDS[c]}"][mode]["size"]; tag = c
                if size == 0:
                    row[f"{name} | {mode}"] = dict(row[f"none | {mode}"], state=tag)
                    continue
                xt = tp.apply_transfer(xf, M.transfer(c, size))
            q, Sr, rr = fit(y, xt, s=s)
            D = effect_size(Sr, Sf, y)
            rng = np.random.default_rng(zlib.crc32(f"{vn}|{sc}|{name}|{mode}".encode()))
            dn = effect_noise(lambda yy: fit(yy, xt, s=s), fit_f, y + rr, rf, NB, rng)
            row[f"{name} | {mode}"] = dict(D=D, size=float(size), state=tag, resid=rms_pct(Sr, y), lead=lead(q, mode),
                                          noise_q50=float(np.median(dn)), noise_q95=float(np.quantile(dn, 0.95)),
                                          noise_max=float(dn.max()))
    return job, row


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    part = OUT / "sim68_partial.json"
    res = json.loads(part.read_text()) if part.exists() else {}
    jobs = [("yard", vn) for vn in M.REFS8] + [(vn, sc) for sc in ("5 min (H3)", "4 min (H1)", "2 min (H2)", "30 s (H1)")
                                                for vn in M.REFS8]
    jobs = [j for j in jobs if " | ".join(j) not in res]
    with RunLog("sim68", f"round 29: effect size against fitting noise ({NB} phase-randomised draws per fit) and per-reference tracing pairs; all 80 node states") as log:
        print(f"{len(res)} done, {len(jobs)} to run", flush=True)
        with Pool(4) as pool:
            for job, row in pool.imap_unordered(run, jobs):
                res[" | ".join(job)] = row
                part.write_text(json.dumps(res, indent=1, default=float))
                print(" | ".join(job), "done", flush=True)
        T = {vn: {m: max(d[m] for d in res[f"yard | {vn}"].values()) for m in ("alone", "drift")} for vn in M.REFS8}
        lines = ["Tracing pairs per reference (largest D, alone / drift; full-fit residual of the pairs):"]
        for vn in M.REFS8:
            rs = [d["resid_full"] for d in res[f"yard | {vn}"].values()]
            lines.append(f"  {vn:22s} {T[vn]['alone']:.2f} / {T[vn]['drift']:.2f}; residual {min(rs):.1f}-{max(rs):.1f} %")
        lines.append("")
        summ = {}
        for sc in SCANS:
            rows = {vn: res[f"{vn} | {sc}"] for vn in M.REFS8}
            fr = [r["full"]["resid"] for r in rows.values()]
            lines.append(f"{sc}: full-fit residual {min(fr):.1f}-{max(fr):.1f} %, lead {min(r['full']['lead'] for r in rows.values()):+.3f} "
                         f"to {max(r['full']['lead'] for r in rows.values()):+.3f} ps")
            for k in rows[M.REFS8[0]]:
                if k == "full":
                    continue
                mode = k.split(" | ")[1]
                thr = {vn: float(np.hypot(r[k]["noise_q95"], T[vn][mode])) for vn, r in rows.items()}
                ex = [vn for vn, r in rows.items() if r[k]["D"] > thr[vn]]
                ds = [r[k]["D"] for r in rows.values()]; th = list(thr.values())
                summ[f"{k} | {sc}"] = dict(excluded=len(ex), of=8, D=[min(ds), max(ds)], threshold=[min(th), max(th)],
                                           not_excluded=[vn for vn in M.REFS8 if vn not in ex])
                lines.append(f"  {k:40s} excluded in {len(ex)}/8; D {min(ds):5.2f}-{max(ds):5.2f} %, threshold {min(th):.2f}-{max(th):.2f} %; "
                             f"residual {min(r[k]['resid'] for r in rows.values()):.1f}-{max(r[k]['resid'] for r in rows.values()):.1f} %; "
                             f"lead {min(r[k]['lead'] for r in rows.values()):+.3f} to {max(r[k]['lead'] for r in rows.values()):+.3f} ps; "
                             f"sizes {sorted(set(round(r[k]['size'], 3) for r in rows.values()))}")
        print("\n".join(lines), flush=True)
        log.outcome = "; ".join(f"{k.split(' | ')[0][:8]} {k.split(' | ')[1][:5]} {k.split(' | ')[2][:5]} {v['excluded']}/8"
                                for k, v in summ.items())[:250]
    (OUT / "sim68.json").write_text(json.dumps({"cases": res, "tracing": T, "summary": summ}, indent=1, default=float))
    (OUT / "sim68.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
