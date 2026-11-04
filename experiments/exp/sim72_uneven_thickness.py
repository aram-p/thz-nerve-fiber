"""Sim 72 - Review round 33: a slice whose thickness varies across the beam.

Round 33 found that every slab candidate so far had one uniform thickness, while the beam (~5 mm) averages
the transmission over regions of different thickness. The detector sees the coherent mean, so a spread of
thickness loses amplitude as exp(-(phase spread)^2 / 2): a loss rising with frequency plus a lag, the shape
the paper said no bulk cause gives. Here the relative thickness change v varies across the beam as a Gaussian
of mean m (-0.3 to 0.51) and rms s (0 to 0.4), cut at a thickness of 5 % (thznerve.physics.slice_change.uneven,
over the uniform-thickness transfers of Sim 64), for all twelve references (floors 0, 0.5, 1.0, 1.5, 2.07,
2.5 dB/THz; a <= 1 and a + b <= 1) at 2, 4, 5 min, and at 30 s for floors 0-0.5. Alone and with a drift.
Screen with single-start fits, refit the best three multi-start; noise, threshold and verdict as Sims 69-71
(noise_scale, 100 draws, level 1 - 0.05/3, any of three scans, threshold_verdict; borderline if the verdict
flips in > 10 % of bootstrap resamples or the largest ratio lies within 0.87-1.15). Best sizes on a grid edge
are flagged. Writes results/sim72/sim72.{txt,json}; resumes from results/sim72/sim72_partial.json.
"""
import json
import zlib
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from thznerve.physics.slice_change import uneven
from thznerve.runlog import RunLog
from thznerve.waveform_fit import effect_noise, effect_size, multistart, noise_scale, threshold_verdict

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim72"
NB = 100
LEVEL = 1 - 0.05 / 3
LATE = ("2 min (H2)", "4 min (H1)", "5 min (H3)")
EARLY = "30 s (H1)"
MS = np.round(np.arange(-0.30, 0.511, 0.03), 3)
SS = np.array([0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.25, 0.3, 0.35, 0.4])
VG = np.round(np.arange(-0.95, 1.801, 0.01), 3)


def load(name, file):
    import importlib.util
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / "exp" / file)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


S70 = load("s70", "sim70_floors_and_models.py")
M, tp, SCANS = S70.M, S70.tp, S70.SCANS
REFS12 = S70.REFS_NEW + M.REFS8
J70 = json.loads((ROOT / "results" / "sim70" / "sim70.json").read_text())
TR = J70["tracing"]
_T = []


def table():
    if not _T:
        _T.append(np.array([M.transfer("thicker path (fraction)", float(v)) for v in VG]))
    return _T[0]


def edge(m, s):
    return bool(m <= MS[0] or m >= MS[-1] or s >= SS[-1])


def run(job):
    vn, sc = job
    T = table()
    k0 = float(vn.split()[1].rstrip(",")); ab = 1.0 if "a + b" in vn else None
    one = lambda y, x, **kw: tp.fit(y, x, a_max=1.0, ab_max=ab, **kw)
    fit = multistart(one)
    y, x0 = SCANS[sc]
    xf = tp.apply_transfer(x0, tp.loss_transfer(k0, "zero"))
    qf, Sf, rf = fit(y, xf, loss="minphase")
    fit_f = lambda yy: fit(yy, xf, loss="minphase")
    scale = noise_scale(fit_f, y + rf, rf, 10, np.random.default_rng(zlib.crc32(f"72scale|{vn}|{sc}".encode())))
    cache = {(float(m), float(s)): tp.apply_transfer(xf, uneven(T, VG, m, s)) for m in MS for s in SS}
    row = {"full": dict(scale=scale)}
    for mode in ("alone", "drift"):
        sv = 0.0 if mode == "alone" else None
        scr = sorted((one(y, xx, s=sv)[1], g) for g, xx in cache.items())[:3]
        best = min((fit(y, cache[g], s=sv)[1], g) for _, g in scr)
        xt = cache[best[1]]
        q, Sr, rr = fit(y, xt, s=sv)
        rng = np.random.default_rng(zlib.crc32(f"72|{vn}|{sc}|{mode}".encode()))
        dn = effect_noise(lambda yy: fit(yy, xt, s=sv), fit_f, y + rr, scale * rf, NB, rng)
        row[f"uneven | {mode}"] = dict(D=effect_size(Sr, Sf, y), size=list(best[1]), edge=edge(*best[1]),
                                       q=[float(v) for v in q], draws=[float(v) for v in dn])
    return job, row


def borderline(v):
    return v["flip"] > 0.1 or 0.87 < v["largest"] < 1.15


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    part = OUT / "sim72_partial.json"
    res = json.loads(part.read_text()) if part.exists() else {}
    jobs = [(vn, sc) for vn in S70.REFS_NEW for sc in LATE + (EARLY,)] + [(vn, sc) for vn in M.REFS8 for sc in LATE]
    jobs = [j for j in jobs if " | ".join(j) not in res]
    with RunLog("sim72", "round 33: slice thickness varying across the beam (Gaussian, mean -0.3-0.51, rms 0-0.4), all twelve references, 2/4/5 min and 30 s at floors 0-0.5; Sim 71 rule") as log:
        print(f"{len(res)} done, {len(jobs)} to run", flush=True)
        with Pool(4) as pool:
            for job, row in pool.imap_unordered(run, jobs):
                res[" | ".join(job)] = row
                part.write_text(json.dumps(res, default=float))
                print(" | ".join(job), "done", {k: (round(v["D"], 2), v["size"]) for k, v in row.items() if k != "full"}, flush=True)
        rng = np.random.default_rng(72)
        lines, summ = [], {}
        for refs, tag in ((S70.REFS_NEW, "floors 0-0.5"), (M.REFS8, "floors 1.0-2.5")):
            for mode in ("alone", "drift"):
                k = f"uneven | {mode}"; per = {}
                for vn in refs:
                    rows = [res[f"{vn} | {sc}"][k] for sc in LATE]
                    v = threshold_verdict([r["D"] for r in rows], [r["draws"] for r in rows], TR[vn][mode], LEVEL, 200, rng)
                    per[vn] = dict(v, borderline=borderline(v), sizes=[r["size"] for r in rows], edge=[r["edge"] for r in rows],
                                   D=[r["D"] for r in rows])
                L = [p["largest"] for p in per.values()]
                s = summ[f"{tag} | {k}"] = dict(excluded=sum(p["excluded"] for p in per.values()),
                                                borderline=sum(p["borderline"] for p in per.values()),
                                                firm_excluded=sum(p["excluded"] and not p["borderline"] for p in per.values()),
                                                largest=[min(L), max(L)], per=per)
                lines.append(f"{tag:14s} {k:15s} excluded {s['excluded']}/{len(refs)} (firm {s['firm_excluded']}, borderline {s['borderline']}); "
                             f"largest D/thr {min(L):.2f}-{max(L):.2f}")
                for vn, p in per.items():
                    lines.append(f"      {vn:22s} " + " / ".join(f"{x:.2f}" for x in p["ratio"]) + "  D " +
                                 " / ".join(f"{x:.1f}" for x in p["D"]) + f"  flip {p['flip']:.2f}{'  BL' if p['borderline'] else ''}"
                                 f"  sizes (m, s) {p['sizes']}{'  GRID EDGE' if any(p['edge']) else ''}")
        lines.append("")
        for vn in S70.REFS_NEW:
            r0 = J70["cases"][f"{vn} | {EARLY}"]; r = res[f"{vn} | {EARLY}"]
            for mode in ("alone", "drift"):
                v = threshold_verdict([r[f"uneven | {mode}"]["D"]], [r[f"uneven | {mode}"]["draws"]], TR[vn][mode], LEVEL, 200, rng)
                summ[f"30 s | {vn} | uneven | {mode}"] = dict(v, size=r[f"uneven | {mode}"]["size"], D=r[f"uneven | {mode}"]["D"])
                lines.append(f"30 s {vn:22s} uneven {mode}: D {r[f'uneven | {mode}']['D']:.2f} % (no change {r0[f'none | {mode}']['D']:.2f} %), "
                             f"D/threshold {v['largest']:.2f}, flip {v['flip']:.2f}, size (m, s) {r[f'uneven | {mode}']['size']}")
        print("\n".join(lines), flush=True)
        log.outcome = "; ".join(f"{k.split(' | ')[0][7:]} {k.split(' | ')[2]} {v['excluded']} ({v['firm_excluded']} firm)"
                                for k, v in summ.items() if not k.startswith("30 s"))[:250]
    (OUT / "sim72.json").write_text(json.dumps({"cases": res, "summary": summ, "level": LEVEL}, default=float))
    (OUT / "sim72.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
