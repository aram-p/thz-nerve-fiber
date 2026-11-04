"""Sim 71 - Review round 32: warming with the temperature change shared freely between the two relaxations of
water, gas with a free bubble radius, the thicker path on a wider grid, and a stricter borderline rule.

Round 32 found: (1) warming models B (slow relaxation only) and C (fast only) are not a bracket on the fit: a mix
w*B + (1 - w)*C, each scaled to the measured 1 THz change, fits better and passes in one reference; (2) Sim 70
tested gas at r = 30 um only, and r = 50 um fits better; (3) per-scan thresholds carry ~11 % Monte Carlo error, so
ratios near 1 are not firm; (4) the thicker path with a drift sits at the 30 % grid edge in four references.
Here, for all twelve references (floors 0, 0.5, 1.0, 1.5, 2.07, 2.5 dB/THz; a <= 1 and a + b <= 1) at 2, 4, 5 min:
  warming  water eps = W0 + w k_B (B(20 + dT) - B(20)) + (1 - w) k_C (C(20 + dT) - C(20)), w 0-1 in 0.2 steps,
           dT 0-40 K (Sim 70 models B, C and their Millon scaling; every w matches the 1 THz change);
  gas      swelling slice, r 20, 30, 50, 75, 100 um, v 0-0.2 (the dilute range);
  thicker  path 0-51 % (floors 1.0-2.5 only; floors 0-0.5 done in Sim 70).
Each grid is screened with single-start fits; the best three are refitted multi-start and the best kept. Noise,
threshold and verdict as Sim 69/70 (noise_scale, 100 draws, level 1 - 0.05/3, any of three scans,
threshold_verdict). A verdict is called borderline if it flips in > 10 % of bootstrap resamples or its largest
ratio lies within 0.87-1.15 (the ~11 % Monte Carlo error of one threshold, round 32).
Also: the floors 0-0.5 30 s test recomputed from Sim 70 at the level 1 - 0.05/3.
Writes results/sim71/sim71.{txt,json}; resumes from results/sim71/sim71_partial.json.
"""
import json
import zlib
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from thznerve.physics.slice_change import D_UM, W0, eps_tissue, ratio
from thznerve.runlog import RunLog
from thznerve.waveform_fit import effect_noise, effect_size, multistart, noise_scale, threshold_verdict

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim71"
NB = 100
LEVEL = 1 - 0.05 / 3
LATE = ("2 min (H2)", "4 min (H1)", "5 min (H3)")


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
W_GRID = np.round(np.linspace(0, 1, 6), 2)
T_GRID = S70.WGRID
RADII = (20, 30, 50, 75, 100)
V_GRID = np.linspace(0, 0.2, 11)
THICK = np.linspace(0, 0.51, 18)


def warm_mix(w, dT):
    kB, kC = S70.K_PART["B"], S70.K_PART["C"]
    B, C = S70.liebe_part, S70.liebe_part

    def water(f):
        return (W0(f) + w * kB * (B(20.0 + dT, "B")(f) - B(20.0, "B")(f))
                + (1 - w) * kC * (C(20.0 + dT, "C")(f) - C(20.0, "C")(f)))
    return lambda f: (eps_tissue(f, water), D_UM)


def transfer(kind, a, b):
    if (kind == "warming" and b == 0) or (kind != "warming" and b == 0):
        return None
    if kind == "thicker":
        return M.transfer("thicker path (fraction)", b)
    if kind == "gas":
        return M.transfer(f"gas, slice swells, r {a} um", b)
    r = np.ones(len(M.fr), complex)
    r[M.band] = ratio(M.fr[M.band] * 1e12, warm_mix(a, b)); r[M.fr > 4.0] = r[M.band][-1]
    return r


def grid(kind):
    if kind == "warming":
        return [(float(w), float(t)) for w in W_GRID for t in T_GRID]
    if kind == "gas":
        return [(r, float(v)) for r in RADII for v in V_GRID]
    return [(0, float(v)) for v in THICK]


def run(job):
    vn, sc = job
    k0 = float(vn.split()[1].rstrip(",")); ab = 1.0 if "a + b" in vn else None
    one = lambda y, x, **kw: tp.fit(y, x, a_max=1.0, ab_max=ab, **kw)
    fit = multistart(one)
    y, x0 = SCANS[sc]
    xf = tp.apply_transfer(x0, tp.loss_transfer(k0, "zero"))
    qf, Sf, rf = fit(y, xf, loss="minphase")
    fit_f = lambda yy: fit(yy, xf, loss="minphase")
    scale = noise_scale(fit_f, y + rf, rf, 10, np.random.default_rng(zlib.crc32(f"71scale|{vn}|{sc}".encode())))
    row = {"full": dict(scale=scale)}
    kinds = ["warming", "gas"] + (["thicker"] if vn in M.REFS8 else [])
    cache = {}
    for kind in kinds:
        for g in grid(kind):
            if g not in cache.setdefault(kind, {}):
                t = transfer(kind, *g)
                cache[kind][g] = xf if t is None else tp.apply_transfer(xf, t)
        for mode in ("alone", "drift"):
            s = 0.0 if mode == "alone" else None
            scr = sorted((one(y, xx, s=s)[1], g) for g, xx in cache[kind].items())[:3]
            best = min((fit(y, cache[kind][g], s=s)[1], g) for _, g in scr)
            xt = cache[kind][best[1]]
            _, Sr, rr = fit(y, xt, s=s)
            rng = np.random.default_rng(zlib.crc32(f"71|{vn}|{sc}|{kind}|{mode}".encode()))
            dn = effect_noise(lambda yy: fit(yy, xt, s=s), fit_f, y + rr, scale * rf, NB, rng)
            row[f"{kind} | {mode}"] = dict(D=effect_size(Sr, Sf, y), size=list(best[1]), draws=[float(v) for v in dn])
    return job, row


def borderline(v):
    return v["flip"] > 0.1 or 0.87 < v["largest"] < 1.15


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    part = OUT / "sim71_partial.json"
    res = json.loads(part.read_text()) if part.exists() else {}
    jobs = [(vn, sc) for vn in REFS12 for sc in LATE if f"{vn} | {sc}" not in res]
    with RunLog("sim71", "round 32: warming with the split between slow and fast relaxation free, gas radius free (20-100 um, v <= 0.2), thicker path to 51 %; borderline if flip > 10 % or ratio 0.87-1.15") as log:
        print(f"{len(res)} done, {len(jobs)} to run", flush=True)
        with Pool(4) as pool:
            for job, row in pool.imap_unordered(run, jobs):
                res[" | ".join(job)] = row
                part.write_text(json.dumps(res, default=float))
                print(" | ".join(job), "done", flush=True)
        rng = np.random.default_rng(71)
        lines, summ = [], {}
        for refs, tag in ((S70.REFS_NEW, "floors 0-0.5"), (M.REFS8, "floors 1.0-2.5")):
            for k in [k for k in res[f"{refs[0]} | {LATE[0]}"] if k != "full"]:
                mode = k.split(" | ")[1]; per = {}
                for vn in refs:
                    rows = [res[f"{vn} | {sc}"][k] for sc in LATE]
                    v = threshold_verdict([r["D"] for r in rows], [r["draws"] for r in rows], TR[vn][mode], LEVEL, 200, rng)
                    per[vn] = dict(v, borderline=borderline(v), sizes=[r["size"] for r in rows], D=[r["D"] for r in rows])
                L = [p["largest"] for p in per.values()]
                s = summ[f"{tag} | {k}"] = dict(excluded=sum(p["excluded"] for p in per.values()),
                                                borderline=sum(p["borderline"] for p in per.values()),
                                                firm_excluded=sum(p["excluded"] and not p["borderline"] for p in per.values()),
                                                largest=[min(L), max(L)], per=per)
                lines.append(f"{tag:14s} {k:18s} excluded {s['excluded']}/{len(refs)} (firm {s['firm_excluded']}, borderline {s['borderline']}); "
                             f"largest D/thr {min(L):.2f}-{max(L):.2f}")
                for vn, p in per.items():
                    lines.append(f"      {vn:22s} " + " / ".join(f"{x:.2f}" for x in p["ratio"]) + "  D " +
                                 " / ".join(f"{x:.1f}" for x in p["D"]) + f"  flip {p['flip']:.2f}{'  BL' if p['borderline'] else ''}  sizes {p['sizes']}")
        lines.append("")
        for vn in S70.REFS_NEW:
            r = J70["cases"][f"{vn} | 30 s (H1)"]
            for mode in ("alone", "drift"):
                v = threshold_verdict([r[f"none | {mode}"]["D"]], [r[f"none | {mode}"]["draws"]], TR[vn][mode], LEVEL, 200, rng)
                summ[f"30 s | {vn} | none | {mode}"] = v
                lines.append(f"30 s {vn:22s} none {mode}: D/threshold at {LEVEL:.4f} {v['largest']:.2f}, flip {v['flip']:.2f}")
        print("\n".join(lines), flush=True)
        log.outcome = "; ".join(f"{k.split(' | ')[0][7:]} {k.split(' | ')[1]} {k.split(' | ')[2][:5]} {v['excluded']} ({v['firm_excluded']} firm)"
                                for k, v in summ.items() if not k.startswith("30 s"))[:250]
    (OUT / "sim71.json").write_text(json.dumps({"cases": res, "summary": summ, "level": LEVEL}, default=float))
    (OUT / "sim71.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
