"""Sim 70 - Review round 31: the Table 2 verdicts at loss floors 0 and 0.5 dB/THz, gas limited to the dilute
range, and warming with three temperature models.

Round 31 found: (1) the floors 1.0-2.5 dB/THz were chosen so that the 30 s scan shows no change, and at floor 0
the thicker path and gas with a drift are excluded; (2) gas was searched up to v = 0.4, beyond the dilute
coherent-medium model (v <= 0.2); (3) the warming verdict rests on Liebe's temperature dependence above 1 THz.
Here, with the Sim 69 rule (noise scaled by the full fit's absorption, `noise_scale`; 100 draws; per-scan level
1 - 0.05/3; any of 2, 4, 5 min; bootstrap flip > 10 % = borderline; `threshold_verdict`):
  new references  floor 0 and 0.5 dB/THz, a <= 1 and a + b <= 1; tracing pairs built like each (as Sim 68);
                  candidates: no change (also at 30 s), thicker path 0-51 %, gas r 30 um v 0-0.2, ions 0-60 S/m,
                  axoplasm 0-150 S/m, and the three warming models below; node states are not refitted (they are
                  excluded by what an aqueous gap can hold, Sec. 4.6).
  the eight references of Table 2: gas v <= 0.2 and warming models B and C.
  warming models  each scaled to the measured 4 -> 21 C absorption change at 1 THz (Millon et al. 2022):
                  A Liebe (1991) as in Table 2; B only the slow relaxation warms (eps_s, tau1; eps_1, eps_inf,
                  tau2 held at 20 C); C only the fast relaxation warms (tau2 as Liebe's, the rest held at 20 C).
                  B and C bracket how the temperature change is shared between the two relaxations above 1 THz.
  node numbers    for every node state at every extreme: loss slope (dB/THz, 0.2-1.7 THz) and delay at 0.6 THz.
Writes results/sim70/sim70.{txt,json}; resumes from results/sim70/sim70_partial.json.
"""
import json
import zlib
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from thznerve.physics.materials import DoubleDebye, alpha_per_cm, water_liebe
from thznerve.physics.slice_change import D_UM, W0, eps_tissue
from thznerve.runlog import RunLog
from thznerve.waveform_fit import effect_noise, effect_size, multistart, noise_scale, threshold_verdict

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim70"
NB = 100
LEVEL = 1 - 0.05 / 3
LATE = ("2 min (H2)", "4 min (H1)", "5 min (H3)")


def load(name, file):
    import importlib.util
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / "exp" / file)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


M = load("s64", "sim64_size_fitted_tests.py")
S68 = load("s68", "sim68_calibrated_effect.py")
tp, SCANS = M.tp, M.SCANS
REFS_NEW = [f"floor {k:g}, " + ab for k in (0.0, 0.5) for ab in ("a <= 1", "a + b <= 1")]
T68 = json.loads((ROOT / "results" / "sim68" / "sim68.json").read_text())["tracing"]


def liebe_part(T, part):
    full, ref = water_liebe(T), water_liebe(20.0)
    if part == "A":
        return full
    if part == "B":
        return DoubleDebye(full.eps_s, ref.eps_1, ref.eps_inf, full.tau1_s, ref.tau2_s)
    return DoubleDebye(ref.eps_s, ref.eps_1, ref.eps_inf, ref.tau1_s, full.tau2_s)


def k_part(part):
    d = alpha_per_cm(liebe_part(21.0, part)(1e12), 1e12) - alpha_per_cm(liebe_part(4.0, part)(1e12), 1e12)
    return 46.0 / d


K_PART = {p: float(k_part(p)) for p in "ABC"}
WARM = {p: f"warming {p}" for p in "ABC"}


def warm_fn(part):
    k = K_PART[part]
    return lambda v: (lambda f: (eps_tissue(f, lambda ff: W0(ff) + k * (liebe_part(20.0 + v, part)(ff) - liebe_part(20.0, part)(ff))), D_UM))


def transfer(name, v):
    if name.startswith("warming"):
        r = np.ones(len(M.fr), complex)
        if v != 0:
            from thznerve.physics.slice_change import ratio
            r[M.band] = ratio(M.fr[M.band] * 1e12, warm_fn(name[-1])(v)); r[M.fr > 4.0] = r[M.band][-1]
        return r
    return M.transfer(name, v)


GAS = "gas, slice swells, r 30 um"
WGRID = np.array([0, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 40.0])
GRIDS_NEW = {
    "thicker path (fraction)": np.linspace(0, 0.51, 18),
    GAS: np.linspace(0, 0.2, 11),
    "ions: all tissue water (S/m)": np.array([0, 1, 2, 3, 4, 5, 7.5, 10, 15, 20, 30, 40, 60.0]),
    "axoplasm conductivity (S/m)": np.array([0, 2.5, 5, 10, 20, 30, 50, 75, 100, 150.0]),
    WARM["A"]: WGRID, WARM["B"]: WGRID, WARM["C"]: WGRID,
}
GRIDS_8 = {GAS: GRIDS_NEW[GAS], WARM["B"]: WGRID, WARM["C"]: WGRID}


def run(job):
    vn, sc = job
    if vn.startswith("yard"):
        return job, S68.yardstick(sc)
    k0 = float(vn.split()[1].rstrip(",")); ab = 1.0 if "a + b" in vn else None
    fit = multistart(lambda y, x, **kw: tp.fit(y, x, a_max=1.0, ab_max=ab, **kw))
    y, x0 = SCANS[sc]
    xf = tp.apply_transfer(x0, tp.loss_transfer(k0, "zero"))
    qf, Sf, rf = fit(y, xf, loss="minphase")
    fit_f = lambda yy: fit(yy, xf, loss="minphase")
    scale = noise_scale(fit_f, y + rf, rf, 10, np.random.default_rng(zlib.crc32(f"70scale|{vn}|{sc}".encode())))
    row = {"full": dict(scale=scale, lead=float(-(qf[2] + qf[3])))}
    if vn in REFS_NEW:
        grids = {} if sc.startswith("30 s") else GRIDS_NEW
        cands = [("none", None)] + list(grids)
    else:
        grids, cands = GRIDS_8, list(GRIDS_8)
    for c in cands:
        c = c[0] if isinstance(c, tuple) else c
        for mode in ("alone", "drift"):
            s = 0.0 if mode == "alone" else None
            if c == "none":
                size, xt = 0.0, xf
            else:
                best = None
                for v in grids[c]:
                    xx = xf if v == 0 else tp.apply_transfer(xf, transfer(c, v))
                    S = fit(y, xx, s=s)[1]
                    if best is None or S < best[0]:
                        best = (S, float(v), xx)
                _, size, xt = best
            _, Sr, rr = fit(y, xt, s=s)
            rng = np.random.default_rng(zlib.crc32(f"70|{vn}|{sc}|{c}|{mode}".encode()))
            dn = effect_noise(lambda yy: fit(yy, xt, s=s), fit_f, y + rr, scale * rf, NB, rng)
            row[f"{c} | {mode}"] = dict(D=effect_size(Sr, Sf, y), size=size, draws=[float(v) for v in dn])
    return job, row


def node_numbers():
    out = []
    sel = (M.fr >= 0.2) & (M.fr <= 1.7)
    for name in M.NODES:
        t = M.transfer("node " + name, 1.0)
        loss = -20 * np.log10(np.abs(t[sel]))
        out.append(dict(state=name, slope=float(np.polyfit(M.fr[sel], loss, 1)[0]),
                        delay=float(-np.interp(0.6, M.fr, np.unwrap(np.angle(t))))))
    return out


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    part = OUT / "sim70_partial.json"
    res = json.loads(part.read_text()) if part.exists() else {}
    jobs = [(f"yard {vn}", vn) for vn in REFS_NEW] + [(vn, sc) for vn in REFS_NEW for sc in LATE + ("30 s (H1)",)] + \
           [(vn, sc) for vn in M.REFS8 for sc in LATE]
    jobs = [j for j in jobs if " | ".join(j) not in res]
    with RunLog("sim70", f"round 31: floors 0 and 0.5 dB/THz, gas v <= 0.2, warming models A/B/C (k {K_PART['A']:.2f}/{K_PART['B']:.2f}/{K_PART['C']:.2f}); Sim 69 rule") as log:
        print(f"k: {K_PART}; {len(res)} done, {len(jobs)} to run", flush=True)
        with Pool(4) as pool:
            for job, row in pool.imap_unordered(run, jobs):
                res[" | ".join(job)] = row
                part.write_text(json.dumps(res, default=float))
                print(" | ".join(job), "done", flush=True)
        T = dict(T68)
        for vn in REFS_NEW:
            T[vn] = {m: max(d[m] for d in res[f"yard {vn} | {vn}"].values()) for m in ("alone", "drift")}
        lines = [f"Warming models scaled to Millon at 1 THz: k A {K_PART['A']:.3f}, B {K_PART['B']:.3f}, C {K_PART['C']:.3f}",
                 "Tracing D (alone / drift): " + "; ".join(f"{vn} {T[vn]['alone']:.2f} / {T[vn]['drift']:.2f}" for vn in REFS_NEW), ""]
        summ = {}
        rng = np.random.default_rng(70)
        for refs, tag in ((REFS_NEW, "floors 0-0.5"), (M.REFS8, "floors 1.0-2.5")):
            keys = [k for k in res[f"{refs[0]} | {LATE[0]}"] if k != "full"]
            for k in keys:
                mode = k.split(" | ")[1]
                per = {}
                for vn in refs:
                    rows = [res[f"{vn} | {sc}"][k] for sc in LATE]
                    v = threshold_verdict([r["D"] for r in rows], [r["draws"] for r in rows], T[vn][mode], LEVEL, 200, rng)
                    per[vn] = dict(v, sizes=[r["size"] for r in rows], D=[r["D"] for r in rows])
                L = [p["largest"] for p in per.values()]
                summ[f"{tag} | {k}"] = dict(excluded=sum(p["excluded"] for p in per.values()),
                                            borderline=sum(p["flip"] > 0.1 for p in per.values()), largest=[min(L), max(L)], per=per)
                s = summ[f"{tag} | {k}"]
                lines.append(f"{tag:14s} {k:36s} excluded {s['excluded']}/{len(refs)} ({s['borderline']} bl); largest D/thr {min(L):.2f}-{max(L):.2f}")
                for vn, p in per.items():
                    lines.append(f"      {vn:22s} " + " / ".join(f"{x:.2f}" for x in p["ratio"]) + f"  D " +
                                 " / ".join(f"{x:.1f}" for x in p["D"]) + f"  flip {p['flip']:.2f}  sizes {p['sizes']}")
        lines.append("")
        for vn in REFS_NEW:
            r = res[f"{vn} | 30 s (H1)"]
            for mode in ("alone", "drift"):
                v = threshold_verdict([r[f"none | {mode}"]["D"]], [r[f"none | {mode}"]["draws"]], T[vn][mode], 0.95, 200, rng)
                summ[f"30 s | {vn} | none | {mode}"] = v
                lines.append(f"30 s {vn:22s} none {mode}: D {r[f'none | {mode}']['D']:.2f} %, D/threshold (single scan, 0.95) {v['largest']:.2f}, flip {v['flip']:.2f}")
        nn = node_numbers()
        md = max(nn, key=lambda x: x["delay"]); ms = max(nn, key=lambda x: x["slope"])
        lines += ["", f"Node states at every extreme ({len(nn)}): largest delay at 0.6 THz {md['delay']:.4f} rad ({md['state']}, slope {md['slope']:.2f} dB/THz); "
                  f"steepest slope {ms['slope']:.2f} dB/THz ({ms['state']}, delay {ms['delay']:+.3f} rad)",
                  "  slope >= 4 dB/THz: " + ", ".join(f"{x['state']} {x['slope']:.2f} ({x['delay']:+.3f} rad)" for x in nn if x["slope"] >= 4)]
        print("\n".join(lines), flush=True)
        log.outcome = "; ".join(f"{k.split(' | ')[0][7:]} {k.split(' | ')[1][:9]} {k.split(' | ')[2][:5]} {v['excluded']}" for k, v in summ.items()
                                if not k.startswith("30 s"))[:250]
    (OUT / "sim70.json").write_text(json.dumps({"cases": res, "summary": summ, "tracing": T, "k": K_PART, "nodes": nn, "level": LEVEL}, default=float))
    (OUT / "sim70.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
