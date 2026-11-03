"""Sim 59 - Review round 25: do the rejections of Sim 58 come from what was held fixed?

Sim 58 capped node states at the every-extreme size (lambda <= 1), gas at v <= 0.2, and fixed the timing at
s = 0. Here, on eight references (floors 1.0, 1.5, 2.07, 2.5 dB/THz; increment held above 2 THz; a <= 1 or
a + b <= 1), each candidate is tested
  (A) alone with a wider size range: node states lambda 0-10 (best three screened at lambda = 1), gas in a
      swelling slice r 30/50 um v 0-0.4, ions 0-120 S/m, axoplasm 0-300 S/m, warming (Liebe x2.31) 0-24 K,
      thicker path 0-30 %; test F(2, N_eff - 6) at the Sim 57 2-dof thresholds;
  (B) with a timing drift: the candidate plus a free shift s of both pulses, no loss term, against the full
      fit (candidate, loss increment, s free); test F(1, N_eff - 6) at the Sim 57 shift-test thresholds.
Writes results/sim59/sim59.txt and results/sim59/sim59.json.
"""
import json
from pathlib import Path

import numpy as np

from thznerve.physics.scattering import foldy_eps
from thznerve.physics.slice_change import D_UM, eps_tissue, k_millon, mechanisms, ratio
from thznerve.runlog import RunLog
from thznerve.waveform_fit import TwoPulse, f_test, load_scans, n_eff, node_transfer

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim59"; OUT.mkdir(parents=True, exist_ok=True)
lines = []


def say(s=""):
    print(s, flush=True); lines.append(s)


tt, K, SCANS = load_scans(ROOT)
LATE = ("2 min (H2)", "4 min (H1)", "5 min (H3)")
S57 = json.loads((ROOT / "results" / "sim57" / "sim57.json").read_text())["A"]
THR2 = {sc: S57[sc]["thr2"] for sc in LATE}; THR1 = {sc: S57[sc]["thr1"] for sc in LATE}
MECH = mechanisms(k_millon())
HX = [n for n in MECH if n.startswith("heating, Liebe x")][0]
GRIDS = {"thicker path (fraction)": np.linspace(0.0, 0.3, 16),
         "ions: all tissue water (S/m)": np.linspace(0.0, 120.0, 13),
         "axoplasm conductivity (S/m)": np.linspace(0.0, 300.0, 13),
         HX: np.linspace(0.0, 24.0, 13),
         "gas, slice swells, r 30 um": np.linspace(0.0, 0.4, 11),
         "gas, slice swells, r 50 um": np.linspace(0.0, 0.4, 11)}
LAMS = np.r_[np.linspace(0.0, 2.0, 9), np.linspace(2.5, 10.0, 16)]
NODES = json.loads((ROOT / "results" / "sim30" / "sim36_signedband.json").read_text())["states"]
fr = TwoPulse(tt, K).fr; band = (fr > 0) & (fr <= 4.0)
cache = {}


def slab_ratio(name, v):
    key = (name, round(float(v), 9))
    if key not in cache:
        r = np.ones(len(fr), complex)
        if v != 0:
            if name.startswith("gas"):
                rad = int(name.split("r ")[1].split()[0])
                r[band] = ratio(fr[band] * 1e12, (lambda f: (foldy_eps(f, eps_tissue(f), 1.0, rad, v), D_UM / (1 - v))))
            else:
                r[band] = ratio(fr[band] * 1e12, MECH[name][0](v))
            r[fr > 4.0] = r[band][-1]
        cache[key] = r
    return cache[key]


def best(fits):
    return min(fits, key=lambda t: t[1])


out = {"thr2": THR2, "thr1": THR1}
with RunLog("sim59", "round 25: candidates with wider sizes (nodes lambda <= 10, gas v <= 0.4) and with a free timing drift, 8 references") as log:
    tp = TwoPulse(tt, K, f_loss=2.0)
    for k0 in (1.0, 1.5, 2.07, 2.5):
        floor = tp.loss_transfer(k0, "zero")
        for ab in (None, 1.0):
            vn = f"floor {k0:g}, " + ("a + b <= 1" if ab else "a <= 1")
            fit = lambda y, x, **kw: tp.fit(y, x, a_max=1.0, ab_max=ab, **kw)
            out[vn] = {}
            for scan in LATE:
                y, x0 = SCANS[scan]
                xf = tp.apply_transfer(x0, floor)
                pr, s_ref, res = fit(y, xf, loss="minphase")
                ne, _ = n_eff(res)
                q0, _, _ = fit(y, xf, s=0.0)
                qd, _, _ = fit(y, xf)
                scr = sorted(((fit(y, tp.apply_transfer(xf, node_transfer(st, fr)), s=0.0, p0=q0)[1], st["state"]) for st in NODES))
                cands = dict(GRIDS)
                for n in [n for _, n in scr[:3]]:
                    cands["node " + n] = LAMS
                rs = {}
                for name, grid in cands.items():
                    if name.startswith("node "):
                        st = next(x for x in NODES if x["state"] == name[5:])
                        tr = lambda v, st=st: node_transfer(st, fr, lam=v)
                    else:
                        tr = lambda v, name=name: slab_ratio(name, v)
                    Sc, Sd, Sf, Ds = [], [], [], []
                    pc, pdd, pf = list(q0), list(qd), list(pr)
                    for v in grid:
                        xt = tp.apply_transfer(xf, tr(v))
                        qc, mc, _ = best((fit(y, xt, s=0.0, p0=pc), fit(y, xt, s=0.0, p0=q0)))
                        qq, md, _ = best((fit(y, xt, p0=pdd), fit(y, xt, p0=qd)))
                        qf, mf, _ = best((fit(y, xt, loss="minphase", p0=pf), fit(y, xt, loss="minphase", p0=pr)))
                        pc, pdd, pf = list(qc), list(qq), list(qf)
                        Sc.append(mc); Sd.append(md); Sf.append(mf); Ds.append(float(qq[2]))
                    S_full = min(min(Sf), s_ref)
                    ic, idd = int(np.argmin(Sc)), int(np.argmin(Sd))
                    rs[name] = dict(p_alone=f_test(Sc[ic], S_full, ne, 6, 2)[1], best_alone=float(grid[ic]),
                                    p_drift=f_test(Sd[idd], S_full, ne, 6, 1)[1], best_drift=float(grid[idd]), drift_ps=Ds[idd])
                out[vn][scan] = rs
                say(f"{vn:22s} {scan[:5]}: " + "; ".join(
                    f"{('node ' + n[5:][:10]) if n.startswith('node ') else n[:12]} alone {d['p_alone']:.3f}@{d['best_alone']:.3g}, "
                    f"drift {d['p_drift']:.3f}@{d['best_drift']:.3g},{d['drift_ps']:+.3f}ps" for n, d in rs.items()))
        (OUT / "sim59.json").write_text(json.dumps(out, indent=1, default=float))
    say()
    refs = [v for v in out if v.startswith("floor")]
    names = list(GRIDS) + ["node states (best of the three screened)"]
    summ = {}
    for c in names:
        row = {}
        for mode, thr in (("alone", THR2), ("drift", THR1)):
            cnt, pm, sz, dr = [], [], [], []
            for sc in LATE:
                if c.startswith("node states"):
                    ds = [max((x for n, x in out[v][sc].items() if n.startswith("node ")), key=lambda x: x[f"p_{mode}"]) for v in refs]
                else:
                    ds = [out[v][sc][c] for v in refs]
                ok = [d for d in ds if d[f"p_{mode}"] >= thr[sc]]
                cnt.append(len(ok)); pm.append(max(d[f"p_{mode}"] for d in ds))
                sz.append([min(d[f"best_{mode}"] for d in ok), max(d[f"best_{mode}"] for d in ok)] if ok else None)
                if mode == "drift":
                    dr.append([min(d["drift_ps"] for d in ok), max(d["drift_ps"] for d in ok)] if ok else None)
            row[mode] = dict(passes=cnt, p_max=pm, sizes=sz, drift_ps=dr or None)
            say(f"{c:40s} {mode:5s} passing at 2/4/5 min (of {len(refs)}): {cnt}; largest p {', '.join(f'{q:.3f}' for q in pm)}; "
                f"sizes {sz}" + (f"; drift {dr} ps" if mode == "drift" else ""))
        summ[c] = row
    out["summary"] = summ
    log.outcome = "; ".join(f"{c[:12]} alone {d['alone']['passes']} drift {d['drift']['passes']}" for c, d in summ.items())[:250]
(OUT / "sim59.json").write_text(json.dumps(out, indent=1, default=float))
(OUT / "sim59.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
