"""Sim 58 - Review round 24: the candidate tests of Sims 53-54 over every reference that passes the 30 s control.

Sim 57 showed that the nominal F tests are close to calibrated once the bootstrap is built on the null fit
(per-scan 2-dof thresholds p < 0.066 / 0.030 / 0.025 at 2 / 4 / 5 min). Sims 53-54 tested the candidates only
against the six references with the 2.07 dB/THz floor; Sim 55 found that floors 1.0, 1.5, 2.07 and 2.5 pass the
30 s control. Here the test of Sim 53 (candidate replaces the causal increment and the shift, F(2, N_eff - 6))
runs on the floors 1.0, 1.5 and 2.5, with the increment held above 1.5, 2 or 4 THz and a <= 1 or a + b <= 1
(18 references; with Sim 53, 24). Candidates: thicker path, ions, axoplasm, warming (two versions), the best
three node states at every extreme, and gas in a slice that swells (r 30 and 50 um).
Writes results/sim58/sim58.txt and results/sim58/sim58.json.
"""
import json
from pathlib import Path

import numpy as np

from thznerve.physics.scattering import foldy_eps
from thznerve.physics.slice_change import D_UM, eps_tissue, k_millon, mechanisms, ratio
from thznerve.runlog import RunLog
from thznerve.waveform_fit import TwoPulse, f_test, load_scans, n_eff, node_transfer

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim58"; OUT.mkdir(parents=True, exist_ok=True)
lines = []


def say(s=""):
    print(s, flush=True); lines.append(s)


tt, K, SCANS = load_scans(ROOT)
LATE = ("2 min (H2)", "4 min (H1)", "5 min (H3)")
THR = {sc: json.loads((ROOT / "results" / "sim57" / "sim57.json").read_text())["A"][sc]["thr2"] for sc in LATE}
MECH = mechanisms(k_millon())
HX = [n for n in MECH if n.startswith("heating, Liebe x")][0]
GRIDS = {"thicker path (fraction)": np.linspace(0.0, 0.2, 11),
         "ions: all tissue water (S/m)": np.linspace(0.0, 60.0, 13),
         "axoplasm conductivity (S/m)": np.linspace(0.0, 150.0, 11),
         "heating, Liebe (dT, K)": np.linspace(0.0, 30.0, 11),
         HX: np.linspace(0.0, 12.0, 13),
         "gas, slice swells, r 30 um": np.linspace(0.0, 0.2, 11),
         "gas, slice swells, r 50 um": np.linspace(0.0, 0.2, 11)}
LAMS = np.linspace(0.0, 1.0, 11)
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


out = {"thr2": THR}
with RunLog("sim58", "round 24: Sim 53 candidate tests on floors 1.0/1.5/2.5 (18 references), judged at the Sim 57 thresholds") as log:
    for k0 in (1.0, 1.5, 2.5):
        for fl in (1.5, 2.0, 4.0):
            tp = TwoPulse(tt, K, f_loss=fl)
            floor = tp.loss_transfer(k0, "zero")
            for ab in (None, 1.0):
                vn = f"floor {k0:g}, held above {fl:g} THz, " + ("a + b <= 1" if ab else "a <= 1")
                fit = lambda y, x, **kw: tp.fit(y, x, a_max=1.0, ab_max=ab, **kw)
                out[vn] = {}
                for scan in LATE:
                    y, x0 = SCANS[scan]
                    xf = tp.apply_transfer(x0, floor)
                    pr, s_ref, res = fit(y, xf, loss="minphase")
                    ne, _ = n_eff(res)
                    q0, s00, _ = fit(y, xf, s=0.0)
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
                        S_c, S_f, prev_c, prev_f = [], [], list(q0), list(pr)
                        for v in grid:
                            xt = tp.apply_transfer(xf, tr(v))
                            qc, mc, _ = min((fit(y, xt, s=0.0, p0=prev_c), fit(y, xt, s=0.0, p0=q0)), key=lambda t: t[1])
                            qf, mf, _ = min((fit(y, xt, loss="minphase", p0=prev_f), fit(y, xt, loss="minphase", p0=pr)), key=lambda t: t[1])
                            prev_c, prev_f = list(qc), list(qf)
                            S_c.append(mc); S_f.append(mf)
                        S_full = min(min(S_f), s_ref)
                        ib = int(np.argmin(S_c))
                        rs[name] = dict(p=f_test(S_c[ib], S_full, ne, 6, 2)[1], best=float(grid[ib]))
                    out[vn][scan] = dict(kappa=float(pr[4]), p_none=f_test(s00, s_ref, ne, 5, 2)[1], cands=rs)
                    say(f"{vn:40s} {scan[:5]}: " + "; ".join(
                        f"{('node ' + n[5:][:12]) if n.startswith('node ') else n[:14]} {d['p']:.3f}@{d['best']:.3g}" for n, d in rs.items()))

    (OUT / "sim58.json").write_text(json.dumps(out, indent=1, default=float))
    say()
    s53 = json.loads((ROOT / "results" / "sim53" / "sim53.json").read_text())
    s54 = json.loads((ROOT / "results" / "sim54" / "sim54.json").read_text())
    refs = {vn: out[vn] for vn in out if vn.startswith("floor")}
    for v in (v for v in s53 if v.startswith("increment")):
        refs["floor 2.07, " + v[len("increment "):]] = {sc: dict(cands={**s53[v][sc]["cands"], **{
            k: s54[v][sc]["cands"][k] for k in ("gas, slice swells, r 30 um", "gas, slice swells, r 50 um")}}) for sc in LATE}
    names = [n for n in GRIDS] + ["node states (best of the three screened)"]
    summ = {}
    for c in names:
        row, pmax = [], []
        for sc in LATE:
            if c.startswith("node states"):
                ps = [max(x["p"] for n, x in refs[v][sc]["cands"].items() if n.startswith("node ")) for v in refs]
            else:
                ps = [refs[v][sc]["cands"][c]["p"] for v in refs]
            row.append(int(sum(q >= THR[sc] for q in ps))); pmax.append(float(max(ps)))
        summ[c] = dict(passes=row, p_max=pmax, n=len(refs))
        say(f"{c:42s} references passing at 2/4/5 min (of {len(refs)}): {row}; largest p {', '.join(f'{q:.3f}' for q in pmax)}")
    out["summary"] = summ
    log.outcome = "; ".join(f"{c[:14]} {d['passes']}" for c, d in summ.items())[:250]
(OUT / "sim58.json").write_text(json.dumps(out, indent=1, default=float))
(OUT / "sim58.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
