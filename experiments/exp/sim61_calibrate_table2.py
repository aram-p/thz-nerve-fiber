"""Sim 61 - Review round 27: every Table 2 verdict that matters, calibrated on its own null with the same statistic
on the data and on the synthetic scans (thznerve.waveform_fit.null_calibrated_test).

Round 27 found that Sim 60 compared an observed p whose full fit could drop the candidate with synthetic p whose
full fit could not, and that the passing candidates were never calibrated. Here, for the eight references of
Sim 59 (floors 1.0, 1.5, 2.07, 2.5 dB/THz; increment held above 2 THz; a <= 1 or a + b <= 1) and each case:
  null        the candidate at its best size from Sim 59 for that test, no loss term; s = 0 (alone) or s free
              (with a drift);
  alternative the better of: the candidate at that size with the causal loss increment and s free, or no
              candidate with the increment and s free;
  statistic   F(n_fixed, N_eff - 6) with N_eff from the alternative's residual (n_fixed 2 alone, 1 with a drift);
  calibration up to 100 synthetic scans: null model + the alternative's residual in circular blocks of one
              correlation length; p_cal = (1 + #{p_b <= p}) / (1 + n); the candidate passes if p_cal >= 0.05.
              The first 45 cases ran all 100 draws; the next 50 drew up to 100 and the last 105 up to 60
              (argv[1]; the 2 h job limit), in batches of 20, stopping once the 99 % interval of the
              fraction lies wholly above or below 0.05. A draw with N_eff <= 6 has p = 1 (no power).
Cases: thicker path and the best node state (alone and with a drift) at 2, 4, 5 min; gas in a swelling slice
r 30/50 um with a drift; ions, axoplasm and warming (Liebe x2.31) alone at 2, 4, 5 min (4 and 5 min added in a
second run) and with a drift at 5 min, warming with a drift at 4 min. Fits start from several points (multi-start) on the data, from the data's solution on the
synthetic scans (thznerve.waveform_fit.multistart). Four worker processes; resumes from a
partial results/sim61/sim61.json. Writes results/sim61/sim61.txt and results/sim61/sim61.json.
"""
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from thznerve.physics.scattering import foldy_eps
from thznerve.physics.slice_change import D_UM, eps_tissue, k_millon, mechanisms, ratio
from thznerve.runlog import RunLog
from thznerve.waveform_fit import TwoPulse, load_scans, multistart, node_transfer, null_calibrated_test

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim61"
NB = int(sys.argv[1]) if len(sys.argv) > 1 else 100
LIMIT = int(sys.argv[2]) if len(sys.argv) > 2 else None    # cases per run, so each run fits the 2 h job limit
LATE = ("2 min (H2)", "4 min (H1)", "5 min (H3)")
tt, K, SCANS = load_scans(ROOT)
S59 = json.loads((ROOT / "results" / "sim59" / "sim59.json").read_text())
MECH = mechanisms(k_millon())
HX = [n for n in MECH if n.startswith("heating, Liebe x")][0]
NODES = {x["state"]: x for x in json.loads((ROOT / "results" / "sim30" / "sim36_signedband.json").read_text())["states"]}
tp = TwoPulse(tt, K, f_loss=2.0)
fr = tp.fr; band = (fr > 0) & (fr <= 4.0)


def transfer(name, v):
    if name.startswith("node "):
        return node_transfer(NODES[name[5:]], fr, lam=v)
    r = np.ones(len(fr), complex)
    if v != 0:
        if name.startswith("gas"):
            rad = int(name.split("r ")[1].split()[0])
            r[band] = ratio(fr[band] * 1e12, (lambda f: (foldy_eps(f, eps_tissue(f), 1.0, rad, v), D_UM / (1 - v))))
        else:
            r[band] = ratio(fr[band] * 1e12, MECH[name][0](v))
        r[fr > 4.0] = r[band][-1]
    return r


def cases():
    out = []
    for vn in (v for v in S59 if v.startswith("floor")):
        for sc in LATE:
            node = {m: max((n for n in S59[vn][sc] if n.startswith("node ")), key=lambda n: S59[vn][sc][n][f"p_{m}"])
                    for m in ("alone", "drift")}
            for m in ("alone", "drift"):
                out += [(vn, sc, "thicker path (fraction)", m), (vn, sc, node[m], m)]
            out += [(vn, sc, "gas, slice swells, r 30 um", "drift"), (vn, sc, "gas, slice swells, r 50 um", "drift")]
            for c in ("ions: all tissue water (S/m)", "axoplasm conductivity (S/m)", HX):
                out.append((vn, sc, c, "alone"))
                if sc.startswith("5 min") or (sc.startswith("4 min") and c == HX):
                    out.append((vn, sc, c, "drift"))
    return out


def run(case):
    vn, sc, c, mode = case
    k0 = float(vn.split()[1].rstrip(",")); ab = 1.0 if "a + b" in vn else None
    fit = multistart(lambda y, x, **kw: tp.fit(y, x, a_max=1.0, ab_max=ab, **kw))
    y, x0 = SCANS[sc]
    xf = tp.apply_transfer(x0, tp.loss_transfer(k0, "zero"))
    v = S59[vn][sc][c][f"best_{mode}"]
    xt = tp.apply_transfer(xf, transfer(c, v))
    d = null_calibrated_test(fit, y, xt, [xt, xf], 0.0 if mode == "alone" else None, 2 if mode == "alone" else 1,
                             nb=NB, batch=20, rng=np.random.default_rng(61), model=lambda q, x, s: tp.model(q, x, s=s))
    return case, dict(size=v, **d)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    lines, res = [], {}
    part = OUT / "sim61.json"
    if part.exists():                                   # resume: keep the cases already done
        old = json.loads(part.read_text())
        res = old.get("cases", old)
    with RunLog("sim61", f"round 27: Table 2 verdicts calibrated on their own null, matched statistic, {NB} block-bootstrap scans, 8 references, 4 workers") as log:
        cs = [c for c in cases() if " | ".join(c) not in res][:LIMIT]
        print(f"{len(res)} cases done, {len(cs)} to run", flush=True)
        with Pool(4) as pool:
            for (vn, sc, c, mode), d in pool.imap_unordered(run, cs):
                key = f"{vn} | {sc} | {c} | {mode}"
                res[key] = d
                s = (f"{vn:22s} {sc[:5]} {c[:30]:30s} {mode:5s} size {d['size']:g}: nominal p {d['p_nominal']:.4f}, calibrated "
                     f"p {d['p_cal']:.3f} ({'pass' if d['p_cal'] >= 0.05 else 'fail'}); null 5 % threshold {d['thr']:.3f}, "
                     f"null rejects {d['rej_nominal']:.2f} at nominal 0.05")
                print(s, flush=True); lines.append(s)
                (OUT / "sim61.json").write_text(json.dumps(res, indent=1, default=float))
        lines = [f"{k}: size {d['size']:g}, nominal p {d['p_nominal']:.4f}, calibrated p {d['p_cal']:.3f} "
                 f"({'pass' if d['p_cal'] >= 0.05 else 'fail'}); null 5 % threshold {d['thr']:.3f}, null rejects "
                 f"{d['rej_nominal']:.2f} at nominal 0.05" for k, d in sorted(res.items())] + [""]
        summ = {}
        for key, d in res.items():
            vn, sc, c, mode = key.split(" | ")
            g = "node states (best screened)" if c.startswith("node ") else c
            row = summ.setdefault(f"{g} | {mode}", {s: [0, 0, 1.0, 0.0] for s in LATE})
            r = row[sc]; r[1] += 1; r[0] += d["p_cal"] >= 0.05
            r[2] = min(r[2], d["p_cal"]); r[3] = max(r[3], d["p_cal"])
        for k, row in summ.items():
            s = f"{k:48s} " + "; ".join(f"{sc[:5]} {r[0]}/{r[1]} pass (calibrated p {r[2]:.3f}-{r[3]:.3f})" for sc, r in row.items() if r[1])
            lines.append(s); print(s, flush=True)
        log.outcome = "; ".join(f"{k[:20]} " + "/".join(str(r[0]) for r in row.values() if r[1]) for k, row in summ.items())[:250]
    (OUT / "sim61.json").write_text(json.dumps({"cases": res, "summary": summ}, indent=1, default=float))
    (OUT / "sim61.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
