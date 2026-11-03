"""Sim 60 - Review round 26, items 219 and 222.

(A) The drift test of Sim 59 (candidate + free shift s, no loss increment, against the full fit with the
increment; F(1, N_eff - 6)) was judged at the shift-test thresholds of Sim 57, calibrated for another null.
Here it is calibrated on its own null: for each of the eight Sim 59 references, at 5 min (ions, axoplasm,
warming) and 4 min (warming), the candidate is fixed at its best drift size, synthetic scans are the
candidate + drift fit plus the full fit's residual resampled in circular blocks of one correlation length
(200 each), and the observed p is converted to a calibrated p (fraction of null p at or below it).
(B) Where the second pulse arrives: lead over the empty-cuvette pulse, -(s + tau) ps, in the reference fit
and in each candidate's best passing fit (alone and with a drift) at 2, 4 and 5 min, against the 0.015 ps to
which the no-voltage minima of the three reports agree.
Writes results/sim60/sim60.txt and results/sim60/sim60.json.
"""
import json
from pathlib import Path

import numpy as np

from thznerve.physics.scattering import foldy_eps
from thznerve.physics.slice_change import D_UM, eps_tissue, k_millon, mechanisms, ratio
from thznerve.runlog import RunLog
from thznerve.waveform_fit import TwoPulse, f_test, load_scans, n_eff, node_transfer

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim60"; OUT.mkdir(parents=True, exist_ok=True)
import sys
NB = int(sys.argv[1]) if len(sys.argv) > 1 else 200
lines = []


def say(s=""):
    print(s, flush=True); lines.append(s)


tt, K, SCANS = load_scans(ROOT)
LATE = ("2 min (H2)", "4 min (H1)", "5 min (H3)")
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


def blocks(rng, res, blk):
    n = len(res)
    starts = rng.integers(0, n, size=n // blk + 1)
    return np.concatenate([np.take(res, np.arange(s0, s0 + blk), mode="wrap") for s0 in starts])[:n]


out = {"A": {}, "B": {}}
with RunLog("sim60", "round 26: drift test calibrated on its own null (ions, axoplasm, warming; 8 references); second-pulse lead per candidate") as log:
    refs = [v for v in S59 if v.startswith("floor")]
    cases = [(c, "5 min (H3)") for c in ("ions: all tissue water (S/m)", "axoplasm conductivity (S/m)", HX)] + [(HX, "4 min (H1)")]
    for vn in refs:
        k0 = float(vn.split()[1].rstrip(",")); ab = 1.0 if "a + b" in vn else None
        floor = tp.loss_transfer(k0, "zero")
        fit = lambda y, x, **kw: tp.fit(y, x, a_max=1.0, ab_max=ab, **kw)
        for c, sc in cases:
            y, x0 = SCANS[sc]
            xf = tp.apply_transfer(x0, floor)
            pr, s_ref, res = fit(y, xf, loss="minphase")
            ne, corr = n_eff(res)
            v = S59[vn][sc][c]["best_drift"]
            xt = tp.apply_transfer(xf, transfer(c, v))
            qd, md, _ = min((fit(y, xt), fit(y, xt, p0=[pr[0], pr[1], pr[2], pr[3]])), key=lambda t: t[1])
            qf, mf, _ = min((fit(y, xt, loss="minphase", p0=pr), fit(y, xt, loss="minphase", p0=list(qd) + [1.0])), key=lambda t: t[1])
            p_obs = f_test(md, min(mf, s_ref), ne, 6, 1)[1]
            base = tp.model(qd, xt)
            rng = np.random.default_rng(60)
            pb = []
            for _ in range(NB):
                yb = base + blocks(rng, res, max(corr, 2))
                _, mbf, rbf = min((fit(yb, xt, loss="minphase", p0=qf), fit(yb, xt, loss="minphase", p0=list(qd) + [1.0])), key=lambda t: t[1])
                _, mbd, _ = fit(yb, xt, p0=qd)
                pb.append(f_test(mbd, mbf, n_eff(rbf)[0], 6, 1)[1])
            pb = np.array(pb)
            d = dict(size=v, p_nominal=p_obs, p_cal=float(np.mean(pb <= p_obs)), thr=float(np.quantile(pb, 0.05)),
                     rej_nominal=float(np.mean(pb < 0.05)))
            out["A"][f"{vn} | {c} | {sc}"] = d
            say(f"{vn:22s} {c[:22]:22s} {sc[:5]} size {v:g}: nominal p {p_obs:.4f}, calibrated p {d['p_cal']:.3f} "
                f"(null rejects {d['rej_nominal']:.2f} at nominal 0.05; 5 % threshold {d['thr']:.3f})")
        (OUT / "sim60.json").write_text(json.dumps(out, indent=1, default=float))

    say()
    # (B) second-pulse lead, a <= 1 and a + b <= 1 with the 2.07 floor
    for vn in ("floor 2.07, a <= 1", "floor 2.07, a + b <= 1"):
        ab = 1.0 if "a + b" in vn else None
        floor = tp.loss_transfer(2.07, "zero")
        fit = lambda y, x, **kw: tp.fit(y, x, a_max=1.0, ab_max=ab, **kw)
        for sc in LATE:
            y, x0 = SCANS[sc]
            xf = tp.apply_transfer(x0, floor)
            pr, _, _ = fit(y, xf, loss="minphase")
            row = {"reference (loss increment, s free)": -(pr[2] + pr[3])}
            node = max((n for n in S59[vn][sc] if n.startswith("node ")), key=lambda n: S59[vn][sc][n]["p_alone"])
            for c in ("thicker path (fraction)", "gas, slice swells, r 30 um", "gas, slice swells, r 50 um", node):
                for mode in ("alone", "drift"):
                    v = S59[vn][sc][c][f"best_{mode}"]
                    xt = tp.apply_transfer(xf, transfer(c, v))
                    q, _, _ = fit(y, xt, s=0.0) if mode == "alone" else fit(y, xt)
                    lead = -q[2] if mode == "alone" else -(q[2] + q[3])
                    row[f"{'node (best)' if c.startswith('node ') else c[:26]}, {mode} ({v:g}; p {S59[vn][sc][c][f'p_{mode}']:.3f})"] = float(lead)
            out["B"][f"{vn} | {sc}"] = row
            say(f"{vn} {sc[:5]} second pulse leads the empty cuvette by (ps): " + "; ".join(f"{k} {x:+.3f}" for k, x in row.items()))
    log.outcome = "; ".join(f"{k.split(' | ')[1][:10]} {k.split(' | ')[2][:5]} cal {d['p_cal']:.3f}" for k, d in list(out["A"].items())[:8])[:250]
(OUT / "sim60.json").write_text(json.dumps(out, indent=1, default=float))
(OUT / "sim60.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
