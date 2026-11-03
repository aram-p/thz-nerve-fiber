"""Sim 56 - Do gas bubbles and a thicker tissue path together supply the change?

Sim 54 (corrected Mie): gas bubbles in a swelling slice give the steep loss (5-8 dB/THz at r 30-50 um) but
almost no delay; a thicker tissue path (Sim 53) gives the delay but a flat loss. Here both act at once:
spheres of gas (r = 30 um, the largest that sits well inside a 100 um slice) at volume fraction v in tissue
whose thickness is D (1 + t) / (1 - v), t the extra tissue path. A post-hoc two-cause model: the test
charges for both sizes (full model 7 parameters; supply_test n_size = 2). References: the 30 s floor times a
causal increment held above 2 THz, a <= 1 and a + b <= 1. p is given with the nominal F test and judged also
against the calibrated 2-dof threshold of Sim 55b. Writes results/sim56/sim56.txt and results/sim56/sim56.json.
"""
import json
from pathlib import Path

import numpy as np

from thznerve.physics.scattering import foldy_eps
from thznerve.physics.slice_change import D_UM, eps_tissue, ratio
from thznerve.runlog import RunLog
from thznerve.waveform_fit import TwoPulse, load_scans, n_eff, supply_test

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim56"; OUT.mkdir(parents=True, exist_ok=True)
lines = []


def say(s=""):
    print(s, flush=True); lines.append(s)


tt, K, SCANS = load_scans(ROOT)
cal = json.loads((ROOT / "results" / "sim55" / "sim55b.json").read_text())
THR2 = min(cal["block x1"]["thr2"], cal["block x2"]["thr2"])
RAD = 30
GRID = [(round(t, 3), round(v, 3)) for t in np.arange(0, 0.161, 0.02) for v in np.arange(0, 0.201, 0.04)]
tp = TwoPulse(tt, K, f_loss=2.0)
KAPPA0 = float(TwoPulse(tt, K).fit(*SCANS["30 s (H1)"], loss="zero", a_max=1.0)[0][4])
fr = tp.fr; band = (fr > 0) & (fr <= 4.0)
FLOOR = tp.loss_transfer(KAPPA0, "zero")
cache = {}


def tr(tv):
    if tv not in cache:
        t, v = tv
        r = np.ones(len(fr), complex)
        if t or v:
            d = D_UM * (1 + t) / (1 - v)
            r[band] = ratio(fr[band] * 1e12, lambda f: (foldy_eps(f, eps_tissue(f), 1.0, RAD, v), d))
            r[fr > 4.0] = r[band][-1]
        cache[tv] = r
    return cache[tv]


out = {"thr2": THR2}
with RunLog("sim56", "gas (r 30 um, v 0-0.2) plus a thicker tissue path (0-16 %), two sizes charged, increment above 2 THz") as log:
    say(f"calibrated 2-dof threshold (Sim 55b): p < {THR2:.4f}")
    for ab in (None, 1.0):
        fit = lambda y, x, **kw: tp.fit(y, x, a_max=1.0, ab_max=ab, **kw)
        for scan in ("2 min (H2)", "4 min (H1)", "5 min (H3)"):
            y, x0 = SCANS[scan]
            xf = tp.apply_transfer(x0, FLOOR)
            pr, s_ref, res = fit(y, xf, loss="minphase")
            ne, _ = n_eff(res)
            q0, _, _ = fit(y, xf, s=0.0)
            d = supply_test(tp, fit, y, xf, tr, GRID, q0, pr, s_ref, ne, n_size=2)
            d["calibrated"] = "pass" if d["p"] >= THR2 else "fail"
            out[f"{'a + b <= 1' if ab else 'a <= 1'}, {scan}"] = d
            ok = d["pass_range"]
            say(f"{'a + b <= 1' if ab else 'a <= 1':10s} {scan}: nominal {d['verdict']:4s} p = {d['p']:.3f} ({d['p_own_neff']:.3f}); "
                f"calibrated {d['calibrated']}; best (t, v) = {tuple(d['best'])}; combinations passing nominal: "
                + (f"t {min(x[0] for x in ok):.2f}-{max(x[0] for x in ok):.2f}, v {min(x[1] for x in ok):.2f}-{max(x[1] for x in ok):.2f}" if ok else "none"))
    log.outcome = "; ".join(f"{k}: {v['verdict']}/{v['calibrated']}" for k, v in out.items() if k != "thr2")[:200]
(OUT / "sim56.json").write_text(json.dumps(out, indent=1, default=float))
(OUT / "sim56.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
