"""Sim 55 - Review round 23: does the waveform result depend on the 30 s floor, and are the F tests calibrated?

(A) The reference of Sims 52-54 fixes a zero-phase loss floor at the value fitted to the 30 s scan (H1) and
applies it to H2 and H3, whose unbiased traces come from other figures. Here the floor kappa0 runs from 0
(a purely causal loss) to 2.5 dB/THz, the causal increment is held above 1.5, 2 or 4 THz, and a <= 1 or
a + b <= 1. For each reference: the 30 s control (s = 0 not rejected), the increment at 2-5 min, the
increment over the 30 s value (the bias-induced loss slope, dB/THz), the net lag at 0.6 THz, and whether
s = 0 and "no change at all" are rejected.
(B) Calibration: synthetic 'no change' scans are the reference fit to the 30 s scan plus its residual
resampled in circular blocks (block = 2 correlation lengths); each is tested as the real scans are. The
fraction rejected at nominal 0.05 gives the false-rejection rate of the shift test (F(1, N_eff - 5)) and
of the no-change test (F(2, N_eff - 5)).
Writes results/sim55/sim55.txt and results/sim55/sim55.json.
"""
import json
from pathlib import Path

import numpy as np

from thznerve.runlog import RunLog
from thznerve.waveform_fit import TwoPulse, f_test, load_scans, n_eff

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim55"; OUT.mkdir(parents=True, exist_ok=True)
lines = []


def say(s=""):
    print(s, flush=True); lines.append(s)


tt, K, SCANS = load_scans(ROOT)
F0, P_PASS = 0.6, 0.05
LATE = ("2 min (H2)", "4 min (H1)", "5 min (H3)")
out = {"A": {}, "B": {}}
with RunLog("sim55", "round 23: floor 0-2.5 dB/THz x hold x bound; block-bootstrap calibration of the shift and no-change F tests") as log:
    for k0 in (0.0, 0.5, 1.0, 1.5, 2.07, 2.5):
        for fl in (1.5, 2.0, 4.0):
            tp = TwoPulse(tt, K, f_loss=fl)
            i06 = int(np.argmin(np.abs(tp.fr - F0)))
            floor = tp.loss_transfer(k0, "zero")
            for ab in (None, 1.0):
                fit = lambda y, x, **kw: tp.fit(y, x, a_max=1.0, ab_max=ab, loss=kw.pop("loss", None), **kw)
                vn = f"floor {k0:g}, held above {fl:g} THz, " + ("a + b <= 1" if ab else "a <= 1")
                rs = {}
                for scan, (y, x0) in SCANS.items():
                    xf = tp.apply_transfer(x0, floor)
                    pr, m, res = fit(y, xf, loss="minphase")
                    ne, _ = n_eff(res)
                    _, m_s0, _ = fit(y, xf, s=0.0, loss="minphase")
                    _, m_none, _ = fit(y, xf, s=0.0)
                    inc = tp.loss_transfer(pr[4], "minphase")
                    rs[scan] = dict(kappa=float(pr[4]), s=float(pr[2]), misfit=m, n_eff=float(ne),
                                    p_s0=f_test(m_s0, m, ne, 5, 1)[1], p_none=f_test(m_none, m, ne, 5, 2)[1],
                                    lag=float(2 * np.pi * F0 * pr[2] - np.angle(inc[i06])))
                k30 = rs["30 s (H1)"]["kappa"]
                for sc in LATE:
                    rs[sc]["d_kappa"] = rs[sc]["kappa"] - k30
                out["A"][vn] = rs
                say(f"{vn:42s} 30 s: kappa {k30:.2f}, s = 0 p {rs['30 s (H1)']['p_s0']:.3f} | " + "; ".join(
                    f"{sc[:5]} dk {rs[sc]['d_kappa']:+.2f} lag {rs[sc]['lag']:.3f} p(s=0) {rs[sc]['p_s0']:.3f} "
                    f"p(none) {rs[sc]['p_none']:.4f} misfit {rs[sc]['misfit']:.2f}" for sc in LATE))
    valid = [v for v, r in out["A"].items() if r["30 s (H1)"]["p_s0"] >= P_PASS]
    allr = list(out["A"])
    for grp, vs in (("every reference", allr), ("references passing the 30 s control", valid)):
        dk = [out["A"][v][sc]["d_kappa"] for v in vs for sc in LATE]
        lag = [out["A"][v][sc]["lag"] for v in vs for sc in LATE]
        pn = [out["A"][v][sc]["p_none"] for v in vs for sc in LATE]
        per = {sc: (min(out["A"][v][sc]["d_kappa"] for v in vs), max(out["A"][v][sc]["d_kappa"] for v in vs)) for sc in LATE}
        say(f"{grp} ({len(vs)}): increment over 30 s " + ", ".join(f"{sc[:5]} {a:.2f}-{b:.2f}" for sc, (a, b) in per.items())
            + f" dB/THz; lag {min(lag):.3f}-{max(lag):.3f} rad; no change rejected at p <= {max(pn):.4f}")
        out["A"][f"summary: {grp}"] = dict(n=len(vs), d_kappa=per, lag=[min(lag), max(lag)], p_none_max=max(pn))
    say(f"floors passing the 30 s control: {sorted({float(v.split()[1].rstrip(',')) for v in valid})}")

    # (B) calibration on the 30 s scan, reference floor 2.07 held above 2 THz, a <= 1
    tp = TwoPulse(tt, K, f_loss=2.0)
    floor = tp.loss_transfer(2.07, "zero")
    fit = lambda y, x, **kw: tp.fit(y, x, a_max=1.0, loss=kw.pop("loss", None), **kw)
    y, x0 = SCANS["30 s (H1)"]
    xf = tp.apply_transfer(x0, floor)
    pr, m, res = fit(y, xf, loss="minphase")
    base = tp.model(pr, xf, loss="minphase")              # the fitted no-change scan (res = model - y)
    _, corr = n_eff(res)
    blk = max(2 * corr, 2)
    rng = np.random.default_rng(55)
    n = len(res); NB = 200
    rej_s, rej_none = 0, 0
    for _ in range(NB):
        starts = rng.integers(0, n, size=n // blk + 1)
        r = np.concatenate([np.take(res, np.arange(s0, s0 + blk), mode="wrap") for s0 in starts])[:n]
        yb = base + r
        pb, mb, resb = fit(yb, xf, loss="minphase", p0=pr)
        neb, _ = n_eff(resb)
        _, mb0, _ = fit(yb, xf, s=0.0, loss="minphase")
        _, mbn, _ = fit(yb, xf, s=0.0)
        rej_s += f_test(mb0, mb, neb, 5, 1)[1] < P_PASS
        rej_none += f_test(mbn, mb, neb, 5, 2)[1] < P_PASS
    out["B"] = dict(n_boot=NB, block=int(blk), false_rej_shift=rej_s / NB, false_rej_none=rej_none / NB)
    say(f"calibration ({NB} block-bootstrap 'no change' scans, block {blk} samples): false rejection at nominal 0.05: "
        f"shift test {rej_s / NB:.3f}, no-change test {rej_none / NB:.3f}")
    log.outcome = lines[-1][:200]
(OUT / "sim55.json").write_text(json.dumps(out, indent=1))
(OUT / "sim55.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
