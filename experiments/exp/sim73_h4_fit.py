"""Sim 73b - Fit the H4 waveforms (Sim 73a traces) with the paper's two-pulse model.

y = a [x0 with a tissue change](t - s) + b K(t - s - tau), x0 the same sample with no voltage, K the H4 empty-cuvette
pulse (Fig. 3a), on the one time axis that H4's figures share (the three unbiased samples peak within 0.05 ps).
tau = 0 means the second pulse arrives when a pulse through the empty cuvette does. Candidates as in Table 2:
no change, shift only, full fit (causal loss, s free), uniform thicker or thinner path (-0.5 to +0.5), and uneven
thickness (mean -0.3 to 0.51, rms 0-0.4; thznerve.physics.slice_change.uneven), each with a <= 1, alone and with a
drift. Tracing spread: every fit redone with each trace shifted by -1, 0, +1 pixel (9 combinations per sample).
Also the frequency of the deepest notch of each 6 min spectrum, of the fitted models, and 1 / (2 |dt|), the
first zero of two equal pulses dt apart. Writes results/sim73/sim73.{txt,json}.
"""
import itertools
import json
from pathlib import Path

import numpy as np

from thznerve.physics.slice_change import uneven
from thznerve.runlog import RunLog
from thznerve.waveform_fit import TwoPulse, effect_size, multistart

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim73"


def load(name, file):
    import importlib.util
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / "exp" / file)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


A = load("s73a", "sim73_h4_trace.py")
S64 = load("s64", "sim64_size_fitted_tests.py")
TT = A.TT
MS = np.round(np.arange(-0.30, 0.511, 0.03), 3)
SS = np.array([0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.25, 0.3, 0.35, 0.4])
VG = np.round(np.arange(-0.95, 1.801, 0.01), 3)
UNI = np.round(np.arange(-0.5, 0.51, 0.02), 3)


def run(W, valid, t2_max=None):
    tp = TwoPulse(TT, W["empty"], f_loss=2.0)
    # the thickness transfers live on Sim 64's frequency grid; interpolate onto this one
    fr64 = S64.fr
    T = np.array([np.interp(tp.fr, fr64, S64.transfer("thicker path (fraction)", float(v)).real) +
                  1j * np.interp(tp.fr, fr64, S64.transfer("thicker path (fraction)", float(v)).imag) for v in VG])
    one = lambda y, x, **kw: tp.fit(y, x, a_max=1.0, t2_max=t2_max, **kw)
    fit = multistart(one)
    out = {}
    for s in ("s1", "s2"):
        y, x0 = np.where(valid, W[f"{s}_6min"], 0.0), np.where(valid, W[f"{s}_none"], 0.0)
        qf, Sf, _ = fit(y, x0, loss="minphase")
        r = {"full": dict(q=list(qf), resid=float(np.sqrt(Sf / np.sum(y ** 2))))}
        for name, grid in (("uniform", [(float(v), 0.0) for v in UNI]), ("uneven", [(float(m), float(sd)) for m in MS for sd in SS]), ("none", [(0.0, 0.0)])):
            xs = {g: tp.apply_transfer(x0, uneven(T, VG, *g)) if g != (0.0, 0.0) else x0 for g in grid}
            for mode in ("alone", "drift"):
                sv = 0.0 if mode == "alone" else None
                scr = sorted((one(y, xx, s=sv)[1], g) for g, xx in xs.items())[:3]
                S, g = min((fit(y, xs[g], s=sv)[1], g) for _, g in scr)
                q = fit(y, xs[g], s=sv)[0]
                r[f"{name} | {mode}"] = dict(D=effect_size(S, Sf, y), size=list(g), q=[float(v) for v in q],
                                             resid=float(np.sqrt(S / np.sum(y ** 2))))
        out[s] = r
    return out


def fmt(q):
    """(a, b, tau) with s fixed at 0; (a, b, s, tau) with a drift; (a, b, s, tau, kappa) for the full fit."""
    names = {3: ("a", "b", "tau"), 4: ("a", "b", "s", "tau"), 5: ("a", "b", "s", "tau", "kappa")}[len(q)]
    return ", ".join(f"{n} {v:+.3f}" for n, v in zip(names, q))


def notch(y, valid, lo=0.2, hi=1.5):
    n = 8192; f = np.fft.rfftfreq(n, TT[1] - TT[0]); Y = 20 * np.log10(np.abs(np.fft.rfft(np.where(valid, y, 0), n)) + 1e-12)
    k = (f > lo) & (f < hi)
    return float(f[k][np.argmin(Y[k])])


if __name__ == "__main__":
    lines = []
    say = lambda s="": (print(s, flush=True), lines.append(s))
    with RunLog("sim73", "H4 (Reports NAS RA 2024) traced and fitted: two-pulse model with the H4 empty pulse; no change, full fit, uniform and uneven thickness; +-1 pixel tracing spread") as log:
        W, valid = A.traces()
        base = run(W, valid)
        spread = []
        names = ["empty", "s1_none", "s1_6min", "s2_none", "s2_6min"]
        for sh in itertools.product((-1, 1), repeat=3):  # empty, unbiased, biased shifted by +-1 pixel
            d = {"empty": sh[0], "s1_none": sh[1], "s2_none": sh[1], "s1_6min": sh[2], "s2_6min": sh[2]}
            spread.append(run(*A.traces(d)))
        for s in ("s1", "s2"):
            b = base[s]
            say(f"H4 sample {s[1]}, 6 min at 100 V (residual = rms misfit / rms trace; D = effect size against the full fit, %)")
            for k, v in b.items():
                sp = [x[s][k] for x in spread]
                q = v["q"]
                par = fmt(q)
                lo, hi = min(x["resid"] for x in sp), max(x["resid"] for x in sp)
                say(f"  {k:16s} residual {100 * v['resid']:5.1f} % ({100 * lo:.1f}-{100 * hi:.1f} over tracing shifts)" +
                    (f"  D {v['D']:5.1f} % ({min(x['D'] for x in sp):.1f}-{max(x['D'] for x in sp):.1f})  size {v['size']}" if "D" in v else "") + f"  {par}")
            say()
        for s in ("s1", "s2"):
            q = base[s]["none | drift"]["q"]
            t_t = TT[np.argmin(W[f"{s}_none"])] + q[2]; t_k = TT[np.argmin(W["empty"])] + q[2] + q[3]
            say(f"  {s}: tissue pulse minimum at {t_t:.3f} ps, second pulse minimum at {t_k:.3f} ps (no change, drift fit); "
                f"first zero expected near {1 / (2 * abs(t_t - t_k)):.2f} THz; deepest notch of the 6 min spectrum (0.2-1.5 THz) at "
                f"{notch(W[f'{s}_6min'], valid):.2f} THz")
            tp = TwoPulse(TT, W["empty"], f_loss=2.0)
            for k in ("none | drift", "full"):
                m = tp.model(base[s][k]["q"], np.where(valid, W[f"{s}_none"], 0.0), s=None, loss="minphase" if k == "full" else None)
                say(f"      deepest notch of the fitted model ({k}): {notch(m, valid):.2f} THz")
        log.outcome = "; ".join(f"{s} b {base[s]['full']['q'][1]:.2f} tau {base[s]['full']['q'][3]:+.2f} resid {100 * base[s]['full']['resid']:.0f} %" for s in ("s1", "s2"))
    (OUT / "sim73.json").write_text(json.dumps({"base": base, "spread": spread}, default=float))
    (OUT / "sim73.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
