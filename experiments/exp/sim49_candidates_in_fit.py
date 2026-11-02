"""Sim 49 - Review round 19: each candidate change of the slice applied inside the waveform fit.

Sim 46 compared candidates with the required delay at 0.6 THz only. Here each candidate's whole-band change
of the slab transfer (slab model of Sec. 5.5/5.7: warming, water gain into thickness, ions in the tissue
water, axoplasm conduction, a thicker path; `thznerve.physics.slice_change`) is applied to the unbiased pulse
B, and the biased pulse is fitted as a B' + b K(t - tau) with a, b, tau free and no shift (s = 0), as Sim 48
does for node states. For each candidate the size that fits best is found; reported are that size, the
share of the misfit reduction a free shift s buys that it recovers, an F test of its fit against s free
(1 and N_eff - 4 dof), the fitted a, and its loss and delay at 0.6 THz. Each fit is done with a free and
with a <= 1: a scales the tissue pulse after the candidate's own change, and a > 1 would need the passive
slice to transmit more than the candidate allows at every frequency. The best node states of Sim 48 (every
extreme) are refitted the same way. Writes results/sim49/sim49.txt and results/sim49/sim49.json.
"""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

from thznerve.physics.slice_change import k_millon, mechanisms, ratio
from thznerve.waveform_fit import TwoPulse, f_test, load_scans, n_eff, node_transfer

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim49"; OUT.mkdir(parents=True, exist_ok=True)
lines = []


def say(s=""):
    print(s, flush=True); lines.append(s)


tt, K, SCANS = load_scans(ROOT)
tp = TwoPulse(tt, K)
F_MAX = 4.0                                              # THz: transfer held constant above (pulse content is noise there)
fr = tp.fr; band = (fr > 0) & (fr <= F_MAX)
MECH = mechanisms(k_millon())
RANGES = {"heating, Liebe (dT, K)": (0, 120), "water content (fractional change)": (0, 0.4),
          "ions: all tissue water (S/m)": (0, 200), "axoplasm conductivity (S/m)": (0, 800),
          "thicker path (fraction)": (0, 0.4)}
RANGES[[n for n in MECH if n.startswith("heating, Liebe x")][0]] = (0, 40)


def full_ratio(fn):
    r = np.ones(len(fr), complex)
    r[band] = ratio(fr[band] * 1e12, fn)
    r[fr > F_MAX] = r[band][-1]
    return r


cache = {}


def get_ratio(name, v):
    key = (name, round(float(v), 6))
    if key not in cache:
        cache[key] = full_ratio(MECH[name][0](v))
    return cache[key]


NODES = {st["state"]: st for st in json.loads((ROOT / "results" / "sim30" / "sim36_signedband.json").read_text())["states"]}
NODE_PICK = ("joint:joint_-2.7", "epsfine:-2.5_0.1", "epsnode:-3_0.3", "sigma:1e+06", "epsnode:-10_30.0")
out = {}
say(f"Sim 49: candidate changes of the slice inside the waveform fit (s = 0; a, b, tau free, then a <= 1); transfer above {F_MAX} THz held")
for scan in ("2 min (H2)", "4 min (H1)", "5 min (H3)"):
    y, x0 = SCANS[scan]
    p, s_free, res = tp.fit(y, x0)
    ne, _ = n_eff(res)
    _, s_zero, _ = tp.fit(y, x0, s=0.0)
    gain = s_zero - s_free
    say(f"  {scan}: s = {p[2]:+.4f} ps, a = {p[0]:.3f}; s free removes {gain / s_zero:.1%} of the s = 0 misfit (N_eff {ne:.1f})")
    out[scan] = {}
    for amax in (np.inf, 1.0):
        say(f"   a {'free' if amax == np.inf else '<= 1'}:")
        for name, (lo, hi) in RANGES.items():
            mis = lambda v: tp.fit(y, tp.apply_transfer(x0, get_ratio(name, v)), s=0.0, a_max=amax)[1]
            grid = np.linspace(lo, hi, 17)
            vals = [mis(v) for v in grid]
            i = int(np.argmin(vals))
            a_, b_ = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
            v = minimize_scalar(mis, bounds=(a_, b_), method="bounded", options={"xatol": (hi - lo) * 1e-3}).x if a_ < b_ else grid[i]
            q, sn, _ = tp.fit(y, tp.apply_transfer(x0, get_ratio(name, v)), s=0.0, a_max=amax)
            F, pv = f_test(sn, s_free, ne)
            r06 = ratio([0.6e12], MECH[name][0](v))[0]
            rec = (s_zero - sn) / gain
            edge = " (end of range)" if i == len(grid) - 1 else ""
            out[scan][f"{name} | a {'free' if amax == np.inf else '<= 1'}"] = dict(
                value=float(v), recovered=float(rec), F=float(F), p=pv, a=float(q[0]), b=float(q[1]),
                loss06_db=float(-20 * np.log10(abs(r06))), delay06_rad=float(-np.angle(r06)))
            say(f"    {name:36s} best {v:8.3g}{edge}: recovers {rec:4.0%}, F = {F:5.1f}, p = {pv:.3f}; a = {q[0]:.3f}; "
                f"at 0.6 THz {-20 * np.log10(abs(r06)):+.2f} dB, {-np.angle(r06):+.3f} rad")
        for st in NODE_PICK:
            r = node_transfer(NODES[st], fr)
            q, sn, _ = tp.fit(y, tp.apply_transfer(x0, r), s=0.0, a_max=amax)
            F, pv = f_test(sn, s_free, ne)
            rec = (s_zero - sn) / gain
            out[scan][f"node {st} | a {'free' if amax == np.inf else '<= 1'}"] = dict(recovered=float(rec), F=float(F), p=pv, a=float(q[0]))
            say(f"    node state {st:28s} every extreme: recovers {rec:4.0%}, F = {F:5.1f}, p = {pv:.3f}; a = {q[0]:.3f}")
(OUT / "sim49.json").write_text(json.dumps(out, indent=1))
(OUT / "sim49.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
