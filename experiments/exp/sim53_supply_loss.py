"""Sim 53 - Can any candidate supply the bias-induced change by itself?

Sim 52: with the loss present at 30 s taken as a zero-phase floor (kappa0, before the reported onset) and the
growth after it as a causal (minimum-phase) loss increment, the biased waveforms need no shift (s = 0 at
p = 0.04-0.89) and fit better than with a zero-phase loss. The bias-induced change is then a smooth causal
extra loss of the tissue pulse, whose own phase carries the lag at 0.6 THz. Here each candidate must supply
that change itself: the reference is floor x causal increment (kappa) x shift (s), and the candidate C_v
replaces the increment and the shift,
    full:       a [C_v L_inc(kappa) F B](t - s) + b K(t - s - tau)      (v, kappa, s, a, b, tau free)
    candidate:  a [C_v F B](t)              + b K(t - tau)              (kappa = s = 0)
F the zero-phase floor. The candidate supplies the change if the two constraints are not rejected
(F(2, N_eff - 6), N_eff from the full fit's and the candidate's own residual; 'edge' if they disagree).
Sizes that pass, the candidate's loss and delay at 0.6 THz, and its loss slope 0.3-1.5 THz against the
fitted increment are reported. The increment is held constant above 1.5, 2 or 4 THz; a <= 1, and as a
check a + b <= 1. Writes results/sim53/sim53.txt and results/sim53/sim53.json.
"""
import json
from pathlib import Path

import numpy as np

from thznerve.physics.slice_change import k_millon, mechanisms, ratio
from thznerve.waveform_fit import TwoPulse, f_test, load_scans, n_eff, node_transfer

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim53"; OUT.mkdir(parents=True, exist_ok=True)
lines = []


def say(s=""):
    print(s, flush=True); lines.append(s)


tt, K, SCANS = load_scans(ROOT)
F0, P_PASS = 0.6, 0.05
MECH = mechanisms(k_millon())
HX = [n for n in MECH if n.startswith("heating, Liebe x")][0]
GRIDS = {"thicker path (fraction)": np.linspace(0.0, 0.2, 11),
         "ions: all tissue water (S/m)": np.linspace(0.0, 60.0, 13),
         "axoplasm conductivity (S/m)": np.linspace(0.0, 150.0, 11),
         "heating, Liebe (dT, K)": np.linspace(0.0, 30.0, 11),
         HX: np.linspace(0.0, 12.0, 13)}
LAMS = np.linspace(0.0, 1.0, 11)
NODES = json.loads((ROOT / "results" / "sim30" / "sim36_signedband.json").read_text())["states"]
tp0 = TwoPulse(tt, K)
KAPPA0 = float(tp0.fit(*SCANS["30 s (H1)"], loss="zero", a_max=1.0)[0][4])
fr = tp0.fr; band = (fr > 0) & (fr <= 4.0); i06 = int(np.argmin(np.abs(fr - F0)))
FLOOR = tp0.loss_transfer(KAPPA0, "zero")
sl = (fr >= 0.3) & (fr <= 1.5)


def slope_db(r):
    """Loss slope (dB per THz of field amplitude) of a transfer over 0.3-1.5 THz."""
    return float(np.polyfit(fr[sl], -20 * np.log10(np.abs(r[sl])), 1)[0])


cache = {}


def slab_ratio(name, v):
    key = (name, round(float(v), 9))
    if key not in cache:
        r = np.ones(len(fr), complex)
        if v != 0:
            r[band] = ratio(fr[band] * 1e12, MECH[name][0](v))
            r[fr > 4.0] = r[band][-1]
        cache[key] = r
    return cache[key]


def verdict(p1, p2):
    return "pass" if min(p1, p2) >= P_PASS else ("fail" if max(p1, p2) < P_PASS else "edge")


say(f"Floor (30 s, zero phase): {KAPPA0:.2f} dB/THz")
out = {"kappa0": KAPPA0}
for fl in (1.5, 2.0, 4.0):
    tp = TwoPulse(tt, K, f_loss=fl)
    for ab in (None, 1.0):
        vn = f"increment held above {fl:g} THz, " + ("a + b <= 1" if ab else "a <= 1")
        say(f"=== {vn} ==="); out[vn] = {}
        fit = lambda y, x, **kw: tp.fit(y, x, a_max=1.0, ab_max=ab, **kw)
        for scan, (y, x0) in SCANS.items():
            xf = tp.apply_transfer(x0, FLOOR)
            pr, s_ref, res = fit(y, xf, loss="minphase")
            ne, _ = n_eff(res)
            inc = tp.loss_transfer(pr[4], "minphase")
            q0, s00, _ = fit(y, xf, s=0.0)                                   # no change at all
            p_none = f_test(s00, s_ref, ne, 5, 2)[1]
            rs = dict(kappa=float(pr[4]), s=float(pr[2]), misfit=s_ref, n_eff=float(ne), p_no_change=p_none,
                      inc_loss06_db=float(-20 * np.log10(abs(inc[i06]))), inc_delay06_rad=float(-np.angle(inc[i06])),
                      inc_slope=slope_db(inc), cands={})
            out[vn][scan] = rs
            say(f"{scan}: increment {pr[4]:.2f} dB/THz (at 0.6 THz {rs['inc_loss06_db']:.2f} dB, {rs['inc_delay06_rad']:.3f} rad), "
                f"s {pr[2]:+.4f} ps, misfit {s_ref:.3g}, N_eff {ne:.1f}; no change at all: p = {p_none:.4f}")
            scr = sorted(((fit(y, tp.apply_transfer(xf, node_transfer(st, fr)), s=0.0, p0=q0)[1], st["state"]) for st in NODES))
            cands = dict(GRIDS)
            for n in [n for _, n in scr[:3]]:
                cands["node " + n] = LAMS
            for name, grid in cands.items():
                if name.startswith("node "):
                    st = next(x for x in NODES if x["state"] == name[5:])
                    tr = lambda v, st=st: node_transfer(st, fr, lam=v)
                else:
                    tr = lambda v, name=name: slab_ratio(name, v)
                S_c, S_f, R_c, prev_c, prev_f = [], [], [], list(q0), list(pr)
                for v in grid:
                    xt = tp.apply_transfer(xf, tr(v))
                    qc, mc, rc = min((fit(y, xt, s=0.0, p0=prev_c), fit(y, xt, s=0.0, p0=q0)), key=lambda t: t[1])
                    qf, mf, _ = min((fit(y, xt, loss="minphase", p0=prev_f), fit(y, xt, loss="minphase", p0=pr)), key=lambda t: t[1])
                    prev_c, prev_f = list(qc), list(qf)
                    S_c.append(mc); S_f.append(mf); R_c.append(rc)
                S_full = min(min(S_f), s_ref)
                ib = int(np.argmin(S_c))
                ne_own, _ = n_eff(R_c[ib])
                p1 = f_test(S_c[ib], S_full, ne, 6, 2)[1]; p1b = f_test(S_c[ib], S_full, ne_own, 6, 2)[1]
                vd = verdict(p1, p1b) if p_none < P_PASS else "n/a"
                pv = [f_test(m, S_full, ne, 6, 3)[1] for m in S_c]
                rng = [float(v) for v, pp in zip(grid, pv) if pp >= P_PASS]
                r = tr(grid[ib])
                rs["cands"][name] = dict(verdict=vd, p=p1, p_own_neff=p1b, best=float(grid[ib]),
                                         pass_range=[min(rng), max(rng)] if rng else None,
                                         loss06_db=float(-20 * np.log10(abs(r[i06]))), delay06_rad=float(-np.angle(r[i06])),
                                         slope=slope_db(r))
                say(f"  {name:36s} {vd:4s} p = {p1:.3f} ({p1b:.3f}); best {grid[ib]:.3g}: at 0.6 THz "
                    f"{rs['cands'][name]['loss06_db']:+.2f} dB, {rs['cands'][name]['delay06_rad']:+.3f} rad, slope "
                    f"{rs['cands'][name]['slope']:+.2f} dB/THz; sizes passing " + (f"{min(rng):.3g} to {max(rng):.3g}" if rng else "none"))
        say()
say("=== Summary (2, 4, 5 min; 6 references) ===")
names = sorted({c for v in out if v != "kappa0" for sc in out[v].values() for c in sc["cands"] if not c.startswith("node ")})
summ = {}
for c in names + ["node (best three)"]:
    vs = []
    for v in out:
        if v == "kappa0":
            continue
        for sc, d in out[v].items():
            if sc.startswith("30 s"):
                continue
            if c.startswith("node"):
                cs = [x["verdict"] for n, x in d["cands"].items() if n.startswith("node ")]
                vs.append("pass" if "pass" in cs else ("edge" if "edge" in cs else ("fail" if "fail" in cs else "n/a")))
            else:
                vs.append(d["cands"][c]["verdict"])
    summ[c] = {k: vs.count(k) for k in ("pass", "edge", "fail", "n/a")}
    say(f"  {c:36s} " + ", ".join(f"{k} {summ[c][k]}" for k in ("pass", "edge", "fail", "n/a")) + f" of {len(vs)}")
out["summary"] = summ
(OUT / "sim53.json").write_text(json.dumps(out, indent=1))
(OUT / "sim53.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
