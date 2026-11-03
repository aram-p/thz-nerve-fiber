"""Sim 52 - Review round 22: candidate causes against reference fits that pass the 30 s control.

Round 22 showed that the causal (minimum-phase) references of Sim 51 find a shift at 30 s, before any lag
exists (s = 0 rejected at p = 0.008), because the loss term's own phase (0.3-1.1 rad at 0.6 THz) exceeds the
lag; at 2-5 min they then test candidates against an artefact. The 30 s scan is used here as a control: a
reference fit is valid only if it does not reject s = 0 at 30 s. Besides no loss term and a zero-phase loss,
a mixed loss is tried: a zero-phase floor fixed at the 30 s value (the loss present before the reported
onset, which carries no phase) times a causal increment held constant above 1.5, 2 or 4 THz; each with
a <= 1 or a + b <= 1. Candidates are tested at all four scans (at 30 s none should be needed) by the nested
test of Sim 51; a candidate test is reported as 'n/a' where the reference itself does not reject s = 0.
States an aqueous gap can hold are now those with eps' between lipid and water (2 to 5) and eps'' <= 3
(water's own ~3 at 0.6 THz plus ions up to ~100 S/m), or the conductivity scan's lowest value. The largest
delay at 0.6 THz any node state gives with every extreme is reported against the lag.
Writes results/sim52/sim52.txt and results/sim52/sim52.json.
"""
import json
from pathlib import Path

import numpy as np

from thznerve.physics.slice_change import k_millon, mechanisms, ratio
from thznerve.waveform_fit import TwoPulse, f_test, load_scans, n_eff, node_transfer

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim52"; OUT.mkdir(parents=True, exist_ok=True)
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
VARIANTS = [(kind, fl, ab) for kind, fl in ((None, 4.0), ("zero", 4.0), ("mixed", 1.5), ("mixed", 2.0), ("mixed", 4.0))
            for ab in (None, 1.0)]


def vname(kind, fl, ab):
    k = {None: "no loss term", "zero": "zero-phase loss", "mixed": f"30 s floor + causal increment held above {fl:g} THz"}[kind]
    return k + (", a + b <= 1" if ab else ", a <= 1")


def aqueous(state):
    fam, k = state.split(":")
    if fam == "sigma":
        return float(k) <= 25.0
    if fam != "epsnode":
        return False
    e1, e2 = (float(x) for x in k.split("_")[:2])
    return 2.0 <= e1 <= 5.0 and e2 <= 3.0


_tp = TwoPulse(tt, K)
KAPPA0 = float(_tp.fit(*SCANS["30 s (H1)"], loss="zero", a_max=1.0)[0][4])     # loss before the onset, zero phase
tps = {fl: TwoPulse(tt, K, f_loss=fl, kappa0=KAPPA0) for fl in (1.5, 2.0, 4.0)}
fr = tps[4.0].fr; band = (fr > 0) & (fr <= 4.0); i06 = int(np.argmin(np.abs(fr - F0)))
NODE_DELAY = max((float(-np.angle(node_transfer(st, fr)[i06])), st["state"]) for st in NODES)
AQ = [st["state"] for st in NODES if aqueous(st["state"])]
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


out = {}
for kind, fl, ab in VARIANTS:
    tp = tps[fl]; vn = vname(kind, fl, ab); out[vn] = {}
    npar = 4 + (1 if kind else 0)                             # a, b, s, tau, [kappa]
    fit = lambda y, x, **kw: tp.fit(y, x, loss=kind, a_max=1.0, ab_max=ab, **kw)
    say(f"=== {vn} ===")
    for scan, (y, x0) in SCANS.items():
        pr, s_ref, res = fit(y, x0)
        ne, _ = n_eff(res)
        q_s0, s0, _ = fit(y, x0, s=0.0)
        p_zero = f_test(s0, s_ref, ne, npar)[1]
        kap = float(pr[4]) if kind else 0.0
        own = -np.angle(tp.loss_transfer(kap, "minphase")[i06]) if kind else 0.0    # min-phase delay of the loss
        extra = own if kind == "mixed" else 0.0
        sg = np.arange(pr[2] - 0.14, pr[2] + 0.1401, 0.004)
        ok = []
        for sv in sg:
            q, ss, _ = fit(y, x0, s=sv, p0=[pr[0], pr[1], pr[3]] + ([pr[4]] if kind else []))
            if f_test(ss, s_ref, ne, npar)[1] >= P_PASS:
                ok.append(2 * np.pi * F0 * sv + (-np.angle(tp.loss_transfer(q[-1], "minphase")[i06]) if kind == "mixed" else 0.0))
        lag = 2 * np.pi * F0 * pr[2] + extra
        rs = dict(misfit=s_ref, n_eff=float(ne), kappa=kap, s=float(pr[2]), a=float(pr[0]), b=float(pr[1]),
                  p_zero=p_zero, lag=float(lag), lag_range=[float(min(ok)), float(max(ok))] if ok else None,
                  minphase_delay_of_loss=float(own), cands={})
        out[vn][scan] = rs
        say(f"{scan}: misfit {s_ref:.3g}, N_eff {ne:.1f}, kappa {kap:.1f} dB/THz (its minimum-phase delay at 0.6 THz "
            f"{own:.2f} rad), s {pr[2]:+.4f} ps, a {pr[0]:.3f}, b {pr[1]:.3f}; s = 0: p = {p_zero:.4f}; lag at 0.6 THz "
            f"{lag:.3f} rad, 95 % " + (f"{min(ok):.3f} to {max(ok):.3f}" if ok else "none"))
        # node states: screen at every extreme with s = 0, keep the best two and the best two aqueous
        scr = sorted(((fit(y, tp.apply_transfer(x0, node_transfer(st, fr)), s=0.0, p0=q_s0)[1], st["state"]) for st in NODES))
        pick = [n for _, n in scr[:2]] + [n for _, n in scr if aqueous(n)][:2]
        cands = dict(GRIDS)
        for n in pick:
            cands["node " + n] = LAMS
        for name, grid in cands.items():
            if name.startswith("node "):
                st = next(x for x in NODES if x["state"] == name[5:])
                tr = lambda v: node_transfer(st, fr, lam=v)
            else:
                tr = lambda v, name=name: slab_ratio(name, v)
            S0, Ss, R0 = [], [], []
            prev0, prevs = list(q_s0), list(pr)
            for v in grid:
                xt = tp.apply_transfer(x0, tr(v))
                f0 = [fit(y, xt, s=0.0, p0=prev0), fit(y, xt, s=0.0, p0=q_s0)]
                q0, m0, r0 = min(f0, key=lambda t: t[1])
                fs = [fit(y, xt, p0=prevs), fit(y, xt, p0=pr), fit(y, xt, p0=[q0[0], q0[1], 0.0] + list(q0[2:]))]
                qs, ms, _ = min(fs, key=lambda t: t[1])
                prev0, prevs = list(q0), list(qs)
                S0.append(m0); Ss.append(ms); R0.append(r0)
            S0, Ss = np.array(S0), np.array(Ss)
            S_full = min(Ss.min(), s_ref)
            ib = int(np.argmin(S0))
            ne_own, _ = n_eff(R0[ib])
            nf = npar + 1
            p1 = f_test(S0[ib], S_full, ne, nf)[1]; p1b = f_test(S0[ib], S_full, ne_own, nf)[1]
            p2 = [f_test(m, S_full, ne, nf, 2)[1] for m in S0]
            rng = [float(v) for v, pp in zip(grid, p2) if pp >= P_PASS]
            vd = verdict(p1, p1b) if p_zero < P_PASS else "n/a"
            rs["cands"][name] = dict(verdict=vd, p=p1, p_own_neff=p1b, best=float(grid[ib]),
                                     pass_range=[min(rng), max(rng)] if rng else None,
                                     loss06_db=float(-20 * np.log10(abs(tr(grid[ib])[i06]))),
                                     delay06_rad=float(-np.angle(tr(grid[ib])[i06])))
            say(f"  {name:36s} {vd:4s} p = {p1:.3f} ({p1b:.3f} with its own N_eff); best {grid[ib]:.3g}; sizes passing "
                + (f"{min(rng):.3g} to {max(rng):.3g}" if rng else "none"))
    say()

# what holds in every reference that passes the 30 s control
VALID = [v for v in out if out[v]["30 s (H1)"]["p_zero"] >= P_PASS]
say(f"Loss floor (30 s, zero phase): {KAPPA0:.2f} dB/THz. Largest delay at 0.6 THz of any node state with every extreme: "
    f"{NODE_DELAY[0]:.3f} rad ({NODE_DELAY[1]}). Aqueous states: {', '.join(AQ)}")
say(f"References passing the 30 s control: {len(VALID)} of {len(out)}: " + "; ".join(VALID))
say("=== Over the valid references (2, 4, 5 min) ===")
names = sorted({c for v in out.values() for sc in v.values() for c in sc["cands"]
                if not c.startswith("node ")})
summ = {}
for c in names:
    vs = [out[v][sc]["cands"][c]["verdict"] for v in VALID for sc in out[v] if c in out[v][sc]["cands"] and not sc.startswith("30 s")]
    summ[c] = {k: vs.count(k) for k in ("pass", "edge", "fail", "n/a")}
    say(f"  {c:36s} pass {summ[c]['pass']:2d}, edge {summ[c]['edge']:2d}, fail {summ[c]['fail']:2d}, n/a {summ[c]['n/a']:2d} of {len(vs)}")
for grp, test in (("node states, every extreme (best two)", lambda n: not aqueous(n[5:])),
                  ("node states an aqueous gap can hold (best two)", lambda n: aqueous(n[5:]))):
    vs = []
    for v in VALID:
        for sc in out[v]:
            cs = [d["verdict"] for n, d in out[v][sc]["cands"].items() if n.startswith("node ") and test(n)]
            if cs and not sc.startswith("30 s"):
                vs.append("pass" if "pass" in cs else ("edge" if "edge" in cs else ("fail" if "fail" in cs else "n/a")))
    summ[grp] = {k: vs.count(k) for k in ("pass", "edge", "fail", "n/a")}
    say(f"  {grp:36s} pass {summ[grp]['pass']:2d}, edge {summ[grp]['edge']:2d}, fail {summ[grp]['fail']:2d}, n/a {summ[grp]['n/a']:2d} of {len(vs)}")
lags = [(out[v][sc]["lag"], out[v][sc]["lag_range"]) for v in VALID for sc in out[v] if not sc.startswith("30 s")]
say(f"  lag at 0.6 THz, 2-5 min: best {min(l for l, _ in lags):.3f} to {max(l for l, _ in lags):.3f} rad; smallest 95 % "
    f"{min(r[0] for _, r in lags if r):.3f} rad; at 30 s best {min(out[v]['30 s (H1)']['lag'] for v in out):+.3f} to "
    f"{max(out[v]['30 s (H1)']['lag'] for v in out):+.3f} rad")
out["summary"] = dict(verdicts=summ, valid=VALID, kappa0=KAPPA0, node_delay_max=NODE_DELAY, aqueous=AQ)
(OUT / "sim52.json").write_text(json.dumps(out, indent=1))
(OUT / "sim52.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
