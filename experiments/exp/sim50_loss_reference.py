"""Sim 50 - Review round 20: candidate causes against a reference fit with a smooth loss term, one rule.

Round 20 found that the s-free fit leaves a structured misfit: a smooth extra loss of the tissue pulse
(5-10 dB/THz) halves it, so the F tests of Sims 48-49 compared candidates with an understated noise level, and
a loss with no delay nearly passed. Here the reference is
    biased(t) = a [L_kappa * B](t - s) + b K(t - s - tau),
L_kappa a smooth loss of kappa dB per THz above 0.2 THz, either zero-phase (a smoothing, e.g. of the traced curve) or causal
(minimum phase, which adds its own delay). Every candidate is fitted with s = 0 and kappa free (a <= 1, b, tau
free); a <= 1 in every fit, the reference included (the loss term is 1 at 0.2 THz, so the
tissue pulse stays at or below the unbiased one over the band); one rule throughout: a candidate reproduces the lag if its F test against the reference (1 and
N_eff - 5 dof) gives p >= 0.05. The same rule applied to a pure shift gives the 95 % interval of s and the
smallest lag at 0.6 THz that the traces allow. For each candidate the sizes that pass are reported (not only
the best one), with the share of the s = 0 misfit reduction its loss alone and its delay alone recover. Node
states (Sim 36 band export) are scaled down from every extreme with their delay scaled by the delay-based path
factor (2.2 at all frequencies) and the smallest scale at which any state passes is found.
Writes results/sim50/sim50.txt and results/sim50/sim50.json.
"""
import json
from pathlib import Path

import numpy as np

from thznerve.physics.slice_change import k_millon, mechanisms, ratio
from thznerve.waveform_fit import TwoPulse, f_test, load_scans, n_eff, node_transfer

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim50"; OUT.mkdir(parents=True, exist_ok=True)
lines = []


def say(s=""):
    print(s, flush=True); lines.append(s)


tt, K, SCANS = load_scans(ROOT)
tp = TwoPulse(tt, K)
fr = tp.fr; band = (fr > 0) & (fr <= 4.0)
F0 = 0.6; i06 = int(np.argmin(np.abs(fr - F0)))
P_PASS = 0.05
NPAR = 5
MECH = mechanisms(k_millon())
GRIDS = {"thicker path (fraction)": np.linspace(0.0, 0.2, 21),
         "water content (fractional change)": np.linspace(0.0, 0.24, 25),
         "ions: all tissue water (S/m)": np.linspace(0.0, 60.0, 31),
         "axoplasm conductivity (S/m)": np.linspace(0.0, 150.0, 31),
         "heating, Liebe (dT, K)": np.linspace(0.0, 30.0, 16)}
GRIDS[[n for n in MECH if n.startswith("heating, Liebe x")][0]] = np.linspace(0.0, 12.0, 13)
NODES = json.loads((ROOT / "results" / "sim30" / "sim36_signedband.json").read_text())["states"]


def aqueous(state):
    fam, k = state.split(":")
    if fam == "sigma":
        return float(k) <= 25.0
    if fam == "joint":
        return False
    e1, e2 = (float(x) for x in k.split("_")[:2])
    return e1 >= 1.0 or e2 >= 3 * abs(e1)


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


def lag06(s_ps, kap, kind):
    """Lag of the tissue pulse at 0.6 THz (rad): the shift plus, for the causal loss, its own phase."""
    extra = -np.angle(tp.loss_transfer(kap, kind)[i06]) if kind == "minphase" else 0.0
    return 2 * np.pi * F0 * s_ps + extra


def passing(vals, ps):
    ok = [v for v, p in zip(vals, ps) if p >= P_PASS]
    return (min(ok), max(ok)) if ok else None


out = {}
for kind in ("zero", "minphase"):
    say(f"=== reference with a {'zero-phase' if kind == 'zero' else 'causal (minimum-phase)'} smooth loss of the tissue pulse ===")
    out[kind] = {}
    for scan, (y, x0) in SCANS.items():
        p_nl, s_nl, _ = tp.fit(y, x0)                                  # Sims 47-49 reference, for comparison
        pr, s_ref, res = tp.fit(y, x0, loss=kind, a_max=1.0)
        ne, _ = n_eff(res)
        q0, s0, _ = tp.fit(y, x0, s=0.0, loss=kind, a_max=1.0)
        F0_, p0_ = f_test(s0, s_ref, ne, NPAR)
        gain = s0 - s_ref
        # 95 % interval of s by the same rule
        sg = np.arange(pr[2] - 0.12, pr[2] + 0.1201, 0.004)
        ps, kap_at = [], []
        prev = [pr[0], pr[1], pr[3], pr[4]]
        for sv in sg:
            q, ss, _ = tp.fit(y, x0, s=sv, loss=kind, p0=prev, a_max=1.0)
            ps.append(f_test(ss, s_ref, ne, NPAR)[1]); kap_at.append(q[-1])
        rng_s = passing(sg, ps)
        k_lo = kap_at[list(sg).index(rng_s[0])] if rng_s else np.nan
        lag_best, lag_lo = lag06(pr[2], pr[4], kind), (lag06(rng_s[0], k_lo, kind) if rng_s else np.nan)
        say(f"{scan}: reference misfit {s_ref:.3g} (without the loss term {s_nl:.3g}), kappa {pr[4]:.1f} dB/THz, "
            f"s {pr[2]:+.4f} ps, a {pr[0]:.3f}, N_eff {ne:.1f}; s = 0: p = {p0_:.4f}; 95 % interval of s "
            f"{rng_s[0]:+.3f} to {rng_s[1]:+.3f} ps; lag at 0.6 THz {lag_best:.3f} rad, smallest passing {lag_lo:.3f} rad")
        res_scan = dict(misfit_ref=s_ref, misfit_noloss=s_nl, kappa=float(pr[4]), s=float(pr[2]), a=float(pr[0]), n_eff=float(ne),
                        p_zero=p0_, s_range=list(map(float, rng_s)) if rng_s else None, lag_best=float(lag_best),
                        lag_lo=float(lag_lo), cands={}, nodes={})
        out[kind][scan] = res_scan
        if scan.startswith("30 s"):
            continue
        # slab candidates
        for name, grid in GRIDS.items():
            pv, fits = [], []
            prev = None
            for v in grid:
                xt = tp.apply_transfer(x0, slab_ratio(name, v))
                q, ss, _ = tp.fit(y, xt, s=0.0, loss=kind, a_max=1.0) if prev is None else \
                    min((tp.fit(y, xt, s=0.0, loss=kind, a_max=1.0, p0=prev), tp.fit(y, xt, s=0.0, loss=kind, a_max=1.0)),
                        key=lambda t: t[1])
                prev = list(q); pv.append(f_test(ss, s_ref, ne, NPAR)[1]); fits.append((ss, q))
            ib = int(np.argmin([f[0] for f in fits])); vb = grid[ib]
            rng_v = passing(grid, pv)
            r = slab_ratio(name, vb)
            rec = lambda rr: (s0 - tp.fit(y, tp.apply_transfer(x0, rr), s=0.0, loss=kind, a_max=1.0)[1]) / gain
            rec_all, rec_mag, rec_ph = rec(r), rec(np.abs(r)), rec(np.exp(1j * np.angle(r)))
            r06 = r[i06]
            res_scan["cands"][name] = dict(best=float(vb), p_best=float(pv[ib]), pass_range=list(map(float, rng_v)) if rng_v else None,
                                           a_best=float(fits[ib][1][0]), recovered=rec_all, from_loss=rec_mag, from_delay=rec_ph,
                                           loss06_db=float(-20 * np.log10(abs(r06))), delay06_rad=float(-np.angle(r06)))
            say(f"  {name:36s} best {vb:7.3g} (p = {pv[ib]:.3f}, a = {fits[ib][1][0]:.3f}); passes "
                + (f"{rng_v[0]:.3g} to {rng_v[1]:.3g}" if rng_v else "nowhere")
                + f"; recovers {rec_all:.0%} (its loss alone {rec_mag:.0%}, its delay alone {rec_ph:.0%})")
        # node states at every extreme, then the smallest scale at which any passes
        pn = []
        for st in NODES:
            q, ss, _ = tp.fit(y, tp.apply_transfer(x0, node_transfer(st, fr)), s=0.0, loss=kind, a_max=1.0)
            pn.append((f_test(ss, s_ref, ne, NPAR)[1], st["state"], ss))
        pn.sort(reverse=True)
        best_aq = max((t for t in pn if aqueous(t[1])), key=lambda t: t[0])
        res_scan["nodes"] = dict(best=pn[0][1], p_best=pn[0][0], n_pass=sum(t[0] >= P_PASS for t in pn),
                                 best_aqueous=best_aq[1], p_best_aqueous=best_aq[0])
        say(f"  node states, every extreme (delay with the path factor 2.2): best {pn[0][1]} p = {pn[0][0]:.3f}, "
            f"recovers {(s0 - pn[0][2]) / gain:.0%}; {res_scan['nodes']['n_pass']} of {len(pn)} pass; best aqueous "
            f"{best_aq[1]} p = {best_aq[0]:.3f}")
        top = [t[1] for t in pn[:3]]
        lam_min = None
        for lam in np.arange(0.1, 1.001, 0.05):
            for name_ in top:
                st = next(x for x in NODES if x["state"] == name_)
                ss = tp.fit(y, tp.apply_transfer(x0, node_transfer(st, fr, lam)), s=0.0, loss=kind, a_max=1.0)[1]
                if f_test(ss, s_ref, ne, NPAR)[1] >= P_PASS:
                    lam_min = (float(lam), name_); break
            if lam_min:
                break
        res_scan["nodes"]["lam_min"] = lam_min
        say(f"  node states: smallest fraction of every extreme at which one passes: "
            + (f"{lam_min[0]:.2f} ({lam_min[1]})" if lam_min else "none up to 1"))
        # with the delay scaled like the loss (3.1 at 2-3 THz), as in Sims 48-49
        st = next(x for x in NODES if x["state"] == pn[0][1])
        ss = tp.fit(y, tp.apply_transfer(x0, node_transfer(st, fr, delay_factor="loss")), s=0.0, loss=kind, a_max=1.0)[1]
        res_scan["nodes"]["p_best_loss_factor"] = f_test(ss, s_ref, ne, NPAR)[1]
        say(f"  same state with the delay scaled like the loss (3.1 at 2-3 THz): p = {res_scan['nodes']['p_best_loss_factor']:.3f}")
    say()
(OUT / "sim50.json").write_text(json.dumps(out, indent=1))
(OUT / "sim50.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
