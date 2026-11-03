"""Sim 54 - Can scattering supply the bias-induced change of the tissue pulse?

Sim 53: the change is a causal loss over the 30 s floor rising 2.1-7.5 dB/THz at 2-5 min, and no passive
change tested (thicker path, warming, conduction, node states) gives a loss that steep. Scattering by
structure on the wavelength scale does. Candidate: spheres of radius r at volume fraction v in the tissue,
as a coherent (Foldy / van de Hulst) effective medium with the Mie forward amplitude
(`thznerve.physics.scattering`):
  - gas bubbles (eps = 1), electrolysis gas, either displacing tissue (slab thickness fixed) or added to it
    (the slice swells: thickness D / (1 - v));
  - lipid-like fragments (eps 2.4 - 0.1j) in the tissue water, a slice disrupted by the current.
Radii 15-120 um, v 0-0.2 (Foldy is a dilute-suspension theory; v <= 0.2 is its edge). Each candidate is
tested as in Sim 53 against the six references of a zero-phase 30 s floor times a causal increment
(`waveform_fit.supply_test`). Writes results/sim54/sim54.txt and results/sim54/sim54.json.
"""
import json
from pathlib import Path

import numpy as np

from thznerve.physics.scattering import foldy_eps
from thznerve.physics.slice_change import D_UM, eps_tissue, ratio
from thznerve.waveform_fit import TwoPulse, f_test, load_scans, n_eff, supply_test

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim54"; OUT.mkdir(parents=True, exist_ok=True)
lines = []


def say(s=""):
    print(s, flush=True); lines.append(s)


tt, K, SCANS = load_scans(ROOT)
F0, P_PASS = 0.6, 0.05
GRID = np.linspace(0.0, 0.2, 11)
KINDS = {"gas, slice swells": (lambda f: 1.0, True),
         "gas, displacing tissue": (lambda f: 1.0, False),
         "lipid-like fragments": (lambda f: 2.4 - 0.1j, False)}
RADII = (15, 30, 50, 80, 120)
tp0 = TwoPulse(tt, K)
KAPPA0 = float(tp0.fit(*SCANS["30 s (H1)"], loss="zero", a_max=1.0)[0][4])
fr = tp0.fr; band = (fr > 0) & (fr <= 4.0); i06 = int(np.argmin(np.abs(fr - F0)))
FLOOR = tp0.loss_transfer(KAPPA0, "zero")
sl = (fr >= 0.3) & (fr <= 1.5)


def slope_db(r):
    return float(np.polyfit(fr[sl], -20 * np.log10(np.abs(r[sl])), 1)[0])


cache = {}


def scat_ratio(kind, rad, v):
    key = (kind, rad, round(float(v), 9))
    if key not in cache:
        r = np.ones(len(fr), complex)
        if v != 0:
            es, swell = KINDS[kind]
            d = D_UM / (1 - v) if swell else D_UM
            r[band] = ratio(fr[band] * 1e12, lambda f: (foldy_eps(f, eps_tissue(f), es(f), rad, v), d))
            r[fr > 4.0] = r[band][-1]
        cache[key] = r
    return cache[key]


say(f"Floor (30 s, zero phase): {KAPPA0:.2f} dB/THz")
say("Forward: candidate at v = 0.1 and 0.2 (loss and delay at 0.6 THz, loss slope 0.3-1.5 THz)")
for kind in KINDS:
    for rad in RADII:
        for v in (0.1, 0.2):
            r = scat_ratio(kind, rad, v)
            say(f"  {kind:24s} r {rad:3d} um v {v:.1f}: {-20 * np.log10(abs(r[i06])):+.2f} dB, "
                f"{-np.angle(r[i06]):+.3f} rad, slope {slope_db(r):+.2f} dB/THz")
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
            q0, s00, _ = fit(y, xf, s=0.0)
            p_none = f_test(s00, s_ref, ne, 5, 2)[1]
            inc = tp.loss_transfer(pr[4], "minphase")
            rs = dict(kappa=float(pr[4]), s=float(pr[2]), misfit=s_ref, n_eff=float(ne), p_no_change=p_none,
                      inc_slope=slope_db(inc), cands={})
            out[vn][scan] = rs
            say(f"{scan}: increment {pr[4]:.2f} dB/THz, s {pr[2]:+.4f} ps, misfit {s_ref:.3g}, N_eff {ne:.1f}; "
                f"no change at all: p = {p_none:.4f}")
            for kind in KINDS:
                for rad in RADII:
                    name = f"{kind}, r {rad} um"
                    tr = lambda v, kind=kind, rad=rad: scat_ratio(kind, rad, v)
                    d = supply_test(tp, fit, y, xf, tr, GRID, q0, pr, s_ref, ne, P_PASS)
                    if p_none >= P_PASS:
                        d["verdict"] = "n/a"
                    r = tr(d["best"])
                    d.update(loss06_db=float(-20 * np.log10(abs(r[i06]))), delay06_rad=float(-np.angle(r[i06])),
                             slope=slope_db(r))
                    rs["cands"][name] = d
                    pr_ = d["pass_range"]
                    say(f"  {name:34s} {d['verdict']:4s} p = {d['p']:.3f} ({d['p_own_neff']:.3f}); best v {d['best']:.2f}: "
                        f"{d['loss06_db']:+.2f} dB, {d['delay06_rad']:+.3f} rad at 0.6 THz, slope {d['slope']:+.2f} dB/THz; "
                        + (f"v passing {pr_[0]:.2f} to {pr_[1]:.2f}" if pr_ else "no v passes"))
        say()
say("=== Summary (2, 4, 5 min; 6 references) ===")
names = sorted({c for v in out if v != "kappa0" for sc in out[v].values() for c in sc["cands"]})
summ = {}
for c in names:
    vs = [out[v][sc]["cands"][c]["verdict"] for v in out if v != "kappa0" for sc in out[v] if not sc.startswith("30 s")]
    summ[c] = {k: vs.count(k) for k in ("pass", "edge", "fail", "n/a")}
    say(f"  {c:34s} " + ", ".join(f"{k} {summ[c][k]}" for k in ("pass", "edge", "fail", "n/a")) + f" of {len(vs)}")
out["summary"] = summ
(OUT / "sim54.json").write_text(json.dumps(out, indent=1))
(OUT / "sim54.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
