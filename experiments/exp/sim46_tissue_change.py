"""Sim 46 - Review rounds 12-13: which change of the tissue itself gives the T -> T^q the waveform fits need?

The two-term fits (Sim 42/43) need the unbiased tissue transfer T raised to q = 1.09 (H1, 4 min) to
1.18 (H3): at each frequency an extra loss (q-1)*|T|^2_dB AND an extra delay (q-1)*(-arg T), with
arg T = -2 pi f 0.611 ps - 0.22 rad from Sim 43 (numpy FFT, delay = -arg). Each candidate change of
the slice is applied to the effective-medium white-matter slab (E-par, packing 0.6, on PTFE, 220 um
optical thickness from Sim 45, Kindt-Schmuttenmaer water as in the node bounds), and its change of
loss and delay is compared with the requirement at 0.6, 1.0 and 1.4 THz:
  heating: the temperature change of Liebe et al. (1991) water added to the base water, as given and
    scaled so that its absorption change from 4 to 21 C at 1 THz equals the measured 46 cm^-1
    (Millon et al. 2022; literature 50 +- 9); Liebe is fitted below 1 THz and 30 C, so neither is
    trusted beyond that, and the two bracket the model dependence;
  water uptake or loss (water fraction 0.77, Looyenga with dry matter; the added volume goes into
    thickness only, as in Sim 34), split into its thickness part and its permittivity part, and
    with isotropic swelling (thickness grows as the cube root of the volume);
  extra conductivity of all tissue water (ions); extra conductivity of the axoplasm only (the
    conducting-fibre model; quasi-static mixing, within 3 % of full wave at <= 1 THz up to 400 S/m,
    Sec. 5.2); and a geometrically thicker path (T^q up to the slab model's own error).
For each, the parameter that gives the required 0.6 THz loss is found and the rest is predicted.
Every transfer is taken relative to air of the same thickness, so a thicker slice displaces air.
Writes results/sim46/sim46.txt and paper/figures/fig_tissue_change.{png,pdf}.
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq, least_squares

from thznerve.physics.materials import alpha_per_cm, water_liebe
from thznerve.physics.slice_change import (D_UM, WF, eps_tissue, k_millon, mechanisms, t_slab,
                                          water_frac)
from thznerve.physics.tissue_change import LOSS_DB, delay_phase, min_required_delay

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim46"; OUT.mkdir(parents=True, exist_ok=True)
FIG = ROOT / "paper" / "figures"
lines = []


def say(s=""):
    print(s, flush=True); lines.append(s)


F = np.array([0.6, 1.0, 1.4]) * 1e12
MEAS_LOSS = np.array([LOSS_DB[f] for f in F])     # unbiased |T|^2 (dB below empty), Sim 43 (c) band medians
DELAY_T = delay_phase(F)                          # tissue delay phase -arg T, Sim 43 (d)

# Liebe's temperature change, scaled to the measured 4 -> 21 C absorption change at 1 THz
_dA = alpha_per_cm(water_liebe(21.0)(1e12), 1e12) - alpha_per_cm(water_liebe(4.0)(1e12), 1e12)
K_MILLON = k_millon()


def change(fn):
    """(extra loss dB, extra delay rad) at F for a function f -> (eps, d_um)."""
    out = []
    for f in F:
        t0 = t_slab(f, eps_tissue(f), D_UM)
        t1 = t_slab(f, *fn(f))
        out.append((-20 * np.log10(abs(t1 / t0)), -np.angle(t1 / t0)))   # exp(+jwt): delay = -arg
    return np.array(out)


MECH = mechanisms(K_MILLON)


def _base_delay():
    """Unwrapped delay phase of the base slab at F (from 0.02 THz up)."""
    fr = np.linspace(0.02e12, F.max(), 400)
    ph = -np.unwrap(np.angle([t_slab(f, eps_tissue(f), D_UM) for f in fr]))
    return np.interp(F, fr, ph)


say(f"Base slab {D_UM:.0f} um, KS water: |T|^2 {', '.join(f'{-20 * np.log10(abs(t_slab(f, eps_tissue(f), D_UM))):.1f}' for f in F)} dB, "
    f"delay phase {', '.join(f'{v:.2f}' for v in _base_delay())} rad at 0.6/1.0/1.4 THz "
    f"(measured {', '.join(f'{v:.1f}' for v in MEAS_LOSS)} dB; {', '.join(f'{v:.2f}' for v in DELAY_T)} rad)")
say(f"Liebe water, 4 -> 21 C absorption change at 1 THz: {_dA:.1f} /cm (measured 46 +- 1.5, Millon et al. 2022): "
    f"scaled version x{K_MILLON:.2f}")
say()
REQ = {}
for q in (1.09, 1.18):
    req = np.column_stack([(q - 1) * MEAS_LOSS, (q - 1) * DELAY_T]); REQ[q] = req
    say(f"q = {q}: required extra loss / delay at 0.6, 1.0, 1.4 THz: "
        + "; ".join(f"{l:.2f} dB, {p:+.3f} rad" for l, p in req))
    roots_at = {}
    for name, (mk, (lo, hi)) in MECH.items():
        g = lambda v: change(mk(v))[0, 0] - req[0, 0]
        grid = np.linspace(lo, hi, 201); vals = [g(v) for v in grid]
        roots = [brentq(g, grid[i], grid[i + 1]) for i in range(200) if vals[i] * vals[i + 1] < 0]
        if not roots:
            say(f"  {name:36s}: never gives the required 0.6 THz loss over {lo}..{hi} "
                f"(range {min(vals) + req[0, 0]:+.2f} to {max(vals) + req[0, 0]:+.2f} dB)"); continue
        roots_at[name] = roots[0]
        for v in roots:
            c = change(mk(v))
            rr = np.concatenate([(c[:, 0] - req[:, 0]) / req[:, 0], (c[:, 1] - req[:, 1]) / req[:, 1]])
            say(f"  {name:36s} = {v:+.3g}: " + "; ".join(f"{l:+.2f} dB ({l / r[0]:.2f}x), {p:+.3f} rad ({p / r[1]:+.2f}x)"
                                                      for (l, p), r in zip(c, req)) + f"; rms misfit {np.sqrt(np.mean(rr ** 2)):.0%}")
    x = roots_at["water content (fractional change)"]
    for mode, lab in (("geom", "its thickness change only"), ("eps", "its permittivity change only"),
                      ("isotropic", "isotropic swelling instead")):
        c = change(water_frac(x, mode))
        say(f"    water gain {x:+.3f}, {lab:30s}: " + "; ".join(f"{l / r[0]:.2f}x loss, {p / r[1]:+.2f}x delay" for (l, p), r in zip(c, req)))
    say(f"    water gain {x:+.3f}: added thickness {D_UM * WF * x:.1f} um, water {D_UM * WF * x * 1e-3 * 400:.1f} mg over 20 x 20 mm")
    best = []
    names = [n for n in MECH if not n.startswith("thicker")]   # the thicker path is excluded from pairs
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            if names[i].startswith("heating") and names[j].startswith("heating"):
                continue
            (mi, (li, hi_)), (mj, (lj, hj)) = MECH[names[i]], MECH[names[j]]

            def res(v):
                fi, fj = mi(v[0]), mj(v[1])
                c = change(lambda f: (fi(f)[0] + fj(f)[0] - eps_tissue(f), fi(f)[1] * fj(f)[1] / D_UM))
                return np.concatenate([(c[:, 0] - req[:, 0]) / req[:, 0], (c[:, 1] - req[:, 1]) / req[:, 1]])
            s = least_squares(res, [0.05 * (li + hi_), 0.05 * (lj + hj)], bounds=([li, lj], [hi_, hj]))
            best.append((np.sqrt(np.mean(s.fun ** 2)), names[i], names[j], s.x))
    best.sort(key=lambda b: b[0])
    for r, a, b, xv in best[:4]:
        say(f"  pair {a} = {xv[0]:+.3g} with {b} = {xv[1]:+.3g}: rms relative misfit {r:.0%} over loss and delay at 3 frequencies")
    say()
# (b) The fits leave the amplitude A of the tissue term free, so only the delay is fixed; extra loss L
# at 0.6 THz must be paid for by a larger A. H1 4 min (Sim 42): A = 0.778, b = 0.164, q = 1.091.
# a <= 1 bounds L <= (q-1)|T|^2 + 20 log10(1/A); a <= 1 - b (the component takes area from the tissue path) is tighter.
A_H1, B_H1, Q_H1 = 0.778, 0.164, 1.091
req = np.column_stack(required := [(Q_H1 - 1) * MEAS_LOSS, (Q_H1 - 1) * DELAY_T])
L_A1 = req[0, 0] + 20 * np.log10(1 / A_H1); L_AB = req[0, 0] + 20 * np.log10((1 - B_H1) / A_H1)
say(f"(b) Delay-matched, H1 4 min (q = {Q_H1}, A = {A_H1}, b = {B_H1}): required delay "
    + ", ".join(f"{v:.3f}" for v in req[:, 1]) + f" rad at 0.6/1.0/1.4 THz; extra loss at 0.6 THz allowed "
    f"<= {L_A1:.2f} dB (a <= 1) or <= {L_AB:.2f} dB (a <= 1 - b); loss rise 0.6->1.0/1.4 THz needed "
    f"{req[1, 0] - req[0, 0]:+.2f}/{req[2, 0] - req[0, 0]:+.2f} dB")
DM = dict(MECH)
for name, (mk, (lo, hi)) in DM.items():
    g = lambda v: change(mk(v))[0, 1] - req[0, 1]
    grid = np.linspace(lo, hi, 201); vals = [g(v) for v in grid]
    roots = [brentq(g, grid[i], grid[i + 1]) for i in range(200) if vals[i] * vals[i + 1] < 0]
    if not roots:
        say(f"  {name:36s}: never gives the 0.6 THz delay over {lo}..{hi} (max {max(vals) + req[0, 1]:+.3f} rad)"); continue
    v = roots[0]; c = change(mk(v))
    ok = "allowed" if c[0, 0] <= L_A1 else "NOT allowed (a > 1)"
    say(f"  {name:36s} = {v:+.3g}: loss at 0.6 THz {c[0, 0]:+.2f} dB ({ok}); delay at 1.0/1.4 THz "
        f"{c[1, 1] / req[1, 1]:.2f}x/{c[2, 1] / req[2, 1]:.2f}x; loss rise {c[1, 0] - c[0, 0]:+.2f}/{c[2, 0] - c[0, 0]:+.2f} dB")
s_ps = req[0, 1] / (2 * np.pi * F[0]) * 1e12
say(f"  {'time shift of the tissue term (drift)':36s} = {s_ps:+.3f} ps: loss 0 dB (allowed); delay at 1.0/1.4 THz "
    f"{2 * np.pi * F[1] * s_ps * 1e-12 / req[1, 1]:.2f}x/{2 * np.pi * F[2] * s_ps * 1e-12 / req[2, 1]:.2f}x; loss rise 0 dB")
say()
(OUT / "sim46.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

# Figure 9: extra loss vs extra delay of the slice at 0.6 and 1.0 THz
D_MIN = min_required_delay(F)
say(f"Lowest extra delay the traces allow (Sim 47, one-sided 95 %, 2-5 min): {', '.join(f'{v:.3f}' for v in D_MIN)} rad at 0.6/1.0/1.4 THz")
for name, (mk, (lo, hi)) in MECH.items():
    g = lambda v: change(mk(v))[0, 1] - D_MIN[0]
    grid = np.linspace(lo, hi, 201); vals = [g(v) for v in grid]
    roots = [brentq(g, grid[i], grid[i + 1]) for i in range(200) if vals[i] * vals[i + 1] < 0]
    if roots:
        c = change(mk(roots[0]))
        say(f"  {name:36s} = {roots[0]:+.3g} gives that delay at 0.6 THz with {c[0, 0]:+.2f} dB "
            f"({'allowed' if c[0, 0] <= L_A1 else 'NOT allowed'}, <= {L_A1:.2f})")
# per scan: each scan's own 95 % lower bound on the delay (Sim 47) against its own loss allowance
# (H1 4 min: A = 0.778, q = 1.091; H3 5 min: A = 0.689, q = 1.109 from Sim 39 over 0.15-1.70 THz; H2 has no spectral fit)
A_H3, Q_H3 = 0.689, 1.109
L_A3 = (Q_H3 - 1) * MEAS_LOSS[0] + 20 * np.log10(1 / A_H3)
for scan, allow in (("2 min", None), ("4 min", L_A1), ("5 min", L_A3)):
    dmin = min_required_delay(F[:1], (scan,))[0]
    say(f"Per scan, {scan}: 95 % lower bound {dmin:.3f} rad at 0.6 THz; allowed extra loss "
        + (f"{allow:.2f} dB" if allow is not None else "not fitted"))
    for name in ("ions: all tissue water (S/m)", "axoplasm conductivity (S/m)", "heating, Liebe (dT, K)",
                 [n for n in MECH if n.startswith("heating, Liebe x")][0], "thicker path (fraction)"):
        mk, (lo, hi) = MECH[name]
        g = lambda v: change(mk(v))[0, 1] - dmin
        grid = np.linspace(lo, hi, 201); vals = [g(v) for v in grid]
        roots = [brentq(g, grid[i], grid[i + 1]) for i in range(200) if vals[i] * vals[i + 1] < 0]
        if roots:
            c = change(mk(roots[0]))
            verdict = "" if allow is None else (" allowed" if c[0, 0] <= allow else " NOT allowed")
            say(f"  {name:36s} = {roots[0]:+.3g} with {c[0, 0]:+.2f} dB{verdict}")
(OUT / "sim46.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
fig, ax = plt.subplots(1, 2, figsize=(11, 5.2))
node = json.loads((ROOT / "results" / "sim30" / "sim36_signed06.json").read_text())["states"]
hk = [n for n in MECH if n.startswith("heating, Liebe x")][0]
for k, axi in enumerate(ax):
    ceil = req[k, 0] + 20 * np.log10(1 / A_H1)
    axi.fill_between([0, ceil], D_MIN[k], REQ[1.18][k, 1], color="0.85", lw=0,
                     label="what the fits allow: extra delay (95 % lower bound to q = 1.18), extra loss up to a = 1 (H1)")
    axi.plot([0], [2 * np.pi * F[k] * s_ps * 1e-12], "kx", ms=7, mew=1.5, label=f"time shift of the tissue term ({s_ps:.3f} ps)")
    for name, col, key, vals in (("heating, Liebe (0-80 K)", "C3", "heating, Liebe (dT, K)", np.linspace(0, 80, 21)),
                                 (f"heating, Liebe x{K_MILLON:.1f} (0-20 K)", "C4", hk, np.linspace(0, 20, 17)),
                                 ("water gain into thickness (0-16 %)", "C0", "water content (fractional change)", np.linspace(0, 0.16, 17)),
                                 ("same tissue thicker (0-16 %)", "C5", "thicker path (fraction)", np.linspace(0, 0.16, 17)),
                                 ("ions in tissue water (0-70 S/m)", "C2", "ions: all tissue water (S/m)", np.linspace(0, 70, 15)),
                                 ("axoplasm conduction (0-180 S/m)", "C1", "axoplasm conductivity (S/m)", np.linspace(0, 180, 15))):
        c = np.array([change(MECH[key][0](v))[k] for v in vals])
        axi.plot(c[:, 0], c[:, 1], "o-", color=col, ms=2.5, lw=1, label=name)
    if k == 0:
        axi.plot([s["loss_db"] for s in node], [s["delay_rad"] for s in node], "k.", ms=3,
                 label=f"node states, every extreme ({len(node)}; 0.6 THz only)")
    axi.axhline(0, color="k", lw=0.5); axi.axvline(0, color="k", lw=0.5)
    axi.set_xlabel("extra loss of the slice (dB)"); axi.set_title(f"{F[k] / 1e12:.1f} THz", fontsize=10)
ax[0].set_ylabel("extra delay of the slice (rad)")
h, l = ax[0].get_legend_handles_labels()
fig.legend(h, l, fontsize=7.5, loc="lower center", ncol=3, frameon=False)
fig.tight_layout(rect=(0, 0.17, 1, 1))
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"fig_tissue_change.{ext}", dpi=200)
