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

from thznerve.physics.effmedium import tissue_eps
from thznerve.physics.materials import (EPS0, PROTEIN, WATER_KS1996, alpha_per_cm, axoplasm, myelin,
                                        water_liebe)
from thznerve.physics.tissue_change import LOSS_DB, delay_phase
from thznerve.physics.tmm import slab_rta

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim46"; OUT.mkdir(parents=True, exist_ok=True)
FIG = ROOT / "paper" / "figures"
lines = []


def say(s=""):
    print(s, flush=True); lines.append(s)


D_UM, P, N_PTFE, WF = 220.0, 0.6, 1.43, 0.77
F = np.array([0.6, 1.0, 1.4]) * 1e12
MEAS_LOSS = np.array([LOSS_DB[f] for f in F])     # unbiased |T|^2 (dB below empty), Sim 43 (c) band medians
DELAY_T = delay_phase(F)                          # tissue delay phase -arg T, Sim 43 (d)
W0 = WATER_KS1996

# Liebe's temperature change, scaled to the measured 4 -> 21 C absorption change at 1 THz
_dA = alpha_per_cm(water_liebe(21.0)(1e12), 1e12) - alpha_per_cm(water_liebe(4.0)(1e12), 1e12)
K_MILLON = 46.0 / _dA


class HeatedWater:
    """Base water plus k times Liebe's change of eps from 20 C to 20 + dT C."""

    def __init__(self, dT, k=1.0):
        self.dT, self.k = dT, k

    def __call__(self, f):
        return W0(f) + self.k * (water_liebe(20.0 + self.dT)(f) - water_liebe(20.0)(f))


def eps_tissue(f, water=W0, sig_ax=0.0, sig_w=0.0):
    ax = axoplasm(water)(f) - 1j * sig_ax / (2 * np.pi * f * EPS0)
    e, _ = tissue_eps(ax, myelin(water)(f), water(f), 0.7, P)
    return e - 1j * WF * sig_w / (2 * np.pi * f * EPS0)


def t_slab(f, e, d_um):
    """Transfer relative to the same thickness of air on PTFE: a thicker slice displaces air, as a
    cuvette window moving outwards would (extra phase k0 (n - 1) dd)."""
    return slab_rta(f, 1.0, [(complex(e), d_um * 1e-6)], N_PTFE)[4] / slab_rta(f, 1.0, [(1.0 + 0j, d_um * 1e-6)], N_PTFE)[4]


def change(fn):
    """(extra loss dB, extra delay rad) at F for a function f -> (eps, d_um)."""
    out = []
    for f in F:
        t0 = t_slab(f, eps_tissue(f), D_UM)
        t1 = t_slab(f, *fn(f))
        out.append((-20 * np.log10(abs(t1 / t0)), -np.angle(t1 / t0)))   # exp(+jwt): delay = -arg
    return np.array(out)


def _mix_delta(f, x):
    """Change of tissue eps when its water content changes by the fraction x (dry matter fixed)."""
    dry = PROTEIN(f)
    mix = lambda wf: (wf * W0(f) ** (1 / 3) + (1 - wf) * dry ** (1 / 3)) ** 3
    wf1 = WF * (1 + x) / (1 + WF * x)
    return mix(wf1) - mix(WF)


def water_frac(x, mode="thickness"):
    """Water content changed by the fraction x (x < 0: loss). mode: 'thickness' (added volume into
    thickness), 'eps' (permittivity only, thickness fixed), 'geom' (thickness only, eps fixed),
    'isotropic' (volume spread in 3D: thickness x cube root)."""
    def fn(f):
        de = 0.0 if mode == "geom" else _mix_delta(f, x)
        d = {"thickness": D_UM * (1 + WF * x), "geom": D_UM * (1 + WF * x), "eps": D_UM,
             "isotropic": D_UM * (1 + WF * x) ** (1 / 3)}[mode]
        return eps_tissue(f) + de, d
    return fn


MECH = {
    "heating, Liebe (dT, K)": (lambda v: (lambda f: (eps_tissue(f, HeatedWater(v)), D_UM)), (0.0, 80.0)),
    f"heating, Liebe x{K_MILLON:.2f} (dT, K)": (lambda v: (lambda f: (eps_tissue(f, HeatedWater(v, K_MILLON)), D_UM)), (0.0, 60.0)),
    "water content (fractional change)": (water_frac, (-0.5, 0.5)),
    "ions: all tissue water (S/m)": (lambda v: (lambda f: (eps_tissue(f, sig_w=v), D_UM)), (0.0, 500.0)),
    "axoplasm conductivity (S/m)": (lambda v: (lambda f: (eps_tissue(f, sig_ax=v), D_UM)), (0.0, 2000.0)),
    "thicker path (fraction)": (lambda v: (lambda f: (eps_tissue(f), D_UM * (1 + v))), (0.0, 0.5)),
}

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
(OUT / "sim46.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

# Figure 9: extra loss vs extra delay of the slice at 0.6 and 1.0 THz
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
node = json.loads((ROOT / "results" / "sim30" / "sim36_signed06.json").read_text())["states"]
hk = [n for n in MECH if n.startswith("heating, Liebe x")][0]
for k, axi in enumerate(ax):
    axi.plot([REQ[1.09][k, 0], REQ[1.18][k, 0]], [REQ[1.09][k, 1], REQ[1.18][k, 1]],
             color="0.6", lw=7, solid_capstyle="butt", label="required: T^q, q = 1.09-1.18")
    qq = np.linspace(0, 0.2, 21)
    axi.plot(qq * MEAS_LOSS[k], qq * DELAY_T[k], color="0.5", lw=0.8, ls=":")
    for name, col, key, vals in (("heating, Liebe (0-40 K)", "C3", "heating, Liebe (dT, K)", np.linspace(0, 40, 17)),
                                 (f"heating, Liebe x{K_MILLON:.1f} (0-20 K)", "C4", hk, np.linspace(0, 20, 17)),
                                 ("water gain (0-16 %)", "C0", "water content (fractional change)", np.linspace(0, 0.16, 17)),
                                 ("ions in tissue water (0-18 S/m)", "C2", "ions: all tissue water (S/m)", np.linspace(0, 18, 10)),
                                 ("axoplasm conduction (0-46 S/m)", "C1", "axoplasm conductivity (S/m)", np.linspace(0, 46, 10))):
        c = np.array([change(MECH[key][0](v))[k] for v in vals])
        axi.plot(c[:, 0], c[:, 1], "o-", color=col, ms=2.5, lw=1, label=name)
    if k == 0:
        axi.plot([s["loss_db"] for s in node], [s["delay_rad"] for s in node], "k.", ms=3,
                 label="any node state, every extreme (116)")
    axi.axhline(0, color="k", lw=0.5); axi.axvline(0, color="k", lw=0.5)
    axi.set_xlabel("extra loss of the slice (dB)"); axi.set_title(f"{F[k] / 1e12:.1f} THz", fontsize=10)
ax[0].set_ylabel("extra delay of the slice (rad)")
ax[0].legend(fontsize=7, loc="upper left")
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"fig_tissue_change.{ext}", dpi=200)
