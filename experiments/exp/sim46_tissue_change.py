"""Sim 46 - Review round 12: which change of the tissue itself gives the T -> T^q the waveform fits need?

The two-term fits (Sim 42/43) need the unbiased tissue transfer T raised to q = 1.09 (H1, 4 min) to
1.18 (H3): at each frequency an extra loss (q-1)*|T|^2_dB AND an extra delay (q-1)*arg T. Each
candidate change of the slice is applied to the effective-medium white-matter slab (E-par, packing
0.6, on PTFE, 220 um optical thickness from Sim 45), and its change of loss and delay is compared
with the requirement at 0.6, 1.0 and 1.4 THz:
  heating (Liebe water at 20 + dT C in every water-bearing compartment); water uptake or loss (water
  fraction 0.77, Looyenga with dry matter, thickness following the water volume as in Sim 34);
  extra conductivity of all tissue water (ions); extra conductivity of the axoplasm only (the
  conducting-fibre model; quasi-static mixing, within 3 % of full wave at <= 1 THz up to 400 S/m,
  Sec. 5.2); and a geometrically thicker path (reproduces T^q by construction).
For each, the parameter that gives the required 0.6 THz loss is found and the rest is predicted.
Every transfer is taken relative to air of the same thickness, so a thicker slice displaces air.
Writes results/sim46/sim46.txt.
"""
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

from thznerve.physics.effmedium import tissue_eps
from thznerve.physics.materials import EPS0, PROTEIN, Looyenga, axoplasm, myelin, water_liebe
from thznerve.physics.tmm import slab_rta

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim46"; OUT.mkdir(parents=True, exist_ok=True)
lines = []


def say(s=""):
    print(s, flush=True); lines.append(s)


D_UM, P, N_PTFE, WF = 220.0, 0.6, 1.43, 0.77
F = np.array([0.6, 1.0, 1.4]) * 1e12
MEAS_LOSS = np.array([12.1, 14.9, 19.5])          # unbiased |T|^2 (dB below empty), Sim 43 (c) band medians
ARG_T = 2 * np.pi * F * 0.611e-12 - 0.22          # tissue phase delay, Sim 43 (d)
W0 = water_liebe(20.0)


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


def water_frac(x):
    """Fractional change x of the water content (x < 0: loss); thickness follows the water volume."""
    def fn(f):
        dry = Looyenga(W0, PROTEIN, 0.0)(f)
        mix = lambda wf: (wf * W0(f) ** (1 / 3) + (1 - wf) * dry ** (1 / 3)) ** 3
        wf1 = WF * (1 + x) / (1 + WF * x)
        return eps_tissue(f) + mix(wf1) - mix(WF), D_UM * (1 + WF * x)
    return fn


MECH = {
    "heating (dT, K)": (lambda v: (lambda f: (eps_tissue(f, water_liebe(20.0 + v)), D_UM)), (0.0, 60.0)),
    "water content (fractional change)": (water_frac, (-0.5, 0.5)),
    "ions: all tissue water (S/m)": (lambda v: (lambda f: (eps_tissue(f, sig_w=v), D_UM)), (0.0, 500.0)),
    "axoplasm conductivity (S/m)": (lambda v: (lambda f: (eps_tissue(f, sig_ax=v), D_UM)), (0.0, 2000.0)),
    "thicker path (fraction)": (lambda v: (lambda f: (eps_tissue(f), D_UM * (1 + v))), (0.0, 0.5)),
}

for q in (1.09, 1.18):
    req = np.column_stack([(q - 1) * MEAS_LOSS, (q - 1) * ARG_T])
    say(f"q = {q}: required extra loss / delay at 0.6, 1.0, 1.4 THz: "
        + "; ".join(f"{l:.2f} dB, {p:+.3f} rad" for l, p in req))
    for name, (mk, (lo, hi)) in MECH.items():
        g = lambda v: change(mk(v))[0, 0] - req[0, 0]
        grid = np.linspace(lo, hi, 201); vals = [g(v) for v in grid]
        roots = [brentq(g, grid[i], grid[i + 1]) for i in range(200) if vals[i] * vals[i + 1] < 0]
        if not roots:
            say(f"  {name:36s}: never gives the required 0.6 THz loss over {lo}..{hi} "
                f"(range {min(vals) + req[0, 0]:+.2f} to {max(vals) + req[0, 0]:+.2f} dB)"); continue
        for v in roots:
            c = change(mk(v))
            say(f"  {name:36s} = {v:+.3g}: " + "; ".join(f"{l:+.2f} dB ({l / r[0]:.2f}x), {p:+.3f} rad ({p / r[1]:+.2f}x)"
                                                      for (l, p), r in zip(c, req)))
    # best pair of mechanisms (least squares on loss in dB and delay in rad*10, all three frequencies)
    best = []
    names = list(MECH)
    for i in range(len(names) - 1):            # the thicker path is excluded from pairs (it is T^q itself)
        for j in range(i + 1, len(names) - 1):
            (mi, (li, hi_)), (mj, (lj, hj)) = MECH[names[i]], MECH[names[j]]
            def res(v):
                fi, fj = mi(v[0]), mj(v[1])
                c = change(lambda f: (fi(f)[0] + fj(f)[0] - eps_tissue(f), fi(f)[1] * fj(f)[1] / D_UM))
                return np.concatenate([(c[:, 0] - req[:, 0]) / req[:, 0], (c[:, 1] - req[:, 1]) / req[:, 1]])
            s = least_squares(res, [0.25 * (li + hi_), 0.25 * (lj + hj)], bounds=([li, lj], [hi_, hj]))
            best.append((np.sqrt(np.mean(s.fun ** 2)), names[i], names[j], s.x))
    best.sort(key=lambda b: b[0])
    for r, a, b, x in best[:3]:
        say(f"  pair {a} = {x[0]:+.3g} with {b} = {x[1]:+.3g}: rms relative misfit {r:.0%} over loss and delay at 3 frequencies")
    say()
(OUT / "sim46.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

# Figure 9: extra loss vs extra delay of the slice at 0.6 and 1.0 THz
import json
import matplotlib.pyplot as plt

FIG = ROOT / "paper" / "figures"
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
node = json.loads((ROOT / "results" / "sim30" / "sim36_signed06.json").read_text())["states"]
for k, axi in enumerate(ax):
    for q0, q1 in ((1.09, 1.18),):
        axi.plot([(q0 - 1) * MEAS_LOSS[k], (q1 - 1) * MEAS_LOSS[k]], [(q0 - 1) * ARG_T[k], (q1 - 1) * ARG_T[k]],
                 color="0.6", lw=7, solid_capstyle="butt", label="required: T^q, q = 1.09-1.18")
    qq = np.linspace(0, 0.2, 21)
    axi.plot((qq) * MEAS_LOSS[k], qq * ARG_T[k], color="0.5", lw=0.8, ls=":")
    for name, col, vals in (("heating (0-30 K)", "C3", np.linspace(0, 30, 16)),
                            ("water gain (0-16 %)", "C0", np.linspace(0, 0.16, 17)),
                            ("ions in tissue water (0-18 S/m)", "C2", np.linspace(0, 18, 10)),
                            ("axoplasm conduction (0-46 S/m)", "C1", np.linspace(0, 46, 10))):
        mk = {"heating": MECH["heating (dT, K)"][0], "water": MECH["water content (fractional change)"][0],
              "ions": MECH["ions: all tissue water (S/m)"][0], "axoplasm": MECH["axoplasm conductivity (S/m)"][0]}[name.split()[0]]
        c = np.array([change(mk(v))[k] for v in vals])
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
