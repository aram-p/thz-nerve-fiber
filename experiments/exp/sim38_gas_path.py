"""Sim 38 - Can a gas path through the slice (electrolysis bubbles, drying) reproduce the reported
voltage-induced spectrum, including the drifting dip, with physical tissue optics?

Sim 29 used the authors' own index (n = 3.6-3.9 at d = 100 um), which no water-based tissue can
have. Here the tissue is the effective medium of Sims 34/36 (E || fibres, n ~ 2.0-2.2). Two numbers
in the reports point to an optical thickness larger than the stated 100 um:
  - their phase delay (n_rep - 1) d with n_rep = 3.6-3.9 equals (n_phys - 1) d_eff for d_eff ~ 2.4x;
  - their loss, 30 dB at 2 THz per "0.1 mm", is 1.9-2.6x pure water.
Model: a fraction x of the beam crosses a gas column of thickness t inside the slice of thickness d
(t <= d); the rest crosses tissue. Relative to the no-voltage spectrum,
    T_on / T_off = (1 - x) + x exp(+j k0 (n_t - 1) t)      (n_t complex, exp(+jwt): loss in the gas path is lower)
Fresnel factors at the PTFE windows are ignored (they change the depth, not the notch frequency).
Two parameters (x, t) are FITTED to points read off the published Fig. 1b of the IRPhE 2024 report
(no voltage vs 5 min at 100 V, E || fibres); +-2 dB reading uncertainty. Stated as a fit.
Then, without refitting: where would the dip sit at earlier times (smaller t), and the second dip?

Writes results/sim38/sim38.txt and paper/figures/fig_gas_path.png
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

from thznerve.physics.effmedium import tissue_eps
from thznerve.physics.materials import WATER_BROADBAND, WATER_KS1996, axoplasm, myelin, refractive_index

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim38"; OUT.mkdir(parents=True, exist_ok=True)
FIG = ROOT / "paper" / "figures"
C0 = 299792458.0
lines = []


def say(s=""):
    print(s, flush=True); lines.append(s)


def n_tissue(f, water=WATER_KS1996):
    e_par, _ = tissue_eps(axoplasm(water)(f), myelin(water)(f), water(f), 0.7, 0.6)
    return refractive_index(e_par)          # n' - j n''


# change of power spectrum (dB), 5 min at 100 V minus no voltage, read off Fig. 1b (paper 3)
F_DIG = np.array([0.3, 0.45, 0.6, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0])
DB_DIG = np.array([-2.0, -6.0, -11.0, -2.0, 2.0, 4.0, 8.0, 12.0, 15.0, 17.0])


def ratio_db(f_thz, x, t_um, water=WATER_KS1996):
    f = np.asarray(f_thz) * 1e12
    n = n_tissue(f, water)
    k0 = 2 * np.pi * f / C0
    g = (1 - x) + x * np.exp(1j * k0 * (n.real - 1) * t_um * 1e-6) * np.exp(k0 * (-n.imag) * t_um * 1e-6)
    return 20 * np.log10(np.abs(g))


say("Effective optical thickness implied by the reports (E || fibres tissue model)")
for fi, n_rep in ((0.3, 3.9), (1.0, 3.85), (2.0, 3.6)):
    n = n_tissue(fi * 1e12)
    say(f"  {fi} THz: reported n {n_rep} at d = 100 um -> delay equals physical n {n.real:.2f} at d = "
        f"{(n_rep - 1) / (n.real - 1) * 100:.0f} um")
for w, name in ((WATER_KS1996, "Kindt-Schmuttenmaer"), (WATER_BROADBAND, "broadband")):
    n = refractive_index(w(2e12)); a = 2 * 2 * np.pi * 2e12 / C0 * (-n.imag)
    say(f"  2 THz loss 30 dB vs pure water ({name}) {10/np.log(10)*a*100e-6:.1f} dB per 100 um -> "
        f"{30/(10/np.log(10)*a*100e-6)*100:.0f} um of water")

fits = {}
for w, name in ((WATER_KS1996, "KS1996"), (WATER_BROADBAND, "broadband")):
    res = least_squares(lambda q: (ratio_db(F_DIG, q[0], q[1], w) - DB_DIG) / 2.0, x0=[0.15, 200.0],
                        bounds=([0.0, 10.0], [1.0, 600.0]))
    x, t = res.x
    rms = np.sqrt(np.mean((ratio_db(F_DIG, x, t, w) - DB_DIG) ** 2))
    fits[name] = (x, t)
    say(f"\nFit ({name} water): beam fraction through gas x = {x:.3f}, gas column t = {t:.0f} um, "
        f"rms misfit {rms:.1f} dB over {len(F_DIG)} points (reading uncertainty ~2 dB)")
    ff = np.linspace(0.2, 3.0, 2801)
    y = ratio_db(ff, x, t, w)
    mins = [ff[i] for i in range(1, len(ff) - 1) if y[i] < y[i - 1] and y[i] < y[i + 1] and y[i] < -3]
    say(f"  dips at {np.round(mins, 2)} THz; change at 2/2.5/3 THz {np.round(ratio_db([2, 2.5, 3], x, t, w), 1)} dB")
    for tt in (0.6 * t, 0.8 * t, t):
        y = ratio_db(ff, x, tt, w)
        m = [ff[i] for i in range(1, len(ff) - 1) if y[i] < y[i - 1] and y[i] < y[i + 1]]
        say(f"  gas column {tt:.0f} um (growing in time): dips at {np.round(m[:3], 2)} THz")

say("\nReading: the first dip falls as the gas column grows (reported 0.7 -> 0.5 THz) and high")
say("frequencies rise because the gas path skips tissue loss. A dip at 0.6 THz needs a column of")
say("~200 um, i.e. an optical thickness ~2x the stated 100 um, which the reported n and 2 THz loss")
say("also imply. This is a two-parameter fit to digitised points, not a prediction; the raw data,")
say("the slice mass and an image of the cuvette during bias would test it.")

x, t = fits["KS1996"]
ff = np.linspace(0.2, 3.0, 600)
fig, ax = plt.subplots(figsize=(6, 3.6))
ax.errorbar(F_DIG, DB_DIG, yerr=2, fmt="ko", ms=4, label="reported, 5 min at 100 V (read off Fig. 1b)")
ax.plot(ff, ratio_db(ff, x, t), "C3", label=f"gas path: x = {x:.2f}, t = {t:.0f} µm (fit)")
for tt, ls in ((0.6 * t, ":"), (0.8 * t, "--")):
    ax.plot(ff, ratio_db(ff, x, tt), "C3", ls=ls, lw=0.8, label=f"same x, t = {tt:.0f} µm (earlier)")
ax.axhline(0, color="0.6", lw=0.6)
ax.set_xlabel("Frequency (THz)"); ax.set_ylabel("Δ power spectrum (dB)")
ax.legend(fontsize=7); fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"fig_gas_path.{ext}", dpi=200)
(OUT / "sim38.txt").write_text("\n".join(lines), encoding="utf-8")
