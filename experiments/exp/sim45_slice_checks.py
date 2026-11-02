"""Sim 45 - Review round 11 checks on the measured slice (no COMSOL).

(a) Optical thickness of the unbiased slice from its measured loss relative to the empty cuvette
    (paper 1 Fig. 3b/4b band medians, Sim 43) and from its delay (0.611 ps, linear phase with offset
    over 0.2-1.4 THz, Sim 43), for an effective-medium white-matter slab on PTFE with the
    Kindt-Schmuttenmaer (the model's) and the broadband water. The delay is inverted with the same
    linear-with-offset fit applied to the model's transmission phase, not with the phase index.
(b) Size of an opening that passes the undelayed component: for mode-matched detection the field
    fraction through an aperture equals the fraction of beam power through it. Gaussian beam 5 mm
    across (H1), read as 1/e^2 or as FWHM diameter; on-axis opening; and how b would change with
    frequency if the beam waist at the sample scales with wavelength.
Writes results/sim45/sim45.txt.
"""
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from thznerve.physics.effmedium import tissue_eps
from thznerve.physics.materials import WATER_BROADBAND, WATER_KS1996, axoplasm, myelin
from thznerve.physics.tmm import slab_rta

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim45"; OUT.mkdir(parents=True, exist_ok=True)
lines = []


def say(s=""):
    print(s, flush=True); lines.append(s)


N_PTFE = 1.43
MEAS_DB = {0.6: -12.1, 1.0: -14.9, 1.5: -20.1, 2.0: -29.2}       # Sim 43 (c), no voltage / empty cuvette
DELAY_PS, OFFSET_RAD = 0.611, -0.22                                # Sim 43 (d)


def t_rel(f, water, d_um, p=0.6):
    e, _ = tissue_eps(axoplasm(water)(f), myelin(water)(f), water(f), 0.7, p)
    t1 = slab_rta(f, 1.0, [(complex(e), d_um * 1e-6)], N_PTFE)[4]
    t0 = slab_rta(f, 1.0, [(1.0 + 0j, d_um * 1e-6)], N_PTFE)[4]
    return t1 / t0


say("(a) Slice optical thickness (effective-medium white matter, packing 0.6, on PTFE)")
for wname, w in (("Kindt-Schmuttenmaer", WATER_KS1996), ("broadband", WATER_BROADBAND)):
    out = []
    for f, db in MEAS_DB.items():
        g = lambda d: 20 * np.log10(abs(t_rel(f * 1e12, w, d))) - db
        out.append(f"{f} THz {brentq(g, 10, 1000):.0f} um")
    fr = np.linspace(0.2, 1.4, 49) * 1e12

    def fit(d):
        ph = np.unwrap(np.angle([t_rel(f, w, d) for f in fr]))   # exp(+jwt): delay = -dphi/domega
        c = np.polyfit(2 * np.pi * fr, -ph, 1)
        return c[0] * 1e12, -c[1]
    d_delay = brentq(lambda d: fit(d)[0] - DELAY_PS, 20, 1000)
    say(f"  {wname:20s} water: from loss {', '.join(out)}; from the delay (linear fit with offset, 0.2-1.4 THz) "
        f"{d_delay:.0f} um, model offset {fit(d_delay)[1]:+.2f} rad (measured {OFFSET_RAD:+.2f})")

say()
say("(b) Opening through the full thickness on the beam axis, Gaussian beam 5 mm across (H1)")
say("    field fraction b = fraction of beam power through the opening (mode-matched detection)")
for lab, w0 in (("1/e^2 diameter 5 mm", 2.5), ("FWHM 5 mm", 5.0 / np.sqrt(2 * np.log(2)))):
    dia = {b: 2 * np.sqrt(-np.log(1 - b) * w0 ** 2 / 2) for b in (0.04, 0.17, 0.22)}
    say(f"  {lab}: opening diameter for b = 0.04 / 0.17 / 0.22: " + " / ".join(f"{v:.2f}" for v in dia.values()) + " mm")
    for b in (0.17, 0.22):
        r2 = (dia[b] / 2) ** 2
        bf = {f: 1 - np.exp(-2 * r2 / (w0 * 0.6 / f) ** 2) for f in (1.0, 1.5, 2.0, 3.0)}
        say(f"    sized for b = {b} at 0.6 THz, waist ~ wavelength: b at 1/1.5/2/3 THz = "
            + " / ".join(f"{v:.2f}" for v in bf.values()))
say("  Measured: b 0.17 (4 min) and 0.22 (5 min) from the pulse fit, 0.20-0.21 and 0.22-0.27 flat over 1.5-3 THz.")
say("  An opening spread across the beam (area fraction b) gives b independent of frequency.")
(OUT / "sim45.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
