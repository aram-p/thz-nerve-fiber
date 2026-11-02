"""Sim 34 - Tissue-level predictions to compare with Hovhannisyan & Makaryan (spinal-cord white
matter slices, 70-100 um, fibres in the slice plane, THz at normal incidence, 100 V DC across
18 mm perpendicular to the fibres; features at 0.6 and 2 THz for E || fibres with DC on).

(a) Polarisation: transmission of a 100 um slab of aligned fibres for E || and E perp to the
    fibres (effective medium, thznerve.physics.effmedium) relative to the bare PTFE window.
(b) Candidate DC-voltage mechanisms, each as a change in slab transmission (dB):
    - node conductivity: physiological channel opening (~1e-4 S/m) and the d_sigma needed for
      9 dB, using the full-wave node field factor eta (Sim 30);
    - Joule heating of the slice by +5 K and +15 K (Liebe water(T) in every water-bearing domain);
    - water loss of 10 % and 30 % of the slab's water (thickness shrinks by the water volume).
    Each mechanism is smooth in frequency; the question is whether any can give narrow
    features at 0.6 and 2 THz.

Writes results/sim34/*.csv and paper/figures/fig_tissue_voltage.png
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from thznerve.physics.effmedium import tissue_eps
from thznerve.physics.materials import WATER_KS1996, Looyenga, PROTEIN, DRY_LIPID, water_liebe
from thznerve.physics.samplebound import DB_PER_NEPER_POWER, node_volume_fraction, required_delta_sigma
from thznerve.physics.tmm import slab_rta
from thznerve.runlog import RunLog

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim34"; OUT.mkdir(parents=True, exist_ok=True)
FIG = ROOT / "paper" / "figures"

f = np.linspace(0.1, 3.0, 291) * 1e12
d = 100e-6
G, P = 0.7, 0.6           # g-ratio, fibre area fraction (white matter)
N_PTFE = 1.43             # PTFE window, n ~ 1.43 and low loss at THz


def tissue(water):
    ax = Looyenga(water, PROTEIN, 0.85)
    my = Looyenga(water, DRY_LIPID, 0.40)
    return tissue_eps(ax(f), my(f), water(f), G, P)


def slab_db(eps, dd=d):
    """|T|^2 of a tissue film on PTFE relative to bare PTFE, dB. The 1.2 mm window is taken as a
    half-space (its own Fabry-Perot fringes cancel in a TDS ratio after windowing)."""
    out = []
    for fi, e in zip(f, eps):
        _, T, _, _, _ = slab_rta(fi, 1.0, [(complex(e), dd)], N_PTFE)
        _, T0, _, _, _ = slab_rta(fi, 1.0, [], N_PTFE)
        out.append(10 * np.log10(T / T0))
    return np.array(out)


def smooth(x):
    return x


with RunLog("sim34", "tissue-level polarisation + DC-voltage mechanisms (analytic, effective medium)") as log:
    w37 = water_liebe(37.0)
    e_par, e_perp = tissue(WATER_KS1996)
    db_par, db_perp = smooth(slab_db(e_par)), smooth(slab_db(e_perp))

    # node conductivity: eta from Sim 30 (E || fibre), node fraction for the Sim 30 fibre
    s30 = np.genfromtxt(ROOT / "results/sim30/base_par.csv", delimiter=",", names=True)
    eta = np.interp(f / 1e12, s30["f_THz"], s30["eta"])
    fn = node_volume_fraction(3.5, 5.0, 1.0, 1000.0, P)
    n_t = np.sqrt(e_par).real
    ds_9db = required_delta_sigma(f, 9.0, eta, fn, n_t, d)
    k0 = 2 * np.pi * f / 2.99792458e8
    d_eps_phys = 1e-4 / (2 * np.pi * f * 8.8541878128e-12)
    db_phys = DB_PER_NEPER_POWER * k0 * d_eps_phys * eta * fn / n_t * d

    # heating: water at T0 = 20 C (lab) vs +5 K and +15 K, Liebe model (valid < 1 THz)
    heat = {}
    for dT in (5.0, 15.0):
        e0 = tissue(water_liebe(20.0))[0]
        e1 = tissue(water_liebe(20.0 + dT))[0]
        heat[dT] = smooth(slab_db(e1)) - smooth(slab_db(e0))

    # water loss: remove a fraction x of the slab's water volume, thickness shrinks accordingly
    wf = P * (G ** 2 * 0.85 + (1 - G ** 2) * 0.40) + (1 - P)   # water volume fraction of tissue
    dry = {}
    for x in (0.1, 0.3):
        wf_new = wf * (1 - x) / (1 - wf * x)
        # Looyenga between bulk tissue water and its dry matter at the new water fraction
        dry_eps = Looyenga(WATER_KS1996, PROTEIN, 0.0)(f)  # dry matter ~ protein/lipid eps 2.3-2.5
        e_new = (wf_new * WATER_KS1996(f) ** (1 / 3) + (1 - wf_new) * dry_eps ** (1 / 3)) ** 3
        e_old = (wf * WATER_KS1996(f) ** (1 / 3) + (1 - wf) * dry_eps ** (1 / 3)) ** 3
        dry[x] = smooth(slab_db(e_new, d * (1 - wf * x))) - smooth(slab_db(e_old))

    i06, i2 = np.argmin(abs(f - 0.6e12)), np.argmin(abs(f - 2e12))
    rows = np.column_stack([f / 1e12, db_par, db_perp, ds_9db, db_phys, heat[5.0], heat[15.0], dry[0.1], dry[0.3]])
    np.savetxt(OUT / "tissue_voltage.csv", rows, delimiter=",",
               header="f_THz,dB_par,dB_perp,dsigma_node_for_9dB,dB_node_channels,dB_heat5K,dB_heat15K,dB_dry10,dB_dry30")
    lines = [
        f"water volume fraction of model tissue: {wf:.2f}",
        f"slab 100um rel. PTFE: E|| {db_par[i06]:.1f} dB @0.6, {db_par[i2]:.1f} dB @2 THz; Eperp {db_perp[i06]:.1f}, {db_perp[i2]:.1f}",
        f"dichroism (par - perp): {db_par[i06] - db_perp[i06]:+.2f} dB @0.6, {db_par[i2] - db_perp[i2]:+.2f} dB @2",
        f"node d_sigma for 9 dB: {ds_9db[i06]:.2e} S/m @0.6, {ds_9db[i2]:.2e} @2 THz (node fraction {fn:.2e})",
        f"physiological channel opening (1e-4 S/m): {db_phys[i06]:.1e} dB @0.6, {db_phys[i2]:.1e} dB @2",
        f"heating +5K: {heat[5.0][i06]:+.2f} dB @0.6; +15K: {heat[15.0][i06]:+.2f} dB @0.6 (Liebe, <1 THz)",
        f"water loss 10%: {dry[0.1][i06]:+.2f} / {dry[0.1][i2]:+.2f} dB; 30%: {dry[0.3][i06]:+.2f} / {dry[0.3][i2]:+.2f} dB (0.6 / 2 THz)",
    ]
    (OUT / "summary.txt").write_text("\n".join(lines))
    print("\n".join(lines))
    log.outcome = "; ".join(lines[1:5])[:250]

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.5))
    ax[0].plot(f / 1e12, db_par, label="E ∥ fibres")
    ax[0].plot(f / 1e12, db_perp, "--", label="E ⊥ fibres")
    ax[0].set_xlabel("Frequency (THz)"); ax[0].set_ylabel("|T|² rel. PTFE window (dB)")
    ax[0].set_title("(a) 100 µm aligned-fibre slab"); ax[0].legend(fontsize=7)
    ax[1].semilogy(f / 1e12, np.abs(db_phys), label="node channels open (1e-4 S/m)")
    lo = f <= 1e12                                   # Liebe water(T) valid below 1 THz
    ax[1].semilogy(f[lo] / 1e12, np.abs(heat[5.0][lo]), label="heating +5 K")
    ax[1].semilogy(f[lo] / 1e12, np.abs(heat[15.0][lo]), label="heating +15 K")
    ax[1].semilogy(f / 1e12, np.abs(dry[0.1]), label="10 % water loss")
    ax[1].semilogy(f / 1e12, np.abs(dry[0.3]), label="30 % water loss")
    ax[1].axhline(9, color="k", lw=0.8, ls=":", label="reported feature ~9 dB")
    ax[1].plot([], [], color="0.8", lw=0.6, label="H₂O vapour lines")
    for x in (0.557, 1.67, 1.72, 1.92, 2.16, 2.20):
        ax[1].axvline(x, color="0.8", lw=0.6)
    ax[1].set_xlabel("Frequency (THz)"); ax[1].set_ylabel("|Δ transmission| (dB)")
    ax[1].set_title("(b) DC-voltage mechanisms"); ax[1].legend(fontsize=6)
    ax[1].set_ylim(1e-13, 1e2)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"fig_tissue_voltage.{ext}", dpi=200)
