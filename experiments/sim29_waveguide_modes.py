"""Sim 29 — Can a myelinated fibre guide THz? 2-D mode analysis.

The THz-signalling hypothesis treats the myelin sheath as a dielectric
waveguide that carries a THz wave along the axon, the way an optical fibre
carries light. This needs two things:

  1. a guided mode: field concentrated in or around the myelin, with
     Re(n_eff) between the water and myelin indices;
  2. low enough loss that the mode survives an internode (~100 µm), given
     that water absorbs THz strongly.

This sim solves the cross-section (axon water core, myelin ring, water
outside) as a 2-D mode analysis in EWFD at each frequency and reports, per
mode: the complex effective index n_eff, the 1/e power propagation length
L = 1 / (2·k0·|Im n_eff|), the fraction of power surviving one internode,
and the fraction of the guided power carried inside the myelin and axon.

The outer boundary is PEC on a box much wider than the fibre. Modes with
little power in the fibre are box modes, not guided modes; the confinement
column separates them.

Run: ``uv run python experiments/sim29_waveguide_modes.py``
2-D solves are fast; about 2-5 minutes in total.
"""

from __future__ import annotations

import csv
import time
from pathlib import Path

import matplotlib.pyplot as plt
import mph
import numpy as np

from thznerve.model.materials import (
    MYELIN_EPS_EXPR, MYELIN_EPS_REAL, _create_material, debye_water_epsilon, water_eps_expr,
)
from thznerve.plots.style import DOUBLE_COL_MM, OKABE_ITO, apply_thesis_style, save_figure

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "results" / "sim29"

AXON_R_UM = 5.0
MYELIN_R_UM = 7.0
BOX_HALF_WIDTH_UM = 60.0
INTERNODE_UM = 100.0
N_MODES = 8
MESH_MAX_UM = 4.0
MESH_FIBRE_UM = 0.5
FREQ_THZ = np.linspace(0.2, 2.1, 12)
C0 = 299_792_458.0

REGIONS = {  # label -> (x, y) inside the region
    "axon": (0.0, 0.0),
    "myelin": (0.5 * (AXON_R_UM + MYELIN_R_UM), 0.0),
    "water": (0.5 * (MYELIN_R_UM + BOX_HALF_WIDTH_UM), 0.0),
}


def _build(client, name: str, f_hz: float):
    model = client.create(name)
    java = model.java
    comp = java.component().create("comp1", True)
    geom = comp.geom().create("geom1", 2)
    geom.lengthUnit("um")

    sq = geom.feature().create("box", "Square")
    sq.set("size", str(2 * BOX_HALF_WIDTH_UM))
    sq.set("base", "center")
    for tag, r in (("c_myelin", MYELIN_R_UM), ("c_axon", AXON_R_UM)):
        c = geom.feature().create(tag, "Circle")
        c.set("r", str(r))
    for label, (x, y) in REGIONS.items():
        sel = geom.feature().create(f"sel_{label}", "DiskSelection")
        sel.set("entitydim", "2")
        sel.set("posx", str(x))
        sel.set("posy", str(y))
        sel.set("r", "0.1")
        sel.set("condition", "intersects")
    geom.run()

    water = water_eps_expr(f_hz)
    for label, expr in (("axon", water), ("myelin", MYELIN_EPS_EXPR), ("water", water)):
        mat = _create_material(comp, f"mat_{label}", label, expr)
        mat.selection().named(f"geom1_sel_{label}")

    comp.physics().create("ewfd", "ElectromagneticWavesFrequencyDomain", "geom1")

    mesh = comp.mesh().create("mesh1")
    mesh.feature("size").set("custom", "on")
    mesh.feature("size").set("hmax", f"{MESH_MAX_UM}[um]")
    tri = mesh.feature().create("ftri1", "FreeTri")
    fine = tri.feature().create("size_fibre", "Size")
    fine.selection().geom("geom1", 2)
    fine.selection().named("geom1_sel_myelin")
    fine.set("custom", "on")
    fine.set("hmax", f"{MESH_FIBRE_UM}[um]")
    fine.set("hmaxactive", True)
    mesh.run()

    study = java.study().create("std1")
    step = study.feature().create("mode", "ModeAnalysis")
    step.set("modeFreq", str(f_hz))
    step.set("neigs", str(N_MODES))
    # Search near the myelin index: a guided mode must sit between water and myelin.
    step.set("shift", f"{np.sqrt(MYELIN_EPS_REAL):.4f}")
    study.run()
    return model


def _integrate_per_mode(model, expr: str, selection: str | None) -> np.ndarray:
    numerical = model.java.result().numerical()
    tag = "int_tmp"
    if tag in [str(t) for t in numerical.tags()]:
        numerical.remove(tag)
    iv = numerical.create(tag, "IntSurface")
    iv.set("data", "dset1")
    iv.set("expr", [expr])
    if selection is None:
        iv.selection().all()
    else:
        iv.selection().named(selection)
    return np.array([row[0] for row in iv.getReal()], dtype=float)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    apply_thesis_style(width_mm=DOUBLE_COL_MM, aspect=0.5)
    t0 = time.monotonic()
    client = mph.start()

    rows = []
    for i, f_thz in enumerate(FREQ_THZ, start=1):
        f_hz = float(f_thz) * 1e12
        ts = time.monotonic()
        model = _build(client, f"sim29_{i}", f_hz)
        neff = np.atleast_1d(np.asarray(model.evaluate("ewfd.neff"), dtype=complex))
        p_total = _integrate_per_mode(model, "ewfd.Poavz", None)
        p_myelin = _integrate_per_mode(model, "ewfd.Poavz", "geom1_sel_myelin")
        p_axon = _integrate_per_mode(model, "ewfd.Poavz", "geom1_sel_axon")
        client.remove(model)

        k0 = 2 * np.pi * f_hz / C0
        n_water = np.sqrt(complex(debye_water_epsilon(f_hz)))
        for m, n in enumerate(neff):
            alpha = 2 * k0 * abs(n.imag)  # power attenuation, 1/m
            length_um = 1e6 / alpha if alpha > 0 else float("inf")
            survive = float(np.exp(-alpha * INTERNODE_UM * 1e-6))
            fibre_share = (p_myelin[m] + p_axon[m]) / p_total[m] if p_total[m] else float("nan")
            rows.append([f_thz, m, n.real, n.imag, n_water.real, length_um, survive,
                         p_myelin[m] / p_total[m] if p_total[m] else float("nan"), fibre_share])
        best = max((r for r in rows if r[0] == f_thz), key=lambda r: r[8])
        print(f"  [{i:>2}/{len(FREQ_THZ)}] f={f_thz:5.3f} THz  best mode: n_eff={best[2]:.3f}"
              f"{best[3]:+.3f}i  L={best[5]:7.1f} µm  internode survival={100 * best[6]:5.1f}%"
              f"  in fibre={100 * best[8]:5.1f}%  ({time.monotonic() - ts:4.1f}s)")

    client.clear()

    csv_path = OUT_DIR / "waveguide_modes.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["frequency_THz", "mode", "neff_real", "neff_imag", "n_water_real",
                    "propagation_length_um", "internode_survival", "myelin_share",
                    "fibre_share"])
        w.writerows(rows)

    # Per frequency, keep the mode with the most power inside the fibre.
    best_rows = [max((r for r in rows if r[0] == f), key=lambda r: r[8]) for f in FREQ_THZ]
    b = np.array(best_rows, dtype=float)

    fig, (ax_n, ax_l) = plt.subplots(1, 2)
    ax_n.plot(b[:, 0], b[:, 2], "o-", ms=3, color=OKABE_ITO[0], label="Re n_eff (most confined)")
    ax_n.plot(b[:, 0], b[:, 4], "--", color=OKABE_ITO[2], label="Re n_water")
    ax_n.axhline(np.sqrt(MYELIN_EPS_REAL), ls=":", color=OKABE_ITO[1], label="n_myelin")
    ax_n.set_xlabel("Frequency (THz)")
    ax_n.set_ylabel("Effective index")
    ax_n.set_title("Is there a mode between water and myelin?")
    ax_n.legend(fontsize=7, framealpha=0.92)
    ax_n.grid(True, alpha=0.3)

    ax_l.semilogy(b[:, 0], b[:, 5], "o-", ms=3, color=OKABE_ITO[0], label="1/e power length")
    ax_l.axhline(INTERNODE_UM, ls=":", color="red", label="internode (100 µm)")
    ax_l.set_xlabel("Frequency (THz)")
    ax_l.set_ylabel("Propagation length (µm)")
    ax_l.set_title("Does it survive one internode?")
    ax_l.legend(fontsize=7, framealpha=0.92)
    ax_l.grid(True, which="both", alpha=0.3)

    fig.suptitle("Sim 29 — THz guided modes of a myelinated fibre")
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim29_waveguide_modes")
    print(f"wrote {pdf}\nwrote {png}\nwrote {csv_path}")
    print(f"TOTAL: {time.monotonic() - t0:.1f}s")


if __name__ == "__main__":
    main()
