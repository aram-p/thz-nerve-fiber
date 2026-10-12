"""Sim 28 — Differential absorption, voltage off vs voltage on.

Paper 1 sees the 0.6 and 2 THz features only with the DC voltage applied.
In this model the voltage enters only through the node conductivity σ
(open ion channels), so the experimental observable is the *change* in
absorption between σ = 0 (voltage off) and σ > 0 (voltage on), not the
field at the node.

Uses the corrected set-up from Sim 27: water-matched background wave and a
60 µm water buffer at each end, so the end-caps sit in homogeneous water.
For each frequency and σ it integrates the EM loss density ``ewfd.Qe``
over every domain and reports:

  * absorbed fraction A(f, σ) = P_abs / P_inc for the whole cell, with P_inc
    taken at z = 0; P_abs includes the water buffers, so A is only a
    relative measure and ΔA is the quantity to compare with experiment;
  * ΔA(f, σ) = A(f, σ) − A(f, 0), the voltage-on minus voltage-off signal;
  * the node's own share of P_abs.

If the grating/ion-channel hypothesis is right, ΔA should peak near 0.6 and
2 THz. A featureless ΔA means the model, as built, does not explain why the
experiment needs the voltage.

Run: ``uv run python experiments/sim28_voltage_on_off.py``
~4 σ × 20 frequencies; budget about 20-30 minutes.
"""

from __future__ import annotations

import csv
import time
from pathlib import Path

import matplotlib.pyplot as plt
import mph
import numpy as np

from thznerve.model.geometry import DOMAIN_LABELS, GeometryParams, build_geometry, selection_tag
from thznerve.model.materials import MaterialParams, apply_materials, debye_water_epsilon
from thznerve.model.mesh import MeshParams, build_mesh
from thznerve.model.study import setup_physics, setup_study, solve_study
from thznerve.plots.style import DOUBLE_COL_MM, OKABE_ITO, apply_thesis_style, save_figure

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "results" / "sim28"

GEOM = GeometryParams(
    axon_radius_um=5, myelin_radius_um=7,
    node_length_um=40, internode_length_um=100, external_half_width_um=20,
    end_padding_um=60,
)
MESH = MeshParams(max_element_size_um=15.0, node_element_size_um=3.0)
SIGMAS_S_PER_M = [0.0, 10.0, 100.0, 1000.0]  # 0 = voltage off
FREQ_THZ = np.linspace(0.2, 2.1, 20)
ETA_0 = 376.730313668  # vacuum wave impedance, Ω


def _integrate_qe(model, label: str) -> float:
    numerical = model.java.result().numerical()
    tag = f"intvol_{label}"
    if tag in [str(t) for t in numerical.tags()]:
        numerical.remove(tag)
    iv = numerical.create(tag, "IntVolume")
    iv.set("data", "dset1")
    iv.set("expr", ["ewfd.Qe"])
    iv.selection().named(selection_tag(label))
    return float(list(iv.computeResult()[0][0])[0])


def _incident_power_w(f_hz: float) -> float:
    """Power of the unit-amplitude background wave crossing the cell at z = 0."""
    n = np.sqrt(complex(debye_water_epsilon(f_hz)))
    area_m2 = (2 * GEOM.external_half_width_um * 1e-6) ** 2
    return n.real / (2 * ETA_0) * area_m2


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    apply_thesis_style(width_mm=DOUBLE_COL_MM, aspect=0.5)
    t0 = time.monotonic()
    client = mph.start()

    rows = []
    for sigma in SIGMAS_S_PER_M:
        for i, f_thz in enumerate(FREQ_THZ, start=1):
            f_hz = float(f_thz) * 1e12
            ts = time.monotonic()
            model = client.create(f"sim28_s{sigma:g}_{i}")
            build_geometry(model, GEOM)
            apply_materials(model, MaterialParams(node_sigma_S_per_m=sigma), freq_hz=f_hz)
            setup_physics(model, GEOM, background_eps=complex(debye_water_epsilon(f_hz)))
            build_mesh(model, MESH)
            setup_study(model, f_hz)
            solve_study(model)
            powers = {label: _integrate_qe(model, label) for label in DOMAIN_LABELS}
            client.remove(model)

            p_total = sum(powers.values())
            absorbed = p_total / _incident_power_w(f_hz)
            node_share = powers["node"] / p_total if p_total else float("nan")
            rows.append([sigma, f_thz, absorbed, node_share, *powers.values()])
            print(f"  σ={sigma:>6g}  f={f_thz:5.3f} THz  A={absorbed:.4f}  "
                  f"node share={100 * node_share:5.2f}%  ({time.monotonic() - ts:4.1f}s)")

    client.clear()

    csv_path = OUT_DIR / "voltage_on_off.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sigma_S_per_m", "frequency_THz", "absorbed_fraction", "node_share",
                    *[f"P_{label}" for label in DOMAIN_LABELS]])
        w.writerows(rows)

    arr = np.array(rows, dtype=float)
    a_off = arr[arr[:, 0] == 0.0][:, 2]

    fig, (ax_a, ax_d) = plt.subplots(1, 2)
    for color, sigma in zip(OKABE_ITO, SIGMAS_S_PER_M):
        sel = arr[arr[:, 0] == sigma]
        ax_a.plot(sel[:, 1], sel[:, 2], "o-", ms=3, color=color, label=f"σ = {sigma:g} S/m")
        if sigma > 0:
            ax_d.plot(sel[:, 1], sel[:, 2] - a_off, "o-", ms=3, color=color,
                      label=f"σ = {sigma:g} S/m")
    for ax in (ax_a, ax_d):
        ax.axvline(0.6, color="red", ls=":", lw=0.7, alpha=0.6)
        ax.axvline(2.0, color="red", ls=":", lw=0.7, alpha=0.6)
        ax.set_xlabel("Frequency (THz)")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7, framealpha=0.92)
    ax_a.set_ylabel("Absorbed fraction of incident power")
    ax_a.set_title("Absorption per cell")
    ax_d.set_ylabel("ΔA = A(σ) − A(0)")
    ax_d.set_title("Voltage on minus voltage off")
    fig.suptitle("Sim 28 — differential absorption from node conductivity")
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim28_voltage_on_off")
    print(f"wrote {pdf}\nwrote {png}\nwrote {csv_path}")
    print(f"TOTAL: {time.monotonic() - t0:.1f}s")


if __name__ == "__main__":
    main()
