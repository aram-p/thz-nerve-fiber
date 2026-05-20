"""Sim 23 — Power dissipation per domain at σ_node = 100 S/m (peak coupling).

Sim 17 ran at σ_node = 1 S/m. The wide-σ sweep in Sim 2 v11 found that
the node-field response peaks at σ ≈ 100 S/m (resonant absorption
regime, just before the metallic limit). This sim re-runs the
power-dissipation-per-domain analysis at that peak σ, with the
water_eps_expr literal-complex fix in apply_materials.

Same frequency grid as Sim 17, same observable (ewfd.Qe integrated
per labelled domain).
"""

from __future__ import annotations

import csv
import time
from pathlib import Path

import matplotlib.pyplot as plt
import mph
import numpy as np

from thznerve.io.provenance import get_git_sha
from thznerve.model.geometry import DOMAIN_LABELS, GeometryParams, build_geometry, selection_tag
from thznerve.model.materials import MaterialParams, apply_materials
from thznerve.model.mesh import MeshParams, build_mesh
from thznerve.model.study import setup_physics, setup_study, solve_study
from thznerve.plots.style import DOUBLE_COL_MM, OKABE_ITO, apply_thesis_style, save_figure

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "results" / "sim23"

GEOM = GeometryParams(
    axon_radius_um=5, myelin_radius_um=7,
    node_length_um=40, internode_length_um=100, external_half_width_um=20,
)
MESH = MeshParams()
SIGMA_S_PER_M = 100.0  # peak coupling σ from sim 2
FREQ_THZ = np.array([0.15, 0.30, 0.45, 0.55, 0.60, 0.65, 0.70, 0.85, 1.00, 1.30, 1.60, 1.90, 2.10])


def _integrate_qe_per_domain(model, labels=DOMAIN_LABELS):
    java = model.java
    numerical = java.result().numerical()
    out = {}
    for label in labels:
        tag = f"intvol_{label}"
        if tag in [str(t) for t in numerical.tags()]:
            numerical.remove(tag)
        iv = numerical.create(tag, "IntVolume")
        iv.set("data", "dset1")
        iv.set("expr", ["ewfd.Qe"])
        iv.selection().named(selection_tag(label))
        result = iv.computeResult()
        try:
            real = float(list(result[0][0])[0])
        except Exception:
            real = float("nan")
        out[label] = real
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    apply_thesis_style(width_mm=DOUBLE_COL_MM, aspect=0.5)
    t0 = time.monotonic()

    client = mph.start()
    git_sha = get_git_sha()
    rows = []

    for i, f_thz in enumerate(FREQ_THZ, start=1):
        f_hz = float(f_thz) * 1e12
        model = client.create(f"sim23_freq_{i:02d}")
        build_geometry(model, GEOM)
        apply_materials(
            model, MaterialParams(node_sigma_S_per_m=SIGMA_S_PER_M),
            freq_hz=f_hz,
        )
        setup_physics(model, GEOM)
        n_elem = build_mesh(model, MESH)
        setup_study(model, f_hz)
        ts = time.monotonic()
        solve_study(model)
        powers = _integrate_qe_per_domain(model)
        dt = time.monotonic() - ts
        client.remove(model)

        total = sum(v for v in powers.values() if np.isfinite(v))
        row = [f_thz, powers["axon"], powers["myelin_proximal"], powers["node"],
               powers["myelin_distal"], powers["external"], total, dt]
        rows.append(row)
        print(f"  [{i:>2}/{len(FREQ_THZ)}] f={f_thz:5.3f} THz  "
              f"P_node={powers['node']:.3e}  P_total={total:.3e}  ({dt:4.1f}s)")

    client.clear()

    csv_path = OUT_DIR / "power_dissipation.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["frequency_THz", "P_axon", "P_myelin_proximal",
                    "P_node", "P_myelin_distal", "P_external", "P_total", "solve_s"])
        w.writerows(rows)

    arr = np.array(rows)
    f_thz_arr = arr[:, 0]
    p_axon, p_myp, p_node, p_myd, p_ext, p_total = arr[:, 1:7].T

    fig, (ax_abs, ax_rel) = plt.subplots(1, 2)

    ax_abs.plot(f_thz_arr, p_total, "k-", lw=1.8, label="total")
    ax_abs.plot(f_thz_arr, p_axon + p_ext, "-", color=OKABE_ITO[2], lw=1.2,
                label="water (axon + external)")
    ax_abs.plot(f_thz_arr, p_myp + p_myd, "-", color=OKABE_ITO[1], lw=1.2,
                label="myelin (prox + distal)")
    ax_abs.plot(f_thz_arr, p_node, "-", color=OKABE_ITO[6], lw=1.8,
                label=fr"node (σ = {SIGMA_S_PER_M:g} S/m)")
    ax_abs.set_xlabel("Frequency (THz)")
    ax_abs.set_ylabel("Absorbed power per cell (W)")
    ax_abs.set_title(fr"Power dissipation per domain, $\sigma_\mathrm{{node}}$ = {SIGMA_S_PER_M:g} S/m")
    ax_abs.set_yscale("log")
    ax_abs.axvline(0.6, color="red", ls=":", lw=0.7, alpha=0.6)
    ax_abs.axvline(2.0, color="red", ls=":", lw=0.7, alpha=0.6)
    ax_abs.legend(loc="best", fontsize=7, framealpha=0.92)
    ax_abs.grid(True, which="both", alpha=0.3)

    frac_node = p_node / np.where(p_total > 0, p_total, np.nan)
    ax_rel.plot(f_thz_arr, frac_node * 100, "o-", color=OKABE_ITO[6], lw=1.5, ms=4)
    ax_rel.set_xlabel("Frequency (THz)")
    ax_rel.set_ylabel("Fraction of total absorbed power in node (%)")
    ax_rel.set_title("Node share of absorption (σ = 100 S/m)")
    ax_rel.axvline(0.6, color="red", ls=":", lw=0.7, alpha=0.6)
    ax_rel.axvline(2.0, color="red", ls=":", lw=0.7, alpha=0.6)
    ax_rel.grid(True, alpha=0.3)

    fig.suptitle(fr"Sim 23 — Power dissipation at $\sigma_\mathrm{{node}}$ = {SIGMA_S_PER_M:g} S/m (peak σ from Sim 2)")
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim23_power_dissipation_sigma100")
    print(f"wrote {pdf}\nwrote {png}\nwrote {csv_path}")
    print(f"TOTAL: {time.monotonic() - t0:.1f}s")


if __name__ == "__main__":
    main()
