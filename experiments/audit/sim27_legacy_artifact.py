"""Sim 27 - Phase 0 audit of the legacy 3D pipeline (Sims 1-26).

Hypothesis: the legacy 0.6 / 2 THz "resonances" are artefacts of the
scattered-field set-up (background field exp(-i*k0*z) with the *vacuum* k0 in a
water-filled box, i.e. a volume source k0^2 (eps_w - 1) E_b in every water
element, plus end-cap reflections), not fibre physics.

Decisive test: run the legacy pipeline unchanged, (A) with the fibre, (B) with
the fibre replaced by water (no scatterer at all), (C) no fibre and a longer box
(internode 100 -> 130 um, box 240 -> 300 um). Physics says B and C must give
|E(z)| = |exp(-i k_w z)| (a decaying plane wave) with no spectral peaks. If B/C
show peaks, or A's peaks move with box length, they are artefacts.

Also overlays the closed-form 1D solution of the legacy equations for B/C.
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import mph
import numpy as np

from thznerve.model.geometry import GeometryParams, build_geometry, total_length_um
from thznerve.model.materials import (
    MaterialParams, apply_materials, debye_water_epsilon, water_eps_expr,
)
from thznerve.model.mesh import MeshParams, build_mesh
from thznerve.model.study import _evaluate_at_points, setup_physics, setup_study
from thznerve.runlog import RunLog

OUT = Path(__file__).resolve().parents[2] / "results" / "sim27"
C0 = 299792458.0
FREQ_THZ = np.round(np.arange(0.2, 2.31, 0.1), 3)


def legacy_1d(f_hz: float, z_um: np.ndarray, L_um: float) -> np.ndarray:
    """Closed-form total field of the legacy formulation in a homogeneous water box.

    E_sc'' + k^2 E_sc = -(k^2 - k0^2) e^{-j k0 z}, first-order absorbing BCs at
    z=0 and z=L (COMSOL e^{jwt} convention, k = k0 sqrt(eps_w), Im k < 0).
    E_tot = B e^{-jkz} + C e^{jkz}, B=(k0+k)/2k, C=-(k0-k) e^{-j k0 L} e^{-jkL}/2k.
    """
    k0 = 2 * np.pi * f_hz / C0
    eps = complex(debye_water_epsilon(f_hz))
    k = k0 * np.sqrt(eps)
    if k.imag > 0:
        k = np.conj(k)
    z, L = z_um * 1e-6, L_um * 1e-6
    B = (k0 + k) / (2 * k)
    C = -(k0 - k) * np.exp(-1j * k0 * L) * np.exp(-1j * k * L) / (2 * k)
    return B * np.exp(-1j * k * z) + C * np.exp(1j * k * z)


def set_all_water(model, f_hz: float) -> None:
    comp = model.java.component("comp1")
    w = water_eps_expr(f_hz)
    for tag in ("mat_myelin_proximal", "mat_myelin_distal", "mat_node"):
        comp.material(tag).propertyGroup("def").set(
            "relpermittivity", [w, "0", "0", "0", w, "0", "0", "0", w])


def run_config(client, name: str, geom: GeometryParams, fibre: bool) -> list[tuple]:
    rows = []
    with RunLog("sim27", f"legacy pipeline config {name}") as log:
        model = client.create(f"sim27_{name}")
        build_geometry(model, geom)
        apply_materials(model, MaterialParams(node_sigma_S_per_m=0.0), freq_hz=FREQ_THZ[0] * 1e12)
        setup_physics(model, geom)
        build_mesh(model, MeshParams())
        setup_study(model, FREQ_THZ[0] * 1e12)
        L = total_length_um(geom)
        zlo, zhi = geom.internode_length_um, geom.internode_length_um + geom.node_length_um
        z_ax = np.linspace(0.5, L - 0.5, 121)
        ax_pts = np.column_stack([np.zeros_like(z_ax), np.zeros_like(z_ax), z_ax])
        r_ann = 0.5 * (geom.axon_radius_um + geom.myelin_radius_um)
        zg = np.linspace(zlo + 1, zhi - 1, 12)
        ang = np.linspace(0, 2 * np.pi, 6, endpoint=False)
        ann_pts = np.array([(r_ann * np.cos(a), r_ann * np.sin(a), z) for z in zg for a in ang])
        for f_thz in FREQ_THZ:
            f_hz = float(f_thz) * 1e12
            apply_materials(model, MaterialParams(node_sigma_S_per_m=0.0), freq_hz=f_hz)
            if not fibre:
                set_all_water(model, f_hz)
            model.java.study("std1").feature("freq").set("plist", str(f_hz))
            ts = time.monotonic()
            model.java.study("std1").run()
            e_ax = _evaluate_at_points(model, "ewfd.normE", ax_pts)
            e_ann = _evaluate_at_points(model, "ewfd.normE", ann_pts)
            m = (z_ax >= zlo) & (z_ax <= zhi)
            e_1d = np.abs(legacy_1d(f_hz, z_ax, L))
            k = 2 * np.pi * f_hz / C0 * np.sqrt(complex(debye_water_epsilon(f_hz)))
            e_true = np.abs(np.exp(-1j * np.conj(k) * z_ax * 1e-6)) if k.imag > 0 else np.abs(np.exp(-1j * k * z_ax * 1e-6))
            rows.append((name, f_thz, L, float(e_ax[m].max()), float(e_ann.mean()),
                         float(e_1d[m].max()), float(e_true[m].max()),
                         float(np.max(np.abs(e_ax - e_1d))), time.monotonic() - ts))
            print(f"{name} f={f_thz:.2f} ax={rows[-1][3]:.3f} ann={rows[-1][4]:.3f} "
                  f"1d={rows[-1][5]:.3f} true={rows[-1][6]:.3f} |FEM-1D|max={rows[-1][7]:.3f}",
                  flush=True)
            np.save(OUT / f"{name}_{f_thz:.2f}_axial.npy", np.vstack([z_ax, e_ax, e_1d, e_true]))
        log.dof = len(model.java.sol("sol1").getU())
        log.outcome = f"{len(rows)} freqs ok"
        client.remove(model)
    return rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mph.option("session", "stand-alone")
    client = mph.start(cores=4)
    base = dict(axon_radius_um=5, myelin_radius_um=7, node_length_um=40, external_half_width_um=20)
    configs = [
        ("A_fibre_L240", GeometryParams(internode_length_um=100, **base), True),
        ("B_water_L240", GeometryParams(internode_length_um=100, **base), False),
        ("C_water_L300", GeometryParams(internode_length_um=130, **base), False),
        ("D_fibre_L300", GeometryParams(internode_length_um=130, **base), True),
    ]
    only = sys.argv[1:]
    allrows = []
    for name, g, fib in configs:
        if only and name not in only:
            continue
        allrows += run_config(client, name, g, fib)
        with (OUT / "summary.csv").open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["config", "f_THz", "L_um", "ax_peak_node", "ann_mean_node",
                        "legacy1d_peak_node", "true_planewave_node", "max_abs_FEM_minus_1D", "solve_s"])
            w.writerows(allrows)


if __name__ == "__main__":
    main()
