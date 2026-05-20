"""Sim 18 — Frequency sweep with E polarised parallel to the fibre.

This is the configuration paper 1 measured: the THz E-field oriented
parallel to the nerve fibres (not perpendicular). The previous sims
all used E ⟂ fibre because the default background field
``exp(-i k0 z) x̂`` has E ∥ x̂ while the fibre is along z.

Fix: rotate the *wave propagation* direction from z to x while keeping
the fibre along z, and set E along z. The new background field is
``exp(-i k0 x) ẑ`` — wave propagates along x, E parallel to fibre.

Caveats
-------
* The unit cell's x-extent is only 2·hw = 40 µm, much less than
  λ_water/2 ≈ 75 µm at 1 THz. The wave samples one "slab thickness"
  but does not develop a long propagation path — consistent with
  paper 1's 100 µm cuvette geometry.
* Periodic BCs run along y (lateral grating direction) and along z
  (fibre continuity).

26 frequencies, ~5 min compute.
"""

from __future__ import annotations

import csv
import time
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import mph
import numpy as np

from thznerve.io.hdf5 import write_result
from thznerve.io.provenance import get_git_sha
from thznerve.model.geometry import (
    GeometryParams, build_geometry, total_length_um,
)
from thznerve.model.materials import MaterialParams, apply_materials
from thznerve.model.mesh import MeshParams, build_mesh
from thznerve.model.study import (
    _add_box_selection, _evaluate_at_points, setup_study, solve_study,
)
from thznerve.plots.style import DOUBLE_COL_MM, apply_thesis_style, save_figure

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "results" / "sim18"

GEOM = GeometryParams(
    axon_radius_um=5, myelin_radius_um=7,
    node_length_um=40, internode_length_um=100, external_half_width_um=20,
)
MAT = MaterialParams(node_sigma_S_per_m=0.0)
MESH = MeshParams(max_element_size_um=30.0, node_element_size_um=5.0)
FREQ_THZ: np.ndarray = np.linspace(0.1, 2.0, 26)


def setup_physics_e_parallel_fibre(model: Any, geom_params: GeometryParams) -> None:
    """EWFD scattered-field: wave propagates along x, E along z (fibre axis).

    Inlet/outlet at x = ±hw (Scattering BC), periodic on y = ±hw and on
    z = 0 / z = L.
    """

    java = model.java
    comp_tag = str(java.component().tags()[0])
    comp = java.component(comp_tag)
    geom_tag = str(comp.geom().tags()[0])
    geom = java.component(comp_tag).geom(geom_tag)

    # Drop any existing EWFD
    for t in [str(x) for x in comp.physics().tags()]:
        comp.physics().remove(t)

    phys = comp.physics().create("ewfd", "ElectromagneticWavesFrequencyDomain", geom_tag)
    phys.label("EWFD — E ∥ fibre")
    phys.prop("BackgroundField").set("SolveFor", "scatteredField")
    # E along z (fibre), wave propagates along x.
    phys.prop("BackgroundField").set("Eb", ["0", "0", "exp(-i*ewfd.k0*x)"])

    L = total_length_um(geom_params)
    hw = geom_params.external_half_width_um
    eps = 0.5
    pad = hw * 0.5

    # Remove old endcap/periodic selections if re-running.
    for t in (
        "sel_inlet_x", "sel_outlet_x",
        "sel_y_minus", "sel_y_plus", "sel_y_pair",
        "sel_z_minus", "sel_z_plus", "sel_z_pair",
    ):
        if t in [str(x) for x in geom.feature().tags()]:
            geom.feature().remove(t)

    # Endcap faces (inlet/outlet): x = ±hw
    _add_box_selection(
        geom, "sel_inlet_x", "Inlet face (x = -hw)",
        entity_dim=2,
        xmin=-hw - eps, xmax=-hw + eps,
        ymin=-hw - pad, ymax=hw + pad,
        zmin=-eps,        zmax=L + eps,
    )
    _add_box_selection(
        geom, "sel_outlet_x", "Outlet face (x = +hw)",
        entity_dim=2,
        xmin=hw - eps,  xmax=hw + eps,
        ymin=-hw - pad, ymax=hw + pad,
        zmin=-eps,      zmax=L + eps,
    )

    # Lateral periodic pairs: y = ±hw, z = 0 / z = L
    _add_box_selection(geom, "sel_y_minus", "y = -hw",
                       entity_dim=2,
                       xmin=-hw - eps, xmax=hw + eps,
                       ymin=-hw - eps, ymax=-hw + eps,
                       zmin=-eps, zmax=L + eps)
    _add_box_selection(geom, "sel_y_plus", "y = +hw",
                       entity_dim=2,
                       xmin=-hw - eps, xmax=hw + eps,
                       ymin=hw - eps,  ymax=hw + eps,
                       zmin=-eps, zmax=L + eps)
    _add_box_selection(geom, "sel_z_minus", "z = 0 (fibre end)",
                       entity_dim=2,
                       xmin=-hw - eps, xmax=hw + eps,
                       ymin=-hw - eps, ymax=hw + eps,
                       zmin=-eps, zmax=eps)
    _add_box_selection(geom, "sel_z_plus", "z = L (fibre end)",
                       entity_dim=2,
                       xmin=-hw - eps, xmax=hw + eps,
                       ymin=-hw - eps, ymax=hw + eps,
                       zmin=L - eps, zmax=L + eps)

    for tag, inputs, label in (
        ("sel_y_pair", ["sel_y_minus", "sel_y_plus"], "y-periodic"),
        ("sel_z_pair", ["sel_z_minus", "sel_z_plus"], "z-periodic (fibre continuity)"),
    ):
        u = geom.feature().create(tag, "UnionSelection")
        u.set("entitydim", "2")
        u.set("input", inputs)
        u.label(label)

    geom.run()

    # Apply features
    for i, (sel_name, label) in enumerate(
        [("sel_inlet_x", "Inlet x=-hw"), ("sel_outlet_x", "Outlet x=+hw")], start=1
    ):
        sbc = phys.feature().create(f"sbc{i}", "Scattering", 2)
        sbc.selection().named(f"{geom_tag}_{sel_name}")
        sbc.label(label)
    for i, (sel_name, label) in enumerate(
        [("sel_y_pair", "Periodic y"), ("sel_z_pair", "Periodic z (fibre continuity)")], start=1
    ):
        pc = phys.feature().create(f"pc{i}", "PeriodicCondition", 2)
        pc.selection().named(f"{geom_tag}_{sel_name}")
        pc.set("PeriodicType", "Continuity")
        pc.label(label)


def _extract_axial_along_fibre(model: Any, geom_params: GeometryParams,
                                n_points: int = 200,
                                ) -> tuple[np.ndarray, np.ndarray]:
    """Axial sampling along the fibre (z-axis, at x=0, y=0)."""
    L = total_length_um(geom_params)
    z = np.linspace(0.5, L - 0.5, n_points)
    pts = np.column_stack([np.zeros_like(z), np.zeros_like(z), z])
    e = _evaluate_at_points(model, "ewfd.normE", pts)
    return z, e


def _extract_node_annulus_perp(model: Any, geom_params: GeometryParams,
                                n_z: int = 30) -> np.ndarray:
    """Sample inside the node annulus.

    In this rotated config, the node annulus is at r in (y, z plane)
    around the fibre axis (=z-axis). At z=z_node_mid, r=node_mid (=6 µm),
    sample at (0, r, z_node_mid_z) and rotations around the axis.

    Wait — fibre is along z, the annulus is around the fibre axis in
    the (x, y) plane. r=6 µm at varying (x, y) angles, at z in node range.
    """
    L_inter = geom_params.internode_length_um
    L_node = geom_params.node_length_um
    z_grid_node = np.linspace(L_inter + 1, L_inter + L_node - 1, n_z)
    r_node = 0.5 * (geom_params.axon_radius_um + geom_params.myelin_radius_um)

    # Sample at multiple (x, y) angles for each z (the annulus is a circle)
    angles = np.linspace(0, 2 * np.pi, 8, endpoint=False)
    pts_list = []
    for z in z_grid_node:
        for ang in angles:
            pts_list.append((r_node * np.cos(ang), r_node * np.sin(ang), z))
    pts = np.array(pts_list)
    return _evaluate_at_points(model, "ewfd.normE", pts)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    apply_thesis_style(width_mm=DOUBLE_COL_MM, aspect=0.5)
    t0 = time.monotonic()

    client = mph.start()
    model = client.create("sim18_e_par_fibre")
    build_geometry(model, GEOM)
    apply_materials(model, MAT)
    setup_physics_e_parallel_fibre(model, GEOM)
    n_elem = build_mesh(model, MESH)
    setup_study(model, FREQ_THZ[0] * 1e12)
    print(f"setup: {n_elem} elements, {time.monotonic() - t0:.1f}s")

    git_sha = get_git_sha()
    rows = []

    for i, f_thz in enumerate(FREQ_THZ, start=1):
        f_hz = float(f_thz) * 1e12
        java = model.java
        java.study(str(java.study().tags()[0])).feature("freq").set("plist", str(f_hz))
        ts = time.monotonic()
        java.study(str(java.study().tags()[0])).run()
        z_axis, e_axis = _extract_axial_along_fibre(model, GEOM)
        e_ann = _extract_node_annulus_perp(model, GEOM)
        dt = time.monotonic() - ts

        L_inter = GEOM.internode_length_um
        L_node = GEOM.node_length_um
        node_mask = (z_axis >= L_inter) & (z_axis <= L_inter + L_node)
        peak_node_axis = float(np.max(e_axis[node_mask])) if node_mask.any() else float("nan")
        peak_node_annulus = float(np.max(e_ann))
        mean_node_annulus = float(np.mean(e_ann))
        peak_global_axis = float(np.max(e_axis))

        h5_path = OUT_DIR / f"freq_{i:03d}_{f_thz:.3f}THz.h5"
        write_result(
            h5_path,
            scalars={
                "frequency_THz": f_thz,
                "peak_E_node_axis": peak_node_axis,
                "peak_E_node_annulus": peak_node_annulus,
                "mean_E_node_annulus": mean_node_annulus,
                "peak_E_global_axis": peak_global_axis,
                "solve_seconds": dt,
            },
            axial_profile=(z_axis, e_axis),
            axial_slice=(np.array([0.0]), z_axis, e_axis[None, :]),
            metadata={"git_sha": git_sha, "sim": "sim18_e_parallel_fibre",
                      "params": {**GEOM.model_dump(), **MAT.model_dump(),
                                 "freq_THz": f_thz},
                      "comsol_version": "6.3", "config": "E parallel to fibre"},
        )
        rows.append((f_thz, peak_node_axis, peak_node_annulus,
                     mean_node_annulus, peak_global_axis, dt))
        print(f"  [{i:>2}/{len(FREQ_THZ)}] f={f_thz:5.3f} THz  "
              f"|E|peak_axis={peak_node_axis:5.3f}  "
              f"|E|peak_ann={peak_node_annulus:5.3f}  ({dt:4.1f}s)")

    client.clear()

    csv_path = OUT_DIR / "spectrum.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["frequency_THz", "peak_E_node_axis", "peak_E_node_annulus",
                    "mean_E_node_annulus", "peak_E_global_axis", "solve_s"])
        w.writerows(rows)

    arr = np.array(rows)
    f_thz_arr = arr[:, 0]
    peak_axis = arr[:, 1]
    peak_ann = arr[:, 2]
    mean_ann = arr[:, 3]

    fig, (ax_a, ax_b) = plt.subplots(1, 2, sharex=True)
    ax_a.plot(f_thz_arr, peak_axis, "o-", lw=1.4, ms=4, label="peak |E| on fibre axis (in node z-range)")
    ax_a.plot(f_thz_arr, mean_ann, "^:", lw=1.0, ms=3, alpha=0.7, label="mean |E| in node annulus")
    ax_a.axvline(0.6, color="red", ls=":", lw=0.8, alpha=0.7)
    ax_a.axvline(2.0, color="red", ls=":", lw=0.8, alpha=0.7)
    ax_a.set_xlabel("Frequency (THz)")
    ax_a.set_ylabel("|E| (normalised, scattered field)")
    ax_a.set_title("Axial / mean-annulus")
    ax_a.legend(loc="best", framealpha=0.92, fontsize=7)
    ax_a.grid(True, alpha=0.3)

    ax_b.plot(f_thz_arr, peak_ann, "s-", lw=1.4, ms=4, color="C2",
              label="peak |E| in node annulus")
    ax_b.axvline(0.6, color="red", ls=":", lw=0.8, alpha=0.7)
    ax_b.axvline(2.0, color="red", ls=":", lw=0.8, alpha=0.7)
    ax_b.set_xlabel("Frequency (THz)")
    ax_b.set_ylabel("|E|")
    ax_b.set_title("Peak |E| in node annulus")
    ax_b.legend(loc="best", framealpha=0.92, fontsize=7)
    ax_b.grid(True, alpha=0.3)

    fig.suptitle("Sim 18 — Frequency sweep with E ∥ fibre (wave along x, E along z)")
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim18_e_parallel_fibre")
    print(f"wrote {pdf}\nwrote {png}\nwrote {csv_path}")
    print(f"TOTAL: {time.monotonic() - t0:.1f}s")


if __name__ == "__main__":
    main()
