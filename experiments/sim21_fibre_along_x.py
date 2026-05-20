"""Sim 21 — Cleanly rotated fibre: cylinders along x, wave along z, E along x.

Sim 18 rotated the *background field* while keeping cylinders along z; the
result was a 40 µm propagation distance, much less than λ_water/2. Sim 21
does the rotation properly: cylinders are constructed with ``axis = "x"``,
the unit cell is *long* along x (240 µm = fibre length, periodic) and y
(grating period, periodic), and *deep enough* along z (200 µm) to let the
incident wave develop on either side of the fibre.

Periodic BCs:
* x = 0 / x = L (continuous fibre periodicity)
* y = ±perp_hw (grating periodicity perpendicular to E and k)

Scattering BCs:
* z = ±wave_hw (inlet / outlet faces, perpendicular to wave propagation)

Background field: ``Eb = [exp(-i*ewfd.k0*z), 0, 0]`` — E along x (fibre axis),
k along z. This is paper 1's resonance condition.
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
from thznerve.model.materials import (
    MYELIN_EPS_EXPR, MaterialParams, water_eps_expr, node_eps_expr, _create_material,
)
from thznerve.model.mesh import MeshParams
from thznerve.model.study import _evaluate_at_points, setup_study
from thznerve.plots.style import DOUBLE_COL_MM, apply_thesis_style, save_figure

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "results" / "sim21"

# Geometry (independent of canonical GeometryParams because the canonical
# layout assumes fibre-along-z; here it's fibre-along-x with separate
# perp / wave-direction extents).
AXON_R = 5.0
MYELIN_R = 7.0
NODE_L = 40.0
INTERNODE_L = 100.0
PERP_HW = 20.0       # y half-width (perpendicular to k AND fibre)
WAVE_HW = 100.0      # z half-width (wave direction) — ~λ_water/2 at 0.6 THz
L_FIBRE = 2 * INTERNODE_L + NODE_L  # 240 µm

FREQ_THZ = np.linspace(0.1, 2.0, 26)
SIGMA = 0.0


def _add_cyl_x(geom, tag, *, radius, length, x_pos):
    cyl = geom.feature().create(tag, "Cylinder")
    cyl.set("r", str(radius))
    cyl.set("h", str(length))
    cyl.set("axis", ["1", "0", "0"])
    cyl.setIndex("pos", str(x_pos), 0)
    cyl.setIndex("pos", "0", 1)
    cyl.setIndex("pos", "0", 2)


def _add_block(geom, tag, *, size, pos):
    blk = geom.feature().create(tag, "Block")
    for i, s in enumerate(size):
        blk.setIndex("size", str(s), i)
    for i, p in enumerate(pos):
        blk.setIndex("pos", str(p), i)


def _add_diff(geom, tag, label, *, inputs, subtract):
    d = geom.feature().create(tag, "Difference")
    d.selection("input").set(inputs)
    d.selection("input2").set(subtract)
    d.label(label)


def _add_ball(geom, tag, label, *, xyz, radius=0.5):
    s = geom.feature().create(tag, "BallSelection")
    s.set("entitydim", "3")
    s.set("posx", str(xyz[0])); s.set("posy", str(xyz[1])); s.set("posz", str(xyz[2]))
    s.set("r", str(radius))
    s.set("condition", "intersects")
    s.label(label)


def _add_box_sel(geom, tag, label, dim, x, y, z):
    s = geom.feature().create(tag, "BoxSelection")
    s.set("entitydim", str(dim))
    s.set("xmin", str(x[0])); s.set("xmax", str(x[1]))
    s.set("ymin", str(y[0])); s.set("ymax", str(y[1]))
    s.set("zmin", str(z[0])); s.set("zmax", str(z[1]))
    s.set("condition", "inside")
    s.label(label)


def build_geometry_x(model):
    java = model.java
    java.component().create("comp1", True)
    java.component("comp1").geom().create("geom1", 3)
    geom = java.component("comp1").geom("geom1")
    geom.lengthUnit("um")

    # Inner axon cylinder along x (full fibre)
    _add_cyl_x(geom, "cyl_axon", radius=AXON_R, length=L_FIBRE, x_pos=0.0)
    geom.feature("cyl_axon").label("Axon")

    # Four annular sheath segments built as Difference(outer, inner)
    # Proximal sheath  x ∈ [0, INTERNODE_L]
    _add_cyl_x(geom, "cyl_myp_out", radius=MYELIN_R, length=INTERNODE_L, x_pos=0.0)
    _add_cyl_x(geom, "cyl_myp_in",  radius=AXON_R,  length=INTERNODE_L, x_pos=0.0)
    _add_diff(geom, "myp", "Myelin proximal",
              inputs=["cyl_myp_out"], subtract=["cyl_myp_in"])

    # Node   x ∈ [INTERNODE_L, INTERNODE_L + NODE_L]
    _add_cyl_x(geom, "cyl_nd_out", radius=MYELIN_R, length=NODE_L, x_pos=INTERNODE_L)
    _add_cyl_x(geom, "cyl_nd_in",  radius=AXON_R,  length=NODE_L, x_pos=INTERNODE_L)
    _add_diff(geom, "nd", "Node of Ranvier",
              inputs=["cyl_nd_out"], subtract=["cyl_nd_in"])

    # Distal sheath   x ∈ [INTERNODE_L + NODE_L, L_FIBRE]
    x_dist = INTERNODE_L + NODE_L
    _add_cyl_x(geom, "cyl_myd_out", radius=MYELIN_R, length=INTERNODE_L, x_pos=x_dist)
    _add_cyl_x(geom, "cyl_myd_in",  radius=AXON_R,  length=INTERNODE_L, x_pos=x_dist)
    _add_diff(geom, "myd", "Myelin distal",
              inputs=["cyl_myd_out"], subtract=["cyl_myd_in"])

    # External: box minus the full-length myelin_r cylinder
    _add_block(geom, "box_ext",
               size=(L_FIBRE, 2 * PERP_HW, 2 * WAVE_HW),
               pos=(0.0, -PERP_HW, -WAVE_HW))
    _add_cyl_x(geom, "cyl_ext_in", radius=MYELIN_R, length=L_FIBRE, x_pos=0.0)
    _add_diff(geom, "ext", "External medium",
              inputs=["box_ext"], subtract=["cyl_ext_in"])

    # Ball selections per labeled domain (for material assignment).
    # Domain centres in (x, y, z): fibre is along x, so "radial" is sqrt(y²+z²).
    x_centres = {
        "axon": L_FIBRE / 2,
        "myelin_proximal": INTERNODE_L / 2,
        "node": INTERNODE_L + NODE_L / 2,
        "myelin_distal": INTERNODE_L + NODE_L + INTERNODE_L / 2,
        "external": L_FIBRE / 2,
    }
    r_mid_sheath = (AXON_R + MYELIN_R) / 2  # = 6
    r_mid_ext = (MYELIN_R + max(PERP_HW, WAVE_HW)) / 2
    centres = {
        "axon": (x_centres["axon"], 0.0, 0.0),
        "myelin_proximal": (x_centres["myelin_proximal"], r_mid_sheath, 0.0),
        "node": (x_centres["node"], r_mid_sheath, 0.0),
        "myelin_distal": (x_centres["myelin_distal"], r_mid_sheath, 0.0),
        "external": (x_centres["external"], r_mid_ext, 0.0),
    }
    for label, xyz in centres.items():
        _add_ball(geom, f"sel_{label}", f"Sel {label}", xyz=xyz)

    # Lateral periodic pair selections and inlet/outlet face selections
    eps = 0.5
    # Inlet at z = -wave_hw, outlet at z = +wave_hw
    _add_box_sel(geom, "sel_inlet_z", "Inlet (z = -wave_hw)", dim=2,
                 x=(-eps, L_FIBRE + eps),
                 y=(-PERP_HW - eps, PERP_HW + eps),
                 z=(-WAVE_HW - eps, -WAVE_HW + eps))
    _add_box_sel(geom, "sel_outlet_z", "Outlet (z = +wave_hw)", dim=2,
                 x=(-eps, L_FIBRE + eps),
                 y=(-PERP_HW - eps, PERP_HW + eps),
                 z=(WAVE_HW - eps, WAVE_HW + eps))
    # Periodic pairs:  x = 0 / x = L (fibre continuity), y = ±perp_hw (grating)
    _add_box_sel(geom, "sel_x_minus", "x = 0", dim=2,
                 x=(-eps, eps), y=(-PERP_HW - eps, PERP_HW + eps),
                 z=(-WAVE_HW - eps, WAVE_HW + eps))
    _add_box_sel(geom, "sel_x_plus", "x = L", dim=2,
                 x=(L_FIBRE - eps, L_FIBRE + eps),
                 y=(-PERP_HW - eps, PERP_HW + eps),
                 z=(-WAVE_HW - eps, WAVE_HW + eps))
    _add_box_sel(geom, "sel_y_minus", "y = -perp_hw", dim=2,
                 x=(-eps, L_FIBRE + eps),
                 y=(-PERP_HW - eps, -PERP_HW + eps),
                 z=(-WAVE_HW - eps, WAVE_HW + eps))
    _add_box_sel(geom, "sel_y_plus", "y = +perp_hw", dim=2,
                 x=(-eps, L_FIBRE + eps),
                 y=(PERP_HW - eps, PERP_HW + eps),
                 z=(-WAVE_HW - eps, WAVE_HW + eps))
    for tag, inputs, label in (
        ("sel_x_pair", ["sel_x_minus", "sel_x_plus"], "x-periodic (fibre)"),
        ("sel_y_pair", ["sel_y_minus", "sel_y_plus"], "y-periodic (grating)"),
    ):
        u = geom.feature().create(tag, "UnionSelection")
        u.set("entitydim", "2")
        u.set("input", inputs)
        u.label(label)

    geom.run()
    return int(geom.getNDomains())


def apply_materials_x(model, sigma, freq_hz):
    java = model.java
    comp_tag = str(java.component().tags()[0])
    comp = java.component(comp_tag)
    geom_tag = str(comp.geom().tags()[0])

    def sel(label):
        return f"{geom_tag}_sel_{label}"

    water_expr = water_eps_expr(freq_hz=freq_hz)
    for tag, label, expr, dom in [
        ("mat_water", "Water (axon)", water_expr, "axon"),
        ("mat_water_ext", "Water (ext)", water_expr, "external"),
        ("mat_myp", "Myelin proximal", MYELIN_EPS_EXPR, "myelin_proximal"),
        ("mat_myd", "Myelin distal",  MYELIN_EPS_EXPR, "myelin_distal"),
        ("mat_node", "Node (water + σ)",
         node_eps_expr(sigma, freq_hz=freq_hz), "node"),
    ]:
        m = _create_material(comp, tag, label, expr)
        m.selection().named(sel(dom))


def setup_physics_x(model):
    java = model.java
    comp_tag = str(java.component().tags()[0])
    comp = java.component(comp_tag)
    geom_tag = str(comp.geom().tags()[0])
    for t in [str(x) for x in comp.physics().tags()]:
        comp.physics().remove(t)
    phys = comp.physics().create("ewfd", "ElectromagneticWavesFrequencyDomain", geom_tag)
    phys.label("EWFD — fibre along x")
    phys.prop("BackgroundField").set("SolveFor", "scatteredField")
    # E along x (fibre axis), k along z (perpendicular)
    phys.prop("BackgroundField").set("Eb", ["exp(-i*ewfd.k0*z)", "0", "0"])
    # Scattering BCs on z = ±wave_hw
    for i, sel_name in enumerate(["sel_inlet_z", "sel_outlet_z"], start=1):
        sbc = phys.feature().create(f"sbc{i}", "Scattering", 2)
        sbc.selection().named(f"{geom_tag}_{sel_name}")
        sbc.label(f"Scattering ({sel_name})")
    # Periodic BCs on x-pair (fibre continuity) and y-pair (grating)
    for i, sel_name in enumerate(["sel_x_pair", "sel_y_pair"], start=1):
        pc = phys.feature().create(f"pc{i}", "PeriodicCondition", 2)
        pc.selection().named(f"{geom_tag}_{sel_name}")
        pc.set("PeriodicType", "Continuity")
        pc.label(f"Periodic ({sel_name})")


def build_mesh_x(model):
    java = model.java
    comp_tag = str(java.component().tags()[0])
    comp = java.component(comp_tag)
    geom_tag = str(comp.geom().tags()[0])
    if comp.mesh().tags().length == 0:
        comp.mesh().create("mesh1", geom_tag)
    mesh = comp.mesh(str(comp.mesh().tags()[0]))
    for t in [str(x) for x in mesh.feature().tags()]:
        if t == "size":
            continue
        mesh.feature().remove(t)
    mesh.feature("size").set("custom", "on")
    mesh.feature("size").set("hmax", "30[um]")
    mesh.feature("size").set("hmin", "1[um]")
    ftet = mesh.feature().create("ftet1", "FreeTet")
    ftet.label("Free Tetrahedral")
    size_node = ftet.feature().create("size_node", "Size")
    size_node.label("Node refinement")
    size_node.selection().geom(geom_tag, 3)
    size_node.selection().named(f"{geom_tag}_sel_node")
    size_node.set("custom", "on")
    size_node.set("hmax", "5[um]")
    size_node.set("hmaxactive", True)
    mesh.run()
    return int(mesh.getNumElem())


def sample_annulus_around_x_axis(model, n_z=20):
    """Sample |E| in the node annulus.

    With fibre along x, the annulus is at r = 6 µm in the (y, z) plane,
    around the x-axis. Sample at multiple azimuthal angles and x positions
    inside the node x-range.
    """
    r = (AXON_R + MYELIN_R) / 2
    x_grid = np.linspace(INTERNODE_L + 1, INTERNODE_L + NODE_L - 1, n_z)
    angles = np.linspace(0, 2 * np.pi, 8, endpoint=False)
    pts = []
    for x in x_grid:
        for a in angles:
            pts.append((x, r * np.cos(a), r * np.sin(a)))
    return _evaluate_at_points(model, "ewfd.normE", np.array(pts))


def sample_axis_through_fibre(model, n=200):
    x = np.linspace(0.5, L_FIBRE - 0.5, n)
    pts = np.column_stack([x, np.zeros_like(x), np.zeros_like(x)])
    return x, _evaluate_at_points(model, "ewfd.normE", pts)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    apply_thesis_style(width_mm=DOUBLE_COL_MM, aspect=0.5)
    t0 = time.monotonic()

    client = mph.start()
    git_sha = get_git_sha()
    rows = []

    for i, f_thz in enumerate(FREQ_THZ, start=1):
        f_hz = float(f_thz) * 1e12
        model = client.create(f"sim21_f{i:02d}")
        n_dom = build_geometry_x(model)
        assert n_dom == 5, f"expected 5 domains, got {n_dom}"
        apply_materials_x(model, SIGMA, f_hz)
        setup_physics_x(model)
        n_elem = build_mesh_x(model)
        setup_study(model, f_hz)
        ts = time.monotonic()
        from thznerve.model.study import solve_study
        solve_study(model)
        x_ax, e_ax = sample_axis_through_fibre(model)
        e_ann = sample_annulus_around_x_axis(model)
        dt = time.monotonic() - ts
        client.remove(model)

        # Mask the node x-range on the axial profile
        x_node_lo = INTERNODE_L
        x_node_hi = INTERNODE_L + NODE_L
        mask = (x_ax >= x_node_lo) & (x_ax <= x_node_hi)
        peak_axis = float(np.max(e_ax[mask])) if mask.any() else float("nan")
        peak_ann = float(np.max(e_ann))
        mean_ann = float(np.mean(e_ann))
        global_axis = float(np.max(e_ax))

        rows.append((f_thz, peak_axis, peak_ann, mean_ann, global_axis, dt, n_elem))
        print(f"  [{i:>2}/{len(FREQ_THZ)}] f={f_thz:5.3f} THz  "
              f"|E|_node_axis={peak_axis:5.3f}  |E|_node_ann={peak_ann:5.3f}  "
              f"({n_elem:5d} elem, {dt:4.1f}s)")

    client.clear()
    csv_path = OUT_DIR / "spectrum.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["frequency_THz", "peak_E_node_axis", "peak_E_node_annulus",
                    "mean_E_node_annulus", "peak_E_global_axis", "solve_s", "n_elem"])
        w.writerows(rows)

    arr = np.array(rows)
    f_arr = arr[:, 0]
    p_ax = arr[:, 1]
    p_an = arr[:, 2]
    m_an = arr[:, 3]

    fig, (ax_a, ax_b) = plt.subplots(1, 2, sharex=True)
    ax_a.plot(f_arr, p_ax, "o-", color="C0", lw=1.4, ms=4, label="axial (along fibre)")
    ax_a.plot(f_arr, p_an, "s--", color="C3", lw=1.4, ms=4, alpha=0.75,
              label="annular (around fibre)")
    ax_a.plot(f_arr, m_an, "^:", color="C2", lw=1.0, ms=3, alpha=0.7,
              label="mean annular")
    ax_a.axvline(0.6, color="red", ls=":", lw=0.7, alpha=0.7)
    ax_a.axvline(2.0, color="red", ls=":", lw=0.7, alpha=0.7)
    ax_a.set_xlabel("Frequency (THz)")
    ax_a.set_ylabel("|E| (normalised)")
    ax_a.set_title("Sim 21 — fibre along x, E ∥ fibre, wave along z")
    ax_a.grid(True, alpha=0.3)
    ax_a.legend(loc="best", fontsize=7, framealpha=0.92)

    # Right: compare to Sim 1
    sim1_csv = REPO_ROOT / "results" / "sim1" / "spectrum.csv"
    if sim1_csv.exists():
        with sim1_csv.open() as f:
            r = csv.reader(f); next(r)
            d = np.array([[float(x) for x in row] for row in r])
        ax_b.plot(d[:, 0], d[:, 1], "o-", color="0.4", lw=1.2, ms=3.5,
                  label="Sim 1 (E ⊥ fibre, fibre along z, axial)")
    ax_b.plot(f_arr, p_an, "s-", color="C3", lw=1.4, ms=4,
              label="Sim 21 (E ∥ fibre, fibre along x, annular)")
    ax_b.axvline(0.6, color="red", ls=":", lw=0.7, alpha=0.7)
    ax_b.axvline(2.0, color="red", ls=":", lw=0.7, alpha=0.7)
    ax_b.set_xlabel("Frequency (THz)")
    ax_b.set_ylabel("|E|")
    ax_b.set_title("Sim 21 vs Sim 1 — polarisation comparison")
    ax_b.grid(True, alpha=0.3)
    ax_b.legend(loc="best", fontsize=7, framealpha=0.92)

    fig.suptitle("Sim 21 — clean fibre reorientation (cylinder rotated, not just BG field)")
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim21_fibre_along_x")
    print(f"wrote {pdf}\nwrote {png}\nwrote {csv_path}")
    print(f"TOTAL: {time.monotonic() - t0:.1f}s")


if __name__ == "__main__":
    main()
