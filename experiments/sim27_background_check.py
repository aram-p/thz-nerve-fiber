"""Sim 27 — Is the 0.6 THz peak physical or a set-up artefact?

Sims 1-26 drive the scattered-field solve with a *vacuum* plane wave
E_bg = exp(-i·k0·z) x̂, but the whole cell is water. That background is not
a solution of the water-filled problem, so the scattered field is sourced
everywhere in the water, not only by the fibre. The end-caps also cut
through the fibre, so the scattering BCs sit on an inhomogeneous face.
Together these can create standing waves set by the 240 µm cell length,
which would explain why the peak appears even at σ = 0.

Four set-ups, same frequency sweep, σ = 0:

  A  vacuum background, no padding        (the original Sims 1-26 set-up)
  B  water-matched background, no padding
  C  water-matched background, 60 µm water buffer at each end
  D  water-matched background, 120 µm water buffer at each end

Each set-up is also solved with the fibre replaced by water ("reference"),
and the fibre result is reported as |E|_node / |E|_reference. A physical
fibre resonance keeps its frequency across B-D; a cell-length artefact
moves when the buffer length changes and flattens after normalisation.

Run: ``uv run python experiments/sim27_background_check.py``
~4 × 2 × 20 solves; budget about 30-40 minutes.
"""

from __future__ import annotations

import csv
import time
from pathlib import Path

import matplotlib.pyplot as plt
import mph
import numpy as np

from thznerve.model.geometry import GeometryParams, build_geometry, selection_tag
from thznerve.model.materials import (
    MaterialParams, _create_material, apply_materials, debye_water_epsilon, water_eps_expr,
)
from thznerve.model.mesh import MeshParams, build_mesh
from thznerve.model.study import _evaluate_at_points, setup_physics, setup_study, solve_study
from thznerve.plots.style import DOUBLE_COL_MM, OKABE_ITO, apply_thesis_style, save_figure

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "results" / "sim27"

BASE = dict(
    axon_radius_um=5, myelin_radius_um=7,
    node_length_um=40, internode_length_um=100, external_half_width_um=20,
)
SETUPS = {
    "A_vacuum_bg": dict(matched=False, pad=0.0),
    "B_water_bg": dict(matched=True, pad=0.0),
    "C_water_bg_pad60": dict(matched=True, pad=60.0),
    "D_water_bg_pad120": dict(matched=True, pad=120.0),
}
MAT = MaterialParams(node_sigma_S_per_m=0.0)
MESH = MeshParams(max_element_size_um=15.0, node_element_size_um=3.0)
FREQ_THZ = np.linspace(0.2, 2.1, 20)


def _replace_fibre_with_water(model, freq_hz: float) -> None:
    """Reference run: myelin and node become water, leaving a homogeneous cell."""
    comp = model.java.component(str(model.java.component().tags()[0]))
    expr = water_eps_expr(freq_hz)
    for label in ("myelin_proximal", "node", "myelin_distal"):
        mat = _create_material(comp, f"mat_{label}", f"Water in place of {label}", expr)
        mat.selection().named(selection_tag(label))


def _node_points(geom: GeometryParams) -> np.ndarray:
    z_lo = geom.internode_length_um
    z_hi = z_lo + geom.node_length_um
    r = 0.5 * (geom.axon_radius_um + geom.myelin_radius_um)
    zs = np.linspace(z_lo + 1, z_hi - 1, 16)
    angles = np.linspace(0, 2 * np.pi, 6, endpoint=False)
    return np.array([(r * np.cos(a), r * np.sin(a), z) for z in zs for a in angles])


def _solve(client, name: str, geom: GeometryParams, f_hz: float, *, matched: bool,
           reference: bool) -> float:
    model = client.create(name)
    build_geometry(model, geom)
    apply_materials(model, MAT, freq_hz=f_hz)
    if reference:
        _replace_fibre_with_water(model, f_hz)
    eps_bg = complex(debye_water_epsilon(f_hz)) if matched else None
    setup_physics(model, geom, background_eps=eps_bg)
    build_mesh(model, MESH)
    setup_study(model, f_hz)
    solve_study(model)
    value = float(np.mean(_evaluate_at_points(model, "ewfd.normE", _node_points(geom))))
    client.remove(model)
    return value


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    apply_thesis_style(width_mm=DOUBLE_COL_MM, aspect=0.5)
    t0 = time.monotonic()
    client = mph.start()

    rows = []
    for key, cfg in SETUPS.items():
        geom = GeometryParams(**BASE, end_padding_um=cfg["pad"])
        for i, f_thz in enumerate(FREQ_THZ, start=1):
            f_hz = float(f_thz) * 1e12
            ts = time.monotonic()
            e_fibre = _solve(client, f"sim27_{key}_{i}", geom, f_hz,
                             matched=cfg["matched"], reference=False)
            e_ref = _solve(client, f"sim27_{key}_{i}_ref", geom, f_hz,
                           matched=cfg["matched"], reference=True)
            rows.append([key, f_thz, e_fibre, e_ref, e_fibre / e_ref])
            print(f"  {key:<18} f={f_thz:5.3f} THz  |E|node={e_fibre:.3f}  "
                  f"ref={e_ref:.3f}  ratio={e_fibre / e_ref:.3f}  ({time.monotonic() - ts:4.1f}s)")

    client.clear()

    csv_path = OUT_DIR / "background_check.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["setup", "frequency_THz", "E_node_fibre", "E_node_reference", "ratio"])
        w.writerows(rows)

    fig, (ax_raw, ax_norm) = plt.subplots(1, 2)
    for color, key in zip(OKABE_ITO, SETUPS):
        sel = [r for r in rows if r[0] == key]
        f = np.array([r[1] for r in sel])
        ax_raw.plot(f, [r[2] for r in sel], "o-", ms=3, color=color, label=key)
        ax_norm.plot(f, [r[4] for r in sel], "o-", ms=3, color=color, label=key)
    for ax in (ax_raw, ax_norm):
        ax.axvline(0.6, color="red", ls=":", lw=0.7, alpha=0.6)
        ax.axvline(2.0, color="red", ls=":", lw=0.7, alpha=0.6)
        ax.set_xlabel("Frequency (THz)")
        ax.grid(True, alpha=0.3)
    ax_raw.set_ylabel("Mean |E| in node annulus")
    ax_raw.set_title("Raw node field")
    ax_norm.set_ylabel("|E| fibre / |E| water-only reference")
    ax_norm.set_title("Normalised to homogeneous cell")
    ax_norm.legend(fontsize=7, framealpha=0.92)
    fig.suptitle("Sim 27 — background-field and end-cap check")
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim27_background_check")
    print(f"wrote {pdf}\nwrote {png}\nwrote {csv_path}")
    print(f"TOTAL: {time.monotonic() - t0:.1f}s")


if __name__ == "__main__":
    main()
