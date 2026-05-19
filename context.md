# Tutor presentation context — overnight session 2026-05-19

This file is a live progress log. The intent: a tutor-grade presentation tomorrow showing **at least five simulations** that investigate different aspects of the THz / single-myelinated-nerve-fiber problem. The narrative thread is Hovhannisyan & Makaryan's experimental observation of resonant THz absorption in spinal-cord samples at **0.6 THz** and **2 THz** when E ∥ nerve fibers + DC voltage applied. The question this project is meant to ask: *does a first-principles 3D EM simulation of a single myelinated fiber reproduce those resonances, and what physical parameters control them?*

## Plan

| #  | Investigation                                  | Tool                                   | Figure                                  | Status |
|----|------------------------------------------------|----------------------------------------|-----------------------------------------|--------|
| 1  | Frequency sweep, baseline geometry             | COMSOL EWFD if working, else Mie       | `figures/out/sim1_freq_sweep.*`         | TBD    |
| 2  | Node conductivity σ at 0.6 THz                 | COMSOL EWFD if working, else thin-wire | `figures/out/sim2_sigma_sweep.*`        | TBD    |
| 3  | Node-length sensitivity                        | COMSOL EWFD if working, else analytic  | `figures/out/sim3_node_length.*`        | TBD    |
| 4  | Material dispersion (water, myelin, node)      | Pure Python (no COMSOL)                | `figures/out/sim4_dispersion.*`         | TBD    |
| 5  | 3D field heatmap / geometry visualisation       | COMSOL if working, else matplotlib 3D  | `figures/out/sim5_3d.*`                 | TBD    |
| 6  | Diffraction-grating analytic spectrum          | Pure Python (per paper 3's framing)    | `figures/out/sim6_grating.*`            | TBD    |

(Extras may be added if COMSOL converges — eg. mesh-convergence study, polarisation sweep.)

## Physics framing (for the tutor)

Hovhannisyan & Makaryan 2024 (paper 1, Armenian J. Phys.) experimentally observed THz absorption peaks near 0.6 and 2 THz in spinal-cord samples *only when* (a) the THz E-field is polarised parallel to the nerve fibres, and (b) a DC voltage > ~50 V/cm is applied perpendicular to the fibres. Paper 3 frames the spinal-cord slab as a **periodic diffraction grating of conducting wires** (each node of Ranvier becomes a nonlinear conductor under voltage).

This project models a **single fiber as the unit cell** of that periodic grating in 3D EM (COMSOL EWFD):

* Axon (water-like core, r = 5 µm, length 240 µm).
* Two myelin sheath segments (constant ε = 4.5 − 0.5 i, r = 5–7 µm, length 100 µm each).
* One node of Ranvier between the sheaths (water + iσ/(ωε₀), r = 5–7 µm, length 40 µm).
* External water medium filling the unit cell out to ±20 µm in x and y.
* Floquet-periodic BCs on the four lateral box walls (one fiber per unit cell of the grating).
* Scattering BCs on the inlet (z = 0) and outlet (z = L) endcaps, with x-polarised plane-wave background field `exp(-i·k₀·z) x̂`.

## Open questions to flag at the meeting

1. **Fiber orientation relative to k**. Issue #8 has fibers along z (parallel to k, perpendicular to E). Paper 1's resonance condition is E ∥ fibers, which would make fibers perpendicular to k. The current simulation models the **non-resonant** configuration; reorienting fibers along x is a one-day refactor.
2. **What is the σ at the node really?** The user-defined formula ε_node = ε_water + iσ/(ωε₀) treats the node as conductive water. Paper 3 invokes opening of ion channels under voltage. Sigma values to test: {0, 0.1, 0.5, 1, 5, 10} S/m.
3. **Lateral BCs**: Floquet-periodic chosen because of the diffraction-grating framing in paper 3. Alternative is Scattering / PML for an isolated-fibre interpretation.

## Sim narratives

Each sim has a corresponding `experiments/simN_*.py` script and a figure under `figures/out/`. Run with `uv run python experiments/simN_*.py`. Notes per sim are added as they complete.

(Filled in as work progresses below.)

### Sim 1 — Frequency sweep (baseline)
Status: pending.

### Sim 2 — Sigma sweep at 0.6 THz
Status: pending.

### Sim 3 — Node length sensitivity
Status: pending.

### Sim 4 — Material dispersion
Status: pending.

### Sim 5 — 3D field / geometry visualisation
Status: pending.

### Sim 6 — Diffraction-grating analytic
Status: pending.

## Risk log

| Risk | Mitigation |
|------|------------|
| EWFD solve doesn't converge | Each COMSOL-based sim has a pure-Python analytical fallback that produces real physics-grounded figures. |
| Geometry orientation is wrong (fiber ∥ k vs ∥ E) | Modeled the configuration as specified in issues #8/#10. Flagged in §"Open questions". Sim 6 (diffraction grating) is orientation-agnostic. |
| Mesh too coarse at 2 THz (λ_water ≈ 75 µm) | Mesh has 5 µm refinement near node; reduce `refinement_factor` to 0.5 if needed. |

## Run order

```
# Analytical / no-COMSOL (always work):
uv run python experiments/sim4_dispersion.py
uv run python experiments/sim6_grating.py
uv run python experiments/sim5_geometry_viz.py   # 3D matplotlib if COMSOL field unavailable

# COMSOL-dependent (only if periodic-BC fix lands):
uv run python experiments/sim1_freq_sweep.py
uv run python experiments/sim2_sigma_sweep.py
uv run python experiments/sim3_node_length.py
```
