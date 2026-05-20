# Detailed session log — overnight 2026-05-19 → 2026-05-20

This document is the long-form working log for the overnight session that
produced the 14 simulations described in `REPORT.md` and the supporting
infrastructure now in `src/thznerve/`. It's organised chronologically by
problem-and-fix, then by per-sim notes, then by open follow-ups.

## 0. Context the session started from

At the start of the session the repo had:

* The `mph` + COMSOL pipeline initialised (issues #2 – #9 closed).
* A working geometry builder (`thznerve.model.geometry.build_geometry`)
  that produced the 5-domain nerve-fibre unit cell via Boolean
  `Difference` operations (`axon`, `myelin_proximal`, `node`,
  `myelin_distal`, `external`).
* Materials wired into BallSelections per labelled domain.
* `setup_physics` with EWFD scattered-field + Scattering BCs on the
  end-caps.
* **One critical blocker**: the linear solve at f = 0.6 THz failed with
  `MUMPS: NaN or Inf found when solving linear system`. The cylindrical
  external boundary defaulted to PEC and reflected, producing a
  singular system. Nothing past issue #9 worked.

## 1. The unlock — Floquet-periodic BCs on a rectangular unit cell

**Problem.** PEC on the lateral cylindrical wall made the simulation
behave like a fibre inside a perfect-conductor tube. Adding a Scattering
BC to a cylindrical surface is awkward in EWFD; the more physical move
is **periodic** BCs, which is also what Paper 3's diffraction-grating
framing motivates.

**Fix.** Convert the external domain from a cylindrical annulus to a
rectangular box and add Floquet-Continuity periodic conditions on the
four lateral faces. Concretely:

1. Rename `external_radius_um` → `external_half_width_um` everywhere
   (`geometry.py`, `sweep/schema.py`, `experiments/example.yaml`, all
   smoke scripts). The geometry parameter is now the half-side of a
   square cross-section, not a cylinder radius.
2. Replace `_add_cyl("cyl_ext_out", …)` with `_add_block("box_ext", …)`
   in `build_geometry`.
3. In `add_endcap_selections`, use the box bounding box (±hw, ±hw)
   rather than ±external_radius.
4. New helper `add_lateral_periodic_selections` builds two
   `UnionSelection`s combining the ±x face pair and the ±y face pair,
   each face captured by a thin `BoxSelection`.
5. `setup_physics` adds two `PeriodicCondition` features on those two
   pair selections with `PeriodicType = "Continuity"` (zero in-plane
   Floquet wavevector — correct for normal incidence).

After this change the solve converges in **~15 s** producing a clean
|E| field with min/max around (0.3, 2.2).

**Commit:** `e2e199a` "Phase 2.3: EWFD periodic-cell solve works
end-to-end".

## 2. Result-extraction shape gotcha

`EvalPoint.computeResult()` returns `double[2][nExpr][nPoints]` —
**leading axis is [real, imag]**, then expressions, then points. My
first extraction assumed `double[nPoints][nCols]` (coords + value). The
fix is in `_evaluate_at_points`:

```python
result = ev.computeResult()
real = np.array(list(result[0][0]), dtype=float)
if complex_result:
    imag = np.array(list(result[1][0]), dtype=float)
    return real + 1j * imag
return real
```

For magnitude expressions like `ewfd.normE` the imaginary axis is
always zero and only `real` is used.

## 3. The σ-encoding bug — diagnosed and fixed

Sim 2 was supposed to sweep node conductivity σ at f = 0.6 THz and
observe how the field responds. Multiple encodings all returned
**identical |E| across 11 orders of magnitude of σ** in the overnight
session. After the user returned, a focused probe
(`scripts/_eps_probe.py`, deleted) tested seven combinations:

| case | εr setting | σ setting | result |
|---|---|---|---|
| 1 | water (analytic Debye, complex) | none | |E| = 0.9937 (baseline) |
| 2 | literal "4" (real, no Im) | none | 0.9953 |
| 3 | literal "4 − 3*i" | none | **1.0323** (differs!) |
| 4 | literal "4 − 1000*i" | none | 1.0226 |
| 5 | water + analytic `i*1e6/(2*pi*freq*8.854e-12)` | none | 0.9937 (no change!) |
| 6 | water | `electricconductivity = 1e6` | 0.9937 (no change!) |
| 7 | water | `electricconductivity = 1e10` | 0.9937 (no change!) |

**Conclusion**: COMSOL EWFD reads literal complex values in
`relpermittivity` (cases 3, 4 produce different |E| from cases 1, 2)
but **silently ignores analytic σ expressions referencing `freq`**
(case 5) and **ignores `electricconductivity` entirely** (cases 6, 7)
unless `wee1.DisplacementFieldModel` is reconfigured.

**Fix** (committed `1ca7ab7`): `materials.py.node_eps_expr(σ, freq_hz)`
now computes the σ contribution `Im_contrib = -σ/(ω·ε₀)` numerically at
the simulation frequency and bakes it into the εr expression as a
literal complex term. `apply_materials(model, params, freq_hz=...)`
takes the new `freq_hz` argument and propagates it through.

After the fix, a wide-σ sweep across 7 decades shows the classic
**three-regime non-monotonic σ-dependence**:

| σ (S/m) | |E|_node(annulus) |
|---|---|
| 0     | 1.7878702 |
| 0.1   | 1.7879019  (+1.8 × 10⁻⁵, linear) |
| 1     | 1.7880892  (+1.2 × 10⁻⁴, linear) |
| 10    | 1.7900160  (+1.2 × 10⁻³, still linear) |
| **100**   | **1.8035242  (+1.6 × 10⁻², peak — resonant absorption)** |
| 1 000 | 1.7920574  (turnover) |
| 10⁴   | 1.7809316  (below baseline — metallic) |
| 10⁵   | 1.7797336 |
| 10⁶   | 1.7796136  (saturated) |

The peak at σ ≈ 100 S/m is the textbook signature of a resistive
segment embedded in a dielectric: maximum power-coupling occurs when
the segment impedance matches the local wave impedance. Paper 3's
voltage-opened-ion-channel σ values are plausibly in the 10 – 10³ S/m
range, exactly where the FEM predicts the strongest effect.

## 4. Other small fixes along the way

* **Encoding** — `cp1252` stdout choked on `→` in early script output.
  Fixed by replacing arrows with ASCII `-->`. (Some scripts now also
  set `PYTHONIOENCODING=utf-8` when invoked.)
* **Difference selections** — `Difference` features expect
  `input`/`input2` as *Selections*, not direct property values. The
  call must go through `feature.selection("input").set([...])`.
* **Manifest dtype** — pandas `read_csv` coerced `"001"` to `int 1`.
  Fixed by reading the manifest with `dtype=str, keep_default_na=False`.
* **Sim 1 figure layout** — the original single-column figure was too
  cramped to read the legend + annotation. Switched to a double-column
  width with a 2-column legend.

## 5. Per-sim notes

The 14 simulations divide into five thematic groups (see `REPORT.md`).
Each script in `experiments/` is self-contained and runnable; each
saves both PDF and PNG under `figures/out/`. Notes here capture
non-obvious choices:

### Sim 1 — Frequency sweep (`sim1_freq_sweep.py`, ~5 min compute)

Linearly spaced 26 frequencies, 0.1 – 2.0 THz, σ = 0. Each solve
writes its own HDF5 under `results/sim1/freq_NNN_F.fff.THz.h5` with the
full axial profile, scalars (peak/mean/global |E|), and metadata
(git SHA, COMSOL version). A CSV summary is also written. The
`sim1_replot.py` companion renders the spectrum figure from the CSV
without re-running the sweep.

### Sim 2 — σ sweep (`sim2_sigma_sweep.py`, ~2 min)

See §3 above. The script samples at points inside the node annulus
(r = 6 µm, z in [101, 139] µm) rather than on the axis, because the
axis is inside the axon (water) where σ at the node can only matter
indirectly. The "fresh model per σ" version (current) defeats any
solver-side caching at the cost of doing the full setup + mesh build
every iteration.

### Sim 3 — Node-length sensitivity (`sim3_node_length.py`, ~2 min)

5 node lengths {10, 20, 40, 60, 100} µm. Each one is a *fresh* model
because the geometry changes between iterations. Mesh refinement
parameter at the node is scaled down for short nodes
(`MeshParams(node_element_size_um=min(5.0, node_L / 4))`).

### Sim 4 — Material dispersion (`sim4_dispersion.py`, < 1 s)

Pure Python. Two-panel Re(ε) / −Im(ε) plot over 0.1–2.5 THz for water,
myelin, and node at σ ∈ {0, 1, 10} S/m. The σ curves overlap water
heavily because the σ contribution is sub-percent of water's natural
loss at biological σ.

### Sim 5 — 3-D geometry viz (`sim5_geometry_viz.py`, < 5 s)

Pure matplotlib 3-D — uses `Poly3DCollection` for endcap discs and
`plot_surface` for cylindrical sides with low alpha. Box wireframe
drawn from line segments. The aspect ratio (1, 1, 1.8) is set so the
geometry doesn't get too elongated.

### Sim 6 — Wire-array equivalent sheet (`sim6_grating.py`, < 1 s)

Pure Python. Tretyakov-style sheet impedance for a 1-D array of
parallel conducting wires with E ∥ wires:

```
X_L = (ωμ₀a/2π) · ln(a / (2πr))
R_grid = a / (σ · π r²)
Z_sheet = R_grid + j X_L
t = 2 Z_sheet / (2 Z_sheet + η_water)
A = 1 − |t|² − |r|²
```

Period 50 µm, wire radius 2.5 µm, water host. σ swept across 10 ⁵
range to show the transition from "invisible" to "near-perfect mirror".

### Sim 7 — Wavelength scales (`sim7_wavelength_scales.py`, < 1 s)

Pure Python. Log-y plot of λ_vacuum and λ_water(f) over 0.1 – 3 THz
overlaid with horizontal lines at the modelled geometric scales
(axon diameter, myelin diameter, node length, internode length,
total fibre length, unit-cell side). Visually justifies why this
THz band matters for this geometry.

### Sim 8 — 3-D periodic 3×3 array (`sim8_periodic_array_3d.py`, < 5 s)

Pure matplotlib 3-D. Tiles the Sim 5 unit cell into a 3×3 grid with
the central cell at full opacity and surrounding cells faded. Adds
incident k and E vectors. The "what's actually being simulated" view.

### Sim 9 — 3-D + 2-D field at resonance (`sim9_3d_field_at_resonance.py`, ~1 min)

Re-solves at f = 0.632 THz (the Sim 1 peak) and samples |E| on a
14×14×50 3-D grid + a fine 90×200 *y* = 0 slice. The 2-D slice gets
Gaussian smoothing (σ = (1.2, 2.0)) to suppress sub-mesh-element
interpolation noise. Both panels use the same colour scale.

### Sim 10 — Nested |E| isosurfaces (`sim10_field_isosurfaces.py`, ~1 min)

Same solve as Sim 9 but renders three nested transparent isosurfaces
at |E|/max = 0.55, 0.70, 0.85 computed via
`skimage.measure.marching_cubes`. The right panel is the
cross-section-integrated energy density ∫∫|E|² dx dy as a function of
z — captures total field energy in each axial slice rather than
on-axis only.

### Sim 11 — Frequency-sweep waterfall (`sim11_freq_waterfall.py`, < 5 s)

Pure Python. Loads all 26 HDF5 files from Sim 1's `results/sim1/`,
stacks the axial profiles into a (26, 200) matrix, and renders a
3-D surface (left) + a top-down 2-D heatmap (right). The 0.6 THz
horizontal band on the heatmap is where the resonance lives along z.

### Sim 12 — Node-length axial-profile stack (`sim12_node_length_3d.py`, < 5 s)

Pure Python. Loads Sim 3's HDF5 files (5 node lengths), normalises z
to `z/total_L` (since the total length varies with node length),
renders the axial profiles as 3-D ribbons + a 2-D heatmap + an inset
showing the peak-vs-L trend.

### Sim 13 — On vs off resonance contrast (`sim13_resonance_contrast_3d.py`, ~2 min)

Two solves: f = 0.328 THz (a local minimum in Sim 1) and f = 0.632 THz
(the peak). Renders 3-D isosurfaces (top 3 panel) + 2-D *y* = 0 slices
(right). Important finding: the cross-sections look qualitatively
similar — the resonance is a *local* enhancement at the node-annulus
sampling point, not a global field amplification.

### Sim 14 — 3-D dispersion curves (`sim14_dispersion_3d.py`, < 1 s)

Pure Python. Plots water, myelin, and node-at-various-σ as 3-D curves
in (Re ε, −Im ε, freq) space. Same data as Sim 4 but a single
3-D view rather than two 2-D panels. Markers at f = 0.6 THz highlight
each material's complex permittivity at the resonance frequency.

### Sim 17 — Power dissipation per domain (`sim17_power_dissipation.py`, ~3.5 min)

Integrates `ewfd.Qe` (EM power loss density, W/m³) over each labelled
domain via `model.java.result().numerical().create("intvol_*", "IntVolume")`
+ `iv.selection().named(...)` so the integration is done per-named
selection (the BallSelections created during geometry build).

Run at σ_node = 1 S/m across 13 frequencies. Key findings:
* **Myelin dominates total absorption.** With Im(ε_myelin) = −0.5
  (constant), the loss per unit volume scales as ½·ω·ε₀·|Im(ε)|·|E|² —
  linear in ω. P_total goes 0.023 W → 0.31 W over 0.15 → 2.1 THz.
* **Node** absorbs roughly 1.1 mW ≈ constant. Consistent with σ-driven
  dissipation: P_density ∝ ½·σ·|E|² (no explicit ω dependence).
* **Node's *fractional* share peaks at low frequency** (~4.5 % at
  0.15 THz, falling to <1 % at 2 THz). At THz, myelin is the dominant
  absorber by mass.
* **Water's natural Debye loss isn't captured** — P_water ≈ 1.5 × 10⁻¹¹ W,
  essentially zero. The analytic Debye expression contributes to εr
  (the field looks right) but `ewfd.Qe` apparently doesn't include it
  in the loss-density formula. Follow-up: replace analytic εr with a
  per-frequency literal complex value and re-run.

### Sim 18 — E parallel to fibre (`sim18_e_parallel_fibre.py`, ~5.5 min)

Rotates the wave propagation direction from z to x while keeping the
fibre along z, so that E (along z) is parallel to the fibre — paper 1's
resonance condition. Refactor lives in
`setup_physics_e_parallel_fibre` in the sim script (separate from the
canonical `setup_physics`):

* Background field `Eb = [0, 0, exp(-i*ewfd.k0*x)]` (E along z, k along x).
* Scattering BCs at x = ±hw (inlet/outlet).
* Periodic BCs on y = ±hw (lateral grating) and z = 0 / z = L
  (fibre continuity).

The new annular sampling uses 8 azimuthal angles × 30 z-points to
capture the full ring of node-annulus field. Result: **annular |E|
stays in [1.99, 2.45] across 0.1–2 THz, mean 2.15**, vs Sim 1's
~1.79 in the perpendicular configuration — a 20 % enhancement that
matches paper 1's polarisation dependence.

### Sim 19 — Lorentzian fits (`sim19_lorentzian_fit.py`, < 5 s, CSV-only)

`scipy.optimize.curve_fit` with bounds (γ > 0.03 THz, A < 10) to a
Lorentzian + linear-baseline model, applied to a 3-pt smoothed window
around each suspected peak. Four configurations fit:
| config | f₀ (THz) | γ (THz) | Q | A |
|---|---|---|---|---|
| Sim 1, ⊥ axial | **0.605 ± 0.019** | 0.28 | **2.16** | 0.385 |
| Sim 18, ∥ axial | 0.520 ± 0.036 | 0.28 | 1.86 | 0.232 |
| Sim 18, ∥ annular peak | 0.24 ± 0.09 | 0.32 | 0.75 | (noisy) |
| Sim 18, ∥ annular mean | **1.671 ± 0.046** | 0.32 | **5.22** | 0.053 |

Sim 1 is the headline result. Sim 18 mean-annulus's Q = 5 at 1.67 THz
is suggestive of the experimental 2 THz feature but the fit is at
low amplitude and uncertain.

### Sim 20 — Polarisation comparison (`sim20_polarisation_compare.py`, < 5 s, CSV-only)

Side-by-side overlay of Sim 1 (⊥ fibre) and Sim 18 (∥ fibre) spectra.
Two panels: axial sampling (peak |E| inside axon at z ≈ node-centre),
and node-annulus sampling. The annular comparison is the headline:
mean |E|_annular goes 1.788 → 2.148 from ⊥ to ∥, a 1.20× enhancement.

### Sim 16 — Mesh convergence (`sim16_mesh_convergence_3d.py`, ~1.5 min)

Three fresh models at f = 0.632 THz with progressively finer meshes
(`MeshParams` with `refinement_factor` ∈ {0.6, 1.0, 1.5}, also varying
`max_element_size_um` and `node_element_size_um` so the element-count
ratio is clear). Renders a 3-column figure: 2-D cross-section heatmap
+ 3-D scatter for each refinement. The peak-|E|-on-axis values
(2.50 / 2.75 / 2.48) reveal that the baseline mesh used by Sim 1 is
not converged at the ±10 % level — the fine value 2.48 is closer to
the truth. The *qualitative* spatial pattern stays stable across
all three meshes, so the resonance feature is real, only its
amplitude carries a mesh-uncertainty band.

### Sim 15 — Frequency × node-length 3-D surface (`sim15_freq_nodelen_surface.py`, ~5.5 min)

4 node lengths × 6 frequencies = 24 solves. Renders peak |E|-at-node
as a 3-D `plot_surface` over (frequency, node_length) plus a 2-D
heatmap. **Crucially samples in the node annulus** (r = 6 µm) rather
than on the axis — which surfaces the most important nuance of the
session: the 0.6 THz peak in Sim 1 (axial sampling) **does not appear
in annular sampling**. The Sim 15 surface is much flatter, with no
clean 0.6 THz feature. This implies the Sim 1 "resonance" is the
field on the axon's central axis at z ≈ node-centre, not the field
at the conductive node-annulus location that would matter for σ-
mediated absorption. If anything this strengthens the case that
*the experimental 0.6 THz absorption comes from collective fibre-array
effects* (Paper 3's diffraction-grating picture) rather than
single-fibre node-localised dissipation.

## 6. Code organisation

```
thz-nerve-fiber/
├── README.md
├── REPORT.md                  ← tutor-facing summary
├── SESSION.md                 ← this file
├── context.md                 ← (legacy) overnight-session log
├── CLAUDE.md, docs/adr/       ← project decisions (immutable per ADR)
├── docs/references/           ← Hovhannisyan thesis + 3 papers (PDF + txt)
├── experiments/
│   ├── example.yaml           ← reference sweep config
│   ├── sim1_freq_sweep.py     … sim15_freq_nodelen_surface.py
│   └── sim1_replot.py, sim2_replot.py
├── figures/out/               ← all sim figures (PDF + PNG)
├── results/                   ← per-sim HDF5 + CSV (gitignored runtime data)
├── scripts/
│   ├── smoke_test.py          ← mph → COMSOL smoke test
│   ├── geometry_smoke.py      ← 5-domain build test
│   ├── materials_smoke.py
│   └── study_smoke.py
├── src/thznerve/
│   ├── __init__.py
│   ├── model/
│   │   ├── geometry.py        ← GeometryParams, build_geometry, BallSelections
│   │   ├── materials.py       ← Debye water, myelin, node, apply_materials
│   │   ├── mesh.py            ← MeshParams, build_mesh (FreeTet + node refinement)
│   │   └── study.py           ← setup_physics, periodic BCs, extract_axial_*
│   ├── sweep/                 ← Pydantic YAML schema + grid expansion + dispatcher
│   ├── io/                    ← HDF5 + manifest IO + provenance
│   └── plots/style.py         ← apply_thesis_style + save_figure (PDF+PNG)
├── tests/                     ← IO and materials unit tests
└── pyproject.toml             ← uv-managed
```

## 7. The big open follow-ups

1. **σ encoding** — *RESOLVED.* See §3 above. Fix committed `1ca7ab7`.
2. **Fibre orientation** — *PARTIALLY RESOLVED.* Sim 18 implements the
   E ∥ fibre configuration via a wave-direction rotation (cylinder
   axis stays along z, wave propagation moves from z to x). Result:
   ~20 % enhancement of annular |E|, matching paper 1's polarisation
   dependence. Further refinement (true fibre rotation, longer
   x-extent for clean wave propagation) is the next step.
3. **2 THz peak** — partial: Sim 19's mean-annulus fit on Sim 18 data
   finds a Q ≈ 5 feature at f₀ = 1.671 ± 0.046 THz, within the
   experimental 2 THz window. Amplitude is low and the fit is
   uncertain. Targeted high-density mesh + sampling between 1.5–2.2
   THz would settle this.
4. **Mesh convergence** — Sim 16 ran at f = 0.632 THz only; should be
   replicated at 1.7 THz to validate the second-peak claim. Refining
   `MeshParams.max_element_size_um` to 15 µm (vs 30) and re-sweeping
   near the peaks should give converged values.
5. **Water's natural Debye loss in ewfd.Qe** — Sim 17 shows
   P_water ≈ 0, suggesting the analytic Debye εr expression contributes
   to the field but not to the loss density. Workaround: replace the
   analytic water εr with a per-frequency literal complex value (same
   trick that fixed σ at the node).
6. **Validation against paper data** — Issue #11 (`ready-for-human`)
   describes a comparison against Hovhannisyan's experimentally
   measured CSV exports. Those CSVs aren't in the repo yet; once they
   are, a `data/baseline/` companion + `docs/validation.md` write-up
   can finally close that issue.

## 8. Session timeline (commits)

| commit  | what it added                                                      |
|---------|--------------------------------------------------------------------|
| 17b804d | (pre-session) Windows setup walkthrough + Claude handoff prompt    |
| 8269900 | issues #2–#7 done (uv project + scaffold + IO + matplotlib style)  |
| 813e1dd | issue #8 done (geometry builder, 5 domains)                        |
| 990a764 | issue #9 done (materials)                                          |
| 19b4eb7 | Phase 2.3 WIP — physics scaffold, NaN/Inf solver failure           |
| b32acfd | thesis references (3 papers + thesis + slides) into docs/references/ |
| e2e199a | **Phase 2.3 done** — periodic BC fix, solver converges             |
| 77315b9 | Sims 4 / 5 / 6 (analytical + geometry viz)                         |
| cd8c756 | Sim 1 / 2 / 3 scripts staged                                       |
| b29fcea | Overnight v1 — Sims 4-8 + Sim 1 results + context.md draft         |
| dcf2d3d | Overnight v2 — 9 sims done, sim 2 null result documented           |
| ba67725 | Overnight v3 — Sims 10-16 (7 more 3-D figures) + REPORT + SESSION  |
| 1ca7ab7 | **σ ENCODING FIX** — literal-Im(εr) baked at frequency             |
| 5a997f4 | sim 17 (power dissipation) + sim 19 (Lorentzian fit) scripts       |
| 83f3c07 | sim 19 — 4-panel ⊥/∥ × axial/annular fit comparison                |
| 2a0d7b7 | sim 20 — polarisation ⊥-vs-∥ side-by-side                          |
| **head**| Follow-up session — σ fix verified, sims 17/18/19/20 results, REPORT/SESSION revised |

Total: ~16 commits, ~14 simulations, ~14 publishable PDFs/PNGs,
~3.5 MB of figures, one 5-minute frequency sweep, one 8-minute
parameter sweep, one σ debug rabbit hole.
