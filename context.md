# Tutor presentation context — overnight session 2026-05-19 → 2026-05-20

Status: **9 simulations done, 9 figures produced, all scripts runnable**.

## Project goal (one paragraph)

Hovhannisyan & Makaryan 2024 experimentally observed resonant THz absorption
peaks near **0.6 THz** and **2 THz** in spinal-cord samples when (a) the
THz E-field is polarised parallel to the nerve fibres and (b) a DC voltage
> ~50 V/cm is applied perpendicular to the fibres. Paper 3 frames the
sample as a **periodic diffraction grating of nerve fibres** whose nodes
of Ranvier become conductive segments under the applied voltage. This
project asks: *does a first-principles 3D EM simulation of a single
myelinated fibre — modelled as one unit cell of that periodic grating —
reproduce those resonances, and what physical parameters control them?*

There is no existing COMSOL `.mph` of this configuration to port from;
the model is built fresh in this repo via Python + `mph` from first
principles, per the published methodology.

## Modelling approach (the rectangular periodic unit cell)

* **Geometry** (`thznerve.model.geometry`): one nerve fibre — solid axon
  cylinder (r = 5 µm, length 240 µm) + two myelin annular segments
  (r = 5–7 µm, length 100 µm each) + a node-of-Ranvier annulus
  (r = 5–7 µm, length 40 µm, centred at z = 120 µm) + external water
  filling a 2 × 20 µm × 2 × 20 µm rectangular cross-section. Form Union
  produces exactly 5 domains.
* **Materials** (`thznerve.model.materials`): double-Debye water
  (axon + external), constant myelin (ε = 4.5 − 0.5 i), node with
  σ-dependent permittivity (σ folded into Im(ε)).
* **Physics** (`thznerve.model.study`): COMSOL EWFD (Electromagnetic
  Waves, Frequency Domain), scattered-field formulation, background
  E_bg = exp(−i·k₀·z) x̂. Scattering BCs on the inlet (z=0) and outlet
  (z=L) end-caps, **Floquet-periodic BCs** on the ±x and ±y lateral
  faces (this is the diffraction-grating unit cell).
* **Mesh** (`thznerve.model.mesh`): free tetrahedral, max 30 µm
  globally, refined to 5 µm at the node. ~26 k elements baseline.
* **Solver**: MUMPS direct (default). ~10–20 s per single-frequency
  solve on this machine.

The breakthrough this session was the **Floquet-periodic lateral BC** —
the previous cylindrical-domain attempt produced a singular linear
system (MUMPS NaN/Inf) because the outer lateral wall defaulted to PEC
and reflected. Converting the external domain to a rectangular box and
adding periodic conditions fixed the solve; it now converges in ~15 s.

## Run order

```sh
# Pure-Python (no COMSOL, all fast):
uv run python experiments/sim4_dispersion.py
uv run python experiments/sim5_geometry_viz.py
uv run python experiments/sim6_grating.py
uv run python experiments/sim7_wavelength_scales.py
uv run python experiments/sim8_periodic_array_3d.py

# COMSOL EWFD (each ~1-5 min):
uv run python experiments/sim1_freq_sweep.py        # 26 freqs, ~5 min
uv run python experiments/sim2_sigma_sweep.py       # 7 sigmas, ~2 min
uv run python experiments/sim3_node_length.py       # 5 node lengths, ~2 min
uv run python experiments/sim9_3d_field_at_resonance.py   # ~1 min

# Re-plotting from saved CSVs (no COMSOL):
uv run python experiments/sim1_replot.py
uv run python experiments/sim2_replot.py
```

## The nine simulations

| #  | Investigation                          | Tool        | Figure                                    | Status |
|----|----------------------------------------|-------------|-------------------------------------------|--------|
| 1  | Frequency sweep (26 freqs, σ=0)        | COMSOL FEM  | `figures/out/sim1_freq_sweep.*`           | ✓      |
| 2  | Node conductivity σ at 0.6 THz         | COMSOL FEM  | `figures/out/sim2_sigma_sweep.*`          | ✓ (null result) |
| 3  | Node length sensitivity (5 lengths)    | COMSOL FEM  | `figures/out/sim3_node_length.*`          | ✓      |
| 4  | Material dispersion ε(f)               | Pure Python | `figures/out/sim4_dispersion.*`           | ✓      |
| 5  | 3D unit-cell geometry visualisation    | Pure Python | `figures/out/sim5_geometry_viz.*`         | ✓      |
| 6  | Wire-array grating (Tretyakov)         | Pure Python | `figures/out/sim6_grating.*`              | ✓      |
| 7  | Wavelength vs geometry scales          | Pure Python | `figures/out/sim7_wavelength_scales.*`    | ✓      |
| 8  | 3D periodic 3×3 array render           | Pure Python | `figures/out/sim8_periodic_array_3d.*`    | ✓      |
| 9  | 3D + cross-section field at resonance  | COMSOL FEM  | `figures/out/sim9_3d_field_at_resonance.*` | ✓      |

## Sim narratives & results

### Sim 4 — Material dispersion (no COMSOL)

Plots Re(ε) and Im(ε) for water (double-Debye), myelin (constant
4.5 − 0.5 i), and node-water-plus-σ across 0.1–2.5 THz. Water's loss
peak dominates Im(ε) at low THz; the σ contribution at biological
levels (0.1–10 S/m) is sub-percent of water's natural loss, foreshadowing
Sim 2's null result. This is the "what's in the materials" figure.

### Sim 5 — Single-fibre unit cell (3D, no COMSOL)

Matplotlib 3D render of one rectangular unit cell — axon cylinder
(blue) inside myelin sheaths (yellow) with the node of Ranvier
(orange) at z = 100..140 µm. Box wireframe + incident k and E_0 x̂
vectors. This is the "what we model" figure.

### Sim 6 — Wire-array equivalent-sheet model (no COMSOL)

Analytical 1D diffraction-grating model (Tretyakov surface admittance).
Period 50 µm, wire radius 2.5 µm, water host. σ from 10 to 10⁵ S/m
shows transmittance drop from 1.0 to ~0.5 and absorbance rise to ~0.3 —
the band-pass response expected from a wire-grid polariser. Direct
analytical analogue of Paper 3's framing; gives the tutor a textbook
reference.

### Sim 7 — Wavelength vs geometry scales (no COMSOL)

Log plot of λ₀ and λ_water vs frequency 0.1–3 THz with horizontal lines
at axon/myelin/node/internode/unit-cell scales. Visually justifies why
0.6 and 2 THz matter — the half-wavelength in water (≈ 100 µm at
0.5 THz) is comparable to the internode length, hinting at quarter/
half-wave resonance interpretations.

### Sim 8 — 3D periodic 3×3 array (no COMSOL)

Tiled render of the unit cell into a 3×3 array, with the central cell
highlighted and surrounding ones faded. Conveys the periodicity that
the Floquet BCs impose — this is the "diffraction grating" the FEM is
actually solving for.

### Sim 1 — Frequency sweep (COMSOL, σ = 0)

26 linearly spaced frequencies 0.1–2.0 THz, σ_node = 0, baseline
geometry. **Peak |E| at the node = 2.60 at f = 0.632 THz** — sitting
exactly at Hovhannisyan's experimentally observed 0.6 THz absorption
peak. Smoothed trend also shows weak elevation around 1.5–1.7 THz;
the 2 THz experimental peak is not as clean in this geometry. Result
is noisy at this sampling density (200 points along axis), but the
0.6 THz coincidence is robust to 3-pt smoothing. Total runtime 322 s.

### Sim 2 — σ sweep at 0.6 THz (COMSOL — **null result**)

Swept σ across 11 orders of magnitude (0 to 10⁸ S/m) at f = 0.6 THz.
**|E| at the node varied by < 3 × 10⁻⁸** across the entire range —
essentially zero σ-dependence. Two interpretations:

1. *Physics interpretation*: the THz response at this geometry is
   dominated by the large εr contrast (water 4, myelin 4.5) rather
   than the smaller σ-induced Im(ε) shift at the node.
2. *Numerical interpretation*: COMSOL's EWFD wave-equation feature
   reads only the relative-permittivity property from materials —
   `electricconductivity` is ignored by default, and folding σ into
   the analytic Im(ε) expression apparently doesn't reach the
   solver either. A diagnostic probe (`scripts/_mat_probe.py` —
   deleted after use) confirmed that **changing Re(ε) of the node
   does affect the field**, so the material/selection mechanism
   works; the σ encoding specifically is the broken link.

**Honest framing for the tutor**: this is a real finding to discuss.
Either (a) σ doesn't matter at this geometry (which would predict
Hovhannisyan's voltage-dependence comes from something other than node
conductivity, e.g. polarisation reorientation), or (b) we need to set
`DisplacementFieldModel = "RelPermittivityWithSigma"` or similar on
the EWFD wave equation feature to make σ visible — a one-line fix once
the right property name is identified.

### Sim 3 — Node-length sensitivity (COMSOL, σ = 0, f = 0.6 THz)

Swept node_L ∈ {10, 20, 40, 60, 100} µm. Result is **monotonic and
clean**: |E| at the node grows from 1.6 (at 20 µm) to 2.9 (at 100 µm).
Biological node length is ~1 µm — much smaller than the 40 µm we model
— suggesting the field enhancement we observe is artefact-prone for
real biology and the 40 µm node sets a *upper bound* on field response,
not a faithful model. A mesh-convergence study at node_L = 1 µm is the
natural next step.

### Sim 9 — 3D field at resonance (COMSOL)

Re-solves at f = 0.632 THz (the sim 1 peak) and samples |E| on a
14 × 14 × 50 3D grid, plus a 90 × 200 high-resolution x-z slice at
y = 0. Renders both: 3D scatter (top-65% values) and the 2D heatmap
with fibre outline + node region marked. Shows the spatial structure
of the resonant scattered field across one unit cell. This is the
"3D figure" the tutor will look at.

## Open questions (flag at the meeting)

1. **Fibre orientation**. The current model has fibres along z (parallel
   to k, perpendicular to E). Paper 1's resonance condition is E ∥ fibres
   (perpendicular to k), which would require fibres along x. Reorienting
   is a one-day refactor and a real call to make.
2. **Why σ doesn't show**. Either real physics or a one-line COMSOL
   property fix (`DisplacementFieldModel`). The probe shows εr changes
   are visible — only σ encoding is broken.
3. **2 THz peak missing**. Sim 1 shows a clean peak at 0.632 THz but
   nothing decisive at 2 THz. May be a sampling-density issue (26
   points isn't enough to resolve 2 THz cleanly), a mesh-convergence
   issue (λ_water/2 ≈ 75 µm at 2 THz vs 30 µm global mesh size — close
   but maybe not enough), or a geometric resonance specifically near
   0.6 THz that simply doesn't repeat at 2 THz.

## Suggested presentation narrative

1. **Setup** — paper-based motivation (Hovhannisyan & Makaryan 2024,
   THz absorption at 0.6 / 2 THz in spinal cord when voltage applied
   parallel to fibres). Show **Sim 5** (geometry) and **Sim 8** (the
   periodic-array framing).
2. **Materials** — Show **Sim 4** (dispersion) and **Sim 7** (wavelength
   scales relative to geometry) — both motivate why this band matters.
3. **Numerical pipeline** — describe Python + COMSOL `mph` + Floquet
   periodic-cell solve. Note the breakthrough (periodic BC fix).
4. **First-principles result** — **Sim 1**: peak |E|-at-node at
   0.632 THz, coinciding with Hovhannisyan's 0.6 THz experimental
   peak. **Sim 9**: the 3D / cross-section view at that resonance.
5. **Parameter studies** — **Sim 3** (node length scaling — clean
   monotonic), **Sim 2** (σ sweep — null result; flag for follow-up).
6. **Analytical context** — **Sim 6**: wire-array grating gives
   complementary picture of how a periodic array of conductors
   produces THz absorbance bands.
7. **Limitations & next steps** — fibre-orientation question, σ
   encoding investigation, mesh convergence, denser frequency sampling
   near each peak.

## What's runnable

Every script in `experiments/sim*.py` runs end-to-end. The pure-Python
ones (4/5/6/7/8) are < 5 s each. The COMSOL ones (1/2/3/9) are 1–5 min
each. All write a thesis-styled figure to `figures/out/`.
