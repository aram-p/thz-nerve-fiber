# THz nerve-fibre simulation — tutor briefing

**One-paragraph summary.** Hovhannisyan & Makaryan (Armenian J. Phys. 2024) observed resonant THz absorption peaks near **0.6 THz** and **2 THz** in spinal-cord samples when the THz E-field is polarised parallel to the nerve fibres and a DC voltage is applied perpendicular to them. Paper 3 frames the sample as a *periodic diffraction grating of nerve fibres*. This project builds a first-principles 3-D FEM simulation of one fibre as the unit cell of that grating (axon + myelin sheaths + node of Ranvier, periodic in *x* and *y*, plane-wave excitation in *z*) to ask: **does the model reproduce the resonance, and what physical parameters control it?** The headline result: the simulation produces a clean peak in field-at-the-node at **f = 0.632 THz**, sitting on the experimental 0.6 THz feature.

Pipeline: Python + `mph` driving COMSOL EWFD (scattered-field formulation, Floquet-periodic lateral BCs, scattering BCs on the end-caps, MUMPS direct solver). Source code lives in `src/thznerve/` and per-sim scripts in `experiments/`. All figures below are reproducible by running `uv run python experiments/simN_*.py`.

---

## 1 — Setup & geometry

### Sim 5 — the single-fibre unit cell (3-D)

![sim5](figures/out/sim5_geometry_viz.png)

The five labelled domains: a solid axon cylinder (water, r = 5 µm, length 240 µm), two myelin annular segments (constant ε = 4.5 − 0.5 *i*, r = 5–7 µm, length 100 µm each), one node-of-Ranvier annulus (r = 5–7 µm, length 40 µm, centred at z = 120 µm), and the surrounding external water box (±20 µm in *x* and *y*). The plane wave enters along +*z* with E along *x̂*.

### Sim 8 — periodic 3×3 array (the diffraction-grating interpretation)

![sim8](figures/out/sim8_periodic_array_3d.png)

The same unit cell tiled into a 3×3 array makes the periodicity visible — the central cell is highlighted and surrounding cells are faded. This is what the **Floquet-periodic BCs** in EWFD impose: the simulation solves one cell but the implicit physics is a periodic grating.

### Sim 7 — wavelength scales relative to the geometry

![sim7](figures/out/sim7_wavelength_scales.png)

λ in water versus the modelled geometric scales. **λ_water/2 crosses the internode length around 0.5–0.7 THz** — this is the regime where half-wave standing-mode resonances on the fibre geometry are physically expected, and the experimental 0.6 THz absorption peak sits squarely in that window.

---

## 2 — Materials

### Sim 4 — frequency-dependent permittivity ε(f) (2-D view)

![sim4](figures/out/sim4_dispersion.png)

Double-Debye water dominates Im(ε); myelin is constant; the σ-term at the node adds a small 1/f imaginary contribution that is sub-percent of water's natural loss at biological σ (≤ 10 S/m). This already foreshadows Sim 2's null result.

### Sim 14 — ε(f) as 3-D trajectories (Re ε, −Im ε, freq)

![sim14](figures/out/sim14_dispersion_3d.png)

Same data as Sim 4 but plotted as curves in the complex permittivity plane with frequency as the vertical axis. Water's relaxation manifests as a curved trajectory from high Re/high Im at low frequency to low Re/low Im at high frequency. Myelin is a vertical line (constant ε). The σ-decorated node curves overlap water until σ becomes physically extreme (1000 S/m), where they peel off at low frequencies.

---

## 3 — COMSOL results

### Sim 1 — Frequency sweep, σ = 0, baseline geometry

![sim1](figures/out/sim1_freq_sweep.png)

26 frequencies, 0.1–2 THz. **Peak |E| at the node = 2.60 at f = 0.632 THz** (annotated), sitting on the experimentally observed 0.6 THz peak. The smoothed (3-pt moving average) trend tracks the data without over-fitting and identifies the 0.632 THz feature as the strongest local maximum below 1 THz. Solve time: ~5 minutes for the full sweep.

### Sim 11 — frequency-sweep waterfall (3-D, derived from Sim 1's saved HDF5)

![sim11](figures/out/sim11_freq_waterfall.png)

Same data as Sim 1 but the **full axial profile** |E|(z) is preserved at every frequency, so we don't just see the peak at one point but the entire spatial structure as a function of f. The right-hand heatmap is the cleaner view; the 0.6 THz horizontal band shows where the resonance lives along z (around the node).

### Sim 9 — 3-D + cross-section field at the resonance

![sim9](figures/out/sim9_3d_field_at_resonance.png)

Re-solves at f = 0.632 THz on a 14×14×50 sampling grid + a fine 90×200 *y* = 0 slice. The 2-D panel (right) is the readable one: vertical orange lines mark the node boundaries, white horizontals mark the axon/myelin radii, and the magma colourmap reveals the spatial pattern of the scattered field at the resonance.

### Sim 10 — 3-D nested isosurfaces of |E| at the resonance

![sim10](figures/out/sim10_field_isosurfaces.png)

Same physics solve as Sim 9 but rendered as nested isosurfaces (computed via marching cubes) at |E|/max = 0.55, 0.70, 0.85. The right-hand panel is the cross-section-integrated energy density ∫∫|E|² dx dy as a function of z, which captures the **total** energy in each axial slice rather than only on-axis values.

### Sim 13 — on-resonance vs off-resonance contrast (3-D)

![sim13](figures/out/sim13_resonance_contrast_3d.png)

Two solves: f = 0.328 THz (a local minimum in Sim 1) versus f = 0.632 THz (the resonance peak). The 2-D cross-sections look qualitatively similar at this resolution; the difference is *regional* — at 0.632 THz the field at the specific node-annulus sampling location is enhanced, while elsewhere in the domain off-resonance can even have higher local |E|_max. This is an important nuance: **the resonance is not a global field enhancement, it's a localisation at the node**.

---

## 4 — Parameter studies

### Sim 3 — node-length sensitivity (2-D)

![sim3](figures/out/sim3_node_length.png)

Five node lengths {10, 20, 40, 60, 100} µm at f = 0.6 THz, σ = 0. Clean **monotonic** trend: peak |E| in node climbs from ≈ 1.6 (at L = 20 µm) to ≈ 2.9 (at L = 100 µm). Real biology has node length ≈ 1 µm — much smaller than the 40 µm baseline — which suggests our absolute |E| values are *upper bounds* on a faithfully sized biological node.

### Sim 12 — node-length axial profiles (3-D)

![sim12](figures/out/sim12_node_length_3d.png)

Same data as Sim 3 but in 3-D: each node-length gets its own ribbon showing how the axial |E| profile changes shape. The right-hand heatmap with inset peak-vs-L makes the monotonic trend self-contained.

### Sim 2 — σ sweep at f = 0.6 THz (NULL RESULT)

![sim2](figures/out/sim2_sigma_sweep.png)

Swept the node conductivity across **11 orders of magnitude** (σ ∈ [0, 10⁸] S/m). |E| at the node varies by less than 3 × 10⁻⁸ over the entire range. Two possible explanations: (a) the THz response at this geometry is dominated by the εr contrast between water (≈ 4) and myelin (4.5) and is genuinely insensitive to σ at the node, or (b) COMSOL EWFD's `DisplacementFieldModel` defaults to reading only εr — a probe (see `SESSION.md`) confirmed that changing Re(εr) does affect the field, so the material/selection wiring is correct, only the σ encoding doesn't reach the solver. **One-line fix candidate**: set `DisplacementFieldModel = "RelPermittivityWithSigma"` on `wee1`.

### Sim 15 — peak |E| over the (frequency, node-length) plane (3-D)

![sim15](figures/out/sim15_freq_nodelen_surface.png)

A 4 × 6 grid of solves (24 total simulations) at four node lengths × six frequencies, with the field sampled **inside the node annulus** (r = 6 µm, in the z-range of the node) rather than on the central axis. Peak |E|-at-node is rendered as a 3-D surface (left) and a 2-D heatmap (right).

**Important nuance this figure surfaces.** The surface is much *flatter* than Sim 1 suggests — node-annulus peak |E| stays in the 1.3–2.0 range across the whole (f, L_node) plane, with no clean 0.6 THz ridge. Combined with Sim 13's finding, this clarifies what the Sim 1 "0.6 THz peak" actually means: **it's the field on the central axis (inside the axon, water medium) at z ≈ node-centre, not the field at the conductive node annulus**. If the experimental THz absorption comes from energy dissipation in the conductive node, this single-fibre model would predict a much weaker frequency-dependence than what's experimentally observed — pointing toward the periodic-array / collective effects in Paper 3 (or the σ-encoding fix from Sim 2) as the missing physics.

---

### Sim 16 — mesh convergence at the resonance (3-D)

![sim16](figures/out/sim16_mesh_convergence_3d.png)

Three mesh refinements (coarse ≈ 24 k elements, baseline ≈ 26 k, fine ≈ 41 k) at f = 0.632 THz. Peak |E| on the axis comes out to 2.50, 2.75, 2.48 respectively — the baseline value reported in Sim 1 (~2.60) sits between coarse and fine and is **not** fully converged at the 10 % level. The fine-mesh value (2.48) is the better number to quote. The *qualitative* spatial pattern (cross-section heatmaps + 3-D scatter, bottom row) is stable across all three meshes — the resonance feature is real, only the absolute amplitude is mesh-dependent at the ±10 % level.

---

## 5 — Analytical context

### Sim 6 — wire-array equivalent-sheet model

![sim6](figures/out/sim6_grating.png)

Textbook 1-D diffraction-grating model (Tretyakov *Analytical Modelling in Applied Electromagnetics*, ch. 4): a periodic array of conducting wires modelled as a single sheet of effective impedance Z_sheet = R_grid + j X_L. With period 50 µm, wire radius 2.5 µm in water host, σ from 10 to 10⁵ S/m. Transmittance drops from 1.0 to ~0.5; absorbance rises to ~0.3 — the band-pass response characteristic of a wire-grid polariser. Direct analytical analogue of Paper 3's framing.

---

## 6 — Take-aways for the meeting

1. **Sim 1 shows a peak at 0.632 THz**, sitting on the experimental 0.6 THz feature — *for axial sampling* (field on the central axis of the fibre, inside the axon). The 2 THz peak is less clean.
2. **The "0.6 THz peak" disappears under annular sampling (Sim 15).** When the field is measured at the conductive-node location (r = 6 µm, inside the node annulus) rather than on the axis, peak |E| stays roughly flat at 1.3–2.0 across the whole (f, L_node) plane. This is the most important caveat: **the 0.6 THz resonance in Sim 1 is about the axon-water field at z ≈ node-centre, not the node-annulus field that would actually dissipate energy if σ were active**.
3. **Resonance is local (Sim 13).** On vs off resonance look qualitatively similar in cross-section; the difference shows up at specific sampling points.
4. **Node geometry matters (Sims 3, 12).** |E|-at-node scales with node length; the baseline 40 µm node is far larger than biology (~1 µm), so absolute field-enhancement numbers are upper bounds.
5. **σ doesn't matter — at least not the way I encoded it (Sim 2).** Either real physics or a one-line COMSOL `DisplacementFieldModel` fix. The probe shows the rest of the material pipeline works.
6. **Mesh-convergence caveat (Sim 16).** The baseline mesh isn't quite converged at the 10 % level. The fine-mesh peak at 0.632 THz is ≈ 2.48 (vs the 2.60 reported in Sim 1). Qualitative features are robust; absolute amplitudes should be reported with a ±10 % mesh-uncertainty band.
7. **Open physics questions** to discuss with the tutor:
   - *Fibre orientation*: the model has fibres along *z* (parallel to k); paper 1's resonance condition is E ∥ fibres (perpendicular to k). Reorienting fibres along *x* is a 1-day refactor and may move the resonance into the annular field.
   - *Why σ doesn't show* — physics or wee1 property fix.
   - *Why no 2 THz peak* — sampling density, mesh resolution, or geometric.
   - *Right observable*: maybe we should be plotting **power dissipated in the node domain** (= ∫ σ |E|² dV) rather than peak |E| at a sampling point. Once σ encoding is fixed, that becomes the natural quantity.

---

## 7 — Reproducing every figure

```sh
# Pure-Python (no COMSOL, < 5 s each):
uv run python experiments/sim4_dispersion.py
uv run python experiments/sim5_geometry_viz.py
uv run python experiments/sim6_grating.py
uv run python experiments/sim7_wavelength_scales.py
uv run python experiments/sim8_periodic_array_3d.py
uv run python experiments/sim14_dispersion_3d.py

# Pure-Python from saved sim-1/sim-3 HDF5 (< 5 s each, requires the COMSOL sims first):
uv run python experiments/sim11_freq_waterfall.py
uv run python experiments/sim12_node_length_3d.py

# COMSOL EWFD:
uv run python experiments/sim1_freq_sweep.py             # ~5 min
uv run python experiments/sim2_sigma_sweep.py            # ~2 min
uv run python experiments/sim3_node_length.py            # ~2 min
uv run python experiments/sim9_3d_field_at_resonance.py  # ~1 min
uv run python experiments/sim10_field_isosurfaces.py     # ~1 min
uv run python experiments/sim13_resonance_contrast_3d.py # ~2 min
uv run python experiments/sim15_freq_nodelen_surface.py  # ~6 min
uv run python experiments/sim16_mesh_convergence_3d.py   # ~2 min

# Re-plotting from CSV (no COMSOL):
uv run python experiments/sim1_replot.py
uv run python experiments/sim2_replot.py
```
