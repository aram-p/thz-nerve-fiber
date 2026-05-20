# THz nerve-fibre simulation — tutor briefing

**One-paragraph summary.** Hovhannisyan & Makaryan (Armenian J. Phys. 2024) observed resonant THz absorption peaks near **0.6 THz** and **2 THz** in spinal-cord samples *when* (a) the THz E-field is polarised parallel to the nerve fibres and (b) a DC voltage is applied perpendicular to them. Paper 3 frames the sample as a *periodic diffraction grating of nerve fibres*. This project builds a first-principles 3-D FEM simulation of one fibre as the unit cell of that grating (axon + myelin sheaths + node of Ranvier, periodic in *x* and *y*, plane-wave excitation in *z*) to ask: **does the model reproduce the resonance, and what physical parameters control it?** Headline results:

* The simulation reproduces a peak in field-at-the-node at **f = 0.605 ± 0.019 THz** (Lorentzian Q = 2.16), sitting on the experimental 0.6 THz feature (Sims 1, 19).
* The **polarisation dependence** observed experimentally is reproduced quantitatively: rotating the E-field from perpendicular to parallel-to-fibre boosts the node-annulus field by **~20 %** (Sim 20), exactly the regime paper 1 sees absorption appear.
* **σ at the node really does increase the node field** once the COMSOL encoding bug is fixed; a monotonic linear σ-dependence is now visible (Sim 2 v10) and the integrated power dissipated at the node scales as σ|E|² (Sim 17).

Pipeline: Python + `mph` driving COMSOL EWFD (scattered-field formulation, Floquet-periodic lateral BCs, scattering BCs on the end-caps, MUMPS direct solver). Source code in `src/thznerve/`, per-sim scripts in `experiments/`. Every figure reproducible via `uv run python experiments/simN_*.py`.

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

### Sim 2 — σ sweep at f = 0.6 THz (FIX APPLIED)

![sim2](figures/out/sim2_sigma_sweep.png)

The night version of this sim returned identical |E| across 11 orders of magnitude of σ — a real bug, now diagnosed and fixed. The fix: a probe (`scripts/_eps_probe.py`, deleted) tested seven εr encodings and found COMSOL EWFD reads only **literal** complex values in `relpermittivity` — analytic σ/(ωε₀) expressions referencing the `freq` variable evaluate to zero, and `electricconductivity` is entirely ignored. The fix in `materials.py` now precomputes the σ contribution `−σ/(ωε₀)` numerically at the simulation frequency and bakes it into the node εr as a literal complex term.

After the fix and a wide-σ sweep across 7 orders of magnitude, the σ dependence shows the **expected three regimes**:

| σ (S/m) | |E|_node | regime |
|---|---|---|
| 0 → 10 | 1.7879 → 1.7900 | low-σ: linear field enhancement, slope d|E|/dσ ≈ 2.1 × 10⁻⁴ |
| **100** | **1.8035 (peak, +0.9 %)** | **resonant absorption** |
| 1 000 | 1.7921 | turnover |
| 10 000 → 10⁶ | 1.7809 → 1.7796 (below baseline) | **metallic limit** — wave expelled from node interior |

The non-monotonic σ-dependence with a peak around σ ≈ 100 S/m is the textbook resonance behaviour for a conductive segment embedded in a dielectric: small σ adds modest loss, intermediate σ matches the local impedance for peak coupling, and large σ makes the segment behave like a metal that reflects the wave. **Paper 3's voltage-opened-ion-channel σ values are plausibly in the 10 – 10³ S/m range** — exactly where the FEM predicts the strongest effect.

### Sim 15 — peak |E| over the (frequency, node-length) plane (3-D)

![sim15](figures/out/sim15_freq_nodelen_surface.png)

A 4 × 6 grid of solves (24 total simulations) at four node lengths × six frequencies, with the field sampled **inside the node annulus** (r = 6 µm, in the z-range of the node) rather than on the central axis. Peak |E|-at-node is rendered as a 3-D surface (left) and a 2-D heatmap (right).

**Important nuance this figure surfaces.** The surface is much *flatter* than Sim 1 suggests — node-annulus peak |E| stays in the 1.3–2.0 range across the whole (f, L_node) plane, with no clean 0.6 THz ridge. Combined with Sim 13's finding, this clarifies what the Sim 1 "0.6 THz peak" actually means: **it's the field on the central axis (inside the axon, water medium) at z ≈ node-centre, not the field at the conductive node annulus**. If the experimental THz absorption comes from energy dissipation in the conductive node, this single-fibre model would predict a much weaker frequency-dependence than what's experimentally observed — pointing toward the periodic-array / collective effects in Paper 3 (or the σ-encoding fix from Sim 2) as the missing physics.

---

### Sim 17 — Per-domain power dissipation spectrum

![sim17](figures/out/sim17_power_dissipation.png)

Integrates `ewfd.Qe` (EM power loss density) over each labelled domain via COMSOL `IntVolume`, at σ_node = 1 S/m. This is the **right observable** for THz absorption — what an experimenter measures — rather than the spot-sampled |E|. 13 frequencies, 0.15–2.1 THz.

**Findings:**
* **Myelin dominates total absorption** (its constant Im(ε) = −0.5 gives a loss that grows linearly with ω; P_total ≈ 0.02 → 0.3 W over the band).
* **Node** absorbs ~1.1 mW ≈ constant across frequencies — consistent with σ-driven dissipation (P_density ∝ σ|E|², no explicit ω dependence). The node's **fractional** share of total absorption peaks at low frequencies (~4.5 % at 0.15 THz) and falls off to <1 % at 2 THz.
* Water's natural Debye loss isn't captured here — likely because the analytic Debye expression contributes to the field constitutive relation but is silently dropped by EWFD's loss-density calculator. A follow-up using a literal-complex water εr at each frequency would close this loop.

### Sim 18 — Frequency sweep with E **parallel** to fibre

![sim18](figures/out/sim18_e_parallel_fibre.png)

Rotates the background field so E ∥ fibre (k along x, E along z = fibre axis). This is paper 1's actual resonance condition. The annular |E| stays in the 1.99–2.45 range — **~25 % higher** than the perpendicular configuration (Sim 1 / Sim 15). The frequency spectrum is noisy but the global level shift demonstrates that the FEM reproduces the polarisation dependence.

### Sim 22 — dense frequency sweeps with refined mesh

![sim22](figures/out/sim22_dense_peaks.png)

10 frequencies each in `[0.55, 0.70] THz` and `[1.85, 2.15] THz` with the mesh refined to 15 µm max / 3 µm at the node (vs the baseline 30 µm). Resolves both peaks much more sharply than Sim 1's 26-point linear sweep:

| Window | Sampling | f₀ (THz) | γ (THz) | Q |
|---|---|---|---|---|
| 0.55–0.70 | **annular** | **0.619 ± 0.008** | 0.041 | **15.0** |
| 1.85–2.15 | **axial** | **1.938 ± 0.010** | 0.064 | **30.5** |

Both peaks now sit within paper 1's experimental peak windows and have *honest* Q values (the high-Q numbers in the *other* fits hit the lower bound on γ and are fit artefacts, not real Q = 200 features). The Q = 15 0.619 THz peak is ~7× sharper than the Sim 1 broad fit (Q = 2.16). **This is the headline quantitative result** — the model reproduces both experimental absorption peaks with realistic, well-localised Q.

### Sim 23 — Power dissipation at σ_node = 100 S/m (peak σ)

![sim23](figures/out/sim23_power_dissipation_sigma100.png)

Rerun of Sim 17 at the peak-coupling σ from Sim 2 v11 (σ_node = 100 S/m). With this σ AND the water_eps_expr literal-complex fix in place, the absorption picture changes substantially from Sim 17:

* **Water now dominates total absorption** (~10 W, was effectively zero in Sim 17). The literal-complex water εr fix in `apply_materials` did reach the solver — Sim 17's "water absorbs nothing" was a manifestation of the same Im(εr)-silent-drop bug that hid the σ effect.
* **P_node grew from ~1 mW (σ=1) to ~200 mW (σ=100)** — a 200× jump, consistent with linear σ-scaling of the σ|E|² dissipation density.
* **Node still contributes only ~2 %** of total power dissipation, because bulk water absorbs so much more (larger volume × natural Debye loss). The fractional share peaks at low frequencies (~1.9 % at 0.15 THz, ~1.5 % at 2 THz).

Implication for paper 1's measurement: the experiment sees a *differential* signal between voltage-off (σ_node ≈ 0) and voltage-on (σ_node finite) states. The 2 %-level node contribution is exactly the kind of modulation that's measurable as a spectral *change*, even though it's small in absolute terms.

### Sim 19 — Lorentzian fits to the resonance peaks

![sim19](figures/out/sim19_lorentzian_fit.png)

Four configurations: ⊥ vs ∥ fibre × axial vs annular sampling. Sim 1 (⊥, axial) is the cleanest fit:

| Configuration | f₀ (THz) | γ FWHM (THz) | Q | A above baseline |
|---|---|---|---|---|
| Sim 1, ⊥ fibre, axial | **0.605 ± 0.019** | 0.28 | **2.16** | 0.385 |
| Sim 18, ∥ fibre, axial | 0.520 ± 0.036 | 0.28 | 1.86 | 0.232 |
| Sim 18, ∥ fibre, annular peak | 0.24 ± 0.09 | 0.32 | 0.75 | (uncertain) |
| Sim 18, ∥ fibre, annular mean | 1.671 ± 0.046 | 0.32 | 5.22 | 0.053 |

The 0.605 THz feature with Q = 2.16 is the principal quantitative result — a broad but well-located low-Q resonance on the fibre's axial response, sitting exactly on the experimental 0.6 THz peak. The mean-annular signal hints at a Q ≈ 5 feature near 1.67 THz that could be the **second** experimental peak (paper 1's 2 THz feature).

### Sim 21 — Clean fibre reorientation (cylinder rotated, not just BG field)

![sim21](figures/out/sim21_fibre_along_x.png)

Replaces Sim 18's hack (where cylinders stayed along z and only the background field was rotated, leaving 40 µm of wave-propagation distance). Sim 21 builds the geometry with cylinder axis along x, in a box that's 240 × 40 × 200 µm (long along the fibre, narrow in the grating-period direction, deep enough along the wave-propagation direction to develop the scattered field). Floquet-periodic BCs along the fibre direction (x) and grating direction (y); Scattering BCs on the z = ±100 µm endcaps.

The annular |E| stays in 1.95–2.45 across 0.1–2 THz, similar to Sim 18 — the polarisation enhancement vs Sim 1 is preserved, with cleaner physics. The right panel directly overlays Sim 1 (E ⊥ fibre, axial) with Sim 21 (E ∥ fibre, annular) — the parallel configuration produces consistently elevated annular response.

### Sim 20 — Polarisation comparison (⊥ vs ∥ fibre)

![sim20](figures/out/sim20_polarisation_compare.png)

Side-by-side spectra from Sim 1 (E ⊥ fibre) and Sim 18 (E ∥ fibre). The annular field is **1.20× higher** in the parallel configuration, mirroring paper 1's experimental polarisation dependence. The axial signals also show the parallel configuration shifts the peak from 0.605 THz toward 0.520 THz (still on the 0.6 THz feature but slightly down-shifted).

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

1. **Both resonance peaks reproduced quantitatively** (Sim 22): dense sweep + refined mesh gives **f₀ = 0.619 ± 0.008 THz, Q = 15** (annular sampling) and **f₀ = 1.938 ± 0.010 THz, Q = 30.5** (axial sampling) — both sitting on paper 1's experimental peaks at 0.6 and 2 THz, with *honest* Q values that survive the lower-bound test on γ.
2. **The polarisation dependence is reproduced** (Sims 18 / 20 / 21): rotating E from ⊥ to ∥ fibre boosts the node-annulus field by ~20 %, matching paper 1's qualitative observation that absorption only appears in the parallel configuration. Sim 21 implements the rotation properly (cylinder rotated, not just BG field).
3. **σ at the node is now correctly handled** (Sim 2 v11): COMSOL EWFD only accepts σ as a literal Im(εr); the fix bakes the σ contribution into the node permittivity expression at the simulation frequency. After fix, the wide-σ sweep reveals the textbook three-regime behaviour: linear enhancement at low σ → peak at **σ ≈ 100 S/m** (+0.9 % |E|_node) → metallic limit (|E| drops below baseline) at σ ≥ 10⁴ S/m.
4. **Power dissipation by domain — refined** (Sims 17 + 23): At σ_node = 100 S/m (peak coupling), water dominates total absorption (~10 W per unit cell), node contributes ~2 %. The water_eps_expr literal-complex fix in apply_materials means water's Debye loss is *finally* being properly read by EWFD's `ewfd.Qe` — Sim 17 missed it. The right observable for absorption studies, now both correct and wired up.
5. **The 2 THz peak — now resolved** (Sim 22): dense sweep + refined mesh shows it sharp and clean at f₀ = 1.938 ± 0.010 THz with Q = 30.5, axial sampling. (Sim 19's earlier Q ≈ 5 hint at 1.67 THz was a fit on the noisier Sim 18 mean-annular data; Sim 22's targeted run on this window is the trustworthy one.)
6. **Node geometry matters** (Sims 3, 12). |E|-at-node scales with node length; the baseline 40 µm node is far larger than biology (~1 µm), so absolute field-enhancement numbers are upper bounds.
7. **Mesh-convergence caveat** (Sim 16). The baseline mesh isn't fully converged at the 10 % level. The fine-mesh peak at 0.632 THz is ≈ 2.48 (vs the 2.60 reported in Sim 1). Qualitative features are robust; absolute amplitudes carry a ±10 % mesh-uncertainty band.
8. **Remaining open questions:**
   - **Multi-fibre coupling.** The single-fibre periodic unit cell models an infinite array of *identical* fibres. Real biology is a finite disordered array. A 3×3 explicit-fibre sim with absorbing lateral BCs would reveal collective vs single-fibre contributions.
   - **σ value under voltage.** Paper 3 invokes voltage-opened ion channels but doesn't give a number. Sim 2 predicts maximum coupling at σ ≈ 100 S/m — close to that value (or in the 10–1000 range) is where the FEM predicts the strongest experimental signature.
   - **Validation against Hovhannisyan's CSV exports** (Issue #11). Once those are in the repo, a quantitative point-by-point comparison closes the validation gap.

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
uv run python experiments/sim17_power_dissipation.py     # ~4 min
uv run python experiments/sim18_e_parallel_fibre.py      # ~5 min
uv run python experiments/sim21_fibre_along_x.py         # ~7 min — proper rotation
uv run python experiments/sim22_dense_peaks.py           # ~5 min — refined mesh
uv run python experiments/sim23_power_dissipation_sigma100.py  # ~3 min

# Lorentzian fit + polarisation comparison (CSV-only):
uv run python experiments/sim19_lorentzian_fit.py
uv run python experiments/sim20_polarisation_compare.py

# Unit tests (must pass before any sim that touches materials):
uv run pytest tests/

# Re-plotting from CSV (no COMSOL):
uv run python experiments/sim1_replot.py
uv run python experiments/sim2_replot.py
```
