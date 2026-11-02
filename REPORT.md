# THz nerve-fibre simulation — report

## Phase 0 audit (2026-11-02): Sims 1–26 are invalid

**Verdict.** The 0.6 THz and 2 THz "resonances", their Q values, the half-wave mode
identification (Sim 24), the polarisation enhancement (Sims 18/20/21) and the σ results
(Sims 2, 17, 23) of the legacy 3D pipeline are numerical artefacts. None of them is evidence
about the fibre. Everything below the "Legacy report" line is kept only as a record.

Evidence (Sim 27, `experiments/audit/`, figure `paper/figures/fig_legacy_artefact.png`):

1. **Wrong background field.** The scattered-field background was `exp(-i*k0*z)` with the
   vacuum k0 in a water-filled box. That is not a solution in water, so every water element
   carries a volume source k0²(ε_w−1)E_b. The on-axis |E| at the node came out 1.7–2.5 V/m,
   where a plane wave entering water would have decayed to 0.3–0.6 V/m (panel a).
2. **The fibre has no effect on the reported metric.** In the solved legacy model, the on-axis
   |E| at the node centre is 1.6466 V/m with lipid myelin, with the myelin replaced by water,
   and with the node set to σ = 1e10 S/m (`scripts/probe/probe_legacy_rematerial.py`).
   Replacing the whole fibre by water reproduced config A to four digits. The material edits did
   reach the solver (the field 1 µm outside the axon changed), so the axis metric is blind to
   the fibre.
3. **The "spectrum" moves with the box.** Changing only the box/fibre length 240 → 300 µm moves
   the local maxima (0.5 → 0.6 THz and others); the 22-point spectra are uncorrelated noise
   between 1.5 and 2.6 V/m (panel a).
4. **Unresolved mesh.** The legacy |E|(z) on the axis jumps by 2× between neighbouring sample
   points (panel b): "peak |E| over the node window" was the maximum of mesh noise. Sim 22 also
   baked the water ε at 0.55 THz into a 1.85–2.15 THz sweep.
5. **Root cause: every legacy solve ran at 10²³–10²⁴ Hz.** The Frequency step defaulted to `punit = THz` and the code passed Hz numbers (`plist = "600000000000.0"` → 6e23 Hz, verified with `scripts/probe/probe_legacy_freq.py`). At that frequency the wavelength is far below the mesh, so the output is noise, and analytic σ/(ωε0) and Debye terms in `freq` collapse to 0 or ε∞. The "EWFD ignores σ / analytic εr" diagnosis was this unit bug and is withdrawn.
   The new pipeline puts all loss into complex εr and tests it against analytic results.
6. **Energy.** Sim 23 reported ~10 W absorbed per unit cell for a 1 V/m incident wave
   (incident power ~10⁻¹² W): energy balance was violated by ~12 orders of magnitude.

What replaces it: a benchmarked pipeline (`src/thznerve/physics`, `src/thznerve/fem`) with
the correct lossy background field, PMLs, a null test, a Mie/TMM benchmark and an energy balance
on every solve. See the iteration log below.

## Iteration log

| # | hypothesis | evidence | verdict |
|---|---|---|---|
| 27 | legacy peaks are artefacts | fibre ≡ no fibre on axis; peaks move with box length; mesh noise | **confirmed: Sims 1–26 invalid** |
| 28 | the new pipeline is exact on problems with known answers | B1 water slab vs transfer matrix: R+T+A−1 ≤ 8e-6, T error 9e-5 (30 freqs, 18 s). B2 2D layered cylinder vs Mie series, TM and TE: absorbed power per layer ≤ 8e-6, flux balance ≤ 5e-7. B3 axisymmetric harmonic model (m = 0..3) vs per-order series, both polarisations: ≤ 2e-4, balance ≤ 2e-5; null test (fibre = water) ≤ 1e-7. Bugs found on the way: flux on the PML interface reads half; axisym integrals need 2πr weighting; axisym PML must be Cylindrical (`experiments/bench/sim28_benchmarks.py`) | **passed: tools trusted** |
| 29 | reported voltage-induced features come from a partial air path, not fibre physics | authors' own n, d give notches at 0.52–0.58 THz (d = 100 µm) and 0.74/2.2 THz (d = 70 µm); the other 2π unwrap branch passes through their n = 7.8/4.2/2.9 points; high-f gain +9/+18/+17 dB vs +8/+11/+18 (`paper/figures/fig_twopath.png`) | consistent, hypothesis only (needs raw data) |
| 30 | H30: the node shows no resonance in 0.1–3 THz; its absorption follows the quasi-static law P = ½ωε0ε″_w|E0|²V_n·F | Realistic fibre (axon a = 3.5 µm, b = 5 µm, g = 0.7, node 1 µm, internode 1000 µm), cited materials, 30 frequencies, both polarisations, m = 0..3 (last order ≤ 5e-4), energy balance ≤ 2e-5. Node field factor η = ⟨|E|²⟩_node/|E0|²: E∥ 0.32 → 0.73, E⊥ 1.75 → 1.18, both smooth and monotonic. Independent quasi-static solver (no fitting): E∥ within 0.1–6 %, E⊥ within 0.5–0.8 % at every frequency. Node share of fibre absorption per period 0.03–0.12 %; node cross-section 0.1–0.8 µm² (E∥) | **confirmed: no resonance; mechanism quasi-static** (convergence: Sim 31) |
| 31 | Sim 30 is converged and independent of the box and boundary | Mesh, three levels (113k/267k/731k DOF, Lint = 200 µm; 4 frequencies 0.1–3 THz): P_node changes L1→L3 ≤ 0.51 % (E∥), ≤ 0.12 % (E⊥); L2→L3 ≤ 0.12 %; P_fibre < 1e-4. Internode 200 vs 1000 µm: η equal to 1.3e-5. Water radius 100 → 150 µm: 7e-5. PML 60 → 100 µm: 1e-6. First-order scattering BC instead of PML: E∥ 1e-4, E⊥ wrong by 2.9× at 0.1 THz at R = 100 µm and by 16 % (0.1 THz) and 5 % (0.3 THz) at R = 300 µm, E∥ exact: the first-order condition is exact only for m = 0 and errs as m²/(k_w R)² (k_w R = 0.7 → 2.1). The PML answer is the one that matched the exact series at 0.1 THz in benchmark B3 and the independent Laplace solver to 0.8 % (`experiments/prod/sim31_convergence.py`) | **passed: < 0.6 % mesh, < 1e-4 boundary (with PML)** |
| 32 | node absorption scales with geometry as the quasi-static law says; no size gives a resonance | Full wave at 0.1/0.6/1/2/3 THz, both polarisations, 10 geometries: node length 0.5–5 µm, myelin 0.5–3 µm, axon radius 1–7.5 µm (g = 0.7), internode 200–1500 µm. All spectra smooth. Full wave / extrapolated Laplace factor: 0.998–1.036 for fibres up to b = 6.5 µm; 1.089 at b = 10.7 µm and 3 THz, where k_w·b = 1.3 and retardation starts. Internode length changes P_node by < 1e-4. Closed law (three constants per polarisation, the transition length scaling with |ε_w/ε_my|, fitted to the Laplace solver only, max error 2.1 % against it; `experiments/prod/fit_nodelaw.py`) vs full wave: ≤ 2.2 % for b ≤ 5 µm, 4.3 % at b = 6.5 µm, 9.8 % at b = 10.7 µm and 3 THz. a = 1 µm at 1 mm internode needed ~1.2 M DOF (8 GB, paged, killed); run at 200 µm internode, justified by the internode independence (`experiments/prod/sim32_scaling.py`, `fig_scaling.py`, `paper/figures/fig_scaling.png`) | **confirmed: the law predicts the scaling; no geometric resonance** |
| 33 | (a) the result is robust to material data; (b) some node state could change THz absorption enough to explain a voltage effect | (a) ε′ and ε″ of axon, myelin and water (water includes the node) each scaled by 0.8 and 1.2: P_node changes −24 % to +34 %, smooth at all frequencies; the closed law predicts each change to ≤ 5 %. Node constriction (axon 3.5 → 2 µm inside the node): full wave vs Laplace solver ≤ 1.2 %. Axolemma (5 nm bilayer, analytic layered cylinder): ≤ 0.2 %. Dense packing (60 % area fraction, periodic Laplace cell): E⊥ shell factor −10 to −29 %, E∥ unchanged. (b) Node gap conductivity +10⁻⁴ … 3×10⁶ S/m, full wave, 11 values, both polarisations: P_node rises up to 43× (E⊥) and then falls as the conducting gap screens its own field. Largest change of a 100 µm slab (1 mm internodes): 0.14 dB (E⊥, 3 THz, 3×10⁴ S/m); for E∥ ≤ 0.002 dB at 0.6 THz and ≤ 0.006 dB at 2 THz. Mesh with the node refined 2×: ≤ 3e-4 at 10⁵–3×10⁶ S/m (skin depth ≥ 0.17 µm). The Laplace solver fails above ~10⁴ S/m (2–100× off at 10⁵), so only full wave is used there (`experiments/prod/sim33_robust_nodestate.py`) | **(a) robust; (b) killed: no node conductivity, physiological or not, gives more than 0.14 dB** |
| 34 | some DC-voltage-driven change of the tissue could produce the reported ~9 dB features at 0.6 and 2 THz with E ∥ fibres | Effective-medium slab (100 µm aligned white matter on PTFE, water fraction 0.77): E∥ −6.2/−9.6 dB, E⊥ −5.8/−9.2 dB at 0.6/2 THz, smooth; dichroism only −0.35 dB. Node conductivity: fully opened Na channels (~1e-4 S/m) give 1e-9 dB; no conductivity up to 3e6 S/m gives more than 0.006 dB for E∥ (Sim 33; an earlier linear estimate of "Δσ ≈ 5e5 S/m for 9 dB" assumed a fixed node field and is withdrawn). DC Joule heating by the bias itself: 2.6 MW/m³ (σ⊥ = 0.083 S/m, Ranck & BeMent 1965), 0.64 K/s adiabatic, enough to evaporate ~9 % of the slice water per minute. Heating +5/+15 K: −0.13/−0.43 dB at 0.6 THz. Water loss 10/30 %: +0.6/+1.5 dB at 0.6 THz, +0.9/+2.8 dB at 2 THz. Every mechanism is smooth in frequency; none produces a narrow feature (`experiments/exp/sim34_tissue_voltage.py`, `paper/figures/fig_tissue_voltage.png`) | **killed: no physiological node state or bulk DC effect gives 9 dB features**; heating and drying give smooth tilts below 3 dB but do not explain why only DC across the fibres acts (Sim 37) |
| 35 | the myelin sheath guides THz along the fibre and focuses it into the node (Liu 2019 picture) | THz index of myelin 1.67–2.25 is below water's 1.91–3.42 at every frequency in 0.1–3 THz (test `test_physics_analytic.py` asserts it; `paper/figures/fig_model.png` c): the sheath cannot guide by total internal reflection. The axoplasm (n 1.86–3.12) is above myelin, so axon plus sheath is a core-cladding pair, but its V number is 0.05–0.18 for a 3.5 µm axon (≤ 0.39 at 7.5 µm, 0.1–3 THz): the fundamental mode spreads far beyond the sheath into water, whose index is higher still, and leaks. Full wave: the node field factor for E∥ is 0.32–0.73 (< 1, the node is shielded, not fed); the node power does not depend on internode length (1.3e-5), so no guided energy arrives from the internode. Water absorption length 35–118 µm (power) is also shorter than any internode | **killed** for 0.1–3 THz (the mid-IR guiding claim is outside this band) |
| 36 | the node bound of Sims 33–34 survives the independent review: other observables, worst-case packing and internode, dense frequencies, thick fibres, lamellar myelin, broadband water, any passive node permittivity | `experiments/prod/sim36_review.py`; every solve energy-balanced ≤ 2.3e-4. (a) Node σ 0–3e6 S/m (12 values), three tissue observables per 100 µm slice: extra fibre absorption, extra absorption inside r = 25 µm (fibre plus the near-field water), and the change of the fibre's induced dipole (coherent transmission; dipole checked against the coated-cylinder value to ≤ 0.6 %). Largest change at internode 1000 µm, fibre packing 0.5: E∥ ≤ 0.006 dB at 0.6 and 2 THz (0.008 dB at 3 THz), E⊥ 0.18 dB (3 THz, 3e4 S/m). The slice effect scales as packing/internode, so internode 200 µm and packing 0.6 is 6× that: E∥ ≤ 0.05 dB, E⊥ ≤ 1.1 dB; the 50–100 µm internodes of the smallest CNS fibres add another 2–4×. Internode 1000 vs 200 µm at 3e4 S/m: node power and dipole change ≤ 7e-4. (b) 16 frequencies 0.1–3 THz at 3e4 and 1e6 S/m, b = 5 and 10.7 µm: smooth; the 10.7 µm fibre has broad maxima (E⊥ at 2.2 THz, width > 2 THz). Slice change at internode 200 µm: E∥ ≤ 0.03 dB; E⊥ ≤ 1.3 dB (3e4) and 3.5 dB (1e6, metallic node), about 0.5 dB at the realistic 1.5 mm internode of such fibres. (c) Lamellar myelin (uniaxial Wiener tensor, radial harmonic, axial arithmetic mean): node power +1–3 % above 0.4 THz, +26–29 % at 0.1 THz; with interlamellar water at half the bulk loss, −6 % to +33 %. Broadband water (Liebe 1991 Lorentz bands at 5.1 and 18.2 THz added to Kindt–Schmuttenmaer): node field factor η changes ≤ 6 %; node power follows water loss (×2 at 3 THz) as the law says. (d) Node with any passive ε = ε′ − jε″, ε′ −30…30, ε″ 0.01…30, at 0.6 and 2 THz: a negative-ε ring resonance near ε′ ≈ −1.5 raises η to ~100 for E∥, but the water around it caps the effect: ≤ 0.08 dB (≤ 0.5 dB at internode 200 µm); mesh check at the peak 3 % in η, 4 % in dB. At ε″ = 0.01 the contrast lies in the range where the gap's 90° corners are singular, so those points are not mesh-converged; the claim is made for ε″ ≥ 0.1. Negative ε′ with ε″ ≪ 1 is not possible in an aqueous gap at THz (ion collisions make ε″ ≫ |ε′|) | **bound holds**: no node state, passive or conducting, changes a 100 µm slice by more than ~0.1 dB for E∥ at 1 mm internodes, ≤ 0.5 dB at 200 µm, ≤ 0.8 dB after the dense ε grid of review round 5 (0.74 dB at 3 THz) (≤ 1 dB E⊥; 3.5 dB only for a metallic node in a 21 µm fibre with 200 µm internodes), against 9–20 dB reported |
| 37 | the experimenters' own model: DC opens node channels, fibres become (nonlinear) conducting wires and the slice a diffraction grating | Analytic (`experiments/exp/sim37_wire_hypothesis.py`, `results/sim37/sim37.txt`): at 0.6/2 THz the node membrane's capacitive admittance (ωC_m = 4e10–1.3e11 S/m²) exceeds that of all Na channels open (3e4 S/m², MRG) by 1.3e6–4.2e6, so channel state cannot change THz current paths; the axoplasm is continuous through the node in any state. Full wave, extra conductivity in the whole axoplasm (`sim36_review.py axon`; whole-axon mesh refined 4× changes results ≤ 7e-4 up to 1e5 S/m, `scripts/probe/probe_axon_mesh.py`), converted to an exact 100 µm slab on PTFE from the fibre dipole: 9 dB for E∥ needs +470–565 S/m at fibre packing 0.5 and 290–980 S/m over packing 0.7–0.3 (0.6–2 THz) (row 40), against 1.4 S/m physiological and ~25 S/m for saturated NaCl. Quasi-static Wiener mixing matches full wave to ≤ 3 % at ≤ 1 THz and ≤ 10 % at 3 THz up to 400 S/m and fails above ~1e4 S/m, where the axon's inductance matters (per-length crossover 0.9–2.9e4 S/m); so at the conductivity needed the fibres are resistive, not wires. Such a conducting axon would be E∥-selective (6.6 vs 2.7 dB at 0.6 THz), so polarisation alone does not discriminate. The tissue is post-mortem; the bias puts ~70 mV across a 10 µm fibre (not the ~1 V the authors estimate). Fibre spacing 2–23 µm: the first diffraction order needs > 6 THz in tissue, > 13 THz in air, and white matter is not periodic. Nonlinearity: Na gating ≥ 10 µs, ion drift in 1 ps ≤ 5e-16 m, TDS fields 1e3–1e4 V/m are 1e-3–1e-4 of the Wien-effect onset; the DC bias breaks inversion symmetry but χ(2)_eff = 3χ(3)E_dc gives ≤ 1e-11 of the linear response; their fitted quadratic current equals the linear one at 13 % of pulse peak. Side effects of the bias: ~0.7–1.1 mA through the slice, 2.6–3.7 MW/m³ Joule heat (0.6–0.9 K/s), and 25–42 mm³ of H₂ in 5 min at the electrodes, 7.5 mm from the beam (an air path through 20 % of the beam needs ~0.1 mm³, so the volume neither supports nor excludes one). But DC along the fibres (σ∥ ≈ 7σ⊥) heats more and reportedly does nothing, so heating alone does not explain the DC-direction selectivity. Reported baseline 30 dB per 0.1 mm at 2 THz is 1.9–2.6× pure water (15.8 / 11.6 dB, broadband / Kindt–Schmuttenmaer); nerve tissue absorbs within ~10 % of water (Berry 2003) | **killed**: open channels cannot make fibres conduct at THz; 9 dB needs 200–700× physiological axoplasm conductivity; no grating orders; no ps nonlinearity. DC-direction and polarisation selectivity remain unexplained |
| 38 | the reported spectrum, including the drifting dip, comes from a gas path through the slice (electrolysis gas or drying), not from tissue physics | `experiments/exp/sim38_gas_path.py`, `paper/figures/fig_gas_path.png`. Physical tissue optics (E∥ effective medium, n 1.9–2.3). Three independent reads of the reports give an optical thickness of 190–300 µm, not the stated 100 µm: the reported n = 3.6–3.9 at 100 µm equals the delay of physical n at 230–300 µm; the 2 THz loss (30 dB) equals 190–260 µm of water. Model: a fraction x of the beam crosses a gas column of thickness t. Fitted to 10 points read off the published 5-min spectrum (±2 dB): x = 0.11–0.22, t = 225–240 µm (water model), rms 3.1–3.5 dB. The fitted column puts dips at 0.55–0.64 THz and 1.98–2.16 THz, the two reported frequencies, and raises 2–3 THz by 7–24 dB (reported +12–17 dB). Growing t moves the first dip down (1.1–1.2 → 0.75–0.86 → 0.55–0.64 THz for t at 60/80/100 %), as reported (0.7 → 0.5 THz). The bias makes ~40–60 mm³ of electrolysis gas in 5 min (Sim 37), the order of the slice volume | **withdrawn** (review round 4): the points were read by eye and biased by up to 5 dB above 1.2 THz; on pixel-traced points the fit gives χ² = 32–52 on 8 degrees of freedom, and with the column held to the stated 100 µm it collapses (x → 0, rms 11 dB). There is no dip near 2 THz in this spectrum (red − blue rises monotonically to 2.5 THz). The thickness estimates are two, not three, and disagree. Replaced by the model-free waveform test of Sim 39 |
| 39 | the voltage-induced change in paper 3 comes from a beam component that does not cross the tissue; tested on the published waveforms with no tissue model | `experiments/exp/sim39_waveform_bypass.py`, `src/thznerve/digitise.py` (tested), `results/sim39/sim39.txt`, `paper/figures/fig_bypass.png`. Paper 3 page 2 rendered at 300 dpi; axes calibrated on gridlines (≤ 0.03 dB, ≤ 0.006 ps); curves traced by colour; tracing jitter 0.013 (trace units). The traced waveforms reproduce the published spectra within 1.5 dB up to 1.7 THz. (a) Biased = 0.587 × unbiased + 0.215 × empty cuvette (τ = −0.02 ps): residual 28 % of the trace rms; moving each trace by ±1 pixel (27 cases) gives 25–30 % with b = 0.214–0.216 and τ = −0.04…−0.01 ps, so the residual is ~5× the pixel-level tracing noise and the sum is not an exact description. (b) One-path changes of the tissue (scaled, shifted, thicker or thinner) leave 65–94 %: a component at the empty-cuvette delay is required (residual 31 % at τ = 0 vs 51–73 % at ±0.1 ps, reviewer). It may cross up to ~20 % of the tissue's optical thickness (29–32 % residual, not distinguishable from air). (c) With that component the tissue path must also change: amplitude 0.59 (a + b = 0.8), or ~11 % more optical thickness (0.69 × tissue path ×1.11 + 0.21 × empty: residual 20 %, dip 0.57 THz −10.7 dB vs reported 0.56 THz −11.4 dB; a fit, so a consistency check; depth ill-conditioned). Over 0.15–1.0 / 1.3 / 1.7 THz the same fit leaves 9 / 17 / 20 % with q = 1.18 / 1.13 / 1.11 and b = 0.212 / 0.206 / 0.207. Traced pulse energies: biased 1.03, no voltage 1.36 (−1.2 dB, the reported 'no significant change of average power'). (d) Fig. 1b alone: biased − empty cuvette is flat at −11.3…−13.2 dB over 1.5–3 THz, an undelayed component of 0.22–0.27 of the empty-cuvette field that dominates where tissue loss is 30–50 dB; model-free and matches (a). No dip near 2 THz. (e) Drift: a pure bypass gives a deep dip only at 0.74–0.76 THz (x 0.10–0.20); with the tissue path 10–20 % thicker the dip moves to 0.56–0.70 THz, so the reported 0.7 → 0.5 THz drift is compatible with a bypass plus a slowly thickening tissue path (qualitative). (f) The traced tissue phase gives n = 2.8–3.1 at the stated 100 µm (model 1.95–2.2) | **supported, cause unknown**: about a fifth of the beam reaches the detector at the empty-cuvette delay, and the tissue path changes too; together they give the dip's frequency and depth. What opens the path (gas, slice pulling away from a window, beam or sample movement) and what changes the tissue path are not identified. Contradicted by the reported E∥-only and DC-direction selectivity, if those controls were done as described. Paper 1: row 42 |
| 40 | bracket the conductivity a 9 dB fibre-conduction change needs | `sim36_review.py axon2`: full-wave points added at +600 and +800 S/m (energy balance ≤ 1.6e-5); exact 100 µm slab on PTFE from the fibre dipole, interpolated in log σ. 9 dB (E∥) needs 470–565 S/m at fibre packing 0.5 and 290–980 S/m over packing 0.7–0.3 (0.6–2 THz); 0.1 THz needs 600–1300 S/m. Physiological axoplasm 1.4 S/m; the axoplasm's own THz water loss ωε0ε″ is 82–115 S/m | **confirmed**: 200–700× physiological, 2.5–12× the water loss already present |
| 41 | the passive node-permittivity bound holds when its peak is sampled densely, meshed finer and the material data are moved ±20 % | `sim36_review.py epsdense / epsdense_mesh / epsrobust`, `results/sim30/sim36_epsdense.json`, `sim36_epsrobust.json`; energy balance ≤ 2.4e-5. (a) ε′ −2…−1 (10 values) at ε″ = 0.1, 0.6/1/2/3 THz: E∥ peak near ε′ −1.45…−1.5, 0.014 / 0.028 / 0.072 / 0.123 dB per 100 µm at internode 1000 µm, packing 0.5 (×6 for 200 µm, 0.6: 0.09 / 0.17 / 0.43 / 0.74 dB); E⊥ ≤ 0.018 (×6 0.11). (b) Three meshes at the peak (node element 0.25 / 0.125 / 0.0625 of the gap): 0.123 / 0.107 / 0.097 dB at 3 THz, 0.071 / 0.066 / 0.062 at 2 THz, within 3 % at 0.6–1 THz; falls monotonically (near-singular corners at this contrast), so the coarse value is an upper bound. (c) ε′, ε″ of axon, myelin, water each ×0.8 / ×1.2 at 3 THz, node ε′ rescanned: peak 0.087–0.148 dB, worst with myelin ε′ ×1.2 (peak moves to ε′ ≈ −1.7) | **bound holds**: ≤ 0.9 dB per 100 µm at internode ≥ 200 µm and packing ≤ 0.6 including ±20 % materials (0.74 dB nominal), 1.5–3.6 dB for the smallest CNS fibres, ≥ 2.5× below 9 dB |
| 42 | paper 1 (the source of the drift and the 2 THz feature) shows the same undelayed component as paper 3, growing with bias time | `experiments/exp/sim42_paper1_waveforms.py`, `results/sim42/sim42.txt`, `paper/figures/fig_paper1_bypass.png`; `thznerve.digitise.trace_track` (tested: follows a curve through crossing labels, stray pixels and stretches hidden under other curves). Paper 1 page 3, Figs. 3a (empty cuvette) and 4a (no voltage, 30 s and 4 min at 100 V), same axes and units: the no-voltage trace in both agrees to 0.017 rms. Spreads over ±1-pixel shifts. (a) 30 s: biased = 0.87 × no voltage + 0.04 (0.034–0.045) × empty, residual 10 % (one path 21 %); traced dip shallow (−4 dB) at 0.78 THz, where a pure bypass puts it. (b) 4 min: b = 0.17 (0.171–0.175) at τ = −0.04 ps (tissue delay ~0.6 ps); one-path fits 83–84 %; air-only 26 %; with the tissue path 1.09× optically thicker 16 %, dip 0.63 THz −20.5 dB vs traced 0.62 THz −19.1 dB. (c) Index from the traced phase at the stated 70 µm: 3.9 (model 1.95–2.2), as in paper 3. (d) The 2 THz feature is above what the 75 px/ps traces resolve and among the unpurged water-vapour lines of Fig. 4b: not tested | **supported**: the undelayed component replicates in the second report and grows 4 % → 17 % (→ 21 % at 5 min in paper 3); both reports need the tissue path ~10 % thicker (1.09, 1.11); the dip moves down as both grow. Fitted, not predicted. The 2 THz feature and the reported saturation after ~6 min are not tested. First pass ran on truncated traces (tracker lost two curves); caught on the overlay and rerun |


### Review after 5 iterations (Sims 27–31)
- Biggest time sink: Sim 30's first mesh (850k DOF, 5.4 GB, 17 min per polarisation), because the
  1000 µm myelin was meshed at the node resolution. Fix: refine only |z| < 10 µm around the node,
  coarser fibre mesh elsewhere (432k DOF, 3.3 GB, 104 s for 4 frequencies); same η to 4 digits.
- Repeated mistake: COMSOL defaults that silently change answers (`punit`, unrevolved axisym
  integrals, Poav averaged on interfaces). Each now has a benchmark that fails loudly (Sim 28), and
  the comsol skill lists them.


### Review after 5 iterations (Sims 32–35)
- Biggest time sink: meshes that silently exceeded RAM. Sim 32 t0 (1.26 M DOF) peaked at 7.8 GB and
  ran 3× slower; a0 paged and was killed. Fix: `axisym.check_memory` predicts DOF (7.0 × mesh
  elements) and peak memory from the mesh alone and refuses > 6 GB before solving.
- Repeated mistake: trusting a cheap model outside the range where it was checked. The linear Δσ
  estimate (Sim 34) and the Laplace solver above 10⁴ S/m both failed against full wave. Rule now:
  check the cheap model against full wave at the extreme of the range before sweeping with it.
- The closed law's first form had no committed fit script and a stated accuracy (7.6 %) that did
  not reproduce (9.7 %); the refit with a contrast-dependent transition is scripted and tested.


### Review after 5 iterations (Sims 36–40)
- Biggest time sink: Sim 38, a two-parameter tissue model fitted to points read off a figure by
  eye. The points were biased by up to 5 dB, the fitted column (225–240 µm) did not fit a 100 µm
  cuvette, and the claim cost a full review round before it was withdrawn. Fix: published curves
  are now traced from the PDF by `thznerve.digitise` (tested), and the published waveforms are tested
  directly (Sim 39) before any physical model is fitted to them.
- Repeated mistake: a claim stated more broadly than its evidence ("both dips", "reproduces"),
  then narrowed by the reviewer. Rule now: each claim about the experiment names the figure, the
  band where the comparison is valid, and the fraction of the effect it accounts for.
- Small losses: wrong output folder for a log file (0.2 min), a mesh check on the wrong region
  (10 min, Sim 36), a default COMSOL input that rejects tensors (2 min, in the comsol skill).

---

# Legacy report (Sims 1–26, invalidated by the Phase 0 audit above)


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

### Sim 24 — Mode identification (the peaks are half-wave standing modes)

![sim24](figures/out/sim24_mode_identification.png)

Predicts each peak frequency from a purely geometric half-wave-standing-mode argument: λ_water(f) / 2 = (segment length). Computed at the same Debye water n(f) used in the FEM, gives:

| Segment | Length (µm) | Predicted f (THz) | Observed f (THz) | Error |
|---|---|---|---|---|
| Node only | 40 | **2.043** | 1.938 | **+5.4 %** |
| Internode + half-node | 120 | **0.589** | 0.619 | **−4.8 %** |

Both peaks fall within 5 % of the half-wave prediction. **The two resonances are physically explained — not just observed — as half-wave standing modes on specific fibre segments**. The lower peak's *annular* sampling and the upper peak's *axial* sampling are also consistent with where each mode lives.

### Sim 25 — Visual confirmation of the two modes (3-D + cross-section)

![sim25](figures/out/sim25_field_at_both_peaks.png)

Solves at both peak frequencies (0.619 and 1.938 THz) with the refined mesh and renders |E| isosurfaces + y = 0 cross-section heatmaps side-by-side. The cross-section panel shows visibly **tighter wave fringes** at the high peak (λ_water/2 ≈ 38 µm at 1.94 THz) and a **broader spatial pattern** at the low peak (λ_water/2 ≈ 120 µm) — direct visual confirmation of the mode-identification claim.

### Sim 26 — Mesh convergence at both peaks (rigour)

![sim26](figures/out/sim26_convergence_both_peaks.png)

Three mesh refinements × two peak frequencies. Peak |E| varies by ~10 % across the coarse → baseline → fine mesh transition at both frequencies; the qualitative feature (peaks at 2.0–2.5 |E|) survives. **Quote peak amplitudes with a ±10 % mesh-uncertainty band**; peak *positions* are mesh-converged to better than 0.5 %.

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

1. **Both resonance peaks reproduced AND explained from first principles** (Sims 22 + 24): dense sweep + refined mesh gives **f₀ = 0.619 ± 0.008 THz, Q = 15** (annular sampling) and **f₀ = 1.938 ± 0.010 THz, Q = 30.5** (axial sampling) — both sitting on paper 1's experimental peaks at 0.6 and 2 THz. The peak frequencies are **predicted to within ±5 % by half-wave standing modes** on geometric segments: 40 µm node → 2.04 THz prediction; 120 µm internode + half-node → 0.59 THz prediction. Both peak amplitudes carry a ±10 % mesh-uncertainty band (Sim 26).
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
