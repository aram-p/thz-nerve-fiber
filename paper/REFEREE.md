# Referee report against our own claim, with responses

**Claim under review.** In 0.1–3 THz a single myelinated fibre in water has no resonant absorption
at the node of Ranvier. The node absorbs as a quasi-static, sub-wavelength gap: P_node =
½ωε0ε″_w|E0|²V_n·F, where F (0.2–2) depends only on node length / myelin thickness, the g-ratio and
the permittivities (F is defined by the absorbed power and computed from a Laplace problem). No node
state, and no fibre state the experimenters proposed (open channels turning fibres into nonlinear
conducting wires forming a grating), changes tissue transmission by an amount close to the reported
~9 dB, so the 0.6/2 THz features reported under DC bias (Hovhannisyan, Makaryan, Martirosyan 2024)
do not come from the linear electromagnetic response of nodes or fibres. The published waveforms
of the three reports, which share one set of reference scans, instead contain a component at the
empty-cuvette delay (4, 12, 17, 22 % after 0.5, 2, 4, 5 min), plus a ~10 % thicker tissue path
(Sims 39, 42, 43); the same two terms, unfitted to the spectra, place a transient dip at 2.24 THz
after 30 s (measured 2.00 THz); their cause is unknown and they conflict with the reported polarisation and bias-direction
selectivity.

Reviewer stance: a sceptical THz-biophysics and computational-EM referee looking for a fatal flaw.
Each objection is followed by the response and where the evidence is. "Fixed" means a simulation,
citation or narrower claim now answers it.

---

### 1. "A null result for a sub-wavelength gap in lossy water is obvious; nothing new."
**Response.** The opposite has been published: resonant or guided node behaviour in the THz/IR band
(Liu 2019; Liu 2020; Guo 2022; Song 2026), and resonant absorption measured at 0.6/2 THz
(Hovhannisyan & Makaryan 2024). None of these used measured dispersive water loss, a verified
solver and realistic internode lengths together. New here: (i) a benchmarked full-wave calculation
(transfer matrix, Mie series per azimuthal order, null test, energy balance on every solve);
(ii) a closed-form law for node absorption with its accuracy against full wave; (iii) a
quantitative bound from full wave: no node conductivity between 10⁻⁴ and 3 × 10⁶ S/m and no passive
node permittivity changes a 100 µm slab by more than the bounds of item 43 (round 4) for E∥, because a conducting gap screens its own field (item 16); fully
open Na channels (~10⁻⁴ S/m) give 10⁻⁹ dB; (iv) a test of the experimenters' own wire/grating model
(item 14). **Fixed** by the prior-work
table (`paper/PLAN.md`) and by narrowing the claim to 0.1–3 THz.

### 2. "Your materials are guesses: Looyenga mixtures for myelin and axoplasm, no measured THz data for myelin; molecular bands are missing from a Debye model."
**Response.** (a) Every material ε′ and ε″ was scaled by 0.8 and 1.2, one at a time
(Sim 33 robust): P_node moves by −24 % to +34 %, smoothly at every frequency; no feature
appears, and the closed law predicts each change to ≤ 5 %. (b) The bound in (1) covers any
node-localised loss, resonant or not: adding Δε″ = Δσ/(ωε0) up to ~10⁵ (Δσ = 3 × 10⁶ S/m) to the
gap, far beyond any molecular band (bulk water's ε″ is 2–7), never gives more than 0.1 dB for E∥ (Sim 36). A node
volume fraction of 3 × 10⁻⁴ cannot carry a 9 dB change. (c) Bulk nerve tissue measured by THz-TPI absorbs like
water (Berry 2003: 246 vs 225 cm⁻¹, featureless 0.5–1.5 THz). **Fixed.**

### 3. "The node gap is not water: there is the axolemma, a constricted axon, perinodal matrix and glial processes."
**Response.** Membrane: a 5 nm bilayer (ε = 2.3) at the node changes the gap field by ≤ 0.04 % and
the shell field by ≤ 0.2 % (exact layered-cylinder series, `thznerve.physics.packing`). Axon
constriction to 57 % of the internodal radius: smooth spectra, and the Laplace solver with the
constriction matches full wave to ≤ 1.2 %, so the law's machinery covers it. Gap contents: any change
of the gap material is covered by the conductivity sweep 10⁻⁴–10⁶ S/m (Sim 33 state) and by (2b).
**Fixed.**

### 4. "One fibre in infinite water is not white matter; dense packing changes local fields."
**Response.** For E ∥ fibres the field in the myelin is the macroscopic tangential field whatever
the packing, so the thin-gap limit F0 = |ε_my/ε_w|² is unchanged. For E ⊥ fibres a periodic Laplace
cell at 60 % area fraction lowers the shell field factor by 10–29 % (0.1–3 THz; Richardson-
extrapolated finite volume, tested against the isolated series in the dilute limit). Along one fibre,
nodes do not interact (internode 200–1500 µm changes P_node by < 10⁻⁴). These are O(1)
corrections; the gap to the reported effect is 10⁹. The tissue-level transmission uses an effective
medium (Wiener / Maxwell Garnett). **Fixed** (limitation stated).

### 5. "Mesh, truncation and boundaries could hide a narrow resonance."
**Response.** 30 frequencies 0.1–3 THz on the base fibre plus 5 frequencies on 10 geometries. Mesh
L1→L3 (113k→731k DOF): P_node ≤ 0.51 %. Azimuthal orders m ≤ 3: last order ≤ 1.5 × 10⁻³ of the
fibre absorption. PML 60→100 µm: 10⁻⁶; water radius 100→150 µm: 7 × 10⁻⁵. A first-order scattering
boundary instead of the PML is wrong for E ⊥ near the fibre (2.9× at R = 100 µm, 16 % at
R = 300 µm, 0.1 THz) as expected from its m²/(k_wR)² error; the PML matches the exact series
(benchmark B3) at the same frequency. Energy balance on every solve ≤ 3 × 10⁻⁵. A narrow resonance
missed between 30 points would need Q > 30 in a medium whose loss tangent is ~0.5–2: water damps
any mode in 35–118 µm (power). **Fixed.**

### 6. "The closed law is curve fitting."
**Response.** Both limits are exact (thin gap: normal-D continuity for E ∥, the coated-cylinder
shell series for E ⊥; long gap: bulk water, F = 1). Three constants per polarisation (two shape constants and the exponent q with which the
transition length scales as |ε_w/ε_my|^q) are fitted to the Laplace solver over x = 0.2–8,
g = 0.6–0.8, 0.1–3 THz, not to full wave and not to any experiment;
max error 2.1 % against the Laplace solver, tested off the fit grid. Against full wave the law
holds to ≤ 2.2 % for fibre radius ≤ 5 µm, 4.3 % at 6.5 µm and 9.8 % at 10.7 µm and 3 THz, where
k_w·b = 1.3 and retardation begins. It also predicts the ±20 % material changes to ≤ 5 %. The Laplace solver itself
matches full wave to 0.2–3.6 % (8.9 % at b = 10.7 µm, 3 THz). **Fixed** (accuracy stated).

### 7. "Normal incidence with k ⊥ fibre is a special case."
**Response.** It is the experimental geometry (slices with fibres in the plane, THz at normal
incidence). The node response depends only on the local field, and no length in the problem
approaches the wavelength in water (k_w·Ln ≤ 0.12), so oblique incidence cannot create a
resonance. **Narrowed** to the experimental geometry.

### 8. "You cannot exclude active or nonlinear processes (e.g. the gain proposed by Liu 2020)."
**Response.** Correct. The model is linear and passive; THz-TDS fields are weak (no heating, no
nonlinearity). The claim is restricted to passive absorption. **Narrowed.**

### 9. "Your experimental comparison is speculative: you have no raw data."
**Response.** Agreed for the instrumental explanations, which are presented as testable hypotheses:
(a) water-vapour lines at 0.557 and 1.92–2.20 THz; (b) a partial air path in the sample beam
(Sim 29); (c) DC Joule heating: 100 V over 18 mm through white matter (σ⊥ = 0.083 S/m, Ranck & BeMent
1965) deposits 2.6 MW/m³, 0.64 K/s without cooling, enough to evaporate ~9 % of the slice water per
minute or raise it by ~27 K at steady state; the timing is compatible with the reported ~5 min onset
and ~1 h recovery. Heating and drying change 100 µm slab transmission smoothly by ≤ 3 dB (Sim 34),
and heating cannot explain why only DC across the fibres acts (item 15). What is not hypothetical:
no node or fibre mechanism reaches the reported magnitude. **Narrowed**; the paper states
what matches and what does not (below).

### 10. "The experiment shows a polarisation dependence; your model does not explain it."
**Response.** Correct, and stated. The model predicts weak dichroism of an aligned slab
(E∥ − E⊥ = −0.35 dB at 0.6 and 2 THz) and no polarisation-selective node effect. DC heating
depends on the current direction relative to the fibres (σ∥/σ⊥ ≈ 7), not on the THz polarisation.
The reported polarisation selectivity is therefore **unexplained** by any mechanism we tested.

### 11. "The reported tissue index n = 3.6–3.9 disagrees with yours; maybe your tissue model is wrong."
**Response.** Water-rich tissue cannot have n above water's (1.9–2.3 above 0.5 THz). Our effective
medium gives n ≈ 1.95–2.2, matching Berry 2003 (nerve 1.95 ± 0.46). The reported values are
inconsistent with any water-based tissue, which supports an extraction artefact (Sim 29). Our
tissue absorption is ~25 % below Berry's nerve value; the node bound is insensitive to this.
**Fixed** (limitation stated).

### 12. "Temperature: you used 20–25 °C water; tissue in vivo is at 37 °C."
**Response.** Liebe water at 37 °C changes α(1 THz) by 2 % and n(0.5 THz) by 3 %, inside the ±20 %
robustness test. **Fixed.**

### 13. "Your own tissue-level estimate used a linear Δσ formula; is the node bound an extrapolation?"
**Response.** It was, and it was wrong: the linear formula assumed a fixed node field and gave
"Δσ ≈ 5 × 10⁵ S/m for 9 dB". The full-wave sweep (11 conductivities, mesh-checked at the largest)
shows the field is screened and the effect saturates (Sim 33; three observables in Sim 36, item 16). The bound now rests on full wave
only; the linear helper is documented as valid for small Δσ only. **Fixed.**

---

## Second round (independent reviewer, 2026-11-02)

### 14. FATAL: "You refute a hypothesis the experimenters never proposed. They propose that open node channels turn fibres into wires with nonlinear conductivity and the slice into a diffraction grating."
**Response.** Correct; the paper is reframed and tests both hypotheses (abstract, §1, §5). Their
model fails on four independent counts (Sim 37; full wave in Sim 36 'axon'):
(a) At 0.6/2 THz the node membrane's capacitive admittance (ωC_m = 4 × 10¹⁰–1.3 × 10¹¹ S/m²)
exceeds that of all Na channels open (3 × 10⁴ S/m², MRG) by 10⁶. Channel state cannot change THz
current paths, and the axoplasm is continuous through the node in any state.
(b) Full wave with extra conductivity in the whole axoplasm, converted to an exact 100 µm slab:
9 dB for E∥ needs +470–565 S/m at packing 0.5 (290–980 S/m over 0.3–0.7; item 44), against 1.4 S/m physiological and ~25 S/m for saturated NaCl.
Quasi-static mixing matches full wave to 3 % (≤ 1 THz) and 10 % (3 THz) up to 400 S/m; fibres act
as wires (inductive regime) only above ~10⁴ S/m. (A conducting axon would be E∥-selective at
0.6 THz, so polarisation does not discriminate; item 25.)
(c) Fibre spacing is 2–23 µm; the first diffraction order needs > 6 THz in tissue, > 13 THz in air,
and white matter is not periodic.
(d) No ps-fast nonlinearity: Na gating ≥ 10 µs; ion drift during the pulse ≤ 5 × 10⁻¹⁶ m; TDS fields
(10³–10⁴ V/m) are 10⁻³–10⁻⁴ of the field where electrolytes turn nonlinear; their fitted quadratic
current equals the linear one at 13 % of pulse peak. Their filter-before/after test does not prove
generation: a linear cascade is order-independent only without reflections between the elements and
with identical beam paths, and the sample changes between the two scans. **Fixed.**

### 15. MAJOR: "Joule heating predicts the opposite DC-direction dependence."
**Response.** Correct. Current along the fibres (σ∥ ≈ 7σ⊥) would heat ~7× more, yet the first report
sees no change for that direction. The paper no longer offers heating as the explanation. It
quantifies the bias current's side effects (~1 mA, 2.6–3.7 MW/m³, electrolysis gas comparable to the
sample volume in 5 min) and lists measurements that would decide (current and mass logging, N₂
purge, reference through a heated or dried slice). DC-direction and polarisation selectivity are
stated as **unexplained**.

### 16. MAJOR: "The bound depends on internode and packing, and uses absorbed power only; near-field changes of the water's own absorption and coherent changes are missed."
**Response.** Sim 36 repeats the node σ sweep (12 values, 0–3 × 10⁶ S/m) with three observables:
fibre absorption, absorption inside r = 25 µm (fibre plus near-field water), and the change of the
fibre's induced dipole, which gives the coherent change of the slice permittivity (dipole checked
against the coated-cylinder value to ≤ 0.6 %). At internode 1000 µm and packing 0.5: E∥ ≤ 0.006 dB at
0.6/2 THz (0.008 dB at 3 THz), E⊥ ≤ 0.18 dB. The slice effect scales as packing/internode; 200 µm
and packing 0.6 is 6× that: E∥ ≤ 0.05 dB, E⊥ ≤ 1.1 dB. Internode independence rechecked at
3 × 10⁴ S/m: ≤ 7 × 10⁻⁴. **Fixed.**

### 17. MAJOR: "'No resonance' covers only conductors; a node with a resonant constituent could resonate."
**Response.** Sim 36 sweeps any passive node permittivity, ε′ −30…30, ε″ 0.01…30, at 0.6 and 2 THz.
A negative-ε ring resonance exists near ε′ ≈ −1.5 (E∥, η ≈ 100), but the lossy water around the gap
caps the slice change at ≤ 0.08 dB (≤ 0.5 dB at internode 200 µm; dense grid in item 51: ≤ 0.8 dB; with mesh and ±20 % materials, item 59: ≤ 0.9 dB); mesh check 4 %. Negative ε′ with
ε″ ≪ 1 is not possible in an aqueous gap at THz. **Fixed**; claim stated for passive media with
ε″ ≥ 0.1 (item 31).

### 18. MAJOR: "The reported baseline exceeds pure water; the 9 dB may sit at the dynamic-range limit."
**Response.** Agreed and stated: 30 dB per 0.1 mm at 2 THz is 1.9–2.6× pure water (15.8 / 11.6 dB). Our
comparisons use changes, not absolute values; the raw data (IEEE DataPort, doi:10.21227/h17h-tp94)
would settle it. **Narrowed.**

### 19. MAJOR: "Isotropic Looyenga myelin is wrong; myelin is lamellar, and confined water differs from bulk."
**Response.** Uniaxial Wiener tensor (radial harmonic, axial and azimuthal arithmetic mean of water
and lipid): node power +1–3 % above 0.4 THz, +26–29 % at 0.1 THz; with interlamellar water at half
the bulk loss, −6 % to +33 %. Spectra stay smooth. **Fixed.**

### 20. MINOR: "Kindt–Schmuttenmaer water is used beyond its 2 THz fit."
**Response.** Adding Liebe's broadband Lorentz bands (5.1 and 18.2 THz) doubles the water loss at
3 THz; the node field factor changes ≤ 6 % and node power follows water loss, as the law says.
**Fixed.**

### 21. MINOR: "Thick fibres and high σ need a dense frequency sweep."
**Response.** 16 frequencies at 3 × 10⁴ and 10⁶ S/m, b = 5 and 10.7 µm: smooth; the 10.7 µm fibre
has only broad maxima (width > 2 THz). Slice change E∥ ≤ 0.03 dB; E⊥ up to 3.5 dB only with a
metallic node and a 200 µm internode, about 0.5 dB at the realistic 1.5 mm. **Fixed.**

### 22. MINOR: wording.
**Response.** Waveguide argument corrected: axoplasm (n 1.9–3.1) is above myelin, but the axon-core
V number is ≤ 0.4 and water's index is higher still, so any mode leaks (REPORT row 35, Fig. 1c).
Electrode gaps stated per report (18 mm, 15 mm). The E∥ bound is quoted in the abstract. **Fixed.**

---

## Third round (independent reviewer, 2026-11-02; verdict before fixes: major revision, no fatal flaw; key numbers recomputed and confirmed)

### 23. MAJOR: "The abstract overstates the bounds (0.1 dB vs 0.5 dB at 200 µm; E⊥ up to 3.5 dB)."
**Response.** Abstract now quotes ≤ 0.5 dB for E∥ (all node states, internode ≥ 200 µm); REPORT
row 36 gives every case: 0.1 dB at 1 mm, 0.5 dB at 200 µm, E⊥ ≤ 1 dB, 3.5 dB only for a metallic
node in a 21 µm fibre with 200 µm internodes. The smallest CNS fibres (50–100 µm internodes) add
2–4×. (Superseded by item 43.) **Fixed.**

### 24. MAJOR: "'Does not come from nodes or fibres' is broader than a linear passive model can show."
**Response.** Narrowed to "the linear electromagnetic response of nodes and fibres". DC-driven
structural changes (swelling, electro-osmotic water flow, pH fronts) are listed as untested except
through their bulk effects (heating, water loss, gas), §6. **Narrowed.**

### 25. MAJOR: "The E⊥ argument against the wire model fails at 0.6 THz."
**Response.** Correct (6.6 vs 2.7 dB at 0.6 THz for +400 S/m). Dropped; the case rests on the
conductivity required (200–700× physiological, item 44) and the missing mechanism (membrane admittance).
**Fixed.**

### 26. MAJOR: "The 9 dB target has no traceable source."
**Response.** Read off the published spectra: the 0.6 THz dip is 9–12 dB below the no-voltage
spectrum after 5 min (paper 3, Fig. 1b, 100 µm) and 15–20 dB after 4 min (paper 1, Fig. 4b, 70 µm);
power spectra in dB. The digitised points (±2 dB) are in `experiments/exp/sim38_gas_path.py` (replaced by pixel tracing, item 37).
The raw data on IEEE DataPort (doi:10.21227/h17h-tp94) need a subscription we do not have; the
paper says so and offers the analysis code. **Fixed.**

### 27. MAJOR: "The dip drifts (0.7 → 0.5 THz); fixed water-vapour lines cannot explain it."
**Response.** Correct; the vapour-line attribution is withdrawn as the explanation. Sim 38 tests a
gas path through the slice with physical tissue optics. The reported n (3.6–3.9 at a stated 100 µm)
and the 2 THz loss independently imply an optical thickness of 190–300 µm. A two-parameter fit
(beam fraction 0.11–0.22, column 225–240 µm) to the published 5-min spectrum (rms 3.1–3.5 dB vs
±2 dB reading error) puts dips at 0.55–0.64 and 2.0–2.2 THz, both reported frequencies, raises
2–3 THz as reported, and moves the first dip down as the column grows. The bias electrolyses
~40–60 mm³ of gas in 5 min, the order of the slice volume. Stated as a fit and a hypothesis; it
does not explain the E∥ or DC-direction selectivity. **Narrowed / new result.** (Withdrawn in round 4, items 36–41.)

### 28. MAJOR: "Novelty: the null result follows from k_w·b ≤ 0.6, Re ε > 0 and n_myelin < n_water; prior work omits Zangari 2018, Kumar 2016, Xiang 2021."
**Response.** Added to the prior-work table. The one-paragraph argument does not give the node's
absorbed power, its scaling law, the full-wave bound on node states including negative-ε ring
resonances (which do exist, η ≈ 100, item 17), the quantitative test of the experimenters' model,
or the test of the published waveforms (item 41); those are the contributions. The introduction now
motivates the closed-form law by node-level THz dosimetry. If the editor prefers, the paper fits a
regular article of moderate length. **Fixed.**

### 29. MINOR: "The axon mesh check refined the wrong region."
**Response.** The production check refined only the node region; a second check refined the whole
axon 4× (`scripts/probe/probe_axon_mesh.py`, 0.44 µm elements vs 0.92 µm skin depth at 10⁵ S/m,
3 THz): changes ≤ 7 × 10⁻⁴ up to 10⁵ S/m. REPORT row 37 now cites it. **Fixed.**

### 30. MINOR: "Why do the absorbed-power and dipole observables differ?"
**Response.** At low σ the absorbed power per incident intensity uses ½n′|E|²/η0 and misses the
coupling to host loss; their ratio is |n_t|²/n_t′² (1.30 at 0.1 THz, 1.04 at 2 THz). The coherent
observable −k0 Im(Δε/n) is the right one for transmission; above ~10⁴ S/m Re Δε dominates and the
first-order form fails, so the axon results are now converted with an exact slab on PTFE (`slab_rta`,
column dB_exact_slab). Text added to §2.4. **Fixed.**

### 31. MINOR: "Negative-ε node with ε″ = 0.01 has contrast in the corner-singular range."
**Response.** Agreed; those points are not mesh-converged (η 100 → 103). The claim is made for
ε″ ≥ 0.1, which covers any aqueous material. **Narrowed.**

### 32. MINOR: "Wiener 'within 3 %' holds only to 1 THz; 'wire-like fibres need 500 S/m' contradicts your own crossover."
**Response.** Wording fixed: 3 % at ≤ 1 THz, 10 % at 3 THz; "a 9 dB change from fibre conduction
would need ~500 S/m"; at that conductivity fibres are resistive, not wires. **Fixed.**

### 33. MINOR: "Tissue can absorb more than water (Berry 2003)."
**Response.** Fixed: nerve tissue absorbs within ~10 % of water; the reported loss is 1.9–2.6× water
(broadband / Kindt–Schmuttenmaer). **Fixed.**

### 34. MINOR: "The DC bias breaks inversion symmetry."
**Response.** Estimated: χ(2)_eff = 3χ(3)E_dc with water's Kerr-derived χ(3) ≈ 2 × 10⁻¹⁹ m²/V² gives
≤ 10⁻¹¹ of the linear response (Sim 37). **Fixed.**

### 35. MINOR: biology and geometry.
**Response.** Added: the tissue is post-mortem; the bias puts ~70 mV across a 10 µm fibre, not ~1 V;
the selectivity is stated in the electrode frame (E_THz ⊥ E_DC, parallel to electrodes 1–2); finite
fibre segments in a 70–100 µm slice only lower the E∥ response, so 500 S/m is a lower bound. **Fixed.**

## Fourth round (independent reviewer, 2026-11-02; node and fibre bounds recomputed and confirmed; gas-path account judged fatal as written)

### 36. FATAL: "The gas column does not reproduce a dip at 2 THz: in paper 3 the 5-min minus no-voltage spectrum has at most a 0.2–2 dB ripple there, and paper 1's 2 THz feature drifts the other way."
**Response.** Correct. Pixel tracing of paper 3 Fig. 1b (Sim 39) shows red − blue rising from
+9.2 dB at 1.5 THz to +18.3 dB at 2.5 THz with no local drop above 0.4 dB: there is no dip near
2 THz in that spectrum. The claim "both dip frequencies" is withdrawn everywhere, and Sim 38 is
withdrawn (REPORT row 38). The paper no longer offers an explanation for the 2 THz features of
paper 1. **Withdrawn.**

### 37. MAJOR: "The hand-read points are biased; the fit is poor."
**Response.** Correct: the eye-read points were 2–5 dB low above 1.2 THz. A reproducible tracer
now renders the PDF page at 300 dpi, calibrates both axes on the gridlines (misfit ≤ 0.03 dB,
≤ 0.006 ps) and traces each curve by colour (`src/thznerve/digitise.py`, tested;
`results/sim39/digitised.npz`). Its waveforms reproduce the published tissue insertion loss to
0.8 dB rms over 0.3–2 THz. On the traced points Sim 38's tissue gas column gives χ² = 32–52 on 8
degrees of freedom. **Fixed (fit withdrawn).**

### 38. MAJOR: "The reported polarisation and DC-direction controls contradict an air path: a hole in the slice would act on any polarisation and any bias direction."
**Response.** Agreed. They are stated as contradictions of the bypass reading, not as open points:
if the controls were measured as described, the bypass does not explain the selectivity, and
nothing in this paper does. **Narrowed.**

### 39. MAJOR: "The gas-volume argument is empty: the gas forms at the electrodes 7.5 mm away, a tungsten anode corrodes, and the air path needs ~1 mm³ against tens produced."
**Response.** Agreed. Sim 37 now gives 25–42 mm³ of H₂ at the cathode (up to 1.5× with O₂) and says
the volume neither supports nor excludes an air path (20 % of a ~2 mm spot in 0.1 mm is ~0.1 mm³).
The cause is now "an air path of unknown origin" (gas, the slice drying or pulling away from a
window, or the beam or sample moving). **Fixed.**

### 40. MAJOR: "A 225–240 µm column in a 100 µm cuvette is implausible; with t ≤ 100 µm the fit collapses."
**Response.** Confirmed (x → 0, rms 11 dB). The model-free test below needs no thickness: it uses the
measured tissue pulse, whose delay already contains whatever thickness and index the slice has.
The reported index (3.6–3.9 at 100 µm) and the 2 THz loss remain inconsistent with water-based
tissue at 100 µm; we report that as the authors' inconsistency, not as an optical thickness. **Fixed.**

### 41. MAJOR: "Fit the published waveforms directly: red = a·blue + b·black(t − τ)."
**Response.** Done (Sim 39, Fig. 7). a = 0.587, b = 0.215, τ = −0.02 ps, residual 28 % of the
trace's rms (a + b = 1 forced: x = 0.24, 33 %; blue also shifted: 23 %). One-path changes fit much
worse: the no-voltage pulse scaled and shifted 77–94 %, the empty-cuvette pulse scaled and
shifted 65 %. So the 5-min trace contains ~20 % of the empty-cuvette pulse at the empty-cuvette
delay: part of the beam no longer crosses tissue. In frequency, the fitted sum puts the dip at
0.61 THz (reported 0.56) with 6.8 of the reported 11.4 dB, and follows the published change within
1.5 dB rms over 0.3–1.7 THz; above 1.7 THz the traced 7-ps window no longer reproduces the published
spectra (> 1.5 dB), so the fit is not tested there. a + b ≈ 0.8: the remaining tissue path also loses ~25 % in amplitude. Stated as
"consistent with an air path through ~20 % of the beam", not as a reproduction. **New result.**

### 42. MINOR: "Fresnel factors and the PTFE etalon were ignored."
**Response.** They cancel in the waveform test, which compares measured pulses through the same
cuvette. For the withdrawn tissue model they changed the ratio by ≤ 0.3 dB (reviewer's estimate,
not used further). **Moot.**

### 43. MAJOR: "One bound everywhere: 0.1 dB in the table, 0.5 dB in the abstract, 1–2 dB for small fibres."
**Response.** One statement now, in abstract, PLAN §4.6 and the table: for E∥, internode ≥ 200 µm
and fibre packing ≤ 0.6, node conductivity changes a 100 µm slice by ≤ 0.05 dB and a passive node
permittivity (ε″ ≥ 0.1) by ≤ 0.9 dB (items 51, 59), the latter only near a negative-ε ring resonance; the smallest
CNS fibres (internodes 50–100 µm) can reach 1.5–3 dB, 3.6 dB with ±20 % material data (items 57, 59). All are ≥ 2.5× below the smallest reported change. **Fixed.**

### 44. MINOR: "Bracket the 9 dB conductivity with full-wave points; give its packing dependence; '350× physiological' needs a second reference."
**Response.** Full-wave points added at 600 and 800 S/m (`sim36_review.py axon2`). 9 dB (E∥, exact slab) needs 470–565 S/m at fibre packing 0.5 and 290–980 S/m over packing 0.7–0.3 (0.6–2 THz); at 0.1 THz 600–1300 S/m.
The second reference is water itself: the axoplasm's THz water loss already equals ωε0ε″ = 82–115 S/m,
so the change needed is 3–8× the conductivity that water's own absorption represents. **Fixed.**

### 45. MINOR: "The table still says the wire model would change E⊥; items 26, 27 and 35 assert without numbers."
**Response.** Table row fixed (E⊥ changes 2.7 dB vs E∥ 6.6 dB at +400 S/m, 0.6 THz; no
discrimination). Item 26 now cites the traced points (Sim 39). Item 27's account is withdrawn
(items 36, 41). Item 35: a finite fibre segment of length L in a slice has a smaller axial dipole
than the infinite one: its axial depolarisation factor N ≈ (a/L)²(ln(2L/a) − 1) is ≤ 7 × 10⁻³ for
an axon radius of 3.5 µm and L ≥ 70 µm, which lowers the dipole by ≤ 2 % at the conductivities
needed (|Δε| ≈ 15 at 0.6 THz). Treating fibres as infinite cannot underestimate the conductivity
needed. **Fixed.**

## Fifth round (independent reviewer, 2026-11-02; verdict: major revision, no fatal flaw; tracing, node bounds and the 9 dB conductivity recomputed and confirmed)

### 46. MAJOR: "The two-path fit is statistically rejected, so 'equals' is wrong; a = 0.59 < 1 − b needs an unexplained loss in the tissue path; a slightly thicker tissue path plus the bypass fits better."
**Response.** Agreed and adopted. The residual is 8× the tracing jitter (0.013, measured before
the pulse), so the sum is not an exact description, and the abstract now says the waveform
"contains" the component. Sim 39 now fits in the frequency domain over 0.15–1.7 THz with the
measured tissue transfer T = B/K raised to a free power: 0.69 × (tissue path 1.11× thicker) +
0.21 × empty cuvette leaves 20 % (vs 27 % for a·B + b·K) and puts the dip at 0.57 THz, −10.7 dB
(reported 0.56 THz, −11.4 dB); a + b = 0.9. So the data need both an undelayed component and a
change of the tissue path; both are stated. **Fixed.**

### 47. MAJOR: "The drift claim (0.76 → 0.50 THz) is in no result file, and a pure bypass cannot give a deep dip below ~0.7 THz."
**Response.** Correct; the claim was computed in a scratch run and not committed. Sim 39 now
prints the sweep: a pure bypass (x 0.05–0.35) dips deeply (−8 to −21 dB) only at 0.74–0.76 THz;
with the tissue path 10 % or 20 % thicker, the deep dip moves to 0.61–0.70 or 0.56–0.65 THz. The
table now says this: the reported 0.7 → 0.5 THz drift is compatible with a bypass plus a slowly
thickening tissue path, qualitatively, not fitted. **Fixed.**

### 48. MAJOR: "'Without crossing tissue' is not identified: a path through 10–20 % of the tissue fits nearly as well."
**Response.** Agreed (Sim 39: 29 % and 32 % against 27 % for air). Now stated as "a component at the
empty-cuvette delay (±0.02 ps) that crosses at most ~20 % of the tissue's optical thickness".
**Narrowed.**

### 49. MAJOR: "The thickness anomaly is used selectively: the traced phase gives n ≈ 2.8–3.1 at 100 µm, and the bypass dip frequency depends on that delay; if the slice is ~200 µm every per-slice bound doubles."
**Response.** Agreed. Sim 39 prints the traced index (3.04/3.10/3.03/2.84 at 0.3/0.6/1/1.5 THz).
PLAN §5.6 now states the anomaly once and its consequence: at ~200 µm the node bounds double
(≤ 0.1 dB for conductivity, ≤ twice the passive-ε bound) and the 9 dB conductivity halves
(~150–500 S/m), still far from the reported change. The waveform test uses the measured tissue
pulse, so its dip frequency does not assume any thickness. **Fixed.**

### 50. MINOR: "Fig. 1b alone tests the bypass up to 3 THz: biased − empty is flat at −11 to −13 dB over 1.5–3 THz."
**Response.** Added (Sim 39, REPORT row 39, PLAN §5.8): an undelayed component of 0.22–0.27 of the
empty-cuvette field dominates where tissue loss is 30–50 dB, model-free and consistent with the
waveform fit (0.21–0.22). **Fixed.**

### 51. MINOR: "The passive-ε bound is under-sampled: the ring resonance lies between grid points and only 0.6/2 THz were run."
**Response.** Dense grid run (`sim36_review.py epsdense`): ε′ = −2.0…−1.0 in 10 steps (0.05 apart
near the peak) at ε″ = 0.1 (the worst case the claim covers; the peak falls as ε″ grows), at 0.6,
1, 2 and 3 THz, with a parabola through the top three points. The peak grows with frequency and sits at ε′ = −1.45 to −1.8: E∥ 0.014/0.028/0.072/0.123 dB at 0.6/1/2/3 THz (internode 1 mm, packing 0.5), so ≤ 0.43 dB at 0.6–2 THz and 0.74 dB at 3 THz for internode 200 µm and packing 0.6; E⊥ ≤ 0.11 dB. The reviewer was right that 0.5 dB was exceeded (at 3 THz); the bound is now ≤ 0.8 dB everywhere it is quoted. **Fixed.**

### 52. MINOR: stale numbers (items 1, 6, 14b, 23, 25; PLAN §4.6 E⊥; abstract law accuracy).
**Response.** Item 1 now points to the single bound of item 43; item 6 states three constants per
polarisation; item 14b gives 470–565 S/m (290–980 over packing); item 23 is marked superseded;
item 25 gives 200–700×; PLAN §4.6 adds E⊥ 3.5 dB for a metallic node in a 21 µm fibre; the
abstract gives the law's 2 % for fibres up to 10 µm in diameter (9.8 % at 21 µm, 3 THz, stated in
§4.2). **Fixed.**

### 53. MINOR: novelty.
**Response.** Framed as in item 28: the contribution is the quantitative bounds, the closed law,
the test of the experimenters' model and the analysis of their published waveforms, not the
qualitative expectation. **Unchanged.**

## Sixth round (independent reviewer, 2026-11-02; verdict: major revision, no fatal flaw; every abstract number traced to a result file; Sim 39 reproduced)

### 54. MAJOR: "'8× tracing jitter' is too strong: on the pulse a 1-pixel timing error gives 0.02–0.09 rms; the residual also depends on the band."
**Response.** Agreed. Sim 39 now redoes the fit with each trace moved by −1, 0, +1 pixel (27 cases):
the red trace changes by up to 0.023 rms, the residual stays at 25–30 % (~5× the pixel-level
noise), and b = 0.214–0.216, τ = −0.04 to −0.01 ps. The four-parameter fit leaves 9 %, 17 % and 20 %
over 0.15–1.0, 1.3 and 1.7 THz, with b = 0.206–0.212 in every band. The text now quotes "~5×" and
the band. **Fixed.**

### 55. MAJOR: "The dip at 0.57 THz is a fit, not a prediction; depth is ill-conditioned. Lead with b ≈ 0.21 at τ = 0."
**Response.** Adopted. The abstract and §5.8 now lead with the component (b and τ, robust across
bands and tracing shifts) and call the dip a consistency check whose depth is ill-conditioned.
**Fixed.**

### 56. MAJOR: "Paper 1's waveforms are untested; the drift and the 2 THz feature come from paper 1."
**Response.** Done (Sim 42, Fig. 8). `trace_track` was added to `thznerve.digitise` (tested) for
curves crossed by labels and arrows or hidden under other curves. Paper 1 Figs. 3a (empty cuvette)
and 4a (no voltage, 30 s, 4 min) share axes and units: the no-voltage trace in both agrees to 0.017
rms. Spreads below are over ±1-pixel tracing shifts. At 30 s, b = 0.04 (0.034–0.045); the residual
is 10 % against 21 % for one path, and the traced dip is shallow at 0.78 THz, where a pure bypass
puts it. At 4 min, b = 0.17 (0.171–0.175) at τ = −0.04 ps (−0.065…−0.012; the tissue delay is
~0.6 ps); one-path fits leave 83–84 %. The frequency-domain fit with the tissue path 1.09× thicker
leaves 16 % and puts the dip at 0.63 THz, −20.5 dB; the traced dip is at 0.62 THz, −19.1 dB. So the
component appears in the second report too, grows with bias time, and needs the same ~10 % thicker
tissue path as paper 3 (1.11×). The dip moves down as both grow, the reported 0.7 → 0.5 THz
direction. The paper 1 index from the traced phase is 3.9 at the stated 70 µm, the same anomaly as
paper 3. The 2 THz feature lies above the band the 75 px/ps traces resolve, among the unpurged
water-vapour lines of paper 1 Fig. 4b, so it is stated as not tested. **New result.**

### 57. MAJOR: "The small-fibre and 200 µm-thickness bounds are stale after the dense ε grid."
**Response.** Corrected everywhere: small CNS fibres 2–4× the 0.74 dB, i.e. 1.5–3 dB, ≥ 3× below
9 dB (3.6 dB and 2.5× with ±20 % materials, item 59); at a doubled slice thickness the passive bound is 1.8 dB. **Fixed.**

### 58. MINOR: wording ("~20 % bypassing", 28 % residual, "three constants", 2 %, packing ranges).
**Response.** Prior-work row and limitations updated; the abstract no longer counts constants and
quotes 2.2 %; the node bounds are stated for packing ≤ 0.6 and the 9 dB conductivity for packing
0.3–0.7, both explicitly. **Fixed.**

### 59. MINOR: "±20 % materials were tested on node power only; the ε-peak mesh check was at 0.6/2 THz, the bound comes from 3 THz."
**Response.** Both added. Mesh at the peak (node ε′ = −1.45, ε″ = 0.1, E∥), node element size 0.25 / 0.125 / 0.0625 of the gap: 0.123 / 0.107 / 0.097 dB at 3 THz, 0.071 / 0.066 / 0.062 dB at 2 THz, within 3 % at 0.6–1 THz, E⊥ within 8 %. The value falls monotonically with refinement (the gap's 90° corners are near-singular at this contrast), so the coarse-mesh 0.123 dB used for the bound is an upper bound. Materials (`sim36_review.py epsrobust`): ε′ and ε″ of axon, myelin and water (host included) each ×0.8 and ×1.2, node ε′ scanned −2…−1 at 3 THz so the peak can move: it moves with myelin ε′ and reaches at most 0.148 dB (myelin ε′ ×1.2), against 0.123 dB nominal; ×6 for internode 200 µm and packing 0.6 gives 0.9 dB. The passive bound is now stated as 0.9 dB including ±20 % materials; for the smallest CNS fibres up to 3.6 dB, still 2.5× below 9 dB. The conductivity bound (0.05 dB, 180× margin) was not rerun. **Fixed.**

### 60. MINOR: "Two reported features lack an explicit match / no match: no further change after ~6 min (paper 1); no significant change of average power (paper 3)."
**Response.** Both added to the table. Traced pulse energy 1.03 vs 1.36 (−1.2 dB), consistent with
"no significant change"; saturation after 6 min cannot be tested without a later trace. **Fixed.**

### 61. MINOR: novelty; lead with the contributions.
**Response.** The abstract now leads with the closed law and the bounds, and closes with the
waveform analysis of both reports. **Fixed.**

## Seventh round (independent reviewer, 2026-11-03; verdict: major revision; every abstract number recomputed and confirmed, scripts reproduce, tracer follows the curves)

### 62. FATAL: "The three reports are one experiment: H1's and H3's empty-cuvette and no-voltage waveforms are the same scans. The paper counts them as independent replications."
**Response.** Confirmed and adopted (Sim 43): H1 vs H3 empty cuvette 2.2 % rms, minima 0.005 ps
apart; no-voltage H1/H2/H3 1.8–3.7 %. §5.5 now treats them as one sample in one session with biased
scans at 0.5, 2, 4 and 5 min, and notes that stated thickness, spectrometer and electrode gap
differ between reports sharing these scans. "Both reports replicate" is removed everywhere. **Fixed.**

### 63. MAJOR: "Trace H2's 2-min waveform; cite the fourth report; get the DataPort data."
**Response.** H2 traced after the other fits: b = 0.120 (0.118–0.121), between 30 s and 4 min
(linear interpolation predicts 0.10), residual 25 %, one path 63 % (Sim 43). The fourth report
(Reports NAS RA 124(2):16–23) is cited as H4, not seen. The DataPort dataset (doi:10.21227/h17h-tp94)
needs a subscription; stated in §5.6. **Fixed / not available.**

### 64. MAJOR: "The 2 THz feature of paper 1 is untested."
**Response.** Tested (Sim 43). H1's spectra traced (Figs. 3b, 4b; the arrow and label in 4b corrupt
0.7–1.3 THz by up to 11 dB, so the no-voltage curve comes from 3b); the 30-s spectrum near 2 THz is
2.4 dB above the traced floor. With |T| from the spectra, the phase of T from the waveforms and the
fitted a, b, the two terms reach antiphase at 2.24 THz where at 30 s they are nearly equal: a −40 dB
dip at 2.25 THz, which vanishes by 4 min as in the published spectrum. Measured 2.00 THz, −43.5 dB:
12 % off. Stated as a partial match with the possible causes (super-linear tissue phase above
1.4 THz, water-vapour lines). **Fixed (partial match, stated).**

### 65. MAJOR: "'Reproduce' overclaims; q depends on the band."
**Response.** "Consistent with" throughout; q = 1.09–1.18 with its band dependence given. **Fixed.**

### 66. MAJOR: "The drift: at 30 s there is no dip at 0.78 THz to drift from."
**Response.** Agreed: the 30-s ratio is a broad −2.9 to −4.1 dB trough over 0.56–0.78 THz. The drift
is now described as a dip that deepens and moves from ~0.7 THz (q = 1) to 0.62–0.65 THz as b and q
grow. **Fixed.**

### 67. MAJOR: "τ < 0 in every fit is unphysical for any path through the cuvette."
**Response.** −0.026 to −0.043 ps (2–3 pixels); explained as delay-stage drift between scans, and
"undelayed" defined as within ~0.1 ps of the empty-cuvette delay (a path through ≤ 20 % of the
tissue fits as well). **Fixed.**

### 68. MAJOR: "§2.4 says the exact slab is used throughout; the bounds are first order."
**Response.** Correct. §2.4 now says the node bounds are first order and checks them against the
exact slab at their largest value (0.1223 vs 0.1227 dB); the exact slab is used for the
conducting-axon case. **Fixed.**

### 69. MINOR: "Spinal white matter packs at 0.7–0.8; the 5-point ε′ grid in the ±20 % check may miss the peak."
**Response.** Packing 0.7–0.8 scales the bound to 1.0–1.2 dB, stated in §4.6; a peak missed by the
5-point grid would have to be 7× the dense grid's to reach 9 dB. **Fixed.**

### 70. MINOR: inconsistencies (18 mm for all reports; author lists; REPORT row 36 "≤ 0.8 dB"; row 39 "frequency and depth"; Fig. 6 uses withdrawn eye-read points; Fig. 7 shows the withdrawn gas column; KeyError in `sim36_sigma.log`).
**Response.** 15–18 mm; "Hovhannisyan and co-workers"; REPORT rows 36 and 39 corrected; the
two-path figure dropped; the gas-column curve removed from `fig_bypass`; the log error comes from a
post-processing call after the values were written (values confirmed from the JSON). **Fixed.**

### 71. MINOR: missing citations (lipid-bilayer dielectric dispersion; salt and THz water absorption); title should say "below 3 THz".
**Response.** Zhang et al. 2022 and Xu et al. 2007 cited; title changed. **Fixed.**

---

## What matches and what does not (experiment)

| Reported (Hovhannisyan, Makaryan, Martirosyan 2024) | This work |
|---|---|
| Dip near 0.6 THz (9–20 dB) | No node or fibre resonance anywhere in 0.1–3 THz. The shared session's biased waveforms contain a component at the empty-cuvette delay (4, 12, 17, 22 % at 0.5, 2, 4, 5 min); no one-path tissue change fits. With the tissue path also 1.09–1.18× thicker, the dip falls at 0.57 THz (5 min, reported 0.56) and 0.63–0.65 THz (4 min, published 0.65): consistent, q fitted (Sims 39, 42, 43) |
| 2 THz dip (30 s, paper 1) | Without fitting the spectra the same two terms give a −40 dB antiphase dip at 2.24 THz after 30 s that vanishes by 4 min, as published; measured 2.00 THz, −43.5 dB (12 % off). Unpurged water-vapour lines lie at 1.92–2.20 THz (Sim 43) |
| Dip drifts 0.7 → 0.5 THz in minutes | Water-vapour lines cannot drift. At 30 s the traced ratio is a broad −3 to −4 dB trough over 0.56–0.78 THz; by 4 min the dip is at 0.62–0.65 THz as b reaches 0.17 and q 1.09 (Sims 42, 43). Matches the direction; fitted, not predicted |
| High frequencies rise by 17–18 dB (2–2.5 THz, traced) | The bypass skips tissue loss: the waveform fit follows the rise to 1.7 THz within 1.5 dB rms (not testable above with the traced window); drying alone ≤ 3 dB (Sim 34) |
| ~9 dB change under 100 V DC | Node (E∥, internode ≥ 200 µm, packing ≤ 0.6): ≤ 0.05 dB for any conductivity, ≤ 0.9 dB for any passive ε (ε″ ≥ 0.1, ±20 % materials); channels 10⁻⁹ dB. Fibre conduction: needs 290–980 S/m axoplasm (200–700× physiological). Heating/drying: smooth, ≤ 3 dB |
| Effect only with E ∥ fibres | Not explained, and contradicts the bypass; model dichroism −0.35 dB; a conducting axon is itself E∥-selective (6.6 vs 2.7 dB at 0.6 THz), so this does not discriminate the wire model |
| Effect only with DC across the fibres | Not explained, and contradicts the bypass; heating predicts the opposite (σ∥ ≈ 7σ⊥) |
| Generated high frequencies (nonlinearity) | No ps nonlinearity at TDS fields (10⁻³ of onset); drying raises high-frequency transmission linearly |
| Sample acts as a grating of wires | No diffraction orders below 6 THz; fibres are not periodic |
| Onset ~5 min, recovery ~1 h; no further change after ~6 min (paper 1) | The undelayed component grows 4 → 12 → 17 → 22 % (0.5, 2, 4, 5 min); saturation after ~6 min not testable (no later trace). Timing compatible with heating/drying (0.64 K/s adiabatic), but see DC direction |
| No significant change of average power (paper 3) | Traced pulse energy 1.03 (biased) vs 1.36 (no voltage): −1.2 dB, consistent with a small change of a heavily attenuated signal |
| Tissue n = 3.6–3.9 | Inconsistent with water-based tissue (model and Berry 2003: 1.95–2.2) |

## Remaining (non-fatal) limitations stated in the paper
- Single fibre plus effective medium; dense-packing corrections are O(1) (item 4).
- No raw experimental data; the undelayed component rests on traced published waveforms from one
  measurement session (fits leave 16–28 %, several times tracing noise; q depends on the band) and its cause is unknown. Tests: image the cuvette during bias; record DC current and slice
  mass; measure slice thickness; repeat the E⊥ and DC-direction controls on the same slice.
- Linear passive EM response only; DC-driven structural changes enter only through bulk effects.
- Passive linear model only; band limited to 0.1–3 THz.
