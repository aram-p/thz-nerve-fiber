# Referee report against our own claim, with responses

**Claim under review.** In 0.1–3 THz a single myelinated fibre in water has no resonant absorption
at the node of Ranvier. The node absorbs as a quasi-static, sub-wavelength gap: P_node =
½ωε0ε″_w|E0|²V_n·F, where F (0.2–2) depends only on node length / myelin thickness, the g-ratio and
the permittivities (F is defined by the absorbed power and computed from a Laplace problem). No node
state, and no fibre state the experimenters proposed (open channels turning fibres into nonlinear
conducting wires forming a grating), changes tissue transmission by an amount close to the reported
~9 dB, so the 0.6/2 THz features reported under DC bias (Hovhannisyan, Makaryan, Martirosyan 2024)
do not come from the linear electromagnetic response of nodes or fibres. The published waveforms
instead contain ~21 % of the empty-cuvette pulse at its own delay, plus a change of the tissue path
(Sim 39); their cause is unknown and they conflict with the reported polarisation and bias-direction
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
caps the slice change at ≤ 0.08 dB (≤ 0.5 dB at internode 200 µm; dense grid in item 51: ≤ 0.8 dB); mesh check 4 %. Negative ε′ with
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
permittivity (ε″ ≥ 0.1) by ≤ 0.8 dB (item 51), the latter only near a negative-ε ring resonance; the smallest
CNS fibres (internodes 50–100 µm) can reach 1–2 dB. All are ≥ 4× below the smallest reported change. **Fixed.**

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

---

## What matches and what does not (experiment)

| Reported (Hovhannisyan, Makaryan, Martirosyan 2024) | This work |
|---|---|
| Dip near 0.6 THz (9–20 dB), changes near 2 THz | No node or fibre resonance anywhere in 0.1–3 THz. Paper 3's biased waveform contains ~21 % of the empty-cuvette pulse at its own delay (no one-path tissue change fits); with the tissue path ~10 % thicker the dip falls at 0.57 THz, −10.7 dB (reported 0.56, −11.4; Sim 39). Paper 3 shows no dip near 2 THz; paper 1's 2 THz features not explained |
| Dip drifts 0.7 → 0.5 THz in minutes | Water-vapour lines cannot drift. A pure bypass dips deeply only at 0.74–0.76 THz; with the tissue path 10–20 % thicker the dip moves to 0.56–0.70 THz (Sim 39, qualitative) |
| High frequencies rise by 17–18 dB (2–2.5 THz, traced) | The bypass skips tissue loss: the waveform fit follows the rise to 1.7 THz within 1.5 dB rms (not testable above with the traced window); drying alone ≤ 3 dB (Sim 34) |
| ~9 dB change under 100 V DC | Node (E∥, internode ≥ 200 µm, packing ≤ 0.6): ≤ 0.05 dB for any conductivity, ≤ 0.8 dB for any passive ε (ε″ ≥ 0.1); channels 10⁻⁹ dB. Fibre conduction: needs 290–980 S/m axoplasm (200–700× physiological). Heating/drying: smooth, ≤ 3 dB |
| Effect only with E ∥ fibres | Not explained, and contradicts the bypass; model dichroism −0.35 dB; a conducting axon is itself E∥-selective (6.6 vs 2.7 dB at 0.6 THz), so this does not discriminate the wire model |
| Effect only with DC across the fibres | Not explained, and contradicts the bypass; heating predicts the opposite (σ∥ ≈ 7σ⊥) |
| Generated high frequencies (nonlinearity) | No ps nonlinearity at TDS fields (10⁻³ of onset); drying raises high-frequency transmission linearly |
| Sample acts as a grating of wires | No diffraction orders below 6 THz; fibres are not periodic |
| Onset ~5 min, recovery ~1 h | Timing compatible with heating/drying (0.64 K/s adiabatic), but see DC direction |
| Tissue n = 3.6–3.9 | Inconsistent with water-based tissue (model and Berry 2003: 1.95–2.2) |

## Remaining (non-fatal) limitations stated in the paper
- Single fibre plus effective medium; dense-packing corrections are O(1) (item 4).
- No raw experimental data; the bypass reading rests on traced published waveforms (residual 28 %)
  and its cause is unknown. Tests: image the cuvette during bias; record DC current and slice
  mass; measure slice thickness; repeat the E⊥ and DC-direction controls on the same slice.
- Linear passive EM response only; DC-driven structural changes enter only through bulk effects.
- Passive linear model only; band limited to 0.1–3 THz.
