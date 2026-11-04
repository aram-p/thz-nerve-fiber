# Referee report against our own claim, with responses

**Claim under review.** In 0.1–3 THz a single myelinated fibre in water has no resonant absorption
at the node of Ranvier. The node absorbs as a quasi-static, sub-wavelength gap: P_node =
½ωε0ε″_w|E0|²V_n·F, where F (0.2–2) depends only on node length / myelin thickness, the g-ratio and
the permittivities (F is defined by the absorbed power and computed from a Laplace problem). No node
state, and no fibre state the experimenters proposed (open channels turning fibres into nonlinear
conducting wires forming a grating), changes tissue transmission for E∥ by the reported ~9 dB or by the
≥ 0.15 rad of extra delay the 0.6 THz feature needs (lowest over ±1-pixel tracing; ≤ 0.06 rad with every extreme combined), so the 0.6/2 THz features reported under DC bias (Hovhannisyan, Makaryan, Martirosyan 2024)
do not come from the linear electromagnetic response of nodes or fibres. The published waveforms
of the three reports, which share reference scans, instead contain a component at the
empty-cuvette delay (4, 12, 17, 22 % in four biased scans after 0.5, 2, 4, 5 min), while the
tissue term arrives 0.04–0.07 ps later (Sims 39, 42, 43), consistent with the 0.6 THz dip; that delay
is a 7–18 % thicker slice, warming or a stage drift; node states cannot give it, nor can ions or
conduction within the loss the fits allow (Sim 46); the 2 THz dip is not
reproduced; their cause is unknown and they conflict with the reported polarisation and bias-direction
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
**Response.** (The 2 THz feature, "not tested" here, was tested in Sims 43–44: items 64, 72, 83.) Done (Sim 42, Fig. 8). `trace_track` was added to `thznerve.digitise` (tested) for
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

## Eighth round (independent reviewer, 2026-11-03; verdict: major revision, no fatal flaw; abstract numbers traced to result files; simulation bounds confirmed)

### 72. MAJOR: "The 2 THz 'prediction' misses; 'within 1.6–3.4 dB above 1.85 THz' is false (5.6 dB at 1.85–2.15 THz); model and data have separate dips; the miss is ~9× the phase-fit error and dispersion makes it 18 %; a transient dip in that region is generic."
**Response.** Agreed on every point. §5.6 now reports the 2 THz dip as not reproduced, gives the
band errors (+5.6 dB at 1.85–2.15 THz, −3.3 dB at 2.3–3 THz), states that the two dips are separate
(model −30 dB at 2.00 THz against −43.5 dB; published −32 dB at 2.20 THz against the model's −40 dB
at 2.25 THz), the ~8× phase-error margin, the dispersive shift to ~2.35 THz, and why a transient dip
there is generic. "Predict" and "accounts for" are removed from abstract, introduction, §5.7 and
conclusion. **Fixed (claim withdrawn).**

### 73. MAJOR: "The 2 THz dip is ~10 dB above the noise floor."
**Response.** Confirmed on H1 Fig. 3b and stated with Jepsen and Fischer 2005 and Withayachumnankul
and Naftaly 2014. **Fixed.** (Superseded by Sim 44 and item 83: traced per column, the dip's bottom
is −71 dB, 2–5 dB above the noise at −73 to −76 dB; the "−63 dB, ~10 dB above" figures were the
reviewer's reading of the binned curve.)

### 74. MAJOR: "With every extreme the paper allows multiplied, the passive-ε bound reaches 9.5 dB at 3 THz."
**Response.** Confirmed and now stated per frequency in §4.6: ×11 for internode 50 µm, packing 0.8
and a 200 µm slice gives 1.1/2.1/5.5/9.5 dB at 0.6/1/2/3 THz (materials ×1.2 included). The reported
changes are 9–20 dB at 0.6 THz (and 3–4 dB after 30 s), ≥ 2.6× above; 14 dB at 2 THz, 2.5× above;
near 3 THz the reported change is a transmission rise of up to 17–18 dB, while every passive node
state adds loss there, and the largest rise any node state gives is 0.008 dB per slice at 1 mm
internodes (0.5 dB with all extremes; checked over all ε grids). The conductivity bound with all
extremes is 0.5 dB. **Fixed.**

### 75. MAJOR: "Shared reference scans do not prove one continuous bias run; the b series is not a time series; the reuse undermines every number from the reports."
**Response.** Adopted: "four biased scans from reports that share reference scans"; whether they come
from one bias run cannot be checked; every number from the reports, including the 9–20 dB, rests on
one sample. **Fixed.**

### 76. MAJOR: "The abstract ties the thicker path to the 2.2 THz cancellation (q = 1 at 30 s); the Fig. 8 caption says nothing is fitted though the 4-min panel uses q = 1.09."
**Response.** Abstract rewritten (the 2 THz claim is withdrawn); the caption now says a, b, τ and
q = 1.09 come from the waveform fits and nothing is fitted to the spectra. **Fixed.**

### 77. MINOR: "The ε mesh check was at two ε′ values; the peak may move under refinement."
**Response.** ε′ re-scanned on the 2× finer mesh (8 values, 0.6/1/2/3 THz; `sim36_review.py
epsdense_hscan`): peaks 0.0138/0.0267/0.0662/0.1117 dB against coarse 0.0143/0.0279/0.0716/0.1227;
the 3 THz peak moves to ε′ −1.35 but stays 9 % below the coarse peak, so the coarse values bound it.
**Fixed.**

### 78. MINOR: figures (Fig. 5 cited for the passive-ε scan; Fig. 2 caption says "benchmarks").
**Response.** Fig. 5 is cited for conductivity only; the Fig. 2 caption now names the legacy artefact.
**Fixed.**

### 79. MINOR: "H3's 0.56 THz dip sits on the 0.557 THz water-vapour line."
**Response.** Stated in §5.6 with van Exter et al. 1989; its exact frequency does not discriminate.
**Fixed.**

### 80. MINOR: "The 30-s change (3–4 dB) is within the small-fibre bound."
**Response.** The per-frequency bound at 0.6 THz is 1.1 dB with every extreme, so the 3–4 dB is
≥ 2.6× above it. **Fixed.**

### 81. MINOR: numbers ("3.4 dB" is 3.3; 0.64 K/s; 300–1000 S/m omits 0.1 THz; "closed law" has fitted constants).
**Response.** 3.3 dB; 0.6–0.9 K/s in both files; 300–1300 S/m; the abstract says the law is exact in
its limits with three constants fitted to a Laplace solver. **Fixed.**

### 82. MINOR: prior work (Wilmink and Grundt 2011; Romanenko et al. 2017; van Exter et al. 1989; Jepsen and Fischer 2005); title undersells the waveform analysis.
**Response.** All cited; title ends "re-analysis of published spinal-cord THz waveforms". **Fixed.**

## Ninth round (independent reviewer, 2026-11-03; verdict: major revision, no fatal flaw; the combined bounds, water model, component sizes and 30-s q fit recomputed and confirmed)

### 83. MAJOR: "A dip whose bottom is at the floor has an unreliable depth; that does not make it noise. Dips found in only one curve are all ≤ 7 dB except this 15 dB one; the vapour-depth comparison proves nothing."
**Response.** Agreed. §5.6 now says the 1.99 THz dip is a real feature of the 30-s scan, twice as deep
as any other single-curve dip, with its depth a lower bound set by dynamic range (bottom 2–5 dB above
the noise) and its origin unknown; the vapour-depth argument is dropped. The binned 3–14 dB deficit
over 1.85–2.10 THz is stated. **Fixed.**

### 84. MAJOR: "'Every passive node state adds loss at 3 THz' is untested: only the ε′ = −2…−1 strip was run there; the quoted 0.008 dB rise is a loss."
**Response.** The full ε grid (ε′ −30…30, ε″ 0.3/3/30; 33 states) was run at 3 THz, E∥
(`sim36_review.py round9`, `results/sim30/round9.log`): the largest rise is 0.0028 dB at
ε = −1 − 0.3j (0.18 dB with all extremes), the largest loss 0.049 dB, below the strip's 0.123 dB.
§4.6 now quotes these signed values. **Fixed.**

### 85. MAJOR: records contradict the text (items 56, 73; REPORT rows 42–43).
**Response.** Items 56 and 73 annotated; REPORT rows 42–44 updated. **Fixed.**

### 86. MINOR: "Shared reference scans" is not shown: one scan traced twice differs by 3.4 %, the reports by 1.8–3.7 %.
**Response.** Now "identical within tracing noise: one sample, possibly the same scans". **Fixed.**

### 87. MINOR: the 2.5× margin holds only for E∥ (E⊥ with all extremes: 2.7/7.4/11.5 dB).
**Response.** Confirmed from `sim36_sigma.json` (E⊥ dipole 0.043/0.116/0.180 dB × 64). §1 and §4.6
now state the E∥ premise, the configuration of every published scan, and give the E⊥ numbers.
**Fixed.**

### 88. MINOR: "The 50 µm internode and packing 0.8 scalings are extrapolated."
**Response.** Solved at internode 50 µm directly (ε′ −1.45 and water, 0.6–3 THz, both polarisations):
the slice change matches the scaled 200 µm value within 0.3 %. Packing 0.8 is stated as assumed
linear (E∥ unchanged at 0.6, §4.5). **Fixed.**

### 89. MINOR: "The 0.6 THz fine-mesh peak is at the scan edge (ε′ −1.8)."
**Response.** ε′ −2.0 run on the fine mesh: 0.0134 dB against 0.0138 at −1.8, so the peak is bracketed
and below the coarse 0.0143. **Fixed.**

### 90. MINOR: figures (Fig. 7 trace ends held flat and used in the fit; Fig. 1b cites Jepsen et al. 2007, not in the references).
**Response.** Sim 42 now fits only where every curve was traced (−362.02 to −357.20 ps): b = 0.039
(30 s) and 0.173 (4 min), residual 9 % and 27 %; traced 30-s dip −4.5 dB at 0.72 THz, 4-min dip
−17.2 dB at 0.62 THz; Sim 43 rerun on the corrected traces (band errors 5.5 and 3.6 dB). Jepsen,
Møller and Merbold 2007 cited in §2.2. **Fixed.**

### 91. MINOR: wording (H4 counted among the reports that measured; "closed-form law" has fitted constants).
**Response.** "Three reports (H1–H3); a fourth was not available"; "a quasi-static law with exact
limits" in the title, §4.2 and the conclusion. **Fixed.**


## Tenth round (independent reviewer, 2026-11-03; verdict: major revision, no fatal flaw; abstract numbers, law accuracy, mesh and PML convergence, node bounds and b, τ, residuals confirmed against the result files)

### 92. MAJOR: the abstract's "no tissue-only change fits" contradicts §5.6 (a component crossing up to ~20 % of the tissue fits as well; a slice thinned over part of the beam is a tissue change).
**Response.** Now "no uniform change of the tissue path fits". **Fixed.**

### 93. MAJOR: the proposed causes are partly impossible and never sized (a slice–window gap is in series; a slice shrinking at its edges cannot clear a central 5 mm beam on a 20 × 20 mm sample).
**Response.** Both dropped. §5.6 now sizes the opening: for a Gaussian beam 5 mm across (1/e² or
FWHM), b = 0.17–0.22 needs an on-axis opening through the full thickness 1.8–3 mm across (0.7–1.2 mm
at 30 s): a hole, a tear or a gas bubble (a 2 mm bubble in a 0.1 mm layer is 0.3 mm³, against
25–42 mm³ of electrolytic gas). An on-axis opening passes more at high frequency, as the high-band
levels show (item 96). **Fixed.**

### 94. MAJOR: §4.6 compares the bounds with dip depths, which under the paper's own model are ill-conditioned interference nulls; the tissue must change T by q − 1 = 0.09–0.18, i.e. 1.1–2.2 dB and 0.2–0.4 rad at 0.6 THz, the size of the all-extremes passive bound.
**Response.** Agreed; this is the right comparison. `tissue_db` now also returns the first-order
phase change. Over every node state with ε″ ≥ 0.1 (conductivity, ε grid, dense ε′ scans) the
largest |Δφ| at 0.6 THz is 1.2 × 10⁻³ rad per 100 µm at internode 1000 µm and packing 0.5
(1.45 × 10⁻³ with the combined materials of item 98). Under nominal parameters (200 µm, 0.6,
230 µm path) node states give ≤ 0.3 dB and ≤ 0.02 rad, ≥ 4× and ≥ 9× short. With every extreme
combined the loss bound (1.5 dB) is inside the required 1.1–2.2 dB, but the phase bound (0.11 rad)
is 1.8× below the required 0.19–0.37 rad. §4.6, the abstract and the conclusion now say this,
including that loss alone does not exclude the node at the extremes. **Fixed.**

### 95. MAJOR: the headline bounds rest on parameters the paper itself rejects (100 µm slice, packing ≤ 0.6); the conclusion drops the E∥ qualifier and the all-extremes values.
**Response.** Abstract and conclusion now carry the E∥ qualifier, and the conclusion gives the
all-extremes values (1.5 dB, 0.11 rad at 0.6 THz) with the comparison of item 94. The slice factor
is now ×2.3 (item 99). **Fixed.**

### 96. MAJOR: uncertainties understated: ±1-pixel ranges omit model error; at 4 min b is 0.173 (time domain), 0.164 (spectral fit), 0.21 (high-band level), and the 4-min model lies 1.2–2 dB below the spectrum above 1.85 THz.
**Response.** The §5.6 table now gives b from both tracings and both fit forms, plus the spectral and
high-band estimates; the stated uncertainty is the spread between methods (up to ~25 %). The text
says the high-band levels exceed the time-domain b by 15–25 % and the 4-min model lies 1.2–2.0 dB
below the published spectrum above 1.85 THz: the component is larger at high frequency than the
pulse fit gives. **Fixed.**

### 97. MAJOR: novelty is not argued in the manuscript; the prior-work table is only in PLAN.md (Zangari row still "closed-form"); microdosimetry is not cited; both exact limits follow from boundary conditions.
**Response.** §1 now has "Prior work and what is new" with the table (Table 0), cites RF/microwave
microdosimetry (Gimsa and Wachner 1999; Kotnik and Miklavčič 2000), says the limits follow from the
continuity of normal D and tangential E, and states what is new: the tested interpolation for the
node in THz water, the node-state bounds including the negative-ε ring resonance, and the waveform
re-analysis. PLAN.md row fixed. **Fixed.**

### 98. MINOR: the ±20 % materials test varied one material at a time; three changes raise the peak together.
**Response.** Run (`sim36_review.py round10`, `results/sim30/round10.log`; axon ε′ ×1.2, myelin ε′
×1.2, myelin ε″ ×0.8; ε′ −2.7…−1.0 at 0.6 and 3 THz, E∥): the peak rises ×1.45 (0.6 THz, now at
ε′ ≈ −2.0, bracketed) and ×1.43 (3 THz), above the ×1.2 used. Bounds recomputed with ×1.45: passive
0.1/0.2/0.6/1.1 dB at 200/0.6, 1.5/3.0/7.7/13 dB with every extreme. **Fixed.**

### 99. MINOR: §4.3 "ten geometries" overstates (seven node geometries, three g-ratios, two outside the fit range; E⊥ error 2.25 % at g = 0.875); §5.5 presents the high loss and the long delay as two anomalies; "none gives a narrow 9 dB feature" (§5.4) contradicts partial coverage.
**Response.** §4.3 lists the variants and states that g = 0.875 and 0.54 are extrapolations; 2.2 %
→ 2.3 % throughout. §5.5 now reads the measured loss as a thickness: an effective-medium slab on
PTFE (broadband water) matching −12.1/−14.9/−20.1/−29.2 dB at 0.6/1/1.5/2 THz is 190/180/195/230 µm
thick, agreeing with the 150–190 µm from the delay; the slice factor is ×2.3. §5.4 is restricted to
uniform changes. **Fixed.**

### 100. MINOR: the drift explanation of τ < 0 is not modelled; a delay-stage drift shifts both terms.
**Response.** Sim 43 now fits a common shift of the whole biased scan plus τ of the component. With
τ = 0 the residuals rise (10/29/34/31 %); with τ free the tissue term moves later by
−0.002/+0.043/+0.053/+0.066 ps (0/7/9/11 % of the tissue delay, growing with bias time), the
residuals fall to 9/18/20/23 %, and the component stays 0.02–0.04 ps ahead. The drift explanation is
withdrawn; the lead is stated as unexplained. The tissue-term shift reproduces the thicker path of
the frequency-domain q fit independently. **Fixed.**

### 101. MINOR: figures and §5.7 (Fig. 8 shades 1.92–2.2 THz as "H₂O lines"; Fig. 6–7 legends differ from the table; Fig. 7 shows held trace ends; Fig. 4 says "closed law"; no traced spectrum reaches 0.5 THz).
**Response.** Fig. 8 now marks only the 30-s dip (1.95–2.05 THz); Fig. 7 draws only the traced
window; Fig. 4 says "interpolation law"; the table's b ranges now include the legend values (two
tracings); §5.7 states that no traced dip lies below 0.56 THz. **Fixed.**


## Eleventh round (independent reviewer, 2026-11-03; verdict: major revision, no fatal flaw; round-10 numbers recomputed with `tissue_db`: 0.19–0.37 rad, 1.1–2.2 dB, ×12.3, ×1.45 confirmed)

### 102. MAJOR: the phase bound ignores sign; the q fit needs extra delay, and states near the ring resonance advance the phase; the positive lobe was still rising at the scan edge (ε′ −2.7); no committed script computes phase over every state.
**Response.** Agreed. The combined-materials scan was extended to ε′ = −10 (−10, −7, −5, −4, −3.5, −3).
`report_round10` now computes the signed first-order loss and phase at 0.6 THz over all 116 states
with ε″ ≥ 0.1 (conductivity scan, coarse and fine ε grids, dense ε′ scans, combined-materials scan;
×1.45 on the nominal-material scans). The largest extra delay is 9.0e-4 rad per 100 µm at internode
1 mm, packing 0.5, i.e. 0.06 rad with every extreme (≥ 3× short of 0.19 rad); the largest-loss state
advances the phase; no state gives the required loss and delay together. §4.6, abstract and
conclusion now say this. **Fixed.**

### 103. MAJOR: the opening sizing is off (1.53–1.76 mm for a 1/e² beam, 2.59–2.99 mm FWHM), assumes mode-matched detection without saying so, and the high-band levels argue against a single on-axis opening (sized at 0.6 THz it would pass 0.86–0.94 at 2 THz if the waist scales with λ).
**Response.** Sim 45 (`experiments/exp/sim45_slice_checks.py`, `results/sim45/sim45.txt`) computes
it: 1.5–3.0 mm (0.7–1.2 mm at 30 s), b at 2 THz 0.87–0.94 for one on-axis opening. §5.6 states the
mode-matching assumption and now says a single on-axis opening is disfavoured by the flat
0.20–0.27 high-band level, which openings spread across the beam fit better. **Fixed.**

### 104. MINOR: the slice factor mixes water models (bounds use Kindt–Schmuttenmaer); with that water the loss gives ~306 µm at 2 THz; the delay inversion should use the same linear-with-offset fit on the model (≈ 220 µm); the model's phase offset disagrees with the measured one.
**Response.** Sim 45: loss gives 198/198/232/314 µm (KS) and 191/181/194/232 µm (broadband) at
0.6/1/1.5/2 THz; the delay, fitted the same way on the model, gives 216–220 µm; model offset −0.5 rad
against −0.22 measured, stated as unexplained. The bounds now use ×2.2 at 0.6–1 THz and ×3.1 at
2–3 THz: 1.5/2.7/10/18 dB passive, ≤ 1.2 dB conductivity, E⊥ 3.0/11.5/18 dB. The 2 THz margin at the
extremes is lost (10 vs ≥ 14 dB); §4.6 says so and gives the nominal margin (7×) and the width
argument. **Fixed.**

### 105. MINOR: with τ free the common-shift fit is a tissue-delay model, not a drift test; "2–4× the alignment error" conflicts with the 0.015 ps agreement of the minima.
**Response.** §5.6 now calls only the τ = 0 fit the drift test (residuals rise: drift does not explain
the lead) and the τ-free fit a tissue-delay model; ±1-pixel ranges added (±0.013 ps); the lead,
0.019–0.044 ps, is 1.2–3× the 0.015 ps agreement and is not interpreted. **Fixed.**

### 106. MINOR: abstract and conclusion ("at most half" is 0.58; 300–1300 S/m is the 100 µm figure; "no uniform change" was tested only as scaled, shifted or thickened copies); §4.6 compares bounds with dip depths after calling depth the wrong comparison, and at 30 s no tissue change is needed.
**Response.** Abstract: "none gives the extra delay … (≤ 0.06 of ≥ 0.19 rad)", "130–1300 S/m", "no
scaled, delayed or thickened copy of the tissue pulse fits". §4.6 drops the 4.5 dB comparison (q = 1
at 30 s needs no tissue change) and treats 2 THz as above. **Fixed.**

## Twelfth round (independent reviewer, 2026-11-03; verdict: major revision, no fatal flaw)

### 107. MAJOR: the fibre-conduction result is still judged by dip depth (130–1300 S/m for 9 dB) after §4.6 says depth is the wrong comparison; it should face the change of T the fits need.
**Response.** Sim 46 (`experiments/exp/sim46_tissue_change.py`, `results/sim46/sim46.txt`) applies
extra axoplasm conductivity to the 220 µm effective-medium slab and compares loss and delay with
(q − 1)·|T|² and (q − 1)·arg T at 0.6/1.0/1.4 THz. 23–46 S/m gives the 0.6 THz loss but 0.15 of the
delay at 0.6 THz and 0.05–0.08 at 1–1.4 THz. Abstract, §5.2 and conclusion now say this; the 9 dB
figures stay in §5.2 as the answer to the authors' own claim. **Fixed.**

### 108. MAJOR: the thicker tissue path (q) is fitted but no bulk change of the slice is tested against it; heating and drying were only compared with dip depth.
**Response.** Sim 46 tests heating (Liebe water), water gain and loss (thickness following the water
volume), ions in all tissue water and axoplasm conduction, each set to give the 0.6 THz loss.
Heating gives 0.5–0.8 of the delay at 0.6 THz and 0.06–0.34 of the loss at 1–1.4 THz; water loss
has the wrong sign; ions and conduction give loss without delay. An 8–16 % water gain gives loss and
delay within 0.84–1.15× at all three frequencies, as close as a geometrically thicker slab
(0.88–1.32×). Pairs of mechanisms barely improve on water gain alone (9–10 % against ~10 % rms
misfit), and water gain dominates every good pair. §5.4, §5.6 (what the change could be: water uptake, e.g. electro-osmotic flow or
swelling), §5.7 and Fig. 9 added. **Fixed.**

### 109. MAJOR: "optically ~220 µm thick" is never squared with the 100 µm cuvette stated in H2/H3.
**Response.** §5.5 now states the conflict: either the gap is about twice the stated one, or the
tissue absorbs ~2× water and has n ≈ 2.8–3.9, which water-based tissue cannot have. We cannot decide;
the ×2.2/×3.1 factors are used as a conservative scaling of the bounds, not as a measured thickness.
Slice and cuvette thickness is added to the decisive measurements (§5.7, §6). **Fixed.**

### 110. MINOR: "~8× what the phase-fit error allows" for the 2 THz offset overstates: the model dip rests on extrapolating the traced phase beyond its band, and the measured level is near the two-term floor.
**Response.** Sim 43 now refits the tissue phase over 0.2–1.0 to 0.2–1.6 THz: the 30-s antiphase
moves over 2.10–2.32 THz (measured 1.99), and the published 1.95–2.10 THz level is 3.2 dB above the
deepest the two terms can reach. §5.6, §5.7 and the abstract now say "not reproduced", and §5.6 "not
excluded". **Fixed.**

### 111. MINOR: §2.4 calls ×6 "the worst case" while §4.6 adds further extremes; Fig. 5b (internode 1 mm) is out of step with the bounds; §5.6 gives the wrong reason for withdrawing the gas column.
**Response.** §2.4 says "nominal values" and points to §4.6; the Fig. 5 caption states internode
1 mm and points to Fig. 9 for the binding comparison; the gas column is now withdrawn for its poor
fit on the traced points (χ² = 32–52 on 8, item 40). **Fixed.**

## Thirteenth round (independent reviewer, 2026-11-03; verdict: major revision, no fatal flaw; node part holds)

### 112. MAJOR: Sim 46 flips the sign of the measured phase offset. Sim 43 fits arg T = −2πfτ − 0.22, so the delay phase at 0.6 THz is 2.52 rad, not 2.08; the required extra delay is 0.23–0.45 rad, not 0.19–0.37.
**Response.** Agreed and fixed in `sim46_tissue_change.py` (`DELAY_T = 2πf·0.611 ps + 0.22`), in the
threshold of `sim36_review.py report_round10` (`results/sim30/round13.log`: still no state with both)
and in Sim 45's offset comparison (measured +0.22 rad against the model's −0.56 rad: 0.8 rad, not 0.3,
stated as unexplained in §5.5). The node conclusion strengthens: 0.06 rad is 3.3× short of the
0.21 rad needed even without the offset. Abstract, §4.6, conclusion, Fig. 9 updated. **Fixed.**

### 113. MAJOR: "heating does not" rests on Liebe water used beyond its range (fitted to 0–30 °C below 1 THz), which gives 20 cm⁻¹ for the 4–21 °C absorption change at 1 THz against 46 measured.
**Response.** Confirmed (Liebe gives 19.9 cm⁻¹; Millon et al. 2022, AIP Adv. 12:115319, measure
46 ± 1.5). Sim 46 now heats with Liebe's change of ε added to the Kindt–Schmuttenmaer base, as given
and scaled ×2.31 to the measured value. Heating fails with both: 0.6 THz loss at 9–29 K, but
0.08–0.60 of the loss at 1.4 THz and rms misfit 37–56 % over loss and delay at three frequencies
(water gain 10 %, thicker slab 7 %). §5.4 states both models and their range. The reviewer's
10 °C-base case is not used: the slices are measured at room temperature. **Fixed.**

### 114. MAJOR: the water gain is a thicker slab under another name: its thickness part gives 0.77 of the loss, its permittivity part 0.22; isotropic swelling gives 0.47.
**Response.** Agreed; Sim 46 now prints the split (thickness only 0.67–0.86×, permittivity only
0.07–0.22×, isotropic 0.32–0.47×). The claim is rewritten everywhere: what the fits need is more
tissue of the same kind in the beam (18–37 µm more of a ~220 µm path); no change of the tissue's
permittivity (heating, ions, conduction, water content at fixed thickness) gives it. "Only a water
gain" is withdrawn. **Fixed.**

### 115. MAJOR: the water gain is not shown to be plausible; §5.4's Joule heat would evaporate water, no water source is named, and the slice fills a 100 µm cuvette.
**Response.** Agreed. §5.6 no longer proposes a cause: it lists what would give a longer path
(swelling only in thickness, tissue moving into the beam) and says the closed cuvette and the Joule
heat argue against both, so no cause is shown to be plausible. Thickness during bias is added to the
decisive measurements (§5.7). **Fixed by withdrawal.**

### 116. MINOR: Sim 46's reference slab used Liebe water at 20 °C while the bounds and Sim 45 use Kindt–Schmuttenmaer; its |T|² and delay differ from the measured by 5–15 %.
**Response.** The base is now Kindt–Schmuttenmaer: |T|² 13.3/16.5/18.7 dB against 12.1/14.9/19.5
measured (within 1.6 dB, stated in §5.4). Sim 46 assumes the first reading of §5.5 (a ~220 µm
water-based slice); under the second (a 100 µm non-water slice) no material model exists to test, and
§5.5 says the paper cannot choose. **Partly fixed; the limit is stated.**

### 117. MINOR: "~10 %" (abstract, conclusion) vs 9–18 % (§5.6, §5.7); the 2.5–5 mg uses a 0.1 mm slice while the calculation uses 220 µm; 46 S/m is twice saturated NaCl.
**Response.** "9–18 %" throughout; the mass is now 18–37 µm more tissue, ~6–11 mg over 20 × 20 mm;
"1–2× saturated NaCl". **Fixed.**

## Fourteenth round (independent reviewer, 2026-11-03; verdict: major revision, no fatal flaw; node part holds)

### 118. MAJOR: q is a timing shift (q − 1 ≈ s/0.611 ps: 0.087 vs 1.091, 0.108 vs 1.109); the τ-free fit is the same model as a delay-stage drift, so "more tissue in the beam" is not established.
**Response.** Agreed. §5.6 now gives both readings of the fit: a later tissue term, or a stage drift
between the shared reference and the biased scans, in which case the component leads the empty
cuvette by 0.07–0.08 ps (5× the alignment error), which would need a path faster than the empty air
gap by 21–25 µm of optical path. We know no such path, but the traces cannot separate the readings.
Abstract, §5.4, §5.7, §6 (new limit) and conclusion now say "a 7–11 % thicker slice or a delay drift";
a second empty-cuvette scan after each biased one is added to the decisive measurements. **Fixed.**

### 119. MAJOR: the extra loss is not measured, since the free amplitude A absorbs it; with the loss offset free heating fits at 17–25 %, and the delay shape discriminating it is within the slab model's own error.
**Response.** Agreed. §4.6 now states that only the delay is fixed and bounds the extra loss by A ≤ 1:
≤ 3.3 dB at 0.6 THz (H1; ≤ 1.7 dB if the component takes its area from the tissue path). Sim 46 (b)
now sets each change to the delay and checks the loss: ions and conduction need 7.5 dB (excluded),
Liebe heating 76 K and 7.0 dB (excluded); Liebe ×2.31 heating (9 K, 1.1 dB) gives 1.4–1.5× the delay at
1–1.4 THz, which we report without excluding it; a thicker slab, a water gain into thickness and a
time shift fit. "The fits need more tissue, not different tissue" is withdrawn. **Fixed.**

### 120. MAJOR: the abstract and conclusion state a result that holds only under §5.5's first reading and needs the window gap to grow.
**Response.** The abstract and conclusion no longer name a tissue mechanism; they say what the data fix
(a later tissue term) and what cannot give it (node states, ions, conduction, using the loss bound,
which does not depend on the reading of §5.5 for the node). §5.6 keeps the cuvette argument. **Fixed.**

### 121. MINOR: node margin overstated if sign or reference were wrong (largest |Δ ln T| 0.175); 116 states are 80 distinct.
**Response.** `report_round10` deduplicates (80 states) and prints the largest change of ln T in any
direction: 0.175, below 0.21; §4.6 states both. **Fixed.**

### 122. MINOR: the abstract quotes 1.1 dB at 100 µm but not the 2 THz bound with every extreme (10 dB vs ≥ 14 dB).
**Response.** Abstract now gives "(10 dB at 2 THz with every extreme)". **Fixed.**

### 123. MINOR: mass 6–11 mg is the water-gain figure, not 18–37 µm of tissue.
**Response.** Now 15–24 µm (7–11 % of ~220 µm), ~6–10 mg over 20 × 20 mm. **Fixed.**

### 124. MINOR: Fig. 9 omits the thicker-path curve.
**Response.** Added, with the loss ceiling (A ≤ 1) and the time-shift point. **Fixed.**

### 125. MINOR: the 3 THz combined-materials scan scatters by ~0.02 dB off peak (12 % of 0.176 dB).
**Response.** §4.6 now gives ×1.43 ± 0.16 at 3 THz; the bounds use ×1.45; at 3 THz the measured change
is a rise, so no conclusion depends on it. **Fixed.**

## Fifteenth round (independent reviewer, whole paper, 2026-11-03; verdict: major revision, no fatal flaw; law fit, law errors, mesh, ×70 scaling and bounds recomputed and confirmed; no prior full-wave node study below 3 THz found)

### 126. MAJOR: heating with the measured temperature coefficient (9 K, 1.1 dB) is within the allowed loss and is set aside only by a frequency shape that comes from the T^q model itself.
**Response.** Agreed. Heating is now one of the readings: abstract, §5.4, §5.6, §5.7, §6 and the
conclusion say "a thicker slice, warming or a stage drift". §5.4 states that the T^q frequency
shape cannot exclude it, and that at the tracing floor 5 K (scaled) or 17 K (Liebe) suffice with
0.6–1.0 dB (Sim 46). **Fixed.**

### 127. MAJOR: the 0.21 rad floor is not robust to the ±1-pixel spread (4 min: 0.039–0.066 ps, 0.15–0.25 rad); the any-direction bound (0.175) does not clear it.
**Response.** Agreed. The floor is now the lowest over ±1 pixel at 4–5 min, 0.039 ps = 0.15 rad
(`thznerve.physics.tissue_change.min_required_delay`, tested); node states give ≤ 0.06 (2.3×). The
any-direction claim is withdrawn; §4.6 instead shows the sign behaves physically (metallic node
delays +0.042 rad, ε′ = 1 advances −0.039 rad) and says the exclusion rests on it. Ions and
conduction need 5.2 dB for the floor, still above 3.3 dB. **Fixed.**

### 128. MINOR: tracing resolution overstated; source bitmaps are 213–244 ppi, so ±1 rendered pixel is less than one native pixel.
**Response.** §5.6 now states the native resolution (0.016 ps H1, 0.012 ps H3) and that the
±1-pixel ranges use rendered pixels. **Fixed.**

### 129. MINOR: "7–11 %" and "9–18 %" for the same quantity.
**Response.** §5.6 defines both (q, frequency domain, band-dependent: 9–18 %; time-domain shift:
7–11 %); abstract and conclusion use 7–18 %. **Fixed.**

### 130. MINOR: "0.04–0.07 ps later" holds from 2 min only; "τ negative in every fit" fails at 30 s over ±1 pixel; 30-s b range should be 0.034–0.044.
**Response.** "from 2 min"; "every central fit … at 30 s the ±1-pixel range reaches +0.001 ps"; table
corrected. **Fixed.**

### 131. MINOR: the abstract says node states fail on loss; they fail on delay.
**Response.** "node states cannot, nor can ions or axoplasm conduction beyond the loss the fits
allow". **Fixed.**

### 132. MINOR: last-order ≤ 5 × 10⁻⁴ and balance ≤ 2.4 × 10⁻⁵ are overstated (9.0 × 10⁻³ thickest fibre E⊥, 1.5 × 10⁻³ thick myelin; balance 4.7–5.2 × 10⁻⁵).
**Response.** Rechecked over 433 production CSVs: last order ≤ 5 × 10⁻⁴ in 350, ≤ 1.5 × 10⁻³ thick
myelin, ≤ 9 × 10⁻³ thickest fibre; balance ≤ 5 × 10⁻⁵. §2.3 and Table 1 corrected. **Fixed.**

### 133. MINOR: k_w·Ln < 0.15 holds only for Ln = 1 µm (0.60 at 5 µm, 3 THz).
**Response.** §4.2 and §4.4 now give k_w·Ln ≤ 0.6 (0.44 at the longest measured node) and rest the
no-resonance point on the law holding at Ln = 5 µm within 2 %. **Fixed.**

### 134. MINOR: committed `fit_nodelaw.log` is stale (two-constant fit, 9.7 %).
**Response.** Rerun and committed: par c0, p, q = 1.066, 0.934, 0.742, max 2.1 %; perp 0.999, 0.826,
−0.784, max 2.0 %. **Fixed.**

### 135. MINOR: Fig. 9 legend says 116 states, hides the time-shift marker, has no node states at 1 THz, stops the ion/conduction curves at 46 S/m, and its grey band fixes the loss.
**Response.** Redrawn: legend below the axes with the state count read from the file (80), node
states labelled "0.6 THz only", curves to 70 / 180 S/m and Liebe heating to 80 K, grey region = the
delay the traces allow up to q = 1.18 and the loss up to amplitude 1. **Fixed.**

## Sixteenth round (independent reviewer, whole paper, 2026-11-03; verdict: major revision, no fatal flaw; round-15 fixes recomputed and confirmed; would move to minor if the argument is rewritten and the 2 THz node test replaced)

### 136. MAJOR: the paper is written for referees, not readers; the abstract uses undefined terms ("empty-cuvette-delay component", "tissue term", "every extreme"), and §4.6 and §5.6 read like response letters.
**Response.** Abstract, §4.6 and §5.6 rewritten. Two terms are defined once and used throughout: the
*tissue pulse* (the unbiased pulse) and the *second pulse* (beam at the empty-cuvette delay); the fit
is written as one equation, biased(t) = a·B(t − s) + b·K(t − τ), with q − 1 ≈ s/0.611 ps stated
there. "Every extreme" is defined where the factors are multiplied (§4.6, Bound). §4.6 now runs: what
the fits need (delay, loss bound, resolution floor), node states give the wrong sign, 2 THz. §5.6 runs:
data, model, second pulse, later tissue pulse (with the drift reading), 0.6 THz dip, 2 THz dip,
what the second pulse could be, what could make the tissue pulse later. **Fixed.**

### 137. MAJOR: the 2 THz node argument compares a loss bound with the dip depth and invokes a "material line", contradicting the paper's own interference reading.
**Response.** Agreed; replaced by the phase test (Sim 43 (e)). Each of 65 distinct node states, with
every extreme, was applied to T in the 30-s two-pulse model as a change constant over 1.8–2.4 THz, for
each of six phase bands. Three state–band pairs move the minimum within 0.03 THz of 1.993 THz, at most
−34.6 dB deep (published −43.5 dB); none does both. The reviewer's state (ε = −2.0 − 0.1j) reaches
1.969 THz at −32.8 dB: its 5.4 dB loss spoils the cancellation, as the reviewer expected, and it
changes the 0.6 THz delay by −0.005 rad. The "material line" sentence is gone. **Fixed.**

### 138. MINOR: the loss margin against ions and conduction is thin and only H1 is reported (H3 allows 4.6 dB against 5.2 dB needed at the floor); lead with plausibility.
**Response.** §4.6 gives both bounds (3.3 dB H1, 4.6 dB H3). §5.2 and §5.4 now lead with the
conductivity: 111 S/m axoplasm and 42 S/m tissue water at the floor are 4× and 1.7× saturated NaCl;
the loss comes second. **Fixed.**

### 139. MINOR: the delay floor is a ±1-rendered-pixel sensitivity, not a confidence bound.
**Response.** Stated as a resolution limit that ignores the model error (residuals 18–23 %); with
native pixels it is about 0.13 rad, 0.11 rad including the 2-min scan; nodes are 1.7–2.1× short of
those. **Fixed.**

### 140. MINOR: abstract and conclusion list stage drift as an equal reading and give "7–18 %" although 4.8 % suffices at the floor.
**Response.** Drift is now "disfavoured" in the abstract's source sections, §5.6, §5.7, §6 and the
conclusion (it needs a path 21–25 µm faster than air); the abstract names thicker slice and warming
only. Range "5–18 %". **Fixed.**

### 141. MINOR: the node sign check is not traceable (`sim36_signed06.json` deleted; `round14.log` prints the retired 0.21 rad threshold).
**Response.** `report_round10` rerun; `results/sim30/round16.log` (committed) prints the floor
(0.15 rad), the sign check (+0.0423 / −0.0393 rad) and the 2 THz state list written to
`results/sim30/sim36_signed2thz.json`. **Fixed.**

### 142. MINOR: nits (Sims 27–44; "2–4 native pixels"; conductivity bound against 9 dB depth; 200–230 µm vs broadband 181–194 µm).
**Response.** Sims 27–46; "−0.026 to −0.042 ps, up to 2.6 native pixels"; §6 notes that the delay
test covers the conductivity scan; §5.5 gives the broadband-water 181–194 µm. **Fixed.**

## Seventeenth round (independent reviewer, whole paper, 2026-11-03; verdict: major revision, no fatal flaw; law, slice bounds, node delays, loss bounds, ion and conduction requirements recomputed and confirmed, the latter also by an independent transfer-matrix calculation; no new prior work beyond Guo et al. 2022)

### 143. MAJOR: abstract and conclusion say no node state "moves the model's dip there", but 3 state–band pairs do move it within 0.03 THz of 1.993 THz; they fail on depth, a criterion §5.6 calls not robust.
**Response.** Abstract and conclusion now say no node state makes the model's dip at 2 THz as deep
as measured. §5.7 ends "Node states fail at 2 THz on depth, not on position" and gives the reason
depth is a fair test there: the measured 2 THz depth is a lower bound set by dynamic range, unlike
the 0.6 THz depth. Sim 43 (e) now also scans b over 0.034–0.044: the deepest placed minimum is
−37.5 dB (b = 0.034), against −43.5 dB, reproducing the reviewer's number. **Fixed.**

### 144. MAJOR: stage drift is judged by a double standard; in the preferred reading the second pulse also leads the empty-cuvette pulse (0.019–0.044 ps), so both readings need a path faster than air.
**Response.** Agreed. Drift is restored as an equal reading in the abstract, §5.6, §5.7, §5.8, §6 and
the conclusion. §5.6 states that both readings have the second pulse lead (6–13 µm of optical path
if the tissue pulse is later, 21–25 µm if the stage drifted), that we cannot explain either lead,
and that the traces cannot separate them; a drift needs no change of the tissue. **Fixed.**

### 145. MAJOR (readability): §4.6 uses T^q, s, the second pulse, the 30-s model and the 220 µm thickness before §5.5–5.6 define them; title 33 words; abstract cryptic.
**Response.** All waveform-based tests moved to a new §5.7 "What could make the tissue pulse later",
after §5.6 defines the model: what the fits need, node states (sign), node states at 2 THz, ions and
conduction, warming, more tissue in the beam, stage drift, the reported selectivity. §4.6 keeps only
the loss bounds and points forward; §5.2 and §5.4 point to §5.7. Title shortened to 19 words. Abstract
rewritten without "4–22 % of it", "Together they" and the bare "(≤ 0.06 of ≥ 0.15 rad)". **Fixed.**

### 146. MAJOR (criterion 5): every cause offered acts for both polarisations and heating predicts the opposite DC direction; the abstract does not say this contradicts two of the four reported observations.
**Response.** Stated in the abstract ("Each remaining cause would act for both polarisations, against
the reported selectivity (described without data)"), the introduction, §5.7 (last paragraph, which
also calls the two-pulse model a decomposition of the waveforms, not an explanation of the effect),
§5.8 and the conclusion. **Fixed.**

### 147. MINOR: the abstract slice bound drops "internode ≥ 200 µm, packing ≤ 0.6".
**Response.** Added. **Fixed.**

### 148. MINOR: the abstract quotes the 0.15 rad floor; native pixels give 0.13 rad, 0.11 rad with the 2-min scan; the margin is 1.7×.
**Response.** Abstract and conclusion use ≥ 0.11 rad; §5.7 gives 1.7× against the lowest floor and
2.3× against 0.15 rad. **Fixed.**

### 149. MINOR: "the any-direction change of 0.175 would not clear the floor" reads backwards.
**Response.** Reworded: 0.175 exceeds the floor, so a test of size alone would not exclude node
states; the exclusion rests on the sign. **Fixed.**

### 150. MINOR: H1 describes the sample "on 1.2 mm thick PTFE substrate"; if open-faced, the thickness conflict and the case against thickening weaken.
**Response.** Checked: H1 p. 2 says samples "are placed inside a cuvette made from PTFE plates with a
thickness of 1.2 mm" and refers to the "empty cuvette"; H2 gives a 100 µm PTFE cuvette. §5.7 now
says this and that H1's Fig. 1b shows the slice on one plate, so whether the gap is fixed is not
stated. **Fixed (clarified).**

### 151. MINOR: the 3.2 dB margin is the 1.95–2.10 THz band mean; at the 2.00 THz bin the published −43.5 dB lies 2.8 dB below the floor at nominal b and is reached only for b ≤ 0.036.
**Response.** Confirmed in Sim 43 (floor −48.0/−40.7/−37.1 dB at b = 0.034/0.039/0.044; reached for
b ≤ 0.037). §5.6 now gives both the band mean and the bin. **Fixed.**

### 152. MINOR: the passive-ε peak is not mesh-converged (0.123 → 0.107 → 0.097 dB, 21 %); the < 3 % bar is met only for node power.
**Response.** Added to §6 as a limit; the coarse value is used as an upper bound because it falls
monotonically with refinement. **Fixed (stated).**

## Eighteenth round (independent reviewer, whole paper, 2026-11-03; verdict: major revision, no fatal flaw; law, slice bounds, delay requirement, loss allowances, waveform fits and bias side effects recomputed and confirmed)

### 153. MAJOR: the node-state exclusion on delay is a pixel limit, not a statistical one; fixing s at the node value raises the misfit only ×1.27–1.38 and the residual is correlated over ~0.2 ps (13–14 independent points).
**Response.** Agreed that the ±1-pixel floor was the wrong yardstick: s is measured between curves on
the same axes, so axis registration cancels and the residual sets the error. New Sim 47 gives a
proper interval two ways. (1) Profile of the misfit with N_eff from the residual autocorrelation
(correlation length 0.31–0.34 ps, N_eff 16–17): 95 % intervals 0.023–0.065, 0.025–0.081 and
0.020–0.108 ps at 2, 4 and 5 min. (2) Moving-block bootstrap of the residual (400 resamples): s =
0.043 ± 0.008, 0.053 ± 0.011, 0.066 ± 0.014 ps, one-sided 95 % lower bounds 0.028/0.035/0.040 ps
(0.10/0.13/0.15 rad at 0.6 THz). The node limit (0.016 ps) is rejected at 2.4/2.5/2.1σ per scan,
which matches the reviewer's 1.5–2σ order, and at 3.3σ over the three independent scans
(p = 9 × 10⁻⁴); s = 0 at 5.3σ. The paper now states these levels instead of "cannot", and §6 notes
that the intervals treat the misfit as correlated noise. `tissue_change.SHIFT_PS` now holds the
bootstrap bounds (tested against `sim47.json`). **Fixed.**

### 154. MAJOR: the abstract says the second pulse "bypasses the tissue", but it arrives 0.019–0.044 ps before the empty-cuvette pulse, which such a pulse cannot do.
**Response.** Abstract: "a second pulse near the empty-cuvette delay … arrives 0.02–0.04 ps early,
unexplained. The waveforms are decomposed, not explained." §5.6 defines the second pulse by its delay
only, and "What the second pulse could be" now calls the bypass reading a hypothesis for this
reason. **Fixed.**

### 155. MAJOR: novelty and readability: both limits of the law are textbook; the cited resonance papers work at 10–40 THz, so "no resonance below 3 THz" contradicts none of them; the abstract carries ~20 numbers.
**Response.** The introduction now says both: the limits are textbook electrostatics and the new part
is a closed-form interpolation checked against full wave over realistic geometry, and a null result
below 3 THz does not contradict the 10–40 THz models; it bears on the node reading of H1–H3 and on
dosimetry. The abstract is rebuilt as three numbered findings with fewer numbers. **Fixed.**

### 156. MINOR: "both g-ratios outside the fit range" are within 2.3 %; only g = 0.875 is (g = 0.54 has b = 6.5 µm, 4.35 %).
**Response.** Corrected in §4.2. **Fixed.**

### 157. MINOR: "≥ 111 S/m to give that delay" belongs to the 0.147 rad floor; at the 0.11 rad bound it is ~71–83 S/m with ~3.3 dB, inside H3's allowance.
**Response.** Sim 46 now tests each scan at its own 95 % lower bound against its own loss
allowance: 2 min 0.104 rad, 31 S/m (tissue water) or 81 S/m (axoplasm), 3.8 dB (no spectral fit, so
no allowance); 4 min 0.130 rad, 38 or 99.5 S/m, 4.65 dB against 3.28 dB (H1); 5 min 0.150 rad, 43 or
113 S/m, 5.27 dB against 4.55 dB (H3). The loss exclusion holds where an allowance exists; the
conductivity (1.2–1.7× and 3.2–4.5× saturated NaCl) is stated first. Conclusion uses ≥ 81 S/m.
**Fixed.**

### 158. MINOR: "5–18 % thicker slice": 18 % is the optical-path factor q − 1; Sim 46 gives 16.6 % thickness for q = 1.18.
**Response.** "3–17 %" (3.4 % at the 2-min lower bound, 16.6 % at q = 1.18). **Fixed.**

### 159. MINOR: the abstract slice bound drops "ε″ ≥ 0.1".
**Response.** Added in the abstract and the conclusion. **Fixed.**

### 160. MINOR: H1's arrow covers the 30-s curve at 1.96–1.99 THz, hiding the dip's low-frequency flank on which the 0.03 THz position test depends.
**Response.** Stated in §5.6 (centre uncertain by up to ~0.03 THz). The node test at 2 THz fails on
depth, which the arrow does not affect. **Fixed (stated).**

## Nineteenth round (independent reviewer, whole paper, 2026-11-03; verdict: major revision, no fatal flaw; recommended change: state the waveform result as a phase lag, test node states with their whole-band change inside the fit with a degrees-of-freedom correction, and take the σ values out of the abstract and conclusion)

### 161. MAJOR: the σ values of Sim 47 are model-dependent, have no degrees-of-freedom correction, and the joint test treats scans that share the empty-cuvette pulse as independent; lower bounds came from the bootstrap and significances from the profile.
**Response.** Agreed. Sim 47 now uses an F test with N_eff − 4 degrees of freedom: a free s beats
s = 0 at p = 0.007/0.008/0.03 (2.7/2.6/2.2σ) per scan, and the profile gives one-sided 95 % lower
bounds of 0.021/0.025/0.020 ps, which now go into `tissue_change.SHIFT_PS` (tested against
`sim47.json`). The joint test is dropped: the scans share K and one session. The abstract and the
conclusion carry no σ values. **Fixed.**

### 162. MAJOR: restate the result as the tissue-pulse phase lag at 0.6 THz, which is robust to the model form; s itself changes sign when a phase term is added.
**Response.** New Sim 48 (A) fits two forms (pure shift; shift plus a Hilbert-transform phase term)
over three windows (whole trace; t0 − 0.6 to t0 + 0.8 ps; without t0 ± 0.4 ps), with block-bootstrap
spreads. The lag at 0.6 THz is 0.16–0.31, 0.20–0.40 and 0.25–0.56 rad at 2, 4 and 5 min (−0.01 to
+0.04 rad at 30 s); every one-sided 95 % lower bound is ≥ 0.11 rad. The range is wider than the
reviewer's 0.16–0.28 rad because leaving out the main pulse raises it; the paper quotes the full
range. §5.6, the abstract and the conclusion now state the lag. **Fixed.**

### 163. MAJOR (from item 6): node states delay up to 0.40 rad at 2 THz, so a test at 0.6 THz alone understates them; apply each state's ΔlnT(f) inside the fit; the maximum at 0.6 THz is 0.063 rad, not 0.06.
**Response.** This reverses the node-state exclusion. `sim36_review.py report_band` exports every
node state's loss and delay at every frequency it was solved at; `waveform_fit.node_transfer`
interpolates them (linear in f, zero at f = 0, the every-extreme factor applied per frequency) and
Sim 48 (B) applies them to the unbiased pulse inside the fit. With every extreme, ε = −2.7 − 0.1j
(combined materials) recovers 96/78/72 % of what a free shift buys and is not rejected (p =
0.54/0.17/0.21): its delay rises from 0.06 rad at 0.6 THz to 0.57 rad at 3 THz. States an aqueous gap
can hold (extra σ ≤ saturated NaCl, ε′ ≥ 1, or ε″ ≥ 3|ε′|) recover ≤ 45 % (p = 0.03/0.03/0.06); at
internode 200 µm, packing 0.6 and the 220 µm path no state recovers more than 22 % (p ≤ 0.04). §5.7
now says node states reproduce the lag only with every extreme combined and a permittivity an aqueous
gap cannot hold; "wrong sign" and "3.3σ" are gone. The 0.063 rad value is used in Sim 47. **Fixed;
the claim is weaker.**

### 164. MAJOR (from item 5): at the corrected lower bound, tissue-water ions need 21–30 S/m, at or below brine, inside the H3 allowance; reword against physiological saline.
**Response.** Confirmed, and the in-fit test goes further. New Sim 49 applies every slab candidate
of Sim 46 over the whole band inside the fit, with the tissue-pulse amplitude a ≤ 1. Ions in the
tissue water at 23/30/41 S/m or axoplasm conduction at 60/78/106 S/m reproduce the lag (92–133 % of
what the shift buys); with a free the fit prefers a = 1.6–1.9, a gain. The paper now says the
waveforms allow conduction, and that 23–41 S/m is 0.9–1.6× saturated NaCl and 15–27× physiological
saline, which is not plausible after 2–5 min. A 6.1/7.5/9.3 % thicker path or 7–11 % water gain
into thickness fits better than a pure shift (a = 0.69–0.83); warming recovers 39–66 % (p =
0.05–0.10). Sim 46 was rerun with the new bounds (Fig. 9). **Fixed; the ion exclusion is now a
plausibility argument.**

### 165. MAJOR (from item 7): the 2 THz depth argument is used both ways, and the model is 5.5 dB off in that band; soften.
**Response.** §5.7 and the match table now say the dip is not reproduced with or without node
states, and that this does not exclude them. **Fixed.**

### 166. MINOR: the abstract sentence on node states is ungrammatical; "REFEREE item 40" in the body; Sim numbers belong in a supplement.
**Response.** The abstract is rewritten (249 words, set with `thznerve.manuscript.set_abstract`);
the REFEREE reference is removed. Sim numbers stay as pointers to the repository for now and will
move to the supplement at submission. **Fixed in part.**

### 167. MINOR: the novelty is modest.
**Response.** Agreed and stated in the introduction since round 17; the in-fit test of every
candidate cause is added to the list of what is new. **Stated.**

## Twentieth round (independent reviewer, whole paper, 2026-11-03; verdict: major revision, no fatal flaw; recommended change: one rule for every candidate, passing ranges instead of best fits, and a reference that does not understate the noise)

### 168. MAJOR: the s-free reference leaves a structured misfit; a smooth extra loss of the tissue pulse halves it, so the F tests of Sims 48–49 used an understated noise level and were lenient.
**Response.** Agreed. New Sim 50 adds a smooth loss L_κ (κ dB/THz above 0.2 THz, 1 at 0.2 THz) to the
reference, zero-phase, with a causal (minimum-phase) version as a check. κ grows with bias time
(2.1/5.6/6.5/10.0 dB/THz, recomputed with the loss held at 1 below 0.2 THz) and the misfit falls from 4.3/6.5/10.7 to 1.6/3.4/5.5 at 2/4/5 min. Every
candidate is now fitted with s = 0 and κ free. With this reference, conduction and node states at every
extreme, which passed in Sims 48–49, fail (items 170, 171). The causal reference finds a shift at 30 s
(p = 0.004, against 0.53 zero-phase) and depends on how the loss is continued above the band, so it is
a check only; both are reported (Sec. 5.7, Sec. 6). **Fixed.**

### 169. MAJOR: two methods give one lower bound (bootstrap for bounds, profile for significance); best-fit sizes are quoted where the claim needs the range that passes.
**Response.** One rule throughout: a candidate, or a fixed shift, passes if its fit is not
significantly worse than the reference (F test, N_eff − 5 dof, p ≥ 0.05). The 95 % intervals of s
(0.020–0.060, 0.024–0.080, 0.028–0.112 ps) and the smallest lags (0.08/0.09/0.11 rad) come from the
same rule. Sec. 5.7 quotes the passing ranges (thicker path 3–8/3–11/3–15 %, warming 6–24 K or 3–10 K
scaled) and the split of each candidate's recovery into its loss and its delay. **Fixed.**

### 170. MAJOR: the node-state claim "only with every extreme" is false: in Sim 48 the best state passes from 0.6–0.75 of every extreme.
**Response.** Against the Sim 50 reference no node state passes at every extreme (best ε = −3.0 −
0.1j, p = 0.050/0.024/0.026; aqueous states p ≤ 0.014), so no scaled-down state passes either. Against
the causal reference, ε′ ≈ −1.45 passes at 2 min from 0.35 of every extreme and at 4 min from 0.70 (rerun), and
fails at 5 min; the text says so and names the ring resonance, which an aqueous gap cannot hold.
**Fixed.**

### 171. MAJOR: the node delay was scaled by the loss-based factor 3.1 at 2–3 THz; the delay-based path factor is 2.2.
**Response.** `node_transfer` now scales the delay by 2.2 at every frequency by default
(`delay_factor="path"`, tested); Sims 48–49 keep the old factor and say so. In Sim 50 the best state
gives the same p with either factor (0.050 vs 0.051). **Fixed.**

### 172. MAJOR: Sim 49 passes tissue-water ions at 23–41 S/m with a ≤ 1; the text excludes conduction on conductivity alone.
**Response.** Against the Sim 50 reference no conductivity passes (best 12–24 S/m tissue water,
30–65 S/m axoplasm, p = 0.008–0.030); the added loss works against the fit (its loss alone recovers
−19 to −24 % of what the shift buys). It fails against the causal reference too (p ≤ 0.010). Sec. 5.7
states that Sim 49 passed it only without the loss term. **Fixed.**

### 173. MINOR: Fig. 9 and Sim 46 label single-frequency candidates "NOT allowed", contradicting the in-fit test.
**Response.** The grey band of Fig. 9 is now the lag from Sim 50; Sim 46 says "within/above the 0.6 THz
allowance" and that the binding test is the in-fit one; the caption points to Sec. 5.7 (Sims 48–50).
**Fixed.**

### 174. MINOR: the mass estimate (6–10 mg) uses the old 15–24 µm.
**Response.** Now 3–15 % of the ~220 µm optical path, 7–33 µm, ~1–13 mg over 20 × 20 mm depending on
a 100 or 220 µm slice. **Fixed.**

### 175. MINOR: a duplicated sentence in Sec. 5.6; "Sims 27–47"; stale rows in the match table; the abstract lag reads as a floor.
**Response.** Sentence removed; "Sims 27–50"; the 0.6 THz, ~9 dB and E∥ rows rewritten from Sim 50;
the abstract gives the best-fit lag (0.15–0.26 rad) and the smallest the traces allow (0.08 rad).
**Fixed.**

## Twenty-first round (independent reviewer, whole paper, 2026-11-03; verdict: major revision, no fatal flaw; node physics holds, the ranking of causes in Sec. 5.7 does not)

### 176. MAJOR: a zero-phase loss on a change of the sample breaks causality, and that choice decides the verdicts; the causal fit was demoted on an argument that fails (its 30 s lag is −0.004 rad).
**Response.** Agreed. The loss belongs to the tissue pulse (put on the whole trace it fits worse:
2.73/5.90/10.2 against 1.62/3.40/5.42, recomputed), so a causal reading cannot be set aside. New Sim 51
repeats every test for ten reference fits (no loss term; zero-phase; causal held above 1.5, 2 or
4 THz; each with a ≤ 1 or a + b ≤ 1), given as Table 2. No candidate passes or fails in all of them.
The abstract, Sec. 5.7, Sec. 6 and the conclusion now claim only what holds in every fit: the tissue
pulse lags at 0.6 THz by 0.15–0.28 rad at 2–5 min (−0.02 to 0.00 rad at 30 s) and loses more at high
frequency as the bias goes on; its cause is not identified. **Fixed.**

### 177. MAJOR: the causal check depends on where the loss stops rising.
**Response.** Table 2 includes the loss held above 1.5, 2 and 4 THz. With 1.5 THz and a + b ≤ 1 the
reference needs no shift (s = 0 at p = 0.10/0.11/0.045) and everything passes; with 2 or 4 THz the
slab causes fail and ring-resonance node states pass (with a + b ≤ 1, ions pass at 4 min). The text
said so. Superseded by item 185: those causal references fail the 30 s control. **Fixed.**

### 178. MAJOR: candidates were compared with the reference by a non-nested test, and the node verdict sat on the threshold.
**Response.** Sim 51 uses a nested test: the candidate's size v and the shift s are both free, and
the candidate suffices if s = 0 is not rejected (F(1, N_eff − n), n including v). The test is done
with N_eff from the full fit's residual and from the candidate's own; a verdict is pass or fail only
if both agree, else edge. Node states are now excluded by the forward bounds (every extreme and, for the
ring-resonance states only, a permittivity an aqueous gap cannot hold), not by p (Sec. 5.7). **Fixed.**

### 179. MAJOR: the bound a ≤ 1 decides verdicts; with a ≤ 1 − b warming fails; the loss term was a gain below 0.2 THz.
**Response.** Table 2 gives both bounds. The loss term is now 1 below 0.2 THz (`TwoPulse.loss_transfer`,
tested). Sec. 6 states that fits with a at its bound are outside the F distribution's assumptions.
**Fixed.**

### 180. MINOR: second-pulse lead, drift-reading lead, path figures, residuals and N_eff come from the fit without the loss term.
**Response.** Sec. 5.6 gives both: lead 0.019–0.044 ps without and 0.033–0.078 ps with the loss term;
drift reading 0.07–0.10 ps; paths 6–23 µm and 21–30 µm; residuals 9/18/20/23 % and 7/11/15/16 %;
N_eff 12–20. **Fixed.**

### 181. MINOR: "warming by a few kelvin" leaves out the range with Liebe as given.
**Response.** The abstract no longer ranks causes; Sec. 5.7 gives up to 12 K scaled and 30 K as given.
**Fixed.**

### 182. MINOR: the fitted loss is never explained.
**Response.** Sec. 5.6 and Sec. 6 state that it grows with bias time and belongs to the tissue pulse;
its origin is open. **Fixed.**

### 183. MINOR: the 0.6 THz passive bound is 0.13 dB, not 0.12; stale rows in this file's match table; internal notes and Sim numbers; data availability needs a URL.
**Response.** 0.13 dB; match-table rows rewritten from Sim 51. Sim numbers and the header note are
kept in this working copy for traceability and removed in the submission copy, which will carry the
repository DOI. **Partly fixed (submission copy).**

### 184. MINOR: novelty and fit; the waveform verdicts changed over rounds 18–20.
**Response.** The waveform section now claims only what is robust, and the abstract says the
waveforms are decomposed, not explained. The node result (law and slice bounds) carries the paper.
**Answered.**

## Twenty-second round (independent reviewer, whole paper, 2026-11-03; verdict: major revision, no fatal flaw; the causal references of Table 2 cannot test a candidate)

### 185. MAJOR: the six causal rows of Table 2 cannot test any candidate: at 30 s, before any change, they find a shift (p = 0.008), and the loss term's own delay (0.29–1.10 rad) exceeds the lag.
**Response.** Agreed, and used as a control. New Sim 52 adds a reference with a zero-phase floor fixed
at the 30 s loss (2.07 dB/THz) times a causal increment held above 1.5, 2 or 4 THz; it passes the 30 s
control (p = 0.41–0.70) in all six forms, and with it no shift is needed in 16 of 18 fits at 2–5 min.
It also fits better than a shift with a zero-phase loss. New Sim 53 then asks whether each candidate
supplies the change itself (increment and shift both zero, F(2, N_eff − 6)), with the 30 s scan as a
control. Table 2 is now Sim 53: no candidate supplies the change at all three bias times; a thicker
path passes at 2 min only. The finding is now positive: the change is a causal loss rising with
frequency (2.1–3.1, 2.8–4.5, 5.0–7.5 dB/THz at 2/4/5 min), steeper than any candidate's (≤ 1.6 dB/THz).
A second pulse crossing part of the tissue is covered by the a + b ≤ 1 rows, not modelled separately.
**Fixed.**

### 186. MAJOR: no node state supplies the lag; the paper cites the loss bound instead of the delay bound; passing states advance the pulse; ring states break Kramers–Kronig.
**Response.** Sec. 5.7 now states the delay bound: 0.063 rad at 0.6 THz with every extreme against a
lag of ≥ 0.15 rad, and a loss slope ≤ 1.6 dB/THz against ≥ 2.1 needed; in Sim 53 node states are edge
in every fit and pass in none. It says a frequency-independent ε is not causal, so those states are
bounds, not models. The 2 THz margin with the ×1.45 materials factor is stated as such in Sec. 4.6.
**Fixed** (the 2 THz factor is assumed, not solved: **partly**).

### 187. MAJOR: the conduction exclusion uses the 9 dB figure (130–1300 S/m mixes scaled and unscaled values).
**Response.** The forward bound is now 130–450 S/m at 0.6–2 THz on the measured optical thickness.
Inside the fits conduction no longer passes anywhere (Sim 53: edge or fail in all 18), and its loss is
nearly flat with frequency (0.2–0.9 dB/THz). **Fixed.**

### 188. MAJOR: stage drift is missing from the abstract and conclusion.
**Response.** With the causal increment no shift is needed, and a drift alone cannot make the tissue
pulse weaker at high frequency (Sec. 5.6). The drift reading is kept in Sec. 5.6 for the shift
description, with its lead (0.07–0.10 ps). **Fixed.**

### 189. MAJOR: stale cause claims ("9–18 % optically longer", "later tissue pulse"); statistics not labelled by reference.
**Response.** Deleted. The text speaks of a change of the tissue pulse; every number in Sec. 5.6 names
the description and Sim it comes from. **Fixed.**

### 190. MAJOR: the bar is only partly met: limits are textbook, 0.6 THz is fitted, 2 THz and selectivity are not explained, one traced sample.
**Response.** The abstract leads with the law and the slice bounds and states what the waveforms do
not explain (2 THz, selectivity, the second pulse's lead). Table 2 stays in the main text because it
now carries a positive result (no passive change supplies the measured loss). **Answered.**

### 191. MINOR: "70× smaller at the nominal geometry" is wrong.
**Response.** Now 11–16× smaller at internode 200 µm and packing 0.6. **Fixed.**

### 192. MINOR: "every passing state needs every extreme" is false; some passing ranges include zero.
**Response.** The sentence is gone; Table 2 (Sim 53) gives the sizes that pass, and none includes
zero, since no change at all is rejected (p ≤ 0.002). **Fixed.**

### 193. MINOR: "only node states at the ring resonance pass at every scan" is wrong.
**Response.** Removed with the old Table 2. **Fixed.**

### 194. MINOR: conduction passes from 5 S/m, not 20.
**Response.** Superseded: against the Sim 53 reference no conductivity passes. **Fixed.**

### 195. MINOR: E⊥ every-extreme values should use ×1.45 for materials.
**Response.** 3.6/13.9/21.7 dB. **Fixed.**

### 196. MINOR: the aqueous test (ε″ ≥ 3|ε′|) is too lenient.
**Response.** Aqueous now means σ ≤ 25 S/m, or 2 ≤ ε′ ≤ 5 with ε″ ≤ 3; none passes. **Fixed.**

### 197. MINOR: the loss is already ~2 dB/THz at 30 s; state the growth over it.
**Response.** The 30 s loss is the floor; every number is the increment over it. **Fixed.**

### 198. MINOR: stale responses (168, 170, 177, 178, REPORT row 55, "undelayed component", Fig. 9 band).
**Response.** 168: κ 2.1/5.6/6.5/10.0. 170 and REPORT row 55: 0.35/0.70. 177: ions pass at 4 min with
a + b ≤ 1, now superseded by 185. 178: the aqueous-gap argument names the ring states only. The
manuscript no longer says "undelayed". Fig. 9 now shows the lag over the ten Sim 52 references and the
Sim 53 loss increment. **Fixed.**

## Twenty-third round (independent reviewer, whole paper, 2026-11-03; verdict: major revision, no fatal flaw; the gas result came from a Mie error)

### 199. MAJOR: the gas-bubble fit came from the Mie approximation: the size parameter used Re(k) with a complex relative index, giving the gas interior a complex index; with x = k_h r the delay at r = 30 µm, v = 0.2 falls from 0.069 to 0.019 rad and no fit passes.
**Response.** Confirmed and fixed: `foldy_eps` now uses the complex host wavenumber (x = k_h r, m x = k_s r,
normalised by k_h³), with new tests (a sphere of the host's material changes nothing; the lossy result
joins the lossless one). Sims 54 and 54b were re-run. Gas in a swelling slice gives the steep loss
(3.3–8.5 dB/THz at r 30–80 µm) but ≤ 0.04 rad of delay; at nominal p it passes no fit. The gas claim is
withdrawn from the abstract, Sec. 5.7 and the conclusion. Sim 56 adds a thicker path to the gas (both
sizes charged): every scan passes, but zero gas passes too, so gas is not needed. **Fixed.**

### 200. MAJOR: the reference was tuned on its own control (the 30 s floor fitted to the 30 s scan), the H1 floor is applied to H2 and H3, and verdicts flipped with each reference.
**Response.** Agreed. Sim 55 varies the floor from 0 to 2.5 dB/THz (36 references with the hold
frequency and the bound). The loss rise above the 30 s value (1.7–2.8, 2.4–4.2, 4.4–7.1 dB/THz at
2/4/5 min) and the lag (0.16–0.28 rad) hold in all 36; only the split between loss and shift depends
on the floor. The paper now claims only these, and that the cause is not identified. **Fixed.**

### 201. MAJOR: the data do not pin down "causal"; the free shift cancels most of the increment's minimum-phase delay, and the hold frequency is a hidden parameter.
**Response.** The paper now describes the change as a loss rising with frequency plus a net lag, says
that whether the loss is causal and how it continues above ~1.7 THz decide only the split into loss and
shift, and covers the hold frequency (1.5/2/4 THz) in the 36 references. **Fixed.**

### 202. MAJOR: the Sim 54b loss slope in Table 2 was wrong.
**Response.** Table 2 rebuilt from the corrected Sims 54 (slopes at the best sizes, 3.3–8.5 dB/THz);
the Sim 54b row is removed with the gas claim. **Fixed.**

### 203. MAJOR: swelling to 1.67× is impossible in a 100 µm cuvette; 60–100 µm bubbles make a monolayer; boiling is ~10× slower with heat flow into PTFE; electrolysis gas forms at the electrodes.
**Response.** Sec. 5.7 states that 80 µm bubbles do not fit a 100 µm slice, that a fixed-gap cuvette
would not let it swell, that the electrolytic gas forms 7.5 mm from the beam; the boiling estimate is
removed and the warming paragraph notes heat flow into the plates within ~13 s. **Fixed.**

### 204. MINOR: "fit every scan" overreaches; "no change of thickness" contradicts the thicker-path passes.
**Response.** Both removed; the thicker path passes every calibrated test and is listed as such.
**Fixed.**

### 205. MINOR: the F tests are uncalibrated (systematic misfit, parameters at bounds, many tests, the Sim 54b grid widened after seeing Sim 54).
**Response.** Sim 55b calibrates them on 600 block-bootstrap 'no change' versions of the 30 s scan: the
nominal 5 % tests reject 19–38 % of the time; true 5 % needs p < 0.0027 (2-dof) and p < 0.0006 (shift).
At those thresholds the change is still significant in each scan (calibrated p 0.005–0.042), the shift is
never needed, and almost every candidate passes (Table 2 now shows calibrated passes). Hence the
waveforms do not identify the cause; node states and conduction rest on the forward bounds. Sec. 6
states the bootstrap's few blocks, the bounds and the number of tests. **Fixed.**

### 206. MINOR: N_eff 12–20 in Sec. 6 against 13–20 in Sec. 5.6.
**Response.** Sec. 6 no longer quotes it; Sec. 5.6 gives 13–20. **Fixed.**

---

## Twenty-fourth round (independent reviewer, whole paper, 2026-11-03; verdict: major revision, no fatal flaw; the bootstrap of Sim 55b was not built on no-change data)

### 207. MAJOR: the Sim 55b bootstrap drew its 'no change' scans from the full fit, which still holds a shift and a loss, so the rejections it counted were true detections; with null-based scans the nominal tests reject 2–4 %, and Table 2 and the claim that the waveforms cannot identify the cause fall.
**Response.** Agreed and confirmed. Sim 57 rebuilds the bootstrap per scan on the null fit of each test (s = 0 for the shift test; s = 0 and no loss for the no-change test) plus the full fit's residual in blocks of one and two correlation lengths, 300 each. The nominal 5 % tests reject 1–8 %; the 2-dof thresholds are p < 0.099/0.066/0.030/0.025 at 30 s/2/4/5 min. No change is rejected at calibrated p = 0.013/0.005/≤ 0.002 at 2/4/5 min and not at 30 s (0.45). Sim 58 reruns the candidate tests on every reference that passes the 30 s control (24) at these thresholds: node states with every extreme, ions, axoplasm conduction, warming (both versions) and gas in a swelling slice are rejected in all 24 at every scan; a thicker tissue path (6/8/10 %) passes 20/18/9. Table 2, the abstract, Sec. 1, Sec. 5.7, Sec. 5.8, Sec. 6 and Sec. 7 are rewritten; Sim 55b's thresholds are withdrawn in the text. **Fixed.**

### 208. MAJOR: the abstract and Sec. 7 said warming passes the calibrated tests; Table 2 showed it failing at 5 min.
**Response.** Superseded by item 207: warming is rejected in all 24 references; the text now says so. **Fixed.**

### 209. MINOR: candidates were tested against six references only, and 12 of the 36 references fail their own 30 s control.
**Response.** Sim 58 uses all 24 references that pass the 30 s control (floors 1.0–2.5 dB/THz). Sec. 5.6 gives the loss rise over those 24 (4.7–7.1 dB/THz at 5 min; 4.4–7.1 over all 36). **Fixed.**

### 210. MINOR: report the effect size against an empirical null, not only p values.
**Response.** Sim 57 (C): one traced unbiased pulse fitted with another (H1, H2, H3) gives a loss slope ≤ 0.12 dB/THz, a lag ≤ 0.03 rad and a change term of 0.2–2.8 % of the trace rms, against 3.1–7.4 dB/THz, 0.17–0.26 rad and 19–22 % at 2–5 min. In Sec. 5.6 and, in short, in the abstract. **Fixed.**

### 211. MINOR: leftovers: the Fig. 9 caption cites Sims 52–53 only; Sec. 4.6 says the waveforms need a delay rather than a loss; Sec. 1 says 'at the empty-cuvette delay'; the 30 s lag; the correlation length.
**Response.** The caption states which fits its bands come from and points to Table 2 (Sims 53, 54, 57, 58); Sec. 4.6 says the waveforms also need a lag; Sec. 1 says within ~0.1 ps of the empty-cuvette delay; the 30 s lag is −0.02 to 0.00 rad; the correlation length is 0.3–0.5 ps. **Fixed.**

### 212. MINOR: the 2 THz forward bound with every extreme (10 dB) is thin against ≥ 14 dB, and the ×1.45 materials factor is assumed at 2 THz.
**Response.** Sec. 4.6 already rests the exclusion on the states that come closest needing ε′ between −3.5 and −2.4 with ε″ = 0.1, which an aqueous gap cannot hold; node states now also fail inside the waveform fits at every scan. **Partly fixed** (the combined-materials case is not solved at 2 THz).

### 213. MINOR: novelty: both limits of the law are textbook; the forensic part rests on traced figures from one session.
**Response.** Sec. 1 states what is new (the interpolation between the limits, the slice bounds with every extreme, the waveform decomposition and the candidate tests). The single-session and tracing limits are in Sec. 6. **Stated.**

---

## Twenty-fifth round (independent reviewer, whole paper, 2026-11-03; verdict: major revision, no fatal flaw; the rejections came from capped sizes and a fixed timing)

### 214. MAJOR: node states fail the fits only because λ is capped at 1; with λ free they pass at λ = 2–4.5, so they are excluded once (forward bounds), not twice.
**Response.** Confirmed by Sim 59 (eight references, λ to 10): node states pass alone in 8 of 8 at every scan at λ = 2, 2.5–3.5 and 3.5–4.5 (2/4/5 min). Sec. 5.7 now says the shape of their change fits and its size does not; they are excluded by the forward bounds (2–4.5× every geometric extreme at once, 22–72× the change at internode 200 µm and packing 0.6). Abstract, Sec. 1, Table 2 and Sec. 7 changed. **Fixed.**

### 215. MAJOR: "only a thicker path fits" assumes no timing drift; with s free gas in a swelling slice passes.
**Response.** Confirmed by Sim 59 (candidate plus a free shift of both pulses, F(1, N_eff − 6) at the shift-test thresholds): gas in a swelling slice passes 8/6/8 (r 30 µm) and 8/5/6 (r 50 µm) of 8 with a drift of +0.03 to +0.06 ps growing with bias time and v 0.08–0.40. Ions, conduction and warming still fail at 4–5 min. Table 2 has a column for each test; the claim is now that the waveforms exclude the water's chemistry and temperature and need a lag, and that a thicker path, gas with a drift and node states beyond every extreme all fit; sizes decide. **Fixed.**

### 216. MAJOR: the passing candidate does not carry the headline loss slope; the slope depends on the amplitude bound.
**Response.** The abstract now leads with the lag and gives the loss as 1.7–7.1 dB/THz "depending on the fit"; Sec. 5.7 says the thicker path fits through its delay and the water loss that comes with it, and the argument that ions and conduction are "too flat" is now "too flat and too little lag", which Sim 59 supports (they fail also at twice the size). **Fixed.**

### 217. MINOR: bootstrap details: one start for the full fit; the no-change threshold reused for candidate tests with a fitted size and screened node states; the 5 % quantile from 600 draws is uncertain by ~±0.015.
**Response.** Stated in Sec. 6; verdicts near a threshold (the thicker path at 5 min) are called marginal. Alone, ions, conduction and warming (p ≤ 0.015) are clear of the thresholds; with a drift their rejection at 5 min is marginal (item 219). **Stated.**

### 218. MINOR: ions p 0.008 for 0.007; 3.1–7.4 against 1.7–7.1 unnamed; sim56.txt still says 'calibrated' against the withdrawn threshold; gas r 120 µm left out of Sim 58.
**Response.** Fixed in the text; Sec. 5.6 names 3.1–7.4 as the increment over the floor; sim56.txt carries a note that all its rows also pass at the Sim 57 thresholds; bubbles of 120 µm radius are larger than the 100 µm slice and are not physical here (Sec. 5.7 already says 80 µm does not fit). **Fixed.**

---

## Twenty-sixth round (independent reviewer, whole paper, 2026-11-03; verdict: no fatal flaw, major revision; three wordings beyond Sim 59)

### 219. MAJOR: "ions, conduction and warming cannot produce it, alone or with a timing drift" goes beyond Sim 59: with a drift they pass 2/1/0 and 2/2/0 of 8 and fail only at 5 min (p 0.013–0.014 against 0.023), and the drift test reuses the shift-test threshold for a different null.
**Response.** Agreed. The abstract, Sec. 1, Sec. 5.7 and Sec. 7 now say they fail alone at every scan and with a drift only at 5 min, marginally; Sec. 6 states that the drift threshold was calibrated for another null. Update (Sim 60): the drift test calibrated on its own null (candidate at its best size plus a drift, 200 block-bootstrap scans per case, eight references) rejects ions and conduction at 5 min at p ≤ 0.005 and warming at p ≤ 0.02 (≤ 0.03 at 4 min). Superseded by item 224: that calibration computed a different statistic in the synthetic scans; with the matched statistic (Sim 61) ions and conduction pass with a drift at 5 min in six of eight references, at sizes beyond saturated brine. **Fixed (see 224).**

### 220. MAJOR: the 1.7–7.1 dB/THz is not a property of the data: the thicker path passes with 0.4–0.7 dB/THz; "too flat" cannot be why ions fail; at 2 min a drift alone passes in 1 of 8.
**Response.** The abstract no longer gives the loss as a number; it says the loss and the lag trade off between fits. Sec. 5.6 says so and names the lag as the robust quantity. Sec. 5.7 says the candidates' change of loss and delay with frequency does not match, and that at 2 min a drift alone passes in one of eight. **Fixed.**

### 221. MINOR: node states with a drift pass at λ 0.75–5.5; alone at 2 min the range is 2–2.5.
**Response.** Table 2, Sec. 5.7 and Sec. 7 give both ranges ("near or beyond every geometric extreme"). The forward-bound exclusion stands (every extreme is 11–16 times the change at internode 200 µm and packing 0.6). **Fixed.**

### 222. MINOR: a positive drift makes the second pulse lead the empty cuvette further; a constraint on its arrival would test gas with a drift.
**Response.** Sim 60 (B): in every candidate's best fit the second pulse leads the empty cuvette, by 0.02–0.04 ps with a thicker path, 0.03–0.07 ps with gas (0.03–0.06 with a drift, not 0.07–0.10), 0.03–0.05 ps with node states, and 0.03–0.04 ps in the reference fit. The lead does not single out a cause; it stays unexplained (Sec. 5.6, Sec. 7). **Answered.**

### 223. MINOR: inconsistencies: 2 % vs 2.3 %; the quasi-static solver at 6.4 %; the every-extreme value at 1 THz; the selectivity row; stale REFEREE entries.
**Response.** Sec. 4.4 says 2.3 %; Sec. 4.1 gives 0.1–6.4 % (largest at 0.1 THz); the every-extreme values recomputed from `results/sim30/sim36_signedband.json` are 1.5/2.7/10/18 dB (were 1.5/2.9/10/17); the Sec. 5.8 selectivity row names only the causes that fit; response 217 and the match-table rows on the bypass and on timing are updated. The figure file names keep 'bypass' (file names only). **Fixed.**

## Twenty-seventh round (independent reviewer, whole paper, 2026-11-03; verdict: no fatal flaw, major revision; the Sim 60 calibration used a different statistic and only failing candidates were calibrated)

### 224. MAJOR: the Sim 60 bootstrap computes a different statistic from the observed one: the observed full fit may drop the candidate (min(mf, s_ref)), the synthetic full fits may not.
**Response.** Agreed. A tested helper (`null_calibrated_test` in `thznerve.waveform_fit`) now computes one statistic on the data and on every synthetic scan: the null at its fitted size, against the better of the candidate with the loss increment and no candidate with the increment, N_eff from the winning fit's residual, the same multi-start fits on both (Sim 61). Calibrated this way the null rejects 10–25 % at nominal 0.05 (not 1–7 %), so the nominal tests are lenient, not conservative. Ions and conduction with a drift at 5 min then pass in six of eight references (calibrated p 0.040–0.148), at 10–40 S/m added to all tissue water or 25–100 S/m added to the axoplasm. Sec. 5.7 now excludes them by size (1.4 S/m physiological, ~25 S/m saturated NaCl) and by the forward bound (130–450 S/m), not by the waveforms; the abstract, Sec. 1, 5.6, 5.8, 6 and 7 say so. The reviewer's own rerun (thresholds 0.039–0.085, calibrated p ≤ 0.020) took N_eff from the candidate fit's residual when the reference fit won and fitted the synthetic scans from one start; both differ from the data. **Fixed; claim narrowed.**

### 225. MAJOR: only the failing candidates were calibrated; the passing ones were judged at thresholds built for other nulls.
**Response.** Sim 61 calibrates every Table 2 cell that is used: thicker path and node states, alone and with a drift, at 2, 4, 5 min; gas r 30/50 µm with a drift; ions, conduction and warming alone at 2, 4, 5 min and with a drift at 5 min (warming also 4 min); 248 cases, eight references, 60–100 synthetic scans each. Results: thicker path 6/6/4 alone, 8/8/8 with a drift; node states 8/7/8 and 8/7/5; gas 8/6/8 and 8/6/6 with a drift; ions and conduction 1/6/6 alone; warming 1/3/0 alone and 0 of 8 with a drift at 5 min. Sim 63 recalibrates the no-change and shift-alone tests the same way: no change is rejected at 2, 4, 5 min in all eight references (calibrated p 0.016–0.049) and in none at 30 s; a shift alone is rejected in 22 of 24. Table 2 is rebuilt from Sim 61; the bounded-size column judged at uncalibrated thresholds is removed, and gas alone, not calibrated, is no longer called rejected. **Fixed.**

### 226. MINOR: warming p numbers in Sec. 5.7 are wrong (Sim 58: up to 0.046 at 2 min).
**Response.** The paragraph now gives the Sim 61 calibrated results (fails at 5 min in all eight; passes at 4 K in one to three at 2 and 4 min). **Fixed.**

### 227. MINOR: the largest node delay is 0.063 rad (ε = −3.0 − 0.1j), not 0.061.
**Response.** Corrected in Sec. 5.7 and response 186. **Fixed.**

### 228. MINOR: Table 1 is tighter than today's rerun.
**Response.** Table 1 now gives the rerun values: cylinder flux balance ≤ 7.5 × 10⁻⁷, null test ≤ 1.4 × 10⁻⁷. **Fixed.**

### 229. MINOR: the second-pulse lead covers 2 of 8 references and includes rejected fits; 0.02–0.04 vs 0.02–0.08 ps.
**Response.** Sim 62 takes every fit that passes the calibrated tests of Sim 61 (all eight references, every candidate): the lead is 0.030–0.058, 0.031–0.049 and 0.020–0.039 ps at 2, 4 and 5 min, varying more between scans than between candidates. Sec. 5.6 gives these ranges; Sec. 5.6 (what the second pulse could be) and Sec. 7 say 0.02–0.06 ps at 2–5 min (0.07–0.08 ps at 30 s, reference fits). **Fixed.**

### 230. MINOR: the fits depend on their starting point (same case p 0.049 in Sim 59, 0.066 in Sim 60).
**Response.** A tested helper (`multistart`) takes the best of the default start, the given start and the nested model's solution; Sims 61–63 use it on the data and on every synthetic scan. **Fixed.**

### 231. MINOR: stale text: "1.5/10/17 dB"; response 219's "nominal test was conservative".
**Response.** The match table now says 1.5/10/18 dB at 0.6/2/3 THz; response 219 points to 224. **Fixed.**

## Twenty-eighth round (independent reviewer, whole paper, 2026-11-03; verdict: no fatal flaw, major revision; the bootstrap doubles N_eff, the size is held fixed, Monte Carlo error decides verdicts)

### 232. MAJOR: the block bootstrap gives synthetic scans about twice the data's N_eff (median ratio 2.03 in Sim 61).
**Response.** Agreed. Sim 64 replaced it with a phase-randomised surrogate of the residual (same amplitude spectrum, so the same autocorrelation), with the data's N_eff used on every synthetic scan (tested helpers `phase_surrogate` and `size_fitted_test`). Before using it we tested it on data with no change (Sim 65): one traced unbiased pulse fitted with another from a different report, six ordered pairs, a ≤ 1 and a free. It rejects "no change" or "shift alone" in 13 of 24 comparisons, "no change" in 10 of 12. Every model leaves a systematic misfit, so with enough points any two traces differ significantly; no significance test on the traced figures keeps its nominal rate. All p-value verdicts (Sims 55b, 57, 60, 61, 63, 64) are withdrawn. Fits are now judged by effect size: D, the rms of what a restricted fit leaves unexplained against the fit with the loss increment and a free shift, in % of the trace rms, against the largest D between traced unbiased pulses (1.63 % no change, 1.18 % shift alone; Sim 66). Sec. 5.6, Table 2, Sec. 5.7, 5.8, 6, 7, 1 and the abstract are rewritten. **Fixed; method replaced.**

### 233. MAJOR: holding the size fixed in the synthetic scans favours passing; Sec. 6 says the opposite.
**Response.** Agreed; Sec. 6's sentence was wrong and is removed. Sim 64 fitted the size inside the statistic, over the grid, on the data and every synthetic scan. With significance tests withdrawn (item 232), Sim 66 takes each candidate at its best size on the grid, alone and with a drift, over a plausible and a wide range, and reports what it leaves. **Fixed.**

### 234. MAJOR: Monte Carlo error decides verdicts (71 of 248 cases with p_cal 0.03–0.07; warming with a drift is the shift-alone test under another name).
**Response.** Sim 64 used up to 400 draws with a 95 % Clopper–Pearson interval and an undecided verdict; Sim 65 then showed the null itself is not valid. Effect size has no Monte Carlo error. Table 2 now says warming's best size with a drift is zero, so its drift column is the shift alone. **Fixed.**

### 235. MAJOR: the grid is too coarse for the conductivity exclusion; electrolysis protons (~55–85 mM, ~2–3 S/m) undercut the saturated-NaCl ceiling; the 9 dB forward bound contradicts Sec. 5.7's own point that the dip depth is not the quantity to bound.
**Response.** Sim 66 searches 0–5 S/m in 1 S/m steps (plausible, covering the protons) and 0–60 S/m in the tissue water and 0–150 S/m in the axoplasm (wide). Conduction leaves 6.9–22 % at 0–5 S/m and 6.8–18.7 % at its best wide size, so it does not reproduce the change at any size; the exclusion no longer rests on a ceiling. Sec. 5.7 gives the proton estimate, and the forward-bound sentence is removed. **Fixed.**

### 236. MAJOR: Sec. 5.7 "How the candidates were tested" describes the withdrawn test.
**Response.** Rewritten for Sim 66: references (eight), grids (plausible and wide), alone and with a drift, criterion and yardstick; Table 2's caption repeats them. **Fixed.**

### 237. MINOR: "calibrating as above can only raise those p" is false.
**Response.** The sentence is removed with the significance tests. **Fixed.**

### 238. MINOR: "10–25 %" (Sim 61 shows 0.00–0.35), "~6–12 blocks" (about 20), "36 fitting assumptions" attached to calibrated tests.
**Response.** Sec. 5.6 now describes the bootstrap only as withdrawn, with a null rejection of 3–23 % at nominal 5 % (5th–95th percentile); the block count is removed; "36 fitting assumptions" now goes only with the loss and the lag (Sim 55), and the effect sizes with the eight references of Sim 66. **Fixed.**

## Twenty-ninth round (independent reviewer, whole paper, 2026-11-03; verdict: fatal flaw in Sec. 5.7, the effect-size yardstick ignores fitting noise; node physics sound; major revision)

### 239. FATAL: the yardstick (unbiased pairs, N_eff 100–230) is on the wrong scale for biased traces (N_eff 13–20, full fits leaving 6–15 %); a true cause would often exceed it, so "no tested cause reproduces the change" does not follow.
**Response.** Agreed; the claim is withdrawn. Sim 68 (tested helpers `effect_size`, `effect_noise`) builds a threshold per restricted fit: the largest D over 30 synthetic traces (the restricted fit's curve plus a phase-randomised copy of the full fit's residual, both models refitted, the candidate's size held at its value on the data, which favours it), added in quadrature to the tracing D built like each reference. A candidate is excluded in a reference if D exceeds the threshold at one of the three scans (about 9 % false exclusion for a true cause). Result: no change and a shift alone are excluded in all eight references (2.5–3.6 and 1.8–3.7 times); warming, ions to 60 S/m and axoplasm conduction to 150 S/m in all eight, alone and with a drift (1.5–3.7 times); gas in a swelling slice alone in all eight. Not excluded: a thicker path (alone in 7 of 8, with a drift in 8 of 8 but by as little as 1.07 times), gas with a drift (6 of 8; 0.85–0.96 in the other two) and node states (6 and 4 of 8, only at 1–6 times every extreme). The abstract, Sec. 1, 5.6, Table 2, 5.7, 5.8, 6 and 7 now say the waveforms exclude warming and conduction, cannot decide between a thicker path and gas with a drift, and accept node states only at sizes the forward bounds exclude. "The yardstick is generous" is removed. **Fixed; claim narrowed.**

### 240. MAJOR: the full model also misses the data (10–15 % at 2–5 min, against 1.7–3.4 % between tracings), yet the abstract says each waveform "fits as".
**Response.** Sec. 5.6 now gives the full fit's residual (10.0–15.2 % at 2–5 min, 6.5–7.1 % at 30 s) against 1.8–9.6 % between unbiased tracings built like the references, and calls the decomposition approximate; the abstract says "fits approximately"; Sec. 6 repeats it. **Fixed.**

### 241. MAJOR: the node-state row uses three states chosen at λ = 1 alone; rescreening all 80 gives D 3.62 with a drift at 2 min.
**Response.** Sim 68 screens all 80 states over λ = 0.5–6 for each reference, scan and mode and refits the best three from several starts: D with a drift 3.6–4.8 % at 2 min, at λ 1–6. Node states are not excluded by the waveforms in two to four references; Sec. 5.7 now excludes them by the forward bounds only. **Fixed.**

### 242. MINOR: the second-pulse lead (Sim 62) covers only fits that passed the withdrawn test; the Sim 57 slope ≤ 0.12 dB/THz holds only with a ≤ 1.
**Response.** The lead is now taken from the best fit of every candidate in Sim 68: 0.020–0.069 ps at 2–5 min (full fits 0.027–0.044 ps); Sec. 5.6 and Sec. 7 say 0.02–0.07 ps. The slope sentence says "with a ≤ 1 (up to 0.58 dB/THz with a free)". **Fixed.**

### 243. MINOR: numbers: 6.9–22.0 vs 22.62; plausible node λ at the edge; two loss ranges (Sim 57 vs Sim 55); conduction optima at 5 S/m.
**Response.** Table 2 and the conduction paragraph are rebuilt from Sim 68, which gives the best sizes found (5–40 S/m in the tissue water, 10–100 S/m in the axoplasm; node λ 1–6). Sec. 5.6 labels the Sim 57 range as the increment over the floor in the eight references, distinct from the Sim 55 rise above the 30 s value over 36. **Fixed.**

### 244. MINOR: the 30 s scan exceeds the yardstick in all eight references.
**Response.** Against the calibrated threshold, no change at 30 s is within it in all eight (0.3–0.9 times); a shift alone exceeds it in two (up to 1.3 times), the part of the change already present at 30 s. Sec. 5.6 says so. **Fixed.**

### 245. MINOR: the tracing pairs are built differently from the references (no floor, no a + b ≤ 1).
**Response.** Sim 68 builds them like each reference (floor applied, a ≤ 1, a + b ≤ 1 where used): largest D 2.6–4.3 % without a shift, ≤ 0.08 % with one. **Fixed.**

## Thirtieth round (independent reviewer, whole paper, 2026-11-03; verdict: no fatal flaw, major revision; the threshold rule differs between script and text, the noise is too small, the node paragraph contradicts the transfers used)

### 246. MAJOR: the Sim 68 script judges at the 95th percentile of 30 draws while the paper uses the largest; Table 2 is in no file.
**Response.** Sim 69 applies one rule in the script and the text: per scan the 1 − 0.05/3 quantile of 100 draws, combined with the tracing D in quadrature; excluded if D exceeds it at one of 2, 4, 5 min. Table 2 is now printed by `sim69_threshold_rule.py` into `results/sim69/sim69.txt`. The rule is a tested helper (`threshold_verdict`), which reproduces every Sim 69 verdict. **Fixed.**

### 247. MAJOR: the full fit's residual understates the noise (5 parameters against N_eff 13–20), and the largest of 30 draws is itself noisy; the false exclusion rate is 13–22 %, not 9 %.
**Response.** Agreed. Sim 69 scales the residual by the factor at which noise injected into the full fit's own curve leaves a residual as large as on the data: 1.07–1.24 (sqrt(N_eff/(N_eff − 5)) gives 1.16–1.29), tested helper `noise_scale`. 100 draws per fit, per-scan level 1 − 0.05/3 so about 5 % per reference over three scans; counted on the draws themselves, 0–6 %. Each verdict is bootstrapped over the draws (200 resamples); a flip rate above 10 % is reported as borderline (Table 2, in brackets). **Fixed.**

### 248. MAJOR: the margins for warming and conduction are overstated ("1.5 times or more").
**Response.** After the correction they are 1.33–3.19 (ions with a drift lowest, 1.33), none borderline, in all eight references; Sec. 5.7, 5.8 and 6 now say 1.3. The exclusions stand. **Fixed.**

### 249. MAJOR: the node paragraph gives a 1.6 dB/THz slope and ε′ −3.5 to −2.4, against 4.6–5.2 dB/THz and −1.0 to −5 in the transfers the fits use; λ = 1 is the forward bound itself, so the bound cannot exclude it.
**Response.** Agreed; recomputed over 0.2–1.7 THz at every extreme. The steepest slope is 5.2 dB/THz, from ε′ −1.7 to −1.4 with ε″ = 0.1, which advance the pulse (−0.08 to −0.11 rad at 0.6 THz); the largest delay, 0.063 rad (ε = −3.0 − 0.1j), comes with 1.1 dB/THz. No state gives both. Sim 69's best states are ε′ −1.3 to −1.0 with a drift (within the threshold in all eight references, 0.83–0.98), −5 to −2.7 alone, or a 1000 S/m node. The exclusion now rests on what an aqueous gap can hold (ε″ ≫ |ε′|, Sec. 4.6), not on the waveforms or the λ range; Sec. 1, 5.7, 5.8 and 7 say so. **Fixed.**

### 250. MINOR: "thicker path excluded in every reference with a drift" is not robust (flips in 31 % of resamples).
**Response.** With the corrected noise it is within the threshold in three references with a drift (0.82–0.95) and one alone (0.96, borderline); the text calls it not excluded. **Fixed.**

### 251. MINOR: leftovers: "Sims 27–67", Fig. 9 caption (Sim 66), a 0.4–0.7 dB/THz slope found nowhere, "bulk effects" excluded while gas is not, 0.96 vs 0.97, the What-matches row citing the withdrawn 1.63 %.
**Response.** All corrected (Sims 27–69; Sim 69; phrase removed; "warming and extra conduction"; the gas numbers come from Sim 69, 0.75–1.00 within and 1.02–1.26 above; the row now cites 2.2–3.3 times the threshold). **Fixed.**

### 252. MINOR: disclose that D compares fits that are not nested and that sizes are free per scan.
**Response.** Both are stated in Sec. 5.6 and Sec. 6, with the note that both favour the candidates. **Fixed.**

## Thirty-first round (independent reviewer, whole paper, 2026-11-04; verdict: no fatal flaw, major revision; the 30 s control and the surviving candidates depend on the loss floor, gas searched beyond its model, warming rests on Liebe above 1 THz)

### 253. MAJOR: the floors 1.0–2.5 dB/THz were chosen so the 30 s scan shows no change, so the 30 s result is not a control; at floor 0 the thicker path and gas with a drift are excluded and no change is already rejected at 30 s.
**Response.** Agreed. Sim 70 adds four references with floors of 0 and 0.5 dB/THz (a ≤ 1, a + b ≤ 1; tracing pairs built like each), Sim 69 rule. There, no change at 30 s exceeds its threshold by 1.2–2.2 times, and every bulk candidate is excluded in all four references (thicker path with a drift 1.21–1.86, grid extended to 51 %; no change at 2–5 min 3.4–4.0). Sec. 5.6 now says the eight references read the 30 s loss as a reference artefact, so their 30 s result holds by construction, and that the tracing pairs argue against that reading. Table 2 has a column for floors 0–0.5; the abstract, Sec. 1, 5.7, 5.8, 6 and 7 say the thicker path remains only if the 30 s loss is an artefact of the reference. **Fixed.**

### 254. MAJOR: gas with a drift is "not excluded" only at v up to 0.4, beyond the dilute coherent-medium model (v ≤ 0.2).
**Response.** Agreed. Sim 70 searches gas only to v = 0.2: with a drift it is excluded in all eight references (1.08–2.00) and in the four at floors 0–0.5 (1.56–2.09); its best sizes sit at the model's limit, so the text says dense gas is not tested. **Fixed.**

### 255. MAJOR: the warming exclusion rests on Liebe's temperature dependence above 1 THz.
**Response.** Ellison (2007) coefficients could not be obtained, so Sim 70 brackets how warming acts above 1 THz with three versions, each scaled to the measured 1 THz change (Millon): A, Liebe's full change (×2.31); B, only the slow relaxation warms (×1.04); C, only the fast one (×3.26), which raises the loss slope. All three are excluded in every reference, alone and with a drift: 1.2–3.2 times the threshold at floors 1.0–2.5 and 1.6–3.8 at floors 0–0.5. **Fixed.**

### 256. MINOR: "cannot produce the change" reads as "absent"; bias heating is present and combinations were not tested.
**Response.** Now "cannot alone produce" / "excluded as sole causes"; the warming paragraph says warming is likely present and combinations are not tested (also Sec. 6). **Fixed.**

### 257. MINOR: numbers: second-pulse lead from two sources; 11–96 and 11–16× vs 11.7–16.5; node slope and delay in no results file; the abstract joins 36 assumptions to the eight-reference ratio.
**Response.** The lead sentence gives both sources (Sim 55; Sim 68: 0.027–0.044 and 0.072–0.084 ps); 12–99 and 12–17×; Sim 70 prints the node numbers to `results/sim70/sim70.txt` (largest delay 0.064 rad at ε = −3.0 − 0.1j with 1.10 dB/THz; steepest slope 5.24 dB/THz at −1.7 − 0.1j, delay −0.084 rad); the abstract now says twelve references for the ratio. **Fixed.**

### 258. MINOR: Sec. 4.6 "ion collisions make ε″ ≫ |ε′|" is false for water itself.
**Response.** Replaced by the ionic Drude increment: its loss exceeds its negative ε′ by 1/(ωτ) ≈ 15 at 0.6 THz (τ ≈ 17 fs); pulling ε′ ≈ 4.5 below −1 needs ~3 × 10³ S/m and brings ε″ ≳ 80. **Fixed.**

### 259. MINOR: the thicker path with a drift reaches 30 %; the mass covers only 6–9 %.
**Response.** Sec. 5.7 adds 30–66 µm, ~12–26 mg at 30 %; at floors 0–0.5 the grid runs to 51 % (best 12–36 %). **Fixed.**

## Thirty-second round (independent reviewer, whole paper, 2026-11-04; verdict: fatal flaw in the warming exclusion as worded; node part sound; major revision)

### 260. FATAL: warming models B and C are not a bracket on the fit; a mixture w B + (1 − w) C, each matched to the 1 THz change, fits better and is within the threshold in one reference at floors 1.0–2.5.
**Response.** Agreed; the claim is withdrawn. Sim 71 fits w (0–1) and the temperature rise (0–40 K) in all twelve references. With floors of 0–0.5 the mixture is firmly excluded (1.49–1.89). With floors of 1.0–2.5 it is not firmly excluded alone (within the threshold in one reference, borderline in four; best w 0.6–1.0, 8–24 K), and is excluded with a drift (8/8, one borderline). The abstract, Sec. 1, 5.7, 5.8 and 7 now say warming is excluded only if the 30 s loss belongs to the slice. Rønne et al. (1997) is cited as the measurement that would fix w. **Fixed.**

### 261. MAJOR: dilute gas tested at one radius (30 µm); 50 µm fits better.
**Response.** Sim 71 frees the radius (20, 30, 50, 75, 100 µm; v ≤ 0.2). Alone, gas is firmly excluded everywhere (1.39–2.13; 1.85–6.67 at floors 0–0.5). With a drift and floors of 1.0–2.5 it is within the threshold in three references and borderline in four (best r 20–50 µm, v at the 0.2 limit); at floors 0–0.5 it is firmly excluded (1.33–2.41). The text now lists gas with a drift among the causes not excluded at floors 1.0–2.5. **Fixed.**

### 262. MINOR: "none borderline" overstates certainty: one threshold carries ~11 % Monte Carlo error.
**Response.** A verdict is now also borderline if its largest ratio lies within 0.87–1.15; "firmly excluded" means excluded and not borderline. Table 2 is rebuilt with this rule for every row by `experiments/exp/table2.py` (results/table2.md). Conduction (≥ 1.33) and warming models A, B, C (≥ 1.20) stay firmly excluded. **Fixed.**

### 263. MINOR: the noise model assumes stationary noise; envelope-shaped noise moves ratios in both directions.
**Response.** Stated in Sec. 5.6 and Sec. 6 as a further reason not to read ratios near 1 as firm. **Fixed.**

### 264. MINOR: the floors 0–0.5 30 s test used the 0.95 level.
**Response.** Recomputed at 1 − 0.05/3 (Sim 71): no change 1.12–1.47 (borderline at floor 0.5), a shift alone 1.37–1.92 (firm). Sec. 5.6 gives these. **Fixed.**

### 265. MINOR: the thicker path with a drift sits at the 30 % grid edge in four references.
**Response.** Sim 71 searches to 51 %: best 6–33 %, inside the grid; with a drift within the threshold in two references and borderline in three. **Fixed.**

### 266. MINOR: stale text (Sec. 6 "1.3 times or more"; "Sims 27–69"; node delay 0.063 vs 0.0640; the What-matches rows).
**Response.** All corrected (Sims 27–71; 0.064 rad, Sim 70; Sec. 6 rewritten; the rows below rewritten). **Fixed.**

### 267. MINOR: the abstract should say the paper's own tracing comparisons favour reading the 30 s loss as part of the slice, under which no tested cause fits and the change precedes the reported onset.
**Response.** The abstract, Sec. 1 and Sec. 7 now say so. **Fixed.** (Withdrawn in round 33, item 268.)

## Thirty-third round (independent reviewer, whole paper, 2026-11-04; verdict: fatal flaw, a slice of uneven thickness was never tested and fits; major revision)

### 268. FATAL: every slab candidate has one uniform thickness; a slice whose thickness varies across the beam reproduces the change at floors 0–0.5, so "no tested bulk cause can" rests on omitting the cheapest cause.
**Response.** Agreed; the claim is withdrawn. Sim 72 adds uneven thickness (tested helper `uneven`: the coherent mean of the uniform-slab transfer over a Gaussian spread of thickness, mean −30 to +51 %, rms 0–40 %; for a pure phase it gives e^{iφm − φ²σ²/2}) in all twelve references, 2/4/5 min and 30 s, with the Sim 71 rule. It is excluded in none of the eight references with floors of 1.0–2.5 (0.51–0.83 of the threshold) and is within or at the threshold in all four with floors of 0–0.5 (0.67–1.11, one borderline); without a drift its best sizes are a mean 6–18 % thicker with an rms of 4–25 %. At 30 s an rms of 8–16 % fits within the threshold (0.24–0.71 with a drift). The reviewer also pointed to H4, which we had listed as not seen: it is open on arar.sci.am. Sim 73 traces it: one sample, probably the H1–H3 one (unbiased pulse within 5.6 % rms), carries after 6 min a second pulse of b = 0.24 at the empty-cuvette time, and the low-frequency resonance H4 reports for it is the cancellation of the two pulses (0.56 THz; model 0.56–0.59); the other sample fits a path 10–15 % thinner. The abstract, Sec. 1, 5.6, 5.7, 5.8, 6 and 7 now say that the waveforms do not single out a cause, that a slice whose geometry in the beam changes fits every reference with no change of the tissue's permittivity, and that H4 shows such a change directly. Thickness mapping is added to the decisive measurements. **Fixed.**

### 269. MAJOR: verdicts disagree where one model contains the other (mixed warming alone not excluded, with a drift excluded; warming A with a drift is the shift row reprinted).
**Response.** Table 2 now judges each candidate in each reference by the better of its fits alone and with a drift (`table2.py`), with one row per candidate and a column counting where the drift fit was better; "no change" and "shift only" stay as the two rows that test the change and its shape. **Fixed.**

### 270. MAJOR: the abstract lost "alone" (settled in item 256).
**Response.** Restored: "Extra conduction alone cannot produce it". **Fixed.**

### 271. MINOR: a repeated phrase in the E∥ row; Sec. 6 residuals from two sources; "36 fitting assumptions" joined to the twelve-reference ratio; a second DataPort entry; H4.
**Response.** Row rewritten; Sec. 6 now quotes the full fit (10–15 % at 2–5 min, Sec. 5.6); the abstract sentence no longer mentions the 36 assumptions; both DataPort pages are the same dataset (doi:10.21227/h17h-tp94) and require a subscription; H4 is traced (Sim 73) and cited with its DOI. **Fixed.**


## Thirty-fourth round (independent reviewer, whole paper, 2026-11-04; verdict: fatal flaw, the conclusion overstates what the uneven-thickness fit and H4 show; major revision)

### 272. FATAL: the uneven-thickness model is the mean over infinitely many patches; one slice under a 5 mm beam does not realise it, so "a slice of uneven thickness can, in every reference" and "the data point to a slice changing shape" rest on a reparameterised loss plus delay.
**Response.** Agreed. Sim 74 lays random thickness fields (the fitted sizes, correlation length 0.25–2 mm) under a mode-matched beam 5 mm across, fixed or scaling with wavelength, 200 slices each. One slice departs from the mean by an rms of 0.7–4.6 dB at 0.6 THz and 2.0–13.8 dB at 1.5 THz for patches of 0.5–2 mm, as large as the mean change itself (−1.6 to −4.0 and −8 to −13 dB); only patches of ≲ 0.25 mm average out, and there a patch is no wider than a wavelength, where the local-slab transfer fails. Sec. 5.7 states this after the method; Sec. 6 adds it as a limitation. The abstract, Sec. 1, 5.7, 5.8 and 7 now say that the waveforms exclude extra conduction alone but do not identify the cause, list uneven thickness with the other causes not firmly excluded, and drop "the data point to a slice changing shape" and "the simplest reading". **Fixed.**

### 273. FATAL: H4's "6 min" scan of sample 1 is the H3 5 min scan or a near copy; it is not independent evidence, and "continuing the sequence" and "shows directly" count it twice.
**Response.** Confirmed by Sim 75: after the best shift and scale it matches H3's 5 min scan within 10.1 %, against 27–47 % between any other two biased scans, and its unbiased pulse matches those of H1–H3 within 3.9–5.0 % (they match each other within 1.8–3.3 %). Sec. 5.6 now says sample 1 is the sample of H1–H3 and its biased scan H3's or one taken soon after, not independent evidence; the H4 sentence is removed from the abstract; "continuing the sequence" and "shows directly" are deleted; the onset row of Sec. 5.8 says it adds no time point. Sample 2 is described as the one new sample. **Fixed.**

### 274. MAJOR: the lead of the second pulse (0.02–0.08 ps), called unexplained, is within how well the empty-cuvette pulse is reproduced between reports.
**Response.** Confirmed by Sim 75: H1's and H4's empty-cuvette pulses differ by 9 % rms and 0.025 ps; with H4's as K the fits without a loss term put the second pulse 0.000–0.015 ps early alone (one pixel 0.015 ps) and from 0.015 ps early to 0.011 ps late with a drift, at 30 s as at 2–5 min, b changing by ≤ 0.004. Sec. 5.6 says so; "which we cannot explain" is removed from Sec. 5.6 and Sec. 7. **Fixed.**

### 275. MINOR: Sec. 5.6 l. 493 "At 2–5 min it leaves 21.6–24.2 %" now reads as the uneven fit.
**Response.** Rewritten: "With these floors no change leaves 21.6–24.2 % at 2–5 min". **Fixed.**

### 276. MINOR: node λ given as 2–6 and as 1–6; uneven thickness called "every reference" though excluded (borderline) in one with floors 0–0.5 and at the grid edges; one wording rule for all rows.
**Response.** Sec. 5.7 gives 2–6 for the better fit of each reference (Table 2) and 1–6 over all fits; uneven thickness is "not firmly excluded in any reference", with the fourth floor-0–0.5 reference "excluded but borderline" and the drift fits at the grid edges (mean −30 %, rms 40 %) named. Every surviving candidate is now described as "not firmly excluded". **Fixed.**

### 277. MINOR: b at 5 min given as 0.21 and as 22 %; a stray "0.6,."; H4 sample 2's full fit cannot add gain, so its D = 0 means nothing.
**Response.** 21 % throughout (fits give 0.209–0.219); typo fixed; Sec. 5.6 says the full fit does no better than no change for sample 2, so no effect size is formed. **Fixed.**


## Thirty-fifth round (independent reviewer, whole paper, 2026-11-04; verdict: no fatal flaw; major revision, six claims go further than the result files)

### 278. MAJOR: "the fits without a drift already pass" is false for uneven thickness: firmly excluded alone in two references (floors 0 and 1.0, a + b ≤ 1; 1.16–1.17), borderline in a third (1.14); there only grid-edge drift fits (rms 40 %) pass.
**Response.** Agreed. Sec. 5.7 now says that without a drift it passes in nine of twelve references, gives the two exclusions and the borderline, says that in those two only the grid-edge drift fit (a mean 12–27 % thinner, rms 40 %) passes, the slice Sim 74 shows a single slice realises least, and counts those two verdicts as unresolved; Sec. 5.7's thicker-path paragraph, Sec. 6 and Sec. 7 follow. REPORT row 76 ("not needed") is corrected. **Fixed.**

### 279. MAJOR: the abstract's "each biased waveform fits …" does not hold for H4's sample 2, the only independent sample.
**Response.** The abstract now says "the biased waveforms of one sample (three reports)" and adds that a second sample transmits more and earlier, with no second pulse. **Fixed.**

### 280. MAJOR: sample 2's "16 % without a change" is a local minimum; unbounded τ lets the empty pulse sit at the tissue delay; "earlier and stronger" is another reading than a thinner path.
**Response.** `TwoPulse.fit` takes a bound on the second pulse's arrival (t2_max, tested). Sim 76 refits H4 with |s + τ| ≤ 0.2 ps over the base tracing and the eight tracing shifts: no change now leaves 15.8 % (15.8–15.9), a path 10–16 % thinner 3.7 %, and a tissue pulse 1.2 times stronger and 0.06 ps earlier with no change of shape 6.4 % (6.1–6.6). Sample 1 is unchanged. Sec. 5.6 gives all three and calls the thinner path the best of these readings, not the only one; Sec. 5.7 and 7 say the sample "transmits more and earlier, as a thinner path would". **Fixed.**

### 281. MAJOR: "excludes extra conduction" contradicts the 1000 S/m node, which passes within the threshold; the abstract lists node states only as negative permittivities.
**Response.** Abstract, Sec. 1, 5.7 and 5.8 now say bulk conduction (tissue water, axoplasm); the abstract and Sec. 7 name node states "as nearly lossless negative permittivities or 1000 S/m conductors". **Fixed.**

### 282. MAJOR: the H4 resonance check is circular (the traced spectrum against a model fitted to the same waveform), counts H3's scan twice, and never quotes H4's frequency.
**Response.** H4 gives no frequency ("deep resonance in the low-frequency region", p. 20); its spectrum (Fig. 4b) has a broad minimum near 0.55 THz and its index (Fig. 5a) jumps at 0.60 THz. Sec. 5.6 now quotes these, calls the match a consistency check on what is H3's scan or one close to it, and notes the 0.557 THz water line; the abstract sentence is removed; Sec. 1, 5.8 and 7 say "consistent with the cancellation". **Fixed.**

### 283. MINOR: Sim 75 tested the lead only without a loss term.
**Response.** Sim 76 repeats it with the causal loss term: with H4's empty pulse the second pulse is within 0.009 ps of the empty-cuvette delay at 2–5 min (H1's: 0.020–0.039 ps early) and 0.011–0.033 ps early at 30 s (b = 0.04). Sec. 5.6 and 7 now limit the claim to 2–5 min and say the candidate fits were not refitted. While doing this a bug was found in `multistart`: with s fixed and a loss term its nested start ignored the fixed s and could return the drift fit; no earlier sim used that path (they call `fit` directly). Fixed and tested. **Fixed.**

### 284. MINOR: inconsistent numbers (30 s rms 8–16 vs drift 8–12; patches ≲ 0.25 mm at 1.5 THz; 1.3–2.5 vs 1.33–2.38; 1.21–2.33 vs 1.76; 12–99 with λ 2–6; "as they match each other").
**Response.** 8–12 % with a drift; patches ≲ 0.25 mm average out at 0.6 THz but still depart by 1.2–6.9 dB at 1.5 THz; 1.3–2.4; 1.21–1.76; 12–99 tied to λ 1–6 over all fits; "close to how well they match each other". **Fixed.**


## Thirty-sixth round (independent reviewer, whole paper, 2026-11-04; verdict: no fatal flaw; minor revision)

The reviewer checked the `multistart` fix and the t2_max bound and found both correct, with no earlier result depending on the old behaviour.

### 285. MINOR: the data prefer H1's empty pulse, the one that gives the lead, so "within how well that pulse is reproduced" is too strong.
**Response.** Confirmed by Sim 77: H1's pulse fits all four H1–H3 scans better (full fit 6.3/9.9/11.9/13.5 % against 6.7/10.2/13.1/14.7 %; no change with s = 0 likewise). Sec. 5.6 and 7 now say the lead is comparable to the difference between two reports' empty-cuvette pulses and that whether it is real is not decided. **Fixed.**

### 286. MINOR: sample 2's fits use a small negative second pulse at the 0.2 ps bound.
**Response.** Sim 77 refits sample 2 with no second pulse: thinner path (8–16 %) 5.4 % with a drift and 11.9 % without, a pulse 1.25 times stronger and 0.05 ps earlier 8.0 %, no change 21.2 %. Sec. 5.6 quotes these and the negative b of the bounded fits. **Fixed.**

### 287. MINOR: H4's own statements on sample 2 (rise from 0.4 THz, 5 dB near 2 THz, index falling over the band, resonance supposed below the band, "all samples show similar regularities") and its "up to 20 dB" for sample 1 against the traced +15 dB are not reported.
**Response.** Sec. 5.6 quotes them and says sample 2, the one independent test, shows neither the dip nor a second pulse; Sec. 5.8 has a sample-2 row. **Fixed.**

### 288. MINOR: one slice's departure "as large as the mean change itself" holds only for large patches.
**Response.** "a quarter to all of the mean change, depending on patch size" in Sec. 5.7 and 6. **Fixed.**

### 289. MINOR: the lead paragraph mixes Sim 55's floored references with Sim 76's unfloored one; 0.009 vs 0.008.
**Response.** The comparison now names its reference (no floor, a ≤ 1), gives H1's pulse in the same reference (0.020–0.039 ps at 2–5 min, 0.039–0.064 at 30 s), and reads 0.008 ps early to 0.009 ps late (s free and s = 0). **Fixed.**

### 290. MINOR: node states "accept them with a drift in every reference" though they were refitted only at floors 1.0–2.5.
**Response.** "in all eight references with floors of 1.0–2.5 dB/THz (node states were not refitted at floors 0–0.5)". The 2–6 in Table 2 is the best size of the better fit in each reference and the 1–6 in the text the range over all fits, as the text says. **Fixed.**

### 291. MINOR: REPORT row 76 and the closing table here still carry withdrawn claims.
**Response.** Both rewritten (row 76 verdict; the dip row: 21 %, H4 without a frequency, bulk conduction, the surviving causes, λ 1–6). **Fixed.**


## Thirty-seventh round (independent reviewer, whole paper, 2026-11-04; verdict: no fatal flaw; minor revision, one major wording item)

### 292. MAJOR: Sec. 7 says "the published waveforms … decompose into" a lagging, lossier tissue pulse plus a second pulse, but H4's sample 2, the one independent sample, has no second pulse and a stronger, earlier tissue pulse.
**Response.** Agreed. Sec. 7 now says "the H1–H3 waveforms, all from one session and one sample"; Sec. 1, 5.6, 5.8 and 7 say that the independent sample's tissue pulse changes in the opposite direction, with no second pulse. The abstract already reports it ("transmits more and earlier, with no second pulse"). **Fixed.**

### 293. MINOR: "a quarter to all of the mean change" understates the most uneven case (1.1–2.9 times at rms 0.40).
**Response.** "a quarter to three times the mean change, depending on patch size and spread" (Sec. 5.7, 6). **Fixed.**

### 294. MINOR: edge fits are not confined to two references; floor 0.5 (a + b ≤ 1) is the same case; a second borderline without a drift (0.91) is missed; Table 2's best sizes are grid limits.
**Response.** Sec. 5.7 now gives nine passing without a drift (one borderline, 0.91) and three excluded (two firmly, one borderline), counts three verdicts as unresolved, and says the drift fits reach the grid edge in ten of twelve references; Sec. 5.7's thicker-path paragraph, Sec. 6 and 7 follow; Table 2's caption marks mean −30 % and rms 40 % as grid limits. A wider grid can only lower D, so "not firmly excluded" stands. **Fixed.**

### 295. MINOR: the H4 "6 min" scan differs from H3's 5 min scan by 2–3 times the tracing noise, so "adds no time point" is too strong.
**Response.** Sec. 5.8 now says it is a later scan of the same sample or H3's, near saturation (b = 0.24 with H4's empty pulse), consistent with H1's "nothing changes after ~6 min", and not an independent sample. **Fixed.**

### 296. MINOR: sample 2's thinner-path fit needs a drift opposite to thinning, is unstable without it (7.6–16.4 %), and with no H4 threshold 5.4 % against 8.0 % does not rank the readings.
**Response.** Sec. 5.6, 5.7 and 5.8 quote the drift (+0.067 ps) and the spread and say the thinner path and a stronger, earlier pulse fit about equally and are not ranked. **Fixed.**

### 297. MINOR: H4 sample 1's second pulse "within 0.07 ps" uses τ without s.
**Response.** "within 0.02 ps of the empty-cuvette delay (s + τ)": 0.007–0.021 ps over the six fits without a loss term (Sim 73). **Fixed.**

### 298. MINOR: Table 0 omits H4.
**Response.** The row reads "H1–H4 (2024)". **Fixed.**

---

## What matches and what does not (experiment)

| Reported (Hovhannisyan, Makaryan, Martirosyan 2024) | This work |
|---|---|
| Dip near 0.6 THz (9–20 dB) | No node or fibre resonance anywhere in 0.1–3 THz. The shared session's biased waveforms contain a second pulse near the empty-cuvette delay (4, 12, 17, 21 % at 0.5, 2, 4, 5 min) and a tissue pulse with an extra loss rising with frequency and a lag of 0.16–0.28 rad at 0.6 THz, over 36 reference fits (Sim 55); without the change the fits miss at 2–5 min by 17.9–23.0 % of the trace rms, 2.2–3.3 times what fitting noise and tracing allow (Sim 69). The dip then falls at 0.57 THz (5 min, reported 0.56) and 0.63–0.65 THz (4 min, published 0.65): consistent, q fitted. H4's "deep resonance" of sample 1, given without a frequency (spectrum minimum near 0.55 THz, index jump at 0.60 THz), is consistent with the same cancellation (model 0.56–0.59; Sim 73) on a scan that matches H3's 5 min scan (Sim 75); sample 2 shows no dip and no second pulse. Against fitting noise and tracing (Sims 69–72), extra bulk conduction alone up to 60/150 S/m is firmly excluded in all twelve references (1.3–2.4× the threshold); uneven or thicker tissue, gas, warming and node states each pass under some reading, uneven thickness only as a mean over many patches (Sims 72, 74); node states pass only at 1–6× every extreme and as states an aqueous gap cannot hold (Table 2) |
| 2 THz dip (30 s, paper 1) | Not reproduced, not excluded. A real feature of one scan: 15 dB at 1.99 THz, twice any other single-curve dip, its bottom (−71 dB) 2–5 dB above the noise, so the depth is a lower bound; origin unknown. With the traced phase extrapolated the two pulses put their antiphase at 2.10–2.32 THz, and the data are only 3.2 dB above the two-pulse floor (Sims 43, 44). No node state (65, every extreme, six phase bands, b = 0.034–0.044) makes the model's minimum within 0.03 THz of 1.99 THz deeper than −37.5 dB (published ≤ −43.5 dB; Sim 43 (e)) |
| Dip drifts 0.7 → 0.5 THz in minutes | Water-vapour lines cannot drift. At 30 s the traced ratio is a broad −3 to −4 dB trough over 0.56–0.78 THz; by 4 min the dip is at 0.62–0.65 THz as b reaches 0.17 and q 1.09 (Sims 42, 43). Matches the direction; fitted, not predicted |
| High frequencies rise by 17–18 dB (2–2.5 THz, traced) | In the two-pulse fit the second pulse skips the tissue loss; the fit follows the rise to 1.7 THz within 1.5 dB rms (not testable above with the traced window); drying alone ≤ 3 dB (Sim 34) |
| ~9 dB change under 100 V DC | Node (E∥, internode ≥ 200 µm, packing ≤ 0.6): ≤ 0.05 dB for any conductivity, ≤ 0.1/0.6/1.1 dB at 0.6/2/3 THz for any passive ε (ε″ ≥ 0.1, combined ±20 % materials); with internode 50 µm, packing 0.8 and the measured slice's optical thickness (×2.2–3.1) 1.5/10/18 dB at 0.6/2/3 THz. Inside the waveform fits node states (all 80) fit only at 1–6 times their every-extreme size (Sim 69) and as states an aqueous gap cannot hold. At 2 THz the extreme bound (10 dB) approaches the ≥ 14 dB change (nominal margin 7×). Near 3 THz the reported change is a rise (node states ≤ 0.4 dB); channels 10⁻⁹ dB. Fibre conduction: 130–450 S/m for 9 dB at 0.6–2 THz on the slice's optical thickness; inside the waveform fits it is excluded in every reference at any size up to 60 S/m in the tissue water or 150 S/m in the axoplasm (Table 2, Sims 69–70). Heating/drying: smooth, ≤ 3 dB |
| Effect only with E ∥ fibres | Described without data (H1 p. 18, H3, H4), so not testable on the published material. Contradicted by the causes the waveforms leave (uneven or thicker tissue, gas with a drift, warming); warming and conduction in the tissue water, which the waveforms exclude, also act for both polarisations); model dichroism −0.35 dB; a conducting axon is itself E∥-selective (6.6 vs 2.7 dB at 0.6 THz) but needs 130–450 S/m for 9 dB |
| Effect only with DC across the fibres | Described without data. Not explained; heating predicts the opposite (σ∥ ≈ 7σ⊥) |
| Generated high frequencies (nonlinearity) | No ps nonlinearity at TDS fields (10⁻³ of onset); drying raises high-frequency transmission linearly |
| Sample acts as a grating of wires | No diffraction orders below 6 THz; fibres are not periodic |
| Onset ~5 min, recovery ~1 h; no further change after ~6 min (paper 1) | The second pulse grows 4 → 12 → 17 → 21 % (0.5, 2, 4, 5 min); H4's sample 1 at "6 min" matches H3's 5 min scan within 10 % (Sim 75), so it adds no time point; saturation not testable. Conduction alone is excluded by the waveforms; no other tested cause is firmly excluded (Table 2) |
| No significant change of average power (paper 3) | Traced pulse energy 1.03 (biased) vs 1.36 (no voltage): −1.2 dB, consistent with a small change of a heavily attenuated signal |
| Tissue n = 3.6–3.9 | Inconsistent with water-based tissue (model and Berry 2003: 1.95–2.2) |

## Remaining (non-fatal) limitations stated in the paper
- Single fibre plus effective medium; dense-packing corrections are O(1) (item 4).
- No raw experimental data; the analysis rests on traced published waveforms (H1–H3 one session, H4
  one figure per sample; the full fit leaves 10–15 %, several times tracing noise; q depends on the band). The waveforms do not identify the cause; uneven thickness fits only as a mean over many patches that one slice under the beam need not realise (Sim 74), and H4 adds one independent sample. Tests: image the cuvette during bias; record DC current and slice
  mass; measure slice thickness; repeat the E⊥ and DC-direction controls on the same slice.
- Linear passive EM response only; DC-driven structural changes enter only through bulk effects.
- Passive linear model only; band limited to 0.1–3 THz.
