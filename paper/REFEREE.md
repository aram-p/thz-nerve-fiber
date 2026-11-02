# Referee report against our own claim, with responses

**Claim under review.** In 0.1–3 THz a single myelinated fibre in water has no resonant absorption
at the node of Ranvier. The node absorbs as a quasi-static, sub-wavelength gap: P_node =
½ωε0ε″_w|E0|²V_n·F, where F (0.2–2) depends only on node length / myelin thickness, the g-ratio and
the permittivities (F is defined by the absorbed power and computed from a Laplace problem). No node
state, and no fibre state the experimenters proposed (open channels turning fibres into nonlinear
conducting wires forming a grating), changes tissue transmission by an amount close to the reported
~9 dB, so the 0.6/2 THz features reported under DC bias (Hovhannisyan, Makaryan, Martirosyan 2024)
do not come from nodes or fibres.

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
node permittivity changes a 100 µm slab by more than 0.1 dB for E∥ (≤ 1 dB for E⊥) even at
worst-case packing and internode, because a conducting gap screens its own field (item 16); fully
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
shell series for E ⊥; long gap: bulk water, F = 1). Two shape constants per polarisation are fitted
to the Laplace solver over x = 0.2–8, g = 0.6–0.8, 0.1–3 THz, not to full wave and not to any
experiment: three constants per polarisation (the transition length scales as |ε_w/ε_my|^q),
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
(b) Full wave with extra conductivity in the whole axoplasm: 9 dB for E∥ needs +460–520 S/m,
against 1.4 S/m physiological and ~25 S/m for saturated NaCl. Quasi-static mixing matches full wave
to 3 % up to 400 S/m; fibres act as wires (inductive regime) only above ~10⁴ S/m. The change is
broadband and also large for E⊥ (6 dB at 2 THz), which was not observed.
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
against the coated-cylinder value to ≤ 0.6 %). At internode 1000 µm and packing 0.5: E∥ ≤ 0.006 dB,
E⊥ ≤ 0.18 dB. The slice effect scales as packing/internode; the worst realistic case (200 µm,
packing 0.6) is 6× that: E∥ ≤ 0.035 dB, E⊥ ≤ 1.1 dB. Internode independence rechecked at
3 × 10⁴ S/m: ≤ 7 × 10⁻⁴. **Fixed.**

### 17. MAJOR: "'No resonance' covers only conductors; a node with a resonant constituent could resonate."
**Response.** Sim 36 sweeps any passive node permittivity, ε′ −30…30, ε″ 0.01…30, at 0.6 and 2 THz.
A negative-ε ring resonance exists near ε′ ≈ −1.5 (E∥, η ≈ 100), but the lossy water around the gap
caps the slice change at ≤ 0.08 dB (≤ 0.5 dB at internode 200 µm); mesh check 4 %. Negative ε′ with
ε″ ≪ 1 is not possible in an aqueous gap at THz. **Fixed**; claim stated for passive media.

### 18. MAJOR: "The reported baseline exceeds pure water; the 9 dB may sit at the dynamic-range limit."
**Response.** Agreed and stated: 30 dB per 0.1 mm at 2 THz is 2.6× pure water (11.6 dB). Our
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

## What matches and what does not (experiment)

| Reported (Hovhannisyan, Makaryan, Martirosyan 2024) | This work |
|---|---|
| Features at 0.6 and 2 THz | No node or fibre resonance anywhere in 0.1–3 THz; both frequencies sit on water-vapour line groups (0.557; 1.92–2.20 THz) |
| ~9 dB change under 100 V DC | Node: ≤ 0.1 dB (E∥) for any conductivity or passive ε; channels 10⁻⁹ dB. Fibres as wires: needs ~500 S/m axoplasm (350× physiological). Heating/drying: smooth, ≤ 3 dB |
| Effect only with E ∥ fibres | Not explained; model dichroism −0.35 dB; the wire model would change E⊥ too |
| Effect only with DC across the fibres | Not explained; heating predicts the opposite (σ∥ ≈ 7σ⊥) |
| Generated high frequencies (nonlinearity) | No ps nonlinearity at TDS fields (10⁻³ of onset); drying raises high-frequency transmission linearly |
| Sample acts as a grating of wires | No diffraction orders below 6 THz; fibres are not periodic |
| Onset ~5 min, recovery ~1 h | Timing compatible with heating/drying (0.64 K/s adiabatic), but see DC direction |
| Tissue n = 3.6–3.9 | Inconsistent with water-based tissue (model and Berry 2003: 1.95–2.2) |

## Remaining (non-fatal) limitations stated in the paper
- Single fibre plus effective medium; dense-packing corrections are O(1) (item 4).
- No raw experimental data; the instrumental explanations are hypotheses with stated tests
  (purge with N₂; record the DC current and slice mass; reference scan with the same humidity).
- Passive linear model only; band limited to 0.1–3 THz.
