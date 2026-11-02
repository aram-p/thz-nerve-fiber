# Referee report against our own claim, with responses

**Claim under review.** In 0.1–3 THz a single myelinated fibre in water has no resonant absorption
at the node of Ranvier. The node absorbs as a quasi-static, sub-wavelength gap: P_node =
½ωε0ε″_w|E0|²V_n·F, where F (0.2–2) depends only on node length / myelin thickness, the g-ratio and
the permittivities. No physiological node state changes tissue transmission measurably, so the
0.6/2 THz features reported under DC bias (Hovhannisyan & Makaryan 2024) do not come from nodes.

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
quantitative bound from full wave: no node conductivity between 10⁻⁴ and 3 × 10⁶ S/m changes a
100 µm slab by more than 0.14 dB (E∥: ≤ 0.006 dB at 0.6 and 2 THz), because a conducting gap
screens its own field; fully open Na channels (~10⁻⁴ S/m) give 10⁻⁹ dB. **Fixed** by the prior-work
table (`paper/PLAN.md`) and by narrowing the claim to 0.1–3 THz.

### 2. "Your materials are guesses: Looyenga mixtures for myelin and axoplasm, no measured THz data for myelin; molecular bands are missing from a Debye model."
**Response.** (a) Every material ε′ and ε″ was scaled by 0.8 and 1.2, one at a time
(Sim 33 robust): P_node moves by −24 % to +34 %, smoothly at every frequency; no feature
appears, and the closed law predicts each change to ≤ 5 %. (b) The bound in (1) covers any
node-localised loss, resonant or not: adding Δε″ = Δσ/(ωε0) up to ~10⁵ (Δσ = 3 × 10⁶ S/m) to the
gap, far beyond any molecular band (bulk water's ε″ is 2–7), never gives more than 0.14 dB. A node
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
minute or raise it by ~27 K at steady state. That matches the reported ~5 min onset and ~1 h
recovery. Heating and drying change 100 µm slab transmission smoothly by ≤ 3 dB (Sim 34). What is
not hypothetical: no node mechanism reaches the reported magnitude. **Narrowed**; the paper states
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
shows the field is screened and the effect saturates at ≤ 0.14 dB. The bound now rests on full wave
only; the linear helper is documented as valid for small Δσ only. **Fixed.**

---

## What matches and what does not (experiment)

| Reported (Hovhannisyan & Makaryan 2024) | This work |
|---|---|
| Features at 0.6 and 2 THz | No node or fibre resonance anywhere in 0.1–3 THz; both frequencies sit on water-vapour line groups (0.557; 1.92–2.20 THz) |
| ~9 dB change under 100 V DC | Node: ≤ 0.006 dB (E∥) for any conductivity up to 3 × 10⁶ S/m; channels 10⁻⁹ dB. Heating/drying: smooth, ≤ 3 dB |
| Effect only with E ∥ fibres | Not explained; model dichroism −0.35 dB |
| Onset ~5 min, recovery ~1 h | Consistent with DC Joule heating and drying (0.64 K/s adiabatic) |
| Tissue n = 3.6–3.9 | Inconsistent with water-based tissue (model and Berry 2003: 1.95–2.2) |

## Remaining (non-fatal) limitations stated in the paper
- Single fibre plus effective medium; dense-packing corrections are O(1) (item 4).
- No raw experimental data; the instrumental explanations are hypotheses with stated tests
  (purge with N₂; record the DC current and slice mass; reference scan with the same humidity).
- Passive linear model only; band limited to 0.1–3 THz.
