# Paper plan

**Working title.** No resonant terahertz absorption at the node of Ranvier: a benchmarked full-wave
study of a single myelinated fibre in water (0.1–3 THz)

**Target journal.** *Journal of Infrared, Millimeter, and Terahertz Waves* (Springer, Scopus-indexed).
It publishes THz biophotonics and computational EM, and a quantitative negative result about a
reported THz bio-effect fits its scope. Fallback: *Journal of Biophotonics* (Wiley), which would want
the tissue-level and experimental-design sections up front.

## Abstract (draft, 250 words max)

Resonant terahertz absorption at the node of Ranvier has been proposed to explain features near
0.6 and 2 THz reported in spinal-cord white matter under DC bias. We test this with a benchmarked
full-wave model of a single myelinated fibre in water with
measured water dispersion. An axisymmetric finite-element model, expanded in azimuthal orders, treats
a periodic chain of nodes with realistic node (0.5–5 µm), internode (0.2–1.5 mm), axon and myelin
sizes. It reproduces transfer-matrix and Mie-series solutions to 2 × 10⁻⁴ and converges to
0.5 %. Between 0.1 and 3 THz the node shows no resonance for either
polarisation. Its absorption is quasi-static, P = ½ωε0ε″|E0|²V·F, where F (0.2–2) depends only on
node length over myelin thickness, the g-ratio and the permittivity contrast; a three-constant
closed form for F matches full wave to about 2 % for fibres up to 10 µm in diameter, and the
internode length does not enter. Myelin's refractive index lies below water's, so the sheath cannot
guide terahertz waves into the node. These results hold under ±20 % material changes, node
constriction, the axolemma and dense packing. No node conductivity from physiological values to
3 × 10⁶ S/m changes the transmission of a 100 µm slice by more than 0.14 dB, against the reported
9 dB. Single nodes therefore cannot produce the reported features. The DC bias itself deposits
about 2.6 MW m⁻³ in the slice; heating and drying, together with water-vapour lines at 0.557 and
1.9–2.2 THz, are testable alternative explanations.

## Section outline

1. **Introduction.** Reported resonant THz absorption at 0.6 and 2 THz in spinal-cord white matter
   under DC bias (Hovhannisyan & Makaryan 2024); proposals that myelin guides THz/IR and that nodes
   resonate or amplify (Liu 2019, 2020; Guo 2022; Song 2026). Question: can a single myelinated
   fibre, with real sizes and measured THz water dispersion, absorb resonantly at the node?
2. **Model.**
   2.1 Geometry: axon a, myelin b (g = a/b), node gap Ln, internode Lint; ranges from MRG and
   Arancibia-Carcamo 2017. 2.2 Materials: double-Debye water (Kindt & Schmuttenmaer 1996; Liebe 1991
   for temperature), Looyenga myelin (40 % water) and axoplasm (85 % water); cited ranges, ±20 %.
   2.3 Full-wave method: 2D-axisymmetric scattered-field FEM with the exact lossy-host background,
   azimuthal harmonics m = 0..3, periodic node chain, cylindrical PML. 2.4 Observables: node absorbed
   power, cross-section, node field factor η, share of fibre absorption.
3. **Verification.** Slab vs transfer matrix; coated cylinder vs Mie series (TM/TE); axisymmetric
   per-order vs series; null test; energy balance on every solve; mesh convergence; boundary choice.
   (Fig. 2, Table 1.)
4. **Results.**
   4.1 Spectra 0.1–3 THz, both polarisations: smooth, no resonance (Fig. 3).
   4.2 Quasi-static law: P_node = ½ωε0ε″|E0|²V_n F, F from a Laplace problem; agreement with full wave.
   4.3 Scaling with Ln, myelin thickness, axon radius, internode length (Fig. 4).
   4.4 Why no resonance: sub-wavelength node (k_w Ln < 0.15), no Fröhlich pole (all Re ε > 0),
   myelin index below water (anti-guiding), water damping (field 1/e length 70–235 µm, below most internodes).
   4.5 Robustness and node states: ±20 % materials, constriction, membrane, dense packing; node
   conductivity 10⁻⁴–3×10⁶ S/m: the gap screens itself, ≤ 0.14 dB (Fig. 5).
5. **Comparison with experiment.** Tissue effective medium for E ∥ / E ⊥; DC mechanisms (channels,
   heating, water loss); what matches and what does not; instrumental explanations to test
   (water-vapour lines, partial air path) (Figs. 5–6).
6. **Discussion and limits.** Single fibre and dilute tissue; quasi-static local field in dense
   packing; no active gain (cannot test Liu 2020); band limited to 0.1–3 THz (Song 2026 is 10–40 THz).
7. **Conclusion.**
   Data and code availability: repository with every script, benchmark and run log.

## Prior work: closest papers and what is new

| Paper | What they did | Band | Water loss / dispersion | Verification | What this work adds |
|---|---|---|---|---|---|
| Liu G et al. 2019, *Adv. Funct. Mater.* 29, 1807862 | SR-FTIR: myelin index ~2× surroundings; sheath as dielectric waveguide | mid-IR | measured contrast, mid-IR only | none reported | In 0.1–3 THz myelin's index (1.6–1.8) is below water's (1.9–3.4): the sheath anti-guides; full-wave shows no field focusing into the node |
| Liu Y et al. 2020, *Sci. China PMA* 63, 274211 | Node as an active amplifier (three-level gain) of guided THz/IR | mid-IR–THz | not stated | none | Passive, energy-conserving calculation: what the node does without assumed gain; their gain hypothesis is outside our scope (stated) |
| Guo L et al. 2022, *Front. Bioeng. Biotechnol.* 10, 1055232 | Analytic 3-layer fibre modes + mode matching at the node, vs CST | 0–30 THz (node 10–30) | constant complex ε, no Debye | vs CST only | 0.1–3 THz with double-Debye water, realistic internodes, absorbed power per node, Mie benchmark, convergence |
| Song M et al. 2026, *Bioengineering* 13, 955 | COMSOL/HFSS: node as impedance transformer, peaks vs node length | 10–40 THz | lossless fluids, ε = 1.44 | none | With water's real ε (4.1−2.1j at 1 THz) the node is quasi-static and non-resonant in 0.1–3 THz; scaling law instead of resonance |
| Hovhannisyan & Makaryan 2024, *Armenian J. Phys.* 17, 17; IRPhE 2024 | THz-TDS of cord slices, peaks at 0.6 and 2 THz with DC bias, E ∥ fibres | 0.2–4 THz | n/a (experiment) | n/a | Full-wave bound: no node conductivity up to 3×10⁶ S/m changes a 100 µm slab by more than 0.14 dB (E∥ ≤ 0.006 dB; channels 10⁻⁹ dB); DC Joule heating (0.6 K/s) and water-vapour lines as testable explanations |

## Figures

1. `figures/fig_model.png`: unit cell, ε(f) of water/axoplasm/myelin, n, absorption length, k_w·Ln.
2. `figures/fig_legacy_artefact.png`: legacy pipeline artefact (supplement) and benchmarks.
3. `figures/fig_node_spectrum.png`: η, node cross-section, node share vs f, both polarisations,
   full wave vs quasi-static.
4. `figures/fig_scaling.png`: η vs node length / myelin thickness with the closed law; full wave vs
   Laplace for all 10 geometries; P_node vs internode length.
5. `figures/fig_tissue_voltage.png`: tissue transmission E ∥ / E ⊥; DC mechanisms (channels, any node
   conductivity, heating, drying) vs the reported 9 dB.
6. `figures/fig_twopath.png`: partial-air-path model vs reported features (discussion).
