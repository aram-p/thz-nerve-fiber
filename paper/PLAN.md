# Paper plan

**Working title.** Can nodes of Ranvier explain voltage-dependent terahertz absorption in white
matter? A verified full-wave study of a myelinated fibre in water (0.1–3 THz)

**Target journal.** *Journal of Infrared, Millimeter, and Terahertz Waves* (Springer, Scopus-indexed).
It publishes THz biophotonics and computational EM, and a quantitative negative result about a
reported THz bio-effect fits its scope. Fallback: *Journal of Biophotonics* (Wiley), which would want
the tissue-level and experimental-design sections up front.

## Abstract (draft, 250 words max)

Terahertz transmission of spinal-cord white matter has been reported to change by about 9 dB near
0.6 and 2 THz when 100 V DC is applied across the slice, with the THz field along the fibres. Two
explanations exist: resonant absorption at the node of Ranvier, and DC-opened node channels that
turn fibres into nonlinear conducting wires forming a diffraction grating. We test both with a
verified full-wave model of a myelinated fibre in water. An
axisymmetric finite-element model treats a periodic node chain at realistic node (0.5–5 µm),
internode (0.2–1.5 mm), axon and myelin sizes; it reproduces transfer-matrix and Mie solutions to
2 × 10⁻⁴ and converges to 0.5 %. Between 0.1 and 3 THz the node shows no resonance. Its absorption
is quasi-static: a three-constant closed form in node length over myelin thickness, g-ratio and
permittivity contrast matches full wave to 2 % for fibres up to 10 µm in diameter, independent of
internode length. The sheath does not guide terahertz waves into the node. No node conductivity
(10⁻⁴–3 × 10⁶ S/m) or passive node permittivity changes a 100 µm slice by more than 0.1 dB for the
field along the fibres. Open channels add a million times less admittance than the node membrane's
capacitance already passes; wire-like fibres would need an axoplasm conductivity near 500 S/m,
350 times physiological; fibre spacing excludes diffraction below 6 THz. The reported effect does
not come from nodes or fibres; we quantify the bias current's side effects and decisive tests.

## Section outline

1. **Introduction.** Reported voltage-dependent THz absorption at 0.6 and 2 THz in spinal-cord white
   matter (Hovhannisyan, Makaryan, Martirosyan 2024, three reports) and the authors' model: DC
   opens node channels, fibres become nonlinear conducting wires, the slice a diffraction grating.
   Separate proposals that myelin guides THz/IR and that nodes resonate or amplify (Liu 2019, 2020;
   Guo 2022; Song 2026). Questions: does a single myelinated fibre with real sizes and measured
   water dispersion absorb resonantly at the node, and can any node or fibre state produce the
   reported change?
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
   myelin index below water and axon-core V number ≤ 0.4 (leaky, no guiding), water damping (field 1/e length 70–235 µm, below most internodes).
   4.5 Robustness: ±20 % materials, lamellar (uniaxial) myelin, broadband water, constriction,
   membrane, dense packing, thick fibres, dense frequency sweeps.
   4.6 Node states: conductivity 10⁻⁴–3×10⁶ S/m (the gap screens itself) and any passive ε (a
   negative-ε ring resonance is capped by the surrounding water's loss); slice change from three
   observables (absorbed power, near-field absorption, induced dipole): E∥ ≤ 0.1 dB, E⊥ ≤ 1 dB at
   worst-case packing and internode (Fig. 5).
5. **The experimenters' model and the experiment.** 5.1 Membrane admittance at THz vs open
   channels. 5.2 Fibres as wires: full-wave axoplasm conductivity needed for 9 dB (~500 S/m),
   quasi-static mixing vs wire (inductive) regime. 5.3 No grating orders below 6 THz. 5.4 No
   ps-fast nonlinearity at TDS fields; why the filter-order test does not prove nonlinearity
   (cascade reflections, beam path, sample changing between scans). 5.5 Tissue effective medium,
   E ∥ / E ⊥; DC side effects (≈1 mA, 3 MW m⁻³, electrolysis, drying) and their limits (DC along the
   fibres heats more yet reportedly does nothing). 5.6 Baseline 30 dB/0.1 mm at 2 THz vs 11.6 dB for
   water: dynamic range. 5.7 What matches, what does not, and decisive measurements (current and
   mass logging, N₂ purge, reference through a heated or dried slice, raw data on IEEE DataPort)
   (Figs. 5–6).
6. **Discussion and limits.** Single fibre and dilute tissue; quasi-static local field in dense
   packing; no active gain (cannot test Liu 2020); band limited to 0.1–3 THz (Song 2026 is 10–40 THz).
7. **Conclusion.**
   Data and code availability: repository with every script, benchmark and run log.

## Prior work: closest papers and what is new

| Paper | What they did | Band | Water loss / dispersion | Verification | What this work adds |
|---|---|---|---|---|---|
| Liu G et al. 2019, *Adv. Funct. Mater.* 29, 1807862 | SR-FTIR: myelin index ~2× surroundings; sheath as dielectric waveguide | mid-IR | measured contrast, mid-IR only | none reported | In 0.1–3 THz myelin's index (1.7–2.2) is below water's (1.9–3.4) and the axon-core V number is ≤ 0.4, so any guided mode leaks into the water; full-wave shows no field focusing into the node |
| Liu Y et al. 2020, *Sci. China PMA* 63, 274211 | Node as an active amplifier (three-level gain) of guided THz/IR | mid-IR–THz | not stated | none | Passive, energy-conserving calculation: what the node does without assumed gain; their gain hypothesis is outside our scope (stated) |
| Guo L et al. 2022, *Front. Bioeng. Biotechnol.* 10, 1055232 | Analytic 3-layer fibre modes + mode matching at the node, vs CST | 0–30 THz (node 10–30) | constant complex ε, no Debye | vs CST only | 0.1–3 THz with double-Debye water, realistic internodes, absorbed power per node, Mie benchmark, convergence |
| Song M et al. 2026, *Bioengineering* 13, 955 | COMSOL/HFSS: node as impedance transformer, peaks vs node length | 10–40 THz | lossless fluids, ε = 1.44 | none | With water's real ε (4.1−2.1j at 1 THz) the node is quasi-static and non-resonant in 0.1–3 THz; scaling law instead of resonance |
| Hovhannisyan, Makaryan (Martirosyan) 2024, *Armenian J. Phys.* 17, 17; IET Conf. Proc. (IRPhE 2024) | THz-TDS of cord slices, features at 0.6 and 2 THz with DC bias, E ∥ fibres; model: open node channels make fibres nonlinear wires, slice a grating | 0.2–4 THz | n/a (experiment) | n/a | Tests their model and node resonance: open channels are 10⁶× below membrane capacitive admittance; wires need ~500 S/m axoplasm (full wave); no grating orders below 6 THz; no node state gives > 0.1 dB (E∥); DC side effects quantified |

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
