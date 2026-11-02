# Paper plan

**Working title.** Can nodes of Ranvier explain voltage-dependent terahertz absorption in white
matter? A verified full-wave study of a myelinated fibre in water (0.1–3 THz)

**Target journal.** *Journal of Infrared, Millimeter, and Terahertz Waves* (Springer, Scopus-indexed).
It publishes THz biophotonics and computational EM, and a quantitative negative result about a
reported THz bio-effect fits its scope. Fallback: *Journal of Biophotonics* (Wiley), which would want
the tissue-level and experimental-design sections up front.

## Abstract (draft, 250 words max)

Spinal-cord slices under 100 V DC bias reportedly develop a 9–20 dB transmission dip near 0.6 THz
and changes near 2 THz, with the THz field along the fibres. Proposed causes are node-of-Ranvier
resonance and DC-opened channels turning fibres into nonlinear wires forming a grating. We test both
with a verified full-wave model of a myelinated fibre in water: a periodic node chain at realistic
sizes (internode 0.2–1.5 mm), reproducing transfer-matrix and Mie solutions to 2 × 10⁻⁴. In 0.1–3
THz the node does not resonate. Its absorption is quasi-static: a three-constant closed form in node
length over myelin thickness, g-ratio and permittivity contrast gives it to 2 % for fibres up to 10
µm in diameter. For internodes ≥ 200 µm and fibre packing ≤ 0.6, no node conductivity changes a 100
µm slice by more than 0.05 dB (field along fibres), and no passive node permittivity by more than
0.8 dB. Open channels add a million times less admittance than the node membrane already passes, 9
dB from fibre conduction needs 300–1000 S/m extra axoplasm conductivity (physiological: 1.4), and
fibre spacing excludes diffraction. Linear electromagnetic responses of nodes and fibres cannot
produce the dip. Instead, the published biased waveform contains 21 % of the empty-cuvette pulse at
its own delay, which no tissue-path change alone reproduces; with the tissue path also ~10 %
optically thicker, the dip falls at 0.57 THz (reported 0.56) at nearly full depth. Its cause, and
the reported polarisation and bias-direction selectivity, remain open.

## Section outline

1. **Introduction.** Motivation: THz exposure and neuromodulation need node-level dosimetry, and
   nodes have been proposed as IR/THz antennas (Zangari 2018) and fibres as optical or THz channels
   (Kumar 2016; Liu 2019; Xiang 2021). Reported voltage-dependent dip of 9–20 dB near 0.6 THz and changes near 2 THz
   in spinal-cord white matter (Hovhannisyan, Makaryan, Martirosyan 2024, three reports; depths read
   off paper 3 Fig. 1b and paper 1 Fig. 4b) and the authors' model: DC
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
   power, cross-section, node field factor η, share of fibre absorption; for tissue, the fibre's
   induced dipole → effective-permittivity change → exact slab transmission (absorbed power alone
   misses the host-loss coupling, ratio |n|²/n′²; first order fails when Re Δε dominates).
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
   observables (absorbed power, near-field absorption, induced dipole). For E∥, internode ≥ 200 µm,
   packing ≤ 0.6: conductivity ≤ 0.05 dB, passive ε (ε″ ≥ 0.1) ≤ 0.8 dB (negative-ε ring resonance
   only, dense ε′ grid at 0.6/1/2/3 THz); smallest CNS fibres (internode 50–100 µm) up to 1–2 dB; E⊥ ≤ 1.1 dB, 3.5 dB only for a metallic
   node in a 21 µm fibre (Fig. 5).
5. **The experimenters' model and the experiment.** 5.1 Membrane admittance at THz vs open
   channels. 5.2 Fibres as conductors: full-wave axoplasm conductivity needed for 9 dB (470–565 S/m at fibre packing 0.5 and 290–980 S/m over packing 0.7–0.3 (0.6–2 THz),
   exact slab; 0.1 THz 600–1300 S/m), quasi-static mixing vs wire (inductive) regime above ~10⁴ S/m. 5.3 No grating orders below 6 THz. 5.4 No
   ps-fast nonlinearity at TDS fields; why the filter-order test does not prove nonlinearity
   (cascade reflections, beam path, sample changing between scans). 5.5 Tissue effective medium,
   E ∥ / E ⊥; DC side effects (≈1 mA, 3 MW m⁻³, electrolysis, drying) and their limits (DC along the
   fibres heats more yet reportedly does nothing); selectivity stated in the electrode frame
   (E_THz ⊥ E_DC, parallel to electrodes 1–2); post-mortem tissue; ~70 mV across a fibre.
   5.6 Baseline 30 dB/0.1 mm at 2 THz vs 11.6–15.8 dB for water, and a tissue delay that implies
   n = 2.8–3.1 at the stated 100 µm (traced phase, Sim 39) vs 1.95–2.2 for water-based tissue: the
   slice is optically ~2× thicker than stated, or not water-based. If it is ~200 µm thick, every
   per-slice bound here doubles (node ≤ 0.1 dB conductivity / ≤ 1 dB passive ε) and the 9 dB
   conductivity halves (~150–500 S/m); the conclusions stand; dynamic range. 5.7 What matches, what does not, and decisive measurements (current and
   mass logging, N₂ purge, reference through a heated or dried slice, raw data on IEEE DataPort)
   5.8 The published waveforms (traced from the PDF, Fig. 7): the biased pulse contains ~21 % of the
   empty-cuvette pulse at its own delay; no one-path change of the tissue fits (65–94 % residual);
   Fig. 1b alone shows it as a flat −12 dB (biased − empty) over 1.5–3 THz. With a tissue path ~10 %
   optically thicker the sum gives the dip at 0.57 THz, −10.7 dB (residual 20 %, 8× tracing jitter, so
   not exact). Drift: pure bypass dips sit at 0.74–0.76 THz; a thickening tissue path moves them to
   0.56–0.70 THz. Causes unknown; the reported selectivity contradicts them. A tissue gas-column fit
   was tried and withdrawn (needs a 225–240 µm column in a 100 µm cuvette).
   (Figs. 5–7).
6. **Discussion and limits.** Linear passive EM response only: DC-driven structural changes of
   fibres (swelling, electro-osmotic water flow, pH fronts) are tested only through their bulk
   effect (heating, water loss, gas). Single fibre plus effective medium; quasi-static local field
   in dense packing; finite fibre segments in a 70–100 µm slice lower the E∥ response; no active gain (cannot test Liu 2020); band limited to 0.1–3 THz (Song 2026 is 10–40 THz).
7. **Conclusion.**
   Data and code availability: repository with every script, benchmark and run log.

## Prior work: closest papers and what is new

| Paper | What they did | Band | Water loss / dispersion | Verification | What this work adds |
|---|---|---|---|---|---|
| Liu G et al. 2019, *Adv. Funct. Mater.* 29, 1807862 | SR-FTIR: myelin index ~2× surroundings; sheath as dielectric waveguide | mid-IR | measured contrast, mid-IR only | none reported | In 0.1–3 THz myelin's index (1.7–2.2) is below water's (1.9–3.4) and the axon-core V number is ≤ 0.4, so any guided mode leaks into the water; full-wave shows no field focusing into the node |
| Liu Y et al. 2020, *Sci. China PMA* 63, 274211 | Node as an active amplifier (three-level gain) of guided THz/IR | mid-IR–THz | not stated | none | Passive, energy-conserving calculation: what the node does without assumed gain; their gain hypothesis is outside our scope (stated) |
| Guo L et al. 2022, *Front. Bioeng. Biotechnol.* 10, 1055232 | Analytic 3-layer fibre modes + mode matching at the node, vs CST | 0–30 THz (node 10–30) | constant complex ε, no Debye | vs CST only | 0.1–3 THz with double-Debye water, realistic internodes, absorbed power per node, Mie benchmark, convergence |
| Song M et al. 2026, *Bioengineering* 13, 955 | COMSOL/HFSS: node as impedance transformer, peaks vs node length | 10–40 THz | lossless fluids, ε = 1.44 | none | With water's real ε (4.1−2.1j at 1 THz) the node is quasi-static and non-resonant in 0.1–3 THz; scaling law instead of resonance |
| Zangari A, Micheli D, Galeazzi R 2018, *Sci. Rep.* 8, 539 | Node of Ranvier as an array of bio-nanoantennas, FEM | IR (tens–hundreds of THz) | lossless or constant ε | none reported | In 0.1–3 THz, with water's dispersive loss, the node is a quasi-static gap with no antenna resonance; closed-form absorption law |
| Kumar S et al. 2016, *Sci. Rep.* 6, 36508 | Myelinated axons as optical waveguides for biophotons | visible–near-IR | n/a | FDTD, no benchmark | Shows the waveguide picture does not carry over to THz: myelin index below water, axon-core V ≤ 0.4 |
| Xiang Z, Tang C, Chang C, Liu G 2021, *Nano Res.* 14, 590 | THz/IR signal transmission along neurons (ordered surface water) | THz–mid-IR | not dispersive water | none | Passive THz response of a myelinated fibre with measured water, verified solver |
| Hovhannisyan, Makaryan (Martirosyan) 2024, *Armenian J. Phys.* 17, 17; IET Conf. Proc. (IRPhE 2024) | THz-TDS of cord slices, features at 0.6 and 2 THz with DC bias, E ∥ fibres; model: open node channels make fibres nonlinear wires, slice a grating | 0.2–4 THz | n/a (experiment) | n/a | Tests their model and node resonance: open channels are 10⁶× below membrane capacitive admittance; 9 dB needs ~500 S/m axoplasm (full wave); no grating orders below 6 THz; no node state gives > 0.8 dB (E∥); their published waveforms show ~20 % of the beam bypassing the tissue |

## Figures

1. `figures/fig_model.png`: unit cell, ε(f) of water/axoplasm/myelin, n, absorption length, k_w·Ln.
2. `figures/fig_legacy_artefact.png`: legacy pipeline artefact (supplement) and benchmarks.
3. `figures/fig_node_spectrum.png`: η, node cross-section, node share vs f, both polarisations,
   full wave vs quasi-static.
4. `figures/fig_scaling.png`: η vs node length / myelin thickness with the closed law; full wave vs
   Laplace for all 10 geometries; P_node vs internode length.
5. `figures/fig_tissue_voltage.png`: tissue transmission E ∥ / E ⊥; DC mechanisms (channels, any node
   conductivity, heating, drying) vs the reported 9 dB.
6. `figures/fig_twopath.png`: partial-air-path model with the authors' own index (supplement).
7. `figures/fig_bypass.png`: traced paper 3 waveforms with the two-path fit; reported change of the
   spectrum vs the bypass sum and the withdrawn tissue gas-column fit.
