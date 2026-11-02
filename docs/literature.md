# Literature: THz (0.1-3 THz) interaction with a single myelinated fibre and the node of Ranvier

Compiled 2026-11-02 from web sources (publisher pages, PMC/Europe PMC, Crossref metadata, author PDFs).
Each entry says how it was checked: **full text** (read the paper), **abstract** (abstract or metadata only), or **secondary** (number quoted by another paper, not checked against the original). Numbers marked **our calc** are arithmetic done here from the cited inputs, not published values.

Planned claim being tested: in 0.1-3 THz a single myelinated fibre in water/saline shows no resonant absorption at the node; node absorption is non-resonant, set by water loss and quasi-static field ratios; myelin cannot guide THz because its THz index (~1.5-1.8) is below water's (~2-2.6); ion-channel conductance changes are orders of magnitude too small to change THz absorption.

---

## 0. Bottom line

- **No paper found that reports a node-of-Ranvier resonance inside 0.1-3 THz from a full-wave model.** Every resonance or amplification claim at the node sits at 10-40 THz (Song 2026), 10-30 THz (Guo 2022, and that one is a reflection minimum, not a resonance), the mid-IR (Liu 2019/2020, Zhang 2022) or the optical/NIR (Kumar 2016, Zangari 2018).
- **No paper found that already states the planned negative result** (no node resonance at 0.1-3 THz with dispersive water, absorbed power per node, convergence, Mie check). The negative result is not pre-empted.
- **Partly pre-empted:** (a) "myelin index below water at low THz": Zhang et al. 2022 (MD simulation of a DMPC bilayer) already state that the head-group permittivity falls below that of water below ~20 THz, so waves leak out of the sheath there; their spectra start at 1 THz. (b) "ion changes are too small to change THz absorption": Xu et al. 2007 measured that physiological salt (~100 mM) does not measurably change water's 0.15-3.4 THz absorption, refuting Masson et al. 2006.
- **Directly contradicted by our claim:** the Yerevan experiments (Hovhannisyan & Makaryan 2024, local copies in `docs/references/paper_1-3`) report resonant absorption at 0.6 and 2 THz in pig/sheep spinal-cord white matter under 100 V DC across 18 mm. Section 5 lists water-vapour lines and DC-field artefacts (Joule heating, electro-osmotic dehydration) that could produce such features.
- **Works that assume what we dispute:** Guo 2022 and Song 2026 both put myelin's THz permittivity above the surrounding fluids (Song: myelin 2, axoplasm 1.69, extracellular 1.44, all lossless). Song's extracellular value is far below water's measured eps' ~4-5 with eps'' ~2 at 1 THz (Section 2a).

---

## 1. Closest prior work

### 1.1 Liu Y, Wu K, Liu C, Cui G, Chang C, Liu G (2020)
"Amplification of terahertz/infrared field at the nodes of Ranvier for myelinated nerve." *Sci. China Phys. Mech. Astron.* **63**(7), 274211. DOI [10.1007/s11433-019-1530-2](https://doi.org/10.1007/s11433-019-1530-2). Checked: abstract and Crossref metadata; full text paywalled and not read.

- **What they did:** Proposed that myelinated fibres carry mid-IR/THz signals as dielectric waveguides, and that the signal is replenished at each node by a "quasi-quantum", three-energy-level amplifying biomolecular ensemble (a periodic relay).
- **Band:** "mid-infrared to terahertz"; the abstract gives no specific band.
- **Geometry, materials, method:** not visible in the abstract. The waveguide premise comes from Liu 2019 (1.2). The model is a gain model, not a passive EM scattering calculation.
- **Resonance at node?** Claims amplification at the node, not a passive resonance. **Water loss / energy balance / convergence:** cannot be confirmed from the abstract.
- **What our study adds:** a passive, energy-conserving full-wave calculation with dispersive water at 0.1-3 THz. It shows what the node does without an assumed gain medium (non-resonant absorption, power per node). It cannot test their active-gain hypothesis, so that limit should be stated.

### 1.2 Liu G, Chang C, Qiao Z, Wu K, Zhu Z, Cui G, Peng W, Tang Y, Li J, Fan C (2019)
"Myelin sheath as a dielectric waveguide for signal propagation in the mid-infrared to terahertz spectral range." *Adv. Funct. Mater.* **29**(7), 1807862. DOI [10.1002/adfm.201807862](https://doi.org/10.1002/adfm.201807862). Checked: abstract. A conference version is Chang C, Qiao Z, Wu KJ, Huang Y, IEEE ICOPS 2020, p. 103, DOI [10.1109/ICOPS37625.2020.9717635](https://doi.org/10.1109/ICOPS37625.2020.9717635).

- **What they did:** Synchrotron-radiation FTIR microspectroscopy on myelinated nerve (frog sciatic, per Zhang 2022's description of it). They report that "within a certain mid-IR to THz spectral range" myelin's refractive index is ~2x that of the outer medium or the axon. They then computed the radial energy distribution and argue that a ~2 um sheath confines IR energy over mm distances, with energy supplied at the nodes.
- **Band:** mid-IR, where SR-FTIR works well. The abstract does not name a band inside 0.1-3 THz. The "~2x index" claim should not be carried into 0.1-3 THz without data there.
- **Resonance at node?** No passive resonance claimed; amplification is asserted. **Water loss:** the index contrast is measured, but how loss enters the propagation model is not visible in the abstract.
- **What our study adds:** in 0.1-3 THz water's index (2.0-3.4, Section 2a) exceeds any plausible myelin index (~1.6-1.8 including its ~40% water, Section 2b). The waveguide condition fails there, so the sheath is anti-guiding. A Mie benchmark of a coated cylinder (or sphere) in water plus mode solving would show this directly.

### 1.3 Guo L, Xu D, Wang K, Sun Y, Zhang Q, Ning H, Lu C, Wang S, Gong Y (2022)
"Electromagnetic characteristics of in vivo nerve fibers at the terahertz-far-infrared band." *Front. Bioeng. Biotechnol.* **10**, 1055232. DOI [10.3389/fbioe.2022.1055232](https://doi.org/10.3389/fbioe.2022.1055232). Checked: full text (numbers in figures not legible).

- **What they did:** Analytical dispersion relation of a 3-layer cylindrical fibre (axoplasm / myelin / extracellular fluid; Bessel/Neumann expansion; guided and radiation modes). Added a mode-matching method ("RNMMA") for the node, with S-parameters cross-checked against CST.
- **Band:** dispersion curves 0-30 THz; modes, energy and attenuation at 30 THz; node S-parameters 10-30 THz.
- **Geometry:** node gap 2 um, myelin 2 um at the node model, 100 um axon length; layer radii from Kwon 2017 (values not legible). Internode length not stated.
- **Materials:** myelin indices taken from Liu 2019 (high index, "low loss"). Fluids complex, with an imaginary part plus an ohmic conductivity term. No explicit Debye dispersion.
- **Resonance at node?** No. S11 has a minimum at 12.5 THz and S21 stays > -5 dB over 10-30 THz ("low-loss forward transmission").
- **Water loss:** yes, qualitatively. **Energy balance / mesh convergence:** not reported.
- **What our study adds:** the 0.1-3 THz band they skip; a double-Debye water model instead of a single complex constant; realistic internode lengths (0.2-1.5 mm); absorbed power per node rather than S-parameters of a guided mode; and a test of their premise that myelin is the high-index layer at these frequencies.

### 1.4 Song M, Xu Z, Xu J (2026)
"Simulation of neural signals' spectrum and transmission from the perspective of impedance-matching effect of the Ranvier node." *Bioengineering* **13**(8), 955. DOI [10.3390/bioengineering13080955](https://doi.org/10.3390/bioengineering13080955). Checked: full text (PMC13509889). **The newest and closest FEM work.**

- **What they did:** COMSOL 6.2 transient EM (vector potential) plus HFSS S21 sweeps of a node + 1 mm myelinated segment, driven by 1 pA current ports on the node. They treat the node as a coaxial-line impedance transformer (quarter-wave argument) and fit an equivalent RLC.
- **Band:** 100 kHz, 10, 20, 27 and 40 THz; S21 over 1-100 THz. **Nothing in 0.1-3 THz.**
- **Geometry:** axon 3 um diameter, myelin 1 um, node length swept 1.0-2.5 um (to 4.5 um at 10 THz), nodal membrane 10 nm, 97.5 um extracellular layer.
- **Materials:** at THz, constant real permittivities: intracellular 1.69, myelin 2, extracellular 1.44, with no loss given. At 100 kHz: intra 136 / 0.487 S/m, myelin 10.3 / 0, extra 136 / 1.454 S/m.
- **Resonance at node?** Yes, in their sense: amplitude peaks vs node length (10 THz: 3.5 um; 20 THz: 2 um; 40 THz: 1.75 um), and the RLC estimate gives ~26.6-28.8 THz. They say water absorbs strongly at THz but argue it does not change the impedance discontinuity. The full-wave model has no THz loss term.
- **Energy balance / convergence:** not reported (physics-controlled mesh only).
- **What our study adds:** this is the paper our results most directly contradict. Lossless fluids with eps = 1.44 are inconsistent with water (eps ~ 4.1 - 2.1j at 1 THz). With dispersive water the extracellular medium becomes the high-index, lossy layer, and any impedance-matching peak would be damped. Show absorbed power vs node length at 0.1-3 THz with converged meshes and a power balance.

### 1.5 Zangari A, Micheli D, Galeazzi R, Tozzi A (2018)
"Node of Ranvier as an array of bio-nanoantennas for infrared communication in nerve tissue." *Sci. Rep.* **8**, 539. DOI [10.1038/s41598-017-18866-x](https://doi.org/10.1038/s41598-017-18866-x). Checked: full text.

- **What they did:** COMSOL FEM model of Na-channel currents at the node as a phased array of dipole nano-antennas (40 dipoles, ~13 nm long, 27-50 nm spacing, densities 1400/um^2 for rat and 400-700/um^2 for rabbit). They also ran beam-envelope mode analysis of the myelinated segment and a full-structure run.
- **Band:** 300-2500 nm (120-1000 THz). Not THz.
- **Geometry:** array model axon 0.5 um; full model core 1.026 um, myelin 0.234 um, fibre 1.674 um; paranode 0.5 um; node length quoted 0.7-1.4 um (mean 1.08 um).
- **Materials:** refractive indices from Kumar 2016; waveguide conductivity 0 (lossless). Water absorption neglected in the simulation, then used afterwards to restrict the useful band to ~700-1600 nm. Mesh max element lambda/10. No convergence study.
- **Resonance at node?** No; the claim is directivity of an antenna array.
- **What our study adds:** shows that at 0.1-3 THz the node is ~1/100-1/1000 of a wavelength in water. It is a quasi-static scatterer, not an antenna, so field enhancement is set by permittivity ratios. Absorption is set by water's eps''.

### 1.6 Kumar S, Boone K, Tuszynski J, Barclay P, Simon C (2016)
"Possible existence of optical communication channels in the brain." *Sci. Rep.* **6**, 36508. DOI [10.1038/srep36508](https://doi.org/10.1038/srep36508). Checked: full text.

- **What they did:** Lumerical FDTD and eigenmode modelling of myelinated axons as optical waveguides, including paranodes, nodes, bends, cross-section changes and crosstalk.
- **Band:** 0.4-1.3 um (230-750 THz).
- **Geometry:** axon radius 0.6-3 um, g-ratio 0.6, node ~1 um bare axon, paranode via "p-ratio" ~5, internode 100-150x axon diameter.
- **Materials:** n = 1.44 (myelin), 1.38 (axon), 1.34 (interstitial fluid); absorption argued to be negligible (white matter 0.07-0.3 mm^-1 at 0.4-1.1 um).
- **Resonance at node?** No; the node and paranode are loss sites with fairly high transmission. **Convergence:** not reported in the parts read.
- **What our study adds:** the optical index ordering (myelin > water) does not carry over to THz. Water's index rises from 1.33 (optical) to 2-3.4 (0.1-3 THz) while lipid stays ~1.5.

### 1.7 Xiang Z, Tang C, Chang C, Liu G (2021) and Xiang Z, Tang C, Yan L, Chang C, Liu G (2020)
- "A new viewpoint and model of neural signal generation and transmission: Signal transmission on unmyelinated neurons." *Nano Research* **14**(3), 590-600 (2021). DOI [10.1007/s12274-020-3016-1](https://doi.org/10.1007/s12274-020-3016-1). Checked: abstract.
- "... Signal transmission on myelinated neuron." *Chin. Phys. B* **29**(10), 108701 (2020). DOI [10.1088/1674-1056/abaee9](https://doi.org/10.1088/1674-1056/abaee9). Checked: abstract.
- **What they did:** A field/circuit model of neural signal generation in which, on short time scales, the signal shows THz and IR electromagnetic oscillations under a slowly varying envelope. In the myelinated case the transmembrane ion current at each node is an energy supplement, and the sheath length has an upper limit l_T (mm scale) beyond which signals fail. Temperature extremes impede conduction.
- **Resonance at node?** The node is an energy source in their model, not a resonator. No full-wave scattering, no water dispersion visible in the abstracts.
- **What our study adds:** an external-illumination scattering and absorption calculation with measured material data. It is a different question from their intrinsic-source model, but it bounds how strongly any node-localised THz field couples to the surrounding water.

### 1.8 Zhang Z, Li Y, Xiang Z, Huang Y, Wang R, Chang C (2022)
"Dielectric dispersion characteristics of the phospholipid bilayer with subnanometer resolution from terahertz to mid-infrared." *Front. Bioeng. Biotechnol.* **10**, 984880. DOI [10.3389/fbioe.2022.984880](https://doi.org/10.3389/fbioe.2022.984880). Checked: full text (values only in figures).

- **What they did:** Atomistic MD (GROMACS, CHARMM36, TIP3P; 128 DMPC + 6968 water; 310 K; 1 ns) giving spatially resolved (0.2 nm slices) dielectric spectra over 1-100 THz via flux-flux correlation plus Kramers-Kronig.
- **Main claim:** head-group eps exceeds tail eps over most of 1-100 THz (head eps' up to ~7 at some slices). **Below ~20 THz (and above ~85 THz) the head-region permittivity is below that of water, so waves would leak from the sheath.** 1-55 THz is "unsuitable" (high eps''); 55-85 THz is proposed as the propagation window.
- **Relevance:** this partly pre-empts our "myelin index below water" claim, but only from a 1 ns MD of one lipid using the TIP3P water model (known to misrepresent water's THz response). It has no data below 1 THz. Our study can cite it as consistent and add measured water dispersion plus a full fibre geometry.

### 1.9 Other relevant items (not "closest" but useful)
- **Masson J-B, Sauviat M-P, Martin J-L, Gallot G (2006).** "Ionic contrast terahertz near-field imaging of axonal water fluxes." *PNAS* **103**(13), 4808-4812. DOI [10.1073/pnas.0510945103](https://doi.org/10.1073/pnas.0510945103). Checked: abstract. THz near-field (200 um aperture) imaging of single axons in saline, with 3D simulation of the axon-aperture system. Claimed THz sensitivity to ion concentration. Later corrected for a unit error in Fig. 1, and the ionic-contrast premise was disputed by Xu 2007 (below). This is the only THz-on-single-axon experiment found.
- **Xu J, Plaxco KW, Allen SJ, Bjarnason JE, Brown ER (2007).** "0.15-3.72 THz absorption of aqueous salts and saline solutions." *Appl. Phys. Lett.* **90**, 031908. DOI [10.1063/1.2430781](https://doi.org/10.1063/1.2430781). Checked: full text. NaCl molar extinction is ~2 orders of magnitude below Masson's. "Physiological concentration of a salt (~100 mM) will not significantly increase the terahertz (0.15-3.42 THz) absorption of the solution above that of pure water." Below 0.4 THz, 5 M NaCl is indistinguishable from dilute buffer. **Direct support for our ion-channel claim.**
- **Berry E et al. (2003)** (Section 2c): broadband THz absorption of freshly excised human peripheral **nerve** is 246 +/- 27 cm^-1 vs water 225 +/- 21 cm^-1 (0.5-1.5 THz). Nerve tissue absorbs like water. Supports "absorption is set by water".
- **Omidi M, Zibaii MI, Granpayeh N (2022).** "Simulation of nerve fiber based on anti-resonant reflecting optical waveguide." *Sci. Rep.* **12**, 19356. DOI [10.1038/s41598-022-23580-4](https://doi.org/10.1038/s41598-022-23580-4). COMSOL optical-band study (n = 1.38/1.44/1.34). Not THz.
- **Romanenko S, Siegel PH, Wagenaar DA, Pikov V (2014).** "Effects of millimeter wave irradiation and equivalent thermal heating on the activity of individual neurons in the leech ganglion." *J. Neurophysiol.* **112**(10), 2423-2431. DOI [10.1152/jn.00357.2014](https://doi.org/10.1152/jn.00357.2014). 60 GHz effects on neurons are largely thermal; useful framing for "THz effect = water heating".
- **Hovhannisyan T, Makaryan A (2024).** "Investigation of nerve fiber network in THz spectrum range." *Armenian J. Phys.* **17**(1), 17-20. DOI 10.54503/18291171-2024.17.1-17. Plus an IRPhE 2024 conference paper (Hovhannisyan, Makaryan, Martirosyan) and the dataset doi:10.21227/h17h-tp94 (IEEE DataPort, 2026). Local copies: `docs/references/paper_1.txt`, `paper_3.txt`. Menlo TDS, 0.2-4 THz; 70-100 um white-matter slices of pig/sheep thoracic cord on 1.2 mm PTFE; tungsten electrodes 18 mm apart; 100 V DC. Report resonant absorption at 0.6 and 2 THz only when the fibres are parallel to the THz E-field and the DC field is perpendicular to the fibres. The effect grows over ~5 min and recovers in ~1 h. The DataPort page gives a threshold of > 5 V/mm. **This is the experimental result our negative computation addresses.**

### 1.10 Per-paper checklist

| Paper | Band | Node resonance reported? | Water loss in model | Dispersive water | Energy balance | Mesh convergence | Mie/analytic benchmark |
|---|---|---|---|---|---|---|---|
| Liu 2020 SCPMA | mid-IR-THz (unspecified) | "amplification" (active gain) | unknown (abstract) | unknown | no (gain model) | unknown | no |
| Liu 2019 AFM | mid-IR (SR-FTIR) | no; relay amplification | measured contrast | n/a | no | unknown | no |
| Guo 2022 | 0-30 THz (node 10-30 THz) | no (S11 min 12.5 THz) | yes, constant complex | no | no | no | vs CST only |
| Song 2026 | 10-40 THz | yes (node-length peaks; RLC ~27 THz) | no (lossless at THz) | no | no | no | no |
| Zangari 2018 | 120-1000 THz | no (antenna directivity) | no (sigma = 0) | no | no | lambda/10 rule only | no |
| Kumar 2016 | 230-750 THz | no | argued negligible | no | not stated | not stated | no |
| Xiang 2020/2021 | THz-IR (intrinsic) | node = source | not visible | not visible | no | no | no |
| Zhang 2022 MD | 1-100 THz | n/a (material) | n/a | MD water (TIP3P) | n/a | n/a | n/a |
| **This study** | **0.1-3 THz** | **tests for it** | **yes** | **double-Debye** | **yes** | **yes** | **Mie** |

---

## 2. THz material data

### 2a. Liquid water, double-Debye

Model: eps(w) = eps_inf + (eps_s - eps_1)/(1 + i w tau_1) + (eps_1 - eps_inf)/(1 + i w tau_2).

| Source | T | eps_s | eps_1 (2nd step) | eps_inf | tau_1 (ps) | tau_2 (ps) | Valid range | Check |
|---|---|---|---|---|---|---|---|---|
| Kindt & Schmuttenmaer 1996, *J. Phys. Chem.* **100**(24), 10373-10379, DOI [10.1021/jp960141g](https://doi.org/10.1021/jp960141g) | 25 C | 78.36 | 4.93 | 3.48 | 8.24 | 0.180 | ~0.1-2 THz | secondary (Table 2 of Zhou et al., PCCP 2021, 23, 4426, DOI 10.1039/d0cp06707e) |
| Jepsen, Moller, Merbold 2007, *Opt. Express* **15**, 14717, DOI 10.1364/OE.15.014717 | ~room | 78.36 | 5.16 | 3.49 | 7.89 | 0.181 | THz | secondary (same table) |
| Pickwell et al. 2004, *Appl. Phys. Lett.* **84**(12), 2190-2192, DOI [10.1063/1.1688448](https://doi.org/10.1063/1.1688448) | ~20 C | 78.8 | 6.6 | 4.1 | 10.06 in table (10.6 often quoted; verify) | 0.180 | THz (tissue model) | secondary |
| Liebe, Hufford, Manabe 1991, *Int. J. Infrared Millim. Waves* **12**(7), 659-675, DOI 10.1007/BF01008897 | formula | 77.66 - 103.3 th | 0.0671 eps_s | 3.52 + 7.52 th | from gamma_1 | gamma_2 = 39.8 gamma_1 | f < 1 THz, about -4 to 40 C | **full text** |
| -> Liebe evaluated (our calc) | 20 C | 80.07 | 5.37 | 3.34 | 9.39 | 0.236 | | |
| -> Liebe evaluated (our calc) | 25 C | 78.30 | 5.25 | 3.47 | 8.24 | 0.207 | | |
| -> Liebe evaluated (our calc) | 37 C | 74.28 | 4.98 | 3.77 | 6.28 | 0.158 | | |

- Liebe details (full text): th = 1 - 300/(273.15 + T[C]); gamma_1 = 20.20 + 146.4 th + 316 th^2 GHz; tau = 1/(2 pi gamma). Broadband extension to 30 THz adds two Lorentzians: f_1 = 5.11 THz, gamma_1 = 4.46 THz, A_1 = 25.03 THz^2; f_2 = 18.2 THz, gamma_2 = 15.4 THz, A_2 = 282.4 THz^2 (room temperature only).
- **Ronne C, Thrane L, Astrand P-O, Wallqvist A, Mikkelsen KV, Keiding SR (1997)**, *J. Chem. Phys.* **107**(14), 5319-5331, DOI [10.1063/1.474242](https://doi.org/10.1063/1.474242). THz reflection spectroscopy 0.1-2.0 THz, 271-367 K; Debye (ps) plus fast (fs) term vs temperature. Fit values not retrieved (paywalled). Use for the temperature dependence of tau_D.
- **Yada H, Nagai M, Tanaka K (2008)**, *Chem. Phys. Lett.* **464**(4-6), 166-170, DOI [10.1016/j.cplett.2008.09.015](https://doi.org/10.1016/j.cplett.2008.09.015). TDS-ATR 0.2-3.5 THz, 270-362 K. Four components (slow Debye, fast relaxation, 5 THz stretch, libration). The fast component has ~2% of the slow strength, a sub-ps time, and almost no T dependence. Numeric values not retrieved.
- **Ellison WJ (2007)**, *J. Phys. Chem. Ref. Data* **36**(1), 1-18, DOI [10.1063/1.2360986](https://doi.org/10.1063/1.2360986). 0-25 THz, 0-100 C: three Debye relaxations plus two resonances. At 25 C the relaxation frequencies are 18.56 GHz, 167.83 GHz and 1.944 THz, and the resonances 4.03 and 14.48 THz (secondary, from abstract mirror). **Recommended reference model for 0.1-3 THz at 37 C** (free PDF at srd.nist.gov/jpcrdreprint/1.555632.pdf).

Water optical constants (our calc; intensity absorption alpha = 2 w k / c):

| f (THz) | Kindt 25 C: eps' / eps'' | n | alpha (cm^-1) | Liebe 37 C: n | alpha (cm^-1) |
|---|---|---|---|---|---|
| 0.1 | 7.55 / 13.83 | 3.41 | 85 | 3.75 | 93 |
| 0.5 | 4.69 / 3.45 | 2.29 | 158 | 2.37 | 176 |
| 1.0 | 4.14 / 2.14 | 2.10 | 213 | 2.17 | 228 |
| 2.0 | 3.72 / 1.25 | 1.96 | 267 | (beyond Liebe validity) | |
| 3.0 | 3.60 / 0.87 | 1.91 | 285 | | |

Check: Berry 2003 measured broadband deionized water at 2.04 +/- 0.07 (n) and 225 +/- 21 cm^-1 (0.5-1.5 THz), consistent with the table. Effective THz conductivity of water at 1 THz: sigma_eff = w eps0 eps'' ~ 119 S/m (our calc).

**Saline:** Xu et al. 2007 (above). ~100 mM NaCl does not measurably change 0.15-3.4 THz absorption, so saline ~ water in this band. The DC/low-frequency conductivity of saline/CSF is ~1.5-1.8 S/m: Baumann SB et al. 1997, "The electrical conductivity of human cerebrospinal fluid at body temperature", *IEEE TBME* **44**(3), 220-223, DOI [10.1109/10.554770](https://doi.org/10.1109/10.554770) (1.79 S/m at 37 C; value recalled, not re-read). At THz the ionic Drude term is negligible next to the water relaxation.

### 2b. Myelin, lipid, phospholipid bilayer

- **No direct measurement of isolated myelin's complex permittivity in 0.1-3 THz was found.** The ~2x index of Liu 2019 is a mid-IR SR-FTIR result.
- **Lipid / adipose (bulk proxy):** Berry et al. 2003 (below): human adipose broadband n = 1.50 +/- 0.47, mu = 89 +/- 23 cm^-1 (0.5-1.5 THz, fresh, ~18-20 C; this adipose still contains some water). Microwave fat value n ~ 1.87 at 23.6 GHz (quoted in Berry from Gabriel's review). Truong BCQ, Tuan HD, Fitzgerald AJ, Wallace VP, Nguyen HT (2015), "A dielectric model of human breast tissue in terahertz regime", *IEEE TBME* **62**(2), 699-707, DOI [10.1109/TBME.2014.2364025](https://doi.org/10.1109/TBME.2014.2364025), gives non-Debye + Debye fits for breast/adipose (parameters not retrieved).
- **Phospholipid THz-TDS:** Lin Y, Wu X, Wang K, Shang S, Gong Y, Zhao H, Wu D, Zhang P, Lu X (2023), "Spectral characteristics and functional responses of phospholipid bilayers in the terahertz band", *Int. J. Mol. Sci.* **24**(8), 7111, DOI [10.3390/ijms24087111](https://doi.org/10.3390/ijms24087111). Checked: abstract. DSPE, DPPC, SPH and lecithin bilayers, 0.5-4 THz. The only distinct feature is a ~3.0-3.1 THz tail collective mode (DPPC/DSPE). No sharp lipid resonance inside 0.5-2.5 THz; n/alpha values not retrieved.
- **MD bilayer:** Zhang 2022 (1.8): head eps' up to ~7 in places, tail lower, 1-100 THz; head eps below water below ~20 THz.
- **Composition:** myelin in situ is ~40% water; dry mass is 70-85% lipid and 15-30% protein. Source: Quarles RH, Macklin WB, Morell P, "Myelin formation, structure and biochemistry", in Siegel GJ et al. (eds.), *Basic Neurochemistry*, 6th ed., Lippincott-Raven, 1999 (NCBI Bookshelf NBK28221); updated chapter Rasband MN, Macklin WB, "Myelin structure and biochemistry", *Basic Neurochemistry* 8th ed., 2012, pp. 180-199, DOI 10.1016/B978-0-12-374947-5.00010-9. Checked: secondary/snippet.
- **Estimated myelin THz index (our calc):** Looyenga mixing of water (eps = 4.14 - 2.14j at 1 THz) with dry lipid/protein (eps ~ 2.3 - 0.05j, n ~ 1.52):

| water fraction | eps (1 THz) | n | alpha (cm^-1) |
|---|---|---|---|
| 0.2 | 2.65 - 0.35j | 1.63 | 45 |
| 0.3 | 2.82 - 0.52j | 1.69 | 64 |
| 0.4 | 3.01 - 0.70j | 1.75 | 85 |
| 0.5 | 3.19 - 0.90j | 1.80 | 105 |

  This supports the planned 1.5-1.8 range (1.75 at 40% water) and is always below water's 2.1. Myelin is lamellar, so its permittivity is uniaxial (Wiener bounds differ for E parallel vs perpendicular to the lamellae). Worth a sensitivity run.

### 2c. White matter, nerve, spinal cord, brain tissue

- **Berry E, Fitzgerald AJ, Zinov'ev NN, Walker GC, Homer-Vanniasinkam S, Sudworth CD, Miles RE, Chamberlain JM, Smith MA (2003).** "Optical properties of tissue measured using terahertz pulsed imaging." *Proc. SPIE* **5030**, 459-470. DOI [10.1117/12.479993](https://doi.org/10.1117/12.479993). Checked: **full text** (author version). Transmission TPI, 0.5-1.5 THz, fresh human tissue within 48 h, TPX cells 50-200 um, room 18-20 C, **not N2-purged**. Broadband n and mu: water 2.04 +/- 0.07 / 225 +/- 21; **nerve 1.95 +/- 0.46 / 246 +/- 27 cm^-1**; adipose 1.50 +/- 0.47 / 89 +/- 23; muscle 2.00 +/- 0.35 / 164 +/- 17; skin 1.73 +/- 0.29 / 121 +/- 18; vein 1.58 / 110; artery 1.86 / 151. The authors note a 1.2 THz feature from **unpurged water vapour** in the reference path (relevant to Section 4).
- **Oh SJ, Kim S-H, Ji YB, Jeong K, Park Y, Yang J, Park DW, Noh SK, Kang S-G, Huh Y-M, Son J-H, Suh J-S (2014).** "Study of freshly excised brain tissues using terahertz imaging." *Biomed. Opt. Express* **5**(8), 2837-2842. DOI [10.1364/BOE.5.002837](https://doi.org/10.1364/BOE.5.002837). Checked: full text. Rat, reflection imaging 0.3-1.3 THz. Grey and white matter are distinguishable because of myelin. "Absorbance and reflectivity of [myelin] are significantly lower than those of water." No n/alpha numbers.
- **Yamaguchi S, Fukushi Y, Kubota O, Itsuji T, Ouchi T, Yamamoto S (2016).** "Brain tumor imaging of rat fresh tissue using terahertz spectroscopy." *Sci. Rep.* **6**, 30124. DOI [10.1038/srep30124](https://doi.org/10.1038/srep30124). Checked: full text. Reflection TDS 0.8-1.5 THz; sections in quartz, dry-air purge, measured within ~10 min. Tumour > normal n and alpha (~5% more water per the companion *Phys. Med. Biol.* 2016 paper). Normal values only in figures.
- **Gavdush AA et al. (2021).** "Terahertz dielectric spectroscopy of human brain gliomas and intact tissues ex vivo: double-Debye and double-overdamped-oscillator models of dielectric response." *Biomed. Opt. Express* **12**(1), 69. DOI [10.1364/BOE.411025](https://doi.org/10.1364/BOE.411025). Checked: abstract. **Gives fitted double-Debye parameters for intact human brain tissue** (table not retrieved; get it from the full text). Also uses the sum rule to estimate water content.
- **Gavdush AA et al. (2019).** *J. Biomed. Opt.* **24**(2), 027001. DOI [10.1117/1.JBO.24.2.027001](https://doi.org/10.1117/1.JBO.24.2.027001). Gelatin embedding preserves the fresh-tissue THz response. Notes grey/white contrast attributed to myelin.
- **Kucheryavenko AS et al. (2021).** "Terahertz dielectric spectroscopy and solid immersion microscopy of ex vivo glioma model 101.8: brain tissue heterogeneity." *Biomed. Opt. Express* **12**(8), 5272-5289. DOI [10.1364/BOE.432758](https://doi.org/10.1364/BOE.432758). Rat; 0.2-1.1 THz; "white matter features slightly higher refractive index than the gray one"; PE film used to prevent (de)hydration.
- **Zou Y et al. (2017).** "Terahertz spectroscopic diagnosis of myelin deficit brain in mice and rhesus monkey with chemometric techniques." *Sci. Rep.* **7**, 5176. DOI [10.1038/s41598-017-05554-z](https://doi.org/10.1038/s41598-017-05554-z). Paraffin-embedded (dehydrated) tissue, so not usable for wet n/alpha. Myelin-deficient vs normal differ at 0.4-1.0 THz.
- **Smolyanskaya OA et al. (2018).** "Terahertz biophotonics as a tool for studies of dielectric and spectral properties of biological tissues and liquids." *Prog. Quantum Electron.* **62**, 1-77. DOI [10.1016/j.pquantelec.2018.10.001](https://doi.org/10.1016/j.pquantelec.2018.10.001). Review with tissue tables (not read here).
- **Spinal cord:** no peer-reviewed THz n/alpha for spinal cord was found except the Hovhannisyan/Makaryan work (1.9), whose raw data are on IEEE DataPort (subscription).

### 2d. Cytoplasm / axoplasm / protein solutions

- No axoplasm-specific THz data found. Axoplasm is ~80-85% water plus protein, so model it as water with a small protein-volume reduction of eps''.
- **Shiraga K, Ogawa Y, Suzuki T, Kondo N, Irisawa A, Imamura M (2014).** "Characterization of dielectric responses of human cancer cells in the terahertz region." *J. Infrared Millim. Terahertz Waves* **35**(5), 493-502. DOI [10.1007/s10762-014-0067-y](https://doi.org/10.1007/s10762-014-0067-y). TDS-ATR of living cell monolayers; intracellular water dynamics, with differences best seen below 1 THz (values not retrieved).
- **Xu J, Plaxco KW, Allen SJ (2006)**, *Protein Sci.* **15**, 1175 and *J. Chem. Phys.* **124**, 036101: protein and buffer THz absorption. Dilute buffers are indistinguishable from water over 0.3-3.72 THz (cited in Xu 2007).
- Low-frequency axoplasm resistivity: 70 Ohm cm (1.43 S/m) in the MRG and Arancibia-Carcamo models (Section 3).

### 2e. Myelin water content
- ~40% water in situ (Basic Neurochemistry; see 2b). Lipid 70-85% of dry weight; human frontal-lobe myelin lipid 78-81% of dry weight (O'Brien & Sampson, *J. Lipid Res.* 1965, secondary).
- Note: the MRI "myelin water fraction" (~10-15% of white-matter water) is a different quantity: water trapped between lamellae as a share of all tissue water.

### Material summary for the model (0.1-3 THz)

| Domain | Recommended model | n at 1 THz | alpha at 1 THz (cm^-1) | Source |
|---|---|---|---|---|
| Water / ECS / saline, 37 C | Ellison 2007, or Liebe 1991 (<1 THz), or Kindt 1996 (25 C) | 2.1-2.2 | 213-228 | 2a |
| Axoplasm | water x (0.8-0.9 water fraction) + protein eps ~2.5 | ~2.0 | ~180-200 (estimate) | 2d |
| Myelin | Looyenga water(40%) + lipid; sweep 20-50% | 1.6-1.8 | 45-105 | 2b |
| Nodal membrane (5-10 nm) | eps ~ 2-3 lipid + conductance term (negligible) | ~1.5 | n/a | 3 |
| Bulk nerve tissue (check) | measured | 1.95 +/- 0.46 | 246 +/- 27 (0.5-1.5 THz avg) | Berry 2003 |

---

## 3. Anatomy and nodal electrical parameters

| Quantity | Value | Source (check) |
|---|---|---|
| Node length | 1.08 +/- 0.02 um (EM, range 0.7-1.4 um); 1.02 um confocal (0.5-2.2 um); cortex 1.50 um (0.43-3.72 um) | Arancibia-Carcamo IL, Ford MC, Cossell L, Ishida K, Tohyama K, Attwell D (2017), "Node of Ranvier length as a potential regulator of myelinated axon conduction speed", *eLife* **6**, e23329, DOI [10.7554/eLife.23329](https://doi.org/10.7554/eLife.23329) (full text) |
| Node length in standard PNS model | 1 um | McIntyre CC, Richardson AG, Grill WM (2002), "Modeling the excitability of mammalian nerve fibers: influence of afterpotentials on the recovery cycle", *J. Neurophysiol.* **87**(2), 995-1006, DOI [10.1152/jn.00353.2001](https://doi.org/10.1152/jn.00353.2001) ("MRG"; ModelDB 3810 code read) |
| Node diameter vs axon | CNS: node 86-88% of internodal axon; MRG: nodeD 1.9-5.5 um vs axonD 3.4-12.7 um (node ~0.45-0.55 x axon) | Arancibia-Carcamo 2017; MRG code |
| Periaxonal space | 15 nm internodal (CNS model); MRG 2 nm at node/MYSA, 4 nm FLUT/STIN | Arancibia-Carcamo 2017; MRG code |
| Paranode | MRG MYSA 3 um, FLUT 35-60 um; Zangari 0.5 um; Kumar p-ratio ~5 | MRG; Zangari 2018; Kumar 2016 |
| Internode length | MRG deltax 500, 750, 1000, 1150, 1250, 1350, 1400, 1450, 1500 um for fibre D = 5.7, 7.3, 8.7, 10, 11.5, 12.8, 14, 15, 16 um (~88-115 x D); 100-150 x axon diameter; CNS cortex 27-154 um (mean 83 um), optic nerve 139 um | MRG code; Kumar 2016 (secondary to Friede/others); Arancibia-Carcamo 2017 |
| g-ratio | PNS ~0.6 (Rushton 1951 optimum); CNS ~0.77 (theory); MRG 0.605-0.791; rat optic nerve 0.78, cortex 0.81 | Rushton WAH (1951), *J. Physiol.* **115**(1), 101-122, DOI 10.1113/jphysiol.1951.sp004655; Chomiak T, Hu B (2009), *PLoS ONE* **4**(11), e7754, DOI [10.1371/journal.pone.0007754](https://doi.org/10.1371/journal.pone.0007754); Friede RL, Bischhausen R (1982), *Brain Res.* **235**(2), 335-350, DOI 10.1016/0006-8993(82)91012-5 |
| Myelin thickness | (D - d)/2: MRG 5.7 um fibre 1.15 um; 10 um 1.55 um; 16 um 1.65 um (our calc); CNS small axons ~0.1 um. Lamellae 80-150 (MRG), period 15.6 nm | MRG; Arancibia-Carcamo |
| Pig spinal-cord WM axons | average axon diameter ~3-5 um (fixed porcine cord, d-PFG MRI vs histology) | Komlosh ME et al. (2013), "Mapping average axon diameters in porcine spinal cord white matter and rat corpus callosum using d-PFG MRI", *NeuroImage* **78**, 210-216, DOI [10.1016/j.neuroimage.2013.03.074](https://doi.org/10.1016/j.neuroimage.2013.03.074) (abstract + secondary) |
| Cross-species SC morphometry | database of axon diameter, myelin thickness and density (human, cat, rat, mouse, others) | Saliani A, Perraud B, Duval T, Stikov N, Rossignol S, Cohen-Adad J (2017), "Axon and myelin morphology in animal and human spinal cord", *Front. Neuroanat.* **11**, 129, DOI [10.3389/fnana.2017.00129](https://doi.org/10.3389/fnana.2017.00129) (abstract) |
| Sheep spinal cord | no morphometry found | - |
| Na channels per node | rat sciatic: N0 = 21,000 per node, gamma = 14.5 pS (20 C) | Neumcke B, Stampfli R (1982), *J. Physiol.* **329**, 163-184, DOI [10.1113/jphysiol.1982.sp014296](https://doi.org/10.1113/jphysiol.1982.sp014296) (abstract) |
| Na channel density | rabbit sciatic: ~12,000/um^2 by saxitoxin binding (~700,000 per node), later revised to ~3000-5000/um^2; rat estimate 21,000 / (pi x 3 um x 1 um) ~ 2,200/um^2 (our calc); Zangari uses 1400 (rat) and 400-700 (rabbit) /um^2 | Ritchie JM, Rogart RB (1977), *PNAS* **74**(1), 211-215, DOI [10.1073/pnas.74.1.211](https://doi.org/10.1073/pnas.74.1.211) (abstract); revision from a Ritchie review (secondary) |
| Single-channel conductance | 14.5 pS rat node; 6.4 +/- 0.9 pS frog node (Sigworth 1980, *J. Physiol.* **307**, 97-129); 7.9 pS frog (Conti 1976); ~18 pS rat muscle (Sigworth F, Neher E (1980), *Nature* **287**, 447-449, DOI 10.1038/287447a0) | as listed (abstracts) |
| Nodal capacitance | MRG node 2 uF/cm^2; Arancibia model 0.9 uF/cm^2; generic ~1 uF/cm^2 | MRG code; Arancibia-Carcamo |
| Nodal conductances (MRG, 37 C) | gNaf 3.0 S/cm^2 (max), gNap 0.01, gKs 0.08, gLeak 0.007 S/cm^2 (rest), E_leak -90 mV, V_rest -80 mV | MRG AXNODE.mod |
| Nodal conductances (CNS model) | gNa 3 S/cm^2, gKs 80 mS/cm^2, gNap 5 mS/cm^2, node leak 80 mS/cm^2 | Arancibia-Carcamo 2017 |
| Peak open conductance per node | 21,000 x 14.5 pS = 0.30 uS, i.e. ~3 S/cm^2 over ~9.4 um^2 (our calc; matches gNaf) | - |
| Myelin electrical | 0.1 uF/cm^2 and 0.001 S/cm^2 per lamella membrane | MRG code |

**Ion-channel effect at THz (our calc):** peak nodal membrane conductance 3 S/cm^2 = 3x10^4 S/m^2. Spread over a 5 nm membrane this is a volume conductivity of 1.5x10^-4 S/m, against water's 119 S/m effective THz conductivity at 1 THz: a ratio of ~10^-6. Channel gating also takes us-ms, while the THz period is ~1 ps. Ion permeation (~10^7-10^8 ions/s per channel) is 4-5 orders slower than the field oscillation, so the DC channel conductance does not even enter the THz response as a conductivity. Xu et al. 2007 shows that changing bulk ion concentration by ~100 mM does not change THz absorption measurably.

---

## 4. Water-vapour absorption lines, 0.4-2.2 THz

Line centres are from the JPL Molecular Spectroscopy catalog, H2(16)O, tag 18003 (Pickett HM, Poynter RL, Cohen EA, Delitsky ML, Pearson JC, Muller HSP (1998), "Submillimeter, millimeter, and microwave spectral line catalog", *J. Quant. Spectrosc. Radiat. Transf.* **60**(5), 883-890, DOI [10.1016/S0022-4073(98)00091-0](https://doi.org/10.1016/S0022-4073(98)00091-0); file spec.jpl.nasa.gov/ftp/pub/catalog/c018003.cat, read 2026-11-02). logI = log10 integrated intensity (nm^2 MHz, 300 K); larger is stronger. HITRAN line positions agree to well below TDS resolution (e.g. HITRAN 556.936 GHz for 1_10-0_01).

All strong lines (logI > -1.2):

| f (THz) | cm^-1 | logI | transition J'Ka'Kc' - J"Ka"Kc" | In the requested list? |
|---|---|---|---|---|
| 0.556936 | 18.58 | -0.82 | 1 1 0 - 1 0 1 | yes (0.557) |
| 0.752033 | 25.09 | -1.00 | 2 1 1 - 2 0 2 | yes |
| 0.987927 | 32.95 | -1.14 | 2 0 2 - 1 1 1 | yes |
| 1.097365 | 36.60 | -0.32 | 3 1 2 - 3 0 3 | yes |
| 1.113343 | 37.14 | -0.84 | 1 1 1 - 0 0 0 | yes |
| 1.153127 | 38.46 | -1.13 | 3 1 2 - 2 2 1 | **missing; add** |
| 1.162912 | 38.79 | -0.28 | 3 2 1 - 3 1 2 | yes |
| 1.207639 | 40.28 | -0.79 | 4 2 2 - 4 1 3 | yes |
| 1.228789 | 40.99 | -0.86 | 2 2 0 - 2 1 1 | yes |
| 1.410618 | 47.05 | -0.38 | 5 2 3 - 5 1 4 | yes |
| 1.602219 | 53.44 | -0.69 | 4 1 3 - 4 0 4 | yes |
| 1.661008 | 55.41 | -0.28 | 2 2 1 - 2 1 2 | yes |
| 1.669905 | 55.70 | +0.10 | 2 1 2 - 1 0 1 | yes (strongest in band) |
| 1.716770 | 57.27 | +0.07 | 3 0 3 - 2 1 2 | yes |
| 1.762043 | 58.78 | -1.04 | 6 3 3 - 6 2 4 | yes |
| 1.794789 | 59.87 | -0.96 | 6 2 4 - 6 1 5 | near 1.797 |
| 1.797159 | 59.95 | -0.83 | 7 3 4 - 7 2 5 | yes |
| 1.867749 | 62.30 | -0.34 | 5 3 2 - 5 2 3 | yes |
| 1.918485 | 63.99 | -1.06 | 5 2 3 - 4 3 2 | pair with 1.919 |
| 1.919360 | 64.02 | -0.58 | 3 2 2 - 3 1 3 | yes |
| 2.040477 | 68.06 | -0.67 | 4 3 1 - 4 2 2 | yes |
| 2.074432 | 69.20 | -0.84 | 4 1 3 - 3 2 2 | **missing; add** |
| 2.164132 | 72.19 | -0.11 | 3 1 3 - 2 0 2 | **missing; add (strong)** |
| 2.196346 | 73.26 | -0.21 | 3 3 0 - 3 2 1 | **missing; add (strong)** |

Weaker lines in the band (logI -2 to -2.9): 0.448, 0.621, 0.916, 0.970, 1.158, 1.278, 1.296, 1.322, 1.441, 1.542, 1.714, 1.766, 1.881, 1.885, 2.016, 2.162 THz.

All 18 requested frequencies are confirmed.

TDS references:
- van Exter M, Fattinger Ch, Grischkowsky D (1989), "Terahertz time-domain spectroscopy of water vapor", *Opt. Lett.* **14**(20), 1128-1130, DOI [10.1364/OL.14.001128](https://doi.org/10.1364/OL.14.001128). First THz-TDS water-vapour spectrum, ~0.2-1.45 THz.
- Xin X, Altan H, Saint A, Matten D, Alfano RR (2006), "Terahertz absorption spectrum of para and ortho water vapors at different humidities at room temperature", *J. Appl. Phys.* **100**(9), 094905, DOI [10.1063/1.2357412](https://doi.org/10.1063/1.2357412).

Relevance to the Yerevan "0.6 and 2 THz resonances": 0.557 THz is the strongest low-frequency line, and 1.92-2.20 THz holds a dense group of strong lines (1.919, 2.040, 2.074, 2.164, 2.196). Berry 2003 saw a spurious 1.2 THz feature from unpurged air in the reference. A humidity difference between the reference and sample scans (minutes apart, sample evaporating under DC heating) would put features at exactly these places. Purging, or checking that the features line up with the table above, is the test.

---

## 5. DC conductivity of white matter; Joule heating, electro-osmosis, dehydration; TDS hydration artefacts

### 5a. White-matter / spinal-cord conductivity (low frequency)

| Tissue | sigma parallel | sigma perpendicular | Source (check) |
|---|---|---|---|
| Cat dorsal columns (spinal WM) | 138-212 Ohm cm, i.e. 0.47-0.72 S/m | 1211 Ohm cm, i.e. 0.083 S/m | Ranck JB Jr, BeMent SL (1965), "The specific impedance of the dorsal columns of cat: an anisotropic medium", *Exp. Neurol.* **11**(4), 451-463, DOI [10.1016/0014-4886(65)90059-2](https://doi.org/10.1016/0014-4886(65)90059-2) (abstract) |
| Cerebral WM (cat internal capsule) | ~1.1-1.2 S/m | ~0.12 S/m | Nicholson P (1965), "Specific impedance of cerebral white matter", *Exp. Neurol.* **13**(4), 386-401, DOI [10.1016/0014-4886(65)90126-3](https://doi.org/10.1016/0014-4886(65)90126-3); compiled in Geddes LA, Baker LE (1967), *Med. Biol. Eng.* **5**(3), 271-293, DOI [10.1007/BF02474537](https://doi.org/10.1007/BF02474537) (values secondary, < 1-20 kHz) |
| Spinal-cord stimulation models | 0.6 S/m | 0.083 S/m | (secondary, quoting Ranck & BeMent) |
| CSF / saline (reference) | 1.79 S/m at 37 C | isotropic | Baumann et al. 1997 |

Anisotropy is ~5-10x. Excised, saline-soaked slices will sit nearer the saline value.

### 5b. Joule heating under the Yerevan DC protocol (our calc)
100 V over 18 mm gives E ~ 5.6 kV/m if the whole drop is across the tissue (an upper bound; electrode polarisation and contact resistance take part of it). Power density p = sigma E^2. With rho c ~ 3.7x10^6 J m^-3 K^-1, the adiabatic heating rates are:

| sigma (S/m) | p (W/m^3) | dT/dt adiabatic |
|---|---|---|
| 0.083 (perpendicular WM) | 2.6x10^6 | 0.7 K/s |
| 0.6 (parallel WM) | 1.9x10^7 | 5 K/s |
| 1.79 (saline film) | 5.5x10^7 | 15 K/s |

Even with strong conductive loss to the PTFE cell, a few-minute application can heat and dry a 70-100 um slice. The reported time course (develops over ~5 min, recovers over ~1 h) fits heating-driven drying and rehydration. Water's THz constants depend weakly on T (Liebe, our calc: 20 to 37 C moves alpha(1 THz) from 223 to 228 cm^-1 and n(0.5 THz) from 2.31 to 2.37). Loss of water changes them strongly.

Background literature:
- Davalos RV, Rubinsky B (2008), "Temperature considerations during irreversible electroporation", *Int. J. Heat Mass Transf.* **51**(23-24), 5617-5622, DOI [10.1016/j.ijheatmasstransfer.2008.04.046](https://doi.org/10.1016/j.ijheatmasstransfer.2008.04.046). Joule heating in tissue under applied fields.
- Romanenko et al. 2014 (1.9): thermal equivalence of mm-wave neuronal effects.

### 5c. Electro-osmosis and electrolysis in tissue under DC
- Guy Y, Sandberg M, Weber SG (2008), "Determination of zeta-potential in rat organotypic hippocampal cultures", *Biophys. J.* **94**(11), 4561-4569, DOI [10.1529/biophysj.107.112722](https://doi.org/10.1529/biophysj.107.112722). Brain tissue has a negative zeta potential (later -22.8 +/- 0.8 mV, tortuosity 2.24), so a DC field drives electro-osmotic water flow through it.
- Pikal MJ (2001), "The role of electroosmotic flow in transdermal iontophoresis", *Adv. Drug Deliv. Rev.* **46**(1-3), 281-305, DOI [10.1016/S0169-409X(00)00138-1](https://doi.org/10.1016/S0169-409X(00)00138-1). Electro-osmotic volume flow in tissue under DC (toward the cathode for negatively charged tissue).
- Nilsson E, von Euler H, Berendson J, Thorne A, Wersall P, Naslund I, Lagerstedt A-S, Narfstrom K (2000), "Electrochemical treatment of tumours", *Bioelectrochemistry* **51**(1), 1-11, DOI [10.1016/S0302-4598(99)00073-2](https://doi.org/10.1016/S0302-4598(99)00073-2). Review of DC tissue treatment: electrolysis at electrodes (acid and Cl2 at the anode, alkali and H2 at the cathode) and electro-osmotic water transport (anode side dries, cathode side swells; described in the EChT literature, secondary here).
- With 100 V on tungsten electrodes in saline-wet tissue, electrolysis and gas evolution are expected. pH fronts and water redistribution across an 18 mm gap would change the THz path's water content over minutes.

### 5d. THz-TDS artefacts from hydration
- Png GM, Choi JW, Ng BW-H, Mickan SP, Abbott D, Zhang X-C (2008), "The impact of hydration changes in fresh bio-tissue on THz spectroscopic measurements", *Phys. Med. Biol.* **53**(13), 3501-3517, DOI [10.1088/0031-9155/53/13/007](https://doi.org/10.1088/0031-9155/53/13/007). Residual hydration "can vastly alter" extracted THz constants. Drying changes sample thickness and texture, which corrupts n and alpha. Proposes lyophilisation.
- Berry et al. 2003: an unpurged spectrometer shows water-vapour features (1.2 THz) through reference/sample mismatch.
- Yamaguchi 2016 (dry-air purge plus quartz cover), Kucheryavenko 2021 (PE film), Gavdush 2019 (gelatin embedding): standard controls against drying during measurement.
- A control for any field-induced spectral change: run the same time course with the field applied in a sham geometry (or with an ion-free control sample) under N2 purge, logging cell temperature and sample mass/thickness.

---

## 6. What a study with full water dispersion, realistic sizes, absorbed power per node, mesh convergence and a Mie benchmark adds (combined)

1. **Band:** first full-wave treatment of 0.1-3 THz. All prior node/waveguide models sit at >= 10 THz or in the optical.
2. **Materials:** dispersive water (double-Debye/Ellison) in axoplasm and extracellular space, instead of constant real permittivities (Song 2026) or one complex constant (Guo 2022). Myelin from a water/lipid mixture rather than the mid-IR "2x" index.
3. **Geometry:** node ~1 um, internode 0.2-1.5 mm, axon 1-15 um, myelin 0.5-3 um (MRG/Arancibia ranges), rather than 100 um to 1 mm segments with single sizes.
4. **Observable:** absorbed power per node (and node vs internode vs bath) under plane-wave illumination, the quantity a THz-TDS absorption experiment measures. Prior papers report S-parameters or field amplitude from internal current sources.
5. **Verification:** mesh convergence and an energy balance (absorbed + scattered = extinguished). Neither appears in any prior node/THz paper. A Mie (coated-cylinder/sphere) benchmark validates the solver in the same lossy, dispersive regime.
6. **Physics:** shows that at 0.1-3 THz the node is deeply sub-wavelength, the sheath is anti-guiding (n_myelin < n_water), and nodal absorption follows quasi-static field ratios times water's eps''. Ion-channel conductance enters at ~10^-6 of water's loss. This gives a quantitative null against the reported 0.6/2 THz resonances, which Section 4-5 traces to plausible water-vapour and DC-hydration artefacts.

## 7. Gaps not closed here
- Liu 2020 full text (geometry, band, any numbers) is paywalled and was not read.
- Gavdush 2021 double-Debye parameters for intact brain are in a table that was not retrieved; Ronne 1997 and Yada 2008 fit values were not retrieved.
- No direct 0.1-3 THz measurement of myelin or of spinal-cord white matter n/alpha (outside the Yerevan data) was found.
- Nicholson 1965 values and Baumann 1997 CSF value are from secondary quotes or recall, not the original text.
