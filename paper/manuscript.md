# Below 3 THz the node of Ranvier absorbs quasi-statically: an interpolation law with exact limits, bounds on voltage-induced changes, and a re-analysis of published spinal-cord THz waveforms

Aram Poghossian

Target journal: *Journal of Infrared, Millimeter, and Terahertz Waves* (Springer).
Figures are in `paper/figures/`; every number below traces to a script and a result file listed in
`REPORT.md` (Sims 27–44) and to `results/run_log.csv`.

## Abstract

Biased (100 V DC) spinal-cord slices reportedly show a 9–20 dB THz transmission dip drifting from
0.7 to 0.5 THz, attributed to resonant nodes of Ranvier or DC-opened channels. We test both with a
full-wave model of a myelinated fibre in water, benchmarked to 2 × 10⁻⁴. In 0.1–3 THz the node does
not resonate. Its absorbed power follows a quasi-static law (exact limits, three constants) in node
length over myelin thickness, g-ratio and permittivity contrast, within 2.3 % of full wave for
fibres up to 10 µm, independent of internode length. For the THz field along the fibres (every
published scan), no node state changes a 100 µm slice by more than 1.1 dB, and with every extreme
combined none gives the extra delay the measured change needs (≤ 0.06 of ≥ 0.23 rad at 0.6 THz).
Open channels add 10⁻⁶ of the membrane admittance. The three reports share reference waveforms; in
their four biased scans the pulse gains a component at the empty-cuvette delay, 4–22 % of the empty-
cuvette pulse after 0.5–5 min, plus a 9–18 % longer optical path through the tissue. Together they
are consistent with the 0.6 THz dip. That needs more tissue in the beam: no change of the tissue's
permittivity (heating, ions, axoplasm conduction of 23–46 S/m, water content) gives its loss and
delay together. A 2.0 THz dip in one scan is not reproduced. The undelayed component's source and
the reported polarisation and bias-direction selectivity (shown without data) remain open.

**Keywords:** terahertz; node of Ranvier; myelinated axon; finite-element method; dosimetry;
terahertz time-domain spectroscopy; spinal cord.

## 1. Introduction

Terahertz (THz) radiation is studied as a probe of nerve tissue and as a possible neuromodulator
(reviews: Wilmink and Grundt 2011; Romanenko et al. 2017).
Both uses need to know where in a nerve fibre THz power is absorbed. Two pictures of the myelinated
fibre suggest that the answer may be unusual. In the first, the myelin sheath is a dielectric
waveguide in the mid-infrared to THz range and the nodes of Ranvier are where the guided field
leaks out or is amplified (Liu G et al. 2019; Liu Y et al. 2020; Xiang et al. 2020, 2021). In the
second, the node is an antenna or impedance transformer with resonances set by its length (Zangari
et al. 2018; Guo et al. 2022; Song et al. 2026). Kumar et al. (2016) proposed optical channels along
myelinated axons on the same grounds. These models work at 10 THz and above or in the optical range,
and most treat the surrounding fluid as lossless or with a constant permittivity.

Experiments now report a THz effect in this tissue. Hovhannisyan and co-workers (three reports,
2024, refs. H1–H3; a fourth, H4, was not available to us) measured 70–100 µm slices of pig and sheep thoracic spinal-cord white
matter by THz time-domain spectroscopy (TDS) with 100 V DC across 15–18 mm. They report an absorption
dip of 9–20 dB near 0.6 THz that drifts from 0.7 to 0.5 THz over minutes, changes near 2 THz,
onset in ~5 min and recovery in ~1 h, and an effect that appears only when the THz field is parallel
to the fibres and the DC field across them. They explain it by DC-opened sodium channels at the
nodes, which would turn fibres into nonlinear conducting wires and the slice into a diffraction
grating.

Below 3 THz water dominates: its permittivity is 4.1 − 2.1j at 1 THz, and its power absorption
length is 35–118 µm. Whether a node can resonate inside such a host, and whether any node or fibre
state could change a 100 µm slice by 9 dB, are questions a full-wave calculation can answer. We ask:

1. Does a single myelinated fibre with realistic dimensions, in water with measured dispersion,
   absorb resonantly at the node in 0.1–3 THz?
2. If not, what law sets the absorbed power?
3. Can any node state, or the fibre state the experimenters propose, change the transmission of a
   100 µm slice by the reported amount?
4. What do the published waveforms themselves show?

The answers are: no resonance; a quasi-static law, exact in its limits; no for the THz field
along the fibres, the configuration of every published scan (with every extreme combined, no node
state gives the extra delay the 0.6 THz feature needs); and an
undelayed beam component plus a thicker tissue path, which are consistent with the 0.6 THz dip but
not with the 2 THz dip. Section 2 describes the model, Section 3 its verification, Section 4 the node results,
Section 5 the experimenters' model and the published data, and Section 6 the limits.

**Prior work and what is new.** Quasi-static field factors of cell and membrane compartments are
standard in radio-frequency and microwave microdosimetry (Gimsa and Wachner 1999; Kotnik and
Miklavčič 2000), and both limits of the law below follow from the boundary conditions (continuity of
normal D and of tangential E). What is new is (i) the interpolation between the limits, with three
constants fitted once to a Laplace solver and tested against full wave over node length, myelin
thickness, axon radius and frequency, for the node in water with measured THz dispersion; (ii)
bounds on what any node state can do to a slice, including the negative-ε ring resonance; and (iii)
a re-analysis of the published waveforms. Table 0 lists the closest work.

| Paper | What they did | Band | Water | This work |
|---|---|---|---|---|
| Liu G et al. 2019 | myelin index ~2× surroundings; sheath as dielectric waveguide (SR-FTIR) | mid-IR | measured, mid-IR | below 3 THz myelin's index is below water's and the axon-core V number ≤ 0.4: guided modes leak |
| Liu Y et al. 2020 | node as an active (three-level gain) amplifier | mid-IR–THz | not stated | passive response only; gain is outside scope |
| Guo et al. 2022 | three-layer fibre modes, mode matching at the node, vs CST | 0–30 THz | constant complex ε | double-Debye water, realistic internodes, node power, Mie benchmark, convergence |
| Song et al. 2026 | node as impedance transformer, peaks vs node length | 10–40 THz | lossless, ε = 1.44 | with water's ε the node is quasi-static and non-resonant below 3 THz |
| Zangari et al. 2018 | node as an array of nano-antennas (FEM) | IR | lossless or constant | no antenna resonance below 3 THz; interpolation law with exact limits |
| Kumar et al. 2016 | myelinated axons as optical waveguides | visible–NIR | n/a | the waveguide picture does not carry over to THz |
| H1–H3 (2024) | THz-TDS of biased cord slices; dips at 0.6 and 2 THz; open node channels | 0.2–4 THz | experiment | tests their model and node resonance; re-analyses their waveforms |

## 2. Model

### 2.1 Geometry

The fibre is an axon of radius a, a myelin sheath to radius b (g-ratio g = a/b) and a node gap of
length Ln repeated with internode length Lint (Fig. 1a). The reference fibre has a = 3.5 µm,
b = 5 µm (g = 0.7), Ln = 1 µm and Lint = 1000 µm, within the ranges of the MRG model (McIntyre
et al. 2002: internodes 500–1500 µm for 5.7–16 µm fibres) and of measured node lengths
(Arancibia-Carcamo et al. 2017: 1.08 ± 0.02 µm, range 0.5–3.7 µm). The scaling study (Sec. 4.3)
covers node length 0.5–5 µm, myelin thickness 0.5–3 µm, axon radius 1–7.5 µm and internode
200–1500 µm. Node constriction (axon radius 3.5 → 2 µm inside the node) and a 5 nm axolemma are
treated as variants.

### 2.2 Materials

Water follows the double-Debye model of Kindt and Schmuttenmaer (1996); Liebe et al. (1991) gives
the temperature dependence and two Lorentz bands (5.1 and 18.2 THz) for a broadband check; the
double-Debye fit of Jepsen et al. (2007) is shown for comparison in Fig. 1b. Myelin
is a Looyenga mixture of 40 % water with dry lipid (ε ≈ 2.3 − 0.05j), and axoplasm one of 85 %
water. In 0.1–3 THz this gives indices of 1.91–3.42 (water), 1.86–3.12 (axoplasm) and 1.67–2.25
(myelin) (Fig. 1b, c). Myelin is lamellar; a uniaxial Wiener tensor (radial harmonic, axial
arithmetic mean) is tested as a variant, with and without confined interlamellar water at half the
bulk loss. Every material is also scaled by ±20 % in ε′ and ε″, one at a time, and for the slice bound
with the three changes that raise it combined. All loss enters as a
literal complex permittivity, tested against analytic results.

### 2.3 Full-wave method

The fibre is solved with the finite-element method (COMSOL 6.3, electromagnetic waves frequency
domain) in 2D-axisymmetric form with azimuthal harmonics m = 0…3. A plane wave in the lossy host
is the background field, and the fibre's scattered field is solved for; using the exact lossy-host
plane wave (not the vacuum one) matters, as Sec. 3 shows. The node chain is periodic along the
axis; a cylindrical perfectly matched layer (PML) closes the 100 µm water radius. Both
polarisations are solved: E∥ (field along the fibre, m = 0) and E⊥ (across, m = ±1 dominant). The
last retained order contributes ≤ 5 × 10⁻⁴ of the absorbed power.

### 2.4 Observables

- Node absorbed power P_node = ½ωε0∫ε″|E|²dV over the gap, and the node field factor
  η = ⟨|E|²⟩_node / |E0|², where E0 is the incident field in water at the fibre axis.
- The node's share of fibre absorption per period and its absorption cross-section.
- For tissue: three slice observables per 100 µm, from the change of (i) fibre absorption,
  (ii) absorption inside a 25 µm radius around the fibre (fibre plus near-field water), and
  (iii) the fibre's induced dipole. The dipole change gives an effective-permittivity change Δε =
  Δp/V_cell, from which the slab transmission is computed exactly (transfer matrix on PTFE,
  n = 1.43). Absorbed power alone misses how the host loss couples to the fibre (a factor |n|²/n′²).
  The node bounds use the first-order change −k0 Im(Δε/n)·d; at their largest value the exact slab
  gives the same (0.1223 vs 0.1227 dB). For the large Δε of the conducting-axon case (Sec. 5.2),
  where first order fails, the exact slab is used.
  V_cell scales with packing / internode; values are given for internode 1000 µm and packing 0.5 and
  scaled to nominal values (internode 200 µm, packing 0.6: ×6); Sec. 4.6 adds the extremes.

## 3. Verification

The legacy 3D pipeline of this project produced "resonances" at 0.6 and 2 THz that turned out to be
artefacts (supplementary Fig. S1): its frequency unit defaulted to THz while Hz values were passed,
so every solve ran at 10²³–10²⁴ Hz; its background field was the vacuum plane wave inside water; the
on-axis field it reported was identical with and without the fibre; and its energy balance was off
by 12 orders of magnitude. Everything here was rebuilt and checked against problems with known
answers (Sim 28; Table 1):

| Test | Reference | Agreement |
|---|---|---|
| B1: water slab, 30 frequencies | transfer matrix | R + T + A − 1 ≤ 8 × 10⁻⁶; T error 9 × 10⁻⁵ |
| B2: 2D layered cylinder, TM and TE | Mie series | absorbed power per layer ≤ 8 × 10⁻⁶; flux balance ≤ 5 × 10⁻⁷ |
| B3: axisymmetric model, m = 0…3, both polarisations | per-order series | ≤ 2 × 10⁻⁴; balance ≤ 2 × 10⁻⁵ |
| Null test (fibre = water) | zero scattering | ≤ 10⁻⁷ |
| Energy balance, every production solve | source = absorbed + outflow | ≤ 2.4 × 10⁻⁵ (≤ 2.3 × 10⁻⁴ in all review runs) |
| Mesh, three levels (113k / 267k / 731k DOF) | Richardson | P_node L1 → L3 ≤ 0.51 % (E∥), ≤ 0.12 % (E⊥) |
| Water radius 100 → 150 µm; PML 60 → 100 µm | — | 7 × 10⁻⁵; 10⁻⁶ |
| Internode 200 vs 1000 µm | — | η equal to 1.3 × 10⁻⁵ |

Three bugs were found on the way and fixed: flux read on the PML interface counts half; axisymmetric
integrals need the 2πr weight; the axisymmetric PML must be cylindrical. A first-order scattering
boundary instead of the PML is exact only for m = 0 and errs as m²/(k_w R)² (E⊥ wrong by 2.9× at
0.1 THz for R = 100 µm), so the PML is used; it matches the exact series and an independent Laplace
solver.

## 4. Results: the node

### 4.1 No resonance in 0.1–3 THz

For the reference fibre (30 frequencies, both polarisations) η rises smoothly from 0.32 to 0.73
for E∥ and falls from 1.75 to 1.18 for E⊥ (Fig. 3). The node absorbs 0.03–0.12 % of the fibre's
absorption per period; its cross-section is 0.1–0.8 µm² (E∥). An independent quasi-static solver
with no fitted parameter reproduces η within 0.1–6 % (E∥) and 0.5–0.8 % (E⊥) at every frequency.
For E∥ the node field is below the incident one: the node is shielded by the myelin, not fed by it.

### 4.2 A quasi-static law with exact limits

The node is a sub-wavelength water gap (k_w·Ln < 0.15) in a quasi-static field, so

P_node = ½ ω ε0 ε″_w |E0|² V_n F,  V_n = π(b² − a²) Ln,

where F is the mean |E|² in the gap relative to the incident field, from a Laplace problem with no
length scale. F therefore depends only on x = Ln/(b − a), the g-ratio and the permittivities. It
has two exact limits. For E∥, a thin gap carries the normal displacement of the myelin, F0 =
|ε_my/ε_w|², and a long gap is bulk water, F = 1. For E⊥, a thin gap takes the tangential field of
the myelin shell of an infinite coated cylinder (F0 = S, an exact 2D series), and a long gap again
gives F = 1. Between them

F = F0 + (1 − F0) xᵖ / (xᵖ + c),  c = c0 |ε_w/ε_my|^q,

with three constants per polarisation (E∥: c0 = 1.066, p = 0.934, q = 0.742; E⊥: 0.999, 0.826,
−0.784) fitted once to the Laplace solver over x = 0.2–8, g = 0.6–0.8 and 0.1–3 THz (maximum error
2.1 %). They are not fitted to full-wave results or to any experiment. No length in the law compares
with the wavelength, and the internode length does not enter: there is no resonance condition.

### 4.3 Scaling

Full wave at 0.1/0.6/1/2/3 THz for the reference fibre and ten variants (Fig. 4): node length
0.5/2/5 µm; myelin thickness 0.5 and 3 µm (g = 0.875 and 0.54, outside the fitted 0.6–0.8); axon
radius 1 and 7.5 µm at g = 0.7; internode 200/500/1500 µm. Every spectrum is smooth. Full wave
over the Laplace factor is 0.998–1.036 for fibres up to b = 6.5 µm and 1.089 for b = 10.7 µm at
3 THz, where k_w·b = 1.3 and retardation begins. The law matches full wave within 2.3 % for
b ≤ 5 µm (fibre diameter ≤ 10 µm, including both g-ratios outside the fit range), 4.3 % at b = 6.5 µm and 9.8 % at b = 10.7 µm and 3 THz. The
internode length changes P_node by < 10⁻⁴ (200–1500 µm).

### 4.4 Why there is no resonance

- The node is short: k_w·Ln < 0.15 at 3 THz.
- No Fröhlich pole: every Re ε is positive in 0.1–3 THz, so no localised plasmon-like mode exists.
- No guiding: myelin's index (1.67–2.25) is below water's (1.91–3.42) at every frequency, so the
  sheath cannot guide by total internal reflection (molecular dynamics of a lipid bilayer gives the
  same ordering below ~20 THz; Zhang et al. 2022). Axon plus sheath is a core-cladding pair, but
  its V number is 0.05–0.18 for a 3.5 µm axon (≤ 0.39 at 7.5 µm); the mode spreads into the water,
  whose index is higher still, and leaks. The node power is independent of internode length, so no
  guided energy arrives from the internode.
- Water damping: the field's 1/e length is 70–235 µm, shorter than most internodes.

The waveguide and antenna pictures of Sec. 1 therefore do not carry over to 0.1–3 THz. They may hold
at 10–40 THz or in the infrared, which this work does not test.

### 4.5 Robustness

- **Materials ±20 %:** P_node changes by −24 % to +34 %, smoothly at all frequencies; the
  law predicts each change within 5 %.
- **Lamellar myelin:** node power changes by +1 to +3 % above 0.4 THz and by +26 to +29 % at 0.1 THz;
  with interlamellar water at half the bulk loss, by −6 % to +33 %.
- **Broadband water** (Liebe Lorentz bands added): η changes ≤ 6 %; node power follows the water
  loss (×2 at 3 THz), as the law says.
- **Constriction:** ≤ 1.2 % against the Laplace solver.
- **Axolemma:** ≤ 0.2 %.
- **Dense packing** (60 % area fraction): the E⊥ shell factor falls by 10 to 29 %; E∥ is unchanged.
- **Thick fibres** (b = 10.7 µm) and **dense frequency sweeps** (16 frequencies): broad maxima
  only (E⊥ near 2.2 THz, width > 2 THz).

### 4.6 Node states and the slice bound

The node gap was given any extra conductivity from 10⁻⁴ to 3 × 10⁶ S/m (12 values) and any passive
permittivity ε′ − jε″ with ε′ from −30 to 30 and ε″ from 0.01 to 30 (conductivity in Fig. 5).

**Conductivity.** P_node rises by up to 43× (E⊥) and then falls as the conducting gap screens its
own field. With the node refined 2×, results change by ≤ 3 × 10⁻⁴ up to 3 × 10⁶ S/m. Over the three
slice observables, the largest change at internode 1000 µm and packing 0.5 is ≤ 0.006 dB for E∥ at
0.6 and 2 THz (0.008 dB at 3 THz) and 0.18 dB for E⊥. Scaled to internode 200 µm and packing 0.6,
this gives ≤ 0.05 dB (E∥) and ≤ 1.1 dB (E⊥); E⊥ reaches 3.5 dB only for a metallic node in a 21 µm
fibre with 200 µm internodes, about 0.5 dB at that fibre's realistic 1.5 mm internode. Fully open
sodium channels correspond to ~10⁻⁴ S/m and give 10⁻⁹ dB.

**Passive permittivity.** A negative-ε ring resonance near ε′ ≈ −1.5 raises η to ~100 for E∥, but
the surrounding water's loss caps its effect on the slice. On a dense ε′ grid (−2 to −1, ε″ = 0.1)
at 0.6/1/2/3 THz the E∥ peak is 0.014/0.028/0.072/0.123 dB per 100 µm (×6: 0.09/0.17/0.43/0.74 dB).
Three mesh levels at the peak give 0.123/0.107/0.097 dB at 3 THz: the value falls with refinement
(the gap's corners are near-singular at this contrast). Re-scanning ε′ on the 2× finer mesh gives
peaks of 0.014/0.027/0.066/0.112 dB, each below the coarse peak (at 3 THz it moves to ε′ ≈ −1.35;
at 0.6 THz the fine peak at ε′ = −1.8 is bracketed by −2.0, 0.013 dB), so the coarse values are
upper bounds.
Re-scanning ε′ with each material scaled by ±20 % one at a time moves the peak to at most
0.148 dB (myelin ε′ ×1.2). The three changes that raise it (axon ε′ ×1.2, myelin ε′ ×1.2, myelin
ε″ ×0.8), applied together, raise it ×1.45 at 0.6 THz and ×1.43 at 3 THz (0.021 and 0.176 dB;
the 0.6 THz peak moves to ε′ ≈ −2.0, bracketed by a scan from −2.7 to −1.0). Including this, the
passive bound at internode 200 µm and packing 0.6 is 0.12/0.24/0.63/1.06 dB at 0.6/1/2/3 THz
(×1.45 assumed at 1 and 2 THz). At ε″ = 0.01 the gap's 90° corners are singular and the points do not converge, so the
claim is made for ε″ ≥ 0.1; negative ε′ with ε″ ≪ 1 is not possible in an aqueous gap at THz, where
ion collisions make ε″ ≫ |ε′|.

**Bound.** For E∥, internode ≥ 200 µm and packing ≤ 0.6, no node conductivity changes a 100 µm
slice by more than 0.05 dB, and no passive node permittivity by more than 0.1/0.2/0.6/1.1 dB at
0.6/1/2/3 THz (combined ±20 % materials included). The factors that could raise these multiply:
internodes of 50 µm (the smallest CNS fibres; ×4), packing 0.8 (×1.33, assumed linear; at 0.6 the
E∥ response is unchanged, Sec. 4.5) and the optical thickness of the measured slice (Sec. 5.5;
×2.2 at 0.6–1 THz, ×3.1 at 2–3 THz, the larger of the loss and delay values with the model's
water). Taken together they give 1.5/2.9/10/17 dB for passive ε and ≤ 1.2 dB for conductivity
(×1.45 materials assumed). The internode scaling was checked by solving at 50 µm directly: the
slice change matches the value scaled from 200 µm within 0.3 % at every frequency and in both
polarisations.

The right comparison at 0.6 THz is not the dip depth: in the two-term reading of Sec. 5.6 the dip is
an interference null whose depth depends on how exactly two terms cancel. What the tissue itself
must change is its transfer, T → T^q with q − 1 = 0.09–0.18, which at 0.6 THz (|T|² = −12.1 dB,
delay phase −arg T = 2.5 rad) is 1.1–2.2 dB of extra loss together with 0.23–0.45 rad of extra
delay (0.21–0.41 rad without the phase offset, Sec. 5.5). Over all
116 node states with ε″ ≥ 0.1 (conductivity scan, ε grids, dense ε′ scans, and the combined-materials
scan from ε′ = −10 to −1), the first-order change of the slice at 0.6 THz is mostly a phase advance:
near the ring resonance, where the loss peaks, the phase falls. The largest extra delay any state
gives is 0.06 rad with every extreme combined, 3.3× short of 0.21 rad, and the state with the
largest loss (1.5 dB with every extreme) advances the phase. No state gives both. Under nominal
parameters (internode 200 µm, packing 0.6, 220 µm path) node states give ≤ 0.3 dB, ≥ 4× short.
At 2 THz the extreme bound (10 dB) approaches the 14 dB change of the 30-s scan, itself a lower
bound (Sec. 5.6); at nominal parameters the margin is 7×. Reaching it at the extremes would also
need every node to hold ε′ ≈ −1.4 at 1.95–2.05 THz only, a material line ≤ 0.1 THz wide, since for
any fixed node material the response is smooth over 0.1–3 THz (Sec. 4.5). At 30 s the 0.6–0.7 THz
change needs no tissue change at all (q = 1). Near 3 THz the extreme bound is added loss, while the
reported change there is a rise in transmission of up to 17–18 dB; over the full ε grid at 3 THz
(33 states) the largest rise any node state gives is 0.003 dB per 100 µm at 1 mm internodes,
≤ 0.4 dB with every extreme. The extremes also need every node, in a slice made only of the
smallest fibres, to hold ε′ between −2 and −1.4 with ε″ = 0.1, which an aqueous gap cannot (above).

These bounds are for the THz field along the fibres, the configuration of every published scan
(H1, p. 18). For E⊥ the conductivity bound with all extremes would be 3.0/11.5/18 dB at
0.6/2/3 THz, so the E∥ margin does not carry over to the reported E⊥ controls, which were described
without data (Sec. 5.6).

## 5. The experimenters' model and the experiment

### 5.1 Channels cannot change THz current paths

At 0.6–2 THz the node membrane's capacitive admittance ωC_m is 4 × 10¹⁰ – 1.3 × 10¹¹ S/m². All
sodium channels open (3 × 10⁴ S/m², MRG) add 1.3 × 10⁶ – 4.2 × 10⁶ times less. The axoplasm is
continuous through the node in any state. Channel state therefore cannot change THz current paths.
Ion concentration changes do not help either: at physiological concentration salt changes water's
absorption only slightly over 0.15–3.7 THz (Xu et al. 2007).

### 5.2 Fibres as conductors

The proposed alternative is that the whole fibre conducts. Adding conductivity to the whole
axoplasm (full wave; the axon mesh refined 4× changes results by ≤ 7 × 10⁻⁴ up to 10⁵ S/m) and
converting the fibre dipole to an exact 100 µm slab on PTFE, a 9 dB change for E∥ needs +470 to
+565 S/m at packing 0.5 and 290–980 S/m over packing 0.7–0.3 (0.6–2 THz); 0.1 THz needs
600–1300 S/m. Physiological axoplasm conducts 1.4 S/m; saturated NaCl about 25 S/m; the axoplasm's
own THz water loss ωε0ε″ is already 82–115 S/m. Quasi-static Wiener mixing matches full wave within
3 % at ≤ 1 THz and within 10 % at 3 THz up to 400 S/m and fails above ~10⁴ S/m, where the axon's
inductance matters; at the conductivity needed, the fibres are resistive, not wires. Such a
conducting axon would also be E∥-selective (6.6 vs 2.7 dB at 0.6 THz), so polarisation selectivity
alone does not discriminate between mechanisms.

Measured against what the waveform fits need (Sec. 5.6: T → T^q, extra loss and extra delay
together) rather than against the dip depth, conduction fails differently. In a 220 µm slab
(Sec. 5.5) 23–46 S/m of extra axoplasm conductivity gives the required 1.1–2.2 dB at 0.6 THz
(16–33× physiological, 1–2× saturated NaCl), but only 0.12 of the required delay at 0.6 THz and
0.04–0.07 at 1–1.4 THz (Sim 46, Fig. 9). Conductivity adds loss with almost no phase.

### 5.3 No grating orders, no picosecond nonlinearity

Fibre spacings of 2–23 µm put the first diffraction order above 6 THz in tissue and 13 THz in air,
and white matter is not periodic. Sodium gating takes ≥ 10 µs; ions drift ≤ 5 × 10⁻¹⁶ m in 1 ps;
TDS fields of 10³–10⁴ V/m are 10⁻³–10⁻⁴ of the Wien-effect onset. The DC bias breaks inversion
symmetry, but χ⁽²⁾_eff = 3χ⁽³⁾E_dc gives ≤ 10⁻¹¹ of the linear response. The quadratic current the
authors fit equals the linear one at 13 % of the pulse peak. The filter-order test they use does not
prove a nonlinearity: reflections in the filter cascade, changes of the beam path and a sample that
changes between scans all give the same signature.

### 5.4 Bulk effects of the bias

An effective-medium slab of 100 µm aligned white matter (water fraction 0.77) on PTFE transmits
−6.2/−9.6 dB (E∥) and −5.8/−9.2 dB (E⊥) at 0.6/2 THz; the dichroism is only 0.35 dB. The bias
drives ~0.7–1.1 mA through the slice: 2.6–3.7 MW/m³ of Joule heat (σ⊥ = 0.083 S/m; Ranck and BeMent
1965), 0.6–0.9 K/s adiabatic, enough to evaporate ~9 % of the slice water per minute. Heating by
5/15 K changes 0.6 THz transmission by −0.13/−0.43 dB; losing 10/30 % of the water by +0.6/+1.5 dB
(0.6 THz) and +0.9/+2.8 dB (2 THz). Electrolysis makes 25–42 mm³ of H₂ in 5 min at the electrodes,
7.5 mm from the beam. Applied uniformly, every bulk mechanism is smooth in frequency and none
gives a narrow 9 dB feature (partial coverage is a different case, Sec. 5.6). DC along the fibres (σ∥ ≈ 7σ⊥) heats more and reportedly does nothing, so heating alone
does not explain the DC-direction selectivity. Against the change of T the fits need (Sim 46,
Fig. 9; 220 µm slab, Kindt–Schmuttenmaer water as base, its own loss within 1.6 dB of the measured), no change of the tissue's permittivity fits.
Heating, with Liebe's temperature dependence (fitted below 1 THz and 30 °C) or with it scaled ×2.31
to the measured 4–21 °C absorption change at 1 THz (46 cm⁻¹; Millon et al. 2022), gives the 0.6 THz
loss at 9–29 K but only 0.08–0.60 of the loss at 1.4 THz (rms misfit over loss and delay at 0.6, 1.0
and 1.4 THz 37–56 %). Ions in the tissue water and axoplasm conduction give loss without delay
(0.04–0.12 of it). More or less water at fixed thickness gives 0.17–0.22 of the loss and 0.07–0.11
of the delay. Only a longer path fits: the same tissue 8–17 % thicker gives 0.87–1.11× of the
requirement (7 % rms), and a water gain of 8–17 % whose volume goes entirely into thickness gives
0.83–1.01× (10 % rms), 0.67–0.86 of it from the added thickness; spread isotropically, the same gain
gives 0.32–0.47×. Adding 1–4 K of warming to such a gain brings the misfit to 5 %. What the fits need
is more tissue in the beam, not a different tissue. The tissue is post-mortem, and the bias puts ~70 mV
across a 10 µm fibre, not the ~1 V the authors estimate.

### 5.5 One session, and a slice optically thicker than stated

The three reports (H1: 70 µm slice, 18 mm electrode gap; H2 and H3: 100 µm) show the same
reference waveforms. Traced, the empty-cuvette pulses of H1 and H3 differ by 2.2 % rms and their
minima by 0.005 ps; the no-voltage pulses of H1, H2 and H3 differ by 1.8–3.7 % with minima within
0.015 ps (Sim 43). One scan traced from two figures of H1 differs by 3.4 %, so the reference
waveforms of the three reports are identical within tracing noise: they come from one sample, and
possibly from the same scans. We treat them as one
sample with four biased scans, after 30 s and 4 min (H1), 2 min (H2) and 5 min (H3). Whether these
come from one continuous bias run cannot be checked, so the trend of b with bias time below is a
sequence of four scans, not a measured time series. The stated thickness, spectrometer and electrode
gap differ between reports that share these scans; every number taken from the reports, including
the 9–20 dB, therefore rests on one sample.

The reported baseline loss of 30 dB per 0.1 mm at 2 THz is 1.9–2.6× that of pure water (15.8 or
11.6 dB with the broadband or Kindt–Schmuttenmaer model), while nerve tissue absorbs within ~10 % of
water (Berry et al. 2003). Read as a thickness instead, the measured loss of the unbiased slice
relative to the empty cuvette (−12.1/−14.9/−20.1/−29.2 dB at 0.6/1/1.5/2 THz) is that of an
effective-medium white-matter slab on PTFE 198/198/232/314 µm thick with the model's water
(Kindt–Schmuttenmaer) and 191/181/194/232 µm with broadband water (Sim 45). The shared tissue delay
is 0.61 ps relative to the empty cuvette (linear phase with offset over 0.2–1.4 THz, rms deviation
0.11 rad); the same fit applied to the model slab's transmission phase gives d = 216–220 µm. Loss
and delay agree: the slice is optically ~200–230 µm thick below 1.5 THz, about twice either stated
thickness; at 2 THz the model's water needs 314 µm, because it underestimates water loss there.
The offset of the model's delay phase (−0.5 rad) differs from the measured one (+0.22 rad) by
0.8 rad, which we do not explain. This conflicts with the cuvette: H2 gives a 100 µm PTFE cuvette and H3 samples 20 × 20 × 0.1 mm.
Either the gap is about twice the stated one, or the tissue absorbs ~2× as much as water and has
n ≈ 2.8–3.9 (the authors report 3.6–3.9), which water-based tissue cannot have (Berry et al. 2003;
Sec. 2.2). We cannot decide which. The per-slice bounds of Sec. 4.6 are therefore scaled by this
optical thickness (×2.2 at 0.6–1 THz, ×3.1 at 2–3 THz) as a conservative device, not as a measured
thickness, and the 9 dB conductivity of Sec. 5.2 falls to ~130–450 S/m on that scaling.

### 5.6 The published waveforms

The papers give no raw data (the deposited dataset, IEEE DataPort doi:10.21227/h17h-tp94, requires
a subscription and was not available to us), so the published waveforms and spectra were traced
from the PDFs: pages rendered at 300 dpi, time axes calibrated on gridlines to ≤ 0.006 ps, curves
followed by colour with a tracker that crosses labels and occlusions (tested on synthetic images;
every trace overlaid on its source in `results/sim42` and `results/sim43`). The traced H3 pulses
reproduce H3's published spectra within 1.5 dB up to 1.7 THz. Spectra were traced as the upper
envelope of each published curve and binned to 0.05 THz.

**An undelayed component.** Each biased pulse was fitted as a·(no voltage) + b·(empty)(t − τ),
with the shared empty-cuvette pulse (Figs. 6 and 7):

| Bias time (report) | b, time domain | b, other estimates | τ (ps) | residual | tissue term shifted too: shift (ps), residual | best one-path fit |
|---|---|---|---|---|---|---|
| 30 s (H1) | 0.039–0.040 | — | −0.042 | 9 % | −0.002, 9 % | 20 % |
| 2 min (H2) | 0.113–0.120 | — | −0.039 | 25 % | +0.043, 18 % | 63 % |
| 4 min (H1) | 0.166–0.174 | 0.164 (spectral fit with q); 0.20–0.21 (1.5–3 THz level) | −0.039 | 27 % | +0.053, 20 % | 84 % |
| 5 min (H3) | 0.209–0.217 | 0.206–0.212 (spectral fit with q); 0.22–0.27 (1.5–3 THz level) | −0.026 | 28 % | +0.066, 23 % | 77 % |

The time-domain range covers two independent tracings (Sims 39/42 and 43), ±1-pixel shifts
(≤ 0.006) and both fit forms; the spread between methods, up to ~25 % at 4–5 min, is the
uncertainty of b and is several times the tracing noise. The tissue pulse arrives 0.61 ps after the empty
cuvette; the second component sits within 0.05 ps of the empty-cuvette delay and grows steadily with
bias time. The 2-min scan of H2 was traced after the others had been fitted and falls between them
(0.12, against 0.10 from linear interpolation). Two spectra show a component of this size without fitting. In H3's Fig. 1b, biased minus empty
cuvette is flat at −11.3 to −13.2 dB over 1.5–3 THz, a field fraction of 0.22–0.27. In H1's
Fig. 4b, the 4-min spectrum is flat at −13.7 dB relative to the empty cuvette over 1.5–3 THz (0.21),
where the unbiased tissue transmits −20 to −32 dB. These high-band levels exceed the time-domain b
by 15–25 %, and above 1.85 THz the 4-min model (Fig. 8) lies 1.2–2.0 dB below the published
spectrum: the component is larger at high frequency than the pulse fit, dominated by 0.3–1 THz,
gives.

τ is slightly negative in every fit (−0.026 to −0.043 ps, 2–3 pixels), where an air path gives
τ = 0 and any tissue in the path gives τ > 0. A drift of the delay stage would shift both terms. Refitted with one common shift of the whole
biased scan and τ = 0, the residuals rise to 10/29/34/31 %, so drift does not explain the lead.
Letting the tissue term and the component shift independently (a tissue-delay model, not a drift
test), the tissue term moves later by −0.002/+0.043/+0.053/+0.066 ps (±0.013 over ±1-pixel shifts;
0/7/9/11 % of the tissue delay, growing with bias time), residuals fall to 9/18/20/23 %, and the
component stays 0.019–0.044 ps ahead of the empty-cuvette pulse (±0.013). The lead is 1.2–3× the
0.015 ps to which the no-voltage minima of the three reports agree, so we do not interpret it.
A component that crosses up to ~20 % of the tissue's optical thickness fits H3 as
well as one that crosses none (residuals 29–32 % against 28 %), so "undelayed" here means "within
~0.1 ps of the empty-cuvette delay".

**A thicker tissue path.** The residuals of 25–28 % at 2–5 min are several times the spread from
±1-pixel shifts (~5× in H3), so the two-term sum is incomplete. In the frequency domain, with the
measured tissue transfer T = B/K (B no voltage, K empty) raised to a power q, the model
A·K·T^q + b·K·e^{−iωτ} leaves 17 % (H1, 4 min, q = 1.09) and 9–20 % (H3, q = 1.18 to 1.11 for bands
ending at 1.0 to 1.7 THz), with b stable (H3: 0.206–0.212). The tissue path is 9–18 % optically
thicker under bias; the value of q depends on the band used. The time-domain fit with the tissue
term free to shift finds the same thing independently: an extra tissue delay of 0, 7, 9 and 11 %
after 0.5, 2, 4 and 5 min (q = 1.00 at 30 s and 1.09 at 4 min in the frequency domain).

**The 0.6 THz dip.** This is a consistency check, not a prediction: q and b are fitted to the same
traces, and q is what places the dip (H1 at 4 min: 0.71 THz with q = 1, 0.62–0.65 THz with
q = 1.09–1.10). With the fitted values the dip falls at 0.57 THz in H3 (reported 0.56 THz) and at
0.63–0.65 THz in H1 after 4 min (published spectrum 0.65 THz). H3's 0.56 THz coincides with the
strongest water-vapour line (0.557 THz; van Exter et al. 1989), which the unpurged spectra show, so
its exact frequency does not discriminate. Dip depth is ill-conditioned: it
depends on how exactly the two terms cancel. A dip is deep where the terms are in antiphase and
equal in size. The tissue delay puts the first antiphase at 0.71 THz, and the size match there
(|a·T| ≈ 0.16 at 0.6–0.7 THz) is reached as b grows to ~0.17, so the 0.6–0.7 THz dip deepens with
bias time. At 30 s the traced ratio to the unbiased spectrum is a shallow trough, deepest
(−4.5 dB) at 0.72 THz, where q = 1 puts the first antiphase. In these data the reported drift from
0.7 to 0.5 THz is a dip that deepens and moves from 0.72 THz (30 s) to 0.62–0.65 THz (4 min) as b
and q grow.

**The 2 THz dip is not reproduced.** The published 30-s spectrum of H1 has a narrow dip at
2.00 THz, 43.5 dB below the empty cuvette and 14 dB below the unbiased spectrum. Traced per pixel
column (Sim 44), it is centred at 1.99 THz, 15 dB below the curve's own level within ±0.05 THz, and
its bottom (−71 dB absolute) is 2–5 dB above the noise (the unbiased sample spectrum is flat noise at
−73 to −76 dB above 3.6 THz; H1, Fig. 3b), so its depth is a lower bound set by dynamic range
(Jepsen and Fischer 2005; Withayachumnankul and Naftaly 2014). In the binned spectrum the 30-s curve
lies 3–14 dB below the unbiased one over 1.85–2.10 THz. The published curves share 17–19 narrow
water-vapour dips over 1.5–2.9 THz (Pickett et al. 1998); dips that appear in only one curve are
≤ 7 dB deep, except this one. It is therefore a real feature of the 30-s scan, of unknown origin;
with one scan, whether it comes from the sample cannot be decided. We still tested whether the two terms
produce it. Taking |T|
from the published unbiased spectrum, the phase of T from the traced waveforms (linear, rms
0.11 rad over 0.2–1.4 THz, extrapolated), and a, b, τ from the waveform fits, nothing is fitted to
the biased spectra (q = 1 at 30 s; q = 1.09 at 4 min, fitted to the waveform). The model follows the
band means within 0.4 dB at 0.5–0.8 THz, but over 1.85–2.15 THz it is 5.5 dB above the published
30-s spectrum and over 2.3–3 THz 3.6 dB below it (Fig. 8). Its second antiphase falls at 2.24 THz,
where at 30 s the two terms are nearly equal (b = 0.04, |a·T| ≈ 0.03), giving a dip at 2.25 THz that
is gone by 4 min, as the measured one is. But the measured and model dips are separate features: at
2.00 THz the model is ~13 dB above the data (about −30 against −43.5 dB), and near 2.25 THz the
published spectrum has only a weaker minimum (−32 dB at 2.20 THz). The model's antiphase rests on
extrapolating the tissue phase beyond the traced band: fitting it over 0.2–1.0 to 0.2–1.6 THz moves
the 30-s antiphase over 2.10–2.32 THz (arg T at 2 THz spread 0.8 rad), still 0.11–0.33 THz above the
measured dip, and the published 1.95–2.10 THz level is only 3.2 dB above the deepest the two terms
can reach with any phase (Sim 43). The two terms therefore cannot be excluded as its source; they do
not place it at 1.99 THz with the traced phase. Any two terms delayed by 0.61 ps reach antiphase near three times
the first-antiphase frequency, and the depth needs only b ≈ |a·T|, so a transient dip in this
region is generic. We count the 2 THz feature as not explained, but not excluded.

**Not tested.** H1's statement that nothing changes after ~6 min, and the recovery over ~1 h, cannot
be tested without later traces. A gas column inside the tissue was fitted and withdrawn: on the
traced points it fitted poorly (χ² = 32–52 on 8 points; REFEREE item 40).

**What the component could be.** It is beam that reaches the detector without crossing the
tissue. A gap between slice and window is in series with the tissue and cannot do this, and a slice
shrinking at its edges cannot clear the beam (H1: beam ~5 mm across; H3: 20 × 20 mm samples;
electrodes 15–18 mm apart). With mode-matched detection the component's field fraction equals the
fraction of beam power through the opening (Sim 45). For a Gaussian beam 5 mm across (1/e² or
FWHM diameter), b = 0.17–0.22 then needs an opening on the beam axis 1.5–3.0 mm across, through the
full thickness (0.7–1.2 mm at 30 s): a hole or tear from drying, or a gas bubble. A 2 mm bubble in a
0.1 mm layer holds 0.3 mm³, against 25–42 mm³ of electrolytic gas in 5 min (Sec. 5.4). If the beam
waist scales with wavelength, one such opening sized at 0.6 THz would pass 0.87–0.94 of the field at
2 THz; the measured high-band fraction is 0.20–0.27 and flat over 1.5–3 THz. A single on-axis
opening is therefore disfavoured; openings spread across the beam (cracks, partial dewetting), whose
area fraction does not depend on frequency, fit the flat level better. The tissue change needs more tissue of the same kind in the beam
(Sec. 5.4; Fig. 9): 18–37 µm more of a ~220 µm optical path. A slice swelling only in thickness
(water drawn in under the bias, e.g. by electro-osmotic flow) or tissue moving into the beam would do
it, but the slice is held in a 100 µm cuvette and the Joule heat favours evaporation (Sec. 5.4), so
we cannot show any such cause to be plausible. The data cannot identify the cause. If the polarisation and DC-direction controls were done
as described, they argue against any such cause, and they remain unexplained. They are described in
words only (H1, p. 18; H3), with no waveforms or spectra, so they cannot be tested on the published
data.

### 5.7 What matches and what does not

| Reported | This work |
|---|---|
| Dip near 0.6 THz (9–20 dB) | No node or fibre resonance in 0.1–3 THz. The undelayed component plus a 9–18 % thicker tissue path is consistent with it (0.57–0.65 THz against 0.56–0.65 measured); q is fitted. The longer path needs more tissue in the beam; no change of the tissue's permittivity (heating, ions, conduction, water content) or node state gives it |
| Drift 0.7 → 0.5 THz | The dip deepens and moves from ~0.7 to 0.62–0.65 THz as b and q grow; water-vapour lines cannot drift. No traced spectrum has a dip below 0.56 THz, so the 0.5 THz endpoint is not in the published figures |
| 2 THz dip (30 s) | Not reproduced, not excluded. A real feature of one scan (15 dB, twice any other single-curve dip), its depth limited by dynamic range; the two terms reach that depth but, with the traced phase extrapolated, at 2.10–2.32 THz rather than 1.99 THz |
| Only E∥ fibres, only DC across fibres | Not explained; bulk heating and drying would act for both DC directions |
| Onset ~5 min, recovery ~1 h | Component grows 4 → 12 → 17 → 22 % (0.5, 2, 4, 5 min); recovery not traced |
| No significant change of average power | Traced energy −1.2 dB (H3), consistent |

Decisive measurements are cheap: logging current and slice mass (18–37 µm more water-based tissue
in the beam is ~6–11 mg over 20 × 20 mm), measuring the cuvette gap and the slice thickness during bias, an N₂ purge, a reference through a
heated or dried slice without bias, a measurement with the slice clamped between two windows, and
release of the raw data.

## 6. Discussion and limits

- **Linear, passive response.** DC-driven structural changes of fibres (swelling, electro-osmotic
  water flow, pH fronts) are tested only through their bulk effect (heating, water loss, gas).
  Active gain at the node (Liu Y et al. 2020) is outside a passive, energy-conserving calculation.
- **Single fibre plus effective medium.** The slice bounds rest on the induced dipole of one fibre
  and a quasi-static local field in dense packing. Finite fibre segments in a 70–100 µm slice lower
  the E∥ response, so the infinite-fibre values are upper bounds.
- **Band.** 0.1–3 THz. Song et al. (2026) work at 10–40 THz; mid-IR guiding (Liu G et al. 2019) is
  not tested.
- **No raw experimental data.** The waveform analysis rests on traced published figures from one
  measurement session; its fits leave 9–28 % (18–23 % with the tissue term free to shift; several times the tracing noise), the tissue-path
  factor q depends on the fitted band, and the cause of the undelayed component is unknown.
- **Slice thickness.** The measured loss and delay give an optical thickness of ~220 µm against
  the stated 100 µm cuvette (Sec. 5.5). We cannot tell a wider gap from tissue that is not
  water-based; the bounds use the optical thickness as a conservative scaling.
- **Conductivity bound** was not rerun with ±20 % materials; with the passive-ε material factor
  (×1.45) and every extreme combined it is 1.2 dB, 7× below 9 dB.

The result is a negative one for the resonance hypothesis and a positive one for dosimetry: node
absorption in this band is fixed by water loss, node volume and a geometric factor that the
law gives within a few percent, with no free parameter fitted to full wave.

## 7. Conclusion

In 0.1–3 THz the node of Ranvier is a quasi-static water gap. Its absorbed power follows
P_node = ½ωε0ε″_w|E0|²V_nF with F from an interpolation law, exact in its limits, in node length over myelin thickness, g-ratio
and permittivity contrast, within 2.3 % of full wave for fibres up to 10 µm. For the THz field along the fibres, as in every
published scan, no node state changes a 100 µm slice by more than 1.1 dB at internodes ≥ 200 µm
and packing ≤ 0.6; with every extreme combined (50 µm internodes, packing 0.8, the ~220 µm optical
thickness of the measured slice) node states reach 1.5 dB at 0.6 THz but give at most 0.06 rad of
the ≥ 0.23 rad of extra delay the 0.6 THz feature needs, and no state gives loss and delay together. The conducting-fibre model proposed for the Yerevan
experiments gives the required loss at 23–46 S/m but almost none of the required delay.
The published waveforms, all from one session, instead contain a component at the empty-cuvette
delay that grows over the four biased scans, together with a 9–18 % longer optical path through the tissue, which needs more tissue in the
beam: heating, ions, conduction and water content do not give its loss and delay. These are
consistent with the 0.6 THz dip but do not reproduce the 2 THz dip, a feature of one scan whose
origin cannot be decided from the published data. Their physical cause, and the reported selectivity, need the raw data or the
control measurements listed in Sec. 5.7.

## Data and code availability

Every script, benchmark, result file and run log is in the project repository: solver
`src/thznerve/fem`, materials and laws `src/thznerve/physics`, tracing `src/thznerve/digitise.py`,
runs `experiments/`, results `results/`, figures `paper/figures/`, and the full iteration log with
withdrawn results in `REPORT.md`.

## Figures

1. `fig_model.png` — unit cell; ε(f) of water, axoplasm and myelin; index; absorption length; k_w·Ln.
2. `fig_legacy_artefact.png` — the legacy pipeline artefact (supplementary Fig. S1).
3. `fig_node_spectrum.png` — η, node cross-section and share vs frequency, both polarisations, full
   wave vs quasi-static.
4. `fig_scaling.png` — η vs node length / myelin thickness with the law; full wave vs Laplace
   for all geometries; P_node vs internode length.
5. `fig_tissue_voltage.png` — slice transmission E∥/E⊥; node conductivity (internode 1 mm) and DC
   bulk effects vs the reported 9 dB; the binding comparison, with the change of T, is Fig. 9.
6. `fig_bypass.png` — traced H3 waveforms with the two-term fit; spectra and fits.
7. `fig_paper1_bypass.png` — traced H1 waveforms (30 s and 4 min at 100 V) with the two-term fits.
8. `fig_paper1_spectra.png` — H1's published spectra relative to the empty cuvette and the two-term
   model with a, b, τ (and q = 1.09 at 4 min) from the waveform fits; nothing fitted to the spectra.
9. `fig_tissue_change.png` — extra loss vs extra delay of the slice at 0.6 and 1.0 THz: the change of
   T the fits need (q = 1.09–1.18) against every node state with every extreme, heating (two water models), water gain,
   ions and axoplasm conduction (Sim 46).

## References

- Arancibia-Carcamo IL, Ford MC, Cossell L, Ishida K, Tohyama K, Attwell D (2017) Node of Ranvier
  length as a potential regulator of myelinated axon conduction speed. eLife 6:e23329.
  doi:10.7554/eLife.23329
- Berry E, Fitzgerald AJ, Zinov'ev NN, Walker GC, Homer-Vanniasinkam S, Sudworth CD, Miles RE,
  Chamberlain JM, Smith MA (2003) Optical properties of tissue measured using terahertz pulsed
  imaging. Proc SPIE 5030:459–470. doi:10.1117/12.479993
- Gimsa J, Wachner D (1999) A polarization model overcoming the geometric restrictions of the Laplace
  solution for spheroidal cells: obtaining new equations for field-induced forces and transmembrane
  potential. Biophys J 77(3):1316–1326. doi:10.1016/S0006-3495(99)76981-X
- Guo L, Xu D, Wang K, Sun Y, Zhang Q, Ning H, Lu C, Wang S, Gong Y (2022) Electromagnetic
  characteristics of in vivo nerve fibers at the terahertz-far-infrared band. Front Bioeng
  Biotechnol 10:1055232. doi:10.3389/fbioe.2022.1055232
- (H1) Hovhannisyan T, Makaryan A (2024) Investigation of nerve fiber network in THz spectrum
  range. Armenian J Phys 17(1):17–20. doi:10.54503/18291171-2024.17.1-17
- (H2) Hovhannisyan TN (2024) Nonlinear properties of biological nerve fiber system in THz spectrum
  range. Reports NAS RA 124(3–4). doi:10.54503/0321-1339-2024.124.3-4-7
- (H3) Hovhannisyan T, Makaryan A, Martirosyan R (2024) Investigation of nonlinear properties of
  animal spinal cord in THz range. IRPhE 2024, IET Conf Proc. doi:10.1049/icp.2025.0984
- (H4) Hovhannisyan T, Makaryan A, Martirosyan R (2024) Study of the parameters of the spinal cord
  neural system in the terahertz frequency range. Reports NAS RA 124(2):16–23 (not seen; cited by
  H2 and H3)
- Kindt JT, Schmuttenmaer CA (1996) Far-infrared dielectric properties of polar liquids probed by
  femtosecond terahertz pulse spectroscopy. J Phys Chem 100(24):10373–10379. doi:10.1021/jp960141g
- Kotnik T, Miklavčič D (2000) Theoretical evaluation of the distributed power dissipation in
  biological cells exposed to electric fields. Bioelectromagnetics 21(5):385–394.
  doi:10.1002/1521-186X(200007)21:5<385::AID-BEM7>3.0.CO;2-F
- Kumar S, Boone K, Tuszynski J, Barclay P, Simon C (2016) Possible existence of optical
  communication channels in the brain. Sci Rep 6:36508. doi:10.1038/srep36508
- Jepsen PU, Møller U, Merbold H (2007) Investigation of aqueous alcohol and sugar solutions with
  reflection terahertz time-domain spectroscopy. Opt Express 15(22):14717–14737.
  doi:10.1364/OE.15.014717
- Jepsen PU, Fischer BM (2005) Dynamic range in terahertz time-domain transmission and reflection
  spectroscopy. Opt Lett 30(1):29–31. doi:10.1364/OL.30.000029
- Liebe HJ, Hufford GA, Manabe T (1991) A model for the complex permittivity of water at frequencies
  below 1 THz. Int J Infrared Millim Waves 12(7):659–675. doi:10.1007/BF01008897
- Liu G, Chang C, Qiao Z, Wu K, Zhu Z, Cui G, Peng W, Tang Y, Li J, Fan C (2019) Myelin sheath as a
  dielectric waveguide for signal propagation in the mid-infrared to terahertz spectral range. Adv
  Funct Mater 29(7):1807862. doi:10.1002/adfm.201807862
- Liu Y, Wu K, Liu C, Cui G, Chang C, Liu G (2020) Amplification of terahertz/infrared field at the
  nodes of Ranvier for myelinated nerve. Sci China Phys Mech Astron 63(7):274211.
  doi:10.1007/s11433-019-1530-2
- McIntyre CC, Richardson AG, Grill WM (2002) Modeling the excitability of mammalian nerve fibers:
  influence of afterpotentials on the recovery cycle. J Neurophysiol 87(2):995–1006.
  doi:10.1152/jn.00353.2001
- Millon C, Schmidt J, Ramos S, van Dam EP, Buchmann A, Saraceno C, Novelli F (2022)
  Temperature-independent non-linear terahertz transmission by liquid water. AIP Adv
  12(11):115319. doi:10.1063/5.0120417
- Pickett HM, Poynter RL, Cohen EA, Delitsky ML, Pearson JC, Müller HSP (1998) Submillimeter,
  millimeter, and microwave spectral line catalog. J Quant Spectrosc Radiat Transf 60(5):883–890.
  doi:10.1016/S0022-4073(98)00091-0
- Ranck JB Jr, BeMent SL (1965) The specific impedance of the dorsal columns of cat: an anisotropic
  medium. Exp Neurol 11(4):451–463. doi:10.1016/0014-4886(65)90059-2
- Romanenko S, Begley R, Harvey AR, Hool L, Wallace VP (2017) The interaction between electromagnetic
  fields at megahertz, gigahertz and terahertz frequencies with cells, tissues and organisms: risks
  and potential. J R Soc Interface 14:20170585. doi:10.1098/rsif.2017.0585
- Song M, Xu Z, Xu J (2026) Simulation of neural signals' spectrum and transmission from the
  perspective of impedance-matching effect of the Ranvier node. Bioengineering 13(8):955.
  doi:10.3390/bioengineering13080955
- van Exter M, Fattinger C, Grischkowsky D (1989) Terahertz time-domain spectroscopy of water
  vapor. Opt Lett 14(20):1128–1130. doi:10.1364/OL.14.001128
- Wilmink GJ, Grundt JE (2011) Invited review article: current state of research on biological
  effects of terahertz radiation. J Infrared Millim Terahertz Waves 32:1074–1122.
  doi:10.1007/s10762-011-9794-5
- Withayachumnankul W, Naftaly M (2014) Fundamentals of measurement in terahertz time-domain
  spectroscopy. J Infrared Millim Terahertz Waves 35:610–637. doi:10.1007/s10762-013-0042-z
- Xu J, Plaxco KW, Allen SJ, Bjarnason JE, Brown ER (2007) 0.15–3.72 THz absorption of aqueous salts
  and saline solutions. Appl Phys Lett 90:031908. doi:10.1063/1.2430781
- Xiang Z, Tang C, Yan L, Chang C, Liu G (2020) A new viewpoint and model of neural signal generation
  and transmission: signal transmission on myelinated neuron. Chin Phys B 29(10):108701.
  doi:10.1088/1674-1056/abaee9
- Xiang Z, Tang C, Chang C, Liu G (2021) A new viewpoint and model of neural signal generation and
  transmission: signal transmission on unmyelinated neurons. Nano Res 14(3):590–600.
  doi:10.1007/s12274-020-3016-1
- Zhang Z, Li Y, Xiang Z, Huang Y, Wang R, Chang C (2022) Dielectric dispersion characteristics of
  the phospholipid bilayer with subnanometer resolution from terahertz to mid-infrared. Front Bioeng
  Biotechnol 10:984880. doi:10.3389/fbioe.2022.984880
- Zangari A, Micheli D, Galeazzi R, Tozzi A (2018) Node of Ranvier as an array of bio-nanoantennas
  for infrared communication in nerve tissue. Sci Rep 8:539. doi:10.1038/s41598-017-18866-x
