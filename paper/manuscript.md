# No terahertz resonance at the node of Ranvier: a quasi-static absorption law and a re-analysis of published spinal-cord waveforms

Aram Poghossian

Target journal: *Journal of Infrared, Millimeter, and Terahertz Waves* (Springer).
Figures are in `paper/figures/`; every number below traces to a script and a result file listed in
`REPORT.md` (Sims 27–77) and to `results/run_log.csv`.

## Abstract

DC-biased spinal-cord slices reportedly show a 9–20 dB terahertz transmission dip drifting from 0.7
to 0.5 THz, attributed to resonant nodes of Ranvier or DC-opened channels. We test this with a
full-wave model of a myelinated fibre in water, benchmarked to 2 × 10⁻⁴. First, in 0.1–3 THz the
node does not resonate: its absorbed power follows a three-constant quasi-static law with exact
limits, within 2.3 % of full wave for fibres up to 10 µm. Second, for the field along the fibres, as
in every published scan, no node state (ε″ ≥ 0.1) changes a 100 µm slice by more than 1.1 dB at
internodes ≥ 200 µm and packing ≤ 0.6. Third, the biased waveforms of one sample (three reports) fit approximately as the unbiased pulse
whose tissue part lags by 0.16–0.28 rad at 0.6 THz and loses more at high frequency,
plus a second pulse near the empty-cuvette delay growing to 21 %. Without that change, fits miss by
2.2–4.0 times what fitting noise and tracing allow. It fits the 0.6 THz dip, not a 2.0 THz dip. A second sample (a fourth report) transmits more and earlier, with no second pulse. Extra bulk
conduction cannot alone produce the change. The waveforms do not identify its cause: uneven or
thicker tissue, gas, warming and node states each fit under some reading, the node states only as
nearly lossless negative permittivities or 1000 S/m conductors, which an aqueous gap cannot hold. The
reported polarisation selectivity stays unexplained.

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

Experiments now report a THz effect in this tissue. Hovhannisyan and co-workers (four reports,
2024, refs. H1–H4) measured 70–100 µm slices of pig and sheep thoracic spinal-cord white
matter by THz time-domain spectroscopy (TDS) with 100 V DC across 15–18 mm. They report an absorption
dip of 9–20 dB near 0.6 THz that drifts from 0.7 to 0.5 THz over minutes, changes near 2 THz,
onset in ~5 min and recovery in ~1 h, and an effect that appears only when the THz field is parallel
to the fibres and the DC field across them. They explain it by DC-opened sodium channels at the
nodes, which would turn fibres into nonlinear conducting wires and the slice into a diffraction
grating.

Below 3 THz water dominates: its permittivity is 4.1 − 2.1j at 1 THz, and its power absorption
length is 35–118 µm. Whether a node can resonate inside such a host, and whether any node or fibre
state could change a 100 µm slice by 9 dB, are questions a full-wave calculation can answer. The resonance models above work at 10 THz and above,
so a null result below 3 THz does not contradict them; it bears on the node reading of H1–H3 and on
dosimetry. The two limits of the absorption law below are textbook electrostatics; what is new is a
closed-form interpolation between them, checked against full wave over realistic fibre geometry. We
ask:

1. Does a single myelinated fibre with realistic dimensions, in water with measured dispersion,
   absorb resonantly at the node in 0.1–3 THz?
2. If not, what law sets the absorbed power?
3. Can any node state, or the fibre state the experimenters propose, change the transmission of a
   100 µm slice by the reported amount?
4. What do the published waveforms themselves show?

The answers are: no resonance; a quasi-static law, exact in its limits; no for the THz field
along the fibres, the configuration of every published scan; and a second pulse within ~0.1 ps of the
empty-cuvette delay plus a tissue pulse that loses more at high frequency, consistent with the 0.6 THz dip but not with the 2 THz
dip. The change of the tissue pulse, a loss rising with frequency and a lag, is robust to the
fitting assumptions. Judged against fitting noise and tracing, the waveforms exclude extra conduction in the tissue water or axoplasm as its sole cause but do not identify the cause: a thicker or unevenly thick slice, gas
in the slice, warming and node states each fit under some reading. Node states are excluded instead by
the forward bounds and because an aqueous gap cannot hold the states that fit. A fourth report (H4),
traced here, adds one new sample, whose tissue pulse after bias is stronger and earlier, the
opposite of the change in H1–H3, with no second pulse; the resonance it reports for the sample of H1–H3 is consistent with the cancellation between the
two pulses, on a scan that matches H3's. Section 2 describes the model, Section 3 its verification, Section 4 the node results, Section 5 the experimenters' model and the published data, and Section 6 the limits.

**Prior work and what is new.** Quasi-static field factors of cell and membrane compartments are
standard in radio-frequency and microwave microdosimetry (Gimsa and Wachner 1999; Kotnik and
Miklavčič 2000), and both limits of the law below follow from the boundary conditions (continuity of
normal D and of tangential E). What is new is (i) the interpolation between the limits, with three
constants fitted once to a Laplace solver and tested against full wave over node length, myelin
thickness, axon radius and frequency, for the node in water with measured THz dispersion; (ii)
bounds on what any node state can do to a slice, including the negative-ε ring resonance; and (iii)
a re-analysis of the published waveforms, with each candidate cause applied inside the waveform fit. Table 0 lists the closest work.

| Paper | What they did | Band | Water | This work |
|---|---|---|---|---|
| Liu G et al. 2019 | myelin index ~2× surroundings; sheath as dielectric waveguide (SR-FTIR) | mid-IR | measured, mid-IR | below 3 THz myelin's index is below water's and the axon-core V number ≤ 0.4: guided modes leak |
| Liu Y et al. 2020 | node as an active (three-level gain) amplifier | mid-IR–THz | not stated | passive response only; gain is outside scope |
| Guo et al. 2022 | three-layer fibre modes, mode matching at the node, vs CST | 0–30 THz | constant complex ε | double-Debye water, realistic internodes, node power, Mie benchmark, convergence |
| Song et al. 2026 | node as impedance transformer, peaks vs node length | 10–40 THz | lossless, ε = 1.44 | with water's ε the node is quasi-static and non-resonant below 3 THz |
| Zangari et al. 2018 | node as an array of nano-antennas (FEM) | IR | lossless or constant | no antenna resonance below 3 THz; interpolation law with exact limits |
| Kumar et al. 2016 | myelinated axons as optical waveguides | visible–NIR | n/a | the waveguide picture does not carry over to THz |
| H1–H4 (2024) | THz-TDS of biased cord slices; dips at 0.6 and 2 THz; open node channels | 0.2–4 THz | experiment | tests their model and node resonance; re-analyses their waveforms |

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
last retained order contributes ≤ 5 × 10⁻⁴ of the absorbed power in 350 of 433 production runs,
≤ 1.5 × 10⁻³ for the thickest myelin and ≤ 9 × 10⁻³ for the thickest fibre (a = 7.5 µm, E⊥).

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
| B2: 2D layered cylinder, TM and TE | Mie series | absorbed power per layer ≤ 8 × 10⁻⁶; flux balance ≤ 7.5 × 10⁻⁷ |
| B3: axisymmetric model, m = 0…3, both polarisations | per-order series | ≤ 2 × 10⁻⁴; balance ≤ 2 × 10⁻⁵ |
| Null test (fibre = water) | zero scattering | ≤ 1.4 × 10⁻⁷ |
| Energy balance, every production solve | source = absorbed + outflow | ≤ 5 × 10⁻⁵ (≤ 2.3 × 10⁻⁴ in all review runs) |
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
with no fitted parameter reproduces η within 0.1–6.4 % (E∥, largest at 0.1 THz) and 0.5–0.8 % (E⊥) at every frequency.
For E∥ the node field is below the incident one: the node is shielded by the myelin, not fed by it.

### 4.2 A quasi-static law with exact limits

The node is a sub-wavelength water gap (k_w·Ln ≤ 0.6 for Ln ≤ 5 µm at 3 THz; < 0.15 for Ln = 1 µm)
in a quasi-static field, so

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
b ≤ 5 µm (fibre diameter ≤ 10 µm, including g = 0.875, outside the fit range; the other g-ratio outside it, 0.54, has b = 6.5 µm), 4.3 % at b = 6.5 µm and 9.8 % at b = 10.7 µm and 3 THz. The
internode length changes P_node by < 10⁻⁴ (200–1500 µm).

### 4.4 Why there is no resonance

- The node is short: k_w·Ln ≤ 0.6 at 3 THz even for Ln = 5 µm (0.44 at the longest measured 3.7 µm),
  and the law, which has no retardation, holds there within 2.3 % (Sec. 4.3).
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
ε″ ×0.8), applied together, raise it ×1.45 at 0.6 THz and ×1.43 at 3 THz (0.021 and 0.176 dB; at 3 THz the off-peak points
scatter by ~0.02 dB, so ×1.43 ± 0.16;
the 0.6 THz peak moves to ε′ ≈ −2.0, bracketed by a scan from −2.7 to −1.0). Including this, the
passive bound at internode 200 µm and packing 0.6 is 0.13/0.24/0.63/1.06 dB at 0.6/1/2/3 THz
(×1.45 assumed at 1 and 2 THz). At ε″ = 0.01 the gap's 90° corners are singular and the points do not converge, so the
claim is made for ε″ ≥ 0.1; negative ε′ with ε″ ≪ 1 is not possible in an aqueous gap at THz. Ions
add a Drude term whose loss exceeds its negative ε′ by 1/(ωτ) ≈ 15 at 0.6 THz (collision time τ ≈ 17 fs);
pulling water's ε′ ≈ 4.5 below −1 needs ~3 × 10³ S/m, over 100 times saturated NaCl, and brings
ε″ ≳ 80 with it.

**Bound.** For E∥, internode ≥ 200 µm and packing ≤ 0.6, no node conductivity changes a 100 µm
slice by more than 0.05 dB, and no passive node permittivity by more than 0.1/0.2/0.6/1.1 dB at
0.6/1/2/3 THz (combined ±20 % materials included). The factors that could raise these multiply:
internodes of 50 µm (the smallest CNS fibres; ×4), packing 0.8 (×1.33, assumed linear; at 0.6 the
E∥ response is unchanged, Sec. 4.5) and the optical thickness of the measured slice (Sec. 5.5;
×2.2 at 0.6–1 THz, ×3.1 at 2–3 THz, the larger of the loss and delay values with the model's
water). Taken together ("every extreme" below) they give 1.5/2.7/10/18 dB for passive ε and ≤ 1.2 dB for conductivity
(×1.45 materials assumed). The internode scaling was checked by solving at 50 µm directly: the
slice change matches the value scaled from 200 µm within 0.3 % at every frequency and in both
polarisations.

These bounds are on loss. The published waveforms also need a lag of the tissue pulse (Sec. 5.6);
Sec. 5.7 tests node states inside the waveform fits.

Near 3 THz the reported change is a rise in transmission of up to 17–18 dB; over the full ε grid at
3 THz (33 states) the largest rise any node state gives is 0.003 dB per 100 µm at 1 mm internodes,
≤ 0.4 dB with every extreme. The extremes also need every node, in a slice made only of the smallest
fibres, to hold ε′ between −2 and −1.4 with ε″ = 0.1, which an aqueous gap cannot (above).

These bounds are for the THz field along the fibres, the configuration of every published scan
(H1, p. 18). For E⊥ the conductivity bound with all extremes (materials ×1.45) would be 3.6/13.9/21.7 dB at
0.6/2/3 THz, so the E∥ margin does not carry over to the reported E⊥ controls, which were described
without data (Sec. 5.7).

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

Sec. 5.7 tests conduction inside the published waveform fits.

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
does not explain the DC-direction selectivity. Sec. 5.7 tests warming inside the published waveform fits. The tissue is post-mortem, and the bias puts ~70 mV
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
and delay agree: the slice is optically ~200–230 µm thick below 1.5 THz (181–194 µm in loss with
broadband water), about twice either stated
thickness; at 2 THz the model's water needs 314 µm, because it underestimates water loss there.
The offset of the model's delay phase (−0.5 rad) differs from the measured one (+0.22 rad) by
0.8 rad, which we do not explain. This conflicts with the cuvette: H2 gives a 100 µm PTFE cuvette and H3 samples 20 × 20 × 0.1 mm.
Either the gap is about twice the stated one, or the tissue absorbs ~2× as much as water and has
n ≈ 2.8–3.9 (the authors report 3.6–3.9), which water-based tissue cannot have (Berry et al. 2003;
Sec. 2.2). We cannot decide which. The per-slice bounds of Sec. 4.6 are therefore scaled by this
optical thickness (×2.2 at 0.6–1 THz, ×3.1 at 2–3 THz) as a conservative device, not as a measured
thickness, and the 9 dB conductivity of Sec. 5.2 falls to ~130–450 S/m on that scaling.

### 5.6 The published waveforms

**Data.** The papers give no raw data (the deposited dataset, IEEE DataPort doi:10.21227/h17h-tp94,
requires a subscription and was not available to us), so the published waveforms and spectra were
traced from the PDFs: pages rendered at 300 dpi from embedded bitmaps of 213–244 ppi (one native
pixel is 0.016 ps in H1 and 0.012 ps in H3; the ±1-pixel ranges below use rendered pixels and are
therefore narrower than one native pixel), time axes calibrated on gridlines to ≤ 0.006 ps, curves
followed by colour with a tracker that crosses labels and occlusions (tested on synthetic images;
every trace overlaid on its source in `results/sim42` and `results/sim43`). The traced H3 pulses
reproduce H3's published spectra within 1.5 dB up to 1.7 THz. Spectra were traced as the upper
envelope of each published curve and binned to 0.05 THz.

**Model.** Each biased pulse is fitted as

  biased(t) = a·B(t − s) + b·K(t − τ),

where B, the unbiased (no-voltage) pulse, is the *tissue pulse* and K is the empty-cuvette pulse.
The second term, the *second pulse*, arrives near the empty-cuvette delay, as beam that has not
crossed the tissue would (it arrives slightly earlier; see below). The tissue pulse arrives 0.61 ps after the empty cuvette. In the
frequency domain the same model is A·K·T^q + b·K·e^{−iωτ}, with T = B/K the measured tissue transfer;
a later tissue pulse (s > 0) and a tissue path optically longer by a factor q are equivalent here,
q − 1 ≈ s/0.611 ps. Fits (Figs. 6 and 7):

| Bias time (report) | b, time domain | b, other estimates | τ (ps), s = 0 | residual, s = 0 | s (ps), residual | best one-path fit |
|---|---|---|---|---|---|---|
| 30 s (H1) | 0.034–0.044 | — | −0.042 | 9 % | −0.002, 9 % | 20 % |
| 2 min (H2) | 0.113–0.120 | — | −0.039 | 25 % | +0.043, 18 % | 63 % |
| 4 min (H1) | 0.166–0.174 | 0.164 (spectral fit with q); 0.20–0.21 (1.5–3 THz level) | −0.039 | 27 % | +0.053, 20 % | 84 % |
| 5 min (H3) | 0.209–0.217 | 0.206–0.212 (spectral fit with q); 0.22–0.27 (1.5–3 THz level) | −0.026 | 28 % | +0.066, 23 % | 77 % |

**The second pulse.** b grows with bias time, from 0.04 at 30 s to 0.21 at 5 min. The time-domain
range covers two independent tracings (Sims 39/42 and 43), ±1-pixel shifts (≤ 0.006) and both fit
forms; the spread between methods, up to ~25 % at 4–5 min, is the uncertainty of b and is several
times the tracing noise. The 2-min scan of H2 was traced after the others had been fitted and falls
between them (0.12, against 0.10 from linear interpolation). Two spectra show a second pulse of this
size without fitting. In H3's Fig. 1b, biased minus empty cuvette is flat at −11.3 to −13.2 dB over
1.5–3 THz, a field fraction of 0.22–0.27. In H1's Fig. 4b, the 4-min spectrum is flat at −13.7 dB
relative to the empty cuvette over 1.5–3 THz (0.21), where the unbiased tissue transmits −20 to
−32 dB. These levels exceed the time-domain b by 15–25 %, and above 1.85 THz the 4-min model
(Fig. 8) lies 1.2–2.0 dB below the published spectrum: the second pulse is larger at high frequency
than the pulse fit, dominated by 0.3–1 THz, gives. With s = 0, τ is slightly negative in every
central fit (−0.026 to −0.042 ps, up to 2.6 native pixels; at 30 s the ±1-pixel range reaches
+0.001 ps), where an air path gives τ = 0 and any tissue in the path gives τ > 0. A second pulse that
crosses up to ~20 % of the tissue's optical thickness fits H3 as well as one that crosses none
(residuals 29–32 % against 28 %), so "at the empty-cuvette delay" means within ~0.1 ps of it.

**The change of the tissue pulse.** With s = 0 and no loss term the residuals at 2–5 min (25–28 %)
are several times the spread from ±1-pixel shifts (~5× in H3), so the tissue pulse itself changes
under bias. A shift alone (the tissue pulse later by −0.002/+0.043/+0.053/+0.066 ps at 0.5/2/4/5 min;
residuals 9/18/20/23 %) leaves a misfit that a smooth extra loss rising with frequency removes. That
loss belongs to the tissue pulse: put on the whole biased trace, as a smoothing would act, it fits
worse (Sim 50). The 30 s scan, traced before the reported onset, already carries ~2 dB per THz of it.
How the change splits into a loss and a shift depends on what the traces cannot fix: whether the loss
is causal, how it continues beyond the traced band (up to ~1.7 THz), and whether the 30 s part is a
tracing artefact (zero-phase) or a property of the slice. Sim 55 covers these with 36 reference fits:
a zero-phase floor of 0–2.5 dB/THz, a causal (minimum-phase) increment held constant above 1.5, 2 or
4 THz, a free shift, and a ≤ 1 or a + b ≤ 1. In the 24 with floors of 1.0–2.5 dB/THz, which absorb
most of the 30 s change, the loss rises above its 30 s value by 1.7–2.8, 2.4–4.2 and
4.7–7.1 dB/THz at 2, 4 and 5 min (4.4–7.1 at 5 min over all 36), and the tissue pulse lags at 0.6 THz
by 0.16–0.28 rad. The loss and the lag trade off: a candidate that supplies more delay needs less loss
(Sec. 5.7), so the lag is the robust quantity (−0.02 to 0.00 rad at 30 s); the residuals are 7/10/13/14 % and N_eff = 13–20.

**Effect size, not significance.** The misfit is systematic rather than noise (residual correlated
over 0.3–0.5 ps, N_eff = 13–20 per trace), and significance tests on traced figures proved unusable.
Four calibrations of the F tests on their own null were tried and withdrawn. Sim 55b drew its
synthetic scans from the full fit, which still holds the change. Sims 57 and 60 fitted the synthetic
scans from fewer starting points than the data. In Sims 61 and 63 a block bootstrap doubled N_eff in
the synthetic scans, so the nominal 5 % tests rejected 3–23 % of them (5th–95th percentile). Sim 64 used a null that
keeps the residual's amplitude spectrum (phase-randomised), with the candidate's size fitted inside
the statistic. That null was then tested on pairs with no tissue change: one traced unbiased pulse
fitted with another from a different report. It rejected "no change" in 10 of 12 such comparisons
(Sim 65). Every model leaves a systematic misfit, so with enough points any two traces differ
significantly. Fits are therefore judged by effect size, against what fitting noise and tracing
give. D is the rms of what a restricted fit leaves unexplained, relative to the full fit (the causal
loss increment and a free shift), in % of the trace rms; for linear models it is the rms difference
between the two fitted waveforms. Two things make D > 0 without any physical difference. Fitting
noise: the full fit has more freedom, and with N_eff = 13–20 it absorbs part of the misfit. For each
restricted fit, phase-randomised copies of the full fit's residual (same autocorrelation) are added to
the restricted fit's curve and both models are refitted, 100 times; the candidate's size is held at its
value on the data, which favours the candidate. The full fit absorbs part of any noise added to it, so
its residual underestimates the noise; the copies are scaled up by 1.07–1.24, the factor at which
noise injected into the full fit's own curve leaves a residual as large as on the data (10 draws per
case; the formula sqrt(N_eff/(N_eff − 5)) gives 1.16–1.29) (Sim 69). Tracing: one traced unbiased
pulse fitted with another from a different report, built like each reference, gives D up to
2.6–4.3 % without a shift and ≤ 0.08 % with one (Sim 68). In each scan the threshold adds the
98.3rd percentile of the noise values (1 − 0.05/3) and the tracing value in quadrature. A restricted
fit is excluded in a reference if its D exceeds the threshold at one of the scans at 2, 4 and 5 min,
so a true cause is excluded by chance in about 5 % of references; counted on the noise draws
themselves, 0–6 %. A verdict is called borderline if it flips in more than 10 % of bootstrap resamples
of the 100 draws or its ratio lies within 0.87–1.15, the ~11 % by which one threshold moves when
the draws are repeated with other seeds; "firmly excluded" means excluded and not borderline. The
noise model spreads the misfit evenly in time, while at 4–5 min it is concentrated in part of the
window; noise shaped to the residual's envelope moves single ratios in both directions, a further
reason not to read ratios near 1 as firm. D compares models that are not nested (a candidate's computed change against a
fitted loss increment), and each candidate's size is chosen per scan on the data; both favour the
candidates. Over the eight references of Sec. 5.7 (floors 1.0–2.5 dB/THz, increment held above
2 THz, a ≤ 1 or a + b ≤ 1), no change leaves 17.9–23.0 % at 2–5 min, up to 2.2–3.3 times its
threshold in every reference, none borderline, and 1.8–4.0 % at 30 s, within it (0.3–0.8 times). A
shift alone leaves 7.7–16.8 %, up to 1.9–3.2 times its threshold, so the tissue pulse changes in
shape, not only in time; at 30 s a shift alone is within its threshold in seven of eight references
and at 1.03 times it in the eighth, as expected by chance. These eight references read the 30 s
loss as an artefact of the reference (their floors of 1.0–2.5 dB/THz absorb it), so in them the 30 s
scan shows no change by construction and is not a control. The tracing pairs argue against that
reading (unbiased tracings differ by ≤ 0.12 dB/THz, below). With floors of 0 and 0.5 dB/THz, which read
the 30 s loss as part of the slice (four references, Sim 70), no change leaves 5.0–6.3 % at 30 s,
1.1–1.5 times its threshold, and a shift alone 1.4–1.9 times, firmly (Sim 71): the slice has then
changed by 30 s. A slice whose thickness varies across the beam by an rms of 8–12 % fits the 30 s scan within its threshold in all four with a drift (0.24–0.71 times; Sim 72), so the waveforms do not show
what this early change is. With these floors no change leaves 21.6–24.2 % at 2–5 min, 3.4–4.0 times
the threshold. The decomposition is approximate: the
full fit leaves 10.0–15.2 % of the trace rms at 2–5 min (6.5–7.1 % at 30 s), against 1.8–9.6 %
between unbiased tracings, so part of each biased trace is captured by no model here. One traced
unbiased pulse fitted with another gives a loss slope ≤ 0.12 dB/THz and a lag ≤ 0.03 rad with a ≤ 1
(up to 0.58 dB/THz with a free), against a loss increment of 3.1–7.4 dB/THz over the floor and
0.17–0.26 rad at 2–5 min in the eight references (Sim 57). Whether the tissue pulse also shifts
beyond the delay of its causal loss is not decided. The curves are traced from the same figure on the
same axes, so the axis calibration cancels; the scans share the empty-cuvette pulse and one session.

The lag is also found without any loss term: 0.16–0.56 rad over model forms (a pure shift, or a shift
plus a frequency-independent phase) and fitting windows (Sim 48). The frequency-domain fit finds the
same independently: q = 1.00 at 30 s, 1.09 at 4 min (H1, residual 17 %) and 1.18 to 1.11 at 5 min (H3,
bands ending at 1.0 to 1.7 THz; residuals 9–20 %), with b stable (H3: 0.206–0.212); q depends on the
band used.

The second pulse leads the empty-cuvette pulse in every fit: by 0.019–0.044 ps with a shift
and no loss term (±0.013), by 0.023–0.044 ps at 2–5 min and 0.074–0.081 ps at 30 s with the loss
term (Sim 55; 0.027–0.044 and 0.072–0.084 ps in the eight references of Sec. 5.7, Sim 68), 1.5–5× the 0.015 ps to which the no-voltage minima of the three reports agree. It also leads
in the best fit of every candidate of Sec. 5.7, alone and with a drift (Sim 68): by 0.020–0.069 ps at
2–5 min. The lead varies more between scans than between candidates, so it does not single out a
cause. A shift of the
tissue pulse also reads as a drift of the delay stage by s between the shared reference scans and each biased
scan (one common shift of the whole biased scan with τ = 0 raises the residuals to 10/29/34/31 %, but
with both shifts free it is the same fit, the second pulse then leading by 0.07–0.10 ps); a drift
alone cannot make the tissue pulse weaker at high frequency. A real lead would need a path faster than the empty air gap. It is comparable to the difference
between the empty-cuvette pulses of two reports: H1's and H4's differ by 9 % rms and 0.025 ps, and
with H4's as K the fits without a loss term put the second pulse 0.000–0.015 ps early alone (one pixel
is 0.015 ps) and 0.015 ps early to 0.011 ps late with a drift, b changing by ≤ 0.004 (Sim 75). In the
fits with the causal loss term (no floor, a ≤ 1) it is 0.008 ps early to 0.009 ps late at 2–5 min with
H4's pulse, against 0.020–0.039 ps early with H1's; at 30 s, where b is 0.04, 0.011–0.033 ps early
against 0.039–0.064 ps (Sim 76). But H1's pulse fits the H1–H3 scans better in all four (full-fit
residual 6.3/9.9/11.9/13.5 % against 6.7/10.2/13.1/14.7 %; Sim 77), and H1's and H3's empty-cuvette
minima agree to 0.005 ps. Whether the lead is real is therefore not decided by the traced figures. The
candidate fits of Sec. 5.7, whose lead reaches 0.069 ps with H1's pulse, were not refitted with H4's. Sec. 5.7 tests what could cause the
change.

**The 0.6 THz dip.** This is a consistency check, not a prediction: q and b are fitted to the same
traces, and q is what places the dip (H1 at 4 min: 0.71 THz with q = 1, 0.62–0.65 THz with
q = 1.09–1.10). With the fitted values the dip falls at 0.57 THz in H3 (reported 0.56 THz) and at
0.63–0.65 THz in H1 after 4 min (published spectrum 0.65 THz). H3's 0.56 THz coincides with the
strongest water-vapour line (0.557 THz; van Exter et al. 1989), which the unpurged spectra show, so
its exact frequency does not discriminate. A dip is deep where the two pulses are in antiphase and
equal in size, so its depth depends on how exactly they cancel and is not a robust observable. The
tissue delay puts the first antiphase at 0.71 THz, and the size match there (|a·T| ≈ 0.16 at
0.6–0.7 THz) is reached as b grows to ~0.17, so the dip deepens with bias time. At 30 s the traced
ratio to the unbiased spectrum is a shallow trough, deepest (−4.5 dB) at 0.72 THz, where q = 1 puts
the first antiphase. In these data the reported drift from 0.7 to 0.5 THz is a dip that deepens and
moves from 0.72 THz (30 s) to 0.62–0.65 THz (4 min) as b and q grow.

**The 2 THz dip.** The published 30-s spectrum of H1 has a narrow dip at 2.00 THz, 43.5 dB below the
empty cuvette and 14 dB below the unbiased spectrum. Traced per pixel column (Sim 44), it is centred
at 1.99 THz, 15 dB below the curve's own level within ±0.05 THz, and its bottom (−71 dB absolute) is
2–5 dB above the noise (the unbiased sample spectrum is flat noise at −73 to −76 dB above 3.6 THz;
H1, Fig. 3b), so its depth is a lower bound set by dynamic range (Jepsen and Fischer 2005;
Withayachumnankul and Naftaly 2014). In the binned spectrum the 30-s curve lies 3–14 dB below the
unbiased one over 1.85–2.10 THz. The published curves share 17–19 narrow water-vapour dips over
1.5–2.9 THz (Pickett et al. 1998); dips that appear in only one curve are ≤ 7 dB deep, except this
one. It is therefore a real feature of the 30-s scan; with one scan, whether it comes from the
sample cannot be decided. H1's own arrow covers the 30-s curve at 1.96–1.99 THz (Sim 44), so the
dip's low-frequency flank is hidden and its centre is uncertain by up to ~0.03 THz.

The two-pulse model was tested against it with nothing fitted to the biased spectra: |T| from the
published unbiased spectrum, the phase of T from the traced waveforms (linear, rms 0.11 rad over
0.2–1.4 THz, extrapolated), a, b, τ from the waveform fits, q = 1 at 30 s and 1.09 at 4 min. The
model follows the band means within 0.4 dB at 0.5–0.8 THz, but over 1.85–2.15 THz it is 5.5 dB
above the published 30-s spectrum and over 2.3–3 THz 3.6 dB below it (Fig. 8). Any two pulses
0.61 ps apart reach antiphase again near three times the first-antiphase frequency, and at 30 s the
two are nearly equal there (b = 0.04, |a·T| ≈ 0.03), so the model has a dip at 2.25 THz that is gone
by 4 min, as the measured one is. But at 2.00 THz the model is ~13 dB above the data (about −30
against −43.5 dB), and near 2.25 THz the published spectrum has only a weaker minimum (−32 dB at
2.20 THz). The model's dip rests on extrapolating the tissue phase beyond the traced band: fitting it
over 0.2–1.0 to 0.2–1.6 THz moves the dip over 2.10–2.32 THz (arg T at 2 THz spread 0.8 rad), still
0.11–0.33 THz above the measured one. The published 1.95–2.10 THz band mean is only 3.2 dB above the
deepest the two pulses can reach with any phase; at the 2.00 THz bin the published −43.5 dB is
reached only for b ≤ 0.037, the low end of the traced range (Sim 43). No node state, with every
extreme, moves the model's dip to 1.99 THz at the measured depth (Sec. 5.7). The 2 THz dip is therefore not
reproduced; the two-pulse model is not excluded as its source, because its position there rests on
an extrapolated phase.

**Not tested.** H1's statement that nothing changes after ~6 min, and the recovery over ~1 h, cannot
be tested without later traces. A gas column inside the tissue was fitted and withdrawn: on the
traced points it fitted poorly (χ² = 32–52 on 8 points).

**What the second pulse could be.** Its delay suggests beam that reaches the detector without
crossing the tissue; its lead over the empty-cuvette pulse is comparable to the difference between two reports'
empty-cuvette pulses and is not decided (above). The reading below is a hypothesis. A gap between slice and window is in series with the tissue and cannot do this, and a slice
shrinking at its edges cannot clear the beam (H1: beam ~5 mm across; H3: 20 × 20 mm samples;
electrodes 15–18 mm apart). With mode-matched detection the second pulse's field fraction equals the
fraction of beam power through the opening (Sim 45). For a Gaussian beam 5 mm across (1/e² or
FWHM diameter), b = 0.17–0.22 then needs an opening on the beam axis 1.5–3.0 mm across, through the
full thickness (0.7–1.2 mm at 30 s): a hole or tear from drying, or a gas bubble. A 2 mm bubble in a
0.1 mm layer holds 0.3 mm³, against 25–42 mm³ of electrolytic gas in 5 min (Sec. 5.4). If the beam
waist scales with wavelength, one such opening sized at 0.6 THz would pass 0.87–0.94 of the field at
2 THz; the measured high-band fraction is 0.20–0.27 and flat over 1.5–3 THz. A single on-axis
opening is therefore disfavoured; openings spread across the beam (cracks, partial dewetting), whose
area fraction does not depend on frequency, fit the flat level better.

**A fourth report (H4).** H4 (Hovhannisyan, Makaryan and Martirosyan 2024), found during review,
shows an empty cuvette and an unbiased sample (its Fig. 3a) and two samples before and after 6 min at
100 V (Fig. 4a, 4v), on the −363 to −357 ps axis of H1. Traced like the others (Sim 73; axis misfit
≤ 0.005 ps, one pixel 0.015 ps), the three unbiased pulses peak within 0.05 ps of each other, 0.63–0.68
ps after the empty cuvette, so the figures share one time axis. After the best shift and scale,
sample 1's unbiased pulse matches those of H1–H3 within 3.9–5.0 % rms, close to how well they match each other (1.8–3.3 %), and its "6 min" scan matches H3's 5 min scan within 10.1 %, against 27–47 % between any
other two biased scans (Sim 75). Sample 1 is therefore the sample of H1–H3, and its biased scan is
H3's or one taken soon after it, so it is not independent evidence. After bias its main minimum falls
at −360.845 ps, the empty-cuvette time (−360.850 ps). Fitted with the two-pulse model and H4's own
empty pulse, its second pulse is b = 0.24 within 0.02 ps of the empty-cuvette delay (s + τ) in the fits
without a loss term; the tissue pulse changes too (residual 24 % without a change, 17 % with uneven
thickness, 15 % with the causal loss term). H4 reports a "deep resonance in the low-frequency region" for this sample without giving its
frequency; its spectrum (Fig. 4b) has a broad minimum near 0.55 THz and its refractive index (Fig. 5a)
jumps at 0.60 THz. The spectrum of the traced waveform has its deepest notch at 0.56 THz, and the
fitted two-pulse models put the cancellation at 0.56–0.59 THz. As for H3, this is a consistency check,
not a prediction, made on what is H3's scan or one close to it, and 0.557 THz is a water-vapour line. Sample 2 is new: its unbiased pulse differs from sample 1's by 35–38 %. After bias it has no second pulse and transmits more and earlier: fits that allow a second pulse
give it a small negative amplitude (b = −0.01 to −0.03, at the 0.2 ps bound on its arrival; Sim 76).
Without one, a tissue path 8–16 % thinner fits with a 5.4 % residual with a drift (11.9 % without), a
tissue pulse 1.25 times stronger and 0.05 ps earlier with no change of shape leaves 8.0 %, and no
change 21.2 % (Sim 77). The thinner-path fit needs a drift of +0.067 ps, opposite to the advance
that thinning gives, and without the drift it spans 7.6–16.4 % over the tracing shifts; with no
noise threshold for H4, the thinner path fits about as well as a stronger, earlier pulse and the two
are not ranked. Either way the tissue pulse changes in the opposite direction to H1–H3's. The full fit,
whose loss term cannot add gain, does no better than no change here, so no effect size is formed.
Both samples transmit more after bias (traced: +15 and +6 dB at 2 THz; H4 gives "up to 20 dB" at 2 THz
for sample 1 and, for sample 2, a rise starting at 0.4 THz and 5 dB near 2 THz), which beam crossing
less tissue gives. For sample 2 H4 reports no dip: its refractive index falls over the whole band, and
the authors suppose its resonance lies below the band, while stating that all samples show similar
regularities. Sample 2, the one independent test, thus shows neither the dip nor a second pulse. H4 has one figure per sample and no second tracing of an
unbiased pair, so no noise threshold is built for it.

### 5.7 What could cause the change

**How the candidates were tested.** In the two-pulse model the 0.6 THz dip is a cancellation, so its
depth is not the quantity to bound; what the fits need is the change of the tissue pulse of Sec. 5.6.
Each candidate's computed change of the slice over the whole band (size v; Sec. 5.5 slab of 220 µm,
Kindt–Schmuttenmaer water, or the node states of Sec. 4.6) multiplies the spectrum of the unbiased
pulse and must replace the loss increment and the shift. The fit with the candidate at its best size
on a grid and no loss term, with s = 0 (alone) or s free (with a timing drift of both pulses), is
compared with the fit with the loss increment and s free by its effect size D, against the threshold
of Sec. 5.6 (fitting noise and tracing). References: the eight with a zero-phase floor of 1.0, 1.5,
2.07 or 2.5 dB/THz, a causal increment held above 2 THz, and a ≤ 1 or a + b ≤ 1 (Sim 59), and four
with floors of 0 and 0.5 dB/THz (Sim 70). Sizes are searched over plausible and wide ranges (Sim 66);
Table 2 uses the best over both, except gas, searched only where its dilute model holds (v ≤ 0.2). Node states are all
80 at every extreme (Sec. 4.6), each scaled by λ = 0.5–6 of that change, the best (state, λ) chosen on
the data for each reference, scan and mode (Sim 68). Uneven thickness averages the transfer of a
uniform slab over a Gaussian spread of thickness across the beam (mean change m from −30 to +51 %, rms
σ from 0 to 40 %, cut at 5 % of the thickness), the coherent mean that the detector sees (Sim 72); for
a pure phase φ the mean is e^{iφm − φ²σ²/2}, a loss rising with frequency (tested). This is the mean
over many independent patches, and under a beam ~5 mm across one slice holds few. For random
thickness fields with the fitted sizes and correlation lengths of 0.5–2 mm, one slice departs from the
mean by an rms of 0.7–4.6 dB at 0.6 THz and 2.0–13.8 dB at 1.5 THz, a quarter to three times the mean change (−1.6 to −4.0 and −8 to −13 dB) depending on patch size and
spread (Sim 74). Only patches of ≲ 0.25 mm average out at 0.6 THz (0.4–1.0 dB; at 1.5 THz still 1.2–6.9 dB),
and there a patch is no wider than a wavelength, where the local-slab transfer used here fails. The
uneven model therefore shows what a spread of thickness gives on average, not what a given slice
does; in the fits it is a loss rising with frequency plus a lag, with two free sizes. Scattering is modelled as spheres of radius r at volume fraction v ≤ 0.2 in the tissue, a coherent
medium (Foldy 1945) with the Mie forward amplitude in the absorbing host (complex size parameter); it
reduces to Maxwell Garnett for small spheres and reproduces the Mie extinction of Bohren and Huffman
(1983), both tested (Sim 54).

Table 2. What each candidate leaves unexplained (Sims 69–72; for uneven thickness the best sizes at
mean −30 % and rms 40 % are grid limits; built by `experiments/exp/table2.py`
into `results/table2.md` with one rule for every row). Each candidate is judged in each reference by
the better of its two fits, alone (s = 0) or with a timing drift s of both pulses (the drift fit
contains the other, so a cause is excluded only if both are); "fits with a drift" counts the references
where the drift fit was better. "No change" and "shift only" are the two claims that the tissue pulse
changes and that it changes shape. D: effect size (% of the trace rms) of the fit with the candidate at
its best size and no loss term, against the fit with the causal loss increment and a free shift, over
the eight reference fits with floors of 1.0–2.5 dB/THz at 2 / 4 / 5 min. D ÷ threshold: the largest over the three
scans in each reference, against fitting noise and tracing added in quadrature (Sec. 5.6); range over
references. Excluded in: references in which it exceeds 1; in brackets, borderline verdicts among all
eight (flip in > 10 % of bootstrap resamples, or ratio within 0.87–1.15). Last column: the same for
the four references with floors of 0 and 0.5 dB/THz. Warming A: Liebe's full temperature change; B:
only the slow relaxation warms; C: only the fast one; shared: w B + (1 − w) C, w and the temperature
rise both fitted; each scaled to the measured absorption change at 1 THz. Gas that displaces tissue
and lipid-like fragments (r 15–80 µm) are best at zero size in every reference (Sims 53–54), i.e. no
change.

| Candidate (sizes searched) | Best sizes | Fits with a drift (of 8) | D (%), 2 / 4 / 5 min | Largest D ÷ threshold | Excluded in (of 8): all (borderline) | Floors 0–0.5: D ÷ threshold; excluded (of 4), borderline |
|---|---|---|---|---|---|---|
| No change | – | 0 | 17.9–20.9 / 19.6–22.5 / 20.4–23.0 | 2.18–3.29 | 8 (0) | 3.41–4.02; 4, 0 |
| Shift only | – | 8 | 7.7–12.5 / 9.6–14.2 / 13.1–16.8 | 1.94–3.19 | 8 (0) | 2.58–3.22; 4, 0 |
| Thicker tissue path, uniform (0–51 %) | 6–33 % | 5 | 3.7–8.5 / 4.4–10.0 / 7.1–12.4 | 0.74–1.40 | 6 (4) | 1.21–1.76; 4, 0 |
| Uneven thickness across the beam (mean −30 to +51 %, rms 0–40 %) | mean -30 to 18 %, rms 16–40 % | 4 | 2.2–3.8 / 3.1–5.0 / 0.0–5.3 | 0.51–0.83 | 0 (0) | 0.67–1.11; 1, 1 |
| Gas in a slice that swells (r 20–100 µm, v 0–0.2) | r 20–50 µm, v 0.16–0.20 | 8 | 1.9–5.5 / 5.5–9.1 / 7.9–10.6 | 0.84–1.61 | 5 (4) | 1.33–2.41; 4, 0 |
| Node states, all 80 (λ 0.5–6 of every extreme) | λ 2–6 | 4 | 3.6–6.4 / 5.9–8.7 / 6.3–9.8 | 0.78–0.95 | 0 (5) | not refitted |
| Ions in all tissue water (0–60 S/m) | 5–40 S/m | 7 | 6.9–14.7 / 8.2–15.8 / 10.8–15.4 | 1.33–2.15 | 8 (0) | 2.05–2.21; 4, 0 |
| Axoplasm conduction (0–150 S/m) | 10–100 S/m | 7 | 6.8–14.3 / 8.1–15.5 / 10.9–15.1 | 1.39–2.36 | 8 (0) | 1.88–2.38; 4, 0 |
| Warming A, Liebe ×2.31 | 0–6 K | 1 | 11.4–15.7 / 13.1–17.8 / 15.8–20.1 | 1.49–2.52 | 8 (0) | 2.42–3.29; 4, 0 |
| Warming B, slow relaxation (×1.04) | 3–24 K | 8 | 5.4–9.7 / 6.3–10.5 / 9.2–12.5 | 1.20–1.59 | 8 (0) | 1.60–2.10; 4, 0 |
| Warming C, fast relaxation (×3.26) | 4–6 K | 0 | 8.7–12.2 / 11.5–14.8 / 14.0–17.0 | 1.22–2.20 | 8 (0) | 2.05–2.53; 4, 0 |
| Warming, slow and fast shared freely (w B + (1 − w) C) | w 0.6–1.0, 8–24 K | 0 | 5.5–9.8 / 7.8–11.6 / 10.3–13.2 | 0.94–1.37 | 7 (4) | 1.49–1.74; 4, 0 |

**What the table shows.** The tissue pulse changes (no change exceeds the threshold by 2.2–4.0
times in all twelve references), and not only in time (a shift only, 1.9–3.2 times). Extra conduction in the bulk alone is firmly excluded: ions in the tissue water over 0–60 S/m, which covers the ~2–3 S/m of
electrolysis protons, and conduction in the axoplasm over 0–150 S/m exceed the threshold in every
reference by 1.3–2.4 times. Warming as Liebe's model has it (A), or with only one relaxation warming
(B, C), is firmly excluded too (1.2–3.3 times). By the better of its two fits a slice whose thickness varies across the beam is
within the threshold in all eight references with floors of 1.0–2.5 dB/THz (0.51–0.83 times) and,
with floors of 0–0.5, within it in three and excluded but borderline in the fourth (0.67–1.11).
Without a drift it passes in nine of the twelve references (one borderline, 0.91), with best sizes a
mean 6–18 % thicker and an rms of 4–25 % of the thickness (9–55 µm on the 220 µm optical path); it is
excluded in three, all with a + b ≤ 1: firmly at floors 0 and 1.0 (1.16–1.17 times) and borderline at
0.5 (1.14). In those three the drift fit passes only at the edge of the grid, a mean 12–27 % thinner
with an rms of 40 %, the most uneven slice tested and the one a single slice realises least (Sim 74);
we count these three verdicts as unresolved. The drift fits reach the grid edge in ten of the twelve
references, so the best sizes listed in Table 2 for this row are partly grid limits, not fitted values. It is the mean over many patches whose limits are given above. A uniformly thicker path, gas with a drift and warming
of both relaxations are not firmly excluded with floors of 1.0–2.5 dB/THz and are firmly excluded with
floors of 0–0.5. Node states are not firmly excluded in any reference (0.78–0.95, five borderline), at
2–6 times their every-extreme change in the better fit of each reference (1–6 over all fits, which is 12–99 times the change at internode 200 µm and packing 0.6), and only as states an aqueous gap cannot hold
(below). The waveforms therefore exclude extra bulk conduction alone but do not identify the cause: each
remaining candidate supplies a loss rising with frequency and a lag through one or two fitted sizes,
and the fits cannot tell them apart. Combinations of candidates are not tested. Earlier versions of this table rested on
significance tests that proved unusable (Sec. 5.6; Sims 55b, 57, 60, 61, 63, 64) and on a tracing
yardstick that ignored fitting noise (Sim 66); their verdicts are withdrawn. Sim 68 used too little
noise and judged at a different percentile in its script and text; Sim 69 corrected both. Sims 70–71
limited gas to its dilute range with the radius free and added warming models B, C and their mixture.
Until Sim 72 every slab candidate had one uniform thickness; the claim, based on them, that no tested
bulk cause produces the change when the 30 s loss belongs to the slice is withdrawn.

**Node states.** With every extreme the largest delay any node state gives at 0.6 THz is 0.064 rad
(ε = −3.0 − 0.1j, combined materials; loss slope 1.1 dB/THz over 0.2–1.7 THz; Sim 70) against a lag of
0.16–0.28 rad. The steepest loss slope, 5.2 dB/THz, comes from states with ε′ between −1.7 and −1.4
and ε″ = 0.1, and these advance the pulse (−0.08 to −0.11 rad at 0.6 THz). No state gives both the
lag and the loss slope. Inside the fits the best of all 80 states (Table 2) need 1–6 times every extreme:
ε′ from −1.3 to −1.0 with a drift, from −5 to −2.7 alone, or a node conducting at 1000 S/m. At
internode 200 µm and packing 0.6 the change is 12–17× smaller still. Every one of these states is a
nearly lossless negative permittivity (ε″ = 0.1) or a conductivity 500 times that of saline, which an
aqueous gap cannot hold: ions large enough to make ε′ negative bring ε″ ≳ 80 (Sec. 4.6). A frequency-independent ε is not
causal (Kramers–Kronig), so these states are bounds on what a gap can do, not models of one. Node
states are thus excluded by what a gap can hold and by the forward bounds, not by the waveforms,
which accept them with a drift in all eight references with floors of 1.0–2.5 dB/THz (node states
were not refitted at floors 0–0.5) at 1–6 times every extreme. The sign of the first-order change behaves as it should: a metallic node, more polarisable
than the water it replaces, delays (+0.042 rad at 0.6 THz), and a node with ε′ = 1 advances
(−0.039 rad).

**Node states at 2 THz.** With every extreme the 2 THz loss bound (10 dB, Sec. 4.6) approaches the
14 dB by which the 30-s scan falls below the unbiased one at its 1.99 THz dip (Sec. 5.6). In the
two-pulse model a narrow dip is a cancellation, so the test is whether a node state moves the
model's dip to 1.99 THz with the measured depth, itself a lower bound (Sec. 5.6). Each of 65 distinct
node states, with every extreme, was applied to T in the 30-s model as a change held constant over
1.8–2.4 THz, for each of six bands used to fit the tissue phase and for b = 0.034–0.044 (Sim 43).
Without a node change the model's lowest point in 1.85–2.15 THz lies at 2.08–2.15 THz, −32 to
−39 dB (b = 0.039). Some states move it to within 0.03 THz of 1.993 THz, but at most −37.5 dB deep
(b = 0.034; −34.6 dB at b = 0.039) against the published −43.5 dB. The state that moves it closest
(ε = −2.0 − 0.1j; 1.969 THz, −32.8 dB) adds 5.4 dB of loss at 2 THz, which spoils the cancellation.
The dip is not reproduced with or without node states; since the model is 5.5 dB off in this band
without them, this does not exclude them.

**Ions and conduction.** Extra conductivity in the tissue water or the axoplasm gives a loss nearly
flat in frequency. Over 0–60 S/m added to all tissue water, which covers the protons of electrolysis
(0.7–1.1 mA for 300 s is 2.2–3.4 µmol of H⁺ in a ~40 µL slice, 55–85 mM or ~2–3 S/m if it stays
there), and 0–150 S/m added to the axoplasm, its best fit is excluded in every reference, by 1.3–2.4 times the threshold (Table 2). Its best sizes, 5–40 S/m in the tissue water and 10–100 S/m in the
axoplasm, compare with ~1.4 S/m physiological and ~25 S/m for saturated NaCl.

**Warming.** Liebe's temperature dependence of water is fitted below 1 THz and 30 °C, so four
versions cover how warming acts above 1 THz, each scaled to the measured 4–21 °C absorption change at
1 THz (46 cm⁻¹; Millon et al. 2022): A, Liebe's full change (×2.31); B, only the slow relaxation warms
(ε_s and τ1; ×1.04); C, only the fast relaxation warms (τ2; ×3.26); and a mixture w B + (1 − w) C with
w fitted (Sims 70–71). A lowers the loss slope over 0.3–1.5 THz, the opposite of the change; C raises
it. A, B and C are firmly excluded in every reference, alone and with a drift (1.2–3.8 times the
threshold; Table 2). The fitted mixture is excluded in all four references with floors of 0–0.5
(1.5–1.9 times) but not firmly with floors of 1.0–2.5: alone it is within the threshold in one
reference and borderline in four, at w = 0.6–1.0 (the slow relaxation warming at 0.6–1.0 and the fast at
0–1.3 times Liebe's rates) and 8–24 K. The measured temperature dependence of water across 0.2–2 THz (e.g. Rønne et al.
1997) would fix w; we did not have it in tabulated form. Liebe as given needs 9 K for the change that
×2.31 gives at 4 K (Sim 53). The bias heats the slice by 0.6–0.9 K/s if no heat is lost (Sec. 5.4);
heat flows into the 1.2 mm PTFE plates within ~13 s, so the real rate is lower, and 8–24 K within
2–5 min is possible. Warming is therefore excluded as the sole cause only if the 30 s loss belongs to
the slice.

**More tissue, or uneven tissue, in the beam.** With floors of 1.0–2.5 dB/THz a path 6–9 % thicker (6–33 % with a drift,
searched to 51 %) is excluded alone in all eight references (two borderline) and with a drift is within
the threshold in two and borderline in three; with floors of 0–0.5 it is firmly excluded (1.21–1.76 times; Table 2). So the waveforms exclude it only if the 30 s loss belongs to the slice. On the
~220 µm optical path 6–9 % is 13–20 µm, ~3–8 mg over 20 × 20 mm depending on a 100 or 220 µm slice;
30 % is 30–66 µm, ~12–26 mg. It
could be a slice swelling (water drawn in under the bias, e.g. by electro-osmotic flow) or tissue
moving into the beam. The slice is held in a PTFE cuvette (H1, p. 2; H2 gives a 100 µm gap, while H1's
Fig. 1b shows the slice on one 1.2 mm plate, so whether the gap is fixed is not stated). A slice whose
thickness varies across the beam is not firmly excluded in nine references and unresolved in three (Table 2), as a mean over many
patches (Sim 74). A slice that dries, swells or tears unevenly under the bias could do this, and
the openings that pass the second pulse would be its zero-thickness extreme. H4's sample 2, the one new sample, transmits more and earlier after bias; a uniformly thinner path
(8–16 %) fits about as well as a stronger, earlier pulse (Sims 76–77).

**The reported selectivity.** More or uneven tissue in the beam, warming, a stage drift, scattering by gas or
disrupted tissue and conduction in the tissue water would each act for both THz polarisations, and DC along the fibres (σ∥ ≈ 7σ⊥) heats more than DC
across them. The reports state that the effect appears only for the THz field along the fibres and
only for DC across them (H1, p. 18; H3). If those controls were done as described, they contradict
these causes, including those the waveforms do not exclude: uneven or thicker tissue, gas with a
drift and warming.
The E∥-selective ones are excluded: a conducting axoplasm by the waveforms at every size up to
150 S/m, node states by the forward bounds. The selectivity therefore remains unexplained. The controls are described in words only, with no waveforms or spectra, so they cannot be tested on the published data. The two-pulse model is
therefore a decomposition of the published waveforms, not an explanation of the reported effect.

### 5.8 What matches and what does not

| Reported | This work |
|---|---|
| Dip near 0.6 THz (9–20 dB) | No node or fibre resonance in 0.1–3 THz. The second pulse plus a changed tissue pulse (lagging 0.16–0.28 rad at 0.6 THz, with a loss rising with frequency that trades off against the lag) is consistent with it (0.57–0.65 THz against 0.56–0.65 measured); q is fitted. H4's "deep resonance" of sample 1, given without a frequency (spectrum minimum near 0.55 THz, index jump at 0.60 THz), is consistent with the same cancellation (model 0.56–0.59 THz), on a scan that matches H3's 5 min scan. Judged against fitting noise and tracing, the change firmly excludes extra conduction in the tissue water or axoplasm alone in every reference (1.3–2.4× the threshold); uneven or thicker tissue, gas, warming and node states each pass under some reading, so the waveforms do not identify the cause; node states pass only at 1–6× every extreme and as states an aqueous gap cannot hold (Table 2) |
| Drift 0.7 → 0.5 THz | The dip deepens and moves from ~0.7 to 0.62–0.65 THz as b and q grow; water-vapour lines cannot drift. No traced spectrum has a dip below 0.56 THz, so the 0.5 THz endpoint is not in the published figures |
| 2 THz dip (30 s) | Not reproduced, with or without node states. A real feature of one scan (15 dB, twice any other single-curve dip), its depth limited by dynamic range; the two pulses reach that depth but, with the traced phase extrapolated, at 2.10–2.32 THz rather than 1.99 THz. With node states (65, every extreme, six phase bands, b = 0.034–0.044) the model's minimum within 0.03 THz of 1.99 THz is at most −37.5 dB (measured ≤ −43.5 dB) |
| Only E∥ fibres, only DC across fibres | Not explained. The causes the waveforms leave (uneven or thicker tissue, gas with a stage drift, warming) act for both polarisations, and DC along the fibres would heat more. Described without data |
| H4 sample 2 (the one independent sample): rise from 0.4 THz, 5 dB near 2 THz; index falling over the band; resonance supposed below the band | No dip and no second pulse; tissue pulse stronger and earlier after bias, the opposite of H1–H3; a path 8–16 % thinner (5.4 %, with a drift of the opposite sign) and a stronger, earlier pulse (8.0 %) fit about equally, no change 21.2 % (Sim 77) |
| Onset ~5 min, recovery ~1 h | Second pulse grows 4 → 12 → 17 → 21 % (0.5, 2, 4, 5 min); H4's sample 1 at "6 min" (24 %) matches H3's 5 min scan within 10 % (2–3 times the tracing noise), so it is a later scan of the same
sample or H3's, near saturation (b = 0.24 with H4's empty pulse), consistent with H1's statement that
nothing changes after ~6 min, and not an independent sample; recovery not traced |
| No significant change of average power | Traced energy −1.2 dB (H3), consistent |

Decisive measurements are cheap: logging current and slice mass (6–9 % more tissue in the beam is
~3–8 mg over 20 × 20 mm, depending on whether the slice is 100 or 220 µm thick), imaging the slice
and mapping its thickness across the beam for gas, openings or disruption, a second empty-cuvette scan after each biased one (stage drift), measuring the cuvette gap and the slice thickness during bias, an N₂ purge, a reference through a
heated or dried slice without bias, a measurement with the slice clamped between two windows, and
release of the raw data.

## 6. Discussion and limits

- **Linear, passive response.** DC-driven structural changes of fibres (swelling, electro-osmotic
  water flow, pH fronts) are tested only through their bulk effect (heating, water loss, gas) and
  through a slice that covers the beam partly or unevenly.
  Active gain at the node (Liu Y et al. 2020) is outside a passive, energy-conserving calculation.
- **Single fibre plus effective medium.** The slice bounds rest on the induced dipole of one fibre
  and a quasi-static local field in dense packing. Finite fibre segments in a 70–100 µm slice lower
  the E∥ response, so the infinite-fibre values are upper bounds.
- **Band.** 0.1–3 THz. Song et al. (2026) work at 10–40 THz; mid-IR guiding (Liu G et al. 2019) is
  not tested.
- **No raw experimental data.** The waveform analysis rests on traced published figures, H1–H3
  from one measurement session and H4 from one figure per sample; the full two-pulse fit leaves
  10–15 % of the trace rms at 2–5 min (Sec. 5.6), several times the tracing noise, and the
  tissue-path factor q depends on the fitted band. H4's sample 1 is the sample of H1–H3 and its
  biased scan matches H3's 5 min scan, so H4 adds one independent sample.
- **Uneven thickness** is modelled as the mean over many independent patches. One slice under a 5 mm
  beam departs from that mean by a quarter to three times the change, depending on patch size, and only
patches of ≲ 0.25 mm average out at 0.6 THz, where
  the local-slab transfer fails (Sim 74), so its fit shows what a spread of thickness gives on
  average, not what a given slice does.
- **Statistics of the waveform fits.** The misfit is systematic, and no significance test we built
  kept its nominal rate on data with no change (Sims 61–65; Sec. 5.6). Fits are judged by effect
  size against a threshold from fitting noise (100 phase-randomised draws per fit, scaled for what
  the full fit absorbs) and from six pairs of traced unbiased pulses; tracings of the biased curves
  could differ more if their shapes are harder to follow. A true cause would be excluded by chance in
  about 5 % of references, and one threshold carries ~11 % Monte Carlo error, so ratios within
  0.87–1.15 are called borderline; the exclusion of conduction holds in all twelve references at 1.3
  times the threshold or more. D compares models that are not nested, candidate sizes are fitted per
  scan, and the noise is taken as spread evenly in time. Which uniform bulk causes are excluded depends on
  whether the 30 s loss is read as an artefact of the reference or as part of the slice (Sec. 5.6;
  the tracing pairs favour the second); uneven thickness is not firmly excluded under either reading (three verdicts unresolved). Gas is tested only where the dilute model holds (v ≤ 0.2), the
  temperature dependence of water above 1 THz is bracketed rather than measured, and combinations of
  causes are not tested. The full fit itself leaves 10–15 % at 2–5 min, so the decomposition is
  approximate. The candidates are searched on grids of one or two sizes, so a cause with another
  frequency dependence is not excluded.
- **Loss term.** Whether the loss is causal and how it continues above the traced band (~1.7 THz)
  changes how the change splits into a loss and a shift, not the loss slope or the lag (Sim 55).
  Scattering is modelled as a dilute coherent medium (Foldy) of spheres in an absorbing host.
- **Node states between solved frequencies.** Each state's change is interpolated linearly in
  frequency between the frequencies it was solved at; the combined-materials states were solved at
  0.6 and 3 THz only.
- **Mesh at the passive-ε peak.** Node power converges within 0.5 %, but the slice change at the
  ring-resonance peak falls with refinement (0.123/0.107/0.097 dB at 3 THz, 21 %) because the gap's
  corners are near-singular there; the coarse value is used as an upper bound (Sec. 4.6).
- **Slice thickness.** The measured loss and delay give an optical thickness of ~220 µm against
  the stated 100 µm cuvette (Sec. 5.5). We cannot tell a wider gap from tissue that is not
  water-based; the bounds use the optical thickness as a conservative scaling.
- **Conductivity bound** was not rerun with ±20 % materials; with the passive-ε material factor
  (×1.45) and every extreme combined it is 1.2 dB, 7× below the 9 dB dip depth.

The result is a negative one for the resonance hypothesis and a positive one for dosimetry: node
absorption in this band is fixed by water loss, node volume and a geometric factor that the
law gives within a few percent, with no free parameter fitted to full wave.

## 7. Conclusion

In 0.1–3 THz the node of Ranvier is a quasi-static water gap. Its absorbed power follows P_node =
½ωε0ε″_w|E0|²V_nF with F from an interpolation law, exact in its limits, in node length over myelin
thickness, g-ratio and permittivity contrast, within 2.3 % of full wave for fibres up to 10 µm. For
the THz field along the fibres, as in every published scan, no node state (ε″ ≥ 0.1) changes a 100
µm slice by more than 1.1 dB at internodes ≥ 200 µm and packing ≤ 0.6. The H1–H3 waveforms, all from one session and one sample, decompose into the unbiased pulse with an extra loss that rises with frequency
and bias time (1.7–7.1 dB per THz above its 30 s value in the reference fits, less where a candidate adds delay), lagging by 0.16–0.28 rad at
0.6 THz, plus a second pulse near the empty-cuvette delay that grows over the four biased scans. This
holds over 36 fitting assumptions, and without the change twelve reference fits miss by 2.2–4.0
times what fitting noise and tracing allow. It is consistent with the 0.6 THz dip but not with the 2 THz dip.
Significance tests on traced figures did not keep their nominal rate on data with no change, so the
causes are judged by effect size against fitting noise and tracing. Extra conduction alone, up to
60 S/m in the tissue water or 150 S/m in the axoplasm, is excluded in every reference fit. The
waveforms do not identify the cause. A slice whose thickness varies across the beam is not firmly excluded in nine of the twelve reference fits and unresolved in three, and only as a mean over many patches that one slice under the beam need
not realise; a uniformly thicker path, gas and warming pass only under one reading of the 30 s scan.
Node states pass only at 1–6 times every geometric extreme and as nearly lossless negative
permittivities or conductors at 1000 S/m, which an aqueous gap cannot hold, so they are excluded by physics rather than by the
waveforms. A fourth report shows the sample of H1–H3 in a scan that matches an earlier one, with its reported resonance consistent with the cancellation between the two pulses, and one new sample whose tissue pulse after bias is stronger and earlier, the opposite of the lag and
extra loss of H1–H3, with no second pulse. The second pulse's lead over the empty-cuvette pulse is comparable to the difference between two
reports' empty-cuvette pulses, so whether it is real is not decided. Every remaining cause would act for both THz polarisations, against the
reported selectivity, which is described without data. Settling the cause needs the raw data or the
control measurements listed in Sec. 5.8.

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
6. `fig_bypass.png` — traced H3 waveforms with the two-pulse fit; spectra and fits.
7. `fig_paper1_bypass.png` — traced H1 waveforms (30 s and 4 min at 100 V) with the two-pulse fits.
8. `fig_paper1_spectra.png` — H1's published spectra relative to the empty cuvette and the two-pulse
   model with a, b, τ (and q = 1.09 at 4 min) from the waveform fits; nothing fitted to the spectra.
9. `fig_tissue_change.png` — extra loss vs extra delay of the slice at 0.6 and 1.0 THz. Horizontal
   bands: the lag of the tissue pulse at 2–5 min over the reference fits with the 2.07 dB/THz floor
   (Sim 52; dark: best fits; light: down to the 95 % lower edge; 0.16–0.28 rad over all 36 references
   of Sim 55), scaled with frequency. Vertical band: the causal loss increment of the tissue pulse
   over the 30 s floor at 2–5 min (Sim 53). Cross: a time shift of the tissue
   pulse. Curves: heating (two water models), water gain into thickness, a thicker slab, ions, axoplasm
   conduction; dots: node states with every extreme, computed at 0.6 THz only (Sim 46). The bands pool
   three scans and two frequencies; a cause must match one scan over the whole band at once, which is
   the test of Table 2 (Sim 69).

## References

- Arancibia-Carcamo IL, Ford MC, Cossell L, Ishida K, Tohyama K, Attwell D (2017) Node of Ranvier
  length as a potential regulator of myelinated axon conduction speed. eLife 6:e23329.
  doi:10.7554/eLife.23329
- Berry E, Fitzgerald AJ, Zinov'ev NN, Walker GC, Homer-Vanniasinkam S, Sudworth CD, Miles RE,
  Chamberlain JM, Smith MA (2003) Optical properties of tissue measured using terahertz pulsed
  imaging. Proc SPIE 5030:459–470. doi:10.1117/12.479993
- Bohren CF, Huffman DR (1983) Absorption and Scattering of Light by Small Particles. Wiley, New York
- Foldy LL (1945) The multiple scattering of waves. I. General theory of isotropic scattering by
  randomly distributed scatterers. Phys Rev 67:107–119
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
  neural system in the terahertz frequency range. Reports NAS RA 124(2):16–23.
  doi:10.54503/0321-1339-2024.124.2-16
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
- Rønne C, Thrane L, Åstrand P-O, Wallqvist A, Mikkelsen KV, Keiding SR (1997) Investigation of the
  temperature dependence of dielectric relaxation in liquid water by THz reflection spectroscopy and
  molecular dynamics simulation. J Chem Phys 107(14):5319–5331. doi:10.1063/1.474242
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
