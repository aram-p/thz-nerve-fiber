#!/usr/bin/env python3
"""Generate thesis report DOCX — English text + all simulation figures."""

from docx import Document
from docx.shared import Pt, Inches, RGBColor, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH
import os

FIGURES_DIR = "/Users/ap/Downloads/Telegram Desktop/thz-nerve-fiber/figures/out"
OUT_PATH = "/Users/ap/Dev/thz-nerve-fiber/thesis_report.docx"

FIGURES = [
    ("sim1_freq_sweep.png",
     "Figure 1. Frequency sweep of |E|-field at the Node of Ranvier (0.1–2.0 THz, σ = 0). "
     "Peak |E| = 2.60 at f = 0.632 THz, coinciding with the experimentally observed 0.6 THz resonance."),
    ("sim2_sigma_sweep.png",
     "Figure 2. Node conductivity σ sweep at f = 0.6 THz across seven orders of magnitude. "
     "Three regimes: linear enhancement (σ < 10 S/m), peak coupling at σ ≈ 100 S/m, "
     "metallic limit (σ ≥ 10⁴ S/m)."),
    ("sim3_node_length.png",
     "Figure 3. Node-length sensitivity — peak |E| at the node vs node length (10–100 µm) "
     "at f = 0.6 THz, σ = 0. Monotonic increase from ≈ 1.6 to ≈ 2.9."),
    ("sim4_dispersion.png",
     "Figure 4. Frequency-dependent dielectric permittivity ε(f): double-Debye water, "
     "constant myelin (ε = 4.5 − 0.5i), and node material with added σ-term (0.1–3 THz)."),
    ("sim5_geometry_viz.png",
     "Figure 5. 3D unit-cell geometry: solid axon (r = 5 µm, L = 240 µm), two myelin annular "
     "sheaths (r = 5–7 µm, L = 100 µm each), node of Ranvier annulus (r = 5–7 µm, L = 40 µm), "
     "surrounding water box (±20 µm). Incident plane wave E₀ x̂ propagates along +z."),
    ("sim6_grating.png",
     "Figure 6. Analytical wire-array diffraction grating model (Tretyakov surface admittance). "
     "Period 50 µm, wire radius 2.5 µm in water host. "
     "Transmittance, reflectance, and absorbance for σ = 10–10⁵ S/m."),
    ("sim7_wavelength_scales.png",
     "Figure 7. Free-space wavelength λ₀ and in-water wavelength λ_water vs frequency (0.1–3 THz), "
     "with horizontal lines at axon, myelin, node, internode, and unit-cell geometric scales. "
     "λ_water/2 crosses the internode length near 0.6 THz."),
    ("sim8_periodic_array_3d.png",
     "Figure 8. 3D visualisation of a 3×3 periodic array of myelinated nerve-fiber unit cells. "
     "Central cell highlighted; surrounding cells faded. Illustrates the diffraction-grating "
     "geometry imposed by Floquet-periodic boundary conditions."),
    ("sim9_3d_field_at_resonance.png",
     "Figure 9. 3D scatter plot and cross-section heatmap of |E|-field at the resonance "
     "(f = 0.632 THz). Orange lines mark node boundaries; white horizontals mark axon/myelin radii."),
    ("sim10_field_isosurfaces.png",
     "Figure 10. Nested isosurfaces of |E|/max at thresholds 0.55, 0.70, and 0.85. "
     "Right panel: cross-section-integrated energy density ∫∫|E|² dx dy vs axial position z."),
    ("sim11_freq_waterfall.png",
     "Figure 11. Frequency-sweep waterfall — full axial profile |E|(z) at each frequency (0.1–2.0 THz). "
     "The heatmap shows the 0.632 THz band as a distinct horizontal feature confirming "
     "field localisation at the node."),
    ("sim12_node_length_3d.png",
     "Figure 12. 3D axial-profile ribbons for five node lengths (10–100 µm) at f = 0.6 THz. "
     "Inset: peak |E|-at-node vs node length — clean monotonic scaling."),
    ("sim13_resonance_contrast_3d.png",
     "Figure 13. Cross-section field maps: off-resonance (f = 0.328 THz) vs on-resonance "
     "(f = 0.632 THz). The resonance is a localised enhancement at the node, "
     "not a global field amplification."),
    ("sim14_dispersion_3d.png",
     "Figure 14. 3D dispersion trajectories of Re(ε) and −Im(ε) with frequency as vertical axis. "
     "Water traces a curved relaxation arc; myelin is a vertical constant line."),
    ("sim15_freq_nodelen_surface.png",
     "Figure 15. Peak |E|-at-node annulus over a 4×6 grid of (node-length, frequency). "
     "3D surface (left) and 2D heatmap (right). Surface is nearly flat in frequency, "
     "confirming the resonance is on-axis, not in the conductive annulus."),
    ("sim16_mesh_convergence_3d.png",
     "Figure 16. Mesh convergence at f = 0.632 THz: coarse (~24 k), baseline (~26 k), "
     "and fine (~41 k) elements. Peak |E| varies by ≈ ±10 %; "
     "qualitative resonance pattern is stable across all three refinements."),
    ("sim17_power_dissipation.png",
     "Figure 17. Per-domain EM power dissipation (W) vs frequency (0.15–2.1 THz) "
     "at σ_node = 1 S/m. Myelin dominates total absorption; "
     "node contribution is ≈ 1–4.5 % and roughly frequency-independent."),
    ("sim18_e_parallel_fibre.png",
     "Figure 18. Frequency sweep with E-field polarised parallel to the fiber axis "
     "(k along x̂, E along ẑ). Annular |E| is ≈ 25 % higher than in the perpendicular "
     "configuration, reproducing the experimentally observed polarisation dependence."),
    ("sim19_lorentzian_fit.png",
     "Figure 19. Lorentzian fits to resonance peaks (⊥ vs ∥ fibre × axial vs annular sampling). "
     "Principal result: f₀ = 0.605 ± 0.019 THz, FWHM = 0.28 THz, Q = 2.16 (⊥ fibre, axial)."),
    ("sim20_polarisation_compare.png",
     "Figure 20. Direct overlay of frequency spectra for E ⊥ fibre (Sim 1) and E ∥ fibre (Sim 18). "
     "Parallel configuration raises the annular field by ≈ 1.20×, consistent with absorption "
     "appearing only in the parallel orientation."),
    ("sim21_fibre_along_x.png",
     "Figure 21. Frequency sweep with fiber geometry properly rotated so the cylinder axis "
     "lies along x̂. Annular |E| stays in 1.95–2.45 across the band — "
     "polarisation enhancement is preserved."),
    ("sim22_dense_peaks.png",
     "Figure 22. Dense frequency sweeps with refined mesh (15 µm / 3 µm at node) near both peaks. "
     "Lower peak: f₀ = 0.619 ± 0.008 THz, Q = 15. Upper peak: f₀ = 1.938 ± 0.010 THz, Q = 30.5. "
     "Headline quantitative result — both experimental peaks reproduced."),
    ("sim23_power_dissipation_sigma100.png",
     "Figure 23. Per-domain power dissipation at peak-coupling conductivity σ_node = 100 S/m. "
     "Water dominates (~10 W); node contribution ≈ 200 mW (200× increase vs σ = 1 S/m). "
     "Node fractional share peaks at low frequencies (~1.9 %)."),
    ("sim24_mode_identification.png",
     "Figure 24. Mode identification — half-wave standing-mode predictions from "
     "λ_water(f)/2 = segment length. "
     "Node (40 µm) → f_pred = 2.043 THz vs f_obs = 1.938 THz (−5 %); "
     "internode + half-node (120 µm) → f_pred = 0.589 THz vs f_obs = 0.619 THz (+5 %)."),
    ("sim25_field_at_both_peaks.png",
     "Figure 25. Side-by-side 3D isosurface renders and cross-section heatmaps at both peak "
     "frequencies (f = 0.619 THz and f = 1.938 THz). Tighter wave fringes at the high peak "
     "confirm the shorter half-wave standing mode on the 40 µm node segment."),
    ("sim26_convergence_both_peaks.png",
     "Figure 26. Mesh convergence at both resonance peaks (three refinements × two frequencies). "
     "Peak |E| varies by ~10 % across mesh levels; peak positions converge to better than 0.5 %."),
]


def font(run, size=12, bold=False, italic=False):
    run.font.name = "Times New Roman"
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = RGBColor(0, 0, 0)


def heading(doc, text, size=13):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(14)
    p.paragraph_format.space_after = Pt(4)
    r = p.add_run(text)
    font(r, size=size, bold=True)
    return p


def body(doc, text):
    p = doc.add_paragraph()
    p.paragraph_format.first_line_indent = Cm(1.25)
    p.paragraph_format.space_after = Pt(6)
    p.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    r = p.add_run(text)
    font(r)
    return p


def figure(doc, fname, caption):
    path = os.path.join(FIGURES_DIR, fname)
    if not os.path.exists(path):
        print(f"  missing: {fname}")
        return
    p = doc.add_paragraph()
    p.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before = Pt(10)
    p.paragraph_format.space_after = Pt(2)
    p.add_run().add_picture(path, width=Inches(5.5))
    cap = doc.add_paragraph()
    cap.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.CENTER
    cap.paragraph_format.space_after = Pt(14)
    r = cap.add_run(caption)
    font(r, size=10, italic=True)


# ── Build ─────────────────────────────────────────────────────────────────
doc = Document()
sec = doc.sections[0]
sec.top_margin = sec.bottom_margin = Cm(2.5)
sec.left_margin = Cm(3.0)
sec.right_margin = Cm(1.5)


def centered(doc, text, size=12, bold=False, after=6):
    p = doc.add_paragraph()
    p.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_after = Pt(after)
    font(p.add_run(text), size=size, bold=bold)


centered(doc, "MASTER'S THESIS", 14, bold=True, after=4)
centered(doc,
    "Computational Physics Study of THz Electromagnetic Wave Interaction\n"
    "with a Single Myelinated Nerve Fiber — Resonant Absorption at the Node of Ranvier",
    13, bold=True, after=4)
centered(doc, "Student: Aram Poghosyan", 12, after=20)


# ══════════════════════════════════════════════════════════════════════════
heading(doc, "1.  Thesis Work Plan (Schedule)")
# ══════════════════════════════════════════════════════════════════════════

body(doc,
    "The thesis is divided into four sequential phases, each with defined deliverables "
    "and success criteria.")

body(doc,
    "Phase I — Literature review (months 1–2). Review the three experimental papers by "
    "Hovhannisyan & Makaryan (Armenian Journal of Physics 2024; NAS Armenia Reports 2024; "
    "IRPhE 2024), covering: the time-domain THz spectroscopy methodology, the observed "
    "resonant absorption peaks at 0.6 THz and 2 THz in spinal-cord samples, the polarisation "
    "and voltage-dependence of those peaks, and the proposed diffraction-grating physical model. "
    "Supplement with background reading on myelinated nerve fiber structure (Hodgkin–Huxley), "
    "double-Debye water dispersion, and COMSOL FEM best practices for EWFD problems.")

body(doc,
    "Phase II — Model construction and software pipeline (months 3–4). Build the 3D "
    "geometry in COMSOL via Python/mph: solid axon cylinder (r = 5 µm, L = 240 µm), two "
    "myelin annular sheaths (r = 5–7 µm, L = 100 µm), one node of Ranvier annulus "
    "(r = 5–7 µm, L = 40 µm), and a surrounding water box (±20 µm) for Floquet-periodic "
    "lateral BCs. Assign double-Debye water, constant myelin (ε = 4.5 − 0.5i), and "
    "σ-dependent node permittivity. Wire up the EWFD scattered-field study with scattering BCs "
    "on the z-axis end-caps and MUMPS direct solver. Set up the YAML sweep configuration "
    "system and HDF5/manifest storage layer.")

body(doc,
    "Phase III — Simulation campaign (months 5–8). Execute 26 distinct simulation scripts: "
    "frequency sweep (0.1–2.0 THz, 26 points), node conductivity sweep (σ = 0 to 10⁸ S/m), "
    "node-length sensitivity (10–100 µm), polarisation comparison (E ⊥ and E ∥ fibre), "
    "dense focused sweeps near each resonance peak with refined mesh, 3D field visualisations, "
    "per-domain power dissipation spectra, and mesh convergence study.")

body(doc,
    "Phase IV — Analysis, interpretation, and writing (months 9–12). Fit Lorentzian profiles "
    "to identify peak frequencies and Q-factors. Apply the half-wave standing-mode model to "
    "explain both resonances from geometry alone. Compare quantitatively with Hovhannisyan's "
    "experimental peak positions. Compile the thesis manuscript and prepare the defence.")


# ══════════════════════════════════════════════════════════════════════════
heading(doc, "2.  Literature Review")
# ══════════════════════════════════════════════════════════════════════════

body(doc,
    "The thesis is grounded in three publications from the Hovhannisyan–Makaryan research "
    "group at the Institute of Physics, Yerevan State University:")

body(doc,
    "Paper 1 — T. Hovhannisyan, A. Makaryan, 'Investigation of Nerve Fiber Network in THz "
    "Spectrum Range', Armenian Journal of Physics, 2024, vol. 17, issue 1, pp. 17–20 "
    "(DOI: 10.54503/18291171-2024.17.1-17). This is the primary experimental reference. "
    "70 µm sections of white matter from pig and sheep spinal cords were placed in PTFE "
    "cuvettes with tungsten electrodes and interrogated by a Menlo Systems LAC 1550 THz "
    "time-domain spectrometer (0.2–4 THz). Key finding: resonant absorption peaks appear at "
    "0.6 THz and 2 THz exclusively when (a) the THz E-field is polarised parallel to the "
    "nerve fibres and (b) a DC voltage > ~50 V/cm is applied perpendicular to the fibres. "
    "No spectral change is observed for any other electrode/polarisation configuration.")

body(doc,
    "Paper 2 — T.N. Hovhannisyan, 'Nonlinear Properties of Biological Nerve Fiber System in "
    "THz Spectrum Range', Reports of NAS of Armenia (Doklady), 2024, vol. 124, no. 3–4 "
    "(DOI: 10.54503/0321-1339-2024.124.3-4-7). Documents the time evolution of the effect: "
    "after voltage is applied, the resonant frequency drifts toward lower values over ~5 minutes, "
    "and the high-frequency spectral components grow significantly. The effect is reversible "
    "within ~1 hour after voltage removal. Mathematical modelling via polynomial nonlinear "
    "transformation provides qualitative agreement with the observed waveform distortion.")

body(doc,
    "Paper 3 — T. Hovhannisyan, A. Makaryan, R. Martirosyan, 'Investigation of Nonlinear "
    "Properties of Animal Spinal Cord in THz Range', IRPhE 2024, Yerevan, September 2024 "
    "(DOI: 10.1049/icp.2025.0984). Proposes the physical mechanism: voltage above a threshold "
    "opens ion channels at the Nodes of Ranvier, converting the fibres into wires with finite "
    "conductivity and the whole sample into a periodic diffraction grating. Analytical modelling "
    "using the Tretyakov surface-admittance wire-array model reproduces the experimental "
    "absorption features qualitatively.")

body(doc,
    "Supporting literature covers: Hodgkin–Huxley nerve fiber electrophysiology; double-Debye "
    "dielectric relaxation model for liquid water in the THz band (two relaxation times: "
    "τ₁ = 9.36 ps, τ₂ = 0.30 ps, ε_inf = 3.17); COMSOL Multiphysics EWFD physics formulation "
    "and scattered-field approach; Tretyakov, Analytical Modelling in Applied Electromagnetics "
    "(Artech House, 2003), ch. 4 on wire-grid filters; and prior computational THz-nerve "
    "waveguide studies predicting that myelinated fibres can support guided THz modes.")


# ══════════════════════════════════════════════════════════════════════════
heading(doc, "3.  Research Plan")
# ══════════════════════════════════════════════════════════════════════════

body(doc,
    "Central hypothesis: the resonant THz absorption peaks at 0.6 THz and 2 THz reported by "
    "Hovhannisyan & Makaryan (2024) can be reproduced and explained from first principles by "
    "a 3D finite-element electromagnetic simulation of a single myelinated nerve fiber modelled "
    "as one unit cell of a periodic diffraction grating.")

body(doc,
    "Objective 1 — Reproduce the 0.6 THz resonance. Run a 26-point frequency sweep "
    "(0.1–2.0 THz, σ = 0) and verify that the peak field-at-node position matches the "
    "experimental value to within the frequency sampling resolution. Target: peak at "
    "0.6 ± 0.1 THz. Result: peak at f = 0.632 THz (Sim 1); confirmed at f₀ = 0.619 ± 0.008 THz "
    "with Q = 15 by dense sweep with refined mesh (Sim 22).")

body(doc,
    "Objective 2 — Reproduce the 2 THz resonance. Identify a second spectral feature near "
    "2 THz. Result: Sim 22 resolves it at f₀ = 1.938 ± 0.010 THz with Q = 30.5.")

body(doc,
    "Objective 3 — Explain the resonances from geometry. Predict both peak frequencies using "
    "the half-wave standing-mode formula λ_water(f)/2 = segment length, evaluated at the "
    "Debye water refractive index. Predicted: 0.589 THz (internode + half-node, 120 µm) "
    "and 2.043 THz (node, 40 µm) — both within 5 % of the FEM values (Sim 24).")

body(doc,
    "Objective 4 — Reproduce the polarisation dependence. Compare E ⊥ fibre vs E ∥ fibre "
    "configurations. Result: rotating to parallel raises the node-annulus field by ~20–25 %, "
    "consistent with Paper 1's observation that absorption appears only in the parallel "
    "orientation (Sims 18, 20, 21).")

body(doc,
    "Objective 5 — Assess the role of node conductivity σ. Sweep σ across 11 orders of "
    "magnitude at f = 0.6 THz. Result: the correct COMSOL encoding (literal complex "
    "permittivity) reveals the textbook three-regime behaviour — linear enhancement at low σ, "
    "peak coupling at σ ≈ 100 S/m, metallic suppression at high σ (Sim 2). The per-domain "
    "power dissipation at σ = 100 S/m shows the node contributes ~2 % of total absorption, "
    "a modulation detectable as a differential signal between voltage-off and voltage-on states "
    "(Sim 23).")

body(doc,
    "Objective 6 — Establish numerical uncertainty. Mesh convergence studies (Sims 16, 26) "
    "show peak amplitudes vary by ±10 % across coarse/baseline/fine meshes; peak positions "
    "are converged to better than 0.5 %.")


# ══════════════════════════════════════════════════════════════════════════
heading(doc, "4.  Research Methods")
# ══════════════════════════════════════════════════════════════════════════

body(doc,
    "Primary numerical method: Finite Element Method (FEM) implemented in COMSOL Multiphysics "
    "using the Electromagnetic Waves, Frequency Domain (EWFD) physics interface. FEM was "
    "chosen over FDTD because COMSOL's EWFD natively supports the scattered-field formulation, "
    "Floquet-periodic lateral boundary conditions, and arbitrary frequency-dependent material "
    "expressions — all required for this problem. FDTD would require a custom periodic-BC "
    "implementation and is less well-suited to frequency-domain sweeps with dispersive media.")

body(doc,
    "Physics formulation. The EWFD solve is set up in scattered-field mode: the background "
    "incident field E_bg = exp(−ik₀z) x̂ is subtracted analytically, and the FEM computes only "
    "the scattered field. Floquet-periodic BCs on the ±x and ±y lateral faces enforce the "
    "diffraction-grating periodicity (unit cell 40 µm × 40 µm in cross-section). Scattering "
    "BCs on the z = 0 and z = L end-caps absorb outgoing waves. The MUMPS direct solver is "
    "used throughout (LU factorisation); each frequency-domain solve takes 10–20 s on the "
    "simulation machine.")

body(doc,
    "Material models. Water (axon interior and external domain): double-Debye permittivity "
    "ε(f) = 3.17 + 73.5/(1 + i·2πf·9.36 ps) + 1.73/(1 + i·2πf·0.30 ps). Myelin sheaths: "
    "constant complex permittivity ε = 4.5 − 0.5i. Node of Ranvier: "
    "ε_node = ε_water(f) + iσ/(2πf·ε₀), with σ baked in as a literal complex value at each "
    "simulation frequency (required by COMSOL EWFD — analytic σ/(ωε₀) expressions referencing "
    "the 'freq' variable evaluate to zero and are silently ignored).")

body(doc,
    "Mesh. Free tetrahedral mesh, maximum element size 30 µm globally, refined to 5 µm at "
    "the node of Ranvier. Baseline element count ~26 k; fine mesh ~41 k. Element size was "
    "chosen to give at least 3 elements per λ_water/2 at 2 THz (λ_water ≈ 75 µm at 2 THz, "
    "so 30 µm is close to the resolution limit — hence the dedicated convergence study).")

body(doc,
    "Software pipeline. COMSOL is driven programmatically from Python using the mph library "
    "(Python binding for the COMSOL Java API). The model geometry, materials, mesh, and study "
    "are defined entirely in Python and rebuilt from scratch for each run — no hand-edited "
    ".mph files. Simulation parameters are declared in YAML experiment files and validated "
    "by Pydantic schemas; a single run.py dispatcher loads any YAML, expands the parameter "
    "grid, and drives the sweep. Results are stored as HDF5 files (one per simulation) "
    "containing scalar metrics, 1D axial profiles |E|(z), and 2D axial slices |E|(x,z), "
    "with provenance metadata (git SHA, COMSOL version, config hash) embedded in each file. "
    "A central manifest.csv indexes all runs for fast filtering.")

body(doc,
    "Analytical verification. Two analytical models cross-check the FEM results. (1) The "
    "Tretyakov wire-array equivalent-sheet model (Sim 6) confirms the grating produces "
    "THz-band absorption with the expected band-pass character as σ increases, providing a "
    "textbook reference for the FEM diffraction-grating setup. (2) The half-wave standing-mode "
    "model (Sim 24) predicts both resonance frequencies to within 5 % from the known fiber "
    "geometry and the Debye water refractive index alone, providing a physical explanation "
    "for the peak positions without any free parameters.")

body(doc,
    "Statistical uncertainty. Lorentzian fits (Sim 19, Sim 22) extract peak frequency f₀, "
    "FWHM linewidth γ, and quality factor Q = f₀/γ from the simulated spectra. "
    "Fit uncertainties are reported at 1σ. Mesh-convergence bounds (±10 % on amplitude, "
    "< 0.5 % on frequency) are quoted alongside all peak values.")


# ══════════════════════════════════════════════════════════════════════════
heading(doc, "Figures")
# ══════════════════════════════════════════════════════════════════════════

for fname, cap in FIGURES:
    figure(doc, fname, cap)

doc.save(OUT_PATH)
print(f"Saved: {OUT_PATH}")
