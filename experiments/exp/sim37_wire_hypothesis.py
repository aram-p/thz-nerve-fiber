"""Sim 37 - order-of-magnitude tests of the experimenters' own hypothesis (Hovhannisyan, Makaryan,
Martirosyan 2024, IRPhE; Armenian J. Phys. 2024): above a DC threshold the node channels open,
fibres become "wires with nonlinear conductivity" and the slice a "diffraction grating" of wires,
which generates combination frequencies.

Each claim is tested with textbook constants; the full-wave version of (2) is Sim 36 'axon'.
  (1) Do open channels change THz current paths?  membrane capacitive admittance vs channel conductance
  (2) What axial conductivity would 9 dB need?   Wiener-parallel (exact for E || fibres, quasi-static)
  (3) When does a fibre act as a wire (inductive)? resistance vs inductance per unit length
  (4) Can the slice diffract?                     fibre spacing vs wavelength in tissue
  (5) Can anything respond nonlinearly in 1 ps?    gating time, ion drift, Wien-effect field
  (6) What else does 100 V DC do?                 current, Joule heat, electrolysis gas
  (7) Is a 9 dB change at 2 THz measurable?       reported baseline vs pure water

Writes results/sim37/sim37.txt
"""
from pathlib import Path

import numpy as np

from thznerve.physics.effmedium import tissue_eps
from thznerve.physics.materials import EPS0, ETA0, MU0, WATER_KS1996, axoplasm, myelin, refractive_index

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim37"; OUT.mkdir(parents=True, exist_ok=True)
C0 = 299792458.0
lines = []


def say(s=""):
    print(s, flush=True); lines.append(s)


f = np.array([0.6, 2.0]) * 1e12
w = 2 * np.pi * f
G, P, B_UM = 0.7, 0.6, 5.0
a = G * B_UM * 1e-6

# (1) membrane
Cm = 1e-2          # F/m^2 (1 uF/cm^2; MRG nodal 2 uF/cm^2)
gNa = 3e4          # S/m^2 (MRG node gNa 3 S/cm^2, all channels open)
say("(1) Node membrane at THz")
for fi, wi in zip(f, w):
    say(f"  {fi/1e12:.1f} THz: capacitive admittance wCm = {wi*Cm:.2e} S/m^2 vs all Na channels open {gNa:.0e} S/m^2"
        f" -> ratio {wi*Cm/gNa:.1e}")
say("  The membrane is already a short circuit for THz currents; opening channels cannot change THz current paths.")

# (2) axial conductivity for 9 dB in 100 um
e_par, _ = tissue_eps(axoplasm()(f), myelin()(f), WATER_KS1996(f), G, P)
n_t = refractive_index(e_par)
d_alpha = 9 / (10 / np.log(10)) / 100e-6                      # 1/m, power
ds_tissue = d_alpha * n_t.real / ETA0                         # small-change limit: d_alpha = d_sigma eta0 / n
f_ax = P * G ** 2
say("\n(2) Extra conductivity needed for a 9 dB change of a 100 um slice (E || fibres)")
for i, fi in enumerate(f):
    eff = w[i] * EPS0 * (-axoplasm()(f[i]).imag)
    say(f"  {fi/1e12:.1f} THz: tissue +{ds_tissue[i]:.0f} S/m; in the axoplasm alone (volume share {f_ax:.2f}) "
        f"+{ds_tissue[i]/f_ax:.0f} S/m; axoplasm now: ionic 1.4 S/m, water loss already w*eps0*eps'' = {eff:.0f} S/m")
say("  Physiological axoplasm 1.4 S/m; saturated NaCl ~25 S/m. Ion channels change the radial membrane path, not the")
say("  axial conductivity of the axoplasm, which is continuous through the node in any state.")

# (3) wire regime: resistance vs inductance per unit length
Lam = np.sqrt(np.pi * (B_UM * 1e-6) ** 2 / P)                 # centre spacing, square packing
L = MU0 / (2 * np.pi) * np.log(Lam / a)
say(f"\n(3) Wire regime: spacing {Lam*1e6:.1f} um, axon radius {a*1e6:.1f} um, inductance {L:.2e} H/m")
for fi, wi in zip(f, w):
    s_x = 1 / (wi * L * np.pi * a ** 2)
    say(f"  {fi/1e12:.1f} THz: inductance matters only above sigma ~ {s_x:.1e} S/m "
        f"({s_x/1.4:.0e} x physiological); below it the fibres just mix in (Wiener), no wire-grid physics")

# (4) diffraction
say("\n(4) Grating orders")
for b_um in (1.0, 5.0, 10.0):
    lam_c = np.sqrt(np.pi * (b_um * 1e-6) ** 2 / P)
    say(f"  fibre radius {b_um:4.1f} um: spacing {lam_c*1e6:5.1f} um; first order propagates in tissue (n~2.1) only above "
        f"{C0/(lam_c*2.1)/1e12:5.1f} THz, in air above {C0/lam_c/1e12:5.1f} THz")
say("  Fibres in white matter are not periodic and their spacing is far below the wavelength: no diffraction orders in 0.1-3 THz.")

# (5) nonlinearity
E_thz = np.array([1e3, 1e4])       # V/m: 10-100 V/cm, typical fibre-coupled photoconductive TDS at focus
mu_na = 5.2e-8                     # m^2/(V s)
say("\n(5) Nonlinear response to a ~2 ps pulse")
say("  Na-channel gating: activation ~0.1 ms, fastest kinetics ~10 us: 1e7 x slower than the pulse.")
for E in E_thz:
    say(f"  E = {E:.0e} V/m: ion drift in 1 ps = {mu_na*E*1e-12:.0e} m; field / Wien-effect onset (1e7 V/m) = {E/1e7:.0e}")
say("  The fitted J = a1 E + a2 E^2 + ... with a2/a1 = 7.5 means the quadratic current equals the linear one at 13 % of")
say("  the pulse peak; bulk electrolytes and water need >~1e7-1e8 V/m for any measurable THz nonlinearity.")
# the DC bias breaks inversion symmetry: chi2_eff = 3 chi3 E_dc. Water's static Kerr constant
# B = 9.4e-14 m/V^2 (589 nm, orientational) gives chi3 ~ 8 n lambda B / 3 ~ 2e-19 m^2/V^2
chi3 = 8 * 1.33 * 589e-9 * 9.4e-14 / 3
E_dc = 100 / 15e-3
say(f"  The DC bias breaks inversion symmetry: chi2_eff E_THz / chi1 = 3 chi3 E_dc E_THz / chi1 <= "
    f"{3 * chi3 * E_dc * 1e4 / 3.0:.0e} (water chi3 ~ {chi3:.0e} m^2/V^2 from its Kerr constant).")
say("  The tissue is post-mortem (slaughterhouse); its channels are unlikely to gate at all.")
say(f"  Their '1 V per fibre' estimate: the bias drops only ~2 E r = {2 * E_dc * 5e-6 * 1e3:.0f} mV across a 10 um fibre.")

# (6) DC side effects
sig_perp = 0.083                   # S/m, white matter across fibres (Ranck & BeMent 1965)
F_C = 96485.0; Vm = 0.0245         # C/mol; m^3/mol gas at 25 C
say("\n(6) 100 V DC across the slice (field across the fibres)")
for gap_mm, t_um, name in ((18, 70, "paper 1"), (15, 100, "paper 3")):
    E = 100 / (gap_mm * 1e-3)
    J = sig_perp * E
    A = t_um * 1e-6 * 20e-3        # slice thickness x 20 mm electrode length (paper 3: 20x20 mm sample)
    I = J * A
    q = sig_perp * E ** 2
    gas = I / (2 * F_C) * Vm * 1.5 * 300 * 1e9   # H2 + O2 (1.5 x H2), mm^3 in 5 min
    say(f"  {name}: E = {E:.0f} V/m, J = {J:.0f} A/m^2, I ~ {I*1e3:.2f} mA, Joule {q/1e6:.1f} MW/m^3 "
        f"({q/4.0e6:.2f} K/s), electrolysis gas ~{gas:.0f} mm^3 in 5 min (sample ~{20*20*t_um/1e3:.0f} mm^3)")
say("  Current along the fibres (sigma_par ~ 0.6 S/m) would heat ~7x more; the papers report no effect for that")
say("  direction, so heating alone does not explain the DC-direction selectivity either.")

# (7) baseline
aw = 2 * w * (-refractive_index(WATER_KS1996(f)).imag) / C0
say("\n(7) Reported baseline vs pure water, 100 um")
for i, fi in enumerate(f):
    say(f"  {fi/1e12:.1f} THz: pure water {10/np.log(10)*aw[i]*100e-6:.1f} dB")
from thznerve.physics.materials import WATER_BROADBAND
ab = 2 * w * (-refractive_index(WATER_BROADBAND(f)).imag) / C0
say(f"  broadband water (Liebe Lorentz bands): {10/np.log(10)*ab[1]*100e-6:.1f} dB at 2 THz")
say("  Reported: ~5 dB at 0.25 THz rising to ~30 dB at 2 THz for 0.1 mm (paper 3): 1.9-2.6x pure water at 2 THz.")
say("  Measured nerve tissue absorbs within ~10 % of water (Berry 2003: 246 vs 225 /cm), so the slice is optically")
say("  thicker than stated or carries extra loss (scatter, reflection); see Sim 38.")

(OUT / "sim37.txt").write_text("\n".join(lines), encoding="utf-8")
