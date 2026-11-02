"""Sim 36 - checks requested by the independent hostile review.

modes (one COMSOL process per mode, run one after another):
  sigma   node conductivity 0..3e6 S/m at Lint = 200 um, 5 frequencies, both polarisations, with
          three tissue observables: (a) extra absorption in the fibre, (b) extra absorption inside
          the flux contour r = Rf (fibre + surrounding water), (c) change of the fibre's induced
          dipole -> change of the tissue's effective permittivity (coherent transmission).
          Plus the worst d_sigma at Lint = 1000 um (internode independence at high d_sigma).
  dense   16 frequencies 0.1-3 THz at d_sigma = 3e4 and 1e6 S/m, b = 5 um and b = 10.7 um
          (a resonance of a conducting node ring would need a dense sweep to be seen).
  aniso   lamellar (uniaxial Wiener) myelin, and lamellar myelin with confined interlamellar water
          (water eps'' halved), base fibre, 10 frequencies.
  epsnode node material with arbitrary passive eps = eps' - j eps'' (eps' from -30 to 30, incl.
          negative values that allow a plasmon-like ring resonance), 0.6 and 2 THz.

Tissue conversion: one fibre period (one node) occupies V_cell = V_node / f_n of tissue, with
f_n = p (1 - g^2) Ln / (Lint + Ln). A change dp of the fibre dipole (|E0| = 1 V/m) changes the
tissue permittivity by d_eps = dp / (eps0-free units) / V_cell, and the power attenuation by
d_alpha = -k0 Im(d_eps / n_t). The dB change of a 100 um slab scales as p (1 - g^2) Ln / Lint;
results are given for Lint = 1000 um, p = 0.5 and rescaled in the summary.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

from thznerve.physics.effmedium import tissue_eps
from thznerve.physics.materials import (DRY_LIPID, ETA0, WATER_BROADBAND, WATER_KS1996, ConstantEps, PlusSigma, Scaled, Uniaxial,
                                        Wiener, axoplasm, myelin, myelin_lamellar, refractive_index)
from thznerve.physics.samplebound import DB_PER_NEPER_POWER, node_volume_fraction
from thznerve.physics.tmm import slab_rta

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("sim30", ROOT / "experiments/prod/sim30_node_spectrum.py")
sim30 = importlib.util.module_from_spec(spec); spec.loader.exec_module(sim30)
OUT = sim30.OUT
F5 = np.array([0.1, 0.6, 1.0, 2.0, 3.0]) * 1e12
MESH = dict(h_water=5.0, h_frac=0.25, h_far_frac=0.5)
A, B, LN = 3.5, 5.0, 1.0
D_SLAB = 100e-6
C0 = 2.99792458e8


def save(name, out):
    (OUT / f"sim36_{name}.json").write_text(json.dumps(out))


def tissue_db(r_new, r_ref, f, pol, Lint_t=1000.0, p=0.5, a=A, b=B, Ln=LN):
    """dB change of a 100 um slab from three observables (fibre absorption, absorption inside
    the flux contour, dipole). Returns dict of arrays."""
    f = np.asarray(f)
    Vn = np.pi * (b ** 2 - a ** 2) * Ln * 1e-18
    V_cell = Vn / node_volume_fraction(a, b, Ln, Lint_t, p)
    e_par, e_perp = tissue_eps(axoplasm()(f), myelin()(f), WATER_KS1996(f), a / b, p)
    n_t = refractive_index(e_par if pol == "par" else e_perp)
    I0 = 0.5 * n_t.real / ETA0
    k0 = 2 * np.pi * f / C0
    out = {}
    for key in ("P_fibre", "flux_in"):
        dP = np.array(r_new[pol][key]) - np.array(r_ref[pol][key])
        out[key] = DB_PER_NEPER_POWER * dP / (V_cell * I0) * D_SLAB
    dp = (np.array(r_new[pol]["p_re"]) + 1j * np.array(r_new[pol]["p_im"])) \
        - (np.array(r_ref[pol]["p_re"]) + 1j * np.array(r_ref[pol]["p_im"]))
    d_eps = dp / V_cell
    out["dipole"] = DB_PER_NEPER_POWER * (-k0 * np.imag(d_eps / n_t)) * D_SLAB
    out["phase_rad"] = k0 * np.real(d_eps / (2 * n_t)) * D_SLAB      # extra phase delay of the slab, first order
    return out


def sigma():
    out = {}
    for ds in (0.0, 1e-4, 1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6, 3e6):
        node = PlusSigma(WATER_KS1996, ds) if ds else WATER_KS1996
        out[f"{ds:g}"] = sim30.run("s36", tag=f"s36_sig_{ds:g}", freqs=list(F5), Lint=200.0,
                                   eps_override=dict(node=node), **MESH)
        save("sigma", out)
    for ds in (0.0, 3e4):                     # internode independence at the worst d_sigma
        node = PlusSigma(WATER_KS1996, ds) if ds else WATER_KS1996
        out[f"L1000_{ds:g}"] = sim30.run("s36", tag=f"s36_sigL1000_{ds:g}", freqs=list(F5), Lint=1000.0,
                                         eps_override=dict(node=node), **MESH)
        save("sigma", out)
    report_sigma(out)


def report_sigma(out=None):
    out = out or json.loads((OUT / "sim36_sigma.json").read_text())
    ref = out["0"]
    rows, summ = [], {}
    for k, r in out.items():
        if k.startswith("L1000"):
            continue
        for pol in ("par", "perp"):
            t = tissue_db(r, ref, F5, pol)
            for i, fi in enumerate(F5):
                rows.append([float(k), 0 if pol == "par" else 1, fi / 1e12, t["P_fibre"][i], t["flux_in"][i], t["dipole"][i]])
    rows = np.array(rows)
    np.savetxt(OUT / "sim36_sigma.csv", rows, delimiter=",",
               header="d_sigma,pol(0=par),f_THz,dB_fibre_abs,dB_flux_in,dB_dipole (100 um slab, Lint 1000 um, p 0.5)")
    for pol in (0, 1):
        for fi in (0.6, 2.0, 3.0):
            sel = rows[(rows[:, 1] == pol) & np.isclose(rows[:, 2], fi)]
            m = np.abs(sel[:, 3:]).max(axis=0)
            j = np.argmax(np.abs(sel[:, 5]))
            line = (f"{'par ' if pol == 0 else 'perp'} {fi} THz: max |dB| over d_sigma: fibre-abs {m[0]:.4f}, "
                    f"flux {m[1]:.4f}, dipole {m[2]:.4f} (dipole max at {sel[j, 0]:.0e} S/m)")
            print(line, flush=True); summ[f"{pol}_{fi}"] = m.tolist()
    # internode independence at 3e4
    for pol in ("par", "perp"):
        a = np.array(out["L1000_30000"][pol]["P_node"]) / np.array(out["30000"][pol]["P_node"]) - 1
        b = (np.array(out["L1000_30000"][pol]["p_re"]) - np.array(out["L1000_0"][pol]["p_re"])) / \
            (np.array(out["30000"][pol]["p_re"]) - np.array(out["0"][pol]["p_re"])) - 1
        print(f"Lint 1000 vs 200 at 3e4 S/m {pol}: P_node {np.round(a, 5)}, d(p_re) {np.round(b, 4)}", flush=True)
    worst = np.abs(rows[:, 3:]).max()
    scale = (1000 / 200) * (0.6 / 0.5)
    print(f"worst |dB| any observable, any pol, any f: {worst:.3f} at Lint 1000/p 0.5; "
          f"x{scale:.0f} = {worst * scale:.2f} dB at Lint 200 um, p 0.6", flush=True)
    (OUT / "sim36_sigma_summary.json").write_text(json.dumps(summ))


def dense():
    out = {}
    f16 = np.linspace(0.1, 3.0, 16) * 1e12
    for (a, b), Lint in (((3.5, 5.0), 200.0), ((7.5, 7.5 / 0.7), 200.0)):
        for ds in (0.0, 3e4, 1e6):
            if ds == 0.0 and a == 3.5:
                continue                       # base spectrum at 30 frequencies exists (Sim 30)
            node = PlusSigma(WATER_KS1996, ds) if ds else WATER_KS1996
            tag = f"s36_dense_a{a:g}_{ds:g}"
            out[tag] = sim30.run("s36", tag=tag, a=a, b=b, Lint=Lint, freqs=list(f16),
                                 eps_override=dict(node=node), **MESH)
            for pol in ("par", "perp"):
                y = np.array(out[tag][pol]["P_node"])
                d2 = np.diff(np.log(y), 2)
                print(f"{tag} {pol}: P_node/P_node(0.1) {np.round(y / y[0], 2)}; max |2nd diff log| {np.abs(d2).max():.3f}", flush=True)
            save("dense", out)


def aniso():
    out = {}
    f10 = np.linspace(0.1, 3.0, 10) * 1e12
    conf = Scaled(WATER_KS1996, im=0.5)
    variants = {
        "iso": None,
        "lamellar": myelin_lamellar(),
        "lamellar_confined": Uniaxial(Wiener(conf, DRY_LIPID, 0.4, "series"), Wiener(conf, DRY_LIPID, 0.4, "par"),
                                      Wiener(conf, DRY_LIPID, 0.4, "par")),
    }
    W = WATER_BROADBAND
    for k, m in variants.items():
        ov = {} if m is None else dict(myelin=m)
        out[k] = sim30.run("s36", tag=f"s36_aniso_{k}", freqs=list(f10), Lint=200.0, eps_override=ov, **MESH)
        save("aniso", out)
    out["water_broadband"] = sim30.run("s36", tag="s36_aniso_wbb", freqs=list(f10), Lint=200.0, **MESH,
                                       eps_override=dict(water=W, node=W, axon=axoplasm(W), myelin=myelin(W)))
    save("aniso", out)
    print("water_broadband eta par", np.round(out["water_broadband"]["par"]["eta"], 3),
          "iso", np.round(out["iso"]["par"]["eta"], 3), flush=True)
    print("water_broadband eta perp", np.round(out["water_broadband"]["perp"]["eta"], 3),
          "iso", np.round(out["iso"]["perp"]["eta"], 3), flush=True)
    for k in ("lamellar", "lamellar_confined", "water_broadband"):
        for pol in ("par", "perp"):
            rel = np.array(out[k][pol]["P_node"]) / np.array(out["iso"][pol]["P_node"])
            print(f"{k} {pol}: P_node / isotropic Looyenga = {np.round(rel, 3)}", flush=True)


def epsnode():
    out = {}
    f2 = np.array([0.6, 2.0]) * 1e12
    ref = sim30.run("s36", tag="s36_epsn_ref", freqs=list(f2), Lint=200.0, **MESH)
    out["ref"] = ref
    for er in (-30, -10, -5, -3, -2, -1, 1, 2, 5, 10, 30):
        for ei in (0.3, 3.0, 30.0):
            k = f"{er}_{ei}"
            out[k] = sim30.run("s36", tag=f"s36_epsn_{k}", freqs=list(f2), Lint=200.0,
                               eps_override=dict(node=ConstantEps(complex(er, -ei))), **MESH)
            save("epsnode", out)
    best = {}
    for k, r in out.items():
        if k == "ref":
            continue
        for pol in ("par", "perp"):
            t = tissue_db(r, ref, f2, pol)
            for i, fi in enumerate(f2):
                v = max(abs(t["flux_in"][i]), abs(t["dipole"][i]))
                key = (pol, fi / 1e12)
                if v > best.get(key, (0, ""))[0]:
                    best[key] = (v, k, r[pol]["eta"][i])
    for key, (v, k, eta) in sorted(best.items()):
        print(f"{key}: max |dB| over passive eps_node = {v:.3f} at eps_node = {k.replace('_', ' - j')} (eta {eta:.2f})", flush=True)


def axon():
    """The experimenters' hypothesis: open channels turn whole fibres into conducting wires.
    Extra conductivity in the whole axoplasm, full wave; tissue dB from the dipole and the
    absorbed power per period (V_cell = period * pi b^2 / p: the effect is per unit length)."""
    out = {}
    for ds, hf in ((0.0, 0.25), (1e1, 0.25), (1e2, 0.25), (4e2, 0.25), (1e3, 0.25), (1e4, 0.25),
                   (3e4, 0.25), (1e5, 0.25), (1e5, 0.125)):
        ax = PlusSigma(axoplasm(), ds) if ds else axoplasm()
        k = f"{ds:g}" + ("" if hf == 0.25 else f"_h{hf}")
        out[k] = sim30.run("s36", tag=f"s36_axon_{k}", freqs=list(F5), Lint=200.0, eps_override=dict(axon=ax),
                           **(MESH | dict(h_frac=hf)))
        save("axon", out)
    report_axon(out)


def axon2():
    """Review round 4: bracket the 9 dB threshold (interpolated between 400 and 1000 S/m) with
    full-wave points at 600 and 800 S/m, then give the threshold against packing."""
    out = json.loads((OUT / "sim36_axon.json").read_text())
    for ds in (6e2, 8e2):
        k = f"{ds:g}"
        if k not in out:
            out[k] = sim30.run("s36", tag=f"s36_axon_{k}", freqs=list(F5), Lint=200.0,
                               eps_override=dict(axon=PlusSigma(axoplasm(), ds)), **(MESH | dict(h_frac=0.25)))
            save("axon", out)
    report_axon(out)
    axon_threshold(out)


def axon_threshold(out=None):
    """Extra axoplasm conductivity for a 9 dB change of a 100 um slab (E || fibres, exact slab on
    PTFE), interpolated in log(sigma) between full-wave points, for fibre packing 0.3/0.5/0.7."""
    out = out or json.loads((OUT / "sim36_axon.json").read_text())
    ref = out["0"]
    keys = sorted((k for k in out if k != "0" and "_" not in k), key=float)
    sig = np.array([float(k) for k in keys])
    for p in (0.3, 0.5, 0.7):
        V_cell = (200.0 + LN) * np.pi * B ** 2 * 1e-18 / p
        e_par, _ = tissue_eps(axoplasm()(F5), myelin()(F5), WATER_KS1996(F5), A / B, p)
        db = np.zeros((len(keys), len(F5)))
        for j, k in enumerate(keys):
            dp = (np.array(out[k]["par"]["p_re"]) + 1j * np.array(out[k]["par"]["p_im"])) \
                - (np.array(ref["par"]["p_re"]) + 1j * np.array(ref["par"]["p_im"]))
            for i, fi in enumerate(F5):
                db[j, i] = -10 * np.log10(slab_rta(fi, 1.0, [(complex(e_par[i] + dp[i] / V_cell), D_SLAB)], 1.43)[1]
                                          / slab_rta(fi, 1.0, [(complex(e_par[i]), D_SLAB)], 1.43)[1])
        thr = [float(np.exp(np.interp(9.0, db[:, i], np.log(sig)))) if db[:, i].max() > 9 else np.nan
               for i in range(len(F5))]
        print(f"packing {p}: dB at sigma {sig.tolist()}:", flush=True)
        for i, fi in enumerate(F5):
            print(f"   {fi/1e12:.1f} THz: {np.round(db[:, i], 2).tolist()}  -> 9 dB at {thr[i]:.0f} S/m", flush=True)


def report_axon(out=None, p=0.5):
    out = out or json.loads((OUT / "sim36_axon.json").read_text())
    ref = out["0"]
    V_cell = (200.0 + LN) * np.pi * B ** 2 * 1e-18 / p
    k0 = 2 * np.pi * F5 / C0
    rows = []
    for k, r in out.items():
        if k == "0":
            continue
        for pol in ("par", "perp"):
            e_par, e_perp = tissue_eps(axoplasm()(F5), myelin()(F5), WATER_KS1996(F5), A / B, p)
            n_t = refractive_index(e_par if pol == "par" else e_perp)
            I0 = 0.5 * n_t.real / ETA0
            dP = np.array(r[pol]["P_fibre"]) - np.array(ref[pol]["P_fibre"])
            dp = (np.array(r[pol]["p_re"]) + 1j * np.array(r[pol]["p_im"])) - (np.array(ref[pol]["p_re"]) + 1j * np.array(ref[pol]["p_im"]))
            db_abs = DB_PER_NEPER_POWER * dP / (V_cell * I0) * D_SLAB
            db_dip = DB_PER_NEPER_POWER * (-k0 * np.imag(dp / V_cell / n_t)) * D_SLAB
            # Wiener parallel (E || fibres): d_eps_tissue = f_ax * d_eps_ax exactly in the quasi-static limit
            f_ax = p * (A / B) ** 2
            d_eps_ax = PlusSigma(axoplasm(), float(k.split("_")[0]))(F5) - axoplasm()(F5)
            db_wiener = DB_PER_NEPER_POWER * (-k0 * np.imag(f_ax * d_eps_ax / n_t)) * D_SLAB
            # exact 100 um slab on PTFE (as Sim 34) with eps_t + d_eps from the dipole: includes the
            # change of reflection, which the first-order -k0 Im(d_eps/n) leaves out
            e_t = e_par if pol == "par" else e_perp
            db_slab = np.array([10 * np.log10(slab_rta(fi, 1.0, [(complex(e_t[i] + dp[i] / V_cell), D_SLAB)], 1.43)[1]
                                              / slab_rta(fi, 1.0, [(complex(e_t[i]), D_SLAB)], 1.43)[1])
                                for i, fi in enumerate(F5)])
            print(f"   exact slab dB {np.round(-db_slab, 3)}", flush=True)
            for i, fi in enumerate(F5):
                rows.append([float(k.split("_")[0]), "_h" in k, 0 if pol == "par" else 1, fi / 1e12,
                             db_abs[i], db_dip[i], db_wiener[i] if pol == "par" else np.nan, -db_slab[i]])
            print(f"axon +{k} S/m {pol}: dB(abs) {np.round(db_abs, 3)} dB(dipole) {np.round(db_dip, 3)}"
                  + (f" Wiener {np.round(db_wiener, 3)}" if pol == "par" else ""), flush=True)
    np.savetxt(OUT / "sim36_axon.csv", np.array(rows, dtype=float), delimiter=",",
               header="d_sigma_axon,fine_mesh,pol(0=par),f_THz,dB_abs,dB_dipole,dB_wiener_par,dB_exact_slab (100 um slab, p 0.5)")


def epsfine():
    """Refine the E-par node resonance found near eps_node = -2 (water loss should cap it),
    with a mesh check at the strongest point."""
    out = json.loads((OUT / "sim36_epsnode.json").read_text())
    f2 = np.array([0.6, 2.0]) * 1e12
    for er in (-2.5, -2.0, -1.6, -1.3):
        for ei in (0.01, 0.1):
            k = f"{er}_{ei}"
            out[k] = sim30.run("s36", tag=f"s36_epsn_{k}", freqs=list(f2), Lint=200.0,
                               eps_override=dict(node=ConstantEps(complex(er, -ei))), **MESH)
            (OUT / "sim36_epsfine.json").write_text(json.dumps(out))
    best = max((k for k in out if k != "ref"), key=lambda k: max(out[k]["par"]["eta"]))
    er, ei = map(float, best.split("_"))
    out[best + "_h"] = sim30.run("s36", tag=f"s36_epsn_{best}_h", freqs=list(f2), Lint=200.0,
                                 eps_override=dict(node=ConstantEps(complex(er, -ei))), **(MESH | dict(h_frac=0.125)))
    (OUT / "sim36_epsfine.json").write_text(json.dumps(out))
    ref = out["ref"]
    rows = []
    for k, r in out.items():
        if k == "ref":
            continue
        for pol in ("par", "perp"):
            t = tissue_db(r, ref, f2, pol)
            for i in range(2):
                rows.append((max(abs(t["flux_in"][i]), abs(t["dipole"][i]), abs(t["P_fibre"][i])), k, pol, f2[i] / 1e12, r[pol]["eta"][i]))
    rows.sort(reverse=True)
    for v, k, pol, fi, eta in rows[:8]:
        print(f"eps_node {k}: {pol} {fi} THz max |dB| {v:.3f} (eta {eta:.1f})", flush=True)
    print(f"mesh check {best}: eta {out[best]['par']['eta']} vs fine {out[best + '_h']['par']['eta']}", flush=True)


def epsdense():
    """Review round 5: the passive-eps bound was sampled only at 0.6/2 THz and with the ring
    resonance between grid points (-1.6, -1.3). Dense eps' grid at eps'' = 0.1 (the worst case the
    claim covers: the peak falls as eps'' grows) at 0.6/1/2/3 THz; parabola through the top three."""
    f4 = np.array([0.6, 1.0, 2.0, 3.0]) * 1e12
    path = OUT / "sim36_epsdense.json"
    out = json.loads(path.read_text()) if path.exists() else {}
    if "ref" not in out:
        out["ref"] = sim30.run("s36", tag="s36_epsd_ref", freqs=list(f4), Lint=200.0, **MESH)
    for er in (-2.0, -1.8, -1.6, -1.5, -1.45, -1.4, -1.35, -1.3, -1.2, -1.0):
        k = f"{er}_0.1"
        if k not in out:
            out[k] = sim30.run("s36", tag=f"s36_epsd_{k}", freqs=list(f4), Lint=200.0,
                               eps_override=dict(node=ConstantEps(complex(er, -0.1))), **MESH)
            path.write_text(json.dumps(out))
    report_epsdense(out)


def epsdense_mesh():
    """Review round 6: mesh check of the passive-eps bound where it peaks (3 THz, eps' -1.45/-1.5)."""
    path = OUT / "sim36_epsdense.json"
    out = json.loads(path.read_text())
    f4 = np.array([0.6, 1.0, 2.0, 3.0]) * 1e12
    n145 = dict(eps_override=dict(node=ConstantEps(complex(-1.45, -0.1))))
    for key, kw, hf in (("ref_h", {}, 0.125), ("-1.45_0.1_h", n145, 0.125),
                        ("-1.5_0.1_h", dict(eps_override=dict(node=ConstantEps(complex(-1.5, -0.1)))), 0.125),
                        ("ref_q", {}, 0.0625), ("-1.45_0.1_q", n145, 0.0625)):
        if key not in out:
            out[key] = sim30.run("s36", tag=f"s36_epsd_{key}", freqs=list(f4), Lint=200.0, **kw,
                                 **(MESH | dict(h_frac=hf)))
            path.write_text(json.dumps(out))
    for k in ("-1.45_0.1", "-1.5_0.1"):
        for pol in ("par", "perp"):
            base = tissue_db(out[k], out["ref"], f4, pol)
            fine = tissue_db(out[k + "_h"], out["ref_h"], f4, pol)
            b = np.max(np.abs([base[o] for o in ("P_fibre", "flux_in", "dipole")]), axis=0)
            fn = np.max(np.abs([fine[o] for o in ("P_fibre", "flux_in", "dipole")]), axis=0)
            print(f"{k} {pol}: max |dB| base {np.round(b, 4)} fine {np.round(fn, 4)} change {np.round(fn / b - 1, 3)}; "
                  f"eta base {np.round(out[k][pol]['eta'], 1)} fine {np.round(out[k + '_h'][pol]['eta'], 1)}", flush=True)
    for pol in ("par", "perp"):
        lv = [tissue_db(out["-1.45_0.1" + s_], out["ref" + s_], f4, pol) for s_ in ("", "_h", "_q")]
        db = [np.max(np.abs([d[o] for o in ("P_fibre", "flux_in", "dipole")]), axis=0) for d in lv]
        print(f"-1.45_0.1 {pol} h_frac 0.25/0.125/0.0625: max |dB| {np.round(db[0], 4)} / {np.round(db[1], 4)} / "
              f"{np.round(db[2], 4)}; eta {np.round(out['-1.45_0.1_q'][pol]['eta'], 1)} at 0.0625", flush=True)


def epsdense_hscan():
    """Review round 8: refinement lowered eps' -1.5 more than -1.45, so the peak may move. Rescan
    eps' at h_frac 0.125 (the peak is converging downwards; see epsdense_mesh) and compare the peak."""
    path = OUT / "sim36_epsdense.json"
    out = json.loads(path.read_text())
    f4 = np.array([0.6, 1.0, 2.0, 3.0]) * 1e12
    for er in (-1.2, -1.3, -1.35, -1.4, -1.6, -1.8):
        key = f"{er}_0.1_h"
        if key not in out:
            out[key] = sim30.run("s36", tag=f"s36_epsd_{key}", freqs=list(f4), Lint=200.0,
                                 eps_override=dict(node=ConstantEps(complex(er, -0.1))), **(MESH | dict(h_frac=0.125)))
            path.write_text(json.dumps(out))
    report_epsdense_h(out)


def report_epsdense_h(out=None):
    out = out or json.loads((OUT / "sim36_epsdense.json").read_text())
    f4 = np.array([0.6, 1.0, 2.0, 3.0]) * 1e12
    keys = sorted((k for k in out if k.endswith("_0.1_h")), key=lambda k: float(k.split("_")[0]))
    er = np.array([float(k.split("_")[0]) for k in keys])
    for pol in ("par",):
        for lvl, ks, ref in (("h_frac 0.25", [k[:-2] for k in keys], out["ref"]), ("h_frac 0.125", keys, out["ref_h"])):
            t = [tissue_db(out[k], ref, f4, pol) for k in ks]
            db = np.array([np.max(np.abs([d[o] for o in ("P_fibre", "flux_in", "dipole")]), axis=0) for d in t])
            sgn = np.array([[d[o][i] for o in ("P_fibre", "flux_in", "dipole")] for d in t for i in range(4)]).reshape(len(ks), 4, 3)
            for i, fi in enumerate(f4):
                j = int(np.argmax(db[:, i]))
                print(f"{pol} {lvl} {fi/1e12:.1f} THz: peak |dB| {db[j, i]:.4f} at eps' {er[j]} "
                      f"(P_fibre, flux_in, dipole signed: {np.round(sgn[j, i], 4)}); all {np.round(db[:, i], 4)}", flush=True)


def round9():
    """Review round 9: (a) the full passive-eps grid at 3 THz (E-par), signed, to find the largest
    transmission rise any node state gives; (b) eps' -2.0 on the 2x finer mesh (the 0.6 THz fine
    peak sat at the scan edge); (c) internode 50 um at the eps peak and in water, to test the
    packing/internode scaling the combined bound extrapolates to."""
    f1 = [3.0e12]
    path = OUT / "sim36_round9.json"
    out = json.loads(path.read_text()) if path.exists() else {}
    if "g3_ref" not in out:
        out["g3_ref"] = sim30.run("s36", tag="s36_r9_g3_ref", freqs=f1, Lint=200.0, pols=("par",), **MESH)
    for er in (-30, -10, -5, -3, -2, -1, 1, 2, 5, 10, 30):
        for ei in (0.3, 3.0, 30.0):
            k = f"g3_{er}_{ei}"
            if k not in out:
                out[k] = sim30.run("s36", tag=f"s36_r9_{k}", freqs=f1, Lint=200.0, pols=("par",),
                                   eps_override=dict(node=ConstantEps(complex(er, -ei))), **MESH)
                path.write_text(json.dumps(out))
    dense = json.loads((OUT / "sim36_epsdense.json").read_text())
    f4 = [0.6e12, 1.0e12, 2.0e12, 3.0e12]
    if "-2.0_0.1_h" not in dense:
        dense["-2.0_0.1_h"] = sim30.run("s36", tag="s36_epsd_-2.0_0.1_h", freqs=f4, Lint=200.0,
                                        eps_override=dict(node=ConstantEps(complex(-2.0, -0.1))),
                                        **(MESH | dict(h_frac=0.125)))
        (OUT / "sim36_epsdense.json").write_text(json.dumps(dense))
    for k, ov in (("L50_ref", None), ("L50_-1.45_0.1", dict(node=ConstantEps(complex(-1.45, -0.1))))):
        if k not in out:
            out[k] = sim30.run("s36", tag=f"s36_r9_{k}", freqs=f4, Lint=50.0, eps_override=ov, **MESH)
            path.write_text(json.dumps(out))
    path.write_text(json.dumps(out))
    report_round9(out)


def report_round9(out=None):
    out = out or json.loads((OUT / "sim36_round9.json").read_text())
    f1 = np.array([3.0e12]); f4 = np.array([0.6, 1.0, 2.0, 3.0]) * 1e12
    obs = ("P_fibre", "flux_in", "dipole")
    rows = []
    for k, r in out.items():
        if k.startswith("g3_") and k != "g3_ref":
            t = tissue_db(r, out["g3_ref"], f1, "par")
            rows.append((k[3:], [t[o][0] for o in obs]))
    rise = min(rows, key=lambda x: min(x[1])); loss = max(rows, key=lambda x: max(x[1]))
    print(f"(a) 3 THz E-par, full eps grid ({len(rows)} states), signed dB per 100 um at internode 1000 um, packing 0.5 "
          f"(positive = more loss): largest rise {min(rise[1]):+.4f} dB at eps {rise[0]} {np.round(rise[1], 4)}; "
          f"largest loss {max(loss[1]):+.4f} dB at eps {loss[0]}", flush=True)
    dense = json.loads((OUT / "sim36_epsdense.json").read_text())
    for pol in ("par",):
        t = tissue_db(dense["-2.0_0.1_h"], dense["ref_h"], f4, pol)
        print(f"(b) eps' -2.0 fine mesh {pol}: max |dB| {np.round(np.max(np.abs([t[o] for o in obs]), axis=0), 4)} at 0.6/1/2/3 THz", flush=True)
    for pol in ("par", "perp"):
        t50 = tissue_db(out["L50_-1.45_0.1"], out["L50_ref"], f4, pol, Lint_t=50.0, p=0.5)
        d50 = np.max(np.abs([t50[o] for o in obs]), axis=0)
        t200 = tissue_db(dense["-1.45_0.1"], dense["ref"], f4, pol, Lint_t=50.0, p=0.5)
        d200 = np.max(np.abs([t200[o] for o in obs]), axis=0)
        print(f"(c) {pol}, slice at internode 50 um, packing 0.5: solved at 50 um {np.round(d50, 4)} dB vs scaled from the "
              f"200 um solve {np.round(d200, 4)} dB (ratio {np.round(d50 / d200, 3)}); eta 50 um "
              f"{np.round(out['L50_-1.45_0.1'][pol]['eta'], 1)} vs 200 um {np.round(dense['-1.45_0.1'][pol]['eta'], 1)}", flush=True)


def epsrobust():
    """Review round 6: the +-20 % material test (Sim 33) covered node power only. Here the passive-eps
    slice bound at its peak (3 THz, E-par) with eps' and eps'' of axon, myelin, water (water includes
    the host around the node) each scaled by 0.8 / 1.2; eps'_node scanned so the peak can move.
    dB uses the unscaled tissue normalisation (V_cell, I0, n_t), so only the fibre response changes."""
    f1 = [3.0e12]
    path = OUT / "sim36_epsrobust.json"
    out = json.loads(path.read_text()) if path.exists() else {}
    base = dict(axon=axoplasm(), myelin=myelin(), water=WATER_KS1996)
    ers = (-2.0, -1.7, -1.45, -1.2, -1.0)
    for dom in ("axon", "myelin", "water"):
        for part in ("re", "im"):
            for s in (0.8, 1.2):
                ov = {dom: Scaled(base[dom], **{part: s})}
                for er in (None,) + ers:
                    k = f"{dom}_{part}{s}_{er}"
                    if k in out:
                        continue
                    o = dict(ov)
                    if er is not None:
                        o["node"] = ConstantEps(complex(er, -0.1))
                    elif dom == "water":
                        o["node"] = ov["water"]
                    out[k] = sim30.run("s36", tag=f"s36_epsr_{k}", freqs=f1, Lint=200.0, eps_override=o,
                                       pols=("par",), **MESH)
                    path.write_text(json.dumps(out))
    report_epsrobust(out)


def report_epsrobust(out=None):
    out = out or json.loads((OUT / "sim36_epsrobust.json").read_text())
    f1 = np.array([3.0e12])
    ers = np.array([-2.0, -1.7, -1.45, -1.2, -1.0])
    rows = []
    for dom in ("axon", "myelin", "water"):
        for part in ("re", "im"):
            for s in (0.8, 1.2):
                ref = out[f"{dom}_{part}{s}_None"]
                db = np.array([max(abs(tissue_db(out[f"{dom}_{part}{s}_{er}"], ref, f1, "par")[o][0])
                                   for o in ("P_fibre", "flux_in", "dipole")) for er in ers])
                j = int(np.argmax(db))
                rows.append((f"{dom} eps{part} x{s}", db[j], ers[j]))
                print(f"{dom} eps{part} x{s}: max |dB| {db[j]:.3f} at eps' {ers[j]} (all: {np.round(db, 3)}); "
                      f"x6 (internode 200 um, packing 0.6) {6 * db[j]:.2f} dB", flush=True)
    worst = max(rows, key=lambda r: r[1])
    print(f"worst: {worst[0]} {worst[1]:.3f} dB (x6 {6 * worst[1]:.2f} dB); unscaled peak 0.123 dB (x6 0.74)", flush=True)


def report_epsdense(out=None):
    out = out or json.loads((OUT / "sim36_epsdense.json").read_text())
    f4 = np.array([0.6, 1.0, 2.0, 3.0]) * 1e12
    ref = out["ref"]
    keys = sorted((k for k in out if not k.startswith("ref") and k.count("_") == 1), key=lambda k: float(k.split("_")[0]))
    er = np.array([float(k.split("_")[0]) for k in keys])
    for pol in ("par", "perp"):
        db = np.array([[max(abs(v) for v in vals) for vals in zip(*(tissue_db(out[k], ref, f4, pol)[o] for o in ("P_fibre", "flux_in", "dipole")))]
                       for k in keys])
        eta = np.array([out[k][pol]["eta"] for k in keys])
        for i, fi in enumerate(f4):
            j = int(np.argmax(db[:, i])); jj = np.clip([j - 1, j, j + 1], 0, len(er) - 1)
            c = np.polyfit(er[jj], db[jj, i], 2) if len(set(jj)) == 3 else None
            peak = (c[2] - c[1] ** 2 / (4 * c[0])) if c is not None and c[0] < 0 else db[j, i]
            print(f"{pol} {fi/1e12:.1f} THz: max |dB| {db[j, i]:.3f} at eps' {er[j]} (eta {eta[j, i]:.1f}); "
                  f"parabola peak {peak:.3f} dB; x6 (internode 200 um, packing 0.6) {6 * peak:.2f} dB", flush=True)


def report_dense(p=0.6, Lint_t=200.0):
    """Tissue dB of the dense sweeps (worst case Lint 200 um, p 0.6)."""
    out = json.loads((OUT / "sim36_dense.json").read_text())
    f16 = np.linspace(0.1, 3.0, 16) * 1e12
    ref = out["s36_dense_a7.5_0"]
    for k in ("s36_dense_a7.5_30000", "s36_dense_a7.5_1e+06"):
        for pol in ("par", "perp"):
            t = tissue_db(out[k], ref, f16, pol, Lint_t=Lint_t, p=p, a=7.5, b=7.5 / 0.7)
            m = np.max(np.abs(np.vstack([t["P_fibre"], t["flux_in"], t["dipole"]])), axis=0)
            print(f"{k} {pol}: max |dB| over f {m.max():.3f} at {f16[np.argmax(m)]/1e12:.2f} THz", flush=True)


JOINT = dict(axon=Scaled(axoplasm(), re=1.2), myelin=Scaled(myelin(), re=1.2, im=0.8))


def round10():
    """Review round 10: the +-20 % materials test varied one material at a time; three changes raise the
    passive-eps peak (axon eps' x1.2, myelin eps' x1.2, myelin eps'' x0.8). Combine them and rescan
    eps'_node at 0.6 and 3 THz (E-par), reporting loss and phase."""
    f2 = [0.6e12, 3.0e12]
    path = OUT / "sim36_round10.json"
    out = json.loads(path.read_text()) if path.exists() else {}
    for er in (None, -10.0, -7.0, -5.0, -4.0, -3.5, -3.0, -2.7, -2.4, -2.2, -2.0, -1.7, -1.45, -1.3, -1.2, -1.0):
        k = f"joint_{er}"
        if k not in out:
            o = dict(JOINT)
            if er is not None:
                o["node"] = ConstantEps(complex(er, -0.1))
            out[k] = sim30.run("s36", tag=f"s36_r10_{k}", freqs=f2, Lint=200.0, eps_override=o, pols=("par",), **MESH)
            path.write_text(json.dumps(out))
    report_round10(out)


def report_round10(out=None):
    """Joint-materials peak (0.6 and 3 THz), then the signed change of slice loss and phase at 0.6 THz over
    every node state with eps'' >= 0.1 (conductivity scan, coarse and fine eps grids, dense eps' scans,
    joint-materials scan). The q fit of Sec. 5.6 needs extra delay (positive phase) and extra loss together."""
    out = out or json.loads((OUT / "sim36_round10.json").read_text())
    f2 = np.array([0.6e12, 3.0e12]); f4 = np.array([0.6, 1.0, 2.0, 3.0]) * 1e12
    dense = json.loads((OUT / "sim36_epsdense.json").read_text())
    obs = ("P_fibre", "flux_in", "dipole")
    dk = [k for k in dense if not k.startswith("ref") and k.count("_") == 1]
    jk = sorted((k for k in out if k != "joint_None"), key=lambda k: float(k[6:]))
    er = np.array([float(k[6:]) for k in jk])
    for i, (fi, j) in enumerate(zip(f2, (0, 3))):      # j: index of 0.6 / 3 THz in the dense scan
        nd = [tissue_db(dense[k], dense["ref"], f4, "par") for k in dk]
        nom = max(max(abs(t[o][j]) for o in obs) for t in nd)
        tj = [tissue_db(out[k], out["joint_None"], f2, "par") for k in jk]
        db = np.array([max(abs(t[o][i]) for o in obs) for t in tj])
        m = int(np.argmax(db)); jj = np.clip([m - 1, m, m + 1], 0, len(er) - 1)
        c = np.polyfit(er[jj], db[jj], 2) if len(set(jj)) == 3 else None
        peak = c[2] - c[1] ** 2 / (4 * c[0]) if c is not None and c[0] < 0 else db[m]
        print(f"{fi/1e12:.1f} THz E-par, joint worst materials (axon eps' x1.2, myelin eps' x1.2, eps'' x0.8), per 100 um at "
              f"internode 1000 um, packing 0.5: |dB| {dict(zip(er.tolist(), np.round(db, 4).tolist()))}; parabola peak "
              f"{peak:.4f} (nominal materials {nom:.4f}: x{peak/nom:.2f})", flush=True)
    # signed loss and phase at 0.6 THz over every state with eps'' >= 0.1
    src = {"sigma": ("sim36_sigma.json", "0", [0.1, 0.6, 1.0, 2.0, 3.0]), "epsnode": ("sim36_epsnode.json", "ref", [0.6, 2.0]),
           "epsfine": ("sim36_epsfine.json", "ref", [0.6, 2.0]), "epsdense": ("sim36_epsdense.json", "ref", [0.6, 1.0, 2.0, 3.0]),
           "joint": ("sim36_round10.json", "joint_None", [0.6, 3.0])}
    rows = []
    for name, (fn, rk, fl) in src.items():
        d = json.loads((OUT / fn).read_text()); fa = np.array(fl) * 1e12; ii = fl.index(0.6)
        for k, r in d.items():
            if k == rk or k.startswith(("ref", "L")) or "_0.01" in k:
                continue
            ref = d.get(f"ref_{k.split('_')[-1]}") if k.endswith(("_h", "_q")) else d[rk]
            if ref is None:
                continue
            t = tissue_db(r, ref, fa, "par")
            rows.append((f"{name}:{k}", t["dipole"][ii], max(t[o][ii] for o in obs), t["phase_rad"][ii], name == "joint"))
    mat = lambda r: 1.0 if r[4] else 1.45            # joint scan already holds the combined materials
    ext = 6 * 4 * 4 / 3 * 2.2                         # internode 50 um, packing 0.8, 220 um optical path (Sim 45)
    pos = max(rows, key=lambda r: r[3] * mat(r)); neg = min(rows, key=lambda r: r[3] * mat(r))
    los = max(rows, key=lambda r: r[2] * mat(r))
    print(f"0.6 THz, {len(rows)} node states (eps'' >= 0.1), first-order dipole values per 100 um at internode 1000 um, packing 0.5 "
          f"(x1.45 for non-joint scans):", flush=True)
    for lab, r in (("largest extra delay (+phase)", pos), ("largest advance (-phase)", neg), ("largest loss", los)):
        print(f"  {lab}: {r[0]} phase {r[3] * mat(r):+.2e} rad, loss {r[2] * mat(r):+.4f} dB; every extreme (x{ext:.0f}): "
              f"{ext * r[3] * mat(r):+.3f} rad, {ext * r[2] * mat(r):+.2f} dB", flush=True)
    (OUT / "sim36_signed06.json").write_text(json.dumps(
        {"note": "0.6 THz, E-par, every node state with eps'' >= 0.1; first-order change of a slice with every extreme "
                 "(internode 50 um, packing 0.8, 220 um path) and combined materials", "factor": ext,
         "states": [{"state": r[0], "loss_db": ext * mat(r) * r[2], "delay_rad": ext * mat(r) * r[3]} for r in rows]}, indent=1))
    both = [r for r in rows if ext * mat(r) * r[3] >= 0.19 and ext * mat(r) * r[2] >= 1.1]
    print(f"  states giving both >= 0.19 rad extra delay and >= 1.1 dB extra loss with every extreme: {len(both)}", flush=True)

if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    {"sigma": sigma, "report_sigma": report_sigma, "dense": dense, "aniso": aniso, "epsnode": epsnode, "axon": axon, "axon2": axon2, "axon_threshold": axon_threshold, "report_axon": report_axon, "epsfine": epsfine, "epsdense": epsdense, "epsdense_mesh": epsdense_mesh, "epsrobust": epsrobust, "report_epsrobust": report_epsrobust, "report_epsdense": report_epsdense, "report_dense": report_dense, "epsdense_hscan": epsdense_hscan, "report_epsdense_h": report_epsdense_h, "round9": round9, "report_round9": report_round9, "round10": round10, "report_round10": report_round10}[sys.argv[1]]()
