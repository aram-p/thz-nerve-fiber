"""Sim 28 - Phase 0 benchmarks for the new pipeline.

B1  water slab in air (FEM) vs transfer-matrix: R, T, A and R+T+A = 1.
B2  layered infinite cylinder, 2D cross-section (FEM), TM and TE, vs exact series:
    absorption per layer + flux balance through r = R.
B3  axisymmetric harmonic model without a node, m = 0..3, both polarisations,
    vs the series per azimuthal order; null test (fibre = water -> no scattering).

Usage: uv run python -u experiments/bench/sim28_benchmarks.py [B1] [B2] [B3]
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

from thznerve.fem import axisym, cyl2d, slab
from thznerve.fem.common import client, complex_literal, debye_expr
from thznerve.physics import mie_cyl, tmm
from thznerve.physics.materials import C0, ETA0, WATER_KS1996, refractive_index
from thznerve.runlog import RunLog

OUT = Path(__file__).resolve().parents[2] / "results" / "sim28"
WATER = WATER_KS1996
EPS_MY = 2.5 - 0.3j   # placeholder myelin for benchmarks only
EW = debye_expr(WATER)
EM = complex_literal(EPS_MY)


def b1():
    f = np.linspace(0.1, 3.0, 30) * 1e12
    d = 100.0
    with RunLog("sim28-B1", "water slab 100um vs TMM, 30 freqs") as log:
        t0 = time.monotonic()
        model, o = slab.build_and_solve(EW, d, f, h_um=3.0)
        arr = o["arr"]
        print("EvalPoint shape", arr.shape, flush=True)
        # arr layout: [expr][freq*points] or [freq][...]: reshape robustly
        # EvalPoint.getReal(): rows = expr-major over points (expr0@p0, expr0@p1, ...), cols = freqs
        a = arr.reshape(4, 2, len(f))
        res = []
        for i, fi in enumerate(f):
            k0 = 2 * np.pi * fi / C0
            esc_l = complex(a[0, 0, i], a[1, 0, i])
            et_r = complex(a[2, 1, i], a[3, 1, i])
            r = esc_l * np.exp(-1j * k0 * o["xl"] * 1e-6)
            t = et_r * np.exp(1j * k0 * (o["xr"] - d) * 1e-6)
            A = o["P"][i] / (0.5 / ETA0 * o["height"] * 1e-6)
            R_a, T_a, A_a, r_a, t_a = tmm.slab_rta(fi, 1.0, [(complex(WATER(fi)), d * 1e-6)], 1.0)
            res.append([fi / 1e12, abs(r) ** 2, abs(t) ** 2, A, R_a, T_a, A_a])
        res = np.array(res)
        bal = res[:, 1] + res[:, 2] + res[:, 3] - 1
        errT = np.max(np.abs(res[:, 2] - res[:, 5]) / res[:, 5])
        errR = np.max(np.abs(res[:, 1] - res[:, 4]) / res[:, 4])
        log.dof = o["dof"]
        log.outcome = f"max|R+T+A-1|={np.max(np.abs(bal)):.2e} relerr T={errT:.2e} R={errR:.2e}"
        print("B1", log.outcome, f"{time.monotonic()-t0:.0f}s", flush=True)
        np.savetxt(OUT / "B1_slab.csv", res, delimiter=",",
                   header="f_THz,R_fem,T_fem,A_fem,R_tmm,T_tmm,A_tmm")
        client().remove(model)
    return dict(balance=float(np.max(np.abs(bal))), errT=float(errT), errR=float(errR))


def b2():
    f = np.array([0.1, 0.3, 0.6, 1.0, 1.5, 2.0, 2.5, 3.0]) * 1e12
    radii = [4.0, 5.5]
    out = {}
    for pol in ("TM", "TE"):
        with RunLog("sim28-B2", f"2D layered cylinder {pol} vs series, 8 freqs") as log:
            spec = cyl2d.Cyl2DSpec(radii_um=radii, eps_exprs=[EW, EM, EW], pol=pol,
                                   R_um=150, pml_um=150, h_water_um=4.0, h_layer_frac=0.2,
                                   freqs_hz=list(f))
            model, o = cyl2d.build_and_solve(spec, name=f"b2{pol}")
            rows = []
            for i, fi in enumerate(f):
                ew = complex(WATER(fi))
                s = mie_cyl.solve(fi, np.array(radii) * 1e-6, [ew, EPS_MY, ew], pol)
                P_an = s.absorbed_per_layer()   # W/m for |E0| = 1
                P_fe = o["P_layers"][:2, i]
                bal = (o["flux_in"][i] - o["P_layers"][:, i].sum()) / o["flux_in"][i]
                rows.append([fi / 1e12, *P_fe, *P_an, bal])
            rows = np.array(rows)
            err = np.abs(rows[:, 1:3] - rows[:, 3:5]) / np.abs(rows[:, 3:5])
            log.dof = o["dof"]
            log.outcome = f"max rel err axon={err[:,0].max():.2e} myelin={err[:,1].max():.2e} balance={np.abs(rows[:,5]).max():.2e}"
            print("B2", pol, log.outcome, flush=True)
            np.savetxt(OUT / f"B2_{pol}.csv", rows, delimiter=",",
                       header="f_THz,P_axon_fem,P_my_fem,P_axon_an,P_my_an,flux_balance")
            out[pol] = dict(err_axon=float(err[:, 0].max()), err_my=float(err[:, 1].max()),
                            balance=float(np.abs(rows[:, 5]).max()))
            client().remove(model)
    return out


def b3():
    f = np.array([0.1, 0.6, 1.5, 3.0]) * 1e12
    a, b = 4.0, 5.5
    out = {}
    for null in (False, True):
        eps = dict(axon=EW, myelin=EW if null else EM, node=EW if null else EM, water=EW)
        spec = axisym.AxiSpec(a_um=a, b_um=b, Ln_um=1.0, Lint_um=60.0, R_um=120, pml_um=80,
                              eps=eps, h_water_um=4.0, freqs_hz=list(f))
        for pol in ("par", "perp"):
            with RunLog("sim28-B3", f"axisym no-node {'NULL ' if null else ''}{pol} m=0..3 vs series") as log:
                model = axisym.build(spec, name=f"b3{pol}{int(null)}")
                L = (spec.Lint_um + spec.Ln_um) * 1e-6
                rows = []
                for m in range(0, 4):
                    o = axisym.solve_mode(model, m, pol)
                    for i, fi in enumerate(f):
                        ew = complex(WATER(fi))
                        s = mie_cyl.solve(fi, np.array([a, b]) * 1e-6, [ew, ew if null else EPS_MY, ew],
                                          "TM" if pol == "par" else "TE", m_max=8)
                        po = s.absorbed_per_layer(per_order=True)
                        im = list(s.ms).index(m)
                        fe_ax = o["P_axon"][i] / L
                        fe_sh = (o["P_myelin"][i] + o["P_node"][i]) / L
                        tot = o["P_axon"][i] + o["P_myelin"][i] + o["P_node"][i] + o["P_win"][i]
                        bal = (o["flux_in"][i] - tot) / o["flux_in"][i]
                        rows.append([m, fi / 1e12, fe_ax, fe_sh, po[0, im], po[1, im], bal])
                        print(f"B3 null={null} {pol} m={m} f={fi/1e12:.1f} axon fem/an={fe_ax/po[0,im]:.5f} "
                              f"shell fem/an={fe_sh/po[1,im]:.5f} bal={bal:.2e}", flush=True)
                rows = np.array(rows)
                big = rows[:, 4] > 1e-3 * rows[:, 4].max()
                e_ax = np.max(np.abs(rows[big, 2] / rows[big, 4] - 1))
                bigs = rows[:, 5] > 1e-3 * rows[:, 5].max()
                e_sh = np.max(np.abs(rows[bigs, 3] / rows[bigs, 5] - 1))
                log.dof = o["dof"]
                log.outcome = f"axon err={e_ax:.2e} shell err={e_sh:.2e} balance={np.abs(rows[:,6]).max():.2e}"
                print("B3", null, pol, log.outcome, flush=True)
                np.savetxt(OUT / f"B3_{pol}_null{int(null)}.csv", rows, delimiter=",",
                           header="m,f_THz,P_axon_fem,P_shell_fem,P_axon_an,P_shell_an,balance")
                out[f"{pol}_null{int(null)}"] = dict(e_ax=float(e_ax), e_sh=float(e_sh),
                                                     balance=float(np.abs(rows[:, 6]).max()))
                client().remove(model)
    return out


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    which = sys.argv[1:] or ["B1", "B2", "B3"]
    summary = {}
    for w in which:
        summary[w] = {"B1": b1, "B2": b2, "B3": b3}[w]()
        (OUT / f"summary_{w}.json").write_text(json.dumps(summary[w], indent=1))
    print(json.dumps(summary, indent=1))
