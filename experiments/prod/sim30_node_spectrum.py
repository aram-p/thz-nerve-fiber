"""Sim 30 - Does the node of Ranvier absorb THz resonantly? Realistic single fibre, 0.1-3 THz.

Hypothesis (H30): the node shows no resonance in 0.1-3 THz. Its absorbed power follows the
quasi-static law P_node = 1/2 w eps0 eps''_node |E0|^2 V_node F with F ~ 0.5 (E || fibre),
smooth in frequency.

Model: 2D-axisymmetric fibre, periodic node chain (period Lint + Ln), azimuthal orders
m = 0..m_max summed, cited materials (water Kindt & Schmuttenmaer 1996; myelin = Looyenga
40 % water + lipid; axoplasm = Looyenga 85 % water + protein; node gap = water).
Control: same fibre with the node gap filled by myelin ("no node").

Observables per frequency and polarisation:
  P_node (W per node for |E0| = 1 V/m), sigma_node = P_node / I0 (m^2),
  eta = P_node / (1/2 w eps0 eps''_w |E0|^2 V_node) = <|E|^2>_node / |E0|^2,
  share = P_node / P_fibre (per period), dP = P_fibre(node) - P_fibre(no node).

Usage: uv run python -u experiments/prod/sim30_node_spectrum.py <case> [key=value ...]
  cases: base | conv (mesh refinement) | nonode (control)
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

from thznerve.fem import axisym
from thznerve.fem.common import client
from thznerve.physics.materials import C0, EPS0, ETA0, WATER_KS1996, axoplasm, myelin, refractive_index
from thznerve.runlog import RunLog

OUT = Path(__file__).resolve().parents[2] / "results" / "sim30"
WATER = WATER_KS1996
MY = myelin(WATER)
AX = axoplasm(WATER)


def run(case: str, a=3.5, b=5.0, Ln=1.0, Lint=1000.0, R=100.0, pml=60.0, h_water=5.0, h_frac=0.25,
        m_max=3, fmin=0.1, fmax=3.0, nf=30, node="water", tag=None):
    f = np.linspace(fmin, fmax, int(nf)) * 1e12
    node_m = {"water": WATER, "myelin": MY}[node]
    eps = dict(axon=AX.comsol(), myelin=MY.comsol(), node=node_m.comsol(), water=WATER.comsol())
    spec = axisym.AxiSpec(a_um=a, b_um=b, Ln_um=Ln, Lint_um=Lint, R_um=R, pml_um=pml, eps=eps,
                          h_water_um=h_water, h_fibre_frac=h_frac, freqs_hz=list(f))
    tag = tag or case
    res = {}
    for pol in ("par", "perp"):
        with RunLog("sim30", f"{tag} {pol} a={a} b={b} Ln={Ln} Lint={Lint} h={h_water}/{h_frac} m<={m_max} nf={nf}") as log:
            t0 = time.monotonic()
            model = axisym.build(spec, name=f"s30{tag}{pol}")
            o = axisym.solve_pol(model, pol, m_max)
            client().remove(model)
            w = 2 * np.pi * f
            ew = WATER(f)
            nw = refractive_index(ew).real
            I0 = 0.5 * nw / ETA0                                   # W/m^2 for |E0| = 1 V/m
            Vn = np.pi * (b ** 2 - a ** 2) * Ln * 1e-18            # m^3
            q_bulk = 0.5 * w * EPS0 * (-np.imag(node_m(f)))         # W/m^3 at |E| = 1
            P_fib = o["P_axon"] + o["P_myelin"] + o["P_node"]
            r = dict(f_THz=f / 1e12, P_node=o["P_node"], P_axon=o["P_axon"], P_myelin=o["P_myelin"],
                     P_fibre=P_fib, sigma_node=o["P_node"] / I0, eta=o["P_node"] / (q_bulk * Vn),
                     share=o["P_node"] / P_fib, balance=o["balance"], last_m=o["last_order_share"])
            log.dof = o["dof"]
            log.outcome = (f"eta {r['eta'].min():.3f}-{r['eta'].max():.3f} share {r['share'].min():.2e}-"
                           f"{r['share'].max():.2e} bal<={r['balance'].max():.1e} lastm<={r['last_m'].max():.1e}")
            print(f"sim30 {tag} {pol}: {log.outcome} ({time.monotonic() - t0:.0f}s)", flush=True)
            cols = list(r)
            np.savetxt(OUT / f"{tag}_{pol}.csv", np.column_stack([r[c] for c in cols]), delimiter=",",
                       header=",".join(cols))
            res[pol] = {k: v.tolist() for k, v in r.items()}
    (OUT / f"{tag}.json").write_text(json.dumps(dict(params=dict(a=a, b=b, Ln=Ln, Lint=Lint, R=R, pml=pml,
                                                                 h_water=h_water, h_frac=h_frac, m_max=m_max,
                                                                 node=node), res=res), indent=0))
    return res


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    case = sys.argv[1] if len(sys.argv) > 1 else "base"
    kw = {}
    for a in sys.argv[2:]:
        k, v = a.split("=")
        kw[k] = v if k in ("node", "tag") else float(v)
    if "m_max" in kw:
        kw["m_max"] = int(kw["m_max"])
    if case == "nonode":
        kw.setdefault("node", "myelin")
    run(case, **kw)
