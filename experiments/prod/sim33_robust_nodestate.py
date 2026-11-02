"""Sim 33 - (a) robustness of node absorption to +/-20 % material changes;
(b) which node states could change THz absorption enough to explain a voltage effect.

(a) Each of axon, myelin, water: eps' and eps'' scaled by 0.8 and 1.2 (one at a time).
    Pass: no new spectral feature (eta/F_qs stays smooth and ~1), observables move
    smoothly by O(20 %).
(b) Node gap = water + extra conductivity d_sigma (0 ... 1e6 S/m). Reports P_node(d_sigma),
    the change of absorption per period, and the transmission change it would give a 100 um
    tissue slab (thznerve.physics.samplebound), versus the 9-13 dB reported features.
    Physiological reference: open Na channels ~3 S/cm^2 over a 5-10 nm membrane is
    ~1e-4 S/m as a volume conductivity (docs/literature.md section 3).

Usage: uv run python -u experiments/prod/sim33_robust_nodestate.py [robust|state]
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

from thznerve.physics.materials import WATER_KS1996, PlusSigma, Scaled, axoplasm, myelin, refractive_index
from thznerve.physics.quasistatic import solve_node_qs
from thznerve.physics.materials import ETA0
from thznerve.physics.samplebound import DB_PER_NEPER_POWER, node_volume_fraction

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("sim30", ROOT / "experiments/prod/sim30_node_spectrum.py")
sim30 = importlib.util.module_from_spec(spec); spec.loader.exec_module(sim30)

F = np.array([0.1, 0.6, 1.0, 2.0, 3.0]) * 1e12
MESH = dict(h_water=5.0, h_frac=0.25, h_far_frac=0.5)     # level L1 (Sim 31)
BASE = dict(axon=axoplasm(), myelin=myelin(), node=WATER_KS1996, water=WATER_KS1996)


def robust():
    out = {}
    ref = sim30.run("rob", tag="rob_ref", freqs=list(F), **MESH)
    out["ref"] = ref
    for dom in ("axon", "myelin", "water"):
        for part in ("re", "im"):
            for s in (0.8, 1.2):
                ov = {dom: Scaled(BASE[dom], **{part: s})}
                if dom == "water":                      # the node gap is extracellular water too
                    ov["node"] = ov["water"]
                tag = f"rob_{dom}_{part}{s}"
                r = sim30.run("rob", tag=tag, freqs=list(F), eps_override=ov, **MESH)
                out[tag] = r
                for pol in ("par", "perp"):
                    rel = np.array(r[pol]["P_node"]) / np.array(ref[pol]["P_node"])
                    print(f"{tag} {pol} P_node/ref = {np.round(rel, 3)}", flush=True)
                (sim30.OUT / "sim33_robust.json").write_text(json.dumps(out))
    return out


def state():
    out = {}
    fn = node_volume_fraction(3.5, 5.0, 1.0, 1000.0)
    for ds in (0.0, 1e-4, 1e2, 1e3, 1e4, 1e5, 1e6):
        node = PlusSigma(WATER_KS1996, ds) if ds else WATER_KS1996
        r = sim30.run("state", tag=f"state_{ds:g}", freqs=list(F), eps_override=dict(node=node), **MESH)
        out[f"{ds:g}"] = r
        (sim30.OUT / "sim33_state.json").write_text(json.dumps(out))
    ref = out["0"]
    nw = refractive_index(WATER_KS1996(F)).real
    summary = {}
    for k, r in out.items():
        for pol in ("par", "perp"):
            dP = np.array(r[pol]["P_fibre"]) - np.array(ref[pol]["P_fibre"])   # W per period
            # extra absorption per tissue volume / incident intensity = extra attenuation (1/m)
            V_cell = np.pi * (5.0 ** 2 - 3.5 ** 2) * 1.0 * 1e-18 / fn          # m^3 per fibre period
            I0 = 0.5 * nw / ETA0
            ddb = DB_PER_NEPER_POWER * dP / (V_cell * I0) * 100e-6
            P_node0 = np.array(ref[pol]["P_node"])
            summary[f"{k}_{pol}"] = dict(P_node_rel=(np.array(r[pol]["P_node"]) / P_node0).tolist(),
                                         dB_100um=ddb.tolist())
            print(f"d_sigma={k} S/m {pol}: P_node x{np.round(np.array(r[pol]['P_node']) / P_node0, 3)} "
                  f"slab dB {np.array2string(ddb, precision=2)}", flush=True)
    (sim30.OUT / "sim33_state_summary.json").write_text(json.dumps(summary, indent=0))
    return out


if __name__ == "__main__":
    sim30.OUT.mkdir(parents=True, exist_ok=True)
    {"robust": robust, "state": state}[sys.argv[1] if len(sys.argv) > 1 else "state"]()
