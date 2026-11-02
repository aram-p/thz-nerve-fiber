"""Sim 31 - Mesh convergence and boundary robustness of the Sim 30 node observables.

Levels refine every element size by ~sqrt(2) and 2. Pass when the node absorption, eta and
the fibre absorption change by < 3 % between the two finest levels. Also varies the outer
boundary (PML vs scattering BC) and the water radius R.

Usage: uv run python -u experiments/prod/sim31_convergence.py [mesh|short|bc]
"""

from __future__ import annotations

import json
import sys
import importlib.util
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("sim30", ROOT / "experiments/prod/sim30_node_spectrum.py")
sim30 = importlib.util.module_from_spec(spec); spec.loader.exec_module(sim30)

COMMON = dict(fmin=0.1, fmax=3.0, nf=4)
LEVELS = {
    "L1": dict(h_water=5.0, h_frac=0.25, h_far_frac=0.5),
    "L2": dict(h_water=3.5, h_frac=0.125, h_far_frac=0.35),
    "L3": dict(h_water=2.5, h_frac=0.0625, h_far_frac=0.25),
}
BCS = {
    "R150": dict(R=150.0),
    "pml100": dict(pml=100.0),
    "sbc": dict(outer="sbc"),
}


def compare(results, ref):
    out = {}
    for k, r in results.items():
        if k == ref:
            continue
        d = {}
        for pol in ("par", "perp"):
            for q in ("P_node", "eta", "P_fibre"):
                a, b = np.array(r[pol][q]), np.array(results[ref][pol][q])
                d[f"{pol}_{q}"] = float(np.max(np.abs(a / b - 1)))
        out[f"{k} vs {ref}"] = d
        print(k, "vs", ref, {kk: f"{v:.2e}" for kk, v in d.items()}, flush=True)
    return out


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "mesh"
    sim30.OUT.mkdir(parents=True, exist_ok=True)
    res = {}
    if which == "short":
        # finest levels do not fit in 8 GB at Lint = 1000 um; use Lint = 200 um (Sim 32 shows
        # node observables do not depend on Lint) plus a node-only refinement L4n.
        lv = LEVELS | {"L4n": dict(h_water=2.5, h_frac=0.03125, h_far_frac=0.25)}
        for k, v in lv.items():
            res[k] = sim30.run("conv", tag=f"convS_{k}", Lint=200.0, **COMMON, **v)
        cmp = {}
        for k in ("L1", "L2", "L4n"):
            cmp |= compare({k: res[k], "L3": res["L3"]}, "L3")
    elif which == "mesh":
        for k, lv in LEVELS.items():
            res[k] = sim30.run("conv", tag=f"conv_{k}", **COMMON, **lv)
        cmp = compare({k: res[k] for k in ("L1", "L3")}, "L3") | compare({k: res[k] for k in ("L2", "L3")}, "L3")
    else:
        res["ref"] = sim30.run("bc", tag="bc_ref", **COMMON, **LEVELS["L1"])
        for k, b in BCS.items():
            res[k] = sim30.run("bc", tag=f"bc_{k}", **COMMON, **LEVELS["L1"], **b)
        cmp = compare(res, "ref")
    (sim30.OUT / f"sim31_{which}.json").write_text(json.dumps(cmp, indent=1))
