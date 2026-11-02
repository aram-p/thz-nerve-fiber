"""Sim 32 - Scaling of node absorption with node length, myelin thickness, axon radius and
internode length, full wave vs the quasi-static law (mechanism test).

Prediction (quasi-static, no free parameter): P_node = 1/2 w eps0 eps''_w |E0|^2 V_node F,
with F from thznerve.physics.quasistatic for the same geometry and materials. F depends on
Ln/t, a/b and the permittivity ratios only; it does not depend on the internode length
once Lint >> b (no inter-node coupling, water damps any standing wave).
A resonance would show up as eta/F_qs departing from 1 at some frequency or size.

Usage: uv run python -u experiments/prod/sim32_scaling.py [Ln|t|a|Lint|all]
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

from thznerve.physics.materials import WATER_KS1996, axoplasm, myelin
from thznerve.physics.quasistatic import solve_node_qs

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("sim30", ROOT / "experiments/prod/sim30_node_spectrum.py")
sim30 = importlib.util.module_from_spec(spec); spec.loader.exec_module(sim30)

F_LIST = [0.1, 0.6, 1.0, 2.0, 3.0]
BASE = dict(a=3.5, b=5.0, Ln=1.0, Lint=1000.0)
SWEEPS = {
    "Ln": [dict(Ln=x) for x in (0.5, 2.0, 5.0)],
    "t": [dict(b=3.5 + t) for t in (0.5, 3.0)],
    "a": [dict(a=1.0, b=1.0 / 0.7, Lint=200.0), dict(a=7.5, b=7.5 / 0.7, Lint=200.0)],   # Lint-independent (Sim 31), keeps a=1 under the RAM budget
    "Lint": [dict(Lint=x) for x in (200.0, 500.0, 1500.0)],
}


def qs(p, f_hz):
    out = {"par": [], "perp": []}
    for fi in f_hz:
        e = [complex(axoplasm()(fi)), complex(myelin()(fi)), complex(WATER_KS1996(fi)), complex(WATER_KS1996(fi))]
        out["par"].append(solve_node_qs(p["a"], p["b"], p["Ln"], *e).F_node)
        out["perp"].append(solve_node_qs(p["a"], p["b"], p["Ln"], *e, field="transverse").F_node)
    return out


def one(tag, over):
    p = BASE | over
    f = np.array(F_LIST) * 1e12
    res = sim30.run("scal", tag=f"scal_{tag}", fmin=0, fmax=0, nf=0, freqs=list(f),
                    h_water=5.0, h_frac=0.25, h_far_frac=0.5, **p)          # level L1 (Sim 31)
    q = qs(p, f)
    row = dict(params=p)
    for pol in ("par", "perp"):
        eta = np.array(res[pol]["eta"])
        row[pol] = dict(eta=eta.tolist(), F_qs=q[pol], ratio=(eta / np.array(q[pol])).tolist(),
                        P_node=res[pol]["P_node"], share=res[pol]["share"])
        print(f"{tag} {pol} eta/F_qs = {np.round(eta / np.array(q[pol]), 3)}", flush=True)
    return row


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    keys = list(SWEEPS) if which == "all" else [which]
    out = {}
    for k in keys:
        for i, over in enumerate(SWEEPS[k]):
            tag = f"{k}{i}"
            out[tag] = one(tag, over)
            (sim30.OUT / f"sim32_{which}.json").write_text(json.dumps(out, indent=0))
