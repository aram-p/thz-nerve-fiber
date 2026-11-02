"""Fit the two geometric constants (c, p) of the closed node law to the extrapolated Laplace
solver. The Laplace problem has no length scale, so F depends only on x = Ln/t, g = a/b and the
permittivities; the axon radius is fixed at 3.5 um.

Grid: x 0.2-8 (10 values), g 0.6/0.7/0.8, f 0.1-3 THz (6), both polarisations (360 solves, ~10 min).
Cached in results/sim30/nodelaw_grid.json; delete it to recompute.
Usage: uv run python -u experiments/prod/fit_nodelaw.py
"""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from thznerve.physics import nodelaw
from thznerve.physics.materials import WATER_KS1996, axoplasm, myelin
from thznerve.physics.quasistatic import node_factor

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "results" / "sim30" / "nodelaw_grid.json"
X = [0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0]
G = [0.6, 0.7, 0.8]
F = [0.1e12, 0.3e12, 0.6e12, 1e12, 2e12, 3e12]
A = 3.5

if CACHE.exists():
    rows = json.loads(CACHE.read_text())
else:
    rows = []
    for g in G:
        b = A / g
        for x in X:
            Ln = x * (b - A)
            for f in F:
                e = [complex(axoplasm()(f)), complex(myelin()(f)), complex(WATER_KS1996(f))]
                for pol, fld in (("par", "axial"), ("perp", "transverse")):
                    v = node_factor(A, b, Ln, e[0], e[1], e[2], e[2], field=fld)
                    rows.append(dict(g=g, x=x, f=f, pol=pol, F=v))
            print(f"g={g} x={x} done", flush=True)
    CACHE.write_text(json.dumps(rows))


def law(cp, r, pol):
    f = r["f"]
    e = [complex(axoplasm()(f)), complex(myelin()(f)), complex(WATER_KS1996(f))]
    b = A / r["g"]
    F0 = abs(e[1] / e[2]) ** 2 if pol == "par" else nodelaw.shell_factor(*e, r["g"])
    xp = r["x"] ** cp[1]
    return F0 + (1 - F0) * xp / (xp + cp[0])


for pol in ("par", "perp"):
    rr = [r for r in rows if r["pol"] == pol]
    fit = least_squares(lambda cp: [law(cp, r, pol) / r["F"] - 1 for r in rr], x0=[1.0, 1.0],
                        bounds=([0.01, 0.2], [20, 3]))
    err = np.abs(fit.fun)
    print(f"{pol}: c = {fit.x[0]:.3f}, p = {fit.x[1]:.3f}; max error {err.max():.3f}, rms {np.sqrt((err ** 2).mean()):.3f}")
