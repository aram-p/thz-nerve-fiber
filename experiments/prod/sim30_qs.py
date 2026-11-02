"""Quasi-static node factor F(f) for the Sim 30 geometry, both polarisations, same materials."""
import time
from pathlib import Path

import numpy as np

from thznerve.physics.materials import WATER_KS1996, axoplasm, myelin
from thznerve.physics.quasistatic import node_factor

OUT = Path(__file__).resolve().parents[2] / "results" / "sim30"
OUT.mkdir(parents=True, exist_ok=True)
f = np.linspace(0.1, 3.0, 30) * 1e12
rows = []
t0 = time.monotonic()
for fi in f:
    e = [complex(axoplasm()(fi)), complex(myelin()(fi)), complex(WATER_KS1996(fi)), complex(WATER_KS1996(fi))]
    Fp = node_factor(3.5, 5.0, 1.0, *e)
    Ft = node_factor(3.5, 5.0, 1.0, *e, field="transverse")
    rows.append([fi / 1e12, Fp, Ft])
    print(f"{fi/1e12:.1f} THz F_par={Fp:.4f} F_perp={Ft:.4f}", flush=True)
np.savetxt(OUT / "qs.csv", np.array(rows), delimiter=",", header="f_THz,F_par,F_perp")
print(f"{time.monotonic()-t0:.0f}s")
