"""Check the full-wave fibre dipole against the infinite coated-cylinder quasi-static value."""
import importlib.util
from pathlib import Path
import numpy as np
from thznerve.physics.effmedium import coated_cylinder_eps
from thznerve.physics.materials import WATER_KS1996, axoplasm, myelin

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("sim30", ROOT / "experiments/prod/sim30_node_spectrum.py")
sim30 = importlib.util.module_from_spec(spec); spec.loader.exec_module(sim30)
f = np.array([0.1e12, 1e12])
a, b, Ln, Lint = 3.5, 5.0, 1.0, 200.0
r = sim30.run("probe", tag="probe_dipole", Lint=Lint, freqs=list(f))
for i, fi in enumerate(f):
    ea, em, ew = complex(axoplasm()(fi)), complex(myelin()(fi)), complex(WATER_KS1996(fi))
    L = Lint * 1e-6                      # fibre length without node (node adds 1 um of water)
    pz = ((ea - ew) * np.pi * a ** 2 + (em - ew) * np.pi * (b ** 2 - a ** 2)) * 1e-12 * L
    ec = coated_cylinder_eps(ea, em, a / b)
    px = np.pi * (b * 1e-6) ** 2 * 2 * ew * (ec - ew) / (ec + ew) * L
    for pol, an in (("par", pz), ("perp", px)):
        fe = r[pol]["p_re"][i] + 1j * r[pol]["p_im"][i]
        print(f"{fi/1e12} THz {pol}: FEM {fe:.4e} analytic {an:.4e} ratio {fe/an:.4f}", flush=True)
