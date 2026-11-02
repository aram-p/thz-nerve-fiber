"""Mesh check for Sim 36 'axon' at high axoplasm conductivity: refine the whole axon, not only the node."""
import importlib.util, json
from pathlib import Path
import numpy as np
from thznerve.physics.materials import PlusSigma, axoplasm
ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("sim30", ROOT / "experiments/prod/sim30_node_spectrum.py")
sim30 = importlib.util.module_from_spec(spec); spec.loader.exec_module(sim30)
F5 = list(np.array([0.1, 0.6, 1.0, 2.0, 3.0]) * 1e12)
base = json.loads((sim30.OUT / "sim36_axon.json").read_text())
for ds in (1e3, 1e4, 1e5):
    r = sim30.run("s36", tag=f"s36_axon_{ds:g}_far", freqs=F5, Lint=200.0, h_water=5.0, h_frac=0.25, h_far_frac=0.125,
                  eps_override=dict(axon=PlusSigma(axoplasm(), ds)))
    for pol in ("par", "perp"):
        b = base[f"{ds:g}"][pol]
        rel_p = np.array(r[pol]["P_fibre"]) / np.array(b["P_fibre"]) - 1
        rel_d = (np.array(r[pol]["p_re"]) + 1j * np.array(r[pol]["p_im"])) / (np.array(b["p_re"]) + 1j * np.array(b["p_im"])) - 1
        print(f"axon +{ds:g} {pol}: fine/base - 1: P_fibre {np.round(rel_p, 4)}, |dipole| {np.round(np.abs(rel_d), 4)}", flush=True)
