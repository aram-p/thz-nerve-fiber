"""Find which feature needs an isotropic permittivity when myelin is a (rr, pp, zz) tensor."""
import importlib.util
from pathlib import Path
from thznerve.fem import axisym
from thznerve.physics.materials import myelin_lamellar, WATER_KS1996, axoplasm
ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("sim30", ROOT / "experiments/prod/sim30_node_spectrum.py")
sim30 = importlib.util.module_from_spec(spec); spec.loader.exec_module(sim30)
eps = dict(axon=axoplasm().comsol(), myelin=myelin_lamellar().comsol(), node=WATER_KS1996.comsol(), water=WATER_KS1996.comsol())
s = axisym.AxiSpec(a_um=3.5, b_um=5.0, Ln_um=1.0, Lint_um=50.0, R_um=100.0, pml_um=60.0, eps=eps, outer="pml",
                   h_water_um=8.0, h_fibre_frac=0.5, h_far_frac=1.0, freqs_hz=[1e12])
model = axisym.build(s, name="probeaniso")
j = model.java; comp = j.component("comp1"); ph = comp.physics("ewfd")
for t in ph.feature().tags():
    f = ph.feature(t); print(t, f.getType())
    for p in map(str, f.properties()):
        if any(k in p.lower() for k in ("model", "epsilon", "n_", "wavelength")):
            try: print("   ", p, f.getString(p))
            except Exception: pass
pml = comp.coordSystem("pml1")
for p in map(str, pml.properties()):
    try: print(" pml", p, pml.getString(p))
    except Exception: pass
def tryit(label):
    try:
        axisym.solve_mode(model, 0, "par"); print("OK", label)
        return True
    except Exception as e:
        print("FAIL", label, str(e).splitlines()[-1][:150])
tryit("as built")
ph.feature("wee1").set("DisplacementFieldModel", "RelativePermittivity")
tryit("RelativePermittivity")
