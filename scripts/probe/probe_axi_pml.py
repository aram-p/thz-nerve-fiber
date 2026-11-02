from thznerve.fem import axisym
from thznerve.fem.common import client
spec = axisym.AxiSpec(Lint_um=20, R_um=40, pml_um=20, eps=dict(axon="1.5", myelin="2", node="2", water="4-1i"), h_water_um=4)
model = axisym.build(spec, "pp")
pml = model.java.component("comp1").coordSystem("pml1")
for k in pml.properties():
    try:
        print(k, pml.getString(k))
    except Exception:
        print(k, "<array>")
for k in ("PMLType", "ScalingType", "typeOfScaling"):
    try:
        print(k, "allowed", list(pml.getAllowedPropertyValues(k)))
    except Exception as e:
        pass
