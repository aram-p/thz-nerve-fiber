"""Calibrate DOF per mesh element for the axisym model (no solve)."""
import mph
from thznerve.fem import axisym
from thznerve.physics.materials import WATER_KS1996, axoplasm, myelin
EPS = {k: m.comsol() for k, m in dict(axon=axoplasm(), myelin=myelin(), node=WATER_KS1996, water=WATER_KS1996).items()}
mph.option("session", "stand-alone"); client = mph.start(cores=4)
for kw in (dict(), dict(b=4.0), dict(Lint=200.0), dict(a=1.0, b=1/0.7), dict(a=1.0, b=1/0.7, Lint=200.0), dict(a=7.5, b=7.5/0.7, Lint=200.0), dict(Lint=1500.0)):
    p = dict(a=3.5, b=5.0, Ln=1.0, Lint=1000.0) | kw
    spec = axisym.AxiSpec(a_um=p["a"], b_um=p["b"], Ln_um=p["Ln"], Lint_um=p["Lint"], R_um=100.0, eps=EPS, freqs_hz=[1e12])
    m = axisym.build(spec, name="cal")
    n = m.java.component("comp1").mesh("mesh1").getNumElem()
    print(kw, n, flush=True)
    client.remove(m)
