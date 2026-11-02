"""Why does the B2 flux balance read 1? Print selection entities and flux values."""
import numpy as np

from thznerve.fem import cyl2d
from thznerve.fem.common import complex_literal, debye_expr, integrate
from thznerve.physics.materials import WATER_KS1996

EW = debye_expr(WATER_KS1996)
spec = cyl2d.Cyl2DSpec(radii_um=[4.0, 5.5], eps_exprs=[EW, complex_literal(2.5 - 0.3j), EW], pol="TM",
                       R_um=150, pml_um=150, h_water_um=4.0, h_layer_frac=0.2, freqs_hz=[1e12])
model, o = cyl2d.build_and_solve(spec, name="pb2")
comp = model.java.component("comp1")
for t in ("geom1_bR", "geom1_bRin", "geom1_bflux", "geom1_d0", "geom1_d1", "geom1_d2", "geom1_d3"):
    print(t, list(comp.selection(t).entities()))
print("P_layers", o["P_layers"].ravel(), "flux", o["flux_in"])
for e in ["-(ewfd.Poavx*x+ewfd.Poavy*y)/sqrt(x^2+y^2)", "ewfd.Poavx", "1"]:
    print(e, integrate(model, "IntLine", "geom1_bflux", [e], tag="t" + str(abs(hash(e)) % 1000))[0])
