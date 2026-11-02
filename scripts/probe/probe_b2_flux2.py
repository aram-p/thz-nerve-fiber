"""B2 flux off by exactly 2: compare ewfd.Poav with a hand-built Poynting vector."""
from thznerve.fem import cyl2d
from thznerve.fem.common import complex_literal, debye_expr, integrate
from thznerve.physics.materials import WATER_KS1996

EW = debye_expr(WATER_KS1996)
for pol in ("TM", "TE"):
    spec = cyl2d.Cyl2DSpec(radii_um=[4.0, 5.5], eps_exprs=[EW, complex_literal(2.5 - 0.3j), EW], pol=pol,
                           R_um=150, pml_um=150, h_water_um=4.0, h_layer_frac=0.2, freqs_hz=[1e12])
    model, o = cyl2d.build_and_solve(spec, name="pb2" + pol)
    print(pol, "sumP", o["P_layers"].sum(), "flux", o["flux_in"])
    ex = {
        "poav": "-(ewfd.Poavx*x+ewfd.Poavy*y)/sqrt(x^2+y^2)",
        "hand": "-0.5*real((ewfd.Ey*conj(ewfd.Hz)-ewfd.Ez*conj(ewfd.Hy))*x+(ewfd.Ez*conj(ewfd.Hx)-ewfd.Ex*conj(ewfd.Hz))*y)/sqrt(x^2+y^2)",
        "relbg": "-(ewfd.Poavx*x+ewfd.Poavy*y)/sqrt(x^2+y^2)",
    }
    for k, e in ex.items():
        print(" ", k, integrate(model, "IntLine", "geom1_bflux", [e], tag="t" + k)[0])
    # inner flux at r=R evaluated slightly inside water only (up side / down side)
    for k, e in (("up", "-(up(ewfd.Poavx)*x+up(ewfd.Poavy)*y)/sqrt(x^2+y^2)"),
                 ("down", "-(down(ewfd.Poavx)*x+down(ewfd.Poavy)*y)/sqrt(x^2+y^2)")):
        print(" ", k, integrate(model, "IntLine", "geom1_bflux", [e], tag="u" + k)[0])
