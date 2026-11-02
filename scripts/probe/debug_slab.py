import numpy as np
from thznerve.fem import slab
from thznerve.fem.common import debye_expr
from thznerve.physics.materials import WATER_KS1996
model, o = slab.build_and_solve(debye_expr(WATER_KS1996), 100.0, [1e12], h_um=3.0)
j = model.java
print("eps in slab expr:", j.component("comp1").variable("var1").get("eps_s"))
xs = np.linspace(-290, 390, 35)
ds = j.result().dataset().create("cl", "CutPoint2D")
ds.set("pointx", [str(x) for x in xs]); ds.set("pointy", ["1"] * len(xs))
ep = j.result().numerical().create("e2", "EvalPoint"); ep.set("data", "cl")
ep.set("expr", ["abs(ewfd.Ez)", "abs(ewfd.relEz)", "real(ewfd.epsilonrzz)", "imag(ewfd.epsilonrzz)", "abs(ewfd.Ex)+abs(ewfd.Ey)"])
a = np.array(ep.getReal()).reshape(5, len(xs))
for i, x in enumerate(xs):
    print(f"{x:7.1f} |E|={a[0,i]:.4f} |Esc|={a[1,i]:.4f} eps={a[2,i]:.3f}{a[3,i]:+.3f}j inplane={a[4,i]:.2e}")
print("dof", o["dof"])
ep2 = j.result().numerical().create("e3", "EvalPoint"); ep2.set("data", "cl")
ep2.set("expr", ["freq", "ewfd.k0", "ewfd.freq"])
print("freq,k0:", np.array(ep2.getReal())[:, 0])
st = j.study("std1").feature("freq")
print("plist", st.getString("plist"), "punit", st.getString("punit") if "punit" in [str(x) for x in st.properties()] else "n/a")
print([str(x) for x in st.properties()][:40])
