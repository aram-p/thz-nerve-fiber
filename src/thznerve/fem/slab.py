"""Benchmark 1: water slab in air, normal incidence, vs the transfer-matrix method.

2D, out-of-plane Ez, x is the propagation direction, periodic in y. Scattered
field with the exact background exp(-j k0 x) in air; PMLs at both x ends.
Extracts r, t from the field and A from the volume integral of the absorbed
power density, so R + T + A = 1 is an independent energy check.
"""

from __future__ import annotations

import numpy as np

from thznerve.fem.common import client, integrate, n_dof, set_eps
from thznerve.physics.materials import C0, ETA0


def build_and_solve(eps_expr: str, d_um: float, freqs_hz, air_um: float = 150.0, pml_um: float = 150.0,
                    h_um: float = 4.0, height_um: float = 2.0, name: str = "slab"):
    c = client()
    model = c.create(name)
    j = model.java
    comp = j.component().create("comp1", True)
    g = comp.geom().create("geom1", 2)
    g.lengthUnit("um")
    xs = [-air_um - pml_um, -air_um, 0.0, d_um, d_um + air_um, d_um + air_um + pml_um]
    for i in range(5):
        r = g.create(f"r{i}", "Rectangle")
        r.set("pos", [str(xs[i]), "0"]); r.set("size", [str(xs[i + 1] - xs[i]), str(height_um)])
    g.run("fin")
    tags = ["pmlL", "airL", "slab", "airR", "pmlR"]
    for i, t in enumerate(tags):
        s = g.create(t, "BoxSelection")
        s.set("entitydim", "2")
        s.set("xmin", str(xs[i] - 1e-3)); s.set("xmax", str(xs[i + 1] + 1e-3))
        s.set("ymin", "-1e-3"); s.set("ymax", str(height_um + 1e-3)); s.set("condition", "inside")
    for t, y in (("yb", 0.0), ("yt", height_um)):
        s = g.create(t, "BoxSelection")
        s.set("entitydim", "1"); s.set("xmin", str(xs[0] - 1)); s.set("xmax", str(xs[-1] + 1))
        s.set("ymin", str(y - 1e-3)); s.set("ymax", str(y + 1e-3)); s.set("condition", "inside")
    u = g.create("yper", "UnionSelection"); u.set("entitydim", "1"); u.set("input", ["yb", "yt"])
    g.run()
    comp.variable().create("var1").set("eps_s", eps_expr)
    set_eps(comp, "m_slab", "geom1_slab", "eps_s")
    for t in ("pmlL", "airL", "airR", "pmlR"):
        set_eps(comp, f"m_{t}", f"geom1_{t}", "1")
    for t in ("pmlL", "pmlR"):
        p = comp.coordSystem().create(f"c_{t}", "PML")
        p.selection().named(f"geom1_{t}")
    ph = comp.physics().create("ewfd", "ElectromagneticWavesFrequencyDomain", "geom1")
    ph.prop("components").set("components", "outofplane")
    ph.prop("BackgroundField").set("SolveFor", "scatteredField")
    ph.prop("BackgroundField").set("Eb", ["0", "0", "exp(-i*ewfd.k0*x)"])
    pc = ph.create("pc1", "PeriodicCondition", 1)
    pc.selection().named("geom1_yper")
    mesh = comp.mesh().create("mesh1")
    mp = mesh.create("map1", "Map")
    sz = mesh.feature("size"); sz.set("custom", "on"); sz.set("hmax", str(h_um))
    mesh.run()
    std = j.study().create("std1")
    std.create("freq", "Frequency").set("plist", " ".join(f"{f:.10g}" for f in freqs_hz))
    std.run()

    xl, xr = -air_um / 2, d_um + air_um / 2
    y0 = height_um / 2
    # point evaluation of fields at all frequencies
    num = j.result().numerical()
    ds = j.result().dataset().create("cpt", "CutPoint2D")
    ds.set("pointx", [str(xl), str(xr)]); ds.set("pointy", [str(y0), str(y0)])
    ep = num.create("ep", "EvalPoint"); ep.set("data", "cpt")
    ep.set("expr", ["real(ewfd.relEz)", "imag(ewfd.relEz)", "real(ewfd.Ez)", "imag(ewfd.Ez)"])
    arr = np.array(ep.getReal(), dtype=float)  # [n_expr*?][...]
    q = "0.5*2*pi*freq*epsilon0_const*(-imag(eps_s))*abs(ewfd.Ez)^2"
    P = integrate(model, "IntSurface", "geom1_slab", [q], tag="iq")[0]
    return model, dict(arr=arr, P=P, xl=xl, xr=xr, d=d_um, height=height_um, dof=n_dof(model))
