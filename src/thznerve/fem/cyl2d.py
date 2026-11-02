"""2D cross-section of an infinite layered cylinder in water, normal incidence.

Scattered-field EWFD with the *correct* background plane wave in the (lossy)
host, exp(-j k_w x) with k_w = k0 sqrt(eps_host), and a cylindrical PML.
Benchmarked against the exact series in thznerve.physics.mie_cyl.

pol "TM": E along the cylinder axis (out-of-plane Ez).
pol "TE": E perpendicular to the axis (in-plane, E along y).
Energy balance: inward flux of the total field through the water/PML circle
r = R equals the power absorbed inside it.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from thznerve.fem.common import ball_sel, client, integrate, n_dof, set_eps
from thznerve.physics.materials import C0, ETA0


@dataclass
class Cyl2DSpec:
    radii_um: list[float]          # interface radii, ascending (core ... last shell)
    eps_exprs: list[str]           # COMSOL eps expressions, core ... host (len = len(radii)+1)
    pol: str = "TM"
    R_um: float = 150.0            # water radius (PML inner radius)
    pml_um: float = 100.0
    h_water_um: float = 6.0        # max element size in water / PML
    h_layer_frac: float = 0.25     # max element size in each layer as fraction of its thickness
    freqs_hz: list[float] = field(default_factory=lambda: [1e12])


def build_and_solve(spec: Cyl2DSpec, name: str = "cyl2d"):
    c = client()
    model = c.create(name)
    j = model.java
    comp = j.component().create("comp1", True)
    geom = comp.geom().create("geom1", 2)
    geom.lengthUnit("um")
    radii = list(spec.radii_um) + [spec.R_um, spec.R_um + spec.pml_um]
    for i, r in enumerate(radii):
        cc = geom.create(f"c{i}", "Circle")
        cc.set("r", str(r))
    geom.run("fin")
    # domain selections: one ball per annulus
    mids = [0.0] + [0.5 * (radii[i] + radii[i + 1]) for i in range(len(radii) - 1)]
    for i, rm in enumerate(mids):
        ball_sel(geom, f"d{i}", rm, 0.0, r=1e-3, dim=2)
    # boundary circle r = R (water/PML interface) for the flux balance
    ibR = len(spec.radii_um)
    s = geom.create("bR", "BallSelection")
    s.set("entitydim", "1"); s.set("posx", "0"); s.set("posy", "0")
    s.set("r", str(spec.R_um * 1.001)); s.set("condition", "inside")
    s2 = geom.create("bRin", "BallSelection")
    s2.set("entitydim", "1"); s2.set("posx", "0"); s2.set("posy", "0")
    s2.set("r", str(spec.R_um * 0.999)); s2.set("condition", "inside")
    d = geom.create("bflux", "DifferenceSelection")
    d.set("entitydim", "1"); d.set("add", ["bR"]); d.set("subtract", ["bRin"])
    geom.run()

    n_lay = len(spec.eps_exprs)
    host = spec.eps_exprs[-1]
    j.param().set("dummy", "0")
    comp.variable().create("var1")
    v = comp.variable("var1")
    v.set("epsh", host)
    v.set("kw", "ewfd.k0*sqrt(epsh)")
    for i, e in enumerate(spec.eps_exprs):
        v.set(f"eps{i}", e)
    # materials: layers 0..n_lay-1, host also in PML domain
    for i in range(n_lay):
        set_eps(comp, f"mat{i}", f"geom1_d{i}", f"eps{i}")
    set_eps(comp, "matpml", f"geom1_d{n_lay}", "epsh")

    pml = comp.coordSystem().create("pml1", "PML")
    pml.selection().named(f"geom1_d{n_lay}")
    pml.set("ScalingType", "Cylindrical")

    ph = comp.physics().create("ewfd", "ElectromagneticWavesFrequencyDomain", "geom1")
    if spec.pol == "TM":
        ph.prop("components").set("components", "outofplane")
        eb = ["0", "0", "exp(-i*kw*x)"]
    else:
        ph.prop("components").set("components", "inplane")
        eb = ["0", "exp(-i*kw*x)", "0"]
    ph.prop("BackgroundField").set("SolveFor", "scatteredField")
    ph.prop("BackgroundField").set("Eb", eb)

    mesh = comp.mesh().create("mesh1")
    ft = mesh.create("ftri", "FreeTri")
    ft.selection().all()
    sz = mesh.feature("size")
    sz.set("custom", "on"); sz.set("hmax", str(spec.h_water_um)); sz.set("hmin", "0.01")
    rr = [0.0] + list(spec.radii_um)
    for i in range(len(spec.radii_um)):
        t = rr[i + 1] - rr[i]
        fs = ft.create(f"s{i}", "Size")
        fs.selection().geom("geom1", 2)
        fs.selection().named(f"geom1_d{i}")
        fs.set("custom", "on"); fs.set("hmaxactive", True)
        fs.set("hmax", str(max(t * spec.h_layer_frac, 0.02)))
    mesh.run()

    std = j.study().create("std1")
    st = std.create("freq", "Frequency")
    st.set("plist", " ".join(f"{f:.10g}" for f in spec.freqs_hz))
    std.run()

    # absorbed power per unit length in each layer and water annulus (W/m)
    q = "0.5*2*pi*freq*epsilon0_const*(-imag(eps{i}))*(abs(ewfd.Ex)^2+abs(ewfd.Ey)^2+abs(ewfd.Ez)^2)"
    P = np.vstack([integrate(model, "IntSurface", f"geom1_d{i}", [q.format(i=i)])[0] for i in range(n_lay)])
    # inward total-field flux through r = R (W/m)
    flux = integrate(model, "IntLine", "geom1_bflux", ["-(ewfd.Poavx*x+ewfd.Poavy*y)/sqrt(x^2+y^2)"])[0]
    out = dict(P_layers=P, flux_in=flux, dof=n_dof(model), freqs=np.array(spec.freqs_hz))
    return model, out


def intensity_host(eps_host: complex) -> float:
    n = np.sqrt(complex(eps_host))
    return 0.5 * abs(n.real) / ETA0
