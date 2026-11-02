"""2D-axisymmetric myelinated fibre with a periodic chain of nodes (normal incidence).

Geometry (r, z), one period Lambda = L_internode + L_node, node centred at z = 0:
    axon    0 < r < a,  |z| < Lambda/2
    myelin  a < r < b,  Ln/2 < |z| < Lambda/2
    node    a < r < b,  |z| < Ln/2           (node material: ECF-like)
    water   b < r < R
    PML     R < r < R + t_pml                 (or a scattering BC at r = R)
Periodic (continuity) at z = +-Lambda/2: an infinite fibre with nodes every
Lambda, illuminated at normal incidence (k perpendicular to the fibre).

The incident plane wave (unit amplitude at the axis, in the lossy host) is split
into azimuthal harmonics exp(-j m phi). For each m the scattered-field
background is the exact m-th harmonic of the plane wave in water:
    E parallel to fibre ("par"):  Eb = (0, 0, J_m(k r))
    E perpendicular ("perp"):     Eb = ((J_{m-1}+J_{m+1})/2, -j (J_{m-1}-J_{m+1})/2, 0)
Harmonics are orthogonal, so absorbed powers add: P = sum_m P_m, with
P_{-m} = P_m by mirror symmetry.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from thznerve.fem.common import client, integrate, n_dof, set_eps

DOMAINS = ("axon", "myelin", "node", "water", "pml")


@dataclass
class AxiSpec:
    a_um: float = 4.0            # axon radius
    b_um: float = 5.5            # myelin outer radius
    Ln_um: float = 1.0           # node length
    Lint_um: float = 400.0       # internode length (period = Lint + Ln)
    R_um: float = 120.0          # water radius
    pml_um: float = 60.0
    outer: str = "pml"           # "pml" or "sbc" (cylindrical-wave scattering BC at r = R)
    eps: dict = field(default_factory=dict)   # COMSOL expressions: axon, myelin, node, water
    h_water_um: float = 5.0
    h_fibre_frac: float = 0.25   # element size in fibre <= frac * min(t, Ln)
    order: int = 2               # shape-function order
    freqs_hz: list = field(default_factory=lambda: [1e12])


def _box(geom, tag, dim, x0, x1, y0, y1):
    s = geom.create(tag, "BoxSelection")
    s.set("entitydim", str(dim))
    s.set("xmin", str(x0)); s.set("xmax", str(x1))
    s.set("ymin", str(y0)); s.set("ymax", str(y1))
    s.set("condition", "inside")
    return s


def _rect(geom, tag, r0, r1, z0, z1):
    rc = geom.create(tag, "Rectangle")
    rc.set("pos", [str(r0), str(z0)])
    rc.set("size", [str(r1 - r0), str(z1 - z0)])


def build(spec: AxiSpec, name: str = "axi"):
    c = client()
    model = c.create(name)
    j = model.java
    comp = j.component().create("comp1", True)
    geom = comp.geom().create("geom1", 2)
    geom.axisymmetric(True)
    geom.lengthUnit("um")
    a, b, Ln, R = spec.a_um, spec.b_um, spec.Ln_um, spec.R_um
    P = spec.Lint_um + Ln
    h = P / 2
    rmax = R + (spec.pml_um if spec.outer == "pml" else 0.0)
    _rect(geom, "r_ax", 0, a, -h, h)
    _rect(geom, "r_my1", a, b, -h, -Ln / 2)
    _rect(geom, "r_nd", a, b, -Ln / 2, Ln / 2)
    _rect(geom, "r_my2", a, b, Ln / 2, h)
    _rect(geom, "r_w", b, R, -h, h)
    if spec.outer == "pml":
        _rect(geom, "r_pml", R, rmax, -h, h)
    geom.run("fin")
    e = 1e-3 * min(Ln, b - a)
    _box(geom, "s_axon", 2, -e, a + e, -h - e, h + e)
    _box(geom, "s_node", 2, a - e, b + e, -Ln / 2 - e, Ln / 2 + e)
    _box(geom, "s_my1", 2, a - e, b + e, -h - e, -Ln / 2 + e)
    _box(geom, "s_my2", 2, a - e, b + e, Ln / 2 - e, h + e)
    u = geom.create("s_myelin", "UnionSelection")
    u.set("entitydim", "2"); u.set("input", ["s_my1", "s_my2"])
    _box(geom, "s_water", 2, b - e, R + e, -h - e, h + e)
    if spec.outer == "pml":
        _box(geom, "s_pml", 2, R - e, rmax + e, -h - e, h + e)
    _box(geom, "s_rR", 1, R - e, R + e, -h - e, h + e)       # flux surface r = R
    _box(geom, "s_zlo", 1, -e, rmax + e, -h - e, -h + e)
    _box(geom, "s_zhi", 1, -e, rmax + e, h - e, h + e)
    up = geom.create("s_per", "UnionSelection")
    up.set("entitydim", "1"); up.set("input", ["s_zlo", "s_zhi"])
    geom.run()

    j.param().set("m", "0")
    var = comp.variable().create("var1")
    for d in ("axon", "myelin", "node", "water"):
        var.set(f"eps_{d}", spec.eps[d])
    var.set("kw", "ewfd.k0*sqrt(eps_water)")
    for d in ("axon", "myelin", "node", "water"):
        set_eps(comp, f"mat_{d}", f"geom1_s_{d}", f"eps_{d}")
    if spec.outer == "pml":
        set_eps(comp, "mat_pml", "geom1_s_pml", "eps_water")
        pml = comp.coordSystem().create("pml1", "PML")
        pml.selection().named("geom1_s_pml")

    ph = comp.physics().create("ewfd", "ElectromagneticWavesFrequencyDomain", "geom1")
    ph.prop("outofplanewavenumber").set("mFloquet", "m")
    ph.prop("BackgroundField").set("SolveFor", "scatteredField")
    ph.prop("ShapeProperty").set("order_electricfield", str(spec.order))
    pc = ph.create("pc1", "PeriodicCondition", 1)
    pc.selection().named("geom1_s_per")
    if spec.outer == "sbc":
        sb = ph.create("sbc1", "Scattering", 1)
        sb.selection().named("geom1_s_rR")
        try:
            sb.set("WaveType", "CylindricalWave")
        except Exception:
            pass

    mesh = comp.mesh().create("mesh1")
    ft = mesh.create("ftri", "FreeTri")
    sz = mesh.feature("size")
    sz.set("custom", "on"); sz.set("hmax", str(spec.h_water_um)); sz.set("hmin", str(1e-3))
    sz.set("hgrad", "1.2")
    hf = spec.h_fibre_frac * min(b - a, Ln)
    for tag, sel, hh in (("sf", "geom1_s_node", hf), ("sm", "geom1_s_myelin", spec.h_fibre_frac * (b - a)),
                         ("sa", "geom1_s_axon", max(spec.h_fibre_frac * a, hf))):
        fs = ft.create(tag, "Size")
        fs.selection().geom("geom1", 2)
        fs.selection().named(sel)
        fs.set("custom", "on"); fs.set("hmaxactive", True); fs.set("hmax", str(hh))
    mesh.run()

    std = j.study().create("std1")
    st = std.create("freq", "Frequency")
    st.set("plist", " ".join(f"{f:.10g}" for f in spec.freqs_hz))
    return model


def set_background(model, pol: str):
    ph = model.java.component("comp1").physics("ewfd")
    if pol == "par":
        eb = ["0", "0", "besselj(m,kw*r)"]
    else:
        eb = ["(besselj(m-1,kw*r)+besselj(m+1,kw*r))/2",
              "-i*(besselj(m-1,kw*r)-besselj(m+1,kw*r))/2", "0"]
    ph.prop("BackgroundField").set("Eb", eb)


QEXPR = ("0.5*2*pi*freq*epsilon0_const*(-imag(eps_{d}))"
         "*(abs(ewfd.Er)^2+abs(ewfd.Ephi)^2+abs(ewfd.Ez)^2)")


def solve_mode(model, m: int, pol: str):
    """Solve one azimuthal order over the frequency list. Returns dict of arrays (W per period)."""
    model.java.param().set("m", str(m))
    set_background(model, pol)
    model.java.study("std1").run()
    out = {}
    for d in ("axon", "myelin", "node", "water"):
        r = integrate(model, "IntSurface", f"geom1_s_{d}", [QEXPR.format(d=d)], tag="iq")
        out[f"P_{d}"] = r[0]
    out["flux_in"] = integrate(model, "IntLine", "geom1_s_rR", ["-ewfd.Poavr"], tag="ifl")[0]
    out["dof"] = n_dof(model)
    return out
