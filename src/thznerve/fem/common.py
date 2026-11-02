"""Shared COMSOL helpers (mph stand-alone, exp(+jwt))."""

from __future__ import annotations

from typing import Any, Sequence

import numpy as np

from thznerve.physics.materials import DoubleDebye, PlusSigma

_CLIENT = None


def client(cores: int = 4):
    """One stand-alone COMSOL session per Python process."""
    global _CLIENT
    if _CLIENT is None:
        import mph

        mph.option("session", "stand-alone")
        _CLIENT = mph.start(cores=cores)
    return _CLIENT


def debye_expr(m: DoubleDebye, extra_sigma: float = 0.0) -> str:
    """COMSOL analytic expression of a double-Debye model in `freq` (+ optional sigma)."""
    return (PlusSigma(m, extra_sigma) if extra_sigma else m).comsol()


def complex_literal(z: complex) -> str:
    return f"({z.real:.10g}+({z.imag:.10g})*i)"


def ball_sel(geom, tag: str, x: float, y: float, r: float = 1e-3, dim: int = 2, z: float | None = None):
    s = geom.create(tag, "BallSelection")
    s.set("entitydim", str(dim))
    s.set("posx", str(x)); s.set("posy", str(y))
    if z is not None:
        s.set("posz", str(z))
    s.set("r", str(r))
    s.set("condition", "intersects")
    return s


def set_eps(comp, tag: str, sel: str, expr: str):
    mats = comp.material()
    if tag in [str(t) for t in mats.tags()]:
        mats.remove(tag)
    m = mats.create(tag, "Common")
    m.propertyGroup("def").set("relpermittivity", [expr, "0", "0", "0", expr, "0", "0", "0", expr])
    m.propertyGroup("def").set("relpermeability", ["1", "0", "0", "0", "1", "0", "0", "0", "1"])
    m.propertyGroup("def").set("electricconductivity", ["0", "0", "0", "0", "0", "0", "0", "0", "0"])
    m.selection().named(sel)
    return m


def integrate(model, kind: str, sel: str | None, exprs: Sequence[str], tag: str = "intx",
              dataset: str = "dset1", revolve: bool = False) -> np.ndarray:
    """Integrate expressions over a named selection for every solution. Returns (n_expr, n_sol).

    revolve=True (2D axisymmetric): weight by 2*pi*r, i.e. IntSurface -> volume integral and
    IntLine -> surface integral. COMSOL leaves this off by default, so an axisymmetric
    integral without it is per radian per metre and off by 1/(2*pi*r_mean).
    """
    num = model.java.result().numerical()
    if tag in [str(t) for t in num.tags()]:
        num.remove(tag)
    f = num.create(tag, kind)
    f.set("data", dataset)
    if revolve:
        f.set({"IntSurface": "intvolume", "IntLine": "intsurface"}[kind], "on")
    if sel is not None:
        f.selection().named(sel)
    f.set("expr", list(exprs))
    res = np.array(f.getReal(), dtype=float)
    # getReal: [n_expr][n_sol] for a frequency list (one inner solution per freq)
    return res


def n_dof(model) -> int:
    return int(len(model.java.sol("sol1").getU()))
