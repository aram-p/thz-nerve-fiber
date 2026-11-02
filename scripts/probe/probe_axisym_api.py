"""Probe COMSOL API names for 2D axisymmetric EWFD (mode number, background field, PML)."""
import time
import mph

mph.option("session", "stand-alone")
c = mph.start(cores=4)
m = c.create("probe")
j = m.java
j.component().create("comp1", True)
g = j.component("comp1").geom().create("geom1", 2)
g.axisymmetric(True)
g.lengthUnit("um")
r = g.create("r1", "Rectangle"); r.set("size", ["50", "20"]); g.run()
ph = j.component("comp1").physics().create("ewfd", "ElectromagneticWavesFrequencyDomain", "geom1")
for pn in ph.prop().tags() if hasattr(ph.prop(), "tags") else []:
    print("prop", pn)
for name in ["components", "BackgroundField", "ShapeProperty", "AzimuthalModeNumber", "outofplanewavenumber", "mFloquet"]:
    try:
        p = ph.prop(name)
        print(name, "->", [str(x) for x in p.properties()])
        for k in p.properties():
            try:
                print("   ", k, "=", [str(v) for v in p.getStringArray(k)])
            except Exception:
                try:
                    print("   ", k, "=", p.getString(k))
                except Exception as e:
                    print("   ", k, "?", e)
    except Exception as e:
        print(name, "ERR", str(e)[:100])
print("physics props:", [str(x) for x in ph.properties()])
for k in ph.properties():
    try:
        print("  ", k, "=", ph.getString(k))
    except Exception:
        pass
wee = ph.feature("wee1")
print("wee1 props:", [str(x) for x in wee.properties()])
