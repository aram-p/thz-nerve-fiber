"""Why does the legacy 3D model ignore the fibre? Inspect selections and eps in the solved model."""
import numpy as np
import mph
from thznerve.model.geometry import GeometryParams, build_geometry
from thznerve.model.materials import MaterialParams, apply_materials
from thznerve.model.mesh import MeshParams, build_mesh
from thznerve.model.study import setup_physics, setup_study, _evaluate_at_points

mph.option("session", "stand-alone")
c = mph.start(cores=4)
m = c.create("legacy_probe")
g = GeometryParams(axon_radius_um=5, myelin_radius_um=7, node_length_um=40, internode_length_um=100, external_half_width_um=20)
nd = build_geometry(m, g)
print("domains", nd)
apply_materials(m, MaterialParams(node_sigma_S_per_m=0.0), freq_hz=0.6e12)
comp = m.java.component("comp1")
for t in comp.material().tags():
    mat = comp.material(str(t))
    print(str(t), "entities", list(mat.selection().entities()), "eps", list(mat.propertyGroup("def").getStringArray("relpermittivity"))[:1])
for lab in ["axon", "myelin_proximal", "node", "myelin_distal", "external"]:
    print(lab, list(comp.selection(f"geom1_sel_{lab}").entities()))
setup_physics(m, g)
build_mesh(m, MeshParams())
setup_study(m, 0.6e12)
m.java.study("std1").run()
pts = np.array([[0, 0, 50.0], [6, 0, 50.0], [6, 0, 120.0], [15, 0, 120.0]])
for e in ["ewfd.epsilonrxx", "ewfd.normE", "abs(ewfd.relEx)", "abs(ewfd.Ex)"]:
    print(e, _evaluate_at_points(m, e, pts, complex_result=True))
ph = comp.physics("ewfd")
print("Eb", list(ph.prop("BackgroundField").getStringArray("Eb")), ph.prop("BackgroundField").getString("SolveFor"))
for f in ph.feature().tags():
    print("feature", str(f), ph.feature(str(f)).getType(), list(ph.feature(str(f)).selection().entities()))
