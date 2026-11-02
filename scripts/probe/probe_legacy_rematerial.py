"""Do material edits after setup reach the solver in the legacy pipeline?"""
import numpy as np
import mph
from thznerve.model.geometry import GeometryParams, build_geometry
from thznerve.model.materials import MaterialParams, apply_materials, water_eps_expr
from thznerve.model.mesh import MeshParams, build_mesh
from thznerve.model.study import setup_physics, setup_study, _evaluate_at_points

mph.option("session", "stand-alone")
c = mph.start(cores=4)
m = c.create("p")
g = GeometryParams(axon_radius_um=5, myelin_radius_um=7, node_length_um=40, internode_length_um=100, external_half_width_um=20)
build_geometry(m, g)
apply_materials(m, MaterialParams(node_sigma_S_per_m=0.0), freq_hz=0.6e12)
setup_physics(m, g); build_mesh(m, MeshParams()); setup_study(m, 0.6e12)
pts = np.array([[0, 0, 50.0], [6, 0, 50.0], [0, 0, 120.0]])
def show(tag):
    m.java.study("std1").run()
    print(tag, "eps_myelin", _evaluate_at_points(m, "ewfd.epsilonrxx", pts[1:2], complex_result=True),
          "|E|", np.round(_evaluate_at_points(m, "ewfd.normE", pts), 4))
show("1 initial")
comp = m.java.component("comp1")
w = water_eps_expr(0.6e12)
for t in ("mat_myelin_proximal", "mat_myelin_distal", "mat_node"):
    comp.material(t).propertyGroup("def").set("relpermittivity", [w, "0", "0", "0", w, "0", "0", "0", w])
show("2 propertyGroup.set water")
apply_materials(m, MaterialParams(node_sigma_S_per_m=0.0), freq_hz=0.6e12)
show("3 apply_materials (remove+create)")
for t in ("mat_myelin_proximal", "mat_myelin_distal", "mat_node"):
    comp.material(t).propertyGroup("def").set("relpermittivity", [w, "0", "0", "0", w, "0", "0", "0", w])
show("4 set water again")
comp.material("mat_node").propertyGroup("def").set("electricconductivity", ["1e10", "0", "0", "0", "1e10", "0", "0", "0", "1e10"])
show("5 node sigma=1e10 via electricconductivity")
