"""What frequency did the legacy pipeline actually solve at?"""
import numpy as np, mph
from thznerve.model.geometry import GeometryParams, build_geometry
from thznerve.model.materials import MaterialParams, apply_materials
from thznerve.model.mesh import MeshParams, build_mesh
from thznerve.model.study import setup_physics, setup_study, _evaluate_at_points
mph.option("session", "stand-alone"); c = mph.start(cores=4); m = c.create("p")
g = GeometryParams(axon_radius_um=5, myelin_radius_um=7, node_length_um=40, internode_length_um=100, external_half_width_um=20)
build_geometry(m, g); apply_materials(m, MaterialParams(node_sigma_S_per_m=0.0), freq_hz=0.6e12)
setup_physics(m, g); build_mesh(m, MeshParams()); setup_study(m, 0.6e12)
st = m.java.study("std1").feature("freq")
print("legacy plist", st.getString("plist"), "punit", st.getString("punit"))
m.java.study("std1").run()
print("solved freq (Hz):", _evaluate_at_points(m, "freq", np.array([[0, 0, 50.0]])))
