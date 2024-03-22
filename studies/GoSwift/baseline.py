import trimesh
from stl import mesh
import os

os.remove("studies/GoSwift/meshes/test_sw.stl")

input_stl = "studies/GoSwift/meshes/unit.stl"
output_file = "studies/GoSwift/meshes/test_sw.stl"

scale_factor = 1 / 304.8
body_mesh = mesh.Mesh.from_file(input_stl)
body_mesh.vectors *= scale_factor

body_mesh.save(output_file)

print("Done")
