import json
from off_body import *
#from deform_tri_module import *


file = "studies/GoSwift/input_files/test.json"

input = json.loads(open(file).read())

Mach = input["flow"]["freestream_mach_number"]
print(Mach)

N_points = 1000
#Mach = 1.5
#gamma = 1.4
r_over_l = 3
ref_length = 200 ###27.432 # for X-59   ########update for n+2
altitude = 50000
PROP_R = r_over_l*ref_length*3.28084  ##### ??? 3.2


#body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/sears_haack_9800.stl')
#body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/100_30_SH.stl')
#body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/pod(SH).stl')
#body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/pod_vsp_smooth.stl')
#body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/test.stl')
#body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/pod_vsp_bump.stl')
#body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/SH_half.stl')
#body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/n+2.stl')
#body_mesh = mesh.Mesh.from_file('studies/GoSwift/meshes/delta_wing.stl')
#body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/low_boom.stl')
#body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/low_boom2.stl')
#body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/biswis.stl')
#body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/JWB_0AOA.stl')
body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/test_sw.stl')
#------------------------------------------------------------------------------------------------
#body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/pod_smooth.stl')
#body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/pod_bump.stl')
find_mins(body_mesh)
off_body(N_points,Mach,r_over_l)  #N_points, Mach Number
#run_machline('studies/GoSwift/input_files/pod_smooth.json')
#run_machline('studies/GoSwift/input_files/pod_bump.json')
#--------------------------------------------------------------------------------------------------
#run_machline('studies/GoSwift/input_files/SH_input.json')
#run_machline('studies/GoSwift/input_files/SH_bump_input.json')
#run_machline('studies/GoSwift/input_files/n+2_input.json')
#run_machline('studies/GoSwift/input_files/input.json')
#run_machline('studies/GoSwift/input_files/low_boom_input.json')
#run_machline('studies/GoSwift/input_files/jaxa_input.json')
run_machline('studies/GoSwift/input_files/test.json')


pressures(243.61, .00036392, 968.08, 1548.928, Mach ) #### run at 50000 ft      , 1452.12 @ 1.5,       1548.928 @ 1.6
#pressures(295.07, 14170 ) #### run at 1400 m 5219

#undefx = []
#def = 
#plt.plot(xundef,pundef)
#plt.plot(xdef,pdef)
#plt.show()