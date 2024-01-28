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
ref_length = 20 ###27.432 # for X-59   ## figureout what the correct refference length is.....
altitude = 50000
PROP_R = r_over_l*ref_length #*3.28084  ##### ??? 3.2

num_azimuth = 0
'''anything over 1 says roll rate is higher than -50.38'''

baseline_stl = "studies/GoSwift/meshes/test_sw.stl"
baseline_tri = "studies/GoSwift/meshes/test_deformed.tri"

body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/test_sw.stl')

#find_mins(body_mesh)
#off_body(N_points,Mach,r_over_l)  #N_points, Mach Number

body_mesh = mesh.Mesh.from_file(baseline_stl)
minx, miny, minz, maxx, maxy, maxz, y_pos, z_pos = find_mins(body_mesh)
#x0, x, xf,y, L, points = off_body(N_points, Mach, r_over_l)
#print("L: ", L)
#off_body(N_points,Mach,r_over_l)
angles = off_body_sheet(N_points, Mach, r_over_l, num_azimuth)
print(angles)


stl_to_tri(baseline_stl, baseline_tri)
#run_machline('studies/GoSwift/input_files/test.json')
#run_machline('studies/GoSwift/input_files/sheet_input.json')     


xg, p = pressures(243.61, .00036392, 968.08, 1548.928, Mach, angles) #### run at 50000 ft      , 1452.12 @ 1.5,       1548.928 @ 1.6

print("len p: ", len(p))
print("len_xg: ", len(xg))
#pressures(295.07, 14170 ) #### run at 1400 m 5219

g_sig, loudness = boom(Mach, r_over_l, ref_length, altitude, N_points, angles)

#print(len(p))
#plt.plot(xg,p)
#plt.show()

print("Loudness: ", loudness)

for i in range(len(angles)):
    plt.plot(g_sig[i][:,0],g_sig[i][:,1])
plt.title("Ground Signature")
plt.ylabel("dP/P")
plt.xlabel("X (ft)")
plt.show()

#plt.plot(g_sig[1][:,0],g_sig[1][:,1])
#plt.title("Ground Signature")
#plt.ylabel("dP/P")
#plt.xlabel("X (ft)")
#plt.show()
#
#plt.plot(g_sig[2][:,0],g_sig[2][:,1])
#plt.title("Ground Signature")
#plt.ylabel("dP/P")
#plt.xlabel("X (ft)")
#plt.show()