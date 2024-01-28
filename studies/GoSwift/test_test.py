import json
import time
import pickle
#import numpy as np
from stl import mesh
from off_body import *
from deform_tri_module import *
start_time = time.time()

## Input Parameters
input_file = "studies/GoSwift/input_files/test.json"



input = json.loads(open(input_file).read())
MACH = input["flow"]["freestream_mach_number"]
N_points = 1000
r_over_l = 3
ref_length = 20 # 154 ft for N+2,  20 ft for SH
altitude = 50000
PROP_R = r_over_l*ref_length #*3.28084  ##### add if mesh is in meters

num_azimuth = 0 ### update in sBOOM function call


## Initialize storage list for FFD Box, Nearfield Sig, Ground Sig, and Loudness
ffd_box = []
x_loc = []
nearfield_sig = []
ground_sig = []
ground_x = []
peound_p = []
loudness = []


## Import Undeformed Test Geometry
baseline_stl = "studies/GoSwift/meshes/test_sw.stl"
baseline_tri = "studies/GoSwift/meshes/test_sw.tri"
#stl_deformed = "studies/GoSwift/meshes/test_sw.stl"


## Create CSV file for off body points
body_mesh = mesh.Mesh.from_file(baseline_stl)
minx, miny, minz, maxx, maxy, maxz, y_pos, z_pos = find_mins(body_mesh)
#off_body(N_points,MACH,r_over_l) 
angles = off_body_sheet(N_points, MACH, r_over_l, num_azimuth)

length_body = maxx - minx

## Convert baseline .stl file to .tri file for FFD
stl_to_tri(baseline_stl, baseline_tri)

# Initialize FFD box info lists
length = []
origin = []
num_points = []
delta_z = []
delta_index = []

## Define initial FFD data for loop     (pre defined locations? or algorithm?)
ffd_lengths0 = (1,1,2)
ffd_origin0 = (0,0,0)
ffd_num_points = (3,3,3) # Constant?
ffd_delta_z0 = (0) 
ffd_delta_index = (1,1,1)  # Constant 

lengths, origins, bumps = (4,int(length_body),4)
#lengths, origins, bumps = (2,int(length_body - ffd_lengths[0]),2)
"'####### Find a way to not check Zero deformation except for the first time #######'"
# for i in range (lengths of mesh - ffd_lengths[0]) ## moves along entire mesh keeping entire ffd box on mesh

for i in range(lengths): # lengths 5
    ffd_lengths = (ffd_lengths0[0] + i, ffd_lengths0[1], ffd_lengths0[2])
    #for j in range(origins): # origins 5
    for j in range(origins): #(int(length_body - ffd_lengths[0])): ## moves along entire mesh keeping entire ffd box on mesh
        ffd_origin = (ffd_origin0[0] + j, ffd_origin0[1]- (ffd_lengths[1]/2), ffd_origin0[2]- (ffd_lengths[2]))
        for k in range(bumps): # delta z 4
            ffd_delta_z = (ffd_delta_z0 - (k/2))


            ## Deform baseline mesh with new FFD box parameters
            vert_coords, tri_verts, comp_num = load_tri_mesh_points(baseline_tri)
            deformed_mesh_points = deform_mesh_points(ffd_lengths, ffd_origin, ffd_num_points, ffd_delta_z, ffd_delta_index, vert_coords,0,True)
            write_tri_file("studies/GoSwift/meshes/test_deformed.tri",deformed_mesh_points,tri_verts,comp_num)

            ## Run MachLine with new deformed geometry
            run_machline('studies/GoSwift/input_files/test.json')

            ## Post processing of bump location and pressure data
            ffd_box.append([ffd_lengths, ffd_origin, ffd_num_points, ffd_delta_z, ffd_delta_index])
            xg,p = pressures(243.61, .00036392, 968.08, 1548.928, MACH ) #### run at 50000 ft      , 1452.12 @ 1.5,       1548.928 @ 1.6
            x_loc.append(xg)

            



            nearfield_sig.append(p)

            #for n in range(len(points)/num_azimuth - 1)   ## loop through lines in pressure file....????
            #   open sheet file 
            #    for row in sheet
            #      find and define one signature for angle
            #    -----skip header
            #   define azimuth angle associated with line
            #   asign azimuth and run sBOOM
            ### figure out how to accurately save each azimuth angle signature data
            ### (header = [azimuth , [signature]] )

            ## Run sBOOM
            g_sig, noise_level = boom(MACH, r_over_l, ref_length, altitude) # add azimuth here....
            ground_sig.append(g_sig)
            loudness.append(noise_level)

    print("Iteration: ", (i*j*k) + 1 , "/", lengths*origins*bumps)

print()
print("Writing Test Case Data to Files....")

### write all data to seperate pickle files for later data handling
with open('studies/Goswift/results/temp_ffd_box.pkl', 'wb') as pickle_file:
    pickle.dump(ffd_box, pickle_file)

with open('studies/Goswift/results/temp_x_loc.pkl', 'wb') as pickle_file:
    pickle.dump(x_loc, pickle_file)

with open('studies/Goswift/results/temp_nearfield.pkl', 'wb') as pickle_file:
    pickle.dump(nearfield_sig, pickle_file)

with open('studies/Goswift/results/temp_ground.pkl', 'wb') as pickle_file:
    pickle.dump(ground_sig, pickle_file)

with open('studies/Goswift/results/temp_loudness.pkl', 'wb') as pickle_file:
    pickle.dump(loudness, pickle_file)

print()
print("------------------------------------------")
print("Test Case Executed Successfully!")
print("Number of Bump Configurations: ", lengths*origins*bumps)
end_time = time.time()
execution_time = end_time - start_time
print("Execution Time: ", execution_time / 60 , " min")
print("------------------------------------------")
print()