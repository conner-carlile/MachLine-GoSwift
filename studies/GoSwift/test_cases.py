import json
import time
import pickle
#import numpy as np
from stl import mesh
from off_body import *
from deform_tri_module import *
from ambiance import Atmosphere
import trimesh
from scipy.signal import savgol_filter
from exponentalSmoothing import exponential_smoothing
start_time = time.time()

## Input Parameters
input_file = "studies/GoSwift/input_files/test.json"

altitude = 40000 #ft
Mach = 1.6
r_over_l = 1.5 ##0.1 for pod and wedge test ## 3 body lengths for wedge  #1.5 ## 1 F-15 body length for wedge     #2.94 #2.3185 
N_points = 900 # @ 1 F-15 body length
gamma = 1.4
num_azimuth = 0
angle_of_attack = 0 #degrees
solver = "QRUP"
formulation = "dirichlet-morino"

ffd_lengths0 =     (1, 1, 1)
ffd_origin0 =      (0,0,0)
ffd_num_points =   (3,3,3)
ffd_delta_z0 =     (0)
ffd_delta_index =  (1,1,1)

## Atmospheric Parameters
atmosphere = Atmosphere(altitude*0.3048)
p_static = atmosphere.pressure[0] * 0.020885434273039   # Pa to lbf/ft^2 conversion factor
density = atmosphere.density[0] * 0.00194032            # kg/m^3 to slugs/ft^3 conversion factor
speed_of_sound = atmosphere.speed_of_sound[0] * 3.28084 # m/s to ft/s conversion factor
v_inf = Mach*speed_of_sound

print("Mach: ", Mach)
print("R/L: ", r_over_l)
print("Altitude: ", altitude)
print("Static Pressure: ", p_static)
print("Density: ", density)
print("Gamma: ", gamma)
print("Speed of Sound: ", speed_of_sound)
print("V_inf: ", v_inf)
print()

#p_static = 393.13 #lbf/ft^2
#density = .00058728
#gamma = 1.4
#speed_of_sound = 968.08 #ft/s
#v_inf = 1548.928

##AOA Calculations
x = np.cos(np.radians(angle_of_attack))*v_inf
y = 0.0
z = np.sin(np.radians(angle_of_attack))*v_inf

set_json(input_file, [x,y,z], gamma, Mach, formulation, solver)


## Initialize storage list for FFD Box, Nearfield Sig, Ground Sig, and Loudness
ffd_box = []
x_loc = []
nearfield_sig = []
ground_sig = []
ground_x = []
peound_p = []
loudness = []


## Import Undeformed Test Geometry
baseline_stl = "studies/GoSwift/meshes/test.stl"
baseline_tri = "studies/GoSwift/meshes/test.tri"


## Create CSV file for off body points
body_mesh = mesh.Mesh.from_file(baseline_stl)
minx, miny, minz, maxx, maxy, maxz, y_pos, z_pos = find_mins(body_mesh) 
angles = off_body_sheet(N_points, Mach, r_over_l, num_azimuth)  ### angle stuff

length_body = maxx - minx
width_body = maxy - miny

## refrence length for sBOOM
ref_length = length_body
PROP_R = r_over_l*ref_length #*3.28084  ##### add if mesh is in meters


## Convert baseline .stl file to .tri file for FFD
stl_to_tri(baseline_stl, baseline_tri)


## Initialize FFD box info lists
length = []
origin = []
num_points = []
delta_z = []
delta_index = []

## Define initial FFD data for loop    

### pod by itself
#ffd_lengths0 = (60,width_body,width_body)
#ffd_origin0 = (60,0,-z_pos)#(length_body/2 - 60/2,0,-z_pos) #(0,0,0)
#ffd_num_points = (3,3,3) # Constant?
#ffd_delta_z0 = (0) 
#ffd_delta_index = (1,1,1)  # Constant 

### Wedge and pod start at front of body
#ffd_lengths0 =     (1, 1.5, 1)
##ffd_lengths0 =    
#ffd_origin0 =      (16, 11.4723-ffd_lengths0[1]/2, -4.5)
#ffd_num_points =   (3,3,3)
#ffd_delta_z0 =     (-0.5)
#ffd_delta_index =  (1,1,1)


##ffd_lengths0 =     (3, 2.5, 1) ## pod and wedge
##ffd_origin0 =      (24, 11.4723-ffd_lengths0[1]/2, -4.5) ## pod and wedge
#ffd_lengths0 =     (1, 1, 1)
#ffd_origin0 =      (0,0,0)
#ffd_num_points =   (3,3,3)
#ffd_delta_z0 =     (0)
#ffd_delta_index =  (1,1,1)


#print("origin: ", ffd_origin0)
print("length: ", length_body)
print("width: ", maxy-miny)

#lengths, origins, bumps = (10,int(length_body/12),4)
lengths, origins, bumps = (1,1,1)
#lengths, origins, bumps = (4, int(length_body - ffd_lengths0[0]), 4) ## moves along entire mesh keeping entire ffd box on mesh (not really true)
"'####### Find a way to not check Zero deformation except for the first time #######'"

## Test Case Loops
iteration = 0
skips = 0
for i in range(lengths):
    ffd_lengths = (ffd_lengths0[0] + (i*2), ffd_lengths0[1], ffd_lengths0[2])
    for j in range(origins):
        #ffd_origin = (ffd_origin0[0] + (j), ffd_origin0[1]- (ffd_lengths[1]/2), ffd_origin0[2]- (ffd_lengths[2]))
        ffd_origin = (ffd_origin0[0] + (j), ffd_origin0[1], ffd_origin0[2])
        for k in range(bumps):
            ffd_delta_z = (ffd_delta_z0 - (k*3.5)) #####!!!!!!!!  was 0.5
            
            ## Skip iteration if FFD box is too big for body
            if ffd_lengths[0] > length_body - ffd_origin[0]:
                skips += 1
                print("box location and/or size is too big to fit on the body")
                print("skipping iteration")
                continue
            
            ## Deform baseline mesh with new FFD box parameters
            vert_coords, tri_verts, comp_num = load_tri_mesh_points(baseline_tri)
            deformed_mesh_points = deform_mesh_points(ffd_lengths, ffd_origin, ffd_num_points, ffd_delta_z, ffd_delta_index, vert_coords,0,True)
            write_tri_file("studies/GoSwift/meshes/test_deformed.tri",deformed_mesh_points,tri_verts,comp_num)

            ## Run MachLine with new deformed geometry
            report = run_machline('studies/GoSwift/input_files/test.json', run = True)
            
            ## Skips folowing calculations if MachLine failed to run
            if report == None:
                print()
                print("MachLine Failed to Run.....      Skipping Iteration")
                skips += 1
                continue
            
            ## Post processing of bump location and pressure data
            ffd_box.append([ffd_lengths, ffd_origin, ffd_num_points, ffd_delta_z, ffd_delta_index])
            xg,p = pressures(gamma, p_static, density, speed_of_sound, v_inf, Mach, angles)

            ## Add smoothing functuion here (Savitzky-Golay or Exponential Smoothing)
            #plt.plot(xg,p)
            #p = savgol_filter(p, 50, 3)
            #plt.plot(xg,p)

            x_loc.append(xg)
            nearfield_sig.append(p)

            ## Run sBOOM
            g_sig, noise_level = boom(Mach, r_over_l, ref_length, altitude, N_points, angles) ## add loop here to update SBoom input flags based on mesh units 
            ground_sig.append(g_sig)
            loudness.append(noise_level)

            iteration += 1
            print()
            print("Iteration: ", iteration, "/", lengths*origins*bumps)
            print()

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
print("Number of Bump Configurations: ", (lengths*origins*bumps))
print("Skipped Cases: ", skips)
end_time = time.time()
execution_time = end_time - start_time
print("Execution Time: ", execution_time / 60 , " min")
print("------------------------------------------")
print()
