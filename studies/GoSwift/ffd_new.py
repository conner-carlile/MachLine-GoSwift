import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from scipy import interpolate
from deform_tri_module import *
import os
import shutil
#from scrape_FIU_deformation import scrape_deformation_data
import time

'This is designed to optimize the FFD_delta_z value to give a specific bump height'

#location of undeformed geometry .tri file
file = "studies/GoSwift/meshes/test.tri"
vert_coords, tri_verts, comp_num = load_tri_mesh_points(file)


origin:  (5941.16668, 0, 1250.0356)
length = 202.0
width = 86.42333

#ffd_lengths1 =       (6.777472119,2,7)
#ffd_origin1 =        (50.38572558, -ffd_lengths1[1]/2, -ffd_lengths1[2]) ## pod aand wedge
#ffd_num_points1 =    (3,3,3) ## I keep this constant
#ffd_delta_z1 =       (.1)    #0-3
#ffd_delta_index1 =   (1,1,1)  #Any other vector 
#target_difference =   0.8954960416 /12

#ffd_lengths1 =       (10.76224884,2,7)
#ffd_origin1 =        (79.18922309, -ffd_lengths1[1]/2, -ffd_lengths1[2]) ## pod aand wedge
#ffd_num_points1 =    (3,3,3) ## I keep this constant
#ffd_delta_z1 =       (1)    #0-3
#ffd_delta_index1 =   (1,1,1)  #Any other vector
#target_difference =  -4.607475831 /12

ffd_lengths1 =       (4.787476729,2,7)
ffd_origin1 =        (158.3793047, -ffd_lengths1[1]/2, -ffd_lengths1[2]) ## pod aand wedge
ffd_num_points1 =    (3,3,3) ## I keep this constant
ffd_delta_z1 =       (-0.08)    #0-3
ffd_delta_index1 =   (1,1,1)  #Any other vector 
target_difference =   4.917918652 /12


def optimize_deformation(ffd_lengths, ffd_origin, ffd_num_points, ffd_delta_z, ffd_delta_index, baseline_points, target_distance):
    def calculate_deformed_points(ffd_delta_z):
        # Deform the mesh using the current ffd_delta_z
        return deform_mesh_points(ffd_lengths, ffd_origin, ffd_num_points, ffd_delta_z, ffd_delta_index, baseline_points, ffd_delta_z, False)

    def objective(x):
        # Calculate the deformed points using the deformation value x
        deformed_points = calculate_deformed_points(x)

        # Calculate the distances between the deformed points
        distances = [( baseline_points[i][2] - deformed_points[i][2]) for i in range(len(deformed_points))]
        
        ## Working
        #max_distance = np.max(distances) if target_distance >= 0 else np.min(distances)
        if target_distance >= 0:
            max_distance = np.max(distances)
        else:
            max_distance = np.min(distances)
        
        print("Max distance between points:", max_distance*12)
        # Return the difference from the target distance
        return abs(max_distance - target_distance)

    # Perform optimization
    #result = optimize.differential_evolution(objective, bounds=[(-100.0, 100.0)], strategy='best1bin')
    result = optimize.differential_evolution(
        objective,
        bounds=[(-10.0, 10.0)],
        strategy='best1bin',
        tol=1e-10,
        polish=True,
        maxiter=500
    )
    print("Optimization result:", result)
    if result.success:
        # Obtain the final deformed points using the optimized ffd_delta_z
        optimized_ffd_delta_z = result.x[0]  # Get the optimized value for ffd_delta_z
        final_deformed_points = calculate_deformed_points(optimized_ffd_delta_z)  # Calculate points only once
        max_distance = objective(optimized_ffd_delta_z)

        return final_deformed_points, optimized_ffd_delta_z, max_distance  # Return both deformed points and the optimized value
    else:
        print("Optimization failed:", result.message)
        return None

#target_difference = 2.914687854/12.0
deformed_mesh_points,optimized_ffd_delta_z, max_distance = optimize_deformation(ffd_lengths1, ffd_origin1, ffd_num_points1, ffd_delta_z1, ffd_delta_index1, vert_coords, target_difference)
print("Optimized ffd_delta_z:", optimized_ffd_delta_z)
print("Convergence distance:", max_distance)
write_tri_file("studies/GoSwift/meshes/optimization_test.tri",deformed_mesh_points,tri_verts,comp_num)

