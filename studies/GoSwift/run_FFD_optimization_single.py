from deform_tri_module import *
from FFD_optimization_wrapper_NEW import *
import numpy as np


final_delta_Z =[]
x_origin = []
z_origin = []

baseline_filename = 'c608.tri'  
original_mesh_points, tri_verts, component_num = load_tri_mesh_points(tri_filename = baseline_filename)

case_ID = 999
print('Case: ', case_ID)
print('Mach, AoA : ', 1.4, 2.15)

'''SET FOR TARGET DEFORMATIONS'''
bump_x_locations = [472.44, 708.66, 984.25]
target_heights = [2.0, 2.0, 2.0]

# bump_x_locations = [X1[i], X2[i], X3[i]]
# target_heights = [H1[i], H2[i], H3[i]]

        
'''SHOULD BE FIXED VALUES (X,Z), MAYBE BASED ON BUMP LOCATIONS (Y)'''
# bump_x_lengths = [60.0, 60.0, 60.0]
# bump_y_lengths = [40.0, 40.0, 40.0]
# bump_z_lengths = [10.0, 10.0, 10.0]

bump_x_lengths = [118.11, 118.11, 118.11]
bump_y_lengths = [28.0, 28.0, 28.0]
bump_z_lengths = [10.0, 10.0, 10.0]

'''IDEALLY ARE FIXED VALUES TO SIMPLIFY INPUTS'''
num_x_nodes = [3, 3, 3]
num_y_nodes = [3, 3, 3]
num_z_nodes = [3, 3, 3]
bump_nodes = [num_x_nodes, num_y_nodes, num_z_nodes]

test_class = FFD_deform_opt(baseline_mesh = original_mesh_points,
                            baseline_tri_verts = tri_verts,
                            baseline_component_num = component_num,
                            set_case_number = case_ID,
                            set_target_heights = target_heights, 
                            set_x_locations = bump_x_locations, 
                            set_x_lengths = bump_x_lengths, 
                            set_y_lengths = bump_y_lengths,
                            set_z_lengths = bump_z_lengths, 
                            set_bump_nodes = bump_nodes,
                            write_delta_mesh = True)

check_return = test_class.run_FFD_optimization()


final_delta_Z.append(check_return)
x_origin.append(test_class._origin_x_locations)
z_origin.append(test_class._origin_z_locations)