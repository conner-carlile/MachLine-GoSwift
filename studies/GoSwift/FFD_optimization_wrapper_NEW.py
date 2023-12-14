import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from scipy import interpolate
from deform_tri_module import *
import shutil
from scrape_FIU_deformation import scrape_deformation_data
import time

class FFD_deform_opt:

    def __init__(self, baseline_mesh, baseline_tri_verts, baseline_component_num,
                 set_case_number, set_target_heights, set_x_locations,
                 set_x_lengths, set_y_lengths, set_z_lengths, set_bump_nodes,
                 write_delta_mesh):
        
        self._case_number = set_case_number
        self._target_heights = set_target_heights
        self._FFD_x_lengths = set_x_lengths
        self._FFD_y_lengths = set_y_lengths
        self._FFD_z_lengths = set_z_lengths
        self.write_delta_mesh = write_delta_mesh
        
        self._keel_points = np.genfromtxt('c608_keel_line_points.txt')
        self._keel_interp = interpolate.CubicSpline(self._keel_points[:,0], self._keel_points[:,1], bc_type=('natural','natural'), extrapolate = False)
        
        mid_FFD_z = []
        self._origin_x_locations = []
        self._origin_y_locations = []
        self._origin_z_locations = []
        
        '''Set the x, y, and z origins based on middle point location and FFD lengths'''
        for i in range(len(set_x_locations)):
            mid_FFD_z.append(self._keel_interp(set_x_locations[i]) - 0.25) # offset the z distance of the FFD mid point by 0.25 inches down
            self._origin_x_locations.append(set_x_locations[i] - (self._FFD_x_lengths[i]/2))
            self._origin_z_locations.append(mid_FFD_z[i] - (self._FFD_z_lengths[i]/2))
            self._origin_y_locations.append(-(self._FFD_y_lengths[i]/2))

        self._num_x_nodes = set_bump_nodes[0]
        self._num_y_nodes = set_bump_nodes[1]
        self._num_z_nodes = set_bump_nodes[2]
        
        self._original_mesh_points = np.copy(baseline_mesh)   
        self._tri_verts = np.copy(baseline_tri_verts) 
        self._component_num = np.copy(baseline_component_num)
        self._final_deformed_mesh_points = np.copy(self._original_mesh_points)

    def objective_function(self, delta_z):
                
        deformed_mesh_points, index_deformed_points = deform_mesh_points(ffd_lengths = self.current_lengths, ffd_origin = self.current_origin, 
                                                                         ffd_num_points = self.current_num_points, ffd_delta_z = delta_z, 
                                                                         ffd_delta_index = self.current_delta_node_index, mesh_points = self._current_mesh_points)
        
        check_current_mesh_points = self._current_mesh_points[index_deformed_points]
        check_deformed_mesh_points = deformed_mesh_points[index_deformed_points]
        delta_mesh_points = check_deformed_mesh_points[:,2] - check_current_mesh_points[:,2]

        max_delta = abs(max(abs(delta_mesh_points)) - (abs(self.current_target_height)))
                
        return max_delta
    
    def deform_final_mesh(self, delta_z):

        self._final_deformed_mesh_points, self._final_index_deformed_points = deform_mesh_points(ffd_lengths = self.current_lengths, 
                                                                                                 ffd_origin = self.current_origin, 
                                                                                                 ffd_num_points = self.current_num_points, 
                                                                                                 ffd_delta_z = delta_z, 
                                                                                                 ffd_delta_index = self.current_delta_node_index, 
                                                                                                 mesh_points = self._final_deformed_mesh_points)
   
    def save_final_mesh(self, write_delta_mesh = False):
        # write_tri_file('case' + str(case_ID_list[j]).zfill(3), deformed_mesh_points = final_unrotated_mesh_points, tri_verts = tri_verts, comp_num = component_num)
        write_tri_file('case' + str(self._case_number).zfill(3) + 'A', deformed_mesh_points = self._final_deformed_mesh_points, tri_verts = self._tri_verts, comp_num = self._component_num)
        write_tecplot_file('case' + str(self._case_number).zfill(3) + 'A', deformed_mesh_points = self._final_deformed_mesh_points, tri_verts = self._tri_verts, comp_num = self._component_num)
        
        if write_delta_mesh == True:
            delta_mesh = np.copy(self._original_mesh_points)
            delta_mesh[:,2] = self._final_deformed_mesh_points[:,2] - self._original_mesh_points[:,2]
            write_tecplot_delta_z_file('case' + str(self._case_number).zfill(3) + 'A_delta_mesh', deformed_mesh_points = self._final_deformed_mesh_points, delta_mesh_points = delta_mesh, tri_verts = self._tri_verts, comp_num = self._component_num)
            
    def run_FFD_optimization(self):
        
        delta_z_FFD = []
        
        for i in range(len(self._num_x_nodes)):
                            
            self.current_lengths = [self._FFD_x_lengths[i], self._FFD_y_lengths[i], self._FFD_z_lengths[i]]
            self.current_origin = [self._origin_x_locations[i], self._origin_y_locations[i], self._origin_z_locations[i]]
            self.current_num_points = [self._num_x_nodes[i], self._num_y_nodes[i], self._num_z_nodes[i]]
            self.current_target_height = self._target_heights[i]
            self.current_delta_node_index = [2,2,2]

            print('FFD Origin (x,y,z):', self.current_origin)
            print('FFD Lengths (x,y,z):', self.current_lengths)
            print('FFD Num Points (x,y,z):', self.current_num_points)
            print('FFD Mid Node Index:', self.current_delta_node_index)
            
            current_deformed_mesh_points, index_current_points = deform_mesh_points(ffd_lengths = self.current_lengths, 
                                                                                    ffd_origin = self.current_origin, 
                                                                                    ffd_num_points = self.current_num_points, 
                                                                                    ffd_delta_z = 0.0, 
                                                                                    ffd_delta_index = self.current_delta_node_index,
                                                                                    mesh_points = self._original_mesh_points,
                                                                                    deformation_num = i + 1,
                                                                                    vtk_box_flag = True)
            
            self._current_mesh_points = np.copy(self._original_mesh_points[index_current_points])
            
            max_bound = 14.0
            target_sign = np.sign(self.current_target_height)
            
            if target_sign > 0.0:
                bounds = ((0.0, max_bound),)
            elif target_sign < 0.0:
                bounds = ((-max_bound, 0.0),)
                     
            init_params = [1.0*target_sign]
            print('Starting Optimization')
            optimum = optimize.minimize(self.objective_function, x0 = init_params, method='SLSQP', bounds = bounds, options={ 'disp': True, 'ftol':1e-10,'maxiter': 50, 'eps': 1.4901161193847656e-05})
        
            print('Optimum Delta_z_FFD: ', optimum.x)
            delta_z_FFD.append(optimum.x)
            self.deform_final_mesh(optimum.x)
            
        self.save_final_mesh(write_delta_mesh = self.write_delta_mesh)
        
        return delta_z_FFD