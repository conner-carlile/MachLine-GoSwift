import numpy as np
#import pygem
from pygem.ffd import FFD
import pandas as pd
import pyevtk.hl as vtk
from scipy import interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def load_tri_mesh_points(tri_filename):
    
    nVerts, nTris = pd.read_table(tri_filename, nrows = 1, delim_whitespace = True)
    nVerts = int(nVerts)
    nTris = int(nTris)
    
    vert_coords = pd.read_table(tri_filename, nrows = nVerts, header = None, index_col = False, skiprows = 1, delim_whitespace = True).to_numpy()
    tri_verts = pd.read_table(tri_filename, nrows = nTris, skiprows = nVerts, delim_whitespace = True).to_numpy()
    comp_num = pd.read_table(tri_filename, nrows = nTris, skiprows = nVerts + nTris, delim_whitespace = True).to_numpy()
    
    return vert_coords, tri_verts, comp_num

def deform_mesh_points(ffd_lengths, ffd_origin, ffd_num_points, ffd_delta_z, ffd_delta_index, mesh_points, deformation_num = 0, vtk_box_flag = False):
    
    # Fixed distribution of control points
    numx, numy, numz = ffd_num_points
    
    z_deform = np.zeros((numx, numy, numz))
    
    ffd = FFD([numx, numy, numz])
    

    # Ensure indices are within bounds *******************************************
    valid_indices = (
        min(max(ffd_delta_index[0], 0), numx - 1),
        min(max(ffd_delta_index[1], 0), numy - 1),
        min(max(ffd_delta_index[2], 0), numz - 1)
        )
    
    z_deform[valid_indices[0], valid_indices[1], valid_indices[2]] = ffd_delta_z    
    #***************************************
    #z_deform[ffd_delta_index[0],ffd_delta_index[1],ffd_delta_index[2]] = ffd_delta_z
    
    ffd.box_length = ffd_lengths
    ffd.box_origin = ffd_origin
    ffd.array_mu_z = z_deform
    
    #print(deform_mesh_points)

    #deformed_mesh_points, index_deformed_points = ffd(mesh_points)     #original#
    deformed_mesh_points = ffd(mesh_points)
    
    # ffd.write_parameters(filename = 'test_params.txt')
    
    if vtk_box_flag == True:
        control_points = ffd.control_points()
        #print("control points:", control_points)
        #fig = plt.figure(100)
        #ax = fig.add_subplot(111, projection='3d')
        #ax.scatter(control_points[0], control_points[1], control_points[2], s=50, c='red')
        #plt.show()
        c_p_b = ffd.control_points(deformed=True)  ## True origionaly
        vtk.gridToVTK('./cont_p_base', c_p_b[:,0], c_p_b[:,1], c_p_b[:,2])
        vtk.gridToVTK('./FFD_control_points_' + str(deformation_num), control_points[:,0], control_points[:,1], control_points[:,2])
        
    
    return deformed_mesh_points#, index_deformed_points

def get_keel_line_normal(keel_points, x_center):
    
    keel_interp = interpolate.CubicSpline(keel_points[:,0], keel_points[:,1], bc_type=('natural','natural'), extrapolate = False)
    dx = 1.0 #inches    
    x1 = x_center - dx
    z1 = keel_interp(x1)
    x2 = x_center + dx
    z2 = keel_interp(x2)
    dx = (x2-x1)
    dz = (z2-z1)
    
    D = np.sqrt(dx**2 + dz**2)
    N = [dz/D,-dx/D]
    theta = np.arctan((0 - dz)/(-1 - (-dx)))
    
    return N, theta

def Euler2Quat(ea):

    '''converts euler angles to quaternion using Eq. (11.7.8)

    Parameters
    ----------
    ea: array
        euler angles

    Returns
    -------
    list of quaternion values
    '''

    ea05 = ea[0]*0.5
    ea15 = ea[1]*0.5
    ea25 = ea[2]*0.5

    CP = np.cos(ea05)
    CT = np.cos(ea15)
    CS = np.cos(ea25)

    SP = np.sin(ea05)
    ST = np.sin(ea15)
    SS = np.sin(ea25)
    c1 = CP*CT
    c2 = SP*ST
    c3 = SP*CT
    c4 = CP*ST

    return [(c1*CS + c2*SS), (c3*CS - c4*SS),
            (c4*CS + c3*SS), (c1*SS - c2*CS)]

def NormQuat2(q):

    '''normalizes quaternion using Eq. (11.10.5) from Phillips

    Parameters
    ----------
    q: array
        quaternion values

    Returns
    -------
    list of normalized quaternion values
    '''
    q0, q1, q2, q3 = q

    norm_den = np.sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3)

    return [q0/norm_den, q1/norm_den, q2/norm_den,
            q3/norm_den]

def rotatePoint(v, e):

    '''converts earth fixed vector to body fixed vector using reduced form 
    of Eq. (11.6.8)

    Parameters
    ----------
    v: array
        vector to be converted
    e: array
        quaternion

    Returns
    -------
    list of body fixed vector components
    '''

    v0, v1, v2 = v
    e0, e1, e2, e3 = e

    ve00 = v0*e0
    ve01 = v0*e1
    ve02 = v0*e2
    ve03 = v0*e3
    ve10 = v1*e0
    ve11 = v1*e1
    ve12 = v1*e2
    ve13 = v1*e3
    ve20 = v2*e0
    ve21 = v2*e1
    ve22 = v2*e2
    ve23 = v2*e3

    return [(e0*(ve00 + ve13 - ve22) - 
                  e1*(-ve01 - ve12 - ve23) -
                  e2*(ve02 - ve11 + ve20) +
                  e3*(-ve03 + ve10 + ve21)),
                (e0*(-ve03 + ve10 + ve21) + 
                 e1*(ve02 - ve11 + ve20) -
                 e2*(-ve01 - ve12 - ve23) -
                 e3*(ve00 + ve13 - ve22)),
                (e0*(ve02 - ve11 + ve20) - 
                 e1*(-ve03 + ve10 + ve21) +
                 e2*(ve00 + ve13 - ve22) -
                 e3*(-ve01 - ve12 - ve23))]

def rotate_mesh_points(x_rotation, y_rotation, z_rotation, mesh_points):
    
    '''converts earth fixed vector to body fixed vector using reduced form 
    of Eq. (11.6.8)

    Parameters
    ----------
    v: array
        vector to be converted
    e: array
        quaternion

    Returns
    -------
    list of body fixed vector components
    '''
    
    quat_angles = Euler2Quat([x_rotation, y_rotation, z_rotation])
    quat_angles = NormQuat2(quat_angles)
    rotated_mesh_points = np.zeros([len(mesh_points[:,0]), len(mesh_points[0,:])])
    
    for i in range(len(mesh_points[:,0])):
        rotated_mesh_points[i,:] = rotatePoint(v = mesh_points[i,:], e = quat_angles)
        
    return rotated_mesh_points

def scale_mesh_points(x_scale, y_scale, z_scale, mesh_points):
    
    scaled_mesh_points = np.zeros([len(mesh_points[:,0]), len(mesh_points[0,:])])
    
    for i in range(len(mesh_points[:,0])):
        scaled_mesh_points[i,:] = np.asarray([mesh_points[i,0]*x_scale, mesh_points[i,1]*y_scale, mesh_points[i,2]*z_scale,])
        
    return scaled_mesh_points

# def write_FFD_2_vtk(case_ID, control_points):
#         vtk.gridToVTK('./' + str(case_ID), control_points[:,0], control_points[:,1], control_points[:,2])

def write_tri_file(new_filename, deformed_mesh_points, tri_verts, comp_num):
    
    print('Starting write-to-file (.tri):')
    nVerts = len(deformed_mesh_points[:,0])
    nTris = len(tri_verts[:,0])
    
    #with open(new_filename + '.tri', 'w') as export_handle:
    with open(new_filename, 'w') as export_handle:

        
        # Write header
        print("{0:<18}{1:<18}".format(nVerts, nTris), file=export_handle)
            
        for vertex in deformed_mesh_points:
            print("{0:<18.10}{1:<18.10}{2:<18.10}".format(*vertex), file=export_handle)
        
        for tri_int in tri_verts:
            print("{0:<12}{1:<12}{2:<12}".format(*tri_int), file=export_handle)
            
        for comp_num in comp_num:
            print("{0:<4d}".format(int(comp_num[0])), file=export_handle)



def write_tecplot_file(new_filename, deformed_mesh_points, tri_verts, comp_num):
    
    print('Starting write-to-file (Tecplot):')
    nVerts = len(deformed_mesh_points[:,0])
    nTris = len(tri_verts[:,0])
    
    with open(new_filename + '.dat', 'w') as export_handle:
        
        # Write header
        print('VARIABLES = X Y Z', file=export_handle)
        print('ZONE T="TEMP",N=', nVerts,',E=', nTris,',DATAPACKING=POINT, ZONETYPE=FETRIANGLE', file=export_handle)

        for vertex in deformed_mesh_points:
            print("{0:<18.10}{1:<18.10}{2:<18.10}".format(*vertex), file=export_handle)
        
        for tri_int in tri_verts:
            print("{0:<12}{1:<12}{2:<12}".format(*tri_int), file=export_handle)
            
def write_tecplot_delta_z_file(new_filename, deformed_mesh_points, delta_mesh_points, tri_verts, comp_num):
    
    print('Starting write-to-file (Tecplot):')
    nVerts = len(deformed_mesh_points[:,0])
    nTris = len(tri_verts[:,0])
    
    with open(new_filename + '.dat', 'w') as export_handle:
        
        # Write header
        print('VARIABLES = X Y Z DELTA_Z', file=export_handle)
        print('ZONE T="TEMP",N=', nVerts,',E=', nTris,',DATAPACKING=POINT, ZONETYPE=FETRIANGLE', file=export_handle)

        for i, vertex in enumerate(deformed_mesh_points):
            print("{0:<18.10}{1:<18.10}{2:<18.10}{3:<18.10}".format(vertex[0], vertex[1], vertex[2], delta_mesh_points[i,2]), file=export_handle)
        
        for tri_int in tri_verts:
            print("{0:<12}{1:<12}{2:<12}".format(*tri_int), file=export_handle)



if __name__ == "__main__":
    #location of undeformed geometry .tri file
    file = "studies/GoSwift/meshes/test_sw.tri"

    #Define ffd box parameters
    #[(12, 32.67716, 32.67716), (25, -16.338579177856445, -48.983498), (3, 3, 3), -1.25, (1, 1, 1)]
    #(12, 32.67716, 32.67716) , (24, 0, -16.30634) , -0.25
#[(24, 32.66058, 32.66058), (36, 0, -16.33746), (3, 3, 3), -0.25, (1, 1, 1)]
    #[(12, 32.66058, 32.66058), (58, 0, -16.33746), (3, 3, 3), -0.25, (1, 1, 1)]
    #[(5, 2.722535, 2.722535), (8.638409614562988, 0, -1.361479), (3, 3, 3), 0.5, (1, 1, 1)]

    origin:  (5941.16668, 0, 1250.0356)
    length:  267.3218
    #width:  32.66058
    width = 22.96588
    #[(2, 1.5, 1), (21.2, 10.7223, -4.5), (3, 3, 3), -1.5, (1, 1, 1)]
    ffd_lengths =       (3, 9.681, 2)
     ## need to shift x and y origin points by half of their length value (origin is corner of box not center)(-2 z to shift box down to place bump on bottom)
    #ffd_origin =        (24, 11.4723-ffd_lengths[1]/2, -4.5) ## pod aand wedge
    ffd_origin =        (.68898-ffd_lengths[0]/2, width/2 -ffd_lengths[1]/2, -ffd_lengths[2]-2) ## pod aand wedge
    #ffd_origin =        (9.36841, 2.722535/2-ffd_lengths[1]/2, -1) ## pod aand wedge 2-2.722535/2
    ffd_num_points =    (3,3,3) ## I keep this constant
    ffd_delta_z =       (20)    #0-3
    ffd_delta_index =   (1,1,1)  #Any other vector 


    vert_coords, tri_verts, comp_num = load_tri_mesh_points(file)
    deformed_mesh_points = deform_mesh_points(ffd_lengths, ffd_origin, ffd_num_points, ffd_delta_z, ffd_delta_index, vert_coords,0,True)
    write_tri_file("studies/GoSwift/meshes/test_deformed.tri",deformed_mesh_points,tri_verts,comp_num)

    ## The following is only to plot the original and deformed geometries and does not effect results

    ##original mesh
    print(vert_coords)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #fig.tight_layout()
    x = vert_coords[:,0]
    y = vert_coords[:,1]
    z = vert_coords[:,2]
    ax.scatter(x,y,z)
    #deformed mesh
    a = deformed_mesh_points[:,0]
    b = deformed_mesh_points[:,1]
    c = deformed_mesh_points[:,2]
    ax.scatter(a,b,c)
    fig.tight_layout()
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_aspect("equal")
    # plot ffd box corners
    ax.scatter(ffd_origin[0],ffd_origin[1],ffd_origin[2], color = "red", s = 20)
    ax.scatter(ffd_origin[0]+ffd_lengths[0],ffd_origin[1],ffd_origin[2], color = "blue", s = 20)
    ax.scatter(ffd_origin[0],ffd_origin[1]+ffd_lengths[1],ffd_origin[2], color = "blue", s = 20)
    ax.scatter(ffd_origin[0],ffd_origin[1],ffd_origin[2]+ffd_lengths[2], color = "blue", s = 20)
    ax.scatter(ffd_origin[0]+ffd_lengths[0],ffd_origin[1]+ffd_lengths[1],ffd_origin[2], color = "blue", s = 20)
    ax.scatter(ffd_origin[0],ffd_origin[1]+ffd_lengths[1],ffd_origin[2]+ffd_lengths[2], color = "blue", s = 20)
    ax.scatter(ffd_origin[0]+ffd_lengths[0],ffd_origin[1],ffd_origin[2]+ffd_lengths[2], color = "blue", s = 20)
    ax.scatter(ffd_origin[0]+ffd_lengths[0],ffd_origin[1]+ffd_lengths[1],ffd_origin[2]+ffd_lengths[2], color = "blue", s = 20)
    plt.show()
           