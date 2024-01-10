import os
import csv
import pandas as pd
from stl import mesh
import subprocess as sp
import sys
import numpy as np
import json
import matplotlib.pyplot as plt
from rapidboom.sboomwrapper import SboomWrapper
import pyldb
#from case_running_functions import run_machline


def find_mins(obj):
    """ Reads in stl mesh file and finds first (x,y,z) point for centerline and length calcs"""

    global minx, miny, minz, maxx, maxy, maxz, y_pos, z_pos
    minx = obj.x.min()
    maxx = obj.x.max()
    miny = obj.y.min()
    maxy = obj.y.max()
    minz = obj.z.min()
    maxz = obj.z.max()

    # Finds y and z location that corresponds to minimum x value
    verts = obj.vectors.reshape(-1,3)
    for row in verts:
        x,y,z = row
        if x == minx:
            y_pos = y
            z_pos = z
    #print("length", maxx - minx)
    return minx, miny, minz, maxx, maxy, maxz, y_pos, z_pos

## keep new y and z just shift x over outside mach cone *****????
def off_body(N_points,M_number,body_lengths):    
    """ Defines and creates csv file of desired off body points"""
    global x0, x, xf,y, L, points
    x_start = minx
    xf = maxx
    #y = (miny + maxy)/2   
    y = y_pos       ########### correct center

    #z = minz
    z = z_pos ############ correct center
    #z = (maxz + minz)/2     # could miss pressure data by having center higher than nose center?????
    
    L = xf-x_start           # length of geometry
    Lu =L + (L/2)            # length of csv line
    
    R = L * body_lengths     # how far away from body 
    zo = z-R
    ds = Lu/N_points 
    

    #define x start and end location with mach cone offset
    mu = np.arcsin(1/M_number)
    x0 = R/np.tan(mu)      


    # define (x,y,z) locations for off body control points
    #assuming x is stream wise direction
    x = x0  ###### - 5
    points = np.zeros((N_points+1, 3), dtype=float)
    i = 0
    j = 0
    for i in range(N_points+1): 
        points[i][j] += (x)
        points[i][j+1] += (y)
        points[i][j+2] += (zo)
        x += ds

    #print(points)
    #print('maxz', maxz)
    #print('minz', minz)
    #print('z', z)
    #print("x",x)
    #print("L", L)
    #print("ds", ds)
    #print("N: ", N_points)
    #print("ymin, ymax: ", miny, maxy)

    with open('studies/Goswift/off_body/off_body.csv','w') as csvfile:
       writer = csv.writer(csvfile)
       for row in points:
           writer.writerow(row)

    return x0, x, xf,y, L, points



#def off_body_sheet(N_points, M_number, body_lengths, num_azimuth):
#    global x0, x, xf, y, L, points
#
#    # Define azimuth angles
#    delta_phi = 90 / (num_azimuth + 1)
#    print(delta_phi)
#    angles = []
#    for i in range((num_azimuth * 2) + 1):
#        angles.append(-delta_phi * (i + 1))
#
#    with open('studies/Goswift/off_body/off_body_sheet.csv', 'a', newline='') as csvfile:
#        writer = csv.writer(csvfile)
#
#        for i in range(len(angles)):
#            x_start = minx
#            xf = maxx
#            y = y_pos  # correct center
#            z = z_pos  # correct center
#
#            L = xf - x_start  # length of geometry
#            Lu = L + (L / 2)  # length of csv line
#
#            R = L * body_lengths  # how far away from the body
#            ds = Lu / N_points
#            #zo = z - R * np.sin((angles[i] - 90) * np.pi / 180)
#            zo = z - np.sin((angles[i] - 90) * np.pi / 180)
#
#            # define x start and end location with mach cone offset
#            mu = np.arcsin(1 / M_number)
#            x0 = R / np.tan(mu)  # will stay constant for all angles
#            #y0 = y - R * np.cos((angles[i] - 90) * np.pi / 180)
#            y0 = y - np.cos((angles[i] - 90) * np.pi / 180)
#
#            # define (x, y, z) locations for off-body control points
#            # assuming x is the streamwise direction
#            points = np.zeros((N_points + 1, 3), dtype=float)
#
#            for j in range(N_points + 1):
#                points[j][0] = x0
#                points[j][1] = y0
#                points[j][2] = zo
#                x0 += ds
#
#            # Write all points for the current azimuth angle
#            for point in points:
#                writer.writerow(point)
#
#        return angles
import csv
import numpy as np

def off_body_sheet(N_points, M_number, body_lengths, num_azimuth):
    global x0, x, xf, y, L, points

    # Define azimuth angles
    delta_phi = 90 / (num_azimuth + 1)
    print(delta_phi)
    angles = []
    for i in range((num_azimuth * 2) + 1):
        angles.append(delta_phi * (i + 1))

    with open('studies/Goswift/off_body/off_body_sheet.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        for i in range(len(angles)):
            x_start = minx
            xf = maxx
            y = y_pos  # correct center
            z = z_pos  # correct center

            L = xf - x_start  # length of geometry
            Lu = L + (L / 2)  # length of csv line

            R = L * body_lengths  # how far away from the body
            ds = Lu / N_points
            #zo = z - R * np.sin((angles[i] - 90) * np.pi / 180)
            zo = z - np.sin((angles[i] - 90) * np.pi / 180)


            # define x start and end location with mach cone offset
            mu = np.arcsin(1 / M_number)
            x0 = R / np.tan(mu)  # will stay constant for all angles
            #y0 = y - R * np.cos((angles[i] - 90) * np.pi / 180)
            y0 = y - np.cos((angles[i] - 90) * np.pi / 180)


            # define (x, y, z) locations for off-body control points
            # assuming x is the streamwise direction
            points = np.zeros((N_points + 1, 3), dtype=float)

            for j in range(N_points + 1):
                points[j][0] = x0
                points[j][1] = y0
                points[j][2] = zo
                x0 += ds

                # Write each point to the CSV file
                writer.writerow(points[j])

    return angles









################################## Run MachLine #############################################
def run_machline(input_filename, delete_input=False, run=True):
    """Runs MachLine with the given input and returns the report if MachLine generated one."""

    # Run
    if run:
        sp.run(["./machline.exe", input_filename])

    ## Get report
    #with open(input_filename, 'r') as input_handle:
    #    input_dict = json.load(input_handle)
    #report_file = input_dict["output"].get("report_file")
    #if report_file is not None:
    #    try:
    #        with open(report_file) as report_handle:
    #            report = json.load(report_handle)
    #    except:
    #        report = None
    #else:
    #    report = None
#
    ## Delete input
    ##if delete_input:
    ##    os.remove(input_filename)
#
    #return report
################################################################################################

def pressures(p_static, density, speed_of_sound, v_inf, Mach):
    """ Calculates pressures from MachLine Calculated Velocities """

    c = speed_of_sound
    gamma = 1.4         ##### ratio of specific heats for air
    #Mach = v_inf / c
    po_inf = p_static*((1 + (gamma-1)/2 * Mach**2) ** (gamma/(gamma-1)))

    #read in off body velocities from MachLine results file
    df = pd.read_csv('studies/Goswift/results/off_body_sample.csv')
    df = df[['V']]
    df.to_csv('studies/Goswift/results/off_body_velocity.csv',index=False)
    csv_file = pd.read_csv('studies/Goswift/results/off_body_velocity.csv')
    #print("Velocities")
    #print(csv_file)
    #print()

    #calculate pressures from velocities
    global p, pressure
    p = []
    xg = []
    for i in range(len(points)-1):
        #xg.append((points[i][0] - x0)/L)  #!
        xg.append((points[i][0]-x0)) #*12)
    #print("x: ", xg)

#========================= working old pressure calcs ===========================================
    with open('studies/Goswift/results/off_body_velocity.csv','r') as file:
        rows = csv.DictReader(file)
        for row in rows:   
            M = abs(float(row['V'])/c)
            pl = ( po_inf/(1 + ((gamma-1)/2)*(M)**2)**(gamma/(gamma-1)) )
            p.append (((pl-p_static)/p_static))
    
    with open('studies/Goswift/results/off_body_pressure.csv','w') as file:
        writer = csv.writer(file)
        #writer.writerow((xg))
        #writer.writerow((p))
        writer.writerow(["x_loc, Pressure"])
        writer.writerows(zip(xg,p))
#================================================================================================

#============================== new pressure ====================================================
    #with open('studies/Goswift/results/off_body_velocity.csv','r') as file:
    #    rows = csv.DictReader(file)
    #    for row in rows: 
    #        v = float(row['V'])
    #        #print(v)
    #        cp = (2/(gamma * Mach**2)) * (((1 + ((gamma-1)/2) * (1 - (v**2/v_inf**2)) * Mach**2))**(gamma/(gamma-1))-1)
    #        print(cp)
    #        pl = ( cp *  ((1/2)*density* v_inf**2) + p_static)
    #        p.append (((pl-p_static)/p_static))
    #        #p.append(cp)
    #with open('studies/Goswift/results/off_body_pressure.csv','w') as file:
    #   writer = csv.writer(file)
    #   writer.writerow((p))
#=================================================================================================

    # Plot pressure distribution
    #plt.figure(1)
    #plt.plot(xg,p)
    #plt.xlabel("X")
    #plt.ylabel("dP/P")
    #plt.title("Near-Field Overpressure")
    #plt.show()

    return xg, p




#======================================== SBoom Stuff ================================================
def boom(altitude, PROP_R, Mach):
    print()
    print("Running sBoom...")
    dp_directory = "studies/GoSwift"
    data = np.genfromtxt('studies/GoSwift/results/off_body_pressure.csv', delimiter=',', skip_header=1)
    #conv = 1/0.0254
    sig = np.asarray([data[:,0], data[:,1]]).T

    MACH = Mach #1.6 # set Mach for current case
    r_over_l = 3
    ref_length = 154 #update for n+1
    altitude = 50000
    PROP_R = r_over_l*ref_length# *3.28084

    _sboom = SboomWrapper('./temp', 'sboom.exe')
    _sboom.set(mach_number=MACH,
                altitude=altitude,
                propagation_start=PROP_R,
                altitude_stop=0,
                output_format=0,
                input_xdim=2,
                azimuthal_angles = 0,
                propagation_points=40000,
                padding_points=8000)

    _sboom.set(signature=sig, input_format=0) # set input signature and format
    sboom_results = _sboom.run(atmosphere_input=None) # run sBOOOM
    g_sig = sboom_results["signal_0"]["ground_sig"]

    # np.savetxt(dp_directory + label + 'g_sig.txt', g_sig)

    noise_level = pyldb.perceivedloudness(g_sig[:, 0],
                                          g_sig[:, 1],)

    #plt.plot(g_sig[:, 0],g_sig[:, 1])
    #plt.title("Ground Signature")
    #plt.show()
    print('PLdB (pressure): ', noise_level)

    return g_sig, noise_level



def stl_to_tri(stl_filename, tri_filename):
    # Load STL file
    mesh_data = mesh.Mesh.from_file(stl_filename)

    n_verts = len(mesh_data.vectors) * 3
    n_tris = len(mesh_data.vectors)

    vertices = mesh_data.vectors.reshape(-1, 3)
    tri_verts = np.arange(1, n_verts + 1).reshape(-1, 3)
    comp_num = [1]

    _write_TRI(n_verts, n_tris, vertices, tri_verts, comp_num, tri_filename)
    print(f"Conversion complete. Output written to {tri_filename}")



def _write_TRI(n_verts, n_tris, vertices, tri_verts, comp_num, tri_filename):
    '''write to .tri file'''

    with open(tri_filename, 'w') as export_handle:
        # Write header
        print("{0:<18}{1:<18}".format(n_verts, n_tris), file=export_handle)

        # Write vertices
        for vertex in vertices:
            print("{0:<18.10f}{1:<18.10f}{2:<18.10f}".format(*vertex), file=export_handle)

        # Write triangle vertices
        for tri_int in tri_verts:
            print("{0:<12}{1:<12}{2:<12}".format(*tri_int), file=export_handle)

        # Write component number
        for comp in comp_num:
            print("{0:<4d}".format(int(comp)), file=export_handle)



def convert_tri_to_stl(tri_filename, stl_filename='output.stl'):
    # Read data from the .tri file
    with open(tri_filename, 'r') as tri_file:
        lines = tri_file.readlines()

    # Extract data from the .tri file
    nVerts, nTris = map(int, lines[0].split())
    vertices = [list(map(float, line.split())) for line in lines[1:nVerts + 1]]
    tri_verts = [list(map(int, line.split())) for line in lines[nVerts + 1:nVerts + nTris + 1]]

    # Write data to the .stl file
    with open(stl_filename, 'w') as stl_file:
        stl_file.write("solid\n")

        for tri_int in tri_verts:
            stl_file.write(" facet normal  {:e} {:e} {:e}\n".format(0.0, 0.0, 0.0))
            stl_file.write("   outer loop\n")

            for vertex_index in tri_int:
                vertex = vertices[vertex_index]
                stl_file.write("     vertex {:e} {:e} {:e}\n".format(*vertex))

            stl_file.write("   endloop\n")
            stl_file.write(" endfacet\n")

        stl_file.write("endsolid\n")


if __name__ == "__main__":

    input_file = "studies/GoSwift/input_files/test.json"
    baseline_stl = "studies/GoSwift/meshes/test_sw.stl"

    input = json.loads(open(input_file).read())
    Mach = input["flow"]["freestream_mach_number"]
    N_points = 1000
    r_over_l = 3
    ref_length = 20 # 154 ft for N+2,  20 ft for SH
    altitude = 50000
    PROP_R = r_over_l*ref_length

    body_mesh = mesh.Mesh.from_file(baseline_stl)
    minx, miny, minz, maxx, maxy, maxz, y_pos, z_pos = find_mins(body_mesh)
    #off_body(N_points,Mach,r_over_l)
    angles = off_body_sheet(N_points, Mach, r_over_l, 3)

    print(angles)