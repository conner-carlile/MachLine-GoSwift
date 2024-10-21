import os
import csv
import pandas as pd
from stl import mesh
import subprocess as sp
import sys
import numpy as np
import json
import pickle
import matplotlib.pyplot as plt
from rapidboom.sboomwrapper import SboomWrapper
import pyldb

def set_json(input_file, freestream_velocity, gamma, mach, formulation,solver):
    """Writes input JSON file for MachLine"""
    
    data = {
        "flow": {
            "freestream_velocity": freestream_velocity,
            "gamma": gamma,
            "freestream_mach_number": mach
        },
        "geometry": {
            "file": "studies/GoSwift/meshes/test_deformed.tri",
            #"mirror_about" : "xz",
            #"spanwise_axis": "+y",
            "wake_model": {
                "wake_present": False,
                "append_wake": False,
                "trefftz_distance": 200.0
            }
        },
        "solver": {
            "formulation": formulation,
            "run_checks": False,
            "matrix_solver": solver
        },
        "post_processing": {
            "pressure_rules": {
                "incompressible": False,
                "isentropic": True,
                "second-order": False
            }
        },
        "output": {
            "verbose": True,
            "body_file": "studies/GoSwift/results/test.vtk",
            #"mirrored_body_file" : "studies/GoSwift/results/test_mirrored.vtk",
            #"wake_file": "studies/GoSwift/results/test_wake.vtk",
            "report_file": "studies/Goswift/results/test.json",
            "offbody_points": {
                "points_file": "studies/Goswift/off_body/off_body_sheet.csv",
                "output_file": "studies/Goswift/results/off_body_sample.csv"
            }
        }
    }
    
    with open(input_file, 'w') as json_file:
        json.dump(data, json_file, indent=4)
    
    return

def find_mins(obj):
    """ Reads in stl mesh file and finds first (x,y,z) point for centerline and length calcs"""

    global minx, miny, minz, maxx, maxy, maxz, y_pos, z_pos, width_body
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
            #y_pos = y
            z_pos = z
            y_pos = 0 #(abs(maxy) - abs(miny))/2
            #z_pos = (maxz + minz)/2
    #print("length", maxx - minx)
    width_body = maxy - miny
    return minx, miny, minz, maxx, maxy, maxz, y_pos, z_pos


def off_body(N_points,M_number,body_lengths):    
    """ Defines and creates csv file of desired off body points (line)"""

    global x0, x, xf,y, L, points
    x_start = minx
    
    xf = maxx
    #y = (miny + maxy)/2   
    y = y_pos      ########### correct center

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
    points = np.zeros((N_points, 3), dtype=float) ##points+1
    i = 0
    j = 0
    for i in range(N_points):  ## +1
        points[i][j] += (x)
        points[i][j+1] += (y)
        points[i][j+2] += (zo)
        x += ds


    with open('studies/Goswift/off_body/off_body_sheet.csv','w') as csvfile:
       writer = csv.writer(csvfile)
       writer.writerow(["x, y, z"])
       for row in points:
           writer.writerow(row)
    print("z_pos:", z_pos)
    return x0, x, xf,y, L, points



def off_body_sheet(N_points, M_number, body_lengths, num_azimuth):
    """ Defines and creates csv file of desired off body points with multiple azimuth angles"""

    global x0, x, xf, y, L, points

    # Define azimuth angles
    delta_phi = 90 / (num_azimuth + 1)

    print("delta phi: ",delta_phi)
    angles = []
    for i in range((num_azimuth * 2) + 1):
        angles.append(delta_phi * (i + 1))
 
    with open('studies/Goswift/off_body/off_body_sheet.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["x,","y,","z"]) ## messes up the csv file when reading in

        for i in range(len(angles)):
            x_start = minx
            #print("x_start: ", x_start)
            xf = maxx
            y = y_pos  # correct center
            #print("y_pos: ", y_pos)
            z = z_pos  # correct center


            L = xf - x_start  # length of geometry
            #Lu = L + (L / 2)  # length of csv line   #!!!!! pod_ensure

            R = L * body_lengths  # how far away from the body
            #ds = Lu / N_points #!!!!! undo pod_ensure
            #zo = z - R * np.sin((angles[i] - 90) * np.pi / 180)
            zo = (-R) * np.sin((angles[i]) * np.pi / 180)# +(1.837082861*12) ##!!!!!!!! change back   pod_ensure
            #print("zo: ", zo)


            # define x start and end location with mach cone offset
            mu = np.arcsin(1 / M_number)

            #x0 = R / np.tan(mu)  # will stay constant for all angles   pod_ensure

            offset = R / np.tan(mu) ###!!!!!!!
            x0 = x_start #+ 3.74 ########!!!!!!!!!!! Change back^^^ for just pod_ensure

            Lu = L + (L) + offset ### was L/2
            ds = Lu / N_points   #####^^^^ pod_ensure

            #print("x0 (trig)", x0)
            #y0 = y - R * np.cos((angles[i] - 90) * np.pi / 180)
            y0 = (R) * np.cos((angles[i]) * np.pi / 180) + y_pos# + width_body/2 #!!!!!!! pod ensure ## add back in???
            #print("y0 (trig)", y0)

            print("Line origin: ", x0,",", y0,",", zo)

            # define (x, y, z) locations for off-body control points
            # assuming x is the streamwise direction
            x = x0
            points = np.zeros((N_points, 3), dtype=float)

            for j in range(N_points): #+ 1
                points[j][0] = x
                points[j][1] = y0
                points[j][2] = zo
                x += ds
                
                # Write each point to the CSV file
                writer.writerow(points[j])
   
    return angles



################################## Run MachLine #############################################
def run_machline(input_filename, delete_output=True, run=True):
    """Runs MachLine with the given input and returns the report if MachLine generated one."""
    
    ## Delete output from previous iteration, ensures no duplicate files are created if MachLine fails to run on the current iteration
    if delete_output:
        input = json.loads(open(input_filename).read())
        output_filename = input["output"]["report_file"]
        if os.path.exists(output_filename):
            os.remove(output_filename)

    ## Run
    if run:
        sp.run(["./machline.exe", input_filename])

    ## Get report
    with open(input_filename, 'r') as input_handle:
        input_dict = json.load(input_handle)
    report_file = input_dict["output"].get("report_file")
    if report_file is not None:
        try:
            with open(report_file) as report_handle:
                report = json.load(report_handle)
        except:
            report = None
    else:
        report = None
    
    ## Delete input
    #if delete_input:
    #    os.remove(input_filename)

    return report
################################################################################################

def pressures(gamma, p_static, density, speed_of_sound, v_inf, Mach, angles): ## delete points
    """ Calculates off body pressures from MachLine calculated velocities """

    c = speed_of_sound
    #gamma = 1.4         ##### ratio of specific heats for air
    po_inf = p_static*((1 + (gamma-1)/2 * Mach**2) ** (gamma/(gamma-1)))

    #read in off body velocities from MachLine results file
    df = pd.read_csv('studies/Goswift/results/off_body_sample.csv')
    df = df[['V']]
    df.to_csv('studies/Goswift/results/off_body_velocity.csv',index=False)
    csv_file = pd.read_csv('studies/Goswift/results/off_body_velocity.csv')
    

    #calculate pressures from velocities
    global p, pressure
    p = []
    xg = []

    #for i in range(len(angles)):   #### added this....... ### Change back
    for i in range(len(angles)):  ## - 1
        for j in range(len(points)):  ## - 1 ## points >> angles
        #xg.append((points[i][0] - x0)/L)  #!
            xg.append((points[j][0]))#-x0)/L) !!!!!!!!! comment for reference -x0/L
    #print("x: ", xg)

    ## Calculate pressures from velocities
    with open('studies/Goswift/results/off_body_velocity.csv','r') as file:
        rows = csv.DictReader(file)
        for row in rows:   
            M = abs(float(row['V'])/c) #### multiply by new calc velocity
            pl = ( po_inf/(1 + ((gamma-1)/2)*(M)**2)**(gamma/(gamma-1)) )
            p.append (((pl-p_static)/p_static))
    
    with open('studies/Goswift/results/off_body_pressure.csv','w',newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["x_loc, Pressure"])
        writer.writerows(zip(xg,p))
    return xg, p


#############  working azimuth angles, just need to add angle constraint ###########
def boom(MACH, r_over_l, ref_length, altitude, num_points, angles):
    """ Runs sBoom and calculates ground signature and perceived loudness """
    
    print()
    print("Running sBoom...")
    data = np.genfromtxt('studies/GoSwift/results/off_body_pressure.csv', delimiter=',', skip_header=1)
    PROP_R = r_over_l*ref_length # *3.28084

    #sig = np.asarray([data[:,0], data[:,1]]).T
    
    sig_list = []
    for i in range(len(angles)):
        start_index = i * num_points
        end_index = (i+1) * num_points
        sig = np.asarray([data[start_index:end_index,0], data[start_index:end_index,1]]).T
        sig_list.append(sig)

    
    g_sig = []
    noise_level = []
    for j in range(len(sig_list)):
        sig = sig_list[j]
        angle = (angles[j]-90)

        _sboom = SboomWrapper('./temp', 'sboom.exe')
        _sboom.set(mach_number=MACH,
                    altitude=altitude,
                    propagation_start=PROP_R,
                    altitude_stop=0,
                    output_format=0,  ########  0 = dp vs ms        ########### 1 = ft vs dp/P
                    input_xdim=0,       ## 1 = inches, 0 = ft
                    num_azimuthal = 1,
                    azimuthal_angles = angle,
                    #propagation_points=40000, #add zero
                    #padding_points=8000
                    propagation_points=20000, #add zero
                    padding_points=5000
                    )

        _sboom.set(signature=sig, input_format=0) # set input signature and format
        sboom_results = _sboom.run(atmosphere_input=None) # run sBOOOM
        g_sig_single = sboom_results["signal_0"]["ground_sig"]

        # np.savetxt(dp_directory + label + 'g_sig.txt', g_sig)

        noise_level_single = pyldb.perceivedloudness(g_sig_single[:, 0],
                                                    g_sig_single[:, 1],)
        g_sig.append(g_sig_single)
        noise_level.append(noise_level_single)

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


def clear_pickle_file(filename):
    # Open in 'wb' mode to overwrite the file with nothing
    with open(filename, 'wb') as pickle_file:
        pass  # Do nothing, just clear the file


def append_to_pickle(filename, data):
    # Check if the file exists and is not empty
    if os.path.exists(filename) and os.path.getsize(filename) > 0:
        try:
            with open(filename, 'rb') as pickle_file:
                existing_data = pickle.load(pickle_file)
        except EOFError:
            # File is empty, proceed with an empty list
            existing_data = []
    else:
        # File does not exist or is empty, proceed with an empty list
        existing_data = []

    # Append the new data
    existing_data.append(data)
    
    # Write the updated data back to the pickle file
    with open(filename, 'wb') as pickle_file:
        pickle.dump(existing_data, pickle_file)
        pickle_file.flush()  # Flush the data to disk
        os.fsync(pickle_file.fileno())  # Ensure the OS writes the file to disk


def heatmap(ffd_box, loudness2):
    # Extract x_location, bump_height, and loudness_2 from ffd_box
    x_locations = [item[1][0] for item in ffd_box]
    bump_heights = [item[3] for item in ffd_box]


    # Reshape data for heatmap grid
    x_unique = list(set(x_locations))
    y_unique = list(set(bump_heights))

    # Create a 2D grid of loudness values
    loudness_grid = np.zeros((len(y_unique), len(x_unique)))

    x_locations = np.array(x_locations)
    bump_heights = np.array(bump_heights)
    loudness2 = np.array(loudness2)

    for i, y in enumerate(y_unique):
        for j, x in enumerate(x_unique):
            mask = (bump_heights == y) & (x_locations == x)
            loudness_grid[i, j] = loudness2[mask]

    # Plot heatmap
    #plt.figure(figsize=(10, 6))
    X, Y = np.meshgrid(x_unique, y_unique)
    plt.pcolormesh(x_unique, y_unique, loudness_grid, shading='auto', cmap='coolwarm')
    plt.colorbar(label='Loudness')
    plt.xlabel('X Location')
    plt.ylabel('Bump Height')
    plt.title('Bump Length = 6ft')
    plt.show()

'''___ Test Run Script / Check ___'''

if __name__ == "__main__":

    input_file = "studies/GoSwift/input_files/test.json"
    baseline_stl = "studies/GoSwift/meshes/test_sw.stl"

    input_file = "studies/GoSwift/input_files/test.json"
    mach = 1.8
    #input = json.loads(open(input_file).read())
    #Mach = input["flow"]["freestream_mach_number"]
    set_json(input_file, [1548.928,0,0], 1.4, mach)

    N_points = 1000
    r_over_l = 3
    ref_length = 22.27682 # 154 ft for N+2,  20 ft for SH
    altitude = 40000
    p_static = 393.13 #lbf/ft^2
    density = .00058728
    speed_of_sound = 968.08
    v_inf = 1548.928
    PROP_R = r_over_l*ref_length #*3.28084  ##### add if mesh is in meters
    num_azimuth = 0

    body_mesh = mesh.Mesh.from_file(baseline_stl)
    minx, miny, minz, maxx, maxy, maxz, y_pos, z_pos = find_mins(body_mesh)
    x0, x, xf,y, L, points = off_body(N_points, mach, r_over_l)
    print("L: ", L)
    print(miny)
    print(maxy)
    #off_body(N_points,Mach,r_over_l)
    angles = off_body_sheet(N_points, mach, r_over_l, num_azimuth)

    #xg,p = pressures(243.61, .00036392, 968.08, 1548.928, Mach ) 
    
    print(angles)