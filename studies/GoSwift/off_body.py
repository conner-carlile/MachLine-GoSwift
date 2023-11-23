import os
import csv
import pandas as pd
from stl import mesh
import subprocess as sp
import sys
import numpy as np
import json
import matplotlib.pyplot as plt
from sboomwrapper import SboomWrapper
#from case_running_functions import run_machline


def find_mins(obj):

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
    print("length", maxx - minx)

# keep new y and z just shift x over outside mach cone
def off_body(N_points,M_number):    
    global x0, x, xf,y, L, points
    x_start = minx
    xf = maxx
    y = (miny + maxy)/2   
    #y = y_pos       ###########

    z = minz
    #z = z_pos ############

    #z = (maxz + minz)/2        # could miss pressure data by having center higher than nose center?????
    
    L = xf-x_start              # length of geometry
    Lu =L + (L/2)               # length of csv line
    
    R = L * 3                   # how far away from body #
    zo = z-R
    ds = Lu/N_points
    

    #define x start and end location with mach cone offset
    mu = np.arcsin(1/M_number)
    x0 = R/np.tan(mu)      #####-.5*L
    
    #print('xo: ', x0)
    #print('new X: ', x)


    # define (x,y,z) locations for off body control points
    #assuming x is stream wise direction
    x = x0 
    points = np.zeros((N_points+1, 3), dtype=float)
    i = 0
    j = 0
    for i in range(N_points+1): 
        points[i][j] += (x)
        points[i][j+1] +=  (y)
        points[i][j+2] +=  (zo)
        x = x + ds

    print(points)
    print('maxz', maxz)
    print('minz', minz)
    print('z', z)
    #print("x",x)
    #print("L", L)
    #print("ds", ds)
    #print("N: ", N_points)
    #print("ymin, ymax: ", miny, maxy)

    with open('studies/Goswift/off_body/off_body.csv','w') as csvfile:
       writer = csv.writer(csvfile)
       for row in points:
           writer.writerow(row)


################################################################################################
def run_machline(input_filename, delete_input=False, run=True):
    """Runs MachLine with the given input and returns the report if MachLine generated one."""

    # Run
    if run:
        sp.run(["./machline.exe", input_filename])

    # Get report
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

    # Delete input
    #if delete_input:
    #    os.remove(input_filename)

    return report
################################################################################################

def pressures(speed_of_sound, p_static):
    c = speed_of_sound
    gamma = 1.4         ##### ratio of specific heats for air
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
    ps = []
    p = []
    xg = []
    for i in range(len(points)-1):
        xg.append((points[i][0] - x0)/L)
    #print("x: ", xg)

    with open('studies/Goswift/results/off_body_velocity.csv','r') as file:
        rows = csv.DictReader(file)
        for row in rows:   
            M = abs(float(row['V'])/c)
            #print(M)
            pl = ( po_inf/(1 + ((gamma-1)/2)*(M)**2)**(gamma/(gamma-1)) )
            p.append (((pl-p_static)/p_static))
        

           
    with open('studies/Goswift/results/off_body_pressure.csv','w') as file:
        writer = csv.writer(file)
        writer.writerow((p))
        #print(file)

    # plot pressures
    #xg = []
    #for i in range(len(points)-1): 
    #    xg.append(((points[i][0]-x0)/L))
    #
    #print(p_boom)
    #plt.plot(xg,p)

    plt.plot(xg,p)
    plt.xlabel("(X-Xo)/L")
    plt.ylabel("(p - pinf)/pinf")
    #print((p))
    plt.show()

    delimiter = ','
    pressure = delimiter.join(str(p))

    df = pd.read_csv("studies/GoSwift/results/off_body_pressure.csv")
    df.to_csv("studies/GoSwift/results/off_body_pressure.dat", sep = ',' )



if __name__ == '__main__':

    N_points = 500
    Mach = 1.5
    #gamma = 1.4
    r_over_l = 3
    ref_length = 200 ###27.432 # for X-59   ########update for n+1
    altitude = 53200
    PROP_R = r_over_l*ref_length*3.28084  ##### ??? 3.2

    #body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/sears_haack_9800.stl')
    #body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/100_30_SH.stl')
    #body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/pod(SH).stl')
    #body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/pod_vsp_smooth.stl')
    #body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/test.stl')
    #body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/pod_vsp_bump.stl')
    #body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/SH_half.stl')
    body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/n+2.stl')
    #body_mesh = mesh.Mesh.from_file('studies/GoSwift/meshes/delta_wing.stl')
    #body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/low_boom.stl')
    #body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/low_boom2.stl')
    #body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/biswis.stl')
    #body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/JWB_0AOA.stl')
    #body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/test_sw.stl')
# --------------------------------------------------------------------------------------------------
    #body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/pod_smooth.stl')
    #body_mesh = mesh.Mesh.from_file('studies/Goswift/meshes/pod_bump.stl')

    find_mins(body_mesh)
    off_body(N_points,Mach)  #N_points, Mach Number

    #run_machline('studies/GoSwift/input_files/pod_smooth.json')
    #run_machline('studies/GoSwift/input_files/pod_bump.json')
#-----------------------------------------------------------------------------------------------------
    #run_machline('studies/GoSwift/input_files/SH_input.json')
    #run_machline('studies/GoSwift/input_files/SH_bump_input.json')

    run_machline('studies/GoSwift/input_files/n+2_input.json')
    #run_machline('studies/GoSwift/input_files/input.json')
    #run_machline('studies/GoSwift/input_files/low_boom_input.json')
    #run_machline('studies/GoSwift/input_files/jaxa_input.json')
    #run_machline('studies/GoSwift/input_files/test.json')

    pressures(968.08, 243.61 ) #### run at 50000 ft
    #pressures(295.07, 14170 ) #### run at 1400 m 5219


    #data =  p #"studies/GoSwift/results/off_body_pressure.dat"  #### data should be of type [x, dp]

    #data = np.genfromtxt( dp_directory)
    #conv = 1/0.0254                     ######dimensiolalize x points so find base mesh units
    #sig = np.asarray([(data[3])/conv, data[4]]).T
    #
#
    #sboom = SboomWrapper('./temp', 'sboom.exe')
    #sboom.set(mach_number=1.5,
    #        altitude=altitude,
    #        propagation_start=PROP_R,
    #        altitude_stop=0,
    #        output_format=0,
    #        input_xdim=2,
    #        azimuthal_angles = 0,
    #        propagation_points=40000,
    #        padding_points=8000)
    #sboom.set(signature=sig, input_format=0) # set input signature and format
    #sboom_results = sboom.run(atmosphere_input=None) # run sBOOOM
    #g_sig = sboom_results["signal_0"]["ground_sig"]
    #print('pass')
    #print(sboom_results)

