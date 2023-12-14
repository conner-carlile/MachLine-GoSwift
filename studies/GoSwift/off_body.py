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


################################## Run MachLine #############################################
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
        xg.append((points[i][0]-x0)*12)
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
    print("x: ", xg)
    print()
    print("dp: ", p)
    plt.figure(1)
    plt.plot(xg,p)
    #plt.ylim(-.015, .015)
    #plt.xlabel("(X-Xo)/L")
    #plt.ylabel("(p - pinf)/pinf")
    plt.xlabel("X")
    plt.ylabel("dP/P")
    plt.title("Near-Field Overpressure")
    plt.show()


    # Future stuff for SBoom
    #delimiter = ','
    #pressure = delimiter.join(str(p))
#
    #df = pd.read_csv("studies/GoSwift/results/off_body_pressure.csv")
    #df.to_csv("studies/GoSwift/results/off_body_pressure.dat", sep = ',' )




#======================================== SBoom Stuff ================================================
#def boom():
    
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

