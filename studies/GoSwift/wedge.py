import json
import time
import pickle
#import numpy as np
from stl import mesh
from off_body import *
from deform_tri_module import *
from scipy.signal import savgol_filter


input_file = "studies/GoSwift/input_files/test.json"

input = json.loads(open(input_file).read())
MACH = input["flow"]["freestream_mach_number"]

#N_points = [1200, 1400, 1800, 2700, 3500] # starting number of points
N_points = [3500, 3500, 3500, 3500, 3500]
r_over_l = [.25, .5, 1, 2, 3] 

ref_length = 13000 #mm 
altitude = 40000
p_static = 393.13 #lbf/ft^2
density = .00058728
speed_of_sound = 968.08
v_inf = 1548.928
#PROP_R = r_over_l*ref_length #*3.28084  ##### add if mesh is in meters
num_azimuth = 0 


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


## Create CSV file for off body points
body_mesh = mesh.Mesh.from_file(baseline_stl)
minx, miny, minz, maxx, maxy, maxz, y_pos, z_pos = find_mins(body_mesh)
print("minx: ", minx)
print("maxx: ", maxx)
print("miny: ", miny)
print("maxy: ", maxy)
print("minz: ", minz)
print("maxz: ", maxz)
print("y_pos: ", y_pos)
print("z_pos: ", z_pos)
#off_body(N_points,MACH,r_over_l) 
#angles = off_body_sheet(N_points, MACH, r_over_l, num_azimuth)

print("width: ", maxy-miny)

def off_body_sheet(N_points, M_number, body_lengths, num_azimuth):
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

        for i in range(len(N_points)):
            x_start = minx
            print("x_start: ", x_start)
            xf = maxx
            y = y_pos  # correct center
            print("y_pos: ", y_pos)
            z = z_pos  # correct center

            L = xf - x_start  # length of geometry
            #Lu = L + (L / 2)  # length of csv line   #!!!!! pod_ensure

            R = L * body_lengths[i]  # how far away from the body
            #ds = Lu / N_points #!!!!! undo pod_ensure
            #zo = z - R * np.sin((angles[i] - 90) * np.pi / 180)
            zo = (-R) * np.sin((angles[0]) * np.pi / 180)# +(1.837082861*12) ##!!!!!!!! change back   pod_ensure
            print("zo: ", zo)


            # define x start and end location with mach cone offset
            mu = np.arcsin(1 / M_number)

            #x0 = R / np.tan(mu)  # will stay constant for all angles   pod_ensure

            #offset = R / np.tan(mu) ###!!!!!!!
            offset = 48710.9843875075
            x0 = x_start #+ 3.74 ########!!!!!!!!!!! Change back^^^ for just pod_ensure

            Lu = L + (L/2) + offset

            print("offset: ", offset)
            ds = Lu / N_points[i]   #####^^^^ pod_ensure

            #print("x0 (trig)", x0)
            #y0 = y - R * np.cos((angles[i] - 90) * np.pi / 180)
            y0 = (R) * np.cos((angles[0]) * np.pi / 180)# + width_body/2 #!!!!!!! pod ensure
            #print("y0 (trig)", y0)


            # define (x, y, z) locations for off-body control points
            # assuming x is the streamwise direction
            x = x0
            points = np.zeros((N_points[i], 3), dtype=float)

            for j in range(N_points[i]): #+ 1
                points[j][0] = x
                points[j][1] = y0
                points[j][2] = zo
                x += ds
                
                # Write each point to the CSV file
                writer.writerow(points[j])

angles = off_body_sheet(N_points, MACH, r_over_l, num_azimuth)  ### angle stuff


length_body = maxx - minx
print(length_body)
width_body = maxy - miny

## Convert baseline .stl file to .tri file for FFD
stl_to_tri(baseline_stl, baseline_tri)


## Initialize FFD box info lists
length = []
origin = []
num_points = []
delta_z = []
delta_index = []


ffd_lengths0 = (1,1,1)
#ffd_origin0 = (5941.16668, -207.018, -1250.0356)#(0,-0.5,-2)  ## need to shift x and y origin points by half of their length value (origin is corner of box not center)(-2 z to shift box down to place bump on bottom)
ffd_origin0 = (0,0,0)#(0,-0.5,-2)  ## need to shift x and y origin points by half of their length value (origin is corner of box not center)(-2 z to shift box down to place bump on bottom)
ffd_num_points = (3,3,1) ## I keep this constant
ffd_delta_z0 = (0)  ## change ffd_deltz_z to ensure correct step size
ffd_delta_index = (1,1,1)

#print("origin: ", ffd_origin0)
#print("length: ", length_body)
#print("width: ", maxy-miny)

#lengths, origins, bumps = (10,int(length_body/12),4)
lengths, origins, bumps = (1,1,1)


vert_coords, tri_verts, comp_num = load_tri_mesh_points(baseline_tri)
deformed_mesh_points = deform_mesh_points(ffd_lengths0, ffd_origin0, ffd_num_points, ffd_delta_z0, ffd_delta_index, vert_coords,0,True)
write_tri_file("studies/GoSwift/meshes/test_deformed.tri",deformed_mesh_points,tri_verts,comp_num)
#report = run_machline('studies/GoSwift/input_files/test.json')


xg,p = pressures(p_static, density, speed_of_sound, v_inf, MACH, angles,points) ## delete "points" from pressure in off body module



#for i in range(len(N_points)):
phat = savgol_filter(p, 100, 3)
#plt.plot(x[start_index:end_index+1], [yi - subtract_value for yi in y[start_index:end_index+1]])

x_new = [xg,xg,xg,xg,xg]

#with open('studies/Goswift/results/wedge_off_body.csv','w',newline='') as file:
#    writer = csv.writer(file)
#    writer.writerow(["pressure"])
#    for item in phat:
#        writer.writerow((item))

pressure = np.array(phat)
np.savetxt("studies/Goswift/results/wedge_off_body_p.csv", pressure, delimiter=",")

print(len(x_new))
print(len(phat))
print("!!!!!!!!!!!")


fig, axs = plt.subplots(figsize=(8, 6))
x = np.linspace(0, 1, len(phat))
plt.plot(x,phat)
#plt.plot(x,p)
plt.show()


x1 = xg[0:N_points[0]]
x2 = xg[0:N_points[1]]
x3 = xg[0:N_points[2]]
x4 = xg[0:N_points[3]]
x5 = xg[0:N_points[4]]


#y1 = [y-0.3250 for y in phat[0:N_points[0]]]
#y2 = [i-0.6500 for i in phat[N_points[0]+1:N_points[1]+N_points[0]+1]]
#y3 = [j-.93000 for j in phat[N_points[0]+N_points[1]+1+1:N_points[2]+N_points[1]+N_points[0]+1+1]]
#y4 = [k-1.26000 for k in phat[N_points[2]+N_points[1]+N_points[0]+1+1:N_points[3]+N_points[2]+N_points[1]+N_points[0]+1+1]]
#y5 = [c-1.539000 for c in phat[N_points[3]+N_points[2]+N_points[1]+N_points[0]:N_points[4]+N_points[3]+N_points[2]+N_points[1]+N_points[0]+1]]
y1 = [y - 0.3250 for y in phat[0:N_points[0]]]
y2 = [i - 0.6500 for i in phat[N_points[0]+1:N_points[0]+N_points[1]+1]]
y3 = [j - 0.93000 for j in phat[N_points[0]+N_points[1]+2:N_points[0]+N_points[1]+N_points[2]+2]]
y4 = [k - 1.26000 for k in phat[N_points[0]+N_points[1]+N_points[2]+3:N_points[0]+N_points[1]+N_points[2]+N_points[3]+3]]
y5 = [c - 1.539000 for c in phat[N_points[0]+N_points[1]+N_points[2]+N_points[3]:N_points[0]+N_points[1]+N_points[2]+N_points[3]+N_points[4]]]

print(N_points[0]+N_points[1]+N_points[2]+N_points[3]+N_points[4])


plt.plot(x1,y1)
plt.plot(x2,y2)
plt.plot(x3,y3)
plt.plot(x4,y4)
plt.plot(x5,y5)

plt.show()


#axs.plot(xg[0:N_points[0]+1],[y-.03250 for y in phat[0:N_points[0]+1]]) #x_loc[0:1000][set],   
#axs.plot(xg[0:N_points[1]],[i-.06500 for i in phat[N_points[0]:N_points[1]+N_points[0]]])
#axs.plot(xg[0:N_points[2]],[j-.13000 for j in phat[N_points[1]+1:N_points[2]+N_points[1]+1]])
#axs.plot(xg[0:N_points[3]],[k-.26000 for k in phat[N_points[2]:N_points[3]+N_points[2]]])
#axs.plot(xg[0:N_points[4]],[c-.39000 for c in phat[N_points[3]+1:N_points[4]+N_points[3]+1]])
#plt.show()

#axs.set_xlabel('X Label')
#axs.set_ylabel('Y Label')
#axs.set_title('Overlayed Subplots')
#axs.legend()

# Show plot
#plt.show()




#x_loc.append(xg)
#nearfield_sig.append(p)
print("Done")



    