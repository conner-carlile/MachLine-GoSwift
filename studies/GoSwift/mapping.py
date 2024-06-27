import pickle
import csv
import numpy as np
from off_body import *
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import PillowWriter
from scipy.signal import savgol_filter
from exponentalSmoothing import exponential_smoothing
from rapidboom.sboomwrapper import SboomWrapper
import pyldb

# Open data files
#with open('studies/Goswift/results/temp_x_loc.pkl', 'rb') as pickle_file:
#    x_loc = pickle.load(pickle_file)
#
#with open('studies/Goswift/results/temp_nearfield.pkl', 'rb') as pickle_file:
#    nearfield_sig  = pickle.load(pickle_file)
#
#with open('studies/Goswift/results/pod_wedge_x_loc.csv', 'w', newline='') as csvfile:
#    writer = csv.writer(csvfile)
#    writer.writerow(['x_loc'])
#    writer.writerows([[x] for x in x_loc[0]])
#
#with open('studies/Goswift/results/pod_wedge_nearfield.csv', 'w', newline='') as csvfile:
#    writer = csv.writer(csvfile)
#    writer.writerow(['nearfield_sig'])
#    writer.writerows([[x] for x in nearfield_sig[0]])


## Read in a CSV file with 2 columns of data x,p
#x_loc = []
with open('studies/Goswift/results/pod_wedge_x_loc.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    next(reader)  # Skip the header row if it exists
    x_loc = [float(row[0]) for row in reader]

with open('studies/Goswift/results/pod_wedge_nearfield.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    next(reader)  # Skip the header row if it exists
    nearfield_sig = [float(row[0]) for row in reader]
print("x_loc: ", x_loc)


sav = savgol_filter(nearfield_sig, 100, 3)
exp = exponential_smoothing(nearfield_sig, 0.05)

#plt.plot(x_loc, nearfield_sig)
#plt.title('Raw Data')
#plt.show()
#plt.plot(x_loc, sav)
#plt.title('Savitzky-Golay Filter')
#plt.show()
#plt.plot(x_loc, exp)
#plt.title('Exponential Smoothing')
#plt.show()

sav_exp = savgol_filter(exp, 100, 3)
exp_sav = exponential_smoothing(sav, 0.05)

#plt.plot(x_loc[0], sav_exp)
#plt.title('Savitzky-Golay Filter on Exponential Smoothing')
#plt.show()
#plt.plot(x_loc[0], exp_sav)
#plt.title('Exponential Smoothing on Savitzky-Golay Filter')
#plt.show()

## Define offbody points and their correct Beta angles
x_2 = [4.198692441,
5.013218237,
11.8256158,
12.6401416,
13.1584762,
15.97229258,
18.04563097,
23.67326374,
28.78256192,
31.30018711,
34.55829029,
36.70567648,
38.85306267,
44.77688664,
47.96094202,
50.18237601]

## Observed beta angles
beta = [0.7343305552,
0.845666837,
0.845666837,
0.845666837,
0.7343305552,
0.9347043834,
0.5951635418,
0.5951635418,
1.006472864,
1.064990429,
0.9347043834,
0.9347043834,
0.9347043834,
1.006472864,
1.113307846,
0.5951635418]

## Manually corrected beta angles
#beta = [0.7343305552,
#0.845666837,
#0.845666837,
#0.845666837,
#0.7343305552,
#0.9347043834,
#0.5951635418,
#0.5951635418,
#0.7,
#0.775,
#0.775,
#0.775,
#0.775,
#0.675,
#0.675,
#0.5951635418]

nearfield_sig = exp_sav


## Loop through nearfield signature and find nodes that are closest to the points in x_2
nodes_x = []
for value2 in x_2:
    closest_value = min(x_loc, key=lambda x:abs(x-value2))
    closest_index = x_loc.index(closest_value)
    nodes_x.append(closest_value)
    print(value2, closest_value, closest_index)
#plt.plot(x_2, nodes_x)
#plt.title("nodes")
#plt.show()


## Interpolate angles based on the nodes
angles_interpolated = np.interp(x_loc, nodes_x , beta)
print("angles_interpolated: ", angles_interpolated)
#plt.plot(x_loc[0], angles_interpolated[0])
#plt.title('Interpolated Angles')
#plt.show()


## Calculate x offset for every point
mu = np.arcsin(1/1.6)
x_offset = []
#R = 1.30675909 # meters
R = 4.287267 # ft

for i in range(len(angles_interpolated)):
    x_offset.append(abs((R/np.tan(angles_interpolated[i])) - (R/np.tan(mu))))

print("x_offset: ", x_offset)
plt.plot(x_loc, x_offset)
plt.title('X Offset')
plt.show()
#print("x_loc[0]: ", x_loc[0])

## Calculate new x values
x_new = []
for i in range(len(x_offset)):
    x_new.append(x_loc[i] - x_offset[i])
    #print(x_loc[0][i], x_offset[i], x_new[i])

plt.plot(x_loc, nearfield_sig)
plt.plot(x_new, nearfield_sig)
plt.title("Original vs Corrected Nearfield Signature")
plt.legend(['Original', 'Corrected'])
plt.xlabel('X (ft)')
plt.ylabel('Pressure (Pa)')
plt.show()
print(x_loc[-1] - x_loc[0], x_new[-1]-x_new[0])

# Write list to a CSV file
#with open('studies/Goswift/results/output_x_new.csv', 'w', newline='') as csvfile:
#    writer = csv.writer(csvfile)
#    writer.writerow(['x_new', 'nearfield_sig'])
#    for i in range(len(x_new)):
#        writer.writerow([x_new[i]-2.788, nearfield_sig[0][i]])
#
#with open('studies/Goswift/results/output_x_loc.csv', 'w', newline='') as csvfile:
#    writer = csv.writer(csvfile)
#    writer.writerow(['x_loc', 'nearfield_sig'])
#    for i in range(len(x_loc[0])):
#        writer.writerow([x_loc[0][i]-2.788, nearfield_sig[0][i]])

# Read in a CSV file with 2 columns of data x,p
data = []
with open('studies/Goswift/results/Untitled spreadsheet - Sheet5.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    next(reader)  # Skip the header row if it exists
    for row in reader:
        x = float(row[0])
        p = float(row[1])
        data.append((x, p))

p_new = []
for i in range(len(nearfield_sig)):
    p_new.append(nearfield_sig[i]*18823.1661 + 18823.1661)


plt.plot([x * 3.281 for x, _ in data], [p for _, p in data])
plt.plot(x_new, p_new)
plt.plot(x_loc, p_new)
plt.legend(['UNS3D', 'MachLine-Corrected', 'MachLine'])
plt.xlabel('X (ft)')
plt.ylabel('Pressure (Pa)')
plt.show()


#data = np.genfromtxt('studies/GoSwift/results/off_body_pressure.csv', delimiter=',', skip_header=1)
r_over_l = 3
ref_length = 42.65 # ft
altitude = 40000
MACH = 1.6

PROP_R = R # *3.28084
#sig = np.asarray([data[:,0], data[:,1]]).T


sig_list = [np.asarray([x_new, nearfield_sig]).T]
#sig_list = np.array([x_new, nearfield_sig[0]]).T
print("sig_list: ", sig_list)
#g_sig = [] ###
#noise_level = [] ###


g_list = [] ###
noise_list = [] ###
rl = []
sig = sig_list[0]
print("sig_list shape: ", sig_list[0].shape)
print("sig_list values: ", sig_list[0])
#angle = (angles[j]-90)
for i in range(1): ###
    g_sig = []
    noise_level = []
    _sboom = SboomWrapper('./temp', 'sboom.exe')
    _sboom.set(mach_number=MACH,
                altitude=altitude,
                propagation_start= PROP_R,
                altitude_stop= altitude-(ref_length*3),  #*(i+1)*.11),
                output_format=0,    ## 0 =  regular time (ms) vs. dp (psf)  ,   1 = X(feet) vs. dp/P
                input_xdim=0,       ## 0 = ft, 1 = in, 2 = m
                #num_azimuthal = 1,
                #azimuthal_angles = 0,
                propagation_points=20000, #add zero
                padding_points=4000
                )
    _sboom.set(signature=sig, input_format=0) # set input signature and format
    sboom_results = _sboom.run(atmosphere_input=None) # run sBOOOM
    #print("sBOOM parameters: ", _sboom._parameters)
    g_sig_single = sboom_results["signal_0"]["ground_sig"]
    # np.savetxt(dp_directory + label + 'g_sig.txt', g_sig)
    noise_level_single = pyldb.perceivedloudness(g_sig_single[:, 0],
                                                g_sig_single[:, 1],)
    rl.append(round((ref_length*(i+1)*.15 / ref_length), 3))
    g_sig.append(g_sig_single)
    noise_level.append(noise_level_single)
    g_list.append(g_sig) ###
    noise_list.append(noise_level) ###
###########################################################

plt.plot(g_sig[0][:,0], g_sig[0][:,1]) ###
#plt.xlim(.007,.055)
plt.show() ###



def animate(i):
    plt.cla()
    plt.plot(g_list[i][0][:, 0], g_list[i][0][:, 1])
    #plt.title("Ground Signature")
    plt.ylabel("dp (psf)")
    plt.xlabel("X (ft)")
    plt.ylim(-100, 100)
    plt.xlim(0, 60)
    plt.legend(["R/L: "+ str(rl[i])], loc='upper right')
    plt.tight_layout()

## Create the animation object
#fig = plt.figure()
#ani = animation.FuncAnimation(fig, animate, frames=len(g_list), interval=300, repeat=True)
##
#gif_writer = PillowWriter(fps=1)  # Set frames per second (fps)
#ani.save("animation.gif", writer=gif_writer)
#
## Display the animation
#plt.show()
