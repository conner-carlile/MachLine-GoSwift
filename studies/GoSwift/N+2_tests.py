import pickle
import csv
import numpy as np
from off_body import *
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.signal import savgol_filter
from exponentalSmoothing import exponential_smoothing

## Open pickle files
with open('studies/Goswift/results/temp_ffd_box.pkl', 'rb') as pickle_file:
    ffd_box = pickle.load(pickle_file)

with open('studies/Goswift/results/temp_x_loc.pkl', 'rb') as pickle_file:
    x_loc = pickle.load(pickle_file)

with open('studies/Goswift/results/temp_nearfield.pkl', 'rb') as pickle_file:
    nearfield_sig  = pickle.load(pickle_file)

with open('studies/Goswift/results/temp_ground.pkl', 'rb') as pickle_file:
    ground_sig = pickle.load(pickle_file)

with open('studies/Goswift/results/temp_loudness.pkl', 'rb') as pickle_file:
    loudness = pickle.load(pickle_file)


def reduce_peaks(data, scale_factor):
    return np.sign(data) * (np.abs(data) ** scale_factor)

print((ffd_box))

## Parameters
L = 202.0
body_lengths = 3
altitude = 40000
M_number = 1.8
R = L * body_lengths 
mu = np.arcsin(1/M_number)
x0 = R/np.tan(mu)  


## Shift and scale x-axis
x_norm = []
for i in range(len(x_loc[0])):
    x_norm.append((x_loc[0][i] - x0)*12-1000) # convert to in., shift by 1000 in.


## split and smooth all case signatures 
sig = []
for i in range(len(nearfield_sig)):
    nearfield_sig[i] = exponential_smoothing(nearfield_sig[i], 0.2)
    nearfield_sig[i] = savgol_filter(nearfield_sig[i], 5, 3)    

    front_sig = nearfield_sig[i][:1348]
    front_x = x_norm[:1348]

    aft_sig = nearfield_sig[i][1348:]
    aft_x = x_norm[1348:]

    aft_smoothed = reduce_peaks(aft_sig, 2) #1.5
    smoothed = np.concatenate((front_sig,aft_smoothed))

    # split again
    front_sig = smoothed[:1365]
    front_x = x_norm[:1365]
    aft_sig = smoothed[1365:]
    aft_x = x_norm[1365:]

    aft_smoothed = exponential_smoothing(aft_sig, 0.03)
    smoothed = np.concatenate((front_sig,aft_smoothed))
    sig.append(smoothed)
print(len(sig))
for i in range(len(sig)):
    plt.plot(x_norm, sig[i])
    plt.title("Near-Field R/L = 3")
    plt.xlabel("X (in)")
    plt.ylabel("dp / P")
    plt.grid(True, linestyle=':')
    plt.xlim(0,4000)
    plt.ylim(-.02,.02)
plt.show()


## Run sBOOM on smoothed signatures
angles = [1]
num_points = len(x_norm)
PROP_R = body_lengths*L
angle = 0

g_sig = []
loudness2 = []
for i in range(len(sig)):
    sig_zip = list(zip(x_norm, sig[i]))
    _sboom = SboomWrapper('./temp', 'sboom.exe')
    _sboom.set(mach_number=M_number,
                altitude=altitude,
                propagation_start=PROP_R,
                altitude_stop=0,
                output_format=0,  ########  0 = dp vs ms        ########### 1 = ft vs dp/P
                input_xdim=1,       ## 1 = inches, 0 = ft
                num_azimuthal = 1,
                azimuthal_angles = angle,
                #propagation_points=40000, #add zero
                #padding_points=8000
                propagation_points=20000, #add zero
                padding_points=4000
                )

    _sboom.set(signature=sig_zip, input_format=0) # set input signature and format
    sboom_results = _sboom.run(atmosphere_input=None) # run sBOOOM
    g_sig_single = sboom_results["signal_0"]["ground_sig"]

    noise_level_single = pyldb.perceivedloudness(g_sig_single[:, 0],
                                                g_sig_single[:, 1],)
    g_sig.append(g_sig_single)
    loudness2.append(noise_level_single)

for i in range(len(g_sig)):
    plt.plot(g_sig[i][:, 0], g_sig[i][:, 1])
    plt.title("Ground Signature")
    plt.xlabel("Time (s)")
    plt.ylabel("dp / P")
    plt.grid(True, linestyle=':')
    #plt.xlim(0, 10)
    #plt.ylim(-.02,.02)
plt.show()
## calculate delta PLDB
baseline_PLDB = 96.25375755910656
delta_PLDB = []
for i in range(len(loudness2)):
    delta_PLDB.append(loudness2[i] - baseline_PLDB)

print("Un-Smoothed, Smoothed, Difference")
for i in range(len(loudness2)):
    print((loudness[i][0], loudness2[i], loudness[i][0] - loudness2[i]))
print(delta_PLDB)


## Find quietest configuration multiple azimuth
min = min(loudness2)
for i in range(len(loudness2)):
    if loudness2[i] == min:
        print("The iteration with the lowest Loudness is: ", i)
        print("The lowest loudness is: ", min)
        break


## plots
loc = np.linspace(0,len(loudness2),len(loudness2))
plt.plot(loc,loudness2)
plt.title("Loudness vs Iteration")
plt.xlabel("Iteration")
plt.ylabel("Loudness (db)")
plt.show()


x = []
for i in range(len(ffd_box)):
    x.append(ffd_box[i][1][0])

plt.scatter(x,loudness2)
plt.title("Loudness vs location")
plt.xlabel("Location Along Body (ft)")
plt.ylabel("Loudness (db)")
plt.show()


y = []
for i in range(len(ffd_box)):
    y.append(ffd_box[i][0][0])

plt.scatter(y,loudness2)
plt.title("Loudness vs Bump Length")
plt.xlabel("Length")
plt.ylabel("Loudness (db)")
plt.show()


z = []
for i in range(len(ffd_box)):
    z.append(ffd_box[i][3])

plt.scatter(z,loudness2)
plt.title("Loudness vs Bump Shape")
plt.xlabel("Delta-z")
plt.ylabel("Loudness (db)")
plt.show()
    
plt.scatter(x,delta_PLDB)
plt.title("Delta PLDB vs Location")
plt.xlabel("Location Along Body (ft)")
plt.ylabel("Delta PLDB")
plt.show()

#save data as csv files
with open('studies/Goswift/results/loudness2.csv', mode='w') as file:
    writer = csv.writer(file)
    writer.writerow(loudness2)
with open('studies/Goswift/results/delta_PLDB.csv', mode='w') as file:
    writer = csv.writer(file)
    writer.writerow(delta_PLDB)
with open('studies/Goswift/results/ffd_box.csv', mode='w') as file:
    writer = csv.writer(file)
    writer.writerow(ffd_box)
#with open('studies/Goswift/results/x_loc.csv', mode='w') as file:
#    writer = csv.writer(file)
#    writer.writerow(x_loc)


## Heatmap
origins = []
for i in range(len(ffd_box)):
    origins.append(ffd_box[i][1][1])
x = np.array(x)

heights = []
for i in range(len(ffd_box)):
    heights.append(ffd_box[i][3])
y = np.array(heights)

loudness2 = np.array(loudness2)

heatmap_data = np.zeros((len(origins), len(heights)))
for i in range(len(origins)):
    for j in range(len(heights)):
        for k in range(len(ffd_box)):
            if origins[k] == origins[i] and heights[k] == heights[j]:
                heatmap_data[i][j] = loudness2[k]

plt.imshow(heatmap_data, cmap='bwr', origin='lower')
plt.colorbar(label='Loudness (db)', cmap='bwr')
plt.title('Loudness Heatmap')
plt.xlabel('Bump Origin')
plt.ylabel('Bump Height')
plt.show()



## plot loudness vs bump location

## plot loudness vs bump length

## plot loudness vs bump height

## plot quietest nearfield signature

## print ffd box parameters for quietest case