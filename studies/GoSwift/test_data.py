import pickle
import csv
import numpy as np
import matplotlib.pyplot as plt
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



## Find quietest configuration multiple azimuth
min = float('inf')
min_zero_azimuth = None

for iteration in loudness:
    center_index = len(iteration) // 2
    center_value = iteration[center_index] ## this defines the center of the list, update to handle variable length lists
    if center_value < min:
        min = center_value
        min_zero_azimuth = iteration

for i in range(len(loudness)):
    if loudness[i] == min_zero_azimuth:
        print("The iteration with the smallest center number is:", i)
        print("The smallest center number is:", min)
        break
#print("The iteration with the smallest center number is:", min_zero_azimuth)
#print("The smallest center number is:", min)




## Plotting nearfield signature and ground signature
plt.plot(x_loc[0],nearfield_sig[0])
plt.title("Nearfield Signature (raw)")
plt.ylabel("dp / P")
plt.xlabel("X (ft)")
plt.xlim(80,115)
#ax[0].set_ylim(-0.1,0.15, 0.5)
plt.show()

for i in range(len(loudness)):
#plt.plot(g_sig[0][:,0],g_sig[0][:,1])
    plt.plot(ground_sig[i][0][:,0],ground_sig[i][0][:,1])
    #plt.plot(ground_sig[i][0][:,0],ground_sig[i][0][:,1])
    #plt.show()
plt.title("Ground Signature")
plt.ylabel("dp (psf)")
plt.xlabel("time (ms)")
plt.show()
#!!!!!!!!!!!!!!!


#Smooth Data     # Apply Savitzky-Golay filter with a 200 window to the entire data set
for i, sig in enumerate(nearfield_sig):
    globals()[f'y{i+1}'] = sig

for i in range(len(nearfield_sig)):
    plt.plot(x_loc[0],globals()[f'y{i+1}'])
plt.show()    
y1 = nearfield_sig[0]
#y2 = nearfield_sig[1]
#y3 = nearfield_sig[2]
#y4 = nearfield_sig[3]
#y5 = nearfield_sig[4]
#y6 = nearfield_sig[5]
#y7 = nearfield_sig[6]


#plt.plot(x_loc[0],y1)
#plt.plot(x_loc[0],y2)
#plt.plot(x_loc[0],y3)
#plt.plot(x_loc[0],y4)
#plt.plot(x_loc[0],y5)
#plt.plot(x_loc[0],y6)
#plt.plot(x_loc[0],y7)

#plt.show()

yhat1 = savgol_filter(y1, 75, 3)
#yhat2 = savgol_filter(y2, 75, 3)
#yhat3 = savgol_filter(y3, 75, 3)
#yhat4 = savgol_filter(y4, 75, 3)
#yhat5 = savgol_filter(y5, 75, 3)
#yhat6 = savgol_filter(y6, 75, 3)
#yhat7 = savgol_filter(y7, 75, 3)


zhst1 = exponential_smoothing(y1, 0.05)
#zhst2 = exponential_smoothing(y2, 0.05)
#zhst3 = exponential_smoothing(y3, 0.05)
#zhst4 = exponential_smoothing(y4, 0.05)
#zhst5 = exponential_smoothing(y5, 0.05)
#zhst6 = exponential_smoothing(y6, 0.05)
#zhst7 = exponential_smoothing(y7, 0.05)

#zhat = savgol_filter(z, 100, 3)


#x_values = np.array(x_loc[0])  # Convert x_loc[0] to a numpy array

#yprime = exponential_smoothing(y, 0.1)
#zprime = exponential_smoothing(z, 0.1)



## Save Data ##
#with open('studies/Goswift/results/pod_and_wedge_raw_1_F15_bodys.csv','w',newline='') as file:
#    writer = csv.writer(file)
#    writer.writerow(["x location (in), Pressure (dp/p)"])
#    writer.writerows(zip(x_loc[0],y))


#!!!!!!!!!!!!!!

plt.plot(x_loc[0],yhat1)
#plt.plot(x_loc[0],yhat2)
#plt.plot(x_loc[0],yhat3)
#plt.plot(x_loc[0],yhat4)
#plt.plot(x_loc[0],yhat5)
#plt.plot(x_loc[0],yhat6)
#plt.plot(x_loc[0],yhat7)

plt.title("sav_gol")

plt.show()

plt.plot(x_loc[0],zhst1)
#plt.plot(x_loc[0],zhst2)
#plt.plot(x_loc[0],zhst3)
#plt.plot(x_loc[0],zhst4)
#plt.plot(x_loc[0],zhst5)
#plt.plot(x_loc[0],zhst6)
#plt.plot(x_loc[0],zhst7)
plt.title("exponential_smoothing")
#plt.legend(["1", "2", "3", "4", "5", "6"])
plt.show()


comp1 = savgol_filter(zhst1, 75, 3)
#comp2 = savgol_filter(zhst2, 75, 3)
#comp3 = savgol_filter(zhst3, 75, 3)
#comp4 = savgol_filter(zhst4, 75, 3)
#comp5 = savgol_filter(zhst5, 75, 3)
#comp6 = savgol_filter(zhst6, 75, 3)
#comp7 = savgol_filter(zhst7, 75, 3)

plt.plot(x_loc[0],comp1)
#plt.plot(x_loc[0],comp2)
#plt.plot(x_loc[0],comp3)
#plt.plot(x_loc[0],comp4)
#plt.plot(x_loc[0],comp5)
#plt.plot(x_loc[0],comp6)
#plt.plot(x_loc[0],comp7)
plt.title("expontnetial_smoothing - sav_gol")
plt.show()

comp7 =     exponential_smoothing(yhat1, 0.1)
#comp8 =     exponential_smoothing(yhat2, 0.1)
#comp9 =     exponential_smoothing(yhat3, 0.1)
#comp10 =    exponential_smoothing(yhat4, 0.1)
#comp11 =    exponential_smoothing(yhat5, 0.1)
#comp12 =    exponential_smoothing(yhat6, 0.1)
#comp13 =    exponential_smoothing(yhat7, 0.1)

#plt.plot(x_loc[0],comp8)
#plt.plot(x_loc[0],comp9)
#plt.plot(x_loc[0],comp10)
#plt.plot(x_loc[0],comp11)
#plt.plot(x_loc[0],comp12)
#plt.plot(x_loc[0],comp13)
plt.plot(x_loc[0],comp7, color = "black")

#plt.xlabel("X (ft)")
#plt.ylabel("dp / P")
#plt.title("Pod and Wedge Nearfield Signature")
#plt.xlim(75,135)
#plt.xlim(230,300)
#plt.savefig("studies/Goswift/results/pod_and_wedge_deformed.png", transparent = True)

plt.title("sav_gol - exponential_smoothing (good for magnitude now)")

plt.show()
#print(ffd_box[-1])

#result = [(a - b)*393.13 for a, b in zip(comp7, comp1)]
#plt.plot(x_loc[0],result)
#plt.xlim(230,300)
#plt.show()

#results = np.subtract(comp7, comp1)
#plt.plot(x_loc[0],results)
#plt.xlim(230,300)
#plt.show()

## compare overpressure for resolution constraints
#new = []
#for i in range(len(comp8)):
#    new.append(((comp8[i]) - (comp7[i]))*393.13)
#    #print(comp7[i],",", comp1[i],",", new[i])
#plt.plot(x_loc[0],new, color = "black")
#plt.xlim(230,300)
#plt.legend(["Deformed - Undeformed"])
#plt.xlabel("X (ft)")
#plt.ylabel("dP, lbf/ft^2")
#plt.show()
##plt.plot(x_loc[0],zhat)
#
##plt.plot(x_loc[0],zhat)
##plt.xlim(20000,42500)
##plt.legend(["Raw Data", "Exponential Smoothing", "Savitzky-Golay Filter"])
#plt.legend(["smothed data"])
##plt.ylim(-.1, .15)
#plt.title("Nearfield Signature")
#plt.ylabel("dp / P")
#plt.xlabel("X (ft)")
#plt.show()
#
#
#loc = np.linspace(0,len(loudness),len(loudness))
#plt.plot(loc,loudness)
#plt.title("Loudness vs Iteration")
#plt.xlabel("Iteration")
#plt.ylabel("Loudness (db)")
#plt.show()
#
#
#x = []
#for i in range(len(ffd_box)):
#    x.append(ffd_box[i][1][0])
#
#plt.scatter(x,loudness)
#plt.title("Loudness vs location")
#plt.xlabel("Location Along Body (ft)")
#plt.ylabel("Loudness (db)")
#plt.show()
#
#
#y = []
#for i in range(len(ffd_box)):
#    y.append(ffd_box[i][0][0])
#
#plt.scatter(y,loudness)
#plt.title("Loudness vs Bump Length")
#plt.xlabel("Length")
#plt.ylabel("Loudness (db)")
#plt.show()
#
#
#z = []
#for i in range(len(ffd_box)):
#    z.append(ffd_box[i][3])
#
#plt.scatter(z,loudness)
#plt.title("Loudness vs Bump Shape")
#plt.xlabel("Delta-z")
#plt.ylabel("Loudness (db)")
#plt.show()
#
##plt.plot(x_loc[0], nearfield_sig[0])
##plt.title("Nearfield Signature")
##plt.ylabel("dP/P")
##plt.xlabel("X (ft)")
##plt.show()
##
###plt.plot(ground_sig[0][:,0],ground_sig[0][:,1])
###plt.title("Ground Signature")
###plt.ylabel("dP/P")
###plt.xlabel("X (ft)")
###plt.show()