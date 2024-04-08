import pickle
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from exponentalSmoothing import exponential_smoothing

## Open pickle files
with open('studies/Goswift/results/temp_ffd_box.pkl', 'rb') as pickle_file:
    ffd_box = pickle.load(pickle_file)
print(ffd_box)

with open('studies/Goswift/results/temp_x_loc.pkl', 'rb') as pickle_file:
    x_loc = pickle.load(pickle_file)

with open('studies/Goswift/results/temp_nearfield.pkl', 'rb') as pickle_file:
    nearfield_sig  = pickle.load(pickle_file)

with open('studies/Goswift/results/temp_ground.pkl', 'rb') as pickle_file:
    ground_sig = pickle.load(pickle_file)

with open('studies/Goswift/results/temp_loudness.pkl', 'rb') as pickle_file:
    loudness = pickle.load(pickle_file)


## Find quietest configuration
#set = 0
#min = min(loudness)
#max = max(loudness)
#print("min: ", min)
#print("max: ", max)
#for i in range(len(loudness)):
#    if loudness[i] == min:
#        set = i
#        print("i for lowest: ", i)
#print(ffd_box)
#
#plt.plot( x_loc[set],nearfield_sig[set]) #x_loc[0:1000][set],   ??????????
#plt.title("Nearfield Signature")
#plt.ylabel("dP")
#plt.xlabel("X (in)")
#plt.show()
#
#
#for i in range(len(loudness)):
#    plt.plot(ground_sig[i][0][:,0],ground_sig[i][0][:,1]) #???????????
#plt.title("Ground Signature")
#plt.ylabel("dP/P")
#plt.xlabel("X (in)")
#plt.show()

print(ffd_box)
set = 0
min = min(loudness)
max = max(loudness)
print("min: ", min)
print("max: ", max)
for i in range(len(loudness)):
    if loudness[i] == min:
        set = i
        print("i for lowest: ", i)

#fig, ax = plt.subplots(1,2, figsize = (14,4.5))
#print(ffd_box[-1])
#plt.plot( xg,p) #x_loc[0:1000][set],
#plt.title("Nearfield Signature")
#plt.ylabel("dP/P")
#plt.xlabel("X (in)")
#plt.show()
#print(x_loc[set])
#print(len(nearfield_sig))
#plt.plot(x_loc[0],nearfield_sig[0])
#plt.plot(x_loc[0],nearfield_sig[1])
#plt.title("Nearfield Signature")
#plt.ylabel("dp / P")
#plt.xlabel("X (in)")

#!!!!!!!!!!!!
#ax[0].plot(x_loc[0],nearfield_sig[0])
#ax[0].set_title("Nearfield Signature")
#ax[0].set_ylabel("dp / P")
#ax[0].set_xlabel("X (in)")
##ax[0].set_ylim(-0.1,0.15, 0.5)
#
##plt.show()
#
#
#for i in range(len(loudness)):
##plt.plot(g_sig[0][:,0],g_sig[0][:,1])
#    ax[1].plot(ground_sig[i][0][:,0],ground_sig[i][0][:,1])
#    #plt.plot(ground_sig[i][0][:,0],ground_sig[i][0][:,1])
#    #plt.show()
#ax[1].set_title("Ground Signature")
#ax[1].set_ylabel("dp (psf)")
#ax[1].set_xlabel("time (ms)")
#
#plt.show()
#!!!!!!!!!!!!!!!


#Smooth Data
#y = nearfield_sig[0]
#yhat = savgol_filter(y, 50, 3)
#yhat[400:len(nearfield_sig[0])] = savgol_filter(y[400:len(nearfield_sig[0])], 200, 3)
        

y = nearfield_sig[0]
#z = nearfield_sig[1]

#print(y)

# Apply Savitzky-Golay filter with a 50 window to the entire data set
yhat = savgol_filter(y, 250, 3)
#zhat = savgol_filter(z, 100, 3)

x_values = np.array(x_loc[0])  # Convert x_loc[0] to a numpy array

yprime = exponential_smoothing(y, 0.1)
#zprime = exponential_smoothing(z, 0.1)

print(len(y))
# Apply Savitzky-Golay filter with a 200 window to the specified range
#yhat[950:1000] = savgol_filter(y[950:1000], 50, 1)


#with open('studies/Goswift/results/pod_and_wedge_raw_1_F15_bodys.csv','w',newline='') as file:
#    writer = csv.writer(file)
#    writer.writerow(["x location (in), Pressure (dp/p)"])
#    writer.writerows(zip(x_loc[0],y))
#
#with open('studies/Goswift/results/pod_and_wedge_sav_gol_filter_1_f_15_bodys.csv','w',newline='') as file:
#    writer = csv.writer(file)
#    writer.writerow(["x location (mm), Pressure (dp/p)"])
#    writer.writerows(zip(x_loc[0],yhat))
#
#with open('studies/Goswift/results/pod_and_wedge_exponential_filter_1_f_15_bodys.csv','w',newline='') as file:
#    writer = csv.writer(file)
#    writer.writerow(["x location (mm), Pressure (dp/p)"])
#    writer.writerows(zip(x_loc[0],yprime))

##with open('studies/Goswift/results/pressure.csv','w',newline='') as file:
##    writer = csv.writer(file)
##    writer.writerow(["Pressure (dp/p)"])
##    for row in yhat:
##        writer.writerow([row])


#for i in range(len(nearfield_sig)):
plt.plot(x_loc[0],y)
#plt.plot(x_loc[0],z)
#plt.plot(x_loc[0], yprime)
plt.plot(x_loc[0],yhat)
#plt.plot(x_loc[0],zhat)
#plt.xlim(20000,42500)
#plt.legend(["Raw Data", "Exponential Smoothing", "Savitzky-Golay Filter"])
plt.legend(["Raw Data", "Savitzky-Golay Filter"])
#plt.ylim(-.1, .15)
plt.title("Nearfield Signature")
plt.ylabel("dp / P")
plt.xlabel("X (mm)")
plt.show()
#
#
##loc = np.linspace(0,len(loudness),len(loudness))
##plt.plot(loc,loudness)
##plt.title("Loudness vs Iteration")
##plt.xlabel("Iteration")
##plt.ylabel("Loudness (db)")
##plt.show()
#
#
#x = []
#for i in range(len(ffd_box)):
#    x.append(ffd_box[i][1][0])
#
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
#plt.plot(x_loc[0], nearfield_sig[0])
#plt.title("Nearfield Signature")
#plt.ylabel("dP/P")
#plt.xlabel("X (ft)")
#plt.show()
#
##plt.plot(ground_sig[0][:,0],ground_sig[0][:,1])
##plt.title("Ground Signature")
##plt.ylabel("dP/P")
##plt.xlabel("X (ft)")
##plt.show()