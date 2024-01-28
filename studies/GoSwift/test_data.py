import pickle
import numpy as np
import matplotlib.pyplot as plt

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


## Find quietest configuration
set = 0
min = min(loudness)
max = max(loudness)
print("min: ", min)
print("max: ", max)
for i in range(len(loudness)):
    if loudness[i] == min:
        set = i
        print("i for lowest: ", i)



plt.plot( x_loc[set],nearfield_sig[set]) #x_loc[0:1000][set],
plt.title("Nearfield Signature")
plt.ylabel("dP/P")
plt.xlabel("X (ft)")
plt.show()

#for i in range(len(loudness)):
plt.plot(ground_sig[set][:,0]/12,ground_sig[set][:,1])
plt.title("Ground Signature")
plt.ylabel("dP/P")
plt.xlabel("X (ft)")
plt.show()


loc = np.linspace(0,len(loudness),len(loudness))
plt.plot(loc,loudness)
plt.title("Loudness vs Iteration")
plt.xlabel("Iteration")
plt.ylabel("Loudness (db)")
plt.show()


x = []
for i in range(len(ffd_box)):
    x.append(ffd_box[i][1][0])


plt.scatter(x,loudness)
plt.title("Loudness vs location")
plt.xlabel("Location Along Body (ft)")
plt.ylabel("Loudness (db)")
plt.show()


y = []
for i in range(len(ffd_box)):
    y.append(ffd_box[i][0][0])

plt.scatter(y,loudness)
plt.title("Loudness vs Bump Length")
plt.xlabel("Length")
plt.ylabel("Loudness (db)")
plt.show()


z = []
for i in range(len(ffd_box)):
    z.append(ffd_box[i][3])

plt.scatter(z,loudness)
plt.title("Loudness vs Bump Shape")
plt.xlabel("Delta-z")
plt.ylabel("Loudness (db)")
plt.show()

#plt.plot(x_loc[0], nearfield_sig[0])
#plt.title("Nearfield Signature")
#plt.ylabel("dP/P")
#plt.xlabel("X (ft)")
#plt.show()
#
#plt.plot(ground_sig[0][:,0],ground_sig[0][:,1])
#plt.title("Ground Signature")
#plt.ylabel("dP/P")
#plt.xlabel("X (ft)")
#plt.show()