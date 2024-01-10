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

print(len(loudness))
## Find quietest configuration
set = 0
min = min(loudness)
print("min: ", min)
for i in range(len(loudness)):
    if loudness[i] == min:
        set = i
        print("i for lowest: ",i)

print("FFD_box: ", ffd_box[set])

plt.plot(x_loc[set], nearfield_sig[set])
plt.title("Nearfield Signature")
plt.ylabel("dP/P")
plt.xlabel("X (ft)")
plt.show()

plt.plot(ground_sig[set][:,0],ground_sig[set][:,1])
plt.title("Ground Signature")
plt.ylabel("dP/P")
plt.xlabel("X (ft)")
plt.show()


x = np.linspace(0,len(loudness),len(loudness))
plt.plot(x,loudness)
plt.title("Iteration-Loudness")
plt.xlabel("Iteration")
plt.ylabel("Loudness (db)")
plt.show()