import numpy as np
import matplotlib.pyplot as plt

def exponential_smoothing(data, alpha):
    """
    Apply exponential smoothing to the data.
    """
    smoothed_data = np.zeros_like(data)
    smoothed_data[0] = data[0]  # Initial value
    for i in range(1, len(data)):
        smoothed_data[i] = alpha * data[i] + (1 - alpha) * smoothed_data[i-1]
    return smoothed_data

# Load data from the file, skipping the first row
#data = np.loadtxt("D:/school/aeroLab/MachLineFiltering/off_body_pressure.csv", delimiter=",", skiprows=1)
#
## Extract x and y values
#x_loc = data[:, 0]
#pressure = data[:, 1]
#
## Define a range of alpha values to sweep through
#alpha_values = [.18]
#
## Initialize an empty list to store smoothed data for each alpha value
#smoothed_data_list = []
#
## Sweep through each alpha value and calculate the smoothed data
#for alpha in alpha_values:
#    smoothed_data = exponential_smoothing(pressure, alpha)
#    smoothed_data_list.append(smoothed_data)
#
## Plot original data and smoothed data for each alpha value
#plt.plot(x_loc, pressure, label='Original Data', color='black', linestyle='dotted')
#for i, alpha in enumerate(alpha_values):
#    plt.plot(x_loc, smoothed_data_list[i], label=f'Alpha = {alpha:.2f}')
#plt.title("Single Exponential Smoothing with Different Alpha Values")
#plt.xlabel("x_loc")
#plt.ylabel("Pressure")
#plt.legend()
#plt.grid(True)
#plt.show()
#
#