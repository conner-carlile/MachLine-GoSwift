import numpy as np
import matplotlib.pyplot as plt
from sboomwrapper import SboomWrapper
#from pyldb import pyldb

dp_directory = "studies/GoSwift/results/off_body_pressure.csv"


# Cart3D format
data = np.genfromtxt(dp_directory , skip_header=32)
conv = 1/0.0254
sig = np.asarray([data[:,3]/conv, data[:,4]]).T

MACH = 1.5 # set Mach for current case
r_over_l = 3
ref_length = 27.432 # for X-59   ########update for n+1
altitude = 53200
PROP_R = r_over_l*ref_length*3.28084

_sboom = SboomWrapper('./temp', 'sboom.exe')
_sboom.set(mach_number=MACH,
            altitude=altitude,
            propagation_start=PROP_R,
            altitude_stop=0,
            output_format=0,
            input_xdim=2,
            azimuthal_angles = 0,
            propagation_points=40000,
            padding_points=8000)

_sboom.set(signature=sig, input_format=0) # set input signature and format
sboom_results = _sboom.run(atmosphere_input=None) # run sBOOOM
g_sig = sboom_results["signal_0"]["ground_sig"]

# np.savetxt(dp_directory + label + 'g_sig.txt', g_sig)

#noise_level = pyldb.perceivedloudness(g_sig[:, 0],
#                                      g_sig[:, 1],)
#
#
#print('PLdB (pressure): ', noise_level)