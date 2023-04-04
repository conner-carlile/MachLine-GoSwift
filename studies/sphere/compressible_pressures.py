import numpy as np
import matplotlib.pyplot as plt

from studies.case_running_functions import write_input_file, run_quad
from studies.paraview_functions import extract_all_data, get_data_column_from_array


RERUN_MACHLINE = True
study_dir = "studies/sphere/"
plot_dir = study_dir + "plots/"


def extract_pressures(result_file):

    # Get data
    headers, data = extract_all_data(result_file, which_data='cell')

    # Get locations and pressures
    C_P = get_data_column_from_array(headers, data, 'C_p_inc')
    locs = np.zeros((len(C_P),3))
    locs[:,0] = get_data_column_from_array(headers, data, 'centroid:0') - 1.0
    locs[:,1] = get_data_column_from_array(headers, data, 'centroid:1')
    locs[:,2] = get_data_column_from_array(headers, data, 'centroid:2')

    return locs, C_P


if __name__=="__main__":

    # Options
    densities = ['ultra_coarse', 'very_coarse', 'coarse', 'medium']
    phis = np.radians([30.0, 45.0, 60.0])
    thetas = np.radians([30.0, 45.0, 60.0])
    Ms = [0.4, 0.5]

    # Loop
    for M in Ms:
        for density in densities:

            for i, phi in enumerate(phis):
                for j, theta in enumerate(thetas):

                    # Initialize input
                    result_file = study_dir + "results/sphere_M_{0}_{1}_{2}_{3}.vtk".format(M, round(np.degrees(phi)), round(np.degrees(theta)), density)
                    report_file = study_dir + "reports/sphere_M_{0}_{1}_{2}_{3}.json".format(M, round(np.degrees(phi)), round(np.degrees(theta)), density)

                    V_inf = [np.cos(phi)*np.cos(theta), np.sin(phi)*np.cos(theta), np.sin(theta)]

                    # Lower-order results
                    input_dict ={
                        "flow" : {
                            "freestream_velocity" : V_inf,
                            "gamma" : 1.4,
                            "freestream_mach_number" : M
                        },
                        "geometry" : {
                            "file" : study_dir + "meshes/sphere_{0}.stl".format(density),
                            "spanwise_axis" : "+y"
                        },
                        "solver" : {
                        },
                        "post_processing" : {
                        },
                        "output" : {
                            "body_file" : result_file,
                            "report_file" : report_file
                        }
                    }

                    # Run
                    input_filename = study_dir + "input.json"
                    write_input_file(input_dict, input_filename)
                    reports = run_quad(input_filename, run=RERUN_MACHLINE)

                    ## Loop through cases
                    #for report, case in zip(reports, ['ML', 'MH', 'SL', 'SH']):

                    #    # Get pressures
                    #    result_file = report["input"]["output"]["body_file"]
                    #    locs, C_P = extract_pressures(result_file)

                    #    # Figure out thetas
                    #    theta_space = np.arccos(np.einsum('ij,j->i', locs, V_inf))
                    #    C_P_anl = 1.0 - 2.25*np.sin(theta_space)**2
                    #    err = abs((C_P - C_P_anl)/C_P_anl)
                    #    plt.figure()
                    #    plt.plot(np.degrees(theta_space), err, 'k.', markersize=1)
                    #    plt.yscale('log')
                    #    plt.xlabel('$\\theta [^\circ]$')
                    #    plt.ylabel('$|(C_P - C_{P_{anl}})/C_{P_{anl}}|$')
                    #    plt.savefig(plot_dir + "compressible_C_P_err_M_{0}_{1}_{2}_{3}_{4}.pdf".format(M, case, round(np.degrees(phi)), round(np.degrees(theta)), density))
                    #    plt.close()