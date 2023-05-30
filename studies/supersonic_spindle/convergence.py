import numpy as np
import matplotlib.pyplot as plt

from studies.case_running_functions import run_quad, get_order_of_convergence, write_input_file, cases, line_styles


RERUN_MACHLINE = False
study_dir = "studies/supersonic_spindle/"
plot_dir = study_dir + "plots/convergence/"


def run_quad_for_mach_and_mesh(M, grid):
    """Runs a case quad for the given Mach number and mesh density."""

    # Parameters
    half_angle = 5.0
    R_G = 287.058
    gamma = 1.4
    T_inf = 300.0
    c_inf = np.sqrt(gamma*R_G*T_inf)
    rho = 1.225
    p_inf = 1.0e5

    # Storage locations
    case_name = "M_{0}_{1}".format(M, grid)
    mesh_file = study_dir + "meshes/ehlers_spindle_{0}.vtk".format(grid)
    results_file = study_dir + "results/"+case_name+".vtk"
    report_file = study_dir + "reports/"+case_name+".json"

    # Write out input file

    # Declare MachLine input
    input_dict = {
        "flow": {
            "freestream_velocity": [1.0, 0.0, 0.0],
            "gamma" : gamma,
            "freestream_mach_number" : M
        },
        "geometry": {
            "file": mesh_file,
            "spanwise_axis" : "+y",
            "wake_model": {
                "wake_present" : False
            }
        },
        "solver": {
            "formulation": "morino"
        },
        "post_processing" : {
            "pressure_rules" : {
                "second-order" : True,
                "isentropic" : True,
                "slender-body" : True,
                "linear" : True
            }
        },
        "output" : {
            "body_file" : results_file,
            "report_file" : report_file
        }
    }

    # Dump
    input_file = study_dir + "input.json"
    write_input_file(input_dict, input_file)

    # Run quad
    reports = run_quad(input_file, run=RERUN_MACHLINE)

    # Pull out forces
    C_F = np.zeros((4,3))
    for i, report in enumerate(reports):
        C_F[i,0] = report["total_forces"]["Cx"]
        C_F[i,1] = report["total_forces"]["Cy"]
        C_F[i,2] = report["total_forces"]["Cz"]

    # Get system dimension and average characteristic length
    N_sys = reports[0]["solver_results"]["system_dimension"]
    l_avg = reports[0]["mesh_info"]["average_characteristic_length"]

    return N_sys, l_avg, C_F


if __name__=="__main__":

    # Study parameters
    grids = ["ultra_coarse", "coarse", "medium", "fine"]#, "ultra_fine"]
    Ms = [1.5, 1.7, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]

    # Loop through parameters
    N_sys = np.zeros(len(grids))
    l_avg = np.zeros(len(grids))
    C_F = np.zeros((len(grids), len(Ms), 4, 3))
    for i, grid in enumerate(grids):
        for j, M in enumerate(Ms):

            N_sys[i], l_avg[i], C_F[i,j] = run_quad_for_mach_and_mesh(M, grid)

    # Set up plotting

    # Get errors
    err = np.abs((C_F[:-1] - C_F[-1])/C_F[-1])

    # Plot C_x errors
    for j, M in enumerate(Ms):

        # Plot
        plt.figure()
        for l, case in enumerate(cases):
            plt.plot(l_avg[:-1], err[:,j,l,0], line_styles[l], label=case)

        # Format
        plt.xlabel('$l_{avg}$')
        plt.ylabel('Fractional Error in $C_x$')
        plt.xscale('log')
        plt.yscale('log')
        if j==0:
            plt.legend()
        plt.savefig(plot_dir+"err_C_x_M_{0}.pdf".format(M))
        plt.savefig(plot_dir+"err_C_x_M_{0}.svg".format(M))
        plt.close()

    #Analyze convergence of Cx
    print()
    print("Cx Convergence Rate")
    print("-------------------")
    slopes = []
    for l, case in enumerate(['ML', 'MH', 'SL', 'SH']):
        slopes.append([])
        for j, M in enumerate(Ms):
            slopes[l].append(get_order_of_convergence(l_avg, C_F[:,j,l,0], truth_from_results=True))

        print("{0}: {1} +/- {2}".format(case, round(np.average(slopes[l]), 4), round(np.std(slopes[l]), 4)))