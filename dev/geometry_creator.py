import copy
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull, convex_hull_plot_2d


def _export_vtk(filename, vertices, panels):
    # Exports the vertices and panel data to a VTK file.

    # Check extension
    if '.vtk' not in filename:
        raise IOError("Filename for VTK export must contain .vtk extension.")

    # Open file
    with open(filename, 'w') as export_handle:
        
        # Write header
        print("# vtk DataFile Version 3.0", file=export_handle)
        print("MachLine geometry file. Generated by MachLine Geometry Creator, USU AeroLab (c) 2022.", file=export_handle)
        print("ASCII", file=export_handle)

        # Write dataset
        print("DATASET POLYDATA", file=export_handle)

        # Write vertices
        print("POINTS {0} float".format(len(vertices)), file=export_handle)
        for vertex in vertices:
            print("{0:<20.12}{1:<20.12}{2:<20.12}".format(*vertex), file=export_handle)

        # Determine polygon list size
        N_panels = panels.shape[0]
        size = N_panels*4

        # Write panel polygons
        print("POLYGONS {0} {1}".format(N_panels, size), file=export_handle)
        for panel in panels:
            print("3 "+" ".join([str(i) for i in panel]), file=export_handle)
            

def _get_proper_points_and_panels_for_closed_right_cone(h, r, N_transverse, N_theta, N_radial):
    # Generates the points and panels for half of a closed right cone, aligned with the x-axis and its base at the origin.
    # resultant mesh should be mirrored about the xy plane.

    # Determine number of vertices
    N_verts = 1 + N_transverse*(N_theta+1) + (N_radial-1)*(N_theta+1) + 1

    # Determine number of panels
    N_panels = (2*(N_transverse-1)+1)*N_theta + (2*(N_radial-1)+1)*N_theta

    # Initialize storage
    vertices = np.zeros((N_verts,3))
    panels = np.zeros((N_panels,3), dtype=int)

    # Apex
    vertices[0,0] = h

    # Generate distribution of locations in x, r, and theta
    X = np.linspace(h, 0.0, N_transverse+1)
    R = np.linspace(r, 0.0, N_radial+1)
    TH = np.linspace(0.0, np.pi, N_theta+1)

    # Calculate radii of side sections
    R_x = (h-X)*r/h

    # Calculate trig functions
    C_TH = np.cos(TH)
    S_TH = np.sin(TH)

    # Generate points going down the side
    for i in range(N_transverse):
        for j in range(N_theta+1):

            # Determine index
            ind = 1 + j + i*(N_theta+1)

            # x coordinate
            vertices[ind,0] = X[i+1]

            # y coordinate
            vertices[ind,1] = C_TH[j]*R_x[i+1]

            # z coordinate
            vertices[ind,2] = S_TH[j]*R_x[i+1]

    # Get end from which the bottom needs to pick up
    end_of_side = copy.copy(ind)

    # Generate points on bottom (from edge to center)
    for i in range(N_radial-1):
        for j in range(N_theta+1):

            # Determine index
            ind = end_of_side+1 + j + i*(N_theta+1)

            # y coordinate
            vertices[ind,1] = C_TH[j]*R[i+1]

            # z coordinate
            vertices[ind,2] = S_TH[j]*R[i+1]

    # Initialize panels

    # First row at tip
    for i in range(N_theta):

        # Set indices
        panels[i,0] = 0
        panels[i,1] = i+1
        panels[i,2] = i+2

    # Subsequent rows up side and on bottom
    # Each loop iteration will create 2 panels
    for i in range(N_transverse-1+N_radial-1):
        for j in range(N_theta):

            # Determine index of first panel
            ind = N_theta + 2*j + 2*i*N_theta

            # Set indices
            panels[ind,0] = 1 + j + i*(N_theta+1)
            panels[ind,1] = 1 + j + (i+1)*(N_theta+1)
            panels[ind,2] = 1 + j+1 + (i+1)*(N_theta+1)

            # Determine index of second panel
            ind += 1

            # Set indices
            panels[ind,0] = 1 + j + i*(N_theta+1)
            panels[ind,1] = 1 + j+1 + (i+1)*(N_theta+1)
            panels[ind,2] = 1 + j+1 + i*(N_theta+1)

    # Final row on bottom
    for j in range(N_theta):

        # Determine index of panel
        ind = N_panels-N_theta+j

        # Set indices
        panels[ind,0] = N_verts-1
        panels[ind,1] = N_verts-N_theta-1+j
        panels[ind,2] = N_verts-N_theta-2+j

    return vertices, panels
            

def _get_proper_points_and_panels_for_closed_right_cone(h, r, N_transverse, N_theta):
    # Generates the points and panels for half of a closed right cone, aligned with the x-axis and its base at the origin.
    # resultant mesh should be mirrored about the xy plane.

    # Determine number of vertices
    N_verts = 1 + N_transverse*(N_theta+1)

    # Determine number of panels
    N_panels = (2*(N_transverse-1)+1)*N_theta

    # Initialize storage
    vertices = np.zeros((N_verts,3))
    panels = np.zeros((N_panels,3), dtype=int)

    # Apex
    vertices[0,0] = h

    # Generate distribution of locations in x, r, and theta
    X = np.linspace(h, 0.0, N_transverse+1)
    TH = np.linspace(0.0, np.pi, N_theta+1)

    # Calculate radii of side sections
    R_x = (h-X)*r/h

    # Calculate trig functions
    C_TH = np.cos(TH)
    S_TH = np.sin(TH)

    # Generate points going down the side
    for i in range(N_transverse):
        for j in range(N_theta+1):

            # Determine index
            ind = 1 + j + i*(N_theta+1)

            # x coordinate
            vertices[ind,0] = X[i+1]

            # y coordinate
            vertices[ind,1] = C_TH[j]*R_x[i+1]

            # z coordinate
            vertices[ind,2] = S_TH[j]*R_x[i+1]

    # Initialize panels

    # First row at tip
    for i in range(N_theta):

        # Set indices
        panels[i,0] = 0
        panels[i,1] = i+1
        panels[i,2] = i+2

    # Subsequent rows up side and on bottom
    # Each loop iteration will create 2 panels
    for i in range(N_transverse-1):
        for j in range(N_theta):

            # Determine index of first panel
            ind = N_theta + 2*j + 2*i*N_theta

            # Set indices
            panels[ind,0] = 1 + j + i*(N_theta+1)
            panels[ind,1] = 1 + j + (i+1)*(N_theta+1)
            panels[ind,2] = 1 + j+1 + (i+1)*(N_theta+1)

            # Determine index of second panel
            ind += 1

            # Set indices
            panels[ind,0] = 1 + j + i*(N_theta+1)
            panels[ind,1] = 1 + j+1 + (i+1)*(N_theta+1)
            panels[ind,2] = 1 + j+1 + i*(N_theta+1)

    return vertices, panels


def generate_proper_right_cone(filename, h, r, N_transverse, N_theta, N_radial=1, close_base=True):
    """Generates a mesh of half of a right cone, aligned with the x-axis and its base at the origin.
    The resultant mesh should be mirrored about the xy plane.
    
    Parameters
    ----------
    filename : str
        Name of the file to write the mesh to. Must have '.vtk' extension.

    h : float
        Height.

    r : float
        Base radius.
    
    N_transverse : int
        Number of sections with which to discretize the sides of the cone in the transverse direction.

    N_theta : int
        Number of sections with which to discretize the sides and bottom of the cone in the angular direction.

    N_radial : int, optional
        Number of sections with which to discretize the bottom of the cone in the radial direction. Defaults to 1.

    close_base : logical, optional
        Whether to close the bottom of the cone. Defaults to True.
    """

    # Get geometry
    if close_base:
        vertices, panels = _get_proper_points_and_panels_for_closed_right_cone(h, r, N_transverse, N_theta, N_radial)
    else:
        vertices, panels = _get_proper_points_and_panels_for_closed_right_cone(h, r, N_transverse, N_theta)

    # Export
    _export_vtk(filename, vertices, panels)


def _get_random_points_on_surface_of_sphere(N, r):
    """Generates a set of points randomly distributed across the surface of a sphere.
    
    Parameters
    ----------
    N : integer
        Number of points.
        
    r : float
        Radius of sphere.
        
    Returns
    -------
    ndarray
        Array of points.
    """

    # Initialize
    points = np.zeros((N,3))

    # Get azimuth angles
    psi = np.random.random(N)*2.0*np.pi

    # Get elevation angles (weighted towards equator, so the distribution remains even)
    theta = np.arccos(1.0-2.0*np.random.random(N))

    # Get coordinates
    points[:,0] = r*np.sin(theta)*np.cos(psi)
    points[:,1] = r*np.sin(theta)*np.sin(psi)
    points[:,2] = r*np.cos(theta)

    return points


def generate_random_sphere(filename, N, r):
    """Generates a random unstructured mesh of a sphere.
    
    Parameters
    ----------
    filename : str
        Name of the file to write the mesh to. Must have '.vtk' extension.

    N : integer
        Number of vertices.
        
    r : float
        Radius of sphere.
        Array of points.
    """
    
    # Get vertices
    verts = _get_random_points_on_surface_of_sphere(N, r)

    # Create convex hull
    hull = ConvexHull(verts)

    # Make sure normal vector points inward
    fixed_simplices = []
    for simplex in hull.simplices:

        # Calculate normal and centroid
        centroid = np.sum(verts[simplex,:], axis=0) / 3.0
        normal = np.cross(verts[simplex[1],:] - verts[simplex[0],:], verts[simplex[2],:] - verts[simplex[1],:])

        # Check the normal points outward
        if np.dot(centroid, normal) > 0.0:
            fixed_simplices.append(simplex)
        else:
            fixed_simplices.append(simplex[::-1])

    # Export mesh
    _export_vtk(filename, verts, np.array(fixed_simplices))


if __name__=="__main__":

    ## Test cone
    #angles = [1.0, 5.0, 10.0, 15.0]
    #for angle in angles:
    #    h = 1.0/np.tan(np.radians(angle))
    #    generate_proper_right_cone('dev/meshes/cone_{0}_deg_coarse.vtk'.format(int(angle)), h, 1.0, 30, 25, close_base=False)

    # Random spheres
    Ns = [500, 1000, 2000, 4000]
    labels = ['coarse', 'medium', 'fine', 'ultra_fine']
    for N, label in zip(Ns, labels):
        for j in range(10):
            generate_random_sphere('dev/meshes/random_spheres/random_sphere_{0}_sample_{1}.vtk'.format(label, j), N, 1.0)