#!/usr/bin/env python

"""
Can the same fault be used in Okada_wrapper and CUTDE triangles?  This script should produce the same outputs
for strain tensor and displacement vector from the same inputs.
"""

# Standard imports
import numpy as np
from okada_wrapper import dc3dwrapper, dc3d0wrapper
import cutde.halfspace as HS
# My own library imports (working to remove for purposes of this example script)
from Tectonic_Utils.geodesy import fault_vector_functions
from Elastic_stresses_py.PyCoulomb import fault_slip_object as fso
from Elastic_stresses_py.PyCoulomb.fault_slip_triangle import fault_slip_triangle as fst


# Global settings for this example
zerolon, zerolat = -115.5, 32.5  # center of cartesian coordinate system
target_lon, target_lat = -115.67, 32.68  # query point of interest
target_depth = 0  # in km
mu = 3e10  # shear modulus
lame1 = 3e10  # first lame parameter

# Mathematical functions that frequently get used in manipulating faults and displacements
def get_poissons_ratio_and_alpha(mu, lame1):
    """
    Return the poisson's ratio from a given mu (shear modulus) and lame1 (lame's first parameter)
    """
    poissons_ratio = lame1 / (2 * (lame1 + mu))
    alpha = (lame1 + mu) / (lame1 + 2 * mu)
    return [poissons_ratio, alpha]

def get_strain_tensor(dUidUj):
    """
    Starts with displacement gradient tensor (3x3 2D array)
    Returns a strain tensor (3x3 2D array).
    """
    rows, cols = np.shape(dUidUj)
    strain_tensor = np.zeros(np.shape(dUidUj))
    for i in range(rows):
        for j in range(cols):
            strain_tensor[i][j] = 0.5 * (dUidUj[i][j] + dUidUj[j][i])
    return strain_tensor

def get_R_from_strike(strike):
    """Compute the rotation matrix into a system with a given fault strike"""
    # Preparing to rotate to a fault-oriented coordinate system.
    theta = strike - 90
    theta = np.deg2rad(theta)
    R = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0],
                  [0, 0, 1]])  # horizontal rotation into strike-aligned coordinates.
    R2 = np.array([[np.cos(-theta), -np.sin(-theta), 0], [np.sin(-theta), np.cos(-theta), 0], [0, 0, 1]])
    return R, R2


# Rectangular fault displacements and strains
def compute_strains_disps_rectangle(source):
    """
    Pseudocode: From a rectangular source and a coordinate system, derive Okada displacements and strain
    """

    print("Calculating displacements by rectangular faults.")
    [xi, yi] = fault_vector_functions.latlon2xy(target_lon, target_lat, zerolon, zerolat)  # coordinate conversion

    z = target_depth
    R, R2 = get_R_from_strike(source.strike)  # rotation matrix from the fault strike
    strike_slip = source.rtlat * -1;  # The dc3d coordinate system has left-lateral positive.
    _poisson_ratio, alpha = get_poissons_ratio_and_alpha(mu, lame1)  # derived from mu and lame1

    # Compute the position relative to the translated, rotated fault.
    translated_pos = np.array([[xi - source.xstart], [yi - source.ystart], [-z]]);   # positions are in km
    xyz = R.dot(translated_pos);
    if source.potency:
        success, u, grad_u = dc3d0wrapper(alpha, [xyz[0], xyz[1], xyz[2]], source.top, source.dipangle,
                                          [source.potency[0], source.potency[1], source.potency[2],
                                           source.potency[3]]);
        grad_u = grad_u * 1e-9;  # DC3D0 Unit correction: potency from N-m results in strain in nanostrain
        u = u * 1e-6;  # Unit correction: potency from N-m results in displacements in microns.
    else:
        success, u, grad_u = dc3dwrapper(alpha, [xyz[0], xyz[1], xyz[2]], source.top, source.dipangle,
                                         [0, source.L], [-source.W, 0],
                                         [strike_slip, source.reverse, source.tensile]);
        grad_u = grad_u * 1e-3;  # DC3D Unit correction. Solve for displacement gradients at certain xyz position

    # Rotate grad_u back into the unprimed coordinates.
    desired_coords_grad_u = np.dot(R2, np.dot(grad_u, R2.T));
    desired_coords_u = R2.dot(np.array([[u[0]], [u[1]], [u[2]]]));
    strain_tensor = get_strain_tensor(desired_coords_grad_u)
    return strain_tensor, desired_coords_u;


def convert_pycoulomb_rectangle_into_two_triangles(source, startlon, startlat):
    """
    Convert one rectangular pycoulomb_fault into two triangular faults. The fault normals are expected to point up.
    """
    [x_all, y_all, _, _] = source.get_fault_four_corners()  # In cartesian position
    top_depth, bottom_depth = source.top, source.bottom
    vertex1 = np.array([x_all[0] * 1000, y_all[0] * 1000, top_depth * 1000])  # in meters
    vertex2 = np.array([x_all[1] * 1000, y_all[1] * 1000, top_depth * 1000])
    vertex3 = np.array([x_all[2] * 1000, y_all[2] * 1000, bottom_depth * 1000])
    vertex4 = np.array([x_all[3] * 1000, y_all[3] * 1000, bottom_depth * 1000])
    first_triangle = fst.TriangleFault(lon=startlon, lat=startlat, segment=source.segment,
                                       tensile=source.tensile, vertex1=vertex1, vertex2=vertex3, vertex3=vertex2,
                                       dip_slip=source.reverse, rtlat_slip=source.rtlat, depth=vertex1[2] / 1000)
    second_triangle = fst.TriangleFault(lon=startlon, lat=startlat, segment=source.segment,
                                        tensile=source.tensile, vertex1=vertex1, vertex2=vertex4, vertex3=vertex3,
                                        dip_slip=source.reverse, rtlat_slip=source.rtlat, depth=vertex1[2] / 1000)
    list_of_two_triangles = [first_triangle, second_triangle]
    return list_of_two_triangles


# CUTDE triangular fault displacements and strains
def compute_strains_disps_triangles(source):
    """
    Compute strains and displacements in cutde
    """
    print("Computing displacements by triangles.")
    fault_triangles = convert_pycoulomb_rectangle_into_two_triangles(source, zerolon, zerolat)
    poisson_ratio, _alpha = get_poissons_ratio_and_alpha(mu, lame1)  # derived from mu and lame1

    [x, y] = fault_vector_functions.latlon2xy(target_lon, target_lat, zerolon, zerolat)  # coordinate conversion

    obsx = [x * 1000]  # calculation works in meters, and uses arrays of position
    obsy = [y * 1000]  # calculation works in meters, and uses arrays of position
    obsz = [target_depth * -1000]   # in meters, negative is down
    pts = np.vstack([obsx, obsy, obsz]).T  # shape: (Npts, 3)

    slip_array = np.array(
        [[-src.rtlat_slip, src.dip_slip, src.tensile] for src in fault_triangles])  # shape:(Ntris, 3)
    fault_pts, fault_tris = fst.extract_mesh_vertices(fault_triangles)
    fault_pts = fst.flip_depth_sign(fault_pts)  # fault_pts shape: N_vertices, 3
    src_tris = fault_pts[fault_tris]  # src_tris shape: (Ntris, 3, 3)
    disp_mat = HS.disp_matrix(obs_pts=pts, tris=src_tris, nu=poisson_ratio)  # disp_mat shape: (Npts, 3, Ntris, 3)
    disp = disp_mat.reshape((-1, np.size(slip_array))).dot(
        slip_array.flatten())  # reshape by len of total slip vector
    disp_grid = disp.reshape((*np.array(obsx).shape, 3))  # disp_grid shape: Npts, 3

    # Get strain
    strain_mat = HS.strain_matrix(obs_pts=pts, tris=src_tris, nu=poisson_ratio)  # strain_mat shape: (Npts, 6, Ntris, 3)
    strain = strain_mat.reshape((-1, np.size(slip_array))).dot(
        slip_array.flatten())  # reshape by len total slip vector
    strain_tensors = strain.reshape((*np.array(obsx).shape, 6))  # strain_tensors shape: Npts, 6
    # strain[:,0] is the xx component of strain, 1 is yy, 2 is zz, 3 is xy, 4 is xz, and 5 is yz.

    modeled_strain_tensors = []
    for i, item in enumerate(strain_tensors):
        comps = strain_tensors[i]
        new_strain_tensor = np.array([[comps[0], comps[3], comps[4]],
                                      [comps[3], comps[1], comps[5]],
                                      [comps[4], comps[5], comps[2]]])
        modeled_strain_tensors.append(new_strain_tensor)

    return modeled_strain_tensors, disp_grid;


if __name__ == "__main__":
    # Generate a source fault in a convenient input format for fault objects
    test_rect = fso.fault_slip_object.FaultSlipObject(lon=zerolon + 0.01, lat=zerolat - 0.2, strike=50, dip=70,
                                                      depth=3, segment=0, length=15, width=5, rake=170, slip=-0.2)

    # Convert to a particular internal format
    pycoulomb_source = fso.fault_slip_object.fault_object_to_coulomb_fault([test_rect], zerolon_system=zerolon,
                                                                           zerolat_system=zerolat)

    strain_tensor, u = compute_strains_disps_rectangle(pycoulomb_source)  # compute Okada_wrapper displacements
    print(strain_tensor, u)

    strain_tensor, u = compute_strains_disps_triangles(pycoulomb_source)  # compute CUTDE displacements
    print(strain_tensor, u)
