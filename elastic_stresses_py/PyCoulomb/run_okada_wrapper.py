import numpy as np
from okada_wrapper import dc3dwrapper, dc3d0wrapper
from Tectonic_Utils.geodesy import fault_vector_functions
from .disp_points_object.disp_points_object import Displacement_points
from . import conversion_math, utilities


def compute_ll_def(inputs, params, disp_points):
    """
    Loop through a list of lon/lat and compute their displacements due to all sources put together.
    """
    if not disp_points:
        return []
    model_disp_points = []
    print("Number of disp_points:", len(disp_points))
    for point in disp_points:
        [xi, yi] = fault_vector_functions.latlon2xy(point.lon, point.lat, inputs.zerolon, inputs.zerolat)
        u_disp, v_disp, w_disp = compute_surface_disp_point(inputs.source_object, params.alpha, xi, yi,
                                                            compute_depth=point.depth)
        model_point = Displacement_points(lon=point.lon, lat=point.lat,
                                          dE_obs=u_disp,
                                          dN_obs=v_disp,
                                          dU_obs=w_disp,
                                          Se_obs=0, Sn_obs=0, Su_obs=0, name=point.name)
        model_disp_points.append(model_point)
    return model_disp_points


def compute_surface_disp_point(sources, alpha, x, y, compute_depth=0):
    """
    A major compute loop for each fault source object at one x/y point.
    x/y in the same coordinate system as the fault object. Computes displacement and strain tensor.

    :param sources: list of fault objects
    :param alpha: float
    :param x: float
    :param y: float
    :param compute_depth: depth of observation. Default depth is at surface of earth
    :returns: three floats
    """
    u_disp, v_disp, w_disp = 0, 0, 0

    for source in sources:
        desired_coords_grad_u, desired_coords_u = compute_strains_stresses_from_one_fault(source, x, y, compute_depth,
                                                                                          alpha)
        # Update the displacements from all sources
        u_disp = u_disp + desired_coords_u[0][0]
        v_disp = v_disp + desired_coords_u[1][0]
        w_disp = w_disp + desired_coords_u[2][0]  # vertical

    return u_disp, v_disp, w_disp


def compute_ll_strain(inputs, params, strain_points):
    """
    Loop through a list of lon/lat and compute their strains due to all sources put together.
    """
    if not strain_points:
        return []
    print("Number of strain_points:", len(strain_points))
    cartesian_strain_points = utilities.convert_ll2xy_disp_points(strain_points, inputs.zerolon, inputs.zerolat)

    strain_tensor_results = []
    # For each coordinate requested.
    for point in cartesian_strain_points:
        strain_tensor = compute_strain_point(inputs.source_object, params.alpha, point.lon, point.lat, point.depth)
        strain_tensor_results.append(strain_tensor)

    return strain_tensor_results


def compute_strain_point(sources, alpha, x, y, compute_depth=0):
    """
    A major compute loop for each fault source object at one x/y point.
    x/y in the same coordinate system as the fault object. Computes strain tensor.

    :param sources: list of fault objects
    :param alpha: float
    :param x: float
    :param y: float
    :param compute_depth: depth of observation. Default depth is at surface of earth
    :returns: 3x3 matrix
    """
    strain_tensor_total = np.zeros((3, 3))
    for source in sources:
        desired_coords_grad_u, desired_coords_u = compute_strains_stresses_from_one_fault(source, x, y, compute_depth,
                                                                                          alpha)
        # Strain tensor math
        strain_tensor = conversion_math.get_strain_tensor(desired_coords_grad_u)
        strain_tensor_total = np.add(strain_tensor, strain_tensor_total)
    return strain_tensor_total


def compute_strains_stresses_from_one_fault(source, x, y, z, alpha):
    """
    The main math of DC3D
    Operates on a source object (e.g., fault),
    and an xyz position in the same cartesian reference frame.
    """
    R = source.R
    R2 = source.R2
    strike_slip = source.rtlat * -1  # The dc3d coordinate system has left-lateral positive.

    # Compute the position relative to the translated, rotated fault.
    translated_pos = np.array(
        [[x - source.xstart], [y - source.ystart], [-z]])
    xyz = R.dot(translated_pos)
    if source.potency:
        success, u, grad_u = dc3d0wrapper(alpha, [xyz[0][0], xyz[1][0], xyz[2][0]], source.top, source.dipangle,
                                          [source.potency[0], source.potency[1], source.potency[2],
                                           source.potency[3]])
        grad_u = grad_u * 1e-9  # DC3D0 Unit correction: potency from N-m results in strain in nanostrain
        u = u * 1e-6  # Unit correction: potency from N-m results in displacements in microns.
    else:
        success, u, grad_u = dc3dwrapper(alpha, [xyz[0][0], xyz[1][0], xyz[2][0]], source.top, source.dipangle,
                                         [0, source.L], [-source.W, 0],
                                         [strike_slip, source.reverse, source.tensile])
        grad_u = grad_u * 1e-3  # DC3D Unit correction.
    # Solve for displacement gradients at certain xyz position

    # Rotate grad_u back into the unprimed coordinates.
    desired_coords_grad_u = np.dot(R2, np.dot(grad_u, R2.T))
    desired_coords_u = R2.dot(np.array([[u[0]], [u[1]], [u[2]]]))
    return desired_coords_grad_u, desired_coords_u
