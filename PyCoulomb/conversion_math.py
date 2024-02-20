# Stress/strain/geometry functions


import numpy as np
from Tectonic_Utils.geodesy import fault_vector_functions


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


def get_stress_tensor(eij, lamda, mu):
    """
    Starts with strain tensor (3x3 2D array)
    Returns a stress tensor (3x3 2D array).
    lamda and mu are Lame parameters
    """
    rows, cols = np.shape(eij)
    stress_tensor = np.zeros(np.shape(eij))
    for i in range(rows):
        for j in range(cols):
            if i == j:
                stress_tensor[i][j] = lamda*(eij[0][0]+eij[1][1]+eij[2][2]) + 2*mu*eij[i][j]
            else:
                stress_tensor[i][j] = 2*mu*eij[i][j]
    return stress_tensor


def get_coulomb_stresses(tau, strike, rake, dip, friction, B):
    """
    Given a stress tensor, receiver strike, receiver rake, and receiver dip
    Resolve the stress changes on the fault plane.

    :param tau: full 3x3 stress tensor
    :param strike: float, in degrees
    :param rake: float, in degrees
    :param dip: float, in degrees
    :param friction: float, coefficient of friction
    :param B: Skepmton's coefficient
    :returns: list of 3 floats, in KPa
    """
    # First compute the geometric vectors associated with the receiver
    strike_unit_vector = fault_vector_functions.get_strike_vector(strike)  # a 3d vector in the horizontal plane.
    dip_unit_vector = fault_vector_functions.get_dip_vector(strike, dip)  # a 3d vector.
    plane_normal = fault_vector_functions.get_plane_normal(strike, dip)  # a 3d vector.

    effective_normal_stress, shear_stress, coulomb_stress = \
        get_coulomb_stresses_internal(tau, strike_unit_vector, rake, dip_unit_vector, plane_normal, friction, B)

    return effective_normal_stress, shear_stress, coulomb_stress


def get_coulomb_stresses_internal(tau, rec_strike_vector, rake, rec_dip_vector, rec_plane_normal, friction, B):
    """
    The math behind Coulomb stresses: given a stress tensor, receiver strike, receiver rake, and receiver dip.
    Resolve the stress changes on the fault plane.

    :param tau: full 3x3 stress tensor
    :param rec_strike_vector: 1d array
    :param rake: float
    :param rec_dip_vector: 1d array
    :param rec_plane_normal: 1d array
    :param friction: float, coefficient of friction
    :param B: float, Skepmton's coefficient
    :returns: list of 3 floats, Return in KPa
    """
    traction_vector = np.dot(tau, rec_plane_normal)

    # The stress that's normal to the receiver fault plane:
    dry_normal_stress = np.dot(rec_plane_normal, traction_vector)  # positive = unclamping (same as Coulomb software)
    effective_normal_stress = dry_normal_stress - (np.trace(tau) / 3.0) * B

    # The shear stress causing strike slip (in the receiver fault plane).
    shear_rtlat = np.dot(rec_strike_vector, traction_vector)

    # The shear stress causing reverse slip (in the receiver fault plane).
    shear_reverse = np.dot(rec_dip_vector, traction_vector)

    # The shear that we want (in the rake direction).
    rake_rad = np.deg2rad(rake)
    R = np.array([[np.cos(rake_rad), -np.sin(rake_rad)], [np.sin(rake_rad), np.cos(rake_rad)]])
    shear_vector = [shear_rtlat, shear_reverse]
    rotated_shear = np.dot(R, shear_vector)
    shear_in_rake_dir = rotated_shear[0]

    # Finally, do unit conversion
    effective_normal_stress = effective_normal_stress/1000.0  # convert to KPa
    shear_stress = shear_in_rake_dir/1000.0

    # The Coulomb Failure Hypothesis
    coulomb_stress = shear_stress + (friction*effective_normal_stress)   # the sign here is important.

    return effective_normal_stress, shear_stress, coulomb_stress


# ----------------------------
# GEOMETRY FUNCTIONS
# ----------------------------

def get_R_from_strike(strike):
    """Compute the rotation matrix into a system with a given fault strike"""
    # Preparing to rotate to a fault-oriented coordinate system.
    theta = strike - 90
    theta = np.deg2rad(theta)
    R = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0],
                  [0, 0, 1]])  # horizontal rotation into strike-aligned coordinates.
    R2 = np.array([[np.cos(-theta), -np.sin(-theta), 0], [np.sin(-theta), np.cos(-theta), 0], [0, 0, 1]])
    return R, R2


def rotate_points(x, y, degrees):
    """Rotate cartesian points into a new orthogonal coordinate system. Implements a rotation matrix. """
    rot_matrix = np.array([[np.cos(np.deg2rad(degrees)), -np.sin(np.deg2rad(degrees))],
                           [np.sin(np.deg2rad(degrees)), np.cos(np.deg2rad(degrees))]])
    unprimed_vector = np.array([[x], [y]])
    xprime, yprime = np.dot(rot_matrix, unprimed_vector)
    return xprime, yprime


def rotate_list_of_points(xlist, ylist, degrees):
    """Calls rotation matrix on a 1d list of points"""
    xprime_list, yprime_list = [], []
    for x, y in zip(xlist, ylist):
        xprime, yprime = rotate_points(x, y, degrees)
        xprime_list.append(xprime)
        yprime_list.append(yprime)
    return xprime_list, yprime_list


def get_geom_attributes_from_receiver_profile(profile):
    """
    Pre-compute geometry for receiver plane
    """
    strike_unit_vector = fault_vector_functions.get_strike_vector(profile.strike)  # 3d vector in horizontal plane.
    dip_unit_vector = fault_vector_functions.get_dip_vector(profile.strike, profile.dip)  # a 3d vector.
    plane_normal = fault_vector_functions.get_plane_normal(profile.strike, profile.dip)  # a 3d vector.
    return strike_unit_vector, dip_unit_vector, plane_normal
