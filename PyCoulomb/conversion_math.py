# Stress and strain functions
# Fault object functions


import numpy as np
from Tectonic_Utils.seismo import moment_calculations
from Tectonic_Utils.geodesy import fault_vector_functions


def get_strain_tensor(dUidUj):
    """
    Starts with displacement gradient tensor (3x3 2D array);
    Returns a strain tensor (3x3 2D array).
    """
    rows, cols = np.shape(dUidUj);
    strain_tensor = np.zeros(np.shape(dUidUj));
    for i in range(rows):
        for j in range(cols):
            strain_tensor[i][j] = 0.5 * (dUidUj[i][j] + dUidUj[j][i]);
    return strain_tensor;

def get_stress_tensor(eij, lamda, mu):
    """
    Starts with strain tensor (3x3 2D array);
    Returns a stress tensor (3x3 2D array).
    lamda and mu are Lame parameters
    """
    rows, cols = np.shape(eij);
    stress_tensor = np.zeros(np.shape(eij));
    for i in range(rows):
        for j in range(cols):
            if i == j:
                stress_tensor[i][j] = lamda*(eij[0][0]+eij[1][1]+eij[2][2]) + 2*mu*eij[i][j];
            else:
                stress_tensor[i][j] = 2*mu*eij[i][j];
    return stress_tensor;

def get_coulomb_stresses(tau, strike, rake, dip, friction, B):
    """
    Given a stress tensor, strike, rake, and dip
    Resolve the stress changes on the fault plane.
    Tau is full 3x3 stress tensor
    Friction is coefficient of friction
    B is Skepmton's coefficient
    Return in KPa
    """

    strike_unit_vector = fault_vector_functions.get_strike_vector(strike);  # a 3d vector in the horizontal plane.
    dip_unit_vector = fault_vector_functions.get_dip_vector(strike, dip);  # a 3d vector.
    plane_normal = fault_vector_functions.get_plane_normal(strike, dip);  # a 3d vector.
    traction_vector = np.dot(tau, plane_normal);

    # The stress that's normal to the receiver fault plane:
    dry_normal_stress = np.dot(plane_normal, traction_vector);  # positive = unclamping (same as Coulomb software)
    effective_normal_stress = dry_normal_stress - np.average(np.diag(tau)) * B;

    # The shear stress causing strike slip (in the receiver fault plane).
    shear_rtlat = np.dot(strike_unit_vector, traction_vector);

    # The shear stress causing reverse slip (in the receiver fault plane).
    shear_reverse = np.dot(dip_unit_vector, traction_vector);

    # The shear that we want (in the rake direction).
    rake_rad = np.deg2rad(rake);
    R = np.array([[np.cos(rake_rad), -np.sin(rake_rad)], [np.sin(rake_rad), np.cos(rake_rad)]]);
    shear_vector = [-shear_rtlat, shear_reverse];  # minus sign here defines right lateral as positive.
    rotated_shear = np.dot(R, shear_vector);
    shear_in_rake_dir = rotated_shear[0];

    # Finally, do unit conversion
    effective_normal_stress = effective_normal_stress/1000.0;  # convert to KPa
    shear_stress = shear_in_rake_dir/1000.0;

    # The Coulomb Failure Hypothesis
    coulomb_stress = shear_stress + (friction*effective_normal_stress);   # the sign here is important.

    return effective_normal_stress, shear_stress, coulomb_stress;


# ----------------------------
# FAULT OBJECT FUNCTIONS
# ----------------------------

def get_fault_center(fault_object):
    """Compute the x-y-z coordinates of the center of a fault patch.
    Index is the i'th fault patch in this fault_object"""
    W = fault_vector_functions.get_downdip_width(fault_object.top, fault_object.bottom, fault_object.dipangle);
    center_z = (fault_object.top+fault_object.bottom)/2.0;
    updip_center_x = (fault_object.xstart+fault_object.xfinish)/2.0;
    updip_center_y = (fault_object.ystart+fault_object.yfinish)/2.0;
    vector_mag = W*np.cos(np.deg2rad(fault_object.dipangle))/2.0;  # how far the middle is displaced
    # downdip from map-view
    center_point = fault_vector_functions.add_vector_to_point(updip_center_x, updip_center_y, vector_mag,
                                                              fault_object.strike+90);
    # strike+90 = downdip direction.
    center = [center_point[0], center_point[1], center_z];
    return center;

def get_fault_four_corners(fault_object):
    """
    Get the four corners of the object, including updip and downdip.
    depth is fault_object.top
    dip is fault_object.dipangle (in case you need it)
    """
    W = fault_vector_functions.get_downdip_width(fault_object.top, fault_object.bottom, fault_object.dipangle);
    strike = fault_object.strike;

    updip_point0 = [fault_object.xstart, fault_object.ystart];
    updip_point1 = [fault_object.xfinish, fault_object.yfinish];
    vector_mag = W*np.cos(np.deg2rad(fault_object.dipangle));  # how far the bottom edge is displaced
    # downdip from map-view
    downdip_point0 = fault_vector_functions.add_vector_to_point(fault_object.xstart, fault_object.ystart, vector_mag,
                                                                strike+90);
    # strike+90 = downdip direction.
    downdip_point1 = fault_vector_functions.add_vector_to_point(fault_object.xfinish, fault_object.yfinish, vector_mag,
                                                                strike+90);

    x_total = [updip_point0[0], updip_point1[0], downdip_point1[0], downdip_point0[0], updip_point0[0]];
    y_total = [updip_point0[1], updip_point1[1], downdip_point1[1], downdip_point0[1], updip_point0[1]];
    x_updip = [updip_point0[0], updip_point1[0]];
    y_updip = [updip_point0[1], updip_point1[1]];
    return [x_total, y_total, x_updip, y_updip];

def get_fault_slip_moment(fault_object, mu):
    """
    From a source fault object, calculate the seismic moment.
    Must be a finite fault, not a point source.
    Not really used yet, but could be useful in the future.
    """
    if fault_object.potency:  # for the case of point source, we can't do the moment calculation
        return None, None;
    W = fault_vector_functions.get_downdip_width(fault_object.top, fault_object.bottom, fault_object.dipangle);
    L = fault_vector_functions.get_strike_length(fault_object.xstart, fault_object.xfinish, fault_object.ystart,
                                                 fault_object.yfinish);
    area = L * W * 1000 * 1000;
    slip = fault_vector_functions.get_vector_magnitude([fault_object.rtlat, fault_object.reverse]);
    seismic_moment = mu * area * slip;
    moment_magnitude = moment_calculations.mw_from_moment(seismic_moment);
    return seismic_moment, moment_magnitude;
