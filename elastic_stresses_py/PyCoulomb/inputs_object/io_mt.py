# Input a moment tensor (math not done)
import numpy as np
from ..pyc_fault_object import Faults_object
from Tectonic_Utils.geodesy import fault_vector_functions
from Tectonic_Utils.seismo import moment_calculations, MT_calculations


def get_DC_potency(rake, momentmagnitude, mu):
    """
    Given the basic double couple parameters,
    Return the four-vector used in Okada DC3D0.
    Pot1 = strike-slip moment of DC / mu
    Pot2 = dip-slip moment of DC / mu
    Pot3 = tensile = M_TENSILE / lambda
    Pot4 = inflation = M_ISO / mu
    In a more general case, we would use a different MT format to handle non-DC parts.
    Right now, it only handles DC focal mechanisms.
    Moment in newton meters
    """
    total_moment = moment_calculations.moment_from_mw(momentmagnitude)
    dc_moment = total_moment * 1.00

    strike_slip_fraction, dip_slip_fraction = fault_vector_functions.get_rtlat_dip_slip(1.0, rake)
    print("strike_slip fraction: ", strike_slip_fraction, " / 1.0")
    print("dip_slip fraction: ", dip_slip_fraction, " / 1.0")
    strike_slip_fraction = -1 * strike_slip_fraction  # DC3D0 wants left lateral slip.
    p1 = dc_moment * strike_slip_fraction / mu
    p2 = dc_moment * dip_slip_fraction / mu
    # In the double-couple case, this is zero.
    p3 = 0
    p4 = 0
    return [p1, p2, p3, p4]


def get_mag_from_dc_potency(potency, mu, rake):
    """
    The opposite of the function above.
    Potency is a four-vector in newton-meters
    Mu in Pascals
    """
    ss_moment = np.abs(potency[0]) * mu
    ds_moment = np.abs(potency[1]) * mu
    strike_slip_fraction, dip_slip_fraction = fault_vector_functions.get_rtlat_dip_slip(1.0, rake)
    if ss_moment >= ds_moment:
        dc_moment = ss_moment / abs(strike_slip_fraction)
    else:  # in case the strike slip moment is zero
        dc_moment = ds_moment / abs(dip_slip_fraction)
    Mw = moment_calculations.mw_from_moment(dc_moment)
    return Mw


def compute_params_for_MT_source(MT, rake, lon, lat, zerolon, zerolat, mu, lame1):
    """
    Given information about point sources from moment tensors,
    Return the right components that get packaged into input_obj.
    """
    [xcenter, ycenter] = fault_vector_functions.latlon2xy(lon, lat, zerolon, zerolat)
    potency = get_MT_potency(MT, rake, mu, lame1)
    comment = ''
    return [xcenter, ycenter, potency, comment]


def get_MT_potency(MT, rake, mu, lame1):
    """
    An unfinished function, since the computation is not trivial.
    Return the four-vector used in Okada DC3D0 from the full six-component moment tensor.
    Pot1 = strike-slip moment of DC / mu
    Pot2 = dip-slip moment of DC / mu
    Pot3 = tensile = M_TENSILE / lambda
    Pot4 = inflation = M_ISO / mu
    Moment in newton meters
    """
    # compute the DC moment, ISO moment, and tensile moment.
    iso, clvd, dc = MT_calculations.decompose_iso_dc_clvd(MT)
    dc_moment = dc[0][0]
    tensile_moment = 0  # THIS IS WHAT OKADA NEEDS.  WILL COMPUTE SEPARATELY WITH MINSON ET AL 2007. NOT DONE YET.
    iso_moment = 0  # SAME. NOT DONE YET.

    strike_slip_fraction, dip_slip_fraction = fault_vector_functions.get_rtlat_dip_slip(1.0, rake)
    strike_slip_fraction = -1 * strike_slip_fraction  # DC3D0 wants left lateral slip.

    p1 = dc_moment * strike_slip_fraction / mu
    p2 = dc_moment * dip_slip_fraction / mu
    p3 = tensile_moment / lame1
    p4 = iso_moment / mu
    return [p1, p2, p3, p4]


def get_MT_source(line, zerolon, zerolat, lame1, mu):
    """
    Create a source object from a six-component moment tensor solution
    """
    [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp, strike, dip, rake, lon, lat, depth] = read_moment_tensor_source_line(line)
    MT = MT_calculations.get_MT(Mrr, Mtt, Mpp, Mrt, Mrp, Mtp)
    [x, y, potency, comment] = compute_params_for_MT_source(MT, rake, lon, lat, zerolon, zerolat, mu, lame1)
    one_source_object = Faults_object(xstart=x, xfinish=x, ystart=y, yfinish=y, rtlat=0, reverse=0, tensile=0,
                                      potency=potency, strike=strike, dipangle=dip, zerolon=zerolon, zerolat=zerolat,
                                      rake=rake, top=depth, bottom=depth, comment=comment)
    return one_source_object


def read_moment_tensor_source_line(line):
    """Format: Mrr Mtt Mpp Mrt Mrp Mtp strike dip rake lon lat depth_km mu lambda"""
    [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp] = [float(i) for i in line.split()[1:7]]
    [strike, dip, rake] = [float(i) for i in line.split()[7:10]]
    [lon, lat, depth] = [float(i) for i in line.split()[10:13]]
    return [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp, strike, dip, rake, lon, lat, depth]


def compute_params_for_point_source(rake, magnitude, lon, lat, zerolon, zerolat, mu):
    """
    Helper function for setting up faults from 'Source_FM' format
    Return the right components that get packaged into input_obj.
    """
    [xcenter, ycenter] = fault_vector_functions.latlon2xy(lon, lat, zerolon, zerolat)
    potency = get_DC_potency(rake, magnitude, mu)
    comment = ''
    return [xcenter, ycenter, potency, comment]
