# Definitions of objects used in this project

import collections
from . import conversion_math
from .inputs_object import input_obj
from Tectonic_Utils.geodesy import fault_vector_functions
from .disp_points_object import disp_points_object


Input_object = input_obj.Input_object

Out_object = collections.namedtuple('Out_object', [
    'x', 'y',
    'x2d', 'y2d', 'u_disp', 'v_disp', 'w_disp',
    'zerolon', 'zerolat',
    'model_disp_points', 'strains',
    'source_object', 'receiver_object',
    'receiver_normal', 'receiver_shear', 'receiver_coulomb', 'receiver_profile']);


class Receiver_Horiz_Profile:
    def __init__(self, depth_km, strike, dip, rake, centerlon, centerlat, lon1d, lat1d, width, length, inc, shape):
        self.depth_km = depth_km;  # km
        self.strike = strike;  # degrees
        self.dip = dip;  # degrees
        self.rake = rake;  # degrees
        self.centerlon = centerlon;  # degrees
        self.centerlat = centerlat;  # degrees
        self.lon1d = lon1d;
        self.lat1d = lat1d;
        self.width = width;
        self.length = length;
        self.inc = inc;
        self.shape = shape;


class Mogi_Source:
    def __init__(self, xstart, ystart, zerolon, zerolat, depth, dV):
        self.xstart = xstart;  # km
        self.ystart = ystart;  # km
        self.zerolon = zerolon;  # degrees
        self.zerolat = zerolat;  # degrees
        self.depth = depth;  # km
        self.dV = dV;  # cubic meters


# Displacement_points are individual disp_point elements, can be put into lists of elements.
Displacement_points = disp_points_object.Displacement_points;


Faults_object = collections.namedtuple('Faults_object', [
    'xstart', 'xfinish',
    'ystart', 'yfinish',
    'Kode', 'zerolon', 'zerolat',
    'rtlat', 'reverse', 'tensile', 'potency',
    'strike', 'dipangle', 'rake',
    'top', 'bottom', 'comment',
    'R', 'R2', 'W', 'L', 'strike_unit_vector', 'dip_unit_vector', 'plane_normal']);
"""
PyCoulomb faults object. 
rtlat, reverse, tensile are in units of meters.
Having both "rake" and "rtlat/reverse" in this object seems a little redundant, but it's used for receiver rake.
See constructor construct_pycoulomb_fault() for the creating of this object.  
We usually don't call this object directly.
"""

def construct_pycoulomb_fault(xstart, xfinish, ystart, yfinish, rtlat, reverse, tensile, potency, strike, dipangle,
                              rake, zerolon, zerolat, top, bottom, Kode=100, comment=None):
    """
    Like a constructor for this named tuple.
    The reason we have BOTH rake and rtlat/reverse is that rake can be specified for receivers, which have slip=0.
    You pass the required input parameters, and the remaining parameters are filled in.
    While keeping the interface simple,
    we are pre-computing some geometry parameters ONCE for each fault, for performance reasons.
    """
    if bottom < top:
        raise ValueError("Error! Provided bad fault- top depth (%f km) below bottom depth (%f km)" % (top, bottom) )
    if dipangle > 90 or dipangle < 0:
        raise ValueError("Error! Provided bad dip of %s (should be between 0 and 90) " % dipangle);
    if strike > 360:
        raise ValueError("Error! Provided bad strike of %s (should be between 0 and 360) " % strike);
    R, R2 = conversion_math.get_R_from_strike(strike);
    L = fault_vector_functions.get_strike_length(xstart, xfinish, ystart, yfinish);
    W = fault_vector_functions.get_downdip_width(top, bottom, dipangle);
    strike_unit_vector = fault_vector_functions.get_strike_vector(strike);  # 3d vector in horizontal plane.
    dip_unit_vector = fault_vector_functions.get_dip_vector(strike, dipangle);  # a 3d vector.
    plane_normal = fault_vector_functions.get_plane_normal(strike, dipangle);  # a 3d vector.

    one_source = Faults_object(xstart=xstart, xfinish=xfinish, ystart=ystart, yfinish=yfinish,
                               Kode=Kode, rtlat=rtlat, reverse=reverse, tensile=tensile,
                               potency=potency, strike=strike, zerolon=zerolon, zerolat=zerolat,
                               dipangle=dipangle, rake=rake, top=top, bottom=bottom, comment=comment,
                               R=R, R2=R2, W=W, L=L,
                               strike_unit_vector=strike_unit_vector, dip_unit_vector=dip_unit_vector,
                               plane_normal=plane_normal);
    return one_source;


def modify_fault_object(default_fault, xstart=None, xfinish=None, ystart=None, yfinish=None, rtlat=None,
                        reverse=None, tensile=None, potency=None, strike=None, dipangle=None, rake=None, zerolon=None,
                        zerolat=None, top=None, bottom=None):
    """
    Modify the fields in a Pycoulomb.Faults_object namedtuple.
    """
    xstart = default_fault.xstart if xstart is None else xstart;
    xfinish = default_fault.xfinish if xfinish is None else xfinish;
    ystart = default_fault.ystart if ystart is None else ystart;
    yfinish = default_fault.yfinish if yfinish is None else yfinish;
    rtlat = default_fault.rtlat if rtlat is None else rtlat;
    reverse = default_fault.reverse if reverse is None else reverse;
    tensile = default_fault.tensile if tensile is None else tensile;
    potency = default_fault.potency if potency is None else potency;
    strike = default_fault.strike if strike is None else strike;
    dipangle = default_fault.dipangle if dipangle is None else dipangle;
    rake = default_fault.rake if rake is None else rake;
    zerolon = default_fault.zerolon if zerolon is None else zerolon;
    zerolat = default_fault.zerolat if zerolat is None else zerolat;
    top = default_fault.top if top is None else top;
    bottom = default_fault.bottom if bottom is None else bottom;
    new_fault = construct_pycoulomb_fault(xstart=xstart, xfinish=xfinish, ystart=ystart, yfinish=yfinish, rtlat=rtlat,
                                          reverse=reverse, tensile=tensile, potency=potency, strike=strike,
                                          dipangle=dipangle, rake=rake, zerolon=zerolon, zerolat=zerolat, top=top,
                                          bottom=bottom);
    return new_fault;
