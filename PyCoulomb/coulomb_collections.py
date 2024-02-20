# Definitions of objects used in this project

from . import pyc_fault_object
from .inputs_object import input_obj
from .disp_points_object import disp_points_object


# A grouping of elements used for inputs into the coulomb calculation
Input_object = input_obj.Input_object

# Displacement_points are individual disp_point elements, can be put into lists of elements.
Displacement_points = disp_points_object.Displacement_points

# PyCoulomb faults object
Faults_object = pyc_fault_object.Faults_object


# Output from the Coulomb stress calculation
class Out_object:
    def __init__(self, x, y, x2d, y2d, u_disp, v_disp, w_disp, zerolon, zerolat, model_disp_points, strains,
                 source_object, receiver_object, receiver_normal, receiver_shear, receiver_coulomb, receiver_profile):
        self.x, self.y = x, y
        self.x2d, self.y2d = x2d, y2d
        self.u_disp, self.v_disp, self.w_disp = u_disp, v_disp, w_disp
        self.zerolon, self.zerolat = zerolon, zerolat
        self.model_disp_points = model_disp_points
        self.strains = strains
        self.source_object = source_object
        self.receiver_object = receiver_object
        self.receiver_normal = receiver_normal
        self.receiver_shear = receiver_shear
        self.receiver_coulomb = receiver_coulomb
        self.receiver_profile = receiver_profile


class Receiver_Horiz_Profile:
    def __init__(self, depth_km, strike, dip, rake, centerlon, centerlat, lon1d, lat1d, width, length, inc, shape):
        self.depth_km = depth_km  # km
        self.strike = strike  # degrees
        self.dip = dip  # degrees
        self.rake = rake  # degrees
        self.centerlon = centerlon  # degrees
        self.centerlat = centerlat  # degrees
        self.lon1d = lon1d
        self.lat1d = lat1d
        self.width = width
        self.length = length
        self.inc = inc
        self.shape = shape


class Mogi_Source:
    def __init__(self, xstart, ystart, zerolon, zerolat, depth, dV):
        self.xstart = xstart  # km
        self.ystart = ystart  # km
        self.zerolon = zerolon  # degrees
        self.zerolat = zerolat  # degrees
        self.depth = depth  # km
        self.dV = dV  # cubic meters
