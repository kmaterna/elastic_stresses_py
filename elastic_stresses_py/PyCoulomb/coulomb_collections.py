# Definitions of objects used in this project

import numpy as np
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
                 source_object, receiver_object, receiver_results, receiver_profile):
        self.x, self.y = x, y
        self.x2d, self.y2d = x2d, y2d
        self.u_disp, self.v_disp, self.w_disp = u_disp, v_disp, w_disp
        self.zerolon, self.zerolat = zerolon, zerolat
        self.model_disp_points = model_disp_points
        self.strains = strains
        self.source_object = source_object
        self.receiver_object = receiver_object
        self.receiver_results = receiver_results  # Stress_Results object
        self.receiver_profile = receiver_profile  # Stress_Results object


class Receiver_Horiz_Profile:
    def __init__(self, depth_km, strike, dip, rake, centerlon, centerlat, lon1d, lat1d, width, length, inc, shape):
        self.depth_km = depth_km  # km, float
        self.strike = strike  # degrees, float
        self.dip = dip  # degrees, float
        self.rake = rake  # degrees, float
        self.centerlon = centerlon  # degrees, float
        self.centerlat = centerlat  # degrees, float
        self.lon1d = lon1d   # reshaped version of grid of target points, into 1D array
        self.lat1d = lat1d   # reshaped version of grid of target points, into 1D array
        self.width = width  # km, float, the north-south half-dimension of the profile
        self.length = length  # km, float, the east-west half-dimension of the profile
        self.inc = inc  # km, float
        self.shape = shape  # tuple, (ny, nx)

        if dip > 90 or dip < 0:
            raise ValueError("Error! Provided bad dip of %s (should be between 0 and 90) " % dip)
        if strike > 360:
            raise ValueError("Error! Provided bad strike of %s (should be between 0 and 360) " % strike)
        if len(self.lon1d) != self.shape[0] * self.shape[1]:
            raise ValueError("Error! Provided lon1d array of ", len(self.lon1d), " doesn't match shape ", self.shape)
        if len(self.lat1d) != self.shape[0] * self.shape[1]:
            raise ValueError("Error! Provided lat1d array of ", len(self.lat1d), " doesn't match shape ", self.shape)


class Mogi_Source:
    def __init__(self, xstart, ystart, zerolon, zerolat, depth, dV):
        self.xstart = xstart  # km
        self.ystart = ystart  # km
        self.zerolon = zerolon  # degrees
        self.zerolat = zerolat  # degrees
        self.depth = depth  # km
        self.dV = dV  # cubic meters


class Stress_Results:
    """
    A container to hold results of a bunch of stress calculations on fault patches.
    The sizes of all these arrays will match the length of the receivers
    """
    def __init__(self, eij, sigmaij, shear, dry_normal, pore_pressure, effective_normal, coulomb):
        def _as_1d_array(name, values):
            array = np.asarray(values).reshape(-1)
            if array.ndim != 1:
                raise ValueError(f"{name} must be convertible to a 1D array.")
            return array

        self.shear = _as_1d_array("shear", shear)  # in kPa
        expected_shape = self.shear.shape

        def _validate_shape(name, values):
            array = _as_1d_array(name, values)
            if array.shape != expected_shape:
                raise ValueError(
                    f"{name} must have shape {expected_shape} to match shear, got {array.shape}."
                )
            return array

        self.eij = eij  # units?
        self.sigmaij = sigmaij  # in Pa
        self.dry_normal = _validate_shape("dry_normal", dry_normal)  # in kPa
        self.pore_pressure = _validate_shape("pore_pressure", pore_pressure)  # in kPa
        self.effective_normal = _validate_shape("effective_normal", effective_normal)  # in kPa
        self.coulomb = _validate_shape("coulomb", coulomb)   # in kPa
