from . import conversion_math
from Tectonic_Utils.geodesy import fault_vector_functions
from Tectonic_Utils.seismo import moment_calculations
import numpy as np


# PyCoulomb faults object.
class Faults_object:
    def __init__(self, xstart, xfinish, ystart, yfinish, zerolon, zerolat, strike, dipangle, rake, top, bottom,
                 rtlat=0, reverse=0, tensile=0, potency=(), comment=None, Kode=100):
        """
        We include BOTH rake and rtlat/reverse because rake can be specified for receivers, which have slip=0.
        Pass the required input parameters, and a few geometry parameters are filled in for performance reasons.
        """
        self.xstart, self.xfinish = xstart, xfinish;  # km
        self.ystart, self.yfinish = ystart, yfinish;  # km
        self.zerolon, self.zerolat = zerolon, zerolat;  # degrees
        self.rtlat = rtlat;  # rtlat, reverse, tensile are in units of meters.
        self.reverse = reverse;
        self.tensile = tensile;
        self.potency = potency;
        self.strike = strike;
        self.dipangle = dipangle;
        self.rake = rake;
        self.top = top;
        self.bottom = bottom;
        self.comment = comment;
        self.Kode = Kode;
        if bottom < top:
            raise ValueError("Error! Provided bad fault- top depth (%f km) below bottom depth (%f km)" % (top, bottom) )
        if dipangle > 90 or dipangle < 0:
            raise ValueError("Error! Provided bad dip of %s (should be between 0 and 90) " % dipangle);
        if strike > 360:
            raise ValueError("Error! Provided bad strike of %s (should be between 0 and 360) " % strike);
        self.R, self.R2 = conversion_math.get_R_from_strike(strike);
        self.L = fault_vector_functions.get_strike_length(xstart, xfinish, ystart, yfinish);
        self.W = fault_vector_functions.get_downdip_width(top, bottom, dipangle);
        self.strike_unit_vector = fault_vector_functions.get_strike_vector(self.strike);  # 3d vector in horiz. plane.
        self.dip_unit_vector = fault_vector_functions.get_dip_vector(strike, dipangle);  # 3d vector.
        self.plane_normal = fault_vector_functions.get_plane_normal(strike, dipangle);  # 3d vector.

    def modify_fault_object(self, xstart=None, xfinish=None, ystart=None, yfinish=None, rtlat=None,
                            reverse=None, tensile=None, potency=None, strike=None, dipangle=None, rake=None,
                            zerolon=None, zerolat=None, top=None, bottom=None):
        """
        Modify the fields in a Pycoulomb.Faults_object.
        """
        xstart = self.xstart if xstart is None else xstart;
        xfinish = self.xfinish if xfinish is None else xfinish;
        ystart = self.ystart if ystart is None else ystart;
        yfinish = self.yfinish if yfinish is None else yfinish;
        rtlat = self.rtlat if rtlat is None else rtlat;
        reverse = self.reverse if reverse is None else reverse;
        tensile = self.tensile if tensile is None else tensile;
        potency = self.potency if potency is None else potency;
        strike = self.strike if strike is None else strike;
        dipangle = self.dipangle if dipangle is None else dipangle;
        rake = self.rake if rake is None else rake;
        zerolon = self.zerolon if zerolon is None else zerolon;
        zerolat = self.zerolat if zerolat is None else zerolat;
        top = self.top if top is None else top;
        bottom = self.bottom if bottom is None else bottom;
        new_fault = Faults_object(xstart=xstart, xfinish=xfinish, ystart=ystart, yfinish=yfinish,
                                  rtlat=rtlat, reverse=reverse, tensile=tensile, potency=potency,
                                  strike=strike, dipangle=dipangle, rake=rake, zerolon=zerolon,
                                  zerolat=zerolat, top=top, bottom=bottom);
        return new_fault;

    def get_fault_center(self):
        """
        Compute the x-y-z coordinates of the center of a PyCoulomb fault patch (a namedtuple)
        """
        center_z = (self.top+self.bottom)/2.0;
        updip_center_x = (self.xstart+self.xfinish)/2.0;
        updip_center_y = (self.ystart+self.yfinish)/2.0;
        vector_mag = self.W*np.cos(np.deg2rad(self.dipangle))/2.0;  # how far the middle is displaced
        # downdip from map-view
        center_point = fault_vector_functions.add_vector_to_point(updip_center_x, updip_center_y, vector_mag,
                                                                  self.strike+90);
        # strike+90 = downdip direction.
        center = [center_point[0], center_point[1], center_z];
        return center;

    def get_fault_slip_moment(self, mu):
        """
        From a source fault object, calculate the seismic moment.
        Must be a finite fault, not a point source.
        Not really used yet, but could be useful in the future.

        :param mu: shear modulus, float, in Pa
        :returns: array with two elements, seismic moment (N-m), and moment magnitude
        """
        if self.potency:  # for the case of point source, we can't do the moment calculation
            return None, None;
        area = self.L * self.W * 1000 * 1000;
        slip = fault_vector_functions.get_vector_magnitude([self.rtlat, self.reverse]);
        seismic_moment = mu * area * slip;
        moment_magnitude = moment_calculations.mw_from_moment(seismic_moment);
        return seismic_moment, moment_magnitude;

    def get_fault_four_corners(self, coords="cartesian"):
        """
        Get the four corners of the object, including updip and downdip.
        depth is fault_object.top
        coords can be "cartesian" or "geographic" for lon/lat
        """
        updip_point0 = [self.xstart, self.ystart];
        updip_point1 = [self.xfinish, self.yfinish];
        vector_mag = self.W*np.cos(np.deg2rad(self.dipangle));  # how far the bottom edge is displaced
        # downdip from map-view
        downdip_point0 = fault_vector_functions.add_vector_to_point(self.xstart, self.ystart, vector_mag,
                                                                    self.strike+90);
        # strike+90 = downdip direction.
        downdip_point1 = fault_vector_functions.add_vector_to_point(self.xfinish, self.yfinish, vector_mag,
                                                                    self.strike+90);

        if coords == 'geographic':
            updip_point0 = fault_vector_functions.xy2lonlat_single(updip_point0[0], updip_point0[1],
                                                                   self.zerolon, self.zerolat);
            updip_point1 = fault_vector_functions.xy2lonlat_single(updip_point1[0], updip_point1[1],
                                                                   self.zerolon, self.zerolat);
            downdip_point0 = fault_vector_functions.xy2lonlat_single(downdip_point0[0], downdip_point0[1],
                                                                     self.zerolon, self.zerolat);
            downdip_point1 = fault_vector_functions.xy2lonlat_single(downdip_point1[0], downdip_point1[1],
                                                                     self.zerolon, self.zerolat);

        x_total = [updip_point0[0], updip_point1[0], downdip_point1[0], downdip_point0[0], updip_point0[0]];
        y_total = [updip_point0[1], updip_point1[1], downdip_point1[1], downdip_point0[1], updip_point0[1]];
        x_updip = [updip_point0[0], updip_point1[0]];
        y_updip = [updip_point0[1], updip_point1[1]];

        return [x_total, y_total, x_updip, y_updip];
