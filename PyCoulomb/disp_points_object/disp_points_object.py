
from Tectonic_Utils.geodesy import euler_pole, insar_vector_functions
from Tectonic_Utils.geodesy import utilities as geod_utilities
import numpy as np


class Displacement_points:
    """
    Displacement_points are individual disp_point elements, can be put into lists of elements.
    dE_obs, Se_obs, etc. in meters.
    Meas_type can be 'GNSS', 'leveling', 'tide_gage', 'survey', 'continuous', 'insar', etc.
    starttime and endtime are optional datetime objects.  Lon is between -180 and 180.
    """
    def __init__(self,  lon, lat, dE_obs, dN_obs, dU_obs, Se_obs=0, Sn_obs=0, Su_obs=0, name=None, starttime=None,
                 endtime=None, refframe=None, meas_type=None):
        self.lon = geod_utilities.wrap_lon(lon);
        self.lat = lat;
        self.dE_obs = dE_obs;  # meters
        self.dN_obs = dN_obs;  # meters
        self.dU_obs = dU_obs;  # meters
        self.Se_obs = Se_obs;  # meters
        self.Sn_obs = Sn_obs;  # meters
        self.Su_obs = Su_obs;  # meters
        self.name = name;
        self.starttime = starttime;
        self.endtime = endtime;
        self.refframe = refframe;
        self.meas_type = meas_type;

    # ------------ PREDICATES -------------- #
    def is_within_bbox(self, bbox):
        """
        :param bbox: [lonW, lonE, latS, latN] floats
        :returns: bool
        """
        if bbox[0] <= self.lon <= bbox[1] and bbox[2] <= self.lat <= bbox[3]:
            return 1;
        else:
            return 0;

    def is_meas_type(self, target_meas_type):
        if self.meas_type == target_meas_type:
            return 1;
        else:
            return 0;

    def has_full_data(self):
        if np.isnan(self.dE_obs) or np.isnan(self.dN_obs) or np.isnan(self.dU_obs):
            return 0;
        else:
            return 1;

    # ------------ REGULAR OBJECT FUNCTIONS -------------- #

    def with_east_as(self, east_value):
        """
        Return a new copy of the object after setting the value of dE_obs to a new value.
        Different from set_east_value(), which modifies the same object in-place.
        """
        obj2 = Displacement_points(lon=self.lon, lat=self.lat, dE_obs=east_value, dN_obs=self.dN_obs,
                                   dU_obs=self.dU_obs, Se_obs=self.Se_obs, Sn_obs=self.Sn_obs, Su_obs=self.Su_obs,
                                   name=self.name, starttime=self.starttime, endtime=self.endtime,
                                   refframe=self.refframe, meas_type=self.meas_type);
        return obj2;

    def set_east_value(self, east_value):
        self.dE_obs = east_value;

    def set_north_value(self, north_value):
        self.dN_obs = north_value;

    def set_vert_value(self, vert_value):
        self.dU_obs = vert_value;

    def multiply_by_value(self, value):
        obj2 = Displacement_points(lon=self.lon, lat=self.lat, dE_obs=value * self.dE_obs, dN_obs=value * self.dN_obs,
                                   dU_obs=value * self.dU_obs, Se_obs=self.Se_obs, Sn_obs=self.Sn_obs,
                                   Su_obs=self.Su_obs, name=self.name, starttime=self.starttime, endtime=self.endtime,
                                   refframe=self.refframe, meas_type=self.meas_type);
        return obj2;

    def project_into_los(self, lkv_e, lkv_n, lkv_u):
        """
        :param lkv_e: east look vector component from ground to satellite
        :param lkv_n: north look vector component from ground to satellite
        :param lkv_u: up look vector component from ground to satellite
        :returns: los_defo, float, in meters
        """
        flight_angle, inc_angle = insar_vector_functions.look_vector2flight_incidence_angles(lkv_e, lkv_n, lkv_u);
        los_defo = insar_vector_functions.def3D_into_LOS(self.dE_obs, self.dN_obs, self.dU_obs, flight_angle, inc_angle)
        return los_defo;

    def translate_point_by_euler_pole(self, euler_pole_components):
        """Rotate a point around a given euler pole.  Euler pole contains (ep_lon, ep_lat, omega).
        Euler pole is in degrees and degrees/Ma."""
        [ep_lon, ep_lat, omega] = euler_pole_components;
        [ep_ve, ep_vn, ep_vu] = euler_pole.point_rotation_by_Euler_Pole([self.lon, self.lat], [ep_lon, ep_lat, omega]);
        obj2 = Displacement_points(lon=self.lon, lat=self.lat, dE_obs=self.dE_obs + ep_ve / 1000,
                                   dN_obs=self.dN_obs + ep_vn / 1000, dU_obs=self.dU_obs + ep_vu / 1000,
                                   Se_obs=self.Se_obs, Sn_obs=self.Sn_obs, Su_obs=self.Su_obs, name=self.name,
                                   starttime=self.starttime, endtime=self.endtime, refframe=self.refframe,
                                   meas_type=self.meas_type);
        return obj2;
