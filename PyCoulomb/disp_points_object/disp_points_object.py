
from Tectonic_Utils.geodesy import euler_pole, insar_vector_functions, utilities


class Displacement_points:
    """
    Displacement_points are individual disp_point elements, can be put into lists of elements.
    dE_obs, Se_obs, etc. in meters.
    Meas_type can be 'GNSS', 'leveling', 'tide_gage', 'survey', 'continuous', 'insar', etc.
    If starttime and endtime are used, they should be datetime objects.  Lon is between -180 and 180.
    """
    def __init__(self,  lon, lat, dE_obs, dN_obs, dU_obs, Se_obs, Sn_obs, Su_obs, name=None, starttime=None,
                 endtime=None, refframe=None, meas_type=None):
        self.lon = utilities.wrap_lon(lon);
        self.lat = lat;
        self.dE_obs = dE_obs;
        self.dN_obs = dN_obs;
        self.dU_obs = dU_obs;
        self.Se_obs = Se_obs;
        self.Sn_obs = Sn_obs;
        self.Su_obs = Su_obs;
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

    # ------------ REGULAR FUNCTIONS -------------- #

    def change_east_value(self, east_value):
        obj2 = Displacement_points(lon=self.lon, lat=self.lat, dE_obs=east_value, dN_obs=self.dN_obs,
                                   dU_obs=self.dU_obs, Se_obs=self.Se_obs, Sn_obs=self.Sn_obs, Su_obs=self.Su_obs,
                                   name=self.name, starttime=self.starttime, endtime=self.endtime,
                                   refframe=self.refframe, meas_type=self.meas_type);
        return obj2;

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
        flight_angle, incidence_angle = insar_vector_functions.look_vector2flight_incidence_angles(lkv_e, lkv_n, lkv_u);
        los_defo = insar_vector_functions.def3D_into_LOS(self.dE_obs, self.dN_obs, self.dU_obs, flight_angle,
                                                         incidence_angle);
        return los_defo;

    def translate_point_by_euler_pole(self, euler_pole_components):
        [ep_lon, ep_lat, omega] = euler_pole_components;
        [ep_ve, ep_vn, ep_vu] = euler_pole.point_rotation_by_Euler_Pole([self.lon, self.lat], [ep_lon, ep_lat, omega]);
        obj2 = Displacement_points(lon=self.lon, lat=self.lat, dE_obs=self.dE_obs + ep_ve / 1000,
                                   dN_obs=self.dN_obs + ep_vn / 1000, dU_obs=self.dU_obs + ep_vu / 1000,
                                   Se_obs=self.Se_obs, Sn_obs=self.Sn_obs, Su_obs=self.Su_obs, name=self.name,
                                   starttime=self.starttime, endtime=self.endtime, refframe=self.refframe,
                                   meas_type=self.meas_type);
        return obj2;
