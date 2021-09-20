"""
The functions in this package operate on cc.disp_points objects
"""

import numpy as np
from Elastic_stresses_py.PyCoulomb import coulomb_collections as cc
from Tectonic_Utils.geodesy import euler_pole


def subtract_disp_points(disp_points1, disp_points2):
    """
    Subtract two lists of objects (1 minus 2) for the residuals
    The metadata for object 1 will be retained.
    """
    residuals = [];
    for i in range(len(disp_points1)):
        res1 = cc.Displacement_points(lon=disp_points1[i].lon, lat=disp_points1[i].lat,
                                      dE_obs=disp_points1[i].dE_obs - disp_points2[i].dE_obs,
                                      dN_obs=disp_points1[i].dN_obs - disp_points2[i].dN_obs,
                                      dU_obs=disp_points1[i].dU_obs - disp_points2[i].dU_obs,
                                      Se_obs=np.nan, Sn_obs=np.nan, Su_obs=np.nan, name="",
                                      starttime=disp_points1[i].starttime,
                                      endtime=disp_points1[i].endtime, refframe=disp_points1[i].refframe,
                                      meas_type=disp_points1[i].meas_type);
        residuals.append(res1);
    return residuals;


def add_disp_points(disp_points1, disp_points2):
    """
    add two lists of objects (1 plus 2).
    The metadata for object 1 will be retained.
    """
    sum_disp_points = [];
    for i in range(len(disp_points1)):
        res1 = cc.Displacement_points(lon=disp_points1[i].lon, lat=disp_points1[i].lat,
                                      dE_obs=disp_points1[i].dE_obs + disp_points2[i].dE_obs,
                                      dN_obs=disp_points1[i].dN_obs + disp_points2[i].dN_obs,
                                      dU_obs=disp_points1[i].dU_obs + disp_points2[i].dU_obs,
                                      Se_obs=np.nan, Sn_obs=np.nan, Su_obs=np.nan, name="",
                                      starttime=disp_points1[i].starttime,
                                      endtime=disp_points1[i].endtime, refframe=disp_points1[i].refframe,
                                      meas_type=disp_points1[i].meas_type);
        sum_disp_points.append(res1);
    return sum_disp_points;


def mult_minus_one(disp_points1):
    """
    Flip list of disp_points
    The metadata for object 1 will be retained.
    """
    residuals = [];
    for i in range(len(disp_points1)):
        res1 = cc.Displacement_points(lon=disp_points1[i].lon, lat=disp_points1[i].lat,
                                      dE_obs=-disp_points1[i].dE_obs,
                                      dN_obs=-disp_points1[i].dN_obs,
                                      dU_obs=-disp_points1[i].dU_obs,
                                      Se_obs=np.nan, Sn_obs=np.nan, Su_obs=np.nan, name="",
                                      starttime=disp_points1[i].starttime,
                                      endtime=disp_points1[i].endtime, refframe=disp_points1[i].refframe,
                                      meas_type=disp_points1[i].meas_type);
        residuals.append(res1);
    return residuals;


def station_vel_object_to_disp_points(velfield):
    """
    Convert from StationVel objects from GNSS Python library into disp_points
    """
    disp_points_list = [];
    for item in velfield:
        new_disp_point = cc.Displacement_points(lon=item.elon, lat=item.nlat, dE_obs=item.e/1000, dN_obs=item.n/1000,
                                                dU_obs=item.u/1000, Se_obs=item.se/1000, Sn_obs=item.sn/1000,
                                                Su_obs=item.su/1000, name=item.name);
        disp_points_list.append(new_disp_point);
    return disp_points_list;


def compute_rms(disp_points):
    list_of_values = [];
    for point in disp_points:
        list_of_values.append(point.dE_obs);
        list_of_values.append(point.dN_obs);
        list_of_values.append(point.dU_obs);
    return np.sqrt(np.mean(np.square(list_of_values)));


def translate_by_euler_pole(disp_points_list, euler_pole_components):
    """Rotate by euler pole"""
    new_disp_points_list = [];
    [ep_lon, ep_lat, omega] = euler_pole_components;
    for item in disp_points_list:
        [ep_ve, ep_vn, ep_vu] = euler_pole.point_rotation_by_Euler_Pole([item.lon, item.lat], [ep_lon, ep_lat, omega]);
        res1 = cc.Displacement_points(lon=item.lon, lat=item.lat, dE_obs=item.dE_obs + ep_ve/1000,
                                      dN_obs=item.dN_obs + ep_vn/1000, dU_obs=item.dU_obs + ep_vu/1000,
                                      Se_obs=item.Se_obs, Sn_obs=item.Sn_obs, Su_obs=item.Su_obs, name=item.name,
                                      starttime=item.starttime, endtime=item.endtime, refframe=item.refframe,
                                      meas_type=item.meas_type);
        new_disp_points_list.append(res1);
    return new_disp_points_list;
