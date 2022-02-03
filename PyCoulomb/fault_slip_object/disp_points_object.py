"""
The functions in this package operate on cc.disp_points objects
"""

import numpy as np
from Elastic_stresses_py.PyCoulomb import coulomb_collections as cc
from Tectonic_Utils.geodesy import euler_pole


def subtract_disp_points(disp_points1, disp_points2, target='all'):
    """
    Subtract two lists of objects (1 minus 2) for the residuals
    The metadata and uncertainties for object 1 will be retained.
    If target=='all', a 3-component correction is used.
    If target=='horizontal', only horizontal correction will be applied (disp_points2.dU is treated like 0).
    """
    residuals = [];
    vert_multiplier = 1;
    if target == 'horizontal':
        vert_multiplier = 0;
    for i in range(len(disp_points1)):
        res1 = cc.Displacement_points(lon=disp_points1[i].lon, lat=disp_points1[i].lat,
                                      dE_obs=disp_points1[i].dE_obs - disp_points2[i].dE_obs,
                                      dN_obs=disp_points1[i].dN_obs - disp_points2[i].dN_obs,
                                      dU_obs=disp_points1[i].dU_obs - vert_multiplier*disp_points2[i].dU_obs,
                                      Se_obs=disp_points1[i].Se_obs, Sn_obs=disp_points1[i].Sn_obs,
                                      Su_obs=disp_points1[i].Su_obs, name="",
                                      starttime=disp_points1[i].starttime,
                                      endtime=disp_points1[i].endtime, refframe=disp_points1[i].refframe,
                                      meas_type=disp_points1[i].meas_type);
        residuals.append(res1);
    return residuals;


def add_disp_points(disp_points1, disp_points2):
    """
    add two lists of objects (1 plus 2).
    The metadata and uncertainties for object 1 will be retained.
    """
    sum_disp_points = [];
    for i in range(len(disp_points1)):
        res1 = cc.Displacement_points(lon=disp_points1[i].lon, lat=disp_points1[i].lat,
                                      dE_obs=disp_points1[i].dE_obs + disp_points2[i].dE_obs,
                                      dN_obs=disp_points1[i].dN_obs + disp_points2[i].dN_obs,
                                      dU_obs=disp_points1[i].dU_obs + disp_points2[i].dU_obs,
                                      Se_obs=disp_points1[i].Se_obs, Sn_obs=disp_points1[i].Sn_obs,
                                      Su_obs=disp_points1[i].Su_obs, name="",
                                      starttime=disp_points1[i].starttime,
                                      endtime=disp_points1[i].endtime, refframe=disp_points1[i].refframe,
                                      meas_type=disp_points1[i].meas_type);
        sum_disp_points.append(res1);
    return sum_disp_points;


def mult_minus_one(disp_points1):
    """
    Flip list of disp_points
    The metadata and uncertainties for object 1 will be retained.
    """
    residuals = [];
    for item in disp_points1:
        res1 = cc.Displacement_points(lon=item.lon, lat=item.lat, dE_obs=-item.dE_obs, dN_obs=-item.dN_obs,
                                      dU_obs=-item.dU_obs, Se_obs=item.Se_obs, Sn_obs=item.Sn_obs, Su_obs=item.Su_obs,
                                      name=item.name, starttime=item.starttime, endtime=item.endtime,
                                      refframe=item.refframe, meas_type=item.meas_type);
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


def obs_vs_model_L1_misfit(obs_disp_points, model_disp_points):
    """Implementing one definition of model misfit: L1 norm"""
    all_misfits_m = [];
    all_misfits_norm = [];
    horiz_misfits_m = [];
    horiz_misfits_norm = [];
    for i in range(len(obs_disp_points)):
        E_misfit = abs(obs_disp_points[i].dE_obs - model_disp_points[i].dE_obs)
        N_misfit = abs(obs_disp_points[i].dN_obs - model_disp_points[i].dN_obs)
        U_misfit = abs(obs_disp_points[i].dU_obs - model_disp_points[i].dU_obs)
        norm_E = E_misfit / obs_disp_points[i].Se_obs
        norm_N = N_misfit / obs_disp_points[i].Sn_obs
        norm_U = U_misfit / obs_disp_points[i].Su_obs
        all_misfits_m = all_misfits_m + [E_misfit, N_misfit, U_misfit];
        all_misfits_norm = all_misfits_norm + [norm_E, norm_N, norm_U];
        horiz_misfits_m = horiz_misfits_m + [E_misfit, N_misfit];
        horiz_misfits_norm = horiz_misfits_norm + [norm_E, norm_N];
    avg_misfit_m = np.nanmean(all_misfits_m)
    avg_misfit_norm = np.nanmean(all_misfits_norm)
    avg_horiz_m = np.nanmean(horiz_misfits_m)
    avg_horiz_norm = np.nanmean(horiz_misfits_norm);
    return avg_misfit_m, avg_misfit_norm, avg_horiz_m, avg_horiz_norm;
