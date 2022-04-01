"""
The functions in this package operate on cc.disp_points objects
Coulomb collections Displacement_points:
Displacement_points = collections.namedtuple('Disp_Points', [
    'lon', 'lat',
    'dE_obs', 'dN_obs', 'dU_obs',
    'Se_obs', 'Sn_obs', 'Su_obs',
    'name', 'starttime', 'endtime', 'refframe', 'meas_type'], defaults=(None,) * 13);
Disp_points are now lists of individual disp_point elements
Displacements are in meters
"""

import numpy as np
import matplotlib.path
from Elastic_stresses_py.PyCoulomb import coulomb_collections as cc
from Tectonic_Utils.geodesy import euler_pole


def subtract_disp_points(disp_points1, disp_points2, target='all', tol=0.001):
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
        if abs(disp_points1[i].lon - disp_points2[i].lon) > tol:
            raise ValueError("Error! Lists of disp points not matching.  Cannot subtract.");
        if abs(disp_points1[i].lat - disp_points2[i].lat) > tol:
            raise ValueError("Error! Lists of disp points not matching.  Cannot subtract.");
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


def add_disp_points(disp_points1, disp_points2, tol=0.001):
    """
    add two lists of objects (1 plus 2).
    The metadata and uncertainties for object 1 will be retained.
    """
    sum_disp_points = [];
    for i in range(len(disp_points1)):
        if abs(disp_points1[i].lon - disp_points2[i].lon) > tol:
            raise ValueError("Error! Lists of disp points not matching.  Cannot subtract.");
        if abs(disp_points1[i].lat - disp_points2[i].lat) > tol:
            raise ValueError("Error! Lists of disp points not matching.  Cannot subtract.");
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


def mult_disp_points_by(disp_points1, multiplier=-1):
    """
    Flip list of disp_points
    The metadata and uncertainties for object 1 will be retained.
    """
    residuals = [];
    for item in disp_points1:
        res1 = cc.Displacement_points(lon=item.lon, lat=item.lat,
                                      dE_obs=multiplier*item.dE_obs, dN_obs=multiplier*item.dN_obs,
                                      dU_obs=multiplier*item.dU_obs, Se_obs=item.Se_obs, Sn_obs=item.Sn_obs,
                                      Su_obs=item.Su_obs,
                                      name=item.name, starttime=item.starttime, endtime=item.endtime,
                                      refframe=item.refframe, meas_type=item.meas_type);
        residuals.append(res1);
    return residuals;


def filter_to_meas_type(obs_points, meas_type='continuous'):
    """Filter a list of disp_points into a list of disp_points with a certain value of meas_type"""
    keep_obs_points = [];
    for item in obs_points:
        if item.meas_type == meas_type:
            keep_obs_points.append(item);
    return keep_obs_points;


def filter_to_meas_type_by_second_table(obs_points, second_table_obs_points, meas_type='continuous'):
    """Keep points on first table if they are a certain meas_type in second table"""
    keep_obs_points = [];
    for i, item in enumerate(second_table_obs_points):
        if item.meas_type == meas_type:
            keep_obs_points.append(obs_points[i]);
    return keep_obs_points;


def filter_to_remove_near_fault(obs_points, fault_points, radius_km=5):
    """
    Filter a list of disp_points to remove those that are very close to a fault, such as a creeping fault

    :param obs_points: list of disp_point objects
    :param fault_points: list of 2-tuples, lon-lat vertices of fault trace in creeping section
    :param radius_km: float, radius around fault trace, in km.  Ex: 5 km on each side of the fault.
    """
    far_field_points = [];
    newpath = matplotlib.path.Path(vertices=fault_points);
    radius_deg = 2 * radius_km * (1 / 111.000);  # 2x because we're taking 'radius' km on each side
    for item in obs_points:
        if newpath.contains_point((item.lon, item.lat), radius=radius_deg):
            continue;
        else:
            far_field_points.append(item);
    print("Filtering creep: Returning %d of %d points " % (len(far_field_points), len(obs_points)));
    return far_field_points;


def filter_by_bounding_box(obs_points, bbox):
    """
    Filter a set of points by bounding box

    :param obs_points: list of disp_points
    :param bbox: list [W, E, S, N]
    """
    keep_obs_points = [];
    for item in obs_points:
        if bbox[0] < item.lon < bbox[1]:
            if bbox[2] < item.lat < bbox[3]:
                keep_obs_points.append(item);
    return keep_obs_points;

def filter_by_bounding_box(obs_points, bbox):
    """
    Filter a set of points by bounding box

    :param obs_points: list of disp_points
    :param bbox: list [W, E, S, N]
    """
    keep_obs_points = [];
    for item in obs_points:
        if bbox[0] < item.lon < bbox[1]:
            if bbox[2] < item.lat < bbox[3]:
                keep_obs_points.append(item);
    return keep_obs_points;


def filter_to_exclude_bounding_box(obs_points, bbox):
    """
    Filter a set of points to exculde a particular bounding box (such as removing a small volcanic area)

    :param obs_points: list of disp_points
    :param bbox: list [W, E, S, N]
    """
    keep_obs_points = [];
    for item in obs_points:
        if bbox[0] < item.lon < bbox[1] and bbox[2] < item.lat < bbox[3]:
            continue;   # station is within exculded box. Ignore it.
        else:
            keep_obs_points.append(item);
    print("Filtering to exculde bounding box: Returning %d of %d points " % (len(keep_obs_points), len(obs_points)));
    return keep_obs_points;


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
    all_misfits_m, all_misfits_norm = [], [];
    horiz_misfits_m, horiz_misfits_norm = [], [];
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


def obs_vs_model_L2_misfit(obs_disp_points, model_disp_points):
    """Implementing one definition of model misfit: L1 norm"""
    all_misfits_m, all_misfits_norm = [], [];
    horiz_misfits_m, horiz_misfits_norm = [], [];
    for i in range(len(obs_disp_points)):
        E_misfit = np.square(obs_disp_points[i].dE_obs - model_disp_points[i].dE_obs)
        N_misfit = np.square(obs_disp_points[i].dN_obs - model_disp_points[i].dN_obs)
        U_misfit = np.square(obs_disp_points[i].dU_obs - model_disp_points[i].dU_obs)
        norm_E = E_misfit / np.square(obs_disp_points[i].Se_obs)
        norm_N = N_misfit / np.square(obs_disp_points[i].Sn_obs)
        norm_U = U_misfit / np.square(obs_disp_points[i].Su_obs)
        all_misfits_m = all_misfits_m + [E_misfit, N_misfit, U_misfit];
        all_misfits_norm = all_misfits_norm + [norm_E, norm_N, norm_U];
        horiz_misfits_m = horiz_misfits_m + [E_misfit, N_misfit];
        horiz_misfits_norm = horiz_misfits_norm + [norm_E, norm_N];
    avg_misfit_m = np.nanmean(all_misfits_m)
    avg_misfit_norm = np.nanmean(all_misfits_norm)
    avg_horiz_m = np.nanmean(horiz_misfits_m)
    avg_horiz_norm = np.nanmean(horiz_misfits_norm);
    return avg_misfit_m, avg_misfit_norm, avg_horiz_m, avg_horiz_norm;
