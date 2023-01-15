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

import matplotlib.path
import numpy as np
from .. import coulomb_collections as cc
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
    Add two lists of objects (1 plus 2).
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


def subtract_reference_from_disp_points(disp_points1, reference_disp_point, target='all'):
    """
    Subtract reference point from each disp_point in list.
    The metadata and uncertainties for object 1 will be retained.
    If target=='all', a 3-component correction is used.
    If target=='horizontal', only horizontal correction will be applied (disp_points2.dU is treated like 0).
    """
    new_station_vels = [];
    vert_multiplier = 1;
    if target == 'horizontal':
        vert_multiplier = 0;
    for i in range(len(disp_points1)):
        res1 = cc.Displacement_points(lon=disp_points1[i].lon, lat=disp_points1[i].lat,
                                      dE_obs=disp_points1[i].dE_obs - reference_disp_point.dE_obs,
                                      dN_obs=disp_points1[i].dN_obs - reference_disp_point.dN_obs,
                                      dU_obs=disp_points1[i].dU_obs - vert_multiplier*reference_disp_point.dU_obs,
                                      Se_obs=disp_points1[i].Se_obs, Sn_obs=disp_points1[i].Sn_obs,
                                      Su_obs=disp_points1[i].Su_obs, name="",
                                      starttime=disp_points1[i].starttime,
                                      endtime=disp_points1[i].endtime, refframe=disp_points1[i].refframe,
                                      meas_type=disp_points1[i].meas_type);
        new_station_vels.append(res1);
    return new_station_vels;


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


def generate_grid_of_disp_points(W, E, S, N, xinc, yinc):
    """Create a synthetic grid of disp_point objects, useful for forward-calculations of displacement."""
    points = [];
    x = np.linspace(W, E, int((E - W) / xinc));
    y = np.linspace(S, N, int((N - S) / yinc));
    [x2d, y2d] = np.meshgrid(x, y);
    for ky in range(len(y)):
        for kx in range(len(x)):
            pt1 = cc.Displacement_points(lon=x2d[ky][kx], lat=y2d[ky][kx], dE_obs=0, dN_obs=0, dU_obs=0, Se_obs=0,
                                         Sn_obs=0, Su_obs=0, name='', starttime=None, endtime=None,
                                         refframe=None, meas_type=None);
            points.append(pt1);
    return points;


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
    print("Filtering to within a bounding box: Returning %d of %d points " % (len(keep_obs_points), len(obs_points)));
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


def extract_particular_station_from_list(disp_points_list, station_lon, station_lat, tol=0.001):
    """Pull one particular station out of a list of stations, with error handling"""
    possible_reference_stations = [];
    for item in disp_points_list:
        if abs(item.lon-station_lon) < tol and abs(item.lat-station_lat) < tol:
            possible_reference_stations.append(item);
    if len(possible_reference_stations) == 1:
        return possible_reference_stations[0];   # return just one item
    if len(possible_reference_stations) > 1:
        print(len(possible_reference_stations));
        print([x.lat for x in possible_reference_stations]);
        raise ValueError("Error in referencing this array!  More than 1 station found at the reference location.");
        # in future, we might want to return just one of these possible matches.
    if len(possible_reference_stations) == 0:
        raise ValueError("Error in referencing this array!  Zero stations found at the reference location.");


def station_vel_object_to_disp_points(velfield):
    """
    Convert from StationVel objects from GNSS Python library into disp_points
    This is the opposite of GNSS_TimeSeries_Viewers.gps_tools.vel_functions.disp_points_to_station_vels()
    """
    disp_points_list = [];
    for item in velfield:
        new_disp_point = cc.Displacement_points(lon=item.elon, lat=item.nlat, dE_obs=item.e/1000, dN_obs=item.n/1000,
                                                dU_obs=item.u/1000, Se_obs=item.se/1000, Sn_obs=item.sn/1000,
                                                Su_obs=item.su/1000, name=item.name, starttime=item.first_epoch,
                                                endtime=item.last_epoch, meas_type=item.meas_type,
                                                refframe=item.refframe);
        disp_points_list.append(new_disp_point);
    return disp_points_list;


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


def extract_region_from_disp_points(disp_points_list):
    """Operates on a list of disp_points"""
    lon = np.array([x.lon for x in disp_points_list]);
    lat = np.array([x.lat for x in disp_points_list]);
    region = [np.min(lon), np.max(lon), np.min(lat), np.max(lat)];
    return region;
