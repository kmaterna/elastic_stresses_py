"""
The functions in this package operate on disp_points objects.
Displacement_points = ['lon', 'lat', 'dE_obs'[m], 'dN_obs'[m], 'dU_obs'[m], 'Se_obs'[m], 'Sn_obs'[m], 'Su_obs'[m],
'name', 'starttime', 'endtime', 'refframe', 'meas_type'];
Disp_points are lists of individual disp_point elements.
"""
import matplotlib.path
import numpy as np
from .disp_points_object import Displacement_points


# ------- LIST COMPREHENSIONS AND WRAPPERS ---------- #

def translate_by_euler_pole(disp_points_list, euler_pole_components):
    """Rotate by euler pole."""
    return [item.translate_point_by_euler_pole(euler_pole_components) for item in disp_points_list];


def mult_disp_points_by(disp_points1, multiplier=-1):
    """
    Flip list of disp_points, or multiply by another value. The metadata and uncertainties will be retained.
    """
    return [item.multiply_by_value(multiplier) for item in disp_points1];


def with_easts_as(disp_points_list, east_array):
    """
    :param disp_points_list: list of disp_point_objects.
    :param east_array: list of floats, matching length with disp_points_list.
    """
    if len(disp_points_list) != len(east_array):
        raise ValueError("Error! Length of east_array not equal to length of disp_points_list.");
    return [disp_points_list[i].with_east_as(east_array[i]) for i in range(len(disp_points_list))];


# ------- REGULAR UTILITY FUNCTIONS ---------- #

def subtract_disp_points(disp_points1, disp_points2, target='all', tol=0.001):
    """
    Subtract two lists of objects (1 minus 2) for the residuals.
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
        res1 = Displacement_points(lon=disp_points1[i].lon, lat=disp_points1[i].lat,
                                   dE_obs=disp_points1[i].dE_obs - disp_points2[i].dE_obs,
                                   dN_obs=disp_points1[i].dN_obs - disp_points2[i].dN_obs,
                                   dU_obs=disp_points1[i].dU_obs - vert_multiplier*disp_points2[i].dU_obs,
                                   Se_obs=disp_points1[i].Se_obs, Sn_obs=disp_points1[i].Sn_obs,
                                   Su_obs=disp_points1[i].Su_obs, starttime=disp_points1[i].starttime,
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
        res1 = Displacement_points(lon=disp_points1[i].lon, lat=disp_points1[i].lat,
                                   dE_obs=disp_points1[i].dE_obs + disp_points2[i].dE_obs,
                                   dN_obs=disp_points1[i].dN_obs + disp_points2[i].dN_obs,
                                   dU_obs=disp_points1[i].dU_obs + disp_points2[i].dU_obs,
                                   Se_obs=disp_points1[i].Se_obs, Sn_obs=disp_points1[i].Sn_obs,
                                   Su_obs=disp_points1[i].Su_obs, starttime=disp_points1[i].starttime,
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
        res1 = Displacement_points(lon=disp_points1[i].lon, lat=disp_points1[i].lat,
                                   dE_obs=disp_points1[i].dE_obs - reference_disp_point.dE_obs,
                                   dN_obs=disp_points1[i].dN_obs - reference_disp_point.dN_obs,
                                   dU_obs=disp_points1[i].dU_obs - vert_multiplier*reference_disp_point.dU_obs,
                                   Se_obs=disp_points1[i].Se_obs, Sn_obs=disp_points1[i].Sn_obs,
                                   Su_obs=disp_points1[i].Su_obs, starttime=disp_points1[i].starttime,
                                   endtime=disp_points1[i].endtime, refframe=disp_points1[i].refframe,
                                   meas_type=disp_points1[i].meas_type);
        new_station_vels.append(res1);
    return new_station_vels;


def generate_grid_of_disp_points(W, E, S, N, xinc, yinc):
    """Create a synthetic grid of disp_point objects, useful for forward-calculations of displacement."""
    points = [];
    x = np.linspace(W, E, int((E - W) / xinc));
    y = np.linspace(S, N, int((N - S) / yinc));
    [x2d, y2d] = np.meshgrid(x, y);
    for ky in range(len(y)):
        for kx in range(len(x)):
            points.append(Displacement_points(lon=x2d[ky][kx], lat=y2d[ky][kx], dE_obs=0, dN_obs=0, dU_obs=0));
    return points;


# ------- MAP / FILTER / REDUCE FUNCTIONS ---------- #

def filter_to_meas_type(obs_points, meas_type='continuous'):
    """Filter a list of disp_points into a list of disp_points with a certain value of meas_type."""
    keep_obs_points = [x for x in obs_points if x.is_meas_type(meas_type)];
    return keep_obs_points;


def filter_to_remove_near_fault(obs_points, fault_points, radius_km=5):
    """
    Filter a list of disp_points to remove those that are very close to a fault, such as a creeping fault.

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


def filter_to_remove_nans(obs_points):
    keep_obs_points = [x for x in obs_points if x.has_full_data()];
    return keep_obs_points;


def filter_by_bounding_box(obs_points, bbox):
    """
    Filter a set of points by bounding box.

    :param obs_points: list of disp_points
    :param bbox: list [W, E, S, N]
    """
    keep_obs_points = [];
    for item in obs_points:
        if item.is_within_bbox(bbox):
            keep_obs_points.append(item);
    print("Filtering to within a bounding box: Returning %d of %d points " % (len(keep_obs_points), len(obs_points)));
    return keep_obs_points;


def filter_to_exclude_bounding_box(obs_points, bbox):
    """
    Filter a set of points to exclude a particular bounding box (such as removing a small volcanic area).

    :param obs_points: list of disp_points
    :param bbox: list [W, E, S, N]
    """
    keep_obs_points = [];
    for item in obs_points:
        if item.is_within_bbox(bbox):
            continue;   # station is within excluded box. Ignore it.
        else:
            keep_obs_points.append(item);
    print("Filtering to exclude bounding box: Returning %d of %d points " % (len(keep_obs_points), len(obs_points)));
    return keep_obs_points;


def extract_particular_station_from_list(disp_points_list, station_lon, station_lat, tol=0.001):
    """Pull one particular station out of a list of stations, with error handling."""
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


def extract_region_from_disp_points(disp_points_list):
    """Operates on a list of disp_points. Returns bbox [W, E, S, N]."""
    lon = np.array([x.lon for x in disp_points_list]);
    lat = np.array([x.lat for x in disp_points_list]);
    region = [np.min(lon), np.max(lon), np.min(lat), np.max(lat)];
    return region;


def station_vel_object_to_disp_points(velfield):
    """
    Convert from StationVel objects from GNSS Python library into disp_points.
    This is the opposite of GNSS_TimeSeries_Viewers.gps_tools.vel_functions.disp_points_to_station_vels().
    """
    disp_points_list = [];
    for item in velfield:
        new_disp_point = Displacement_points(lon=item.elon, lat=item.nlat, dE_obs=item.e/1000, dN_obs=item.n/1000,
                                             dU_obs=item.u/1000, Se_obs=item.se/1000, Sn_obs=item.sn/1000,
                                             Su_obs=item.su/1000, name=item.name, starttime=item.first_epoch,
                                             endtime=item.last_epoch, meas_type=item.meas_type, refframe=item.refframe);
        disp_points_list.append(new_disp_point);
    return disp_points_list;
