"""
The functions in this package operate on cc.disp_points objects
"""

import numpy as np
from Elastic_stresses_py.PyCoulomb import coulomb_collections as cc


def subtract_disp_points(disp_points1, disp_points2):
    """
    Subtract two objects (1 minus 2) for the residuals
    """
    residuals = [];
    for i in range(len(disp_points1)):
        res1 = cc.Displacement_points(lon=disp_points1[i].lon, lat=disp_points1[i].lat,
                                      dE_obs=disp_points1[i].dE_obs - disp_points2[i].dE_obs,
                                      dN_obs=disp_points1[i].dN_obs - disp_points2[i].dN_obs,
                                      dU_obs=disp_points1[i].dU_obs - disp_points2[i].dU_obs,
                                      Se_obs=np.nan, Sn_obs=np.nan, Su_obs=np.nan, name="");
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
