"""
Special input functions for disp_point_objects
More specific than the basic ones listed in io_additionals.py
"""

import numpy as np
from .. import coulomb_collections as cc

def read_USGS_file(filename):
    """Read a list of offsets produced by USGS's Jerry Svarc"""
    list_of_disp_points = [];
    [lon, lat, E, N, E_sig, N_sig, U, U_sig] = np.loadtxt(filename, unpack=True, skiprows=1,
                                                          usecols=(0, 1, 2, 3, 4, 5, 8, 9))
    for i in range(len(lon)):
        item = cc.Displacement_points(lon=lon[i], lat=lat[i], dE_obs=E[i]/1000,
                                      dN_obs=N[i]/1000, dU_obs=U[i]/1000,
                                      Se_obs=E_sig[i]/1000, Sn_obs=N_sig[i]/1000, Su_obs=U_sig[i]/1000, name='');
        list_of_disp_points.append(item);
    print("Reading %s and returning %d points " % (filename, len(list_of_disp_points)) );
    return list_of_disp_points;


def read_pycoulomb_displacements(filename):
    lon, lat, disp_x_Okada, disp_y_Okada, disp_z_Okada = np.loadtxt(filename, skiprows=1,
                                                                    usecols=(0, 1, 2, 3, 4), unpack=True);
    disp_points = [];
    for i in range(len(lon)):
        disp_point = cc.Displacement_points(lon=lon[i], lat=lat[i], dE_obs=disp_x_Okada[i], dN_obs=disp_y_Okada[i],
                                            dU_obs=disp_z_Okada[i], Se_obs=np.nan, Sn_obs=np.nan, Su_obs=np.nan,
                                            name="", starttime=None, endtime=None, meas_type=None, refframe=None);
        disp_points.append(disp_point);
    return disp_points;
