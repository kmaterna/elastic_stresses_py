"""
Special input functions for disp_point_objects.
More specific than the basic ones listed in io_additionals.py
"""

import numpy as np
from .disp_points_object import Displacement_points

def read_USGS_file(filename):
    """Read a list of offsets produced by USGS's Jerry Svarc"""
    list_of_disp_points = [];
    [lon, lat, E, N, E_sig, N_sig, U, U_sig] = np.loadtxt(filename, unpack=True, skiprows=1,
                                                          usecols=(0, 1, 2, 3, 4, 5, 8, 9))
    for i in range(len(lon)):
        item = Displacement_points(lon=lon[i], lat=lat[i], dE_obs=E[i]/1000, dN_obs=N[i]/1000, dU_obs=U[i]/1000,
                                   Se_obs=E_sig[i]/1000, Sn_obs=N_sig[i]/1000, Su_obs=U_sig[i]/1000, meas_type='gnss');
        list_of_disp_points.append(item);
    print("Reading %s and returning %d points " % (filename, len(list_of_disp_points)) );
    return list_of_disp_points;


def read_pycoulomb_displacements(filename):
    lon, lat, disp_x_Okada, disp_y_Okada, disp_z_Okada = np.loadtxt(filename, skiprows=1,
                                                                    usecols=(0, 1, 2, 3, 4), unpack=True);
    disp_points = [];
    for i in range(len(lon)):
        disp_point = Displacement_points(lon=lon[i], lat=lat[i], dE_obs=disp_x_Okada[i], dN_obs=disp_y_Okada[i],
                                         dU_obs=disp_z_Okada[i], Se_obs=np.nan, Sn_obs=np.nan, Su_obs=np.nan);
        disp_points.append(disp_point);
    return disp_points;


def read_nota_offsets_file(filename):
    """Read the kalts file created by NOTA for each earthquake."""
    disp_points = [];
    start = 0;
    with open(filename) as ifile:
        for oneline in ifile:
            if ' mm ' in oneline and ' deg ' in oneline:
                start = 1;
                continue;
            if start and len(oneline.split()) > 1:
                oneline = oneline.split();
                lon, lat, dE, dN = float(oneline[0]), float(oneline[1]), float(oneline[2])/1000, float(oneline[3])/1000;
                Se, Sn = float(oneline[4])/1000, float(oneline[5])/1000;
                dU, Su = float(oneline[7])/1000, float(oneline[8])/1000;
                name = oneline[9];
                new_disp = Displacement_points(lon=lon, lat=lat, dE_obs=dE, dN_obs=dN, Se_obs=Se, Sn_obs=Sn,
                                               dU_obs=dU, Su_obs=Su, name=name);
                disp_points.append(new_disp);
    return disp_points;
