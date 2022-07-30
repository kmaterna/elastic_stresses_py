""""
Functions for conversion of faults and slip distributions from PyCoulomb into list of fault_slip_object dictionaries

Some of these faults may come in through .intxt files. Use io_intxt to read them.

Format intxt:
    strike
    rake
    dip
    length(km)
    width(km)
    updip_corner_lon
    updip_corner_lat
    updip_corner_dep(km)
    slip(m)
"""

from .. import coulomb_collections as cc
from Tectonic_Utils.geodesy import fault_vector_functions
import numpy as np


def coulomb_fault_to_fault_dict(source_object):
    """Convert a list of fault objects from Elastic_stresses_py into a list of internal dictionary objects"""
    fault_dict_list = [];
    for src in source_object:
        if src.potency:
            print("ERROR! Cannot convert a point source into a rectangular source. Skipping...");
            continue;
        one_fault = {"strike": src.strike, "dip": src.dipangle, "depth": src.top,
                     "rake": fault_vector_functions.get_rake(rtlat_strike_slip=src.rtlat, dip_slip=src.reverse),
                     "slip": fault_vector_functions.get_total_slip(src.rtlat, src.reverse),
                     "tensile": src.tensile,
                     "length": fault_vector_functions.get_strike_length(src.xstart, src.xfinish, src.ystart,
                                                                        src.yfinish),
                     "width": fault_vector_functions.get_downdip_width(src.top, src.bottom, src.dipangle)};
        lon, lat = fault_vector_functions.xy2lonlat(src.xstart, src.ystart, src.zerolon, src.zerolat);
        one_fault["lon"] = lon;
        one_fault["lat"] = lat
        fault_dict_list.append(one_fault);
    return fault_dict_list;


def fault_dict_to_coulomb_fault(fault_dict_list, zerolon_system=None, zerolat_system=None):
    """
    Convert a list of internal dictionary objects into a list of source objects for Elastic_stresses_py
    By default, the bottom corner of the fault is the center of the coordinate system, but
    Parameters zerolon_system and zerolat_system can be passed in for system with 1+ faults.
    """
    source_object = [];
    if len(fault_dict_list) > 1 and zerolon_system is None:
        print("Warning! You are converting multiple faults to cartesian without defining a system coordinate system!");
    for onefault in fault_dict_list:
        zerolon = onefault['lon'] if not zerolon_system else zerolon_system;
        zerolat = onefault['lat'] if not zerolat_system else zerolat_system;
        _top, bottom = fault_vector_functions.get_top_bottom_from_top(onefault['depth'], onefault['width'],
                                                                      onefault['dip']);
        [startx, starty] = fault_vector_functions.latlon2xy_single(onefault['lon'], onefault['lat'], zerolon, zerolat);

        rtlat, reverse = fault_vector_functions.get_rtlat_dip_slip(onefault['slip'], onefault['rake']);
        xfinish, yfinish = fault_vector_functions.add_vector_to_point(startx, starty, onefault['length'],
                                                                      onefault['strike']);
        one_source = cc.construct_fault_object(xstart=startx, xfinish=xfinish, ystart=starty, yfinish=yfinish, Kode=100,
                                               rtlat=rtlat, reverse=reverse, tensile=onefault['tensile'],
                                               potency=[], strike=onefault['strike'],
                                               dipangle=onefault['dip'], zerolon=zerolon, zerolat=zerolat,
                                               rake=onefault['rake'], top=onefault['depth'], bottom=bottom, comment='');
        source_object.append(one_source);
    return source_object;


def read_pycoulomb_displacements(filename):
    lon, lat, disp_x_Okada, disp_y_Okada, disp_z_Okada = np.loadtxt(filename, skiprows=1,
                                                                    usecols=(0, 1, 2, 3, 4), unpack=True);
    disp_points = [];
    for i in range(len(lon)):
        disp_point = cc.Displacement_points(lon=lon[i], lat=lat[i], dE_obs=disp_x_Okada[i], dN_obs=disp_y_Okada[i],
                                            dU_obs=disp_z_Okada[i], Se_obs=np.nan, Sn_obs=np.nan, Su_obs=np.nan,
                                            name="");
        disp_points.append(disp_point);
    return disp_points;
