"""
Read other fault formats of rectangular patches into fault slip objects

"""
import numpy as np
from Tectonic_Utils.geodesy import fault_vector_functions
from . import io_four_corners


def io_hamling_2017(filename):
    """Read Ian Hamling's 2017 Science paper fault model"""
    print("Reading file %s " % filename);
    fault_dict_list = [];
    [centerlon, centerlat, strike, dip, rake, slip, l_km, top_km, bottom_km, segment] = np.loadtxt(filename,
                                                                                                   unpack=True,
                                                                                                   skiprows=8);
    for i in range(len(centerlon)):
        width = fault_vector_functions.get_downdip_width(top_km[i], bottom_km[i], dip[i]);
        one_fault = {"strike": strike[i], "dip": dip[i], "length": l_km[i], "width": width, "depth": top_km[i],
                     "rake": rake[i], "slip": slip[i], "tensile": 0, "segment": int(segment[i])};
        x_start, y_start = fault_vector_functions.add_vector_to_point(0, 0, one_fault["length"] / 2,
                                                                      one_fault["strike"] - 180);  # in km
        downdip_width_proj = one_fault["width"] * np.cos(np.deg2rad(one_fault["dip"]));
        x_start, y_start = fault_vector_functions.add_vector_to_point(x_start, y_start, downdip_width_proj/2,
                                                                      one_fault["strike"] - 90);
        corner_lon, corner_lat = fault_vector_functions.xy2lonlat(x_start, y_start, centerlon[i], centerlat[i]);

        one_fault["lon"] = corner_lon;
        one_fault["lat"] = corner_lat;
        fault_dict_list.append(one_fault);
    return fault_dict_list;


def io_wallace_sse(filename):
    """
    Read Laura Wallace's slow slip event files for New Zealand
    """
    print("Reading file %s " % filename);
    fault_dict_list = [];
    ifile = open(filename, 'r');
    for line in ifile:
        if line[0] == "#":
            continue;
        if line[0:4] == "> -Z":
            temp = line.split();
            patch_slip_m = float(temp[3])/1000;  # into meters
            patch_tensile_m = float(temp[6])/1000;
            rake = float(temp[7]);
            c1_lon, c1_lat, c1_depth = ifile.readline().split();
            c2_lon, c2_lat, c2_depth = ifile.readline().split();
            c3_lon, c3_lat, c3_depth = ifile.readline().split();
            c4_lon, c4_lat, c4_depth = ifile.readline().split();
            lons = [float(c1_lon), float(c2_lon), float(c3_lon), float(c4_lon)];
            lats = [float(c1_lat), float(c2_lat), float(c3_lat), float(c4_lat)];
            depths = [-float(c1_depth), -float(c2_depth), -float(c3_depth), -float(c4_depth)];
            if np.sum(lons) == 0 and np.sum(lats) == 0 and np.sum(depths) == 0:  # skip pathological case
                continue;
            one_fault = io_four_corners.get_fault_dict_from_four_corners(lons, lats, depths);
            one_fault["slip"] = patch_slip_m;
            one_fault["rake"] = rake
            one_fault["tensile"] = patch_tensile_m;
            fault_dict_list.append(one_fault);
    ifile.close();
    return fault_dict_list;
