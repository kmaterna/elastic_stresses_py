"""
Read other fault formats of rectangular patches into fault slip objects

"""
import numpy as np
from Tectonic_Utils.geodesy import fault_vector_functions


def io_hamling_2017(filename):
    """Read Ian Hamling's 2017 Science paper fault model"""
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
