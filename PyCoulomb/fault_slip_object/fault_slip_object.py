"""
The functions in this package convert between formats of faults/ slip distributions

* Format 1: json format for Slippy
* Format 2: slip distribution format for Slippy
* Format 3: .intxt, Kathryn Materna's format for Elastic_stresses_py.  USE ELASTIC_STRESSES_PY TO READ THAT FORMAT.
* Format 4: static1d, Fred Pollitz's STATIC1D code. USE ELASTIC_STRESSES_PY TO READ THAT FORMAT.
* ANY OTHER FORMATS: make sure you build a conversion function into the internal format.

The internal format here is a dictionary containing:
Utility_Dict:
{
    strike(deg),
    dip(deg),
    length(km),
    width(km),
    lon(back top corner),
    lat(back top corner),
    depth(top, km),
    rake(deg),
    slip(m),
    tensile(m)
}
If the fault is a receiver fault, we put slip = 0
"""

from .io_pycoulomb import fault_dict_to_coulomb_fault
from Tectonic_Utils.geodesy import fault_vector_functions
from Tectonic_Utils.seismo import moment_calculations
from Elastic_stresses_py.PyCoulomb import conversion_math
import numpy as np


def get_four_corners_lon_lat(fault_dict_object):
    """
    Return the lon/lat of all 4 corners of a fault_dict_object
    """
    [source] = fault_dict_to_coulomb_fault([fault_dict_object]);
    [x_total, y_total, _, _] = conversion_math.get_fault_four_corners(source);
    lons, lats = fault_vector_functions.xy2lonlat(x_total, y_total, source.zerolon, source.zerolat);
    return lons, lats;


def get_total_moment(fault_dict_object_list, mu=30e9):
    """
    Return the total moment of a list of slip objects, in fault_dict_object
    Moment in newton-meters
    """
    total_moment = 0;
    for item in fault_dict_object_list:
        A = item["width"] * 1000 * item["length"] * 1000;
        d = item["slip"];
        total_moment += moment_calculations.moment_from_muad(mu, A, d);
    return total_moment;


def get_total_moment_depth_dependent(fault_dict_object_list, depths, mus):
    """Compute total moment using a depth-dependent G calculation"""
    total_moment = 0;
    for item in fault_dict_object_list:
        depth = item["depth"];
        idx = np.abs(depths - depth).argmin()
        G = mus[idx];
        A = item["width"] * 1000 * item["length"] * 1000;
        d = item["slip"];
        total_moment += moment_calculations.moment_from_muad(G, A, d);
    return total_moment;


def add_two_fault_dict_lists(list1, list2):
    """Assuming identical geometry in the two lists"""
    if len(list1) != len(list2):
        raise Exception("Error! Two fault_dict lists are not identical");
    new_list = [];
    for item1, item2 in zip(list1, list2):
        ss_1, ds_1 = fault_vector_functions.get_rtlat_dip_slip(item1["slip"], item1["rake"]);
        ss_2, ds_2 = fault_vector_functions.get_rtlat_dip_slip(item2["slip"], item2["rake"]);
        ss_total = ss_1 + ss_2;  # rtlat
        ds_total = ds_1 + ds_2;  # reverse
        slip_total = fault_vector_functions.get_total_slip(ss_total, ds_total);
        rake_total = fault_vector_functions.get_rake(ss_total, ds_total);
        new_item = {"strike": item1["strike"],
                    "dip": item1["dip"],
                    "length": item1["length"],
                    "width": item1["width"],
                    "lon": item1["lon"],
                    "lat": item1["lat"],
                    "depth": item1["depth"],
                    "tensile": item1["tensile"]+item2["tensile"],
                    "slip": slip_total,
                    "rake": rake_total };
        new_list.append(new_item);
    return new_list;


def change_fault_slip(fault_dict_list, new_slip):
    """Change the fault slip to something different."""
    new_list = [];
    for item in fault_dict_list:
        new_obj = {"strike": item["strike"],
                   "dip": item["dip"],
                   "length": item["length"],
                   "width": item["width"],
                   "lon": item["lon"],
                   "lat": item["lat"],
                   "depth": item["depth"],
                   "tensile": item["tensile"],
                   "slip": new_slip,
                   "rake": item["rake"]};
        new_list.append(new_obj);
    return new_list;
