"""
The functions in this package convert between formats of faults/ slip distributions
For all other formats, make sure you build read/write conversion functions into the internal format and test functions.

Written file formats:
    * Format: geojson format for Slippy input faults
    * Format: slippy format, output format for Elastic_stresses_py and Slippy
    * Format: static1d/visco1d, Fred Pollitz's STATIC1D code.
    * Format: four-corners, used for Wei et al. 2015
Computer memory formats:
    * Format: fault_dict dictionary with all the elements for a slip distribution (INTERNAL FORMAT FOR THIS LIBRARY)
    * Format: slip distribution format for Slippy

The internal format here is a dictionary containing:
Fault_Dict:
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
        rake_total = fault_vector_functions.get_rake(rtlat_strike_slip=ss_total, dip_slip=ds_total);
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


def change_fault_slip(fault_dict_list, new_slip, new_rake=None):
    """
    Set the fault slip on a list of fault_slip_dictionaries to something different.
    Can optionally also set the rake on all fault patches to a constant value; otherwise, leave rake unchanged.
    :param fault_dict_list: list of fault_slip_dictionaries
    :param new_slip: float, in meters
    :param new_rake: float, in degrees
    :returns new_list: a list of fault_slip_dictionaries
    """
    new_list = [];
    for item in fault_dict_list:
        if new_rake is None:
            new_rake = item["rake"];
        new_obj = {"strike": item["strike"],
                   "dip": item["dip"],
                   "length": item["length"],
                   "width": item["width"],
                   "lon": item["lon"],
                   "lat": item["lat"],
                   "depth": item["depth"],
                   "tensile": item["tensile"],
                   "slip": new_slip,
                   "rake": new_rake};
        new_list.append(new_obj);
    return new_list;


def write_gmt_fault_file(fault_dict_list, outfile):
    """
    Write the 4 corners of a fault and its slip values into a multi-segment file for plotting in GMT
    """
    print("Writing file %s " % outfile);
    ofile = open(outfile, 'w');
    for fault in fault_dict_list:
        lons, lats = get_four_corners_lon_lat(fault);
        ofile.write("> -Z"+str(fault['slip'])+"\n");
        ofile.write("%f %f\n" % (lons[0], lats[0]));
        ofile.write("%f %f\n" % (lons[1], lats[1]));
        ofile.write("%f %f\n" % (lons[2], lats[2]));
        ofile.write("%f %f\n" % (lons[3], lats[3]));
        ofile.write("%f %f\n" % (lons[0], lats[0]));
    ofile.close();
    return;
