"""
The functions in this package operate on an internal format of faults/ slip distributions
For all other formats, make sure you build read/write conversion AND TEST functions into internal format.

Written file formats:
    * Format: geojson format for Slippy input faults
    * Format: slippy format, output format for Elastic_stresses_py and Slippy
    * Format: static1d/visco1d, Fred Pollitz's STATIC1D code.
    * Format: four-corners, used for Wei et al. 2015
Computer memory formats:
    * Format: fault_dict dictionary with all the elements for a slip distribution (INTERNAL FORMAT FOR THIS LIBRARY)
    * Format: PyCoulomb format, a named tuple (For Elastic_stresses_py.PyCoulomb)

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
from .. import conversion_math
import numpy as np


def get_four_corners_lon_lat(fault_dict_object):
    """
    Return the lon/lat of all 4 corners of a fault_dict_object
    """
    [source] = fault_dict_to_coulomb_fault([fault_dict_object]);
    [x_total, y_total, _, _] = conversion_math.get_fault_four_corners(source);
    lons, lats = fault_vector_functions.xy2lonlat(x_total, y_total, source.zerolon, source.zerolat);
    return lons, lats;


def get_four_corners_lon_lat_multiple(fault_dict_list):
    """
    Return the lon/lat of all 4 corners of a list of fault_dict_object
    Basically the bounding box for this list of fault_dict_objects
    """
    lons_all, lats_all = [], [];
    for item in fault_dict_list:
        [source] = fault_dict_to_coulomb_fault([item]);
        [x_total, y_total, _, _] = conversion_math.get_fault_four_corners(source);
        lons, lats = fault_vector_functions.xy2lonlat(x_total, y_total, source.zerolon, source.zerolat);
        lons_all = lons_all + lons;
        lats_all = lats_all + lats;
    return np.min(lons_all), np.max(lons_all), np.min(lats_all), np.max(lats_all);


def get_updip_corners_lon_lat(fault_dict_object):
    """
    Return the lon/lat of 2 shallow corners of a fault_dict_object
    """
    [source] = fault_dict_to_coulomb_fault([fault_dict_object]);
    [_, _, x_updip, y_updip] = conversion_math.get_fault_four_corners(source);
    lons, lats = fault_vector_functions.xy2lonlat(x_updip, y_updip, source.zerolon, source.zerolat);
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

def filter_by_depth(fault_dict_list, upper_depth, lower_depth):
    """
    Filter a list of fault_dicts to only those that fall within the depth range [upper_depth, lower_depth]
    :param fault_dict_list: list of fault_slip_dictionaries
    :param upper_depth: float, km
    :param lower_depth: float, km
    :returns new_list: a list of fault_slip_dictionaries
    """
    new_list = [];
    for item in fault_dict_list:
        if upper_depth <= item['depth'] <= lower_depth:
            new_list.append(item);
    return new_list;


def write_gmt_fault_file(fault_dict_list, outfile, colorcode='slip'):
    """
    Write the 4 corners of a fault and its slip values into a multi-segment file for plotting in GMT
    By default, colors patches by slip on fault segment.
    Can color patches by their depth.
    """
    print("Writing file %s " % outfile);
    ofile = open(outfile, 'w');
    for fault in fault_dict_list:
        lons, lats = get_four_corners_lon_lat(fault);
        if colorcode == 'depth':
            color_string = str(fault['depth']);
        else:
            color_string = str(fault['slip'])
        ofile.write("> -Z"+color_string+"\n");
        ofile.write("%f %f\n" % (lons[0], lats[0]));
        ofile.write("%f %f\n" % (lons[1], lats[1]));
        ofile.write("%f %f\n" % (lons[2], lats[2]));
        ofile.write("%f %f\n" % (lons[3], lats[3]));
        ofile.write("%f %f\n" % (lons[0], lats[0]));
    ofile.close();
    return;

def write_gmt_surface_trace(fault_dict_list, outfile):
    """
    Write the 2 updip corners of a fault into a multi-segment file for plotting in GMT
    """
    print("Writing file %s " % outfile);
    ofile = open(outfile, 'w');
    for fault in fault_dict_list:
        lons, lats = get_four_corners_lon_lat(fault);
        ofile.write("> -Z\n");
        ofile.write("%f %f\n" % (lons[0], lats[0]));
        ofile.write("%f %f\n" % (lons[1], lats[1]));
    ofile.close();
    return;

def write_gmt_vertical_fault_file(fault_dict_list, outfile):
    """
    Write the vertical coordinates of planar fault patches (length and depth, in local coords instead of lon/lat)
    and associated slip values into a multi-segment file for plotting in GMT.
    Good for vertical faults.  Plots with depth as a negative number.
    Works for only one planar fault segment.
    """
    print("Writing file %s " % outfile);

    # Get origin: extremal patch at the top. First, find bounding box for top points
    depth_array = [x["depth"] for x in fault_dict_list];
    top_row_patches = [x for x in fault_dict_list if x['depth'] == np.nanmin(depth_array)];
    top_row_lon, top_row_lat = [], [];
    for patch in top_row_patches:
        lon_updip, lat_updip = get_updip_corners_lon_lat(patch);
        top_row_lon = top_row_lon + lon_updip;
        top_row_lat = top_row_lat + lat_updip;  # joining two lists
    bbox = [np.nanmin(top_row_lon), np.nanmax(top_row_lon), np.nanmin(top_row_lat), np.nanmax(top_row_lat)];

    # Find fault corner coordinates that are candidates for extremal points on fault. Choose one for origin.
    origin_ll = [np.nan, np.nan];
    for lon, lat in zip(top_row_lon, top_row_lat):
        if lon in bbox and lat in bbox:
            origin_ll = [lon, lat];  # this should be guaranteed to happen twice, once for each end.
            break;

    ofile = open(outfile, 'w');
    for fault in fault_dict_list:
        [source] = fault_dict_to_coulomb_fault([fault], zerolon_system=origin_ll[0], zerolat_system=origin_ll[1]);
        [_, _, x_updip, y_updip] = conversion_math.get_fault_four_corners(source);
        deeper_offset = fault["width"]*np.sin(np.deg2rad(fault["dip"]));
        [xprime, _] = conversion_math.rotate_list_of_points(x_updip, y_updip, 90+fault["strike"]);
        start_x = xprime[0];
        finish_x = xprime[1];
        [rtlat, _dipslip] = fault_vector_functions.get_rtlat_dip_slip(fault['slip'], fault['rake']);   # decision point
        ofile.write("> -Z"+str(-rtlat)+"\n");  # currently writing left-lateral slip as positive
        ofile.write("%f %f\n" % (start_x, -fault["depth"]));
        ofile.write("%f %f\n" % (finish_x, -fault["depth"]));
        ofile.write("%f %f\n" % (finish_x, -fault["depth"]-deeper_offset));
        ofile.write("%f %f\n" % (start_x, -fault["depth"]-deeper_offset));
        ofile.write("%f %f\n" % (start_x, -fault["depth"]));

    ofile.close();
    return;
