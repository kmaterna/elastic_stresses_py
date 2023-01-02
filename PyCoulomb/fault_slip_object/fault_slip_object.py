"""
The functions in this package operate on an internal format of faults/ slip distributions
For all other formats, make sure you build read/write conversion AND TEST functions into internal format.

    * Format(File): geojson format for Slippy input faults
    * Format(File): slippy format, output format for Elastic_stresses_py and Slippy
    * Format(File): static1d/visco1d, Fred Pollitz's STATIC1D code.
    * Format(File): four-corners, used for Wei et al. 2015
    * Format(File): .fsp format, used by USGS fault distributions and srcmod database
    * Format(Memory): fault_dict dictionary with elements for a slip distribution (INTERNAL FORMAT FOR THIS LIBRARY)
    * Format(Memory): PyCoulomb format, a named tuple (For Elastic_stresses_py.PyCoulomb)
"""

# TODO total library: create branch, turn fso into class and fst into class (then bump version when it works)
from .. import coulomb_collections as cc
from Tectonic_Utils.geodesy import fault_vector_functions, haversine
from Tectonic_Utils.seismo import moment_calculations
from .. import conversion_math
import numpy as np
from typing import TypedDict
import collections.abc

class FaultDict(TypedDict):
    """
    The internal format is a dictionary containing:  FaultDict:
    {
        strike(deg),
        dip(deg),
        length(km),
        width(km),
        lon(back top corner),
        lat(back top corner),
        depth(top, km), positive downwards
        rake(deg),
        slip(m),
        tensile(m)
        segment(int)
    }
    If the fault is a receiver fault, we put slip = 0
    """
    strike: float
    dip: float
    length: float
    width: float
    lon: float
    lat: float
    depth: float
    rake: float
    slip: float
    tensile: float
    segment: int


def get_total_slip(one_fault_dict):
    """Helper function to return the total slip amount of a fault object (always > 0)"""
    return one_fault_dict["slip"];

def get_total_slip_mm(one_fault_dict):
    """Helper function to return the total slip amount of a fault object in mm (always > 0)"""
    return one_fault_dict["slip"]*1000;

def get_rtlat_slip(one_fault_dict):
    """Helper function to return the right lateral slip amount of a fault object"""
    [rtlat, _] = fault_vector_functions.get_rtlat_dip_slip(one_fault_dict['slip'], one_fault_dict['rake']);
    return rtlat;

def get_dip_slip(one_fault_dict):
    """Helper function to return the dip slip amount of a fault object"""
    [_, dipslip] = fault_vector_functions.get_rtlat_dip_slip(one_fault_dict['slip'], one_fault_dict['rake']);
    return dipslip;

def get_fault_depth(one_fault_dict):
    """Helper function to return the depth of a fault object"""
    return one_fault_dict["depth"];

def get_fault_rake(one_fault_dict):
    """Helper function to return the depth of a fault object"""
    return one_fault_dict["rake"];

def get_segment(one_fault_dict):
    """Helper function to return the segment_num of a fault object"""
    return one_fault_dict["segment"];

def get_blank_fault(_one_fault_dict):
    """Helper function to an empty color palette"""
    return "";


def get_four_corners_lon_lat(fault_dict_object):
    """
    Return the lon/lat of all 4 corners of a fault_dict_object
    """
    [source] = fault_dict_to_coulomb_fault([fault_dict_object]);
    [x_total, y_total, _, _] = conversion_math.get_fault_four_corners(source);
    lons, lats = fault_vector_functions.xy2lonlat(x_total, y_total, source.zerolon, source.zerolat);
    return lons, lats;


def get_updip_corners_lon_lat(fault_dict_object):
    """
    Return the lon/lat of 2 shallow corners of a fault_dict_object
    """
    [source] = fault_dict_to_coulomb_fault([fault_dict_object]);
    [_, _, x_updip, y_updip] = conversion_math.get_fault_four_corners(source);
    lons, lats = fault_vector_functions.xy2lonlat(x_updip, y_updip, source.zerolon, source.zerolat);
    return lons, lats;


def get_fault_element_distance(fault_dict1, fault_dict2, threedimensional=True):
    """Return the distance between two fault elements, in 3d, in km"""
    h_distance = haversine.distance([fault_dict1["lat"], fault_dict1["lon"]], [fault_dict2["lat"], fault_dict2["lon"]]);
    if threedimensional:
        depth_distance = fault_dict1['depth'] - fault_dict2['depth'];
    else:
        depth_distance = 0;
    distance_3d = np.sqrt(np.square(h_distance) + np.square(depth_distance))
    return distance_3d;


def get_fault_area(fault_dict1):
    """Return the area of a rectangular fault. Units: square meters"""
    A = fault_dict1["width"] * 1000 * fault_dict1["length"] * 1000;
    return A;


def get_fault_moment(fault_dict1, mu=30e9):
    """Get the moment for a fault patch. Units: Newton-meters"""
    A = get_fault_area(fault_dict1);
    d = fault_dict1["slip"];
    moment = moment_calculations.moment_from_muad(mu, A, d);
    return moment;


def change_fault_slip(fault_dict1, new_slip=None, new_rake=None):
    """
    Set the fault slip on a fault_slip_dictionary to something different.
    Can optionally also set the rake; otherwise, leave rake unchanged.

    :param fault_dict1: fault_slip_dictionary
    :param new_slip: float, in meters
    :param new_rake: float, in degrees
    :returns: fault_slip_object
    """
    desired_slip = fault_dict1["slip"] if new_slip is None else new_slip;
    desired_rake = fault_dict1["rake"] if new_rake is None else new_rake;
    new_obj = FaultDict(strike=fault_dict1["strike"], dip=fault_dict1["dip"], length=fault_dict1["length"],
                        width=fault_dict1["width"], lon=fault_dict1["lon"], lat=fault_dict1["lat"],
                        depth=fault_dict1["depth"], tensile=fault_dict1["tensile"],
                        slip=desired_slip, rake=desired_rake, segment=fault_dict1["segment"]);
    return new_obj;


def add_slip_from_two_faults(fault_dict1, fault_dict2):
    """Assuming identical fault_dicts except for different slip amounts. Add slip amounts. """
    ss_1, ds_1 = fault_vector_functions.get_rtlat_dip_slip(fault_dict1["slip"], fault_dict1["rake"]);
    ss_2, ds_2 = fault_vector_functions.get_rtlat_dip_slip(fault_dict2["slip"], fault_dict2["rake"]);
    ss_total = ss_1 + ss_2;  # rtlat
    ds_total = ds_1 + ds_2;  # reverse
    slip_total = fault_vector_functions.get_total_slip(ss_total, ds_total);
    rake_total = fault_vector_functions.get_rake(rtlat_strike_slip=ss_total, dip_slip=ds_total);
    new_item = FaultDict(strike=fault_dict1["strike"], dip=fault_dict1["dip"], length=fault_dict1["length"],
                         width=fault_dict1["width"],
                         lon=fault_dict1["lon"], lat=fault_dict1["lat"], depth=fault_dict1["depth"],
                         tensile=fault_dict1["tensile"] + fault_dict2["tensile"], slip=slip_total, rake=rake_total,
                         segment=fault_dict1["segment"])
    return new_item;


def is_within_depth_range(fault_dict1, upper_depth, lower_depth):
    """
    Filter fault_dict if it falls within depth range [upper_depth, lower_depth]

    :param fault_dict1: fault_slip_dictionary
    :param upper_depth: float, km
    :param lower_depth: float, km
    :returns: bool
    """
    if upper_depth <= fault_dict1['depth'] <= lower_depth:
        return 1;
    else:
        return 0;


"""
Functions that work on lists of fault items
"""


def get_four_corners_lon_lat_multiple(fault_dict_list):
    """
    Return the lon/lat of all 4 corners of a list of fault_dict_object
    Basically the bounding box for this list of fault_dict_objects
    """
    lons_all, lats_all = [], [];
    for item in fault_dict_list:
        lons, lats = get_four_corners_lon_lat(item);
        lons_all = lons_all + lons;
        lats_all = lats_all + lats;
    return np.min(lons_all), np.max(lons_all), np.min(lats_all), np.max(lats_all);


def get_total_moment(fault_dict_object_list, mu=30e9):
    """
    Return the total moment of a list of slip objects, in fault_dict_object
    Moment in newton-meters
    """
    total_moment = 0;
    for item in fault_dict_object_list:
        total_moment += get_fault_moment(item, mu);
    return total_moment;


def get_total_moment_depth_dependent(fault_dict_object_list, depths, mus):
    """Compute total moment using a depth-dependent G calculation"""
    total_moment = 0;
    for item in fault_dict_object_list:
        depth = get_fault_depth(item);
        idx = np.abs(depths - depth).argmin()
        G = mus[idx];
        A = get_fault_area(item);
        d = get_total_slip(item);
        total_moment += moment_calculations.moment_from_muad(G, A, d);
    return total_moment;


def add_two_fault_dict_lists(list1, list2):
    """Assuming identical geometry in the two lists"""
    if len(list1) != len(list2):
        raise Exception("Error! Two fault_dict lists are not identical");
    new_list = [];
    for item1, item2 in zip(list1, list2):
        new_item = add_slip_from_two_faults(item1, item2);
        new_list.append(new_item);
    return new_list;


def change_fault_slip_list(fault_dict_list, new_slip=None, new_rake=None):
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
        new_obj = change_fault_slip(item, new_slip, new_rake);
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
        if is_within_depth_range(item, upper_depth, lower_depth):
            new_list.append(item);
    return new_list;


def get_how_many_segments(fault_dict_list):
    """
    Reduce fault_dict_list into the number of segments (often 1) and the number of individual fault patches

    :param fault_dict_list: list of fault_slip_dictionaries
    """
    segments = [get_segment(x) for x in fault_dict_list];
    num_segments = len(set(segments));
    num_patches = len(segments);
    return num_segments, num_patches;


def filter_by_segment(fault_dict_list, segment_num=0):
    """
    Filter a list of fault_dicts to only those that have segment=x

    :param fault_dict_list: list of fault_slip_dictionaries
    :param segment_num: int
    :returns new_list: a list of fault_slip_dictionaries
    """
    new_list = [];
    for item in fault_dict_list:
        if get_segment(item) == segment_num:
            new_list.append(item);
    return new_list;


def write_gmt_fault_file(fault_dict_list, outfile, color_mappable=get_blank_fault, verbose=True):
    """
    Write the 4 corners of a fault and its slip values into a multi-segment file for plotting in GMT
    By default, does not provide color on the fault patches
    color_mappable can be 1d array of scalars, or a function of the object's class.
    """
    if verbose:
        print("Writing file %s " % outfile);
    ofile = open(outfile, 'w');
    for i, fault in enumerate(fault_dict_list):
        lons, lats = get_four_corners_lon_lat(fault);
        if isinstance(color_mappable, collections.abc.Sequence):
            color_string = "-Z"+str(color_mappable[i]);  # if separately providing the color array
        else:
            color_string = "-Z"+str(color_mappable(fault));  # call the function that you've provided
        ofile.write("> "+color_string+"\n");
        ofile.write("%f %f\n" % (lons[0], lats[0]));
        ofile.write("%f %f\n" % (lons[1], lats[1]));
        ofile.write("%f %f\n" % (lons[2], lats[2]));
        ofile.write("%f %f\n" % (lons[3], lats[3]));
        ofile.write("%f %f\n" % (lons[0], lats[0]));
    ofile.close();
    return;


def write_gmt_surface_trace(fault_dict_list, outfile, verbose=True):
    """
    Write the 2 updip corners of a rectangular fault into a multi-segment file for plotting in GMT
    """
    if verbose:
        print("Writing file %s " % outfile);
    ofile = open(outfile, 'w');
    for fault in fault_dict_list:
        lons, lats = get_four_corners_lon_lat(fault);
        ofile.write("> -Z\n");
        ofile.write("%f %f\n" % (lons[0], lats[0]));
        ofile.write("%f %f\n" % (lons[1], lats[1]));
    ofile.close();
    return;


def write_gmt_vertical_fault_file(fault_dict_list, outfile, color_mappable=get_rtlat_slip):
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
        start_x, finish_x = xprime[0], xprime[1];
        slip_amount = color_mappable(fault);
        ofile.write("> -Z"+str(-slip_amount)+"\n");  # currently writing left-lateral slip as positive
        ofile.write("%f %f\n" % (start_x, -fault["depth"]));
        ofile.write("%f %f\n" % (finish_x, -fault["depth"]));
        ofile.write("%f %f\n" % (finish_x, -fault["depth"]-deeper_offset));
        ofile.write("%f %f\n" % (start_x, -fault["depth"]-deeper_offset));
        ofile.write("%f %f\n" % (start_x, -fault["depth"]));

    ofile.close();
    return;


def construct_intxt_source_patch(fault_dict):
    """
    Turn a fault_slip_object into a string corresponding to a slip_patch in .intxt file.
    """
    fault_patch_string = "Source_Patch: " + str(fault_dict["strike"]) + " " + str(fault_dict["rake"]) + " ";
    fault_patch_string = fault_patch_string + str(fault_dict["dip"]) + " " + str(fault_dict["length"]) + " ";
    fault_patch_string = fault_patch_string + str(fault_dict["width"]) + " " + str(fault_dict["lon"]) + " ";
    fault_patch_string = fault_patch_string + str(fault_dict["lat"]) + " " + str(fault_dict["depth"]) + " ";
    fault_patch_string = fault_patch_string + str(fault_dict["slip"]) + " " + str(fault_dict["tensile"]) + " ";
    return fault_patch_string;


def coulomb_fault_to_fault_dict(source_object):
    """Convert a list of fault objects from Elastic_stresses_py into a list of internal dictionary objects"""
    if not source_object:
        return [];
    fault_dict_list = [];
    for src in source_object:
        if src.potency:
            print("ERROR! Cannot convert a point source into a rectangular source. Skipping...");
            continue;
        lon, lat = fault_vector_functions.xy2lonlat(src.xstart, src.ystart, src.zerolon, src.zerolat);
        slip = fault_vector_functions.get_total_slip(src.rtlat, src.reverse);
        length = fault_vector_functions.get_strike_length(src.xstart, src.xfinish, src.ystart, src.yfinish);
        width = fault_vector_functions.get_downdip_width(src.top, src.bottom, src.dipangle);
        one_fault = FaultDict(strike=src.strike, dip=src.dipangle, depth=src.top, rake=src.rake,
                              slip=slip, length=length, width=width, tensile=src.tensile, lon=lon, lat=lat, segment=0);
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
        one_source = cc.construct_pycoulomb_fault(xstart=startx, xfinish=xfinish, ystart=starty, yfinish=yfinish,
                                                  Kode=100, rtlat=rtlat, reverse=reverse, tensile=onefault['tensile'],
                                                  potency=[], strike=onefault['strike'], dipangle=onefault['dip'],
                                                  zerolon=zerolon, zerolat=zerolat, rake=onefault['rake'],
                                                  top=onefault['depth'], bottom=bottom, comment='');
        source_object.append(one_source);
    return source_object;
