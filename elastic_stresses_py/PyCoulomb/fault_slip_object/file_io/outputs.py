
import collections.abc
from ... import conversion_math
from ..fault_slip_object import fault_object_to_coulomb_fault
import numpy as np


def get_blank_fault_function(x):
    return x.get_blank_fault()


def get_total_slip(x):
    return x.get_total_slip()


def get_right_lateral_slip(x):
    return x.get_rtlat_slip()


def get_left_lateral_slip(x):
    return -x.get_rtlat_slip()


def write_gmt_fault_file(fault_object_list, outfile, color_mappable=get_blank_fault_function, verbose=True):
    """
    Write the 4 corners of a fault and its slip values into a multi-segment file for plotting in GMT.
    By default, does not provide color on the fault patches.

    :param fault_object_list: list of fault elements
    :param outfile: string
    :param color_mappable: 1d array of scalars, or a function that takes an object of the fault's type.
    :param verbose: bool
    """
    if verbose:
        print("Writing file %s " % outfile)
    ofile = open(outfile, 'w')
    for i, fault in enumerate(fault_object_list):
        lons, lats = fault.get_four_corners_lon_lat()
        if isinstance(color_mappable, collections.abc.Sequence):
            color_string = "-Z"+str(color_mappable[i])  # if separately providing the color array
        else:
            color_string = "-Z"+str(color_mappable(fault))  # call the function that you've provided
        ofile.write("> "+color_string+"\n")
        ofile.write("%f %f\n" % (lons[0], lats[0]))
        ofile.write("%f %f\n" % (lons[1], lats[1]))
        ofile.write("%f %f\n" % (lons[2], lats[2]))
        ofile.write("%f %f\n" % (lons[3], lats[3]))
        ofile.write("%f %f\n" % (lons[0], lats[0]))
    ofile.close()
    return


def write_gmt_surface_trace(fault_object_list, outfile, verbose=True):
    """
    Write the 2 updip corners of a rectangular fault into a multi-segment file for plotting in GMT.

    :param fault_object_list: list of fault elements
    :param outfile: string
    :param verbose: bool
    """
    if verbose:
        print("Writing file %s " % outfile)
    ofile = open(outfile, 'w')
    for fault in fault_object_list:
        lons, lats = fault.get_four_corners_lon_lat()
        ofile.write("> -Z\n")
        ofile.write("%f %f\n" % (lons[0], lats[0]))
        ofile.write("%f %f\n" % (lons[1], lats[1]))
    ofile.close()
    return


def write_gmt_vertical_fault_file(fault_object_list, outfile, color_mappable=get_blank_fault_function,
                                  desired_rotation_strike=None, flip_x=False):
    """
    Write the vertical coordinates of planar fault patches (length and depth, in local coords instead of lon/lat).
    and associated slip values into a multi-segment file for plotting in GMT.
    Good for vertical faults.  Plots with depth as a negative number.
    Works for only one planar fault segment.

    :param fault_object_list: list of fault_slip_objects
    :param outfile: string, filename
    :param color_mappable: function, such as the example functions on the top of this file
    :param desired_rotation_strike: optional, provide a given strike for the coordinate transformation
    :param flip_x: boolean, in case the strike is rotated 180° off the desired strike
    """
    print("Writing file %s " % outfile)

    # Get origin: extremal patch at the top. First, find bounding box for top points
    depth_array = [x.depth for x in fault_object_list]
    top_row_patches = [x for x in fault_object_list if x.depth == np.nanmin(depth_array)]
    top_row_lon, top_row_lat = [], []
    for patch in top_row_patches:
        lon_updip, lat_updip = patch.get_updip_corners_lon_lat()
        top_row_lon = top_row_lon + lon_updip
        top_row_lat = top_row_lat + lat_updip  # joining two lists
    bbox = [np.nanmin(top_row_lon), np.nanmax(top_row_lon), np.nanmin(top_row_lat), np.nanmax(top_row_lat)]

    # Find fault corner coordinates that are candidates for extremal points on fault. Choose one for origin.
    origin_ll = [np.nan, np.nan]
    for lon, lat in zip(top_row_lon, top_row_lat):
        if lon in bbox and lat in bbox:
            origin_ll = [lon, lat]  # this should be guaranteed to happen twice, once for each end.
            break

    ofile = open(outfile, 'w')
    for fault in fault_object_list:
        if desired_rotation_strike is None:
            rot_theta = fault.strike
        else:
            rot_theta = desired_rotation_strike
        if flip_x:
            mult = -1
        else:
            mult = 1
        [source] = fault_object_to_coulomb_fault([fault], zerolon_system=origin_ll[0], zerolat_system=origin_ll[1])
        [_, _, x_updip, y_updip] = source.get_fault_four_corners()
        deeper_offset = fault.width*np.sin(np.deg2rad(fault.dip))
        [xprime, _] = conversion_math.rotate_list_of_points(x_updip, y_updip, 90+rot_theta)
        start_x, finish_x = xprime[0], xprime[1]
        slip_amount = color_mappable(fault)
        ofile.write("> -Z"+str(slip_amount)+"\n")
        ofile.write("%f %f\n" % (mult * start_x, -fault.depth))
        ofile.write("%f %f\n" % (mult * finish_x, -fault.depth))
        ofile.write("%f %f\n" % (mult * finish_x, -fault.depth-deeper_offset))
        ofile.write("%f %f\n" % (mult * start_x, -fault.depth-deeper_offset))
        ofile.write("%f %f\n" % (mult * start_x, -fault.depth))

    ofile.close()
    return
