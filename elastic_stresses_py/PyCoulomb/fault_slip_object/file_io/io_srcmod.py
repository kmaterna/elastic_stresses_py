"""
Functions for SRCMOD IO of faults and slip distributions into list of fault_slip_objects

"""
import numpy as np
from tectonic_utils.geodesy import fault_vector_functions
import tectonic_utils.seismo.moment_calculations as moment_calcs
from .. import fault_slip_object


def read_srcmod_distribution(infile):
    """
    Let's assume that the lon/lat/depth given in the SRCMOD .fsp file is for the top center of the fault patch.
    This function doesn't have a unit test yet.

    :param infile: name of input slip distribution file, defined to be the '.fsp' file
    :type infile: string
    :returns: list of fault slip objects
    :rtype: list
    """
    print("Reading SRCMOD distribution %s " % infile)

    fault_list = []
    overall_strike, overall_dip, total_len_km, total_width_km, nx, nz, segnum = 0, 90, 10, 10, 10, 10, 0   # defaults.
    segment_list, nx_list, nz_list, rake_col = determine_nx_nz_for_multiple_segments(infile)
    nx, nz, segnum = nx_list[0], nz_list[0], segment_list[0]

    ifile = open(infile, 'r')
    for line in ifile:
        temp = line.split()
        if len(temp) == 0:   # skip nonexistent lines
            continue
        if len(temp) <= 3 and line[0] == '%':  # skip short comment lines
            continue
        if len(temp) > 3:  # for real lines of content...
            if line[0:6] == "% Mech":
                overall_strike = float(temp[5])
                overall_dip = float(temp[8])
            if line[0:6] == "% Size":
                total_len_km = float(temp[5])
                total_width_km = float(temp[9])
            if "% SEGMENT #" in line:  # Enter into segment definition for multi-segment file
                segment_number = int(temp[3].strip(":"))
                overall_strike = float(temp[6])
                overall_dip = float(temp[10])
                length_width_line = ifile.readline()  # subsequent line has information about length and width
                total_len_km = float(length_width_line.split()[3])
                total_width_km = float(length_width_line.split()[7])
                segnum = segment_list[segment_number-1]  # one-indexed segments, zero-indexed python
                nx = nx_list[segment_number-1]
                nz = nz_list[segment_number-1]  # for multi-segment files
                print("segnum, strike, dip, len, width, nx, nz")
                print(segnum, overall_strike, overall_dip, total_len_km, total_width_km, nx, nz)
            if line[0] != '%':
                one_fault = read_srcmod_line(line, overall_strike, overall_dip, total_len_km, total_width_km, nx, nz,
                                             segment=segnum, rake_col=rake_col)
                fault_list.append(one_fault)
    ifile.close()
    print("--> Returning %d fault patches" % len(fault_list))
    print("Mw = ", moment_calcs.mw_from_moment(fault_slip_object.get_total_moment(fault_list)))
    return fault_list


def read_srcmod_line(line, overall_strike, overall_dip, total_len_km, total_width_km, nx, nz, segment=0,
                     rake_col=7):
    # read one content line of srcmod files
    temp = line.split()
    lon_top_center, lat_top_center, depth_top_center = float(temp[1]), float(temp[0]), float(temp[4])
    patch_width = total_width_km / nz
    depth_top, _ = fault_vector_functions.get_top_bottom_from_center(depth_top_center, patch_width, overall_dip)
    slip_m, rake = float(temp[5]), float(temp[rake_col])   # RAKE IS COLUMN 6 FOR MULTISEGMENT
    strike, dip = overall_strike, overall_dip
    length = total_len_km / nx
    x_start, y_start = fault_vector_functions.add_vector_to_point(0, 0, length / 2, strike - 180)  # in km
    _downdip_width_proj = patch_width * np.cos(np.deg2rad(overall_dip))
    # x_start, y_start = fault_vector_functions.add_vector_to_point(x_start, y_start, downdip_width_proj/2,
    #                                                               one_fault["strike"] - 90)
    # ^^ offset the fault location for center. Optional/unknown.
    corner_lon, corner_lat = fault_vector_functions.xy2lonlat(x_start, y_start, lon_top_center, lat_top_center)
    one_fault = fault_slip_object.FaultSlipObject(strike=strike, dip=dip, length=length, width=patch_width,
                                                  depth=depth_top_center, rake=rake, slip=slip_m, tensile=0,
                                                  lon=corner_lon, lat=corner_lat, segment=segment)
    return one_fault


def determine_nx_nz_for_multiple_segments(filename):
    """Returns 3 lists of integers that are num_segments long"""
    segment_list, nx_list, nz_list = [], [], []
    single_segment_flag = determine_if_single_segment(filename)
    if single_segment_flag:  # For a one-segment file
        total_patch_depths = np.loadtxt(filename, unpack=True, comments="%", usecols=[4, ])
        nz = int(len(set(total_patch_depths)))
        nz_list = [nz]
        nx_list = [int(len(total_patch_depths) / nz)]
        segment_list = [1]   # for single-segment files
        rake_col = 7  # column for rake data
    else:  # multi-segment rupture: how many nx and nz?
        rake_col = 6  # column for rake data
        ifile = open(filename, 'r')
        for oneline in ifile:
            if "% SEGMENT #" in oneline:
                segment_number = int(oneline.split()[3].strip(":"))
                entering_segment = 1
                depth_list = []
                while entering_segment:
                    newline = ifile.readline()
                    if is_segment_break(newline) and len(depth_list) > 0:   # finishing a loop
                        segment_list.append(segment_number)
                        nz = int(len(set(depth_list)))
                        nx = int(len(depth_list) / nz)
                        nz_list.append(nz)
                        nx_list.append(nx)
                        break
                    if len(newline) > 0 and newline[0] != '%':
                        depth_list.append(float(newline.split()[4]))
                    if len(newline) == 0:
                        break
        ifile.close()
    return segment_list, nx_list, nz_list, rake_col


def is_segment_break(line):
    """Filter for the various ways that .fsp files encode segment breaks in multi-segment files"""
    is_segment_break_line = 0
    if '% -------------------------------' in line:  # has leading space
        is_segment_break_line = 1
    if '%-------------------------------' in line:  # has no leading space
        is_segment_break_line = 1
    return is_segment_break_line


def determine_if_single_segment(filename):
    # Let's determine whether the file is multi-segment or not.
    single_segment_flag = 1
    ifile = open(filename, 'r')
    for line in ifile:
        if "MULTISEGMENT MODEL" in line:
            single_segment_flag = 0
    ifile.close()
    return single_segment_flag
