""""
Functions for SRCMOD IO of faults and slip distributions into list of fault_slip_object dictionaries

"""
import numpy as np
from Tectonic_Utils.geodesy import fault_vector_functions


def read_srcmod_distribution(infile):
    """
    Let's assume that the lon/lat/depth given in the SRCMOD .fsp file is for the top center of the fault patch.
    This function doesn't have a unit test yet.

    :param infile: name of input slip distribution file, defined to be the '.fsp' file
    :type infile: string
    :returns: list of fault dictionaries
    :rtype: list
    """
    print("Reading SRCMOD distribution %s " % infile);
    fault_list = [];
    overall_strike, overall_dip, total_len_km, nx, total_width_km, nz = 0, 90, 10, 10, 10, 10;   # defaults.
    ifile = open(infile, 'r');
    for line in ifile:
        temp = line.split();
        if len(temp) == 0:
            continue;
        if len(temp) <= 3 and line[0] == '%':
            continue;
        if len(temp) > 3:
            if line[0] == '%' and temp[1] == 'Mech':
                overall_strike = float(temp[5]);
                overall_dip = float(temp[8]);
            if line[0] == '%' and temp[1] == 'Size':
                total_len_km = float(temp[5]);
                total_width_km = float(temp[9]);
            if line[0] == '%' and temp[3] == 'Nx':
                nx = int(temp[5]);
                nz = int(temp[8]);
            if line[0] != '%':
                lon_top_center = float(temp[1]);
                lat_top_center = float(temp[0]);
                depth_top_center = float(temp[4]);
                depth_top, _ = fault_vector_functions.get_top_bottom_from_center(depth_top_center, total_width_km/nz,
                                                                                 overall_dip)
                slip_m = float(temp[5]);
                rake = float(temp[7]);
                one_fault = {"strike": overall_strike, "dip": overall_dip, "length": total_len_km/nx,
                             "width": total_width_km/nz,
                             "depth": depth_top_center, "rake": rake, "slip": slip_m, "tensile": 0};
                x_start, y_start = fault_vector_functions.add_vector_to_point(0, 0, one_fault["length"] / 2,
                                                                              one_fault["strike"] - 180);  # in km
                _downdip_width_proj = one_fault["width"]*np.cos(np.deg2rad(overall_dip));
                # x_start, y_start = fault_vector_functions.add_vector_to_point(x_start, y_start, downdip_width_proj/2,
                #                                                               one_fault["strike"] - 90);
                # ^^ offset the fault location for center. Optional/unknown.
                corner_lon, corner_lat = fault_vector_functions.xy2lonlat(x_start, y_start, lon_top_center,
                                                                          lat_top_center);
                one_fault["lon"] = corner_lon;
                one_fault["lat"] = corner_lat;
                fault_list.append(one_fault);
    ifile.close();
    print("  -->Returning %d fault segments" % len(fault_list));
    return fault_list;


def write_srcmod_distribution(_faults_list, outfile):
    """
    :param _faults_list: a list of fault dictionaries
    :param outfile: name of output file.
    """
    print("Not writing file %s " % outfile);
    return;
