""""
Functions for slippy IO of faults and slip distributions into list of fault_slip_object dictionaries

Format slippy: lon lat depth[m] strike[deg] dip[deg] length[m] width[m] left-lateral[m] thrust[m] tensile[m]
"""

import numpy as np
from Tectonic_Utils.geodesy import fault_vector_functions


def read_slippy_distribution(infile, desired_segment=-1):
    """
    Read a file from the Slippy inversion outputs (lon[degrees] lat[degrees] depth[m] strike[degrees] dip[degrees]
    length[m] width[m] left-lateral[m] thrust[m] tensile[m] segment_num) into a list of fault dictionaries.
    Lon/lat usually refer to the center top of the fault, so it must convert the lon/lat to the top left corner.

    :param infile: name of input slip distribution file
    :type infile: string
    :param desired_segment: starting at 0, which fault segment do we want to return? default of -1 means all.
    :type desired_segment: int
    :returns: list of fault dictionaries
    :rtype: list
    """
    print("Reading slippy distribution %s " % infile);
    fault_list = [];
    [lon, lat, depth, strike, dip, length, width, ll_slip,
     thrust_slip, _, num] = np.loadtxt(infile, skiprows=1, unpack=True, dtype={"names": ('lon', 'lat', 'depth',
                                                                                         'strike', 'dip', 'length',
                                                                                         'width', 'ss', 'ds', 'tensile',
                                                                                         'num'),
                                                                               "formats": (float, float, float, float,
                                                                                           float, float, float, float,
                                                                                           float, float, float)});
    for i in range(len(lon)):
        if desired_segment == -1 or desired_segment == num[i]:
            one_fault = {"strike": strike[i], "dip": dip[i], "length": length[i] / 1000, "width": width[i] / 1000,
                         "depth": -depth[i] / 1000};
            center_lon = lon[i];
            center_lat = lat[i];
            x_start, y_start = fault_vector_functions.add_vector_to_point(0, 0, one_fault["length"] / 2,
                                                                          one_fault["strike"] - 180);  # in km
            corner_lon, corner_lat = fault_vector_functions.xy2lonlat(x_start, y_start, center_lon, center_lat);
            one_fault["lon"] = corner_lon;
            one_fault["lat"] = corner_lat;
            one_fault["rake"] = fault_vector_functions.get_rake(ll_slip[i], thrust_slip[i]);
            one_fault["slip"] = fault_vector_functions.get_total_slip(ll_slip[i], thrust_slip[i]);
            one_fault["tensile"] = 0;
            fault_list.append(one_fault);
    return fault_list;


def write_slippy_distribution(faults_list, outfile):
    """
    :param faults_list: a list of fault dictionaries
    :param outfile: name of output file.
    Caveat: can only do one fault segment right now.  That part of the read/write cycle is lossy.
    """
    print("Writing file %s " % outfile);
    ofile = open(outfile, 'w');
    ofile.write("# lon[degrees] lat[degrees] depth[m] strike[degrees] dip[degrees] length[m] width[m] left-lateral[m] "
                "thrust[m] tensile[m] segment_num\n");
    for item in faults_list:
        x_center, y_center = fault_vector_functions.add_vector_to_point(0, 0, item["length"] / 2, item["strike"]);
        center_lon, center_lat = fault_vector_functions.xy2lonlat(x_center, y_center, item["lon"], item["lat"]);
        rtlat_slip, dip_slip = fault_vector_functions.get_rtlat_dip_slip(item["slip"], item["rake"]);
        tensile_slip = 0;
        ofile.write("%f %f %f " % (center_lon, center_lat, item["depth"]*-1000) );
        ofile.write("%f %f %f %f %f %f %f 0 \n" % (item["strike"], item["dip"], item["length"]*1000,
                                                   item["width"]*1000, -1*rtlat_slip, -1*dip_slip, tensile_slip) );
    ofile.close();
    return;
