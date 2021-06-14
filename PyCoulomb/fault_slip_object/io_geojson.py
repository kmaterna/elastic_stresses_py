""""
Functions for slippy IO of faults and slip distributions into list of fault_slip_object dictionaries

Format json: basis1, basis2, length(m), width(m), nlength, nwidth, strike, dip, position [lon, lat, dep], penalty
Format slippy: lon lat depth[m] strike[deg] dip[deg] length[m] width[m] left-lateral[m] thrust[m] tensile[m]
"""

import numpy as np
import json
from Tectonic_Utils.geodesy import fault_vector_functions


def read_faults_json(infile):
    """
    Read all faults from a json file (just geometry; no slip or rake) into a list of fault dictionaries.
    It has to convert from fault center to fault corner.
    Faults read from JSON have zero slip.

    :param infile: name of input json file
    :type infile: string
    :returns: list of fault dictionaries
    :rtype: list
    """
    fault_list = [];
    config_file = open(infile, 'r')
    config = json.load(config_file);
    for key in config.keys():
        one_fault = {"strike": config[key]["strike"], "dip": config[key]["dip"],
                     "length": config[key]["length"] / 1000.0,
                     "width": config[key]["width"] / 1000.0};
        center_lon = config[key]["position"][0];
        center_lat = config[key]["position"][1];
        x_start, y_start = fault_vector_functions.add_vector_to_point(0, 0, one_fault["length"] / 2,
                                                                      one_fault["strike"] - 180);  # in km
        corner_lon, corner_lat = fault_vector_functions.xy2lonlat(x_start, y_start, center_lon, center_lat);  #
        one_fault["lon"] = corner_lon;
        one_fault["lat"] = corner_lat;
        one_fault["depth"] = -config[key]["position"][2] / 1000;
        one_fault["rake"] = 0;
        one_fault["slip"] = 0;
        one_fault["tensile"] = 0;
        fault_list.append(one_fault);
    return fault_list;


def write_faults_json(faults_list, outfile):
    """
    Writes faults to json as receivers with zero slip

    :param faults_list: list of fault dictionaries
    :type faults_list: list
    :param outfile: name of output json file
    :type outfile: string
    """
    output = {};
    faults = {};
    count = 0;
    for fault in faults_list:
        # Convert the fault (which has top left corner) into a fault with top center coordinate
        corner_lon = fault["lon"]
        corner_lat = fault["lat"]
        x_center, y_center = fault_vector_functions.add_vector_to_point(0, 0, fault["length"] / 2, fault["strike"]);
        center_lon, center_lat = fault_vector_functions.xy2lonlat(x_center, y_center, corner_lon, corner_lat);

        count = count + 1;
        label = "fault" + str(count)
        fault["length"] = fault["length"] * 1000;
        fault["width"] = fault["width"] * 1000;
        fault["basis1"] = [1, 0, 0];
        fault["basis2"] = None;
        fault["Nlength"] = 1;
        fault["Nwidth"] = 1;
        fault["penalty"] = 1;
        fault["position"] = [center_lon, center_lat, fault["depth"]];
        fault.pop("lon");
        fault.pop("lat");
        fault.pop("depth");
        faults[label] = fault;
    output["faults"] = faults;
    output["plotter"] = "gmt";
    with open(outfile, 'w') as ofile:
        json.dump(output, ofile, indent=4);
    return;


def read_slippy_distribution(infile):
    """
    Read a file from the Slippy inversion outputs (lon[degrees] lat[degrees] depth[m] strike[degrees] dip[degrees]
    length[m] width[m] left-lateral[m] thrust[m] tensile[m] segment_num) into a list of fault dictionaries.
    Lon/lat usually refer to the center top of the fault, so it must convert the lon/lat to the top left corner.

    :param infile: name of input slip distribution file
    :type infile: string
    :returns: list of fault dictionaries
    :rtype: list
    """
    print("Reading slippy distribution %s " % infile);
    fault_list = [];
    [lon, lat, depth, strike, dip, length, width, ll_slip,
     thrust_slip, _, _] = np.loadtxt(infile, skiprows=1, unpack=True, dtype={"names": ('lon', 'lat', 'depth', 'strike',
                                                                                       'dip', 'length', 'width', 'ss',
                                                                                       'ds', 'tensile', 'num'),
                                                                             "formats": (float, float, float, float,
                                                                                         float, float, float, float,
                                                                                         float, float, float)});
    for i in range(len(lon)):
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
