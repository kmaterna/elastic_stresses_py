"""
Functions for slippy IO of faults and slip distributions into list of fault_slip_object dictionaries
Format json: basis1, basis2, length(m), width(m), nlength, nwidth, strike, dip, position [lon, lat, dep], penalty
"""

import json
from Tectonic_Utils.geodesy import fault_vector_functions
from Elastic_stresses_py.PyCoulomb.fault_slip_object import fault_slip_object


def read_faults_json(infile):
    """
    Read all faults from a json file (just geometry. no slip or rake) into a list of fault slip objects.
    It has to convert from fault center to fault corner.
    Faults read from JSON have zero slip.

    :param infile: name of input json file
    :type infile: string
    :returns: list of fault slip objects
    :rtype: list
    """
    fault_list = []
    print("Reading file %s " % infile)
    config_file = open(infile, 'r')
    config = json.load(config_file)
    for key in config.keys():
        length = config[key]["length"] / 1000.0
        width = config[key]["width"] / 1000.0
        center_lon = config[key]["position"][0]
        center_lat = config[key]["position"][1]
        x_start, y_start = fault_vector_functions.add_vector_to_point(0, 0, length / 2,
                                                                      config[key]["strike"] - 180)  # in km
        corner_lon, corner_lat = fault_vector_functions.xy2lonlat(x_start, y_start, center_lon, center_lat)
        depth = -config[key]["position"][2] / 1000
        one_fault = fault_slip_object.FaultSlipObject(strike=config[key]["strike"], dip=config[key]["dip"], depth=depth,
                                                      length=length, width=width, lon=corner_lon, lat=corner_lat,
                                                      rake=0, slip=0, tensile=0, segment=0)
        fault_list.append(one_fault)
    config_file.close()
    print("--> Returning %d fault patches " % len(fault_list))
    return fault_list


def write_faults_json(faults_list, outfile):
    """
    Writes faults to json as receivers with zero slip
    position is lon/lat/depth in meters (negative means below surface)

    :param faults_list: list of fault dictionaries
    :type faults_list: list
    :param outfile: name of output json file
    :type outfile: string
    """
    output = {}
    for k, fault in enumerate(faults_list):
        # Convert the fault (which has top left corner) into a fault with top center coordinate
        label = "fault" + str(k)
        newfault = {}

        corner_lon = fault.lon
        corner_lat = fault.lat
        x_center, y_center = fault_vector_functions.add_vector_to_point(0, 0, fault.length / 2, fault.strike)
        center_lon, center_lat = fault_vector_functions.xy2lonlat(x_center, y_center, corner_lon, corner_lat)

        newfault["strike"] = fault.strike
        newfault["dip"] = fault.dip
        newfault["length"] = fault.length * 1000
        newfault["width"] = fault.width * 1000
        newfault["basis1"] = [1, 0, 0]
        newfault["basis2"] = None
        newfault["Nlength"] = 1
        newfault["Nwidth"] = 1
        newfault["penalty"] = 1
        newfault["position"] = [center_lon, center_lat, -fault.depth*1000]
        output[label] = newfault
    with open(outfile, 'w') as ofile:
        json.dump(output, ofile, indent=4)
    return
