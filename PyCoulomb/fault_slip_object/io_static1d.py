"""
Functions to read and write STATIC1D input files for fault slip and displacement at GPS points
"""
import numpy as np
from Tectonic_Utils.geodesy import fault_vector_functions
from . import fault_slip_object


def write_static1D_file(fault_dict_list, disp_points, filename):
    """
    Write a slip source and a list of GPS points into an input file for static1D
    """
    print("Writing static1d file %s " % filename);
    ofile = open(filename, 'w');

    # Choosing first fault.  Possibly only able to do one at a time.
    one_fault = fault_dict_list[0];
    top, bottom = fault_vector_functions.get_top_bottom_from_top(one_fault["depth"], one_fault["width"],
                                                                 one_fault["dip"]);
    ofile.write("{0:<5.1f}".format(np.round(bottom, 1)));
    ofile.write("{0:<6.1f}".format(np.round(top, 1)));
    ofile.write("{0:<6.1f}\n".format(np.round(one_fault["dip"], 1)));
    ofile.write("1\n");
    lons, lats = fault_slip_object.get_four_corners_lon_lat(one_fault);
    fault_lon = lons[2];
    fault_lat = lats[2];  # the deeper edge towards the strike direction
    ofile.write(" %f %f %.2f %.2f %.2f %.2f \n" % (fault_lon, fault_lat, one_fault["length"], one_fault["strike"],
                                                   one_fault["rake"], one_fault["slip"]*100) );  # lon/lat etc

    ofile.write("%d\n" % len(disp_points.lon));
    for i in range(len(disp_points.lon)):
        ofile.write("{0:>13f}".format(disp_points.lon[i]) );
        ofile.write("{0:>13f}\n".format(disp_points.lat[i]));
    ofile.close();
    return;


def read_static1D_input_file(_filename):
    """Not written yet"""
    fault_dict_list = [];
    disp_points = [];
    return fault_dict_list, disp_points;


def read_static1D_output_file(output_filename, _input_filename):
    """
    Read the displacements from the output of a Static-1D calculation.
    Input_filename will eventually be implemented to make the displacements and lon/lat into one object.
    """
    print("Reading file %s " % output_filename);
    ifile = open(output_filename);
    xdisp, ydisp, zdisp = [], [], [];
    for line in ifile:
        xdisp.append(float(line[20:33]));
        ydisp.append(float(line[34:46]));
        zdisp.append(float(line[47:59]));
    ifile.close();
    return xdisp, ydisp, zdisp;
