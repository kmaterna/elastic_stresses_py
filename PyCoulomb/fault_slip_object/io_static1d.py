"""
Functions to read and write STATIC1D input files for fault slip and displacement at GPS points
"""
import numpy as np
from Tectonic_Utils.geodesy import fault_vector_functions
from . import fault_slip_object
from Elastic_stresses_py.PyCoulomb import coulomb_collections


def write_static1D_source_file(fault_dict_list, disp_points, filename):
    """
    Write a slip source and a list of GPS points into an input file for static1D
    Constraint: each fault in this fault_dict_list should have the same dip, top depth, and bottom depth
    If you have more than one of these, you should write more than one source file.
    """
    print("Writing static1d file %s " % filename);
    ofile = open(filename, 'w');

    # Setting header information: top depth, bottom depth, and dip
    one_fault = fault_dict_list[0];
    top, bottom = fault_vector_functions.get_top_bottom_from_top(one_fault["depth"], one_fault["width"],
                                                                 one_fault["dip"]);
    ofile.write("{0:<5.1f}".format(np.round(bottom, 1)));
    ofile.write("{0:<6.1f}".format(np.round(top, 1)));
    ofile.write("{0:<6.1f}\n".format(np.round(one_fault["dip"], 1)));

    ofile.write("%d \n" % len(fault_dict_list) );
    for fault in fault_dict_list:
        line = write_fault_slip_line_static1d_visco1d(fault);
        ofile.write(line);

    ofile.write("%d\n" % len(disp_points.lon));
    for i in range(len(disp_points.lon)):
        ofile.write("{0:>13f}".format(disp_points.lat[i]) );
        ofile.write("{0:>13f}\n".format(disp_points.lon[i]));
    ofile.close();
    return;


def read_static1D_source_file(filename, gps_filename=None):
    """
    Read a number of static1d fault segments and gps station locations into objects.
    If the points are inside the same input file, you don't need the second argument
    Object 1: fault dict, in the internal format.
    Object 2: disp points (lats and lons only since this is inputs)
    """
    disp_points = read_static1d_disp_points(filename, gps_filename);
    fault_dict_list = [];
    ifile = open(filename);
    headerline = ifile.readline();   # reading header information, skipping first line.
    lower_depth = float(headerline.split()[0]);  # in km
    upper_depth = float(headerline.split()[1]);  # in km
    dip = float(headerline.split()[2]);  # in degrees
    for line in ifile:
        if len(line.split()) == 6:   # reading one fault segment
            new_fault = read_fault_slip_line_static1d_visco1d(line, upper_depth, lower_depth, dip);
            fault_dict_list.append(new_fault);
    ifile.close();
    return fault_dict_list, disp_points;


def read_fault_slip_line_static1d_visco1d(line, upper_depth, lower_depth, dip):
    """
    read a line from fred's format of faults into my format of faults
    """
    lower_lat_corner = float(line.split()[0]);  # in degrees
    lower_lon_corner = float(line.split()[1]);  # in degrees
    length = float(line.split()[2]);  # in km
    strike = float(line.split()[3]);  # in degrees
    rake = float(line.split()[4]);  # in degrees
    slip = float(line.split()[5]);  # in cm
    downdip_width = fault_vector_functions.get_downdip_width(upper_depth, lower_depth, dip);
    vector_mag = downdip_width * np.cos(np.deg2rad(dip));  # how far the bottom edge is displaced
    upper_corner_along_strike = fault_vector_functions.add_vector_to_point(0, 0, vector_mag, strike - 90);
    upper_corner_back_edge = fault_vector_functions.add_vector_to_point(upper_corner_along_strike[0],
                                                                        upper_corner_along_strike[1],
                                                                        length, -strike);
    fault_lon, fault_lat = fault_vector_functions.xy2lonlat_single(upper_corner_back_edge[0],
                                                                   upper_corner_back_edge[1], lower_lon_corner,
                                                                   lower_lat_corner);
    new_fault = {"strike": strike, "dip": dip, "length": length, "rake": rake, "slip": slip / 100, "tensile": 0,
                 "depth": upper_depth, "width": downdip_width, "lon": fault_lon, "lat": fault_lat};
    return new_fault;


def write_fault_slip_line_static1d_visco1d(one_fault):
    """
    write a line of fred's format of faults from my fault_dictionary
    """
    lons, lats = fault_slip_object.get_four_corners_lon_lat(one_fault);
    fault_lon = lons[2];
    fault_lat = lats[2];  # the deeper edge towards the strike direction
    writestring = (" %f %f %.2f %.2f %.2f %.2f \n" % (fault_lat, fault_lon, one_fault["length"], one_fault["strike"],
                                                      one_fault["rake"], one_fault["slip"]*100) );  # lon/lat etc
    return writestring;


def read_static1d_disp_points(filename1, filename2):
    """
    General function to read static1d lat/lon pairs
    It seems that static1d can work with its gps points located in a single file with the source faults,
    or with a separate file called "latlon.inDEF".
    """
    disp_points = read_disp_points_from_static1d(filename1);
    if len(disp_points.lon) == 0:
        disp_points = read_disp_points_from_static1d(filename2);
    return disp_points;


def read_latloninDEF(gps_filename):
    """
    Read gps station locations from static1d inputs (latlon.inDEF) into a disp_points object.
    """
    print("Reading file %s" % gps_filename);
    [lat, lon] = np.loadtxt(gps_filename, skiprows=1, unpack=True);
    disp_points = coulomb_collections.Displacement_points(lon=lon, lat=lat, dE_obs=None, dN_obs=None, dU_obs=None,
                                                          Se_obs=None, Sn_obs=None, Su_obs=None, name=None);
    print("Returning %d lat/lon pairs " % len(lon));
    return disp_points;


def read_disp_points_from_static1d(filename):
    """
    Read gps station locations from static1d inputs into a disp_points object.
    """
    print("Reading file %s" % filename);
    lat, lon = [], [];
    ifile = open(filename, 'r');
    for line in ifile:
        if len(line.split()) == 2:
            lat.append(float(line.split()[0]));
            lon.append(float(line.split()[1]));
    ifile.close();
    disp_points = coulomb_collections.Displacement_points(lon=lon, lat=lat, dE_obs=None, dN_obs=None, dU_obs=None,
                                                          Se_obs=None, Sn_obs=None, Su_obs=None, name=None);
    print("Returning %d lat/lon pairs " % len(lon));
    return disp_points;


def write_disp_points_static1d(disp_points, filename):
    print("Writing %d points in file %s" % (len(disp_points.lon), filename));
    ofile = open(filename, 'w');
    ofile.write("%d\n" % (len(disp_points.lon)) );
    for i in range(len(disp_points.lon)):
        ofile.write('%f %f\n' % (disp_points.lat[i], disp_points.lon[i]));
    ofile.close();
    return;


def read_static1D_output_file(output_filename, gps_input_filename):
    """
    Read the displacements from the output of a Static-1D calculation into a disp_points object.
    Returns displacements in m.
    """
    print("Reading file %s " % output_filename);
    ifile = open(output_filename);
    xdisp, ydisp, zdisp = [], [], [];
    for line in ifile:
        xdisp.append((1/100)*float(line[20:33]));
        ydisp.append((1/100)*float(line[33:46]));
        zdisp.append((1/100)*float(line[46:59]));
    ifile.close();
    disp_points_only = read_static1d_disp_points(gps_input_filename, None);
    modeled_disp_points = coulomb_collections.Displacement_points(lon=disp_points_only.lon, lat=disp_points_only.lat,
                                                                  dE_obs=xdisp, dN_obs=ydisp, dU_obs=zdisp, Se_obs=None,
                                                                  Sn_obs=None, Su_obs=None, name=None);
    return modeled_disp_points;


def write_visco1D_source_file(fault_dict_list, filename):
    """
    Writing the slightly more complicated visco1d source file. Very similar to static1d source files, but
    a few extra parameters included.
    """
    print("Writing static1d file %s " % filename);
    ofile = open(filename, 'w');
    ofile.write("Source fault\n");   # this format has a header line

    # Setting general information: top depth, bottom depth, and dip
    one_fault = fault_dict_list[0];
    top, bottom = fault_vector_functions.get_top_bottom_from_top(one_fault["depth"], one_fault["width"],
                                                                 one_fault["dip"]);
    ofile.write("{0:<5.1f}".format(np.round(bottom, 1)));
    ofile.write("{0:<6.1f}".format(np.round(top, 1)));
    ofile.write("{0:<6.1f}\n".format(np.round(one_fault["dip"], 1)));
    ofile.write("1900. 1900. 2000. 1000. 500.0\n");  # Hard coding for simplicity:
    # earthquake cycle begins at 1900,
    # sample at year 1990 to 2000, the periodicity is 500 years.
    ofile.write("%d \n" % len(fault_dict_list) );
    for fault in fault_dict_list:
        line = write_fault_slip_line_static1d_visco1d(fault);
        ofile.write(line);
    ofile.write('1\n');  # for velocities instead of displacements
    ofile.write('0\n');  # evaluate at depth specified in earlier runs of green's functions
    ofile.close();
    return;
