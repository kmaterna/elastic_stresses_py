"""
Functions to read and write STATIC1D input files for fault slip and displacement at GPS points
"""
import numpy as np
from Tectonic_Utils.geodesy import fault_vector_functions
from . import fault_slip_object
from Elastic_stresses_py.PyCoulomb import coulomb_collections


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
    ofile.write(" %f %f %.2f %.2f %.2f %.2f \n" % (fault_lat, fault_lon, one_fault["length"], one_fault["strike"],
                                                   one_fault["rake"], one_fault["slip"]*100) );  # lon/lat etc

    ofile.write("%d\n" % len(disp_points.lon));
    for i in range(len(disp_points.lon)):
        ofile.write("{0:>13f}".format(disp_points.lat[i]) );
        ofile.write("{0:>13f}\n".format(disp_points.lon[i]));
    ofile.close();
    return;


def read_static1D_input_files(filename, gps_filename=None):
    """
    Read a number of static1d fault segments and gps station locations into objects.
    If the points are inside the same input file, you don't need the second argument
    Object 1: fault dict, in the internal format.
    Object 2: disp points (lats and lons only since this is inputs)

    lon(back top corner),
    lat(back top corner),
    """
    disp_points = read_static1d_disp_points(filename, gps_filename);
    fault_dict_list = [];
    ifile = open(filename);
    headerline = ifile.readline();   # reading header information, first line.
    lower_depth = float(headerline.split()[0]);  # in km
    upper_depth = float(headerline.split()[1]);  # in km
    dip = float(headerline.split()[2]);  # in degrees
    for line in ifile:
        if len(line.split()) == 6:   # reading one fault segment
            lower_lat_corner = float(line.split()[0]);   # in degrees
            lower_lon_corner = float(line.split()[1]);   # in degrees
            length = float(line.split()[2]);  # in km
            strike = float(line.split()[3]);  # in degrees
            rake = float(line.split()[4]);   # in degrees
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

            new_fault = {"strike": strike, "dip": dip, "length": length, "rake": rake, "slip": slip/100, "tensile": 0,
                         "depth": upper_depth, "width": downdip_width, "lon": fault_lon, "lat": fault_lat};
            fault_dict_list.append(new_fault);

    ifile.close();
    return fault_dict_list, disp_points;


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


def read_visco1D_output_file(output_filename, gps_input_filename):
    print("Reading file %s " % output_filename);
    disp_points_only = read_static1d_disp_points(gps_input_filename, None);
    ifile = open(output_filename);
    xdisp, ydisp, zdisp = [], [], [];
    for line in ifile:
        xdisp.append((1/100)*float(line[20:33]));
        ydisp.append((1/100)*float(line[33:46]));
        zdisp.append((1/100)*float(line[46:59]));
    ifile.close();
    modeled_disp_points = coulomb_collections.Displacement_points(lon=disp_points_only.lon, lat=disp_points_only.lat,
                                                                  dE_obs=xdisp, dN_obs=ydisp, dU_obs=zdisp, Se_obs=None,
                                                                  Sn_obs=None, Su_obs=None, name=None);
    return modeled_disp_points;
