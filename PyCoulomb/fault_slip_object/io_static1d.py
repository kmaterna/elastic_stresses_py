"""
Functions to read and write STATIC1D input files for fault slip and displacement at GPS points
"""
import numpy as np
from Tectonic_Utils.geodesy import fault_vector_functions
from . import fault_slip_object
from Elastic_stresses_py.PyCoulomb import coulomb_collections as cc


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

    ofile.write("%d\n" % len(disp_points));
    for point in disp_points:
        ofile.write("{0:>13f}".format(point.lat) );
        ofile.write("{0:>13f}\n".format(point.lon));
    ofile.close();
    return;


def read_static1D_source_file(filename, gps_filename=None, headerlines=0):
    """
    Read a number of static1d fault segments and gps station locations into objects.
    If the points are inside the same input file, you don't need the second argument
    Object 1: fault dict, in the internal format.
    Object 2: disp points (lats and lons only since this is inputs)
    """
    if gps_filename:
        disp_points = read_static1d_disp_points(filename, gps_filename);
    else:
        disp_points = [];
    fault_dict_list = [];
    ifile = open(filename);
    headerline = ifile.readline();
    for i in range(headerlines):
        headerline = ifile.readline();   # skipping any other header information.
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
    for Visco1D, the slip field is pretty meaningless
    """
    lower_lat_corner = float(line.split()[0]);  # in degrees
    lower_lon_corner = float(line.split()[1]);  # in degrees
    length = float(line.split()[2]);  # in km
    strike = float(line.split()[3]);  # in degrees
    rake = float(line.split()[4]);  # in degrees
    slip = float(line.split()[5]);  # in cm
    downdip_width = fault_vector_functions.get_downdip_width(upper_depth, lower_depth, dip);
    vector_mag = downdip_width * np.cos(np.deg2rad(dip));  # how far bottom edge is displaced from top edge
    upper_corner_along_strike = fault_vector_functions.add_vector_to_point(0, 0, vector_mag, strike - 90);
    upper_corner_back_edge = fault_vector_functions.add_vector_to_point(upper_corner_along_strike[0],
                                                                        upper_corner_along_strike[1],
                                                                        length, strike+180);
    fault_lon, fault_lat = fault_vector_functions.xy2lonlat_single(upper_corner_back_edge[0],
                                                                   upper_corner_back_edge[1], lower_lon_corner,
                                                                   lower_lat_corner);
    new_fault = {"strike": strike, "dip": dip, "length": length, "rake": rake, "slip": slip / 100, "tensile": 0,
                 "depth": upper_depth, "width": downdip_width, "lon": fault_lon, "lat": fault_lat};
    return new_fault;


def write_fault_slip_line_static1d_visco1d(one_fault):
    """
    write a line of Fred's format of faults from internal fault_dictionary
    """
    lons, lats = fault_slip_object.get_four_corners_lon_lat(one_fault);
    fault_lon = lons[2];
    fault_lat = lats[2];  # the deeper edge towards the strike direction
    writestring = (" %f %f %.2f %.2f %.2f %.2f \n" % (fault_lat, fault_lon, one_fault["length"], one_fault["strike"],
                                                      one_fault["rake"], one_fault["slip"]*100) );  # lon/lat etc
    # Slip written in cm
    return writestring;


def read_static1d_disp_points(filename1, filename2):
    """
    General function to read static1d lat/lon pairs
    It seems that static1d can work with its gps points located in a single file with the source faults,
    or with a separate file called "latlon.inDEF".
    """
    disp_points = read_disp_points_from_static1d(filename1);
    if len(disp_points) == 0:
        disp_points = read_disp_points_from_static1d(filename2);
    return disp_points;


def read_latloninDEF(gps_filename):
    """
    Read gps station locations from static1d inputs (latlon.inDEF) into a list of disp_points objects.
    """
    disp_points = [];
    [lat, lon] = np.loadtxt(gps_filename, skiprows=1, unpack=True);
    for i in range(len(lon)):
        disp_point = cc.Displacement_points(lon=lon[i], lat=lat[i], dE_obs=np.nan, dN_obs=np.nan, dU_obs=np.nan,
                                            Se_obs=np.nan, Sn_obs=np.nan, Su_obs=np.nan, name="");
        disp_points.append(disp_point);
    print("Reading file %s... %d lat/lon pairs" % (gps_filename, len(disp_points)));
    return disp_points;


def read_disp_points_from_static1d(filename):
    """
    Read gps station locations from static1d inputs into a disp_points object.
    """
    disp_points = [];
    ifile = open(filename, 'r');
    for line in ifile:
        if len(line.split()) == 2:
            disp_point = cc.Displacement_points(lon=float(line.split()[1]), lat=float(line.split()[0]),
                                                dE_obs=np.nan, dN_obs=np.nan, dU_obs=np.nan,
                                                Se_obs=np.nan, Sn_obs=np.nan, Su_obs=np.nan, name="");
            disp_points.append(disp_point);
    ifile.close();
    print("Reading file %s... %d lat/lon pairs" % (filename, len(disp_points)));
    return disp_points;


def write_disp_points_static1d(disp_points, filename):
    """A very small function for taking cc.displacements_points into static1D format"""
    print("Writing %d points in file %s" % (len(disp_points), filename));
    ofile = open(filename, 'w');
    ofile.write("%d\n" % (len(disp_points)) );
    for point in disp_points:
        ofile.write('%f %f\n' % (point.lat, point.lon));
    ofile.close();
    return;


def write_stationvel_points_static1d(stationvels, filename):
    """A very small function for taking gnss stationvels into static1D format"""
    print("Writing %d points in file %s" % (len(stationvels), filename));
    ofile = open(filename, 'w');
    ofile.write("%d\n" % (len(stationvels)) );
    for point in stationvels:
        ofile.write('%f %f\n' % (point.nlat, point.elon));
    ofile.close();
    return;


def read_static1D_output_file(output_filename, gps_input_filename):
    """
    Read the displacements from the output of a Static-1D calculation into a disp_points object.
    Returns displacements in m.
    """
    ifile = open(output_filename);
    xdisp, ydisp, zdisp, modeled_disp_points = [], [], [], [];
    for line in ifile:
        xdisp.append((1/100)*float(line[20:33]));
        ydisp.append((1/100)*float(line[33:46]));
        zdisp.append((1/100)*float(line[46:59]));
    ifile.close();
    disp_points_only = read_static1d_disp_points(gps_input_filename, None);
    for i in range(len(disp_points_only)):
        modeled_disp_point = cc.Displacement_points(lon=disp_points_only[i].lon, lat=disp_points_only[i].lat,
                                                    dE_obs=xdisp[i], dN_obs=ydisp[i], dU_obs=zdisp[i], Se_obs=np.nan,
                                                    Sn_obs=np.nan, Su_obs=np.nan, name="");
        modeled_disp_points.append(modeled_disp_point);
    print("Reading file %s... %d points" % (output_filename, len(modeled_disp_points)));
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
    # 1000 is parameter vmult.  Value of 1000 is arbitrarily high, representing the high-viscosity limit.
    ofile.write("%d \n" % len(fault_dict_list) );
    for fault in fault_dict_list:
        line = write_fault_slip_line_static1d_visco1d(fault);
        ofile.write(line);
    ofile.write('1\n');  # for velocities instead of displacements
    ofile.write('0\n');  # evaluate at depth specified in earlier runs of green's functions
    ofile.close();
    return;


def read_stat2C_geometry(infile):
    """
    Reading a fault geometry file used for stat2c.f, specifically for the Cascadia subduction zone.
    Returns a list of fault dictionaries in the internal format.
    """
    faultlist = [];
    with open(infile, 'r') as ifile:
        for line in ifile:
            temp = line.split();
            if len(temp) == 2:
                upper_depth, lower_depth = float(temp[0]), float(temp[1]);
            elif len(temp) == 7:
                lower_lat_corner = float(temp[0]);  # in degrees
                lower_lon_corner = float(temp[1]);  # in degrees
                length = float(temp[2]);  # in km
                strike = float(temp[3]);  # in degrees
                rake = float(temp[4]);  # in degrees
                dip = float(temp[6]);  # in degrees
                slip = float(temp[5])/100;  # from cm/yr into m/yr
                downdip_width = fault_vector_functions.get_downdip_width(upper_depth, lower_depth, dip);

                vector_mag = downdip_width * np.cos(np.deg2rad(dip));  # how far the bottom edge is displaced
                upper_corner_along_strike = fault_vector_functions.add_vector_to_point(0, 0, vector_mag, strike - 90);
                upper_corner_back_edge = fault_vector_functions.add_vector_to_point(upper_corner_along_strike[0],
                                                                                    upper_corner_along_strike[1],
                                                                                    length, strike+180);
                fault_lon, fault_lat = fault_vector_functions.xy2lonlat_single(upper_corner_back_edge[0],
                                                                               upper_corner_back_edge[1],
                                                                               lower_lon_corner,
                                                                               lower_lat_corner);

                new_fault = {"strike": strike, "dip": dip, "length": length, "rake": rake, "slip": slip, "tensile": 0,
                             "depth": upper_depth, "width": downdip_width, "lon": fault_lon, "lat": fault_lat};
                faultlist.append(new_fault);

    return faultlist;
