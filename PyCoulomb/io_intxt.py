# The purpose of these functions is to read/write a convenient input file
# Since we're using an input format that involves lon/lat/depth/mag, we have to convert
# We will use Wells and Coppersmith 1994, in the same way that Coulomb does. 
# This input file assumes a fixed rake indicated for sources and receivers

import numpy as np
from . import conversion_math
from . import coulomb_collections
from Tectonic_Utils.seismo import wells_and_coppersmith


def read_intxt(input_file):
    print("Reading source and receiver fault information from file %s " % input_file);
    sources = [];
    receivers = [];

    ifile = open(input_file, 'r');
    for line in ifile:
        temp = line.split();
        if len(temp) == 0:
            continue;
        if temp[0] == 'S:':
            # reading wells and coppersmith convenient format
            if " SS " in line or " N " in line or " R " in line or " ALL " in line:
                [strike, rake, dip, magnitude, faulting_type, fault_lon, fault_lat,
                 fault_depth] = read_source_line_WCconvention(line);
                [xstart, xfinish, ystart, yfinish, Kode, rtlat, reverse, top, bottom,
                 comment] = compute_params_for_WC_source(strike, dip, rake, fault_depth, magnitude, faulting_type,
                                                         fault_lon, fault_lat, zerolon, zerolat);
            else:  # reading the source-slip convenient format
                [strike, rake, dip, L, W, slip, fault_lon, fault_lat, fault_depth] = read_source_line_slip_convention(
                    line);
                [xstart, xfinish, ystart, yfinish, Kode, rtlat, reverse, top, bottom,
                 comment] = compute_params_for_slip_source(strike, dip, rake, fault_depth, L, W, fault_lon, fault_lat,
                                                           slip, zerolon, zerolat);
            one_source_object = coulomb_collections.Faults_object(xstart=xstart, xfinish=xfinish, ystart=ystart,
                                                                  yfinish=yfinish, Kode=Kode, rtlat=rtlat,
                                                                  reverse=reverse, potency=[], strike=strike,
                                                                  dipangle=dip, rake=rake, top=top, bottom=bottom,
                                                                  comment=comment);
            sources.append(one_source_object);
            print("RtLat slip: %f m, Reverse slip: %f m" % (rtlat, reverse));
        elif temp[0] == 'R:':
            [strike, rake, dip, L, W, fault_lon, fault_lat, fault_depth] = read_receiver_line(line);
            [xstart, xfinish, ystart, yfinish, Kode, rtlat, reverse, top, bottom,
             comment] = compute_params_for_slip_source(strike, dip, rake, fault_depth, L, W, fault_lon, fault_lat, 0,
                                                       zerolon, zerolat);
            one_receiver_object = coulomb_collections.Faults_object(xstart=xstart, xfinish=xfinish, ystart=ystart,
                                                                    yfinish=yfinish, Kode=Kode, rtlat=rtlat,
                                                                    reverse=reverse, potency=[], strike=strike,
                                                                    dipangle=dip, rake=rake, top=top, bottom=bottom,
                                                                    comment=comment);
            receivers.append(one_receiver_object);
        elif temp[0] == 'G:':
            [PR1, FRIC, minlon, maxlon, zerolon, minlat, maxlat, zerolat] = read_general_line(line);
        else:
            continue;
    ifile.close();

    # Wrapping up the inputs.
    [start_gridx, finish_gridx, start_gridy, finish_gridy, xinc, yinc] = compute_grid_parameters(minlon, maxlon,
                                                                                                 minlat, maxlat,
                                                                                                 zerolat);
    input_obj = coulomb_collections.Input_object(PR1=PR1, FRIC=FRIC, depth=0, start_gridx=start_gridx,
                                                 finish_gridx=finish_gridx, start_gridy=start_gridy,
                                                 finish_gridy=finish_gridy,
                                                 xinc=xinc, yinc=yinc, minlon=minlon, maxlon=maxlon, zerolon=zerolon,
                                                 minlat=minlat, maxlat=maxlat, zerolat=zerolat,
                                                 receiver_object=receivers, source_object=sources);
    return input_obj;


def read_source_line_WCconvention(line):
    strike = float(line.split()[1]);
    rake = float(line.split()[2]);
    dip = float(line.split()[3]);
    magnitude = float(line.split()[4]);
    faulting_type = line.split()[5];
    fault_center_lon = float(line.split()[6]);
    fault_center_lat = float(line.split()[7]);
    fault_center_dep = float(line.split()[8]);
    return [strike, rake, dip, magnitude, faulting_type, fault_center_lon, fault_center_lat, fault_center_dep];


def read_source_line_slip_convention(line):
    strike = float(line.split()[1]);
    rake = float(line.split()[2]);
    dip = float(line.split()[3]);
    length = float(line.split()[4]);
    width = float(line.split()[5]);
    updip_corner_lon = float(line.split()[6]);
    updip_corner_lat = float(line.split()[7]);
    updip_corner_dep = float(line.split()[8]);
    slip = float(line.split()[9]);
    return [strike, rake, dip, length, width, slip, updip_corner_lon, updip_corner_lat, updip_corner_dep];


def read_receiver_line(line):
    strike = float(line.split()[1]);
    rake = float(line.split()[2]);
    dip = float(line.split()[3]);
    length = float(line.split()[4]);
    width = float(line.split()[5]);
    updip_corner_lon = float(line.split()[6]);
    updip_corner_lat = float(line.split()[7]);
    updip_corner_dep = float(line.split()[8]);
    return [strike, rake, dip, length, width, updip_corner_lon, updip_corner_lat, updip_corner_dep];


def read_general_line(line):
    PR1 = float(line.split()[1]);
    FRIC = float(line.split()[2]);
    lon_min = float(line.split()[3]);
    lon_max = float(line.split()[4]);
    lon_zero = float(line.split()[5]);
    lat_min = float(line.split()[6]);
    lat_max = float(line.split()[7]);
    lat_zero = float(line.split()[8]);
    return [PR1, FRIC, lon_min, lon_max, lon_zero, lat_min, lat_max, lat_zero];


def compute_grid_parameters(minlon, maxlon, minlat, maxlat, zerolat):
    # Compute the grid parameters that we'll be using, based on the size of the map.
    # Assuming 100 increments in both directions.
    deltalon = (maxlon - minlon) * 111.00 * np.cos(np.deg2rad(zerolat));  # in km.
    deltalat = (maxlat - minlat) * 111.00;  # in km.
    start_gridx = -deltalon / 2.0;
    finish_gridx = deltalon / 2.0;
    start_gridy = -deltalat / 2.0;
    finish_gridy = deltalat / 2.0;
    xinc = deltalon / 100.0;
    yinc = deltalat / 100.0;
    return [start_gridx, finish_gridx, start_gridy, finish_gridy, xinc, yinc];


def compute_params_for_WC_source(strike, dip, rake, depth, magnitude, faulting_type, fault_lon, fault_lat, zerolon,
                                 zerolat):
    [xcenter, ycenter] = conversion_math.latlon2xy(fault_lon, fault_lat, zerolon, zerolat);
    L = wells_and_coppersmith.RLD_from_M(magnitude, faulting_type);  # rupture length
    W = wells_and_coppersmith.RW_from_M(magnitude, faulting_type);  # rupture width
    slip = wells_and_coppersmith.rectangular_slip(L * 1000, W * 1000, magnitude);  # must input in meters
    # xstart,ystart=conversion_math.add_vector_to_point(xcenter,ycenter,-L/2,strike[i]);  # if the hypocenter is really the center of the rupture
    # xfinish,yfinish=conversion_math.add_vector_to_point(xcenter,ycenter,L/2,strike[i]);
    xstart, ystart = conversion_math.add_vector_to_point(xcenter, ycenter, 0,
                                                         strike);  # if the hypocenter is on one side of the rupture
    xfinish, yfinish = conversion_math.add_vector_to_point(xcenter, ycenter, L, strike);
    rtlat, reverse = conversion_math.get_rtlat_dip_slip(slip, rake);
    top, bottom = conversion_math.get_top_bottom(depth, W, dip);
    Kode = 100;
    comment = '';
    return [xstart, xfinish, ystart, yfinish, Kode, rtlat, reverse, top, bottom, comment];


def compute_params_for_slip_source(strike, dip, rake, depth, L, W, fault_lon, fault_lat, slip, zerolon, zerolat):
    [xcorner, ycorner] = conversion_math.latlon2xy(fault_lon, fault_lat, zerolon, zerolat);
    xstart, ystart = conversion_math.add_vector_to_point(xcorner, ycorner, 0, strike);
    xfinish, yfinish = conversion_math.add_vector_to_point(xcorner, ycorner, L, strike);
    rtlat, reverse = conversion_math.get_rtlat_dip_slip(slip, rake);
    top, bottom = conversion_math.get_top_bottom_from_top(depth, W, dip);
    Kode = 100;
    comment = '';
    return [xstart, xfinish, ystart, yfinish, Kode, rtlat, reverse, top, bottom, comment];
