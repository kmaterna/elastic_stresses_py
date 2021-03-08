# The purpose of these functions is to read/write a convenient input file
# Since we're using an input format that involves lon/lat/depth/mag, we have to convert
# We will use Wells and Coppersmith 1994, in the same way that Coulomb does.

import numpy as np
from . import coulomb_collections as cc
from Tectonic_Utils.seismo import wells_and_coppersmith
from Tectonic_Utils.seismo import moment_calculations
from Tectonic_Utils.geodesy import fault_vector_functions


def read_intxt(input_file):
    print("Reading source and receiver fault information from file %s " % input_file);
    sources, receivers = [], [];
    [PR1, FRIC, minlon, maxlon, zerolon, minlat, maxlat, zerolat] = input_general_compute_params(input_file);
    ifile = open(input_file, 'r');
    for line in ifile:
        temp = line.split();
        if len(temp) > 0:
            if temp[0] == 'Source_WC:':  # wells and coppersmith convenient format
                one_source_object = get_source_wc(line, zerolon, zerolat);
                sources.append(one_source_object);
                print("RtLat slip: %f m, Reverse slip: %f m" % (one_source_object.rtlat, one_source_object.reverse));
            if temp[0] == 'Source_Patch:':  # source-slip convenient format
                one_source_object = get_source_patch(line, zerolon, zerolat);
                sources.append(one_source_object);
                print("RtLat slip: %f m, Reverse slip: %f m" % (one_source_object.rtlat, one_source_object.reverse));
            if temp[0] == 'Source_FM:':  # point source from focal mechanism
                one_source_object = get_FocalMech_source(line, zerolon, zerolat);
                sources.append(one_source_object);
            if temp[0] == "Source_MT":  # point source from moment tensor
                one_source_object = get_MT_source(line, zerolon, zerolat);
                sources.append(one_source_object);
            if temp[0] == 'Receiver:':
                one_receiver_object = get_receiver_fault(line, zerolon, zerolat);
                receivers.append(one_receiver_object);
    ifile.close();

    # Wrapping up the inputs.
    [start_gridx, end_gridx, start_gridy, end_gridy, xinc, yinc] = compute_grid_params_general(minlon, maxlon, minlat,
                                                                                               maxlat, zerolat);
    input_obj = cc.Input_object(PR1=PR1, FRIC=FRIC, depth=0, start_gridx=start_gridx, finish_gridx=end_gridx,
                                start_gridy=start_gridy, finish_gridy=end_gridy, xinc=xinc, yinc=yinc, minlon=minlon,
                                maxlon=maxlon, zerolon=zerolon, minlat=minlat, maxlat=maxlat, zerolat=zerolat,
                                receiver_object=receivers, source_object=sources);
    return input_obj;

def write_intxt(input_object, output_file):
    ofile = open(output_file, 'w');
    ofile.write("# General: poissons_ratio friction_coef lon_min lon_max lon_zero lat_min lat_max lat_zero\n");
    ofile.write("# Source_Patch: strike rake dip length_km width_km lon lat depth_km slip_m\n");
    ofile.write("# Receiver: strike rake dip length_km width_km lon lat depth_km\n\n");
    ofile.write("General: %f %f %f %f %f %f %f %f \n" % (input_object.PR1, input_object.FRIC, input_object.minlon,
                                                         input_object.maxlon, input_object.zerolon, input_object.minlat,
                                                         input_object.maxlat, input_object.zerolat) );
    for src in input_object.source_object:
        if not src.potency:  # write a finite source
            L = fault_vector_functions.get_strike_length(src.xstart, src.xfinish, src.ystart, src.yfinish);  # in km
            W = fault_vector_functions.get_downdip_width(src.top, src.bottom, src.dipangle);  # in km
            fault_lon, fault_lat = fault_vector_functions.xy2lonlat(src.xstart, src.ystart, input_object.zerolon,
                                                                    input_object.zerolat);
            slip = fault_vector_functions.get_vector_magnitude([src.rtlat, src.reverse]);  # in m
            ofile.write("Source_Patch: %f %f %f %f %f %f %f %f %f\n" % (src.strike, src.rake, src.dipangle, L, W,
                                                                        fault_lon, fault_lat, src.top, slip));
        if src.potency:   # write a focal mechanism source
            continue;   # still working on this.
    for rec in input_object.receiver_object:
        L = fault_vector_functions.get_strike_length(rec.xstart, rec.xfinish, rec.ystart, rec.yfinish);  # in km
        W = fault_vector_functions.get_downdip_width(rec.top, rec.bottom, rec.dipangle);  # in km
        fault_lon, fault_lat = fault_vector_functions.xy2lonlat(rec.xstart, rec.ystart, input_object.zerolon,
                                                                input_object.zerolat);
        ofile.write("Receiver: %f %f %f %f %f %f %f %f \n" % (rec.strike, rec.rake, rec.dipangle, L, W, fault_lon,
                                                              fault_lat, rec.top) );
    ofile.close();
    print("Writing file %s " % output_file);
    return;


# ------------------------------------
# Functions to parse lines of text into objects
# ------------------------------------

def input_general_compute_params(input_file):
    [PR1, FRIC, minlon, maxlon, zerolon, minlat, maxlat, zerolat] = None, None, None, None, None, None, None, None;
    ifile = open(input_file, 'r');
    for line in ifile:
        if len(line.split()) == 0:
            continue;
        if line.split()[0] == 'General:':  # line of general parameters
            [PR1, FRIC, minlon, maxlon, zerolon, minlat, maxlat, zerolat] = read_general_line(line);
    ifile.close();
    assert (FRIC is not None), RuntimeError("Line General: general parameters not read from input file.");
    assert (PR1 is not None), RuntimeError("Line General: general parameters not read from input file.");
    assert (zerolat is not None), RuntimeError("Line General: general parameters not read from input file.");
    return [PR1, FRIC, minlon, maxlon, zerolon, minlat, maxlat, zerolat];

def get_receiver_fault(line, zerolon, zerolat):
    """Create a receiver object from line in text file"""
    [strike, rake, dip, L, W, fault_lon, fault_lat, fault_depth] = read_receiver_line(line);
    [xstart, xfinish, ystart, yfinish, Kode, rtlat, reverse, top, bottom,
     comment] = compute_params_for_slip_source(strike, dip, rake, fault_depth, L, W, fault_lon, fault_lat, 0,
                                               zerolon, zerolat);
    one_receiver_object = cc.Faults_object(xstart=xstart, xfinish=xfinish, ystart=ystart, yfinish=yfinish, Kode=Kode,
                                           rtlat=rtlat, reverse=reverse, potency=[], strike=strike, dipangle=dip,
                                           rake=rake, top=top, bottom=bottom, comment=comment);
    return one_receiver_object;

def get_source_patch(line, zerolon, zerolat):
    """Create a source object from a patch in text file"""
    [strike, rake, dip, L, W, slip, fault_lon, fault_lat, fault_depth] = read_source_line_slip_convention(line);
    [xstart, xfinish, ystart, yfinish, Kode, rtlat, reverse, top, bottom,
     comment] = compute_params_for_slip_source(strike, dip, rake, fault_depth, L, W, fault_lon, fault_lat, slip,
                                               zerolon, zerolat);
    one_source_object = cc.Faults_object(xstart=xstart, xfinish=xfinish, ystart=ystart, yfinish=yfinish,
                                         Kode=Kode, rtlat=rtlat, reverse=reverse, potency=[], strike=strike,
                                         dipangle=dip, rake=rake, top=top, bottom=bottom, comment=comment);
    return one_source_object;

def get_source_wc(line, zerolon, zerolat):
    """Create a source object from wells and coppersmith and text file"""
    [strike, rake, dip, magnitude, faulting_type, fault_lon, fault_lat,
     fault_depth] = read_source_line_WCconvention(line);
    [xstart, xfinish, ystart, yfinish, Kode, rtlat, reverse, top, bottom,
     comment] = compute_params_for_WC_source(strike, dip, rake, fault_depth, magnitude, faulting_type, fault_lon,
                                             fault_lat, zerolon, zerolat);
    one_source_object = cc.Faults_object(xstart=xstart, xfinish=xfinish, ystart=ystart, yfinish=yfinish,
                                         Kode=Kode, rtlat=rtlat, reverse=reverse, potency=[], strike=strike,
                                         dipangle=dip, rake=rake, top=top, bottom=bottom, comment=comment);
    return one_source_object;

def get_FocalMech_source(line, zerolon, zerolat):
    """Create a source object from a point source focal mechanism"""
    [strike, rake, dip, lon, lat, depth, magnitude, mu, lame1] = read_point_source_line(line);
    [x, y, Kode, rtlat, reverse, potency, comment] = compute_params_for_point_source(rake, magnitude, lon, lat,
                                                                                     zerolon, zerolat, mu,
                                                                                     lame1);
    one_source_object = cc.Faults_object(xstart=x, xfinish=x, ystart=y, yfinish=y, Kode=Kode, rtlat=rtlat,
                                         reverse=reverse, potency=potency, strike=strike, dipangle=dip,
                                         rake=rake, top=depth, bottom=depth, comment=comment);
    return one_source_object;

def get_MT_source(line, zerolon, zerolat):
    """Write this next. """
    print(line, zerolon, zerolat);
    _ = read_moment_tensor_source_line(line);
    # step 1: read line
    # step 2: compute necessary things
    # step 3: package into a source
    return 0;


# ------------------------------------
# Helper functions to parse lines of text into variables
# ------------------------------------

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
    return [strike, rake, dip, length, width, updip_corner_lon, updip_corner_lat, updip_corner_dep]

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

def read_point_source_line(line):
    """Format: strike rake dip lon lat depth magnitude mu lamdba """
    strike = float(line.split()[1]);
    rake = float(line.split()[2]);
    dip = float(line.split()[3]);
    lon = float(line.split()[4]);
    lat = float(line.split()[5]);
    depth = float(line.split()[6]);
    magnitude = float(line.split()[7]);
    mu = float(line.split()[8]);
    lame1 = float(line.split()[9]);
    return [strike, rake, dip, lon, lat, depth, magnitude, mu, lame1];

def read_moment_tensor_source_line(line):
    """write this next"""
    print(line);
    return 0;


# ------------------------------------
# Compute functions to generate fault object's quantities
# ------------------------------------

def compute_grid_params_general(minlon, maxlon, minlat, maxlat, zerolat):
    """
    Compute the grid parameters that we'll be using, based on the size of the map.
    Assuming 100 increments in both directions.
    """
    deltalon = (maxlon - minlon) * 111.00 * np.cos(np.deg2rad(zerolat));  # in km.
    deltalat = (maxlat - minlat) * 111.00;  # in km.
    start_gridx = -deltalon / 2.0;
    finish_gridx = deltalon / 2.0;
    start_gridy = -deltalat / 2.0;
    finish_gridy = deltalat / 2.0;
    xinc = deltalon / 100.0;
    yinc = deltalat / 100.0;
    return [start_gridx, finish_gridx, start_gridy, finish_gridy, xinc, yinc];

def compute_params_for_WC_source(strike, dip, rake, depth, mag, faulting_type, fault_lon, fault_lat, zerolon, zerolat):
    [xcenter, ycenter] = fault_vector_functions.latlon2xy(fault_lon, fault_lat, zerolon, zerolat);
    L = wells_and_coppersmith.RLD_from_M(mag, faulting_type);  # rupture length
    W = wells_and_coppersmith.RW_from_M(mag, faulting_type);  # rupture width
    slip = wells_and_coppersmith.rectangular_slip(L * 1000, W * 1000, mag);  # must input in meters
    # if the hypocenter is really the center of the rupture:
    # xstart,ystart=conversion_math.add_vector_to_point(xcenter,ycenter,-L/2,strike[i]);
    # xfinish,yfinish=conversion_math.add_vector_to_point(xcenter,ycenter,L/2,strike[i]);
    # assuming hypocenter is on one side of rupture:
    xstart, ystart = fault_vector_functions.add_vector_to_point(xcenter, ycenter, 0, strike);
    xfinish, yfinish = fault_vector_functions.add_vector_to_point(xcenter, ycenter, L, strike);
    rtlat, reverse = fault_vector_functions.get_rtlat_dip_slip(slip, rake);
    top, bottom = fault_vector_functions.get_top_bottom_from_top(depth, W, dip);
    Kode = 100;
    comment = '';
    return [xstart, xfinish, ystart, yfinish, Kode, rtlat, reverse, top, bottom, comment];

def compute_params_for_slip_source(strike, dip, rake, depth, L, W, fault_lon, fault_lat, slip, zerolon, zerolat):
    [xcorner, ycorner] = fault_vector_functions.latlon2xy(fault_lon, fault_lat, zerolon, zerolat);
    xstart, ystart = fault_vector_functions.add_vector_to_point(xcorner, ycorner, 0, strike);
    xfinish, yfinish = fault_vector_functions.add_vector_to_point(xcorner, ycorner, L, strike);
    rtlat, reverse = fault_vector_functions.get_rtlat_dip_slip(slip, rake);
    top, bottom = fault_vector_functions.get_top_bottom_from_top(depth, W, dip);
    Kode = 100;
    comment = '';
    return [xstart, xfinish, ystart, yfinish, Kode, rtlat, reverse, top, bottom, comment]

def compute_params_for_point_source(rake, magnitude, lon, lat, zerolon, zerolat, mu, lame1):
    """ Given information about point sources from focal mechanisms,
    Return the right components that get packaged into input_obj. """
    [xcenter, ycenter] = fault_vector_functions.latlon2xy(lon, lat, zerolon, zerolat);
    potency = get_DC_potency(rake, magnitude, mu, lame1);
    # Filler variables for the point source case
    Kode = 100;
    comment = '';
    return [xcenter, ycenter, Kode, 0, 0, potency, comment];

def get_DC_potency(rake, momentmagnitude, mu, lame1):
    """
    Given the basic double couple parameters,
    Return the four-vector used in Okada DC3D0.
    Pot1 = strike-slip moment of DC / mu
    Pot2 = dip-slip moment of DC / mu
    Pot3 = inflation = M_ISO / lambda
    Pot4 = tensile = M_lineardipole / mu
    In a more general case, we would use a different MT format to handle non-DC parts.
    Right now, it only handles DC focal mechanisms.
    Moment in newton meters
    """
    total_moment = moment_calculations.moment_from_mw(momentmagnitude);
    strike_slip_fraction, dip_slip_fraction = fault_vector_functions.get_rtlat_dip_slip(1.0, rake);
    print("strike_slip: ", strike_slip_fraction);
    print("dip_slip: ", dip_slip_fraction);
    strike_slip_fraction = -1 * strike_slip_fraction;  # DC3D0 wants left lateral slip.
    p1 = total_moment * strike_slip_fraction / mu;
    p2 = total_moment * dip_slip_fraction / mu;
    # In the double-couple case, this is zero.
    p3 = 0;
    p4 = 0;
    return [p1, p2, p3, p4];
