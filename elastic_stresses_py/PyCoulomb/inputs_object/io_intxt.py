"""
Read/write convenient input files (.intxt and .inzero).
Since we're using an input format that involves lon/lat/depth/mag, we have to convert
using Wells and Coppersmith 1994, in the same way that Coulomb does.
"""

import numpy as np
from . import io_mt
from .. import coulomb_collections as cc
from .. import conversion_math
from .input_obj import Input_object
from ..pyc_fault_object import Faults_object
from Tectonic_Utils.seismo import wells_and_coppersmith
from Tectonic_Utils.geodesy import fault_vector_functions


def read_intxt(input_file, mu=30e9, _lame1=30e9):
    """
    :param input_file: string, filename
    :param mu: optional, mu for focal mechanism sources
    :param _lame1: optional, lame1 never used
    :return: an input object
    """
    print("Reading source and receiver fault information from file %s " % input_file)
    sources, receivers = [], []
    receiver_horiz_profile = None
    [PR1, FRIC, minlon, maxlon, zerolon, minlat, maxlat, zerolat] = get_general_compute_params(input_file)
    [start_x, end_x, start_y, end_y, xinc, yinc] = compute_grid_params_general(minlon, maxlon, minlat, maxlat,
                                                                               zerolon, zerolat)
    ifile = open(input_file, 'r')
    for line in ifile:
        temp = line.split()
        if len(temp) > 0:
            if temp[0] == 'Source_WC:':  # wells and coppersmith convenient format
                one_source_object = get_source_wc(line, zerolon, zerolat)
                sources.append(one_source_object)
            if temp[0] == 'Source_Patch:':  # source-slip convenient format
                one_source_object = get_source_patch(line, zerolon, zerolat)
                sources.append(one_source_object)
            if temp[0] == 'Source_FM:':  # point source from focal mechanism
                one_source_object = get_FocalMech_source(line, zerolon, zerolat, mu)
                sources.append(one_source_object)
            if temp[0] == 'Receiver:':  # receiver fault
                one_receiver_object = get_receiver_fault(line, zerolon, zerolat)
                receivers.append(one_receiver_object)
            if temp[0] == 'Receiver_Horizontal_Profile:':
                receiver_horiz_profile = get_receiver_profile(line)
            if temp[0] == 'Source_Mogi:':  # Mogi Source
                one_mogi_source = get_mogi_source(line, zerolon, zerolat)
                sources.append(one_mogi_source)
    ifile.close()

    # Wrapping up the inputs.
    input_obj = Input_object(PR1=PR1, FRIC=FRIC, depth=0, start_gridx=start_x, finish_gridx=end_x,
                             start_gridy=start_y, finish_gridy=end_y, xinc=xinc, yinc=yinc, minlon=minlon,
                             maxlon=maxlon, zerolon=zerolon, minlat=minlat, maxlat=maxlat, zerolat=zerolat,
                             receiver_object=receivers, source_object=sources,
                             receiver_horiz_profile=receiver_horiz_profile)
    return input_obj


def write_intxt(input_object, output_file, label=None, mu=30e9, lame1=30e9):
    """
    :param input_object: cc.Input_object
    :param output_file: string, filename
    :param label: string, optional header
    :param mu: float, shear modulus
    :param lame1: float, first lame parameter
    """
    ofile = open(output_file, 'w')
    if label:
        ofile.write(label)
    ofile.write("# General: poissons_ratio friction_coef lon_min lon_max lon_zero lat_min lat_max lat_zero\n")
    ofile.write("# Source_Patch: strike rake dip length_km width_km lon lat depth_km slip_m\n")
    ofile.write("# Source_FM: strike rake dip lon lat depth_km magnitude\n")
    ofile.write("# Receiver: strike rake dip length_km width_km lon lat depth_km\n\n")
    real_poissons_ratio, _ = conversion_math.get_poissons_ratio_and_alpha(mu, lame1)

    # WRITE POISSON'S RATIO FROM MU, LAME1
    ofile.write("General: %.3f %.3f %f %f %f %f %f %f \n" % (real_poissons_ratio, input_object.FRIC,
                                                             input_object.minlon, input_object.maxlon,
                                                             input_object.zerolon, input_object.minlat,
                                                             input_object.maxlat, input_object.zerolat))
    for src in input_object.source_object:
        if isinstance(src, cc.Mogi_Source):
            src_lon, src_lat = fault_vector_functions.xy2lonlat(src.xstart, src.ystart, src.zerolon, src.zerolat)
            ofile.write("Source_Mogi: %.3f %.3f %.3f %.3f\n" % (src_lon, src_lat, src.depth, src.dV))
        if isinstance(src, cc.Faults_object):
            if not src.potency:  # write a finite source as a Source_Patch
                L = fault_vector_functions.get_strike_length(src.xstart, src.xfinish, src.ystart, src.yfinish)  # in km
                W = fault_vector_functions.get_downdip_width(src.top, src.bottom, src.dipangle)  # in km
                fault_lon, fault_lat = fault_vector_functions.xy2lonlat(src.xstart, src.ystart,
                                                                        src.zerolon, src.zerolat)
                slip = fault_vector_functions.get_vector_magnitude([src.rtlat, src.reverse])  # in m
                ofile.write("Source_Patch: %.3f %.3f %.3f %f %f %f %f %.3f %f\n" % (src.strike, src.rake, src.dipangle,
                                                                                    L, W, fault_lon, fault_lat,
                                                                                    src.top, slip))
            if src.potency:  # write a DC focal mechanism
                fault_lon, fault_lat = fault_vector_functions.xy2lonlat(src.xstart, src.ystart,
                                                                        src.zerolon, src.zerolat)
                mag = io_mt.get_mag_from_dc_potency(src.potency, mu, src.rake)
                ofile.write("Source_FM: %.3f %.3f %.3f %f %f %.3f %.3f\n" % (src.strike, src.rake, src.dipangle,
                                                                             fault_lon, fault_lat, src.top, mag))
    for rec in input_object.receiver_object:
        L = fault_vector_functions.get_strike_length(rec.xstart, rec.xfinish, rec.ystart, rec.yfinish)  # in km
        W = fault_vector_functions.get_downdip_width(rec.top, rec.bottom, rec.dipangle)  # in km
        fault_lon, fault_lat = fault_vector_functions.xy2lonlat(rec.xstart, rec.ystart, rec.zerolon, rec.zerolat)
        ofile.write("Receiver: %.3f %.3f %.3f %f %f %f %f %.3f \n" % (rec.strike, rec.rake, rec.dipangle, L, W,
                                                                      fault_lon, fault_lat, rec.top))
    if input_object.receiver_horiz_profile:
        rec = input_object.receiver_horiz_profile
        ofile.write("Receiver_Horizontal_Profile: ")
        ofile.write("%.3f %.3f %.3f %.3f " % (rec.depth_km, rec.strike, rec.dip, rec.rake))
        ofile.write("%f %f " % (rec.centerlon, rec.centerlat))
        ofile.write("%f %f %f\n" % (rec.length, rec.width, rec.inc))
    ofile.close()
    print("Writing file %s " % output_file)
    return


# ------------------------------------
# Functions to parse lines of text into objects
# ------------------------------------

def get_general_compute_params(input_file):
    [PR1, FRIC, minlon, maxlon, zerolon, minlat, maxlat, zerolat] = None, None, None, None, None, None, None, None
    ifile = open(input_file, 'r')
    for line in ifile:
        if len(line.split()) == 0:
            continue
        if line.split()[0] == 'General:':  # line of general parameters
            [PR1, FRIC, minlon, maxlon, zerolon, minlat, maxlat, zerolat] = read_general_line(line)
    ifile.close()
    assert (FRIC is not None), RuntimeError("Line General: general parameters not read from input file.")
    assert (PR1 is not None), RuntimeError("Line General: general parameters not read from input file.")
    assert (zerolat is not None), RuntimeError("Line General: general parameters not read from input file.")
    return [PR1, FRIC, minlon, maxlon, zerolon, minlat, maxlat, zerolat]


def get_receiver_fault(line, zerolon, zerolat):
    """
    Create a receiver object from line in text file
    """
    [strike, rake, dip, L, W, fault_lon, fault_lat, fault_depth] = read_receiver_line(line)
    [xstart, xfinish, ystart, yfinish, _, _, top, bottom, comment] = compute_params_for_slip_source(strike, dip, rake,
                                                                                                    fault_depth, L, W,
                                                                                                    fault_lon,
                                                                                                    fault_lat, 0,
                                                                                                    zerolon,
                                                                                                    zerolat)
    one_receiver_object = Faults_object(xstart=xstart, xfinish=xfinish, ystart=ystart, yfinish=yfinish, rtlat=0,
                                        reverse=0, tensile=0, strike=strike, dipangle=dip, zerolon=zerolon,
                                        zerolat=zerolat, rake=rake, top=top, bottom=bottom, comment=comment)
    defensive_programming_faults(one_receiver_object)
    return one_receiver_object


def get_source_patch(line, zerolon, zerolat):
    """
    Create a source object from a patch in text file
    """
    [strike, rake, dip, L, W, slip, tensile, flt_lon, flt_lat, fault_depth] = read_source_line_slip_convention(line)
    [xstart, xfinish, ystart, yfinish, rtlat, reverse, top, bottom,
     comment] = compute_params_for_slip_source(strike, dip, rake, fault_depth, L, W, flt_lon, flt_lat, slip, zerolon,
                                               zerolat)
    one_source_object = Faults_object(xstart=xstart, xfinish=xfinish, ystart=ystart, yfinish=yfinish, rtlat=rtlat,
                                      reverse=reverse, tensile=tensile, strike=strike, zerolon=zerolon, zerolat=zerolat,
                                      dipangle=dip, rake=rake, top=top, bottom=bottom, comment=comment)
    defensive_programming_faults(one_source_object)
    return one_source_object


def get_source_wc(line, zerolon, zerolat):
    """
    Create a source object from wells and coppersmith and text file
    """
    [strike, rake, dip, magnitude, faulting_type, fault_lon, fault_lat,
     fault_depth] = read_source_line_WCconvention(line)
    [xstart, xfinish, ystart, yfinish, rtlat, reverse, top, bottom,
     comment] = compute_params_for_WC_source(strike, dip, rake, fault_depth, magnitude, faulting_type, fault_lon,
                                             fault_lat, zerolon, zerolat)
    one_source_object = Faults_object(xstart=xstart, xfinish=xfinish, ystart=ystart, yfinish=yfinish, rtlat=rtlat,
                                      reverse=reverse, strike=strike, zerolon=zerolon, zerolat=zerolat, dipangle=dip,
                                      rake=rake, top=top, bottom=bottom, comment=comment)
    defensive_programming_faults(one_source_object)
    return one_source_object


def get_FocalMech_source(line, zerolon, zerolat, mu):
    """
    Create a source object from a point source focal mechanism
    """
    [strike, rake, dip, lon, lat, depth, magnitude] = read_point_source_line(line)
    [x, y, potency, comment] = io_mt.compute_params_for_point_source(rake, magnitude, lon, lat, zerolon, zerolat, mu)
    one_source_object = Faults_object(xstart=x, xfinish=x, ystart=y, yfinish=y, rtlat=0, reverse=0, potency=potency,
                                      strike=strike, dipangle=dip, zerolon=zerolon, zerolat=zerolat, rake=rake,
                                      top=depth, bottom=depth, comment=comment)
    defensive_programming_faults(one_source_object)
    return one_source_object


def get_mogi_source(line, zerolon, zerolat):
    [lon, lat, depth, dV] = read_mogi_source_line(line)
    [x_source, y_source] = fault_vector_functions.latlon2xy(lon, lat, zerolon, zerolat)
    one_mogi_source = cc.Mogi_Source(xstart=x_source, ystart=y_source,
                                     zerolon=zerolon, zerolat=zerolat, depth=depth, dV=dV)
    return one_mogi_source


def get_receiver_profile(line):
    """
    Create a horizontal profile based on the provided params.
    'Receiver_Horizontal_Profile: depth strike dip rake centerlon centerlat length_km width_km inc_km'
    """
    [depth, strike, dip, rake, centerlon, centerlat, length, width, inc] = [float(i) for i in line.split()[1:10]]
    startlon = fault_vector_functions.xy2lonlat(-length, 0, centerlon, centerlat)[0]
    endlon = fault_vector_functions.xy2lonlat(length, 0, centerlon, centerlat)[0]
    startlat = fault_vector_functions.xy2lonlat(0, -width, centerlon, centerlat)[1]
    endlat = fault_vector_functions.xy2lonlat(0, width, centerlon, centerlat)[1]
    nlonpts = (2.0 * length + 1.e-4) / inc
    nlatpts = (2.0 * width + 1.e-4) / inc
    lon_inc = (endlon - startlon) / float(nlonpts)
    lat_inc = (endlat - startlat) / float(nlatpts)
    lonpts = np.arange(startlon, endlon, lon_inc)
    latpts = np.arange(startlat, endlat, lat_inc)
    X, Y = np.meshgrid(lonpts, latpts)
    lon1d = np.reshape(X, (np.size(X),))
    lat1d = np.reshape(Y, (np.size(X),))
    horiz_profile = cc.Receiver_Horiz_Profile(depth_km=depth, strike=strike, dip=dip, rake=rake, centerlon=centerlon,
                                              centerlat=centerlat, lon1d=lon1d, lat1d=lat1d, length=length, width=width,
                                              inc=inc, shape=np.shape(X))
    return horiz_profile


def defensive_programming_faults(onefault):
    """
    Assert sanity on faults to make sure they've been read properly
    """
    critical_distance = 1000  # in km
    assert (-360 < onefault.strike < 360), RuntimeError("Invalid strike parameter for fault object")
    assert (-360 < onefault.rake < 360), RuntimeError("Invalid rake parameter for fault object")
    assert (0 < onefault.dipangle < 90), RuntimeError("Invalid dip parameter for fault object")
    assert (abs(onefault.xstart) < critical_distance), RuntimeError("Too distant fault for fault object")
    assert (abs(onefault.xfinish) < critical_distance), RuntimeError("Too distant fault for fault object")
    assert (abs(onefault.ystart) < critical_distance), RuntimeError("Too distant fault for fault object")
    assert (abs(onefault.yfinish) < critical_distance), RuntimeError("Too distant fault for fault object")
    assert (onefault.top >= 0), RuntimeError("Invalid depth parameter for fault object")
    assert (onefault.bottom >= onefault.top), RuntimeError("Bottom of fault above top of fault. ")
    if onefault.rtlat > 0 or onefault.reverse > 0:
        print("RtLat slip: %f m, Reverse slip: %f m" % (onefault.rtlat, onefault.reverse))
    return


# ------------------------------------
# Helper functions to parse lines of text into variables
# ------------------------------------

def read_source_line_WCconvention(line):
    """Format: strike rake dip magnitude faulting_type lon lat depth_km"""
    [strike, rake, dip, magnitude] = [float(i) for i in line.split()[1:5]]
    faulting_type = line.split()[5]
    [fault_center_lon, fault_center_lat, fault_center_dep] = [float(i) for i in line.split()[6:9]]
    return [strike, rake, dip, magnitude, faulting_type, fault_center_lon, fault_center_lat, fault_center_dep]


def read_source_line_slip_convention(line):
    """Format: strike rake dip length_km width_km lon lat depth_km slip_m"""
    line = line.split('#')[0]
    [strike, rake, dip, length, width] = [float(i) for i in line.split()[1:6]]
    [updip_corner_lon, updip_corner_lat, updip_corner_dep] = [float(i) for i in line.split()[6:9]]
    slip = float(line.split()[9])
    if len(line.split()) > 10:
        tensile = float(line.split()[10])
    else:
        tensile = 0
    return [strike, rake, dip, length, width, slip, tensile, updip_corner_lon, updip_corner_lat, updip_corner_dep]


def read_receiver_line(line):
    """Format: strike rake dip length_km width_km lon lat depth_km"""
    [strike, rake, dip, length, width] = [float(i) for i in line.split()[1:6]]
    [updip_corner_lon, updip_corner_lat, updip_corner_dep] = [float(i) for i in line.split()[6:9]]
    return [strike, rake, dip, length, width, updip_corner_lon, updip_corner_lat, updip_corner_dep]


def read_general_line(line):
    """Format: poissons_ratio friction_coef lon_min lon_max lon_zero lat_min lat_max lat_zero"""
    [PR1, FRIC, lon_min, lon_max, lon_zero] = [float(i) for i in line.split()[1:6]]
    [lat_min, lat_max, lat_zero] = [float(i) for i in line.split()[6:9]]
    return [PR1, FRIC, lon_min, lon_max, lon_zero, lat_min, lat_max, lat_zero]


def read_point_source_line(line):
    """Format: strike rake dip lon lat depth magnitude mu lamdba """
    [strike, rake, dip, lon, lat, depth, magnitude] = [float(i) for i in line.split()[1:8]]
    return [strike, rake, dip, lon, lat, depth, magnitude]


def read_mogi_source_line(line):
    """Format: lon lat depth dV"""
    [lon, lat, depth_km, dV] = [float(i) for i in line.split()[1:5]]
    return [lon, lat, depth_km, dV]


# ------------------------------------
# Compute functions to generate fault object's quantities
# ------------------------------------

def compute_grid_params_general(minlon, maxlon, minlat, maxlat, zerolon, zerolat):
    """
    Compute the Cargesian grid parameters (in km) that we'll be using, based on the size of the map.
    For synthetic displacements, assuming 100 increments in both directions.
    """
    deltalon = (maxlon - minlon) * 111.00 * np.cos(np.deg2rad(zerolat))  # in km.
    deltalat = (maxlat - minlat) * 111.00  # in km.
    start_gridx = (minlon - zerolon) * 111.00 * np.cos(np.deg2rad(zerolat))  # in km
    finish_gridx = (maxlon - zerolon) * 111.00 * np.cos(np.deg2rad(zerolat))  # in km.
    start_gridy = (minlat - zerolat) * 111.00  # in km.
    finish_gridy = (maxlat - zerolat) * 111.00  # in km.
    xinc = deltalon / 100.0
    yinc = deltalat / 100.0
    return [start_gridx, finish_gridx, start_gridy, finish_gridy, xinc, yinc]


def compute_params_for_WC_source(strike, dip, rake, depth, mag, faulting_type, fault_lon, fault_lat, zerolon, zerolat,
                                 refpoint='updip_corner'):
    """
    Helper function for setting up faults from 'Source_WC' format
    """
    [xcenter, ycenter] = fault_vector_functions.latlon2xy(fault_lon, fault_lat, zerolon, zerolat)
    L = wells_and_coppersmith.RLD_from_M(mag, faulting_type)  # rupture length
    W = wells_and_coppersmith.RW_from_M(mag, faulting_type)  # rupture width
    slip = wells_and_coppersmith.rectangular_slip(L * 1000, W * 1000, mag)  # must input in meters
    if refpoint == 'center':  # if hypocenter is really the center of the rupture:
        xstart, ystart = fault_vector_functions.add_vector_to_point(xcenter, ycenter, -L / 2, strike)
        top, bottom = fault_vector_functions.get_top_bottom_from_center(depth, W, dip)
        xfinish, yfinish = fault_vector_functions.add_vector_to_point(xcenter, ycenter, L / 2, strike)
    else:  # if hypocenter is top updip corner of rupture:
        xstart, ystart = fault_vector_functions.add_vector_to_point(xcenter, ycenter, 0, strike)
        xfinish, yfinish = fault_vector_functions.add_vector_to_point(xcenter, ycenter, L, strike)
        top, bottom = fault_vector_functions.get_top_bottom_from_top(depth, W, dip)
    rtlat, reverse = fault_vector_functions.get_rtlat_dip_slip(slip, rake)
    comment = ''
    return [xstart, xfinish, ystart, yfinish, rtlat, reverse, top, bottom, comment]


def compute_params_for_slip_source(strike, dip, rake, depth, L, W, fault_lon, fault_lat, slip, zerolon, zerolat,
                                   refpoint='updip_corner'):
    """
    Helper function for setting up faults from 'Source_Patch' and 'Receiver' format
    """
    [xcorner, ycorner] = fault_vector_functions.latlon2xy(fault_lon, fault_lat, zerolon, zerolat)
    if refpoint == 'center':  # if hypocenter is the center of the rupture:
        xstart, ystart = fault_vector_functions.add_vector_to_point(xcorner, ycorner, -L / 2, strike)
        xfinish, yfinish = fault_vector_functions.add_vector_to_point(xcorner, ycorner, L / 2, strike)
        top, bottom = fault_vector_functions.get_top_bottom_from_center(depth, W, dip)
    else:  # if hypocenter is top updip corner of rupture:
        xstart, ystart = fault_vector_functions.add_vector_to_point(xcorner, ycorner, 0, strike)
        xfinish, yfinish = fault_vector_functions.add_vector_to_point(xcorner, ycorner, L, strike)
        top, bottom = fault_vector_functions.get_top_bottom_from_top(depth, W, dip)
    rtlat, reverse = fault_vector_functions.get_rtlat_dip_slip(slip, rake)
    comment = ''
    return [xstart, xfinish, ystart, yfinish, rtlat, reverse, top, bottom, comment]
