"""
Read other fault formats of rectangular patches into fault slip objects

"""
import numpy as np
import scipy.io
from Tectonic_Utils.geodesy import fault_vector_functions
from Elastic_stresses_py.PyCoulomb.fault_slip_object import fault_slip_object
from . import io_four_corners
import json


def io_hamling_2017(filename):
    """Read Ian Hamling's 2017 Science paper fault model"""
    print("Reading file %s " % filename);
    fault_object_list = [];
    [centerlon, centerlat, strike, dip, rake, slip, l_km, top_km, bottom_km, segment] = np.loadtxt(filename,
                                                                                                   unpack=True,
                                                                                                   skiprows=8);
    for i in range(len(centerlon)):
        width = fault_vector_functions.get_downdip_width(top_km[i], bottom_km[i], dip[i]);
        x_start, y_start = fault_vector_functions.add_vector_to_point(0, 0, l_km[i] / 2, strike[i] - 180);  # in km
        downdip_width_proj = width * np.cos(np.deg2rad(dip[i]));
        x_start, y_start = fault_vector_functions.add_vector_to_point(x_start, y_start, downdip_width_proj/2,
                                                                      strike[i] - 90);
        corner_lon, corner_lat = fault_vector_functions.xy2lonlat(x_start, y_start, centerlon[i], centerlat[i]);
        one_fault = fault_slip_object.FaultSlipObject(strike=strike[i], dip=dip[i], length=l_km[i], width=width,
                                                      depth=top_km[i], rake=rake[i], slip=slip[i], tensile=0,
                                                      segment=int(segment[i]), lon=corner_lon, lat=corner_lat);
        fault_object_list.append(one_fault);
    print("--> Returning %d fault patches " % len(fault_object_list));
    return fault_object_list;


def io_wallace_sse(filename):
    """
    Read Laura Wallace's slow slip event files for New Zealand
    """
    print("Reading file %s " % filename);
    fault_object_list = [];
    ifile = open(filename, 'r');
    for line in ifile:
        if line[0] == "#":
            continue;
        if line[0:4] == "> -Z":
            temp = line.split();
            patch_slip_m = float(temp[3])/1000;  # into meters
            patch_tensile_m = float(temp[6])/1000;
            rake = float(temp[7]);
            c1_lon, c1_lat, c1_depth = ifile.readline().split();
            c2_lon, c2_lat, c2_depth = ifile.readline().split();
            c3_lon, c3_lat, c3_depth = ifile.readline().split();
            c4_lon, c4_lat, c4_depth = ifile.readline().split();
            lons = [float(c1_lon), float(c2_lon), float(c3_lon), float(c4_lon)];
            lats = [float(c1_lat), float(c2_lat), float(c3_lat), float(c4_lat)];
            depths = [-float(c1_depth), -float(c2_depth), -float(c3_depth), -float(c4_depth)];
            if np.sum(lons) == 0 and np.sum(lats) == 0 and np.sum(depths) == 0:  # skip pathological case
                continue;
            one_fault = io_four_corners.get_fault_object_from_four_corners(lons, lats, depths);
            new_fault = fault_slip_object.FaultSlipObject(strike=one_fault.strike, dip=one_fault.dip,
                                                          lon=one_fault.lon, lat=one_fault.lat,
                                                          depth=one_fault.depth, length=one_fault.length,
                                                          width=one_fault.width, segment=one_fault.segment,
                                                          slip=patch_slip_m, rake=rake, tensile=patch_tensile_m);
            fault_object_list.append(new_fault);
    ifile.close();
    print("--> Returning %d fault patches " % len(fault_object_list));
    return fault_object_list;


def io_dreger_finite_fault(filename):
    """
    Read Doug Dreger's .mat format for point-source faults
    keys: 'dep', 'dip', 'lat', 'lon', 'rake', 'slip (cm)', 'str'.  Each is 2D array of 1km x 1km patches.
    """
    print("Reading file %s " % filename);
    fault_object_list = [];
    mat_data = scipy.io.loadmat(filename);  # loads into a dictionary.  Can print keys() to see it.
    for i in range(len(mat_data['dep'])):  # for each row in the slip distribution (constant depth)
        for j in range(len(mat_data['dep'][i])):  # for each elemnet inside the row
            new_fault = fault_slip_object.FaultSlipObject(strike=mat_data['str'][i][j], dip=mat_data['dip'][i][j],
                                                          rake=mat_data['rake'][i][j], depth=mat_data['dep'][i][j],
                                                          lon=mat_data['lon'][i][j], lat=mat_data['lat'][i][j],
                                                          slip=mat_data['slip'][i][j]*0.01,
                                                          length=1, width=1, segment=0, tensile=0);
            fault_object_list.append(new_fault);
    print("--> Returning %d fault patches " % len(fault_object_list));
    return fault_object_list;


def io_iceland_met_office(filename):
    """
    Read a Json that contains the iceland met office fault discritization
    """
    fault_object_list = [];
    f = open(filename);
    print("Reading file %s " % filename);
    read_object = json.load(f)['fault1'];
    X = read_object['X']/1000;  # x-coordinate in km
    Y = read_object['Y']/1000   # y-coordinate in km
    reflon = read_object['reflon']
    reflat = read_object['reflat']
    strike_slip, dip_slip = read_object['strike_slip'], read_object['dip_slip'];
    new_fault = fault_slip_object.FaultSlipObject(strike=read_object['strike'], dip=-read_object['dip'],
                                                  depth=read_object['depth']/1000,
                                                  lon=fault_vector_functions.xy2lonlat_single(X, Y, reflon, reflat)[0],
                                                  lat=fault_vector_functions.xy2lonlat_single(X, Y, reflon, reflat)[1],
                                                  slip=fault_vector_functions.get_total_slip(strike_slip, dip_slip),
                                                  rake=fault_vector_functions.get_rake(strike_slip, dip_slip),
                                                  length=read_object['length']/1000, width=read_object['width']/1000,
                                                  segment=0, tensile=0);
    fault_object_list.append(new_fault);
    print("--> Returning %d fault patches " % len(fault_object_list));
    return fault_object_list;
