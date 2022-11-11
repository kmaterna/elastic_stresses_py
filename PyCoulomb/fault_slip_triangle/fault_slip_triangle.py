
from typing import TypedDict
import numpy as np
from Tectonic_Utils.geodesy import fault_vector_functions
from Tectonic_Utils.seismo import moment_calculations
from .. import conversion_math
from ..fault_slip_object.io_pycoulomb import fault_dict_to_coulomb_fault

class FaultDict(TypedDict):
    """
    The internal format is a dictionary for a triangular fault segment:
    {
        vertex1 [x, y, z] in m, z is positive downward
        vertex2 [x, y, z] in m, z is positive downward
        vertex3 [x, y, z] in m, z is positive downward
        lon(coord reference),
        lat(coord reference),
        depth (km),
        rt_lat slip(m),
        dip_slip(m),
        tensile(m),
        segment(int)
    }
    If the fault is a receiver fault, we put slip = 0
    """
    vertex1: np.ndarray
    vertex2: np.ndarray
    vertex3: np.ndarray
    lon: float
    lat: float
    depth: float
    rtlat_slip: float
    dip_slip: float
    tensile: float
    segment: int


def write_gmt_plots_cartesian(triangle_list, outfile):
    # Write triangle edges out to file for GMT plots, in cartesian space
    print("Writing %d triangles to file %s " % (len(triangle_list), outfile) );
    with open(outfile, 'w') as ofile:
        for item in triangle_list:
            total_slip = get_total_slip(item);
            ofile.write("> -Z%f \n" % total_slip);
            ofile.write("%f %f \n" % (item['vertex1'][0], item['vertex1'][1]));
            ofile.write("%f %f \n" % (item['vertex2'][0], item['vertex2'][1]));
            ofile.write("%f %f \n" % (item['vertex3'][0], item['vertex3'][1]));
            ofile.write("%f %f \n" % (item['vertex1'][0], item['vertex1'][1]));
    return;


def write_gmt_plots_geographic(triangle_list, outfile):
    # Write triangle edges out to file for GMT plots, in lon/lat space assuming a cartesian-to-geographic transform
    print("Writing %d triangles to file %s " % (len(triangle_list), outfile) );
    with open(outfile, 'w') as ofile:
        for item in triangle_list:
            total_slip = get_total_slip(item);
            vertex1, vertex2, vertex3 = get_ll_corners(item);
            ofile.write("> -Z%f \n" % total_slip);
            ofile.write("%f %f \n" % (vertex1[0], vertex1[1]));
            ofile.write("%f %f \n" % (vertex2[0], vertex2[1]));
            ofile.write("%f %f \n" % (vertex3[0], vertex3[1]));
            ofile.write("%f %f \n" % (vertex1[0], vertex1[1]));
    return;


def get_ll_corners(triangle_fault):
    # Convert m to km before coordinate conversion
    # Returns 3 tuples (lon, lat)
    lon1, lat1 = fault_vector_functions.xy2lonlat_single(triangle_fault['vertex1'][0] / 1000,
                                                         triangle_fault['vertex1'][1] / 1000,
                                                         triangle_fault['lon'], triangle_fault['lat']);
    lon2, lat2 = fault_vector_functions.xy2lonlat_single(triangle_fault['vertex2'][0] / 1000,
                                                         triangle_fault['vertex2'][1] / 1000,
                                                         triangle_fault['lon'], triangle_fault['lat']);
    lon3, lat3 = fault_vector_functions.xy2lonlat_single(triangle_fault['vertex3'][0] / 1000,
                                                         triangle_fault['vertex3'][1] / 1000,
                                                         triangle_fault['lon'], triangle_fault['lat']);
    vertex1 = [lon1, lat1]; vertex2 = [lon2, lat2]; vertex3 = [lon3, lat3];
    return vertex1, vertex2, vertex3;


def get_total_slip(triangle_fault):
    total_slip = np.sqrt(np.square(triangle_fault['dip_slip']) + np.square(triangle_fault['rtlat_slip']));
    return total_slip;


def check_consistent_reference_frame(triangle_fault_list):
    """Returns True if all the triangles have the coordinate reference system"""
    return_code = 1;  # default true
    reflon, reflat = triangle_fault_list[0]['lon'], triangle_fault_list[0]['lat'];
    for item in triangle_fault_list:
        if item['lon'] != reflon or item['lat'] != reflat:
            return_code = 0;
    if return_code == 0:
        print("Warning! Not all triangles have the same reference lon/lat");
    return return_code;


def convert_rectangle_into_two_triangles(one_fault_dict):
    """
    Convert a fault_dict into two triangular fault dicts. The fault normals are expected to point up.
    """
    [source] = fault_dict_to_coulomb_fault([one_fault_dict]);
    [x_all, y_all, _, _] = conversion_math.get_fault_four_corners(source);  # This is cartesian
    top_depth, bottom_depth = source.top, source.bottom;
    vertex1 = np.array([x_all[0]*1000, y_all[0]*1000, top_depth*1000]);  # in meters
    vertex2 = np.array([x_all[1]*1000, y_all[1]*1000, top_depth*1000]);
    vertex3 = np.array([x_all[2]*1000, y_all[2]*1000, bottom_depth*1000]);
    vertex4 = np.array([x_all[3]*1000, y_all[3]*1000, bottom_depth*1000]);
    first_triangle = FaultDict(lon=one_fault_dict['lon'], lat=one_fault_dict['lat'], segment=one_fault_dict['segment'],
                               tensile=source.tensile, vertex1=vertex1, vertex2=vertex3,
                               vertex3=vertex2, dip_slip=source.reverse, rtlat_slip=source.rtlat,
                               depth=vertex1[2]/1000);
    second_triangle = FaultDict(lon=one_fault_dict['lon'], lat=one_fault_dict['lat'], segment=one_fault_dict['segment'],
                                tensile=source.tensile, vertex1=vertex1, vertex2=vertex4,
                                vertex3=vertex3, dip_slip=source.reverse, rtlat_slip=source.rtlat,
                                depth=vertex1[0]/1000);
    list_of_two_triangles = [first_triangle, second_triangle];
    return list_of_two_triangles;


def change_fault_slip(one_fault_dict, rtlat=None, dipslip=None, tensile=None):
    new_rtlat = one_fault_dict['rtlat_slip'] if rtlat is None else rtlat;
    new_dipslip = one_fault_dict['dip_slip'] if dipslip is None else dipslip;
    new_tensile = one_fault_dict['tensile'] if tensile is None else tensile;
    new_fault_dict = FaultDict(lon=one_fault_dict['lon'], lat=one_fault_dict['lat'], segment=one_fault_dict['segment'],
                               depth=one_fault_dict['depth'], vertex1=one_fault_dict['vertex1'],
                               vertex2=one_fault_dict['vertex2'], vertex3=one_fault_dict['vertex3'],
                               dip_slip=new_dipslip, rtlat_slip=new_rtlat, tensile=new_tensile);
    return new_fault_dict;


def change_reference_loc(one_fault_dict):
    """Move the reference to the location of vertex1"""
    vertex1_ll, vertex2_ll, vertex3_ll = get_ll_corners(one_fault_dict);
    x2, y2 = fault_vector_functions.latlon2xy(vertex2_ll[0], vertex2_ll[1], vertex1_ll[0], vertex1_ll[1]);
    x3, y3 = fault_vector_functions.latlon2xy(vertex3_ll[0], vertex3_ll[1], vertex1_ll[0], vertex1_ll[1]);
    vertex1_newref = np.array([0, 0, one_fault_dict['vertex1'][2]]);
    vertex2_newref = np.array([x2*1000, y2*1000, one_fault_dict['vertex2'][2]]);
    vertex3_newref = np.array([x3*1000, y3*1000, one_fault_dict['vertex3'][2]]);
    new_fault_dict = FaultDict(lon=vertex1_ll[0], lat=vertex1_ll[1], segment=one_fault_dict['segment'],
                               depth=one_fault_dict['depth'],
                               vertex1=vertex1_newref, vertex2=vertex2_newref,
                               vertex3=vertex3_newref, dip_slip=one_fault_dict['dip_slip'],
                               rtlat_slip=one_fault_dict['rtlat_slip'], tensile=one_fault_dict['tensile']);
    return new_fault_dict;


def compute_triangle_centroid(one_fault_dict):
    """Compute cartesian centroid of a triangular fault. Returns a numpy array."""
    x_centroid = np.mean([one_fault_dict['vertex1'][0], one_fault_dict['vertex2'][0], one_fault_dict['vertex3'][0]]);
    y_centroid = np.mean([one_fault_dict['vertex1'][1], one_fault_dict['vertex2'][1], one_fault_dict['vertex3'][1]]);
    z_centroid = np.mean([one_fault_dict['vertex1'][2], one_fault_dict['vertex2'][2], one_fault_dict['vertex3'][2]]);
    return np.array([x_centroid, y_centroid, z_centroid]);


def compute_triangle_area(one_fault_dict):
    """Compute area of a triangle: A = 1/2 (AB x AC), with all cartesian coordinates in meters. Returns in sq-m."""
    AB = np.array([one_fault_dict['vertex2'][0] - one_fault_dict['vertex1'][0],
                   one_fault_dict['vertex2'][1] - one_fault_dict['vertex1'][1],
                   one_fault_dict['vertex2'][2] - one_fault_dict['vertex1'][2]]);
    AC = np.array([one_fault_dict['vertex3'][0] - one_fault_dict['vertex1'][0],
                   one_fault_dict['vertex3'][1] - one_fault_dict['vertex1'][1],
                   one_fault_dict['vertex3'][2] - one_fault_dict['vertex1'][2]]);
    cross_product = np.cross(AB, AC);
    total_vector_mag = np.sqrt(np.square(cross_product[0]) + np.square(cross_product[1]) + np.square(cross_product[2]));
    area_sq_m = (1/2) * total_vector_mag;
    return area_sq_m;


def get_total_moment(fault_dict_object_list, mu=30e9):
    """
    Return total moment of a list of triangle slip objects
    Moment in newton-meters
    """
    total_moment = 0;
    for item in fault_dict_object_list:
        A = compute_triangle_area(item);
        d = get_total_slip(item);
        total_moment += moment_calculations.moment_from_muad(mu, A, d);
    return total_moment;
