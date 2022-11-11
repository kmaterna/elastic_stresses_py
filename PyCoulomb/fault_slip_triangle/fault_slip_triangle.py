
from typing import TypedDict
import numpy as np
from Tectonic_Utils.geodesy import fault_vector_functions
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
    zerolon, zerolat = one_fault_dict['lon'], one_fault_dict['lat']
    [source] = fault_dict_to_coulomb_fault([one_fault_dict]);
    [x_all, y_all, _, _] = conversion_math.get_fault_four_corners(source);  # This is cartesian
    top_depth, bottom_depth = source.top, source.bottom;
    vertex1 = np.array([x_all[0]*1000, y_all[0]*1000, top_depth*1000]);
    vertex2 = np.array([x_all[1]*1000, y_all[1]*1000, top_depth*1000]);
    vertex3 = np.array([x_all[2]*1000, y_all[2]*1000, bottom_depth*1000]);
    vertex4 = np.array([x_all[3]*1000, y_all[3]*1000, bottom_depth*1000]);
    first_triangle = FaultDict(lon=zerolon, lat=zerolat, segment=0, tensile=source.tensile,
                               vertex1=vertex1, vertex2=vertex3,
                               vertex3=vertex2, dip_slip=source.reverse, rtlat_slip=source.rtlat);
    second_triangle = FaultDict(lon=zerolon, lat=zerolat, segment=0, tensile=source.tensile,
                                vertex1=vertex1, vertex2=vertex4,
                                vertex3=vertex3, dip_slip=source.reverse, rtlat_slip=source.rtlat);
    list_of_two_triangles = [first_triangle, second_triangle];

    return list_of_two_triangles;
