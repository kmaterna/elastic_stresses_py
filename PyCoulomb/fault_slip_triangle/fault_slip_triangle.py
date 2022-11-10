
from typing import TypedDict
import numpy as np
from Tectonic_Utils.geodesy import fault_vector_functions

class FaultDict(TypedDict):
    """
    The internal format is a dictionary for a triangular fault segment:
    {
        vertex1 [x, y, z] in m,
        vertex2 [x, y, z] in m,
        vertex3 [x, y, z] in m,
        lon(coord reference),
        lat(coord reference),
        rt_lat slip(m),
        dip_slip(m),
        tensile(m),
        segment(int)
    }
    If the fault is a receiver fault, we put slip = 0
    """
    vertex1: np.array
    vertex2: np.array
    vertex3: np.array
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
            # Convert m to km before coordinate conversion
            lon1, lat1 = fault_vector_functions.xy2lonlat_single(item['vertex1'][0]/1000, item['vertex1'][1]/1000,
                                                                 item['lon'], item['lat']);
            lon2, lat2 = fault_vector_functions.xy2lonlat_single(item['vertex2'][0]/1000, item['vertex2'][1]/1000,
                                                                 item['lon'], item['lat']);
            lon3, lat3 = fault_vector_functions.xy2lonlat_single(item['vertex3'][0]/1000, item['vertex3'][1]/1000,
                                                                 item['lon'], item['lat']);
            ofile.write("> -Z%f \n" % total_slip);
            ofile.write("%f %f \n" % (lon1, lat1));
            ofile.write("%f %f \n" % (lon2, lat2));
            ofile.write("%f %f \n" % (lon3, lat3));
            ofile.write("%f %f \n" % (lon1, lat1));
    return;


def get_total_slip(triangle_fault):
    total_slip = np.sqrt(np.square(triangle_fault['dip_slip']) + np.square(triangle_fault['rtlat_slip']));
    return total_slip;
