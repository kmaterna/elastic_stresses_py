import numpy as np
from .. import fault_slip_triangle


def write_slippy_like_triangles(fault_object_list, filename):
    """
    Write the slippy-like file format for triangular fault slip elements
    """
    print("Writing file %s" % filename);
    ofile = open(filename, 'w');
    ofile.write("# x1 y1 z1 x2 y2 z2 x3 y3 z3[m] lon[deg] lat[deg] depth[km] rtlat_slip[m] dip_slip[m] "
                "tensile[m] segment[int]\n");
    for item in fault_object_list:
        ofile.write("%f %f %f " % (item.vertex1[0], item.vertex1[1], item.vertex1[2]));
        ofile.write("%f %f %f " % (item.vertex2[0], item.vertex2[1], item.vertex2[2]));
        ofile.write("%f %f %f " % (item.vertex3[0], item.vertex3[1], item.vertex3[2]));
        ofile.write("%f %f %f " % (item.lon, item.lat, item.depth));
        ofile.write("%f %f %f %d\n" % (item.rtlat_slip, item.dip_slip, item.tensile, item.segment));
    ofile.close();
    return;


def read_slippy_like_triangles(filename):
    """
    Read the slippy-like file format for triangular fault slip elements
    """
    fault_triangle_list = [];
    print("Reading file %s" % filename);
    [x1, y1, z1, x2, y2, z2, x3, y3, z3] = np.loadtxt(filename, skiprows=1, unpack=True,
                                                      dtype={"names": ('x1', 'y1', 'z1', 'x2',
                                                                       'y2', 'z2', 'x3', 'y3',
                                                                       'z3'),
                                                             "formats": (float, float, float, float,
                                                                         float, float, float, float,
                                                                         float)});

    [lon, lat, depth, rtlat_slip, dip_slip, tensile, segment] = np.loadtxt(filename, skiprows=1, unpack=True,
                                                                           dtype={"names": ('lon', 'lat', 'depth',
                                                                                            'rtlat_slip',
                                                                                            'dip_slip',
                                                                                            'tensile', 'segment'),
                                                                                  "formats": (
                                                                                  float, float, float, float,
                                                                                  float, float, int)});
    for i in range(len(x1)):
        new_triangle = fault_slip_triangle.TriangleFault(vertex1=[x1, y1, z1], vertex2=[x2, y2, z2],
                                                         vertex3=[x3, y3, z3],
                                                         lon=lon, lat=lat, depth=depth, rtlat_slip=rtlat_slip,
                                                         dip_slip=dip_slip, tensile=tensile, segment=segment);
        fault_triangle_list.append(new_triangle);
    print("--> Returning %d fault patches " % len(fault_triangle_list));
    return fault_triangle_list;
