import scipy.io
import numpy as np
from osgeo import osr  # gdal library, works inside pygmt environment
from . import fault_slip_triangle


def convert_points_to_wgs84(x, y):
    """
    Converts points from UTM zone 11 to lat/lon coordinates

    :param x: float, easting (m)
    :param y: float, northing (m)
    :returns: tuple of (lat, lon, depth (m, positive down) )
    """
    # Equivalent Command line API: gdaltransform -s_srs EPSG:32611 -t_srs EPSG:4326
    src = osr.SpatialReference()
    tgt = osr.SpatialReference()
    src.ImportFromEPSG(32611)    # source: UTM Zone 11
    tgt.ImportFromEPSG(4326)   # destination: WGS84 CGS
    transform = osr.CoordinateTransformation(src, tgt);
    newtuple = transform.TransformPoint(x, y);  # TRANSFORM
    return newtuple;


def read_brawley_lohman_2005(filename):
    """
    Read a matlab structure from Rowena Lohman, originally reported in utm zone 11 easting and northing meters
    Matlab structure from Rowena Lohman, from McGuire et al. 2015:
    dict_keys(['__header__', '__version__', '__globals__', 'xfault', 'yfault', 'zfault'])
    xfault = 3 arrays of 250 each in the range of 600K [632840.42547992]  UTM Zone 11
    yfault = 3 arrays of 250 each in the range of 3M [3670484.48912303]  UTM Zone 11
    zfault = 3 arrays of 250 each in the range of ~1000, assuming meters below the surface
    """
    print("Reading file %s " % filename);
    triangle_list = [];
    mat = scipy.io.loadmat(filename)
    (reference_lat, reference_lon, ref_depth) = convert_points_to_wgs84(mat['xfault'][0][0], mat['yfault'][0][0])
    for i in range(len(mat['xfault'][0])):  # collect lon/lat of 3 vertices of the triangles
        first_vertex = np.array([mat['xfault'][0][i]-mat['xfault'][0][0], mat['yfault'][0][i]-mat['yfault'][0][0],
                                mat['zfault'][0][i]]);
        second_vertex = np.array([mat['xfault'][1][i]-mat['xfault'][0][0], mat['yfault'][1][i]-mat['yfault'][0][0],
                                 mat['zfault'][1][i]]);
        third_vertex = np.array([mat['xfault'][2][i]-mat['xfault'][0][0], mat['yfault'][2][i]-mat['yfault'][0][0],
                                mat['zfault'][2][i]]);
        new_triangle = fault_slip_triangle.TriangleFault(lon=reference_lon, lat=reference_lat, dip_slip=0,
                                                         rtlat_slip=0, tensile=0, segment=0, vertex1=first_vertex,
                                                         vertex2=second_vertex, vertex3=third_vertex,
                                                         depth=first_vertex[2]/1000);
        triangle_list.append(new_triangle);
    print("--> Returning %d triangular fault patches" % len(triangle_list));
    return triangle_list;
