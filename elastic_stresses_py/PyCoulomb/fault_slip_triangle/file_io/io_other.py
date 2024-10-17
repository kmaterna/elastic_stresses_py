import scipy.io
import numpy as np
from osgeo import osr  # gdal library, works inside pygmt environment
from .. import fault_slip_triangle
from Tectonic_Utils.geodesy import fault_vector_functions as fvf


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
    transform = osr.CoordinateTransformation(src, tgt)
    newtuple = transform.TransformPoint(x, y)  # TRANSFORM
    return newtuple


def read_brawley_lohman_2005(filename):
    """
    Read a matlab structure from Rowena Lohman, originally reported in utm zone 11 easting and northing meters
    Matlab structure from Rowena Lohman, from McGuire et al. 2015:
    dict_keys(['__header__', '__version__', '__globals__', 'xfault', 'yfault', 'zfault'])
    xfault = 3 arrays of 250 each in the range of 600K [632840.42547992]  UTM Zone 11
    yfault = 3 arrays of 250 each in the range of 3M [3670484.48912303]  UTM Zone 11
    zfault = 3 arrays of 250 each in the range of ~1000, assuming meters below the surface
    """
    print("Reading file %s " % filename)
    triangle_list = []
    mat = scipy.io.loadmat(filename)
    (reference_lat, reference_lon, ref_depth) = convert_points_to_wgs84(mat['xfault'][0][0], mat['yfault'][0][0])
    for i in range(len(mat['xfault'][0])):  # collect lon/lat of 3 vertices of the triangles
        first_vertex = np.array([mat['xfault'][0][i]-mat['xfault'][0][0], mat['yfault'][0][i]-mat['yfault'][0][0],
                                mat['zfault'][0][i]])
        second_vertex = np.array([mat['xfault'][1][i]-mat['xfault'][0][0], mat['yfault'][1][i]-mat['yfault'][0][0],
                                 mat['zfault'][1][i]])
        third_vertex = np.array([mat['xfault'][2][i]-mat['xfault'][0][0], mat['yfault'][2][i]-mat['yfault'][0][0],
                                mat['zfault'][2][i]])
        new_triangle = fault_slip_triangle.TriangleFault(lon=reference_lon, lat=reference_lat, dip_slip=0,
                                                         rtlat_slip=0, tensile=0, segment=0, vertex1=first_vertex,
                                                         vertex2=second_vertex, vertex3=third_vertex,
                                                         depth=float(first_vertex[2])/1000)
        triangle_list.append(new_triangle)
    print("--> Returning %d triangular fault patches" % len(triangle_list))
    return triangle_list


def extract_given_patch_helper(nodes, idx):
    """ nodes = 1859 x 3.  idx = [a b c]."""
    xs = [nodes[idx[0]-1][0], nodes[idx[1]-1][0], nodes[idx[2]-1][0], nodes[idx[0]-1][0]]
    ys = [nodes[idx[0]-1][1], nodes[idx[1]-1][1], nodes[idx[2]-1][1], nodes[idx[0]-1][1]]
    depths = [nodes[idx[0]-1][2], nodes[idx[1]-1][2], nodes[idx[2]-1][2], nodes[idx[0]-1][2]]
    return xs, ys, depths


def read_csz_bartlow_2019(input_file):
    """
    Read matlab file format used for the inversions in Materna et al., 2019.
    Returns a list of triangular fault patches and a list of node points.
    """
    print("Reading file %s " % input_file)
    data_structure = scipy.io.loadmat(input_file)
    nodes = data_structure['nd_ll']  # all the nodes for the entire CSZ, a big array from Canada to MTJ.
    elements = data_structure['el']  # elements

    # Open all the fault patches
    fault_patches = []
    for i, item in enumerate(elements):
        xs, ys, depths = extract_given_patch_helper(nodes, item)
        reflon, reflat = xs[0], ys[0]
        v1x, v1y = fvf.latlon2xy_single(xs[0], ys[0], reflon, reflat)
        v2x, v2y = fvf.latlon2xy_single(xs[1], ys[1], reflon, reflat)
        v3x, v3y = fvf.latlon2xy_single(xs[2], ys[2], reflon, reflat)
        new_ft = fault_slip_triangle.TriangleFault(vertex1=[v1x*1000, v1y*1000, depths[0]*-1000],
                                                   vertex2=[v2x*1000, v2y*1000, depths[1]*-1000],
                                                   vertex3=[v3x*1000, v3y*1000, depths[2]*-1000],
                                                   lon=reflon, lat=reflat, depth=depths[0],
                                                   rtlat_slip=0, dip_slip=1)
        fault_patches.append(new_ft)
    print("--> Returning %d triangular fault patches" % len(fault_patches))
    return fault_patches, nodes


def read_superstition_hills_mesh_2024(triangles, mesh, slip):
    """
    Reading the triangular mesh of the Superstition Hills model from Vavra et al., GRL, 2024

    :param triangles: string, filename
    :param mesh: string, filename
    :param slip: string, filename
    :return: list of fault patches
    """
    print("Reading file %s " % triangles)
    reflon, reflat = -115.70124, 32.93049  # creepmeter location
    iv1, iv2, iv3 = np.loadtxt(triangles, unpack=True, skiprows=1, usecols=(0, 1, 2))
    slip_mm = np.loadtxt(slip, unpack=True, skiprows=1, usecols=(0, ))
    vx, vy, vz = np.loadtxt(mesh, unpack=True, skiprows=2, usecols=(0, 1, 2))  # in km, with negative meaning down

    # Open all the fault patches
    fault_patches = []
    for i in range(len(iv1)):
        index1, index2, index3 = int(iv1[i]), int(iv2[i]), int(iv3[i])
        new_ft = fault_slip_triangle.TriangleFault(vertex1=[vx[index1]*1000, vy[index1]*1000, vz[index1]*-1000],
                                                   vertex2=[vx[index2]*1000, vy[index2]*1000, vz[index2]*-1000],
                                                   vertex3=[vx[index3]*1000, vy[index3]*1000, vz[index3]*-1000],
                                                   lon=reflon, lat=reflat, depth=0,
                                                   rtlat_slip=float(slip_mm[i])*0.001, dip_slip=0)
        fault_patches.append(new_ft)
    print("--> Returning %d triangular fault patches" % len(fault_patches))
    return fault_patches
