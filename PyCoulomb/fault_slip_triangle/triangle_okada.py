"""
Implementing Okada on a fault_slip_triangle object using Ben Thompson's cutde library
"""

import numpy as np
from Tectonic_Utils.geodesy import fault_vector_functions
from ..disp_points_object.disp_points_object import Displacement_points
import cutde.halfspace as HS
from . import fault_slip_triangle
from .. import pyc_fault_object

def convert_rect_sources_into_tris(rect_sources):
    """
    :param rect_sources: list of sources that include PyCoulomb fault rectangles
    :return: list of triangular faults
    """
    tri_faults = []
    for source in rect_sources:
        if isinstance(source, pyc_fault_object.Faults_object):
            two_tris = fault_slip_triangle.convert_pycoulomb_rectangle_into_two_triangles(source, source.zerolon,
                                                                                          source.zerolat)
            tri_faults.append(two_tris[0])
            tri_faults.append(two_tris[1])
    return tri_faults


def compute_ll_def_tris(inputs, params, obs_disp_points, coords='geographic'):
    """
    Similar to the okada_wrapper version in run_dc3d.py.

    :param inputs: pycoulomb Inputs object
    :param params: pycoulomb Params object
    :param obs_disp_points: list of disp_points
    :param coords: string telling us whether disp_points are in km or lon/lat
    :return: list of disp_points
    """
    if not obs_disp_points:
        return []
    if isinstance(obs_disp_points, Displacement_points):
        obs_disp_points = [obs_disp_points]
    tri_faults = convert_rect_sources_into_tris(inputs.source_object)
    modeled_tri_points, _ = compute_disp_points_from_triangles(tri_faults, obs_disp_points, params.nu, coords)
    return modeled_tri_points


def compute_ll_strain_tris(inputs, params, strain_points, coords='geographic'):
    """
    Similar to the compute_ll_strain in PyCoulomb.
    Loop through a list of lon/lat and compute their strains due to all sources put together.
    """
    if not strain_points:
        return []
    if isinstance(strain_points, Displacement_points):
        strain_points = [strain_points]
    tri_faults = convert_rect_sources_into_tris(inputs.source_object)
    _, strain_tensor = compute_disp_points_from_triangles(tri_faults, strain_points, params.nu, coords=coords)
    return strain_tensor


def compute_disp_points_from_triangles(fault_triangles, disp_points, poisson_ratio, coords='geographic'):
    """
    Similar to run_dc3d.compute_ll_def(inputs, alpha, disp_points). Only lon and lat of disp_points will be used.
    Requires all fault_triangles to have the same reference lon/lat

    :param fault_triangles: list
    :param disp_points: list
    :param poisson_ratio: float
    :param coords: string telling us whether disp_points are in km or lon/lat
    :returns: list of disp_points objects, list of strain tensors in 3x3 matrix
    """
    proceed_code = fault_slip_triangle.check_consistent_reference_frame(fault_triangles)
    if not proceed_code:
        raise ValueError("Error! Triangular faults do not have same reference")
    obsx, obsy, obsz, pts = [], [], [], []

    if coords == 'cartesian':
        obsx = [point.lon*1000 for point in disp_points]
        obsy = [point.lat*1000 for point in disp_points]
    else:
        for point in disp_points:
            [x, y] = fault_vector_functions.latlon2xy(point.lon, point.lat, fault_triangles[0].lon,
                                                      fault_triangles[0].lat)
            obsx.append(x*1000)
            obsy.append(y*1000)  # calculation works in meters
    obsz = [point.depth * -1000 for point in disp_points]  # in meters, negative is down
    pts = np.vstack([obsx, obsy, obsz]).T   # shape: (Npts, 3)

    slip_array = np.array([[-src.rtlat_slip, src.dip_slip, src.tensile] for src in fault_triangles])  # shape:(Ntris, 3)
    fault_pts, fault_tris = fault_slip_triangle.extract_mesh_vertices(fault_triangles)
    fault_pts = fault_slip_triangle.flip_depth_sign(fault_pts)  # fault_pts shape: N_vertices, 3
    src_tris = fault_pts[fault_tris]  # src_tris shape: (Ntris, 3, 3)
    disp_mat = HS.disp_matrix(obs_pts=pts, tris=src_tris, nu=poisson_ratio)  # disp_mat shape: (Npts, 3, Ntris, 3)
    disp = disp_mat.reshape((-1, np.size(slip_array))).dot(slip_array.flatten())   # reshape by len of total slip vector
    disp_grid = disp.reshape((*np.array(obsx).shape, 3))  # disp_grid shape: Npts, 3

    # Get strain
    strain_mat = HS.strain_matrix(obs_pts=pts, tris=src_tris, nu=poisson_ratio)  # strain_mat shape: (Npts, 6, Ntris, 3)
    strain = strain_mat.reshape((-1, np.size(slip_array))).dot(slip_array.flatten())  # reshape by len total slip vector
    strain_tensors = strain.reshape((*np.array(obsx).shape, 6))  # strain_tensors shape: Npts, 6
    # strain[:,0] is the xx component of strain, 1 is yy, 2 is zz, 3 is xy, 4 is xz, and 5 is yz.

    # Package the results into usable formats
    modeled_disp_points = []
    for i, item in enumerate(disp_points):
        new_disp_pt = Displacement_points(lon=item.lon, lat=item.lat, dE_obs=disp_grid[i][0],
                                          dN_obs=disp_grid[i][1], dU_obs=disp_grid[i][2], Se_obs=0, Sn_obs=0, Su_obs=0,
                                          depth=item.depth, endtime=item.endtime, starttime=item.starttime,
                                          meas_type=item.meas_type, refframe=item.refframe, name=item.name)
        modeled_disp_points.append(new_disp_pt)

    modeled_strain_tensors = []
    for i, item in enumerate(strain_tensors):
        comps = strain_tensors[i]
        new_strain_tensor = np.array([[comps[0], comps[3], comps[4]],
                                      [comps[3], comps[1], comps[5]],
                                      [comps[4], comps[5], comps[2]]])
        modeled_strain_tensors.append(new_strain_tensor)
    return modeled_disp_points, modeled_strain_tensors
