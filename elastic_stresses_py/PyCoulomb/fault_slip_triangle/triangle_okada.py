"""
Implementing Okada on a fault_slip_triangle object using Ben Thompson's cutde library
"""

import numpy as np
from ..disp_points_object.disp_points_object import Displacement_points
import cutde.halfspace as hs
from . import fault_slip_triangle
from .. import pyc_fault_object, utilities


def convert_rect_sources_into_tris(rect_sources):
    """
    :param rect_sources: list of sources that include PyCoulomb fault rectangles
    :return: list of triangular faults
    """
    tri_faults = []
    for source in rect_sources:
        if isinstance(source, pyc_fault_object.Faults_object):
            if source.is_point_source:
                continue
            sub_tris = fault_slip_triangle.convert_pycoulomb_rectangle_into_three_triangles(source, source.zerolon,
                                                                                            source.zerolat)
            for item in sub_tris:
                tri_faults.append(item)
        if isinstance(source, fault_slip_triangle.TriangleFault):
            tri_faults.append(source)
    return tri_faults


def compute_cartesian_strain_tris(inputs, params, strain_points):
    """
    Loop through a list of lon/lat and compute their strains due to all sources put together.
    Returns list of strain tensors
    """
    tri_faults = convert_rect_sources_into_tris(inputs.source_object)
    strain_tensors = compute_strain_points_from_triangles(tri_faults, strain_points, params.nu)
    return strain_tensors


def compute_cartesian_def_tris(inputs, params, obs_disp_points):
    tri_faults = convert_rect_sources_into_tris(inputs.source_object)
    modeled_tri_points = compute_disp_points_from_triangles(tri_faults, obs_disp_points, params.nu)
    return modeled_tri_points


def compute_disp_points_from_triangles(fault_triangles, disp_points, poisson_ratio):
    """
    Similar to run_dc3d.compute_ll_def(inputs, alpha, disp_points).
    Requires all fault_triangles to have the same reference lon/lat

    :param fault_triangles: list
    :param disp_points: list
    :param poisson_ratio: float
    :returns: list of disp_points objects
    """
    if not disp_points:
        return []
    if len(fault_triangles) == 0:
        return utilities.get_zeros_disp_points(disp_points)
    proceed_code = fault_slip_triangle.check_consistent_reference_frame(fault_triangles)
    if not proceed_code:
        raise ValueError("Error! Triangular faults do not have same reference")

    obsx = [point.lon*1000 for point in disp_points]
    obsy = [point.lat*1000 for point in disp_points]   # calculation works in meters
    obsz = [point.depth * -1000 for point in disp_points]  # in meters, negative is down
    pts = np.vstack([obsx, obsy, obsz]).T   # shape: (Npts, 3)

    slip_array = np.array([[-src.rtlat_slip, src.dip_slip, src.tensile] for src in fault_triangles])  # shape:(Ntris, 3)
    fault_pts, fault_tris = fault_slip_triangle.extract_mesh_vertices(fault_triangles)
    fault_pts = fault_slip_triangle.flip_depth_sign(fault_pts)  # fault_pts shape: N_vertices, 3
    src_tris = fault_pts[fault_tris]  # src_tris shape: (Ntris, 3, 3)
    disp_mat = hs.disp_matrix(obs_pts=pts, tris=src_tris, nu=poisson_ratio)  # disp_mat shape: (Npts, 3, Ntris, 3)
    disp = disp_mat.reshape((-1, np.size(slip_array))).dot(slip_array.flatten())   # reshape by len of total slip vector
    disp_grid = disp.reshape((*np.array(obsx).shape, 3))  # disp_grid shape: Npts, 3

    # Package the results into usable formats
    modeled_disp_points = []
    for i, item in enumerate(disp_points):
        new_disp_pt = Displacement_points(lon=item.lon, lat=item.lat, dE_obs=disp_grid[i][0],
                                          dN_obs=disp_grid[i][1], dU_obs=disp_grid[i][2], Se_obs=0, Sn_obs=0, Su_obs=0,
                                          depth=item.depth, endtime=item.endtime, starttime=item.starttime,
                                          meas_type=item.meas_type, refframe=item.refframe, name=item.name)
        modeled_disp_points.append(new_disp_pt)

    return modeled_disp_points


def compute_strain_points_from_triangles(fault_triangles, strain_points, poisson_ratio):
    """
    Similar to run_dc3d.compute_ll_def(inputs, alpha, disp_points).
    Requires all fault_triangles to have the same reference lon/lat

    :param fault_triangles: list
    :param strain_points: list
    :param poisson_ratio: float
    :returns: list of disp_points objects, list of strain tensors in 3x3 matrix
    """
    if not strain_points:
        return []
    if len(fault_triangles) == 0:
        return utilities.get_zeros_strain_points(strain_points)
    proceed_code = fault_slip_triangle.check_consistent_reference_frame(fault_triangles)
    if not proceed_code:
        raise ValueError("Error! Triangular faults do not have same reference")

    obsx = [point.lon*1000 for point in strain_points]
    obsy = [point.lat*1000 for point in strain_points]   # calculation works in meters
    obsz = [point.depth * -1000 for point in strain_points]  # in meters, negative is down
    pts = np.vstack([obsx, obsy, obsz]).T   # shape: (Npts, 3)

    slip_array = np.array([[-src.rtlat_slip, src.dip_slip, src.tensile] for src in fault_triangles])  # shape:(Ntris, 3)
    fault_pts, fault_tris = fault_slip_triangle.extract_mesh_vertices(fault_triangles)
    fault_pts = fault_slip_triangle.flip_depth_sign(fault_pts)  # fault_pts shape: N_vertices, 3
    src_tris = fault_pts[fault_tris]  # src_tris shape: (Ntris, 3, 3)

    # Get strain
    strain_mat = hs.strain_matrix(obs_pts=pts, tris=src_tris, nu=poisson_ratio)  # strain_mat shape: (Npts, 6, Ntris, 3)
    strain = strain_mat.reshape((-1, np.size(slip_array))).dot(slip_array.flatten())  # reshape by len total slip vector
    strain_tensors = strain.reshape((*np.array(obsx).shape, 6))  # strain_tensors shape: Npts, 6
    # strain[:,0] is the xx component of strain, 1 is yy, 2 is zz, 3 is xy, 4 is xz, and 5 is yz.

    modeled_strain_tensors = []
    for i, item in enumerate(strain_tensors):
        comps = strain_tensors[i]
        new_strain_tensor = np.array([[comps[0], comps[3], comps[4]],
                                      [comps[3], comps[1], comps[5]],
                                      [comps[4], comps[5], comps[2]]])
        modeled_strain_tensors.append(new_strain_tensor)
    return modeled_strain_tensors
