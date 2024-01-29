"""
Implementing Okada on a fault_slip_triangle object using Ben Thompson's cutde library
"""

import numpy as np
from Tectonic_Utils.geodesy import fault_vector_functions
from ..disp_points_object.disp_points_object import Displacement_points
import cutde.halfspace as HS
from . import fault_slip_triangle


def compute_disp_points_from_triangles(fault_triangles, disp_points, poisson_ratio):
    """
    Similar to run_dc3d.compute_ll_def(inputs, alpha, disp_points). Only lon and lat of disp_points will be used.
    Requires all fault_triangles to have the same reference lon/lat

    :param fault_triangles: list
    :param disp_points: list
    :param poisson_ratio: float
    :returns: list of disp_points objects
    """
    proceed_code = fault_slip_triangle.check_consistent_reference_frame(fault_triangles)
    if not proceed_code:
        raise ValueError("Error! Triangular faults do not have same reference")
    obsx, obsy, pts = [], [], []

    for point in disp_points:
        [x, y] = fault_vector_functions.latlon2xy(point.lon, point.lat, fault_triangles[0].lon, fault_triangles[0].lat)
        obsx.append(x*1000)
        obsy.append(y*1000)  # calculation works in meters
    pts = np.vstack([obsx, obsy, np.zeros(np.shape(obsx))]).T   # shape: (Npts, 3)

    slip_array = np.array([[-src.rtlat_slip, src.dip_slip, src.tensile] for src in fault_triangles])  # shape:(Ntris, 3)
    fault_pts, fault_tris = fault_slip_triangle.extract_mesh_vertices(fault_triangles)
    fault_pts = fault_slip_triangle.flip_depth_sign(fault_pts)  # fault_pts shape: N_vertices, 3
    src_tris = fault_pts[fault_tris]  # src_tris shape: (Ntris, 3, 3)
    disp_mat = HS.disp_matrix(obs_pts=pts, tris=src_tris, nu=poisson_ratio)  # disp_mat shape: (Npts, 3, Ntris, 3)
    disp = disp_mat.reshape((-1, np.size(slip_array))).dot(slip_array.flatten())   # reshape by len of total slip vector
    disp_grid = disp.reshape((*np.array(obsx).shape, 3))  # disp_grid shape: Npts, 3

    # Get strain
    strain_mat = HS.strain_matrix(obs_pts=pts, tris=src_tris, nu=poisson_ratio)  # strain_mat shape: (Npts, 6, Ntris, 3)

    modeled_disp_points = []
    for i, item in enumerate(disp_points):
        new_disp_pt = Displacement_points(lon=item.lon, lat=item.lat, dE_obs=disp_grid[i][0],
                                          dN_obs=disp_grid[i][1], dU_obs=disp_grid[i][2], Se_obs=0,
                                          Sn_obs=0, Su_obs=0, endtime=item.endtime, starttime=item.starttime,
                                          meas_type=item.meas_type, refframe=item.refframe, name=item.name)
        modeled_disp_points.append(new_disp_pt)
    return modeled_disp_points
