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
        [xi, yi] = fault_vector_functions.latlon2xy(point.lon, point.lat,
                                                    fault_triangles[0].lon, fault_triangles[0].lat)
        obsx.append(xi*1000)
        obsy.append(yi*1000)  # calculation works in meters
    pts = np.vstack([obsx, obsy, np.zeros(np.shape(obsx))]).T   # shape: (Npts, 3)
    resulting_model = np.zeros(np.shape(pts))

    slip_array = []
    for source in fault_triangles:
        slip_array.append([-source.rtlat_slip, source.dip_slip, source.tensile])
    slip_array = np.array(slip_array)

    for source in fault_triangles:
        fault_pts = np.array([[source.vertex1[0], source.vertex1[1], -source.vertex1[2]],
                              [source.vertex2[0], source.vertex2[1], -source.vertex2[2]],
                              [source.vertex3[0], source.vertex3[1], -source.vertex3[2]]])  # vertex coords
        fault_tris = np.array([[0, 1, 2]], dtype=np.int64)  # triangles, indexing into vertices array
        src_tris = fault_pts[fault_tris]  # src_tris shape (Ntris, 3, 3)
        disp_mat = HS.disp_matrix(obs_pts=pts, tris=src_tris, nu=poisson_ratio)  # disp_mat: shape (Npts, 3, Ntris, 3)
        slip = np.array([[-source.rtlat_slip, source.dip_slip, source.tensile]])  # shape: (1, 3)

        disp = disp_mat.reshape((-1, 3)).dot(slip.flatten())   # 3 here is the length of the slip vector
        disp_grid = disp.reshape((*np.array(obsx).shape, 3))
        resulting_model = np.add(resulting_model, disp_grid)

        # Get strain
        strain_mat = HS.strain_matrix(obs_pts=pts, tris=src_tris, nu=poisson_ratio)

    modeled_disp_points = []
    for i, item in enumerate(disp_points):
        new_disp_pt = Displacement_points(lon=item.lon, lat=item.lat, dE_obs=resulting_model[i][0],
                                          dN_obs=resulting_model[i][1], dU_obs=resulting_model[i][2], Se_obs=0,
                                          Sn_obs=0, Su_obs=0, endtime=item.endtime, starttime=item.starttime,
                                          meas_type=item.meas_type, refframe=item.refframe, name=item.name)
        modeled_disp_points.append(new_disp_pt)
    return modeled_disp_points
