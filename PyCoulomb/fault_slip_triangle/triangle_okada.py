"""
Implementing Okada on a fault_slip_triangle object using Ben Thompson's cutde library
A work in progress
"""

import numpy as np
import cutde.fullspace as HS


def compute_disp_points_from_triangles(fault_triangles, disp_points, poisson_ratio):
    """
    Similar to run_dc3d.compute_ll_def(inputs, alpha, disp_points)
    """
    for source in fault_triangles:
        pts = [];
        # pts = np.array([obsx, obsy, 0 * obsy]).reshape((3, -1)).T.copy()
        fault_pts = np.array([[source['vertex1'][0], source['vertex1'][1], source['vertex1'][2]],
                              [source['vertex2'][0], source['vertex2'][1], source['vertex2'][2]],
                              [source['vertex3'][0], source['vertex3'][1], source['vertex3'][2]]])  # vertex coordinates
        fault_tris = np.array([[0, 1, 2]], dtype=np.int64)  # triangles, indexing into vertices array
        disp_mat = HS.disp_matrix(obs_pts=pts, tris=fault_pts[fault_tris], nu=poisson_ratio);
        slip = np.array([[source['rtlat_slip'], source['dip_slip'], source['tensile']]]);
        # disp = disp_mat.reshape((-1, 6)).dot(slip.flatten())

    return [];
