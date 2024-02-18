from .okada_pt_src import DC3D0
import numpy as np
from ..disp_points_object.disp_points_object import Displacement_points


def compute_cartesian_def_point(inputs, params, disp_points):
    """
    Adds the results of the point-source deformation to the existing disp_points objects

    :param inputs:
    :param params:
    :param disp_points:
    :return:
    """
    if not disp_points:
        return []
    model_disp_points = []
    for point in disp_points:
        disp_u = np.zeros((3,))
        for source in inputs.source_object:
            if source.is_point_source:
                _, new_disp = compute_displacements_strains_point(source, point.lon, point.lat, point.depth,
                                                                  params.alpha)
                disp_u = np.add(disp_u, new_disp)
        model_point = Displacement_points(lon=point.lon, lat=point.lat,
                                          dE_obs=point.dE_obs+disp_u[0],
                                          dN_obs=point.dN_obs+disp_u[1],
                                          dU_obs=point.dU_obs+disp_u[2], name=point.name)
        model_disp_points.append(model_point)
    return model_disp_points


def compute_cartesian_strain_point(inputs, params, strain_points):
    """
    Adds the results of the point-source strain

    :param inputs:
    :param params:
    :param strain_points:
    :return: list of strain tensors
    """
    if not strain_points:
        return []
    strain_tensors = []
    for point in strain_points:
        point_strain_tensor = np.zeros((3, 3))
        for source in inputs.source_object:
            if source.is_point_source:
                new_strain, _ = compute_displacements_strains_point(source, point.lon, point.lat, point.depth,
                                                                    params.alpha)
                point_strain_tensor = np.add(point_strain_tensor, new_strain)
        strain_tensors.append(point_strain_tensor)
    return strain_tensors


def compute_displacements_strains_point(source, x, y, z, alpha):
    """
    The main math of DC3D
    Operates on a source object (e.g., fault),
    and an xyz position in the same cartesian reference frame.
    """
    R = source.R
    R2 = source.R2

    # Compute the position relative to the translated, rotated fault.
    translated_pos = np.array(
        [[x - source.xstart], [y - source.ystart], [-z]])
    xyz = R.dot(translated_pos)
    success, u, grad_u = DC3D0(alpha, xyz[0], xyz[1], xyz[2], source.top, source.dipangle,
                               source.potency[0], source.potency[1], source.potency[2], source.potency[3])
    grad_u = grad_u * 1e-9  # DC3D0 Unit correction: potency from N-m results in strain in nanostrain
    u = u * 1e-6  # Unit correction: potency from N-m results in displacements in microns.
    # Solve for displacement gradients at certain xyz position

    # Rotate grad_u back into the unprimed coordinates.
    desired_coords_grad_u = np.dot(R2, np.dot(grad_u, R2.T))
    desired_coords_u = R2.dot(np.array([[u[0]], [u[1]], [u[2]]]))
    return desired_coords_grad_u, desired_coords_u
