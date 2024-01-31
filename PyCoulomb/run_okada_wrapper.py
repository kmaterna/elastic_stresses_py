import numpy as np
from okada_wrapper import dc3dwrapper, dc3d0wrapper


def compute_strains_stresses_from_one_fault(source, x, y, z, alpha):
    """
    The main math of DC3D
    Operates on a source object (e.g., fault),
    and an xyz position in the same cartesian reference frame.
    """
    R = source.R
    R2 = source.R2
    strike_slip = source.rtlat * -1  # The dc3d coordinate system has left-lateral positive.

    # Compute the position relative to the translated, rotated fault.
    translated_pos = np.array(
        [[x - source.xstart], [y - source.ystart], [-z]])
    xyz = R.dot(translated_pos)
    if source.potency:
        success, u, grad_u = dc3d0wrapper(alpha, [xyz[0], xyz[1], xyz[2]], source.top, source.dipangle,
                                          [source.potency[0], source.potency[1], source.potency[2],
                                           source.potency[3]])
        grad_u = grad_u * 1e-9  # DC3D0 Unit correction: potency from N-m results in strain in nanostrain
        u = u * 1e-6  # Unit correction: potency from N-m results in displacements in microns.
    else:
        success, u, grad_u = dc3dwrapper(alpha, [xyz[0], xyz[1], xyz[2]], source.top, source.dipangle,
                                         [0, source.L], [-source.W, 0],
                                         [strike_slip, source.reverse, source.tensile])
        grad_u = grad_u * 1e-3  # DC3D Unit correction.
    # Solve for displacement gradients at certain xyz position

    # Rotate grad_u back into the unprimed coordinates.
    desired_coords_grad_u = np.dot(R2, np.dot(grad_u, R2.T))
    desired_coords_u = R2.dot(np.array([[u[0]], [u[1]], [u[2]]]))
    return desired_coords_grad_u, desired_coords_u
