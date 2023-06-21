import numpy as np


def compute_surface_disp_point(sources, nu, x, y, compute_depth=0):
    """
    A compute loop for each source object at one x/y point.
    x/y in the same coordinate system as the fault object. Computes displacement.

    :param sources: list of mogi_source objects
    :param nu: float, poisson's ratio
    :param x: float
    :param y: float
    :param compute_depth: depth of observation. Default depth is at surface of earth
    """
    u_disp, v_disp, w_disp = 0, 0, 0;

    for source in sources:
        dx, dy, dz = compute_disps_from_one_mogi(source, x, y, nu, compute_depth);
        # Update the displacements from all sources
        u_disp = u_disp + dx;
        v_disp = v_disp + dy;
        w_disp = w_disp + dz;  # vertical

    return u_disp, v_disp, w_disp;


def compute_disps_from_one_mogi(source, x, y, nu, compute_depth=0):
    x, y, z = 0, 0, 0;
    return x, y, z;
