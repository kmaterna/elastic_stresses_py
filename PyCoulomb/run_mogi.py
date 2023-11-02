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
    :returns: list of 3 floats
    """
    u_disp, v_disp, w_disp = 0, 0, 0;

    for source in sources:
        dx, dy, dz = compute_disps_from_one_mogi(source, x, y, nu, compute_depth);
        # Update the displacements from all sources
        u_disp = u_disp + dx;
        v_disp = v_disp + dy;
        w_disp = w_disp + dz;  # vertical
    return u_disp, v_disp, w_disp;


def compute_disps_from_one_mogi(source, x, y, nu, _compute_depth=0):
    """
    Calculates surface deformation based on point source
    References: Mogi 1958, Segall 2010 p.203

    Mogi source at x=source.xstart, y=source.ystart in km, volume = source.dV in m^3, depth=source.depth in km.
    Compute the displacement at (x, y) in km, at depth=compute_depth in km
    Returns in meters.
    """
    dx = 1000 * (x - source.xstart)
    dy = 1000 * (y - source.ystart)

    # Convert to surface cylindrical coordinates
    th, rho = cart2pol(dx, dy)
    R = np.hypot(source.depth*1000, rho)  # in meters

    # Mogi displacement calculation
    C = ((1-nu) / np.pi) * source.dV
    ur = C * rho / R**3
    uz = C * source.depth*1000 / R**3

    ux, uy = pol2cart(th, ur)
    return ux, uy, uz;


def cart2pol(x1, x2):
    theta = np.arctan2(x2, x1)
    r = np.hypot(x2, x1)
    return theta, r


def pol2cart(theta, r):
    x1 = r * np.cos(theta)
    x2 = r * np.sin(theta)
    return x1, x2
