# Running dc3d, given an input namedtuple.

import numpy as np
from okada_wrapper import dc3dwrapper, dc3d0wrapper
from . import coulomb_collections as cc
from . import conversion_math
from . import run_mogi, utilities
from .disp_points_object.disp_points_object import Displacement_points
from Tectonic_Utils.geodesy import fault_vector_functions


def do_stress_computation(params, inputs, disp_points=(), strain_points=()):
    """
    * Step 0. Split receiver fault into many sub-faults if necessary
    * Step 1. Compute strains and displacements
    * Step 2. Resolve stresses on receiver faults
    """

    print("Beginning stress calculation.");
    print("Number of sources: %d " % len(inputs.source_object));
    print("Number of receivers: %d " % len(inputs.receiver_object));
    subfaulted_inputs = split_subfault_receivers(params, inputs);

    # Computes here.
    [x, y, x2d, y2d, u_disps, v_disps, w_disps] = compute_grid_def(subfaulted_inputs, params);
    model_disp_points = compute_ll_def(subfaulted_inputs, params.alpha, disp_points);
    strain_tensor_results = compute_ll_strain(subfaulted_inputs, params.alpha, strain_points);
    [receiver_normal, receiver_shear, receiver_coulomb] = compute_strains_stresses(params, subfaulted_inputs);
    receiver_profile_results = compute_stresses_horiz_profile(params, subfaulted_inputs);

    MyOutObject = cc.Out_object(x=x, y=y, x2d=x2d, y2d=y2d, u_disp=u_disps, v_disp=v_disps, w_disp=w_disps,
                                strains=strain_tensor_results, model_disp_points=model_disp_points,
                                zerolon=inputs.zerolon, zerolat=inputs.zerolat,
                                source_object=inputs.source_object, receiver_object=subfaulted_inputs.receiver_object,
                                receiver_normal=receiver_normal, receiver_shear=receiver_shear,
                                receiver_coulomb=receiver_coulomb, receiver_profile=receiver_profile_results);
    return MyOutObject;


def split_subfault_receivers(params, inputs):
    strike_split = params.strike_num_receivers;
    dip_split = params.dip_num_receivers;

    if strike_split == 1 and dip_split == 1:
        # If we're not splitting the subfaults...
        subfaulted_receivers = inputs.receiver_object;
        print("Not subdividing receiver faults further.");
    else:
        subfaulted_receivers = [];
        print("Splitting %d receiver faults into %d subfaults each." % (
            len(inputs.receiver_object), strike_split * dip_split));

        for fault in inputs.receiver_object:  # for each receiver...
            # We find the depths corresponding to the tops and bottoms of our new sub-faults
            zsplit_array = get_split_z_array(fault.top, fault.bottom, dip_split);

            for j in range(dip_split):  # First we split it up by dip.
                # Get the new coordinates of the top of the fault plane.
                W = fault_vector_functions.get_downdip_width(fault.top, zsplit_array[j], fault.dipangle);
                vector_mag = W * np.cos(
                    np.deg2rad(fault.dipangle));  # how far the bottom edge is displaced downdip from map-view

                # Get the starting points for the next row of fault subpatches.
                [start_x_top, start_y_top] = fault_vector_functions.add_vector_to_point(fault.xstart, fault.ystart,
                                                                                        vector_mag, fault.strike + 90);
                [finish_x_top, finish_y_top] = fault_vector_functions.add_vector_to_point(fault.xfinish, fault.yfinish,
                                                                                          vector_mag, fault.strike +
                                                                                          90);

                [xsplit_array, ysplit_array] = get_split_x_y_arrays(start_x_top, finish_x_top, start_y_top,
                                                                    finish_y_top, strike_split);

                for k in range(strike_split):
                    single_subfaulted_receiver = cc.construct_pycoulomb_fault(xstart=xsplit_array[k],
                                                                              xfinish=xsplit_array[k + 1],
                                                                              ystart=ysplit_array[k],
                                                                              yfinish=ysplit_array[k + 1],
                                                                              Kode=fault.Kode, rtlat=0, reverse=0,
                                                                              tensile=0, potency=[],
                                                                              strike=fault.strike,
                                                                              dipangle=fault.dipangle,
                                                                              zerolon=inputs.zerolon,
                                                                              zerolat=inputs.zerolat,
                                                                              rake=fault.rake, top=zsplit_array[j],
                                                                              bottom=zsplit_array[j + 1],
                                                                              comment=fault.comment);
                    subfaulted_receivers.append(single_subfaulted_receiver);

    subfaulted_objects = cc.Input_object(PR1=inputs.PR1, FRIC=inputs.FRIC, depth=inputs.depth,
                                         start_gridx=inputs.start_gridx, finish_gridx=inputs.finish_gridx,
                                         start_gridy=inputs.start_gridy, finish_gridy=inputs.finish_gridy,
                                         xinc=inputs.xinc, yinc=inputs.yinc, minlon=inputs.minlon, maxlon=inputs.maxlon,
                                         zerolon=inputs.zerolon, minlat=inputs.minlat, maxlat=inputs.maxlat,
                                         zerolat=inputs.zerolat, source_object=inputs.source_object,
                                         receiver_object=subfaulted_receivers,
                                         receiver_horiz_profile=inputs.receiver_horiz_profile);

    return subfaulted_objects;


def get_split_x_y_arrays(start_x_top, finish_x_top, start_y_top, finish_y_top, strike_split):
    """
    Take the coordinates of the top of a receiver fault plane.
    Generate the list of coordinates that will help split it up along-strike
    strike_slip : int
    """
    if start_x_top == finish_x_top:
        xsplit_array = [start_x_top for _j in range(strike_split + 1)];
    else:
        xincrement = (finish_x_top - start_x_top) / strike_split;
        xsplit_array = np.arange(start_x_top, finish_x_top + xincrement, xincrement);
    # length : xsplit+1. contains all the xlocations that could be used as start-stop points in the top row.
    if start_y_top == finish_y_top:
        ysplit_array = [start_y_top for _j in range(strike_split + 1)];
    else:
        yincrement = (finish_y_top - start_y_top) / strike_split;
        ysplit_array = np.arange(start_y_top, finish_y_top + yincrement, yincrement);
    # length : xsplit+1. contains all the ylocations that could be used as start-stop points in the top row.
    return [xsplit_array, ysplit_array];


def get_split_z_array(top, bottom, dip_split):
    if top == bottom:
        zsplit_array = [top for _j in range(dip_split + 1)];
    else:
        zincrement = abs(top - bottom) / dip_split;
        zsplit_array = np.arange(top, bottom + zincrement, zincrement);
    return zsplit_array;


def compute_grid_def(inputs, params):
    """
    Loop through a grid and compute the displacements at each point from all sources put together.
    """
    x = np.linspace(inputs.start_gridx, inputs.finish_gridx,
                    int((inputs.finish_gridx - inputs.start_gridx) / inputs.xinc));
    y = np.linspace(inputs.start_gridy, inputs.finish_gridy,
                    int((inputs.finish_gridy - inputs.start_gridy) / inputs.yinc));
    [x2d, y2d] = np.meshgrid(x, y);
    u_displacements = np.zeros((len(y), len(x)));
    v_displacements = np.zeros((len(y), len(x)));
    w_displacements = np.zeros((len(y), len(x)));

    if not params.plot_grd_disp:
        return [x, y, x2d, y2d, u_displacements, v_displacements, w_displacements];

    print("Computing synthetic grid of displacements");
    numrows = np.shape(u_displacements)[0]
    numcols = np.shape(u_displacements)[1]
    rectangles, points, mogis = utilities.separate_source_types(inputs.source_object);
    fault_sources = rectangles + points;

    for ky in range(numrows):
        for kx in range(numcols):
            u_disp, v_disp, w_disp, _ = compute_surface_disp_point(fault_sources, params.alpha,
                                                                   x2d[ky][kx], y2d[ky][kx]);
            u_mogi, v_mogi, w_mogi = run_mogi.compute_surface_disp_point(mogis, params.nu, x2d[ky][kx], y2d[ky][kx]);
            u_displacements[ky][kx] = u_disp + u_mogi;
            v_displacements[ky][kx] = v_disp + v_mogi;
            w_displacements[ky][kx] = w_disp + w_mogi;
    return [x, y, x2d, y2d, u_displacements, v_displacements, w_displacements];


def compute_ll_strain(inputs, alpha, strain_points):
    """
    Loop through a list of lon/lat and compute their strains due to all sources put together.
    """
    if not strain_points:
        return [];
    x, y = [], [];
    print("Number of strain_points:", len(strain_points));
    for i in range(len(strain_points)):
        [xi, yi] = fault_vector_functions.latlon2xy(strain_points[i].lon, strain_points[i].lat, inputs.zerolon,
                                                    inputs.zerolat);
        x.append(xi);
        y.append(yi);

    strain_tensor_results = [];
    # For each coordinate requested.
    for k in range(len(x)):
        _, _, _, strain_tensor = compute_surface_disp_point(inputs.source_object, alpha, x[k], y[k]);
        strain_tensor_results.append(strain_tensor);

    return strain_tensor_results;


def compute_ll_def(inputs, alpha, disp_points):
    """
    Loop through a list of lon/lat and compute their displacements due to all sources put together.
    """
    if not disp_points:
        return [];
    model_disp_points = [];
    print("Number of disp_points:", len(disp_points));
    for point in disp_points:
        [xi, yi] = fault_vector_functions.latlon2xy(point.lon, point.lat, inputs.zerolon, inputs.zerolat);
        u_disp, v_disp, w_disp, _ = compute_surface_disp_point(inputs.source_object, alpha, xi, yi);
        model_point = Displacement_points(lon=point.lon, lat=point.lat, dE_obs=u_disp[0], dN_obs=v_disp[0],
                                          dU_obs=w_disp[0], Se_obs=0, Sn_obs=0, Su_obs=0, name=point.name);
        model_disp_points.append(model_point);
    return model_disp_points;


def compute_surface_disp_point(sources, alpha, x, y, compute_depth=0):
    """
    A major compute loop for each source object at one x/y point.
    x/y in the same coordinate system as the fault object. Computes displacement and strain tensor.

    :param sources: list of fault objects
    :param alpha: float
    :param x: float
    :param y: float
    :param compute_depth: depth of observation. Default depth is at surface of earth
    """
    u_disp, v_disp, w_disp = 0, 0, 0;
    strain_tensor_total = np.zeros((3, 3));

    for source in sources:
        desired_coords_grad_u, desired_coords_u = compute_strains_stresses_from_one_fault(source, x, y, compute_depth,
                                                                                          alpha);

        # Strain tensor math
        strain_tensor = conversion_math.get_strain_tensor(desired_coords_grad_u);
        strain_tensor_total = np.add(strain_tensor, strain_tensor_total);

        # Update the displacements from all sources
        u_disp = u_disp + desired_coords_u[0];
        v_disp = v_disp + desired_coords_u[1];
        w_disp = w_disp + desired_coords_u[2];  # vertical

    return u_disp, v_disp, w_disp, strain_tensor_total;


def compute_stresses_horiz_profile(params, inputs):
    """
    Pseudocode: Create a grid of points, and compute all stresses on them, given specified strike/dip/rake.
    horiz_profile = [depth_km, strike, dip, rake, centerlon, centerlat, length_km, width_km, inc_km];

    :param params: named tuple
    :param inputs: named tuple
    :returns: list of 3 lists, representing receiver normal, shear, and coulomb stress results
    """
    if not inputs.receiver_horiz_profile:
        return None;

    # Build a regular grid and iterate through.
    print("Resolving stresses on a horizontal profile.")
    profile = inputs.receiver_horiz_profile;
    receiver_normal, receiver_shear, receiver_coulomb = [], [], [];

    # perf improvement: Compute receiver geometry just once, since it's a profile of fixed geometry
    rec_strike_v, rec_dip_v, rec_plane_normal = conversion_math.get_geom_attributes_from_receiver_profile(profile);

    for i in range(len(inputs.receiver_horiz_profile.lon1d)):
        [xi, yi] = fault_vector_functions.latlon2xy(profile.lon1d[i], profile.lat1d[i], inputs.zerolon, inputs.zerolat);
        normal_sum, shear_sum, coulomb_sum = 0, 0, 0;
        for source in inputs.source_object:
            # A major compute loop for each source object.

            desired_coords_grad_u, _ = compute_strains_stresses_from_one_fault(source, xi, yi, profile.depth_km,
                                                                               params.alpha);

            # Then rotate again into receiver coordinates.
            strain_tensor = conversion_math.get_strain_tensor(desired_coords_grad_u);
            stress_tensor = conversion_math.get_stress_tensor(strain_tensor, params.lame1, params.mu);

            # Then compute shear, normal, and coulomb stresses.
            [normal, shear, coulomb] = conversion_math.get_coulomb_stresses_internal(stress_tensor, rec_strike_v,
                                                                                     profile.rake,
                                                                                     rec_dip_v, rec_plane_normal,
                                                                                     inputs.FRIC, params.B);
            normal_sum = normal_sum + normal;
            shear_sum = shear_sum + shear;
            coulomb_sum = coulomb_sum + coulomb;

        receiver_normal.append(normal_sum);
        receiver_shear.append(shear_sum);
        receiver_coulomb.append(coulomb_sum);

    return receiver_normal, receiver_shear, receiver_coulomb;


def compute_strains_stresses(params, inputs):
    """
    Pseudocode:
    For each receiver, at the center point, sum up the strain and stress for each source.
    Return : source object, receiver object, shear stress, normal stress, and coulomb stress on each receiver.
    """

    # The values we're actually going to output.
    receiver_shear, receiver_normal, receiver_coulomb = [], [], [];
    if not inputs.receiver_object:
        return [receiver_normal, receiver_shear, receiver_coulomb];

    if not params.plot_stress:
        return [receiver_normal, receiver_shear, receiver_coulomb];

    print("Resolving stresses on receiver fault(s).")
    for receiver in inputs.receiver_object:
        centercoords = conversion_math.get_fault_center(receiver);
        normal_sum, shear_sum, coulomb_sum = 0, 0, 0;

        for source in inputs.source_object:
            # A major compute loop for each source object.
            desired_coords_grad_u, _ = compute_strains_stresses_from_one_fault(source, centercoords[0],
                                                                               centercoords[1], centercoords[2],
                                                                               params.alpha);

            # Then rotate again into receiver coordinates.
            strain_tensor = conversion_math.get_strain_tensor(desired_coords_grad_u);
            stress_tensor = conversion_math.get_stress_tensor(strain_tensor, params.lame1, params.mu);

            # Then compute shear, normal, and coulomb stresses.
            [normal, shear, coulomb] = conversion_math.get_coulomb_stresses_internal(stress_tensor,
                                                                                     receiver.strike_unit_vector,
                                                                                     receiver.rake,
                                                                                     receiver.dip_unit_vector,
                                                                                     receiver.plane_normal,
                                                                                     inputs.FRIC, params.B);
            normal_sum = normal_sum + normal;
            shear_sum = shear_sum + shear;
            coulomb_sum = coulomb_sum + coulomb;

        receiver_normal.append(normal_sum);
        receiver_shear.append(shear_sum);
        receiver_coulomb.append(coulomb_sum);

    # return lists of normal, shear, coulomb values for each receiver.
    return [receiver_normal, receiver_shear, receiver_coulomb];


def compute_strains_stresses_from_one_fault(source, x, y, z, alpha):
    """
    The main math of DC3D
    Operates on a source object (e.g., fault),
    and an xyz position in the same cartesian reference frame.
    """
    R = source.R
    R2 = source.R2
    strike_slip = source.rtlat * -1;  # The dc3d coordinate system has left-lateral positive.

    # Compute the position relative to the translated, rotated fault.
    translated_pos = np.array(
        [[x - source.xstart], [y - source.ystart], [-z]]);
    xyz = R.dot(translated_pos);
    if source.potency:
        success, u, grad_u = dc3d0wrapper(alpha, [xyz[0], xyz[1], xyz[2]], source.top, source.dipangle,
                                          [source.potency[0], source.potency[1], source.potency[2],
                                           source.potency[3]]);
        grad_u = grad_u * 1e-9;  # DC3D0 Unit correction: potency from N-m results in strain in nanostrain
        u = u * 1e-6;  # Unit correction: potency from N-m results in displacements in microns.
    else:
        success, u, grad_u = dc3dwrapper(alpha, [xyz[0], xyz[1], xyz[2]], source.top, source.dipangle,
                                         [0, source.L], [-source.W, 0],
                                         [strike_slip, source.reverse, source.tensile]);
        grad_u = grad_u * 1e-3;  # DC3D Unit correction.
    # Solve for displacement gradients at certain xyz position

    # Rotate grad_u back into the unprimed coordinates.
    desired_coords_grad_u = np.dot(R2, np.dot(grad_u, R2.T));
    desired_coords_u = R2.dot(np.array([[u[0]], [u[1]], [u[2]]]));
    return desired_coords_grad_u, desired_coords_u;
