# Running dc3d, given an input object and params object.

import numpy as np
from . import coulomb_collections as cc
from . import run_mogi, pyc_fault_object, conversion_math
from .fault_slip_triangle import triangle_okada
from .point_source_object import point_sources
from .disp_points_object.disp_points_object import Displacement_points
from Tectonic_Utils.geodesy import fault_vector_functions


def do_stress_computation(params, inputs, disp_points=(), strain_points=()):
    """
    * Step 0. Split receiver fault into many sub-faults if necessary
    * Step 1. Compute strains and displacements
    * Step 2. Resolve stresses on receiver faults
    """

    print("Beginning stress calculation.")
    print("Number of sources: %d " % len(inputs.source_object))
    print("Number of receivers: %d " % len(inputs.receiver_object))
    print("Coefficient of friction: %f" % inputs.FRIC)
    subfaulted_inputs = split_subfault_receivers(params, inputs)

    # Computes here.
    [x, y, x2d, y2d, u_disps, v_disps, w_disps] = compute_grid_def(subfaulted_inputs, params)
    model_disp_points = compute_ll_def(subfaulted_inputs, params, disp_points)
    strain_tensor_results = compute_ll_strain(inputs, params, strain_points)
    receiver_normal, receiver_shear, receiver_coulomb = compute_strains_stresses(params, subfaulted_inputs)
    receiver_profile_results = compute_stresses_horiz_profile(params, subfaulted_inputs)

    MyOutObject = cc.Out_object(x=x, y=y, x2d=x2d, y2d=y2d, u_disp=u_disps, v_disp=v_disps, w_disp=w_disps,
                                strains=strain_tensor_results, model_disp_points=model_disp_points,
                                zerolon=inputs.zerolon, zerolat=inputs.zerolat,
                                source_object=inputs.source_object, receiver_object=subfaulted_inputs.receiver_object,
                                receiver_normal=receiver_normal, receiver_shear=receiver_shear,
                                receiver_coulomb=receiver_coulomb, receiver_profile=receiver_profile_results)
    return MyOutObject


def split_subfault_receivers(params, inputs):
    strike_split = params.strike_num_receivers
    dip_split = params.dip_num_receivers

    if strike_split == 1 and dip_split == 1:
        # If we're not splitting the subfaults...
        subfaulted_receivers = inputs.receiver_object
        print("Not subdividing receiver faults further.")
    else:
        subfaulted_receivers = []
        print("Splitting %d receiver faults into %d subfaults each." % (
            len(inputs.receiver_object), strike_split * dip_split))

        for fault in inputs.receiver_object:  # for each receiver...
            # We find the depths corresponding to the tops and bottoms of our new sub-faults
            zsplit_array = get_split_z_array(fault.top, fault.bottom, dip_split)

            for j in range(dip_split):  # First we split it up by dip.
                # Get the new coordinates of the top of the fault plane.
                W = fault_vector_functions.get_downdip_width(fault.top, zsplit_array[j], fault.dipangle)
                vector_mag = W * np.cos(
                    np.deg2rad(fault.dipangle))  # how far the bottom edge is displaced downdip from map-view

                # Get the starting points for the next row of fault subpatches.
                [start_x_top, start_y_top] = fault_vector_functions.add_vector_to_point(fault.xstart, fault.ystart,
                                                                                        vector_mag, fault.strike + 90)
                [finish_x_top, finish_y_top] = fault_vector_functions.add_vector_to_point(fault.xfinish, fault.yfinish,
                                                                                          vector_mag, fault.strike +
                                                                                          90)

                [xsplit_array, ysplit_array] = get_split_x_y_arrays(start_x_top, finish_x_top, start_y_top,
                                                                    finish_y_top, strike_split)

                for k in range(strike_split):
                    single_subfaulted_receiver = pyc_fault_object.Faults_object(xstart=xsplit_array[k],
                                                                                xfinish=xsplit_array[k + 1],
                                                                                ystart=ysplit_array[k],
                                                                                yfinish=ysplit_array[k + 1],
                                                                                Kode=fault.Kode, strike=fault.strike,
                                                                                dipangle=fault.dipangle,
                                                                                zerolon=inputs.zerolon,
                                                                                zerolat=inputs.zerolat,
                                                                                rake=fault.rake, top=zsplit_array[j],
                                                                                bottom=zsplit_array[j + 1],
                                                                                comment=fault.comment)
                    subfaulted_receivers.append(single_subfaulted_receiver)

    subfaulted_objects = cc.Input_object(PR1=inputs.PR1, FRIC=inputs.FRIC, depth=inputs.depth,
                                         start_gridx=inputs.start_gridx, finish_gridx=inputs.finish_gridx,
                                         start_gridy=inputs.start_gridy, finish_gridy=inputs.finish_gridy,
                                         xinc=inputs.xinc, yinc=inputs.yinc, minlon=inputs.minlon, maxlon=inputs.maxlon,
                                         zerolon=inputs.zerolon, minlat=inputs.minlat, maxlat=inputs.maxlat,
                                         zerolat=inputs.zerolat, source_object=inputs.source_object,
                                         receiver_object=subfaulted_receivers,
                                         receiver_horiz_profile=inputs.receiver_horiz_profile)

    return subfaulted_objects


def get_split_x_y_arrays(start_x_top, finish_x_top, start_y_top, finish_y_top, strike_split):
    """
    Take the coordinates of the top of a receiver fault plane.
    Generate the list of coordinates that will help split it up along-strike
    strike_slip : int
    """
    if start_x_top == finish_x_top:
        xsplit_array = [start_x_top for _j in range(strike_split + 1)]
    else:
        xincrement = (finish_x_top - start_x_top) / strike_split
        xsplit_array = np.arange(start_x_top, finish_x_top + xincrement, xincrement)
    # length : xsplit+1. contains all the xlocations that could be used as start-stop points in the top row.
    if start_y_top == finish_y_top:
        ysplit_array = [start_y_top for _j in range(strike_split + 1)]
    else:
        yincrement = (finish_y_top - start_y_top) / strike_split
        ysplit_array = np.arange(start_y_top, finish_y_top + yincrement, yincrement)
    # length : xsplit+1. contains all the ylocations that could be used as start-stop points in the top row.
    return [xsplit_array, ysplit_array]


def get_split_z_array(top, bottom, dip_split):
    if top == bottom:
        zsplit_array = [top for _j in range(dip_split + 1)]
    else:
        zincrement = abs(top - bottom) / dip_split
        zsplit_array = np.arange(top, bottom + zincrement, zincrement)
    return zsplit_array


def compute_ll_strain(inputs, params, strain_points):
    """
    Loop through a list of lon/lat and compute their strains due to all sources put together.
    """
    # NEED TO ADD MOGI-SOURCE STRAINS
    strain_tensors1 = triangle_okada.compute_ll_strain_tris(inputs, params, strain_points)
    strain_tensors2 = point_sources.compute_ll_strain_point(inputs, params, strain_points)
    strain_tensors_total = [np.add(x, y) for x, y in zip(strain_tensors1, strain_tensors2)]
    return strain_tensors_total


def compute_ll_def(inputs, params, disp_points):
    """
    Loop through a list of lon/lat and compute their displacements due to all sources put together.
    """
    model_disp_points = triangle_okada.compute_ll_def_tris(inputs, params, disp_points)
    model_disp_points = run_mogi.compute_ll_def_mogi(inputs, params, model_disp_points)
    model_disp_points = point_sources.compute_ll_def_point(inputs, params, model_disp_points)
    return model_disp_points


def compute_xy_strain(inputs, params, strain_points):
    """ Loop through inputs and compute strain from all sources at given points, in cartesian coords.
    Strain_points is a list. """
    # NEED TO ADD MOGI-SOURCE STRAINS
    strain_tensors1 = triangle_okada.compute_cartesian_strain_tris(inputs, params, strain_points)
    strain_tensors2 = point_sources.compute_cartesian_strain_point(inputs, params, strain_points)
    strain_tensors_total = [np.add(x, y) for x, y in zip(strain_tensors1, strain_tensors2)]
    return strain_tensors_total


def compute_xy_def(inputs, params, disp_points):
    """ Loop through inputs and compute displacements from all sources at given points, in cartesian coords.
    Disp_points is a list."""
    modeled_disp_points = triangle_okada.compute_cartesian_def_tris(inputs, params, disp_points)
    modeled_disp_points = run_mogi.compute_cartesian_def_mogi(inputs, params, modeled_disp_points)
    modeled_disp_points = point_sources.compute_cartesian_def_point(inputs, params, modeled_disp_points)
    return modeled_disp_points


def compute_grid_def(inputs, params):
    """
    Loop through a grid and compute the displacements at each point from all sources put together.
    """
    x = np.linspace(inputs.start_gridx, inputs.finish_gridx,
                    int((inputs.finish_gridx - inputs.start_gridx) / inputs.xinc))
    y = np.linspace(inputs.start_gridy, inputs.finish_gridy,
                    int((inputs.finish_gridy - inputs.start_gridy) / inputs.yinc))
    [x2d, y2d] = np.meshgrid(x, y)
    u_disps, v_disps, w_disps = np.zeros(np.shape(x2d)), np.zeros(np.shape(x2d)), np.zeros(np.shape(x2d))

    if not params.plot_grd_disp:
        return [x, y, x2d, y2d, u_disps, v_disps, w_disps]

    print("Computing synthetic grid of displacements")
    numrows, numcols = np.shape(u_disps)

    disp_points = []
    for ky in range(numrows):
        for kx in range(numcols):
            disp_points.append(Displacement_points(lon=x2d[ky][kx], lat=y2d[ky][kx]))

    modeled_disps = compute_xy_def(inputs, params, disp_points)
    u_displacements = np.array([x.dE_obs for x in modeled_disps]).reshape(np.shape(u_disps))
    v_displacements = np.array([x.dN_obs for x in modeled_disps]).reshape(np.shape(v_disps))
    w_displacements = np.array([x.dU_obs for x in modeled_disps]).reshape(np.shape(w_disps))
    return [x, y, x2d, y2d, u_displacements, v_displacements, w_displacements]


def compute_stresses_horiz_profile(params, inputs):
    """
    Pseudocode: Create a grid of points, and compute all stresses on them, given specified strike/dip/rake.
    horiz_profile = [depth_km, strike, dip, rake, centerlon, centerlat, length_km, width_km, inc_km]

    :param params: named tuple
    :param inputs: named tuple
    :returns: list of 3 lists, representing receiver normal, shear, and coulomb stress results
    """
    if not inputs.receiver_horiz_profile:
        return None

    # Build a regular grid and iterate through.
    print("Resolving stresses on a horizontal profile.")
    profile = inputs.receiver_horiz_profile
    receiver_normal, receiver_shear, receiver_coulomb, strain_points = [], [], [], []

    # perf improvement: Compute receiver geometry just once, since it's a profile of fixed geometry
    rec_strike_v, rec_dip_v, rec_plane_normal = conversion_math.get_geom_attributes_from_receiver_profile(profile)

    for i in range(len(inputs.receiver_horiz_profile.lon1d)):
        [xi, yi] = fault_vector_functions.latlon2xy(profile.lon1d[i], profile.lat1d[i], inputs.zerolon, inputs.zerolat)
        strain_points.append(Displacement_points(lon=xi, lat=yi, depth=profile.depth_km))

    strain_tensors = compute_xy_strain(inputs, params, strain_points)

    for i in range(len(strain_tensors)):
        stress_tensor = conversion_math.get_stress_tensor(strain_tensors[i], params.lame1, params.mu)

        # Then compute shear, normal, and coulomb stresses.
        [normal, shear, coulomb] = conversion_math.get_coulomb_stresses_internal(stress_tensor, rec_strike_v,
                                                                                 profile.rake, rec_dip_v,
                                                                                 rec_plane_normal, inputs.FRIC,
                                                                                 params.B)

        receiver_normal.append(normal)
        receiver_shear.append(shear)
        receiver_coulomb.append(coulomb)

    return receiver_normal, receiver_shear, receiver_coulomb


def compute_strains_stresses(params, inputs):
    """
    Pseudocode:
    For each receiver, at the center point, sum up the strain and stress for each source.
    Return : source object, receiver object, shear stress, normal stress, and coulomb stress on each receiver.
    """

    # The values we're actually going to output.
    receiver_shear, receiver_normal, receiver_coulomb, target_points = [], [], [], []
    if not inputs.receiver_object:
        return [receiver_normal, receiver_shear, receiver_coulomb]
    if not params.plot_stress:
        return [receiver_normal, receiver_shear, receiver_coulomb]

    print("Resolving stresses on receiver fault(s).")
    for rec in inputs.receiver_object:
        centercoords = rec.get_fault_center()  # in cartesian coordinates
        strain_point = Displacement_points(lon=centercoords[0], lat=centercoords[1], depth=centercoords[2])
        target_points.append(strain_point)

    strain_tensors = compute_xy_strain(inputs, params, target_points)

    for rec, strain_tensor in zip(inputs.receiver_object, strain_tensors):
        stress_tensor = conversion_math.get_stress_tensor(strain_tensor, params.lame1, params.mu)

        # Then compute shear, normal, and coulomb stresses.
        [normal, shear, coulomb] = conversion_math.get_coulomb_stresses_internal(stress_tensor, rec.strike_unit_vector,
                                                                                 rec.rake, rec.dip_unit_vector,
                                                                                 rec.plane_normal, inputs.FRIC,
                                                                                 params.B)
        receiver_normal.append(normal)
        receiver_shear.append(shear)
        receiver_coulomb.append(coulomb)

    # return lists of normal, shear, coulomb values for each receiver.
    return receiver_normal, receiver_shear, receiver_coulomb
