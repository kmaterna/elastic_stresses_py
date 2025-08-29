# Running dc3d, given an input object and params object.

import numpy as np
from . import coulomb_collections as cc
from . import run_mogi, conversion_math
from .fault_slip_triangle import triangle_okada
from .point_source_object import point_sources
from .disp_points_object.disp_points_object import Displacement_points
from tectonic_utils.geodesy import fault_vector_functions
from . import utilities


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
    subfaulted_inputs = split_subfault_receivers(inputs, params.strike_num_receivers, params.dip_num_receivers)

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


def split_subfault_receivers(inputs, strike_split, dip_split):
    if strike_split == 1 and dip_split == 1:
        # If we're not splitting the subfaults...
        subfaulted_receivers = inputs.receiver_object
        print("Not subdividing receiver faults further.")
    else:
        subfaulted_receivers = []
        print("Splitting %d receiver faults into %d subfaults each." % (
            len(inputs.receiver_object), strike_split * dip_split))

        for fault in inputs.receiver_object:  # for each receiver...
            subfaults = fault.split_single_fault(strike_split, dip_split, inputs.zerolon, inputs.zerolat)
            subfaulted_receivers.extend(subfaults)  # unpack into new list, equivalent to L1 += L2

    subfaulted_objects = cc.Input_object(PR1=inputs.PR1, FRIC=inputs.FRIC, depth=inputs.depth,
                                         start_gridx=inputs.start_gridx, finish_gridx=inputs.finish_gridx,
                                         start_gridy=inputs.start_gridy, finish_gridy=inputs.finish_gridy,
                                         xinc=inputs.xinc, yinc=inputs.yinc, minlon=inputs.minlon, maxlon=inputs.maxlon,
                                         zerolon=inputs.zerolon, minlat=inputs.minlat, maxlat=inputs.maxlat,
                                         zerolat=inputs.zerolat, source_object=inputs.source_object,
                                         receiver_object=subfaulted_receivers,
                                         receiver_horiz_profile=inputs.receiver_horiz_profile)
    return subfaulted_objects


def compute_ll_strain(inputs, params, strain_points):
    """
    Loop through a list of lon/lat and compute their strains due to all sources put together.
    NEED TO ADD MOGI-SOURCE STRAINS
    """
    strain_points = utilities.convert_ll2xy_disp_points(strain_points, inputs.zerolon, inputs.zerolat)
    strain_tensors_total = compute_xy_strain(inputs, params, strain_points)
    return strain_tensors_total


def compute_ll_def(inputs, params, disp_points):
    """
    Loop through a list of lon/lat and compute their displacements due to all sources put together.
    """
    cart_disp_points = utilities.convert_ll2xy_disp_points(disp_points, inputs.zerolon, inputs.zerolat)  # ll to cart.
    model_disp_points = compute_xy_def(inputs, params, cart_disp_points)
    model_disp_points = utilities.transform_disp_points_ll_by_key_array(model_disp_points, disp_points)  # back to ll
    return model_disp_points


def compute_xy_strain(inputs, params, strain_points):
    """ Loop through inputs and compute strain from all sources at given points, in cartesian coords.
    Strain_points is a list. # NEED TO ADD MOGI-SOURCE STRAINS  """
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

    :param params: object of type configure_calc.Params
    :param inputs: object of type Inputs_object
    :returns: list of 5 lists, representing receiver normal, shear, and coulomb stress results,
     and a list of the full stress tensors and full strain tensors
    """
    if not inputs.receiver_horiz_profile:
        return None

    # Build a regular grid and iterate through.
    print("Resolving stresses on a horizontal profile with points of shape ", inputs.receiver_horiz_profile.shape)
    profile = inputs.receiver_horiz_profile
    receiver_normal = np.empty(np.shape(profile.lon1d))
    receiver_shear = np.empty(np.shape(profile.lon1d))
    receiver_coulomb = np.empty(np.shape(profile.lon1d))
    sigmaij, strain_points = [], []

    # perf improvement: Compute receiver geometry just once, since it's a profile of fixed geometry
    rec_strike_v, rec_dip_v, rec_plane_normal = conversion_math.get_geom_attributes_from_receiver_profile(profile)

    [xi, yi] = fault_vector_functions.latlon2xy(profile.lon1d, profile.lat1d, inputs.zerolon, inputs.zerolat)
    for (x, y) in zip(xi, yi):
        strain_points.append(Displacement_points(lon=x, lat=y, depth=profile.depth_km))

    strain_tensors = np.array(compute_xy_strain(inputs, params, strain_points))  # list of strain tensors
    stress_tensors = conversion_math.stress_tensor_batch(strain_tensors, params.lame1, params.mu)

    for i, stress_tensor in enumerate(stress_tensors):
        # Then compute shear, normal, and coulomb stresses.
        [normal, shear, coulomb] = conversion_math.get_coulomb_stresses_internal(stress_tensor, rec_strike_v,
                                                                                 profile.rake, rec_dip_v,
                                                                                 rec_plane_normal, inputs.FRIC,
                                                                                 params.B)
        receiver_normal[i] = normal
        receiver_shear[i] = shear
        receiver_coulomb[i] = coulomb
        sigmaij.append(stress_tensor)  # optionally, send out to GRD files if rec_full_stress_tensor is set. In Pa.

    return receiver_normal, receiver_shear, receiver_coulomb, sigmaij, strain_tensors


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

    strain_tensors = np.array(compute_xy_strain(inputs, params, target_points))
    stress_tensors = conversion_math.stress_tensor_batch(strain_tensors, params.lame1, params.mu)

    for rec, stress_tensor in zip(inputs.receiver_object, stress_tensors):
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
