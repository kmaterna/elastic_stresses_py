import numpy as np
from subprocess import call
import os
from tectonic_utils.seismo import moment_calculations
from . import fault_slip_object as fso
from . import pyc_fault_object as pycfaults
from . import coulomb_collections as cc
from .disp_points_object.disp_points_object import Displacement_points
from tectonic_utils.geodesy import fault_vector_functions


def define_colorbar_series(plotting_array, vmin=None, vmax=None, tol=0.0005, v_labeling_interval=None):
    """
    Create the min, max, and intervals/labels for a colormap and its associated colorbar.
    Various scenarios of input plotting_arrays are coded here.
    The goal is to avoid issue of "z_high larger than highest z"

    :param plotting_array: 1d array that is being plotted with pygmt colorscale
    :param vmin: imposed lower bound of color scale
    :param vmax: imposed upper bound of color scale
    :param tol: lower tolerance on edges of color scale
    :param v_labeling_interval: hard-code the label interval on the colorbar
    :return cmap_options: [float, float, float] for vmin, vmax, and cmap_interval
    :return cbar_options: [float, float, float] for gmin, gmax, and cbar_interval
    """
    if vmin and vmax:   # easy case where vmin and vmax are specified
        total_interval = vmax - vmin  # guaranteed positive
        bumper = total_interval*0.07
        if bumper < tol:
            bumper = tol
        cmap_options = [vmin-bumper, vmax+bumper, total_interval/100]  # color maps
        cbar_options = [vmin, vmax, np.round((vmax-vmin)/8, 3)]  # color bars
        return [cmap_options, cbar_options]

    if plotting_array is None:    # if no plotting array is passed
        cmap_options, cbar_options = [], []

    elif len(plotting_array) == 0:   # if an empty list
        cmap_options, cbar_options = [], []

    elif len(plotting_array) == 1:   # in the case of a single value sent in
        total_range = 2*tol
        cmap_options = [plotting_array[0]-tol, plotting_array[0]+tol, total_range/10]
        cbar_options = [plotting_array[0]-tol/2, plotting_array[0]+tol/2, total_range/8]

    else:   # in the case of an array, even with the same values
        total_range = np.nanmax(plotting_array) - np.nanmin(plotting_array)
        if total_range < tol:  # in the case of an array made of practically the same values
            total_range = 10*tol
            bumper = total_range * 0.05
            if v_labeling_interval:
                label_int = v_labeling_interval  # using default value
            else:
                label_int = 0.001
            cmap_interval = total_range/100
            cmap_options = [np.nanmin(plotting_array) - bumper, np.nanmax(plotting_array) + bumper, cmap_interval]
            cbar_options = [np.nanmin(plotting_array) - bumper/2, np.nanmax(plotting_array) + bumper / 2, label_int]

        else:   # in the case of a normal array with multiple distinct values
            bumper = total_range*0.05
            if bumper < tol:
                bumper = tol
            if v_labeling_interval:
                label_int = v_labeling_interval  # using default value
            else:
                label_int = np.round(total_range/8, 5)
            cmap_interval = total_range/100
            cmap_options = [np.nanmin(plotting_array)-bumper, np.nanmax(plotting_array)+bumper, cmap_interval]
            cbar_options = [np.nanmin(plotting_array)-bumper/2, np.nanmax(plotting_array)+bumper/2, label_int]

    return [cmap_options, cbar_options]


def produce_vmin_vmax_symmetric(plotting_array, vmin, vmax):
    """
    Determine boundaries of a symmetric colormap object (like stress change), coding all the edge-case logic here

    :param plotting_array: 1d array or None.
    :param vmin: float
    :param vmax: float
    """
    if not plotting_array:
        return -1, 1
    auto_vmin = np.nanmin(plotting_array)  # one number
    auto_vmax = np.nanmax(plotting_array)  # one number
    auto_extreme = np.nanmax(np.abs([auto_vmin, auto_vmax]))
    auto_bounds = [-auto_extreme, auto_extreme]
    if not vmin:
        vmin = auto_bounds[0]
    if not vmax:
        vmax = auto_bounds[1]
    return vmin, vmax


def define_vector_scale_size(model_dE, model_dN):
    """
    Based on modeled displacements, determine an appropriate vector scale for map visualization

    :param model_dE: 1d-array
    :param model_dN: 1d-array
    :returns: float for arrow scale, string for vector text
    """
    max_disp = np.max(np.sqrt(np.square(model_dN) + np.square(model_dE)))
    if max_disp > 0.5:
        scale_arrow = 0.500
        vectext = "50 cm"
    elif max_disp > 0.2:
        scale_arrow = 0.200
        vectext = "20 cm"
    elif max_disp > 0.10:
        scale_arrow = 0.100
        vectext = "10 cm"  # 10 cm, large vectors
    elif max_disp > 0.050:
        scale_arrow = 0.05
        vectext = "5 cm"
    elif max_disp > 0.02:
        scale_arrow = 0.02
        vectext = "20 mm"
    elif max_disp > 0.01:
        scale_arrow = 0.01
        vectext = "10 mm"   # medium vectors
    elif max_disp > 0.005:
        scale_arrow = 0.005
        vectext = "5 mm"
    elif max_disp > 0.002:
        scale_arrow = 0.002
        vectext = "2 mm"  # small vectors
    else:
        scale_arrow = 0.001
        vectext = "1 mm"
    return scale_arrow, vectext


def displacements_to_3_grds(outdir, efiles, nfiles, ufiles, region, inc=0.0005):
    """
    Call gmt surface on each component. efiles, nfiles, and ufiles are tuples of inputs and outputs: (txtfile, grdfile)
    """
    call_gmt_surface(os.path.join(outdir, ufiles[0]), os.path.join(outdir, ufiles[1]), region, inc=inc)
    print("Printing "+ufiles[1])
    call_gmt_surface(os.path.join(outdir, efiles[0]), os.path.join(outdir, efiles[1]), region, inc=inc)
    print("Printing " + efiles[1])
    call_gmt_surface(os.path.join(outdir, nfiles[0]), os.path.join(outdir, nfiles[1]), region, inc=inc)
    print("Printing " + nfiles[1])
    return


def call_gmt_surface(xyzfile, outfile, region, inc):
    """
    Create a grd file from an xyz text file

    :param xyzfile: string
    :param outfile: string
    :param region: list of 4 floats
    :param inc: float
    """
    call(['gmt', 'surface', xyzfile, '-G' + outfile,
          '-R' + str(region[0]) + '/' + str(region[1]) + '/' + str(region[2]) + '/' + str(region[3]), '-I'+str(inc),
          '-rp'], shell=False)
    return


def print_metrics_on_sources(source_object, mu):
    """
    Print overall magnitude of rectangular source objects.
    """
    print("Number of total sources:", len(source_object))
    rect_sources, pt_sources, mogi_sources = separate_source_types(source_object)
    if len(rect_sources) > 0:
        mw = moment_calculations.mw_from_moment(pycfaults.get_faults_slip_moment(rect_sources, mu))
        print("Number of rectangular sources: %d" % len(rect_sources))
        print("Moment Magnitude from Rectangular Fault Patches (assuming G=%.1fGPa): %f" % (mu/1e9, mw))
    if len(pt_sources) > 0:
        print("Number of point sources: %d" % len(pt_sources))
    if len(mogi_sources) > 0:
        print("Number of Mogi sources: %d" % len(mogi_sources))
    return


def separate_source_types(source_object):
    """
    Take a list of source objects and separate it into lists:
    one of rectangular sources, one of point sources, one of mogi sources
    """
    rect_sources, pt_sources, mogi_sources = [], [], []
    for source in source_object:
        if isinstance(source, cc.Mogi_Source):
            mogi_sources.append(source)
        if isinstance(source, cc.Faults_object):
            if source.is_point_source:
                pt_sources.append(source)
            else:
                rect_sources.append(source)
    return rect_sources, pt_sources, mogi_sources


def convert_ll2xy_disp_points(disp_points, zerolon, zerolat):
    """
    :param disp_points: list of disp_points
    :param zerolon: longitude for center of coordinate system
    :param zerolat: latitude for center of coordinate system
    :return: list of disp_points relative to center of coordinate system. All lon/lat/depths are in km.
    """
    cartesian_disp_points = []
    for point in disp_points:
        [xi, yi] = fault_vector_functions.latlon2xy(point.lon, point.lat, zerolon, zerolat)
        model_point = Displacement_points(lon=xi, lat=yi, depth=point.depth, dE_obs=point.dE_obs, dN_obs=point.dN_obs,
                                          dU_obs=point.dU_obs, name=point.name)
        cartesian_disp_points.append(model_point)
    return cartesian_disp_points


def transform_disp_points_ll_by_key_array(cart_disp_points, ll_disp_points):
    result_disp_points = []
    for cart, ll in zip(cart_disp_points, ll_disp_points):
        new_disp_point = Displacement_points(lon=ll.lon, lat=ll.lat, depth=ll.depth, dE_obs=cart.dE_obs,
                                             dN_obs=cart.dN_obs, dU_obs=cart.dU_obs, Se_obs=cart.Se_obs,
                                             Sn_obs=cart.Sn_obs, Su_obs=cart.Su_obs, name=cart.name)
        result_disp_points.append(new_disp_point)
    return result_disp_points


def get_zeros_disp_points(disp_points):
    """
    :param disp_points: list of disp points
    :return: matching list of disp points where all data values will be set to zero
    """
    zero_disp_points = []
    for item in disp_points:
        new_pt = Displacement_points(lon=item.lon, lat=item.lat, depth=item.depth, name=item.name)
        zero_disp_points.append(new_pt)
    return zero_disp_points


def get_zeros_strain_points(strain_points):
    """
    :param strain_points: list of disp points
    :return: list of 3x3 zero matrices, for holding null strain tensor results
    """
    zero_strains = []
    for _x in strain_points:
        new_strain_tensor = np.zeros((3, 3))
        zero_strains.append(new_strain_tensor)
    return zero_strains


def write_fault_edges_to_gmt_file(fault_object, outfile='tmp.txt', color_array=lambda x: x.get_total_slip()):
    """
    This speeds up pygmt plotting.  Color_array can be a function of x, or a numpy array
    """
    rect_sources, _, _ = separate_source_types(fault_object)
    fault_dict_list = fso.fault_slip_object.coulomb_fault_to_fault_object(rect_sources)
    fso.file_io.outputs.write_gmt_fault_file(fault_dict_list, outfile, color_mappable=color_array, verbose=False)
    return
