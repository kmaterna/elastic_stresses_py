
import numpy as np
from subprocess import call
from Tectonic_Utils.seismo import moment_calculations
from Tectonic_Utils.geodesy import fault_vector_functions
import Elastic_stresses_py.PyCoulomb.fault_slip_object as fso


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
        total_interval = vmax - vmin;  # guaranteed positive
        bumper = total_interval*0.07;
        if bumper < tol:
            bumper = tol;
        cmap_options = [vmin-bumper, vmax+bumper, total_interval/100];  # color maps
        cbar_options = [vmin, vmax, np.round((vmax-vmin)/8, 3)];  # color bars
        return [cmap_options, cbar_options];

    if plotting_array is None:    # if no plotting array is passed
        cmap_options, cbar_options = [], [];

    elif len(plotting_array) == 0:   # if an empty list
        cmap_options, cbar_options = [], [];

    elif len(plotting_array) == 1:   # in the case of a single value sent in
        total_range = 2*tol;
        cmap_options = [plotting_array[0]-tol, plotting_array[0]+tol, total_range/10];
        cbar_options = [plotting_array[0]-tol/2, plotting_array[0]+tol/2, total_range/8];

    else:   # in the case of an array, even with the same values
        total_range = np.nanmax(plotting_array) - np.nanmin(plotting_array);
        if total_range < tol:  # in the case of an array made of practically the same values
            total_range = 10*tol;
            bumper = total_range * 0.05;
            if v_labeling_interval:
                label_int = v_labeling_interval;  # using default value
            else:
                label_int = 0.001;
            cmap_interval = total_range/100;
            cmap_options = [np.nanmin(plotting_array) - bumper, np.nanmax(plotting_array) + bumper, cmap_interval];
            cbar_options = [np.nanmin(plotting_array) - bumper/2, np.nanmax(plotting_array) + bumper / 2, label_int];

        else:   # in the case of a normal array with multiple distinct values
            bumper = total_range*0.05;
            if bumper < tol:
                bumper = tol;
            if v_labeling_interval:
                label_int = v_labeling_interval;  # using default value
            else:
                label_int = np.round(total_range/8, 5);
            cmap_interval = total_range/100;
            cmap_options = [np.nanmin(plotting_array)-bumper, np.nanmax(plotting_array)+bumper, cmap_interval];
            cbar_options = [np.nanmin(plotting_array)-bumper/2, np.nanmax(plotting_array)+bumper/2, label_int];

    return [cmap_options, cbar_options];


def produce_vmin_vmax_symmetric(plotting_array, vmin, vmax):
    """
    Determine boundaries of a symmetric colormap object (like stress change), coding all the edge-case logic here
    plotting_array: 1d array or None.
    """
    if not plotting_array:
        return -1, 1;
    auto_vmin = np.min(plotting_array);  # one number
    auto_vmax = np.max(plotting_array);  # one number
    auto_extreme = np.max(np.abs([auto_vmin, auto_vmax]));
    auto_bounds = [-auto_extreme, auto_extreme];
    if not vmin:
        vmin = auto_bounds[0];
    if not vmax:
        vmax = auto_bounds[1];
    return vmin, vmax;


def define_vector_scale_size(model_dE, model_dN):
    """
    Based on the modeled displacements, determine an appropriate vector scale for map visualization
    """
    max_disp = np.max(np.sqrt(np.square(model_dN) + np.square(model_dE)));
    if max_disp > 0.5:
        scale_arrow = 0.500;  vectext = "50 cm";
    elif max_disp > 0.2:
        scale_arrow = 0.200;  vectext = "20 cm";
    elif max_disp > 0.10:
        scale_arrow = 0.100;   vectext = "10 cm";  # 10 cm, large vectors
    elif max_disp > 0.050:
        scale_arrow = 0.05;   vectext = "5 cm";
    elif max_disp > 0.02:
        scale_arrow = 0.02;  vectext = "20 mm"
    elif max_disp > 0.01:
        scale_arrow = 0.01;  vectext = "10 mm"   # medium vectors
    elif max_disp > 0.005:
        scale_arrow = 0.005;  vectext = "5 mm"
    elif max_disp > 0.002:
        scale_arrow = 0.002;  vectext = "2 mm"  # small vectors
    else:
        scale_arrow = 0.001;  vectext = "1 mm"
    return scale_arrow, vectext;


def define_map_region(inputs):
    """
    Define the bounding box for map in pygmt [W, E, S, N] based on sources and receivers, if bigger than coord system
    """
    region = [inputs.minlon, inputs.maxlon, inputs.minlat, inputs.maxlat];
    allfaults = inputs.receiver_object + inputs.source_object
    for item in allfaults:
        lon, lat = fault_vector_functions.xy2lonlat(item.xstart, item.ystart, inputs.zerolon, inputs.zerolat);
        if lon < region[0]:
            region[0] = lon;
        if lon > region[1]:
            region[1] = lon;
        if lat < region[2]:
            region[2] = lat;
        if lat > region[3]:
            region[3] = lat;
    return region;


def displacements_to_3_grds(outdir, efiles, nfiles, ufiles, region, inc=0.0005):
    """
    Call gmt surface on each component. efiles, nfiles, and ufiles are tuples of inputs and outputs: (txtfile, grdfile)
    """
    call_gmt_surface(outdir+'/'+ufiles[0], outdir+'/'+ufiles[1], region, inc=inc);
    print("Printing "+ufiles[1]);
    call_gmt_surface(outdir+'/'+efiles[0], outdir+'/'+efiles[1], region, inc=inc);
    print("Printing " + efiles[1]);
    call_gmt_surface(outdir+'/'+nfiles[0], outdir+'/'+nfiles[1], region, inc=inc);
    print("Printing " + nfiles[1]);
    return;


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
          '-r'], shell=False);
    return;


def print_metrics_on_sources(source_object, mu):
    """
    Print overall magnitude of rectangular source objects.
    """
    print("Number of sources:", len(source_object));
    rect_sources, point_sources = separate_rectangular_from_point_sources(source_object);
    if len(point_sources) > 0:
        print("Number of point sources: %d" % len(point_sources) );
    if len(rect_sources) > 0:
        fault_dict_list = fso.io_pycoulomb.coulomb_fault_to_fault_dict(rect_sources);
        # default mu = 30GPa
        Mw = moment_calculations.mw_from_moment(fso.fault_slip_object.get_total_moment(fault_dict_list, mu=mu));
        print("Moment Magnitude from Rectangular Fault Patches (assuming G=%.1fGPa): %f" % (mu/1e9, Mw));
    return;


def separate_rectangular_from_point_sources(source_object):
    """
    Take a list of source objects and separate it into two lists: one of rectangular sources, one of point sources
    """
    point_sources, rect_sources = [], [];
    for source in source_object:
        if source.potency:
            point_sources.append(source);
        else:
            rect_sources.append(source);
    return rect_sources, point_sources;


def write_fault_edges_to_gmt_file(fault_object, outfile='tmp.txt', color_array=None):
    """
    This speeds up pygmt plotting
    """
    rect_sources, _ = separate_rectangular_from_point_sources(fault_object);
    fault_dict_list = fso.io_pycoulomb.coulomb_fault_to_fault_dict(fault_object);
    if color_array is None:
        color_mappable = fso.fault_slip_object.get_total_slip;  # put the total slip if nothing is provided
    else:
        color_mappable = color_array;
    fso.fault_slip_object.write_gmt_fault_file(fault_dict_list, outfile, color_mappable=color_mappable, verbose=False);
    return;


def check_each_fault_has_same_coord_system(list_of_fault_objects, system_zerolon, system_zerolat):
    """
    Check that all faults have the same coord system
    """
    tol = 0.00001
    for item in list_of_fault_objects:
        if np.abs(item.zerolon-system_zerolon) > tol:
            raise(ValueError, "input or receiver faults lack a good longitude coordinate system.");
        if np.abs(item.zerolat-system_zerolat) > tol:
            raise(ValueError, "input or receiver faults lack a good latitude coordinate system.");
    return;
