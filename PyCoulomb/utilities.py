
import numpy as np
from subprocess import call


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
        cbar_options = [vmin, vmax, np.round((vmax-vmin)/8, 5)];  # color bars
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
