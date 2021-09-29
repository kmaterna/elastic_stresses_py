
import numpy as np


def define_colorbar_series(plotting_array, vmin=None, vmax=None, tol=0.0005, v_labeling_interval=None):
    """
    Create the min, max, and intervals/labels for a colormap and its associated colorbar.
    Various scenarios of input plotting_arrays are coded here.
    The goal is to avoid issue of "z_high larger than highest z"

    :param plotting_array: array that is being plotted with pygmt colorscale
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
        cmap_options = [plotting_array[0]-tol, plotting_array[0]+tol, total_range/100];
        cbar_options = [plotting_array[0], plotting_array[0], total_range/8];

    else:   # in the case of an array, even with the same values
        total_range = np.nanmax(plotting_array) - np.nanmin(plotting_array);
        if total_range < tol:  # in the case of an array made of practically the same values
            total_range = 10*tol;
            bumper = total_range * 0.05;
            if v_labeling_interval:
                label_int = v_labeling_interval;  # using default value
            else:
                label_int = 0.001;
            cmap_options = [np.nanmin(plotting_array) - bumper, np.nanmax(plotting_array) + bumper, total_range/100];
            cbar_options = [np.nanmin(plotting_array) - bumper/2, np.nanmax(plotting_array) + bumper / 2, label_int];

        else:   # in the case of a normal array with multiple distinct values
            bumper = total_range*0.05;
            if bumper < tol:
                bumper = tol;
            if v_labeling_interval:
                label_int = v_labeling_interval;  # using default value
            else:
                label_int = np.round(total_range/8, 5);
            cmap_options = [np.nanmin(plotting_array)-bumper, np.nanmax(plotting_array)+bumper, total_range/100];
            cbar_options = [np.nanmin(plotting_array)-bumper/2, np.nanmax(plotting_array)+bumper/2, label_int];

    return [cmap_options, cbar_options];
