# output_manager

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
from matplotlib.patches import Polygon
from subprocess import call
from . import coulomb_collections as cc
from . import conversion_math, io_inp, pygmt_plots, io_additionals, utilities, configure_calc
from .fault_slip_object import io_pycoulomb
from .fault_slip_object import io_slippy
from Tectonic_Utils.geodesy import fault_vector_functions


def produce_outputs(params, inputs, obs_disp_points, obs_strain_points, out_object):
    """
    Passes the program into the right plotting and mapping routines.
    Reads parameter flags such as plot_stress and plot_grd_disp.
    """
    call(['mkdir', '-p', params.outdir], shell=False);
    if params.config_file:
        configure_calc.write_params_into_config(params, params.outdir+"used_config.txt");  # for record-keeping
    if params.input_file:
        call(['cp', params.input_file, params.outdir], shell=False);  # for record-keeping
    write_output_files(params, out_object, obs_strain_points);
    write_subfaulted_inp(inputs, out_object, params.outdir+"subfaulted.inp");
    pygmt_plots.map_displacement_vectors(params, inputs, obs_disp_points, out_object.model_disp_points,
                                         params.outdir+"vector_plot.png");  # map point displacements
    if params.plot_stress:
        stress_plot(params, out_object, 'shear');  # can give vmin, vmax here if desired.
        stress_plot(params, out_object, 'normal');
        stress_plot(params, out_object, 'coulomb');
        pygmt_plots.map_stress_plot(params, inputs, out_object, 'coulomb');
        pygmt_plots.map_stress_plot(params, inputs, out_object, 'normal');
        pygmt_plots.map_stress_plot(params, inputs, out_object, 'shear');
        stress_cross_section_cartesian(params, out_object, 'coulomb', writefile=params.outdir+'coulomb_xsection.txt');
        stress_cross_section_cartesian(params, out_object, 'normal');
        stress_cross_section_cartesian(params, out_object, 'shear');
    if params.plot_grd_disp:  # create grd files and plot vertical. Can take a while.
        surface_def_plot(params, out_object);  # grid of synthetic points in cartesian space
        write_disp_grd_files(params, inputs);  # based on txt files already written
        pygmt_plots.map_vertical_def(params, inputs, params.outdir+"vertical_map.png");
    if out_object.receiver_profile:
        write_horiz_profile(params, inputs.receiver_horiz_profile, out_object.receiver_profile);
        map_horiz_profile(params, inputs.receiver_horiz_profile, out_object.receiver_profile);
    return;


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


def write_subfaulted_inp(inputs, out_object, outfile):
    subfaulted_inputs = cc.Input_object(PR1=inputs.PR1, FRIC=inputs.FRIC, depth=inputs.depth,
                                        start_gridx=inputs.start_gridx, start_gridy=inputs.start_gridy,
                                        finish_gridx=inputs.finish_gridx, finish_gridy=inputs.finish_gridy,
                                        xinc=inputs.xinc, yinc=inputs.yinc, minlon=inputs.minlon, maxlon=inputs.maxlon,
                                        zerolon=inputs.zerolon, minlat=inputs.minlat, maxlat=inputs.maxlat,
                                        zerolat=inputs.zerolat, source_object=out_object.source_object,
                                        receiver_object=out_object.receiver_object, receiver_horiz_profile=None);
    # write an inp for the subfaulted configuration.
    io_inp.write_inp(subfaulted_inputs, outfile);
    return;


def surface_def_plot(params, out_object):
    """
    Plots surface displacements on the synthetic cartesian grid domain
    """
    print("Making plot of predicted displacement throughout model domain.")
    print("Max vertical displacement is %f m" % (out_object.w_disp.max()));
    # Plot of elastic surface deformation from the given input model.
    plt.figure(figsize=(16, 16))
    plt.pcolormesh(out_object.x2d, out_object.y2d, out_object.w_disp, cmap='jet');
    cb = plt.colorbar();
    cb.set_label('Vertical displacement (meters)', fontsize=22);
    for label in cb.ax.yaxis.get_ticklabels():
        label.set_size(18)
    for i in np.arange(0, len(out_object.y), 5):
        for j in np.arange(0, len(out_object.x), 5):
            plt.quiver(out_object.x2d[i][j], out_object.y2d[i][j], out_object.u_disp[i][j], out_object.v_disp[i][j],
                       units='width', scale=0.2)
    for source in out_object.source_object:
        [x_total, y_total, x_updip, y_updip] = conversion_math.get_fault_four_corners(source);
        plt.plot(x_total, y_total, 'k', linewidth=1);
        plt.plot(x_updip, y_updip, 'g', linewidth=3);
        center = conversion_math.get_fault_center(source);
        plt.plot(center[0], center[1], '.g', markersize=8);
    for rec in out_object.receiver_object:
        [x_total, y_total, x_updip, y_updip] = conversion_math.get_fault_four_corners(rec);
        plt.plot(x_total, y_total, 'b', linewidth=1);
        plt.plot(x_updip, y_updip, 'b', linewidth=3);
        center = conversion_math.get_fault_center(rec);
        plt.plot(center[0], center[1], '.b', markersize=8);
    plt.xlim([out_object.x.min(), out_object.x.max()])
    plt.ylim([out_object.y.min(), out_object.y.max()])
    plt.gca().tick_params(axis='both', which='major', labelsize=18)
    plt.grid();
    plt.axis('equal');
    plt.xlabel('X (km)', fontsize=20)
    plt.ylabel('Y (km)', fontsize=20)
    plt.gca().tick_params(labelsize=16)
    plt.title('Surface Dipslacement', fontsize=28)
    plt.savefig(params.outdir + "Displacement_model_on_grid.png")
    plt.close();
    return;


def stress_plot(params, out_object, stress_type, vmin=None, vmax=None):
    """
    default vmin,vmax are in KPa
    Plots of fault patches in cartesian space, colored by the magnitude of the stress component.
    """

    if not out_object.receiver_object:
        return;

    print("Making plot of %s stress on receiver fault patches: %s. " % (stress_type, params.outdir +
                                                                        'Stresses_' + stress_type + '.png'));

    if stress_type == 'shear':
        stress_component = out_object.receiver_shear;
    elif stress_type == 'normal':
        stress_component = out_object.receiver_normal;
    else:
        stress_component = out_object.receiver_coulomb;

    # Select boundaries of color map.
    vmin, vmax = produce_vmin_vmax_symmetric(stress_component, vmin, vmax);
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=True);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='RdYlBu_r');

    # Figure of stresses.
    plt.figure(figsize=(12, 10));

    for i in range(len(stress_component)):
        [x_total, y_total, _, _] = conversion_math.get_fault_four_corners(out_object.receiver_object[i]);
        xcoords = x_total[0:4];
        ycoords = y_total[0:4];
        fault_vertices = np.column_stack((xcoords, ycoords));
        patch_color = custom_cmap.to_rgba(stress_component[i]);

        mypolygon = Polygon(fault_vertices, color=patch_color, alpha=1.0);
        plt.gca().add_patch(mypolygon);

    custom_cmap.set_array(np.arange(vmin, vmax, 100));
    cb = plt.colorbar(custom_cmap);
    cb.set_label('Kilopascals', fontsize=22);
    for label in cb.ax.yaxis.get_ticklabels():
        label.set_size(18)

    # Adding source and receiver faults
    for source in out_object.source_object:
        [x_total, y_total, x_updip, y_updip] = conversion_math.get_fault_four_corners(source);
        plt.plot(x_total, y_total, 'k', linewidth=1);
        plt.plot(x_updip, y_updip, 'g', linewidth=3);
        center = conversion_math.get_fault_center(source);
        plt.plot(center[0], center[1], '.g', markersize=8);
    for rec in out_object.receiver_object:
        [x_total, y_total, x_updip, y_updip] = conversion_math.get_fault_four_corners(rec);
        plt.plot(x_total, y_total, 'b', linewidth=1);
        plt.plot(x_updip, y_updip, 'b', linewidth=3);
        center = conversion_math.get_fault_center(rec);
        plt.plot(center[0], center[1], '.b', markersize=8);

    plt.grid();
    plt.axis('equal');
    plt.gca().tick_params(axis='both', which='major', labelsize=18)
    plt.title(stress_type + ' stress from source faults', fontsize=22)
    plt.xlim([out_object.x.min(), out_object.x.max()])
    plt.ylim([out_object.y.min(), out_object.y.max()])
    plt.xlabel('X (km)', fontsize=20);
    plt.ylabel('Y (km)', fontsize=20);
    plt.gca().tick_params(labelsize=16)
    plt.savefig(params.outdir + 'Stresses_' + stress_type + '.png');
    plt.close();
    return;


def stress_cross_section_cartesian(params, out_object, stress_type, vmin=None, vmax=None, writefile=None):
    """
    default vmin,vmax are in KPa
    Vertical plots of fault patches in cartesian space, colored by the magnitude of the stress component.
    """

    if not out_object.receiver_object:
        return;

    print("Making plot of %s stress on receiver fault patches: %s. " % (stress_type, params.outdir +
                                                                        'Stresses_cross_section_' +
                                                                        stress_type + '.png'));

    if stress_type == 'shear':
        stress_component = out_object.receiver_shear;
    elif stress_type == 'normal':
        stress_component = out_object.receiver_normal;
    else:
        stress_component = out_object.receiver_coulomb;

    max_depth_array = [x.bottom for x in out_object.receiver_object];
    min_depth_array = [x.top for x in out_object.receiver_object];

    # Select boundaries of color map, usually forcing even distribution around 0
    vmin, vmax = produce_vmin_vmax_symmetric(stress_component, vmin, vmax);
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=True);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='RdYlBu_r');

    # Figure of stresses.
    plt.figure(figsize=(17, 8));
    if writefile:
        print("Writing file %s" % writefile);
        ofile = open(writefile, 'w');
        ofile.close();

    total_y_array = [];
    for i in range(len(stress_component)):
        xcoords = [out_object.receiver_object[i].xstart, out_object.receiver_object[i].xstart,
                   out_object.receiver_object[i].xfinish, out_object.receiver_object[i].xfinish];
        ycoords = [out_object.receiver_object[i].ystart, out_object.receiver_object[i].ystart,
                   out_object.receiver_object[i].yfinish, out_object.receiver_object[i].yfinish];
        zcoords = [out_object.receiver_object[i].top, out_object.receiver_object[i].bottom,
                   out_object.receiver_object[i].bottom, out_object.receiver_object[i].top];
        fault_strike = out_object.receiver_object[i].strike
        x1, y1 = conversion_math.rotate_points(xcoords[0], ycoords[0], fault_strike);
        x2, y2 = conversion_math.rotate_points(xcoords[1], ycoords[1], fault_strike);
        x3, y3 = conversion_math.rotate_points(xcoords[2], ycoords[2], fault_strike);
        x4, y4 = conversion_math.rotate_points(xcoords[3], ycoords[3], fault_strike);
        ycoords = [y1, y2, y3, y4];
        total_y_array = total_y_array + ycoords;
        fault_vertices = np.column_stack((ycoords, zcoords));
        patch_color = custom_cmap.to_rgba(stress_component[i]);

        mypolygon = Polygon(fault_vertices, color=patch_color, alpha=1.0);
        plt.gca().add_patch(mypolygon);

        if writefile:   # the x-axis here is an arbitrary position along the orientation of the receiver fault
            ofile = open(writefile, 'a');
            ofile.write("> -Z%f\n" % stress_component[i] );
            ofile.write("%f %f \n" % (-y1, zcoords[0]))
            ofile.write("%f %f \n" % (-y2, zcoords[1]))
            ofile.write("%f %f \n" % (-y3, zcoords[2]))
            ofile.write("%f %f \n" % (-y4, zcoords[3]))
            ofile.write("%f %f \n" % (-y1, zcoords[0]))
            ofile.close();

    custom_cmap.set_array(np.arange(vmin, vmax, 100));
    cb = plt.colorbar(custom_cmap);
    cb.set_label('Kilopascals', fontsize=22);
    for label in cb.ax.yaxis.get_ticklabels():
        label.set_size(18)

    plt.grid();
    plt.gca().tick_params(axis='both', which='major', labelsize=18)
    plt.title(stress_type + ' stress from source faults', fontsize=22)
    plt.xlim([np.min(total_y_array)-1, np.max(total_y_array)+1]);
    plt.xlabel('Distance along fault profile (km)', fontsize=18);
    plt.ylabel('Depth (km)', fontsize=18);
    plt.ylim([min(min_depth_array), max(max_depth_array)]);
    plt.gca().invert_yaxis();
    plt.gca().invert_xaxis();
    plt.savefig(params.outdir + 'Stresses_cross_section_' + stress_type + '.png');
    plt.close();
    return;

def map_horiz_profile(params, horiz_profile, profile_results):
    """Display a small map of a horizontal profile of stresses. Default colors for now."""

    X = np.reshape(horiz_profile.lon1d, horiz_profile.shape);
    Y = np.reshape(horiz_profile.lat1d, horiz_profile.shape);
    _normal_stress = np.reshape(profile_results[0], horiz_profile.shape);
    _shear_stress = np.reshape(profile_results[1], horiz_profile.shape);
    coulomb_stress = np.reshape(profile_results[2], horiz_profile.shape);

    # Figure of stresses.
    plt.figure(figsize=(17, 8));
    plt.contourf(X, Y, coulomb_stress, cmap='RdYlBu_r');
    plt.title('Coulomb stresses on horizontal profile, fixed strike/dip/rake/depth of '+str(horiz_profile.strike)+', ' +
              str(horiz_profile.dip)+', '+str(horiz_profile.rake)+', '+str(horiz_profile.depth_km));

    cb = plt.colorbar();
    cb.set_label('Kilopascals', fontsize=22);
    for label in cb.ax.yaxis.get_ticklabels():
        label.set_size(18)
    plt.xlabel('Longitude', fontsize=18)
    plt.ylabel('Latitude', fontsize=18)
    plt.gca().tick_params(labelsize=16)
    plt.savefig(params.outdir+'horizontal_profile_stresses.png');
    return;


def write_synthetic_grid_triplets(x, y, x2d, y2d, zerolon, zerolat, u, v, w, outdir):
    """
    Write lists of lon/lat/def for each component of deformation in synthetic grid
    Used for GMT plots
    """
    ofile_w = open(outdir + 'xy_vert_model.txt', 'w');
    ofile_u = open(outdir + 'xy_east_model.txt', 'w');
    ofile_v = open(outdir + 'xy_north_model.txt', 'w');
    for i in np.arange(0, len(y)):
        for j in np.arange(0, len(x)):
            loni, lati = fault_vector_functions.xy2lonlat(x2d[i][j], y2d[i][j], zerolon, zerolat);
            ofile_w.write("%f %f %f\n" % (loni, lati, w[i][j]));
            ofile_u.write("%f %f %f\n" % (loni, lati, u[i][j]));
            ofile_v.write("%f %f %f\n" % (loni, lati, v[i][j]));
    ofile_w.close();
    ofile_u.close();
    ofile_v.close();
    return;


def write_synthetic_grid_full_results(x, y, x2d, y2d, zerolon, zerolat, u, v, w, outdir):
    outfile = outdir + 'disps_model_grid.txt';
    print("Writing synthetic grid of displacements in %s " % outfile)
    ofile = open(outfile, 'w');
    ofile.write("# Format: x y lon lat x_disp[m] y_disp[m] z_disp[m] \n");
    for i in np.arange(0, len(y)):
        for j in np.arange(0, len(x)):
            loni, lati = fault_vector_functions.xy2lonlat(x2d[i][j], y2d[i][j], zerolon, zerolat);
            ofile.write("%f %f %f %f %f %f %f\n" % (x2d[i][j], y2d[i][j], loni, lati, u[i][j], v[i][j], w[i][j]));
    ofile.close();
    return;


def write_output_files(params, out_object, obs_strain_points):
    # Write synethetic displacement output file and lists of lon/lat/def for synthetic grid
    if params.plot_grd_disp:
        write_synthetic_grid_full_results(out_object.x, out_object.y, out_object.x2d, out_object.y2d,
                                          out_object.zerolon,out_object.zerolat, out_object.u_disp, out_object.v_disp,
                                          out_object.w_disp, params.outdir);
        write_synthetic_grid_triplets(out_object.x, out_object.y, out_object.x2d, out_object.y2d, out_object.zerolon,
                                      out_object.zerolat, out_object.u_disp, out_object.v_disp, out_object.w_disp,
                                      params.outdir);

    # Write output file for stresses.
    if out_object.receiver_object:
        fault_dict_list = io_pycoulomb.coulomb_fault_to_fault_dict(out_object.receiver_object);
        io_slippy.write_stress_results_slippy_format(fault_dict_list, out_object.receiver_shear,
                                                     out_object.receiver_normal, out_object.receiver_coulomb,
                                                     params.outdir+'stresses_full.txt');

    # Write output files for GPS displacements and strains at specific lon/lat points (if used)
    io_additionals.write_disp_points_results(out_object.model_disp_points, params.outdir+"ll_disps.txt");
    io_additionals.write_strain_results(obs_strain_points, out_object.strains, params.outdir+'ll_strains.txt');
    io_additionals.write_receiver_traces_gmt(out_object.receiver_object, params.outdir+"receiver_traces.txt");
    return;

def write_horiz_profile(params, horiz_profile, profile_results):
    print("Writing %s " % params.outdir+"stresses_horiz_profile.txt");
    ofile = open(params.outdir+"stresses_horiz_profile.txt", 'w');
    ofile.write("# lon lat depth_km normal_kPa shear_kPa coulomb_kPa\n");
    ofile.write("# strike %f, dip %f, rake %f\n" % (horiz_profile.strike, horiz_profile.dip, horiz_profile.rake) );
    for i in range(len(horiz_profile.lon1d)):
        ofile.write("%f %f %f %f %f %f\n" % (horiz_profile.lon1d[i], horiz_profile.lat1d[i], horiz_profile.depth_km,
                                             profile_results[0][i], profile_results[1][i], profile_results[2][i]) );
    return;


def write_disp_grd_files(params, inputs):
    # Make surfaces of east/north/up deformation for plotting
    region = [inputs.minlon, inputs.maxlon, inputs.minlat, inputs.maxlat];
    inc = (region[1]-region[0]) / 500;  # for gmt mapping, default is 500 points per axis
    utilities.displacements_to_3_grds(params.outdir,
                                      efiles=('xy_east_model.txt', 'east.grd'),
                                      nfiles=('xy_north_model.txt', 'north.grd'),
                                      ufiles=('xy_vert_model.txt', 'vert.grd'), region=region, inc=inc)
    return;
