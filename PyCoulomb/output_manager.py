# output_manager

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
from matplotlib.patches import Polygon
from subprocess import call
from . import coulomb_collections as cc
from . import conversion_math, io_inp, pygmt_plots, io_additionals
from . import fault_slip_object
from Tectonic_Utils.geodesy import fault_vector_functions


def produce_outputs(params, inputs, obs_disp_points, obs_strain_points, out_object):
    call(['mkdir', '-p', params.outdir], shell=False);
    call(['cp', params.config_file, params.outdir], shell=False);
    call(['cp', params.input_file, params.outdir], shell=False);
    write_subfaulted_inp(inputs, out_object, params.outdir+"subfaulted.inp");
    write_output_files(params, out_object, obs_strain_points);
    surface_def_plot(params, out_object);  # a grid of synthetic points
    # stress_plot(params, out_object, 'shear');  # can give vmin, vmax here if desired.
    # stress_plot(params, out_object, 'normal');
    stress_plot(params, out_object, 'coulomb');
    stress_cross_section_cartesian(params, out_object, 'coulomb');
    pygmt_plots.map_stress_plot(params, inputs, out_object, 'coulomb');
    # pygmt_plots.map_stress_plot(params, inputs, out_object, 'normal');
    # pygmt_plots.map_stress_plot(params, inputs, out_object, 'shear');
    pygmt_plots.map_displacement_vectors(params, inputs, obs_disp_points, out_object, params.outdir+"vector_plot.png")
    # pygmt_plots.map_vertical_def(params, inputs, params.outdir+"vert.grd", params.outdir+"vertical_map.png");
    return;


def produce_vmin_vmax_symmetric(plotting_array, vmin, vmax):
    """
    Determine boundaries of a symmetric colormap object (like stress change), coding all the edge-case logic here
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
                                        receiver_object=out_object.receiver_object);
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
    plt.title('Surface Dipslacement', fontsize=28)
    plt.savefig(params.outdir + "Displacement_model_on_grid.eps")
    plt.close();
    return;


def stress_plot(params, out_object, stress_type, vmin=None, vmax=None):
    """
    default vmin,vmax are in KPa
    Plots of fault patches in cartesian space, colored by the magnitude of the stress component.
    """
    print("Making plot of %s stress on receiver fault patches: %s. " % (stress_type, params.outdir +
                                                                        'Stresses_' + stress_type + '.eps'));

    if stress_type == 'shear':
        stress_component = out_object.receiver_shear;
    elif stress_type == 'normal':
        stress_component = out_object.receiver_normal;
    else:
        stress_component = out_object.receiver_coulomb;

    if not out_object.receiver_object:
        return;

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
    plt.savefig(params.outdir + 'Stresses_' + stress_type + '.eps');
    plt.close();
    return;


def stress_cross_section_cartesian(params, out_object, stress_type, vmin=None, vmax=None):
    """
    default vmin,vmax are in KPa
    Vertical plots of fault patches in cartesian space, colored by the magnitude of the stress component.
    """
    print("Making vertical plot of %s stress on receiver fault patches: %s. " % (stress_type, params.outdir +
                                                                                 'Stresses_cross_section_' +
                                                                                 stress_type + '.eps'));

    if stress_type == 'shear':
        stress_component = out_object.receiver_shear;
    elif stress_type == 'normal':
        stress_component = out_object.receiver_normal;
    else:
        stress_component = out_object.receiver_coulomb;

    # Select boundaries of color map, usually forcing even distribution around 0
    vmin, vmax = produce_vmin_vmax_symmetric(stress_component, vmin, vmax);
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=True);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='RdYlBu_r');

    # Figure of stresses.
    plt.figure(figsize=(17, 8));

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
        fault_vertices = np.column_stack((ycoords, zcoords));
        patch_color = custom_cmap.to_rgba(stress_component[i]);

        mypolygon = Polygon(fault_vertices, color=patch_color, alpha=1.0);
        plt.gca().add_patch(mypolygon);

    custom_cmap.set_array(np.arange(vmin, vmax, 100));
    cb = plt.colorbar(custom_cmap);
    cb.set_label('Kilopascals', fontsize=22);
    for label in cb.ax.yaxis.get_ticklabels():
        label.set_size(18)

    plt.grid();
    plt.gca().tick_params(axis='both', which='major', labelsize=18)
    plt.title(stress_type + ' stress from source faults', fontsize=22)
    plt.xlim([out_object.x.min(), out_object.x.max()])
    plt.xlabel('Distance (km)', fontsize=18);
    plt.ylabel('Depth (km)', fontsize=18);
    plt.ylim([0, 10]);
    plt.gca().invert_yaxis();
    plt.gca().invert_xaxis();
    plt.savefig(params.outdir + 'Stresses_cross_section_' + stress_type + '.eps');
    plt.close();
    return;


def write_synthetic_grid_triplets(x, y, x2d, y2d, zerolon, zerolat, u, v, w, outdir):
    """
    Write lists of lon/lat/def for each component of deformation in synthetic grid
    Used for GMT plots
    """
    ofile_w = open(outdir + 'xyz_model.txt', 'w');
    ofile_u = open(outdir + 'xyu_model.txt', 'w');
    ofile_v = open(outdir + 'xyv_model.txt', 'w');
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
    ofile = open(outdir + 'disps_model_grid.txt', 'w');
    ofile.write("# Format: x y lon lat udisp vdisp wdisp (m) \n");
    for i in np.arange(0, len(y)):
        for j in np.arange(0, len(x)):
            loni, lati = fault_vector_functions.xy2lonlat(x2d[i][j], y2d[i][j], zerolon, zerolat);
            ofile.write("%f %f %f %f %f %f %f\n" % (x2d[i][j], y2d[i][j], loni, lati, u[i][j], v[i][j], w[i][j]));
    ofile.close();
    return;


def write_output_files(params, out_object, obs_strain_points):
    # Write synethetic displacement output file and lists of lon/lat/def for synthetic grid
    write_synthetic_grid_full_results(out_object.x, out_object.y, out_object.x2d, out_object.y2d, out_object.zerolon,
                                      out_object.zerolat, out_object.u_disp, out_object.v_disp, out_object.w_disp,
                                      params.outdir);
    write_synthetic_grid_triplets(out_object.x, out_object.y, out_object.x2d, out_object.y2d, out_object.zerolon,
                                  out_object.zerolat, out_object.u_disp, out_object.v_disp, out_object.w_disp,
                                  params.outdir);

    # Write output file for stresses.
    if out_object.receiver_object:
        fault_dict_list = fault_slip_object.io_pycoulomb.coulomb_fault_to_fault_dict(out_object.receiver_object);
        fault_slip_object.io_slippy.write_stress_results_slippy_format(fault_dict_list, out_object.receiver_normal,
                                                                       out_object.receiver_shear,
                                                                       out_object.receiver_coulomb,
                                                                       params.outdir+'stresses_full.txt');

    # Write output files for GPS displacements and strains at specific lon/lat points (if used)
    io_additionals.write_disp_points_results(out_object.model_disp_points, params.outdir+"ll_disps.txt");
    io_additionals.write_strain_results(obs_strain_points, out_object.strains, params.outdir+'ll_strains.txt');
    return;
