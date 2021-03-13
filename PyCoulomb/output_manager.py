# output_manager

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
from matplotlib.patches import Polygon
from subprocess import call
from . import coulomb_collections as cc
from . import conversion_math
from . import io_inp
from . import pygmt_plots
from Tectonic_Utils.geodesy import fault_vector_functions


def produce_outputs(params, inputs, disp_points, out_object):
    call(['mkdir', '-p', params.outdir], shell=False);
    call(['cp', params.config_file, params.outdir], shell=False);
    call(['cp', params.input_file, params.outdir], shell=False);
    write_subfaulted_inp(inputs, out_object, params.outdir+"subfaulted.inp");
    write_output_files(params, out_object, disp_points);
    surface_def_plot(params, out_object);  # a grid of synthetic points
    stress_plot(params, out_object, 'shear');  # can give vmin, vmax here if desired.
    stress_plot(params, out_object, 'normal');
    stress_plot(params, out_object, 'coulomb');
    pygmt_plots.map_stress_plot(params, inputs, out_object, 'coulomb');
    pygmt_plots.map_stress_plot(params, inputs, out_object, 'normal');
    pygmt_plots.map_stress_plot(params, inputs, out_object, 'shear');
    slip_vector_map(params, inputs, disp_points, out_object);
    pygmt_plots.map_vertical_def(params, inputs, params.outdir+"vert.grd", params.outdir+"vertical_map.png");
    return;


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
    print("Making plot of %s stress on receiver fault patches. " % stress_type);

    if stress_type == 'shear':
        stress_component = out_object.receiver_shear;
    elif stress_type == 'normal':
        stress_component = out_object.receiver_normal;
    else:
        stress_component = out_object.receiver_coulomb;

    # Select boundaries of color map.
    if not vmin:
        vmin = np.min(stress_component);
    if not vmax:
        vmax = np.max(stress_component);
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


def slip_vector_map(params, input_object, disp_points, out_object, vmin=None, vmax=None):
    """
    Here we will make a plot of vector displacement from the given GPS_LL file
    We will include vertical deformation too.
    The source faults will be colored by slip values.
    Sources colored by slip would probably be great here.
    Doing this in basic plotting until pygmt has good vector plotting utilities
    CAN I RE-WRITE THIS IN PYGMT NOW?
    """

    if not disp_points:
        return;

    # Vertical cmap (mm)
    if not vmin:
        vmin = -10;
    if not vmax:
        vmax = 10;
    color_boundary_object_vert = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=True);
    vertical_cmap = cm.ScalarMappable(norm=color_boundary_object_vert, cmap='RdYlBu_r');

    # slip cmap
    slip_total = [np.sqrt(source.rtlat ** 2 + source.reverse ** 2) for source in input_object.source_object];
    slipmin = 0
    slipmax = np.max(slip_total) + 0.1;
    color_boundary_object_slip = matplotlib.colors.Normalize(vmin=slipmin, vmax=slipmax, clip=True);
    slip_cmap = cm.ScalarMappable(norm=color_boundary_object_slip, cmap='jet');

    # MAKING THE FITURE
    plt.figure(figsize=(16, 10), dpi=300);
    lonW = input_object.minlon;
    lonE = input_object.maxlon;
    latS = input_object.minlat;
    latN = input_object.maxlat;

    # Drawing Sources
    for source in input_object.source_object:
        slip = np.sqrt(source.rtlat ** 2 + source.reverse ** 2);
        [x_total, y_total, _, _] = conversion_math.get_fault_four_corners(source);
        lons, lats = fault_vector_functions.xy2lonlat(x_total, y_total, input_object.zerolon, input_object.zerolat);
        fault_vertices = np.column_stack((lons[0:4], lats[0:4]));
        patch_color = slip_cmap.to_rgba(slip);
        mypolygon = Polygon(fault_vertices, color=patch_color, alpha=1.0);
        plt.gca().add_patch(mypolygon);
        # plt.plot(lons, lats, linewidth=1,color='purple');
        # plt.plot(lons[0:2], lats[0:2], linewidth=3,color='purple');  # the top edge of the fault patch
        if source.potency:
            plt.plot(lons[0], lats[0], marker='s', markersize=10, color='purple');  # the location of focal mechanism

    # Displacement vectors at GPS stations
    # Smaller scale here means zooming in on the vectors (appropriate for small events)
    max_vector = np.nanmax(
        np.sqrt(np.multiply(out_object.u_ll, out_object.u_ll) + np.multiply(out_object.v_ll, out_object.v_ll)));
    scale = max_vector * 3000;  # 2000 gives a reasonable sized vector for 1m of slip
    if disp_points:
        for i in range(len(disp_points.lon)):
            patch_color = vertical_cmap.to_rgba(1000 * out_object.w_ll[i]);
            plt.plot(disp_points.lon[i], disp_points.lat[i], marker='o', markersize=15, markeredgecolor='black',
                     color=patch_color);
            plt.quiver(disp_points.lon[i], disp_points.lat[i], 1000 * out_object.u_ll[i], 1000 * out_object.v_ll[i],
                       scale=scale, color='black', zorder=10);
            if len(disp_points.lon) == len(disp_points.dE_obs):
                plt.quiver(disp_points.lon[i], disp_points.lat[i], 1000 * disp_points.dE_obs[i],
                           1000 * disp_points.dN_obs[i], scale=scale, color='red', zorder=10);
        # plt.text(disp_points[0][i],disp_points[1][i],disp_points[2][i]);  # If you want to label the GPS stations
        plt.quiver(lonW + 0.02, latS + 0.03, 20.0, 0.0, scale=scale, color='black');
        plt.text(lonW + 0.02, latS + 0.05, "20mm model", color="black", fontsize=20);

    # Drawing colorbars
    cbar = plt.colorbar(vertical_cmap);
    cbar.set_label("Vertical Motion (mm)", fontsize=16);
    cbar.ax.tick_params(labelsize=14);
    cbar2 = plt.colorbar(slip_cmap);
    cbar2.set_label("Fault Slip (m)", fontsize=16);
    cbar2.ax.tick_params(labelsize=14);

    plt.xlim([lonW, lonE]);
    plt.ylim([latS, latN]);
    plt.xlabel("Longitude", fontsize=24);
    plt.ylabel("Latitude", fontsize=24);
    plt.gca().tick_params(axis='both', which='major', labelsize=24)
    plt.savefig(params.outdir + 'slip_vector_map.png');
    plt.close();
    return;


def write_output_files(params, out_object, disp_points):
    # Write synethetic displacement output file
    ofile = open(params.outdir + 'disps_model_grid.txt', 'w');
    ofile.write("# Format: x y lon lat udisp vdisp wdisp (m) \n");
    for i in np.arange(0, len(out_object.y)):
        for j in np.arange(0, len(out_object.x)):
            loni, lati = fault_vector_functions.xy2lonlat(out_object.x2d[i][j], out_object.y2d[i][j],
                                                          out_object.zerolon, out_object.zerolat);
            ofile.write("%f %f %f %f %f %f %f\n" % (
                out_object.x2d[i][j], out_object.y2d[i][j], loni, lati, out_object.u_disp[i][j],
                out_object.v_disp[i][j], out_object.w_disp[i][j]));
    ofile.close();

    # Write lists of lon/lat/def for each component of deformation
    ofile = open(params.outdir + 'xyz_model.txt', 'w');
    for i in np.arange(0, len(out_object.y)):
        for j in np.arange(0, len(out_object.x)):
            loni, lati = fault_vector_functions.xy2lonlat(out_object.x2d[i][j], out_object.y2d[i][j],
                                                          out_object.zerolon, out_object.zerolat);
            ofile.write("%f %f %f\n" % (loni, lati, out_object.w_disp[i][j]));
    ofile.close();
    ofile = open(params.outdir + 'xyu_model.txt', 'w');
    for i in np.arange(0, len(out_object.y)):
        for j in np.arange(0, len(out_object.x)):
            loni, lati = fault_vector_functions.xy2lonlat(out_object.x2d[i][j], out_object.y2d[i][j],
                                                          out_object.zerolon, out_object.zerolat);
            ofile.write("%f %f %f\n" % (loni, lati, out_object.u_disp[i][j]));
    ofile.close();
    ofile = open(params.outdir + 'xyv_model.txt', 'w');
    for i in np.arange(0, len(out_object.y)):
        for j in np.arange(0, len(out_object.x)):
            loni, lati = fault_vector_functions.xy2lonlat(out_object.x2d[i][j], out_object.y2d[i][j],
                                                          out_object.zerolon, out_object.zerolat);
            ofile.write("%f %f %f\n" % (loni, lati, out_object.v_disp[i][j]));
    ofile.close();

    # Write output file for stresses.
    ofile = open(params.outdir + 'stresses.txt', 'w');
    ofile.write("# Format: centerx centery centerz rake normal shear coulomb (kpa)\n");
    for i in range(len(out_object.receiver_object)):
        rec = out_object.receiver_object[i];
        # [x_total, y_total, x_updip, y_updip] = conversion_math.get_fault_four_corners(rec);
        center = conversion_math.get_fault_center(rec);
        ofile.write("%f %f %f %f %f %f %f \n" % (
            center[0], center[1], center[2], rec.rake, out_object.receiver_normal[i], out_object.receiver_shear[i],
            out_object.receiver_coulomb[i]));
    ofile.close();
    print("Write file %s " % params.outdir+"stresses.txt");

    if disp_points:
        ofile = open(params.outdir + 'll_disps.txt', 'w');
        ofile.write("# Format: lon lat u v w (m)\n");
        for i in range(len(out_object.u_ll)):
            ofile.write("%f %f %f %f %f\n" % (
                disp_points.lon[i], disp_points.lat[i], out_object.u_ll[i], out_object.v_ll[i], out_object.w_ll[i]));
        ofile.close();
    return;
