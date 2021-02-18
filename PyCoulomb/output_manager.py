# output_manager

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
from matplotlib.patches import Polygon
import pygmt
from subprocess import call
from . import coulomb_collections
from . import conversion_math
from . import io_inp
from . import io_additionals


def produce_outputs(params, inputs, disp_points, out_object):
    call(['mkdir', '-p', params.outdir], shell=False);
    call(['cp', params.config_file, params.outdir], shell=False);
    call(['cp', params.input_file, params.outdir], shell=False);

    subfaulted_inputs = coulomb_collections.Input_object(PR1=inputs.PR1, FRIC=inputs.FRIC, depth=inputs.depth,
                                                         start_gridx=inputs.start_gridx, start_gridy=inputs.start_gridy,
                                                         finish_gridx=inputs.finish_gridx,
                                                         finish_gridy=inputs.finish_gridy, xinc=inputs.xinc,
                                                         yinc=inputs.yinc, minlon=inputs.minlon, maxlon=inputs.maxlon,
                                                         zerolon=inputs.zerolon, minlat=inputs.minlat,
                                                         maxlat=inputs.maxlat, zerolat=inputs.zerolat,
                                                         source_object=out_object.source_object,
                                                         receiver_object=out_object.receiver_object);
    # make a new object of the subfaulted configuration.
    io_inp.write_inp(params.outdir + 'subfaulted.inp', subfaulted_inputs);
    surface_def_plot(params, out_object);  # THIS WORKS
    stress_plot(params, out_object, 'shear');  # can give vmin, vmax here if desired.
    stress_plot(params, out_object, 'normal');
    stress_plot(params, out_object, 'coulomb');
    map_plot(params, inputs, out_object, 'coulomb');
    map_plot(params, inputs, out_object, 'normal');
    map_plot(params, inputs, out_object, 'shear');
    write_output_files(params, inputs, disp_points, out_object);
    slip_vector_map(params, inputs, disp_points, out_object);
    # side_on_plot(params);
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


def stress_plot(params, out_object, stress_type, vmin="", vmax=""):
    # default vmin,vmax are in KPa
    # Here we will put plots of fault patches, colored by the magnitude of the stress component.

    print("Making plot of %s stress on receiver fault patches. " % stress_type);

    if stress_type == 'shear':
        stress_component = out_object.receiver_shear;
    elif stress_type == 'normal':
        stress_component = out_object.receiver_normal;
    elif stress_type == 'coulomb':
        stress_component = out_object.receiver_coulomb;
    else:
        print("Error! Invalid stress type : %s " % stress_type);

    # Select boundaries of color map.
    if vmin == "":
        vmin = np.min(stress_component);
    if vmax == "":
        vmax = np.max(stress_component);
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=True);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='RdYlBu_r');

    # Figure of stresses.
    plt.figure(figsize=(12, 10));

    for i in range(len(stress_component)):
        [x_total, y_total, x_updip, y_updip] = conversion_math.get_fault_four_corners(out_object.receiver_object[i]);
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


def side_on_plot(params):
    [x, y, z, rake, normal, shear, coulomb] = np.loadtxt(params.outdir + 'stresses.txt', skiprows=1, unpack=True)
    plt.figure(figsize=(10, 6));
    plotting_data = "Normal";
    # plotting_data="Coulomb";
    # plotting_data="Shear";
    vmin = -1;
    vmax = 1;  # kpa
    if plotting_data == "Coulomb":
        plt.scatter(x, z, c=coulomb, s=1450, marker='s', cmap='jet', edgecolor='black', vmin=vmin, vmax=vmax);
    elif plotting_data == "Normal":
        plt.scatter(x, z, c=normal, s=1450, marker='s', cmap='jet', edgecolor='black', vmin=vmin, vmax=vmax);
    else:
        plt.scatter(x, z, c=shear, s=1450, marker='s', cmap='jet', edgecolor='black', vmin=vmin, vmax=vmax);
    plt.ylim(z.min() - 2, z.max() + 2);
    plt.xlim(x.min() - 10, x.max() + 10);
    plt.xlabel('X axis (km)');
    plt.ylabel('Depth (km)')
    plt.gca().invert_yaxis();
    cb = plt.colorbar();
    if len(set(rake)) == 1:
        plt.title(plotting_data + ' stress change on fault planes, rake = %.1f (KPa)' % rake[0], fontsize=20);
    else:
        plt.title(plotting_data + ' stress change for variable rake (KPa)', fontsize=20);
    cb.set_label('Kilopascals', fontsize=18);
    plt.savefig(params.outdir + 'side_view.eps');
    plt.close();
    return;


def map_plot(params, inputs, out_object, stress_component):
    # Using PyGMT
    # Filling in fault patches with colors corresponding to their stress changes
    # Some options:
    if stress_component == 'shear':
        plotting_stress = out_object.receiver_shear;
        label = 'Shear';
    elif stress_component == 'normal':
        plotting_stress = out_object.receiver_normal;
        label = 'Normal';
    else:
        plotting_stress = out_object.receiver_coulomb;  # The default option
        label = 'Coulomb';

    # Make stress bounds for map.
    stress_bounds = [abs(np.min(plotting_stress)), abs(np.max(plotting_stress))];
    stress_bound = np.max(stress_bounds);  # setting the scale to symmetric about zero
    smallest_stress = -stress_bound;  # units: KPa
    largest_stress = stress_bound;  # units: KPa
    smallest_stress = -1;  # units: KPa
    largest_stress = 1;  # units: KPa

    # Make cpt
    pygmt.makecpt(C="jet", T=str(smallest_stress - 0.1) + "/" + str(largest_stress + 0.1) + "/0.05", H="mycpt.cpt",
                  D=True);

    # Make Map
    region = [inputs.minlon, inputs.maxlon, inputs.minlat, inputs.maxlat];
    proj = "M7i"
    fig = pygmt.Figure()
    title = "+t\"" + stress_component + " stress\"";  # must put escaped quotations around the title.
    fig.basemap(region=region, projection=proj, B=title);
    fig.coast(shorelines="1.0p,black", region=region, N="1", projection=proj, B="1.0");  # the boundary.
    fig.coast(region=region, projection=proj, N='2', W='0.5p,black', S='white', L="g-125.5/39.6+c1.5+w50");

    # Draw each source
    eq_lon = [];
    eq_lat = [];
    for source in out_object.source_object:
        source_lon, source_lat = conversion_math.xy2lonlat(source.xstart, source.ystart, inputs.zerolon,
                                                           inputs.zerolat);
        eq_lon.append(source_lon);
        eq_lat.append(source_lat);
        [x_total, y_total, x_updip, y_updip] = conversion_math.get_fault_four_corners(source);
        lons, lats = conversion_math.xy2lonlat(x_total, y_total, inputs.zerolon, inputs.zerolat);
        if source.potency == []:
            fig.plot(x=lons, y=lats, pen="thick,black");  # in case of area sources, just outline them.
        else:
            fig.plot(x=lons, y=lats, style='s0.3c', G="purple", pen="thin,black");  # in case of point sources

    # Draw each receiver, with associated data
    for i in range(len(out_object.receiver_object)):
        rec = out_object.receiver_object[i];
        [x_total, y_total, x_updip, y_updip] = conversion_math.get_fault_four_corners(rec);
        mylon_top, mylat_top = conversion_math.xy2lonlat(x_total[0], y_total[0], inputs.zerolon, inputs.zerolat);
        mylon_bot, mylat_bot = conversion_math.xy2lonlat(x_total[2], y_total[2], inputs.zerolon, inputs.zerolat);
        lons, lats = conversion_math.xy2lonlat(x_total, y_total, inputs.zerolon, inputs.zerolat);
        fig.plot(x=lons, y=lats, Z="f" + str(plotting_stress[i]), pen="thick,black",
                 C="mycpt.cpt");  # coloring by stress value

    # Colorbar annotation
    fig.coast(shorelines="1.0p,black", region=region, projection=proj);  # the boundary.
    fig.colorbar(D="jBr+w3.5i/0.2i+o2.5c/1.5c+h", C="mycpt.cpt", I="0.8",
                 G=str(smallest_stress) + "/" + str(largest_stress - 0.1), B=["x" + str(0.2), "y+L\"KPa\""]);

    # Annotate with earthquake location.
    fig.plot(eq_lon, eq_lat, style='s0.3c', G="purple", pen="thin,black");

    # Annotate with aftershock locations
    if len(params.aftershocks) > 0:
        [lon, lat, _, _, _] = io_additionals.read_aftershock_table(params.aftershocks);
        fig.plot(lon, lat, style='c0.1c', G='black', pen="thin,black");

    fig.savefig(params.outdir + label + '_map.png');
    plt.close();
    return;


def slip_vector_map(params, input_object, disp_points, out_object):
    # Here we will make a plot of vector displacement from the given GPS_LL file
    # We will include vertical deformation too.
    # The source faults will be colored by slip values.
    # Sources colored by slip would probably be great here.
    # Doing this in basic plotting until pygmt has good vector plotting utilities

    # Vertical cmap
    vmin = -10;
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
        [x_total, y_total, x_updip, y_updip] = conversion_math.get_fault_four_corners(source);
        lons, lats = conversion_math.xy2lonlat(x_total, y_total, input_object.zerolon, input_object.zerolat);
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


def write_output_files(params, inputs, disp_points, out_object):
    # Write displacement output file
    ofile = open(params.outdir + 'disps_model_grid.txt', 'w');
    ofile.write("# Format: x y udisp vdisp wdisp (m) \n");
    for i in np.arange(0, len(out_object.y)):
        for j in np.arange(0, len(out_object.x)):
            ofile.write("%f %f %f %f %f\n" % (
                out_object.x2d[i][j], out_object.y2d[i][j], out_object.u_disp[i][j], out_object.v_disp[i][j],
                out_object.w_disp[i][j]));
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
    print("Outputs written to file.")

    if disp_points:
        ofile = open(params.outdir + 'll_disps.txt', 'w');
        ofile.write("# Format: lon lat u v w (m)\n");
        for i in range(len(out_object.u_ll)):
            ofile.write("%f %f %f %f %f\n" % (
                disp_points.lon[i], disp_points.lat[i], out_object.u_ll[i], out_object.v_ll[i], out_object.w_ll[i]));
        ofile.close();
    return;
