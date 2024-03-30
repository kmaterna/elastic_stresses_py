# output_manager
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import shutil
import matplotlib.cm as cm
from matplotlib.patches import Polygon
from . import conversion_math, pygmt_plots, io_additionals, utilities
from .inputs_object import io_intxt, io_inp
from . import fault_slip_object as fso
from Tectonic_Utils.geodesy import fault_vector_functions


def produce_outputs(params, inputs, obs_disp_points, obs_strain_points, out_object):
    """
    Passes the program into the right plotting and mapping routines.
    Reads parameter flags such as plot_stress and plot_grd_disp.
    """
    if not os.path.exists(params.outdir):
        os.makedirs(params.outdir)
    # if params.config_file:
    params.write_params_into_config(os.path.join(params.outdir, "used_config.txt"))  # for record-keeping
    io_intxt.write_intxt(inputs, os.path.join(params.outdir, "used_inputs.txt"), mu=params.mu, lame1=params.lame1)
    if params.input_file:
        shutil.copy(params.input_file, params.outdir)  # for record-keeping

    # Write output files for GPS displacements and strains at specific lon/lat points (if used)
    io_additionals.write_disp_points_results(out_object.model_disp_points, os.path.join(params.outdir, "ll_disps.txt"))
    io_additionals.write_strain_results(obs_strain_points, out_object.strains,
                                        os.path.join(params.outdir, 'll_strains.txt'))
    io_additionals.write_stress_results(obs_strain_points, out_object.strains, params.lame1, params.mu,
                                        os.path.join(params.outdir, 'll_stresses.txt'))
    io_additionals.write_fault_traces_gmt(out_object.receiver_object,
                                          os.path.join(params.outdir, "receiver_traces.txt"))
    io_additionals.write_fault_traces_gmt(out_object.source_object, os.path.join(params.outdir, "source_traces.txt"))
    write_subfaulted_inp(inputs, out_object, os.path.join(params.outdir, "subfaulted.inp"))
    pygmt_plots.map_displacement_vectors(params, inputs, obs_disp_points, out_object.model_disp_points,
                                         os.path.join(params.outdir, "vector_plot.png"))  # map point displacements
    if params.plot_stress:  # write the outputs of stress calculation, if doing a stress calculation
        fault_dict_list = fso.fault_slip_object.coulomb_fault_to_fault_object(out_object.receiver_object)
        fso.file_io.io_slippy.write_stress_results_slippy_format(fault_dict_list, out_object.receiver_shear,
                                                                 out_object.receiver_normal,
                                                                 out_object.receiver_coulomb,
                                                                 os.path.join(params.outdir, 'stresses_full.txt'))
        stress_plot(params, out_object, stress_type='Shear')  # can give vmin, vmax here if desired.
        stress_plot(params, out_object, stress_type='Normal')
        stress_plot(params, out_object, stress_type='Coulomb')
        pygmt_plots.map_stress_plot(params, inputs, out_object, stress_type='Coulomb')
        pygmt_plots.map_stress_plot(params, inputs, out_object, stress_type='Normal')
        pygmt_plots.map_stress_plot(params, inputs, out_object, stress_type='Shear')
        stress_cross_section_cartesian(params, out_object, stress_type='Coulomb',
                                       writefile=os.path.join(params.outdir, 'coulomb_xsection.txt'))
        stress_cross_section_cartesian(params, out_object, stress_type='Normal')
        stress_cross_section_cartesian(params, out_object, stress_type='Shear')
    if params.plot_grd_disp:  # create synthetic grid outputs, grd files, and plot vertical.
        write_synthetic_grid_full_results(out_object, os.path.join(params.outdir, 'disps_model_grid.txt'))
        surface_def_plot(out_object, os.path.join(params.outdir, "Displacement_model_on_grid.png"))  # synthetic grid
        write_synthetic_grid_triplets(out_object, params.outdir, 'xy_east.txt', 'xy_north.txt', 'xy_vert.txt')
        write_disp_grd_files(inputs, params.outdir, 'xy_east.txt', 'xy_north.txt', 'xy_vert.txt')  # from txt files
        pygmt_plots.map_vertical_def(params, inputs, os.path.join(params.outdir, "vertical_map.png"))
    if out_object.receiver_profile:
        write_horiz_profile(inputs.receiver_horiz_profile, out_object.receiver_profile,
                            os.path.join(params.outdir, "stresses_horiz_profile.txt"))
        map_horiz_profile(inputs.receiver_horiz_profile, out_object.receiver_profile,
                          os.path.join(params.outdir, 'horizontal_profile_stresses.png'))
    return


def write_subfaulted_inp(inputs, out_object, outfile):
    # write an inp for the sub-faulted configuration.
    subfaulted_inputs = inputs.modify_inputs_object(source_object=out_object.source_object,
                                                    receiver_object=out_object.receiver_object)
    io_inp.write_inp(subfaulted_inputs, outfile)
    return


def surface_def_plot(out_object, outfile):
    """
    Plots surface displacements on the synthetic cartesian grid domain
    """
    print("Making plot of predicted displacement throughout model domain.")
    print("Max vertical displacement is %f m" % (out_object.w_disp.max()))
    # Plot of elastic surface deformation from the given input model.
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    fig = plt.figure(figsize=(16, 16), dpi=300)
    display_map = plt.pcolormesh(out_object.x2d, out_object.y2d, out_object.w_disp, cmap='jet')
    cb = fig.colorbar(display_map, ax=plt.gca(), location='right')
    cb.set_label('Vertical displacement (meters)', fontsize=22)
    for label in cb.ax.yaxis.get_ticklabels():
        label.set_size(18)
    for i in np.arange(0, len(out_object.y), 5):
        for j in np.arange(0, len(out_object.x), 5):
            plt.quiver(out_object.x2d[i][j], out_object.y2d[i][j], out_object.u_disp[i][j], out_object.v_disp[i][j],
                       units='width', scale=0.2)
    rectangles, points, mogis = utilities.separate_source_types(out_object.source_object)
    fault_sources = rectangles + points
    for source in fault_sources:
        [x_total, y_total, x_updip, y_updip] = source.get_fault_four_corners()
        plt.plot(x_total, y_total, 'k', linewidth=1)
        plt.plot(x_updip, y_updip, 'g', linewidth=3)
        center = source.get_fault_center()
        plt.plot(center[0], center[1], '.g', markersize=8)
    for rec in out_object.receiver_object:
        [x_total, y_total, x_updip, y_updip] = rec.get_fault_four_corners()
        plt.plot(x_total, y_total, 'b', linewidth=1)
        plt.plot(x_updip, y_updip, 'b', linewidth=3)
        center = rec.get_fault_center()
        plt.plot(center[0], center[1], '.b', markersize=8)
    plt.xlim([out_object.x.min(), out_object.x.max()])
    plt.ylim([out_object.y.min(), out_object.y.max()])
    plt.gca().tick_params(axis='both', which='major', labelsize=18)
    plt.grid()
    plt.axis('equal')
    plt.xlabel('X (km)', fontsize=20)
    plt.ylabel('Y (km)', fontsize=20)
    plt.gca().tick_params(labelsize=16)
    plt.title('Surface Displacement', fontsize=28)
    plt.savefig(outfile, facecolor="w")
    plt.close()
    return


def stress_plot(params, out_object, stress_type, vmin=None, vmax=None):
    """
    default vmin,vmax are in KPa
    Plots of fault patches in cartesian space, colored by the magnitude of the stress component.
    """

    if not out_object.receiver_object:
        return

    outfile = os.path.join(params.outdir, 'Stresses_' + stress_type + '.png')
    print("Making plot of %s stress on receiver fault patches: %s. " % (stress_type, outfile))

    if stress_type == 'Shear':
        stress_component = out_object.receiver_shear
    elif stress_type == 'Normal':
        stress_component = out_object.receiver_normal
    else:
        stress_component = out_object.receiver_coulomb

    # Select boundaries of color map.
    vmin, vmax = utilities.produce_vmin_vmax_symmetric(stress_component, vmin, vmax)
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='RdYlBu_r')

    # Figure of stresses.
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    fig = plt.figure(figsize=(12, 10), dpi=300)

    for i in range(len(stress_component)):
        [x_total, y_total, _, _] = out_object.receiver_object[i].get_fault_four_corners()
        xcoords = x_total[0:4]
        ycoords = y_total[0:4]
        fault_vertices = np.column_stack((xcoords, ycoords))
        patch_color = custom_cmap.to_rgba(stress_component[i])

        mypolygon = Polygon(fault_vertices, color=patch_color, alpha=1.0)
        plt.gca().add_patch(mypolygon)

    custom_cmap.set_array(np.arange(vmin, vmax, 100))
    ax1 = plt.gca()
    cb = fig.colorbar(custom_cmap, ax=ax1, location='right')
    cb.set_label(stress_type + ' stress (kPa)', fontsize=22)
    for label in cb.ax.yaxis.get_ticklabels():
        label.set_size(18)

    # Adding source and receiver faults
    for source in out_object.source_object:
        [x_total, y_total, x_updip, y_updip] = source.get_fault_four_corners()
        plt.plot(x_total, y_total, 'k', linewidth=1)
        plt.plot(x_updip, y_updip, 'g', linewidth=3)
        center = source.get_fault_center()
        plt.plot(center[0], center[1], '.g', markersize=8)
    for rec in out_object.receiver_object:
        [x_total, y_total, x_updip, y_updip] = rec.get_fault_four_corners()
        plt.plot(x_total, y_total, 'b', linewidth=1)
        plt.plot(x_updip, y_updip, 'b', linewidth=3)
        center = rec.get_fault_center()
        plt.plot(center[0], center[1], '.b', markersize=8)

    plt.grid()
    plt.axis('equal')
    plt.gca().tick_params(axis='both', which='major', labelsize=18)
    plt.title(stress_type + ' stress from source faults', fontsize=22)
    plt.xlim([out_object.x.min(), out_object.x.max()])
    plt.ylim([out_object.y.min(), out_object.y.max()])
    plt.xlabel('X (km)', fontsize=20)
    plt.ylabel('Y (km)', fontsize=20)
    plt.gca().tick_params(labelsize=16)
    plt.savefig(outfile, facecolor="w")
    plt.close()
    return


def stress_cross_section_cartesian(params, out_object, stress_type, vmin=None, vmax=None, writefile=None):
    """
    Rotate existing receivers into a depth cross-section and plot their stress values.
    vmin,vmax are in kPa
    Vertical plots of fault patches in cartesian space, colored by the magnitude of the stress component.
    """

    if not out_object.receiver_object:
        return
    outfile = os.path.join(params.outdir, 'Stresses_cross_section_' + stress_type + '.png')

    print("Making plot of %s stress on receiver fault patches: %s. " % (stress_type, outfile))

    if stress_type == 'Shear':
        stress_component = out_object.receiver_shear
    elif stress_type == 'Normal':
        stress_component = out_object.receiver_normal
    else:
        stress_component = out_object.receiver_coulomb

    max_depth_array = [x.bottom for x in out_object.receiver_object]
    min_depth_array = [x.top for x in out_object.receiver_object]

    # Select boundaries of color map, usually forcing even distribution around 0
    vmin, vmax = utilities.produce_vmin_vmax_symmetric(stress_component, vmin, vmax)
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='RdYlBu_r')

    # Figure of stresses.
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    fig = plt.figure(figsize=(17, 8), dpi=300)
    if writefile:
        print("Writing file %s" % writefile)
        ofile = open(writefile, 'w')
        ofile.close()

    total_y_array = []
    for i in range(len(stress_component)):
        xcoords = [out_object.receiver_object[i].xstart, out_object.receiver_object[i].xstart,
                   out_object.receiver_object[i].xfinish, out_object.receiver_object[i].xfinish]
        ycoords = [out_object.receiver_object[i].ystart, out_object.receiver_object[i].ystart,
                   out_object.receiver_object[i].yfinish, out_object.receiver_object[i].yfinish]
        zcoords = [out_object.receiver_object[i].top, out_object.receiver_object[i].bottom,
                   out_object.receiver_object[i].bottom, out_object.receiver_object[i].top]
        fault_strike = out_object.receiver_object[i].strike
        x1, y1 = conversion_math.rotate_points(xcoords[0], ycoords[0], fault_strike)
        x2, y2 = conversion_math.rotate_points(xcoords[1], ycoords[1], fault_strike)
        x3, y3 = conversion_math.rotate_points(xcoords[2], ycoords[2], fault_strike)
        x4, y4 = conversion_math.rotate_points(xcoords[3], ycoords[3], fault_strike)
        ycoords = [y1, y2, y3, y4]
        total_y_array = total_y_array + ycoords
        fault_vertices = np.column_stack((ycoords, zcoords))
        patch_color = custom_cmap.to_rgba(stress_component[i])

        mypolygon = Polygon(fault_vertices, color=patch_color, alpha=1.0)
        plt.gca().add_patch(mypolygon)

        if writefile:   # the x-axis here is an arbitrary position along the orientation of the receiver fault
            ofile = open(writefile, 'a')
            ofile.write("> -Z%f\n" % stress_component[i])
            ofile.write("%f %f \n" % (-y1, zcoords[0]))
            ofile.write("%f %f \n" % (-y2, zcoords[1]))
            ofile.write("%f %f \n" % (-y3, zcoords[2]))
            ofile.write("%f %f \n" % (-y4, zcoords[3]))
            ofile.write("%f %f \n" % (-y1, zcoords[0]))
            ofile.close()

    custom_cmap.set_array(np.arange(vmin, vmax, 100))
    cb = fig.colorbar(custom_cmap, ax=plt.gca(), location='right')
    cb.set_label(stress_type + ' Stress (kPa)', fontsize=22)
    for label in cb.ax.yaxis.get_ticklabels():
        label.set_size(18)

    plt.grid()
    plt.gca().tick_params(axis='both', which='major', labelsize=18)
    plt.title(stress_type + ' stress from source faults', fontsize=22)
    plt.xlim([np.min(total_y_array)-1, np.max(total_y_array)+1])
    plt.xlabel('Distance along fault profile (km)', fontsize=18)
    plt.ylabel('Depth (km)', fontsize=18)
    plt.ylim([min(min_depth_array), max(max_depth_array)])
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()
    plt.savefig(outfile, facecolor="w")
    plt.close()
    return


def map_horiz_profile(horiz_profile, profile_results, outfile, vmin=-200, vmax=200, cap_colorbar=1):
    """Display a small map of a horizontal profile of stresses. Default colors for now."""

    x = np.reshape(horiz_profile.lon1d, horiz_profile.shape)
    y = np.reshape(horiz_profile.lat1d, horiz_profile.shape)
    _normal_stress = np.reshape(profile_results[0], horiz_profile.shape)
    _shear_stress = np.reshape(profile_results[1], horiz_profile.shape)
    coulomb_stress = np.reshape(profile_results[2], horiz_profile.shape)

    # Figure of stresses.
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    fig = plt.figure(figsize=(14, 8), dpi=300)
    levels = np.linspace(vmin, vmax, 20)
    if cap_colorbar:  # do we force the higher values to the top of the colorbar?
        coulomb_stress[coulomb_stress > vmax] = vmax
        coulomb_stress[coulomb_stress < vmin] = vmin
    dislay_map = plt.contourf(x, y, coulomb_stress, levels=levels, cmap='RdYlBu_r', vmin=vmin, vmax=vmax)
    plt.title('Coulomb stresses on horizontal profile, fixed strike/dip/rake/depth of '+str(horiz_profile.strike)+', ' +
              str(horiz_profile.dip)+', '+str(horiz_profile.rake)+', '+str(horiz_profile.depth_km))

    cb = fig.colorbar(dislay_map, ax=plt.gca(), location='right')
    cb.set_label('Stress (kPa)', fontsize=22)
    for label in cb.ax.yaxis.get_ticklabels():
        label.set_size(18)
    plt.xlabel('Longitude', fontsize=18)
    plt.ylabel('Latitude', fontsize=18)
    plt.gca().tick_params(labelsize=16)
    print("Saving figure %s " % outfile)
    plt.savefig(outfile, facecolor="w")
    return


def write_synthetic_grid_triplets(out_object, outdir, east_model_file, north_model_file, vert_model_file):
    """
    Write lists of lon/lat/def for each component of deformation in synthetic grid.
    Used for GMT plots.
    """
    print("Writing synthetic grid into files %s etc." % outdir+east_model_file)
    ofile_w = open(outdir + vert_model_file, 'w')
    ofile_u = open(outdir + east_model_file, 'w')
    ofile_v = open(outdir + north_model_file, 'w')
    for i in np.arange(0, len(out_object.y)):
        for j in np.arange(0, len(out_object.x)):
            loni, lati = fault_vector_functions.xy2lonlat(out_object.x2d[i][j], out_object.y2d[i][j],
                                                          out_object.zerolon, out_object.zerolat)
            ofile_w.write("%f %f %f\n" % (loni, lati, out_object.w_disp[i][j]))
            ofile_u.write("%f %f %f\n" % (loni, lati, out_object.u_disp[i][j]))
            ofile_v.write("%f %f %f\n" % (loni, lati, out_object.v_disp[i][j]))
    ofile_w.close()
    ofile_u.close()
    ofile_v.close()
    return


def write_synthetic_grid_full_results(out_object, outfile):
    # Write output of synthetic displacement grid in cartesian and lon/lat coordinates
    print("Writing synthetic grid of displacements in %s " % outfile)
    ofile = open(outfile, 'w')
    ofile.write("# Format: x y lon lat x_disp[m] y_disp[m] z_disp[m] \n")
    for i in np.arange(0, len(out_object.y)):
        for j in np.arange(0, len(out_object.x)):
            loni, lati = fault_vector_functions.xy2lonlat(out_object.x2d[i][j], out_object.y2d[i][j],
                                                          out_object.zerolon, out_object.zerolat)
            ofile.write("%f %f %f %f %f %f %f\n" % (out_object.x2d[i][j], out_object.y2d[i][j], loni, lati,
                                                    out_object.u_disp[i][j], out_object.v_disp[i][j],
                                                    out_object.w_disp[i][j]))
    ofile.close()
    return


def write_horiz_profile(horiz_profile, profile_results, outfile):
    print("Writing %s " % outfile)
    ofile = open(outfile, 'w')
    ofile.write("# lon lat depth_km normal_kPa shear_kPa coulomb_kPa\n")
    ofile.write("# strike %f, dip %f, rake %f\n" % (horiz_profile.strike, horiz_profile.dip, horiz_profile.rake))
    for i in range(len(horiz_profile.lon1d)):
        ofile.write("%f %f %f %f %f %f\n" % (horiz_profile.lon1d[i], horiz_profile.lat1d[i], horiz_profile.depth_km,
                                             profile_results[0][i], profile_results[1][i], profile_results[2][i]))
    return


def write_disp_grd_files(inputs, outdir, east_txt, north_txt, vert_txt):
    # Make surfaces of east/north/up deformation for plotting. Based on three text files in outdir, with x-y-disp
    region = inputs.define_map_region()
    inc = (region[1]-region[0]) / 500  # for gmt mapping, default is 500 points per axis
    utilities.displacements_to_3_grds(outdir, efiles=(east_txt, 'east.grd'), nfiles=(north_txt, 'north.grd'),
                                      ufiles=(vert_txt, 'vert.grd'), region=region, inc=inc)
    return
