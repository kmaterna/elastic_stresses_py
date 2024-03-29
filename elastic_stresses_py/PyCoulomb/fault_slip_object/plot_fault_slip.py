"""Functions to take a slip distribution and make a map. """

import numpy as np
import pygmt
import os
from . import fault_slip_object as fso
from . import file_io
from .. import utilities
from ..disp_points_object import utilities as dpo_utils


def unpack_disp_points(disp_points):
    """
    Unpack all displacement points.
    Also unpack any displacement points that have vertical data
    """
    lon = np.array([x.lon for x in disp_points])
    lat = np.array([x.lat for x in disp_points])
    disp_x = np.array([x.dE_obs for x in disp_points])
    disp_y = np.array([x.dN_obs for x in disp_points])
    disp_z = np.array([x.dU_obs for x in disp_points])
    lon_vert = np.array([x.lon for x in disp_points if ~np.isnan(x.dU_obs)])
    lat_vert = np.array([x.lat for x in disp_points if ~np.isnan(x.dU_obs)])
    disp_z_vert = np.array([x.dU_obs for x in disp_points if ~np.isnan(x.dU_obs)])
    return [lon, lat, disp_x, disp_y, disp_z, lon_vert, lat_vert, disp_z_vert]


def unpack_horiz_disp_points_for_vectors(disp_points):
    """Unpack any displacement points that have horizontal data will be plotted as vectors"""
    lon_horiz = np.array([x.lon for x in disp_points if ~np.isnan(x.dE_obs)])
    lat_horiz = np.array([x.lat for x in disp_points if ~np.isnan(x.dE_obs)])
    disp_x_horiz = np.array([x.dE_obs for x in disp_points if ~np.isnan(x.dE_obs)])
    disp_y_horiz = np.array([x.dN_obs for x in disp_points if ~np.isnan(x.dE_obs)])
    return [lon_horiz, lat_horiz, disp_x_horiz, disp_y_horiz]


def write_patch_edges_for_plotting(fault_dict_list, plotting_function=lambda x: file_io.outputs.get_total_slip(x)):
    """
    Plot the full patches and the surface traces of a collection of rectangular faults.  Colorcode by slip (default).
    plotting_function is a function that operates on a fault_dict object to return a float.
    """
    file_io.outputs.write_gmt_fault_file(fault_dict_list, 'tmp.txt', color_mappable=plotting_function, verbose=False)
    file_io.outputs.write_gmt_surface_trace(fault_dict_list, 'tmp2.txt', verbose=False)
    return "tmp.txt", "tmp2.txt"


def automatically_determine_map_region(region_given=None, fault_dict_list=None, disp_points=None, buffer_deg=0.15):
    """
    Automatically select a map region [W, E, S, N] for a given list of fault segments and disp_point_objects.
    """
    region = region_given
    if not region_given:  # automatically determine region
        if len(fault_dict_list) > 0:
            minlon, maxlon, minlat, maxlat = fso.get_four_corners_lon_lat_multiple(fault_dict_list)
        elif len(disp_points) > 0:
            minlon, maxlon, minlat, maxlat = dpo_utils.extract_region_from_disp_points(disp_points)
        else:
            raise ValueError("Error! Cannot automatically determine map region no fault patches or disp_points.")
        region = [minlon-buffer_deg, maxlon+buffer_deg, minlat-buffer_deg, maxlat+buffer_deg]
    return region


def map_source_slip_distribution(fault_dict_list, outfile, disp_points=(), region=None,
                                 scale_arrow=(1.0, 0.010, "10 mm"),
                                 fault_traces_from_memory=None, fault_traces_from_dict=None,
                                 fault_traces_from_file=None, title="", vmin=None, vmax=None, v_labeling_interval=None,
                                 plot_slip_colorbar=True, vert_disp_units="m", vert_mult=1, map_scale=25,
                                 slip_cbar_opts=(-1, 1, 0.001)):
    """
    Plot a map of slip distribution from fault_dict_list, a general format for slip distributions of rectangular faults.
    In order to use this function with other fault formats, convert to the internal rectangular fault dict first.

    :param fault_dict_list: list of fault_dict objects
    :param outfile: string, name of file
    :param disp_points: list of disp_point objects
    :param region: tuple of 4 numbers, (W, E, S, N)
    :param scale_arrow: tuple of 3 numbers
    :param fault_traces_from_memory: list of tuples with ([lons], [lats]) for plotting multiple fault traces
    :param fault_traces_from_dict: a list of fault_dict objects that will be used for updip fault traces
    :param fault_traces_from_file: string, optional filename with fault traces to be plotted
    :param title: string
    :param vmin: float, bottom of vertical color bar for disp_points (optional)
    :param vmax: float, top of vertical color bar for disp_points (optional)
    :param v_labeling_interval: float for labeling vertical displacement color bar (optional)
    :param plot_slip_colorbar: bool, whether to show a color bar for fault slip
    :param vert_disp_units: string, describing the units of the vertical scale bar
    :param vert_mult: can turn verticals into mm by providing 1000 if you want (default is meters)
    :param map_scale: int, in km
    :param slip_cbar_opts: tuple with (cbar_min, cbar_max, cbar_int) for fault slip cbar
    """
    print("Plotting outfile %s " % outfile)
    proj = "M7i"
    if not region:  # automatically determine the region
        region = automatically_determine_map_region(region, fault_dict_list, disp_points, buffer_deg=0.15)
    fig = pygmt.Figure()
    fig_width_deg = region[1] - region[0]
    fig_height_deg = region[3] - region[2]
    fig.basemap(region=region, projection=proj, frame="+t\"" + title + "\"")
    fig.coast(shorelines="1.0p,black", region=region, borders="1", projection=proj, frame=str(fig_width_deg/5))
    fig.coast(region=region, projection=proj, borders='2', shorelines='0.5p,black', water='lightblue')

    # Draw fault slip for a list of faults and draw a color scale.
    if len(fault_dict_list) > 0:
        fault_colors = [_i.slip for _i in fault_dict_list]
        if slip_cbar_opts:
            cmap_opts, cbar_opts = slip_cbar_opts, slip_cbar_opts  # user-defined color bar
        else:
            [cmap_opts, cbar_opts] = utilities.define_colorbar_series(fault_colors)  # auto-defined color bar
        pygmt.makecpt(cmap="devon", series=str(cmap_opts[0])+"/"+str(cmap_opts[1])+"/"+str(cmap_opts[2]),
                      truncate="0/1", background="o", reverse=True, output="mycpt.cpt")
        file1, file2 = write_patch_edges_for_plotting(fault_dict_list)  # write source patches, colored by slip
        fig.plot(data=file1, pen="0.2p,black", fill="+z", cmap="mycpt.cpt")
        fig.plot(data=file2, pen="0.6p,black")   # shallow edges
        os.remove(file1)
        os.remove(file2)
        if plot_slip_colorbar:
            fig.colorbar(position="jBr+w3.5i/0.2i+o2.5c/1.5c+h", cmap="mycpt.cpt",
                         truncate=str(cbar_opts[0]) + "/" + str(cbar_opts[1]),
                         frame=["x" + str(cbar_opts[2]), "y+L\"Slip(m)\""])
        os.remove("mycpt.cpt")

    # Draw pre-written fault edges from GMT format, with colorbar (useful for triangles or descriptive data)
    if fault_traces_from_file:
        if slip_cbar_opts is None:
            raise ValueError("Error! Requesting to plot slip from GMT file, but no slip_cbar_opts provided.")
        pygmt.makecpt(cmap="devon", truncate="0/1", background="o", reverse=True, output="mycpt.cpt",
                      series=str(slip_cbar_opts[0]) + "/" + str(slip_cbar_opts[1]) + "/" + str(slip_cbar_opts[2]))
        fig.plot(data=fault_traces_from_file, pen="0.2p,black", fill="+z", cmap='mycpt.cpt')
        fig.colorbar(position="jBr+w3.5i/0.2i+o2.5c/1.5c+h", cmap="mycpt.cpt",
                     truncate=str(slip_cbar_opts[0]) + "/" + str(slip_cbar_opts[1]),
                     frame=["x" + str(slip_cbar_opts[2]), "y+L\"Slip(m)\""])
        os.remove("mycpt.cpt")

    # Draw fault traces (lines) on the plot
    if fault_traces_from_memory:
        for item in fault_traces_from_memory:
            fig.plot(x=item[0], y=item[1], pen="thickest,darkred")

    # Draw the updip trace of each fault segment
    if fault_traces_from_dict:
        for item in fault_dict_list:
            lons, lats = item.get_updip_corners_lon_lat()
            fig.plot(x=lons, y=lats, pen="thickest,darkred")

    # Draw vector displacements
    if len(disp_points) > 0:
        [lon, lat, _, _, disp_z, lon_vert, lat_vert, disp_z_vert] = unpack_disp_points(disp_points)
        [lon_horiz, lat_horiz, disp_x_horiz, disp_y_horiz] = unpack_horiz_disp_points_for_vectors(disp_points)
        fig.plot(x=lon, y=lat, style='s0.07i', fill='blue', pen="thin,black")
        if sum(~np.isnan(disp_z)) > 0:  # display vertical data if it's provided
            disp_z_vert = np.multiply(disp_z_vert, vert_mult)
            [v_cmap_opts, v_cbar_opts] = utilities.define_colorbar_series(disp_z_vert, tol=0.0001,
                                                                          v_labeling_interval=v_labeling_interval,
                                                                          vmin=vmin, vmax=vmax)
            series_str = str(v_cmap_opts[0])+"/" + str(v_cmap_opts[1]) + "/" + str(v_cmap_opts[2])
            pygmt.makecpt(cmap="roma", series=series_str, background="o", output="vert.cpt")
            fig.plot(x=lon_vert, y=lat_vert, style='c0.3c', fill=disp_z_vert, cmap='vert.cpt', pen="thin,black")
            truncate_str = str(v_cbar_opts[0]) + "/" + str(v_cbar_opts[1]),
            fig.colorbar(position="JCR+w4.0i+v+o0.7i/0i", cmap="vert.cpt", truncate=truncate_str,
                         frame=["x"+str(v_cbar_opts[2]), "y+L\"Vert Disp("+vert_disp_units+")\""])
            os.remove('vert.cpt')
        if sum(~np.isnan(disp_x_horiz) > 0):
            scale = scale_arrow[0] * (1/scale_arrow[1])  # empirical scaling for convenient display
            fig.plot(x=lon_horiz, y=lat_horiz, style='v0.2c+e+gblack+h0+p1p,black+z'+str(scale),
                     direction=[disp_x_horiz, disp_y_horiz], pen="thin,black")
            fig.plot(x=[region[0]+0.1*fig_width_deg], y=[region[2]+0.1*fig_height_deg],
                     style='v0.2c+e+gblack+h0+p1p,black+z'+str(scale),
                     direction=[[scale_arrow[1]], [0]], pen="thin,black")  # scale vector
            fig.text(x=[region[0]+0.15*fig_width_deg], y=[region[2]+0.15*fig_height_deg], text=scale_arrow[2],
                     font='14p,Helvetica,black')  # scale label

    # Map km scale
    fig.coast(region=region, projection=proj, borders='2', shorelines='0.5p,black', map_scale="jBL+o0.7c/1c+w" +
                                                                                              str(map_scale))
    fig.savefig(outfile)
    return


def plot_data_model_residual(outfile, disp_points, model_disp_points, resid_disp_points, region,
                             scale_arrow, fault_dict_list=(), v_labeling_interval=None, rms=None):
    """
    Plot data, model, residual vectors
    """
    print("Plotting outfile %s " % outfile)
    fig = pygmt.Figure()
    proj = 'M3.2i'
    point_size = 0.15  # for vertical symbol
    numrows, numcols = 1, 3
    with fig.subplot(nrows=numrows, ncols=numcols, figsize=("12i", "12i"), frame="lrtb"):
        # Data
        with fig.set_panel(panel=0):
            fig.basemap(region=region, projection=proj, frame=["WeSn", "2.0"])
            fig.coast(region=region, projection=proj, borders='2', shorelines='1.0p,black', water='lightblue')
            [_, _, _, _, disp_z, lon_vert, lat_vert, disp_z_vert] = unpack_disp_points(disp_points)
            [lon_horiz, lat_horiz, dispx_horiz, dispy_horiz] = unpack_horiz_disp_points_for_vectors(disp_points)

            if sum(~np.isnan(disp_z)) > 0:  # display vertical data if it's provided
                [cmap_opts, cbar_opts] = utilities.define_colorbar_series(disp_z_vert,
                                                                          v_labeling_interval=v_labeling_interval)
                series_str = str(cmap_opts[0]) + "/" + str(cmap_opts[1]) + "/" + str(cmap_opts[2])
                pygmt.makecpt(cmap="roma", series=series_str, background="o", output="vert.cpt")
                fig.plot(x=lon_vert, y=lat_vert, projection=proj, style='c'+str(point_size)+'c', fill=disp_z_vert,
                         cmap='vert.cpt', pen="0.1p,black")

            scale = scale_arrow[0] * (1/scale_arrow[1])  # empirical scaling for convenient display
            fig.plot(x=lon_horiz, y=lat_horiz, projection=proj, style='v0.12c+e+gblack+h0+p0.1p,black+z'+str(scale),
                     direction=[dispx_horiz, dispy_horiz], pen="0.1p,black")
            fig.plot(x=[region[0]+0.30], y=[region[2]+1.05], projection=proj,
                     style='v0.2c+e+gblack+h0+p1p,black+z'+str(scale),
                     direction=[[scale_arrow[1]], [0]],  pen="thin,black")  # scale vector
            fig.text(x=[region[0]+0.85], y=[region[2]+1.25], projection=proj, text=scale_arrow[2]+' data')  # label
            fig.text(position="TL", projection=proj, font="18p,Helvetica,black", fill="white", offset='j0.1/0.1',
                     text='A) Data')

        # Model
        with fig.set_panel(panel=1):
            fig.basemap(region=region, projection=proj, frame=["WeSn", "2.0"])
            fig.coast(region=region, projection=proj, borders='2', shorelines='1.0p,black', water='lightblue')
            [_, _, _, _, disp_z, lon_vert, lat_vert, disp_z_vert] = unpack_disp_points(model_disp_points)
            [lon_horiz, lat_horiz, dispx_horiz, dispy_horiz] = unpack_horiz_disp_points_for_vectors(model_disp_points)

            file1, file2 = write_patch_edges_for_plotting(fault_dict_list,
                                                          plotting_function=lambda x: x.get_blank_fault)
            fig.plot(data=file1, pen="0.2p,black", fill="white")  # annotate the rectangular fault patches
            fig.plot(data=file2, pen="0.6p,black", projection=proj)  # annotate shallow edges of rect. fault patches
            os.remove(file1)
            os.remove(file2)

            if sum(~np.isnan(disp_z)) > 0:  # display vertical data if it's provided
                fig.plot(x=lon_vert, y=lat_vert, projection=proj, style='c'+str(point_size)+'c', fill=disp_z_vert,
                         cmap='vert.cpt', pen="0.1p,black")

            scale = scale_arrow[0] * (1/scale_arrow[1])  # empirical scaling for convenient display
            fig.plot(x=lon_horiz, y=lat_horiz, projection=proj, style='v0.12c+e+gblack+h0+p0.1p,black+z'+str(scale),
                     direction=[dispx_horiz, dispy_horiz], pen="0.1p,black")
            fig.plot(x=[region[0]+0.30], y=[region[2]+1.05], projection=proj,
                     style='v0.2c+e+gblack+h0+p1p,black+z'+str(scale),
                     direction=[[scale_arrow[1]], [0]],  pen="thin,black")  # scale vector
            fig.text(x=[region[0]+0.95], y=[region[2]+1.25], projection=proj, text=scale_arrow[2]+' model')  # label
            fig.text(position="TL", projection=proj, font="18p,Helvetica,black", fill="white", offset='j0.1/0.1',
                     text='B) Model')

        # Residual
        with fig.set_panel(panel=2):
            fig.basemap(region=region, projection=proj, frame=["WeSn", "2.0"])
            fig.coast(region=region, projection=proj, borders='2', shorelines='1.0p,black', water='lightblue')
            [_, _, _, _, disp_z, lon_vert, lat_vert, disp_z_vert] = unpack_disp_points(resid_disp_points)
            [lon_horiz, lat_horiz, dispx_horiz, dispy_horiz] = unpack_horiz_disp_points_for_vectors(resid_disp_points)

            if sum(~np.isnan(disp_z)) > 0:  # display vertical data if it's provided
                fig.plot(x=lon_vert, y=lat_vert, projection=proj, style='c'+str(point_size)+'c', fill=disp_z_vert,
                         cmap='vert.cpt', pen="0.1p,black")
                fig.colorbar(position="JCR+w4.0i+v+o0.4i/0i", projection=proj, cmap="vert.cpt",
                             truncate=str(cbar_opts[0]) + "/" + str(cbar_opts[1]),
                             frame=["x" + str(cbar_opts[2]), "y+L\"Vert Disp(m)\""])

            scale = scale_arrow[0] * (1/scale_arrow[1])  # empirical scaling for convenient display
            fig.plot(x=lon_horiz, y=lat_horiz, projection=proj, style='v0.12c+e+gblack+h0+p0.1p,black+z'+str(scale),
                     direction=[dispx_horiz, dispy_horiz], pen="0.1p,black")
            fig.plot(x=[region[0]+0.30], y=[region[2]+1.05], projection=proj,
                     style='v0.2c+e+gblack+h0+p1p,black+z'+str(scale),
                     direction=[[scale_arrow[1]], [0]],  pen="thin,black")  # scale vector
            fig.text(x=[region[0]+0.85], y=[region[2]+1.25], projection=proj, text=scale_arrow[2]+' res')  # label
            fig.text(position="TL", projection=proj, font="18p,Helvetica,black", fill="white", offset='j0.1/0.1',
                     text='C) Residual')
            if rms:
                fig.text(x=[region[0] + 1.20], y=[region[2] + 0.65], projection=proj,
                         text='rms='+str(np.round(rms, 2)) + ' mm')  # label

    fig.savefig(outfile)
    return
