# output_manager

import numpy as np
import pygmt, os
from . import io_additionals, utilities
from Tectonic_Utils.geodesy import fault_vector_functions


def map_stress_plot(params, inputs, out_object, stress_type, vmin=-1, vmax=1):
    """
    Fill in mapped fault patches with colors corresponding to their stress changes.
    """
    if not out_object.receiver_object:
        return

    if stress_type == 'Shear':
        plotting_stress = out_object.receiver_shear
    elif stress_type == 'Normal':
        plotting_stress = out_object.receiver_normal
    else:
        plotting_stress = out_object.receiver_coulomb  # The default option

    # Make stress bounds for color map.
    [cmap_opts, cbar_opts] = utilities.define_colorbar_series(plotting_stress, vmin=vmin, vmax=vmax)

    # Make cpt
    pygmt.makecpt(cmap="jet", series=str(cmap_opts[0]) + "/" + str(cmap_opts[1]) + "/"+str(cmap_opts[2]),
                  output="mycpt.cpt", background=True)

    # Make Map
    region = inputs.define_map_region()
    proj = "M7i"
    fig = pygmt.Figure()
    title = "+t\"" + stress_type + " stress\""  # must put escaped quotations around the title.
    fig.basemap(region=region, projection=proj, frame=title)
    fig.coast(shorelines="1.0p,black", region=region, borders="1", projection=proj, frame="1.0")  # the boundary.
    fig.coast(region=region, projection=proj, borders='2', shorelines='0.5p,black', water='white',
              map_scale="jBL+c1.5+w50")

    fig = annotate_figure_with_sources(fig, inputs, params.mu)

    # Plot the stress components
    utilities.write_fault_edges_to_gmt_file(out_object.receiver_object, 'tmp.txt', color_array=plotting_stress)
    fig.plot(data='tmp.txt', pen="thick,black", fill="+z", cmap="mycpt.cpt")
    os.remove('tmp.txt')

    # Colorbar annotation
    fig.coast(shorelines="1.0p,black", region=region, projection=proj)  # the boundary.
    fig.colorbar(position="jBr+w3.5i/0.2i+o2.5c/1.5c+h", cmap="mycpt.cpt", shading="0.8",
                 truncate=str(cbar_opts[0]) + "/" + str(cbar_opts[1]), frame=["x" + str(0.2), "y+L\"Stress (kPa)\""])

    # Annotate with aftershock locations
    fig = annotate_figure_with_aftershocks(fig, aftershocks_file=params.aftershocks, style='c0.02c')
    print("Mapping stress in %s " % os.path.join(params.outdir, stress_type + '_map.png'))
    fig.savefig(os.path.join(params.outdir, stress_type + '_map.png'))
    return


def map_vertical_def(params, inputs, outfile, max_down=-0.045, max_up=0.045):
    """
    Simple map of grdfile with subsampled vertical deformation.
    Currently a proof of concept!
    """
    print("Mapping vertical deformation in %s " % outfile)

    proj = 'M4i'
    region = inputs.define_map_region()

    # Build a PyGMT plot
    fig = pygmt.Figure()
    pygmt.makecpt(cmap="roma", series=str(max_down) + "/"+str(max_up)+"/0.001", background="o", output="mycpt.cpt")
    fig.basemap(region=region, projection=proj, frame="+t\"Modeled Vertical Displacement\"")
    fig.grdimage(os.path.join(params.outdir, 'vert.grd'), region=region, cmap="mycpt.cpt")
    fig.coast(region=region, projection=proj, borders='1', shorelines='1.0p,black', water='lightblue',
              map_scale="n0.23/0.06+c" + str(region[2]) + "+w20", frame="1.0")

    fig = annotate_figure_with_sources(fig, inputs, params.mu, dotstyle='s0.01c')
    fig = annotate_figure_with_aftershocks(fig, aftershocks_file=params.aftershocks, style='c0.02c')

    fig.colorbar(position="JCR+w4.0i+v+o0.7i/0i", cmap="mycpt.cpt", truncate=str(max_down)+"/"+str(max_up),
                 frame=["x0.01", "y+L\"Disp(m)\""])
    fig.savefig(outfile)
    return


def map_displacement_vectors(params, inputs, obs_disp_points, model_disp_points, outfile, vmin=None, vmax=None):
    """
    Make a plot of modeled vector displacement from model_disp_points.
    obs_disp_points is a list of disp_points objects that can be used to also plot observations at the same points.
    """
    if len(model_disp_points) == 0:
        return

    proj = 'M4i'
    region = inputs.define_map_region()

    # Unpack modeled displacements, in meters
    model_dE = np.array([x.dE_obs for x in model_disp_points])
    model_dN = np.array([x.dN_obs for x in model_disp_points])
    model_dU = np.array([x.dU_obs for x in model_disp_points])
    model_lon = np.array([x.lon for x in model_disp_points])
    model_lat = np.array([x.lat for x in model_disp_points])

    # Make modeled vertical displacement color map
    [cmap_opts, cbar_opts] = utilities.define_colorbar_series(model_dU, vmin, vmax)
    pygmt.makecpt(cmap="roma", series=str(cmap_opts[0])+"/"+str(cmap_opts[1])+"/"+str(cmap_opts[2]), background="o",
                  output="mycpt.cpt")

    # Automatically define the vector scale bar based on max modeled displacements
    scale_arrow, vectext = utilities.define_vector_scale_size(model_dE, model_dN)
    scale = 2.5 / scale_arrow  # empirical scaling, convenient display

    # Build a PyGMT plot
    fig = pygmt.Figure()
    fig.basemap(region=region, projection=proj, frame="+t\"Coseismic Displacements\"")
    fig.coast(region=region, projection=proj, borders='1', shorelines='1.0p,black', water='lightblue',
              map_scale="n0.4/0.06+c" + str(region[2]) + "+w20", frame="1.0")
    fig.plot(x=model_lon, y=model_lat, style='c0.3c', fill=model_dU, cmap='mycpt.cpt', pen="thin,black")

    fig.plot(x=model_lon, y=model_lat, style='v0.2c+e+gblack+h0+p1p,black+z'+str(scale),
             direction=[model_dE, model_dN], pen="thin,black")
    fig.plot(x=[region[0]+0.30], y=[region[2]+0.05],  style='v0.2c+e+gblack+h0+p1p,black+z'+str(scale),
             direction=[[scale_arrow], [0]],  pen="thin,black")  # scale vector
    fig.text(x=[region[0]+0.45], y=[region[2]+0.15], text=vectext+" model")  # scale label
    # Plot the observations if they exist
    obs_dE = np.array([x.dE_obs for x in obs_disp_points])
    obs_dN = np.array([x.dN_obs for x in obs_disp_points])
    if sum(~np.isnan(obs_dE)) > 0:
        fig.plot(x=model_lon, y=model_lat, style='v0.2c+e+gred+h0+p1p,red+z' + str(scale),
                 direction=[obs_dE, obs_dN], pen="thin,red")
        fig.plot(x=[region[0]+0.30], y=[region[2]+0.35],  style='v0.2c+e+gred+h0+p1p,red+z'+str(scale),
                 direction=[[scale_arrow], [0]],  pen="thin,red")  # scale vector
        fig.text(x=[region[0]+0.50], y=[region[2]+0.45], text=vectext+' obs')  # scale label

    fig = annotate_figure_with_sources(fig, inputs, params.mu)
    fig = annotate_figure_with_aftershocks(fig, aftershocks_file=params.aftershocks)

    fig.colorbar(position="JCR+w4.0i+v+o0.7i/0i", cmap="mycpt.cpt", truncate=str(cbar_opts[0])+"/"+str(cbar_opts[1]),
                 frame=["x"+str(cbar_opts[2]), "y+L\"Vert Disp(m)\""])
    fig.savefig(outfile)
    return


def annotate_figure_with_sources(fig, inputs, mu=30e9, fmscale="0.3c", dotstyle="s0.3c"):
    """
    Draw each source in inputs.source_object, a list of sources.
    Inputs and mu provide additional information about the calculation.

    :param fig: pygmt figure
    :param inputs: coulomb_collection.inputs_object, which has a list of sources in inputs.source_object
    :param mu: shear modulus, float
    :param fmscale: focal mechanism size
    :param dotstyle: psxy style for source dots
    """
    if len(inputs.source_object) == 0:
        return fig

    # Draw dots for EQ sources
    eq_lon, eq_lat = [], []
    for source in inputs.source_object:  # Plotting all types of sources
        source_lon, source_lat = fault_vector_functions.xy2lonlat(source.xstart, source.ystart, inputs.zerolon,
                                                                  inputs.zerolat)
        eq_lon.append(source_lon)
        eq_lat.append(source_lat)
    fig.plot(x=eq_lon, y=eq_lat, style=dotstyle, fill="purple", pen="thin,black")

    rect_sources, _ = utilities.separate_source_types(inputs.source_object)
    for source in rect_sources:  # draw focal mechanisms
        # Draw small sources with focal mechanisms
        if source.area < 0.001 * 0.001:  # If the area is < 1 m^2:
            [x_total, y_total, _, _] = source.get_fault_four_corners()
            lons, lats = fault_vector_functions.xy2lonlat(x_total, y_total, inputs.zerolon, inputs.zerolat)
            _, mag = source.get_fault_slip_moment(mu)
            focal_mechanism = dict(strike=source.strike, dip=source.dipangle, rake=source.rake, magnitude=mag)
            fig.meca(focal_mechanism, scale=fmscale, longitude=lons[0], latitude=lats[0], depth=source.top)
    # draw the fault patches, no special color code
    utilities.write_fault_edges_to_gmt_file(rect_sources, "tmp.txt")
    fig.plot(data='tmp.txt', pen="0.2p,black")
    os.remove('tmp.txt')
    return fig


def annotate_figure_with_aftershocks(fig, aftershocks_file=None, style='c0.02c', pen="0.1p,black", color="black"):
    """
    Annotate a figure with aftershock locations.
    """
    if aftershocks_file:
        [lon, lat, _, _, _] = io_additionals.read_aftershock_table(aftershocks_file)
        fig.plot(x=lon, y=lat, style=style, fill=color, pen=pen)
    return fig
