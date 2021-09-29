# output_manager

import numpy as np
import pygmt
from subprocess import call
from . import conversion_math
from . import io_additionals, utilities
from Tectonic_Utils.geodesy import fault_vector_functions


def map_stress_plot(params, inputs, out_object, stress_component):
    """
    Using PyGMT
    Filling in fault patches with colors corresponding to their stress changes
    """
    if stress_component == 'shear':
        plotting_stress = out_object.receiver_shear;
        label = 'Shear';
    elif stress_component == 'normal':
        plotting_stress = out_object.receiver_normal;
        label = 'Normal';
    else:
        plotting_stress = out_object.receiver_coulomb;  # The default option
        label = 'Coulomb';

    if not out_object.receiver_object:
        return;

    # Make stress bounds for color map.
    vmin, vmax = -1, 1;
    [cmap_opts, cbar_opts] = utilities.define_colorbar_series(plotting_stress, vmin=vmin, vmax=vmax);

    # Make cpt
    pygmt.makecpt(cmap="jet", series=str(cmap_opts[0]) + "/" + str(cmap_opts[1]) + "/"+str(cmap_opts[2]),
                  output="mycpt.cpt", background=True);

    # Make Map
    region = [inputs.minlon, inputs.maxlon, inputs.minlat, inputs.maxlat];
    proj = "M7i"
    fig = pygmt.Figure()
    title = "+t\"" + stress_component + " stress\"";  # must put escaped quotations around the title.
    fig.basemap(region=region, projection=proj, frame=title);
    fig.coast(shorelines="1.0p,black", region=region, borders="1", projection=proj, frame="1.0");  # the boundary.
    fig.coast(region=region, projection=proj, borders='2', shorelines='0.5p,black', water='white',
              map_scale="g-125.5/39.6+c1.5+w50");

    # Draw each source
    eq_lon, eq_lat = [], [];
    for source in out_object.source_object:
        source_lon, source_lat = fault_vector_functions.xy2lonlat(source.xstart, source.ystart, inputs.zerolon,
                                                                  inputs.zerolat);
        eq_lon.append(source_lon);
        eq_lat.append(source_lat);
        [x_total, y_total, _, _] = conversion_math.get_fault_four_corners(source);
        lons, lats = fault_vector_functions.xy2lonlat(x_total, y_total, inputs.zerolon, inputs.zerolat);
        if not source.potency:
            fig.plot(x=lons, y=lats, pen="thick,black");  # in case of area sources, outline them.
        else:
            fig.plot(x=lons, y=lats, style='s0.3c', color="purple", pen="thin,black");  # in case of point sources
    # Annotate with earthquake location.
    fig.plot(eq_lon, eq_lat, style='s0.3c', color="purple", pen="thin,black");

    # Draw each receiver, with associated data
    for i in range(len(out_object.receiver_object)):
        rec = out_object.receiver_object[i];
        [x_total, y_total, _, _] = conversion_math.get_fault_four_corners(rec);
        lons, lats = fault_vector_functions.xy2lonlat(x_total, y_total, inputs.zerolon, inputs.zerolat);
        fig.plot(x=lons, y=lats, zvalue=str(plotting_stress[i]), pen="thick,black", color="+z", cmap="mycpt.cpt");
        # color = stress

    # Colorbar annotation
    fig.coast(shorelines="1.0p,black", region=region, projection=proj);  # the boundary.
    fig.colorbar(position="jBr+w3.5i/0.2i+o2.5c/1.5c+h", cmap="mycpt.cpt", shading="0.8",
                 truncate=str(cbar_opts[0]) + "/" + str(cbar_opts[1]), frame=["x" + str(0.2), "y+L\"KPa\""]);

    # Annotate with aftershock locations
    if params.aftershocks:
        [lon, lat, _, _, _] = io_additionals.read_aftershock_table(params.aftershocks);
        fig.plot(lon, lat, style='c0.1c', color='black', pen="thin,black");

    fig.savefig(params.outdir + label + '_map.png');
    return;


def map_vertical_def(params, inputs, filename, outfile):
    """Simple map of grdfile with subsampled vertical deformation.
    Currently mess, but a proof of concept!
    Takes a grd file created by gmt surface from the xyz file written in this software. """
    print("Mapping file %s " % filename);

    proj = 'M4i'
    region = [inputs.minlon, inputs.maxlon, inputs.minlat, inputs.maxlat];

    # First make surfaces of east/north/up deformation for later plotting
    outdir = params.outdir;
    inc = 0.0005;
    call(['gmt', 'surface', outdir + '/xyz_model.txt', '-G' + outdir + '/vert.grd',
          '-R' + str(region[0]) + '/' + str(region[1]) + '/' + str(region[2]) + '/' + str(region[3]), '-I'+str(inc),
          '-r'], shell=False);
    call(['gmt', 'surface', outdir + '/xyu_model.txt', '-G' + outdir + '/east.grd',
          '-R' + str(region[0]) + '/' + str(region[1]) + '/' + str(region[2]) + '/' + str(region[3]), '-I'+str(inc),
          '-r'], shell=False);
    call(['gmt', 'surface', outdir + '/xyv_model.txt', '-G' + outdir + '/north.grd',
          '-R' + str(region[0]) + '/' + str(region[1]) + '/' + str(region[2]) + '/' + str(region[3]), '-I'+str(inc),
         '-r'], shell=False);

    # Build a PyGMT plot
    fig = pygmt.Figure();
    pygmt.makecpt(cmap="roma", series="-0.045/0.045/0.001", background="o", output="mycpt.cpt");
    fig.basemap(region=region, projection=proj, frame="+t\"Vertical Displacement\"");
    fig.grdimage(filename, region=region, cmap="mycpt.cpt");
    fig.coast(region=region, projection=proj, borders='1', shoreliens='1.0p,black', water='lightblue',
              map_scale="n0.23/0.06+c" + str(region[2]) + "+w20", frame="1.0");
    # Annotate with aftershock locations
    if params.aftershocks:
        [lon, lat, _, _, _] = io_additionals.read_aftershock_table(params.aftershocks);
        fig.plot(lon, lat, style='c0.1c', color='black', pen="thin,black");

    # Draw each source
    eq_lon, eq_lat = [], [];
    for source in inputs.source_object:
        source_lon, source_lat = fault_vector_functions.xy2lonlat(source.xstart, source.ystart, inputs.zerolon,
                                                                  inputs.zerolat);
        eq_lon.append(source_lon);
        eq_lat.append(source_lat);
        [x_total, y_total, _, _] = conversion_math.get_fault_four_corners(source);
        lons, lats = fault_vector_functions.xy2lonlat(x_total, y_total, inputs.zerolon, inputs.zerolat);
        if not source.potency:
            fig.plot(x=lons, y=lats, pen="thick,black");  # in case of area sources, outline them.
        else:
            fig.plot(x=lons, y=lats, style='s0.3c', color="purple", pen="thin,black");  # in case of point sources
    # Annotate with earthquake location.
    fig.plot(eq_lon, eq_lat, style='s0.05c', color="purple", pen="thin,black");
    fig.colorbar(position="JCR+w4.0i+v+o0.7i/0i", cmap="mycpt.cpt", truncate="-0.045/0.045",
                 frame=["x0.01", "y+L\"Disp(m)\""]);
    fig.savefig(outfile);
    return;


def map_displacement_vectors(params, inputs, obs_disp_points, out_object, outfile, vmin=None, vmax=None):
    """
    Make a plot of modeled vector displacement from model_disp_points
    obs_disp_points is an object that can be used to also plot observations at the same points.
    """
    if len(out_object.model_disp_points) == 0:
        return;

    proj = 'M4i'
    region = [inputs.minlon, inputs.maxlon, inputs.minlat, inputs.maxlat];

    # Unpack
    model_dE = np.array([x.dE_obs for x in out_object.model_disp_points]);
    model_dN = np.array([x.dN_obs for x in out_object.model_disp_points]);
    model_dU = np.array([x.dU_obs for x in out_object.model_disp_points]);
    model_lon = np.array([x.lon for x in out_object.model_disp_points]);
    model_lat = np.array([x.lat for x in out_object.model_disp_points]);

    # Make modeled vertical displacement color map
    [cmap_opts, cbar_opts] = utilities.define_colorbar_series(model_dU, vmin, vmax);
    pygmt.makecpt(cmap="roma", series=str(cmap_opts[0])+"/"+str(cmap_opts[1])+"/"+str(cmap_opts[2]), background="o",
                  output="mycpt.cpt");

    # Build a PyGMT plot
    fig = pygmt.Figure();
    fig.basemap(region=region, projection=proj, frame="+t\"Coseismic Displacements\"");
    fig.coast(region=region, projection=proj, borders='1', shorelines='1.0p,black', water='lightblue',
              map_scale="n0.4/0.06+c" + str(region[2]) + "+w20", frame="1.0");
    fig.plot(model_lon, model_lat, style='c0.3c', color=model_dU, cmap='mycpt.cpt', pen="thin,black");

    # Draw vectors and vector scale bar
    scale_factor = 1; scale_arrow = 0.100;   vectext = "10 cm";  # 10 cm, large vectors
    # scale_factor = 100; scale_arrow = 0.010;  vectext = "10 mm";  # 10 mm, small vectors

    scale = (np.max(np.abs(model_dE)) * scale_factor / 0.012);  # empirical scaling, convenient display
    fig.plot(x=model_lon, y=model_lat, style='v0.2c+e+gblack+h0+p1p,black+z'+str(scale),
             direction=[model_dE, model_dN], pen="thin,black");
    fig.plot(x=[region[0]+0.30], y=[region[2]+0.05],  style='v0.2c+e+gblack+h0+p1p,black+z'+str(scale),
             direction=[[scale_arrow], [0]],  pen="thin,black");  # scale vector
    fig.text(x=[region[0]+0.45], y=[region[2]+0.15], text=vectext+" model");  # scale label
    # Plot the observations if they exist
    obs_dE = np.array([x.dE_obs for x in obs_disp_points]);
    obs_dN = np.array([x.dN_obs for x in obs_disp_points]);
    if sum(~np.isnan(obs_dE)) > 0:
        fig.plot(x=model_lon, y=model_lat, style='v0.2c+e+gred+h0+p1p,red+z' + str(scale),
                 direction=[obs_dE, obs_dN], pen="thin,red");
        fig.plot(x=[region[0]+0.30], y=[region[2]+0.35],  style='v0.2c+e+gred+h0+p1p,red+z'+str(scale),
                 direction=[[scale_arrow], [0]],  pen="thin,red");  # scale vector
        fig.text(x=[region[0]+0.50], y=[region[2]+0.45], text=vectext+' obs');  # scale label

    # Draw each source
    eq_lon, eq_lat = [], [];
    for source in inputs.source_object:
        srclon, srclat = fault_vector_functions.xy2lonlat(source.xstart, source.ystart, inputs.zerolon, inputs.zerolat);
        eq_lon.append(srclon);
        eq_lat.append(srclat);
        [x_total, y_total, _, _] = conversion_math.get_fault_four_corners(source);
        lons, lats = fault_vector_functions.xy2lonlat(x_total, y_total, inputs.zerolon, inputs.zerolat);
        if not source.potency:
            fig.plot(x=lons, y=lats, pen="thick,black");  # in case of area sources, outline them.
        else:
            fig.plot(x=lons, y=lats, style='s0.3c', color="purple", pen="thin,black");  # in case of point sources
    fig.plot(eq_lon, eq_lat, style='s0.3c', color="purple", pen="thin,black");

    # Annotate with aftershock locations
    if params.aftershocks:
        [lon, lat, _, _, _] = io_additionals.read_aftershock_table(params.aftershocks);
        fig.plot(lon, lat, style='c0.1c', color='black', pen="thin,black");

    fig.colorbar(position="JCR+w4.0i+v+o0.7i/0i", cmap="mycpt.cpt", truncate=str(cbar_opts[0])+"/"+str(cbar_opts[1]),
                 frame=["x"+str(cbar_opts[2]), "y+L\"Vert Disp(m)\""]);
    fig.savefig(outfile);
    return;
