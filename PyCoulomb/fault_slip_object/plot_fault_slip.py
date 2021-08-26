"""Functions to take a slip distribution and make a map. """

import numpy as np
import pygmt
from Elastic_stresses_py.PyCoulomb.fault_slip_object import fault_slip_object


def map_source_slip_distribution(fault_dict_list, outfile, disp_points=(), region=None,
                                 scale_arrow=(1.0, 0.010, "10 mm"), v_labeling_interval=0.005, fault_traces=None,
                                 title=""):
    """
    Plot a map of slip distribution from fault_dict_list, a general format for slip distributions.
    In order to use this function with other formats, like intxt or slippy, convert to the internal fault dict first.
    """
    print("Plotting outfile %s " % outfile);
    proj = "M7i"
    if not region:  # automatically determine the region
        buffer_deg = 0.15
        lons, lats = fault_slip_object.get_four_corners_lon_lat(fault_dict_list[0]);
        region = [np.min(lons)-buffer_deg, np.max(lons)+buffer_deg, np.min(lats)-buffer_deg, np.max(lats)+buffer_deg];
    fig = pygmt.Figure();
    fig_width_deg = region[1] - region[0];
    fig.basemap(region=region, projection=proj, B="+t\"" + title + "\"");
    fig.coast(shorelines="1.0p,black", region=region, N="1", projection=proj, B=str(fig_width_deg/5));  # the boundary
    fig.coast(region=region, projection=proj, N='2', W='0.5p,black', S='lightblue');

    # Drawing fault slip with a color scale.
    fault_colors = [_i["slip"] for _i in fault_dict_list];
    if len(fault_colors) > 0:
        range_of_slip = np.max(fault_colors) - np.min(fault_colors);
    else:
        fault_colors = [0];  # some default arguments in case of empty fault_dict_list
        range_of_slip = 0;
    if range_of_slip == 0:
        range_of_slip = 0.010;  # very small slip.
    vmin = np.min(fault_colors) - (0.1*range_of_slip+0.001);
    vmax = np.max(fault_colors) + (0.1*range_of_slip+0.001);
    pygmt.makecpt(C="devon", T=str(vmin)+"/"+str(vmax)+"/"+str((vmax-vmin)/100), G="0/1", D="o", I=True,
                  H="mycpt.cpt");  # slip divided into 100
    for item in fault_dict_list:
        lons, lats = fault_slip_object.get_four_corners_lon_lat(item);
        fig.plot(x=lons, y=lats, Z=str(item["slip"]), pen="thick,black", G="+z", C="mycpt.cpt");
        fig.plot(x=lons[0:2], y=lats[0:2], pen="thickest,black", G="+z", C="mycpt.cpt");
    fig.coast(region=region, projection=proj, N='2', W='0.5p,black', L="g-124.75/40.2+c1.5+w50");

    # Optional lines you can draw on the plot
    if fault_traces:
        for item in fault_traces:
            fig.plot(x=item[0], y=item[1], pen="thickest,darkred");

    if len(disp_points) > 0:
        lon = np.array([x.lon for x in disp_points]);
        lat = np.array([x.lat for x in disp_points]);
        disp_x = np.array([x.dE_obs for x in disp_points]);
        disp_y = np.array([x.dN_obs for x in disp_points]);
        disp_z = np.array([x.dU_obs for x in disp_points]);
        lon_vert = np.array([x.lon for x in disp_points if ~np.isnan(x.dU_obs)]);
        lat_vert = np.array([x.lat for x in disp_points if ~np.isnan(x.dU_obs)]);
        disp_z_vert = np.array([x.dU_obs for x in disp_points if ~np.isnan(x.dU_obs)]);
        fig.plot(x=lon, y=lat, style='s0.07i', color='blue', pen="thin,black");
        if sum(~np.isnan(disp_z)) > 0:
            vmin_v = np.nanmin(disp_z);
            vmax_v = np.nanmax(disp_z);
            pygmt.makecpt(C="roma", T=str(vmin_v)+"/"+str(vmax_v)+"/"+str((vmax_v-vmin_v)/100), D="o", H="vert.cpt");
            fig.plot(lon_vert, lat_vert, style='c0.3c', color=disp_z_vert, C='vert.cpt', pen="thin,black");
            fig.colorbar(D="JCR+w4.0i+v+o0.7i/0i", C="vert.cpt", G=str(vmin_v*0.99)+"/"+str(vmax_v*0.99),
                         B=["x"+str(v_labeling_interval), "y+L\"Vert Disp(m)\""]);
            scale = scale_arrow[0] * (1/scale_arrow[1]);  # empirical scaling for convenient display
            fig.plot(x=lon, y=lat, style='v0.2c+e+gblack+h0+p1p,black+z'+str(scale),
                     direction=[disp_x, disp_y], pen="thin,black");
            fig.plot(x=[region[0]+0.30], y=[region[2]+0.05],  style='v0.2c+e+gblack+h0+p1p,black+z'+str(scale),
                     direction=[[scale_arrow[1]], [0]],  pen="thin,black");  # scale vector
            fig.text(x=[region[0]+0.45], y=[region[2]+0.15], text=scale_arrow[2]+' model');  # scale label

    fig.colorbar(D="jBr+w3.5i/0.2i+o2.5c/1.5c+h", C="mycpt.cpt",
                 G=str(vmin+0.001) + "/" + str(vmax-0.001), B=["x" + str(range_of_slip/10), "y+L\"Slip(m)\""]);
    fig.savefig(outfile);
    return;
