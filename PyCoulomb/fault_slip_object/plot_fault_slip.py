"""Functions to take a slip distribution and make a map. """

import numpy as np
import pygmt
from Elastic_stresses_py.PyCoulomb.fault_slip_object import fault_slip_object


def map_source_slip_distribution(fault_dict_list, outfile, disp_points=None, region=None,
                                 scale_arrow=(1.0, 0.010, "10 mm"), v_labeling_interval=0.005):
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
    fig.basemap(region=region, projection=proj, B="+t");
    fig.coast(shorelines="1.0p,black", region=region, N="1", projection=proj, B="1.0");  # the boundary
    fig.coast(region=region, projection=proj, N='2', W='0.5p,black', S='lightblue');

    # Drawing fault slip with a color scale.
    fault_colors = [_i["slip"] for _i in fault_dict_list];
    range_of_slip = np.max(fault_colors) - np.min(fault_colors);
    if range_of_slip == 0:
        range_of_slip = 0.010;  # very small slip.
    vmin = np.min(fault_colors) - 0.1*range_of_slip;
    vmax = np.max(fault_colors) + 0.1*range_of_slip;
    pygmt.makecpt(C="devon", T=str(vmin)+"/"+str(vmax)+"/"+str((vmax-vmin)/100), G="0/1", D="o", I=True,
                  H="mycpt.cpt");  # slip divided into 100
    for item in fault_dict_list:
        lons, lats = fault_slip_object.get_four_corners_lon_lat(item);
        fig.plot(x=lons, y=lats, Z=str(item["slip"]), pen="thick,black", G="+z", C="mycpt.cpt");
        fig.plot(x=lons[0:2], y=lats[0:2], pen="thickest,black", G="+z", C="mycpt.cpt");
    fig.coast(region=region, projection=proj, N='2', W='0.5p,black', L="g-124.75/40.6+c1.5+w50");

    if disp_points:
        fig.plot(x=disp_points.lon, y=disp_points.lat, style='s0.07i', color='blue', pen="thin,black");
        if disp_points.dE_obs is not None:
            disp_z = disp_points.dU_obs;
            disp_x = disp_points.dE_obs;
            disp_y = disp_points.dN_obs;
            vmin_v = np.min(disp_z);
            vmax_v = np.max(disp_z);
            pygmt.makecpt(C="roma", T=str(vmin_v)+"/"+str(vmax_v)+"/"+str((vmax_v-vmin_v)/100), D="o", H="vert.cpt");
            fig.plot(disp_points.lon, disp_points.lat, style='c0.3c', color=disp_z, C='vert.cpt', pen="thin,black");
            fig.colorbar(D="JCR+w4.0i+v+o0.7i/0i", C="vert.cpt", G=str(vmin_v*0.99)+"/"+str(vmax_v*0.99),
                         B=["x"+str(v_labeling_interval), "y+L\"Vert Disp(m)\""]);
            scale = scale_arrow[0] * (1/scale_arrow[1]);  # empirical scaling for convenient display
            fig.plot(x=disp_points.lon, y=disp_points.lat, style='v0.2c+e+gblack+h0+p1p,black+z'+str(scale),
                     direction=[disp_x, disp_y], pen="thin,black");
            fig.plot(x=[region[0]+0.30], y=[region[2]+0.05],  style='v0.2c+e+gblack+h0+p1p,black+z'+str(scale),
                     direction=[[scale_arrow[1]], [0]],  pen="thin,black");  # scale vector
            fig.text(x=[region[0]+0.45], y=[region[2]+0.15], text=scale_arrow[2]+' model');  # scale label

    fig.colorbar(D="jBr+w3.5i/0.2i+o2.5c/1.5c+h", C="mycpt.cpt",
                 G=str(vmin) + "/" + str(vmax), B=["x" + str(range_of_slip/10), "y+L\"Slip(m)\""]);
    fig.savefig(outfile);
    return;
