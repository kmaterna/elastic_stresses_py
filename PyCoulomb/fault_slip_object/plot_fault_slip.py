"""Functions to take a slip distribution and make a map. """

import numpy as np
import pygmt
from Elastic_stresses_py.PyCoulomb.fault_slip_object import fault_slip_object


def map_source_slip_distribution(fault_dict_list, outfile, disp_points=None, region=None):
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
    fig.coast(shorelines="1.0p,black", region=region, N="1", projection=proj, B="1.0");  # the boundary.
    fig.coast(region=region, projection=proj, N='2', W='0.5p,black', S='white');
    fault_colors = [_i["slip"] for _i in fault_dict_list];
    vmin = np.min(fault_colors) - 0.03;
    vmax = np.max(fault_colors) + 0.03;
    pygmt.makecpt(C="roma", T=str(vmin)+"/"+str(vmax)+"/"+str((vmax-vmin)/100), G="0/0.5", D="o", I=True,
                  H="mycpt.cpt");
    for item in fault_dict_list:
        lons, lats = fault_slip_object.get_four_corners_lon_lat(item);
        fig.plot(x=lons, y=lats, Z=str(item["slip"]), pen="thick,black", G="+z", C="mycpt.cpt");
    fig.coast(region=region, projection=proj, N='2', W='0.5p,black', L="g-125.5/39.6+c1.5+w50");
    if disp_points:
        fig.plot(x=disp_points.lon, y=disp_points.lat, style='s0.07i', color='blue', pen="thin,black");
    fig.colorbar(D="jBr+w3.5i/0.2i+o2.5c/1.5c+h", C="mycpt.cpt", I="0.8",
                 G=str(vmin) + "/" + str(vmax), B=["x" + str(0.05), "y+L\"Slip(m)\""]);
    fig.savefig(outfile);
    return;
