"""
Special output functions for disp_point_objects.
More specific than the basic ones listed in io_additionals.py.
"""


def write_disp_points_gmt(disp_points, filename, write_meas_type=False):
    """
    Write disp_points in GMT psvelo format.

    :param disp_points: list of disp_points_objects
    :param filename: string
    :param write_meas_type: bool, similar to 'verbose', to write an extra column with meas_type
    """
    print("Writing %s " % filename);
    ofile = open(filename, 'w');
    ofile.write("# lon lat dE dN dU Se Sn Su name\n");
    for item in disp_points:
        ofile.write("%f %f %f %f %f " % (item.lon, item.lat, item.dE_obs, item.dN_obs, item.dU_obs) );
        ofile.write("%f %f %f %s " % (item.Se_obs, item.Sn_obs, item.Su_obs, item.name));
        if write_meas_type:
            ofile.write("%s " % item.meas_type);
        ofile.write("\n");
    ofile.close();
    return;


def write_disp_points_lon_lat(disp_points, filename):
    """
    :param disp_points: a list of disp_points
    :param filename: string
    """
    print("Writing %s " % filename);
    with open(filename, 'w') as ofile:
        for item in disp_points:
            ofile.write("%f %f \n" % (item.lon, item.lat) );
    return;
