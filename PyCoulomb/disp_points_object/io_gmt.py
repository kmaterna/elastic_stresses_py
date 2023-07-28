"""
Special i/o functions for disp_point_objects.
More specific than the basic ones listed in io_additionals.py.
"""


from .disp_points_object import Displacement_points

def read_disp_points_gmt(filename):
    """
    Read disp_points from GMT psvelo format.
    # lon lat dE dN dU Se Sn Su name\n

    :param filename: string
    """
    disp_pts = [];
    ifile = open(filename, 'r');
    for line in ifile:
        if line.split()[0] == "#":
            continue;
        else:
            temp = line.split();
            new_disp = Displacement_points(lon=float(temp[0]), lat=float(temp[1]), dE_obs=float(temp[2]),
                                           dN_obs=float(temp[3]), Se_obs=float(temp[4]), Sn_obs=0,
                                           dU_obs=0, Su_obs=0, name=temp[8]);
            disp_pts.append(new_disp);
    return disp_pts;


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
