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
    :returns: a list of Displacement_points objects
    """
    print("Reading file %s " % filename)
    disp_pts = []
    with open(filename, 'r') as ifile:
        for line in ifile:
            if line.split()[0] == "#":
                continue
            else:
                temp = line.split()
                new_disp = Displacement_points(lon=float(temp[0]), lat=float(temp[1]), dE_obs=float(temp[2])/1000,
                                               dN_obs=float(temp[3])/1000, dU_obs=float(temp[4])/1000,
                                               Se_obs=float(temp[5])/1000,
                                               Sn_obs=float(temp[6])/1000, Su_obs=float(temp[7])/1000, name=temp[8])
                disp_pts.append(new_disp)
    return disp_pts


def write_disp_points_gmt(disp_points, filename, write_meas_type=False, multiply_by=1):
    """
    Write disp_points in GMT psvelo format.

    :param disp_points: list of disp_points_objects
    :param filename: string
    :param write_meas_type: bool, similar to 'verbose', to write an extra column with meas_type
    :param multiply_by: a unit conversion for displacements and uncertainties. Default 1.
    """
    print("Writing %s " % filename)
    with open(filename, 'w') as ofile:
        ofile.write("# lon lat dE dN dU Se Sn Su name\n")
        for item in disp_points:
            ofile.write("%f %f %f %f %f " % (item.lon, item.lat, multiply_by*item.dE_obs, multiply_by*item.dN_obs,
                                             multiply_by*item.dU_obs))
            ofile.write("%f %f %f %s " % (multiply_by*item.Se_obs, multiply_by*item.Sn_obs, multiply_by*item.Su_obs,
                                          item.name))
            if write_meas_type:
                ofile.write("%s " % item.meas_type)
            ofile.write("\n")
    return


def write_disp_points_lon_lat(disp_points, filename):
    """
    :param disp_points: a list of disp_points
    :param filename: string
    """
    print("Writing %s " % filename)
    with open(filename, 'w') as ofile:
        for item in disp_points:
            ofile.write("%f %f \n" % (item.lon, item.lat))
    return
