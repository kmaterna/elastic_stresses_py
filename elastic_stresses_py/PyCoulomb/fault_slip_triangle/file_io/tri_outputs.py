from ... import conversion_math
from tectonic_utils.geodesy import fault_vector_functions


def return_total_slip(fault_object):
    """ lambda function to get the slip of the object, for plotting."""
    return fault_object.get_total_slip()


def write_gmt_plots_cartesian(triangle_list, outfile, color_mappable=return_total_slip, verbose=True):
    """
    Write triangle edges out to file for GMT plots, in X-Y cartesian space in m.

    :param triangle_list: list of triangle fault elements
    :param outfile: string
    :param color_mappable: a simple function of one fault object that returns the plotting value
    :param verbose: print the status. default True
    """
    if verbose:
        print("Writing %d triangles to file %s " % (len(triangle_list), outfile))
    with open(outfile, 'w') as ofile:
        for item in triangle_list:
            total_slip = color_mappable(item)
            ofile.write("> -Z%f \n" % total_slip)
            ofile.write("%f %f \n" % (item.vertex1[0], item.vertex1[1]))
            ofile.write("%f %f \n" % (item.vertex2[0], item.vertex2[1]))
            ofile.write("%f %f \n" % (item.vertex3[0], item.vertex3[1]))
            ofile.write("%f %f \n" % (item.vertex1[0], item.vertex1[1]))
    return


def write_gmt_plots_geographic(triangle_list, outfile, color_mappable=return_total_slip, verbose=True):
    """
    Write triangle edges out to file for GMT plots, in lon/lat space assuming a cartesian-to-geographic transform

    :param triangle_list: list of triangle fault elements
    :param outfile: string
    :param color_mappable: a simple function of one fault object that returns the plotting value
    :param verbose: print the status. default True
    """
    if verbose:
        print("Writing %d triangles to file %s " % (len(triangle_list), outfile))
    with open(outfile, 'w') as ofile:
        for item in triangle_list:
            slip_for_coloring = color_mappable(item)
            vertex1, vertex2, vertex3 = item.get_ll_corners()
            ofile.write("> -Z%f \n" % slip_for_coloring)
            ofile.write("%f %f \n" % (vertex1[0], vertex1[1]))
            ofile.write("%f %f \n" % (vertex2[0], vertex2[1]))
            ofile.write("%f %f \n" % (vertex3[0], vertex3[1]))
            ofile.write("%f %f \n" % (vertex1[0], vertex1[1]))
    return


def write_colored_triangles(triangle_list, outfile, color_array):
    """
    Write geographic edges of triangles for psxy plotting.

    :param triangle_list: list of triangle fault objects
    :param outfile: string
    :param color_array: 1-d array with the same length as triangle_list
    """
    print("Writing %d triangles to file %s " % (len(triangle_list), outfile))
    if len(triangle_list) != len(color_array):
        raise ValueError("Error! List of fault elements and list of colors have different sizes.")
    with open(outfile, 'w') as ofile:
        for i, item in enumerate(triangle_list):
            vertex1, vertex2, vertex3 = item.get_ll_corners()
            ofile.write("> -Z%f \n" % color_array[i])
            ofile.write("%f %f \n" % (vertex1[0], vertex1[1]))
            ofile.write("%f %f \n" % (vertex2[0], vertex2[1]))
            ofile.write("%f %f \n" % (vertex3[0], vertex3[1]))
            ofile.write("%f %f \n" % (vertex1[0], vertex1[1]))
    return


def write_gmt_vertical_fault_file(fault_object_list, outfile, color_mappable=return_total_slip, strike=45):
    """
    Write the vertical coordinates of triangular fault patches (length and depth, in local coords instead of lon/lat)
    and associated slip values into a multi-segment file for plotting in GMT.
    Good for vertical faults.  Plots with depth as a negative number.
    Works for only planar-esque fault segments.

    :param fault_object_list: list of triangle fault elements
    :param outfile: string
    :param color_mappable: a simple function of one fault object that returns the plotting value
    :param strike: the approximate strike of the fault. default 45
    """
    print("Writing file %s " % outfile)
    origin_lon, origin_lat = fault_object_list[0].lon, fault_object_list[0].lat

    ofile = open(outfile, 'w')
    for fault in fault_object_list:
        slip = color_mappable(fault)
        vertex1, vertex2, vertex3 = fault.get_ll_corners()  # vertex1 = [lon1, lat1]
        x1, y1 = fault_vector_functions.latlon2xy_single(vertex1[0], vertex1[1], origin_lon, origin_lat)
        x2, y2 = fault_vector_functions.latlon2xy_single(vertex2[0], vertex2[1], origin_lon, origin_lat)
        x3, y3 = fault_vector_functions.latlon2xy_single(vertex3[0], vertex3[1], origin_lon, origin_lat)
        [xprime1, _] = conversion_math.rotate_points(x1, y1, 90 + strike)
        [xprime2, _] = conversion_math.rotate_points(x2, y2, 90 + strike)
        [xprime3, _] = conversion_math.rotate_points(x3, y3, 90 + strike)
        ofile.write("> -Z"+str(slip)+"\n")  # whatever slip value the user chooses
        ofile.write("%f %f\n" % (xprime1[0], -fault.vertex1[2]/1000))
        ofile.write("%f %f\n" % (xprime2[0], -fault.vertex2[2]/1000))
        ofile.write("%f %f\n" % (xprime3[0], -fault.vertex3[2]/1000))
        ofile.write("%f %f\n" % (xprime1[0], -fault.vertex1[2]/1000))
    ofile.close()
    return


def write_centroid_with_value(fault_object_list, outfile, color_array):
    """
    Write the lon/lat of the centroid of each triangle, with a color array.

    :param fault_object_list: list of triangular fault objects
    :param outfile: string
    :param color_array: list of floats, same length as fault_object_list
    """
    print("Writing %d triangles to file %s " % (len(fault_object_list), outfile))
    if len(fault_object_list) != len(color_array):
        raise ValueError("Error! List of fault elements and list of colors have different sizes.")
    with open(outfile, 'w') as ofile:
        for i, fault in enumerate(fault_object_list):
            centroid_ll = fault.get_llz_point(fault.compute_triangle_centroid())
            ofile.write("%f %f %f\n" % (centroid_ll[0], centroid_ll[1], color_array[i]))
    return
