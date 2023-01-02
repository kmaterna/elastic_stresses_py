
import numpy as np
from Tectonic_Utils.geodesy import fault_vector_functions, haversine
from Tectonic_Utils.seismo import moment_calculations
from .. import conversion_math
from ..fault_slip_object.io_pycoulomb import fault_dict_to_coulomb_fault

class TriangleFaultDict:
    """
    The internal format is a dictionary for a triangular fault segment:
    {
        vertex1 [x, y, z] in m, z is positive downward
        vertex2 [x, y, z] in m, z is positive downward
        vertex3 [x, y, z] in m, z is positive downward
        lon(reference coord),
        lat(reference coord),
        depth (km),
        rt_lat slip(m),
        dip_slip(m),
        tensile(m),
        segment(int)
    }
    If the fault is a receiver fault, we put slip = 0
    """

    def __init__(self, vertex1, vertex2, vertex3, lon, lat, depth, rtlat_slip=0, dip_slip=0, tensile=0, segment=0):
        self.vertex1 = vertex1;  # np.ndarray
        self.vertex2 = vertex2;  # np.ndarray
        self.vertex3 = vertex3;  # np.ndarray
        self.lon = lon;  # float
        self.lat = lat;  # float
        self.depth = depth;  # float
        self.rtlat_slip = rtlat_slip;  # float
        self.dip_slip = dip_slip;  # float
        self.tensile = tensile;  # float
        self.segment = segment;  # int

    def get_total_slip(self):
        """Helper function to return the total slip amount of a fault object (always > 0)"""
        return np.sqrt(np.square(self.dip_slip) + np.square(self.rtlat_slip));

    def get_total_slip_mm(self):
        """Helper function to return the total slip amount of a fault object in mm (always > 0)"""
        return self.get_total_slip()*1000;

    def get_rtlat_slip(self):
        """Helper function to return the right lateral slip amount of a fault object"""
        return self.rtlat_slip;

    def get_dip_slip(self):
        """Helper function to return the dip slip amount of a fault object"""
        return self.dip_slip;

    def get_depth(self):
        """Helper function to return the depth of a fault object"""
        return self.depth;

    def get_segment(self):
        """Helper function to return the segment_num of a fault object"""
        return self.segment;

    def get_ll_corners(self):
        # Convert m to km before coordinate conversion
        # Returns 3 tuples (lon, lat)
        lon1, lat1 = fault_vector_functions.xy2lonlat_single(self.vertex1[0] / 1000, self.vertex1[1] / 1000,
                                                             self.lon, self.lat);
        lon2, lat2 = fault_vector_functions.xy2lonlat_single(self.vertex2[0] / 1000, self.vertex2[1] / 1000,
                                                             self.lon, self.lat);
        lon3, lat3 = fault_vector_functions.xy2lonlat_single(self.vertex3[0] / 1000, self.vertex3[1] / 1000,
                                                             self.lon, self.lat);
        vertex1 = [lon1, lat1]; vertex2 = [lon2, lat2]; vertex3 = [lon3, lat3];
        return vertex1, vertex2, vertex3;

    def compute_triangle_centroid(self):
        """Compute cartesian centroid of a triangular fault. Returns a numpy array."""
        x_centroid = np.mean([self.vertex1[0], self.vertex2[0], self.vertex3[0]]);
        y_centroid = np.mean([self.vertex1[1], self.vertex2[1], self.vertex3[1]]);
        z_centroid = np.mean([self.vertex1[2], self.vertex2[2], self.vertex3[2]]);
        return np.array([x_centroid, y_centroid, z_centroid]);

    def get_fault_element_distance(self, fault_dict2, threedimensional=True):
        """Return the distance between two triangular fault elements, in 3d, in km"""
        [xc1, yc1, zc1] = self.compute_triangle_centroid();
        [xc2, yc2, zc2] = fault_dict2.compute_triangle_centroid();
        lon_t1, lat_t1 = fault_vector_functions.xy2lonlat_single(xc1/1000, yc1/1000, self.lon, self.lat)
        lon_t2, lat_t2 = fault_vector_functions.xy2lonlat_single(xc2/1000, yc2/1000, fault_dict2.lon, fault_dict2.lat)
        h_distance = haversine.distance([lat_t1, lon_t1], [lat_t2, lon_t2]);
        if threedimensional:
            depth_distance = (zc1 - zc2)/1000;
        else:
            depth_distance = 0;
        distance_3d = np.sqrt(np.square(h_distance) + np.square(depth_distance))
        return distance_3d;

    def change_fault_slip(self, rtlat=None, dipslip=None, tensile=None):
        new_rtlat = self.rtlat_slip if rtlat is None else rtlat;
        new_dipslip = self.dip_slip if dipslip is None else dipslip;
        new_tensile = self.tensile if tensile is None else tensile;
        new_fault_dict = TriangleFaultDict(lon=self.lon, lat=self.lat, segment=self.segment, depth=self.depth,
                                           vertex1=self.vertex1, vertex2=self.vertex2, vertex3=self.vertex3,
                                           dip_slip=new_dipslip, rtlat_slip=new_rtlat, tensile=new_tensile);
        return new_fault_dict;

    def change_reference_loc(self, new_refcoords=None):
        """Move the reference lon/lat.  If none given, moves it to the location of vertex1

        :param new_refcoords: tuple of (lon, lat).
        """
        vertex1_ll, vertex2_ll, vertex3_ll = self.get_ll_corners();
        if new_refcoords is None:
            new_reflon, new_reflat = vertex1_ll[0], vertex1_ll[1];  # default behavior sends the reference to vertex1
        else:
            new_reflon, new_reflat = new_refcoords[0], new_refcoords[1];
        x1, y1 = fault_vector_functions.latlon2xy(vertex1_ll[0], vertex1_ll[1], new_reflon, new_reflat);
        x2, y2 = fault_vector_functions.latlon2xy(vertex2_ll[0], vertex2_ll[1], new_reflon, new_reflat);
        x3, y3 = fault_vector_functions.latlon2xy(vertex3_ll[0], vertex3_ll[1], new_reflon, new_reflat);
        vertex1_newref = np.array([x1*1000, y1*1000, self.vertex1[2]]);
        vertex2_newref = np.array([x2*1000, y2*1000, self.vertex2[2]]);
        vertex3_newref = np.array([x3*1000, y3*1000, self.vertex3[2]]);
        new_fault_dict = TriangleFaultDict(lon=new_reflon, lat=new_reflat, segment=self.segment, depth=self.depth,
                                           vertex1=vertex1_newref, vertex2=vertex2_newref, vertex3=vertex3_newref,
                                           dip_slip=self.dip_slip, rtlat_slip=self.rtlat_slip, tensile=self.tensile);
        return new_fault_dict;

    def get_fault_area(self):
        """Compute area of a triangle: A = 1/2 (AB x AC), with all cartesian coordinates in meters. Units: sq-m."""
        AB = np.array([self.vertex2[0] - self.vertex1[0], self.vertex2[1] - self.vertex1[1],
                       self.vertex2[2] - self.vertex1[2]]);
        AC = np.array([self.vertex3[0] - self.vertex1[0], self.vertex3[1] - self.vertex1[1],
                       self.vertex3[2] - self.vertex1[2]]);
        xprod = np.cross(AB, AC);
        total_vector_mag = np.sqrt(np.square(xprod[0]) + np.square(xprod[1]) + np.square(xprod[2]));
        area_sq_m = (1/2) * total_vector_mag;
        return area_sq_m;

    def get_fault_moment(self, mu=30e9):
        """Get the moment for a triangular fault patch. Units: Newton-meters"""
        A = self.get_fault_area();
        d = self.get_total_slip();
        moment = moment_calculations.moment_from_muad(mu, A, d);
        return moment;


def convert_rectangle_into_two_triangles(one_fault_dict):
    """
    Convert a fault_dict into two triangular fault dicts. The fault normals are expected to point up.
    """
    [source] = fault_dict_to_coulomb_fault([one_fault_dict]);
    [x_all, y_all, _, _] = conversion_math.get_fault_four_corners(source);  # This is cartesian
    top_depth, bottom_depth = source.top, source.bottom;
    vertex1 = np.array([x_all[0]*1000, y_all[0]*1000, top_depth*1000]);  # in meters
    vertex2 = np.array([x_all[1]*1000, y_all[1]*1000, top_depth*1000]);
    vertex3 = np.array([x_all[2]*1000, y_all[2]*1000, bottom_depth*1000]);
    vertex4 = np.array([x_all[3]*1000, y_all[3]*1000, bottom_depth*1000]);
    first_triangle = TriangleFaultDict(lon=one_fault_dict['lon'], lat=one_fault_dict['lat'],
                                       segment=one_fault_dict['segment'],
                                       tensile=source.tensile, vertex1=vertex1, vertex2=vertex3,
                                       vertex3=vertex2, dip_slip=source.reverse, rtlat_slip=source.rtlat,
                                       depth=vertex1[2]/1000);
    second_triangle = TriangleFaultDict(lon=one_fault_dict['lon'], lat=one_fault_dict['lat'],
                                        segment=one_fault_dict['segment'],
                                        tensile=source.tensile, vertex1=vertex1, vertex2=vertex4,
                                        vertex3=vertex3, dip_slip=source.reverse, rtlat_slip=source.rtlat,
                                        depth=vertex1[2]/1000);
    list_of_two_triangles = [first_triangle, second_triangle];
    return list_of_two_triangles;


def return_total_slip(fault_dict):
    """ lambda function to get the slip of the object, for plotting """
    return fault_dict.get_total_slip();


"""
Functions that work on lists of fault items
"""


def get_total_moment(fault_dict_object_list, mu=30e9):
    """
    Return total moment of a list of triangle slip objects
    Moment in newton-meters
    """
    total_moment = 0;
    for item in fault_dict_object_list:
        total_moment += moment_calculations.moment_from_muad(mu, item.get_fault_area(), item.get_total_slip());
    return total_moment;


def check_consistent_reference_frame(triangle_fault_list):
    """Returns True if all the triangles have the coordinate reference system"""
    return_code = 1;  # default true
    reflon, reflat = triangle_fault_list[0].lon, triangle_fault_list[0].lat;
    for item in triangle_fault_list:
        if item.lon != reflon or item.lat != reflat:
            return_code = 0;
    if return_code == 0:
        print("Warning! Not all triangles have the same reference lon/lat");
    return return_code;


def write_gmt_plots_cartesian(triangle_list, outfile, color_mappable=return_total_slip):
    """
    Write triangle edges out to file for GMT plots, in X-Y cartesian space in m
    color_mappable is a simple function of one fault dictionary object that returns the plotting value
    """
    print("Writing %d triangles to file %s " % (len(triangle_list), outfile) );
    with open(outfile, 'w') as ofile:
        for item in triangle_list:
            total_slip = color_mappable(item);
            ofile.write("> -Z%f \n" % total_slip);
            ofile.write("%f %f \n" % (item.vertex1[0], item.vertex1[1]));
            ofile.write("%f %f \n" % (item.vertex2[0], item.vertex2[1]));
            ofile.write("%f %f \n" % (item.vertex3[0], item.vertex3[1]));
            ofile.write("%f %f \n" % (item.vertex1[0], item.vertex1[1]));
    return;


def write_gmt_plots_geographic(triangle_list, outfile, color_mappable=return_total_slip):
    """
    Write triangle edges out to file for GMT plots, in lon/lat space assuming a cartesian-to-geographic transform
    color_mappable is a simple function of one fault dictionary object that returns the plotting value
    """
    print("Writing %d triangles to file %s " % (len(triangle_list), outfile) );
    with open(outfile, 'w') as ofile:
        for item in triangle_list:
            slip_for_coloring = color_mappable(item);
            vertex1, vertex2, vertex3 = item.get_ll_corners();
            ofile.write("> -Z%f \n" % slip_for_coloring);
            ofile.write("%f %f \n" % (vertex1[0], vertex1[1]));
            ofile.write("%f %f \n" % (vertex2[0], vertex2[1]));
            ofile.write("%f %f \n" % (vertex3[0], vertex3[1]));
            ofile.write("%f %f \n" % (vertex1[0], vertex1[1]));
    return;


def write_gmt_vertical_fault_file(fault_dict_list, outfile, color_mappable=return_total_slip, strike=45):
    """
    Write the vertical coordinates of triangular fault patches (length and depth, in local coords instead of lon/lat)
    and associated slip values into a multi-segment file for plotting in GMT.
    Good for vertical faults.  Plots with depth as a negative number.
    Works for only planar-esque fault segments.
    Strike is the approximate strike of the fault
    """
    print("Writing file %s " % outfile);
    origin_lon, origin_lat = fault_dict_list[0].lon, fault_dict_list[0].lat;

    ofile = open(outfile, 'w');
    for fault in fault_dict_list:
        slip = color_mappable(fault);
        vertex1, vertex2, vertex3 = fault.get_ll_corners();  # vertex1 = [lon1, lat1]
        x1, y1 = fault_vector_functions.latlon2xy_single(vertex1[0], vertex1[1], origin_lon, origin_lat);
        x2, y2 = fault_vector_functions.latlon2xy_single(vertex2[0], vertex2[1], origin_lon, origin_lat);
        x3, y3 = fault_vector_functions.latlon2xy_single(vertex3[0], vertex3[1], origin_lon, origin_lat);
        [xprime1, _] = conversion_math.rotate_points(x1, y1, 90+strike);
        [xprime2, _] = conversion_math.rotate_points(x2, y2, 90 + strike);
        [xprime3, _] = conversion_math.rotate_points(x3, y3, 90 + strike);
        ofile.write("> -Z"+str(slip)+"\n");  # whatever slip value the user chooses
        ofile.write("%f %f\n" % (xprime1[0], -fault.vertex1[2]/1000));
        ofile.write("%f %f\n" % (xprime2[0], -fault.vertex2[2]/1000));
        ofile.write("%f %f\n" % (xprime3[0], -fault.vertex3[2]/1000));
        ofile.write("%f %f\n" % (xprime1[0], -fault.vertex1[2]/1000));
    ofile.close();
    return;
