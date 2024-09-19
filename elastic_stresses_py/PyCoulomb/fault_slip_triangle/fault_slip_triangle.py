
import numpy as np
from Tectonic_Utils.geodesy import fault_vector_functions, haversine
from Tectonic_Utils.seismo import moment_calculations
from ..pyc_fault_object import Faults_object
from ..fault_slip_object.fault_slip_object import fault_object_to_coulomb_fault


class TriangleFault:
    """
    The internal format is a class for a triangular fault segment, specified by three vertices and a reference point.
    If the fault is a receiver fault, we put slip = 0.

    :param vertex1: array containing [x, y, z] of vertex 1 in m, z is positive downward. z in meters below surface.
    :param vertex2: array containing [x, y, z] of vertex 2 in m, z is positive downward. z in meters below surface.
    :param vertex3: array containing [x, y, z] of vertex 3 in m, z is positive downward. z in meters below surface.
    :param lon: longitude of reference coordinate
    :type lon: float
    :param lat: latitude of reference coordinate
    :type lat: float
    :param depth: depth of reference coordinate, in km (not used for much)
    :type depth: float
    :param rtlat_slip: amount of right-lateral slip, in meters
    :type rtlat_slip: float
    :param dip_slip: amount of dip slip, in meters
    :type dip_slip: float
    :param tensile: amount of tensile opening, in meters
    :type tensile: float
    :param segment: segment number for this fault
    :type segment: int
    """

    def __init__(self, vertex1, vertex2, vertex3, lon, lat, depth, rtlat_slip=0, dip_slip=0, tensile=0, segment=0):
        self.vertex1 = vertex1  # np.ndarray
        self.vertex2 = vertex2  # np.ndarray
        self.vertex3 = vertex3  # np.ndarray
        self.lon = lon  # float
        self.lat = lat  # float
        self.depth = depth  # float
        self.rtlat_slip = rtlat_slip  # float
        self.dip_slip = dip_slip  # float
        self.tensile = tensile  # float
        self.segment = segment  # int

    def get_total_slip(self) -> float:
        """Helper function to return the total slip amount of a fault object (always > 0)."""
        return np.sqrt(np.square(self.dip_slip) + np.square(self.rtlat_slip))

    def get_total_slip_mm(self) -> float:
        """Helper function to return the total slip amount of a fault object in mm (always > 0)."""
        return self.get_total_slip()*1000

    def get_rtlat_slip(self):
        """Helper function to return the right lateral slip amount of a fault object."""
        return self.rtlat_slip

    def get_dip_slip(self):
        """Helper function to return the dip slip amount of a fault object."""
        return self.dip_slip

    def get_fault_depth(self):
        """Helper function to return the depth of a fault object."""
        return self.depth

    def get_centroid_depth(self) -> float:
        """Helper function to return the computed centroid depth of a fault object, in meters, positive downwards."""
        return float(self.compute_triangle_centroid()[2])

    def get_segment(self) -> int:
        """Helper function to return the segment_num of a fault object."""
        return self.segment

    def get_llz_point(self, meters_pt):
        """
        meters_pt is a length-three tuple. XY are in meters from the reference pixel. Z(depth) is in meters.
        Returns array of [lon, lat, depth_km].
        """
        lon1, lat1 = fault_vector_functions.xy2lonlat_single(meters_pt[0] / 1000, meters_pt[1] / 1000, self.lon,
                                                             self.lat)
        return np.array([lon1, lat1, (meters_pt[2] / 1000)])

    def get_llz_tri_coords(self, tri_coords_pt):
        """tri_coords_pt is a length-two tuple between 0 and 1"""
        target_pt_meters = np.add(self.vertex1, tri_coords_pt[0] * self.get_vector_12())
        target_pt_meters = np.add(target_pt_meters, tri_coords_pt[1] * self.get_vector_13())
        return self.get_llz_point(target_pt_meters)

    def get_ll_corners(self):
        """Return three tuples representing (lon, lat) pairs of all three vertices for mapping."""
        # Convert m to km before coordinate conversion
        v1 = self.get_llz_point(self.vertex1)
        v2 = self.get_llz_point(self.vertex2)
        v3 = self.get_llz_point(self.vertex3)
        vertex1 = [v1[0], v1[1]]
        vertex2 = [v2[0], v2[1]]
        vertex3 = [v3[0], v3[1]]
        return vertex1, vertex2, vertex3

    def compute_triangle_centroid(self):
        """
        Compute cartesian centroid of a triangular fault in m, relative to self.lon and self.lat and the surface.
        Returns a numpy array [x, y, z] in meters
        """
        x_centroid = np.mean([self.vertex1[0], self.vertex2[0], self.vertex3[0]])
        y_centroid = np.mean([self.vertex1[1], self.vertex2[1], self.vertex3[1]])
        z_centroid = np.mean([self.vertex1[2], self.vertex2[2], self.vertex3[2]])
        return np.array([x_centroid, y_centroid, z_centroid])

    def get_fault_element_distance(self, fault_object2, threedimensional=True) -> float:
        """Return the distance between two triangular fault elements, in 3d, in km."""
        [xc1, yc1, zc1] = self.compute_triangle_centroid()
        [xc2, yc2, zc2] = fault_object2.compute_triangle_centroid()
        lon_t1, lat_t1 = fault_vector_functions.xy2lonlat_single(xc1/1000, yc1/1000, self.lon, self.lat)
        lon_t2, lat_t2 = fault_vector_functions.xy2lonlat_single(xc2/1000, yc2/1000, fault_object2.lon,
                                                                 fault_object2.lat)
        h_distance = haversine.distance([lat_t1, lon_t1], [lat_t2, lon_t2])
        if threedimensional:
            depth_distance = (zc1 - zc2)/1000
        else:
            depth_distance = 0
        distance_3d = np.sqrt(np.square(h_distance) + np.square(depth_distance))
        return distance_3d

    def change_fault_slip(self, rtlat=None, dipslip=None, tensile=None):
        new_rtlat = self.rtlat_slip if rtlat is None else rtlat
        new_dipslip = self.dip_slip if dipslip is None else dipslip
        new_tensile = self.tensile if tensile is None else tensile
        new_fault_obj = TriangleFault(lon=self.lon, lat=self.lat, segment=self.segment, depth=self.depth,
                                      vertex1=self.vertex1, vertex2=self.vertex2, vertex3=self.vertex3,
                                      dip_slip=new_dipslip, rtlat_slip=new_rtlat, tensile=new_tensile)
        return new_fault_obj

    def change_reference_loc(self, new_refcoords=None):
        """Move the reference lon/lat.  If none given, moves it to the location of vertex1.

        :param new_refcoords: tuple of (lon, lat).
        """
        vertex1_ll, vertex2_ll, vertex3_ll = self.get_ll_corners()
        if new_refcoords is None:
            new_reflon, new_reflat = vertex1_ll[0], vertex1_ll[1]  # default behavior sends the reference to vertex1
        else:
            new_reflon, new_reflat = new_refcoords[0], new_refcoords[1]
        x1, y1 = fault_vector_functions.latlon2xy(vertex1_ll[0], vertex1_ll[1], new_reflon, new_reflat)
        x2, y2 = fault_vector_functions.latlon2xy(vertex2_ll[0], vertex2_ll[1], new_reflon, new_reflat)
        x3, y3 = fault_vector_functions.latlon2xy(vertex3_ll[0], vertex3_ll[1], new_reflon, new_reflat)
        vertex1_newref = np.array([x1*1000, y1*1000, self.vertex1[2]])
        vertex2_newref = np.array([x2*1000, y2*1000, self.vertex2[2]])
        vertex3_newref = np.array([x3*1000, y3*1000, self.vertex3[2]])
        new_fault_obj = TriangleFault(lon=new_reflon, lat=new_reflat, segment=self.segment, depth=self.depth,
                                      vertex1=vertex1_newref, vertex2=vertex2_newref, vertex3=vertex3_newref,
                                      dip_slip=self.dip_slip, rtlat_slip=self.rtlat_slip, tensile=self.tensile)
        return new_fault_obj

    def get_fault_area(self) -> float:
        """Compute area of a triangle: A = 1/2 (AB x AC), with all cartesian coordinates in meters. Units: sq-m."""
        AB = np.array([self.vertex2[0] - self.vertex1[0], self.vertex2[1] - self.vertex1[1],
                       self.vertex2[2] - self.vertex1[2]])
        AC = np.array([self.vertex3[0] - self.vertex1[0], self.vertex3[1] - self.vertex1[1],
                       self.vertex3[2] - self.vertex1[2]])
        xprod = np.cross(AB, AC)
        total_vector_mag = np.sqrt(np.square(xprod[0]) + np.square(xprod[1]) + np.square(xprod[2]))
        area_sq_m = (1/2) * total_vector_mag
        return area_sq_m

    def get_fault_moment(self, mu=30e9) -> float:
        """Get the moment for a triangular fault patch. Units: Newton-meters."""
        A = self.get_fault_area()
        d = self.get_total_slip()
        moment = moment_calculations.moment_from_muad(mu, A, d)
        return moment

    def get_fault_normal(self):
        """
        Return the upward unit-normal-vector to the plane of the triangle
        """
        vector1 = self.get_vector_12()
        vector2 = self.get_vector_13()
        fault_normal = fault_vector_functions.get_unit_vector(np.cross(vector1, vector2))
        if fault_normal[2] > 0:
            fault_normal = np.multiply(fault_normal, -1)
        return fault_normal

    def get_strike(self) -> float:
        fault_normal = self.get_fault_normal()
        strike = np.rad2deg(np.arctan2(fault_normal[0], fault_normal[1])) - 90
        if strike < 0:
            strike = strike + 360
        return strike

    def get_dip(self) -> float:
        dip = np.rad2deg(np.arccos(np.dot(self.get_fault_normal(), [0, 0, -1])))
        if dip > 90:
            raise ValueError("Error! Dip returned greater than 90 degrees.")
        if dip < 0:
            raise ValueError("Error! Dip returned less than 0 degrees.")
        return dip

    def get_vector_12(self):
        vector1 = np.array([self.vertex2[0] - self.vertex1[0],
                            self.vertex2[1] - self.vertex1[1],
                            self.vertex2[2] - self.vertex1[2]])
        return vector1

    def get_vector_13(self):
        vector2 = np.array([self.vertex3[0] - self.vertex1[0],
                            self.vertex3[1] - self.vertex1[1],
                            self.vertex3[2] - self.vertex1[2]])
        return vector2

    def fault_triangle_to_coulomb_fault(self, zerolon_system=None, zerolat_system=None):
        """
        Convert a triangular fault object into a source object for Elastic_stresses_py.
        It will be a small rectangle centered around the centroid of the triangle
        By default, the bottom corner of the fault is the center of the coordinate system, but
        Parameters zerolon_system and zerolat_system can be passed in.
        """
        zerolon = self.lon if not zerolon_system else zerolon_system
        zerolat = self.lat if not zerolat_system else zerolat_system
        rake = fault_vector_functions.get_rake(self.rtlat_slip, self.dip_slip)

        centroid = self.compute_triangle_centroid()  # in meters
        [refx0, refy0] = fault_vector_functions.latlon2xy_single(self.lon, self.lat, zerolon, zerolat)  # in km
        [xstart, ystart] = refx0 + centroid[0]/1000, refy0 + centroid[1]/1000  # in km with respect to zerolon/lat
        xfinish, yfinish = fault_vector_functions.add_vector_to_point(xstart, ystart, 2, self.get_strike())

        one_source = Faults_object(xstart=xstart, xfinish=xfinish, ystart=ystart, yfinish=yfinish,
                                   rtlat=self.rtlat_slip, reverse=self.dip_slip, tensile=self.tensile,
                                   strike=self.get_strike(), dipangle=self.get_dip(), rake=rake, zerolon=zerolon,
                                   zerolat=zerolat, top=centroid[2]/1000, bottom=centroid[2]/1000+1.2)
        return one_source

    # ------------ PREDICATES -------------- #
    def is_within_depth_range(self, upper_depth, lower_depth) -> bool:
        """
        Filter fault_object if it falls within depth range [upper_depth, lower_depth].

        :param upper_depth: float, km
        :param lower_depth: float, km
        :returns: bool
        """
        centroid_depth_km = self.get_centroid_depth()/1000
        if upper_depth <= centroid_depth_km <= lower_depth:
            return True
        else:
            return False


def convert_rectangle_into_two_triangles(one_fault_obj):
    """
    Convert a rectangular fault_slip_object into two triangular faults. The fault normals are expected to point up.
    """
    [source] = fault_object_to_coulomb_fault([one_fault_obj])
    list_of_two_triangles = convert_pycoulomb_rectangle_into_two_triangles(source, one_fault_obj.lon, one_fault_obj.lat)
    return list_of_two_triangles


def convert_pycoulomb_rectangle_into_two_triangles(source, startlon, startlat):
    """
    Convert one rectangular pycoulomb_fault into two triangular faults. The fault normals are expected to point up.
    """
    [x_all, y_all, _, _] = source.get_fault_four_corners()  # This is cartesian
    top_depth, bottom_depth = source.top, source.bottom
    vertex1 = np.array([x_all[0]*1000, y_all[0]*1000, top_depth*1000])  # in meters
    vertex2 = np.array([x_all[1]*1000, y_all[1]*1000, top_depth*1000])
    vertex3 = np.array([x_all[2]*1000, y_all[2]*1000, bottom_depth*1000])
    vertex4 = np.array([x_all[3]*1000, y_all[3]*1000, bottom_depth*1000])
    first_triangle = TriangleFault(lon=startlon, lat=startlat, segment=source.segment,
                                   tensile=source.tensile, vertex1=vertex1, vertex2=vertex3, vertex3=vertex2,
                                   dip_slip=source.reverse, rtlat_slip=source.rtlat, depth=float(vertex1[2])/1000)
    second_triangle = TriangleFault(lon=startlon, lat=startlat, segment=source.segment,
                                    tensile=source.tensile, vertex1=vertex1, vertex2=vertex4, vertex3=vertex3,
                                    dip_slip=source.reverse, rtlat_slip=source.rtlat, depth=float(vertex1[2])/1000)
    list_of_two_triangles = [first_triangle, second_triangle]
    return list_of_two_triangles


def convert_pycoulomb_rectangle_into_three_triangles(source, startlon, startlat):
    """
    Convert one rectangular pycoulomb_fault into three triangular faults. The fault normals are expected to point up.
    """
    [x_all, y_all, _, _] = source.get_fault_four_corners()  # This is cartesian
    top_depth, bottom_depth = source.top, source.bottom
    vertex1 = np.array([x_all[0]*1000, y_all[0]*1000, top_depth*1000])  # in meters
    vertex2 = np.array([x_all[1]*1000, y_all[1]*1000, top_depth*1000])
    vertex3 = np.array([x_all[2]*1000, y_all[2]*1000, bottom_depth*1000])
    vertex4 = np.array([x_all[3]*1000, y_all[3]*1000, bottom_depth*1000])
    far_midpoint = np.array([1000 * np.mean([x_all[1], x_all[2]]),
                             1000 * np.mean([y_all[1], y_all[2]]),
                             1000 * np.mean([top_depth, bottom_depth])])
    first_triangle = TriangleFault(lon=startlon, lat=startlat, segment=source.segment, tensile=source.tensile,
                                   vertex1=vertex1, vertex2=far_midpoint, vertex3=vertex2,
                                   dip_slip=source.reverse, rtlat_slip=source.rtlat, depth=float(vertex1[2])/1000)
    second_triangle = TriangleFault(lon=startlon, lat=startlat, segment=source.segment, tensile=source.tensile,
                                    vertex1=vertex1, vertex2=vertex4, vertex3=far_midpoint,
                                    dip_slip=source.reverse, rtlat_slip=source.rtlat, depth=float(vertex1[2])/1000)
    third_triangle = TriangleFault(lon=startlon, lat=startlat, segment=source.segment, tensile=source.tensile,
                                   vertex1=vertex4, vertex2=vertex3, vertex3=far_midpoint,
                                   dip_slip=source.reverse, rtlat_slip=source.rtlat, depth=float(vertex1[2])/1000)
    list_of_triangles = [first_triangle, second_triangle, third_triangle]
    return list_of_triangles


"""
Functions that work on lists of fault items
"""


def get_total_moment(fault_object_list, mu=30e9) -> float:
    """
    Return total moment of a list of triangle slip objects. Moment in newton-meters.
    """
    total_moment = 0
    for item in fault_object_list:
        total_moment += moment_calculations.moment_from_muad(mu, item.get_fault_area(), item.get_total_slip())
    return total_moment


def get_total_moment_depth_dependent(fault_object_list, depths, mus) -> float:
    """Compute total moment using a depth-dependent G calculation."""
    total_moment = 0
    for item in fault_object_list:
        depth = item.get_centroid_depth()
        idx = np.abs(depths - depth).argmin()
        G = mus[idx]
        A = item.get_fault_area()
        d = item.get_total_slip()
        total_moment += moment_calculations.moment_from_muad(G, A, d)
    return total_moment


def check_consistent_reference_frame(triangle_fault_list) -> int:
    """Returns True if all the triangles have the coordinate reference system."""
    return_code = 1  # default true
    reflon, reflat = triangle_fault_list[0].lon, triangle_fault_list[0].lat
    for item in triangle_fault_list:
        if item.lon != reflon or item.lat != reflat:
            return_code = 0
    if return_code == 0:
        print("Warning! Not all triangles have the same reference lon/lat")
    return return_code


def fault_triangles_to_coulomb_fault(fault_object_list, zerolon_system=None, zerolat_system=None):
    """
    Convert a list of triangular fault objects into a list of small Coulomb source objects for Elastic_stresses_py.
    Parameters zerolon_system and zerolat_system can be passed in for system with 1+ faults.
    """
    source_object = []
    if len(fault_object_list) > 1 and zerolon_system is None:
        print("Warning! You are converting multiple faults to cartesian without defining a system coordinate system!")
    for onefault in fault_object_list:
        one_source = onefault.fault_triangle_to_coulomb_fault(zerolon_system, zerolat_system)
        source_object.append(one_source)
    return source_object


def extract_mesh_vertices(triangle_fault_list):
    """
    From a list of triangle objects, determine the set of vertices and indices that form the mesh.

    :param triangle_fault_list: list of triangle objects
    :return: mesh_points, index_list numpy arrays
    """
    index = 0
    lookup = {}
    mesh_points, index_list = [], []
    for tri in triangle_fault_list:
        for vertex in [tri.vertex1, tri.vertex2, tri.vertex3]:
            if tuple(vertex) not in lookup.keys():
                lookup[tuple(vertex)] = index
                mesh_points.append(vertex)
                index += 1

    for tri in triangle_fault_list:  # create the list of indices
        id1 = lookup[tuple(tri.vertex1)]
        id2 = lookup[tuple(tri.vertex2)]
        id3 = lookup[tuple(tri.vertex3)]
        index_list.append([id1, id2, id3])

    return np.array(mesh_points), np.array(index_list)


def flip_depth_sign(mesh_points):
    """
    Flip the z-coordinate of a numpy array mesh
    """
    flipped_mesh = []
    for item in mesh_points:
        flipped_mesh.append([item[0], item[1], -item[2]])
    return np.array(flipped_mesh)
