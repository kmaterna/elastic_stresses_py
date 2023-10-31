"""
The functions in this package operate on an internal class for faults/slip distributions and read some common formats.
For all other formats, make sure you build read/write conversion AND TEST functions into internal format.
"""

from .. import coulomb_collections as cc
from Tectonic_Utils.geodesy import fault_vector_functions, haversine
from Tectonic_Utils.seismo import moment_calculations
from .. import conversion_math
import numpy as np

class FaultSlipObject:
    """
    An object for a rectangular fault patch, potentially with slip on it.
    If the fault is a receiver fault, we put slip = 0.

    :param strike: fault strike, in degrees,
    :type strike: float
    :param dip: dip, in degrees, 0<dip<90
    :type dip: float
    :param length: km
    :type length: float
    :param width: km
    :type width: float
    :param lon: longitude of back top corner (see Aki and Richards convention)
    :type lon: float
    :param lat: latitude of back top corner (see Aki and Richards convention)
    :type lat: float
    :param depth: depth of top of fault, in km, positive downwards
    :type depth: float
    :param rake: degrees.  Should be zero if slip is zero.
    :type rake: float
    :param slip: amount of slip, in meters
    :type slip: float
    :param tensile: amount of tensile opening, in meters
    :type tensile: float
    :param segment: segment number for this fault
    :type segment: int
    """

    def __init__(self, strike, dip, length, width, lon, lat, depth, rake, slip, tensile=0, segment=0):
        self.strike = strike;  # float
        self.dip = dip;  # float
        self.length = length;  # float
        self.width = width;  # float
        self.lon = lon;  # float
        self.lat = lat;  # float
        self.depth = depth;  # float
        self.rake = rake;  # float
        self.slip = slip;  # float
        self.tensile = tensile;  # float
        self.segment = segment;  # int
        if self.dip > 90:
            raise ValueError("Error! Dip greater than 90 degrees.");
        if self.dip < 0:
            raise ValueError("Error! Dip less than 0 degrees.");
        if self.length < 0:
            raise ValueError("Error! Length less than 0 km.");
        if self.width < 0:
            raise ValueError("Error! Width less than 0 km.");

    def __iter__(self):
        for attr, value in self.__dict__.items():
            yield attr, value

    def get_total_slip(self) -> float:
        """Helper function to return the total slip amount of a fault object (always > 0)"""
        return self.slip;

    def get_total_slip_mm(self) -> float:
        """Helper function to return the total slip amount of a fault object in mm (always > 0)"""
        return self.slip*1000;

    def get_rtlat_slip(self) -> float:
        """Helper function to return the right lateral slip amount of a fault object"""
        [rtlat, _] = fault_vector_functions.get_rtlat_dip_slip(self.slip, self.rake);
        return rtlat;

    def get_dip_slip(self) -> float:
        """Helper function to return the dip slip amount of a fault object"""
        [_, dipslip] = fault_vector_functions.get_rtlat_dip_slip(self.slip, self.rake);
        return dipslip;

    def get_fault_depth(self) -> float:
        """Helper function to return the depth of a fault object"""
        return self.depth;

    def get_fault_rake(self) -> float:
        """Helper function to return the depth of a fault object"""
        return self.rake;

    def get_segment(self) -> int:
        """Helper function to return the segment_num of a fault object"""
        return self.segment;

    @staticmethod
    def get_blank_fault():
        """Helper function to an empty color palette"""
        return "";

    def get_four_corners_lon_lat(self):
        """
        Return the lon/lat of all 4 corners of a fault_object.
        """
        [source] = fault_object_to_coulomb_fault([self]);
        [x_total, y_total, _, _] = conversion_math.get_fault_four_corners(source);
        lons, lats = fault_vector_functions.xy2lonlat(x_total, y_total, source.zerolon, source.zerolat);
        return lons, lats;

    def get_geographic_center(self):
        """
        Return the geographic center of a rectangular fault patch, in (lon, lat, depth).
        """
        lons, lats = self.get_four_corners_lon_lat();
        center_lon = np.nanmean(lons[0:4]);
        center_lat = np.nanmean(lats[0:4]);
        top, bottom = self.get_top_bottom_depths();
        center_depth = np.nanmean([top, bottom]);
        return center_lon, center_lat, center_depth;

    def get_updip_corners_lon_lat(self):
        """
        Return the lon/lat of 2 shallow corners of a fault_object.
        """
        [source] = fault_object_to_coulomb_fault([self]);
        [_, _, x_updip, y_updip] = conversion_math.get_fault_four_corners(source);
        lons, lats = fault_vector_functions.xy2lonlat(x_updip, y_updip, source.zerolon, source.zerolat);
        return lons, lats;

    def get_fault_element_distance(self, fault_object2, threedimensional=True) -> float:
        """Return the distance between two fault elements, in 3d, in km."""
        h_distance = haversine.distance([self.lat, self.lon], [fault_object2.lat, fault_object2.lon]);
        if threedimensional:
            depth_distance = self.depth - fault_object2.depth;
        else:
            depth_distance = 0;
        distance_3d = np.sqrt(np.square(h_distance) + np.square(depth_distance))
        return distance_3d;

    def get_top_bottom_depths(self):
        top, bottom = fault_vector_functions.get_top_bottom_from_top(self.depth, self.width, self.dip);
        return top, bottom;

    def get_fault_area(self) -> float:
        """Return the area of a rectangular fault. Units: square meters."""
        A = self.width * 1000 * self.length * 1000;
        return A;

    def get_fault_moment(self, mu=30e9) -> float:
        """Get the moment for a fault patch. Units: Newton-meters."""
        A = self.get_fault_area();
        d = self.slip;
        moment = moment_calculations.moment_from_muad(mu, A, d);
        return moment;

    def change_fault_slip(self, new_slip=None, new_rake=None, new_tensile=None):
        """
        Set the fault slip to something different.
        Can optionally also set the rake; otherwise, leave rake unchanged.

        :param new_slip: float, in meters
        :param new_rake: float, in degrees
        :param new_tensile: float, in meters
        :returns: fault_slip_object
        """
        desired_slip = self.slip if new_slip is None else new_slip;
        desired_rake = self.rake if new_rake is None else new_rake;
        desired_tensile = self.tensile if new_tensile is None else new_tensile;
        new_obj = FaultSlipObject(strike=self.strike, dip=self.dip, length=self.length, width=self.width, lon=self.lon,
                                  lat=self.lat, depth=self.depth, tensile=desired_tensile,
                                  slip=desired_slip, rake=desired_rake, segment=self.segment);
        return new_obj;

    # ------------ PREDICATES -------------- #
    def is_within_depth_range(self, upper_depth, lower_depth) -> bool:
        """
        Filter fault_object if it falls within depth range [upper_depth, lower_depth].

        :param upper_depth: float, km
        :param lower_depth: float, km
        :returns: bool
        """
        if upper_depth <= self.depth <= lower_depth:
            return True;
        else:
            return False;

    def is_within_bbox(self, bbox) -> bool:
        """
        Filter fault_object if it falls within bounding box.

        :param bbox: list of four floats, [w, e, s, n]
        :returns: bool
        """
        if bbox[0] <= self.lon <= bbox[1] and bbox[2] <= self.lat <= bbox[3]:
            return True;
        else:
            return False;

    # ------------ REGULAR FUNCTIONS -------------- #
    def construct_intxt_source_patch(self):
        """
        Turn a fault_slip_object into a string corresponding to a slip_patch in .intxt file.
        """
        fault_patch_string = "Source_Patch: " + str(self.strike) + " " + str(self.rake) + " ";
        fault_patch_string = fault_patch_string + str(self.dip) + " " + str(self.length) + " ";
        fault_patch_string = fault_patch_string + str(self.width) + " " + str(self.lon) + " ";
        fault_patch_string = fault_patch_string + str(self.lat) + " " + str(self.depth) + " ";
        fault_patch_string = fault_patch_string + str(self.slip) + " " + str(self.tensile) + " ";
        return fault_patch_string;

    def add_slip_from_two_faults(self, fault_object2):
        """Assuming identical fault_slip_object except for different slip amounts. Add slip amounts. """
        ss_1, ds_1 = fault_vector_functions.get_rtlat_dip_slip(self.slip, self.rake);
        ss_2, ds_2 = fault_vector_functions.get_rtlat_dip_slip(fault_object2.slip, fault_object2.rake);
        ss_total = ss_1 + ss_2;  # rtlat
        ds_total = ds_1 + ds_2;  # reverse
        slip_total = fault_vector_functions.get_total_slip(ss_total, ds_total);
        rake_total = fault_vector_functions.get_rake(rtlat_strike_slip=ss_total, dip_slip=ds_total);
        new_item = FaultSlipObject(strike=self.strike, dip=self.dip, length=self.length, width=self.width,
                                   lon=self.lon, lat=self.lat, depth=self.depth, slip=slip_total, rake=rake_total,
                                   tensile=self.tensile + fault_object2.tensile,
                                   segment=self.segment)
        return new_item;

    def fault_object_to_coulomb_fault(self, zerolon_system=None, zerolat_system=None):
        """
        Convert an internal fault object into a source object for Elastic_stresses_py.
        By default, the bottom corner of the fault is the center of the coordinate system, but
        Parameters zerolon_system and zerolat_system can be passed in.
        """
        zerolon = self.lon if not zerolon_system else zerolon_system;
        zerolat = self.lat if not zerolat_system else zerolat_system;
        _, bottom = self.get_top_bottom_depths()
        [startx, starty] = fault_vector_functions.latlon2xy_single(self.lon, self.lat, zerolon, zerolat);

        rtlat, reverse = fault_vector_functions.get_rtlat_dip_slip(self.slip, self.rake);
        xfinish, yfinish = fault_vector_functions.add_vector_to_point(startx, starty, self.length, self.strike);
        one_source = cc.construct_pycoulomb_fault(xstart=startx, xfinish=xfinish, ystart=starty, yfinish=yfinish,
                                                  Kode=100, rtlat=rtlat, reverse=reverse, tensile=self.tensile,
                                                  potency=[], strike=self.strike, dipangle=self.dip,
                                                  zerolon=zerolon, zerolat=zerolat, rake=self.rake,
                                                  top=self.depth, bottom=bottom, comment='');
        return one_source;

    def subdivide_by_depth(self, num_divisions):
        """
        Subdivide a rectangular fault into an integer number of sub-faults along its width (i.e., down-dip) axis.
        """
        divided_faults_list = [];
        lons, lats = self.get_four_corners_lon_lat();
        lons_range = lons[-2] - lons[-1];
        lats_range = lats[-2] - lats[-1];
        top, bottom = self.get_top_bottom_depths();
        depth_range = bottom - top;
        for i in range(num_divisions):
            new_item = FaultSlipObject(strike=self.strike, dip=self.dip, length=self.length,
                                       width=self.width/num_divisions,
                                       lon=self.lon + lons_range * i / num_divisions,
                                       lat=self.lat + lats_range * i / num_divisions,
                                       depth=self.depth + depth_range * i / num_divisions,
                                       slip=self.slip, rake=self.rake, tensile=self.tensile, segment=self.segment);
            divided_faults_list.append(new_item);
        return divided_faults_list;


"""
Functions that work on lists of fault items
"""


def get_four_corners_lon_lat_multiple(fault_object_list):
    """
    Return global bounding-box lon/lat for a list of fault_objects as (W, E, S, N)
    """
    lons_all, lats_all = [], [];
    for item in fault_object_list:
        lons, lats = item.get_four_corners_lon_lat();
        lons_all = lons_all + lons;
        lats_all = lats_all + lats;
    return np.min(lons_all), np.max(lons_all), np.min(lats_all), np.max(lats_all);


def get_total_moment(fault_object_list, mu=30e9) -> float:
    """
    Return the total moment of a list of slip objects, in fault_object
    Moment in newton-meters
    """
    total_moment = 0;
    for item in fault_object_list:
        total_moment += item.get_fault_moment(mu=mu);
    return total_moment;


def get_total_moment_depth_dependent(fault_object_list, depths, mus) -> float:
    """Compute total moment using a depth-dependent G calculation"""
    total_moment = 0;
    for item in fault_object_list:
        depth = item.get_fault_depth();
        idx = np.abs(depths - depth).argmin()
        G = mus[idx];
        A = item.get_fault_area();
        d = item.get_total_slip();
        total_moment += moment_calculations.moment_from_muad(G, A, d);
    return total_moment;


def add_two_fault_object_lists(list1, list2):
    """Assuming identical geometry in the two lists"""
    if len(list1) != len(list2):
        raise Exception("Error! Two fault_object lists are not identical");
    new_list = [];
    for item1, item2 in zip(list1, list2):
        new_item = item1.add_slip_from_two_faults(item2);
        new_list.append(new_item);
    return new_list;


def change_fault_slip_list(fault_object_list, new_slip=None, new_rake=None):
    """
    Set the fault slip on a list of fault_slip_objects to something different.
    Can optionally also set the rake on all fault patches to a constant value; otherwise, leave rake unchanged.

    :param fault_object_list: list of fault_slip_objects
    :param new_slip: float, in meters
    :param new_rake: float, in degrees
    :returns new_list: a list of fault_slip_objects
    """
    new_list = [];
    for item in fault_object_list:
        new_obj = item.change_fault_slip(new_slip=new_slip, new_rake=new_rake);
        new_list.append(new_obj);
    return new_list;


def filter_by_depth(fault_object_list, upper_depth, lower_depth):
    """
    Filter a list of fault_objects to only those that fall within the depth range [upper_depth, lower_depth].

    :param fault_object_list: list of fault_slip_objects
    :param upper_depth: float, km
    :param lower_depth: float, km
    :returns: a list of fault_slip_objects
    """
    return [x for x in fault_object_list if x.is_within_depth_range(upper_depth, lower_depth)];


def get_how_many_segments(fault_object_list):
    """
    Reduce fault_object_list into the number of segments (often 1) and the number of individual fault patches.

    :param fault_object_list: list of fault_slip_objects
    """
    segments = [x.get_segment() for x in fault_object_list];
    num_segments = len(set(segments));
    num_patches = len(segments);
    return num_segments, num_patches;


def filter_by_segment(fault_object_list, segment_num=0):
    """
    Filter a list of fault_objects to only those that have segment=x.

    :param fault_object_list: list of fault_slip_objects
    :param segment_num: int
    :returns: a list of fault_slip_objects
    """
    return [x for x in fault_object_list if x.get_segment() == segment_num];


def coulomb_fault_to_fault_object(source_object):
    """Convert a list of fault objects from Elastic_stresses_py into a list of internal fault objects."""
    if not source_object:
        return [];
    fault_object_list = [];
    for src in source_object:
        if src.potency:
            print("ERROR! Cannot convert a point source into a rectangular source. Skipping...");
            continue;
        lon, lat = fault_vector_functions.xy2lonlat(src.xstart, src.ystart, src.zerolon, src.zerolat);
        slip = fault_vector_functions.get_total_slip(src.rtlat, src.reverse);
        length = fault_vector_functions.get_strike_length(src.xstart, src.xfinish, src.ystart, src.yfinish);
        width = fault_vector_functions.get_downdip_width(src.top, src.bottom, src.dipangle);
        one_fault = FaultSlipObject(strike=src.strike, dip=src.dipangle, depth=src.top, rake=src.rake,
                                    slip=slip, length=length, width=width, tensile=src.tensile, lon=lon, lat=lat);
        fault_object_list.append(one_fault);
    return fault_object_list;


def fault_object_to_coulomb_fault(fault_object_list, zerolon_system=None, zerolat_system=None):
    """
    Convert a list of internal fault objects into a list of source objects for Elastic_stresses_py.
    By default, the bottom corner of the fault is the center of the coordinate system, but
    Parameters zerolon_system and zerolat_system can be passed in for system with 1+ faults.
    """
    source_object = [];
    if len(fault_object_list) > 1 and zerolon_system is None:
        print("Warning! You are converting multiple faults to cartesian without defining a system coordinate system!");
    for onefault in fault_object_list:
        one_source = onefault.fault_object_to_coulomb_fault(zerolon_system, zerolat_system);
        source_object.append(one_source);
    return source_object;
