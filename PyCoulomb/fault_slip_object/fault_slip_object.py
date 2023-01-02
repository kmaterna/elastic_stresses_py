"""
The functions in this package operate on an internal format of faults/ slip distributions
For all other formats, make sure you build read/write conversion AND TEST functions into internal format.

    * Format(File): geojson format for Slippy input faults
    * Format(File): slippy format, output format for Elastic_stresses_py and Slippy
    * Format(File): static1d/visco1d, Fred Pollitz's STATIC1D code.
    * Format(File): four-corners, used for Wei et al. 2015
    * Format(File): .fsp format, used by USGS fault distributions and srcmod database
    * Format(Memory): fault_object with elements for a slip distribution (INTERNAL FORMAT FOR THIS LIBRARY)
    * Format(Memory): PyCoulomb format, a named tuple (For Elastic_stresses_py.PyCoulomb)
"""

from .. import coulomb_collections as cc
from Tectonic_Utils.geodesy import fault_vector_functions, haversine
from Tectonic_Utils.seismo import moment_calculations
from .. import conversion_math
import numpy as np
import collections.abc

class FaultSlipObject:
    """
    The internal format is a class:
    {
        strike(deg),
        dip(deg),
        length(km),
        width(km),
        lon(back top corner),
        lat(back top corner),
        depth(top, km), positive downwards
        rake(deg),
        slip(m),
        tensile(m)
        segment(int)
    }
    If the fault is a receiver fault, we put slip = 0
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

    def __iter__(self):
        for attr, value in self.__dict__.items():
            yield attr, value

    def get_total_slip(self):
        """Helper function to return the total slip amount of a fault object (always > 0)"""
        return self.slip;

    def get_total_slip_mm(self):
        """Helper function to return the total slip amount of a fault object in mm (always > 0)"""
        return self.slip*1000;

    def get_rtlat_slip(self):
        """Helper function to return the right lateral slip amount of a fault object"""
        [rtlat, _] = fault_vector_functions.get_rtlat_dip_slip(self.slip, self.rake);
        return rtlat;

    def get_dip_slip(self):
        """Helper function to return the dip slip amount of a fault object"""
        [_, dipslip] = fault_vector_functions.get_rtlat_dip_slip(self.slip, self.rake);
        return dipslip;

    def get_fault_depth(self):
        """Helper function to return the depth of a fault object"""
        return self.depth;

    def get_fault_rake(self):
        """Helper function to return the depth of a fault object"""
        return self.rake;

    def get_segment(self):
        """Helper function to return the segment_num of a fault object"""
        return self.segment;

    @staticmethod
    def get_blank_fault():
        """Helper function to an empty color palette"""
        return "";

    def get_four_corners_lon_lat(self):
        """
        Return the lon/lat of all 4 corners of a fault_object
        """
        [source] = fault_object_to_coulomb_fault([self]);
        [x_total, y_total, _, _] = conversion_math.get_fault_four_corners(source);
        lons, lats = fault_vector_functions.xy2lonlat(x_total, y_total, source.zerolon, source.zerolat);
        return lons, lats;

    def get_updip_corners_lon_lat(self):
        """
        Return the lon/lat of 2 shallow corners of a fault_object
        """
        [source] = fault_object_to_coulomb_fault([self]);
        [_, _, x_updip, y_updip] = conversion_math.get_fault_four_corners(source);
        lons, lats = fault_vector_functions.xy2lonlat(x_updip, y_updip, source.zerolon, source.zerolat);
        return lons, lats;

    def get_fault_element_distance(self, fault_dict2, threedimensional=True):
        """Return the distance between two fault elements, in 3d, in km"""
        h_distance = haversine.distance([self.lat, self.lon], [fault_dict2.lat, fault_dict2.lon]);
        if threedimensional:
            depth_distance = self.depth - fault_dict2.depth;
        else:
            depth_distance = 0;
        distance_3d = np.sqrt(np.square(h_distance) + np.square(depth_distance))
        return distance_3d;

    def get_fault_area(self):
        """Return the area of a rectangular fault. Units: square meters"""
        A = self.width * 1000 * self.length * 1000;
        return A;

    def get_fault_moment(self, mu=30e9):
        """Get the moment for a fault patch. Units: Newton-meters"""
        A = self.get_fault_area();
        d = self.slip;
        moment = moment_calculations.moment_from_muad(mu, A, d);
        return moment;

    def change_fault_slip(self, new_slip=None, new_rake=None):
        """
        Set the fault slip to something different.
        Can optionally also set the rake; otherwise, leave rake unchanged.

        :param new_slip: float, in meters
        :param new_rake: float, in degrees
        :returns: fault_slip_object
        """
        desired_slip = self.slip if new_slip is None else new_slip;
        desired_rake = self.rake if new_rake is None else new_rake;
        new_obj = FaultSlipObject(strike=self.strike, dip=self.dip, length=self.length, width=self.width, lon=self.lon,
                                  lat=self.lat, depth=self.depth, tensile=self.tensile,
                                  slip=desired_slip, rake=desired_rake, segment=self.segment);
        return new_obj;

    def is_within_depth_range(self, upper_depth, lower_depth):
        """
        Filter fault_object if it falls within depth range [upper_depth, lower_depth]

        :param upper_depth: float, km
        :param lower_depth: float, km
        :returns: bool
        """
        if upper_depth <= self.depth <= lower_depth:
            return 1;
        else:
            return 0;

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

    def add_slip_from_two_faults(self, fault_dict2):
        """Assuming identical fault_slip_object except for different slip amounts. Add slip amounts. """
        ss_1, ds_1 = fault_vector_functions.get_rtlat_dip_slip(self.slip, self.rake);
        ss_2, ds_2 = fault_vector_functions.get_rtlat_dip_slip(fault_dict2.slip, fault_dict2.rake);
        ss_total = ss_1 + ss_2;  # rtlat
        ds_total = ds_1 + ds_2;  # reverse
        slip_total = fault_vector_functions.get_total_slip(ss_total, ds_total);
        rake_total = fault_vector_functions.get_rake(rtlat_strike_slip=ss_total, dip_slip=ds_total);
        new_item = FaultSlipObject(strike=self.strike, dip=self.dip, length=self.length, width=self.width,
                                   lon=self.lon, lat=self.lat, depth=self.depth, slip=slip_total, rake=rake_total,
                                   tensile=self.tensile + fault_dict2.tensile,
                                   segment=self.segment)
        return new_item;


"""
Functions that work on lists of fault items
"""


def get_four_corners_lon_lat_multiple(fault_object_list):
    """
    Return the lon/lat of all 4 corners of a list of fault_objects
    Basically the bounding box for this list of fault_dict_objects
    """
    lons_all, lats_all = [], [];
    for item in fault_object_list:
        lons, lats = item.get_four_corners_lon_lat();
        lons_all = lons_all + lons;
        lats_all = lats_all + lats;
    return np.min(lons_all), np.max(lons_all), np.min(lats_all), np.max(lats_all);


def get_total_moment(fault_object_list, mu=30e9):
    """
    Return the total moment of a list of slip objects, in fault_dict_object
    Moment in newton-meters
    """
    total_moment = 0;
    for item in fault_object_list:
        total_moment += item.get_fault_moment(mu=mu);
    return total_moment;


def get_total_moment_depth_dependent(fault_object_list, depths, mus):
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
        raise Exception("Error! Two fault_dict lists are not identical");
    new_list = [];
    for item1, item2 in zip(list1, list2):
        new_item = item1.add_slip_from_two_faults(item2);
        new_list.append(new_item);
    return new_list;


def change_fault_slip_list(fault_object_list, new_slip=None, new_rake=None):
    """
    Set the fault slip on a list of fault_slip_dictionaries to something different.
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
    Filter a list of fault_objects to only those that fall within the depth range [upper_depth, lower_depth]

    :param fault_object_list: list of fault_slip_objects
    :param upper_depth: float, km
    :param lower_depth: float, km
    :returns new_list: a list of fault_slip_objects
    """
    new_list = [];
    for item in fault_object_list:
        if item.is_within_depth_range(upper_depth, lower_depth):
            new_list.append(item);
    return new_list;


def get_how_many_segments(fault_object_list):
    """
    Reduce fault_object_list into the number of segments (often 1) and the number of individual fault patches

    :param fault_object_list: list of fault_slip_objects
    """
    segments = [x.get_segment() for x in fault_object_list];
    num_segments = len(set(segments));
    num_patches = len(segments);
    return num_segments, num_patches;


def filter_by_segment(fault_object_list, segment_num=0):
    """
    Filter a list of fault_objects to only those that have segment=x

    :param fault_object_list: list of fault_slip_objects
    :param segment_num: int
    :returns new_list: a list of fault_slip_objects
    """
    new_list = [];
    for item in fault_object_list:
        if item.get_segment() == segment_num:
            new_list.append(item);
    return new_list;


def get_blank_fault_function(x):
    return x.get_blank_fault();


def write_gmt_fault_file(fault_object_list, outfile, color_mappable=get_blank_fault_function, verbose=True):
    """
    Write the 4 corners of a fault and its slip values into a multi-segment file for plotting in GMT
    By default, does not provide color on the fault patches
    color_mappable can be 1d array of scalars, or a function of the object's class.
    """
    if verbose:
        print("Writing file %s " % outfile);
    ofile = open(outfile, 'w');
    for i, fault in enumerate(fault_object_list):
        lons, lats = fault.get_four_corners_lon_lat();
        if isinstance(color_mappable, collections.abc.Sequence):
            color_string = "-Z"+str(color_mappable[i]);  # if separately providing the color array
        else:
            color_string = "-Z"+str(color_mappable(fault));  # call the function that you've provided
        ofile.write("> "+color_string+"\n");
        ofile.write("%f %f\n" % (lons[0], lats[0]));
        ofile.write("%f %f\n" % (lons[1], lats[1]));
        ofile.write("%f %f\n" % (lons[2], lats[2]));
        ofile.write("%f %f\n" % (lons[3], lats[3]));
        ofile.write("%f %f\n" % (lons[0], lats[0]));
    ofile.close();
    return;


def write_gmt_surface_trace(fault_object_list, outfile, verbose=True):
    """
    Write the 2 updip corners of a rectangular fault into a multi-segment file for plotting in GMT
    """
    if verbose:
        print("Writing file %s " % outfile);
    ofile = open(outfile, 'w');
    for fault in fault_object_list:
        lons, lats = fault.get_four_corners_lon_lat();
        ofile.write("> -Z\n");
        ofile.write("%f %f\n" % (lons[0], lats[0]));
        ofile.write("%f %f\n" % (lons[1], lats[1]));
    ofile.close();
    return;


def write_gmt_vertical_fault_file(fault_object_list, outfile, color_mappable=get_blank_fault_function):
    """
    Write the vertical coordinates of planar fault patches (length and depth, in local coords instead of lon/lat)
    and associated slip values into a multi-segment file for plotting in GMT.
    Good for vertical faults.  Plots with depth as a negative number.
    Works for only one planar fault segment.
    """
    print("Writing file %s " % outfile);

    # Get origin: extremal patch at the top. First, find bounding box for top points
    depth_array = [x.depth for x in fault_object_list];
    top_row_patches = [x for x in fault_object_list if x.depth == np.nanmin(depth_array)];
    top_row_lon, top_row_lat = [], [];
    for patch in top_row_patches:
        lon_updip, lat_updip = patch.get_updip_corners_lon_lat();
        top_row_lon = top_row_lon + lon_updip;
        top_row_lat = top_row_lat + lat_updip;  # joining two lists
    bbox = [np.nanmin(top_row_lon), np.nanmax(top_row_lon), np.nanmin(top_row_lat), np.nanmax(top_row_lat)];

    # Find fault corner coordinates that are candidates for extremal points on fault. Choose one for origin.
    origin_ll = [np.nan, np.nan];
    for lon, lat in zip(top_row_lon, top_row_lat):
        if lon in bbox and lat in bbox:
            origin_ll = [lon, lat];  # this should be guaranteed to happen twice, once for each end.
            break;

    ofile = open(outfile, 'w');
    for fault in fault_object_list:
        [source] = fault_object_to_coulomb_fault([fault], zerolon_system=origin_ll[0], zerolat_system=origin_ll[1]);
        [_, _, x_updip, y_updip] = conversion_math.get_fault_four_corners(source);
        deeper_offset = fault.width*np.sin(np.deg2rad(fault.dip));
        [xprime, _] = conversion_math.rotate_list_of_points(x_updip, y_updip, 90+fault.strike);
        start_x, finish_x = xprime[0], xprime[1];
        slip_amount = color_mappable(fault);
        ofile.write("> -Z"+str(-slip_amount)+"\n");  # currently writing left-lateral slip as positive
        ofile.write("%f %f\n" % (start_x, -fault.depth));
        ofile.write("%f %f\n" % (finish_x, -fault.depth));
        ofile.write("%f %f\n" % (finish_x, -fault.depth-deeper_offset));
        ofile.write("%f %f\n" % (start_x, -fault.depth-deeper_offset));
        ofile.write("%f %f\n" % (start_x, -fault.depth));

    ofile.close();
    return;


def coulomb_fault_to_fault_object(source_object):
    """Convert a list of fault objects from Elastic_stresses_py into a list of internal fault objects"""
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
                                    slip=slip, length=length, width=width, tensile=src.tensile, lon=lon, lat=lat,
                                    segment=0);
        fault_object_list.append(one_fault);
    return fault_object_list;


def fault_object_to_coulomb_fault(fault_object_list, zerolon_system=None, zerolat_system=None):
    """
    Convert a list of internal fault objects into a list of source objects for Elastic_stresses_py
    By default, the bottom corner of the fault is the center of the coordinate system, but
    Parameters zerolon_system and zerolat_system can be passed in for system with 1+ faults.
    """
    source_object = [];
    if len(fault_object_list) > 1 and zerolon_system is None:
        print("Warning! You are converting multiple faults to cartesian without defining a system coordinate system!");
    for onefault in fault_object_list:
        zerolon = onefault.lon if not zerolon_system else zerolon_system;
        zerolat = onefault.lat if not zerolat_system else zerolat_system;
        _top, bottom = fault_vector_functions.get_top_bottom_from_top(onefault.depth, onefault.width, onefault.dip);
        [startx, starty] = fault_vector_functions.latlon2xy_single(onefault.lon, onefault.lat, zerolon, zerolat);

        rtlat, reverse = fault_vector_functions.get_rtlat_dip_slip(onefault.slip, onefault.rake);
        xfinish, yfinish = fault_vector_functions.add_vector_to_point(startx, starty, onefault.length,
                                                                      onefault.strike);
        one_source = cc.construct_pycoulomb_fault(xstart=startx, xfinish=xfinish, ystart=starty, yfinish=yfinish,
                                                  Kode=100, rtlat=rtlat, reverse=reverse, tensile=onefault.tensile,
                                                  potency=[], strike=onefault.strike, dipangle=onefault.dip,
                                                  zerolon=zerolon, zerolat=zerolat, rake=onefault.rake,
                                                  top=onefault.depth, bottom=bottom, comment='');
        source_object.append(one_source);
    return source_object;
