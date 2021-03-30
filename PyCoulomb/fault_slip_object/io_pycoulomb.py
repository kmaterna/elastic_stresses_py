""""
Functions for conversion of faults and slip distributions from PyCoulomb into list of fault_slip_object dictionaries

Some of these faults may come in through .intxt files. Use io_intxt to read them.

Format intxt: strike rake dip length(km) width(km) updip_corner_lon updip_corner_lat updip_corner_dep(km) slip
"""
from Elastic_stresses_py.PyCoulomb import coulomb_collections as cc
from Tectonic_Utils.geodesy import fault_vector_functions


def coulomb_fault_to_fault_dict(source_object):
    """Convert a list of fault objects from Elastic_stresses_py into a list of internal dictionary objects"""
    fault_dict_list = [];
    for src in source_object:
        if src.potency:
            print("ERROR! Cannot convert a point source into a rectangular source. Skipping...");
            continue;
        one_fault = {"strike": src.strike, "dip": src.dipangle, "depth": src.top,
                     "rake": fault_vector_functions.get_rake(src.rtlat, src.reverse),
                     "slip": fault_vector_functions.get_total_slip(src.rtlat, src.reverse),
                     "length": fault_vector_functions.get_strike_length(src.xstart, src.xfinish, src.ystart,
                                                                        src.yfinish),
                     "width": fault_vector_functions.get_downdip_width(src.top, src.bottom, src.dipangle)};
        lon, lat = fault_vector_functions.xy2lonlat(src.xstart, src.ystart, src.zerolon, src.zerolat);
        one_fault["lon"] = lon;
        one_fault["lat"] = lat;
        fault_dict_list.append(one_fault);
    return fault_dict_list;


def fault_dict_to_coulomb_fault(fault_dict_list):
    """Convert a list of internal dictionary objects into a list of source objects for Elastic_stresses_py"""
    source_object = [];
    for onefault in fault_dict_list:
        _top, bottom = fault_vector_functions.get_top_bottom_from_top(onefault['depth'], onefault['width'],
                                                                      onefault['dip']);
        rtlat, reverse = fault_vector_functions.get_rtlat_dip_slip(onefault['slip'], onefault['rake']);
        xfinish, yfinish = fault_vector_functions.add_vector_to_point(0, 0, onefault['length'], onefault['strike']);
        one_source = cc.Faults_object(xstart=0, xfinish=xfinish, ystart=0, yfinish=yfinish, Kode=100,
                                      rtlat=rtlat, reverse=reverse, tensile=0,
                                      potency=[], strike=onefault['strike'],
                                      dipangle=onefault['dip'], zerolon=onefault['lon'], zerolat=onefault['lat'],
                                      rake=onefault['rake'], top=onefault['depth'], bottom=bottom, comment='');
        source_object.append(one_source);
    return source_object;
