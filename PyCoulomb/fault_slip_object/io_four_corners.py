""""
Function for conversion of faults and slip distributions
from four corners into list of fault_slip_object dictionaries
"""

import numpy as np
from Tectonic_Utils.geodesy import fault_vector_functions


def read_four_corners_fault_file(filename):
    """
    Read fault file from Shengji Wei, EPSL, 2015.
    Provided: lat/lon/depth of each corner.
    Read into an internal fault dictionary
    """
    print("Reading file %s" % filename);
    lons, lats, depths = [], [], [];
    ifile = open(filename, 'r');
    for line in ifile:
        temp = line.split();
        if len(temp) == 0:
            continue;
        if temp[0][0] == "#":
            continue;
        else:
            lons.append(float(temp[0]))
            lats.append(float(temp[1]))
            depths.append(float(temp[2]))
    ifile.close();

    # Get parameters of fault patch; assuming four vertices are given, and top depth is same between two top vertices.
    x, y = fault_vector_functions.latlon2xy(lons, lats, lons[0], lats[0]);

    # Find the top of the fault plane
    shallow_depth = np.min(depths);
    deep_depth = np.max(depths);
    shallow_depth_idx = np.where(depths == shallow_depth)
    idx0 = shallow_depth_idx[0][0];
    idx1 = shallow_depth_idx[0][1];

    deltax = x[idx1] - x[idx0];
    deltay = y[idx1] - y[idx0];
    strike_initial = fault_vector_functions.get_strike(deltax, deltay)  # STRIKE MIGHT BE 180 DEGREES AWAY
    dip = fault_vector_functions.get_dip_degrees(x[-2], y[-2], depths[-2], x[-1], y[-1], depths[-1]);
    # Assuming the last two entries go from the bottom to the top.

    # MIGHT REVERSE THE STRIKE HERE.
    strike_vector = fault_vector_functions.get_strike_vector(strike_initial);
    dip_vector = fault_vector_functions.get_dip_vector(strike_initial, dip);
    xp = np.cross(dip_vector, strike_vector);  # dip x strike for outward facing normal, by right hand rule.
    if xp[2] < 0:
        print("WARNING: STRIKE AND DIP DON'T OBEY RIGHT HAND RULE. FLIPPING STRIKE.");
        strike = strike_initial+180;
    else:
        strike = strike_initial;
    if strike > 360:
        strike = strike-360;

    length = fault_vector_functions.get_strike_length(x[idx0], x[idx1], y[idx0], y[idx1]);
    width = fault_vector_functions.get_downdip_width(shallow_depth, deep_depth, dip);

    # GET THE TOP CORNER OF THE FAULT PLANE (ONE OF TWO OPTIONS)
    print("  ", strike, "degrees strike");
    print("  ", dip, "degrees dip");
    print("  ", length, "km length");
    print("  ", width, "km width");

    if strike_initial == strike:
        updip_corner_lon = lons[idx0];
        updip_corner_lat = lats[idx0];
        updip_corner_depth = depths[idx0];
    else:
        updip_corner_lon = lons[idx1];
        updip_corner_lat = lats[idx1];
        updip_corner_depth = depths[idx1];

    fault_object = {'strike': strike, 'slip': 0, 'tensile': 0, 'rake': 0, 'length': length, 'width': width, 'dip': dip,
                    'lon': updip_corner_lon, 'lat': updip_corner_lat, 'depth': updip_corner_depth, 'segment': 0};
    return fault_object;
