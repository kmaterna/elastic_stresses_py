import numpy as np
from PyCoulomb import coulomb_collections
from PyCoulomb import conversion_math
from PyCoulomb import io_intxt


def read_intxt(input_file):
    print("Reading source and receiver fault information from file %s " % input_file);
    sources = [];
    receivers = [];
    ifile = open(input_file, 'r');
    for line in ifile:
        temp = line.split();
        if len(temp) == 0:
            continue;
        if temp[0] == 'S:':
            [strike, rake, dip, lon, lat, depth, magnitude, mu, lame1] = read_point_source_line(line);
            [x, y, Kode, rtlat, reverse, potency, comment] = compute_params_for_point_source(strike, dip, rake,
                                                                                             magnitude, lon, lat,
                                                                                             zerolon, zerolat, mu,
                                                                                             lame1);
            one_source_object = coulomb_collections.Faults_object(xstart=x, xfinish=x, ystart=y, yfinish=y, Kode=Kode,
                                                                  rtlat=rtlat, reverse=reverse, potency=potency,
                                                                  strike=strike, dipangle=dip, rake=rake, top=depth,
                                                                  bottom=depth, comment=comment);
            sources.append(one_source_object);
        elif temp[0] == 'R:':
            [strike, rake, dip, L, W, fault_lon, fault_lat, fault_depth] = io_intxt.read_receiver_line(line);
            [xstart, xfinish, ystart, yfinish, Kode, rtlat, reverse, top, bottom,
             comment] = io_intxt.compute_params_for_slip_source(strike, dip, rake, fault_depth, L, W, fault_lon,
                                                                fault_lat, 0, zerolon, zerolat);
            one_receiver_object = coulomb_collections.Faults_object(xstart=xstart, xfinish=xfinish, ystart=ystart,
                                                                    yfinish=yfinish, Kode=Kode, rtlat=rtlat,
                                                                    reverse=reverse, potency=[], strike=strike,
                                                                    dipangle=dip, rake=rake, top=top, bottom=bottom,
                                                                    comment=comment);
            receivers.append(one_receiver_object);
        elif temp[0] == 'G:':
            [PR1, FRIC, minlon, maxlon, zerolon, minlat, maxlat, zerolat] = io_intxt.read_general_line(line);
        else:
            continue;
    ifile.close();

    # Wrapping up the inputs
    [start_gridx, finish_gridx, start_gridy, finish_gridy, xinc, yinc] = io_intxt.compute_grid_parameters(minlon,
                                                                                                          maxlon,
                                                                                                          zerolon,
                                                                                                          minlat,
                                                                                                          maxlat,
                                                                                                          zerolat);
    input_obj = coulomb_collections.Input_object(PR1=PR1, FRIC=FRIC, depth=0, start_gridx=start_gridx,
                                                 finish_gridx=finish_gridx, start_gridy=start_gridy,
                                                 finish_gridy=finish_gridy,
                                                 xinc=xinc, yinc=yinc, minlon=minlon, maxlon=maxlon, zerolon=zerolon,
                                                 minlat=minlat, maxlat=maxlat, zerolat=zerolat,
                                                 receiver_object=receivers, source_object=sources);
    return input_obj;


def read_point_source_line(line):
    # Format: strike rake dip lon lat depth magnitude mu lamdba
    strike = float(line.split()[1]);
    rake = float(line.split()[2]);
    dip = float(line.split()[3]);
    lon = float(line.split()[4]);
    lat = float(line.split()[5]);
    depth = float(line.split()[6]);
    magnitude = float(line.split()[7]);
    mu = float(line.split()[8]);
    lame1 = float(line.split()[9]);
    return [strike, rake, dip, lon, lat, depth, magnitude, mu, lame1];


def compute_params_for_point_source(strike, dipangle, rake, magnitude, lon, lat, zerolon, zerolat, mu, lame1):
    # Given information about point sources from focal mechanisms,
    # Return the right components that get packaged into input_obj.
    [xcenter, ycenter] = conversion_math.latlon2xy(lon, lat, zerolon, zerolat);
    potency = get_DC_potency(rake, magnitude, mu, lame1);
    # Filler variables for the point source case
    Kode = 100;
    comment = '';
    return [xcenter, ycenter, Kode, 0, 0, potency, comment];


def get_DC_potency(rake, momentmagnitude, mu, lame1):
    # Given the basic double couple parameters,
    # Return the four-vector used in Okada DC3D0.
    # Pot1 = strike-slip moment of DC / mu
    # Pot2 = dip-slip moment of DC / mu
    # Pot3 = inflation = M_ISO / lambda
    # Pot4 = tensile = M_lineardipole / mu
    # In a more general case, we would use a different MT format to handle non-DC parts.
    # Right now, it only handles DC focal mechanisms.
    total_moment = moment(momentmagnitude);
    strike_slip_fraction, dip_slip_fraction = conversion_math.get_rtlat_dip_slip(1.0, rake);
    print("strike_slip: ", strike_slip_fraction);
    print("dip_slip: ", dip_slip_fraction);
    strike_slip_fraction = -1 * strike_slip_fraction;  # DC3D0 wants left lateral slip.
    p1 = total_moment * strike_slip_fraction / mu;
    p2 = total_moment * dip_slip_fraction / mu;
    # In the double-couple case, this is zero.
    p3 = 0;
    p4 = 0;
    return [p1, p2, p3, p4];


def moment(Mw):
    # Current definition of moment magnitude, returning moment in newton meters
    exponent = 1.5 * Mw + 1.5 * 10.7;
    moment = np.power(10, exponent);
    moment_newton_meters = moment * 1e-7;
    return moment_newton_meters;
