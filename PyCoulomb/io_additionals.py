"""
Reading aftershock tables and GPS lon/lat pairs
"""

from . import coulomb_collections as cc


def read_aftershock_table(infile):
    print("Reading aftershocks from file %s " % infile);
    lon, lat, time, depth, magnitude = [], [], [], [], [];

    ifile = open(infile);
    for line in ifile:
        temp = line.split();
        if temp[0][0] == '#':
            continue;
        else:
            time.append(temp[0]);
            lon.append(float(temp[3]));
            lat.append(float(temp[2]));
            depth.append(float(temp[4]));
            magnitude.append(float(temp[5]));
    ifile.close();
    return [lon, lat, depth, magnitude, time];


def read_disp_points(infile):
    """
    A file with lon/lat points that we are computing displacements.
    If the observed displacements are given in the additional columns,
    then we add them to the object so we can plot them against the model later.
    A slightly flexible-format read for:
    # lon lat
    or
    # lon lat de dn du se sn su name
    """
    print("Reading displacement points from file %s " % infile);
    lon, lat, names = [], [], [];
    dE_obs, dN_obs, dU_obs = [], [], [];
    Se_obs, Sn_obs, Su_obs = [], [], [];
    ifile = open(infile, 'r');
    for line in ifile:
        temp = line.split();
        if temp[0][0] == '#':
            continue;
        else:
            lon.append(float(temp[0]));
            lat.append(float(temp[1]));
            if len(temp) >= 8:
                names.append(temp[-1]);
                dE_obs.append(float(temp[2]));
                dN_obs.append(float(temp[3]));
                dU_obs.append(float(temp[4]));
                Se_obs.append(float(temp[5]));
                Sn_obs.append(float(temp[6]));
                Su_obs.append(float(temp[7]));
            else:
                names.append("");
    ifile.close();
    disp_points = cc.Displacement_points(lon=lon, lat=lat, dE_obs=dE_obs, dN_obs=dN_obs, dU_obs=dU_obs, Se_obs=Se_obs,
                                         Sn_obs=Sn_obs, Su_obs=Su_obs, name=names);
    return disp_points;
