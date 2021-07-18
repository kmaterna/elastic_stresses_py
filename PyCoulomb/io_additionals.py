"""
Reading aftershock tables and GPS lon/lat pairs
"""

from . import coulomb_collections as cc
import numpy as np


def read_aftershock_table(infile):
    """Simple catalog format: time, lon, lat, depth, magnitude"""
    print("Reading aftershocks from file %s " % infile);
    lon, lat, time, depth, magnitude = [], [], [], [], [];

    ifile = open(infile);
    for line in ifile:
        temp = line.split();
        if temp[0][0] == '#':
            continue;
        else:
            time.append(temp[0]);
            lon.append(float(temp[1]));
            lat.append(float(temp[2]));
            depth.append(float(temp[3]));
            magnitude.append(float(temp[4]));
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
    disp_points_list = [];
    ifile = open(infile, 'r');
    for line in ifile:
        temp = line.split();
        if temp[0][0] == '#':
            continue;
        else:
            lon, lat = float(temp[0]), float(temp[1]);
            dE_obs, dN_obs, dU_obs = np.nan, np.nan, np.nan;
            Se_obs, Sn_obs, Su_obs = np.nan, np.nan, np.nan;
            name = "";
            if len(temp) >= 8:  # if we have longer GPS format with uncertainties
                name = temp[-1];
                dE_obs, dN_obs, dU_obs = float(temp[2]), float(temp[3]), float(temp[4]);
                Se_obs, Sn_obs, Su_obs = float(temp[5]), float(temp[6]), float(temp[7]);
            new_disp_point = cc.Displacement_points(lon=lon, lat=lat, dE_obs=dE_obs, dN_obs=dN_obs, dU_obs=dU_obs,
                                                    Se_obs=Se_obs, Sn_obs=Sn_obs, Su_obs=Su_obs, name=name);
            disp_points_list.append(new_disp_point);
    ifile.close();
    return disp_points_list;


def write_disp_points_results(disp_points, outfile):
    """
    Write the contents of disp_points (dE_obs etc.)
    """
    if len(disp_points) > 0:
        ofile = open(outfile, 'w');
        ofile.write("# Format: lon lat u v w (m)\n");
        for point in disp_points:
            ofile.write("%f %f %f %f %f\n" % (point.lon, point.lat, point.dE_obs, point.dN_obs, point.dU_obs));
        ofile.close();
    return;

def write_strain_results(obs_strain_points, strains, outfile):
    """
    obs_strain_points is an object of format cc.dips_points
    strains is a list of tensors
    """
    if obs_strain_points:
        ofile = open(outfile, 'w');
        ofile.write("# Format: lon lat strain_tensor (microstrain)\n")
        for i in range(len(obs_strain_points.lon)):
            eij = np.multiply(strains[i], 1e6);  # microstrain
            ofile.write("%f %f\n" % (obs_strain_points.lon[i], obs_strain_points.lat[i]));
            ofile.write("%f %f %f\n" % (eij[0][0], eij[0][1], eij[0][2]));
            ofile.write("%f %f %f\n" % (eij[1][0], eij[1][1], eij[1][2]));
            ofile.write("%f %f %f\n" % (eij[2][0], eij[2][1], eij[2][2]));
            ofile.write("\n");
        ofile.close();
    return;
