"""
Reading aftershock tables and GPS lon/lat pairs
"""

from .disp_points_object.disp_points_object import Displacement_points
from . import utilities, conversion_math
import numpy as np
from tectonic_utils.seismo import moment_calculations
from . import pyc_fault_object as pycfaults


def read_aftershock_table(infile):
    """
    Simple catalog format: time, lon, lat, depth, magnitude

    :param infile: string, name of input file
    :returns: list of five lists, representing lon, lat, depth, magnitude (all floats), and time (str)
    """
    print("Reading aftershocks from file %s " % infile)
    lon, lat, time, depth, magnitude = [], [], [], [], []

    ifile = open(infile)
    for line in ifile:
        temp = line.split()
        if temp[0][0] == '#':
            continue
        else:
            time.append(temp[0])
            lon.append(float(temp[1]))
            lat.append(float(temp[2]))
            depth.append(float(temp[3]))
            magnitude.append(float(temp[4]))
    ifile.close()
    return [lon, lat, depth, magnitude, time]


def read_strain_points(infile):
    """
    A file with lon/lat points for which we are computing displacements.
    Format: "lon lat [depth]"

    :param infile: string, filename
    :returns: list of disp_point objects
    """
    print("Reading displacement points from file %s " % infile)
    disp_points_list = []
    with open(infile, 'r') as ifile:
        for line in ifile:
            temp = line.split()
            if temp[0][0] == '#':
                continue
            else:
                lon, lat = float(temp[0]), float(temp[1])
                depth = 0
                if len(temp) == 3:
                    depth = float(temp[2])
                new_disp_point = Displacement_points(lon=lon, lat=lat, depth=depth)
                disp_points_list.append(new_disp_point)
    print("--> Read %d strain points " % len(disp_points_list))
    return disp_points_list


def read_disp_points(infile):
    """
    A file with lon/lat points that we are computing displacements.
    If the observed displacements are given in the additional columns,
    then we add them to the object for later plotting against the model.
    A slightly flexible-format read for:
    - "lon lat"
    - "lon lat name"
    - "lon lat de dn du name"
    - "lon lat de dn du se sn su name"

    :param infile: input filename, string
    :returns: list of disp_point objects
    """
    print("Reading displacement points from file %s " % infile)
    disp_points_list = []
    with open(infile, 'r') as ifile:
        for line in ifile:
            temp = line.split()
            if temp[0][0] == '#':
                continue
            else:
                # Ultimately it might be better to have a different way of determining formats, but for now...
                lon, lat = float(temp[0]), float(temp[1])
                if len(temp) == 2:  # if file is: lon, lat
                    new_disp_point = Displacement_points(lon=lon, lat=lat)
                    disp_points_list.append(new_disp_point)
                elif len(temp) == 3:  # if file is: lon, lat, name
                    new_disp_point = Displacement_points(lon=lon, lat=lat, name=temp[2])
                    disp_points_list.append(new_disp_point)
                elif len(temp) == 5 or len(temp) == 6:  # if file is: lon, lat, de, dn, du [name] without uncertainties
                    dE_obs, dN_obs, dU_obs = float(temp[2]), float(temp[3]), float(temp[4])
                    new_disp_point = Displacement_points(lon=lon, lat=lat, dE_obs=dE_obs, dN_obs=dN_obs, dU_obs=dU_obs)
                    disp_points_list.append(new_disp_point)
                elif len(temp) >= 8:  # if we have longer GPS format with uncertainties
                    name = temp[-1]
                    dE_obs, dN_obs, dU_obs = float(temp[2]), float(temp[3]), float(temp[4])
                    Se_obs, Sn_obs, Su_obs = float(temp[5]), float(temp[6]), float(temp[7])
                    new_disp_point = Displacement_points(lon=lon, lat=lat, dE_obs=dE_obs, dN_obs=dN_obs, dU_obs=dU_obs,
                                                         Se_obs=Se_obs, Sn_obs=Sn_obs, Su_obs=Su_obs, name=name)
                    disp_points_list.append(new_disp_point)
    print("--> Read %d displacement points " % len(disp_points_list))
    return disp_points_list


def write_disp_points_locations(disp_points_list, filename, precision=6):
    """
    Write the simplest format for lon/lat points

    :param disp_points_list: list of disp points
    :param filename: string
    :param precision: number of significant digits to be printed for the point locations
    """
    print("Writing file %s " % filename)
    with open(filename, 'w') as ofile:
        for item in disp_points_list:
            format_string = "%."+str(precision)+"f %."+str(precision)+"f \n"
            ofile.write(format_string % (item.lon, item.lat))
    return


def write_fault_traces_gmt(fault_list, outfile):
    """
    Write the top edge of each fault in GMT map coordinates

    :param fault_list: pycoulomb source object
    :param outfile: string, filename for printing traces
    """
    fault_list, pt_sources, mogis = utilities.separate_source_types(fault_list)
    if not fault_list:
        return
    print("Writing %s" % outfile)
    with open(outfile, 'w') as ofile:
        for fault in fault_list:
            [_, _, x_updip, y_updip] = fault.get_fault_four_corners_geographic()
            ofile.write("> \n")
            ofile.write("%f %f\n" % (x_updip[0], y_updip[0]))
            ofile.write("%f %f\n" % (x_updip[1], y_updip[1]))
    return


def write_disp_points_results(disp_points, outfile):
    """
    Write the contents of disp_points (dE_obs etc.) into a file.  Mirrors the function read_disp_points().

    :param disp_points: list of disp_points objects
    :param outfile: filename, string
    """
    if not disp_points:
        return
    print("Writing %s " % outfile)
    with open(outfile, 'w') as ofile:
        ofile.write("# Format: lon lat east[m] north[m] up[m]\n")
        for point in disp_points:
            ofile.write("%f %f %f %f %f" % (point.lon, point.lat, point.dE_obs, point.dN_obs, point.dU_obs))
            if point.name != "":
                ofile.write(" %s" % point.name)
            ofile.write("\n")
    easts = [x.dE_obs for x in disp_points]
    norths = [x.dN_obs for x in disp_points]
    print("Min, Max disp_points east: %f, %f m" % (np.min(easts), np.max(easts)))
    print("Min, Max disp_points north: %f, %f m" % (np.min(norths), np.max(norths)))
    return


def write_strain_results(obs_strain_points, strains, outfile):
    """
    obs_strain_points is an object of format cc.dips_points
    strains is a list of tensors

    :param obs_strain_points: a list of disp_points
    :param strains: a list of strain tensors (matrices)
    :param outfile: string, name of output file
    """
    if not obs_strain_points:
        return
    print("Writing file %s " % outfile)
    with open(outfile, 'w') as ofile:
        ofile.write("# Format: lon lat depth_km exx exy exz eyy eyz ezz (microstrain)\n")
        for i in range(len(obs_strain_points)):
            eij = np.multiply(strains[i], 1e6)  # microstrain
            ofile.write("%f %f %f " % (obs_strain_points[i].lon, obs_strain_points[i].lat, obs_strain_points[i].depth))
            ofile.write("%f %f %f " % (eij[0][0], eij[0][1], eij[0][2]))
            ofile.write("%f %f %f\n" % (eij[1][1], eij[1][2], eij[2][2]))
    return


def write_stress_results(obs_strain_points, strains, lame1, mu, outfile):
    """
    Write the full stress tensor at a calculation point.
    obs_strain_points is an object of format cc.dips_points
    strains is a list of tensors

    :param obs_strain_points: a list of disp_points
    :param strains: a list of strain tensors (matrices)
    :param lame1: lame's first parameter (Pa)
    :param mu: shear modulus (Pa)
    :param outfile: string, name of output file
    """
    if not obs_strain_points:
        return
    print("Writing file %s " % outfile)
    with open(outfile, 'w') as ofile:
        ofile.write("# Format: lon lat depth_km sigma_xx sigma_xy sigma_xz sigma_yy sigma_yz sigma_zz (kPa)\n")
        for i in range(len(obs_strain_points)):
            eij = strains[i]
            stress_tensor = conversion_math.get_stress_tensor(eij, lame1, mu)
            stress_tensor = np.divide(stress_tensor, 1000)  # convert to kPa
            ofile.write("%f %f %f " % (obs_strain_points[i].lon, obs_strain_points[i].lat, obs_strain_points[i].depth))
            ofile.write("%f %f %f " % (stress_tensor[0][0], stress_tensor[0][1], stress_tensor[0][2]))
            ofile.write("%f %f %f\n" % (stress_tensor[1][1], stress_tensor[1][2], stress_tensor[2][2]))
    return


def write_results_metrics(inputs, metrics_file, mu):
    """
    :param inputs: inputs object
    :param metrics_file: string, filename
    :param mu: shear modulus, in Pa
    """
    with open(metrics_file, 'w') as f:
        f.write("Number of sources: %d \n" % len(inputs.source_object))
        rect_sources, pt_sources, mogi_sources = utilities.separate_source_types(inputs.source_object)
        if len(rect_sources) > 0:
            mw = moment_calculations.mw_from_moment(pycfaults.get_faults_slip_moment(rect_sources, mu))
            f.write(" ->  Number of rectangular sources: %d\n" % len(rect_sources))
            f.write("Moment Magnitude from Rectangular Fault Patches (assuming G=%.1fGPa): %f\n\n" % (mu / 1e9, mw))
        if len(pt_sources) > 0:
            f.write(" ->  Number of point sources: %d\n" % len(pt_sources))
        if len(mogi_sources) > 0:
            f.write(" ->  Number of Mogi sources: %d\n" % len(mogi_sources))
        f.write("Number of receivers: %d \n" % len(inputs.receiver_object))
        f.write("Coefficient of friction: %f\n" % inputs.FRIC)
    return
