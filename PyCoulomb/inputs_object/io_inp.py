"""
Read Coulomb input files in the .inp format. Important parameters are:
* 1. Poisson's Ratio
* 2. Coefficient of Friction
* 3. Depth
* 4. Fault params (x start, x finish, y start, y finish, Kode, rt. lat, reverse, dip angle, top, bottom, comment)
* 5. Grid parameters (start_gridx, finish_gridx, start_gridy, finish_gridy, xinc, yinc)
* 6. Map info (min lon, max lon, zero lon, min lat, max lat, zero lat)
"""

from .. import utilities
from .input_obj import Input_object
from ..pyc_fault_object import Faults_object
from Tectonic_Utils.geodesy import fault_vector_functions
import subprocess


def read_inp(input_file, fixed_rake):
    fixed_rake = float(fixed_rake)
    print("Reading source and receiver information from file %s with fixed_rake = %f " % (input_file, fixed_rake))
    # inp files require fixed rake for receiver faults, since they don't provide a fault-specific one in input file.

    [minlon, maxlon, zerolon, minlat, maxlat, zerolat] = get_map_info(input_file)
    [start_gridx, start_gridy, finish_gridx, finish_gridy, xinc, yinc] = get_grid_parameters(input_file)

    read_faults = 0
    sources, receivers = [], []
    PR1, FRIC, depth = 0, 0, 0

    try:  # if the rake is explicitly defined, as in USGS NEIC FFM Coulomb solution
        grep_result = subprocess.check_output("grep netslip " + input_file, shell=True)
    except subprocess.CalledProcessError:
        grep_result = []
    if grep_result:
        rake_is_explicit = 1
    else:
        rake_is_explicit = 0

    ifile = open(input_file)
    for line in ifile:
        temp = line.split()
        if 'PR1=' in line:
            prspot = temp.index('PR1=')
            PR1 = float(temp[prspot + 1])
            dpspot = temp.index('DEPTH=')
            depth = float(temp[dpspot + 1])
        if 'FRIC=' in line:
            fric_place = temp.index('FRIC=')
            FRIC = float(temp[fric_place + 1])
        if 'xxxxxxxxxx' in line:  # Moving into the fault definitions
            read_faults = 1
            continue
        if read_faults == 1:
            if len(temp) == 0:
                read_faults = 0
                continue  # Getting out of fault definitions
            else:
                slip = abs(float(temp[6])) + abs(float(temp[7]))
                if slip > 0.0000001:  # Here we have a source fault
                    [xstart, ystart, xfinish, yfinish, Kode, rtlat, reverse, strike, dipangle, top, bottom,
                     comment] = read_fault_line(line, rake_is_explicit=rake_is_explicit)
                    rake = fault_vector_functions.get_rake(rtlat_strike_slip=rtlat, dip_slip=reverse)
                    one_source_object = Faults_object(xstart=xstart, xfinish=xfinish, ystart=ystart, yfinish=yfinish,
                                                      Kode=Kode, rtlat=rtlat, reverse=reverse, strike=strike,
                                                      dipangle=dipangle, rake=rake, top=top,
                                                      bottom=bottom, zerolon=zerolon, zerolat=zerolat, comment=comment)
                    sources.append(one_source_object)
                else:  # here we have a receiver fault
                    [xstart, ystart, xfinish, yfinish, Kode, _, _, strike, dipangle, top, bottom,
                     comment] = read_fault_line(line)
                    rake = fixed_rake
                    one_receiver_object = Faults_object(xstart=xstart, xfinish=xfinish, ystart=ystart, yfinish=yfinish,
                                                        Kode=Kode, rtlat=0, reverse=0, strike=strike, dipangle=dipangle,
                                                        rake=rake, top=top, zerolon=zerolon, zerolat=zerolat,
                                                        bottom=bottom, comment=comment)
                    receivers.append(one_receiver_object)
    ifile.close()

    my_inputs = Input_object(PR1=PR1, FRIC=FRIC, depth=depth, start_gridx=start_gridx,
                             finish_gridx=finish_gridx, start_gridy=start_gridy, finish_gridy=finish_gridy,
                             xinc=xinc, yinc=yinc, minlon=minlon, maxlon=maxlon, zerolon=zerolon,
                             minlat=minlat, maxlat=maxlat, zerolat=zerolat,
                             receiver_object=receivers, source_object=sources, receiver_horiz_profile=None)
    return my_inputs


def read_fault_line(line, rake_is_explicit=0):
    """
    By default, we read #   X-start    Y-start     X-fin      Y-fin   Kode  rtlat  reverse  dip angle top  bot
    If required, we read #   X-start    Y-start     X-fin      Y-fin   Kode  rake   netslip  dip angle top  bot
    """
    temp = line.split()
    xstart = float(temp[1])
    ystart = float(temp[2])
    xfinish = float(temp[3])
    yfinish = float(temp[4])
    Kode = int(temp[5])
    if rake_is_explicit:
        rake = float(temp[6])
        netslip = float(temp[7])
        lftlat, reverse = fault_vector_functions.get_leftlat_reverse_slip(netslip, rake)
        rtlat = -lftlat
    else:
        rtlat = float(temp[6])
        reverse = float(temp[7])
    dipangle = float(temp[8])
    top = float(temp[9])
    bottom = float(temp[10])
    comment = " ".join(temp[11:-1])
    strike = fault_vector_functions.get_strike(xfinish - xstart, yfinish - ystart)
    return [xstart, ystart, xfinish, yfinish, Kode, rtlat, reverse, strike, dipangle, top, bottom, comment]


def get_grid_parameters(input_file):
    # Open a Coulomb file and get the grid information
    # Moving into Grid Parameters and Map Information
    read_grid = 0
    start_gridx, start_gridy, finish_gridx, finish_gridy, xinc, yinc = 0, 0, 0, 0, 0, 0
    ifile = open(input_file)
    for line in ifile:
        temp = line.split()
        if 'Grid Parameters' in line:
            read_grid = 1
            continue
        if read_grid == 1:
            if '  1  -' in line:
                start_gridx = float(temp[-1])
            elif '  2  -' in line:
                start_gridy = float(temp[-1])
            elif '  3  -' in line:
                finish_gridx = float(temp[-1])
            elif '  4  -' in line:
                finish_gridy = float(temp[-1])
            elif '  5  -' in line:
                xinc = float(temp[-1])
            elif '  6  -' in line:
                yinc = float(temp[-1])
            else:
                read_grid = 0
                continue
    ifile.close()
    return [start_gridx, start_gridy, finish_gridx, finish_gridy, xinc, yinc]


def get_map_info(input_file):
    # Open a Coulomb file and get the map information
    read_maps = 0
    minlon, maxlon, zerolon, minlat, maxlat, zerolat = 0, 0, 0, 0, 0, 0
    ifile = open(input_file)
    for line in ifile:
        temp = line.split()
        # Reading Map Information
        if 'Map info' in line:
            read_maps = 1
            continue
        if read_maps == 1:
            if '  1  -' in line:
                minlon = float(temp[-1])
            elif '  2  -' in line:
                maxlon = float(temp[-1])
            elif '  3  -' in line:
                zerolon = float(temp[-1])
            elif '  4  -' in line:
                minlat = float(temp[-1])
            elif '  5  -' in line:
                maxlat = float(temp[-1])
            elif '  6  -' in line:
                zerolat = float(temp[-1])
            else:
                read_maps = 0
                continue
    ifile.close()
    return [minlon, maxlon, zerolon, minlat, maxlat, zerolat]


def write_inp(data_obj, outfile):
    rectangles, points, _ = utilities.separate_source_types(data_obj.source_object)
    sources = rectangles + points
    receivers = data_obj.receiver_object
    num_faults = len(sources) + len(receivers)
    ofile = open(outfile, 'w')
    ofile.write('This is a test file for Coulomb 3.0\n')
    ofile.write('This file is prepared by Python\n')
    ofile.write('#reg1=  0  #reg2=  0   #fixed=  %d  sym=  1\n' % num_faults)
    ofile.write(' PR1=       %4.3f     PR2=       .250    DEPTH=        %.1f\n' % (data_obj.PR1, data_obj.depth))
    ofile.write('  E1=   0.800000E+06   E2=   0.800000E+06\n')  # The young's modulus, in bar
    ofile.write('XSYM=       .000     YSYM=       .000\n')  # An old pattern- ignore
    ofile.write('FRIC=       %4.3f\n' % data_obj.FRIC)
    ofile.write('S1DR=    19.0001     S1DP=     -0.0001    S1IN=    100.000     S1GD=   .000000\n')
    ofile.write('S3DR=    89.9999     S3DP=      89.999    S3IN=     30.000     S3GD=   .000000\n')
    ofile.write('S2DR=   109.0001     S2DP=     -0.0001    S2IN=      0.000     S2GD=   .000000\n\n')
    ofile.write(
        '  #   X-start    Y-start     X-fin      Y-fin   Kode  rt.lat    reverse   dip angle     top      bot\n')
    ofile.write(
        'xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx\n')
    for fault in sources:
        ofile.write('  1 % 10.4f % 10.4f % 10.4f % 10.4f %3d % 10.4f % 10.4f % 10.4f % 10.2f % 10.2f %s\n' % (
            fault.xstart, fault.ystart, fault.xfinish, fault.yfinish,
            fault.Kode, fault.rtlat, fault.reverse, fault.dipangle, fault.top, fault.bottom, fault.comment))
    for fault in receivers:
        ofile.write('  1 % 10.4f % 10.4f % 10.4f % 10.4f %3d % 10.4f % 10.4f % 10.4f % 10.2f % 10.2f %s\n' % (
            fault.xstart, fault.ystart, fault.xfinish, fault.yfinish,
            fault.Kode, fault.rtlat, fault.reverse, fault.dipangle, fault.top, fault.bottom, fault.comment))
    ofile.write('\n')
    ofile.write('    Grid Parameters\n')
    ofile.write('  1  ----------------------------  Start-x =    % 10.5f\n' % data_obj.start_gridx)
    ofile.write('  2  ----------------------------  Start-y =    % 10.5f\n' % data_obj.start_gridy)
    ofile.write('  3  --------------------------   Finish-x =    % 10.5f\n' % data_obj.finish_gridx)
    ofile.write('  4  --------------------------   Finish-y =    % 10.5f\n' % data_obj.finish_gridy)
    ofile.write('  5  ------------------------  x-increment =    % 10.5f\n' % data_obj.xinc)
    ofile.write('  6  ------------------------  y-increment =    % 10.5f\n' % data_obj.yinc)
    ofile.write('     Size Parameters\n')
    ofile.write('  1  --------------------------  Plot size =     2.000000\n')
    ofile.write('  2  --------------  Shade/Color increment =     1.000000\n')
    ofile.write('  3  ------  Exaggeration for disp.& dist. =     10000.00\n\n')
    ofile.write(
        'Cross section default\n')
    # Because the python program doesn't do cross sections yet, I'm leaving this hard-coded.
    ofile.write('  1  ----------------------------  Start-x =    -36.00000\n')
    ofile.write('  2  ----------------------------  Start-y =     36.00000\n')
    ofile.write('  3  --------------------------   Finish-x =     38.00000\n')
    ofile.write('  4  --------------------------   Finish-y =    -36.00000\n')
    ofile.write('  5  ------------------  Distant-increment =     1.000000\n')
    ofile.write('  6  ----------------------------  Z-depth =     30.00000\n')
    ofile.write('  7  ------------------------  Z-increment =     1.000000\n')
    ofile.write('     Map infomation\n')
    ofile.write('  1  ---------------------------- min. lon =    % 9.4f \n' % data_obj.minlon)
    ofile.write('  2  ---------------------------- max. lon =    % 9.4f \n' % data_obj.maxlon)
    ofile.write('  3  ---------------------------- zero lon =    % 9.4f \n' % data_obj.zerolon)
    ofile.write('  4  ---------------------------- min. lat =    % 9.4f \n' % data_obj.minlat)
    ofile.write('  5  ---------------------------- max. lat =    % 9.4f \n' % data_obj.maxlat)
    ofile.write('  6  ---------------------------- zero lat =    % 9.4f \n' % data_obj.zerolat)
    ofile.close()
    print("Writing outfile %s " % outfile)
    return
