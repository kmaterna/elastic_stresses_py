# Configures a stress calculation

import os, configparser
from . import coulomb_collections as cc
from . import conversion_math


def configure_stress_calculation(config_file):
    print("Config file: ", config_file);
    assert(os.path.isfile(config_file)), FileNotFoundError("config file "+config_file+" not found.");

    configobj = configparser.ConfigParser();
    configobj.optionxform = str  # make the config file case-sensitive
    configobj.read(config_file);

    # Basic parameters
    exp_name = configobj.get('io-config', 'exp_name');
    input_file = configobj.get('io-config', 'input_file');
    output_dir = configobj.get('io-config', 'output_dir');
    aftershocks = configobj.get('io-config', 'aftershocks') \
        if configobj.has_option('io-config', 'aftershocks') else None;
    gps_file = configobj.get('io-config', 'gps_disp_points') \
        if configobj.has_option('io-config', 'gps_disp_points') else None;
    strain_file = configobj.get('io-config', 'strain_file') \
        if configobj.has_option('io-config', 'strain_file') else None;
    if output_dir[-1] != '/':
        output_dir = output_dir + '/';
    output_dir = output_dir + exp_name + '/';
    if configobj.has_option('io-config', 'plot_stress'):
        plot_stress = configobj.getint('io-config', 'plot_stress');
    else:
        plot_stress = 1;
    if configobj.has_option('io-config', 'plot_grd_disp'):
        plot_grd_disp = configobj.getint('io-config', 'plot_grd_disp');
    else:
        plot_grd_disp = 1;

    # Computation parameters
    strike_num_receivers = configobj.getint('compute-config', 'strike_num_receivers');
    dip_num_receivers = configobj.getint('compute-config', 'dip_num_receivers');
    mu = configobj.getfloat('compute-config', 'mu');
    lame1 = configobj.getfloat('compute-config', 'lame1');  # this is lambda
    B = configobj.getfloat('compute-config', 'B');
    [poissons_ratio, alpha] = conversion_math.get_poissons_ratio_and_alpha(mu, lame1);
    print("For this stress calculation:\n--> Poisson's ratio = %f, alpha = %f" % (poissons_ratio, alpha));
    # alpha = parameter for Okada. It is 2/3 for nu=1/4. See DC3D.f docs.
    fixed_rake = configobj.get('compute-config', 'fixed_rake') \
        if configobj.has_option('compute-config', 'fixed_rake') else None;
    # on receiver faults, we need to specify rake globally if we're using .inp format. No effect for other file formats.
    if '.inp' in input_file:
        assert fixed_rake, ValueError("Must provide fixed_rake for receiver faults in .inp file. ex: 90 (reverse).");

    MyParams = cc.Params(config_file=config_file, input_file=input_file, aftershocks=aftershocks,
                         disp_points_file=gps_file, strain_file=strain_file,
                         strike_num_receivers=strike_num_receivers, fixed_rake=fixed_rake,
                         dip_num_receivers=dip_num_receivers, mu=mu, lame1=lame1, B=B,
                         alpha=alpha, plot_stress=plot_stress, plot_grd_disp=plot_grd_disp,
                         outdir=output_dir);
    print(MyParams);
    return MyParams;


def write_params_into_config(params, outfile):
    """
    Write an output file from the param object. First you must unpack Params into a ConfigParser.
    Useful for Python API work.
    """
    configobj = configparser.ConfigParser()
    configobj["io-config"] = {};
    ioconfig = configobj["io-config"];
    ioconfig["exp_name"] = params.outdir.split('/')[-2]
    ioconfig["input_file"] = params.input_file
    ioconfig["output_dir"] = '/'.join(params.outdir.split('/')[0:-2])+'/'
    ioconfig["plot_stress"] = str(params.plot_stress)
    ioconfig["plot_grd_disp"] = str(params.plot_grd_disp)
    ioconfig["gps_disp_points"] = "" if params.disp_points_file is None else params.disp_points_file
    ioconfig["aftershocks"] = "" if params.aftershocks is None else params.aftershocks
    ioconfig["strain_file"] = "" if params.strain_file is None else params.strain_file
    configobj["compute-config"] = {};
    computeconfig = configobj["compute-config"];
    computeconfig["strike_num_receivers"] = str(params.strike_num_receivers)
    computeconfig["dip_num_receivers"] = str(params.dip_num_receivers)
    computeconfig["mu"] = str(params.mu)
    computeconfig["lame1"] = str(params.lame1)
    computeconfig["B"] = str(params.B)
    computeconfig["fixed_rake"] = str(params.fixed_rake)
    print("Writing file %s " % outfile);
    with open(outfile, 'w') as ofile:
        configobj.write(ofile);
    return;


def configure_default_displacement_params(outdir='output/', plot_stress=1, plot_grd_disp=1, config_file=None,
                                          input_file=None, aftershocks=None, disp_points_file=None, strain_file=None,
                                          strike_num_receivers=1, dip_num_receivers=1, fixed_rake=0,
                                          mu=30e9, lame1=30e9, B=0):
    """
    Build a default Params object used for displacement-only calculations.
    All arguments are optional.
    Displacement from Okada only uses alpha.  Stress uses alpha, mu, lame1, and B.
    Sends the output to an outdir directory.
    """
    if outdir[-1] != '/':
        outdir += '/';
    [pr, alpha] = conversion_math.get_poissons_ratio_and_alpha(mu, lame1);  # derived from provided mu and lame1
    print("For this stress calculation:\n--> Poisson's ratio = %f, alpha = %f" % (pr, alpha));
    MyParams = cc.Params(config_file=config_file, input_file=input_file, aftershocks=aftershocks,
                         disp_points_file=disp_points_file, strain_file=strain_file,
                         strike_num_receivers=strike_num_receivers, fixed_rake=fixed_rake,
                         dip_num_receivers=dip_num_receivers, mu=mu, lame1=lame1, B=B,
                         alpha=alpha, plot_stress=plot_stress, plot_grd_disp=plot_grd_disp, outdir=outdir);
    return MyParams;


def modify_params_object(default_params, config_file=None, input_file=None, aftershocks=None, disp_points_file=None,
                         strain_file=None, strike_num_receivers=None, dip_num_receivers=None, fixed_rake=None,
                         mu=None, lame1=None, B=None, plot_stress=None, plot_grd_disp=None, outdir=None):
    """
    Modify the fields in a Pycoulomb.Params named tuple
    By default, none of the properties will be altered.

    :param default_params: (required) existing named tuple
    :param config_file: optional, string
    :param input_file: optional, string
    :param aftershocks: optional, string
    :param disp_points_file: optional, string
    :param strain_file: optional, string
    :param strike_num_receivers: optional, int
    :param dip_num_receivers: optional, int
    :param fixed_rake: optional, float
    :param mu: optional, float
    :param lame1: optional, float
    :param B: optional, float
    :param plot_stress: optional, int
    :param plot_grd_disp: optional, int
    :param outdir: optional, float
    """
    config_file = default_params.config_file if config_file is None else config_file;
    input_file = default_params.input_file if input_file is None else input_file;
    aftershocks = default_params.aftershocks if aftershocks is None else aftershocks;
    disp_points_file = default_params.disp_points_file if disp_points_file is None else disp_points_file;
    strain_file = default_params.strain_file if strain_file is None else strain_file;
    str_num_receivers = default_params.strike_num_receivers if strike_num_receivers is None else strike_num_receivers;
    dip_num_receivers = default_params.dip_num_receivers if dip_num_receivers is None else dip_num_receivers;
    fixed_rake = default_params.fixed_rake if fixed_rake is None else fixed_rake;
    mu = default_params.mu if mu is None else mu;
    lame1 = default_params.lame1 if lame1 is None else lame1;
    B = default_params.B if B is None else B;
    plot_stress = default_params.plot_stress if plot_stress is None else plot_stress;
    plot_grd_disp = default_params.plot_grd_disp if plot_grd_disp is None else plot_grd_disp;
    outdir = default_params.outdir if outdir is None else outdir;
    if outdir[-1] != '/':
        outdir += '/';
    [pr, alpha] = conversion_math.get_poissons_ratio_and_alpha(mu, lame1);  # derived from provided mu and lame1
    print("For this stress calculation:\n--> Poisson's ratio = %f, alpha = %f" % (pr, alpha));

    MyParams = cc.Params(config_file=config_file, input_file=input_file, aftershocks=aftershocks,
                         disp_points_file=disp_points_file, strain_file=strain_file,
                         strike_num_receivers=str_num_receivers, fixed_rake=fixed_rake,
                         dip_num_receivers=dip_num_receivers, mu=mu, lame1=lame1, B=B,
                         alpha=alpha, plot_stress=plot_stress, plot_grd_disp=plot_grd_disp,
                         outdir=outdir);
    return MyParams;


def configure_default_displacement_input(source_object, zerolon, zerolat, bbox, domainsize=20):
    """
    Build a default Input object for displacement-only calculations.
    Unused parameters get default values or None.
    domainsize: in km. default is 20.
    bbox: [w, e, s, n]
    """
    inputs = cc.Input_object(PR1=0.25, FRIC=0.4, depth=0,
                             start_gridx=-domainsize, finish_gridx=domainsize,
                             start_gridy=-domainsize, finish_gridy=domainsize,
                             xinc=domainsize/20, yinc=domainsize/20,
                             minlon=bbox[0], maxlon=bbox[1], zerolon=zerolon,
                             minlat=bbox[2], maxlat=bbox[3], zerolat=zerolat,
                             source_object=source_object, receiver_object=[],
                             receiver_horiz_profile=None);
    return inputs;


def modify_inputs_object(default_inputs, PR1=None, FRIC=None, depth=None, start_gridx=None, finish_gridx=None,
                         start_gridy=None, finish_gridy=None, xinc=None, yinc=None, minlon=None, maxlon=None,
                         zerolon=None, minlat=None, maxlat=None, zerolat=None, source_object=None, receiver_object=None,
                         receiver_horiz_profile=None):
    """
    Modify the fields in a PyCoulomb.Input_object namedtuple
    """

    PR1 = default_inputs.PR1 if PR1 is None else PR1;
    FRIC = default_inputs.FRIC if FRIC is None else FRIC;
    depth = default_inputs.depth if depth is None else depth;
    start_gridx = default_inputs.start_gridx if start_gridx is None else start_gridx;
    finish_gridx = default_inputs.finish_gridx if finish_gridx is None else finish_gridx;
    start_gridy = default_inputs.start_gridy if start_gridy is None else start_gridy;
    finish_gridy = default_inputs.finish_gridy if finish_gridy is None else finish_gridy;
    xinc = default_inputs.xinc if xinc is None else xinc;
    yinc = default_inputs.yinc if yinc is None else yinc;
    minlon = default_inputs.minlon if minlon is None else minlon;
    maxlon = default_inputs.maxlon if maxlon is None else maxlon;
    zerolon = default_inputs.zerolon if zerolon is None else zerolon;
    minlat = default_inputs.minlat if minlat is None else minlat;
    maxlat = default_inputs.maxlat if maxlat is None else maxlat;
    zerolat = default_inputs.zerolat if zerolat is None else zerolat;
    source_object = default_inputs.source_object if source_object is None else source_object;
    receiver_object = default_inputs.receiver_object if receiver_object is None else receiver_object;
    rec_profile = default_inputs.receiver_horiz_profile if receiver_horiz_profile is None else receiver_horiz_profile;

    modified_inputs = cc.Input_object(PR1=PR1, FRIC=FRIC, depth=depth, start_gridx=start_gridx,
                                      finish_gridx=finish_gridx, start_gridy=start_gridy, finish_gridy=finish_gridy,
                                      xinc=xinc, yinc=yinc, minlon=minlon, maxlon=maxlon, zerolon=zerolon,
                                      minlat=minlat, maxlat=maxlat, zerolat=zerolat,
                                      source_object=source_object, receiver_object=receiver_object,
                                      receiver_horiz_profile=rec_profile);
    return modified_inputs;


def write_valid_config_file(directory):
    config_filename = "my_config.txt";
    configobj = configparser.ConfigParser()
    configobj.optionxform = str   # case-sensitive config options
    configobj["io-config"] = {};
    ioconfig = configobj["io-config"];
    ioconfig["exp_name"] = 'my_experiment';
    ioconfig["input_file"] = 'M6.8_2014.intxt';
    ioconfig["output_dir"] = 'Outputs/';
    ioconfig["plot_stress"] = '1'
    ioconfig["plot_grd_disp"] = '1'
    ioconfig["gps_disp_points"] = 'CA_GPS_ll.txt';
    ioconfig["aftershocks"] = 'CA_aftershocks_2014.txt';
    ioconfig["strain_file"] = ''
    configobj["compute-config"] = {};
    computeconfig = configobj["compute-config"];
    computeconfig["strike_num_receivers"] = '10';
    computeconfig["dip_num_receivers"] = '10';
    computeconfig["mu"] = '30000000000';
    computeconfig["lame1"] = '30000000000';
    computeconfig["B"] = '0';
    computeconfig["fixed_rake"] = '';
    with open(directory+'/'+config_filename, 'w') as configfile:
        configobj.write(configfile)
    print("Writing file %s " % directory+"/"+config_filename);
    return;


def modify_fault_object(default_fault, xstart=None, xfinish=None, ystart=None, yfinish=None, rtlat=None,
                        reverse=None, tensile=None, potency=None, strike=None, dipangle=None, rake=None, zerolon=None,
                        zerolat=None, top=None, bottom=None):
    """
    Modify the fields in a Pycoulomb.Faults_object namedtuple
    """

    xstart = default_fault.xstart if xstart is None else xstart;
    xfinish = default_fault.xfinish if xfinish is None else xfinish;
    ystart = default_fault.ystart if ystart is None else ystart;
    yfinish = default_fault.yfinish if yfinish is None else yfinish;
    rtlat = default_fault.rtlat if rtlat is None else rtlat;
    reverse = default_fault.reverse if reverse is None else reverse;
    tensile = default_fault.tensile if tensile is None else tensile;
    potency = default_fault.potency if potency is None else potency;
    strike = default_fault.strike if strike is None else strike;
    dipangle = default_fault.dipangle if dipangle is None else dipangle;
    rake = default_fault.rake if rake is None else rake;
    zerolon = default_fault.zerolon if zerolon is None else zerolon;
    zerolat = default_fault.zerolat if zerolat is None else zerolat;
    top = default_fault.top if top is None else top;
    bottom = default_fault.bottom if bottom is None else bottom;
    new_fault = cc.construct_pycoulomb_fault(xstart=xstart, xfinish=xfinish, ystart=ystart, yfinish=yfinish,
                                             rtlat=rtlat, reverse=reverse, tensile=tensile,
                                             potency=potency, strike=strike, dipangle=dipangle, rake=rake,
                                             zerolon=zerolon, zerolat=zerolat, top=top, bottom=bottom);
    return new_fault;
