# Configures a stress calculation 

import os
import configparser
from . import coulomb_collections as cc


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
    output_dir = output_dir + exp_name + '/';

    # Computation parameters
    strike_num_receivers = configobj.getint('compute-config', 'strike_num_receivers');
    dip_num_receivers = configobj.getint('compute-config', 'dip_num_receivers');
    mu = configobj.getfloat('compute-config', 'mu');
    lame1 = configobj.getfloat('compute-config', 'lame1');  # this is lambda
    B = configobj.getfloat('compute-config', 'B');
    alpha = (lame1 + mu) / (lame1 + 2 * mu);
    # alpha = parameter for Okada functions. It is 2/3 for simplest case. See DC3D.f documentation.
    fixed_rake = configobj.getfloat('compute-config', 'fixed_rake') \
        if configobj.has_option('compute-config', 'fixed_rake') else None;
    # on receiver faults, we need to specify rake globally if we're using .inp format. No effect for other file formats.
    if '.inp' in input_file:
        assert fixed_rake, ValueError("Must provide fixed_rake for receiver faults in .inp file. ex: 90 (reverse).");

    MyParams = cc.Params(config_file=config_file, input_file=input_file, aftershocks=aftershocks,
                         disp_points_file=gps_file, strain_file=strain_file,
                         strike_num_receivers=strike_num_receivers, fixed_rake=fixed_rake,
                         dip_num_receivers=dip_num_receivers, mu=mu, lame1=lame1, B=B,
                         alpha=alpha, outdir=output_dir);
    print(MyParams);
    return MyParams;


def write_valid_config_file(directory):
    configobj = configparser.ConfigParser()
    configobj["io-config"] = {};
    ioconfig = configobj["io-config"];
    ioconfig["exp_name"] = 'my_experiment';
    ioconfig["input_file"] = 'my_input.intxt';
    ioconfig["output_dir"] = 'Outputs/';
    ioconfig["aftershocks"] = '[optional]';
    ioconfig["gps_disp_points"] = '[optional]';
    ioconfig["strain_file"] = '[optional]';
    configobj["compute-config"] = {};
    computeconfig = configobj["compute-config"];
    computeconfig["strike_num_receivers"] = '10';
    computeconfig["dip_num_receivers"] = '10';
    computeconfig["mu"] = '30000000';
    computeconfig["lame1"] = '30000000';
    computeconfig["B"] = '0';
    computeconfig["fixed_rake"] = '[optional]';
    with open(directory+'/dummy_config.txt', 'w') as configfile:
        configobj.write(configfile)
    print("Writing file %s " % directory+"/dummy_config.txt");
    return;
