# Configures a stress calculation 

import os
import configparser
from . import coulomb_collections as cc


def configure_stress_calculation(config_file):
    print("Config file: ", config_file);
    assert(os.path.isfile(config_file)), FileNotFoundError("config file "+config_file+"not found.");

    configobj = configparser.ConfigParser();
    configobj.optionxform = str  # make the config file case-sensitive
    configobj.read(config_file);

    # Basic parameters
    exp_name = configobj.get('io-config', 'exp_name');
    input_file = configobj.get('io-config', 'input_file');
    output_dir = configobj.get('io-config', 'output_dir');
    aftershocks = configobj.get('io-config', 'aftershocks') if configobj.has_option('io-config', 'aftershocks') else None;
    gps_file = configobj.get('io-config', 'gps_disp_points') if configobj.has_option('io-config', 'gps_disp_points') else None;
    output_dir = output_dir + exp_name + '/';

    # Computation parameters
    strike_num_receivers = configobj.getint('compute-config', 'strike_num_receivers');
    dip_num_receivers = configobj.getint('compute-config', 'dip_num_receivers');
    mu = configobj.getfloat('compute-config', 'mu');
    lame1 = configobj.getfloat('compute-config', 'lame1');  # this is lambda
    alpha = (lame1 + mu) / (lame1 + 2 * mu);
    # alpha = parameter for Okada functions. It is 2/3 for simplest case. See DC3D.f documentation.
    fixed_rake = configobj.getfloat('compute-config', 'fixed_rake');
    # on receiver faults, we need to specify rake globally if we're using .inp format. 90=reverse.
    # No effect if using .inr, .inzero, or .intxt format.

    MyParams = cc.Params(config_file=config_file, input_file=input_file, aftershocks=aftershocks,
                         disp_points_file=gps_file, strike_num_receivers=strike_num_receivers, fixed_rake=fixed_rake,
                         dip_num_receivers=dip_num_receivers, mu=mu, lame1=lame1, alpha=alpha, outdir=output_dir);
    print(MyParams);
    return MyParams;
