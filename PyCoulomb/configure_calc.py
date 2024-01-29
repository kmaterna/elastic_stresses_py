# Configures a stress calculation

import os, configparser
from . import conversion_math


class Params:
    """
    Only need to provide mu and lame1 elastic parameters.
    Displacement calculations from Okada only use alpha.  Stress uses alpha, mu, lame1, and B.
    alpha is Okada alpha parameter. It is 2/3 for nu=1/4. See DC3D.f docs.   nu = poisson's ratio
    """

    def __init__(self, config_file=None, input_file=None, outdir=os.path.join('output', ''),
                 aftershocks=None, disp_points_file=None, strain_file=None,
                 strike_num_receivers=1, dip_num_receivers=1, fixed_rake=None, mu=30e9, lame1=30e9, B=0,
                 plot_stress=1, plot_grd_disp=1):
        self.config_file = config_file  # string, filename
        self.input_file = input_file  # string, filename
        self.outdir = os.path.join(outdir, '')  # string, directory name
        self.aftershocks = aftershocks  # string, filename
        self.disp_points_file = disp_points_file  # string, filename
        self.strain_file = strain_file  # string, filename
        self.strike_num_receivers = strike_num_receivers  # int, number to split along-strike
        self.dip_num_receivers = dip_num_receivers  # int, number to split along-dip
        self.fixed_rake = fixed_rake  # on receivers, we must specify rake globally if we're using .inp format.
        self.mu = mu  # float, shear modulus
        self.lame1 = lame1  # float, first lame parameter
        self.B = B  # float, Skempton's coefficient
        self.nu, self.alpha = conversion_math.get_poissons_ratio_and_alpha(mu, lame1)  # derived from mu and lame1
        self.plot_stress = plot_stress  # boolean
        self.plot_grd_disp = plot_grd_disp  # boolean
        if self.B > 1:
            raise ValueError("Error! Provided Skepmton's coefficient is invalid because 0 <= B <= 1.")

    def print_summary(self):
        print("Configuring with the following Params:")
        print("  -Config file: %s " % self.config_file)
        print("  -Input file: %s " % self.input_file)
        print("  -Disp_points_file: %s " % self.disp_points_file)
        print("  -Outdir: %s " % self.outdir)
        print("  -Elastic moduli: %f, %f " % (self.mu, self.lame1))
        print("  -Poisson's ratio and Alpha: %f, %f" % (self.nu, self.alpha))
        return

    def write_params_into_config(self, outfile):
        """
        Write an output file from the param object. First you must unpack Params into a ConfigParser.
        Useful for Python API work.
        """
        configobj = configparser.ConfigParser()
        configobj["io-config"] = {}
        ioconfig = configobj["io-config"]
        ioconfig["exp_name"] = os.path.split(os.path.split(self.outdir)[0])[1]
        ioconfig["input_file"] = self.input_file
        ioconfig["output_dir"] = os.path.split(self.outdir)[0]
        ioconfig["plot_stress"] = str(self.plot_stress)
        ioconfig["plot_grd_disp"] = str(self.plot_grd_disp)
        ioconfig["gps_disp_points"] = "" if self.disp_points_file is None else self.disp_points_file
        ioconfig["aftershocks"] = "" if self.aftershocks is None else self.aftershocks
        ioconfig["strain_file"] = "" if self.strain_file is None else self.strain_file
        configobj["compute-config"] = {}
        computeconfig = configobj["compute-config"]
        computeconfig["strike_num_receivers"] = str(self.strike_num_receivers)
        computeconfig["dip_num_receivers"] = str(self.dip_num_receivers)
        computeconfig["mu"] = str(self.mu)
        computeconfig["lame1"] = str(self.lame1)
        computeconfig["B"] = str(self.B)
        computeconfig["fixed_rake"] = str(self.fixed_rake)
        print("Writing file %s " % outfile)
        with open(outfile, 'w') as ofile:
            configobj.write(ofile)
        return

    def modify_params_object(self, config_file=None, input_file=None, aftershocks=None, disp_points_file=None,
                             strain_file=None, strike_num_receivers=None, dip_num_receivers=None, fixed_rake=None,
                             mu=None, lame1=None, B=None, plot_stress=None, plot_grd_disp=None, outdir=None):
        """
        Set the fields in a Pycoulomb.Params object. By default, none of the properties will be altered.
        """
        config_file = self.config_file if config_file is None else config_file
        input_file = self.input_file if input_file is None else input_file
        aftershocks = self.aftershocks if aftershocks is None else aftershocks
        disp_points_file = self.disp_points_file if disp_points_file is None else disp_points_file
        strain_file = self.strain_file if strain_file is None else strain_file
        str_num_receivers = self.strike_num_receivers if strike_num_receivers is None else strike_num_receivers
        dip_num_receivers = self.dip_num_receivers if dip_num_receivers is None else dip_num_receivers
        fixed_rake = self.fixed_rake if fixed_rake is None else fixed_rake
        mu = self.mu if mu is None else mu
        lame1 = self.lame1 if lame1 is None else lame1
        B = self.B if B is None else B
        plot_stress = self.plot_stress if plot_stress is None else plot_stress
        plot_grd_disp = self.plot_grd_disp if plot_grd_disp is None else plot_grd_disp
        outdir = self.outdir if outdir is None else os.path.join(outdir, '')
        MyParams = Params(config_file=config_file, input_file=input_file, aftershocks=aftershocks,
                          disp_points_file=disp_points_file, strain_file=strain_file,
                          strike_num_receivers=str_num_receivers, fixed_rake=fixed_rake,
                          dip_num_receivers=dip_num_receivers, mu=mu, lame1=lame1, B=B,
                          plot_stress=plot_stress, plot_grd_disp=plot_grd_disp, outdir=outdir)
        return MyParams


def configure_stress_calculation(config_file):
    assert (os.path.isfile(config_file)), FileNotFoundError("config file " + config_file + " not found.")
    configobj = configparser.ConfigParser()
    configobj.optionxform = str  # make the config file case-sensitive
    configobj.read(config_file)

    # Basic parameters
    exp_name = configobj.get('io-config', 'exp_name')
    input_file = configobj.get('io-config', 'input_file')
    output_dir = configobj.get('io-config', 'output_dir')
    aftershocks = configobj.get('io-config', 'aftershocks') \
        if configobj.has_option('io-config', 'aftershocks') else None
    gps_file = configobj.get('io-config', 'gps_disp_points') \
        if configobj.has_option('io-config', 'gps_disp_points') else None
    strain_file = configobj.get('io-config', 'strain_file') \
        if configobj.has_option('io-config', 'strain_file') else None
    output_dir = os.path.join(output_dir, exp_name, '')
    if configobj.has_option('io-config', 'plot_stress'):
        plot_stress = configobj.getint('io-config', 'plot_stress')
    else:
        plot_stress = 1
    if configobj.has_option('io-config', 'plot_grd_disp'):
        plot_grd_disp = configobj.getint('io-config', 'plot_grd_disp')
    else:
        plot_grd_disp = 1

    # Computation parameters
    strike_num_receivers = configobj.getint('compute-config', 'strike_num_receivers')
    dip_num_receivers = configobj.getint('compute-config', 'dip_num_receivers')
    mu = configobj.getfloat('compute-config', 'mu')
    lame1 = configobj.getfloat('compute-config', 'lame1')  # this is lamda
    B = configobj.getfloat('compute-config', 'B')
    fixed_rake = configobj.get('compute-config', 'fixed_rake') \
        if configobj.has_option('compute-config', 'fixed_rake') else None
    if '.inp' in input_file:
        assert fixed_rake, ValueError("Must provide fixed_rake for receiver faults in .inp file. ex: 90 (reverse).")

    MyParams = Params(config_file=config_file, input_file=input_file, aftershocks=aftershocks,
                      disp_points_file=gps_file, strain_file=strain_file,
                      strike_num_receivers=strike_num_receivers, fixed_rake=fixed_rake,
                      dip_num_receivers=dip_num_receivers, mu=mu, lame1=lame1, B=B,
                      plot_stress=plot_stress, plot_grd_disp=plot_grd_disp, outdir=output_dir)
    MyParams.print_summary()
    return MyParams


def write_valid_config_file(directory):
    config_filename = "my_config.txt"
    configobj = configparser.ConfigParser()
    configobj.optionxform = str  # case-sensitive config options
    configobj["io-config"] = {}
    ioconfig = configobj["io-config"]
    ioconfig["exp_name"] = 'my_experiment'
    ioconfig["input_file"] = 'M6.8_2014.intxt'
    ioconfig["output_dir"] = os.path.join('Outputs', '')
    ioconfig["plot_stress"] = '1'
    ioconfig["plot_grd_disp"] = '1'
    ioconfig["gps_disp_points"] = 'CA_GPS_ll.txt'
    ioconfig["aftershocks"] = 'CA_aftershocks_2014.txt'
    ioconfig["strain_file"] = ''
    configobj["compute-config"] = {}
    computeconfig = configobj["compute-config"]
    computeconfig["strike_num_receivers"] = '10'
    computeconfig["dip_num_receivers"] = '10'
    computeconfig["mu"] = '30000000000'
    computeconfig["lame1"] = '30000000000'
    computeconfig["B"] = '0'
    computeconfig["fixed_rake"] = ''
    with open(os.path.join(directory, config_filename), 'w') as configfile:
        configobj.write(configfile)
    print("Writing file %s " % os.path.join(directory, config_filename))
    return
