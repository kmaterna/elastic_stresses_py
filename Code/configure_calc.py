# Configures a stress calculation 

import argparse, configparser
import coulomb_collections


def welcome_and_parse():
	print("\n\nWelcome to a simple forward modeling tool for calculating elastic displacements and coulomb stresses. ");
	parser = argparse.ArgumentParser(description='Run elastic stress models in Python', epilog='\U0001f600 \U0001f600 \U0001f600 ');
	parser.add_argument('config',type=str,help='name of config file for calculation. Required.')
	args = parser.parse_args()
	print("Config file:",args.config);
	return args;

def configure_stress_calculation(config_file):
	configobj=configparser.ConfigParser();
	configobj.optionxform = str # make the config file case-sensitive
	configobj.read(config_file);

	# Basic parameters
	exp_name=configobj.get('io-config','exp_name');
	title=configobj.get('io-config','title');
	aftershocks=configobj.get('io-config','aftershocks');
	input_file=configobj.get('io-config','input_file');
	gps_disp_points=configobj.get('io-config','gps_disp_points');
	output_dir=configobj.get('io-config','output_dir');
	output_dir = output_dir+exp_name+'/';

	# Computation parameters
	strike_num_receivers = configobj.getint('compute-config','strike_num_receivers');
	dip_num_receivers = configobj.getint('compute-config','dip_num_receivers');
	mu = configobj.getfloat('compute-config','mu');
	lame1 = configobj.getfloat('compute-config','lame1');  # this is lambda
	alpha = (lame1+mu)/(lame1+2*mu);  # a parameter for Okada functions. It is 2/3 for simplest case.  See documentation for DC3D.f
	fixed_rake = configobj.getfloat('compute-config','fixed_rake');
	# on receiver faults, we need to specify the rake globally if we're using .inp format. 90=reverse. 
	# No effect if you're using .inr or .intxt format. 

	MyParams = coulomb_collections.Params(input_file=input_file, aftershocks=aftershocks, disp_points_file=gps_disp_points,
		strike_num_receivers=strike_num_receivers, dip_num_receivers=dip_num_receivers, fixed_rake=fixed_rake, 
		mu=mu, lame1=lame1, alpha=alpha, outdir=output_dir, title=title);
	print(MyParams)
	return MyParams;	