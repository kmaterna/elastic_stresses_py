# Configures a stress calculation 

import collections

Params = collections.namedtuple('Params',
	['input_file','x_num_receivers','y_num_receivers','fixed_rake','mu','lame1','alpha','outdir']);


def configure_stress_calculation():
	input_file = "../Inputs/Example-2KZM_LL.inp";
	outdir = "../Outputs/"
	x_num_receivers = 1;
	y_num_receivers = 1;  # how many sub-faults do you want? 

	mu = 30e9; # 30 GPa for shear modulus
	lame1 = 30e9;  # This is LAMDA, but I'm not using Lamda as a variable name. 

	alpha = 2.0/3.0;  # a parameter for the Okada functions. Check if this is appropriate. It is (lamda+mu)/(lamda+2*mu).  See documentation for DC3D.f
	fixed_rake = 90; # on receiver faults, we need to specify the rake globally if we're using .inp format. 90=reverse. 

	MyParams = Params(input_file=input_file, x_num_receivers=x_num_receivers, y_num_receivers=y_num_receivers, fixed_rake=fixed_rake, mu=mu, lame1=lame1, alpha=alpha, outdir=outdir);
	return MyParams;