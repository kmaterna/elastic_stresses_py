# Read the Coulomb or user-defined input files
# Valid formats: 
# .inp (Coulomb)
# .inr (Coulomb)
# .intxt (My own definition, convenient for an experiment I'm running)

import io_inp
import io_inr


def read_inputs(params):
	print("Reading %s. " % params.input_file);
	if '.inp' in params.input_file:
		input_object = io_inp.read_inp(params.input_file,params.fixed_rake);
	elif '.inr' in params.input_file:
		input_object = io_inr.read_inr(params.input_file);
	elif '.intxt' in params.input_file:
		input_object = [];  # Will write later. 
	else:
		print("Error! Unrecognized type of input file!");
		input_object = [];
	return input_object;