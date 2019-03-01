# Read the Coulomb or user-defined input files
# Valid formats: 
# .inp (Coulomb)
# .inr (Coulomb)
# .intxt (My own definition, convenient for an experiment I'm running)

import io_inp
import io_inr
import io_intxt
import io_additionals


def read_inputs(params):
	print("Reading %s. " % params.input_file);
	if '.inp' in params.input_file:
		input_object = io_inp.read_inp(params.input_file,params.fixed_rake);  # fixed rake format
	elif '.inr' in params.input_file:
		input_object = io_inr.read_inr(params.input_file);  # variable rake format (will write later);
	elif '.intxt' in params.input_file:
		input_object = io_intxt.read_intxt(params.input_file);  # convenient input format (fixed rake)
	else:
		print("Error! Unrecognized type of input file!");
		input_object = [];
	if params.disp_points_file=='':
		print("Not reading any GPS points for displacements.");
		disp_points = [];
	else:
		disp_points = io_additionals.read_disp_points(params.disp_points_file);

	return [input_object, disp_points];