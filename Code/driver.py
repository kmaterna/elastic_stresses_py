# Driver to run Coulomb Stress Calculations in Python
# This code uses Ben Thompson's Python wrapper for Okada's fortran codes (in particular DC3D.f).
# It computes the displacements and strains inside the half-space
# From these components, it resolves stresses on receiver faults, including Coulomb stresses. 
# This toolbox performs functions similar to Coulomb in matlab, although in command line.
# I found myself running the same thing over and over in Coulomb, so I built a command-line analog. 
# The inputs and outputs are not as fancy as Coulomb. 
# It only calculates stress on faults. 
# Does not have functionality to do cross-sections yet. 
# Requirements: Python 3, Ben Thompson's DC3D.f python wrapper on your pythonpath. 
# Convention: positive strike slip is right-lateral. 
# Convention: positive dip slip is reverse. 

import configure_calc
import input_values
import run_dc3d
import output_manager
import sys

def do_calculation():
	params = configure_calc.configure_stress_calculation();
	[inputs, disp_points] = input_values.read_inputs(params);
	out_object = run_dc3d.do_stress_computation(params, inputs, disp_points);
	output_manager.produce_outputs(params, inputs, disp_points, out_object);
	return;





do_calculation();