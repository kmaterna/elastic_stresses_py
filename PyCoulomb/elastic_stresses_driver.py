#!/usr/bin/env python

"""
Driver to run Coulomb Stress Calculations in Python
This code uses a Python wrapper for Okada's fortran codes (in particular DC3D.f).
It computes the displacements and strains inside the elastic half-space
From these components, it resolves stresses on receiver faults.
This toolbox performs functions similar to Coulomb in matlab, although in command line.
I found myself running the same thing over and over in Coulomb, so I built a command-line analog.
The inputs and outputs are not as fancy as Coulomb.
It only calculates stress on faults, not cross sections.
"""

import argparse
from Elastic_stresses_py import PyCoulomb


def welcome_and_parse_runstring():
    print("\n\nWelcome to a simple forward modeling tool for calculating elastic displacements and coulomb stresses. ");
    parser = argparse.ArgumentParser(description='Run elastic stress models in Python',
                                     epilog='\U0001f600 \U0001f600 \U0001f600 ');
    parser.add_argument('config', type=str, help='name of config file for calculation. Required.')
    args = parser.parse_args()
    return args;


def drive_calculation(config_file):
    params = PyCoulomb.configure_calc.configure_stress_calculation(config_file);
    [inputs, disp_points] = PyCoulomb.input_values.read_inputs(params);
    out_object = PyCoulomb.run_dc3d.do_stress_computation(params, inputs, disp_points);
    PyCoulomb.output_manager.produce_outputs(params, inputs, disp_points, out_object);
    return;


if __name__ == "__main__":
    args = welcome_and_parse_runstring();
    drive_calculation(args.config);
