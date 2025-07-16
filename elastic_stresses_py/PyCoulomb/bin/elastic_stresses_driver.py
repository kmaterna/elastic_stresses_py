#!/usr/bin/env python

"""
Driver to run Coulomb Stress Calculations in Python
This code uses a Python wrapper for Okada's fortran codes (in particular DC3D.f).
It computes the displacements and strains inside the elastic half-space
From these components, it resolves stresses on receiver faults.
This toolbox performs functions similar to Coulomb in matlab, although in command line.
I found myself running the same thing over and over in Coulomb, so I built a command-line analog.
The inputs and outputs are not as fancy as Coulomb.
"""

import argparse
import elastic_stresses_py.PyCoulomb as PyCoulomb


def welcome_and_parse_runstring():
    print("\n\nWelcome to a simple forward modeling tool for calculating elastic displacements and coulomb stresses. ")
    parser = argparse.ArgumentParser(description='Run elastic stress models in Python',
                                     epilog='\U0001f600 \U0001f600 \U0001f600 ')
    parser.add_argument('config', type=str,
                        help='name of config file for calculation. Required.')
    parser.add_argument('--input_file', type=str,
                        help='override the input file from config file. Optional.', required=False)
    parser.add_argument('--exp_name', type=str,
                        help='override the experiment name from config file. Optional.', required=False)
    parser.add_argument('--output_dir', type=str,
                        help='override the output directory from config file. Optional.', required=False)
    parser.add_argument('--plot_stress', type=bool,
                        help='override plot_stress from the config file. Optional.', required=False)
    parser.add_argument('--plot_grd_disp', type=bool,
                        help='override plot_grd_disp from the config file. Optional.', required=False)
    parser.add_argument('--rec_full_stress_tensor', type=str,
                        help='override rec_full_stress_tensor from the config file. Optional.', required=False)
    parser.add_argument('--gps_disp_points', type=str,
                        help='override gps_disp_points from the config file. Optional.', required=False)
    parser.add_argument('--strain_file', type=str,
                        help='override strain_file from the config file. Optional.', required=False)
    parser.add_argument('--save_file_type', type=str,
                        help='override save_file_type from the config file. Optional.', required=False)
    args = parser.parse_args()
    return args


def drive_calculation(args):
    params = PyCoulomb.configure_calc.configure_stress_calculation(config_file=args.config, cli_opts=args)
    [inputs, obs_disp_points, obs_strain_points] = PyCoulomb.input_values.read_inputs(params)
    out_object = PyCoulomb.run_dc3d.do_stress_computation(params, inputs, obs_disp_points, obs_strain_points)
    PyCoulomb.output_manager.produce_outputs(params, inputs, obs_disp_points, obs_strain_points, out_object)
    return


def main():
    my_args = welcome_and_parse_runstring()
    drive_calculation(my_args)
    return


if __name__ == "__main__":
    main()
