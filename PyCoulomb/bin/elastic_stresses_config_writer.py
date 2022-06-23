#!/usr/bin/env python

# A simple script to write a valid config template in a particular directory
# Call this from a new experiment directory to get a valid config file template

import argparse
import Elastic_stresses_py.PyCoulomb as PyCoulomb

def welcome_and_parse_runstring():
    print("\nWrite an empty config file for elastic_stresses_py. ");
    parser = argparse.ArgumentParser(description='Write an empty config file into a specified directory',
                                     epilog='\U0001f600 \U0001f600 \U0001f600 ');
    parser.add_argument('directory', type=str, help='name of directory. Required.')
    args = parser.parse_args()
    return args;


if __name__ == "__main__":
    args = welcome_and_parse_runstring();
    PyCoulomb.configure_calc.write_valid_config_file(args.directory);
