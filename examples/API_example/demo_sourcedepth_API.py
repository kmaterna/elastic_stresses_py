#!/usr/bin/env python

from Elastic_stresses_py import PyCoulomb
from Elastic_stresses_py.PyCoulomb import coulomb_collections as cc

# My configuration: special experiment, focal mechanism as a function of depth.  This is my "lab notebook"
config = {
    "exp1": [1.0, "depth1.0"],
    "exp2": [2.0, "depth2.0"],
    "exp3": [3.0, "depth3.0"],
    "exp4": [4.0, "depth4.0"],
    "exp5": [5.0, "depth5.0"],
    "exp6": [6.0, "depth6.0"],
    "exp7": [7.0, "depth7.0"]
}

def swap_output_dir(default_params, dirname):
    # write lines to replace params.outdir with the proper params.outdir
    # PACK MY SUITCASE
    MyParams = PyCoulomb.configure_calc.modify_params_object(default_params, outdir=dirname);
    return MyParams;


def swap_input_depth(default, new_depth):
    # PACK MY SUITCASE. Replace source(old_depth) with source(new_depth)
    new_sources = [];
    inputs = cc.Input_object(PR1=default.PR1, FRIC=default.FRIC, depth=default.depth,
                             start_gridx=default.start_gridx, finish_gridx=default.finish_gridx,
                             start_gridy=default.start_gridy, finish_gridy=default.finish_gridy,
                             xinc=default.xinc, yinc=default.yinc,
                             minlon=default.minlon, maxlon=default.minlon, zerolon=default.zerolon,
                             minlat=default.minlat, maxlat=default.maxlat, zerolat=default.zerolat,
                             source_object=new_sources, receiver_object=default.receiver_object,
                             receiver_horiz_profile=default.receiver_object);
    return inputs;


if __name__ == "__main__":   # CONFIGURE, INPUT, COMPUTE, OUTPUT
    default_params = PyCoulomb.configure_calc.configure_stress_calculation("myconfig.txt");
    [default_inputs, obs_disp_points, _] = PyCoulomb.input_values.read_inputs(default_params);
    for expkey in config.keys():
        exp_config = swap_output_dir(default_params, config[expkey][1])   # replace default outdir with special outdir
        inputs = swap_input_depth(default_inputs, config[expkey][0])
        out_object = PyCoulomb.run_dc3d.do_stress_computation(exp_config, inputs, obs_disp_points, ());
        PyCoulomb.output_manager.produce_outputs(exp_config, inputs, obs_disp_points, (), out_object);
