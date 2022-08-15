#!/usr/bin/env python

from Elastic_stresses_py import PyCoulomb

# My configuration: special experiment, focal mechanism as a function of depth.  This is my "lab notebook".
# First column: depth.  Second column: Directory name.
config = {
    "exp1": [1.0, "depth1.0/"],
    "exp2": [2.0, "depth2.0/"],
    "exp2.5": [2.5, "depth2.5/"],
    "exp3": [3.0, "depth3.0/"],
    "exp4": [4.0, "depth4.0/"],
    "exp5": [5.0, "depth5.0/"],
    "exp6": [6.0, "depth6.0/"],
    "exp7": [7.0, "depth7.0/"]
}

def swap_output_dir(default_params, dirname):
    # write lines to replace params.outdir with the proper params.outdir
    # PACK MY SUITCASE
    old_outdir = '/'.join(default_params.outdir.split('/')[0:-2]);
    MyParams = PyCoulomb.configure_calc.modify_params_object(default_params, outdir=old_outdir+"/"+dirname);
    return MyParams;


def swap_input_depth(default_inputs, new_depth):
    # PACK MY SUITCASE. Replace source(old_depth) with source(new_depth)
    base_source = default_inputs.source_object[0];
    new_point_source = PyCoulomb.configure_calc.modify_fault_object(base_source, top=new_depth, bottom=new_depth);
    inputs = PyCoulomb.configure_calc.modify_inputs_object(default_inputs, source_object=[new_point_source]);
    return inputs;


if __name__ == "__main__":   # CONFIGURE, INPUT, COMPUTE, OUTPUT
    default_params = PyCoulomb.configure_calc.configure_stress_calculation("myconfig.txt");
    [default_inputs, obs_disp_points, _] = PyCoulomb.input_values.read_inputs(default_params);
    for expkey in config.keys():
        exp_config = swap_output_dir(default_params, config[expkey][1])   # replace default outdir with special outdir
        inputs = swap_input_depth(default_inputs, config[expkey][0])
        out_object = PyCoulomb.run_dc3d.do_stress_computation(exp_config, inputs, obs_disp_points, ());
        PyCoulomb.output_manager.produce_outputs(exp_config, inputs, obs_disp_points, (), out_object);
