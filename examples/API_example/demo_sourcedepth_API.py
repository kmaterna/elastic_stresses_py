#!/usr/bin/env python

# My configuration: special experiment, focal mechanism as a function of depth
config = {
	"exp1": [1.0, "depth1.0"], 
	"exp2": [2.0, "depth2.0"], 
	"exp3": [3.0, "depth3.0"], 
	"exp4": [4.0, "depth4.0"], 
	"exp5": [5.0, "depth5.0"], 
	"exp6": [6.0, "depth6.0"], 
	"exp7": [7.0, "depth7.0"]
}

def swap_output_dir(default, dirname):
	# write lines to replace params.outdir with the proper params.outdir
	# PACK MY SUITCASE
    MyParams = cc.Params(config_file=default.config_file, input_file=default.input_file, aftershocks=default.aftershocks,
                         disp_points_file=default.gps_file, strain_file=default.strain_file,
                         strike_num_receivers=default.strike_num_receivers, fixed_rake=default.fixed_rake,
                         dip_num_receivers=default.dip_num_receivers, mu=default.mu, lame1=default.lame1, B=default.B,
                         alpha=default.alpha, plot_stress=default.plot_stress, plot_grd_disp=default.plot_grd_disp,
                         outdir=DIRNAME);
	return better_configured_params;

def swap_input_depth(default, new_depth):
	# PACK MY SUITCASE. Replace source(old_depth) with source(new_depth)
	return;


if __name__ == "__main__":   # CONFIGURE, INPUT, COMPUTE, OUTPUT
    default_params = PyCoulomb.configure_calc.configure_stress_calculation("myconfig.txt");    
    [inputs, obs_disp_points, _] = PyCoulomb.input_values.read_inputs(default_params);
    for expkey in config.keys():
    	exp_config = swap_output_dir(default_params, config[expkey][1])   # replace default outputdir with special outputdir
    	inputs = swap_depth(default_inputs, config[expkey][0])
    	out_object = PyCoulomb.run_dc3d.do_stress_computation(params, inputs, obs_disp_points, ());
    	PyCoulomb.output_manager.produce_outputs(params, inputs, obs_disp_points, (), out_object);



