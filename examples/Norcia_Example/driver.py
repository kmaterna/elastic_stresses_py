#!/usr/bin/env python
# Run this in your pygmt environment

import Elastic_stresses_py.PyCoulomb.fault_slip_object as fso
import Tectonic_Utils.seismo as seismo
import Elastic_stresses_py.PyCoulomb as PyCoulomb

filedict = {"fault_slip_file": "s2016NORCIA01PIZZ.fsp",
            "demo_config": "demo_config.txt",
            "real_config": "real_config.txt",
            "norcia_input": "norcia_inputs.intxt"}


def replace_inputs(demo_obj, fault_dict_list):
    pycoulomb_faults = fso.fault_slip_object.fault_object_to_coulomb_fault(fault_dict_list, demo_obj.zerolon,
                                                                           demo_obj.zerolat)
    input_object = PyCoulomb.configure_calc.modify_inputs_object(demo_obj, source_object=pycoulomb_faults);
    return input_object;


def convert_norcia_to_pycoulomb(filedict):
    # DO ONCE: READ SRCMOD INTO FAULT_SLIP_OBBJECT, WRITE PYCOULOMB INTXT
    italy_faults = fso.file_io.io_srcmod.read_srcmod_distribution(filedict["fault_slip_file"]);
    fso.plot_fault_slip.map_source_slip_distribution(italy_faults, "fault_slip.png", region=[12, 14, 41.7, 43.2]);
    Mw = seismo.moment_calculations.mw_from_moment(fso.fault_slip_object.get_total_moment(italy_faults))
    print("Moment Magnitude: ", Mw);
    params = PyCoulomb.configure_calc.configure_stress_calculation(filedict["demo_config"]);
    inputs = PyCoulomb.io_intxt.read_intxt(params.input_file, params.mu, params.lame1);  # convenient input format
    inputs = replace_inputs(inputs, italy_faults);
    PyCoulomb.io_intxt.write_intxt(inputs, filedict["norcia_input"], mu=params.mu, lame1=params.lame1);
    return;


if __name__ == "__main__":

    convert_norcia_to_pycoulomb(filedict);

    # Steps of elastic_stresses_driver
    params = PyCoulomb.configure_calc.configure_stress_calculation(filedict["config"]);
    [inputs, obs_disp_points, _] = PyCoulomb.input_values.read_inputs(params);
    out_object = PyCoulomb.run_dc3d.do_stress_computation(params, inputs, obs_disp_points);
    PyCoulomb.output_manager.produce_outputs(params, inputs, obs_disp_points, (), out_object);
