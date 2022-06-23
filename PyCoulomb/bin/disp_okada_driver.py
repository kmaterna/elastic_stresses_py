#!/usr/bin/env python
# An example pseudo-driver for elastic calculation. 
# To forward-model displacements from a slip distribution at known locations

from PyCoulomb.fault_slip_object import io_slippy
from PyCoulomb.fault_slip_object import fault_slip_object as fso
from PyCoulomb import run_dc3d, configure_calc, output_manager, io_additionals

filedict = {
    'datafile': 'lon_lat_points.txt',
    'fault_model': 'predicted_slip.txt',
    'lon0_sys': -115.50,
    'lat0_sys': 32.750,
    'bbox': [-115.65, -115.40, 32.65, 32.85] }


if __name__ == "__main__":
    # Input stage: Getting things into proper formats
    fault_model_list = io_slippy.read_slippy_distribution(filedict['fault_model']);
    coulomb_fault_model = fso.fault_dict_to_coulomb_fault(fault_model_list, filedict['lon0_sys'], filedict['lat0_sys']);
    disp_points = io_additionals.read_disp_points(filedict['datafile']);
    params = configure_calc.configure_default_displacement_params();
    inputs = configure_calc.configure_default_displacement_input(coulomb_fault_model, filedict['lon0_sys'],
                                                                 filedict['lat0_sys'], filedict['bbox']);
    # Compute and Output
    outobj = run_dc3d.do_stress_computation(params, inputs, disp_points=disp_points, strain_points=[]);
    output_manager.produce_outputs(params, inputs, disp_points, [], outobj);
