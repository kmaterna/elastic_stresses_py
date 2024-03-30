#!/usr/bin/env python
# An example pseudo-driver for elastic calculation. 
# To forward-model displacements from a slip distribution at known locations

from elastic_stresses_py.PyCoulomb import fault_slip_object as fso
from elastic_stresses_py.PyCoulomb import run_dc3d, configure_calc, output_manager, io_additionals, inputs_object

filedict = {
    'datafile': 'lon_lat_points.txt',
    'fault_model': 'predicted_slip.txt',
    'lon0_sys': -115.50,
    'lat0_sys': 32.750,
    'bbox': [-115.65, -115.40, 32.65, 32.85]}


if __name__ == "__main__":
    # Input stage: Getting things into proper formats
    fault_model_list = fso.file_io.io_slippy.read_slippy_distribution(filedict['fault_model'])
    coulomb_fault_model = fso.fault_slip_object.fault_object_to_coulomb_fault(fault_model_list, filedict['lon0_sys'],
                                                                              filedict['lat0_sys'])
    disp_points = io_additionals.read_disp_points(filedict['datafile'])
    params = configure_calc.Params()  # configure with default values
    inputs = inputs_object.input_obj.configure_default_displacement_input(coulomb_fault_model,
                                                                          zerolon=filedict['lon0_sys'],
                                                                          zerolat=filedict['lat0_sys'],
                                                                          bbox=filedict['bbox'])
    # Compute and Output
    outobj = run_dc3d.do_stress_computation(params, inputs, disp_points=disp_points, strain_points=())
    output_manager.produce_outputs(params, inputs, obs_disp_points=disp_points,
                                   obs_strain_points=(), out_object=outobj)
