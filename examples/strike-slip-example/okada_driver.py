#!/usr/bin/env python

import Elastic_stresses_py.PyCoulomb.fault_slip_object as fso
from Elastic_stresses_py.PyCoulomb import run_dc3d, configure_calc, output_manager, io_additionals, inputs_object


# Definitions
lon0_sys, lat0_sys = -120.5, 36
bbox = (-121.5, -119.5, 35.2, 36.8)
lonlatfile = "Inputs/lon_lats.txt"
source_slip_dist = "Inputs/s2004PARKFI01CUST.fsp"

# Inputs
parkfield_faults = fso.file_io.io_srcmod.read_srcmod_distribution(source_slip_dist)
coulomb_fault_model = fso.fault_slip_object.fault_object_to_coulomb_fault(parkfield_faults, lon0_sys, lat0_sys)
disp_points = io_additionals.read_disp_points(lonlatfile)

# Configure, Compute, Output
params = configure_calc.Params()  # configure with default values
inputs = inputs_object.input_obj.configure_default_displacement_input(coulomb_fault_model, zerolon=lon0_sys,
                                                                      zerolat=lat0_sys, bbox=bbox, domainsize=100)
outobj = run_dc3d.do_stress_computation(params, inputs, disp_points=disp_points, strain_points=[])
output_manager.produce_outputs(params, inputs, disp_points, obs_strain_points=[], out_object=outobj)
