"""
Read the Coulomb or user-defined input files:

* .inp (Coulomb)
* .inr (Coulomb)
* .intxt (My own definition, convenient for slip or magnitude)
* .inzero (My own definition, convenient for point sources)
"""

from .inputs_object import io_inp, io_inr, io_intxt
from . import io_additionals, utilities
from .disp_points_object.disp_points_object import Displacement_points


def read_inputs(params):
    """Driver to read any format of source and receiver faults for calculation."""
    if '.inp' in params.input_file:
        input_object = io_inp.read_inp(params.input_file, params.fixed_rake)  # fixed rake format
    elif '.inr' in params.input_file:
        input_object = io_inr.read_inr(params.input_file)  # variable rake format (will write later)
    elif '.intxt' in params.input_file:
        input_object = io_intxt.read_intxt(params.input_file, params.mu, params.lame1)  # convenient input format
    elif '.inzero' in params.input_file:
        input_object = io_intxt.read_intxt(params.input_file, params.mu, params.lame1)  # point source
    else:
        raise Exception("Error! Unrecognized type of input file!")

    # Read points where displacement or strain will be written out
    disp_points, strain_points = (), ()
    if params.disp_points_file:
        disp_points = io_additionals.read_disp_points(params.disp_points_file)
    if params.strain_file:
        strain_points = io_additionals.read_disp_points(params.strain_file)
    if isinstance(strain_points, Displacement_points):
        strain_points = [strain_points]
    if isinstance(disp_points, Displacement_points):
        disp_points = [disp_points]

    assert input_object.source_object, ValueError("You have not specified any sources.")
    assert input_object.xinc > 0, ValueError("Your cartesian coordinate system's x-increment defaults to zero.")
    assert input_object.yinc > 0, ValueError("Your cartesian coordinate system's y-increment defaults to zero.")
    utilities.print_metrics_on_sources(input_object.source_object, params.mu)
    input_object.check_each_fault_has_same_coord_system()
    return [input_object, disp_points, strain_points]
