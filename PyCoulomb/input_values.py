"""
Read the Coulomb or user-defined input files:

* .inp (Coulomb)
* .inr (Coulomb)
* .intxt (My own definition, convenient for slip or magnitude)
* .inzero (My own definition, convenient for point sources)
"""

from . import io_inp
from . import io_inr
from . import io_intxt
from . import io_additionals
from . import utilities


def read_inputs(params):
    # Read source and receiver faults for calculation
    if '.inp' in params.input_file:
        input_object = io_inp.read_inp(params.input_file, params.fixed_rake);  # fixed rake format
    elif '.inr' in params.input_file:
        input_object = io_inr.read_inr(params.input_file);  # variable rake format (will write later);
    elif '.intxt' in params.input_file:
        input_object = io_intxt.read_intxt(params.input_file);  # convenient input format
    elif '.inzero' in params.input_file:
        input_object = io_intxt.read_intxt(params.input_file);  # a point source in convenient format
    else:
        raise Exception("Error! Unrecognized type of input file!");

    # Read points where displacement or strain will be written out
    disp_points, strain_points = (), ();
    if params.disp_points_file:
        disp_points = io_additionals.read_disp_points(params.disp_points_file);
    if params.strain_file:
        strain_points = io_additionals.read_disp_points(params.strain_file);

    assert input_object.source_object, ValueError("You have not specified any sources.");
    utilities.print_metrics_on_sources(input_object.source_object);
    return [input_object, disp_points, strain_points];
