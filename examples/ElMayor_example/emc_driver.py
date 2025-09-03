#!/usr/bin/env python
# Run this in your elastic_py environment

import elastic_stresses_py.PyCoulomb.fault_slip_object as fso
from elastic_stresses_py.PyCoulomb.fault_slip_object.file_io import outputs
import tectonic_utils.seismo.moment_calculations as seismo_mo

filedict = {"slip_file": "./fialko_emc_dlc_format.txt"}


if __name__ == "__main__":
    fault_dict_list = fso.file_io.io_other.read_fialko_dlc(filedict["slip_file"])  # read into memory
    print("Moment Magnitude: ",
          seismo_mo.mw_from_moment(fso.fault_slip_object.get_total_moment(fault_dict_list)))
    fso.plot_fault_slip.map_source_slip_distribution(fault_dict_list, "fault_slip_fialko.png",
                                                     slip_cbar_opts=(0, 5, 0.5))
    outputs.write_gmt_fault_file(fault_dict_list, "fault_patches_gmt_format.txt")
