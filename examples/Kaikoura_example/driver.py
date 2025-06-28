#!/usr/bin/env python
# Run this in your pygmt environment

import elastic_stresses_py.PyCoulomb.fault_slip_object as fso
import tectonic_utils.seismo.moment_calculations as seismo_mo
import elastic_stresses_py.PyCoulomb as PyCoulomb

filedict = {"usgs_slip_file": "../../_Data/files_MTMOD_WS/Kaikoura_usgs_finite_fault.fsp",
            "hamling": "../../_Data/files_MTMOD_WS/hamling_coseismic_model/hamling_coseismic_crustal_fault_only.dat",
            "wallace_file": "../../_Data/files_MTMOD_WS/W18_SSE_time_evolution_total_slip/kuni_time_002.atr",
            "demo_config": "demo_config.txt",   # contains the setup of the fault system
            "finished_input": "usgs_inputs.intxt"}


def write_pycoulomb_input(filedict, fault_dict_list):
    params = PyCoulomb.configure_calc.configure_stress_calculation(filedict["demo_config"])
    demo_obj = PyCoulomb.inputs_object.io_intxt.read_intxt(params.input_file, params.mu, params.lame1)
    pycoulomb_faults = fso.fault_slip_object.fault_object_to_coulomb_fault(fault_dict_list, demo_obj.zerolon,
                                                                           demo_obj.zerolat)
    inputs = demo_obj.modify_inputs_object(source_object=pycoulomb_faults)
    PyCoulomb.inputs_object.io_intxt.write_intxt(inputs, filedict["finished_input"], mu=params.mu, lame1=params.lame1)
    return


if __name__ == "__main__":
    # TO READ SRCMOD AND WRITE PYCOULOMB INPUTS
    # fault_dict_list = fso.io_other.io_hamling_2017(filedict["hamling"])  # examples of other read functions
    # fault_dict_list = fso.io_other.io_wallace_sse(filedict["wallace_file2"])  # examples of other read functions
    fault_dict_list = fso.file_io.io_srcmod.read_srcmod_distribution(filedict["usgs_slip_file"])  # read into memory
    print("Moment Magnitude: ",
          seismo_mo.mw_from_moment(fso.fault_slip_object.get_total_moment(fault_dict_list)))
    fso.plot_fault_slip.map_source_slip_distribution(fault_dict_list, "fault_slip_usgs.png")

    # # Set up calculation
    write_pycoulomb_input(filedict, fault_dict_list)

    # NOW GO RUN "elastic_stresses_driver.py real_config.txt" (from the terminal)
