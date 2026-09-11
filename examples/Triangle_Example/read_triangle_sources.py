#!/usr/bin/env python

"""
An example script to read and compute triangles from R. Lohman
9/10/2026 K. Materna
"""

from elastic_stresses_py.PyCoulomb.fault_slip_triangle.file_io import io_other
from elastic_stresses_py.PyCoulomb.fault_slip_object import plot_fault_slip
from elastic_stresses_py.PyCoulomb.disp_points_object import disp_points_object
from elastic_stresses_py.PyCoulomb import configure_calc, inputs_object, run_dc3d, output_manager

filedict = {"mesh": "Data/forKatherine.mat",
            "slip_values": "Data/2005_2015inversion_fault_vertices_latlondepthslip.txt",
            'bbox': [-115.75, -115.50, 33.00, 33.30]
            }


def read_slip_values(filename):
    """ Read the slip values from this file, simple function """
    slip_values = []
    with open(filename) as ifile:
        count = 0
        for line in ifile:
            if count == 0:
                slip_values.append(float(line.split()[-1]))
            count += 1
            if line.split()[0] == '>':
                count = 0
    return slip_values


def do_main():
    """ Main function """
    """ Part 1: Read triangles from MAT file and assign slip assuming the MAT is read in same order as slip file """
    fault_slip_triangles = io_other.read_brawley_lohman_2005(filedict['mesh'])  # turns out my function already works
    slip_values = read_slip_values(filedict['slip_values'])
    print(slip_values)
    for i, (item, slip) in enumerate(zip(fault_slip_triangles, slip_values)):
        fault_slip_triangles[i] = item.change_fault_slip(rtlat=-slip)  # change each fault to have associated slip
        # note: rtlat = -slip because the fault is assumed to slip left-laterally.  Change if incorrect.
    lon0_sys = fault_slip_triangles[0].lon
    lat0_sys = fault_slip_triangles[0].lat  # generate the coordinate system for the whole calculation

    """ Part 2: Use triangular mesh source to produce displacements (strains also possible). Calculated at GNSS."""
    params = configure_calc.Params()  # configure with default values
    p507 = disp_points_object.Displacement_points(lon=-115.612, lat=33.200, name='P507')  # from https://geodesy.unr.edu/NGLStationPages/stations/P507.sta
    p506 = disp_points_object.Displacement_points(lon=-115.510, lat=33.081, name='P506')  # from UNR magnet
    disp_points = [p507, p506]  # represents each GNSS station

    # # Compute Displacements and Outputs with the normal higher-level API
    inputs = inputs_object.input_obj.configure_default_displacement_input(source_object=fault_slip_triangles,
                                                                          zerolon=lon0_sys,
                                                                          zerolat=lat0_sys,
                                                                          bbox=filedict['bbox'])
    outobj = run_dc3d.do_stress_computation(params, inputs, disp_points=disp_points, strain_points=())
    output_manager.produce_outputs(params, inputs, obs_disp_points=disp_points,
                                   obs_strain_points=(), out_object=outobj)

    # Extra plot for prettiness
    plot_fault_slip.map_source_slip_distribution(fault_slip_triangles, 'model_prediction.png',
                                                 disp_points=outobj.model_disp_points,
                                                 slip_cbar_opts=(0, .5, 0.05),
                                                 scale_arrow=(1.0, 0.002, "2 mm"),
                                                 region=filedict['bbox'])  # mapping total slip distribution
    return


if __name__ == "__main__":
    do_main()
