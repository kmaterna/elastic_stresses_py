import unittest
import Elastic_stresses_py.PyCoulomb as PyCoulomb
from Elastic_stresses_py.PyCoulomb import fault_slip_triangle as fst
from Elastic_stresses_py.PyCoulomb import fault_slip_object as fso
import numpy as np

class Tests(unittest.TestCase):

    def test_okada_rect_vs_tri(self):
        obs_disp_points = PyCoulomb.io_additionals.read_disp_points("test/example_disp_field.txt");
        zerolon, zerolat = -115.68, 33.1012233;
        test_params = PyCoulomb.configure_calc.configure_default_displacement_params(mu=30e9, lame1=30e9);
        _, pr = PyCoulomb.conversion_math.get_poissons_ratio_and_alpha(mu=30e9, lame1=30e9)
        pr = 0.25;

        test_rectangular_fault = fso.fault_slip_object.FaultSlipObject(lon=zerolon, lat=zerolat, strike=50,
                                                                       dip=70, depth=3, segment=0, length=15, width=5,
                                                                       rake=170, slip=-0.2, tensile=0);
        [pycoulomb_rectangle] = fso.fault_slip_object.fault_object_to_coulomb_fault([test_rectangular_fault]);
        inputs = PyCoulomb.coulomb_collections.Input_object(zerolon=zerolon, zerolat=zerolat, finish_gridy=None,
                                                            finish_gridx=None, depth=None, FRIC=None, PR1=None,
                                                            maxlat=None, maxlon=None, minlat=None, minlon=None,
                                                            receiver_object=None, receiver_horiz_profile=None,
                                                            source_object=[pycoulomb_rectangle],
                                                            start_gridy=None, start_gridx=None, xinc=None, yinc=None);

        tri_faults = fst.fault_slip_triangle.convert_rectangle_into_two_triangles(test_rectangular_fault);
        modeled_tri_points = fst.triangle_okada.compute_disp_points_from_triangles(tri_faults, obs_disp_points, pr);
        modeled_rect_points = PyCoulomb.run_dc3d.compute_ll_def(inputs, test_params, obs_disp_points);

        # # Optional: Visualizing the faults
        # fst.fault_slip_triangle.write_gmt_plots_geographic(tri_faults, "test/temp-outfile.txt");
        # fso.plot_fault_slip.map_source_slip_distribution([test_rectangular_fault], 'test/model_rect_disps.png',
        #                                                  disp_points=modeled_rect_points);
        # fso.plot_fault_slip.map_source_slip_distribution([test_rectangular_fault], 'test/model_tri_disps.png',
        #                                                  disp_points=modeled_tri_points, plot_slip_colorbar=False,
        #                                                  fault_traces_from_file="test/temp-outfile.txt");

        obs_E_triangles = np.array([x.dE_obs for x in modeled_tri_points]);
        obs_E_rectangles = np.array([x.dE_obs for x in modeled_rect_points]);
        np.testing.assert_allclose(obs_E_rectangles, obs_E_triangles, atol=0.0001);
        return;


if __name__ == "__main__":
    unittest.main();
