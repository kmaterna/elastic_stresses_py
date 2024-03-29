import unittest
import elastic_stresses_py.PyCoulomb as PyCoulomb
from elastic_stresses_py.PyCoulomb.pyc_fault_object import Faults_object
import numpy as np
import elastic_stresses_py.PyCoulomb.inputs_object.io_mt as io_mt
from elastic_stresses_py.PyCoulomb.point_source_object import point_sources
from elastic_stresses_py.PyCoulomb import run_okada_wrapper


class Tests(unittest.TestCase):

    def test_pt_source(self):
        obs_disp_points = PyCoulomb.io_additionals.read_disp_points("test/example_disp_field.txt")
        obs_disp_points = [i.with_depth_as(4) for i in obs_disp_points]
        zerolon, zerolat = -115.68, 33.1012233
        depth = 10
        dip = 80.9
        strike = 135
        rake = 70
        test_params = PyCoulomb.configure_calc.Params(mu=30e9, lame1=30e9)
        [x, y, potency, comment] = io_mt.compute_params_for_point_source(rake=rake, magnitude=7.5,
                                                                         lon=-115.4, lat=32.8, zerolon=zerolon,
                                                                         zerolat=zerolat, mu=30e9)
        one_source_object = Faults_object(xstart=x, xfinish=x, ystart=y, yfinish=y, rtlat=0, reverse=0, potency=potency,
                                          strike=strike, dipangle=dip, zerolon=zerolon, zerolat=zerolat, rake=rake,
                                          top=depth, bottom=depth, comment=comment)
        one_tensile = Faults_object(xstart=x, xfinish=x, ystart=y, yfinish=y, rtlat=0, reverse=0,
                                    potency=[0, 0, 1e8, 0], strike=strike, dipangle=dip, zerolon=zerolon,
                                    zerolat=zerolat, rake=rake, top=depth, bottom=depth, comment=comment)
        inputs = PyCoulomb.coulomb_collections.Input_object(zerolon=zerolon, zerolat=zerolat, depth=None,
                                                            maxlat=None, maxlon=None, minlat=None, minlon=None,
                                                            receiver_object=None, source_object=[one_source_object],
                                                            xinc=None, yinc=None)

        # Try a single point source, tensile component
        z = depth
        ow_strain, ow_disp = run_okada_wrapper.compute_strains_stresses_from_one_fault(one_tensile, 5, -3, z, 0.6667)
        ps_strain, ps_disp = point_sources.compute_displacements_strains_point(one_tensile, 5, -3, z, 0.6667)
        np.testing.assert_allclose(ow_disp, ps_disp, atol=1e-5)
        np.testing.assert_allclose(ow_strain, ps_strain, atol=1e-7)

        # The heart of the test: same interface, different guts
        modeled_ow_points = PyCoulomb.run_okada_wrapper.compute_ll_def(inputs, test_params, obs_disp_points)  # o_w
        modeled_new_points = PyCoulomb.run_dc3d.compute_ll_def(inputs, test_params, obs_disp_points)  # chatgpt

        # The comparison part: compare displacements
        obs_E_ow = np.array([x.dE_obs for x in modeled_ow_points])
        obs_E_py = np.array([x.dE_obs for x in modeled_new_points])
        np.testing.assert_allclose(obs_E_ow, obs_E_py, atol=0.0001)
        obs_N_ow = np.array([x.dN_obs for x in modeled_ow_points])
        obs_N_new = np.array([x.dN_obs for x in modeled_new_points])
        np.testing.assert_allclose(obs_N_ow, obs_N_new, atol=0.0001)

        # The strain tensor comparison: compare strains
        modeled_ow_tensors = PyCoulomb.run_okada_wrapper.compute_ll_strain(inputs, test_params, obs_disp_points)  # o_w
        modeled_new_tensors = PyCoulomb.run_dc3d.compute_ll_strain(inputs, test_params, obs_disp_points)  # chatgpt

        for i in range(len(modeled_ow_tensors)):
            # print("okada_wrapper:")
            # print(modeled_ow_tensors[i])
            # print("new translation:")
            # print(modeled_new_tensors[i])
            # print("\n\n")
            np.testing.assert_allclose(modeled_ow_tensors[i], modeled_new_tensors[i], atol=1e-7)

        # Zooming in on the closest-function issue
        ow_strain, ow_disp = run_okada_wrapper.compute_strains_stresses_from_one_fault(one_source_object, 6, -8, 5, alpha=2 / 3)
        ps_strain, ps_disp = point_sources.compute_displacements_strains_point(one_source_object, 6, -8, 5, alpha=2/3)
        print("Do the inner functions produce the same thing?")
        print(ps_strain)
        print(ow_strain)
        np.testing.assert_allclose(ps_strain, ow_strain, atol=1e-7)
        return


if __name__ == "__main__":
    unittest.main()
