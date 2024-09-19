import unittest
import elastic_stresses_py.PyCoulomb as PyCoulomb
from elastic_stresses_py.PyCoulomb import fault_slip_triangle as fst
from elastic_stresses_py.PyCoulomb import fault_slip_object as fso
from elastic_stresses_py.PyCoulomb import run_okada_wrapper, run_dc3d
from elastic_stresses_py.PyCoulomb.disp_points_object.disp_points_object import Displacement_points
import numpy as np


class Tests(unittest.TestCase):

    def test_okada_rect_vs_tri(self):
        obs_disp_points = PyCoulomb.io_additionals.read_disp_points("test/example_disp_field.txt")
        obs_disp_points = [i.with_depth_as(1) for i in obs_disp_points]  # can test different depths
        zerolon, zerolat = -115.68, 33.1012233
        test_params = PyCoulomb.configure_calc.Params(mu=30e9, lame1=30e9)
        test_rect = fso.fault_slip_object.FaultSlipObject(lon=zerolon + 0.01, lat=zerolat - 0.2, strike=50, dip=70,
                                                          depth=3, segment=0, length=15, width=5, rake=170, slip=-0.2)
        test_rect2 = fso.fault_slip_object.FaultSlipObject(lon=zerolon + 0.41, lat=zerolat - 0.4, strike=170, dip=70,
                                                           depth=3, segment=0, length=15, width=5, rake=170, slip=0.4,
                                                           tensile=0.3)

        pycoulomb_sources = fso.fault_slip_object.fault_object_to_coulomb_fault([test_rect, test_rect2],
                                                                                zerolon_system=zerolon,
                                                                                zerolat_system=zerolat)
        inputs = PyCoulomb.coulomb_collections.Input_object(zerolon=zerolon, zerolat=zerolat, depth=None,
                                                            maxlat=None, maxlon=None, minlat=None, minlon=None,
                                                            receiver_object=None,
                                                            source_object=pycoulomb_sources, xinc=None, yinc=None)

        # The heart of the test: same interface, different guts
        rect_points = run_okada_wrapper.compute_ll_def(inputs, test_params, obs_disp_points)  # okada_wrapper
        tri_points = run_dc3d.compute_ll_def(inputs, test_params, obs_disp_points)  # cutde version
        rect_strains = run_okada_wrapper.compute_ll_strain(inputs, test_params, obs_disp_points[0:3])
        tri_strains = run_dc3d.compute_ll_strain(inputs, test_params, obs_disp_points[0:3])

        # # The comparison part
        obs_E_triangles = np.array([x.dE_obs for x in tri_points])
        obs_E_rectangles = np.array([x.dE_obs for x in rect_points])
        np.testing.assert_allclose(obs_E_rectangles, obs_E_triangles, atol=0.0001)
        obs_N_triangles = np.array([x.dN_obs for x in tri_points])
        obs_N_rectangles = np.array([x.dN_obs for x in rect_points])
        np.testing.assert_allclose(obs_N_rectangles, obs_N_triangles, atol=0.0001)

        for rect_strain, tri_strain in zip(rect_strains, tri_strains):
            np.testing.assert_allclose(rect_strain[0][0], tri_strain[0][0], atol=0.0001)
            np.testing.assert_allclose(rect_strain[0][1], tri_strain[0][1], atol=0.0001)
            np.testing.assert_allclose(rect_strain[0][2], tri_strain[0][2], atol=0.0001)
            np.testing.assert_allclose(rect_strain[1][1], tri_strain[1][1], atol=0.0001)
            np.testing.assert_allclose(rect_strain[1][2], tri_strain[1][2], atol=0.0001)
            np.testing.assert_allclose(rect_strain[2][2], tri_strain[2][2], atol=0.0001)
        return

    def test_triangle_geometry_functions(self):
        """ Test the dip and strike calculation for a triangular fault element. """
        tri_fault = fst.fault_slip_triangle.TriangleFault(lon=-121.0, lat=40, depth=5000,
                                                          vertex1=np.array([0, 0, 5000]),
                                                          vertex2=np.array([0, 1000, 5000]),
                                                          vertex3=np.array([1000, 1000, 6000]))
        np.testing.assert_almost_equal(0, tri_fault.get_strike())
        np.testing.assert_almost_equal(45, tri_fault.get_dip())
        return

    def test_self_stress(self):
        """
        Test the stress on its own patch, when source and receiver are the same patch.
        Makes the self-stress pass under triangles.
        """
        zerolon, zerolat = -115.68, 33.1012233
        test_params = PyCoulomb.configure_calc.Params(mu=30e9, lame1=30e9)
        test_rect = fso.fault_slip_object.FaultSlipObject(lon=-115, lat=33, strike=50, dip=70,
                                                          depth=3, length=15, width=5, rake=170, slip=-0.2)
        # receiver = fso.fault_slip_object.FaultSlipObject(lon=-115, lat=33.01, strike=50, dip=70,
        #                                                  depth=3, length=15, width=5, rake=170, slip=0)
        # WHEN THE RECEIVER IS THE SAME FAULT, DOES IT SUCCEED?
        receiver = test_rect.change_fault_slip(new_slip=0)
        pycoulomb_sources = fso.fault_slip_object.fault_object_to_coulomb_fault([test_rect],
                                                                                zerolon_system=zerolon,
                                                                                zerolat_system=zerolat)
        pycoulomb_recs = fso.fault_slip_object.fault_object_to_coulomb_fault([receiver],
                                                                             zerolon_system=zerolon,
                                                                             zerolat_system=zerolat)
        inputs = PyCoulomb.coulomb_collections.Input_object(zerolon=zerolon, zerolat=zerolat, depth=None,
                                                            maxlat=None, maxlon=None, minlat=None, minlon=None,
                                                            receiver_object=pycoulomb_recs,
                                                            source_object=pycoulomb_sources, xinc=None, yinc=None)

        # The heart of the test: same interface, different guts
        centercoords = inputs.receiver_object[0].get_fault_center()  # in cartesian coordinates
        strain_point = [Displacement_points(lon=centercoords[0], lat=centercoords[1], depth=centercoords[2])]
        rect_strains = run_okada_wrapper.compute_strain_point(inputs.source_object, 2/3, x=centercoords[0],
                                                              y=centercoords[1], compute_depth=centercoords[2])  # o_w
        tri_strains = run_dc3d.compute_xy_strain(inputs, test_params, strain_point)[0]  # cutde_version

        print("tri_strains:", tri_strains)
        print("rect_strains:", rect_strains)

        np.testing.assert_almost_equal(tri_strains[0][0], rect_strains[0][0])
        np.testing.assert_almost_equal(tri_strains[0][1], rect_strains[0][1])
        np.testing.assert_almost_equal(tri_strains[2][1], rect_strains[2][1])
        return


if __name__ == "__main__":
    unittest.main()
