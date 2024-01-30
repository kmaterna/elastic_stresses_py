import unittest
import Elastic_stresses_py.PyCoulomb as PyCoulomb
from Elastic_stresses_py.PyCoulomb import fault_slip_triangle as fst
from Elastic_stresses_py.PyCoulomb import fault_slip_object as fso
import numpy as np


class Tests(unittest.TestCase):

    def test_okada_rect_vs_tri(self):
        obs_disp_points = PyCoulomb.io_additionals.read_disp_points("test/example_disp_field.txt")
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
        modeled_tri_points = fst.triangle_okada.compute_ll_def_tris(inputs, test_params, obs_disp_points)
        modeled_rect_points = PyCoulomb.run_dc3d.compute_ll_def(inputs, test_params, obs_disp_points)

        # The comparison part
        obs_E_triangles = np.array([x.dE_obs for x in modeled_tri_points])
        obs_E_rectangles = np.array([x.dE_obs for x in modeled_rect_points])
        np.testing.assert_allclose(obs_E_rectangles, obs_E_triangles, atol=0.0001)
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


if __name__ == "__main__":
    unittest.main()
