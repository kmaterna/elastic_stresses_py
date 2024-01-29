# Test the conversion between four different formats for slip distributions and fault geometry.

import unittest
import Elastic_stresses_py.PyCoulomb.fault_slip_object as fso


example_fault = fso.fault_slip_object.FaultSlipObject(strike=5, dip=75, length=40, width=20, lon=-123.00, lat=40.00,
                                                      depth=3, rake=10, slip=1, tensile=0, segment=0)


class Tests(unittest.TestCase):

    def test_io_geojson(self):
        flush_file = "test/example_geojson.txt"
        fso.file_io.io_geojson.write_faults_json([example_fault], flush_file)
        [returned_fault] = fso.file_io.io_geojson.read_faults_json(flush_file)
        for a, b in zip(example_fault.__iter__(), returned_fault.__iter__()):
            # don't test the lossy part of this conversion (rake, slip, tensile)
            if a[0] in ['lon', 'lat']:
                self.assertAlmostEqual(a[1], b[1], places=3)
            elif a[0] in ['strike', 'dip', 'length', 'width', 'depth']:
                self.assertAlmostEqual(a[1], b[1])
        return

    def test_io_pycoulomb(self):
        pycoulomb_faults = fso.fault_slip_object.fault_object_to_coulomb_fault([example_fault])
        [returned_fault] = fso.fault_slip_object.coulomb_fault_to_fault_object(pycoulomb_faults)
        for a, b in zip(example_fault.__iter__(), returned_fault.__iter__()):
            self.assertAlmostEqual(a[1], b[1])   # print(a[1], b[1])
        return

    def test_io_slippy(self):
        flush_file = "test/example_slippy.txt"
        fso.file_io.io_slippy.write_slippy_distribution([example_fault, example_fault], flush_file)
        [returned_fault, _] = fso.file_io.io_slippy.read_slippy_distribution(flush_file)
        for a, b in zip(example_fault.__iter__(), returned_fault.__iter__()):
            if a[0] in ['lon', 'lat', 'rake']:
                self.assertAlmostEqual(a[1], b[1], places=3)
            elif a[0] in ['slip']:
                self.assertAlmostEqual(a[1], b[1], places=5)
            else:
                self.assertAlmostEqual(a[1], b[1])
        return

    def test_io_static1d(self):
        flush_file = "test/example_static1d.txt"
        fso.file_io.io_static1d.write_static1D_source_file([example_fault, example_fault], [], flush_file)
        returned_fault, _ = fso.file_io.io_static1d.read_static1D_source_file(flush_file)
        returned_fault = returned_fault[0]
        for a, b in zip(example_fault.__iter__(), returned_fault.__iter__()):
            if a[0] in ['width']:
                self.assertLess(abs(example_fault.width-returned_fault.width), 0.2)
            elif a[0] in ['lon', 'lat']:
                self.assertAlmostEqual(a[1], b[1], places=3)
            else:
                self.assertAlmostEqual(a[1], b[1])
        return


if __name__ == "__main__":
    unittest.main()
