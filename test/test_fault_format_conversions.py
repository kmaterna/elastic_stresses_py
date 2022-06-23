# Test the conversion between four different formats for slip distributions and fault geometry.

import unittest
import PyCoulomb.fault_slip_object as fso

example_fault = {'strike': 5, 'dip': 75, 'length': 40, 'width': 20, 'lon': -123.00, 'lat': 40.00,
                 'depth': 3, 'rake': 10, 'slip': 1, 'tensile': 0};


class Tests(unittest.TestCase):

    def test_io_geojson(self):
        flush_file = "test/example_geojson.txt";
        fso.io_geojson.write_faults_json([example_fault], flush_file);
        [returned_fault] = fso.io_geojson.read_faults_json(flush_file);
        for key in example_fault.keys():
            if key in ['rake', 'slip', 'tensile']:
                continue;   # don't test the lossy part of the conversions
            elif key == 'lon':
                self.assertAlmostEqual(example_fault[key], returned_fault[key], places=3);
            else:
                self.assertAlmostEqual(example_fault[key], returned_fault[key]);
        return;

    def test_io_pycoulomb(self):
        pycoulomb_faults = fso.io_pycoulomb.fault_dict_to_coulomb_fault([example_fault]);
        [returned_fault] = fso.io_pycoulomb.coulomb_fault_to_fault_dict(pycoulomb_faults);
        for key in example_fault.keys():
            self.assertAlmostEqual(example_fault[key], returned_fault[key]);
        return;

    def test_io_slippy(self):
        flush_file = "test/example_slippy.txt";
        fso.io_slippy.write_slippy_distribution([example_fault, example_fault], flush_file);
        [returned_fault, _] = fso.io_slippy.read_slippy_distribution(flush_file);
        for key in example_fault.keys():
            if key in ['lon', 'lat', 'rake']:
                self.assertAlmostEqual(example_fault[key], returned_fault[key], places=3);
            elif key in ['slip']:
                self.assertAlmostEqual(example_fault[key], returned_fault[key], places=5);
            else:
                self.assertAlmostEqual(example_fault[key], returned_fault[key]);
        return;

    def test_io_static1d(self):
        flush_file = "test/example_static1d.txt";
        fso.io_static1d.write_static1D_source_file([example_fault, example_fault], [], flush_file);
        returned_fault, _ = fso.io_static1d.read_static1D_source_file(flush_file);
        for key in example_fault.keys():
            if key in ['width']:
                self.assertLess(abs(example_fault[key]-returned_fault[0][key]), 0.2);
            elif key in ['lon', 'lat']:
                self.assertAlmostEqual(example_fault[key], returned_fault[0][key], places=3);
            else:
                self.assertAlmostEqual(example_fault[key], returned_fault[0][key]);
        return;


if __name__ == "__main__":
    unittest.main();
