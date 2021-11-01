# Testing code

import unittest
import Elastic_stresses_py.PyCoulomb as PyCoulomb

class Tests(unittest.TestCase):

    def test_read_config(self):
        config_file = 'examples/example_config.txt'
        myParams = PyCoulomb.configure_calc.configure_stress_calculation(config_file);
        self.assertIsNotNone(myParams);
        return;

    def test_against_Coulomb3_benchmark(self):
        """
        Running a benchmark against Coulomb3.4 and PyCoulomb.
        Input File: examples/example_inputs/simple_receiver_bm.pm, 135 degree rake
        Coulomb Output: examples/example_inputs/Coulomb34_benchmark_outs.csv
        """
        config_file = 'examples/benchmark_config.txt'
        params = PyCoulomb.configure_calc.configure_stress_calculation(config_file);
        [inputs, obs_disp_points, obs_strain_points] = PyCoulomb.input_values.read_inputs(params);
        out_object = PyCoulomb.run_dc3d.do_stress_computation(params, inputs, obs_disp_points, obs_strain_points);
        self.assertAlmostEqual(out_object.receiver_normal[0], -539.037, places=3);
        self.assertAlmostEqual(out_object.receiver_coulomb[0], -230.576, places=3);
        return;


if __name__ == "__main__":
    unittest.main();
