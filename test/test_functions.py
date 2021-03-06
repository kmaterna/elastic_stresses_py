# Testing code

import numpy as np
import unittest
from PyCoulomb import conversion_math
from PyCoulomb import configure_calc

class Tests(unittest.TestCase):

    def test_strike(self):
        strike = conversion_math.get_strike(deltax=1.0, deltay=0.0);
        self.assertEqual(strike, 90);
        strike = conversion_math.get_strike(deltax=-1.0, deltay=0.0);
        self.assertEqual(strike, 270);
        strike = conversion_math.get_strike(deltax=0.0, deltay=-1.0);
        self.assertEqual(strike, 180);
        angle = -160;
        strike = conversion_math.get_strike(deltax=np.cos(np.deg2rad(angle)), deltay=np.sin(np.deg2rad(angle)));
        self.assertEqual(strike, 250);
        return;

    def test_rake(self):
        rake = conversion_math.get_rake(strike_slip=1.0, dip_slip=0.5);
        self.assertAlmostEqual(rake, 26.5650511770779);
        return;

    def test_plane_normal(self):
        plane_normal = conversion_math.get_plane_normal(strike=0, dip=0);
        self.assertAlmostEqual(plane_normal[0], 0);
        self.assertAlmostEqual(plane_normal[1], 0);
        self.assertAlmostEqual(plane_normal[2], 1);
        plane_normal = conversion_math.get_plane_normal(strike=90, dip=89.99);
        self.assertAlmostEqual(plane_normal[0], 0);
        self.assertAlmostEqual(plane_normal[1], -1);
        self.assertAlmostEqual(plane_normal[2], 0.0001745329);
        plane_normal = conversion_math.get_plane_normal(strike=270, dip=89.99);
        self.assertAlmostEqual(plane_normal[0], 0);
        self.assertAlmostEqual(plane_normal[1], 1);
        self.assertAlmostEqual(plane_normal[2], 0.0001745329);
        plane_normal = conversion_math.get_plane_normal(strike=180, dip=1);
        self.assertAlmostEqual(plane_normal[0],  -0.0174524064372835);
        self.assertAlmostEqual(plane_normal[1], 0);
        self.assertAlmostEqual(plane_normal[2], 0.999847695156391);
        return;

    def test_read_config(self):
        config_file = 'examples/example_config.txt'
        myParams = configure_calc.configure_stress_calculation(config_file);
        self.assertIsNotNone(myParams);
        return;


if __name__ == "__main__":
    unittest.main();
