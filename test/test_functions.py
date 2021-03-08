# Testing code

import unittest
from PyCoulomb import configure_calc

class Tests(unittest.TestCase):

    def test_read_config(self):
        config_file = 'examples/example_config.txt'
        myParams = configure_calc.configure_stress_calculation(config_file);
        self.assertIsNotNone(myParams);
        return;


if __name__ == "__main__":
    unittest.main();
