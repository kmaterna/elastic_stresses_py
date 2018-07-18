# Elastic_stresses_py

This code uses Okada's (1992) DC3D function to compute elastic displacements, strains, and stresses in an elastic half-space due to fault slip. It performs a similar type of calculation as Coulomb, and in fact reads Coulomb input files. I wrote this tool to iterate faster than I could using the Coulomb GUI. On test cases, it reproduces the Coulomb outputs (see Examples). 

### Capabilities: ###
* Reads source and receiver faults from .inp formatted Coulomb input files.
* Takes a single receiver fault and splits it into subfaults in the dip- and strike- directions.
* Writes .inp formatted Coulomb files with split subfaults.
* Plots displacements at the surface.
* Plots shear stress, normal stress, and Coulomb stress on receiver faults.
* Produces output tables of quantities.

### Future work: ###
* Convert to map coordinates and plot a map
* Read .inr files
* Read user-defined input files (based on EQ lat/lon, source parameters, etc.)

### Specs: ###
This code uses Python3, numpy, and matplotlib. It requires you to have Ben Thompson's Okada Python wrapper on your pythonpath (https://github.com/tbenthompson/okada_wrapper). The original Okada documentation can be found at http://www.bosai.go.jp/study/application/dc3d/DC3Dhtml_E.html. The elastic parameters mu and lamda are set in configure_calc.py. By convention, right lateral strike slip is positive, and reverse dip slip is positive. Strike is defined from 0 to 360 degrees, clockwise from north; dip is defined from 0 to 90 degrees by the right hand rule.

Example with two vertical strike-slip faults, one source (green) and one receiver (blue):
![CoulombCalc](https://github.com/kmaterna/Elastic_stresses_py/blob/master/Example/Python_Displacement_model.png)

Python-produced Coulomb Stress Changes:
![Python_stresses](https://github.com/kmaterna/Elastic_stresses_py/blob/master/Example/Python_test_case.png)

Coulomb-produced Coulomb Stress Changes:
![Coulomb_stresses](https://github.com/kmaterna/Elastic_stresses_py/blob/master/Example/Coulomb_test_case.png)

