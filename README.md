# Elastic_stresses_py

This code uses Okada's (1992) DC3D function to compute elastic displacements, strains, and stresses in an elastic half-space due to fault slip. It performs a similar type of calculation as Coulomb, and in fact reads Coulomb input files. I wrote this tool to iterate faster than I could using the Coulomb GUI. On test cases, it reproduces the Coulomb outputs (see Examples). 

## Description

### Capabilities: 
* Reads source and receiver faults from .inp formatted Coulomb input files.
* Reads source and receiver faults from .intxt files, a more convenient input format
* Reads point sources and receiver faults from .inzero, a focal mechanism input format
* Takes a single receiver fault and splits it into subfaults in the dip- and strike- directions.
* Computes elastic displacements and stresses due to slip on source faults.
* Writes .inp formatted Coulomb files with split subfaults.
* Plots displacements at the surface.
* Plots shear stress, normal stress, and Coulomb stress on receiver faults.
* Maps the faults and stress values using PyGMT.
* Produces output tables of stresses and displacements.

### Future work: 
* Read .inr files (like Coulomb)
* Reshape input arrays 

### Usage: 
Most of the flow of the program is controlled from config.txt. The elastic parameters mu and lamda, as well as your input/output options, are set in config file. You call the program by 
```bash
python elastic_stresses_driver.py config.txt
```

### New Input Formats (Not Coulomb Format): 
Source Faults (or faults have slip on them) and Receiver Faults (or faults receive stress from slip on source faults) can be specified in several types of more convenient input files beyond the .inp file that Coulomb uses. Each fault is specified by a row in the input file. 
* **Slip Format:** "S: strike rake dip length width lon lat depth slip"
	* For sources in .intxt files, best for slip distributions and fault patches
* **WC Format:** "S: strike rake dip magnitude faulting_type lon lat depth"
	* For sources in .intxt files, best for catalogs using Wells and Coppersmith (1994)
* **FM Format:** "S: strike rake dip lon lat depth magnitude mu lambda"
	* For sources in .inzero files, best for focal mechanisms
* **Receiver Format:** "R: strike rake dip length width lon lat depth"
	* For all files, describes general receiver faults
* **General Format:** "G: poissons_ratio friction_coef lon_min lon_max lon_zero lat_min lat_max lat_zero"
	* For all files, describes coordinate system and domain setup

For all finite sources, lon/lat/depth refer to the back updip corner of the fault plane, the corner where one looks along the strike direction to see the fault's upper edge (start_x and start_y in the Aki and Richards (1980) diagram in Coulomb's documentation). 

For Source Format 2, faulting_type = ["SS","R","N","ALL"] from the Wells and Coppersmith indications of Strike Slip, Reverse, Normal, and All. 

## Notes

### Sign Conventions: 
By convention, right lateral strike slip is positive, and reverse dip slip is positive. Strike is defined from 0 to 360 degrees, clockwise from north; dip is defined from 0 to 90 degrees by the right hand rule. As in Coulomb, positive shear stress is towards failure, and positive normal stress is unclamping. The original Okada documentation can be found at http://www.bosai.go.jp/study/application/dc3d/DC3Dhtml_E.html. 

### Requirements: 
This code uses Python, numpy, matplotlib, and Pygmt (Basemap is being deprecated, so I have switched to Pygmt: https://www.pygmt.org/dev/). This code also requires you to have Ben Thompson's Okada Python wrapper on your pythonpath (https://github.com/tbenthompson/okada_wrapper). It requires a few utility functions (haversine, wells_and_coppersmith) in a separate utilities repository (https://github.com/kmaterna/Utility_Code). 


## Results: 

Elastic_stresses_py reproduces the Coulomb outputs in the simple case of two vertical strike-slip faults, one source (green) and one receiver (blue):
![CoulombCalc](https://github.com/kmaterna/Elastic_stresses_py/blob/master/Example/Python_Displacement_model.png)

Python-produced Coulomb Stress Changes:
![Python_stresses](https://github.com/kmaterna/Elastic_stresses_py/blob/master/Example/Python_test_case.png)

Coulomb-produced Coulomb Stress Changes:
![Coulomb_stresses](https://github.com/kmaterna/Elastic_stresses_py/blob/master/Example/Coulomb_test_case.png)


In a real research application, a computation would look more like this: 
![Ex_Coulomb_stresses](https://github.com/kmaterna/Elastic_stresses_py/blob/master/Example/Coulomb_map.png)
