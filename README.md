# Elastic_stresses_py

This code uses Okada's (1992) DC3D function to compute elastic displacements, strains, and stresses in an elastic half-space due to fault slip. It performs a similar type of calculation as Coulomb, and in fact reads Coulomb input files. I wrote this tool to iterate faster than I could using the Coulomb GUI. On test cases, it reproduces the Coulomb outputs (see Examples). 

## Description

### Code Requirements: 
This code uses Python3, numpy, matplotlib, and Pygmt (https://www.pygmt.org/dev/). This code also requires you to have Ben Thompson's Okada Python wrapper on your pythonpath (https://github.com/tbenthompson/okada_wrapper). It requires a few utility functions in a separate utilities repository (https://github.com/kmaterna/Tectonic_Utils).  To get the utilities library, please ```pip install Tectonic-Utils```. 

### Installation and Usage: 
To install, you can clone this library onto your computer into a location that is on your $PYTHONPATH (in other words, the parent directory holding Elastic_stresses_py/ must be on your $PYTHONPATH).  For convenience, it is nice (but not required) to put the full path to Elastic_stresses_py/PyCoulomb/ on your regular system $PATH so the main executable ```elastic_stresses_driver.py``` can be found.   

Most of the behavior of the program is controlled by a config text file. An example config file is provided in examples/. The elastic parameters mu and lamda, as well as your input/output options, are set in config file. You call the program by 
```bash
elastic_stresses_driver.py config.txt
```

### Capabilities: 
* Reads source and receiver faults from .inp formatted Coulomb input files.
* Reads source and receiver faults from .intxt files, a more convenient input format wheere slip is specified by length/width or by Wells and Coppersmith (1994) scaling rules
* Reads point sources and receiver faults from .inzero, a focal mechanism input format
* Takes a single receiver fault and splits it into subfaults in the dip- and strike- directions.
* Computes elastic displacements and stresses due to slip on source faults.
* Writes .inp formatted Coulomb files with split subfaults.
* Plots displacements at the surface.
* Plots shear stress, normal stress, and Coulomb stress on receiver faults.
* Maps the faults and stress values using PyGMT.
* Produces output tables of stresses and displacements.

### Future work: 
* Generate a valid but empty config file and/or input file
* Read in full moment tensor (not just double couple focal mechanisms)
* Make output vector plots (gps obs vs. model etc.) to pygmt
* Read .inr files (like Coulomb)

### Additional Input Formats beyond Coulomb Format: 
Source Faults (or faults that have slip on them) and Receiver Faults (or faults that receive stress from slip on source faults) can be specified in several human-readable formats beyond the .inp file that Coulomb uses. Each fault is specified by a row in an input text file of extention .intxt or .inzero. Valid rows of the input file include: 
* **General Format:** "General: poissons_ratio friction_coef lon_min lon_max lon_zero lat_min lat_max lat_zero"
	* Describes coordinate system and domain setup. Required. 
* **Receiver Format:** "Receiver: strike rake dip length_km width_km lon lat depth_km"
	* Describes receiver faults
* **Slip Format:** "Source_Patch: strike rake dip length_km width_km lon lat depth_km slip_m"
	* For slip distributions and fault patches.  
* **WC Format:** "Source_WC: strike rake dip magnitude faulting_type lon lat depth_km"
	* For catalogs using Wells and Coppersmith (1994)
* **FM Format:** "Source_FM: strike rake dip lon lat depth_km magnitude mu lambda"
	* For focal mechanisms
* **MT Format:** "Source_MT: Mrr Mtt Mpp Mrt Mrp Mtp strike rake dip lon lat depth_km mu lambda"
    * For full moment tensors, not just double-couple

For all finite sources (i.e., Patch, WC), lon/lat/depth refer to the back updip corner of the fault plane, the corner where one looks along the strike direction to see the fault's upper edge (start_x and start_y in the Aki and Richards (1980) diagram in Coulomb's documentation). 

For WC Source Format, faulting_type = ["SS","R","N","ALL"] from the Wells and Coppersmith indications of Strike Slip, Reverse, Normal, and All. 

An example .intxt file with WC source format might look like: 
```
# General: poissons_ratio friction_coef lon_min lon_max lon_zero lat_min lat_max lat_zero
# Source_Patch: strike rake dip length_km width_km lon lat depth_km slip_m
# Receiver: strike rake dip length_km width_km lon lat depth_km

General: 0.25 0.40 -125.80 -122.60 -124.50 39.30 41.70 40.30 
Source_Patch: 228.0 -2.0 79.0 44.25 11.91 -125.126 40.834 16.40 1.124
Receiver: 355.0 90.0 12.0 140.00 120.00 -124.560 40.300 10.80 
```

## Notes

### Sign Conventions and Units: 
By convention, right lateral strike slip is positive, and reverse dip slip is positive. Strike is defined from 0 to 360 degrees, clockwise from north; dip is defined from 0 to 90 degrees by the right hand rule. As in Coulomb, positive shear stress is towards failure, and positive normal stress is unclamping. The original Okada documentation can be found at http://www.bosai.go.jp/study/application/dc3d/DC3Dhtml_E.html. 

For input files, strike/rake/dip have units of degrees. Length/width/depth have units of km. Slip has units of meters.   


## Results: 

Elastic_stresses_py reproduces the Coulomb outputs in the simple case of two vertical strike-slip faults, one source (green) and one receiver (blue):
![CoulombCalc](https://github.com/kmaterna/Elastic_stresses_py/blob/master/examples/pngs/Python_Displacement_model.png)

Python-produced Coulomb Stress Changes:
![Python_stresses](https://github.com/kmaterna/Elastic_stresses_py/blob/master/examples/pngs/Python_test_case.png)

Coulomb-produced Coulomb Stress Changes:
![Coulomb_stresses](https://github.com/kmaterna/Elastic_stresses_py/blob/master/examples/pngs/Coulomb_test_case.png)


In a real research application, a computation would look more like this: 
![Ex_Coulomb_stresses](https://github.com/kmaterna/Elastic_stresses_py/blob/master/examples/pngs/Coulomb_map.png)
