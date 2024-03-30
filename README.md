# Elastic_stresses_py
[![Python 3.9+](https://img.shields.io/badge/python-3.9-blue.svg)](https://www.python.org/downloads/release/python-390/)
[![DOI](https://zenodo.org/badge/141371162.svg)](https://zenodo.org/badge/latestdoi/141371162)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](https://github.com/kmaterna/Strain_2D/blob/dev/license.md)


This code uses Okada's (1992) formulation to compute elastic displacements, strains, stresses, and Coulomb failure stresses (e.g., King et al., 1994) in a linear, isotropic elastic half-space due to fault slip. It performs a similar type of calculation as Coulomb, and in fact reads Coulomb input files. I wrote this tool to iterate faster than I could using the Coulomb GUI. On test cases (in ```examples/```), it reproduces the Coulomb outputs. 

The equation being solved, following the sign conventions used in Coulomb, is:

$$ \Delta CFS = \Delta\tau + \mu (\Delta\sigma_n - \Delta P) $$ 

where the instantaneous pore pressure change P is determined by the mean stress and the Skempton's coefficient (B):

$$ \Delta P = {B} * {\sigma_{kk} \over 3} $$

e.g., Beeler et al., 2000.


## Installation and Usage

### Requirements:
* Python 3.9+
* gfortran and gcc (For mac/linux, this is done through brew, port, apt-get, whatever works for your system)
* Ben Thompson's Python [cutde](https://github.com/tbenthompson/cutde). 
* Kathryn Materna's [Tectonic_Utils](https://github.com/kmaterna/Tectonic_Utils).  
* Several standard Python libraries such as numpy, matplotlib, and [PYGMT](https://www.pygmt.org/dev/); these are listed in ```environment.yml```.  


### Installation on command line
To install Elastic_stresses_py, first clone this library onto your computer with ```git clone ``` and the address of the git repository (see the green **Code** button). 
The easiest way to gather all the dependencies is to create a new conda environment with ```conda env create -f environment.yml``` in the directory where you've cloned the repository. 
Then, in the new conda environment ```elastic_py``` or similar, run ```pip install . ``` from the directory where you've cloned the repository. 

**NOTE 1:** Mac users switching from Intel architecture to M1 architecture may experience some errors when compiling code on their new architectures for the first time.  
Such errors may appear as ```(mach-o file, but is an incompatible architecture)``` or similar.  If this happens, please re-install Xcode on your new architecture with ```xcode-select --install``` and try again.


### Installation in Anaconda GUI:

* Download the code's zip file from this Github page. Click Code --> Download Zip.  Save it somewhere onto your local disk. 
  
* Create a new Python environment in the Anaconda GUI. Environments Tab --> Import --> Local drive --> navigate to the local Elastic_stresses_py folder and select ```environment.yml``` --> Open.  Select name=espy (or desired) --> Create environment.  This will take a few minutes. 
    
* To install ```Elastic_stresses_py``` itself, we will use the terminal once more. Move to the location where you downloaded the source code and type ```pip install .``` 
  
* After refreshing the Anaconda navigator, you will be able to see all the installed packages in your new environment, including Elastic_stresses_py. 

* On the Anaconda home tab, install Jupyter Notebook in your new environment by clicking the Install button.
  
* Close all your terminals, notebooks, pythons, and anacondas. Then re-open.  

* In the Anaconda Environments tab, click the play button beside the new environment you just created, and select "Open with Jupyter Notebook".  This will open in a browser window and display the contents of your home directory.  Navigate to the location of the Elastic_stresses_py package and then examples/Example_2/. Open the Jupyter nobebook Run_Elastic_stresses_py.ipynb. 

* You can use this Jupyter notebook by clicking each cell (i.e., any body of text or code) and executing it with shift+Enter. Start at the top and work your way down the tutorial. Some code lines will take a moment and produce screen output into the notebook. More information about using Jupyter notebooks can be found [here](https://www.dataquest.io/blog/jupyter-notebook-tutorial/).  


### Citing:
If you use this code in your research, please cite the DOI for the current version:

[![DOI](https://zenodo.org/badge/141371162.svg)](https://zenodo.org/badge/latestdoi/141371162)

Example citation: Materna, K. (2023). Elastic_stresses_py (version 1.0.0). https://github.com/kmaterna/Elastic_stresses_py. Archived at DOI: 10.5281/zenodo.79751979.

Also cite the underlying theoretical work behind elastic dislocations and Coulomb Failure Stress (Okada, 1992; King et al., 1994; and others). 

## Software Usage:
The main executable is ```elastic_stresses_driver```, which takes a config file as the first argument.

Most of the behavior of the program is controlled by the config text file.  From an experiment directory on your system, you can generate a default config file in the current working directory (.) with: 
```bash
elastic_stresses_config_writer .
```
You should change the parameters to your own experiment needs.  An example config file is provided in examples/. The elastic parameters mu and lambda, as well as your input/output options, are set in config file. 

Then, you call the program by passing the config file into the main executable: 
```bash
elastic_stresses_driver my_config.txt
```

*New:* 
To get the hang of the usage, check out our recently added calculation in ```examples/Example_2``` through a Jupyter Notebook, located [here](examples/Example_2/Run_Elastic_stresses_py.ipynb).

### Config Parameters

A complete description of the config file is shown below, with both required and optional arguments. Optional arguments can be omitted or be left blank.

![CoulombCalc](https://github.com/kmaterna/Elastic_stresses_py/blob/master/examples/pngs/annotated_config.png)

#### Important Config Switches
* If ```plot_grd_disp```: Will produce grid of 100x100 synthetic points and their surface displacements, in txt/grd/plots [default: True]
* If ```plot_stress```: Will compute and plot shear stress, normal stress, and Coulomb stress on receiver faults [default: True]

### ```Intxt``` and ```Inzero``` Input Formats: 
Source Faults (or faults that have slip on them) and Receiver Faults (or faults that receive stress from slip on source faults) can be specified in several human-readable formats beyond the .inp file that Coulomb uses. In ```.inzero``` and ```.intxt``` text files, each fault is specified by a row in an input text file of extention .intxt or .inzero. Valid rows of the input file include: 
* **General Format:** Describes coordinate system and domain setup. Required.  
    * "General: poissons_ratio* friction_coef lon_min lon_max lon_zero lat_min lat_max lat_zero" 
* **Receiver Format:** Describes receiver faults 
    * "Receiver: strike rake dip length_km width_km lon lat depth_km"
* **Slip Format:** For slip distributions and fault patches. 
    * "Source_Patch: strike rake dip length_km width_km lon lat depth_km slip_m (opt: tensile_m)"
* **WC Format:** For catalogs using Wells and Coppersmith (1994) 
    * "Source_WC: strike rake dip magnitude faulting_type lon lat depth_km" 
* **FM Format:** For focal mechanisms 
    * "Source_FM: strike rake dip lon lat depth_km magnitude"
* **Horizontal Profile Format:** Specify an orientation and compute stresses on that plane/orientation over an area. Like a horizontal cross-section.
    * "Receiver_Horizontal_Profile: depth_km strike dip rake centerlon centerlat length_km width_km inc_km" 
* **Mogi Source Format:** Specify a location, depth in km, and volume change in meters^3.  **Only implemented for displacement at the moment**, not stresses and strains. 
    * "Source_Mogi: lon lat depth_km dV_m3"  

*PR1 in the intxt/inzero format is never actually used; it's just a placeholder from the old Coulomb format. 
Specifying Poisson's ratio should be done at the computation level through lame1 and mu in the config file, thus altering alpha in the Okada formulation.
On the output stage, the actual poisson's ratio from mu and lame1 will be calculated and written into ```used_inputs.txt```.

For all finite sources (i.e., Patch, WC), lon/lat/depth refer to the back updip corner of the fault plane, the corner where one looks along the strike direction to see the fault's upper edge (start_x and start_y in the Aki and Richards (1980) diagram in Coulomb's documentation). 

For WC Source Format, faulting_type = ["SS","R","N","ALL"] from the Wells and Coppersmith indications of Strike Slip, Reverse, Normal, and All. 

An example file ```input.intxt``` with a "Slip Format" source might look like: 
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

Set domain elastic parameters at the computation level through lame1 and mu in the config file.



### List of Code Capabilities:
* Reads source and receiver faults from .inp formatted Coulomb input files.
* Reads source and receiver faults from .intxt files, a more convenient input format where slip is specified by length/width or by Wells and Coppersmith (1994) scaling rules
* Reads point sources and receiver faults from .inzero, a focal mechanism / moment tensor input format
* Takes a single receiver fault and splits it into subfaults in the dip- and strike- directions.
* Computes elastic displacements and stresses due to slip on source faults.
* Writes .inp formatted Coulomb files with split subfaults.
* Maps the faults and stress values using PyGMT.
* Produces output tables of stresses and displacements.

### Future work: 
* Output computations at depths other than the surface
* Read .inr files (like Coulomb)


## Example 1: Command Line

Running from the command line, from ```examples/Example_1``` in this repository, type:
```bash
elastic_stresses_driver my_config.txt
```
where my_config.txt is a text file in the local directory containing:
```
[io-config]
exp_name = my_experiment
input_file = M6.8_2014.intxt
output_dir = Outputs/
plot_stress = 1
plot_grd_disp = 1
gps_disp_points = CA_GPS_ll.txt
aftershocks = CA_aftershocks_2014.txt
strain_file = 

[compute-config]
strike_num_receivers = 10
dip_num_receivers = 10
mu = 30e9
lame1 = 30e9
B = 0
fixed_rake = 
```

should produce files and plots such as:

![Gallery](https://github.com/kmaterna/Elastic_stresses_py/blob/master/examples/pngs/example_plots.png)

## Example 2: Jupyter Notebook

We can do the same calculation as above in ```examples/Example_2``` through a Jupyter Notebook, located [here](examples/Example_2/Run_Elastic_stresses_py.ipynb).    

This notebook explores modifying the inputs and configuration parameters through the Python API (in addition to the regular command line).  This is more similar to a real research application.  



## Example 3: Real Application
The code can also be used for larger numbers of source faults and receiver faults. 
The largest we have tried is 10,000 sources and 10,000 receivers. That application took about half an hour and the results were compiled into something like this:
 
![NZ](https://github.com/kmaterna/Elastic_stresses_py/blob/master/examples/pngs/nz_example.png)


## References:

* Beeler, N. M., Simpson, R. W., Hickman, S. H., & Lockner, D. A. (2000). Pore fluid pressure, apparent friction, and Coulomb failure. Journal of Geophysical Research: Solid Earth, 105(B11), 25533-25542.

* King, G. C., Stein, R. S., & Lin, J. (1994). Static stress changes and the triggering of earthquakes. Bulletin of the Seismological Society of America, 84(3), 935-953.

* Okada, Y. (1992). Internal deformation due to shear and tensile faults in a half-space. Bulletin of the seismological society of America, 82(2), 1018-1040.

* Wells, D. L., and K. J. Coppersmith (1994). New empirical relationships among magnitude, rupture length, rupture width, rupture area, and surface displacement. Bulletin of the seismological Society of America 84.4: 974-1002.
