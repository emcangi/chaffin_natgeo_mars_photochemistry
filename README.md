# Mars D/H fractionation as a function of atmospheric parameters

---
## Description/Overview

This project was forked from planetarymike/chaffin_natgeo_mars_photochemistry and modified substantially to study the effects of atmospheric temperature and water vapor content on the D/H fractionation on Mars. We submitted a paper using this code to JGR: Plants on 16 July 2020.

The code is a 1D photochemical model of the Martian atmosphere, with the only coordinate being altitude in the atmosphere. The model output is stored as .h5 files, in a matrix of species number density (rows) by altitude (columns). 

The model includes the following physics, chemistry, and computation:

- Chemistry
	- Species: (majors) CO2, O2, O3, O1D, O, CO, 
			   (H-bearing) H, H2, OH, HO2, H2O2, HOCO, H2O (fixed profile)
			   (D-bearing) D, HD, OD, DO2, HDO2, DOCO, HDO (fixed profile)
               (Nonreactive) N2, Ar
               (ions) CO2+ (fixed profile to approximate ionosphere)
    - Current reaction rates 
- Photochemistry 
	- Cross sections
	- Solar UV irradiation appropriate for solar mean conditions
	- Calculation of photodissociation rates
- Vertical transport via eddy diffusion, molecular diffusion
- Upper and lower boundary conditions
	- Outgoing flux
	- Number density
	- Particle velocity
- Coupled network of production and loss equations
	- Chemistry
	- Transport
	- Solution using the Gear method

---
## Requirements

Languages:
- Julia 1.4.1
- Python 3+ (some of the analysis files are in Python because I didn't want to make the effort to fit lines in Julia)

Packages (Julia), used in either all scripts or some scripts:
- PyPlot
- PyCall
- HDF5
- JLD
- LaTeXStrings
- Distributed
- DelimitedFiles
- SparseArrays
- LinearAlgebra
- ProgressMeter
- Printf
- DataFrames
- Analysis.jl (my custom module including functions shared by multiple scripts within this project. Included here at Github)

---
## Script/file call order and usage

Files are called as follows:

(0) Set up files - optional

I used all of these, but anyone trying to use my code should only need to use these if the goal is to make a significant change to the model before using it.

(0a) `HDO_crosssection.py`: Uses data for HDO cross sections from Cheng+1999 and Cheng+2004 to fully extrapolate to missing wavelengths ~121 and 220 nm. Figure S1 in paper.
	
(0b) `OD_crosssection.py` - Smashes together OD cross section data between 115-180 nm from Nee 1985 with OH cross sections from Barfield 1972. Figure S2. NOTE: this file doesn't exist yet, the code is buried in a jupyter notebook. I'll get it in here I promise.

(0c) `SolarCycle.ipynb`: Scales data on solar irradiance by wavelength from Earth orbiters for solar mean conditions to Mars' orbit. 


(1) Model files - required

(1a) `converge_new_file.jl`

This can be used to either (1) converge a new atmosphere with a different altitudinal extent, starting from a converged file for a 0-200 km atmosphere, or (2) use the default 250 km atmosphere to run a simulation for given temperature or water vapor parameters. It can also be used to test different atmospheric D/H ratios and values of escape of O in #/cm^2/s, but I only used those for testing and did not include in the paper.

Note that the code is written so that you can converge an atmosphere of any new max altitude, BUT it assumes the new max altitude is higher than 200 km. Also, the max altitude of 250 km is hard coded in `PARAMETERS.jl`, so even if you choose some other max altitude, the code will probably break because I didn't udpate it to account for user-chosen max altitudes since I only needed to go up to 250 km.

Inputs for running new simulations:
- converged_250km_atmosphere.h5: a pre-converged atmosphere. The water profiles are goofy in the upper atmosphere, which is why this gets re-converged once new water profiles are prescribed within the simulation run.
- cross section data, stored in uvxsect folder. For a list of files, see Photochemical Crosssections/Crosssection Files within this file. 
- `marssolarphotonflux_solarmean.dat` - photon flux (γ/s/cm^2/nm) by wavelength for 0-2400 nm. Files for solar max, solar min also available.
- Temperatures: either specified by calling the file with command line args "temp Tsurf Ttropo Texo" or the mean values
- Water vapor mixing ratio in the lower atmosphere: This is a little weird, you have to try a few options while printing the resulting pr μm of water vapor to see what mixing ratio corresponds to what pr μm, aiming for pr μm of 1, 10, 25, 50, 100. Sorry.
Outputs:
- A converged atmosphere h5 file with a specific name that depends on the type of experiment run (temperature or water vapor). The file contains an array of species densities by altitude.

(1b) `PARAMETERS.jl`

This file contains key global constants and general simulation parameters, as well as the chemistry reaction network. 

(1c) `plot_profiles.jl`

This allows the user to plot the temperature and water vapor profiles that are used in the model. Produces Figures 1-3.

- Individual temperature profiles: can be specified, or can do the standard atmosphere temperature profile.
- 6-panel temperature profiles: "Climate extrema" profiles; profiles with cold/warm temperature at the surface/tropopaues/exobase.
- water profiles: This will plot the H2O and HDO profiles for 1, 10, 25, 50, 100 pr μm total water. 
- 3-panel temperature profiles: This plots all the temperature profiles used when each control temperature is tweaked by a small amount. 

(2) Analysis of model - required scripts

These scripts are required to analyze the model. If you run all these scripts, you will get all the plots contained in the paper. And perhaps some extra information too!

The order here doesn't matter except check_eq.jl should be run on each model output file and return an affirmative result before making any plots.

(2a) `check_eq.jl`: Used to check that the atmospheric escape has reached stoichiometric balance, such that Φ_H + Φ_D = 2Φ_O

(2b) `plot_all_results.jl`: Contains code to make all results Figures (Figures 4 through 7, 10, 11, S4)

(2c) `plot_water_loss.jl`: Makes a plot of water loss as a function of various assumptions about current water inventory and the results for f. Figure 8 in paper.

(2d) `plot_waterloss_results_vs_others.py`: Simple comparison plot, Figure 9 in paper 

(2e) `reproductions_plot.jl`: Makes a plot of results for f when reproducing past studies. Includes past study values for comparison. Figure S3

(2f) `plot_flux_contributions.jl`: Plots the contribution of H, D, H2, and HD to the escape flux of H and D as a function of tropopause and exobase temperatures.
                                   Figure S5. Shows that estimates of non-thermal escape may be a bit off for low temperatures.

(2g) `plot_rxn_rates.jl`: Plots chemical reaction rates for atomic H or D as a function of altitude, as well as flux quantities by altitude, for each experiment.
                          These were created after receiving reviewer comments and are intended to demonstrate that higher temperatures leads to a 
                          significant relative increase in flux of D vs for H.

---

## Notes for using this code

Note that any instructions here are written for Ubuntu Linux. They probably also work on Macs, but I haven't used Macs in a long time. If you're using Windows, sorry, you're on your own... 

### Setting up the main folder for the codebase

Add the folder in which you will store the scripts to your julia startup file. Assuming the folder is /home/username/dh_fractionation/:

```shell
	~/.julia/config
	vim startup.jl
	i
	cd("/home/username/dh_fractionation/")
	:wq
```

### Using the custom Julia modules

In order for the various files to load `Analysis.jl` and `Photochemistry.jl`, you will need to add the root folder to your machine's `$JULIA_LOAD_PATH` variable. The codebase is not set up to specify a definite path for this module. 

This section assumes that you have downloaded the contents of the repo to /home/username/dh_fractionation/.

To see what's on your `JULIA_LOAD_PATH` variable:

```shell
echo $JULIA_LOAD_PATH
```

To add the folder to the `JULIA_LOAD_PATH` variable, add the following line to either your ~/.bashrc or ~/.profile files:

`export JULIA_LOAD_PATH=$JULIA_LOAD_PATH:/home/username/dh_fractionation/CustomJuliaModules`

Within CustomModules you should have the Analysis folder and Photochemistry folders (they are already organized there in the repo)

---

## Support

If you have any questions, please contact the lead author/repo owner at: eryn.cangi@colorado.edu