SPRAI
=====

Yet another method to solve the radiation transfer problem. SPRAI -- Simplex Photon Radiation in the Arepo Implementation -- is a radiation
transfer method based on the SimpleX algorithm. It is essentially a form of short characteristics ray tracing, and has the major advantage compared to the conventional
long characteristics approach that the computational cost is independent of the number of radiation sources. In SPRAI, every cell in the volume is a small source that
keeps information about the photons that pass through the cell with a specific frequency and in a certain direction.

Usage
-----

In the basic setup one needs to enable the chemistry module `SGCHEM` and the SimpleX radiation transfer module `SIMPLEX`.
One can choose between different types of sources `SX_SOURCES=int` and set the desired number of rotations `SX_NUM_ROT=int`.

In the parameter file one has to set the unit number of photons in the system `UnitPhotons_per_s` and the threshold of the photon dumping `MinNumPhotons`.
You will find more options and settings below in the corresponding sections.

The test sources can be either hardcoded (`SX_SOURCES>10`) in the header file `sx_test_sources.h` or provided as a separate binary table (`SX_SOURCES==10`). The name of the source table is specified by the parameter 'TestSrcFile' described below.

Additional Parameters
---------------------

* `UnitPhotons_per_s=float`  
  Unit photon emission from sources in units of [ph/s]. This setting is especially important when `SX_FLOAT_PRECISION` is defined. Usual values are of an order of 1e49 photons/s.
* `MinNumPhotons=float`  
  If the number of photons in the site is lower than the `Number-of-nucleons-in-the-cell * MinNumPhotons`, the photons are lost. This should reduces calculation time and remove small amounts of photons from the grid. Usual value is of the order of 1e-5.
* `TestSrcFile=string`  
  Name of the file with the test sources (i.e. `sources.bin`). This parameter is required if `SX_SOURCES=10`. The format of the binary source file is defined in `sx_photon_sources.c`. It is possible to include also values of all cross-sections and photon energies.
  For more information about the source file structure look at the file 'sx_photon_sources.c' and function 'sx_phot_sources_load()'

Additional Config.sh Options
----------------------------

* `SGCHEM_RT`  
  Main switch for SGChem that has to be turned on when using SPRAI radiation module.
* `SGCHEM_CONSTANT_ALPHAB=float`  
  This options modifies SGChem for test purposes and sets a fixed recombination rate $\alpha_\mathrm{B}$ instead of temperature dependent $\alpha_\mathrm{B}(T)$ that is calculated in the chemistry. Usually one sets $\alpha_\mathrm{B}$=2.59e-13 .
* `SGCHEM_DISABLE_COMPTON_COOLING`  
  Sets the Compton cooling rate in the SGChem cooling function to zero. This should only be used for testing purposes.
* `SIMPLEX`  
   Main switch for SPRAI.
* `SX_CHEMISTRY=int`  
  One uses almost exclusively the SGChem module, but SPRAI can be modified to work also with other chemistry modules. This sets a desired chemistry module:  
  > 3 = SGChem (default)  
  > 4 = FiBY (deprecated)   
* `SX_NDIR=int`  
  Number of direction bins. There are several options listed in the header file: `sx_direction_bins.h`
* `SX_SOURCES=int`  
  Sets the radiation sources:   
  > 4 = Calculated from properties of star particles (particles type 4)  
  > 5 = Calculated from properties of sink particles (particles type 5)    
  > 10 = Uses sources listed in a binary file. The name of the binary file is set in the parameter file using `TestSrcFile`  
  > 11+ = Uses corresponding sources defined in the header file `sx_test_sources.h`  
* `SX_RUNNING_MODE=int`  
  There are two possible running modes of SimpleX but the second type is considered as the default when running SGChem.    
  > 1 = Photons are transported only a constant number of steps (set by flag `SX_NUM_STEPS`) during one run. Unused photons are saved in the gas particles and wait until the next SimpleX run occurs. This method should in principle have advantage that the RT takes only small number of steps and can be done more frequently. Unfortunately, it cannot be used in combination with the SGChem. (deprecated)  
  > 2 = Photons are transported on the gas particle grid until they are completely attenuated. This option requires that the following flags are set: `SX_CLEAR_PHOTONS_BEFORE_RUN` and `SX_SINGLE_SHOT`. (default)
* `SX_NUM_STEPS=int`  
  Number of radiation sub-steps. If set to 0 or undefined photons are evolved until they are completely attenuated. Non-zero number of steps is used only in case of the `SX_RUNNING_MODE=1`. (deprecated)
* `SX_NUM_ROT=int`  
  Number of rotation of the direction base per one RT run. If set to 0 or undefined the directional base is static. This function can be used only in case of the full attenuation with `SX_NUM_STEPS=0` and `SX_RUNNING_MODE=2`. The possibility to increase the rotations was introduced in order to suppress the photon field inhomogeneity caused by the finite number of the directions. 
* `SX_HYDROGEN_ONLY`  
  Sets the helium mass fraction in Arepo and SGChem to zero. It is used only for testing purpose.
* `SX_PP_ZERO_SNAPSHOT`  
  This option is used in combination with the post-processing run. If set, the post-processing run creates also a zero snapshot.
* `SX_SINGLE_SHOT`  
  Controls the source radiation.  
  > defined = Sources radiate all photons only at the first RT step. This is used in combination with the `SX_RUNNING_MODE=2`. (default)  
  > undefined = The total number of photons is divided by the number of steps (`SX_NUM_STEPS`) and radiated during each RT step. (deprecated)  
* `SX_FLOAT_PRECISION`  
  If defined the photon transport uses float variables. If undefined doubles are used. The rate equation in the particular site is calculated always with doubles. The point of this flag is to reduce memory consumption of the SphP structure in the code.
* `SX_ROT_EVERY_STEP` (experimental,deprecated)
  Rotate directional base during every photon step. 
* `SX_CMB_TEMP_LIMIT` (experimental)
  This sets a temperature floor threshold for gas paticles equal to the CMB  temperature. This should avoid some 'zero-temperature' errors in during post-processing. However, there are also other ways how to limit temperature e.g.: setting Arepo parameter 'MinGasTemp' or/and SGChem parameter 'SGChemTemperatureFloor'.
* `SX_SKIP_RADIUS=float` (experimental)
  This flag can be used with the sink particles. When switched on, the attenuation of photons is not calculated within the accretion radius of the sink particles with the mass greater than the value of this parameter. The units of the SX_SKIP_RADIUS are solar masses.
* `SX_RECOMBINE`
  Recombination fix that estimates the number of recombined species in case of high recombination rates in high density regions.
* `SX_RADIATION_PRESSURE`
  Add radiation pressure from each ionization reaction to the gas 

* `SX_DISPLAY_STATS`  
  Dumps out detailed photon statistics into the output file `simplex.txt`. If undefined only basic info is printed.
* `SX_DISPLAY_SOURCES`
  Prints out all source positions and source cells during each RUN start to the output. Use with caution for large sets of sources, because printing may slow the code!
* `SX_DISPLAY_MEMORY`
  Prints out memory usage statistics into the output file `simplex.txt`.
* `SX_DISPLAY_TIMERS`
  Measure a processing time of various SPRAI code parts and prints it into the `simplex.txt` in seconds.
* `SX_DISPLAY_TIMERS_SITE`
  Print out more detailed processing time including also the site calculations.

* `SX_OUTPUT_IMAGE`  
  If this option is set Arepo creates together with the snapshot file also images (slices and projections) of the ionization (rih) and heating (hrih) rate of the Hydrogen.
* `SX_OUTPUT_IMAGE_ALL`  
  This option creates only one Arepo image file with all SimpleX rates (slices only). The Arepo image file data has dimensions (X,Y,R) where X,Y are image pixel dimensions and R is the corresponding rate index as defined in the file `sx_def.h`.
* `SX_OUTPUT_FLUX`
  This option outputs fluxes into the snapshot file for each gas cell and frequency. The fluxes are in the code units so in order to get ph/s one has to multiply the values by the `UnitPhotons_per_s`.
* `POPIII_SNE`
  Enable Pop III Supernovae
* `SINK_SIMPLEX`
  This switch links simplex and sink particles

Authors
-------

  * Ondrej Jaura (ondrej.jaura@uni-heidelberg.de)
  * Simon Glover (glover@uni-heidelberg.de)

Usage Policy and Citation
-------------------------

Please contact the authors before using this code, so that we can
avoid any unnecessary project overlap. Co-authorship on independent publications
produced with the code are not necessary, although we are happy to
actively collaborate with users of the code that are interested in
doing so.

Papers to cite:

  * [SPRAI: coupling of radiative feedback and primordial chemistry in moving mesh hydrodynamics](https://doi.org/10.1093/mnras/stx3356)     
  
