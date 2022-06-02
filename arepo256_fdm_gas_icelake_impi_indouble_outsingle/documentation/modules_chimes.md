
CHIMES 
====== 

The CHIMES module tracks the time-dependent evolution of 157 ions and molecules, and uses the resulting ion and molecule abundances to compute the radiative cooling rates of the gas. This replaces the standard cooling routines in AREPO. 

The details of the CHIMES model are described in `Richings+ 2014a`_ and `Richings+ 2014b`_.  

.. _Richings+ 2014a: http://adsabs.harvard.edu/abs/2014MNRAS.440.3349R 
.. _Richings+ 2014b: http://adsabs.harvard.edu/abs/2014MNRAS.442.2780R 


Usage 
----- 

The CHIMES module can be switched on using the ``CHIMES`` Config option, along with the ``COOLING`` Config option. There are several additional Config options that control the behavior of CHIMES. These are described further below.

CHIMES can also be used in conjunction with the ``USE_SFR`` option, for simulations that include star formation. When the `Springel & Hernquist 2003`_ multiphase ISM model is used, gas cells that are above the density threshold are set to chemical equilibrium, using pre-computed equilibrium abundance tables. The abundances of cells that lie on the equation of state (EOS) are computed using the effective temperature of the EOS. However, when computing the EOS, CHIMES uses the equilibrium abundances in the hot phase to compute the cooling rate of the hot phase, which goes into computing the effective temperature of the EOS.

.. _Springel & Hernquist 2003: http://adsabs.harvard.edu/abs/2003MNRAS.339..289S

Note that there are several Config options in AREPO that are incompatible with CHIMES. These include: ``GFM_COOLING_METAL``, ``GFM_AGN_RADIATION``, ``GFM_UVB_CORRECTIONS``, ``GFM_LAMBDA``, ``FM_RADIATION_FEEDBACK``, ``FM_STOCHASTIC_HII_PHOTOIONIZATION``, ``LOCAL_FEEDBACK``, ``RADCOOL_HOTHALO``, ``MODIFIED_EOS``, ``UVB_SELF_SHIELDING``, ``DELAYED_COOLING``, ``DELAYED_COOLING_TURB``, ``SIMPLE_COOLING``, ``EXPLICIT_COOLING``, ``FLD``, ``PREHEATING``, ``COSMIC_RAYS``, ``GRACKLE``, ``TILE_ICS``. There may also be other Config options that are incompatible with CHIMES.

To run CHIMES, you will also need the CHIMES data files that contain the various reaction rate coefficients, photoionisation cross sections etc. These can be downloaded from the following bitbucket repository:

https://bitbucket.org/richings/chimes-data

The paths to the various data files can then be specified via the parameter file, as described below.

Additional Parameters
---------------------

* ``Chimes_data_path`` Path to the directory containing the CHIMES data files, which you downloaded from the above bitbucket repository. The directory structure of this repository should not be changed - CHIMES will then be able to find the various individual files that it needs within this directory. 

* ``PhotoIon_table_path`` Path to a text file that contains the path to the photoionisation cross sections HDF5 file that will be used for the given UV radiation field. The HDF5 file is specified in this way, via a text file, because CHIMES can take multiple UV radiation fields, in which case the text file would contain multiple paths, one per line for each cross sections HDF5 file, corresponding to each individual UV radiation field. Currently, we only include a single UV radiation field in AREPO, so the text file will only contain a single line, but we may add options for multiple radiation fields in AREPO in the future. Note that, if the ``CHIMES_REDSHIFT_DEPENDENT_UVB`` Config option is enabled, the cross sections table that you specify here will be replaced according to the current redshift. However, you need to specify something here so that CHIMES can be properly intialised. For a cosmological run, it is best to set this to the HM12 cross sections table closest to the initial redshift. 

* ``EqAbundance_table_path`` Path to the HDF5 table containing the pre-computed equilibrium abundances. If you are using the `Springel & Hernquist 2003`_ multiphase ISM model, this table is used to set gas cells above the density threshold to chemical equilibrium. We typically use the ``chimes-data/EqAbundancesTables/HM12/EqAbundances_HM12_noShield.hdf5`` table here, although you can create your own tables as well. If you are not using chemical equilibrium at all, you can use the dummy table: ``chimes-data/EqAbundancesTables/DummyTable.hdf5``. This contains only a handful of values, and so is much faster to read in. 

* ``Thermal_Evolution_On`` Flag to switch on temperature evolution in CHIMES. When this flag is 0, the chemical abundances will be evolved at constant temperature in CHIMES. When it is 1, CHIMES will evolve the temperature along with the chemical abundances. *Typical value: 1*. 

* ``Chemistry_eqm`` Flag to determine whether to folow the time-dependent chemical evolution, or set the abundances to chemical equilibrium. When this flag is 0, the full time-dependent evolution is used. When it is 1, the abundances are set to equilibrium using the pre-computed equilibrium abundance tables. When it is 2, the equilibrium abundances are computed on the fly by setting the chemical rate equations to zero. Note that, if you are using the `Springel & Hernquist 2003`_ multiphase ISM model, gas cells above the density threshold are set to chemical equilibrium, using option 1, regardless of what this parameter is set to. This parameter then determines the behavior of the low-density gas cells. *Typical value: 0*.

* ``StaticMolCooling`` The molecular cooling from CO and H2O depends on line broadening. When this flag is set to 0, we take into account the divergence of the velocity here. When this flag is set to 1, we only include thermal broadening for the CO and H2O cooling. *Typical value: 0*. 

* ``CellSelfShieldingOn`` Flag to determine whether self shielding is included. When this flag is 0, self shielding is not included. When it is 1, it is included, and the column densities of individual species (e.g. HI, H2 etc.) are updated throughout the course of the integration of the rate equations. However, updating the column densities in this way is very slow. When this flag is 2, self shielding is included, but column densities are only updated at the beginning of the hydrodynamic time-step, and not throughout the course of the CHIMES integration. Note that, to switch on self shielding, you will also need to enable either the ``CHIMES_JEANS_SHIELDING`` or the ``CHIMES_SOBOLEV_SHIELDING`` Config option, which will define the shielding length. *Typical value: 2*. 

* ``ShieldingLengthFactor`` The shielding length (either the Jeans length or the Sobolev length) is multiplied by this factor. This allows you to control the normalisation of the shielding length. 

* ``MaxShieldingLength_kpc`` The shielding length (either the Jeans length or the Sobolev length) is capped to be no greater than this maximum value, in kpc. *Typical value: 100*. 

* ``Grain_Temperature`` The temperature of dust grains in Kelvin, as used when computing the formation rate of H2 on dust grains. *Typical value: 10*. 

* ``CrRate`` Cosmic ray ionisation rate of HI. The cosmic ray ionisation rate of all other species are then scaled relative to this parameter. *Typical value: 1.8e-16*. 

* ``max_mol_temperature`` Molecules are excluded above this temperature. *Typical value: 1.0e5*. 

* ``IsotropicPhotonDensity`` This sets the normalisation of the UV radiation field. It is defined as the number density of hydrogen-ionising photons in cgs units. The flux of hydrogen-ionising photons is then the IsotropicPhotonDensity multiplied by the speed of light. If the ``CHIMES_REDSHIFT_DEPENDENT_UVB`` Config option is switched on, this value will be overwritten with the correct value for the given redshift, but it still needs to be set when CHIMES is first initialised. 

* ``relativeTolerance`` This controls the accuracy of the thermo-chemistry integration. CHIMES will sub-cycle each hydro time-step, aiming to achieve a relative error as given by this parameter. *Typical value: 1.0e-4*. 

* ``absoluteTolerance`` This controls the accuracy of the chemistry integration. Species with an abundance much below the absolute tolerance are not taken into account when determining the sub-steps needed to achieve a given relative error. *Typical value: 1.0e-6 (when scale_metal_tolerances is set to 1)*. 

* ``thermalAbsoluteTolerance`` As the ``absoluteTolerance`` parameter, but for the thermal energy. *Typical value: 1.0e-40*. 

* ``scale_metal_tolerances`` In cosmological simulations, the metal element abundances can be arbitrarily small. If an element abundance is zero, its ions will not be included in the chemical network. However, if it is non-zero but lower than the absolute tolerance, this can cause problems for the chemical integration, because none of that element's ions are taken into account when determining the sub-steps. To avoid this problem, set the ``scale_metal_tolerances`` parameter to 1. Then the absolute tolerance of each individual ion and molecule species will be set to the ``absoluteTolerance`` parameter multiplied by that species' corresponding element abundance. If this parameter is set to 0, CHIMES will use a constant absolute tolerance. *Typical value: 1*. 

* ``IncludeCarbon`` Set this flag to 1 to include carbon in the CHIMES network. Set it to 0 to exclude carbon. *Typical value: 1*. 

* ``IncludeNitrogen`` Set this flag to 1 to include nitrogen in the CHIMES network. Set it to 0 to exclude nitrogen. *Typical value: 1*. 

* ``IncludeOxygen`` Set this flag to 1 to include oxygen in the CHIMES network. Set it to 0 to exclude oxygen. *Typical value: 1*. 

* ``IncludeNeon`` Set this flag to 1 to include neon in the CHIMES network. Set it to 0 to exclude neon. *Typical value: 1*. 

* ``IncludeMagnesium`` Set this flag to 1 to include magnesium in the CHIMES network. Set it to 0 to exclude magnesium. *Typical value: 1*. 

* ``IncludeSilicon`` Set this flag to 1 to include silicon in the CHIMES network. Set it to 0 to exclude silicon. *Typical value: 1*. 

* ``IncludeSulphur`` Set this flag to 1 to include sulphur in the CHIMES network. Set it to 0 to exclude sulphur. *Typical value: 1*. 

* ``IncludeCalcium`` Set this flag to 1 to include calcium in the CHIMES network. Set it to 0 to exclude calcium. *Typical value: 1*. 

* ``IncludeIron`` Set this flag to 1 to include iron in the CHIMES network. Set it to 0 to exclude iron. *Typical value: 1*. 

The following parameters are related to an obsolete option in CHIMES that is no longer used. These will be removed from CHIMES in the future, but for now you should just set these to their typical values, so that this option is turned off. 

* ``Reduction_On`` *Typical value: 0*. 

* ``Reduction_N_Ions_Low`` *Typical value: 3*. 

* ``Reduction_N_Ions_Med`` *Typical value: 4*. 

* ``Reduction_N_Ions_High`` *Typical value: 5*. 

* ``macrostep_tolerance`` *Typical value: 1.0*. 

* ``min_macrostep`` *Typical value: 1.0e2*. 

* ``TemperatureEqm_Thresh`` *Typical value: 0.05*. 


Additional Config.sh Options
----------------------------

CHIMES
  Main switch to enable the CHIMES module. 

CHIMES_JEANS_SHIELDING
  The shielding length is set to the Jeans length, multiplied by the ``ChimesShieldingLengthFactor`` parameter. Column densities are then set to the cell density times this shielding length. We limit this to be no greater than 100 kpc. 

CHIMES_SOBOLEV_SHIELDING
  The shielding length is set by a Sobolev-like approximation, with L_shield = rho / grad(rho), multiplied by the ``ChimesShieldingLengthFactor`` parameter. Column densities are then set to the cell density times this shielding length. We limit this to be no greater than 100 kpc. 

CHIMES_PREENRICH_AT_START
  Applies pre-enrichment at the beginning of the simulation, within init(), rather than at some specified time. This ensures that the CHIMES element abundances are consistent with the CHIMES ion and molecule abundances read in from the initial conditions. Requires GFM_PREENRICH. 

CHIMES_PTHREADS=28
  Multi-threading option to improve work load balancing for the chemistry solver. Without this option, the code is parallelised with pure MPI, which means that particles are first distributed between MPI tasks and then, each time-step, a given MPI task will loop through its active particles and apply the chemistry and cooling. However, the computational cost of the chemistry solver can vary dramatically between particles. For example, a particle that is already close to chemical equilibrium will be much cheaper than one that is far from equilibrium. This makes it difficult to balance the work load between MPI tasks. A more efficient approach is to use multi-threading, where a single MPI task on a given node distributes its active particles between multiple threads, which are split between the CPUs on that node. This improves work load balancing because particles are dynamically allocated between the threads - in other words, a given thread will run the chemistry on one particle, and then move on to the next available particle once it is done, and so on until all active particles on that MPI task have been completed. However, since the rest of the AREPO code is not set up to use multi-threading, we cannot simply reduce the number of MPI tasks and use multiple threads in the chemistry and cooling loop. Therefore, with the CHIMES_PTHREADS option, we continue to use pure MPI for the rest of the AREPO code. Then, for the chemistry and cooling routines, all MPI tasks on a given node send the particle data of their active gas particles to one task on that node. That task then creates multiple threads to fill up all the CPUs on the node, using PTHREADS, and these threads apply the chemistry solver to those active particles. Finally, once they are completed, the resulting chemical abundances and temperatures are sent back to their original MPI tasks. 

CHIMES_INITIALISE_IN_EQM
  By default, AREPO will read in the initial ion and molecule abundances from the ICs file. However, when this option is enabled, AREPO will instead compute the initial abundances in chemical equilibrium, by evolving the chemistry of each gas cell for a long period of time. 

CHIMES_ADVECT_ABUNDANCES
  Ion and molecule abundances are advected between gas cells. 

CHIMES_DISABLE_SHIELDING_ON_EOS
  Self shielding is disabled when computing the cooling rate of the hot phase in the multiphase ISM model, which is used to calculate the effective temperature of the equation of state. 

CHIMES_DISABLE_ZCOOL_ON_EOS
  Metal cooling is disabled when computing the cooling rate of the hot phase in the multiphase ISM model, which is used to calculate the effective temperature of the equation of state. 

CHIMES_REDSHIFT_DEPENDENT_UVB
  With this option, it uses the `Haardt & Madau 2012`_ UV background, with photoionisation cross sections computed in redshift bins from z=8 to z=0. These cross section tables can be found in the ``chimes-data/HM12_cross_sections/`` directory. It will read in the tables from the redshifts just before and after the current redshift, and then interpolate the cross sections and the isotropic photon density to the current redshift. At redshifts greater than 8, it uses the redshift 8 tables but sets the isotropic photon density to zero, i.e. it switches off the UV background. 

.. _Haardt & Madau 2012: http://adsabs.harvard.edu/abs/2012ApJ...746..125H 


Authors
-------

* Alex Richings 
* Joop Schaye 
* Ben Oppenheimer 


Usage Policy and Citation
------------------------- 

If you would like to use this module for a new project, please contact the authors. 

Papers to cite: 

* Richings A. J., Schaye J., Oppenheimer B. D., 2014a, MNRAS, 440, 3349 
* Richings A. J., Schaye J., Oppenheimer B. D., 2014b, MNRAS, 442, 2780 
