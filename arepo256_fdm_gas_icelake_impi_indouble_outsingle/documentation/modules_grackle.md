
GRACKLE
=======

This module replaces the standard cooling implementation with functions from the Grackle (v3) library
(see https://grackle.readthedocs.org). Grackle provides:

* Two options for primordial chemistry:

  #. non-equilibrium primordial chemistry network. This is the default behavior when GRACKLE is
     uncommented in Config.sh. The network can be switched on in three stages. By default atomic
     species are enabled (HI, HII, HeI, HeII, HeIII and e). GRACKLE_H2 enables molecular hydrogen
     cooling and chemistry (HM, H2I, H2II). GRACKLE_D enables deuterium (DI, DII, HDI). The fractional
     abundances of these species are tracked and advected.

  #. tabulated H and He cooling rates calculated from Cloudy, activated with GRACKLE_TAB

* tabulated metal cooling rates calculated from Cloudy

* optional photo-heating and photo-ionization from a UV-background


Usage
-----

Grackle must be installed separately (see https://grackle.readthedocs.org) and then linked to Arepo as an 
external library. Arepo currently expects to find Grackle installed at $(HOME)/local/include and $(HOME)/local/lib. Modification of the ``GRACKLE_INCL`` and ``GRACKLE_LIBS`` variables may be required if you have installed Grackle elsewhere (please do not commit changes to this without contacting the authors), as well as taking care of the linking (e.g. by setting your ``LD_LIBRARY_PATH``).

When the code is run (with non-equilibrium chemistry), primordial abundances of the relevant species are
assumed and the correct ionization states based on internal energy provided in the ICs are obtained by 
iterating on cells until converged. This behavior can be overidden with ``GRACKLE_ABUNDANCE_IN_ICS``,
whereupon Arepo will look for the species abundances in the ICs.

When enabled, Grackle completely replaces the standard cooling implementation; this means that any calls
to functions found in cooling.c will fail, with the exception of cooling_only() which has a new definition
in grackle.c. Calling this function evolves the chemistry of the cells, carries out cooling and updates the
pressure of the cells. To reiterate, because GRACKLE replaces the standard cooling implementation it should 
be assumed that any extra modules beyond 'vanilla' Arepo that rely in some way on the cooling routines 
(e.g. star formation) will not be compatible with GRACKLE until they are appropriately modified; however, 
this should not be particularly difficult.

When Grackle is enabled, cell temperatures and fractional abundances of relevant species are outputed in 
snapshots.


Additional Parameters
---------------------

Certain parameters must be provided in the parameterfile:

* ``GrackleOn`` Master switch, 1 for on, 0 for off (in which case existing species will be
  advected as passive scalars, but no chemistry or cooling will be performed).

* ``GrackleRadiativeCooling`` Use radiative cooling

* ``GrackleMetalCooling`` Use metal cooling

* ``GrackleUVB`` UV background, 0 for off, 1 for Faucher-Giguere et al. (2009), 2 for Haardt & Madau (2012)

* ``GrackleSelfShieldingMethod`` 0 for off, see Grackle documentation for details of other options.

* ``GrackleDataFile`` path to tabulated Cloudy datafile

* ``GrackleInitialMetallicity`` If metals are not being tracked in the simulation, this value gives a default
  metallicity to be used by Grackle

* ``GracklePhotoelectricHeatingRate`` Needed if GRACKLE_PHOTOELECTRIC defined. 
  The heating rate in units of (erg cm^-3 s^-1) n^-1, where n is the total hydrogen number density.
  In other words, this is the volumetric heating rate at a hydrogen number density of n = 1 cm^-3.


Additional Config.sh Options
----------------------------

GRACKLE
  master switch, defaults to non-equilibrium atomic mode.

GRACKLE_H2
  Adds molecules.

GRACKLE_D
  Adds Deuterium.

GRACKLE_TAB
  Turn off non-equilibrium cooling and use Grackle in tabulated equilibrium mode.

GRACKLE_ABUNDANCE_IN_ICS
  Use the species that are in the ICs.

GRACKLE_TEMPERATURE_FLOOR
  Obey the floor imposed by MinGasTemp or MinEgySpec.

GRACKLE_VERBOSE
  Turn on internal print statements from the Grackle library. This will likely produce a lot of output to stdout from each task. Grackle's error and warning messages will still be output to stdout and stderr even with this option switched off.

GRACKLE_UNCONVERGED_IGNORE
  On startup (without ``GRACKLE_ABUNDANCE_IN_ICS``), ordinarily we will iteratively attempt to find equilibrium abundances given a particular internal energy for each cell. If abundances fail to converge within a certain number of iterations (currently 10,000 is hard-coded) the code will terminate. Turning on this option will instead allow the code to continue even if the abundances are not converged.


Authors
-------

* Benjamin W. Keller (benjamin.keller@uni-heidelberg.de)
* Matthew C. Smith (msmith@mpia.de)


Usage Policy and Citation
-------------------------

Please contact the authors (of the Arepo module) before carrying out any development work so
coding conflicts can be avoided.

If Grackle is used in a publication, the authors of the library request the following recognition:

If you used the Grackle library in your work, please cite it as “the Grackle chemistry and cooling library 
(Smith et al. 2017).” Also, please add a footnote to 
https://grackle.readthedocs.org/.
