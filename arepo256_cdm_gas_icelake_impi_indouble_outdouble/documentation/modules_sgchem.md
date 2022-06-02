
SGCHEM
======

A set of routines for modelling gas-phase chemistry and cooling.
Several chemical networks are currently implemented:

    Network 1  -- primordial chemistry (see Glover 2015a,b and refs. therein)
    Network 4  -- Simple H chemistry from Glover & Mac Low (2007)
    Network 5  -- Network 4 + approximate CO chemistry from Nelson & Langer (1997)
    Network 15 -- H, C, O network based on Glover & Mac Low (2007), Nelson & Langer (1999);
                  see Glover & Clark (2012) for a full description

Usage
-----

TODO: Networks > 1

Network 1 (primordial chemisitry):
  Six chemical species are gives in the snapshot files, in units of relative abundance (by number) in the cell, with respect to the total number of hydrogen nuclei. (Note that with this convention, the H2 abundance is 0.5 in gas in which all of the hydrogen is in molecular form.) They are stored in the following order: H2, H+, D+, HD, He+, He++.

Additional Parameters
---------------------

TODO: not complete

SGChemConstInitAbundances     
  Can be set to 0 or 1. If set to 1 and starting with flag 0 (from initial conditions), it sets the relative abundances to the values specified below. 

SGChemInitH2Abund             
  Initial abundance of H2. 

  Important for all networks. For cosmological runs: e.g. set to 2.0e-6. Used if started from initial conditions and with parameter flag SGChemConstInitAbundances set to 1.  

SGChemInitHPAbund             
  Initial abundance of H+. 

  Important for all networks. For cosmological runs: e.g. set to 0.0004. Used if started from initial conditions and with parameter flag SGChemConstInitAbundances set to 1.

SGChemInitHePAbund           
  Initial abundance of He+. 

  Important for all networks. For cosmological runs: e.g. set to 0. Used if started from initial conditions and with parameter flag SGChemConstInitAbundances set to 1. 

SGChemInitDIIAbund            
  Initial abundance of D+.

  Important for network 1 (but needs to be set in any case). For cosmological runs: e.g. set to DeutAbund*SGChemInitHPAbund. Used if started from initial conditions and with parameter flag SGChemConstInitAbundances set to 1.

SGChemInitHDAbund             
  Initial abundance of HD.

  Important for network 1 (but needs to be set in any case). For cosmological runs: e.g. set to 0.Used if started from initial conditions and with parameter flag SGChemConstInitAbundances set to 1.

SGChemInitHeIIIAbund          
  Initial abundance of He++.

  Important for network 1 (but needs to be set in any case). For cosmological runs: e.g. set to 0. Used if started from initial conditions and with parameter flag SGChemConstInitAbundances set to 1. 

SGChemInitCPAbund             
  Initial abundance of C+.

  Important for networks 5, 7 and 15 (but needs to be set in any case). Set to 0 otherwise. Used if started from initial conditions and with parameter flag SGChemConstInitAbundances set to 1. 

SGChemInitCOAbund             
  Initial abundance of CO.

  Important for network 15 (but needs to be set in any case). Set to 0 otherwise. Used if started from initial conditions and with parameter flag SGChemConstInitAbundances set to 1. 

SGChemInitCHxAbund             
  Initial abundance of the CHx pseduospecies.

  Important for network 15 (but needs to be set in any case). For cosmological runs: e.g. set to 0. Used if started from initial conditions and with parameter flag SGChemConstInitAbundances set to 1.

SGChemInitOHxAbund             
  Initial abundance of the OHx pseduospecies.

  Important for network 15 (but needs to be set in any case). For cosmological runs: e.g. set to 0. Used if started from initial conditions and with parameter flag SGChemConstInitAbundances set to 1. 

SGChemInitHCOPAbund           
  Initial abundance of HCO+

  Important for network 15 (but needs to be set in any case). For cosmological runs: e.g. set to 0. Used if started from initial conditions and with parameter flag SGChemConstInitAbundances set to 1. 

SGChemInitMPAbund            
  Initial abundance of M+. [M is a representative low ionization potential metal used in network 15].

  Important for network 15 (but needs to be set in any case). For cosmological runs: e.g. set to 0. Used if started from initial conditions and with parameter flag SGChemConstInitAbundances set to 1. 

CarbAbund                     
  Elemental abundance (by number, relative to hydrogen) of carbon.

  Set to 0 for network 1.   

OxyAbund                    
  Elemental abundance of oxygen.

  Set to 0 for network 1.   

MAbund                       
  Elemental abundance of "M" [see above].

  Set to 0 for network 1.   

ZAtom                         
  Total metallicity, relative to solar. Influences strength of high-T atomic cooling rate.
 
  Set to 0 for network 1. 

DeutAbund                     
  Elemental abundance of deuterium.

  Important for network 1. For cosmological runs: e.g. set to 2.6e-5. 

AtomicCoolOption              
  SGCHEM contains several different models for the atomic cooling function at high temperatures. 
  For runs using network 1, set to 2.

H2OpacityOption               
  Determines what approximation is used to model H2 cooling in the optically thick limit. Only important
  in very dense gas [n > 1e8 cm^-3].

  Currently available options are:
    0 -- no opacity correction; use this if the calculation will never reach densities > 1e8 cm^-3
    1 -- power-law fit from Ripamonti & Abel (2004)
    2 -- improved power-law fit from Greif et al (2013)
    3 -- LVG approximation from Clark et al (2011)

InitDustTemp                         
  Initial dust temperature, needs to be set (but ignored e.g. in network 1). 

DustToGasRatio                         
  Dust-to-gas ratio, relative to the value for the local ISM. 

  Can be varied independently of metallicity. For runs using network 1, should generally be set to 0.

UVFieldStrength               
  Strength of UV field.

  For runs using network 1, this is the value of the mean specific intensity at the Lyman limit (13.6 eV) in
  units of [1e-21 erg s^-1 cm^-2 Hz^-1 sr^-1]; i.e. it's equivalent to the J21 parameter used by many high-z
  studies.

  For runs using chemistry networks other than network 1, this is instead the value of the UV field strength
  in Habing (1968) units. For the local ISM, an appropriate default value is 1.7 [per Draine, 1978]

LWBGType                      
  Type of Lyman-Werner background field. 

  Current options: 0 or 1. For 0, a constant field is assumed throughout the run, for option 1 the field is set to 0 before the LWBGStartRedsh and to the value specified by UVFieldStrength afterwards. 

LWBGStartRedsh                
  Start redshift for Lyman-Werner background. 

  Unused if LWBGType is not set to 1.


Additional Config.sh Options
----------------------------

SGCHEM
  Master switch used to enable SGChem framework

CHEMISTRYNETWORK
  Selects chemistry network used in the code. No default: user must explicitly select a network. Currently, valid values are 1, 4, 5 or 15, as described above. [The odd numbering is for historical reasons]

MCMA
  Enables the modified consistent multispecies advection (MCMA) algorithm described in Glover et al (2010), which ensures local conservation of total elemental abundances even when species containing multiple elements (e.g. CO) are advected. Only needed (and implemented) for chemistry network 15.

SGCHEM_VARIABLE_Z
  Allows the metallicity and dust-to-gas ratio to vary between different cells. If this option is not selected, the total gas-phase metallicity, the values of selected elemental abundances (C, O, etc.) and the dust-to-gas ratio are global variables. If this option is selected, these are converted to local variables that can vary from cell to cell. This adds some memory overhead and so is not enabled by default.

SGCHEM_VARIABLE_ISRF
  Allows the strength of the interstellar radiation field to vary spatially. If this option is not selected, the radiation field strength is a global variable (although not a constant - it is still allowed to vary over time). If this option is selected, it becomes a local variable, with a value specified by the find_and_set_local_ISRF routine in sgchem.c. At the moment, this routine assumes that you are running a simulation of a present-day Milky Way type galaxy, but it will be simple to generalize the infrastructure to support other ISRF models.

SGCHEM_VARIABLE_CRION
  Allows the strength of the cosmic ray ionization rate to vary spatially. If this option is not selected, the cosmic ray ionization rate is a global variable. If this option is selected, then it becomes a local variable. If the COSMIC_RAYS option is also selected, the value of the local variable is set by the local cosmic ray energy density. Otherwise, a simple Galactic model is used. See the routine find_and_set_local_CRION in sgchem.c for details.

CHEM_IMAGE
  Enables support for the output of chemical and dust information in the makeimage code.

SGCHEM_ACCRETION_LUMINOSITY
  Enables support for accretion luminosity heating from accreting protostars. Only useful in combination with SINK_PARTICLES. Does not account for dust, so only really useful in simulations run with CHEMISTRYNETWORK=1

SGCHEM_NO_HIGHN_DCHEM
  Disables explicit tracking of deuterium chemistry at high densities (specifically, H nuclei number densities n > 1e8 cm^-3). When selected, the HD/H2 and D+/H+ ratios are assumed to have their equilibrium values. This option is highly recommended for anyone running simulations with primordial gas chemistry (CHEMISTRYNETWORK=1) at very high densities, as it can lead to a very significant speed-up of the chemistry solver. This option has no effect when any of the other chemistry networks are used.

SGCHEM_DUMP_THERMAL_RATES
  Output individual cooling & heating rates for each cell in snapshot. This significantly increases the snapshot size and hence is disabled by default. Currently, only supports chemistry networks 4 and 5.

DEBUG_EVOLVE
DEBUG_RATE_EQ
  These options switch on additional debugging output in the chemistry solver. This can result in immense amounts of data being written to the screen, so you are advised not to turn these options on unless you already know what you are doing.

SGCHEM_NO_MOLECULES
  Switches off all molecular chemistry. Molecular species are still advected, but the chemical source and sink terms are set to zero. May be useful for testing purposes. NOTE: at present, only chemistry network 1 supports this option. Enabling it when running with any other network is useless but should be harmless.
SGCHEM_GONG17_TEST
  Our implementation of the Gong et al. (2017) simplified CO chemistry network (CHEMISTRYNETWORK 16) differs from the one presented in their paper in several respects, in one case because of an error in their implementation and in a number of other cases because we choose to adopt different rate coefficients for some reactions (e.g. H2 formation on dust grains, where we account for the dependence of gas and dust temperatures and Gong et al. do not). However, for testing purposes, it can be useful to run with exactly the same setup as in their tests. This compile-time flag enables this. It should **not** be used for production runs.

OLD_CO_SELF_SHIELDING
  Our default treatment of CO self-shielding and the shielding of CO by H2 is based on Visser et al (2009). However, we have run a number of simulations with a treatment based on the earlier work of Lee et al (1996). This older treatment is selected by this compile-time flag. It is not recommended for new simulations, but is retained in the code for backwards compatibility.

Authors
-------

  * Simon Glover (glover@uni-heidelberg.de)
  * Rowan Smith (rowanjsmith.astro@googlemail.com)
  * Paul Clark (paul.clark@astro.cf.ac.uk)
  * Tilman Hartwig (Tilman.Hartwig@ipmu.jp)  [network 1]

Acknowledgements
----------------
 
The authors thank Anna Schauer, Katharina Wollenberg, Mattis Magg, Robin Tress and Ondrej Jaura for their help with testing various aspects of the SGCHEM module. Further thanks go to Anna Schauer for contributing the time-varying LW background option and for helping with the documentation.


Usage Policy and Citation
-------------------------

Please contact the authors before using this code, so that we can
avoid any unnecessary project overlap. Co-authorship on publications
produced with the code is not required, although we are happy to
actively collaborate with users of the code that are interested in
doing so.

Papers to cite:

  * Please contact the authors for details, as the precise list will depend on which features of the code are being used.

