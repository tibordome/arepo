
MRT
===

Moment based radiative transfer module using the M1 closure relation. Handles both single scattering UV radiation and multisctattering IR radiation fields. 


Usage
-----

The main switch to turn on this option is MRT. The internal units for the UV frequency bins is 1e63 photons while for the IR bin is just the internal energy unit of the code. 


Additional Parameters
---------------------

``RTNUmSubCycles`` Number of RT subcycles for every hydro step. Ideally use a number that is a poer of 2. WOuld no recommend going above 32 subcycles.

``HLLEFile`` The file that contains the eingenvalues of the riemann solution. Needs to be specified if the MRT_RIEMANN_HLLE or MRT_RIEMANN_HLLE_NEW options are used.

``SpectrumTablePath`` Path to the file that contains the input spectrum. Needs to be specified if MRT_LOCAL_FEDDBACK or MRT_STARS option is used.


Additional Config.sh Options
----------------------------

MRT
  Master Switch.

MRT_COMOVING
  Solve the RT equations in the comoving reference frame - needs more testing.

MRT_INIT_IONIZATION
  Set the ionization states to the equilibrium values at the begining of the simulation. If not specified, the code assumes fully neutral gas irrespective of the temperature.

MRT_OUTPUT_FLUX
  Output the photon flux vector.

MRT_TIME_EXTRAPOLATION
  Include the first order time prediction step - needed for second order convergence.

MRT_FLUX_EXTRAPOLATION
  Time extrapolate the photon flux as well - induces too much noise in the solution. Do not use for now.

MRT_COOLING_HEATING
  Include cooling and heating.

MRT_RADIATION_PRESSURE
  Include radition pressure as a soucre term for the momentum conservation equation.

MRT_INCLUDE_HE
  Include Helium chemistry, heating and cooling.

MRT_LSF_GRADIENTS
  Use the least square fit gradient estimates, should always be turned on.

MRT_RIEMANN_ROSUNOV
  Use the Rosunov (GLF) riemann solver.

MRT_RIEMANN_ROSUNOV_NEW
  Use the Rosunov (GLF) riemann solver (accounts mesh motion via advection step).

MRT_RIEMANN_HLLE
  Use the Harten-Lax-van Leer flux function - Preferred riemann solver.

MRT_RIEMANN_HLLE_NEW
  Use the Harten-Lax-van Leer flux function (accounts mesh motion via advection step)

MRT_MULTI_FREQUENCY
  Multi Frequency radiative transfer. Number of frequency bins specified in allvars.h. The frequency range of each frequency bin specified in src/MRT/RT_init.c

MRT_CHEMISTRY_PS2009
  Petkova & Springel 2009 Chemistry. Unstable for large variations in radiation field intensity.

MRT_CHEMISTRY_PS2011
  Petkova & Springel 2011 Chemistry. Unstable for large	variations in radiation	field intensity. Only compatible with H chemistry and single frequency RT.

MRT_COUPLED_THERMOCHEMISTRY
  This option estimates the change in quatities using the PS2009 method, if the change chnage is smaller than 10%, then it keeps the solution, else it solves the coupled chemistry and cooling equations using the CVODE solver from the SUNDIALS library. Needs to be liked with the SUNDIALS library.

MRT_NO_OTSA
  Do not apply On The Spot Approximation (OTSA), i.e., account for the recombination radiation.

MRT_NOCOLLISION_IONIZATION
  No Collisional ionisation. Useful in some test problems.

MRT_SLOWLIGHT
  Reduce the speed of light to the value of CLIGHT set in allvars.h.

MRT_CONSTANT_KAPPA
  Constant speed of light set in RT.c. Mainly used for testing purposes.

MRT_IR
  Include IR (multi scattering) radiative transfer. The dust opacities need to be set in RT.c if MRT_IR_GRAIN_KAPPA is not defined.
  
MRT_IR_ONLY_CHEMISTRY
  Only do absorption terms and not the cooling terms - assume there is no absorption of energy, only the flux is absorbed.

MRT_IR_LTE
  Do gas heating and cooling according to single fluid gas-dust-IR radiation LTE assumption.

MRT_IR_LTE_SEMI_IMPLICIT
  Follows  Rosdahl+15 iterative approach to solve the coupled dust/gas-IR radiation temperatures. Can lead to large number of interative steps (sometimes infinite)

MRT_IR_LTE_GSL
  Use the GSL provided Implicit Bulirsch-Stoer method of Bader and Deuflhard to solve the  dust/gas-IR radiation temperatures (Prefered Method).

MRT_IR_PHOTON_TRAPPING
  Trap unresolved IR photons - not compatible with gradient extrapolations (hardcoded - no need to turn off time and space extrapolations).

MRT_IR_GRAIN_KAPPA
  Use grain opacities calculated from Draine & Lee 1984, Laor & Draine 1993, requires local dust-to-gas ratio from GFM_DUST.

MRT_UV_ONLY_DUST
  Let UV radiation only interact with dust, no H-He thermochemistry.

MRT_NO_UV
  Do not include UV RT (mainly for testing purposes).

MRT_SETUP_SPECIAL_BOUNDARIES
  Setup special boundary conditions (mainly for testing purposes). Need to edit the src/MRT/RT_setup.c file.

MRT_LEVITATION_TEST
  Enable the appropriate options for the levitation test (Section 4.10 of Kannan et al. 2018).

MRT_SOURCES=0
  Enable source treatment for GFM stellar particles and black holes

MRT_REDUCE_OUTPUT
  Reduce output from MRT routine.

MRT_EQUIL_CHEM_COOL
  Use the fiducial equilibrium chemistry + cooling of Arepo. Not as accurate as the non-equilibrium thermochemistry.

MRT_MOLECULAR_COOLING
  Include low temprature molecular cooling, fit from grackle.

MRT_PHOTOELECTRIC_HEATING
  Include photo electric heating, needs an additional frequency bin in the 5.6-11.2 eV range.

MRT_METAL_COOLING
  Add tabulated metal cooling. Needs GFM and GFM_COOLING_METAL to be turned on.

MRT_UVB
  Add contribution from UVB - only z=0 for now.

MRT_UPDATE_AT_END_OF_STEP
  Update the hydro primitive varibales only at the end of hydro loop. This makes the RT more compatible with the hydro, but can lead to unphysical oscillations if SF is turned on.

MRT_STARS
  Include ionizing photons from GFM stellar particles.

MRT_STARS_EXACT_NGB
  Do not do ngb mass-weighting.

MRT_BH
  Include photons from black hole particles.

MRT_BH_EXACT_NGB
  Do not do ngb mass-weighting.

MRT_BH_PULSED
  Pulsed radiation injection, inorder to mimick the quasar duty cycles.

MRT_BH_UV_INJECTION
  Inject photons in UV bin(s) for BH particle only.

MRT_BH_IR_INJECTION
  Inject photons in IR bin for BH particle only.

MRT_BH_BIPOLAR
  Inject photons in a bipolar jet instead of in a isotropic manner.

MRT_BH_BIPOLAR_SET_FLUX
  Set photon flux for bipolar injection.

MRT_BH_OMEGA_WEIGHT
  Weight photon injection by solid angle subtended by the cell.

MRT_LOCAL_FEEDBACK
  Main switch to couple to local feedback module of Chrsitine Simpson.

MRT_CHEM_SG
  Switch to couple to SGchem module of Simon Glover. 



Authors
-------

  * Rahul Kannan (rahul.kannan@cfa.harvard.edu)


Usage Policy and Citation
-------------------------

Papers to cite:
  * Kannan et al. 2018, ArXiv 1804.01987
  
