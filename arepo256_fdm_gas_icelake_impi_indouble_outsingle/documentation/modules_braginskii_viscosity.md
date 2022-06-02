
BRAGINSKII_VISCOSITY
====================

Braginskii MHD (Braginskii 1965) is an extension to ideal MHD which is
relevant when the plasma is weakly collisional and magnetized, i.e., when the
mean-free-path of particle collisions is much larger than the gyroradius but
much smaller than the scale-height. This leads to anisotropic heat diffusion
(see the conduction module by Rahul Kannan) and anisotropic momentum diffusion
(known as Braginskii viscosity, this module) which is directed along the local
magnetic field.

The module is explained in detail in the Braginskii code paper:

Berlok, T., Pakmor, R., and Pfrommer, C., “Braginskii viscosity on an unstructured, moving mesh accelerated with super-time-stepping”, MNRAS, vol. 491, no. 2, pp. 2919–2938, 2020. doi:10.1093/mnras/stz3115.

Usage
-----

We have included a few examples in the examples directory:

*  ``brag_decay_3d`` Decay of a velocity profile directed along the field. See section 4.1 in the Braginskii
code paper.
*  ``brag_2d`` Braginskii vicosity without the MHD solver enabled. See section 4.2 in the Braginskii
code paper.
*  ``fastwave_2d`` Decay of a magnetosonic wave due to Braginskii viscosity. Five versions with different types of timestepping. See section 4.3 in the Braginskii code paper.
*  ``rising_bubble_2d`` Buoyant bubble rising in an isothermal galaxy cluster with Braginskii viscosity.

These tests should be run with the test.sh bash script.


Additional Parameters
---------------------

* ``BragViscosityCoefficient`` A constant viscosity coefficient.

* ``BragViscosityCourant`` The Courant number, see equation 52 in the Braginskii code paper.

* ``BragViscosityMaxRKL2Stages`` The maximum number of RKL2 stages allowed when ``BRAGINSKII_RKL2_SUPER_TIME_STEPPING`` is enabled. Should be an odd integer. Very large values will make the calculation inaccurate.

* ``BragViscosityMaxSubcycles`` The maximum number of subcycling steps allowed when ``BRAGINSKII_VISCOSITY_SUBCYCLE`` is enabled.

* ``BragViscosityMaximumCoefficient`` The Spitzer viscosity coefficient depends on density and temperature as temperature^(5/2)/density. To avoid extremely small time steps, this parameter is used to set a maximum value for the viscosity coefficient when ``BRAGINSKII_SPITZER`` is enabled.


Additional Config.sh Options
----------------------------

Braginskii viscosity currently has the following extra options.

BRAGINSKII_VISCOSITY
  enables the module.

BRAGINSKII_VISCOSITY_GLOBAL_TIMESTEP
  Use the same time step for all cells in the simulation. Activate together with ``FORCE_EQUAL_TIMESTEPS``.

BRAGINSKII_RKL2_SUPER_TIME_STEPPING
  Second order accurate Runge-Kutta Legendre super timestepping. See section 3.2 in the Braginskii code paper.

BRAGINSKII_VISCOSITY_SUBCYCLE
  Subcycle the explicit call.

BRAGINSKII_SPITZER
  Use Spitzer viscosity coefficient (temperature and density dependent). Otherwise a constant value, set in the parameter file, ``BragViscosityCoefficient``, is used.

BRAGINSKII_LIMIT_PRESSURE_ANISOTROPY
  Limit the pressure anisotropy to marginal stability of firehose and mirror instabilities when evaluating the viscous flux.

BRAGINSKII_OUTPUT_PRESSURE_ANISOTROPY
  Output the (full, unlimited) pressure anisotropy (even if ``BRAGINSKII_LIMIT_PRESSURE_ANISOTROPY`` is turned on).

BRAGINSKII_OUTPUT_NUM_SUBSTEPS
  This is a diagnostic output for RKL2 or subcycling.

Authors
-------

  * Thomas Berlok (tberlok@aip.de)
  * Ruediger Pakmor (rpakmor@MPA-Garching.MPG.DE)
  * Christoph Pfrommer (cpfrommer@aip.de)


Usage Policy and Citation
-------------------------

Please contact the author before using this code for a new project.
Co-authorship on a first paper may be required.

Papers to cite:

  Berlok, T., Pakmor, R., and Pfrommer, C., “Braginskii viscosity on an unstructured, moving mesh accelerated with super-time-stepping”, MNRAS, vol. 491, no. 2, pp. 2919–2938, 2020. doi:10.1093/mnras/stz3115.
