NON_IDEAL_MHD
=============

.. warning::

  The documentation of this module is work in progress and still incomplete

Implementation of the non-ideal MHD terms in Arepo. Currently only ohmic resistivity is available.

Usage
-----

The ``NON_IDEAL_MHD`` master switch enables the module, which then allows to choose between different non-ideal MHD terms and time integration strategies. The module works with both the default **MHD_POWELL** and the **MHD_CT** schemes. More details on the usage of the module can be found in the descritpion of the Config.sh options below. Please note that currently only ohmic resistivity has been fully implemented in the code.

Additional Parameters
---------------------

* ``OhmicDiffusionCoefficient`` is the value, in internal code units, of the gas ohmic resistivity (assumed to be constant). It is requested by **OHMIC_DIFFUSION** or **IMPLICIT_OHIMC_DIFFUSION**.

* ``AmbipolarDiffusionCoefficient`` is the value, in internal code units, of the ambipolar diffusion coefficient (assumed to be constant). It is requested by **AMBIPOLAR_DIFFUSION**.

* ``NonidealMHDTimelimit`` is the fraction of the Alfven crossing time of each cell to which the time step has to be limited. It is requested by **NON_IDEAL_MHD_EXPLICIT_LIMIT_TIMESTEP**.

Additional Config.sh Options
----------------------------

NON_IDEAL_MHD
  Master switch needed to enable the non-ideal MHD terms in the code. See the switches below for the choice of the type of non-ideal MHD effect and of its time integration strategy.

OHMIC_DIFFUSION
  Enables the treatment of ohimic resistivity in the code. Currently, only constant resistivity is supported. The switch works with either the **MHD_POWELL** or **MHD_CT** options.

IMPLICIT_OHMIC_DIFFUSION
  If enabled, integrates the resistive ohmic term with an implicit times integration strategy. Otherwise an explicit timestepping techinques, which simply adds the resistive fluxes to the ideal MHD fluxes, is used. The convergence order of the implicit time integration scheme can be controlled with the option **OHM_CRANCK_NICHOLSON**.

OHM_CRANK_NICHOLSON
  If enabled, the code resorts to a second-order Crank-Nicholson scheme for the (implicit) time integration of the resistive terms rather than a simple but more robust backward Euler scheme. Requires **IMPLICIT_OHMIC_DIFFUSION**.

ONLY_OHMIC_DIFFUSION
  The switch is required for test problems modeling only the resistive ohmic terms (such as the diffusion of a magnetic field). If enabled, ideal MHD fluxes are zeroed and only the resistive ohmic fluxes are considered instead.

OHMIC_HEATING
  Only if **IMPLICIT_OHMIC_DIFFUSION** is disabled. Allows for a separate treatment of the Joule heating term for the explicit time integration case. Rather than adding this to the total energy flux, the thermal energy of the gas is updated at the end of the time step instead, in the same way as it is done for implicit time integrations schemes.

NON_IDEAL_MHD_EXPLICIT_LIMIT_TIMESTEP
  Limits the time step of each cell to a fraction ``NonidealMHDTimelimit`` of the Alfven crossing time for the explicit treatment of time integration of the non-ideal MHD terms.

AMBIPOLAR_DIFFUSION
  Enables the treatment of ohimic resistivity in the code. Currently under development.

ONLY_AMBIPOLAR_DIFFUSION
  The switch is required for test problems modeling only the ambipolar diffusion terms. If enabled, ideal MHD fluxes are zeroed and only the ambipolar diffusion fluxes fluxes are considered instead. Currently under development.


Authors
-------

* Federico Marinacci (federico.marinacci@cfa.harvard.edu)


Usage Policy and Citation
-------------------------

Please contact the author before using this code for a new project. Co-authorship on a first paper may be required.

Papers to cite:

  * `Marinacci et al. 2018, MNRAS, 442, 43 <http://adsabs.harvard.edu/abs/2014MNRAS.442...43M>`_
