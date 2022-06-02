
SGS_TURBULENCE
==============

.. warning::

	This module needs further development, do not use!

This module implements a simple model for turbulence below the grid scale. It follows the evolution of the subgrid-scale (SGS) turbulence energy and accounts for the impact on the momentum and energy of the fluid. SGS turbulence exerts an additional (isotropic) pressure and adds a stress tensor to the momentum and energy equations (similar to viscosity). Source terms describe the production of SGS turbulence energy in the turbulent cascade and its dissipation into internal energy. 

This SGS turbulence model follows the work of ``W. Schmidt 2015`` (Living Review) and is described further in the thesis of S. Jacob.


Usage
-----

This module is still under development and it has only been used for test simulations so far. The model is only suited for simple, non-cosmological hydrodynamic simulations with static Cartesian meshes. Moreover, the derivation of the underlying equations assumes 3D. 2D simulations can be used for testing, but the description of SGS turbulence energy is then not physical any more.

The master switch is SGS_TURBULENCE, but individual parts of the model have to be enabled separately. For the full model, the switches SGS_TURBULENCE, SGS_TURBULENCE_RIEMANN_PRESSURE, SGS_TURBULENCE_STRESS_TENSOR, SGS_TURBULENCE_TURBULENT_PRODUCTION and SGS_TURBULENCE_VISCOUS_DISSIPATION have to be enabled. Additionally, a way to compute the effective grid scale from the volume of a cell and a closure scheme have to be chosen. For the different choices and more details on the usage of the switches, see section 'Additional Config.sh Options'.

Additional Parameters
---------------------

Parameters that need to be set in the parameter file:

* ``MinimumSgsTSpecificEnergy`` Sets the minimum SGS turbulence specific energy (in internal code units) that is enforced in all cells.

If the eddy viscosity closure is used (SGS_TURBULENCE_EDDY_VISCOSITY):

* ``EddyViscosityParameterCnu``  This parameter controls the eddy viscosity. If it is set to zero, a default value of ``0.05`` is chosen (``Schmidt et al. 2006, Schmidt 2015``).

* ``EddyViscosityParameterCepsilon`` This parameter controls the viscous dissipation of SGS turbulence energy. If it is set to zero, a default value of ``1.58`` is chosen (``Schmidt et al. 2014`` and references therein).


Additional Config.sh Options
----------------------------

SGS_TURBULENCE_IN_ICS
  Reads the SGS turbulence specific energy from the initial conditions. This switch is not very general. It only works with ICFormat 2 and a selected set of additional fields in the ICs. If this switch is not set, the SGS turbulence energy is set to zero.

SGS_TURBULENCE_VIEW_CELLS_AS_CUBES
	Treats cells as cubes in the computation of the effective filter scale from the volume.

SGS_TURBULENCE_VIEW_CELLS_AS_SPHERES
	Treats cells as spheres in the computation of the effective filter scale from the volume.

SGS_TURBULENCE_EDDY_VISCOSITY
	Enables the eddy viscosity closure (which is the only implemented closure scheme).

SGS_TURBULENCE_CONSTANT_EDDY_VISCOSITY
	Artificially sets the eddy viscosity to a constant value.
	
SGS_TURBULENCE_RIEMANN_PRESSURE
	Enables the isotropic pressure of the SGS turbulence model. This switch requires the MHD and the RIEMANN_HLLD switches. This means that a magnetic field of 0 has to be set in the ICs.

SGS_TURBULENCE_STRESS_TENSOR
	Enables the stress tensor that is part of the SGS turbulence model. For the eddy viscosity closure, this acts like a physical viscosity (see 	``Munoz et al. 2013``). To use the reconstruction of the velocity gradients, the switches SECOND_DERIVATIVES, SLOPE_LIMIT_HESSIANS and RECONSTRUCT_GRADIENTS have to be enabled.

SGS_TURBULENCE_TURBULENT_PRODUCTION
	Enables the production of SGS turbulence through the turbulent cascade.
	
SGS_TURBULENCE_VISCOUS_DISSIPATION
	Enables viscous dissipation of SGS turbulence energy.

SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION
	Writes a log-file for the different production and dissipation mechanisms of SGS turbulence energy.
                                                         
SGS_TURBULENCE_MAKE_PRODUCTION_DISSIPATION_SPECTRUM 
	Makes spectra of the dissipation and production of SGS turbulence energy. This only works for turbulent box simulations with the driving routine of ``Bauer et al. 2012``. The additional switches AB_TURB, GAMMA and POWERSPEC_GRID are required for this.                                            

SGS_TURBULENCE_OUTPUT_SGS_PRESSURE
	Writes the SGS turbulence pressure in addition to the SGS turbulence energy into snapshots.                           

OUTPUT_SGS_T_PRESSURE_GRADIENT
	Writes the SGS turbulence pressure gradient into snapshots.

OUTPUT_DENSTROPHY
	Writes the denstrophy (1/2 * |\nabla \times \rho^{1/2} v |^2) into snapshots.
                                                                                 


Authors
-------

  * Svenja Jacob (svenja.jacob@h-its.org)


Usage Policy and Citation
-------------------------

This module is not yet usable and needs further development. Please contact Volker Springel, if you are interested in working on it.

