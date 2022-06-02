
Initial Conditions
==================

Arepo is a multi-purpose code that supports a number of different types of simulations.

For **examples of actual simulation setups** see the ``examples/`` directory.

Below are comments and tips for setting up simulations of various common types.


Cosmological Simulations
------------------------

One of the main simulation types in Arepo are simulations on an expanding spacetime, where the coordinates are treated as comoving
coordinates.  This mode is actived whenever the parameter flag ``ComovingIntegrationOn`` is set to "1" in the parameter file. In this
case, the code models an expanding spacetime described by the Friedmann equations.


Periodic Cosmological Volumes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One of the standard types of cosmological simulations are cosmological volume simulations. These simulations have a uniform mass resolution
and are set up by using cosmological perturbation theory, converted to small initial position and velocity perturbations in a periodic box
with fixed comoving extent. These simulations can be done with gravity only, in which case the simulation particles are all of a single type
(type > 0) and gravity is recommended to be calculated using a TreePM algorithm. This means that the compile time flags ``SELFGRAVITY`` and
``PMGRID`` are set. The example ``cosmo_box_gravity_only_3d`` is a very small representative of such a simulation.

In these cosmological simulations it is also possible to include gas. This can be done, for example, by reading in gravity only initial
conditions and setting the compile-time option ``GENERATE_GAS_IN_ICS``, which will split every input particle in a
gas cell and a dark matter particle, with the mass ratio following the cosmic average mass and baryon fractions, ``Omega0`` and
``OmegaBaryon``, respectively, which are defined in the parameter file.  Including gas in the simulation, it is also useful to use mesh
regularization options as well as refinement and derefinement. On top of this, it is also possible to use prescriptions for radiative
cooling and star formation using the respective compile-time flags ``COOLING`` and ``USE_SFR``. The example
``cosmo_box_star_formation_3d`` is using such a setup.

To generate such initial conditions, there are many public codes. Easiest recommendations:

* `N-GenIC <https://www.h-its.org/2014/11/05/ngenic-code/>`_ (MPI parallel, Zeldovich approximation).
* `S-GenIC <https://github.com/sbird/S-GenIC>`_ is an extended version of N-GenIC, written by Simeon Bird (OpenMP parallel, 2LPT).
* `MUSIC <https://www.h-its.org/2014/11/05/ngenic-code/>`_ is a more general IC code, which can also generate zooms, written by Oliver Hahn. It includes AREPO format output support, written by Dylan Nelson. (OpenMP parallel, 2LPT).


Cosmological Zoom (Multi-Mass)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Instead of a fully periodic volume, the Lagrangian region belonging to a single final object can be sampled with significantly higher resolution than the
rest of the cosmological volume. In this way, it is possible to accurately resolve the formation of an individual object, while still
correctly modeling the large-scale environmental effects. These simulations can be run, as above, in gravity only mode or including
gas physics. A gravity only example is given as ``cosmo_zoom_gravity_only_3d`` in the examples (note however that the
ICs need to be created separately in this example, as they are slightly too large to include them into the main code repository).
Technically, such a zoom simulation needs slightly different configuration in the particle-mesh algorithm. In particular, it is
often useful to place a second particle-mesh layer which calculates forces only in the high-resolution region.  This can be triggered by
the compile-time flag ``PLACEHIGHRESREGION``.  In addition, the flag ``PM_ZOOM_OPTIMIZED`` will change the employed particle-mesh algorithm
to one that is specifically optimized for this type of simulations.

To generate such initial conditions, there are several codes. Easiest recommendations:

* `MUSIC <https://www.h-its.org/2014/11/05/ngenic-code/>`_ is a more general IC code, which can also generate zooms, written by Oliver Hahn.
* `P-Resim <https://www.tng-project.org/>`_ is an extension of N-Genic, written by Volker Springel. This code is not public, contact the author for usage. An updated version which produces HDF5 ICs based on Illustris[TNG] halos is available, contact Dylan Nelson.




Non-Cosmological Simulations
----------------------------


Apart from comoving integration, Arepo can also handle an ordinary Newtonian spacetime by choosing the parameter option
``ComovingIntegrationOn 0``.  While cosmological simulations usually assume periodic boundary conditions, simulations in 
Newtonian space can also have reflective or inflow/outflow boundary conditions.


Isolated Galaxy or Merger
^^^^^^^^^^^^^^^^^^^^^^^^^

One example of a non-cosmological galaxy simulation is an isolated, self-gravitating object, such as in the
``isolated_galaxy_collisionless_3d`` example. This case only contains collisionless particles, namely the dark matter 
and stellar component of a disk galaxy, and their gravitational interactions are calculated with a tree algorithm only, 
assuming non-periodic forces. Mergers between such idealized galaxies are shown in the ``galaxy_merger_star_formation_3d``
example, which also includes gas and hydrodynaimcs. Since AREPO requires a finite extent of the simulation box as well 
as positive density at every point, isolated galaxy initial conditions must be prepared somewhat carefully.

(0) Generating the initial conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Note that the ``makenewdisk`` of old has been updated as ``GALIC``. This is now nicely documented,
described in a `paper <http://adsabs.harvard.edu/abs/2014MNRAS.444...62Y>`_, and
`publicly available for download <https://www.h-its.org/tap-software-de/galic-code-2/>`_.
The GALIC code implements an iterative method for the construction of N-body galaxy models in
collisionless equilibrium, and can be used to set up isolated galactic disks, or, for example,
idealized binary merger simulations.

(1) Preparing the initial conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A major conceptual difference between GADGET and AREPO is that the
former can more easily deal with vacuum boundary conditions, while in
AREPO there always needs to be a well defined volume in which the gas
calculation takes place.

My suggestion for dealing with this is to place
the merger into a (very) large box, and to fill the volume with a coarse
“background grid” of cells. The box needs to be large enough such that
the hydrodynamics is not affected by the presence of these outer
boundaries, which is the case if they are far enough away (for example
enclosing at least the virial radii of the halos and a buffer zone). In
this case it should also not matter if the hydrodynamic boundary
conditions for the box are assumed to be periodic or reflective, which
is both possible in the code. Gravity should in any case be treated as
if there was no box at all. This means in particular that dark matter
and star particles are allowed to reach locations outside the fiducial
gas-box.

There is a special procedure implemented in AREPO to transform an
existing SPH galaxy collision initial conditions to one you can use
with AREPO. To this end, you set the ::

    ADDBACKGROUNDGRID=N

option in the Makefile. You give it a value N equal to the resolution of
the desired background grid, for example 32. The largest volume of a
cell in the initial conditions will then be ``(BoxSize/ADDBACKGROUNDGRID)^3``.

In the parameter file, you point
``InitCondFile`` to the SPH initial conditions, and all other parameters
you can set to the values you want to use in the merger calculation. In
particular, you should set BoxSize, which will determine the extension
of the fiducial gas-box used by AREPO. Note that the procedure to
transform the initial conditions assumes that your SPH calculation has
its center at ~(0/0/0), and it will first shift all particles by
``BoxSize/2`` along each coordinate axis. In your AREPO calculaion, the
merger remnant will hence be found at ~(``BoxSize/2``, ``BoxSize/2``,
``BoxSize/2``).

Then you start AREPO in the usual way. After a short time,
the code should terminate, after having created a new set of initial
conditions with identical name(s) as ``InitCondFile`` except for an appended
“with-grid”.

(2) Running the simulation code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, you need to disable ``ADDBACKGROUNDGRID`` again and recompile the
code, such that it is now ready to evolve the new initial conditions. In
the parameterfile, change the "``InitCondFile``" variable to now
point to the file that was created in step (1), then start the run.



Three-dimensional simulations
-----------------------------

Simulations may also be completely hydrodynamical, with no
gravitational forces involved. One such example is the setup of the
``Noh_3d`` problem.  This example is an important test problem as it
has an analytic solution, thus it is suitable for code
verification. However, like for the more complex examples, we also
here check against outcomes from previous code versions in order to
ensure consistency.


Two-dimensional simulations
---------------------------

Two-dimensional simulations are often used as simplified examples
which are significantly cheaper to run than their three-dimensional
counterparts, and thus also very useful as test problems. Note that
while the one-dimensional simulation code is in substantive parts
detached from the rest, two-dimensional and three-dimensional
simulations largely use the same routines. Thus these routines can be
efficiently tested with 2d test problems.

One example of such a test is ``gresho_2d``, a stationary vortex
problem for which the pressure gradient balances the centrifugal
forces. The same is true for the Yee vortex, ``yee_2d``, which has the
advantage of being smooth in the hydrodynamic properties and therefore
is better suited to determine the convergence order of the code. Both
these tests are sensitive to accurate hydrodynamical modelling and
gradient estimates.

Other 2d test problems are ``noh_2d``, a converging gas flow in 2d
which introduces a strong shock and is a challenging problem for the
Riemann solver, as well as ``current_sheet_2d``, an MHD test probing
numerical reconnection properties of the code.


One-dimensional simulations
---------------------------

Most of the one-dimensional simulations are test problems for
particular solvers or numerical methods. Note that Arepo for the case
of 1d problems is **NOT** MPI parallel, i.e. all the following
examples can only be calculated using one MPI rank, although 1D 
problems are fully OpenMP parallelized. Furthermore, refinement and 
derefinement are not supported in 1D.

One of the most basic one-dimensional examples is a simple linear wave
propagation test as in the example ``wave_1d``. Such a test is
well-suited to test the convergence order of a scheme. For such
applications it is important to ensure that both input/output and
calculation are done in double precision by using the compile-time
flags ``INPUT_IN_DOUBLEPRECISION``, ``OUTPUT_IN_DOUBLEPRECISION`` and
``DOUBLEPRECISION=1``.

Another very important basic test is a 1d shocktube problem. Such a
shocktube is calculated in the example ``shocktube_1d``, where the
time evolution can be compared to the exact solution. This tests both
the Riemann solver as well as the gradient estimates and associated
limiters. Similarly, ``mhd_shocktube_1d`` tests the time evolution of
an MHD Riemann problem. Note that the verification in the latter case
is only against a high-resolution simulation.

One example of a more complex one-dimensional test problem is the
example ``interacting_blastwaves_1d``. This is checked against a
high-resolution solution calculated on a fixed grid.
