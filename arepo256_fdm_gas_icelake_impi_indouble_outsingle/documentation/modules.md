
Code Modules List
=================

The AREPO code contains a number of physics modules that take the
form of extensions of the code for specific science applications. 
Examples include MHD, non-equilibrium chemistry, tracer particles,
etc.

Usually, these modules have been developed by a single person or a
small group of people, and are proprietary to them. If other people
want to use these code parts for their own project, they have to
check with the corresponding authors first. Depending on the policy
they have defined, such a use may (or may not) be possible and may 
involve a request for co-authorship in publications. 

In any case, **please check this list for a first info about the
usage policies** of different, clearly defined modules. In case of
doubt, please feel free to discuss with the corresponding module
authors. 

**When starting to develop a new module, please add a short 
note here about its purpose and its 'in progress' status. When 
finalizing a module, possibly in conjunction with a first publication, 
please update and expand its documentation as appropriate.**


-----


Geometry
--------

AMR
^^^
Use Cartesian-based adaptive mesh refinement instead of the usual Voronoi approach for the discretization of the hydrodynamics.

* Documentation: :doc:`modules_amr`
* Authors: Andreas Bauer, Ruediger Pakmor
* Status: in progress, experimental only.


-----

Solvers and Equation Systems
----------------------------

MHD
^^^
Ideal MHD implementation. This implementation uses cell-centered magnetic fields and employs the Powell/Dedner scheme for divergence cleaning.

* Documentation: :doc:`modules_mhd`
* Authors: Ruediger Pakmor
* Status: working, production ready.


MHD_CT
^^^^^^
Constrained transport Ideal MHD implementation. More accurate and robust (and experimental) than the Powell/Dedner scheme.

* Documentation: :doc:`modules_mhd_ct`
* Authors: Philip Mocz
* Status: ?


GLOBAL_VISCOSITY
^^^^^^^^^^^^^^^^
Solve Navier-Stokes instead of Euler equations.

* Documentation: TODO.
* Authors: Diego Munoz
* Status: working, production ready?


SPECIAL_RELATIVITY
^^^^^^^^^^^^^^^^^^
Special relativistic hydro solver for Arepo.

* Documentation: :doc:`modules_special_relativity`
* Authors: Andreas Bauswein
* Status: in progress development only.


DG
^^
Discontinuous galerkin method for hydrodynamics instead of usual linear gradient approach.

* Documentation: TODO.
* Authors: Thomas Guillet
* Status: in progress development only?


NON_IDEAL_MHD
^^^^^^^^^^^^^
Enables non-ideal MHD module, which can be used in conjunction with either the Powell divergence cleaning or the constrained transport MHD schemes. Currently, only the ohmic resistive term is available.

* Documentation: :doc:`modules_non_ideal_mhd`
* Authors: Federico Marinacci
* Status: ?


-----

Radiation Transport
-------------------

RT
^^
Radiative transfer module (cone method).

* Documentation: TODO.
* Authors: Margarita Petkova?
* Status: ?


OTVET
^^^^^
Radiative transfer module (optically thin variable eddington tensor method).

* Documentation: TODO.
* Authors: Margarita Petkova?
* Status: ?


FLD
^^^
Radiative transfer module (flux limited diffusion method).

* Documentation: TODO.
* Authors: Margarita Petkova?
* Status: ?

TREECOLV2
^^^^^^^^^
The TreeCol module (obtain column density maps of sky during treewalk)

* Documentation: :doc:`modules_treecolv2`
* Authors: Paul Clark and Simon Glover
* Status: Production ready.

MONOTONE_CONDUCTION
^^^^^^^^^^^^^^^^^^^

* Documentation: :doc:`modules_monotone_conduction`
* Authors: Rahul Kannan
* Status: Production ready.

BRAGINSKII_VISCOSITY
^^^^^^^^^^^^^^^^^^^^

Enables Braginskii viscosity, a special type of viscosity which is relevant
for weakly-collisional and magnetized plasmas, e.g. the intracluster medium
in galaxy clusters.

* Documentation: :doc:`modules_braginskii_viscosity`
* Authors: Thomas Berlok, Ruediger Pakmor and Christoph Pfrommer
* Status: Basic functionality working, comoving integration is work in progress.

MRT
^^^

Moment based radiative transfer using the M1 closure relation. 

* Documentation: :doc:`modules_mrt`
* Authors: Rahul Kannan
* Status: Production ready.

-----


Physical Models (Galaxy/Cosmological)
-------------------------------------

GFM
^^^
The 'galaxy formation model' with several components originally developed for the Illustris simulation. Includes galactic-scale winds, stellar evolution and enrichment.

* Documentation: :doc:`modules_gfm`
* Authors: Mark Vogelsberger, Paul Torrey, Shy Genel (and others)
* Status: working, production ready.


GFM_CHEMTAGS
^^^^^^^^^^^^
Chemical tagging to follow in more detail the origin of mass/metals from different processes.

* Documentation: :doc:`modules_gfm_chemtags`
* Authors: Jill Naiman
* Status: working, production ready.


GFM_DUST
^^^^^^^^
Passive scalar model for the formation and evolution of dust.

* Documentation: :doc:`modules_gfm_dust`
* Authors: Ryan McKinnon
* Status: working, production ready?


DUST_LIVE
^^^^^^^^^
Particle-based framework for dust grain dynamics and size evolution.

* Documentation: :doc:`modules_dust_live`
* Authors: Ryan McKinnon
* Status: in progress development only.


BLACK_HOLES
^^^^^^^^^^^
Treatment of supermassive blackholes, including their seeding, accretion, and feedback. 

* Documentation: :doc:`modules_black_holes`
* Authors: Debora Sijacki, Mark Vogelsberger (and others)
* Status: working, production ready.


GRACKLE
^^^^^^^
Replaces the standard cooling implementation with functions from the Grackle library.

* Documentation: :doc:`modules_grackle`
* Authors: Matthew Smith
* Status: in progress development only.


COSMIC_RAYS
^^^^^^^^^^^
Cosmic Ray implementation.

* Documentation: :doc:`modules_cosmic_rays`
* Authors: Ruediger Pakmor (and others)
* Status: usable, under active development.


FM_SFR
^^^^^^
Stellar feedback and ISM physics module of FM, intended for 'resolved' ISM resolution runs.

* Documentation: :doc:`modules_fm_sfr`
* Authors: Federico Marinacci (and others)
* Status: in progress development only.


SIDM
^^^^
Self-interacting dark matter model.

* Documentation: TODO.
* Authors: Mark Vogelsberger
* Status: working, production ready?


LOCAL_FEEDBACK
^^^^^^^^^^^^^^
Localized supernova feedback model.

* Documentation: :doc:`modules_local_feedback`
* Authors: Christine Simpson
* Status: in progress development only.


SFR_MCS
^^^^^^
Star formation and stellar feedback routines for high resolution galaxy sims.

* Documentation: :doc:`modules_sfr_mcs`
* Authors: Matthew Smith
* Status: usable, under active development.


-----

Physical Models (Other)
-----------------------

SGCHEM
^^^^^^
A set of routines for modelling gas-phase chemistry and cooling. Currently implements a very simple H, C, O network.

* Documentation: :doc:`modules_sgchem`
* Authors: Simon Glover (and others)
* Status: usable, under active development.


CIRCUMSTELLAR
^^^^^^^^^^^^^
Models for circumstellar disks.

* Documentation: TODO.
* Authors: Diego Munoz
* Status: ?


EOS_OPAL
^^^^^^^^
OPAL equation of state for stellar applications.

* Documentation: :doc:`modules_eos_opal`
* Authors: Sebastian Ohlmann
* Status: working, production ready.


-----


Analysis
--------

TRACER_*
^^^^^^^^
Implementation of tracer particles that allow to follow the flow of gas (as well as stars/BHs if relevant) in a 'Lagrangian' way. Also the passive scalar field.

* Documentation: :doc:`modules_tracer`
* Authors: Shy Genel, Dylan Nelson (and others)
* Status: working, production ready.


SHOCK_FINDER 
^^^^^^^^^^^^
Hydrodynamic shock finder implementation. 

* Documentation: :doc:`modules_shock_finder`
* Authors: Kevin Schaal 
* Status: working, production ready.



-----


These are what seem to be 'master' module switches, entirely lacking documentation at present:

Unknown
-------

* AB_TURB
* EOS_DEGENERATE
* NUCLEAR_NETWORK
* SPECIAL_BOUNDARY
* WINDTUNNEL
* RADCOOL
* TGSET
* SPIRAL
* SNE_FEEDBACK
* MODGRAV
* CONDUCTION
