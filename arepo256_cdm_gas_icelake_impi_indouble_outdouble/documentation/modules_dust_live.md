
DUST_LIVE
============

Particle-based framework treating dust as live simulation particles.  Dust
particles represent ensembles of physical grains and are subject to dynamical
forces.  In principle, dust motion can be decoupled from gas motion.
Optionally, DUST_LIVE allows dust particles to represent grains of different
sizes using a discretized grain size distribution that can evolve in time.
Grain size evolution models dust grain collisions with gas-phase metals and
other dust grains and is implemented in a subgrid fashion.

Usage
-----

By default, simulation dust particles move dynamically according to gravity and
an aerodynamic drag force.  The drag force can be solved with a semi-implicit
time integrator, which alleviates the need to resolve the (possibly short) drag
stopping time-scale.  Instead, it makes use of the analytic behavior of drag
over time intervals where the stopping time-scale is roughly constant.

Optionally, dust particles can be assigned grain size distributions that track
the number of grains of different sizes.  Unlike a stellar IMF, the grain size
distribution can vary in time and from dust particle to dust particle.  For a
dust particle, the grain size distribution in bin :math:`i` is parameterized as

.. math:: \frac{\mathrm{d}n(a, t)}{\mathrm{d}a} = \frac{N_i(t)}{a_{i+1}^\mathrm{e} - a_{i}^\mathrm{e}} + s_i(t) \times (a - a_{i}^\mathrm{c})

where :math:`\mathrm{d} n(a, t) / \mathrm{d}a \times \Delta a` is the number of
grains in the size interval :math:`[a, a + \Delta a]` at time :math:`t` on a
dust particle, :math:`N_i(t)` and :math:`s_i(t)` are the number of grains and
slope in bin :math:`i`, :math:`a_{i}^\mathrm{e}` and
:math:`a_{i+1}^\mathrm{e}` are the bin edges, and :math:`a_{i}^\mathrm{c} =
(a_{i}^\mathrm{e} + a_{i+1}^\mathrm{e})/2` is the bin center.

Star particles stochastically spawn dust particles during ISM enrichment
according to tabulated dust yields.  This approach creates the correct amount
of dust in expectation, and the desired mass of simulation dust particles can
be adjusted to balance computational cost and resolution.

Depending on which configuration options are on, various snapshot fields are
written.  These are summarized below.

===================  =================  ====================================  ==============================================================================================
HDF5 Dataset         Physical Factor    Units                                 Description
===================  =================  ====================================  ==============================================================================================
Coordinates          :math:`a / h`      UnitLength                            Same as in snapshots.
Velocities           :math:`\sqrt{a}`   UnitVelocity                          Same as in snapshots.
Masses               :math:`1 / h`      UnitMass                              Same as in snapshots.
Dust_GasHsml         :math:`a / h`      UnitLength                            Smoothing length enclosing gas neighbors around a dust particle.
Dust_GasDensity      :math:`h^2 / a^3`  UnitMass/UnitLength^3                 Local gas density around a dust particle.
Dust_NumGrains       :math:`1`          none                                  Array indicating number of grains in each grain size bin.
Dust_BinSlopes       :math:`1`          1/micron^2                            Array indicating bin slope in each grain size bin.
Dust_MetalFractions  :math:`1`          none                                  Array indicating which fraction of dust particle mass comes from different chemical elements.
Dust_DustHsml        :math:`a / h`      UnitLength                            Smoothing length enclosing dust neighbors around a dust particle.
Dust_DustDensity     :math:`h^2 / a^3`  UnitMass/UnitLength^3                 Local dust density around a dust particle.
Dust_NgbID           none               none                                  ID of closest dust neighbor that a dust particle could derefine into.
Dust_IsWind          none               none                                  Flag indicating whether a dust particle is in wind phase (1) or not (0).
Dust_WindTimeLeft    :math:`1 / h`      UnitLength/UnitVelocity               Amount of time a dust particle has left in wind phase.
===================  =================  ====================================  ==============================================================================================


Additional Parameters
---------------------

The DUST_LIVE module requires several parameters, depending on which
configuration options are active.

* ``DesNumNgbDust`` Number of neighbors (gas or dust) around dust particles
  used in kernel estimates.
  *Default value: 64*.

* ``MaxNumNgbDeviationDust`` Tolerance in number of neighbors for searches
  around dust particles.
  *Default value: varies*.

* ``MinGrainSize`` Minimum size for grain size distribution, in microns.
  *Default value: 0.001*.

* ``MaxGrainSize`` Maximum size for grain size distribution, in microns.
  *Default value: 1.0*.

* ``GrainDensity`` Physical density of a solid dust grain, in g cm^-3.
  *Default value: 2.4*.

* ``MaxBinFracMassChg`` Fractional amount by which a grain size distribution
  bin is allowed to change in mass during a (possibly sub-cycled) time-step.
  *Default value: 0.1*.

* ``GrainDataPath`` Path to HDF5 file containing various initialization data
  (dust nucleosynthetic yields, initial grain size distributions, etc.).
  *Default value: varies*.

* ``DustTargetFrac`` Fraction of initial mass of a star particle that converts
  to a dust particle when stochastic condition is met.
  *Default value: 1.0e-1*.

* ``DustMaxFrac`` Fraction of TargetGasMass, indicating the maximum allowable
  dust mass before refinement occurs.
  *Default value: 1.0*.

* ``DustMinFrac`` Fraction of TargetGasMass, indicating the minimum allowable
  dust mass before derefinement occurs.
  *Default value: 1.0e-2*.

* ``NumDustPerSpawn`` Number of dust particles created when a star particle
  stochastically decides to form dust.  If set to one (recommended), a single
  dust particle contains the full grain size distribution.  If set greater than
  one, multiple dust particles are created with contiguous, equal-mass chunks of
  the grain size distribution.
  *Default value: 1*.

* ``DustSubcycleFac`` Force dust particle time-steps to resolve this factor
  times the grain size evolution time-step.  For example, setting this value to
  2.0 allows dust particle time-steps to be twice as large as what grain size
  evolution time-steps would require.
  *Default value: varies*.


Additional Config.sh Options
----------------------------

DUST_LIVE
  master switch, set to an integer for the particle type of dust in output

DL_STOPPING_TIME_CORRECTION
  calculate drag stopping time-scales using correction for supersonic motion

DL_DRAG_SEMI_IMPLICIT
  highly recommended: use second-order semi-implicit drag integrator instead of
  explicit integrator

DL_GRAIN_BINS=10
  set to an integer N and track grain size distribution broken up into N
  log-spaced size bins

DL_GRAIN_BINS_PIECEWISE_LINEAR
  highly recommended: discretize grain size distribution using piecewise linear
  scheme instead of just piecewise constant scheme

DL_GROWTH
  allow dust grains to accrete gas-phase metals

DL_SPUTTERING
  allow dust grains to shrink and lose mass through thermal sputtering

DL_SNE_DESTRUCTION
  allow dust grains to shrink and lose mass through supernova destruction

DL_SHATTERING
  allow dust grains to collide together and shatter to smaller sizes

DL_COAGULATION
  allow dust grains to collide together and coagulate to larger sizes

DL_SHATTERING_DETAILED_INTEGRALS
  calculate shattering and coagulation rates using full piecewise linear grain
  size discretization instead of piecewise constant approximation; this option
  requires more floating point operations

DL_PRODUCTION
  stochastically spawn dust particles from star particles based on dust yield
  calculations; requires GFM (for star particles) and DL_GRAIN_BINS (since dust
  particles are initialized with grain size information)

DL_REFINEMENT
  split dust particles in two whenever masses exceed a threshold value

DL_DEREFINEMENT
  merge dust particles into a neighbor whenever masses drop below a threshold
  value, with option set to an integer (0=merge into closest dust neighbor with
  mass above derefinement mass limit, 1=merge into closest dust neighbor
  without regard to neighbor mass)

DL_SUBCYCLE
  set to an integer N and perform grain size evolution over N subcycles within
  particle timesteps

DL_WINDS
  give dust particles stochastic kicks to mimic winds; requires GFM_WINDS,
  which will set the method for calculating wind direction

DL_ONLY_HIGHRES_DUST
  in zoom-in runs, only spawn dust particles from high-resolution star
  particles to avoid the formation of very large dust particles in
  low-resolution regions; requires DL_REFINEMENT


Authors
-------

* Ryan McKinnon (ryanmck@mit.edu; permanent: ryanmmckinnon@gmail.com)
* Mark Vogelsberger
* Paul Torrey
* Federico Marinacci
* Rahul Kannan


Usage Policy and Citation
-------------------------

Users are free to adopt DUST_LIVE in their own work.  We ask that any papers
based on this work cite

McKinnon R., Vogelsberger M., Torrey P., Marinacci F., Kannan R., 2018, MNRAS
