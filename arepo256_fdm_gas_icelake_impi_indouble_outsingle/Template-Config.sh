#!/bin/bash            # this line only there to enable syntax highlighting in this file

##################################################
#  Enable/Disable compile-time options as needed #
##################################################

#--------------------------------------- Basic operation mode of the code
NTYPES=6                       # number of particle types
#RANDOM_REALIZATION
#TWODIMS
#AXISYMMETRY                   # This is for axisymmetry in cylindrical coordinates (requires TWODIMS and a stationary mesh)
#ONEDIMS
#ONEDIMS_PARALLEL
#ONEDIMS_SPHERICAL
#LONG_X=10.0
#LONG_Y=2.0
#LONG_Z=10.0
#REFLECTIVE_X=1 #=2            # if set to 2, the boundary is inflow/outflow
#REFLECTIVE_Y=1 #=2
#REFLECTIVE_Z=1 #=2

#COOLING
#UVB_SELF_SHIELDING            # gas is self-shielded from the cosmic background based on its density
#USE_SFR
#QUICK_LYALPHA                 # turns dense and cold gas to stars immediately
#QUICK_LYALPHA_LATETIMEONLY    # cooling and star formation only after a certain time
#SFR_KEEP_CELLS
#GAMMA=1.4
#ISOTHERM_EQS
#USE_ENTROPY_FOR_COLD_FLOWS
#ENTROPY_MACH_THRESHOLD=1.1
#PREHEATING
#NOHYDRO
#NOHYDRO_NOTIMESTEP


#--------------------------------------- MPI/Threading Hybrid
#NUM_THREADS=4                           # use OpenMP, with the given number of threads per MPI task
IMPOSE_PINNING
#IMPOSE_PINNING_OVERRIDE_MODE
#GENERIC_ASYNC                           # enables asynchronous communication scheme


#--------------------------------------- Mesh Type
#AMR
VORONOI


#--------------------------------------- SR/GR
#SPECIAL_RELATIVITY
#SPECIAL_RELATIVITY_HLLC
#SR_HLLC_ZERO_COMPVEL
#GENERAL_RELATIVITY
#METRIC_TYPE=1
#ATMOSPHERE_GENERAL_RELATIVITY=1
#ADIABATIC_GENERAL_RELATIVITY=0


#--------------------------------------- MHD
#MHD
#MHD_CT
#MHD_CT_IC
#MHD_CT_CTR_B_ALT
#MHD_POWELL
#MHD_POWELL_LIMIT_TIMESTEP
#MHD_POWELL_SPLIT
#MHD_POWELL_ENERGYLIMITER
#MHD_DEDNER
#MHD_DEDNER_VARIABLE_SPEED
#MHD_DEDNER_WITHOUT_POWELL_OVERRIDE
#MHD_SEEDFIELD
#MHD_SEEDPSPEC
#MHD_THERMAL_ENERGY_SWITCH


#--------------------------------------- NON-IDEAL MHD
#NON_IDEAL_MHD
#OHMIC_DIFFUSION
#AMBIPOLAR_DIFFUSION
#IMPLICIT_OHMIC_DIFFUSION
#OHM_CRANK_NICHOLSON
#ONLY_OHMIC_DIFFUSION
#ONLY_AMBIPOLAR_DIFFUSION
#OHMIC_HEATING
#NON_IDEAL_MHD_EXPLICIT_LIMIT_TIMESTEP


#--------------------------------------- COSMIC RAYS
#DIFFUSION
#COSMIC_RAYS
#COSMIC_RAYS_EXTRA_DIAGNOSTICS
#COSMIC_RAYS_COOLING
#COSMIC_RAYS_ALFVEN_COOLING
#COSMIC_RAYS_STREAMING
#COSMIC_RAYS_STREAMING_EXPLICIT
#COSMIC_RAYS_STREAMING_FACE_OUTPUT
#COSMIC_RAYS_DIFFUSION
#COSMIC_RAYS_DIFFUSION_CONSTANT_TIMESTEP
#COSMIC_RAYS_DIFFUSION_GLOBAL_TIMESTEP
#COSMIC_RAYS_DIFFUSION_EXPLICIT
#COSMIC_RAYS_DIFFUSION_EXPLICIT_LIMITER
#COSMIC_RAYS_DIFFUSION_FULL_NORMAL_GRADIENT
#COSMIC_RAYS_DIFFUSION_ALWAYS_USE_PRECONDITIONER
#COSMIC_RAYS_DIFFUSION_ANISOTROPIC
#COSMIC_RAYS_DIFFUSION_BOUNDARY_X
#COSMIC_RAYS_DIFFUSION_BOUNDARY_Y
#COSMIC_RAYS_DIFFUSION_BOUNDARY_Z
#COSMIC_RAYS_DIFFUSION_LIMITER
#COSMIC_RAYS_DIFFUSION_OLD
#COSMIC_RAYS_DIFFUSION_RELATIVE_CR_ENERGY_THRESHOLD
#COSMIC_RAYS_SN_INJECTION
#COSMIC_RAYS_SHOCK_ACCELERATION
#COSMIC_RAYS_IN_ICS
#COSMIC_RAYS_MAGNETIC_OBLIQUITY
#OUTPUT_CR_PRESSURE_GRADIENT


#--------------------------------------- Riemann solver
#VARIABLE_GAMMA
#RIEMANN_HLL
#RIEMANN_HLLC
#RIEMANN_ROSUNOV
#RIEMANN_HLLD
#RIEMANN_GAMMA
#TRACER

#AMR_CONNECTIONS
#AMR_GRADIENTS
#AMR_REDUCE_DOMAIN_DECOMPOISTION


#--------------------------------------- Reconstruction
#TVD_SLOPE_LIMITER
#TVD_SLOPE_LIMITER_VANLEER
#TVD_SLOPE_LIMITER_SUPERBEE
#TVD_SLOPE_LIMITER_ALBADA
#TVD_SLOPE_LIMITER_MINBEE
#TVD_SLOPE_LIMITER_MINMOD
#TVD_SLOPE_LIMITER_MC
#GRADIENT_LIMITER_DUFFELL
#DISABLE_TIME_EXTRAPOLATION              # use only when you know exactly what you are doing; activating this option will make your results wrong but can tell you about the behavior of your code
#DISABLE_SPATIAL_EXTRAPOLATION           # use only when you know exactly what you are doing; activating this option will make your results wrong but can tell you about the behavior of your code
#FINITE_VOLUME_EXTRAPOLATION_IN_LABFRAME # do extrapolation in labframe instead of moving frame of the gas
#DISABLE_SPATIAL_RECONSTRUCTION
#NO_SCALAR_GRADIENTS                     # disables time and spatial extrapolation for passive scalar fields
#GRADIENTS_GREEN_GAUSS                   # original (now deprecated) gradient estimate, reduced hydro scheme to first order


#--------------------------------------- Mesh motion and regularization
#VORONOI_STATIC_MESH
#VORONOI_STATIC_MESH_DO_DOMAIN_DECOMPOSITION  # for VORONOI_STATIC_MESH force domain decomposition if there exist non-gas particles
REGULARIZE_MESH_CM_DRIFT
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED
REGULARIZE_MESH_FACE_ANGLE
#REGULARIZE_MESH_LLOYD
#REGULARIZE_MESH_SMOOTH
#OUTPUT_MESH_FACE_ANGLE


#--------------------------------------- Time integration options
#FORCE_EQUAL_TIMESTEPS    # this chooses a variable but global time step
TREE_BASED_TIMESTEPS      # non-local time step criterion (take “signal speed” into account)
#DECOUPLE_TIMESTEPS       # allows different timebins for gravity and hydro. use only WITHOUT FORCE_EQUAL_TIMESTEPS
#MUSCL_HANCOCK            # original (now deprecated) time integration scheme, only first order
#RUNGE_KUTTA_FULL_UPDATE
#PM_TIMESTEP_BASED_ON_TYPES=2+4      # select particle types that should be considered in setting the PM time step
#NO_PMFORCE_IN_SHORT_RANGE_TIMESTEP  # if this is on, PM force is not included in short-range time step criterion
#ENLARGE_DYNAMIC_RANGE_IN_TIME  # This extends the dynamic range of the integer timeline from 32 to 64 bit
#NOSTOP_WHEN_BELOW_MINTIMESTEP
#TIMESTEP_OUTPUT_LIMIT          # Limit time steps to write snaps on time for output lists with huge range
#LEGACY_DISPLACEMENT_CONSTRAINT


#--------------------------------------- Direct Volume Raycast Rendering
#DVR_RENDER=0                             # direct volumetric raycasting, 0->stand-alone, 1->on-the-fly
#DVR_RENDER_SMOOTH                        # smooth field
#DVR_RENDER_ORTHOGONAL                    # orthogonal projection
#DVR_NUM_FIELDS=3                         # can be used to set output to a subset of fields, defaults to 13 (all fields), not allowed to be <3
#DVR_STAY_IN_BOX                          # do not ray-trace beyond simulation volume
#DVR_TESS_ALL
#DVR_TOPHAT


#--------------------------------------- Image generation
#VORONOI_MESHOUTPUT                      # 2D and 3D mesh output
#VORONOI_IMAGES_FOREACHSNAPSHOT
#VORONOI_FREQUENT_IMAGES                 # creates images with frequency 'TimeBetweenImages' given in parameterfile, independent of snapshots
#VORONOI_FIELD_DUMP_PIXELS_X=1536
#VORONOI_FIELD_DUMP_PIXELS_Y=150
#VORONOI_VELOCITY_FIELD_2D
#VORONOI_FIELD_COMPENSATE_VX=4.0
#VORONOI_FIELD_COMPENSATE_VY=0
#VORONOI_NEW_IMAGE
#VORONOI_PROJ_TEMP                       # project T instead of u
#VORONOI_PROJ                            # do projection along any predefined direction
#VORONOI_MULTIPLE_PROJECTIONS            # do face-on and edge-on projections by swapping y and z axes
#VORONOI_NOGRADS                         # add an additional set of images where density gradients are not taken into account
#IMAGE_FOOTERS
#VORONOI_PROJ_TAU                        # Projection using opacities, can be used to compute photospheres for stellar calculations -> also use VORONOI_PROJ and OPACITIES
#VORONOI_PROJ_SUBSTEPS=10
#PROJ_WRITE_RAYS
#OPACITIES                               # include opacities for stellar material, tables have to be supplied
#AURIGA_MOVIE
#HCOUTPUT=1+16+32                        # 2^type for the types used for HCOUTPUT


#--------------------------------------- Refinement and derefinement
REFINEMENT_SPLIT_CELLS
REFINEMENT_MERGE_CELLS
#REFINEMENT_LIMIT_STARFORMING_GAS
#GRADIENTREFINEMENT
#REFINEMENT_SPLIT_MOST_DISTANCE_NEIGHBOUR
#REFINEMENT_MERGE_PAIRS
#REFINEMENT_VOLUME_LIMIT
#REFINEMENT_KEEP_INITIAL_VOLUME
#REFINEMENT_HIGH_RES_GAS
#REFINEMENT_CGM
#REFINEMENT_CGM_USE_R200M
#REFINEMENT_AROUND_BH=0                    # spatial refinement scheme near BHs (0: default, 1: ignore cell shape constraints and always refine)
#REFINEMENT_VOLUME_LIMIT_MASS_LIMIT
#DEREFINE_ONLY_DENSE_GAS
#NODEREFINE_BACKGROUND_GRID
#DEREFINE_GENTLY
#OPTIMIZE_MESH_MEMORY_FOR_REFINEMENT       # deletes the mesh structures not needed for refinement/derefinemet to lower the peak memory consumption
#REFINEMENT_AROUND_DM                      # refine around DM particles according to their softening length (useful for binary systems)
#GMC_REFINEMENT
#JEANS_DEREFINEMENT_DENSITY_THRESHOLD
#NO_TARGET_MASS_CONDITION
#DISC_REFINE_ONLY
#REFINE_ONLY_WITH_TRACER
#ROTATING_HIGHRES_REGION
#TRACK_ROTATING_HIGHRES_REGION
#RAMP_REFINE
#SNE_RAMP_REFINE
#MHD_REFINE_ON_DIVB_FACTOR
#REFINE_ABOVE_WNM_DENSITY
#BH_BASED_CGM_ZOOM
#REFINE_MCTR


#--------------------------------------- Mesh-relaxing or mesh-adding (this will not carry out a simulation)
#MESHRELAX                     # this keeps the mass constant and only regularizes the mesh
#MESHRELAX_DENSITY_IN_INPUT
#ADDBACKGROUNDGRID=16
#AMR_REMAP


#--------------------------------------- Gravity treatment
SELFGRAVITY                    # switch on for self-gravity
#HIERARCHICAL_GRAVITY          # use hierarchical splitting of the time integration of the gravity
#CELL_CENTER_GRAVITY           # uses geometric centers to calculate gravity of cells, only possible with HIERARCHICAL_GRAVITY
#NO_GAS_SELFGRAVITY            # switch off gas self-gravity in tree
#GRAVITY_NOT_PERIODIC          # if gravity is not to be treated periodically
#GRAVITY_TALLBOX               # special switch for making treating gravity in z-extended box, with x/y periodic, and z nonperiodic. LONG_Z may be used but must be an integer.
#ALLOW_DIRECT_SUMMATION
#DIRECT_SUMMATION_THRESHOLD=1000
#EXACT_GRAVITY_FOR_PARTICLE_TYPE=4 # N-squared fashion gravity for a small number of particles of the given type
#NO_SELFGRAVITY_TYPE=1         # exclude particle type from self-gravity (can be used with exact gravity)
#NO_GRAVITY_TYPE=1             # disable computation of gravity on particle type
#EXACT_GRAVITY_REACTION        # include reaction to other particle types when using exact gravity
#EXTERNALGRAVITY               # switch on for external potential
#EXTERNALGY=0.0
#EXTERNALSHEARBOX
#STATIC_ISOTHERMAL_CLUSTER
#EXTERNALSHEARBOX_KSRATE_RANDOM
#EXTERNALSHEARBOX_KSRATE_UPDATE_PARAM
#ENFORCE_JEANS_STABILITY_OF_CELLS_EEOS
ENFORCE_JEANS_STABILITY_OF_CELLS    # this imposes an adaptive floor for the temperature
#EVALPOTENTIAL                  # computes gravitational potential
#EXTERNALSHEETY
#COMPUTE_POTENTIAL_ENERGY
#ACCRETE_ONTO_CENTRAL_POTENTIAL # Allow mass to be accreted onto the central potential (needs CENTRAL_MASS_POTENTIAL)


#--------------------------------------- Gravity softening
#NSOFTTYPES=4                  # Number of different softening values to which particle types can be mapped.
#MULTIPLE_NODE_SOFTENING       # If a tree node is to be used which is softened, this is done with the softenings of its different mass components
#INDIVIDUAL_GRAVITY_SOFTENING=2+4  # bitmask with particle types where the softenig type should be chosen with that of parttype 1 as a reference type
#ADAPTIVE_HYDRO_SOFTENING
#NSOFTTYPES_HYDRO=64           # this is only relevant for ADAPTIVE_HYDRO_SOFTENING can can be set to override default value of 64


#--------------------------------------- TreePM Options
PMGRID=512
#ASMTH=1.25
#RCUT=6.0

#PLACEHIGHRESREGION=2
#ENLARGEREGION=1.1
#GRIDBOOST=2
#ONLY_PM                # only use long-range (PM) force for gravity, without tree force

#FFT_COLUMN_BASED
#PM_ZOOM_OPTIMIZED

#BINS_PS=2000           # number of k-bins to use for power spectrum calculations
#POWERSPEC_FOLDFAC=16.  # folding factor to obtain an estimate of the power spectrum on
                        # very small scales

#NUMPART_PER_TASK_LARGE # Should be activated if eight times the particle load per
                        # processor exceeds 2^31


#--------------------------------------- Things that are always recommended
#AUTO_SWAP_ENDIAN_READIC # Enables automatic ENDIAN swapping for reading ICs
CHUNKING                 # will calculate the gravity force in interleaved blocks.
                         # This can reduce imbalances in case multiple iterations due to
                         # insufficient buffer size need to be done.


#--------------------------------------- Single/Double Precision
DOUBLEPRECISION=1
DOUBLEPRECISION_FFTW
#OUTPUT_IN_DOUBLEPRECISION                # snapshot files will be written in double precision
#INPUT_IN_DOUBLEPRECISION                 # initial conditions are in double precision
#OUTPUT_COORDINATES_IN_DOUBLEPRECISION    # will always output coordinates in double precision
#NGB_TREE_DOUBLEPRECISION                 # if this is enabled, double precision is used for the neighbor node extension


#--------------------------------------- On the fly FoF group finder
#FOF                                # enable FoF output
#FOF_PRIMARY_LINK_TYPES=2           # 2^type for the primary dark matter type
#FOF_SECONDARY_LINK_TYPES=1+16+32   # 2^type for the types linked to nearest primaries
#FOF_SECONDARY_LINK_TARGET_TYPES=   # should normally be set to a list of all dark matter types (in zoom runs), if not set defaults to FOF_PRIMARY_LINK_TYPES
#FOF_GROUP_MIN_LEN=32               # enforce minimum number of particles within a group
#FOF_LINKLENGTH=0.16                # Linking length for FoF (default=0.2)
#FOF_FUZZ_SORT_BY_NEAREST_GROUP=0   # sort fuzz particles by nearest group and generate offset table in catalog (=1 writes nearest group number to snapshot)
#FOF_STOREIDS                       # store IDs in group/subfind catalogue, do not order particles in snapshot files by group order
#USE_AREPO_FOF_WITH_GADGET_FIX      # Needed in order to run FOF with Arepo on Gadget snapshot files, if gas is present and should be linked to the FOFs
#ADD_GROUP_PROPERTIES               # This can be used to calculate additional properties for an already existing group catalogue. These are then added as additional columns to the HDF5 group catalogues.
#ADD_MAGNETIC_GROUP_PROPERTIES

#--------------------------------------- Subfind
#SUBFIND                            # enables substructure finder
#SAVE_HSML_IN_SNAPSHOT              # stores hsml, density, and velocity dispersion values in the snapshot files

#SUBFIND_MEASURE_H2MASS             # special measuremenat option for mass in molecular hydrogen
#SUBFIND_CALC_MORE                  # calculates also the velocity dispersion in the local density estimate (this is automatically enabled by several other options, e.g. SAVE_HSML_IN_SNAPSHOT)
#SUBFIND_EXTENDED_PROPERTIES        # adds calculation of further quantities related to angular momentum in different components


#--------------------------------------- SFR/feedback model
#METALS
#MIN_METALLICITY_ON_STARTUP
#STELLARAGE

#SOFTEREQS
#MODIFIED_EOS
#SLOW_RELAX_TO_EOS
#STEEPER_SFR_FOR_STARBURST
#SF_STELLAR_MASS_TO_GAS_MASS_RATIO

#-------------------------------------- AGN stuff
#BLACK_HOLES               # enables Black-Holes (master switch)
#BH_THERMALFEEDBACK        # quasar-mode: couple a fraction of the BH luminosity into surrounding
#BH_THERMALFEEDBACK_ACC    # quasar-mode: bursty quasar-mode, accumulate thermal energy
#BH_NF_RADIO               # radio-mode model based on Nulsen & Fabian theory
#DRAINGAS=1                # non-stochastic smooth accretion (1: on, 2: on + cell rho, 3: on + gas drained from all cells within hsml)
#BH_EXACT_INTEGRATION      # integrates analytically mass accretion
#BH_BONDI_DEFAULT          # default Bondi prescription
#BH_BONDI_DENSITY          # Bondi -> density dependent
#BH_BONDI_CAPTURE          # Bondi -> capture mass dependent
#BH_BONDI_DISK_VORTICITY   # Bondi -> vorticity dependent
#BH_DO_NOT_PREVENT_MERGERS # When this is enabled, BHs can merge irrespective of their relative velocity
#BH_USE_GASVEL_IN_BONDI    # only when this is enabled, the surrounding gas velocity is used in addition to the sounds speed in the Bondi rate
#BH_USE_ALFVEN_SPEED_IN_BONDI  # when this is enabled the alfven speed is added to the gas sound speed in the Bondi rate and the total gas pressure around the BH includes the magnetic contribution when compared the the reference pressure in BH_PRESSURE_CRITERION (requires MHD)
#MASSIVE_SEEDS             # BH seeds assigned large dynamical mass, such that ideally no repositioning is needed anymore
#MASSIVE_SEEDS_MERGER
#BH_NEW_CENTERING          # an alternative to the BH_FRICTION and REPOSITION_ON_POTMIN switches
#BH_DRAG                   # Drag on black-holes due to accretion: current implementation simply double-accounts for the accretion momentum transfer, and should therefore not be used
#REPOSITION_ON_POTMIN      # repositions hole on potential minimum (requires EVALPOTENTIAL)
#BH_REPOSITION_POTMIN_TRUST_THRESHOLD
#BH_PRESSURE_CRITERION
#BH_RELATIVE_NGB_DEVIATION # Maximum NGB number deviation calculated relative to total number of neighbours
#OUTPUT_BLACK_HOLE_TIMESTEP # outputs the 3 time steps for BH particles
#BH_FRICTION                # Estimates the local DM density around BH and applies a friction force to the relative velocity, meant as a replacement for REPOSITION_ON_POTMIN
#BH_FRICTION_AGGRESSIVE
#BH_HARMONIC_OSCILLATOR_FORCE
BH_INFLOW_RATE

#-------------------------------------- AGN spin evolution and recoil merger kicks
#BH_RECOIL_KICK             # Includes the remnant recoil of a BH merger
#BH_SPIN_EVOLUTION          # When this is enabled, spin evolution of black holes is computed
#BH_SPIN_MODEL=0            # Spin model to be used: 0-Prolonged spin model. 1-Chaotic spin model. 2-Mass dependend model. 3-Self-gravity dependend model.

#-------------------------------------- other AGN stuff
#UNIFIED_FEEDBACK        # activates BH_THERMALFEEDBACK at high Mdot and BH_BUBBLES FEEDBACK al low Mdot (-->OBSOLETE: replaced by BH_NEW_RADIO)
#BH_BUBBLES              # calculate bubble energy directly from the black hole accretion rate (-->OBSOLETE: replaced by BH_NEW_RADIO)
#BH_MAGNETIC_BUBBLES     # inject part of the  bubble energy as magnetic energy
#BH_MAGNETIC_DIPOLAR_BUBBLES #inject part of the bubble energy as magnetic energy, field arranged as a dipole with random orientation
#BH_ADIOS_WIND
#BH_ADIOS_DENS_DEP_EFFICIANCY  # makes the radiative efficiency density dependend
#BH_ADIOS_WIND_WITH_QUASARTHRESHOLD  # use a threshold value ("qusarthrehold") of bondi-rate over Eddington rate to decide about quasar mode vs. adios wind
#BH_ADIOS_WIND_WITH_VARIABLE_QUASARTHRESHOLD  # scales the threshold with black hole mass (with a factor (M_BH/M_ref)^2, where M_ref = 10^8 Msun)
#BH_ADIOS_WIND_DIRECTIONAL  # puts in momentum preferentially along a random direction
#BH_ADIOS_RANDOMIZED        # inputs momentum along alternating random directions
#BH_ADIOS_ONLY_ABOVE_MINIMUM_DENSITY   # disable ADIOS wind if density around blackhole drops below a certain fraction of the star formation density
#BH_CONTINOUS_MODE_SWITCH # calculates fraction of thermal and mechanical feedback energy depending on eddington factor and mass (continously in both quantities)

#-------------------------------------- Black Hole Refinement and Bipolar Options
#REFINEMENT_AROUND_BH_FIXED
#SUPPRESS_SF_IN_REFINEMENT_REGION
#BH_BIPOLAR_FEEDBACK


#---------------------------------------- Passive Tracers
#TRACER_FIELD                        # passive scalar field which is advected in proportion to fluid mass fluxes

#TRACER_MC=3                         # Monte Carlo tracer particles: master switch (value specifies output parttype)
#GENERATE_TRACER_MC_IN_ICS           # add a fixed number (given in the parameter file) of MC tracers to each gas cell in ICs
#TRACER_MC_NUM_FLUID_QUANTITIES=13   # number of fluid quantities to be stored for MC tracers - must match the number in TRACER_MC_STORE_WHAT
#TRACER_MC_STORE_WHAT=1+2+4          # bit mask for quantities to store (see allvars.h for bitmask)
#TRACER_NO_RESET_EACH_SNAP           # do not set tracked fluid quantities to zero after writing each snapshot
#TRACER_MC_CHECKS                    # carries out frequent consistency checks

#TRACER_PARTICLE=2                   # Velocity Field tracer particles: master switch (value specified parttype)
#GENERATE_TRACER_PARTICLE_IN_ICS     # add tracer particles at positions of cell vertices in ICs
#TRACER_PART_NUM_FLUID_QUANTITIES=8  # number of fluid quantities to be stored for velocity tracers - must match the value given to TRACER_PART_STORE_WHAT
#TRACER_PART_STORE_WHAT=1+2+4        # bit mask for quantities to store (see allvars.h for bitmask)

#TRACER_TRAJECTORY
#TRACER_TRAJECTORY_GENERATE
#TRACER_TRAJECTORY_EXTENDED_OUTPUT

#OUTPUT_MCTRNUM                      # write number of MC tracers in each cell to the output snapshots


#-------------------------------------------- Things for special behavior
#READ_DM_AS_GAS
#NO_ID_UNIQUE_CHECK
#RUNNING_SAFETY_FILE            # if file './running' exists, do not start the run
#LOAD_TYPES=1+2+4+16+32
#READ_COORDINATES_IN_DOUBLE
#IDS_OFFSET=1                   # offset for gas particles if created from DM
#TILE_ICS
#COMBINETYPES                   # reads in the IC file types 4+5 as type 3 (useful for doing gas runs of Aquarius ICs)
#MULTIPLE_RESTARTS
#TOLERATE_WRITE_ERROR
#OPTIMIZE_MEMORY_USAGE          # optimize for memory, not for speed. Note: this is dangerous for high dynamic range simulations with mixed precision, since some position variables are singles instead of doubles
#SUBBOX_SNAPSHOTS
#PROCESS_TIMES_OF_OUTPUTLIST
#EXTENDED_GHOST_SEARCH          # This extends the ghost search to the full 3x3 domain instead of the principal domain
#DOUBLE_STENCIL                 # this will ensure that the boundary region of the local mesh is deep enough to have a valid double stencil for all local cells
#TETRA_INDEX_IN_FACE            # adds an index to each entry of VF[] and DC[] to one of the tetrahedra that share this edge
VORONOI_DYNAMIC_UPDATE          # keeps track of mesh connectivity, which speeds up mesh construction
VORONOI_MESH_KEEP_DT_AND_DTC    # keeps DTC and DT in memory, i.e. for anisotropic transport solvers
#FREE_DC_BEFORE_FULL_MESH_CONSTRUCTION # frees DC in memory when full mesh gets constructed. Saves memory but slows down mesh construction on global time steps.
#COFFEE_PROBLEM
#NOH_PROBLEM
#SHIFT_BY_HALF_BOX
#DISABLE_VELOCITY_CSND_SLOPE_LIMITING
NO_MPI_IN_PLACE
NO_ISEND_IRECV_IN_DOMAIN
FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
#USE_MPIALLTOALLV_IN_DOMAINDECOMP
#MPI_HYPERCUBE_ALLGATHERV       # some MPI libraries may use quite a bit of internal storage for MPI_Allgatherv. This uses hypercubes instead as a work-around
#MPISENDRECV_CHECKSUM
#MPISENDRECV_SIZELIMIT
#MPI_MESSAGE_SIZELIMIT_IN_MB=200 # maximum amount of data to transfer in a single MPI call
#MYIBARRIER                     # include src/mpi_utils/myIBarrier.c in compilation (unused?)
#DO_NOT_CREATE_STAR_PARTICLES
#ALLOWEXTRAPARAMS               # do not terminate when encountering unexpected names in
                                # the parameter file
#NOCALLSOFSYSTEM                # do not use the system() function, i.e. disable all
                                # calls to external programs
#FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES  # this can be used to load SPH ICs that contain identical particle coordinates
#VEL_POWERSPEC                  # compiles in a code module that allows, via restart flag 7, the calculation of a gas velocity power spectrum of a snapshot
#VEL_POWERSPEC_BOX
#ADJ_BOX_POWERSPEC              # compiles in a code module that allows, via restart flag 7, the calculation of gas power spectra of a snapshot with an adjustable box (user defined center and size)
#DISABLE_OPTIMIZE_DOMAIN_MAPPING
#DO_NOT_RANDOMIZE_DOMAINCENTER
#RECOMPUTE_POTENTIAL_IN_SNAPSHOT  # needed for postprocess option 18 that can be used to calculate potential values for a snapshot
#COMPUTE_VORONOI_DM_DENSITY_IN_POSTPROC
#ACTIVATE_MINIMUM_OPENING_ANGLE   # this does not open tree nodes under the relative opening criterion any more if their opening angle has dropped below a minimum angle
#USE_DIRECT_IO_FOR_RESTARTS     # Try to use O_DIRECT for low-level read/write operations of restart files to circumvent the linux kernel page caching
#EQUIPARTITION
#PERTURB_VELOCITIES             # continuously perturb velocities when running simulation
#UVB_OFF                        # No UVB
#UVB_START                      # UVB switched on after a redshift supplied in parameterfile

#CUDA                       # enables CUDA support in Arepo
#CUDA_INSTRUMENT            # This enables instrumentation support for the nvidia profiler
#USE_DSDE                   # try to use a dynamic sparse data exchange paradigm to get rid off sparse MPI_Alltoall patterns on large partitions

#DISABLE_MEMORY_MANAGER     # disable AREPO's custom memory manager (use standard
                            # malloc() instead)
#HUGEPAGES                  # use huge pages for memory allocation, through hugetlbfs library
#DETAILEDTIMINGS            # creates individual timings entries for primary/secondary kernels to diagnose work-load balancing

#BITS_PER_DIMENSION=42      # Peano-Hilbert order
#OVERRIDE_PEANOGRID_WARNING

#WALLCLOCK
#NEW_FFT


#--------------------------------------- Compatibility options
#NOTYPEPREFIX_FFTW
#DECALPHA_NOPOW


#--------------------------------------- Output/Input options
#READ_IN_ALL_IC_FIELDS      # Read in all fields that are in the initial condition files
#UPDATE_GRADIENTS_FOR_OUTPUT
#REDUCE_FLUSH
#OUTPUT_REFBHCOUNTER
#OUTPUT_EVERY_STEP
#GODUNOV_STATS
#OUTPUT_CPU_CSV
#OUTPUT_TASK
#OUTPUT_TIMEBIN_HYDRO
#OUTPUT_PRESSURE_GRADIENT
#OUTPUT_DENSITY_GRADIENT
#OUTPUT_VELOCITY_GRADIENT
#OUTPUT_BFIELD_GRADIENT
#OUTPUT_VERTEX_VELOCITY
#OUTPUT_VERTEX_VELOCITY_DIVERGENCE
#OUTPUT_VOLUME
#OUTPUT_CENTER_OF_MASS
#OUTPUT_SURFACE_AREA
#OUTPUT_PRESSURE
#OUTPUTPOTENTIAL
#OUTPUTACCELERATION
#OUTPUTTIMESTEP
#OUTPUT_SOFTENINGS            # output particle softenings
#OUTPUTGRAVINTERACTIONS       # output gravitatational interactions (from the tree) of particles
HAVE_HDF5                     # needed when HDF5 I/O support is desired
#HDF5_FILTERS                  # activate snapshot compression and checksum for HDF5 output
#OUTPUT_XDMF                   #writes an .xmf file for each snapshot, which can be read by visit (with the hdf5 snapshot)
#OUTPUTCOOLRATE                # outputs cooling rate, and conduction rate if enabled
#OUTPUT_HE_IONIZATION_STATE    # outputs the helium ionization state in detail (HeI, HeII, HeIII)
#OUTPUT_DIVVEL                 # output  velocity divergence
#OUTPUT_CURLVEL                 # output  velocity curl
#OUTPUT_COOLHEAT               # output actual energy loss/gain in cooling/heating routine
#OUTPUT_VORTICITY
#OUTPUT_CELL_SPIN
#MEASURE_DISSIPATION_RATE      # measures and outputs dissipation rate. Note: requires USE_ENTROPY_FOR_COLD_FLOWS, even though it will then always use the thermal energy update
#OUTPUT_MACHNUM                # output maximum mach number of a cell
#OUTPUT_ENTROPY
#OUTPUT_CSND
#OUTPUT_AMR_REFFLAG
#HIGH_FREQUENCY_OUTPUT_STARS

#POWERSPECTRUM_ON_THE_FLY


#--------------------------------------- Testing and Debugging options
DEBUG                          # enables core dumps
#NDEBUG                        # disables assertions (specified by the C standard)
#DEBUG_ENABLE_FPU_EXCEPTIONS   # tries to enable FPU exceptions
#DEBUG_REFINE
#RESTART_DEBUG
#VERBOSE                       # reports readjustments of buffer sizes
HOST_MEMORY_REPORTING          # reports the available system memory after start-up by analyzing /proc/meminfo
MEMORY_MANAGER_USE_MPROTECT    # mark memory for unallocated memory blocks as
                               # inaccessible using mprotect() in order to detect memory
                               # access bugs (e.g. buffer overflows)
#VTUNE_INSTRUMENT
#FORCETEST=0.001               # calculates for given fraction of particles direct summation forces to check accuracy of tree force
#FORCETEST_TESTFORCELAW=1      # this enables a special test to measure the effective force law of the code, can be set to 1 or 2

#MAX_VARIATION_TOLERANCE=0.1   # tolerance for performance variation in health tests
#DISABLE_HEALTHTEST            # disables health tests at start-up (healthtest.c)
#DISABLE_HEALTHTEST_FULL_HYPERCUBE  # disables full hypercube health test at
                                    # start-up, which can take a long time for
                                    # large numbers of MPI tasks (healthtest.c)

#PERFORMANCE_TEST_SPARSE_MPI_ALLTOALL
#VORONOI_TEST
#THERMAL_INSTABILITY_TEST
#CHECK_LOCAL_RANK              # additional checks in parallel_sort.c
#MHD_DONT_PRINT_BMAGSUM


#--------------------------------------- Static Disk Potential
#EXTERNALDISKPOTENTIAL
#DISK_MASS_M0=1.0
#DISK_SCALE_R0=1.0


#--------------------------------------- Static NFW Potential
#STATICNFW
#NFW_C=12
#NFW_M200=100.0
#NFW_Eps=0.01
#NFW_DARKFRACTION=0.87
#NFW_h=0.7


#--------------------------------------- Static Isothermal Sphere Potential
#STATICISO
#ISO_M200=100.0
#ISO_R200=160.0
#ISO_Eps=0.1
#ISO_FRACTION=0.9


#--------------------------------------- Static Hernquist Potential
#STATICHQ
#HQ_M200=186.015773
#HQ_C=10.0
#HQ_A=10.0
#HQ_DARKFRACTION=0.9


#--------------------------------------- Growing Disk Potential
#GROWING_DISK_POTENTIAL


#--------------------------------------- Dark energy
#DARKENERGY # Enables Dark Energy
#TIMEDEPDE  # read w(z) from a DE file
#RESCALEVINI # rescale v_ini in read_ic/read_ic_cluster
#EXTERNALHUBBLE # reads the hubble function from the DE file
#TIMEDEPGRAV # resacles H and G according to DE model
#DARKENERGY_DEBUG # enable writing of drift/kick table


#--------------------------------------- Initial conditions options
#LONGIDS
#OFFSET_FOR_NON_CONTIGUOUS_IDS
#GENERATE_GAS_IN_ICS
#SPLIT_PARTICLE_TYPE=4+8
#NTYPES_ICS=6 # number of particle types in ICs, if not NTYPES (only works for 6, and non-HDF5 ICs!)


#--------------------------------------- Simple turbulence test
#VS_TURB
#POWERSPEC_GRID=128

#AB_TURB
#AB_TURB_DECAYING


#--------------------------------------- Degenerate Equation of State
#EOS_NSPECIES=3  # first species is X(H), which is relevant for the EOS
#EOS_DEGENERATE
#EOS_COULOMB_CORRECTIONS
#EOS_COULOMB_CORRECTIONS_SMOOTH
#RELAXOBJECT
#RELAXOBJECT_COOLING
#RELAXOBJECT_COOLING2
#RELAXOBJECT_BINARY
#RELAX_RUNTIME
#INSPIRAL
#PASSIVE_SCALARS=3
#GW_SIGNAL

#--------------------------------------- Passive element tracking
#EOS_PASSIVE

#--------------------------------------- OPAL Equation of State
#EOS_OPAL


#--------------------------------------- Nuclear Network
#NUCLEAR_NETWORK
#NETWORK_NSE
#NETWORK_PARDISO
#NETWORK_SCREENING
#REACLIB1
#NUCLEAR_NETWORK_DETONATE
#NUCLEAR_NETWORK_DETONATE_CORE
#NUCLEAR_NETWORK_DETONATE_POSITION
#NUCLEAR_NETWORK_TIMESTEP_LIMITER
#NUCLEAR_NETWORK_USE_SHOCKFINDER
#NUCLEAR_NETWORK_DISABLE_BURNING_IN_SHOCK
#NUCLEAR_NETWORK_LIMIT_COMPOSITION_CHANGE
#NUCLEAR_NETWORK_ALT_MPI_PARALLELISATION


#--------------------------------------- Radiative transfer options
#RT_ENABLE                   # RT master switch
#RT_COOLING_PHOTOHEATING
#RT_ADVECT                   # enable advection of radiation field
#RT_CGMETHOD                 # enables CG method solution of the RT advection
#RT_SLOWLIGHT                # enable slow light approximation
#RT_N_DIR=2                  # track this number of locally brightest sources (one of them is diffuse field)
#RT_COMBINE_N_DIR_IN_OUTPUT  # writes only a single summed photon/photon-density field into output files
#RT_ALLOW_ABSORBING_CELLS    # if this is set, all cells with ID >= 1000000000 will absorb radiation
#RT_SPREAD_SOURCE
#RT_STELLAR_SOURCES
#RT_HEALPIX_NSIDE=1          # if this is set, a discretization of the solid angle is used instead of brightest source selection
#RT_INCLUDE_HE
#RT_SHORT_CHARACTERISTICS
#RT_OUTPUT_COL_DENS
#SOURCE_PERIODIC

#DO_NOT_MOVE_GAS
#HYDROGEN_ONLY


#--------------------------------------- Calculate Hessian Matrix
#SECOND_DERIVATIVES
#SLOPE_LIMIT_HESSIANS
#RECONSTRUCT_GRADIENTS
#OUTPUT_HESSIAN


#--------------------------------------- Navier-Stokes terms
#GLOBAL_VISCOSITY           # needs dynamic and bulk coefficients
#USE_KINEMATIC_VISCOSITY    # needs only one input parameter
#ALPHA_VISCOSITY=2          # for accretion disks
#LOCAL_VISCOSITY=1          # = 1 Sutherland viscosity/= 2 Spitzer viscosity
#THERMAL_CONDUCTION
#TRACER_DIFFUSION           # requires TRACER_FIELD switched on


#--------------------------------------- Circumstellar Disks
#CIRCUMSTELLAR              # Master switch
#CIRCUMSTELLAR_WBOUNDARIES
#CIRCUMSTELLAR_IRRADIATION
#CIRCUMSTELLAR_REFINEMENTS
#CIRCUMSTELLAR_SINKS
#CIRCUMSTELLAR_SINKS_ALTERNATIVE # alternative calculation of cirumstellar sink, i.e. swallow cell entirely
#CIRCUMSTELLAR_PLANET_GROWTH     # Requires BLACK_HOLES turned on
#GRAVITY_FROM_STARS_PLANETS_ONLY # Requires EXTERNALGRAVITY turned on
#CENTRAL_MASS_POTENTIAL     # Point-mass potential
#BINARY_POTENTIAL           # Fixed star-planet circular orbit
#LOCALLY_ISOTHERM_DISK      # Isothermal Equation of state at each radii.


#--------------------------------------- Special boundaries within domain
#SPECIAL_BOUNDARY          # Main Switch
#COAXIAL_BOUNDARIES        # e.g. Couette flow-type boundaries
#BOUNDARY_FLAG


#--------------------------------------- Windtunnel
#WINDTUNNEL
#WINDTUNNEL_COORD=0                    # sets the coordinate in which the wind blows (0,1,2 for x,y,z)
#WINDTUNNEL_EXTERNAL_SOURCE
#WINDTUNNEL_FIXVARIABLESININJECTIONREGION # enables a region with fixed properties
#WINDTUNNEL_REFINEMENT_VOLUME_LIMIT # Volume refinement option for windtunnel setup. REFINEMENT_VOLUME_LIMIT should also be enabled.
#WINDTUNNEL_READ_IN_BFIELD # Overwrites B-field in injection region with a field specified in a separate file.


#--------------------------------------- Dark Matter Windtunnel
#DM_WINDTUNNEL                         # Master switch
#DM_WINDTUNNEL_EXTERNAL_SOURCE         # Reads injection region parameters from a separate file
#DM_WINDTUNNEL_STARS


#--------------------------------------- Boundaries with optional inflow/outflow
#BOUNDARY_INFLOWOUTFLOW_MINID=10000000   # defines the ID range describing inflow/outflow nozzle of wind-tunnel
#BOUNDARY_INFLOWOUTFLOW_MAXID=20000000
#BOUNDARY_INFLOWOUTFLOW_FLUIDSIDE_MINID=10000000   # defines the ID range describing inflow/outflow nozzle of wind-tunnel
#BOUNDARY_INFLOWOUTFLOW_FLUIDSIDE_MAXID=20000000
#BOUNDARY_REFL_FLUIDSIDE_MINID=30000000  # defines the ID ranges describing a reflective boundary
#BOUNDARY_REFL_FLUIDSIDE_MAXID=30000000  # defines the ID ranges describing a reflective boundary
#BOUNDARY_REFL_SOLIDSIDE_MINID=40000000
#BOUNDARY_REFL_SOLIDSIDE_MAXID=40000000
#BOUNDARY_REFL_ACTS_AS_SOURCE            # makes the boundary act as a source (using the inner values)
#BOUNDARY_STICKY_MINID=50000000          # this-ID range specifies cells that will not be moved, and neighbors of these cells will only do mesh regularization motions
#BOUNDARY_STICKY_MAXID=60000000
#STICKYFLAGS
#OUTPUT_STICKYFLAGS


#--------------------------------------- GFM - Galaxy Formation Module
#GFM                                    #master switch
#GFM_STELLAR_EVOLUTION=0                #stellar evolution: 0->default, 1->no mass loss (beta value changes + MassMetallicity & MassMetals inconsistent internally with cell dynamical mass) 2->call only test routine
#GFM_STELLAR_EVOLUTION_NO_ELEMENTS
#GFM_CONST_IMF=1                        #0 for Chabrier (default), 1 for a pure power-law (requires parameter IMFslope, e.g. -2.35 for Salpeter)
#GFM_VARIABLE_IMF=0                     #0 for a pure power-law that depends on DM-veldisp
#GFM_PREENRICH                          #pre enrich gas at given redshift
#GFM_SET_METALLICITY                    #set the metallicity of gas in solar metallicity units
#GFM_NO_METAL_ENRICHMENT                #disable metal production by SNII and AGB stars
#GFM_EXACT_NUMNGB                       #use direct neighbor count instead of kernel weighted neighbor count
#GFM_WINDS                              #decoupled ISM winds
#GFM_WINDS_VARIABLE=0                   #decoupled ISM winds: 0->scale winds with halo mass, requires FoF, 1->sigma winds
#GFM_WINDS_VARIABLE_HUBBLE              #add an additional H(z)^(-1/3) factor to the wind scaling, such that it scales with halo mass not halo velocity dispersion
#GFM_WINDS_HUBBLESCALING                #scale the wind energy fraction with the Hubble rate, limit the maximum to 1
#GFM_WINDS_MASSSCALING                  #scale the wind energy mass loading with halo mass (equivalent to scaling the wind energy fraction with halo virial radius)
#GFM_WIND_ENERGY_METAL_DEPENDENCE       #this can be used to decrease the wind energy for high metallicity (mimicking higher cooling losses)
#GFM_WIND_ENERGY_METAL_DEPENDENCE_TANH  #this selects an alternative functional form for the transition, requires GFM_WIND_ENERGY_METAL_DEPENDENCE
#GFM_WINDS_STRIPPING                    #wind metal stripping
#GFM_WINDS_THERMAL                      #not only give the wind kinetic energy but also thermal energy
#GFM_WINDS_THERMAL_NEWDEF               #with this switch, the thermal energy is specified as a fraction of the total energy
#GFM_BIPOLAR_WINDS=1                    #decoupled ISM winds: bipolar winds: 0->default, 1->relative to motion of FOF group, 3->parallel to spin of star-forming gas in halo
#GFM_WINDS_LOCAL                        #energy-driven decoupled local sigma winds
#GFM_STELLAR_FEEDBACK                   #local SNIa and AGB energy and momentum feedback
#GFM_PRIMORDIAL_RATES                   #updated coefficients for primordial chemistry and cooling
#GFM_COOLING_METAL                      #metal line cooling
#GFM_UVB_CORRECTIONS                    #reionization energy corrections
#GFM_AGN_RADIATION                      #cooling suppression/heating due to AGN radiation field (proximity effect)
#GFM_STELLAR_PHOTOMETRICS               #calculate stellar magnitudes for different filters based on GALAXEV/BC03
#GFM_OUTPUT_MASK=1+2+4+8+16+32+64+128   #which fields to output (see io_fields.c)
#GFM_CHECKS                             #this checks the consistency of the AuxDataID/PID indices of stars and black holes every time step
#GFM_DISCARD_ENRICHMENT_GRADIENTS       #this disables the gradient extrapolation of the passively advected metallicity scalar variables
#GFM_NORMALIZED_METAL_ADVECTION         #this introduces an additional pseudo element for all untracked metals and normalizes the extrapolated abundance vectors to unity
#GFM_OUTPUT_BIRTH_POS                   #output BirthPos and BirthVel for all star particles
#GFM_CHEMTAGS                           #see documentation/modules_GFM_chemtags
#GFM_WINDS_SAVE_PARTTYPE=2              #save wind particles as separate particle type instead of mixed with 4 (stars)
#GFM_DISCRETE_ENRICHMENT                #allow stars to enrich nearby gas from stellar evolution only above some delta mass fraction threshold
#GFM_NO_NEGATIVE_ELEMENT_MASS_RELEASED  #do not allow that negative yields for each element consume more mass than contained in ejecta elemental composition
#GFM_SPLITFE                            #see documentation/modules_GFM_chemtags
#GFM_SPLITFE_ADDINAGB                   #add in the AGB iron half-half on the two iron SNIa/SNII tags such that the sum of them should be equal to the total iron
#GFM_RPROCESS                           #see documentation/modules_GFM_chemtags, must have GFM_SPLITFE toggled as well
#GFM_LAMBDA                             #output all cooling rates
#GFM_RPROCESS_CHANNELS=10               #alternative to GFM_PROCESS, use many different channels, number is number of independent r-process channels
#GFM_RPROCESS_CHANNELS_NS_KICKS         #include neutron star kicks for NSNS mergers
#GFM_RPROCESS_NSNS=7                    #the number of channels of GFM_RPROCESS_CHANNELS that are NSNS channels, the rest are SN channels
#GFM_SPROCESS                           #adds s-process elements, need to exist in yield tables
#GFM_SNIA_ENERGY_INJECTION              #add thermal energy if Ia's
#GFM_SINGLE_CELL_INJECTION


#--------------------------------------- Dust physics
#GFM_DUST                               #formation and evolution of dust, requires GFM_STELLAR_EVOLUTION
#GFM_DUST_DESTMODE=0                    #dust destruction mode: 0->default (uses supernova rate), 1->constant destruction timescale
#GFM_DUST_SPUTTERING=1                  #sputtering of dust grains by gas-phase metals: 0->first principles calculation, 1->using empirical timescale
#GFM_DUST_COOLING                       #high temperature dust cooling
#GFM_DUST_MRN                           #MRN grain size distribution; otherwise single grain size with size in mu specified in parameter file
#GFM_DUST_CAP                           #cap negative dust masses
#GFM_DUST_ISMGROWTH                     #stop dust from growing outside of star-forming regions
#GFM_DUST_LOWGROWTH                     #use low dust growth as in McKinnon+ 2017 paper


#--------------------------------------- Live dust physics
#DUST_LIVE=3                            #turns on live dust particles, value specifies output parttype; parttype 3 is recommended
#DL_STOPPING_TIME_CORRECTION            #include higher-order corrections to stopping timescale for supersonic flow, makes analytic tests more difficult
#DL_DRAG_SEMI_IMPLICIT                  #make use of drag analytic solution for velocity updates, instead of requiring explicit drag time steps
#DL_DRAG_BACKREACTION                   #allow dust to drag gas using an equal and opposite drag force
#DL_NODRAG                              #do not couple dust to gas through a drag force
#DL_GRAIN_BINS=10                       #track grain size distribution information for dust particles, using the specified number of bins; requires cooling, star formation
#DL_GRAIN_BINS_PIECEWISE_LINEAR         #allow grain size distribution bins to be piecewise linear, not just piecewise constant
#DL_GROWTH                              #enable growth of dust mass by accumulating gas-phase metals
#DL_SPUTTERING                          #loss of dust mass due to thermal sputtering
#DL_SNE_DESTRUCTION                     #loss of dust mass due to supernova shocks
#DL_SHATTERING                          #grain size shattering
#DL_COAGULATION                         #grain size coagulation
#DL_SHATTERING_DETAILED_INTEGRALS       #do not use piecewise constant approximation for integrals used in shattering and coagulation calculations to compute grain collision cross sections, but use slower piecewise linear integrals
#DL_PRODUCTION                          #creation of dust particles from star particles
#DL_REFINEMENT                          #refinement of large dust particles
#DL_DEREFINEMENT                        #derefinement of small dust particles, 0=derefine into closest dust particle with mass above derefinement mass, 1=derefine into closest dust particle
#DL_SUBCYCLE                            #grain size evolution takes places over multiple subcycles during larger particle time steps
#DL_WINDS                               #stochastic prescription giving dust particles kicks to mimic dust winds
#DL_ONLY_HIGHRES_DUST                   #avoid spawning very high mass dust particles; useful for zoom-in runs, where some star particles may have very high mass; requires DL_REFINEMENT for desired maximum dust mass
#DL_RADIATION                           #main switch to allow dust particles and radiation to interact
#DL_RADIATION_PRESSURE                  #radiation computed from MRT imparts momentum on by dust particles
#DL_RADIATION_ABSORPTION                #dust particles provide opacity to absorb photons
#DL_OUTPUT_RT_FLUX                      #print interpolated radiation field flux in snapshots
#DL_THERMAL_IR                          #if tracking IR radiation, include the thermal coupling between dust and IR radiation


#--------------------------------------- Modified Gas Cooling
#RADCOOL                                # Include the effects of local radiation fields on gas cooling rates, requires GFM, if PMGRID not defined then uncomment GRAVITY_NOT_PERIODIC, currently not compatible with PLACEHIGHRESREGION
#RADCOOL_HOTHALO                        # Include radiation field from HOT GAS, works only with RADCOOL option
#RADCOOL_HOTHALO_METAL_BOOST            # Include an additional boost factor to account for the additional luminosity emitted in emission lines
#TEST_COOLING_METAL                     # call only cooling test routine (save cooling function with metal cooling for solar metallicity)
#EXPLICIT_COOLING                       # switch to a 2nd order explicit method for cooling if (u^{n+1} - u_{n}) < tol * u^{n}


#--------------------------------------- SMUGGLE - Star formation and feedback module
#SMUGGLE_SFR                                #turns on star formation (needs USE_SFR)
#SMUGGLE_STAR_FEEDBACK                      #turns on stellar feedback
#SMUGGLE_STAR_FEEDBACK_TIME_LIMITER         #turns on time step limiter for stellar evolution
#SMUGGLE_VARIABLE_EFFICIENCY                #allows for variation of star formation efficiency based on virial parameter
#SMUGGLE_OUTPUT_SF_PROBABILITY              #enables output of the probability of transforming gas cell into star particle
#SMUGGLE_TEST_SFR                           #only calls the SF initialization and saves Kennicutt law (and the gas effective EOS if available)
#SMUGGLE_USE_POLYTROPIC_EQSTATE             #imposes a minimum temperature to star forming gas (through a polytropic equation of state)
#SMUGGLE_COMPUTE_SFR_FROM_H2                #links the SFR to the H2 gas fraction
#SMUGGLE_OUTPUT_STELLAR_FEEDBACK            #outputs SNII number, feedback energy and mass released for stellar particles and log files for feedback (requires GFM_STELLAR_EVOLUTION)
#SMUGGLE_OUTPUT_MOLECULAR_FRACTION          #outputs the H2 gas fraction (requires SMUGGLE_COMPUTE_SFR_FROM_H2 switched on)
#SMUGGLE_OUTPUT_OPTICAL_DEPTH               #outputs the gas optical depth (requires SMUGGLE_COMPUTE_SFR_FROM_H2 switched on)
#SMUGGLE_OUTPUT_VIRIAL_PARAM                #outputs the gas cell virial parameter
#SMUGGLE_RADPRESS_OPT_THIN                  #adds radiative pressure in optically thin approximation. If GFM active only young stars are considered. Needs OTVET
#SMUGGLE_RADPRESS_OPT_THIN_LUMPERMASS       #source emits at a rate proportional to mass (IonizingLumPerSolarMass in parameterfile). Otherwise constant given by IonizingLumPerSolarMass
#SMUGGLE_RADPRESS_OPT_THICK                 #adds radiation pressure feedback using radiative transfer. Needs OTVET active
#SMUGGLE_RADIATION_FEEDBACK                 #inputs momentum to gas particles within stromgren radius, keep cells at 10^4 K and prevents star formation in them.
#SMUGGLE_RADIATION_FEEDBACK_DEBUG           #extra output fields for SMUGGLE_RADIATION_FEEDBACK
#SMUGGLE_MASS_WEIGHT_SN                     #feedback energy weighted by mass instead of volume
#SMUGGLE_OMEGA_WEIGHT_SN                    #feedback energy weighted by solid angle instead of volume
#SMUGGLE_VAR_SN_EFF                         #SN efficiency scales with neighboring gas metallicity
#SMUGGLE_MOLEC_COOLING                      #approx extra molecular cooling contribution addition based on fit to CLOUDY cooling curves
#SMUGGLE_DUST_HEATING_COOLING               #approx extra gas-dust collisional heating cooling (Meijerink & Spaans 2005)
#SMUGGLE_COSMIC_RAY_HEATING                 #approx cosmic rate heating based on formula of Guo & Oh (2008)
#SMUGGLE_PHOTOELECTRIC_HEATING              #approx photoelectric heating based on formula of Wolfire (2003)
#SMUGGLE_SN_COOLING_RADIUS_BOOST            #returns momentum and energy to the ISM accounting for an unresolved energy conserving early ST blast wave phase
#SMUGGLE_DISCRETE_SN                        #SN feedback is done in discrete SN rather then as a continuous injection
#SMUGGLE_AGB_WINDS                          #returns momentum and energy to the ISM accounting OB and AGB stellar winds
#SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION=0   #stochastic photoionization based on ionizing photon budget (0=each time step; 1=one ionization event per star particle)
#SMUGGLE_SUPERBUBBLE_LIMITER                #sets feedback coupling radius according to estimates of superbubble sizes


#-------------------------------------- Conduction
#MONOTONE_CONDUCTION                   # Monotonicity Preserving anisotropic diffusion solver for thermal conduction
#CONDUCTION_ISOTROPIC                  # Isotropic conduction
#CONDUCTION_ANISOTROPIC                # Anisotropic Conduction
#CONDUCTION_CONSTANT                   # Set Conduction coefficient constant
#CONDUCTION_SATURATION                 # Saturation of Conduction coefficient at low densities
#IMPLICIT_TI                           # Implicit time integration , backwards euler scheme -no limitation of time step
#SEMI_IMPLICIT_TI                      # Semi Implicit Time Integration, stable upto ncfl=4
#RESTRICT_KAPPA                        # Set a maximum value for diffusivity, done in order to avoid very small time steps
#MULTIPLE_TIME_STEPPING                # An approximated implicit time integration scheme which works on multiple time steps
#NON_LINEAR_SLOPE_LIMITERS             # ADIITIONAL SLOPE LIMITERS TO LIMIT NUMERICAL DIFFUSION


#-------------------------------------- Braginskii Viscosity
#BRAGINSKII_VISCOSITY                   # Explicit anisotropic viscosity. Implementation reference: https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2919B/abstract
#BRAGINSKII_VISCOSITY_GLOBAL_TIMESTEP   # Activate together with FORCE_EQUAL_TIMESTEPS
#BRAGINSKII_RKL2_SUPER_TIME_STEPPING    # Second order accurate Runge-Kutta Legendre super timestepping. Unfortunately not compatible with local time steps.
#BRAGINSKII_VISCOSITY_SUBCYCLE          # Subcycle the explicit call. Compatible with local time steps.
#BRAGINSKII_SPITZER                     # Use Spitzer viscosity coefficient (temperature and density dependent). Otherwise a constant is used.
#BRAGINSKII_LIMIT_PRESSURE_ANISOTROPY   # Limit the pressure anisotropy to marginal stability of firehose and mirror instabilities when evaluating the viscous flux.
#BRAGINSKII_OUTPUT_PRESSURE_ANISOTROPY  # Output the (full, unlimited) pressure anisotropy (even if BRAGINSKII_LIMIT_PRESSURE_ANISOTROPY is turned on).
#BRAGINSKII_OUTPUT_NUM_SUBSTEPS         # This is a diagnostic output for RKL2 or subcycling.


#-------------------------------------- MRT
#MRT                           # Moment based RT - currently support M1 closure
#MRT_COMOVING                  # Solve the RT equations in the comoving reference frame
#MRT_SUBCYCLE                  # Subcycle the RT step
#MRT_INIT_IONIZATION           # initialize ionization
#MRT_OUTPUT_FLUX               # output photon fluxes
#MRT_TIME_EXTRAPOLATION        # Include the Runge Kutta time extrapolation - needed for second order convergence
#MRT_FLUX_EXTRAPOLATION        # Extrapolate the photon flux as well - right now induces a bit of noise
#MRT_COOLING_HEATING           # Cooling and Heating
#MRT_RADIATION_PRESSURE        # Include radition pressure as soucre term for the momentum conservation equation
#MRT_INCLUDE_HE                # Include Helium in heating and cooling
#MRT_LSF_GRADIENTS             # Use the least square fit gradient estimates
#MRT_RIEMANN_ROSUNOV           # Use the Rosunov (GLF) riemann solver
#MRT_RIEMANN_ROSUNOV_NEW       # Use the Rosunov (GLF) riemann solver (accounts mesh motion via advection step)
#MRT_RIEMANN_HLLE              # Use the Harten-Lax-van Leer flux function
#MRT_RIEMANN_HLLE_NEW          # Use the Harten-Lax-van Leer flux function (accounts mesh motion via advection step)
#MRT_MULTI_FREQUENCY           # Multi Frequency radiative transfer
#MRT_CHEMISTRY_PS2009          # Petkova & Springel 2009 Chemistry
#MRT_CHEMISTRY_PS2011          # Petkova & Springel 2011 Chemistry
#MRT_COUPLED_THERMOCHEMISTRY   # Solve the coupled chemistry and cooling equations
#MRT_NO_OTSA                   # Do not apply On The Spot Approximation (OTSA)
#MRT_NOCOLLISION_IONIZATION    # No Collisional ionisation
#MRT_SLOWLIGHT                 # Reduce speed of light
#MRT_DUAL_LIGHTSPEED           # Use dual speed of light reduced c for porpogation equations full c for the chemistry network
#MRT_CONSTANT_KAPPA            # Constant Kappa
#MRT_IR                        # Include IR radiative transfer - ala Rosdahl+15
#MRT_IR_ONLY_CHEMISTRY         # Only do absorption terms and not the cooling terms - assume there is no absorption of energy, only the flux is absorbed
#MRT_IR_LTE                    # Do gas heating and cooling according to single fluid gas-dust-IR radiation LTE assumption
#MRT_IR_LTE_SEMI_IMPLICIT      # Follows  Rosdahl+15 iterative approach. Can lead to large number of interative steps (sometimes infinite)
#MRT_IR_LTE_GSL                # Use the GSL provided Implicit Bulirsch-Stoer method of Bader and Deuflhard to solve the energy equation (preferred Method)
#MRT_IR_PHOTON_TRAPPING        # Trap unresolved IR photons - not compatible with gradient extrapolations (hardcoded - no need to turn off time and space extrapolations)
#MRT_IR_GRAIN_KAPPA            # Use grain opacities calculated from Draine & Lee 1984, Laor & Draine 1993, requires local dust-to-gas ratio from GFM_DUST
#MRT_UV_ONLY_DUST              # Let UV radiation only interact with dust
#MRT_NO_UV                     # Do not include UV RT (mainly for testing purposes)
#MRT_SETUP_SPECIAL_BOUNDARIES  # Setup special boundary conditions (mainly for testing purposes)
#MRT_LEVITATION_TEST           # Enable external gravity for levitation test
#MRT_SOURCES=0                 # Enable source treatment for GFM stellar particles and black holes
#MRT_REDUCE_OUTPUT             # less output from MRT
#MRT_EQUIL_CHEM_COOL           # Do equilibrium chemistry + cooling
#MRT_MOLECULAR_COOLING         # Molecular cooling  - fit from grackle
#MRT_PHOTOELECTRIC_HEATING     # Photo electric heating
#MRT_METAL_COOLING             # Add tabulated metal cooling
#MRT_UVB                       # Add contribution from UVB - only z=0 for now
#MRT_UPDATE_AT_END_OF_STEP     # update the hydro primitive varibales only at the end of hydro loop

#-------------------------------------- MRT - STARS
#MRT_STARS                     # Include ionizing photons from GFM stellar particles
#MRT_STARS_EXACT_NGB           # no ngb mass-weighting

#-------------------------------------- MRT - AGN
#MRT_BH                        # Include photons from black hole particles
#MRT_BH_EXACT_NGB              # no ngb mass-weighting
#MRT_BH_PULSED                 # pulsed radiation injection
#MRT_BH_UV_INJECTION           # Inject photons in UV bin(s)
#MRT_BH_IR_INJECTION           # Inject photons in IR bin
#MRT_BH_BIPOLAR                # Inject photons in a bipolar manner
#MRT_BH_BIPOLAR_SET_FLUX       # set photon flux for bipolar injection
#MRT_BH_OMEGA_WEIGHT           # Weight photon injection by solid angle subtended by the cell

#-------------------------------------- MRT - LOCAL FEEDBACK
#MRT_LOCAL_FEEDBACK            # Main switch to couple to local feedback module
#MRT_CHEM_SG                   # Switch to couple to SGchem module
#MRT_INJECT_PHOTONS_EVERY_STEP # Inject photons every step regardless of wether the particle is active or not
#MRT_SINGLE_STAR               # Setup to simulate single star feedback


#-------------------------------------- OTVET IMPLEMENTATION
#OTVET                                 #Master switch
#OTVET_CHEMISTRY_PS2009                #uses Petkova&Springel 2009 original chemical network, with CGmethod solving transport+absorption
#OTVET_CHEMISTRY_PS2011                #uses Petkova&Springel 2011 chemical network, with CGmethod solving only transport
#OTVET_NOGRAVITY                       #builds the tree but does not apply the kicks. Needs improvement to not account optionally for self-gravity of the gas, instead of all gravity like now
#OTVET_NOTMOVEGAS                      #Does not move gas according to flow
#OTVET_OUTPUT_ET                       #outputs the Eddington Tensor for all cells
#EDDINGTON_TENSOR_STARS                #activate stars as sources of radiation
#OTVET_MODIFY_EDDINGTON_TENSOR         #fully anisotropic Eddington Tensor
#OTVET_FLUXLIMITER                     #activate flux limited diffusion, expresion from Petkova & Springel (2009)
#OTVET_CHANGEFLUXLIMITER               #change flux limiter formula to Levermore and Pomraning (1981)
#OTVET_FIXTIMESTEP                     #fixed time step, activate together with FORCE_EQUAL_TIMESTEPS
#OTVET_COOLING_HEATING                 #activate cooling and heating following Petkova & Springel (2009)
#OTVET_MULTI_FREQUENCY                 #relaxes the assumption of monochromatic emission
#OTVET_INCLUDE_HE                      #follow also Helium
#OTVET_SCATTER_SOURCE                  #distributes luminosity in an SPH-way in otvet_Ngb_source neighbour gas cells
#OTVET_OUTPUT_SOURCEHSML               #if OTVET_SCATTER_SOURCE is active, this outputs the HSML and density of sources
#OTVET_CHECK_PHOTONCOUNT               #checks explicity photon conservation during OTVET transport. Only makes sense if OTVET_CHEMISTRY_PS2011 isactive.
#OTVET_NOCOLLISION_IONIZATION          #switch off collisional ionization
#OTVET_SILENT                          # --inactive--
#OTVET_MODIFY_EDDINGTON_TENSOR         # --inactive--
#OTVET_TEST_SST
#EDDINGTON_TENSOR_SFR


#-------------------------------------- TG's switches for primordial simulations
#TGSET                                 #some custom settings
#TGCHEM                                #primordial chemistry and cooling network (Greif 2014)
#TGCHEM_TEST                           #primordial chemistry and cooling network test (Greif 2014)
#HEALRAY                               #adaptive ray-tracing (Greif 2014)
#SINKS                                 #sink particles (under construction)
#SINKS_MERGERS                         #enable mergers for sink particles


#-------------------------------------- Axion dark matter
#BECDM
#BECDM_INPUT_PHASE_AS_VX


#-------------------------------------- SIDM - Self-Interacting DM
#SIDM=2                                #activate and set types
#SIDM_CONST_CROSS                      #constant cross section
#SIDM_STATES=2                         #number of DM states (for inelastic models)
#SIDM_REACTIONS=5                      #number of scatter reactions (for inelasitc models)
#SIDM_NO_SCATTER                       #DEBUG: no scattering at all
#SIDM_NO_TIMESTEP                      #DEBUG: do not change time step
#SIDM_NO_KINEMATICS                    #DEBUG: do not change particle velocities, but still run through full scattering process
#SIDM_NO_NGB_SEL                       #DEBUG: take closest particle to scatter with; will select particles in wrong state for multiple states; scatter state check is then turned off
#SIDM_NO_MASSCHANGE                    #DEBUG: no mass change during inelastic scattering
#SIDM_NO_ENERGYCHANGE                  #DEBUG: no energy change during inelastic scattering   NOTE: to get fully inelastic behavior turn on SIDM_NO_MASSCHANGE and SIDM_NO_ENERGYCHANGE

#-------------------------------------- On-the-fly shock finder
#SHOCK_FINDER_BEFORE_OUTPUT             #Use this flag if you want to run the shock finder before a snapshot dump, no additional flags or parameters needed.
#SHOCK_FINDER_ON_THE_FLY                #Run the shock finder at every local time step, no additional flags or parameters needed.


#--------------------------------------- Post-processing shock finder, please read the instructions in shock_finder.h
#SHOCK_FINDER_POST_PROCESSING           #post-processing shock finder
#SHOCK_FINDER_AREPO                     #standard operating mode
#UNLIMITED_GRADIENTS                    #standard option
#ZONE_JUMP_P                            #standard option
#ZONE_JUMP_T                            #standard option
#SHOCK_DIR_GRAD_T                       #standard option
#SHOCK_JUMP_T                           #standard option
#SURFACE_SPHERE_APPROX                  #use this for 2d sims
#SURFACE_ANGLE_APPROX                   #use this for 3d sims
#RESET_WRONG_JUMPS                      #standard option
#RESET_WRONG_RHO_JUMPS                  #standard option
#RESET_WRONG_P_JUMPS                    #standard option
#SKIP_BORDER                            #for non-periodic boundaries of the snapshot/subbox


#--------------------------------------- atomic dark matter (in fluid approximation)
#ATOMIC_DM                              # master switch


#--------------------------------------- Binary stellar systems
#BINARYLOG
#BINARYLOG_FOR_MERGERS                  # Output stats based on passive scalars. Requires BINARYLOG
#SPECIAL_SOFTENINGS
#REDUCE_SOFTENINGS
#ID_RGCORE=1000000000
#DMLOWESTTIMEBIN
#DMFIXED


#--------------------------------------- FLD
#FLD

#FLD_CONES
#FLD_NCONES=12

#FLD_CONST_KAPPA
#FLD_MARSHAK

#FLD_HYPRE
#FLD_HYPRE_IJ1
#FLD_HYPRE_IJ2
#HYPRE_PCG

#FLD_MG
#FLD_MG_GS

#FLD_ANISOTROPIC_CIRCULAR
#FLD_NO_TEMP_UPDATE
#FLD_SILENT

#FLD_TEST_BOUNDARY
#FLD_UPPER_BOUNDARY_MINID=1
#FLD_UPPER_BOUNDARY_MAXID=1000000
#FLD_LOWER_BOUNDARY_MINID=5000000
#FLD_LOWER_BOUNDARY_MAXID=6000000


#--------------------------------------- Calculate quantities in post-processing
#CALCULATE_QUANTITIES_IN_POSTPROCESS
#POWERSPECTRUM_IN_POSTPROCESSING
#POWERSPECTRUM_IN_POSTPROCESSING_ICS


#--------------------------------------- Discontinuous Galerkin (DG)
#DG                                    # master switch
#DG_SET_IC_FROM_AVERAGES               # loads ordinary non-DG ICs and sets initial weights according to density, velocity and internal energy
#DG_TEST_PROBLEM                       # initial conditions are created from file src/dg/test_problems.c
#DG_VERBOSE                            # additional output for debugging
#DG_DEBUG                              # run in debug mode, additional checks are active
#DEGREE_K=1                            # spatial degree of the scheme
#CALC_QUADRATURE_DATA                  # calculate the quadrature data instead of reading it from a table
#RIEMANN_HLLC                          # use the hllc riemann solver instead of the normal one

#MINMOD_B                              # use the slope limiter as a total variaton bounded limiter
#DISCONTINUITY_DETECTION               # limit only when a discontinuity is found
#OUTPUT_DG_DISCONTINUITIES
#OUTPUT_DG_INFLOW_BOUNDARIES
#ANGLE_BOUND                           # use the angle bound methdod
#CHARACTERISTIC_LIMITER                # limit the characteristic variables instead of the conserved variables
#CONSERVED_LIMITER                     # limit the conserved variables
#POSITIVITY_LIMITER                    # keep the cell average values positive
#FIX_MEAN_VALUES                       # reset negative values


#--------------------------------------- Spiral potential as used by Dobbs
#SPIRAL


#--------------------------------------- Supernova Energy or Momentum cons.
#SNE_FEEDBACK
#CLUSTERED_SNE
#INJECT_TRACER_INTO_SN
#SNE_solarring


#--------------------------------------- Deprecated
#CONDUCTION
#CONDUCTION_CRANK_NICOLSON
#CONDUCTION_PCG
#MAX_COND_ITER=400
#COND_ITER_ACCURACY=1.0e-10

#WG15_INTPL
#SUNRISE


#--------------------------------------- Grackle
#GRACKLE                               #master switch
#GRACKLE_H2                            #Turn on H2 cooling and chemistry
#GRACKLE_D                             #Turn on Deuterium cooling and chemistry
#GRACKLE_TAB                           #Run in tabulated mode
#GRACKLE_ABUNDANCE_IN_ICS              #Use abundances in ICs instead of converging on startup
#GRACKLE_PHOTOELECTRIC                 #Use global volumetric heating rate
#GRACKLE_TEMPERATURE_FLOOR             #Obey temperature floor set by MinGasTemp or MinEgySpec (only with GRACKLE_TAB)
#GRACKLE_VERBOSE                       #Turn on Grackle internal print statements
#GRACKLE_UNCONVERGED_IGNORE            #On startup, proceed even if Grackle failed to find an equilibrium within iter limit. 


#--------------------------------------- Modified gravity solver (Private to Ewald Puchwein, Volker Springel and Christian Arnold)
#MODGRAV                               #master switch
#MODGRAV_EFF_MASS                      #use the effective mass scheme to obtain the forces (currently the only method implemented)
#MODGRAV_INTERPOLATE_PHI

#BAROTROPIC
#BARO_CONSTANT_GAMMA_EOS
#BHATTAL98


#--------------------------------------- SINK particles
#SINK_PARTICLES
#DUMP_SINK_PARTICLE_INFO
#SINK_PARTICLES_SKIM_CELL_MASS
#SINK_PARTICLES_VARIABLE_ACC_RADIUS
#SINK_PARTICLES_VARIABLE_CREATION
#SINK_PARTICLES_LIMIT_TIMESTEP
#SINK_PARTICLES_FEEDBACK
#SINK_PHOTOION_FEEDBACK
#SINK_PARTICLES_FEEDBACK_RETURN_MASS
#SINK_FEEDBACK_SINGLE_STAR
#SINK_PARTICLES_REFINEMENT_LIMIT
#SINK_PARTICLES_FORCE_FORMATION
#SINK_PARTICLES_OUTPUT_EVERY_NEW_SINK
#ALLOW_MULTIPLE_SINK_CREATION_PER_TIMESTEP
#DEBUG_SINK_PARTICLES
#SINK_PARTICLE_FREE_FALL_TEST
#SINK_MERGERS
#SINK_MERGERS_DEBUG
#STORE_SINK_PARTICLE_SPIN


#--------------------------------------- TreeColV2
#TREECOLV2
#TREECOLV2_C
#TREECOLV2_CO
#TREECOLV2_H2
#TREECOLV2_VEL
#TREECOLV2_NO_GAS_SELFGRAVITY
#TREECOLV2_DEBUG
#NSIDE=2
#OUTPUTCOL

#TURBULENT_METALDIFFUSION
#TURBULENT_METALDIFFUSION_EXPLICIT


#--------------------------------------- Model for subgrid-scale turbulence
#SGS_TURBULENCE
#SGS_TURBULENCE_IN_ICS
#SGS_TURBULENCE_VIEW_CELLS_AS_CUBES
#SGS_TURBULENCE_VIEW_CELLS_AS_SPHERES
#SGS_TURBULENCE_EDDY_VISCOSITY
#SGS_TURBULENCE_CONSTANT_EDDY_VISCOSITY
#SGS_TURBULENCE_RIEMANN_PRESSURE
#SGS_TURBULENCE_STRESS_TENSOR
#SGS_TURBULENCE_TURBULENT_PRODUCTION
#SGS_TURBULENCE_VISCOUS_DISSIPATION
#SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION
#SGS_TURBULENCE_MAKE_PRODUCTION_DISSIPATION_SPECTRUM
#SGS_TURBULENCE_OUTPUT_SGS_PRESSURE
#OUTPUT_SGS_T_PRESSURE_GRADIENT
#OUTPUT_DENSTROPHY


#--------------------------------------- external Galaxy potential
#GALPOT                               


#--------------------------------------- SGChem chemistry module
#SGCHEM
#CHEMISTRYNETWORK=5
#MCMA
#ABHE
#CHEM_IMAGE
#IMAGE_FOOTERS
#SGCHEM_VARIABLE_Z                     #Allow metallicity and dust-to-gas ratio to vary between different cells
#SGCHEM_VARIABLE_ISRF                  #Allow interstellar radiation field strength to vary spatially
#SGCHEM_VARIABLE_CRION                 #Allow cosmic ray ionization rate to vary spatially
#SGCHEM_TEMPERATURE_FLOOR
#SGCHEM_ACCRETION_LUMINOSITY
#SGCHEM_NO_HIGHN_DCHEM
#SGCHEM_DUMP_THERMAL_RATES
#SGCHEM_NO_MOLECULES
#SGCHEM_NO_COOL
#CHEMCOOL
#NO_CO_COOL
#ONLY_DUST_COOL
#NO_DUST_H2_HEAT
#NO_GAS_DUST_HEAT
#ADIABATIC_DENSITY_THRESHOLD
#_13CO_
#COLLECT_EQB_DATA
#OLD_CO_SELF_SHIELDING
#THERMAL_INFO_DUMP
#DEBUG_SGCHEM
#DEBUG_COOLING_ID
#DEBUG_PARTICLE_ID
#DEBUGDUSTID
#DEBUG_RATE_EQ
#DEBUG_EVOLVE


#--------------------------------------- CHIMES Chemistry module
#CHIMES 
#CHIMES_JEANS_SHIELDING 
#CHIMES_SOBOLEV_SHIELDING
#CHIMES_PREENRICH_AT_START
#CHIMES_PTHREADS
#CHIMES_INITIALISE_IN_EQM
#CHIMES_ADVECT_ABUNDANCES
#CHIMES_DISABLE_SHIELDING_ON_EOS
#CHIMES_DISABLE_ZCOOL_ON_EOS
#CHIMES_REDSHIFT_DEPENDENT_UVB



#--------------------------------------- AGB_WIND
#AGB_WIND

#------------------------------------- SFR_MCS (Matthew C. Smith, private but collaboration encouraged)
#SFR_MCS                         #master switch, requires USE_SFR
#SFR_MCS_SELECT_CRITERIA=0       #Criteria that selects SF gas (see sfr_criteria_announce())
#SFR_MCS_RATE_CRITERIA=0         #SF rate deterination method (see sfr_criteria_announce())
#IMF_SAMPLING_MCS                #Explicitly sample the IMF instead of using Starburst99 properties
#N_STAR_SLOTS=8                  #Number of available slots for massive stars per star particle
#SFR_MCS_LOG                     #Histogram of SF site densities recorded and output
#SFR_MCS_LOG_N=200               #Number of density bins in log10 H/cc
#SFR_MCS_LOG_MIN=0               #log10 min density in H/cc
#SFR_MCS_LOG_MAX=8               #log10 max density in H/cc
#SFR_MCS_BIRTH_RECORDS           #Information about individual star particle birth written to snapshot
#SFR_MCS_LOG_DETAILS             #Information about individual star particle birth written to output file
 
#JEANS_PRESSURE_LIMIT_MCS=8      #Pressure floor, enforce resolution of Jeans length by N cells
#JEANS_MASS_PRESSURE_LIMIT_MCS   #Instead enforce resolution of Jeans mass (instead of length)
#SMAUG_PRESSURE_FLOOR            #Ramses-like pressure floor formulation

#SN_MCS                          #SN master switch, requires SFR_MCS
#SN_MCS_MECHANICAL               #Use mechanical feedback scheme
#SN_MCS_NNGB_MAX=64              #Buffer size for neighbours (and terminate if exceeded)
#SN_MCS_NO_ENERGY                #Mass return etc. but no energy/momentum injection
#SN_MCS_CHANCE=1                 #SN Energy boosted by SN_MCS_CHANCE, occur with prob. 1/SN_MCS_CHANCE
#SN_MCS_LOG                      #Histogram of SN site densities
#SN_MCS_LOG_N=200                #Number of density bins in log10 H/cc
#SN_MCS_LOG_MIN=0                #log10 min density in H/cc
#SN_MCS_LOG_MAX=8                #log10 max density in H/cc
#SN_MCS_LOCATION_RECORDS         #Information about individual SNe written to snapshot
#SN_MCS_LOG_DETAILS              #Information about individual SNe written to output file
#SN_MCS_VARIABLE_EJECTA          #Ejecta mass and metallicity also interpolated from tables
#SN_MCS_PROMPT                   #Only works with IMF_SAMPLING_MCS, SNe explode immediately
#SN_MCS_SINGLE_INJECTION         #One SN event per star particle after fixed delay time.                      
#SN_MCS_INITIAL_DRIVING          #Inject SNe directly from gas, used to settle an isolated simulation

#HII_MCS                         #Photoionization using overlapping Stromgren approximation
#HII_MCS_ANISO                   #Solve Stromgren radius using independent healpix bins
#HII_MCS_LR                      #Long range photoionisation/photoheating, tree based
#HII_MCS_EVERY_SYNCPOINT         #Recalculate Hii regions every timestep
#HII_MCS_DENSCUT                 #Exclude low density gas from HII calculation
#HII_MCS_LOG                     #Turn on log file

#PE_MCS                          #Photoelectric heating, ISRF calculated using tree
#PE_MCS_PRESHIELD=0              #Attenuation at star particle. 0 = none, 1 = Jeans length 
#PE_MCS_POSTSHIELD=0             #Attenuation at gas cell. 0 = none, 1 = Jeans length
#PE_MCS_LINEAR_DGR               #Make dust to gas ratio depend linearly on metallicity instead of broken power law
#PE_MCS_FIXED_EPS                #Use a fixed efficiency supplied in parameterfile instead of varying with density

#TURB_APPROX_MCS                 #ILES model, tracking generation and dissipation of unresolved turbulence
#TURB_APPROX_MCS_GRAD_UNLIM      #Use unlimited velocity gradients
#TURB_APPROX_MCS_RENORM          #Renormalise turb. energy as filter scale (from volume) changes
#TURB_APPROX_MCS_SEED            #On startup, initialise turbulent energy from velocity shear
#TURB_APPROX_MCS_KEEP_ICS        #On startup, use the turbulent energy from the ICs
#TURB_APPROX_MCS_OUTPUT_RATES    #Output related quantities to snapshots

#--------------------------------------- SPRAI (SimpleX Radiation Transfer)
#SIMPLEX
#SX_CHEMISTRY=3
#SX_NDIR=84
#SX_SOURCES=10
#SX_NUM_ROT=5
#SX_HYDROGEN_ONLY
#SX_DISPLAY_STATS
#SX_DISPLAY_TIMERS
#SX_OUTPUT_IMAGE
#SX_OUTPUT_IMAGE_ALL
#SX_OUTPUT_FLUX


#--------------------------------------- MAGNETIC FIELD SEEDING
#BIERMANN_BATTERY                      # activate Biermann (1950) battery
#DURRIVE_BATTERY                       # activate Durrive & Langer (2015)
#GFM_INJECT_B_FROM_SN                  # inject a dipolar magnetic field during SN events
#MAGNETIC_BATTERIES_OUTPUT_GRADIENTS   # output gradients of electron number density, pressure, and momentum transfer rate

#--------------------------------------- STELLAR RADIATIVE TRANSFER
#SOLAR
#SOLAR_RADIATIVE_TRANSFER_DIFF
#SOLAR_RADIATIVE_TRANSFER_EDD
#OUTPUT_QRAD
