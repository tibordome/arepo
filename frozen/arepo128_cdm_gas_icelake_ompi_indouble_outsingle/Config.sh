#!/bin/bash            # this line only there to enable syntax highlighting in this file

PMGRID=128 
INPUT_IN_DOUBLEPRECISION





##################################################
#  Enable/Disable compile-time options as needed #
##################################################

#--------------------------------------- Basic operation mode of code
NTYPES=6                                 # number of particle types

#GENERIC_ASYNC                           # enables asynchronous communication scheme
PERIODIC



MHD
MHD_POWELL
MHD_SEEDFIELD
MHD_POWELL_LIMIT_TIMESTEP

COOLING
UVB_SELF_SHIELDING                # gas is self-shielded from the cosmic background based on its density
USE_SFR


#--------------------------------------- Mesh Type
VORONOI


#--------------------------------------- Riemann solver
RIEMANN_HLLD



#--------------------------------------- Mesh motion and regularization
REGULARIZE_MESH_CM_DRIFT
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED
REGULARIZE_MESH_FACE_ANGLE


#--------------------------------------- Time integration options
TREE_BASED_TIMESTEPS ###         # non-local timestep criterion (take 'signal speed' into account)
###FORCE_EQUAL_TIMESTEPS


#--------------------------------------- Image generation
#VORONOI_IMAGES_FOREACHSNAPSHOT
#VORONOI_FREQUENT_IMAGES                     # creates images with frequency 'TimeBetweenImages' given in parameterfile, independent of snapshots  
#VORONOI_FIELD_DUMP_PIXELS_X=1536
#VORONOI_FIELD_DUMP_PIXELS_Y=150
#VORONOI_VELOCITY_FIELD_2D
#VORONOI_FIELD_COMPENSATE_VX=4.0
#VORONOI_FIELD_COMPENSATE_VY=0
#VORONOI_NEW_IMAGE
#VORONOI_PROJ_TEMP                           #project T instead of u
#VORONOI_PROJ                                # do projection along any predefined direction
#VORONOI_MULTIPLE_PROJECTIONS                # do face-on and edge-on projections by swapping y and z axes




#--------------------------------------- Refinement and derefinement
REFINEMENT_SPLIT_CELLS
REFINEMENT_MERGE_CELLS



#--------------------------------------- Gravity treatment
SELFGRAVITY                       # switch on for self-gravity         
HIERARCHICAL_GRAVITY ###             # use hierarchical splitting of the time integration of the gravity
CELL_CENTER_GRAVITY ###              # uses geometric centers to calculate gravity of cells, only possible with HIERARCHICAL_GRAVITY
ALLOW_DIRECT_SUMMATION
DIRECT_SUMMATION_THRESHOLD=5000 #8000 for L75n1820: can be changed along the way

ENFORCE_JEANS_STABILITY_OF_CELLS_EEOS # was on for Illustris, check for resolution convergence; if below is on, this is on too.
ENFORCE_JEANS_STABILITY_OF_CELLS        # this imposes an adaptive floor for the temperature

EVALPOTENTIAL                     # computes gravitational potential




#--------------------------------------- Gravity softening
NSOFTTYPES=4                      # Number of different softening values to which particle types can be mapped.
MULTIPLE_NODE_SOFTENING           # If a tree node is to be used which is softened, this is done with the softenings of its different mass components
INDIVIDUAL_GRAVITY_SOFTENING=32  # bitmask with particle types where the softenig type should be chosen with that of parttype 1 as a reference type
ADAPTIVE_HYDRO_SOFTENING



#--------------------------------------- TreePM Options
#PMGRID=2048 		#4096 for L75n1820
#ASMTH=1.25
RCUT=5.0

#PLACEHIGHRESREGION=2
#ENLARGEREGION=1.1
#GRIDBOOST=2

FFT_COLUMN_BASED ###
#PM_ZOOM_OPTIMIZED


#--------------------------------------- Things that are always recommended
CHUNKING                     # will calculated the gravity force in interleaved blocks. This can reduce imbalances in case multiple iterations due to insufficient buffer size need to be done
                              


#---------------------------------------- Single/Double Precision
DOUBLEPRECISION=1
DOUBLEPRECISION_FFTW

#OUTPUT_COORDINATES_IN_DOUBLEPRECISION # to be implemented

NGB_TREE_DOUBLEPRECISION  # if this is enabled, double precision is used for the neighbor node extension



#---------------------------------------- On the fly FOF groupfinder
FOF                                    # enable FoF output
FOF_PRIMARY_LINK_TYPES=2               # 2^type for the primary dark matter type
FOF_SECONDARY_LINK_TYPES=1+16+32   # 2^type for the types linked to nearest primaries
#FOF_SECONDARY_LINK_TARGET_TYPES=   # should normally be set to a list of all dark matter types (in zoom runs), if not set defaults to FOF_PRIMARY_LINK_TYPES
#FOF_GROUP_MIN_LEN=32                   # default is 32
#FOF_LINKLENGTH=0.16                    # Linkinglength for FoF (default=0.2)
#FOF_FUZZ_SORT_BY_NEAREST_GROUP=0   # sort fuzz particles by nearest group and generate offset table in catalog (=1 writes nearest group number to snapshot)
#FOF_STOREIDS                           # store IDs in group/subfind catalogue, do not order particles in snapshot files by group order
#USE_AREPO_FOF_WITH_GADGET_FIX          # Needed in order to run FOF with Arepo on Gadget snapshot files, if gas is present and should be linked to the FOFs
#ADD_GROUP_PROPERTIES                   # This can be used to calculate additional properties for an already existing group catalogue. These are then added as additional columns to the HDF5 group catalogues.
#ADD_SO_GROUP_PROPERTIES                # This can be used to calculate additional properties for an already existing group catalogue. These are then added as additional columns to the HDF5 group catalogues.


#---------------------------------------- Subfind
SUBFIND                                # enables substructure finder
SAVE_HSML_IN_SNAPSHOT                  # this will store hsml and density values in the snapshot files


#SUBFIND_MEASURE_H2MASS                 # special measurement option for mass in molecular hydrogen
SUBFIND_CALC_MORE                      # calculates also the velocity dispersion in the local density estimate
SUBFIND_EXTENDED_PROPERTIES            # adds calculation of further quantities related to angular momentum in different components


#--------------------------------------- SFR/feedback model


SOFTEREQS
#MODIFIED_EOS
#SLOW_RELAX_TO_EOS


#-------------------------------------- AGN stuff
###BLACK_HOLES                   # enables Black-Holes (master switch)
###BH_THERMALFEEDBACK            # quasar-mode: couple a fraction of the BH luminosity into surrounding
###DRAINGAS=3                   # non-stochastic smooth accretion (1: on, 2: on + cell rho, 3: on + gas drained from all cells within hsml) In Illustris was 2, but with 3 accretion is more stabel at low res
###BH_EXACT_INTEGRATION          # integrates analytically mass accretion
###BH_BONDI_DEFAULT              # default Bondi prescription
###BH_DO_NOT_PREVENT_MERGERS # When this is enabled, BHs can merge irrespective of their relative velocity
####BH_USE_GASVEL_IN_BONDI        # only when this is enabled, the surrounding gas velocity is used in addition to the sounds speed in the Bondi rate
###BH_USE_ALFVEN_SPEED_IN_BONDI  # when this is enabled the alfven speed is added to the gas sound speed in the Bondi rate and the total gas pressure around the BH includes the magnetic contribution when compared the the reference pressure in BH_PRESSURE_CRITERION (requires MHD)
####BH_FRICTION                   # estimates the local DM density around BH and applies a friction force to the relative velocity, meant as a replacement for REPOSITION_ON_POTMIN
####BH_FRICTION_AGGRESSIVE
###BH_NEW_CENTERING
####BH_DRAG                       # Drag on black-holes due to accretion: current implementation simply double-accounts for the accretion momentum transfer, and should therefore not be used
####REPOSITION_ON_POTMIN          # repositions hole on potential minimum (requires EVALPOTENTIAL)
###BH_PRESSURE_CRITERION
####BH_RELATIVE_NGB_DEVIATION # Maximum NGB number deviation calculated relative to total number of neighbours
####OUTPUT_BLACK_HOLE_TIMESTEP #outputs the 3 time-steps for BH particles


#-------------------------------------- other AGN stuff
#UNIFIED_FEEDBACK            # activates BH_THERMALFEEDBACK at high Mdot and BH_BUBBLES FEEDBACK al low Mdot (-->OBSOLETE: replaced by BH_NEW_RADIO)
#BH_BUBBLES                  # calculate bubble energy directly from the black hole accretion rate (-->OBSOLETE: replaced by BH_NEW_RADIO)
#BH_MAGNETIC_BUBBLES         # inject part of the  bubble energy as magnetic energy
#BH_MAGNETIC_DIPOLAR_BUBBLES #inject part of the bubble energy as magnetic energy, field arranged as a dipole with random orientation
###BH_ADIOS_WIND
###BH_ADIOS_RANDOMIZED
###BH_ADIOS_WIND_WITH_QUASARTHRESHOLD  # use a threshold value ("qusarthrehold") of bondi-rate over Eddington rate to decide about quasar mode vs. adios wind
###BH_ADIOS_WIND_WITH_VARIABLE_QUASARTHRESHOLD  # scales the threshold with black hole mass (with a factor (M_BH/M_ref)^2, where M_ref = 10^8 Msun)
###BH_ADIOS_ONLY_ABOVE_MINIMUM_DENSITY
####BH_ADIOS_DENS_DEP_EFFICIANCY
####BH_CONTINOUS_MODE_SWITCH # calculates fraction of thermal and mechanical feedback energy depending on eddington factor and mass (continously in both quantities)




#---------------------------------------- Passive Tracers
#TRACER_FIELD
#TRACER_PARTICLE=2                     # advect massless tracer particles of type TRACER_PARTICLE with velocity field

###TRACER_MC=3                       # Monte Carlo tracer particles (=1 to enable, >=2 to output as that partType)
#GENERATE_TRACER_PARTICLE_IN_ICS   # add tracer particles at positions of cell vertices in ICs
###GENERATE_TRACER_MC_IN_ICS             # add a fixed number (given in the parameter file) of MC tracers to each gas cell in ICs


#TRACER_MC_NUM_FLUID_QUANTITIES=0                      # number of fluid quantities to be stored for MC tracers - must match the value gien to TRACER_MC_STORE_WHAT
#TRACER_MC_STORE_WHAT=0         # bit mask for quantities to store (see allvars.h for bitmask)


#TRACER_MC_SKIPLOAD=3                       # skip reading this particle type when reading initial conditions from a snapshot file
#FOF_DISABLE_SNAP_REWRITE                   # do not rewrite a snapshot file when RestartFlag==3
#TRACER_TRAJECTORY
#TRACER_TRAJECTORY_GENERATE




#-------------------------------------------- Things for special behaviour
#READ_DM_AS_GAS
#NO_ID_UNIQUE_CHECK
#RUNNING_SAFETY_FILE                    # if file './running' exists, do not start the run
#LOAD_TYPES=1+2+4+16+32
#READ_COORDINATES_IN_DOUBLE
#IDS_OFFSET=1           #offset for gas particles if created from DM
#TILE_ICS
#COMBINETYPES            # reads in the IC file types 4+5 as type 3 (useful for doing gas runs of Aquarius ICs)
#USE_RANDOM_GENERATOR
#MULTIPLE_RESTARTS
#TOLERATE_WRITE_ERROR
#OPTIMIZE_MEMORY_USAGE                       #optimize for memory, not for speed. Note: this is dangerous for high dynamic range simulations with mixed precision, since some position variables are singles instead of doubles
###SUBBOX_SNAPSHOTS ### # note: maybe one only subbox0 of Illustris, better time spacing, smaller number of chunks (to be made a new param separated from the standard snaps)
PROCESS_TIMES_OF_OUTPUTLIST
#EXTENDED_GHOST_SEARCH               # This extends the ghost search to the full 3x3 domain instead of the principal domain
 
#DOUBLE_STENCIL                     # this will ensure that the boundary region of the local mesh is deep enough to have a valid double stencil for all local cells


#TETRA_INDEX_IN_FACE                # adds an index to each entry of VF[] and DC[] to one of the tetrahedra that share this edge


VORONOI_DYNAMIC_UPDATE              # keeps track of mesh connectivity, which speeds up mesh construction
#COFFEE_PROBLEM
#NOH_PROBLEM
#SHIFT_BY_HALF_BOX
#DISABLE_VELOCITY_CSND_SLOPE_LIMITING
#NO_MPI_IN_PLACE
#NO_ISEND_IRECV_IN_DOMAIN
#FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
#USE_MPIALLTOALLV_IN_DOMAINDECOMP
#MPI_HYPERCUBE_ALLGATHERV         # some MPI-libraries may use quite a bit of internal storage for MPI_Allgatherv. This uses hypercubes instead as a work-around
#MPISENDRECV_CHECKSUM
#NOTREERND
ENLARGE_DYNAMIC_RANGE_IN_TIME   # This extends the dynamic range of the integer timeline from 32 to 64 bit Likely for the smaller box IllustrisDwarf
#NOSTOP_WHEN_BELOW_MINTIMESTEP
#DO_NOT_CREATE_STAR_PARTICLES
#DMPIC                              # enable special image code for dark matter simulations   
#ALLOWEXTRAPARAMS
#RADIATIVE_RATES                   # used in non-equilibrium chemistry model
FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES  # this can be used to load SPH ICs that contain identical particle coordinates
#VEL_POWERSPEC                     # compiles in a code module that allows via restart-flag 7 the calculation of a gas velocity power spectrum of a snapshot
#ADJ_BOX_POWERSPEC             # compiles in a code module that allows via restart-flag 7 the calculation of gas power spectra of a snapshot with an adjustable box (user defined center and size)
#DISABLE_OPTIMIZE_DOMAIN_MAPPING


#CUDA                           # enables CUDA support in Arepo
#CUDA_INSTRUMENT                # This enables instrumentation support for the nvidia profiler
#GPU_TREE                       # uses the GPU to calculate the tree grav interactions
#GPU_TREE_CALC_CPU              # does the computation on the CPU instead (for debugging)
#GPU_TREE_VERBOSE               # output more informations


#GPU_PM                         # uses the GPU for the FFTs of the PM force


#USE_DSDE                       # try to use a dynamic sparse data exchange paradigm to get rid off sparse MPI_Alltoall patterns on large partitions
#USE_NBC_FOR_IBARRIER           # use the NBC library to implement non-blocking collectives (only relevant when USE_DSDE is used)


#--------------------------------------- Output/Input options
#UPDATE_GRADIENTS_FOR_OUTPUT
REDUCE_FLUSH
#OUTPUT_REFBHCOUNTER                 
#OUTPUT_EVERY_STEP
#GODUNOV_STATS
OUTPUT_CPU_CSV
#OUTPUT_TASK
#OUTPUT_TIMEBIN_HYDRO
#OUTPUT_PRESSURE_GRADIENT
#OUTPUT_DENSITY_GRADIENT
#OUTPUT_VELOCITY_GRADIENT
#OUTPUT_BFIELD_GRADIENT
#OUTPUT_VERTEX_VELOCITY
#OUTPUT_VERTEX_VELOCITY_DIVERGENCE  # requires CALCULATE_VERTEX_VELOCITY_DIVERGENCE
OUTPUT_CENTER_OF_MASS
#OUTPUT_SURFACE_AREA
#OUTPUT_PRESSURE
OUTPUTPOTENTIAL
#OUTPUTACCELERATION
#OUTPUTTIMESTEP
#OUTPUT_SOFTENINGS                # output particle softenings
#OUTPUTGRAVINTERACTIONS           # output gravitatational interactions (from the tree) of particles
HAVE_HDF5                         # needed when HDF5 I/O support is desired
#PARAMS_IN_SNAP                    # add the compiler flags and parameter file values to every snapshot file (requires HAVE_HDF5)
#HDF5_FILTERS                      # activate snapshot compression and checksum for HDF5 output
#OUTPUT_XDMF                       #writes an .xmf file for each snapshot, which can be read by visit (with the hdf5 snapshot)
#OUTPUTCOOLRATE                    # outputs cooling rate, and conduction rate if enabled
#OUTPUT_DIVVEL                             # output  velocity divergence
#OUTPUT_CURLVEL                     # output  velocity curl
#OUTPUT_COOLHEAT                   # output actual energy loss/gain in cooling/heating routine
#OUTPUT_VORTICITY
#OUTPUT_CELL_SPIN                  
#MEASURE_DISSIPATION_RATE          # measures and outputs dissipation rate. Note: requires USE_ENTROPY_FOR_COLD_FLOWS, even though it will then always use the thermal energy update
#OUTPUT_MACHNUM                    # output maximum mach number of a cell


#--------------------------------------- Testing and Debugging options
#DEBUG                             # enables core-dumps
#DEBUG_ENABLE_FPU_EXCEPTIONS   # tries to enable FPU exceptions
#CHECKSUM_DEBUG
#RESTART_DEBUG
#VERBOSE                           # reports readjustments of buffer sizes
HOST_MEMORY_REPORTING             # reports after start-up the available system memory by analyzing /proc/meminfo
#FORCETEST=0.001                   # calculates for given fraction of particles direct summation forces to check accuracy of tree force
#FORCETEST_TESTFORCELAW=1          # this enables a special test to measure the effective force law of the code, can be set to 1 or 2




#--------------------------------------- Glass making/ 2nd-order initial conditions / Initial conditions options
#SECOND_ORDER_ICS
#LONGIDS ### XXX for 2048+ run only
#OFFSET_FOR_NON_CONTIGUOUS_IDS
###GENERATE_GAS_IN_ICS###
###SPLIT_PARTICLE_TYPE=2 ###
#NTYPES_ICS=6 # number of particle types in ICs, if not NTYPES (only works for 6, and non-HDF5 ICs!)





#-------------------------------------- GFM - Galaxy Formation Module
GFM                                                        #master switch
GFM_STELLAR_EVOLUTION=0                    #stellar evolution: 0->default, 1->no mass loss (beta value changes + MassMetallicity & MassMetals inconsistent internally with cell dynamical mass) 2->call only test routine
GFM_CONST_IMF=0                            #0 for Chabrier (default), 1 for a pure power-law (requires parameter IMFslope, e.g. -2.35 for Salpeter)
#GFM_VARIABLE_IMF=0                         #0 for a pure power-law that depends on DM-veldisp
GFM_PREENRICH                              #pre enrich gas at given redshift
#GFM_EXACT_NUMNGB                           #use direct neighbor count instead of kernel weighted neighbor count
GFM_WINDS                                  #decoupled ISM winds
###GFM_WINDS_VARIABLE=1                   #decoupled ISM winds: 0->scale winds with halo mass, requires FoF, 1->sigma winds
###GFM_WINDS_VARIABLE_HUBBLE                  #add an additional H(z)^(-1/3) factor to the wind scaling, such that it scales with halo mass not halo velocity dispersion
###GFM_WIND_ENERGY_METAL_DEPENDENCE           #this can be used to decrease the wind energy for high metallicity (mimicking higher cooling losses)
###GFM_WINDS_STRIPPING                        #wind metal stripping
###GFM_WINDS_THERMAL_NEWDEF                          #not only give the wind kinetic energy but also thermal energy
GFM_COOLING_METAL                          #metal line cooling
GFM_AGN_RADIATION                          #cooling suppression/heating due to AGN radiation field (proximity effect)
GFM_STELLAR_PHOTOMETRICS                   #calculate stellar magnitudes for different filters based on GALAXEV/BC03
GFM_OUTPUT_MASK=1+2+4+8+16+32+64+256   #which fields to output (search GFM_OUTPUT_MASK in io.c to see which fields the bits encode)
#GFM_DUST                                   #formation and evolution of dust, requires GFM_STELLAR_EVOLUTION
#GFM_DUST_DESTMODE=0                        #dust destruction mode: 0->default (uses supernova rate), 1->constant destruction timescale
#GFM_CHECKS                                 #this checks the consistency of the AuxDataID/PID indices of stars and black holes every timestep
#GFM_DISCARD_ENRICHMENT_GRADIENTS           #this disables the gradient extrapolation of the passively advected metallicity scalar variables
GFM_NORMALIZED_METAL_ADVECTION             #this introduces an additional pseudo element for all untracked metals and normalizes the extrapolated abundance vectors to unity
GFM_OUTPUT_BIRTH_POS

GFM_CHEMTAGS
GFM_SPLITFE
GFM_RPROCESS

GFM_DISCRETE_ENRICHMENT


#-------------------------------------- On-the-fly shock finder


###SHOCK_FINDER_BEFORE_OUTPUT                 #Use this flag if you want to run the shock finder before a snapshot dump, no additional flags or parameters needed.
#SHOCK_FINDER_ON_THE_FLY                    #Run the shock finder at every local timestep, no additional flags or parameters needed.




#--------------------------------------- Post-processing shock finder, please read the instructions in shock_finder.h
#SHOCK_FINDER_POST_PROCESSING               #post-processing shock finder  
#SHOCK_FINDER_AREPO                                    #standard operating mode
#UNLIMITED_GRADIENTS                            #standard option    
#ZONE_JUMP_P                                        #standard option    
#ZONE_JUMP_T                                #standard option
#SHOCK_DIR_GRAD_T                           #standard option
#SHOCK_JUMP_T                               #standard option
#SURFACE_SPHERE_APPROX                      #use this for 2d sims
#SURFACE_ANGLE_APPROX                       #use this for 3d sims
#RESET_WRONG_JUMPS                          #standard option
#SKIP_BORDER                                            #for non-periodic boundaries of the snapshot/subbox
