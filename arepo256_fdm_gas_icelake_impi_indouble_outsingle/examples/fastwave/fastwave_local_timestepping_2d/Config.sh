# examples/fastwave_local_timestepping_2d
NTYPES=6
NSOFTTYPES=2
TWODIMS
VORONOI
DOUBLEPRECISION=1
GAMMA=1.6666666666666667
REGULARIZE_MESH_CM_DRIFT
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED
REGULARIZE_MESH_FACE_ANGLE

LONG_Y=0.2

# generic
VORONOI_DYNAMIC_UPDATE
HAVE_HDF5
DEBUG

MHD
# FORCE_EQUAL_TIMESTEPS
VORONOI_MESH_KEEP_DT_AND_DTC
RIEMANN_HLLD
#MHD_DIVBCLEANING

MESHRELAX_DENSITY_IN_INPUT

OUTPUT_IN_DOUBLEPRECISION                # snapshot files will be written in double precision
INPUT_IN_DOUBLEPRECISION                 # initial conditions are in double precision

# Below this can be removed if we do not use braginskii
BRAGINSKII_VISCOSITY
# BRAGINSKII_VISCOSITY_GLOBAL_TIMESTEP
TETRA_INDEX_IN_FACE
TREE_BASED_TIMESTEPS     # non-local timestep criterion (take 'signal speed' into account)

# BRAGINSKII_RKL2_SUPER_TIME_STEPPING
# VORONOI_STATIC_MESH

# BRAGINSKII_VISCOSITY_CONSTANT_TIMESTEP
# DEBUG
# NOHYDRO