# AREPO Makefile
#   see documentation/getting_started.md
#
# Add new systems in makefiles/systypes.make.

AREPO_ROOT      := $(dir $(lastword $(MAKEFILE_LIST)))
AREPO_ROOT      := $(AREPO_ROOT:/=)

EXEC            := Arepo
LIBRARY         := arepo
CONFIG          := $(AREPO_ROOT)/Config.sh
BUILD_DIR       := $(AREPO_ROOT)/build
SRC_DIR         := $(AREPO_ROOT)/src
SUBMAKEFILE_DIR := $(AREPO_ROOT)/makefiles


#####################
# determine SYSTYPE #
#####################
ifdef SYSTYPE
  SYSTYPE := "$(SYSTYPE)"
  -include $(AREPO_ROOT)/Makefile.systype
else
  include $(AREPO_ROOT)/Makefile.systype
endif

MAKEFILES := $(AREPO_ROOT)/Makefile $(SUBMAKEFILE_DIR)/config.make \
  $(SUBMAKEFILE_DIR)/systypes.make $(SUBMAKEFILE_DIR)/modules.make \
  $(SUBMAKEFILE_DIR)/modules-lib.make
ifneq ($(wildcard $(AREPO_ROOT)/Makefile.systype), )
  MAKEFILES += $(AREPO_ROOT)/Makefile.systype
endif

$(info Build configuration:)
$(info SYSTYPE: $(SYSTYPE))
$(info CONFIG: $(CONFIG))
$(info EXEC: $(EXEC))
$(info )

AR     := ar
PYTHON := python

RESULT     := $(shell CONFIG=$(CONFIG) PYTHON=$(PYTHON) AREPO_ROOT=$(AREPO_ROOT) BUILD_DIR=$(BUILD_DIR) make -f $(SUBMAKEFILE_DIR)/config.make)
CONFIGVARS := $(shell cat $(BUILD_DIR)/arepoconfig.h)
RESULT     := $(shell SRC_DIR=$(SRC_DIR) BUILD_DIR=$(BUILD_DIR) $(AREPO_ROOT)/git_version.sh)

MPICHLIB  = -lmpich
GMPLIB    = -lgmp
GSLLIB    = -lgsl -lgslcblas
MATHLIB   = -lm -lstdc++
HWLOC_LIB = -lhwloc


# Add new systems in makefiles/systypes.make.
# If you add a new system, also add that SYSTYPE to Template-Makefile.systype.
include $(SUBMAKEFILE_DIR)/systypes.make

ifeq ($(SYSTYPE),"dirac")
CC       =  mpiicc -g -Wall -Wno-unknown-pragmas -std=c99 -DH5Dopen_vers=1 -DH5Dcreate_vers=1 -DH5Gopen_vers=1 -DH5Gcreate_vers=1 -DH5Acreate_vers=1
OPTIMIZE = -O3
GSL_INCL =
GSL_LIBS = -lgsl -lgslcblas
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
HDF5INCL = -DH5_USE_16_API
HDF5LIB  = -lhdf5 -lm -lz
LIBS = $(GSL_LIBS) $(FFTW_LIBS)
endif

ifndef LINKER
  LINKER := $(CC)
endif


############################################
# determine the needed object/header files #
############################################

SUBDIRS = . mpi_utils pm gitversion subfind

OBJS = \
  powerspec_vel.o                                                            \
  forcetree.o forcetree_walk.o forcetree_ewald.o                             \
  ngbtree.o ngbtree_walk.o forcetree_optimizebalance.o                       \
  domain.o domain_balance.o domain_box.o domain_counttogo.o                  \
  domain_DC_update.o domain_exchange.o domain_toplevel.o                     \
  domain_rearrange.o domain_sort_kernels.o domain_vars.o domain_checks.o     \
  voronoi_dynamic_update.o voronoi_ghost_search.o logs.o                     \
  timestep_treebased.o                                                       \
  external_disk.o riemann.o riemann_rosunov.o riemann_hll.o                  \
  riemann_hllc.o riemann_hlld.o riemann_gamma.o finite_volume_solver.o       \
  set_vertex_velocities.o do_gravity_hydro.o voronoi.o voronoi_utils.o       \
  voronoi_3d.o voronoi_1d.o voronoi_1d_spherical.o                           \
  voronoi_exchange.o voronoi_check.o                                         \
  voronoi_makeimage.o voronoi_makeimage_new.o coffee.o                       \
  windtunnel.o growing_disk_potential.o voronoi_3d_test.o                    \
  allvars.o gradients.o scalars.o subfind/subfind_density.o hdf5_util.o      \
  pinning.o mpz_extension.o run.o predict.o begrun.o                         \
  global.o pm/pm_periodic2d.o timestep.o init.o restart.o io.o io_fields.o   \
  starformation.o mpi_utils/checksummed_sendrecv.o accel.o read_ic.o         \
  parallel_sort.o system.o allocate.o                         \
  mpi_utils/myalltoall.o density.o noh.o mpi_utils/sizelimited_sendrecv.o    \
  mpi_utils/hypercube_allgatherv.o                                           \
  gravtree.o gravdirect.o grav_softening.o grav_external.o                   \
  gravtree_forcetest.o driftfac.o darkenergy.o                               \
  peano.o pm/pm_periodic.o pm/pm_mpi_fft.o pm/pm_nonperiodic.o longrange.o   \
  mymalloc.o helm_eos.o mpi_utils/healthtest.o                               \
  second_derivatives.o                                                       \
  diffusion_fluxes.o extrapolate_quantities.o boundaries.o                   \
  tracer_particle.o tracer_mc.o tracer_trajectory.o debug.o                  \
  update_primitive_variables.o voronoi_derefinement.o                        \
  voronoi_refinement.o voronoi_gradients.o                                   \
  refinement.o criterion_refinement.o parallel_logs.o                        \
  criterion_derefinement.o mpi_utils/mpi_util.o                              \
  relaxobject.o voronoi_proj.o voronoi_derefinement_pairs.o                  \
  sfr_eEOS.o sfr_quicklyalpha.o parameters.o voronoi_gradients_onedims.o     \
  twopoint.o diffusion_general.o runge_kutta_full.o inspiral.o               \
  ngbtree_search.o voronoi_gradients_lsf.o main.o

INCL = \
  allvars.h proto.h forcetree.h domain.h dtypes.h hdf5_util.h                \
  voronoi.h mesh.h helm_eos.h voronoi_proj.h gitversion/version.h            \
  timer.h timestep.h generic_comm_helpers2.h generic_comm_helpers_async.h    \
  runge_kutta_full.h parallel_logs.h coord_util.h subfind/subfind.h

OBJS    += debug_md5/Md5.o debug_md5/calc_checksum.o
INCL    += debug_md5/Md5.h
SUBDIRS += debug_md5


# Add source and header files for new modules in makefiles/modules.make.
# Add libraries for new modules in makefiles/modules-lib.make (see below).
include $(SUBMAKEFILE_DIR)/modules.make


##################################
# determine the needed libraries #
##################################

# we only need FFTW if PMGRID is turned on
ifeq (PMGRID, $(findstring PMGRID, $(CONFIGVARS)))
  # test for double-precision libraries
  ifeq (DOUBLEPRECISION_FFTW, $(findstring DOUBLEPRECISION_FFTW, $(CONFIGVARS)))
    FFTW_LIB = $(FFTW_LIBS) -lfftw3
  else
    FFTW_LIB = $(FFTW_LIBS) -lfftw3f
  endif
# or if POWERSPEC_GRID is activated
else ifeq (POWERSPEC_GRID, $(findstring POWERSPEC_GRID, $(CONFIGVARS)))
  # test for double-precision libraries
  ifeq (DOUBLEPRECISION_FFTW, $(findstring DOUBLEPRECISION_FFTW, $(CONFIGVARS)))
    FFTW_LIB = $(FFTW_LIBS) -lfftw3
  else
    FFTW_LIB = $(FFTW_LIBS) -lfftw3f
  endif
else
  FFTW_LIB =
endif

ifneq (HAVE_HDF5, $(findstring HAVE_HDF5, $(CONFIGVARS)))
  HDF5LIB  =
endif

ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
  THREAD_LIB =
endif

ifneq (CUDA, $(findstring CUDA, $(CONFIGVARS)))
  CUDA_INCL =
  CUDA_LIBS =
endif

ifneq (IMPOSE_PINNING, $(findstring IMPOSE_PINNING, $(CONFIGVARS)))
  HWLOC_INCL =
  HWLOC_LIB =
endif


# Add source and header files for new modules in makefiles/modules.make (see above).
# Add libraries for new modules in makefiles/modules-lib.make.
include $(SUBMAKEFILE_DIR)/modules-lib.make


############################
# combine compiler options #
############################

CFLAGS = \
  $(OPTIMIZE) $(OPT) $(HDF5INCL) $(GSL_INCL) $(EIGEN_INCL) $(FFTW_INCL)      \
  $(CVODE_INCL) $(CFITSIO_INCL) $(HEALPIX_INCL) $(GMP_INCL) $(MKL_INCL)      \
  $(CUDA_INCL) $(HWLOC_INCL) -I$(BUILD_DIR) $(NBC_INCL) $(HYPRE_INCL)        \
  $(VTUNE_INCL) $(GRACKLE_INCL) $(CHIMESINCL)

CFLAGS_CUDA = \
  $(CUDA_OPTIMIZE) $(OPT) $(GSL_INCL) $(FFTW_INCL) $(HDF5INCL) $(CVODE_INCL) \
  $(CFITSIO_INCL) $(HEALPIX_INCL) $(GMP_INCL) $(MKL_INCL) $(CUDA_INCL)       \
  $(VTUNE_INCL) $(GRACKLE_INCL) $(EIGEN_INCL) -I$(BUILD_DIR)

LIBS = \
  $(MATHLIB) $(HDF5LIB) $(MPICHLIB) $(GSL_LIBS) $(GSLLIB) $(FFTW_LIB)        \
  $(GMP_LIBS) $(GMPLIB) $(CVODE_LIB) $(CFITSIO_LIB) $(HEALPIX_LIB)           \
  $(MKL_LIBS) $(THREAD_LIB) $(CUDA_LIBS) $(HWLOC_LIB) $(NBC_LIB)             \
  $(HYPRE_LIB) $(VTUNE_LIBS) $(GRACKLE_LIBS) $(LAPACK_LIB) $(CHIMESLIBS)

FOPTIONS = $(OPTIMIZE)
FFLAGS = -I$(BUILD_DIR) $(FOPTIONS)


SUBDIRS := $(addprefix $(BUILD_DIR)/, $(SUBDIRS))
OBJS := $(addprefix $(BUILD_DIR)/, $(OBJS)) \
  $(BUILD_DIR)/compile_time_info.o $(BUILD_DIR)/compile_time_info_hdf5.o \
  $(BUILD_DIR)/version.o
INCL := $(addprefix $(SRC_DIR)/, $(INCL)) $(BUILD_DIR)/arepoconfig.h


##################
# create subdirs #
##################
RESULT := $(shell mkdir -p $(SUBDIRS))


#############################################
# create info file for command line options #
#############################################
RESULT := $(shell echo 'static const char *compiler_flags = "$(CC) $(CFLAGS)";' > $(BUILD_DIR)/compiler-command-line-args.h)


###############
# build rules #
###############

all: check build

build: $(EXEC)

clean:
	@echo 'Cleaning all build files...'
	@rm -f $(OBJS) $(EXEC) lib$(LIBRARY).a
	@rm -rf $(BUILD_DIR)/
	@rm -rf doxygen/

build_with_config_check: check_config build

$(EXEC): $(OBJS)
	$(LINKER) $(OPTIMIZE) $(OPT) $(LINKER_OPT) $(OBJS) $(LIBS) -o $@

lib$(LIBRARY).a: $(filter-out $(BUILD_DIR)/main.o, $(OBJS))
	$(AR) -rcs lib$(LIBRARY).a $(OBJS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.F   $(INCL) $(MAKEFILES)
	$(FC) -c $(FFLAGS) $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.F90 $(INCL) $(MAKEFILES)
	$(FC) -c $(FFLAGS) $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90 $(INCL) $(MAKEFILES)
	$(FC) -c $(FFLAGS) $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c   $(INCL) $(MAKEFILES)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cc  $(INCL) $(MAKEFILES)
	$(CPPC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cu  $(INCL) $(MAKEFILES)
	$(NVCC) $(CFLAGS_CUDA) -c $< -o $@

$(BUILD_DIR)/%.o: $(BUILD_DIR)/%.c $(MAKEFILES)
	$(CC) $(CFLAGS) -c $< -o $@


#################
# sanity checks #
#################

# find all source files to check for unused options in Template-Config.sh and
# defines_extra
rwildcard = $(foreach d, $(wildcard $(1:/=)/*), $(call rwildcard, $d, $2) $(filter $(subst *, %, $2), $d))
ALL_SOURCE_FILES := $(call rwildcard, $(SRC_DIR), *.h *.c *.cc *.F *.f90 *.F90)
ALL_SOURCE_FILES_CHECK := $(addsuffix .check, $(patsubst $(SRC_DIR)%, $(BUILD_DIR)%, $(ALL_SOURCE_FILES)))
ALL_SOURCE_FILES_CHECK := $(ALL_SOURCE_FILES_CHECK:.c.check=.o.check)
ALL_SOURCE_FILES_CHECK := $(ALL_SOURCE_FILES_CHECK:.cc.check=.o.check)
ALL_SOURCE_FILES_CHECK := $(ALL_SOURCE_FILES_CHECK:.F.check=.o.check)
ALL_SOURCE_FILES_CHECK := $(ALL_SOURCE_FILES_CHECK:.f90.check=.o.check)
ALL_SOURCE_FILES_CHECK := $(ALL_SOURCE_FILES_CHECK:.F90.check=.o.check)
# create directories in $(BUILD_DIR) to contain the .check files
RESULT := $(shell mkdir -p $(dir $(ALL_SOURCE_FILES_CHECK)))

# collect files relevant to current $(CONFIG) to check (given by $(OBJS),
# $(INCL) and $(MAKEFILES))
SOURCE_FILES_TO_CHECK := $(addsuffix .check, $(OBJS) $(patsubst $(SRC_DIR)%, $(BUILD_DIR)%, $(INCL)))
MAKEFILES_TO_CHECK    := $(addsuffix .check, $(addprefix $(BUILD_DIR)/, $(notdir $(MAKEFILES))))
TO_CHECK := $(SOURCE_FILES_TO_CHECK) $(MAKEFILES_TO_CHECK)

# .check file paths for Config.sh, Template-Config.sh, defines_extra, and
# documentation
CONFIG_CHECK = $(BUILD_DIR)/$(notdir $(CONFIG)).check

TEMPLATE_CHECK = $(BUILD_DIR)/Template-Config.sh.check $(BUILD_DIR)/defines_extra.check $(BUILD_DIR)/Template-Config-duplicates.check

DOCS_CHECK = $(BUILD_DIR)/documentation.check


# rules for checks

check: check_config check_template

check_config: $(CONFIG_CHECK)

check_template: $(TEMPLATE_CHECK)

check_docs: $(DOCS_CHECK)

$(CONFIG_CHECK): $(CONFIG) $(TO_CHECK) $(AREPO_ROOT)/check.py
	@$(PYTHON) $(AREPO_ROOT)/check.py config   $< $@ $(TO_CHECK)

$(BUILD_DIR)/%.o.check: $(SRC_DIR)/%.c $(AREPO_ROOT)/Template-Config.sh $(AREPO_ROOT)/defines_extra $(AREPO_ROOT)/check.py
	@$(PYTHON) $(AREPO_ROOT)/check.py code     $< $@ $(AREPO_ROOT)/Template-Config.sh $(AREPO_ROOT)/defines_extra

$(BUILD_DIR)/%.o.check: $(SRC_DIR)/%.cc $(AREPO_ROOT)/Template-Config.sh $(AREPO_ROOT)/defines_extra $(AREPO_ROOT)/check.py
	@$(PYTHON) $(AREPO_ROOT)/check.py code     $< $@ $(AREPO_ROOT)/Template-Config.sh $(AREPO_ROOT)/defines_extra

$(BUILD_DIR)/%.o.check: $(SRC_DIR)/%.F $(AREPO_ROOT)/Template-Config.sh $(AREPO_ROOT)/defines_extra $(AREPO_ROOT)/check.py
	@$(PYTHON) $(AREPO_ROOT)/check.py code     $< $@ $(AREPO_ROOT)/Template-Config.sh $(AREPO_ROOT)/defines_extra

$(BUILD_DIR)/%.o.check: $(SRC_DIR)/%.f90 $(AREPO_ROOT)/Template-Config.sh $(AREPO_ROOT)/defines_extra $(AREPO_ROOT)/check.py
	@$(PYTHON) $(AREPO_ROOT)/check.py code     $< $@ $(AREPO_ROOT)/Template-Config.sh $(AREPO_ROOT)/defines_extra

$(BUILD_DIR)/%.o.check: $(SRC_DIR)/%.F90 $(AREPO_ROOT)/Template-Config.sh $(AREPO_ROOT)/defines_extra $(AREPO_ROOT)/check.py
	@$(PYTHON) $(AREPO_ROOT)/check.py code     $< $@ $(AREPO_ROOT)/Template-Config.sh $(AREPO_ROOT)/defines_extra

$(BUILD_DIR)/%.h.check: $(SRC_DIR)/%.h $(AREPO_ROOT)/Template-Config.sh $(AREPO_ROOT)/defines_extra $(AREPO_ROOT)/check.py
	@$(PYTHON) $(AREPO_ROOT)/check.py code     $< $@ $(AREPO_ROOT)/Template-Config.sh $(AREPO_ROOT)/defines_extra

$(BUILD_DIR)/%.o.check: $(BUILD_DIR)/%.c $(AREPO_ROOT)/Template-Config.sh $(AREPO_ROOT)/defines_extra $(AREPO_ROOT)/check.py
	@$(PYTHON) $(AREPO_ROOT)/check.py code     $< $@ $(AREPO_ROOT)/Template-Config.sh $(AREPO_ROOT)/defines_extra

$(BUILD_DIR)/%.h.check: $(BUILD_DIR)/%.h $(AREPO_ROOT)/Template-Config.sh $(AREPO_ROOT)/defines_extra $(AREPO_ROOT)/check.py
	@$(PYTHON) $(AREPO_ROOT)/check.py code     $< $@ $(AREPO_ROOT)/Template-Config.sh $(AREPO_ROOT)/defines_extra

$(BUILD_DIR)/Makefile.check $(BUILD_DIR)/Makefile.systype.check: $(AREPO_ROOT)/Makefile $(AREPO_ROOT)/Template-Config.sh $(AREPO_ROOT)/defines_extra $(AREPO_ROOT)/check.py
	@$(PYTHON) $(AREPO_ROOT)/check.py makefile $< $@ $(AREPO_ROOT)/Template-Config.sh $(AREPO_ROOT)/defines_extra

$(BUILD_DIR)/%.make.check: $(SUBMAKEFILE_DIR)/%.make $(AREPO_ROOT)/Template-Config.sh $(AREPO_ROOT)/defines_extra $(AREPO_ROOT)/check.py
	@$(PYTHON) $(AREPO_ROOT)/check.py makefile $< $@ $(AREPO_ROOT)/Template-Config.sh $(AREPO_ROOT)/defines_extra

$(BUILD_DIR)/Template-Config.sh.check: $(AREPO_ROOT)/Template-Config.sh $(AREPO_ROOT)/check.py $(ALL_SOURCE_FILES_CHECK) $(TO_CHECK)
$(BUILD_DIR)/defines_extra.check: $(AREPO_ROOT)/defines_extra $(AREPO_ROOT)/check.py $(ALL_SOURCE_FILES_CHECK) $(TO_CHECK)
$(BUILD_DIR)/Template-Config.sh.check $(BUILD_DIR)/defines_extra.check:
	@$(PYTHON) $(AREPO_ROOT)/check.py template-unused $< $@ $(ALL_SOURCE_FILES_CHECK) $(MAKEFILES_TO_CHECK)

$(BUILD_DIR)/Template-Config-duplicates.check: $(BUILD_DIR)/Template-Config.sh.check $(BUILD_DIR)/defines_extra.check $(AREPO_ROOT)/check.py
	@$(PYTHON) $(AREPO_ROOT)/check.py template-duplicate $(TEMPLATE_CHECK)

$(DOCS_CHECK): $(AREPO_ROOT)/Template-Config.sh $(AREPO_ROOT)/check.py
	@$(PYTHON) $(AREPO_ROOT)/check.py documentation $< $@

