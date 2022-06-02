# AREPO Makefile
#   see documentation/getting_started.md
#
# Add source and header files for new modules in makefiles/modules.make.
# Add libraries for new modules below.

ifneq (FLD_HYPRE, $(findstring FLD_HYPRE, $(CONFIGVARS)))
ifneq (COSMIC_RAYS_DIFFUSION, $(findstring COSMIC_RAYS_DIFFUSION, $(CONFIGVARS)))
ifneq (COSMIC_RAYS_STREAMING, $(findstring COSMIC_RAYS_STREAMING, $(CONFIGVARS)))
ifneq (IMPLICIT_TI, $(findstring IMPLICIT_TI, $(CONFIGVARS)))
ifneq (IMPLICIT_OHMIC_DIFFUSION, $(findstring IMPLICIT_OHMIC_DIFFUSION, $(CONFIGVARS)))
ifneq (SOLAR_RADIATIVE_TRANSFER_DIFF, $(findstring SOLAR_RADIATIVE_TRANSFER_DIFF, $(CONFIGVARS)))
ifneq (SOLAR_RADIATIVE_TRANSFER_EDD, $(findstring SOLAR_RADIATIVE_TRANSFER_EDD, $(CONFIGVARS)))
  HYPRE_INCL =
  HYPRE_LIB =
endif
endif
endif
endif
endif
endif
endif

ifneq (VTUNE_INSTRUMENT, $(findstring VTUNE_INSTRUMENT, $(CONFIGVARS)))
  VTUNE_INCL =
  VTUNE_LIBS =
endif

ifneq (CHIMES, $(findstring CHIMES, $(CONFIGVARS)))
  CHIMESINCL =
  CHIMESLIBS =
endif

ifeq (NETWORK_PARDISO, $(findstring NETWORK_PARDISO, $(CONFIGVARS)))
  MKL_INCL =  -I/cm/shared/apps/intel/composer_xe/2011_sp1.12.361/mkl/include
  MKL_LIBS =  -L/cm/shared/apps/intel/composer_xe/2011_sp1.12.361/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
else
  MKL_INCL =
  MKL_LIBS =
endif

ifeq (EOS_OPAL, $(findstring EOS_OPAL, $(CONFIGVARS)))
  GSL_LIBS += -lz
endif

ifeq (MCMA, $(findstring MCMA, $(CONFIGVARS)))
  LAPACK_LIB = -llapack
endif

ifeq (GRACKLE, $(findstring GRACKLE, $(CONFIGVARS)))
  GRACKLE_INCL = -I$(HOME)/local/include
  GRACKLE_LIBS = -L$(HOME)/local/lib -lgrackle -Xlinker -rpath=$(HOME)/local/lib
else
  GRACKLE_INCL =
  GRACKLE_LIBS =
endif

ifeq (GALPOT, $(findstring GALPOT, $(CONFIGVARS)))
  EIGEN_INCL =
  CPPC = mpic++ -std=c++11
  LINKER   = $(CPPC) -lgfortran
endif

ifeq (SX_CHEMISTRY 4, $(findstring SX_CHEMISTRY 4, $(CONFIGVARS)))
  LINKER   = $(FC)
endif

ifeq (SGCHEM, $(findstring SGCHEM, $(CONFIGVARS)))
  LINKER   = $(FC)
endif

ifeq (TGCHEM, $(findstring TGCHEM, $(CONFIGVARS)))
  CVODE_INCL = -I$(CVODE_DIR)/include
  CVODE_LIB  = -L$(CVODE_DIR)/lib -lsundials_cvode -lsundials_nvecserial -Xlinker -R -Xlinker $(CVODE_DIR)/lib
endif

ifeq (HEALRAY, $(findstring HEALRAY, $(CONFIGVARS)))
  HEALPIX_INCL = -I$(HEALPIX_DIR)/include
  HEALPIX_LIB  = -L$(HEALPIX_DIR)/lib -lchealpix
endif

