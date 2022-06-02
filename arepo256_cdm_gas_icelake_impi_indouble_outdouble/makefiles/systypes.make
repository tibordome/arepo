# AREPO Makefile
#   see documentation/getting_started.md
#
# Add new systems here.
# If you add a new system below, also add that SYSTYPE to Template-Makefile.systype.


ifeq ($(SYSTYPE), "ARCCA")
  MPICHLIB =
  CC       =  mpicc -std=c99
  FC       =  mpif90 -nofor_main
  OPTIMIZE =  -O3
  GSL_INCL =  -I/software/libraries/gsl/1.15/gnu-4.1.2/include/
  GSL_LIBS =  -L/software/libraries/gsl/1.15/gnu-4.1.2/lib/ -lgsl -lgslcblas -lm
  GMP_INCL =  -I/software/libraries/gmp/5.0.2/include/
  GMP_LIBS =  -L/software/libraries/gmp/5.0.2/lib/
  FFTW_INCL=
  FFTW_LIBS=
  HDF5INCL =
  HDF5LIB  =
  HWLOC_INCL= -I/usr/include/
  HWLOC_LIB = -L/usr/local/lib64 -lhwloc -Xlinker -R -Xlinker /usr/local/lib64
  LINKER   = $(FC)
endif


ifeq ($(SYSTYPE), "aurora")
  CC       =  mpicc   # sets the C compiler
  OPTIMIZE =  -std=c11 -g -O3 -Wall -ipo -vec_report0
  GSL_INCL =
  GSL_LIBS =
  FFTW_INCL=
  FFTW_LIBS= -L/home/spb/.local/fftw/lib
  MPICHLIB = -lmpi
  HDF5INCL = -DH5_USE_16_API
  HDF5LIB  = -lhdf5
  AR = xiar
endif


ifeq ($(SYSTYPE), "binac")
  CC         = mpicc
  OPTIMIZE   = -std=c11 -O3 -g -Wall
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -openmp
  else
    OPTIMIZE += -Wno-unknown-pragmas -Wno-unused-function
  endif
  GSL_INCL   = -I$(GSL_INC_DIR)
  GSL_LIBS   = -L$(GSL_LIB_DIR)
  GMP_INCL   = -I$(GMP_DIR)/include
  GMP_LIBS   = -L$(GMP_DIR)/lib
  FFTW_INCL  = -I$(FFTW_INC_DIR)
  FFTW_LIBS  = -L$(FFTW_LIB_DIR)
  MPICHLIB   =
endif

ifeq ($(SYSTYPE), "BInAC")
  # bwfor BInAC in Tuebingen
  CC        = mpicc
  FC        = mpif90
  OPTIMIZE  = -std=c11 -O3 -g -xhost
  GSL_INCL  = -I$(GSL_INC_DIR)
  GSL_LIBS  = -L$(GSL_LIB_DIR)
  FFTW_INCL = -I$(FFTW_INC_DIR)
  FFTW_LIBS = -L$(FFTW_LIB_DIR)
  MKL_INCL  = -I$(MKL_LIB_DIR)
  MKL_LIBS  = -L$(MKL_INC_DIR)
  HDF5INCL  = $(HDF5_INC) -DH5_USE_16_API=1
  HDF5LIB   = $(HDF5_LIB) -lhdf5
  HWLOC_INCL= -I/usr/include
  HWLOC_LIBS= -L/usr/lib64
  GMP_INCL  = -I/usr/include
  GMP_LIBS  = -L/usr/lib64
  MPICHLIB  = -lmpi
  LINKER    = $(FC) -nofor-main
endif


ifeq ($(SYSTYPE), "bocksbeutel")
  CC       =  mpicc   # sets the C compiler
  OPTIMIZE =  -std=c11 -ggdb -O3 -funroll-all-loops -fprefetch-loop-arrays --param prefetch-latency=300 -ftree-vectorize -march=bdver1 -mprefer-avx128 -mcx16 -msahf -mno-movbe -maes -mpclmul -mpopcnt -mabm -mlwp -mno-fma -mfma4 -mxop -mno-bmi -mno-bmi2 -mno-tbm -mavx -mno-avx2 -msse4.2 -msse4.1 -mlzcnt -mno-rdrnd -mno-f16c -mno-fsgsbase --param l1-cache-size=16 --param l1-cache-line-size=64 --param l2-cache-size=2048 -mtune=bdver1
  GSL_INCL =
  GSL_LIBS =
  FFTW_INCL=
  FFTW_LIBS=
  MPICHLIB = -lmpi -lhwloc
  HDF5INCL = -DH5_USE_16_API=1
  HDF5LIB  = -lhdf5
endif


ifeq ($(SYSTYPE), "bwFor")
  CC        = mpicc
  FC        = mpif90
  OPTIMIZE  = -std=c11 -O3 -g -xhost
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -fopenmp
  else
    OPTIMIZE += -Wno-unknown-pragmas -Wno-unused-function
  endif
  GSL_INCL  = -I$(GSL_INC_DIR)
  GSL_LIBS  = -L$(GSL_LIB_DIR)
  GMP_INCL  = -I$(GMPDIR)/include
  GMP_LIBS  = -L$(GMPDIR)/lib -Xlinker -R -Xlinker $(GMPDIR)/lib
  FFTW_INCL = -I$(FFTW_INC_DIR)
  FFTW_LIBS = -L$(FFTW_LIB_DIR)
  MKL_INCL  = -I$(MKL_LIB_DIR)
  MKL_LIBS  = -L$(MKL_INC_DIR)
  HDF5INCL  = -I$(HDF5_INC_DIR) -DH5_USE_16_API
  HDF5LIB   = -L$(HDF5_LIB_DIR) -lhdf5 -lz
  MPICHLIB  = -lmpi
  LINKER    = $(FC) -nofor-main
endif


ifeq ($(SYSTYPE), "bwforcluster")
  CC         = mpicc
  OPTIMIZE   = -std=c11 -O3 -g -Wall
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE   += -openmp
  else
    OPTIMIZE   += -Wno-unknown-pragmas -Wno-unused-function
  endif
  GSL_INCL   = -I$(GSL_INC_DIR)
  GSL_LIBS   = -L$(GSL_LIB_DIR)
  GMP_INCL   = -I$(GMP_DIR)/include
  GMP_LIBS   = -L$(GMP_DIR)/lib
  FFTW_INCL  = -I$(FFTW_INC_DIR)
  FFTW_LIBS  = -L$(FFTW_LIB_DIR)
  MKL_INCL   = -I$(MKL_LIB_DIR)
  MKL_LIBS   = -L$(MKL_INC_DIR)
  HDF5INCL   = -I$(HDF5_INC_DIR) -DH5_USE_16_API
  HDF5LIB    = -L$(HDF5_LIB_DIR) -lhdf5 -lz
  HWLOC_INCL =  -I/home/hd/hd_hd/hd_mm002/libraries/include
  HWLOC_LIBS =  -lhwloc -L/home/hd/hd_hd/hd_mm002/libraries/lib -Xlinker -R -Xlinker /home/hd/hd_hd/hd_mm002/libraries/lib
  MPICHLIB   = -lmpi
endif


ifeq ($(SYSTYPE), "bwForDev")
  CC        = mpicc
  FC        = mpif90
  OPTIMIZE  = -std=c11 -O3 -g # -xhost
  GSL_INCL  = -I$(USER_INC_DIR)
  GSL_LIBS  = -L$(USER_LIB_DIR)
  GMP_INCL  =
  GMP_LIBS  =
  FFTW_INCL =
  FFTW_LIBS =
  MKL_INCL  =
  MKL_LIBS  =
  HWLOC_INCL=
  HWLOC_LIB = -lhwloc
  HDF5INCL  = -DH5_USE_16_API
  HDF5LIB   = -lhdf5 # -lz
  MPICHLIB  = -lmpi
  LINKER    = $(FC) # -nofor-main
endif


ifeq ($(SYSTYPE), "bwfordev")
  CC         =  mpicc     # sets the C compiler
  OPTIMIZE   =  -std=c99 -O3 -g -Wall -Wno-unknown-pragmas -Wno-unused-function
  GMP_INCL   =  -I$(GMP_DIR)/include
  GMP_LIBS   =  -L$(GMP_DIR)/lib
  HDF5INCL   =
  HDF5LIB    =
  HWLOC_INCL =
  HWLOC_LIBS =
  MPICHLIB   =  -lmpi
endif


ifeq ($(SYSTYPE), "Centos5-Gnu")
  CC       =  mpicc   # sets the C compiler
  OPTIMIZE =  -std=c11 -ggdb -O3
  GSL_INCL =
  GSL_LIBS =
  FFTW_INCL= -L/usr/local/lib
  FFTW_LIBS=
  MPICHLIB = -lmpi
  HDF5INCL =
  HDF5LIB  = -lhdf5
endif


ifeq ($(SYSTYPE), "Centos5-Intel")
  CC       =  icc   # sets the C compiler
  OPTIMIZE =  -std=c11 -g -O3 -ipo  -mtune=host -mcpu=host -mp
  GSL_INCL =
  GSL_LIBS =
  FFTW_INCL=
  FFTW_LIBS=
  MPICHLIB = -lmpi
  HDF5INCL =
  HDF5LIB  =
endif


ifeq ($(SYSTYPE), "Cobra")
  # alias for SYSTYPE "MPCDF-OpenMPI" (with some extra settings)
  systype_mpcdf_open_mpi = 1
endif


ifeq ($(SYSTYPE), "Cobra-Intel")
  # alias for SYSTYPE "MPCDF-Intel"
  systype_mpcdf_intel = 1
endif


ifeq ($(SYSTYPE), "Cobra-OpenMPI")
  # alias for SYSTYPE "MPCDF-OpenMPI"
  systype_mpcdf_open_mpi = 1
endif


ifeq ($(SYSTYPE), "CobraOpenMPI")
  # HPC system Cobra at MPCDF, Garching:
  #   https://www.mpcdf.mpg.de/services/supercomputing/cobra
  # module load gsl fftw-serial hdf5-serial
  # need to start job with /u/vrs/Libs/openmpi-4.0.1/bin/mpiexec
  CC       = /u/vrs/Libs/openmpi-4.0.1/bin/mpicc -std=c11
  OPT     += -Wall -Wextra -Wno-unused-parameter
  LINKER   = /u/vrs/Libs/openmpi-4.0.1/bin/mpicc -Xlinker -rpath=$(GSL_HOME)/lib:$(FFTW_HOME)/lib:$(HDF5_HOME)/lib
  OPTIMIZE = -O3
  GSL_INCL = -I$(GSL_HOME)/include/
  GSL_LIBS = -L$(GSL_HOME)/lib/
  FFTW_INCL= -I$(FFTW_HOME)/include/
  FFTW_LIBS= -L$(FFTW_HOME)/lib/
  MPICHLIB =
  HDF5INCL = -DH5_USE_16_API=1 -I$(HDF5_HOME)/include/
  HDF5LIB  = -lhdf5 -L$(HDF5_HOME)/lib/
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -fopenmp
  else
    OPTIMIZE += -Wno-unknown-pragmas
  endif
endif


ifeq ($(SYSTYPE), "cosma")
  CC       =  mpicc
  FC       =  mpif90
  OPTIMIZE =  -std=c99 -O2 -g
  GSL_INCL =
  GSL_LIBS =  -lgsl
  FFTW_INCL=
  FFTW_LIBS=  -lfftw
  HDF5INCL =
  HDF5LIB  =  -lhdf5
  MPICHLIB = -lmpi
  HWLOC_INCL= -I/usr/include
  HWLOC_LIB = $(LDFLAGS) -lhwloc
  LINKER   = $(FC)
endif


ifeq ($(SYSTYPE), "cosma5")
  CC       = $(MPIROOT)/bin/mpicc
  OPTIMIZE = -std=c11 -DUSE_IRECV -g -O3 -ip -DH5_USE_16_API -Wno-unknown-pragmas
  HDF5LIB  = -lhdf5
  MPICHLIB = -lmpi
  EXTRA_LDFLAGS = -shared-intel $(OMPI_LDFLAGS) # only for intel compiler
endif


ifeq ($(SYSTYPE), "Curie")
  CC         = mpicc
  OPTIMIZE   = -std=c11 -g -O2 #-xW -ipo -Wall
  GMP_INCL   = -I/ccc/work/cont005/ra0844/springev/libs/gmp/include
  GMP_LIBS   = -L/ccc/work/cont005/ra0844/springev/libs/gmp/lib  -Xlinker -R -Xlinker /ccc/work/cont005/ra0844/springev/libs/gmp/lib
  GSL_INCL   = -I$(GSL_INC_DIR)
  GSL_LIBS   = -L$(GSL_LIB_DIR)
  FFTW_INCL  = -I$(FFTW2_INC_DIR)
  FFTW_LIBS  = -L$(FFTW2_LIB_DIR)
  HDF5INCL   = -I$(HDF5_INC_DIR) -DH5_USE_16_API
  HDF5LIB    = -L$(HDF5_LIB_DIR) -lhdf5 -lz
  HWLOC_INCL = -I/ccc/work/cont005/ra0844/springev/libs/hwloc/include/
  HWLOC_LIB  = -L/ccc/work/cont005/ra0844/springev/libs/hwloc/lib/ -lhwloc -Xlinker -R -Xlinker /ccc/work/cont005/ra0844/springev/libs/hwloc/lib/
  MPICHLIB   =
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -openmp
  else
    OPTIMIZE += -Wno-unknown-pragmas
  endif
endif


ifeq ($(SYSTYPE), "Darwin")
  CC       = mpicxx   # sets the C compiler
  FC       = mpicxx
  LINKER   = $(FC)
  OPTIMIZE = -std=c++11 -O2 -Wall -Wno-format-security -Wno-unknown-pragmas -Wno-unused-function
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -fopenmp
    MPI_COMPILE_FLAGS = $(shell mpicc --showme:compile)
    #CC       =  /opt/local/bin/gcc-mp-7  $(MPI_COMPILE_FLAGS)      # to replace clang with gcc (mpicc uses clang for some reason)
  endif
  GSL_INCL = -I/opt/local/include
  GSL_LIBS = -L/opt/local/lib
  FFTW_INCL= -I/opt/local/include -I/usr/local/include
  FFTW_LIBS= -L/opt/local/lib -I/usr/local/lib
  HDF5INCL = -I/opt/local/include -DH5_USE_16_API
  HDF5LIB  = -L/opt/local/lib -lhdf5 -lz
  HWLOC_INCL= -I/opt/local/include
  HWLOC_LIB = -L/opt/local/lib -lhwloc
  HYPRE_INCL = -I/Users/varan/codes/libs/hypre-2.10.0b/src/hypre/include
  HYPRE_LIB = -L/Users/varan/codes/libs/hypre-2.10.0b/src/hypre/lib -lHYPRE
  MPICHLIB = -lmpi
endif


ifeq ($(SYSTYPE), "Darwin-mpich")
  CC       = mpicc   # sets the C compiler
  LINKER   = mpicc
  OPTIMIZE = -std=c11 -m64 -ggdb -O3 -Wall -Wno-format-security -Wno-unknown-pragmas
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE+= -fopenmp
  endif
  GSL_INCL = -I/opt/local/include
  GSL_LIBS = -L/opt/local/lib
  FFTW_INCL=
  FFTW_LIBS=
  HDF5INCL = -DH5_USE_16_API
  HDF5LIB  = -lhdf5 -lz
  HWLOC_INCL =
  HWLOC_LIB = -lhwloc
  CUDA_INCL= -I/Developer/NVIDIA/CUDA-5.0/include
  CUDA_LIBS= -Xlinker -rpath /Developer/NVIDIA/CUDA-5.0/lib -L/Developer/NVIDIA/CUDA-5.0/lib -lcudart -lnvToolsExt -framework CUDA
  NVCC     = /Developer/NVIDIA/CUDA-5.0/bin/nvcc
  CUDA_OPTIMIZE = -g -G -O3 -m64 --ptxas-options=-v -Xptxas="-v" --maxrregcount=32 -arch=sm_30 $(filter -I%, $(shell mpicc -show))
endif


ifeq ($(SYSTYPE), "darwinUK")
  CC       =  mpiicc -g -Wall -Wno-unknown-pragmas -std=c99 -DH5Dopen_vers=1 -DH5Dcreate_vers=1 -DH5Gopen_vers=1 -DH5Gcreate_vers=1 -DH5Acreate_vers=1
  OPTIMIZE = -O3
  GSL_INCL =
  GSL_LIBS = -lgsl -lgslcblas
  FFTW_INCL=
  FFTW_LIBS= -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
  MPICHLIB =
  HDF5INCL =
  HDF5LIB  = -lhdf5 -lm -lz
  LIBS     = $(GSL_LIBS) $(FFTW_LIBS)
endif


ifeq ($(SYSTYPE), "draco")
  # alias for SYSTYPE "MPCDF-Intel" (with some extra settings)
  systype_mpcdf_intel = 1
endif


ifeq ($(SYSTYPE), "Draco")
  # alias for SYSTYPE "MPCDF-OpenMPI"
  systype_mpcdf_open_mpi = 1
endif


ifeq ($(SYSTYPE), "Draco-Intel")
  # alias for SYSTYPE "MPCDF-Intel"
  systype_mpcdf_intel = 1
endif


ifeq ($(SYSTYPE), "Draco-OpenMPI")
  # alias for SYSTYPE "MPCDF-OpenMPI"
  systype_mpcdf_open_mpi = 1
endif


ifeq ($(SYSTYPE), "DracoAndHydra")
  # HPC system Draco at MPCDF, Garching:
  #   https://www.mpcdf.mpg.de/services/supercomputing/draco
  CC       =  mpiicc -std=c11
  OPTIMIZE =  -march=native -g -ipo -O3 -Wundef
  GSL_INCL =  -I$(GSL_HOME)/include/
  GSL_LIBS =  -L$(GSL_HOME)/lib/
  FFTW_INCL=  -I$(FFTW_HOME)/include/
  FFTW_LIBS=  -L$(FFTW_HOME)/lib/
  HWLOC_INCL=
  HWLOC_LIB = -lhwloc -L/u/dnelson/.local/lib/
  MPICHLIB =
  HDF5INCL =  -DH5_USE_16_API=1 -I$(HDF5_HOME)/include/
  HDF5LIB  =  -lhdf5 -L$(HDF5_HOME)/lib/ -Xlinker -R -Xlinker $(HDF5_HOME)/lib
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -qopenmp
  else
    OPTIMIZE += -Wno-unknown-pragmas
  endif
endif


ifeq ($(SYSTYPE), "elgato")
  CC       =  mpicc   # sets the C compiler
  OPTIMIZE =  -std=c11 -ggdb -O3 -m64 -Wall -Wno-format-security -Wno-unknown-pragmas -Wno-unused-function
  GMP_INCL = -I$(GMP_INC_DIR)
  GMP_LIBS = -L$(GMP_LIB_DIR) # -Xlinker -R -Xlinker $(GMP_LIB_DIR)
  GSL_INCL = -I$(GSL_INC_DIR)
  GSL_LIBS = -L$(GSL_LIB_DIR)
  HDF5INCL = -I$(HDF5_INC_DIR) -DH5_USE_16_API
  HDF5LIB  = -L$(HDF5_LIB_DIR) -lhdf5 -lz #-Xlinker -R -Xlinker $(HDF5_LIB_DIR)
  FFTW_INCL=
  FFTW_LIBS=
  MPICHLIB = -lmpi
endif


ifeq ($(SYSTYPE), "Flatiron")
  CC       =  mpicc     # sets the C compiler
  OPT      += -DMPICH_IGNORE_CXX_SEEK
  OPTIMIZE =  -std=c99 -O3 -g -Wall -Wno-unused-but-set-variable -Wno-uninitialized -Wno-unknown-pragmas -Wno-unused-function -march=native
  GSL_INCL =
  GSL_LIBS =
  FFTW_INCL=
  FFTW_LIBS=
  MPICHLIB =
  HDF5INCL =  -DH5_USE_16_API
  HDF5LIB  =  -lhdf5 -lz
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -fopenmp
    OPT      += -DIMPOSE_PINNING -DSOCKETS=4 -DMAX_CORES=16
  endif
endif


ifeq ($(SYSTYPE), "ForHLRI")
  # module load compiler/gnu/5 mpi/openmpi/1.10 lib/hdf5/1.8
  CC       = mpicc   # sets the C compiler
  OPTIMIZE = -std=c11 -ggdb -O3 -ftree-vectorize -march=native
  GMP_INCL = -I/opt/gmp/6.0/include
  GMP_LIBS = -L/opt/gmp/6.0/lib
  GSL_INCL =
  GSL_LIBS =
  FFTW_INCL=
  FFTW_LIBS=
  MPICHLIB = -lmpi
  HDF5INCL = -DH5_USE_16_API=1
  HDF5LIB  = -lhdf5
endif


ifeq ($(SYSTYPE), "Freya")
  # alias for SYSTYPE "MPCDF-Intel"
  systype_mpcdf_intel = 1
endif


ifeq ($(SYSTYPE), "Freya-Intel")
  # alias for SYSTYPE "MPCDF-Intel"
  systype_mpcdf_intel = 1
endif


ifeq ($(SYSTYPE), "Freya-OpenMPI")
  # alias for SYSTYPE "MPCDF-OpenMPI"
  systype_mpcdf_open_mpi = 1
endif


ifeq ($(SYSTYPE), "FreyaOpenMPI")
  # Linux cluster FREYA for the Max Planck Institute for Astrophysics (MPA)
  # at MPCDF, Garching:
  #   https://docs.mpcdf.mpg.de/doc/computing/clusters/systems/Astrophysics.html
  # module load gcc fftw-serial hdf5-serial gsl
  CC       = /u/vrs/Libs/openmpi-3.1.2/bin/mpicxx -std=c++11 # sets the C compiler
  CPP      = /u/vrs/Libs/openmpi-3.1.2/bin/mpicxx -std=c++11 # sets the C++ compiler
  OPT     += -Wall -Wextra -Wno-unused-parameter
  LINKER   = /u/vrs/Libs/openmpi-3.1.2/bin/mpicxx -Xlinker -rpath=$(GSL_HOME)/lib:$(FFTW_HOME)/lib:$(HDF5_HOME)/lib
  OPTIMIZE = -O2 -g
  GSL_INCL = -I$(GSL_HOME)/include/
  GSL_LIBS = -L$(GSL_HOME)/lib/
  FFTW_INCL= -I$(FFTW_HOME)/include/
  FFTW_LIBS= -L$(FFTW_HOME)/lib/
  MPICHLIB =
  HDF5INCL = -DH5_USE_16_API=1 -I$(HDF5_HOME)/include/
  HDF5LIB  = -lhdf5 -L$(HDF5_HOME)/lib/
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -fopenmp
  else
    OPTIMIZE += -Wno-unknown-pragmas
  endif
endif


ifeq ($(SYSTYPE), "Gordon")
  # module load fftw/2.1.5 hdf5/1.8.11 gsl/1.15
  CC       = mpicc     # sets the C compiler
  OPT      += -DNOCALLSOFSYSTEM  -DMPICH_IGNORE_CXX_SEEK
  OPTIMIZE = -std=c99 -O3 -g -Wall #-mavx #-xHOST
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -openmp
  else
    OPTIMIZE += -Wno-unknown-pragmas
  endif
  GSL_INCL = -I/opt/gsl/gnu/include -I/opt/gnu/gmp/4.3.2/include
  GSL_LIBS = -L/opt/gsl/gnu/lib -L/opt/gnu/gmp/4.3.2/lib
  FFTW_INCL= -I/opt/fftw/2.1.5/gnu/mvapich2/include
  FFTW_LIBS= -L/opt/fftw/2.1.5/gnu/mvapich2/lib
  MPICHLIB =
  HDF5INCL = -I/opt/hdf5/intel/mvapich2/ib/include -DH5_USE_16_API
  HDF5LIB  = -L/opt/hdf5/intel/mvapich2/ib/lib -lhdf5
endif


ifeq ($(SYSTYPE), "Grendel")
    # There's a LOT of non-standard packages in Grendel, so point them to your personal installation
    CC       = mpicxx
    FC       = gcc
    OPTIMIZE = -std=c++11 -O2 -Wall -Wno-format-security -Wno-unknown-pragmas -Wno-unused-function -g
    GSL_INCL = -I/com/gsl/include -I/home/mark/MESA/mesasdk/include
    GSL_LIBS = -L/com/gsl/lib -L/home/mark/MESA/mesasdk/lib
    HDF5INCL = -I/home/mark/cmodules/hdf5/include -DH5_USE_16_API
    HDF5LIB  = -L/home/mark/cmodules/hdf5/lib -lhdf5
    HYPRE_INCL = -I/home/mark/cmodules/hypre-master/src/hypre/include
    HYPRE_LIB  = -L/home/mark/cmodules/hypre-master/src/hypre/lib -lHYPRE
    MPICHLIB = -lmpi
endif


ifeq ($(SYSTYPE), "Haswell")
  CC         = mpicc
  FC         = mpif90
  OPTIMIZE   = -std=c11 -O3 -msse3 -g -Wall -m64
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -fopenmp
  else
    OPTIMIZE += -Wno-unknown-pragmas -Wno-unused-function
  endif
  GMP_INCL   = -I/hits/basement/tap/sw/libs/gmp-6.1.1/include
  GMP_LIBS   = -L/hits/basement/tap/sw/libs/gmp-6.1.1/lib -Xlinker -R -Xlinker /hits/basement/tap/sw/libs/gmp-6.1.1/lib
  GSL_INCL   = -I/hits/basement/tap/sw/libs/gsl-2.1/include
  GSL_LIBS   = -L/hits/basement/tap/sw/libs/gsl-2.1/lib -Xlinker -R -Xlinker /hits/basement/tap/sw/libs/gsl-2.1/lib
  FFTW_INCL  =
  FFTW_LIBS  =
  MPICHLIB   =
  HDF5INCL   = -I/hits/basement/tap/sw/libs/hdf5-1.8.17/include -DH5_USE_16_API
  HDF5LIB    = -L/hits/basement/tap/sw/libs/hdf5-1.8.17/lib -lhdf5 -Xlinker -R -Xlinker /hits/basement/tap/sw/libs/hdf5-1.8.17/lib
  HWLOC_INCL = -I/hits/basement/tap/sw/libs/hwloc-1.11.3/include
  HWLOC_LIB  = -L/hits/basement/tap/sw/libs/hwloc-1.11.3/lib -lhwloc -Xlinker -R -Xlinker /hits/basement/tap/sw/libs/hwloc-1.11.3/lib
  LINKER     = $(CC) -lgfortran #-nofor-main
  HYPRE_LIB  = -lHYPRE
  VTUNE_INCL = -I/cm/shared/apps/intel/vtune_u13/vtune_amplifier_xe_2015/include
  VTUNE_LIBS = -L/cm/shared/apps/intel/vtune_u13/vtune_amplifier_xe_2015/lib64 -littnotify
endif

ifeq ($(SYSTYPE), "Haswell_pso")
  CC         = mpicc
  FC         = mpif90
  OPTIMIZE   = -std=c11 -O3 -msse3 -g -Wall -m64
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -fopenmp
  else
    OPTIMIZE += -Wno-unknown-pragmas -Wno-unused-function
  endif
  LINKER     = $(CC) -lgfortran #-nofor-main
  # IF ERROR "is directory" when using the cluster-provided OpenMPI, use simply instead:
  # MPICHLIB   = --enable-new-dtags -lmpi 
  # Otherwise:
  MPICHLIB   = $(shell pkg-config --libs ompi) 
  GSL_LIBS   = $(shell gsl-config --libs)
  GSL_INCL   = $(shell gsl-config --cflags)
  HWLOC_LIBS = $(shell pkg-config --libs hwloc)
  HWLOC_INCL = $(shell pkg-config --cflags hwloc)
  GMP_LIBS   = $(shell pkg-config --libs gmp)
  GMP_INCL   = $(shell pkg-config --cflags gmp)
  FFTW_INCL  = 
  FFTW_LIBS  = -lfftw
  HDF5LIB    = -lhdf5 -lz
  HDF5INCL   = -DH5_USE_16_API
  HYPRE_LIB  = -lHYPRE
endif

ifeq ($(SYSTYPE), "hecate")
  CC       = mpicc   # sets the Ccompiler
  OPTIMIZE = -std=c11 -g -O3 -w1 -ipo -vec_report0 -xHOST
  GSL_INCL =
  GSL_LIBS =
  FFTW_INCL=
  FFTW_LIBS=
  MPICHLIB = -lmpi
  HDF5INCL = -DH5_USE_16_API
  HDF5LIB  = -lhdf5
endif


ifeq ($(SYSTYPE),"HiPerGator2")
  CC       = mpicc
  CXX      = mpicxx
  FC       = $(CC)
  OPT     += -DH5_USE_16_API
  OPTIMIZE = -std=c11 -O3 -g -Wall -parallel -ipo -funroll-loops -xcore-avx2
  GSL_INCL = -I$(HPC_GSL_INC)
  GSL_LIBS = -L$(HPC_GSL_LIB)
  FFTW_INCL= -I$(HPC_FFTW_INC)
  FFTW_LIBS= -L$(HPC_FFTW_LIB)
  HDF5LIB  = -lhdf5
  HDF5INCL =
  MPICHLIB =
  ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
    OPTIMIZE  +=  -openmp
  else
    OPTIMIZE  +=  -diag-disable 3180
  endif
endif


ifeq ($(SYSTYPE), "Hornet")
  CC       = cc -std=c11 -g -O2  #-xW -ipo -Wall
  OPTIMIZE =
  GMP_INCL = -I/zhome/academic/HLRS/lha/zahrpakm/libs/include
  GMP_LIBS = -L/zhome/academic/HLRS/lha/zahrpakm/libs/lib  -Xlinker -R -Xlinker /zhome/academic/HLRS/lha/zahrpakm/libs/lib
  GSL_INCL = -I/zhome/academic/HLRS/lha/zahrpakm/libs/include
  GSL_LIBS = -L/zhome/academic/HLRS/lha/zahrpakm/libs/lib  -Xlinker -R -Xlinker /zhome/academic/HLRS/lha/zahrpakm/libs/lib
  HDF5INCL = -I$(HDF5_DIR)/include -DH5_USE_16_API # need to do module load cray-hdf5
  HDF5LIB  = -L$(HDF5_LIB)/lib -lhdf5 -lz
  HWLOC_INCL = -I/zhome/academic/HLRS/lha/zahrpakm/libs/include
  HWLOC_LIB  = -L/zhome/academic/HLRS/lha/zahrpakm/libs/lib -lhwloc
  MPICHLIB =
  OPT      +=  -DNOCALLSOFSYSTEM
endif


ifeq ($(SYSTYPE), "Judge")
  CC       =  mpicc
  LINKER   =  mpicxx
  OPTIMIZE =  -std=c11 -g -O3 -Wall
  GSL_INCL =  -I/homeb/zam/baueras/libs/include
  GSL_LIBS =  -L/homeb/zam/baueras/libs/lib  -Xlinker -R -Xlinker /homeb/zam/baueras/libs/lib
  FFTW_INCL=
  FFTW_LIBS=
  MPICHLIB =
  HDF5INCL = -DH5_USE_16_API
  HDF5LIB  = -lhdf5
  CUDA_INCL= -I/opt/cuda/include/
  CUDA_LIBS= -L/opt/cuda/lib64 -Xlinker -R -Xlinker /opt/cuda/lib64 -lcuda -lcudart
  NVCC     = nvcc
  CUDA_OPTIMIZE = -g -G -O3 --ptxas-options=-v -Xptxas="-v" -arch=sm_20 $(filter -I%, $(shell mpicc -show))
endif


ifeq ($(SYSTYPE), "JuQueen")
  CC       = mpixlc_r -g
  OPTIMIZE = -std=c11 -qstrict -q64 -qcpluscmt -O3 -qipa
  GMP_INCL = -I/bgsys/local/gmp/5.0.5/include
  GMP_LIBS = -L/bgsys/local/gmp/5.0.5/lib
  GSL_INCL = -I/bgsys/local/gsl/1.15_O3g/include
  GSL_LIBS = -L/bgsys/local/gsl/1.15_O3g/lib
  FFTW_INCL= -I/bgsys/local/fftw2/2.1.5/include
  FFTW_LIBS= -L/bgsys/local/fftw2/2.1.5/lib -q64 # -qipa
  MPICHLIB =
  HDF5INCL = -I/bgsys/local/hdf5/include
  HDF5LIB  = -L/bgsys/local/hdf5/lib -lhdf5 -lz
endif


ifeq ($(SYSTYPE), "Leibniz")
  CC       = mpiicc
  OPTIMIZE = -std=c11 -g -O3 -Wall
  GSL_INCL = $(GSL_INC)
  GSL_LIBS = $(GSL_LIB)
  FFTW_INCL= $(FFTW_INC)
  FFTW_LIBS= -L$(FFTW_LIBDIR)
  HWLOC_INCL = $(HWLOC_INC)
  HWLOC_LIB = -lhwloc
  MPICHLIB =
  HYPRE_INCL = -I/opt/hypre/2.11.2/include
  HYPRE_LIB = -L/opt/hypre/2.11.2/lib -lHYPRE
  HDF5INCL = $(HDF5_INC) -DH5_USE_16_API=1
  HDF5LIB  = $(HDF5_LIB) -lhdf5
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE +=  -openmp
  else
    OPTIMIZE += -diag-disable 3180
  endif
endif


ifeq ($(SYSTYPE), "lonestar")
  CC        = mpicc
  OPTIMIZE  = -std=c11 -O2 -m64 -Wno-uninitialized
  GSL_INCL  = -I$(TACC_GSL_INC)
  GSL_LIBS  = -L$(TACC_GSL_LIB)
  FFTW_INCL = -I$(TACC_FFTW2_INC)
  FFTW_LIBS = -L$(TACC_FFTW2_LIB)
  MPICHLIB  =
  HDF5INCL  = -I$(TACC_HDF5_INC)
  HDF5LIB   = -L$(TACC_HDF5_LIB) -lhdf5 -lz
  HWLOC_INCL= -I/home1/00025/tgreif/libs/hwloc/include
  HWLOC_LIB = -L/home1/00025/tgreif/libs/hwloc/lib -lhwloc
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -openmp
  else
    OPTIMIZE += -Wno-unknown-pragmas
  endif
  OPT      += -DNOCALLSOFSYSTEM
endif


ifeq ($(SYSTYPE), "MacOS")
  CC       = mpicc   # sets the C-compiler
  FC       = mpif90  # sets fortran compiler
  OPTIMIZE = -std=c11 -O3 -g -w -Wall
  GSL_INCL =
  GSL_LIBS =
  FFTW_INCL=
  FFTW_LIBS=
  MPICHLIB = -lmpi
  HDF5INCL = -I/opt/local/lib/hdf5-18/include -DH5_USE_16_API
  HDF5LIB  = -L/opt/local/lib/hdf5-18/lib -lhdf5
  LINKER   = $(FC) -Wno-main
endif


ifeq ($(SYSTYPE), "Magny")
  # modules for Magny
  # module add mvapich2/gcc/64/1.6-qlc
  CC       = mpicc
  FC       = mpif90
  OPTIMIZE = -std=c11 -O3 -msse3 -g -Wall -m64
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -fopenmp
  else
    OPTIMIZE += -Wno-unknown-pragmas -Wno-unused-function
  endif
  GMP_INCL = -I/hits/tap/sw/libs/gmp-5.0.5/include
  GMP_LIBS = -L/hits/tap/sw/libs/gmp-5.0.5/lib -Xlinker -R -Xlinker /hits/tap/sw/libs/gmp-5.0.5/lib
  GSL_INCL = -I/hits/tap/sw/libs/gsl-1.15/include
  GSL_LIBS = -L/hits/tap/sw/libs/gsl-1.15/lib -Xlinker -R -Xlinker /hits/tap/sw/libs/gsl-1.15/lib
  FFTW_INCL= -I/hits/tap/sw/libs/fftw-3.3.4/include
  FFTW_LIBS= -L/hits/tap/sw/libs/fftw-3.3.4/lib -Xlinker -R -Xlinker /hits/tap/sw/libs/fftw-3.3.4/lib
  MPICHLIB =
  HDF5INCL = -I/hits/tap/sw/libs/hdf5-1.8.10/include -DH5_USE_16_API
  HDF5LIB  = -L/hits/tap/sw/libs/hdf5-1.8.10/lib -lhdf5 -Xlinker -R -Xlinker /hits/tap/sw/libs/hdf5-1.8.10/lib
  HWLOC_INCL=
  HWLOC_LIB = -lhwloc
  NBC_INCL  = -I/hits/tap/sw/libs/libNBC-1.1.1/include
  NBC_LIB   = -L/hits/tap/sw/libs/libNBC-1.1.1/lib -lnbc
  LINKER    = $(CC) -lgfortran #-nofor-main
  HYPRE_LIB = -lHYPRE
  VTUNE_INCL = -I/cm/shared/apps/intel/vtune_u13/vtune_amplifier_xe_2015/include
  VTUNE_LIBS = -L/cm/shared/apps/intel/vtune_u13/vtune_amplifier_xe_2015/lib64 -littnotify
endif


ifeq ($(SYSTYPE), "Magny-Intel")
  # modules for Magny-Intel
  # module load intel/compiler
  # module load mvapich2/intel/64/1.6-qlc
  CC       = mpiicc
  OPTIMIZE = -std=c11 -O3 -g -Wall -m64
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -openmp
  else
    OPTIMIZE += -Wno-unknown-pragmas
  endif
  GSL_INCL = -I/hits/tap/sw/libs/include
  GSL_LIBS = -L/hits/tap/sw/libs/lib -Xlinker -R -Xlinker /hits/tap/sw/libs/lib
  FFTW_INCL=
  FFTW_LIBS=
  MPICHLIB =
  MATHLIB  = -limf -lm
  HDF5INCL = -DH5_USE_16_API
  HDF5LIB  = -lhdf5
  VTUNE_INCL = -I/cm/shared/apps/intel/vtune_u13/vtune_amplifier_xe_2015/include
  VTUNE_LIBS = -L/cm/shared/apps/intel/vtune_u13/vtune_amplifier_xe_2015/lib64 -littnotify
endif


ifeq ($(SYSTYPE), "Milkyway")
  CC       = mpicc
  FC       = mpif90
  OPTIMIZE = -std=c11 -O3 -g
  MPICHLIB = -lmpi
  LINKER   = $(FC)
endif


ifeq ($(SYSTYPE), "mira")
  CC       = mpixlc   # sets the C compiler
  OPTIMIZE = -std=c11 -g -O2 -qarch=qp -qtune=qp
  GMP_INCL = -I/home/mvogelsb/gmp/include
  GMP_LIBS = -L/home/mvogelsb/gmp/lib
  GSL_INCL = -I/soft/libraries/unsupported/gsl/1.9/xl/include
  GSL_LIBS = -L/soft/libraries/unsupported/gsl/1.9/xl/lib -lgsl
  FFTW_INCL= -I/home/mvogelsb/fftw/include
  FFTW_LIBS= -L/home/mvogelsb/fftw/lib
  HWLOC_INCL=
  MPICHLIB =  -Wl,-dy -L/bgsys/drivers/ppcfloor/gnu-linux/powerpc64-bgq-linux/lib -lstdc++
  HDF5INCL = -I/soft/libraries/hdf5/1.8.10/cnk-xl/current/include -DH5_USE_16_API
  HDF5LIB  = -L/soft/libraries/hdf5/1.8.10/cnk-xl/current/lib -lhdf5 -L/soft/libraries/alcf/current/xl/ZLIB/lib -lz
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -openmp
  endif
endif


ifeq ($(SYSTYPE), "monstrum")
  CC       = mpicc     # sets the C compiler
  OPT     += -DMPICH_IGNORE_CXX_SEEK
  OPTIMIZE = -std=c11 -O3 -g -Wall #-Wno-unused-but-set-variable -Wno-uninitialized -Wno-unknown-pragmas -Wno-unused-function
  GSL_INCL =
  GSL_LIBS =
  FFTW_INCL=
  FFTW_LIBS=
  MPICHLIB =
  HDF5INCL = -DH5_USE_16_API
  HDF5LIB  = -lhdf5 -lz
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -fopenmp
    OPT      += -DIMPOSE_PINNING -DSOCKETS=4 -DMAX_CORES=16
  endif
endif


ifeq ($(SYSTYPE), "MPCDF-Intel")
  # see below for definition
  systype_mpcdf_intel = 1
endif


ifeq ($(SYSTYPE), "MPCDF-OpenMPI")
  # see below for definition
  systype_mpcdf_open_mpi = 1
endif


ifeq ($(SYSTYPE), "Odin")
  CC       = mpiicc
  OPTIMIZE = -std=c11 -g -O2
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -openmp
  else
    OPTIMIZE += -Wno-unknown-pragmas
  endif
  GMP_INCL =
  GMP_LIBS =
  GSL_INCL = -I$(GSL_INCDIR)
  GSL_LIBS = -L$(GSL_LIBDIR) -Xlinker -R -Xlinker  $(GSL_LIBDIR)
  FFTW_INCL= -I$(FFTW_HOME)/include
  FFTW_LIBS= -L$(FFTW_HOME)/lib -Xlinker -R -Xlinker $(FFTW_HOME)/lib
  HDF5INCL = -I$(HDF5_HOME)/include -DH5_USE_16_API
  HDF5LIB  = -L$(HDF5_HOME)/lib -lhdf5 -lz  -Xlinker -R -Xlinker $(HDF5_HOME)/lib
  HWLOC_INCL= -I/u/vrs/Libs/include
  HWLOC_LIB = -L/u/vrs/Libs/lib -lhwloc -Xlinker -R -Xlinker /u/vrs/Libs/lib
  MPICHLIB =
  OPT      += -DNOCALLSOFSYSTEM
endif


ifeq ($(SYSTYPE), "odyssey")
  CC       = mpicc     # sets the C compiler
  OPT     += -DMPICH_IGNORE_CXX_SEEK
  OPTIMIZE = -std=c11 -O3 -g -Wall -Wno-unused-but-set-variable -Wno-uninitialized -Wno-unused-function
  GSL_INCL =
  GSL_LIBS =
  FFTW_INCL= -I/n/home12/mvogelsberger/opt/include/
  FFTW_LIBS= -L/n/home12/mvogelsberger/opt/lib/
  MPICHLIB =
  HDF5INCL = -DH5_USE_16_API
  HDF5LIB  = -lhdf5
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -fopenmp
    OPT      += -DIMPOSE_PINNING -DSOCKETS=4 -DMAX_CORES=16
  else
    OPTIMIZE += -Wno-unknown-pragmas
  endif
endif

ifeq ($(SYSTYPE), "odyssey2")
  CC       = mpicc     # sets the C compiler
  OPT     += -DMPICH_IGNORE_CXX_SEEK
  OPTIMIZE = -std=c11 -O3 -g -Wall -Wno-unused-but-set-variable -Wno-uninitialized -Wno-unknown-pragmas -Wno-unused-function -mprefer-avx128 -march=native
  GSL_INCL =
  GSL_LIBS =
  FFTW_INCL=
  FFTW_LIBS=
  MPICHLIB =
  HDF5INCL = -DH5_USE_16_API
  HDF5LIB  = -lhdf5 -lz
  ifeq (MONOTONE_CONDUCTION, $(findstring MONOTONE_CONDUCTION, $(CONFIGVARS)))
    HYPRE_INCL = -I/n/home13/kannan/hypre-2.10.0b/src/hypre/include
    HYPRE_LIB = -L/n/home13/kannan/hypre-2.10.0b/src/hypre/lib -lHYPRE
  endif
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -fopenmp
    OPT      += -DIMPOSE_PINNING -DSOCKETS=4 -DMAX_CORES=16
  endif
endif


ifeq ($(SYSTYPE), "odyssey-gcc")
  CC       = mpicc     # sets the C compiler
  OPT     += -DMPICH_IGNORE_CXX_SEEK
  #OPTIMIZE = -std=c11 -O3 -g -Wall -Wno-unused-but-set-variable -Wno-uninitialized -Wno-unknown-pragmas -Wno-unused-function -mprefer-avx128 -march=native
  OPTIMIZE = -std=c11 -O3 -g -Wall -Wno-unused-but-set-variable -Wno-uninitialized -Wno-unknown-pragmas -Wno-unused-function -mprefer-avx128
  GSL_INCL =
  GSL_LIBS =
  FFTW_INCL= #-I/n/home13/kannan/fftw-3.3.4/include
  FFTW_LIBS= #-L/n/home13/kannan/fftw-3.3.4/lib
  MPICHLIB =
  HDF5INCL = -DH5_USE_16_API
  HDF5LIB  = -lhdf5 -lz
  ifeq (TGCHEM, $(findstring TGCHEM, $(CONFIGVARS)))
    CVODE_INCL = -I/n/home00/fmarinacci/sundials/install/include/
    CVODE_LIB  = -L/n/home00/fmarinacci/sundials/install/lib -lsundials_cvode -lsundials_nvecserial
  endif
  ifeq (MRT, $(findstring MRT, $(CONFIGVARS)))
    CVODE_INCL = -I/n/home12/mvogelsberger/sundials/instdir/include/
    CVODE_LIB  = -L/n/home12/mvogelsberger/sundials/instdir/lib -Xlinker -R -Xlinker /n/home12/mvogelsberger/sundials/instdir/lib -lsundials_cvode -lsundials_nvecserial
  endif
  ifeq (MONOTONE_CONDUCTION, $(findstring MONOTONE_CONDUCTION, $(CONFIGVARS)))
    HYPRE_INCL = -I/n/home13/kannan/hypre-2.10.0b/src/hypre/include
    HYPRE_LIB  = -L/n/home13/kannan/hypre-2.10.0b/src/hypre/lib -lHYPRE
  endif
  ifeq (IMPLICIT_OHMIC_DIFFUSION, $(findstring IMPLICIT_OHMIC_DIFFUSION, $(CONFIGVARS)))
    HYPRE_INCL = -I/n/home13/kannan/hypre-2.10.0b/src/hypre/include
    HYPRE_LIB  = -L/n/home13/kannan/hypre-2.10.0b/src/hypre/lib -lHYPRE
  endif
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE  += -fopenmp
    OPT       += -DIMPOSE_PINNING -DSOCKETS=4 -DMAX_CORES=16
  endif
endif


ifeq ($(SYSTYPE), "odyssey-intel")
  CC       = mpiicc
  OPT     += -DMPICH_IGNORE_CXX_SEEK
  OPTIMIZE = -std=c99 -O2 -parallel -ipo -funroll-loops
  GSL_INCL =
  GSL_LIBS =
  FFTW_INCL= -I/n/home03/cpopa/fftw3/include
  FFTW_LIBS= -L/n/home03/cpopa/fftw3/lib
  MPICHLIB =
  HDF5INCL = -DH5_USE_16_API
  HDF5LIB  = -lhdf5 -lz
  ifeq (MONOTONE_CONDUCTION, $(findstring MONOTONE_CONDUCTION, $(CONFIGVARS)))
    HYPRE_INCL = -I/n/home13/kannan/hypre-2.10.0b/src/hypre/include
    HYPRE_LIB = -L/n/home13/kannan/hypre-2.10.0b/src/hypre/lib -lHYPRE
  endif
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -fopenmp
    OPT      += -DIMPOSE_PINNING -DSOCKETS=4 -DMAX_CORES=16
  else
    OPTIMIZE += -diag-disable 3180
  endif
endif


ifeq ($(SYSTYPE), "odyssey-opteron")
  CC       = mpiicc
  OPT     +=
  OPTIMIZE = -std=c11 -O2 -m64 -Wno-uninitialized -Wno-unknown-pragmas
  GSL_INCL =
  GSL_LIBS =
  FFTW_INCL=
  FFTW_LIBS=
  MPICHLIB =
  MATHLIB  =
  HDF5INCL = -I/n/home12/mvogelsberger/optmine/include -DH5_USE_16_API
  HDF5LIB  = -L/n/home12/mvogelsberger/optmine/lib -Xlinker -R -Xlinker /n/home12/mvogelsberger/optmine/lib -lhdf5 -lz
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -fopenmp
    OPT      += -DIMPOSE_PINNING -DSOCKETS=4 -DMAX_CORES=16
  else
    OPTIMIZE += -Wno-unknown-pragmas
  endif
endif


ifeq ($(SYSTYPE), "OPA")
  CC       = mpicc   # sets the C compiler
  OPTIMIZE = -std=c11 -ggdb -O3 -ftree-vectorize -march=core2 -mcx16 -msahf -mno-movbe -mno-aes -mno-pclmul -mno-popcnt -mno-abm -mno-lwp -mno-fma -mno-fma4 -mno-xop -mno-bmi -mno-bmi2 -mno-tbm -mno-avx -mno-avx2 -mno-sse4.2 -mno-sse4.1 -mno-lzcnt -mno-rdrnd -mno-f16c -mno-fsgsbase --param l1-cache-size=32 --param l1-cache-line-size=64 --param l2-cache-size=4096 -mtune=core2
  GSL_INCL =
  GSL_LIBS =
  FFTW_INCL=
  FFTW_LIBS=
  MPICHLIB = -lmpi -lhwloc
  HDF5INCL = -DH5_USE_16_API=1
  HDF5LIB  = -lhdf5
endif


ifeq ($(SYSTYPE), "OpenSuse")
  CC       = mpicc -g -O2 -Wall
  OPTIMIZE = -std=c11
  GSL_INCL =
  GSL_LIBS =
  FFTW_INCL=
  FFTW_LIBS=
  MPICHLIB = -L/usr/lib/mpi/gcc/openmpi/lib -Xlinker -R -Xlinker /usr/lib/mpi/gcc/openmpi/lib -lmpi
  HDF5INCL = -DH5_USE_16_API
  HDF5LIB  =
endif


ifeq ($(SYSTYPE), "OpenSuse64")
  CC       = mpicc
  OPTIMIZE = -std=c11 -g -O3 -Wall
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -fopenmp
  else
    OPTIMIZE += -Wno-unknown-pragmas
  endif
  GSL_INCL = -I/sw/tap/include
  GSL_LIBS = -L/sw/tap/lib -Xlinker -R -Xlinker /sw/tap/lib
  FFTW_INCL=
  FFTW_LIBS=
  MPICHLIB =
  HDF5INCL = -DH5_USE_16_API
  HDF5LIB  = -lhdf5 -lhdf5_hl
endif


ifeq ($(SYSTYPE), "OpenSuse64-cuda")
  CC       = mpicc
  LINKER   = mpicxx
  OPTIMIZE = -std=c11 -g -O3 -Wall
  GSL_INCL = -I/sw/tap/include
  GSL_LIBS = -L/sw/tap/lib -Xlinker -R -Xlinker /sw/tap/lib
  FFTW_INCL=
  FFTW_LIBS=
  MPICHLIB =
  HDF5INCL = -DH5_USE_16_API
  HDF5LIB  = -lhdf5
  CUDA_INCL= -I/usr/local/cuda/include
  CUDA_LIBS= -L/usr/local/cuda/lib64 -lcuda -lcudart
  NVCC     = /usr/local/cuda/bin/nvcc
  CUDA_OPTIMIZE = -g -G -O0 --ptxas-options=-v -Xptxas="-v" -arch=sm_20 $(filter -I%, $(shell mpicc -show))
endif


ifeq ($(SYSTYPE), "pascal")
  # Linux cluster pascal at the Max Planck Institute for Astrophysics (MPA), Garching:
  #   https://wwwmpa.mpa-garching.mpg.de/mpa/internal/comp_support/hardware/computers/cs_comp_list-en.html
  # module load gcc mpich fftw-mpich hdf5
  CC       = mpicc
  FC       = mpif90
  OPTIMIZE = -std=c99 -O2 -g
  GSL_INCL =
  GSL_LIBS = -lgsl
  FFTW_INCL=
  FFTW_LIBS= -lfftw3
  HDF5INCL = -I$(HDF5_HOME)/include -DH5_USE_16_API
  HDF5LIB  = -L$(HDF5_HOME)/lib -lhdf5 -Xlinker -R -Xlinker $(HDF5_HOME)/lib
  MPICHLIB = -lmpi
  HWLOC_INCL= -I/usr/include
  HWLOC_LIB = $(LDFLAGS) -lhwloc
  LINKER   = $(FC)
endif


ifeq ($(SYSTYPE), "Peta4")
  CC         = mpicc
  OPTIMIZE   = -std=c11 -g -O3 -xHost -ipo
  GMP_INCL   =
  GMP_LIBS   =
  GSL_INCL   =
  GSL_LIBS   =
  FFTW_INCL  =
  FFTW_LIBS  =
  HDF5INCL   = -DH5_USE_16_API
  HDF5LIB    = -lhdf5 -lz
  HWLOC_INCL =
  HWLOC_LIB  = -lhwloc
  MPICHLIB   = -lmpi
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -openmp
  else
    OPTIMIZE += -Wno-unknown-pragmas
  endif
endif


ifeq ($(SYSTYPE), "PhilipNUC")
  CC       = mpicc
  CPP      = mpicxx   # sets the C++ compiler
  OPTIMIZE = -std=c11 -O3 -Wall -m64 -g -Wno-error=format-security -Wno-unknown-pragmas
  GSL_INCL = -I/usr/include
  GSL_LIBS = -L/usr/lib -lgsl -lm
  FFTW_INCL= -I/usr/include
  FFTW_LIBS= -L/usr/lib -lfftw3
  MPICHLIB = -lmpich
  HDF5INCL = -I/usr/local/include
  HDF5LIB  = -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -lz
  HWLOC_INCL= -I/usr/include/hdf5/serial -DH5_USE_16_API
  HWLOC_LIB = -lhwloc -Xlinker -R -Xlinker /usr/lib
  CPV      = $(CC)
  LINKER = mpicxx
endif

ifeq ($(SYSTYPE), "Puma")
  # compiler and its optimization options
  OPTIMIZE  = -std=c11 -ggdb -O3 -Wall -Wno-format-security -Wno-unknown-pragmas -Wno-unused-function

  # overwrite default:
  OPT      += -I/opt/ohpc/pub/mpi/openmpi3-gnu8/3.1.4/include
  MPICHLIB  = -L/opt/ohpc/pub/mpi/openmpi3-gnu8/3.1.4/lib -lmpi
  GSL_INCL  = -I/opt/ohpc/pub/libs/gnu8/gsl/2.6/include
  GSL_LIBS  = -L/opt/ohpc/pub/libs/gnu8/gsl/2.6/lib -lgsl -lgslcblas
  HWLOC_LIB = -L/opt/ohpc/pub/apps/hwloc/2.2.0/lib -lhwloc
  GMP_LIBS  = -lgmp

  # libraries that are included on demand, depending on Config.sh options
  FFTW_INCL = # in Puma's OS, no path needed
  FFTW_LIBS = # in Puma's OS, no path needed
  HDF5INCL  = -I/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/include -DH5_USE_16_API
  HDF5LIB   = -L/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/lib -lhdf5 -lz
  HWLOC_INCL= -I/opt/ohpc/pub/apps/hwloc/2.2.0/include
endif

ifeq ($(SYSTYPE), "Quest-intel")
  CC       = mpicc   # sets the C compiler
  FC       = $(CC)
  OPTIMIZE = -std=c99 -O2
  ifeq (CHIMES, $(findstring CHIMES, $(CONFIGVARS)))
    CHIMESINCL = -I/home/ajr882/sundials/include
    CHIMESLIBS = -L/home/ajr882/sundials/lib -lsundials_cvode -lsundials_kinsol -lsundials_nvecserial
    ifeq (CHIMES_PTHREADS, $(findstring CHIMES_PTHREADS, $(CONFIGVARS)))
      OPTIMIZE += -pthread
    endif
  endif
  GSL_INCL =
  GSL_LIBS =
  FFTW_INCL=
  FFTW_LIBS=
  MPICHLIB = -lmpi
  HDF5INCL = -DH5_USE_16_API
  HDF5LIB  = -lhdf5
  HWLOC_INCL = -I/software/hwloc/1.4.2/include
  HWLOC_LIB = -L/software/hwloc/1.4.2/lib -lhwloc
endif


ifeq ($(SYSTYPE), "Raijin")
  CC       = mpicc   # sets the C compiler
  OPTIMIZE = -std=c11 -g -O3 -ipo -m64 -Wall -xCORE-AVX2
  GSL_INCL =
  GSL_LIBS =
  FFTW_INCL=
  FFTW_LIBS=
  MPICHLIB = -lmpi
  HDF5INCL = -DH5_USE_16_API
  HDF5LIB  = -lhdf5
endif


ifeq ($(SYSTYPE), "Ranger_intel" )
  AR = xiar
endif


ifeq ($(SYSTYPE), "Raven-Intel")
  # alias for SYSTYPE "MPCDF-Intel" (with some extra settings)
  systype_mpcdf_intel = 1
endif


ifeq ($(SYSTYPE), "Raven-OpenMPI")
  # alias for SYSTYPE "MPCDF-OpenMPI" (with some extra settings)
  systype_mpcdf_open_mpi = 1
endif


ifeq ($(SYSTYPE), "ScientificLinux")
  CC       = mpicc   # sets the C compiler
  OPTIMIZE = -std=c11 -ggdb -O2 -Wno-uninitialized
  GSL_INCL = -L/usr/lib64/
  GSL_LIBS = -L/usr/lib64/
  FFTW_INCL= -I/home/djmunoz/opt/include
  FFTW_LIBS= -L/home/djmunoz/opt/lib
  MPICHLIB = -L/home/djmunoz/opt/lib -lmpi
  HDF5INCL = -L/home/djmunoz/opt/lib/ -DH5_USE_16_API
  HDF5LIB  = -L/home/djmunoz/opt/lib/ -lhdf5
endif


ifeq ($(SYSTYPE), "stampede")
  CC        = mpicc
  OPTIMIZE  = -std=c99 -O2 -xhost -Wno-uninitialized -Wno-unknown-pragmas
  GMP_INCL  =
  GMP_LIBS  =
  GSL_INCL  = -I$(TACC_GSL_INC)
  GSL_LIBS  = -L$(TACC_GSL_LIB)
  FFTW_INCL = -I$(TACC_FFTW3_INC)
  FFTW_LIBS = -L$(TACC_FFTW3_LIB)
  MPICHLIB  =
  HDF5INCL  = -I$(TACC_HDF5_INC) -DH5_USE_16_API
  HDF5LIB   = -L$(TACC_HDF5_LIB) -lhdf5 -lz
  HWLOC_INCL=
  HWLOC_LIB =
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -openmp
  endif
    OPT      += -DNOCALLSOFSYSTEM
endif


ifeq ($(SYSTYPE), "stampede_illustris")
  CC        = mpicc
  OPTIMIZE  = -std=c11 -O2 -m64 -Wno-uninitialized
  GMP_INCL  = -I/home1/01637/mvogelsb/libs/gmp/include
  GMP_LIBS  = -L/home1/01637/mvogelsb/libs/gmp/lib -Xlinker -R -Xlinker /home1/01637/mvogelsb/libs/gmp/lib
  GSL_INCL  = -I$(TACC_GSL_INC)
  GSL_LIBS  = -L$(TACC_GSL_LIB)
  FFTW_INCL = -I$(TACC_FFTW2_INC)
  FFTW_LIBS = -L$(TACC_FFTW2_LIB)
  MPICHLIB  =
  HDF5INCL  = -I/home1/01637/mvogelsb/libs/hdf5/include/
  HDF5LIB   = -L/home1/01637/mvogelsb/libs/hdf5/lib/ -lhdf5 -lz -Xlinker -R -Xlinker /home1/01637/mvogelsb/libs/hdf5/lib/
  HWLOC_INCL= -I/home1/01637/mvogelsb/libs/hwloc/include
  HWLOC_LIB = -L/home1/01637/mvogelsb/libs/hwloc/lib -lhwloc -Xlinker -R -Xlinker /home1/01637/mvogelsb/libs/hwloc/lib
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -openmp
  else
    OPTIMIZE += -Wno-unknown-pragmas
  endif
endif


ifeq ($(SYSTYPE), "stampede2-SKX")
  CC        = mpicc
  OPTIMIZE  = -std=c99 -O0 -g
  GMP_INCL  =
  GMP_LIBS  =
  GSL_INCL  = -I$(TACC_GSL_INC)
  GSL_LIBS  = -L$(TACC_GSL_LIB)
  FFTW_INCL = -I/home1/04121/tg834202/libs_skx/fftw3/include
  FFTW_LIBS = -L/home1/04121/tg834202/libs_skx/fftw3/lib
  MPICHLIB  =
  HDF5INCL  = -I$(TACC_HDF5_INC) -DH5_USE_16_API
  HDF5LIB   = -L$(TACC_HDF5_LIB) -lhdf5 -lz
  HWLOC_INCL=
  HWLOC_LIB =
  ifeq (CHIMES, $(findstring CHIMES, $(CONFIGVARS)))
    CHIMESINCL = -I/home1/04121/tg834202/libs_skx/sundials/include
    CHIMESLIBS = -L/home1/04121/tg834202/libs_skx/sundials/lib -lsundials_cvode -lsundials_kinsol -lsundials_nvecserial
    ifeq (CHIMES_PTHREADS, $(findstring CHIMES_PTHREADS, $(CONFIGVARS)))
      OPTIMIZE += -pthread
    endif
  endif
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -openmp
  endif
  OPT      += -DUSE_MPI_IN_PLACE
endif

# Sunnyvale @CITA; used by Rainer Weinberger
ifeq ($(SYSTYPE), "Sunnyvale")
  CC       = mpicc   # sets the C compiler
  FC       = mpif90
  CP       = mpicpc
  OPTIMIZE = -std=c11 -O3 -Wall
  GSL_INCL = 
  GSL_LIBS = 
  FFTW_INCL= 
  FFTW_LIBS= 
  HWLOC_INCL = 
  HWLOC_LIB = -lhwloc
  MPICHLIB =
  HDF5INCL =  -DH5_USE_16_API
  HDF5LIB  =  -L/opt/hdf5/1.8.20-gcc-7.3.0/lib -lhdf5 -lz
  LINKER   = $(FC)  
endif

ifeq ($(SYSTYPE), "SuperMuc")
  # module load fftw   needed on SuperMUC
  CC         = mpicc   # sets the C compiler
  OPTIMIZE   = -std=c11 -xAVX -g -O3 -Wall
  GSL_INCL   = $(GSL_INC)
  GSL_LIBS   = -L/lrz/sys/libraries/gsl/1.16/lib/
  FFTW_INCL  = $(FFTW_INC)
  FFTW_LIBS  = -L$(FFTW_LIBDIR)
  HWLOC_INCL = $(HWLOC_INC)
  HWLOC_LIB  = -L/lrz/sys/tools/hwloc/1.8.1/lib -lhwloc
  MPICHLIB   =
  HDF5INCL   = $(HDF5_INC) -DH5_USE_16_API=1
  HDF5LIB    = $(HDF5_LIB) -lhdf5
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -openmp
  else
    OPTIMIZE += -diag-disable 3180
  endif
endif


ifeq ($(SYSTYPE), "SuperMucIntel")
  CC       = mpicc   # sets the C compiler
  FC       = mpif90
  OPTIMIZE = -std=c11 -xAVX -g -O3 -Wall
  GSL_INCL = $(GSL_INC)
  GSL_LIBS = -L/lrz/sys/libraries/gsl/1.16/lib/
  FFTW_INCL= $(FFTW_INC)
  FFTW_LIBS= -L$(FFTW_LIBDIR)
  HWLOC_INCL = $(HWLOC_INC)
  HWLOC_LIB = -L/lrz/sys/tools/hwloc/1.8.1/lib -lhwloc
  MPICHLIB =
  HDF5INCL = $(HDF5_INC) -DH5_USE_16_API=1
  HDF5LIB  = -L/lrz/sys/libraries/hdf5/1.8.15/serial/lib -lhdf5_hl_cpp -lhdf5_cpp #-L/lrz/sys/libraries/hdf5/1.8.12/serial_gpfs/lib -lhdf5
  LINKER   = $(FC) -nofor-main
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -openmp
  else
    OPTIMIZE += -diag-disable 3180
  endif
endif

ifeq ($(SYSTYPE), "SuperMucIntelFat")
  CC       = mpicc   # sets the C compiler
  FC       = mpif90
  OPTIMIZE = -std=c11 -g -O3 -Wall
  GMP_INCL = -I/home/hpc/pr53ka/di57jir/lib/include/
  GMP_LIBS = -L/home/hpc/pr53ka/di57jir/lib/lib/
  GSL_INCL = $(GSL_INC)
  GSL_LIBS = -L/lrz/sys/libraries/gsl/1.16/lib/
  FFTW_INCL= $(FFTW_INC)
  FFTW_LIBS= -L$(FFTW_LIBDIR)
  HWLOC_INCL = $(HWLOC_INC)
  HWLOC_LIB = -L/lrz/sys/tools/hwloc/1.8.1/lib -lhwloc
  MPICHLIB =
  HDF5INCL = $(HDF5_INC) -DH5_USE_16_API=1
  HDF5LIB  = -L/lrz/sys/libraries/hdf5/1.8.15/serial/lib -lhdf5_hl_cpp -lhdf5_cpp #-L/lrz/sys/libraries/hdf5/1.8.12/serial_gpfs/lib -\
  lhdf5
  LINKER   = $(FC) -nofor-main
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -openmp
  else
    OPTIMIZE += -diag-disable 3180
  endif
endif


ifeq ($(SYSTYPE), "SuperMuc-NG")
  CC       = mpicc   # sets the C compiler
  FC       = mpif90
  CP       = mpicpc
  OPTIMIZE = -std=c11 -ggdb -O3 -Wall -march=native
  GSL_INCL = $(GSL_INC)
  GSL_LIBS = -L$(GSL_LIBDIR)
  FFTW_INCL= $(FFTW_INC)
  FFTW_LIBS= -L$(FFTW_LIBDIR)
  HWLOC_INCL = $(HWLOC_INC)
  HWLOC_LIB = -L$(HWLOC_LIBDIR) -lhwloc
  MPICHLIB =
  HDF5INCL = $(HDF5_INC) -DH5_USE_16_API
  HDF5LIB  = $(HDF5_SHLIB)
  GMP_INCL = -I/dss/dsshome1/0F/di98fem3/libs/include
  GMP_LIBS = -L/dss/dsshome1/0F/di98fem3/libs/lib
  LINKER   = $(FC)
endif

ifeq ($(SYSTYPE), "SuperMuc-NG-GCC")
  CC       = mpicxx   # sets the C compiler
  FC       = mpif90
  CPPC     = mpicxx   # sets the C compiler
  OPTIMIZE = -g -std=c++11 -O3 -Wall -march=native -Wformat=0 -Wno-unknown-pragmas -Wno-unused-function -Wno-unused -fstack-protector-all -fno-omit-frame-pointer #-fsanitize=address
  GSL_INCL = $(GSL_INC)
  GSL_LIBS = -L$(GSL_LIBDIR)
  FFTW_INCL= $(FFTW_INC)
  FFTW_LIBS= -L$(FFTW_LIBDIR)
  HWLOC_INCL = $(HWLOC_INC)
  HWLOC_LIB = -L$(HWLOC_LIBDIR) -lhwloc
  MPICHLIB = -lbfd -ldl -lz
  HDF5INCL = $(HDF5_INC) -DH5_USE_16_API
  HDF5LIB  = $(HDF5_SHLIB) $(CFITSIO_LIB)
  GMP_INCL = $(GMP_INC)
  GMP_LIBS = $(GMP_SHLIB)
  LINKER   = $(CPPC) -fstack-protector-all -fno-omit-frame-pointer #-fsanitize=address
endif

ifeq ($(SYSTYPE), "SuperMuc-NG-GCC")
  CC       = mpicxx   # sets the C compiler
  FC       = mpif90
  CPPC     = mpicxx   # sets the C compiler
  OPTIMIZE = -g -std=c++11 -O3 -Wall -march=native -Wformat=0 -Wno-unknown-pragmas -Wno-unused-function -Wno-unused -fstack-protector-all -fno-omit-frame-pointer #-fsanitize=address
  GSL_INCL = $(GSL_INC)
  GSL_LIBS = -L$(GSL_LIBDIR)
  FFTW_INCL= $(FFTW_INC)
  FFTW_LIBS= -L$(FFTW_LIBDIR)
  HWLOC_INCL = $(HWLOC_INC)
  HWLOC_LIB = -L$(HWLOC_LIBDIR) -lhwloc
  MPICHLIB = -lbfd -ldl -lz
  HDF5INCL = $(HDF5_INC) -DH5_USE_16_API
  HDF5LIB  = $(HDF5_SHLIB) $(CFITSIO_LIB)
  GMP_INCL = $(GMP_INC)
  GMP_LIBS = $(GMP_SHLIB)
  HYPRE_LIB = -lHYPRE
  LINKER   = $(CPPC) -fstack-protector-all -fno-omit-frame-pointer #-fsanitize=address
endif


ifeq ($(SYSTYPE), "SuperMuc-NG-GCC")
  CC       = mpicxx   # sets the C compiler
  FC       = mpif90
  CPPC     = mpicxx   # sets the C compiler
  OPTIMIZE = -g -std=c++11 -O3 -Wall -march=native -Wformat=0 -Wno-unknown-pragmas -Wno-unused-function -Wno-unused -fstack-protector-all -fno-omit-frame-pointer #-fsanitize=address
  GSL_INCL = $(GSL_INC)
  GSL_LIBS = -L$(GSL_LIBDIR)
  FFTW_INCL= $(FFTW_INC)
  FFTW_LIBS= -L$(FFTW_LIBDIR)
  HWLOC_INCL = $(HWLOC_INC)
  HWLOC_LIB = -L$(HWLOC_LIBDIR) -lhwloc
  MPICHLIB = -lbfd -ldl -lz
  HDF5INCL = $(HDF5_INC) -DH5_USE_16_API
  HDF5LIB  = $(HDF5_SHLIB) $(CFITSIO_LIB)
  GMP_INCL = $(GMP_INC)
  GMP_LIBS = $(GMP_SHLIB)
  LINKER   = $(CPPC) -fstack-protector-all -fno-omit-frame-pointer #-fsanitize=address
endif


ifeq ($(SYSTYPE), "Ubuntu")
  # Ubuntu Linux distribution: https://www.ubuntu.com/
  CC        = mpicc  -std=c11
  CPP       = mpic++ -std=c++11
  FC        = mpifort
  OPT      += -Wall -Wextra -Wno-unused-parameter
  OPTIMIZE  = -g -O2 -flto -march=native
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -fopenmp
  else
    OPTIMIZE += -Wno-unknown-pragmas
  endif
  GSL_INCL  =
  GSL_LIBS  =
  FFTW_INCL =
  FFTW_LIBS =
  MPICHLIB  =
  HDF5INCL  = -I/usr/include/hdf5/serial/ -DH5_USE_16_API
  HDF5LIB   = -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ -lhdf5
  ifeq (MRT_COUPLED_THERMOCHEMISTRY, $(findstring MRT_COUPLED_THERMOCHEMISTRY, $(CONFIGVARS)))
    CVODE_INCL= -I/home/kannan/codes/sundials-3.1.0/builddir/include
    CVODE_LIB = -L/home/kannan/codes/sundials-3.1.0/builddir/lib -lsundials_cvode -lsundials_nvecserial
  endif
  HYPRE_INCL= -I/home/kannan/codes/hypre-2.10.0b/src/hypre/include
  HYPRE_LIB = -L/home/kannan/codes/hypre-2.10.0b/src/hypre/lib -lHYPRE
endif


ifeq ($(SYSTYPE), "UbuntuSJ")
  CC       = mpicc    # sets the C compiler
  OPTIMIZE = -std=c11 -O3 -msse3 -g -Wall -m64 -Wno-format-security
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -fopenmp
  else
    OPTIMIZE += -Wno-unknown-pragmas -Wno-unused-function
  endif
  GSL_INCL =
  GSL_LIBS =
  FFTW_INCL= -lfftw3l -lm
  FFTW_LIBS= -lfftw3l -lm
  MPICHLIB =
  HDF5INCL = -DH5_USE_16_API
  HDF5LIB  = -lhdf5
endif


ifeq ($(SYSTYPE), "VIP")
  CC       = mpcc_r -g # -qflttrap=enable:zerodivide:nanq # sets the C compiler
  OPTIMIZE = -qstrict -q64 -qcpluscmt -O3 -qipa
  GSL_INCL = -I/afs/rzg/u/vrs/Libs/regatta/include
  GSL_LIBS = -L/afs/rzg/u/vrs/Libs/regatta/lib
  FFTW_INCL= -I/afs/rzg/u/vrs/Libs/regatta/include
  FFTW_LIBS= -L/afs/rzg/u/vrs/Libs/regatta/lib -q64 # -qipa
  MPICHLIB =
  HDF5INCL = -I/afs/rzg/u/vrs/Libs/regatta/include
  HDF5LIB  = -L/afs/rzg/u/vrs/Libs/regatta/lib -lhdf5 -lz
  OPT     += -DNOINLINE
endif


ifeq ($(SYSTYPE), "wuerzburg")
  CC       = mpicc   # sets the C compiler
  OPTIMIZE = -std=c11 -ggdb -O3 -march=native -ftree-vectorize
  GSL_INCL =
  GSL_LIBS =
  FFTW_INCL=
  FFTW_LIBS=
  MPICHLIB = -lmpi -lhwloc
  HDF5INCL = -DH5_USE_16_API=1
  HDF5LIB  = -lhdf5
endif


ifeq ($(SYSTYPE), "yeti")
  CC       =  mpicc     # sets the C compiler
  OPT     += -DMPICH_IGNORE_CXX_SEEK
  OPTIMIZE = -std=c99 -O3 -g -Wall -Wno-unused-but-set-variable -Wno-uninitialized -Wno-unknown-pragmas -Wno-unused-function -march=native
  GSL_INCL =
  GSL_LIBS =
  FFTW_INCL= -I/vega/opt/fftw2-2.1.5/include
  FFTW_LIBS= -L/vega/opt/fftw2-2.1.5/lib
  MPICHLIB =
  HDF5INCL = -DH5_USE_16_API
  HDF5LIB  = -lhdf5 -lz
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPTIMIZE += -fopenmp
    OPT      += -DIMPOSE_PINNING -DSOCKETS=4 -DMAX_CORES=16
  endif
endif


ifeq ($(systype_mpcdf_intel), 1)
  # HPC systems Draco at Max Planck Computing and Data Facility (MPCDF),
  # Garching (Freya, Cobra, Raven):
  #   https://docs.mpcdf.mpg.de/doc/computing/clusters/systems/Astrophysics.html
  #   https://www.mpcdf.mpg.de/services/supercomputing/cobra
  #   https://www.mpcdf.mpg.de/services/supercomputing/raven
  # module load intel impi gsl fftw-serial hdf5-serial
  # disabled ICC compiler warnings/remarks:
  #   #424: extra ";" ignored
  #   #981: operands are evaluated in unspecified order
  #  #1418: external function definition with no prior declaration
  #  #1572: floating-point equality and inequality comparisons are unreliable
  #  #2259: non-pointer conversion may lose significant bits
  #  #2282: unrecognized GCC pragma
  # #11074: Inlining inhibited by limit max-size (or max-total-size)
  # #11076: To get full report use -qopt-report=4 -qopt-report-phase ipo
  CC       = mpiicc  -std=c11
  CPP      = mpiicpc -std=c++11
  OPT     += -w3 -Wno-unused-parameter -diag-disable=424,981,1418,1572,2259,2282,11074,11076
  LINKER_OPT += -Xlinker -rpath=$(GSL_HOME)/lib:$(FFTW_HOME)/lib:$(HDF5_HOME)/lib
  OPTIMIZE = -g -O2 -ipo -march=native
  GSL_INCL = -I$(GSL_HOME)/include/
  GSL_LIBS = -L$(GSL_HOME)/lib/
  FFTW_INCL= -I$(FFTW_HOME)/include/
  FFTW_LIBS= -L$(FFTW_HOME)/lib/
  MPICHLIB =
  HDF5INCL = -DH5_USE_16_API=1 -I$(HDF5_HOME)/include/
  HDF5LIB  = -lhdf5 -L$(HDF5_HOME)/lib/
  # OpenMP
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPT += -qopenmp
  else
    OPT += -Wno-unknown-pragmas
  endif
  # system-dependent settings
  ifeq ($(SYSTYPE), "Raven-Intel")
    # Raven requires access to the hwloc library through a module
    LINKER_OPT += -Xlinker -rpath=$(HWLOC_HOME)/lib
    HWLOC_INCL = -I$(HWLOC_HOME)/include/
    HWLOC_LIB  = -L$(HWLOC_HOME)/lib/ -lhwloc
  endif
endif


ifeq ($(systype_mpcdf_open_mpi), 1)
  # HPC systems Draco at Max Planck Computing and Data Facility (MPCDF),
  # Garching (Freya, Cobra, Raven):
  #   https://docs.mpcdf.mpg.de/doc/computing/clusters/systems/Astrophysics.html
  #   https://www.mpcdf.mpg.de/services/supercomputing/cobra
  #   https://www.mpcdf.mpg.de/services/supercomputing/raven
  # module load gcc openmpi gsl fftw-serial hdf5-serial
  CC       = mpicxx -std=c++11
  CPP      = mpicxx -std=c++11
  FC       = mpifort
  OPT     += -Wall -Wextra -Wno-unused-parameter
  LINKER_OPT += -Xlinker -rpath=$(GSL_HOME)/lib:$(FFTW_HOME)/lib:$(HDF5_HOME)/lib
  OPTIMIZE = -g -O2 -flto -march=native
  GSL_INCL = -I$(GSL_HOME)/include/
  GSL_LIBS = -L$(GSL_HOME)/lib/
  FFTW_INCL= -I$(FFTW_HOME)/include/
  FFTW_LIBS= -L$(FFTW_HOME)/lib/
  MPICHLIB = 
  HDF5INCL = -DH5_USE_16_API=1 -I$(HDF5_HOME)/include/
  HDF5LIB  = -lhdf5 -L$(HDF5_HOME)/lib/
  # OpenMP
  ifeq (NUM_THREADS, $(findstring NUM_THREADS, $(CONFIGVARS)))
    OPT += -fopenmp
  else
    OPT += -Wno-unknown-pragmas
  endif
  # system-dependent settings
  ifeq ($(SYSTYPE), "Cobra")
    OPTIMIZE += -O3 -march=skylake-avx512
  endif
  ifeq ($(SYSTYPE), "Raven-OpenMPI")
    # Raven requires access to the hwloc library through a module
    LINKER_OPT += -Xlinker -rpath=$(HWLOC_HOME)/lib
    HWLOC_INCL = -I$(HWLOC_HOME)/include/
    HWLOC_LIB  = -L$(HWLOC_HOME)/lib/ -lhwloc
  endif
endif

