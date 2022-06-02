#!/bin/bash            # this line only there to enable syntax highlighting in this file

## examples/OrszagTang/Config_common_nolabframe.sh
## config file for 2d OrszagTang decaying turbulence MHD problem 
## (OrszagTang_2d and OrszagTang_Boost_2d)


#--------------------------------------- Basic operation mode of code
VORONOI
CHUNKING
TWODIMS                                  # 2d simulation
MESHRELAX_DENSITY_IN_INPUT               # Reads the mass field in the IC as density

#--------------------------------------- MHD options
MHD                                      # enable MHD
MHD_POWELL                               # enable divergence cleaning
RIEMANN_HLLD                             # switch to HLLD Riemann solver

#--------------------------------------- Mesh motion and regularization
REGULARIZE_MESH_CM_DRIFT                 # Mesh regularization; Move mesh generating point towards center of mass to make cells rounder.
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED  # Limit mesh regularization speed by local sound speed
REGULARIZE_MESH_FACE_ANGLE               # Use maximum face angle as roundness criterion in mesh regularization
VORONOI_DYNAMIC_UPDATE                   # keeps track of mesh connectivity, which speeds up mesh construction
VORONOI_MESH_KEEP_DT_AND_DTC             # keep DT and DTC structures, required to output mesh

#--------------------------------------- Time integration options
FORCE_EQUAL_TIMESTEPS                    # global timesteps only

#---------------------------------------- Single/Double Precision
DOUBLEPRECISION=1                        # Mode of double precision: not defined: single; 1: full double precision 2: mixed, 3: mixed, fewer single precisions; unless short of memory, use 1.
INPUT_IN_DOUBLEPRECISION                 # initial conditions are in double precision
OUTPUT_IN_DOUBLEPRECISION                # snapshot files will be written in double precision
OUTPUT_CENTER_OF_MASS                    # output centers of cells
OUTPUT_PRESSURE                          # output pressure of cells

#--------------------------------------- Output/Input options
HAVE_HDF5                                # needed when HDF5 I/O support is desired (recommended)

#--------------------------------------- Testing and Debugging options
DEBUG                                    # enables core-dumps



