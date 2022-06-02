#!/bin/bash            # this line only there to enable syntax highlighting in this file

## examples/Noh_3d/Noh_refinement_3d/Config.sh
## config file for 3d Noh problem

#--------------------------------------- Basic operation mode of code
VORONOI
CHUNKING
REFLECTIVE_X=2                           # in-/outflow boundary conditions in x direction
REFLECTIVE_Y=2                           # in-/outflow boundary conditions in y direction
REFLECTIVE_Z=2                           # in-/outflow boundary conditions in z direction

#--------------------------------------- Mesh motion and regularization
REGULARIZE_MESH_CM_DRIFT                 # Mesh regularization; Move mesh generating point towards center of mass to make cells rounder.
REGULARIZE_MESH_FACE_ANGLE               # Use maximum face angle as roundness criterion in mesh regularization
VORONOI_DYNAMIC_UPDATE                   # keeps track of mesh connectivity, which speeds up mesh construction

#--------------------------------------- Explicit refinement
REFINEMENT_SPLIT_CELLS                   # enable refinement
REFINEMENT_MERGE_CELLS                   # enable derefinement
REFINEMENT_VOLUME_LIMIT                  # enforces minimum and maximum volume for cells by refinement

#--------------------------------------- Time integration options
TREE_BASED_TIMESTEPS                     # non-local timestep criterion (take 'signal speed' into account)

#---------------------------------------- Single/Double Precision
DOUBLEPRECISION=1                        # Mode of double precision: not defined: single; 1: full double precision 2: mixed, 3: mixed, fewer single precisions; unless short of memory, use 1.
INPUT_IN_DOUBLEPRECISION                 # initial conditions are in double precision
OUTPUT_CENTER_OF_MASS                    # output centers of cells

#--------------------------------------- Output/Input options
HAVE_HDF5                                # needed when HDF5 I/O support is desired; should this be standard?

#--------------------------------------- Testing and Debugging options
DEBUG                                    # enables core-dumps, should this be standard?