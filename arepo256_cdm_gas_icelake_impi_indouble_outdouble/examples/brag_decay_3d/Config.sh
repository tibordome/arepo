#!/bin/bash            # this line only there to enable syntax highlighting in this file
MHD
NOHYDRO
NOHYDRO_NOTIMESTEP
VORONOI
RIEMANN_HLLD
VORONOI_STATIC_MESH
TREE_BASED_TIMESTEPS
MESHRELAX_DENSITY_IN_INPUT
CHUNKING
DOUBLEPRECISION=1
OUTPUT_IN_DOUBLEPRECISION
INPUT_IN_DOUBLEPRECISION
TETRA_INDEX_IN_FACE
VORONOI_DYNAMIC_UPDATE
# NO_MPI_IN_PLACE
# NO_ISEND_IRECV_IN_DOMAIN
# FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
# OUTPUT_CENTER_OF_MASS
HAVE_HDF5
DEBUG


FORCE_EQUAL_TIMESTEPS
VORONOI_MESH_KEEP_DT_AND_DTC
# UNLIMITED_GRADIENTS

#-------------------------------------- Braginskii Viscosity
BRAGINSKII_VISCOSITY
#BRAGINSKII_VISCOSITY_CONSTANT_TIMESTEP
BRAGINSKII_VISCOSITY_GLOBAL_TIMESTEP
BRAGINSKII_RKL2_SUPER_TIME_STEPPING

