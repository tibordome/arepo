#include <stdio.h>
#include "arepoconfig.h"
#ifdef HAVE_HDF5
#include <hdf5.h>
#include "hdf5_util.h"

void write_compile_time_options_in_hdf5(hid_t handle)
{
hid_t hdf5_dataspace, hdf5_attribute;
double val;
hid_t atype = H5Tcopy(H5T_C_S1);
H5Tset_size(atype, 1);
hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "PMGRID", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 128;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "PMGRID");
my_H5Aclose(hdf5_attribute, "PMGRID");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "NTYPES", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 6;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "NTYPES");
my_H5Aclose(hdf5_attribute, "NTYPES");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "MHD", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "MHD");
my_H5Aclose(hdf5_attribute, "MHD");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "MHD_POWELL", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "MHD_POWELL");
my_H5Aclose(hdf5_attribute, "MHD_POWELL");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "MHD_SEEDFIELD", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "MHD_SEEDFIELD");
my_H5Aclose(hdf5_attribute, "MHD_SEEDFIELD");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "MHD_POWELL_LIMIT_TIMESTEP", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "MHD_POWELL_LIMIT_TIMESTEP");
my_H5Aclose(hdf5_attribute, "MHD_POWELL_LIMIT_TIMESTEP");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "COOLING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "COOLING");
my_H5Aclose(hdf5_attribute, "COOLING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "UVB_SELF_SHIELDING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "UVB_SELF_SHIELDING");
my_H5Aclose(hdf5_attribute, "UVB_SELF_SHIELDING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "USE_SFR", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "USE_SFR");
my_H5Aclose(hdf5_attribute, "USE_SFR");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "VORONOI", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "VORONOI");
my_H5Aclose(hdf5_attribute, "VORONOI");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "RIEMANN_HLLD", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "RIEMANN_HLLD");
my_H5Aclose(hdf5_attribute, "RIEMANN_HLLD");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "REGULARIZE_MESH_CM_DRIFT", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "REGULARIZE_MESH_CM_DRIFT");
my_H5Aclose(hdf5_attribute, "REGULARIZE_MESH_CM_DRIFT");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED");
my_H5Aclose(hdf5_attribute, "REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "REGULARIZE_MESH_FACE_ANGLE", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "REGULARIZE_MESH_FACE_ANGLE");
my_H5Aclose(hdf5_attribute, "REGULARIZE_MESH_FACE_ANGLE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "TREE_BASED_TIMESTEPS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "TREE_BASED_TIMESTEPS");
my_H5Aclose(hdf5_attribute, "TREE_BASED_TIMESTEPS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "REFINEMENT_SPLIT_CELLS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "REFINEMENT_SPLIT_CELLS");
my_H5Aclose(hdf5_attribute, "REFINEMENT_SPLIT_CELLS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "REFINEMENT_MERGE_CELLS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "REFINEMENT_MERGE_CELLS");
my_H5Aclose(hdf5_attribute, "REFINEMENT_MERGE_CELLS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SELFGRAVITY", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SELFGRAVITY");
my_H5Aclose(hdf5_attribute, "SELFGRAVITY");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "HIERARCHICAL_GRAVITY", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "HIERARCHICAL_GRAVITY");
my_H5Aclose(hdf5_attribute, "HIERARCHICAL_GRAVITY");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "CELL_CENTER_GRAVITY", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "CELL_CENTER_GRAVITY");
my_H5Aclose(hdf5_attribute, "CELL_CENTER_GRAVITY");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "ALLOW_DIRECT_SUMMATION", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "ALLOW_DIRECT_SUMMATION");
my_H5Aclose(hdf5_attribute, "ALLOW_DIRECT_SUMMATION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "DIRECT_SUMMATION_THRESHOLD", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 5000;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "DIRECT_SUMMATION_THRESHOLD");
my_H5Aclose(hdf5_attribute, "DIRECT_SUMMATION_THRESHOLD");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "ENFORCE_JEANS_STABILITY_OF_CELLS_EEOS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "ENFORCE_JEANS_STABILITY_OF_CELLS_EEOS");
my_H5Aclose(hdf5_attribute, "ENFORCE_JEANS_STABILITY_OF_CELLS_EEOS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "ENFORCE_JEANS_STABILITY_OF_CELLS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "ENFORCE_JEANS_STABILITY_OF_CELLS");
my_H5Aclose(hdf5_attribute, "ENFORCE_JEANS_STABILITY_OF_CELLS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "EVALPOTENTIAL", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "EVALPOTENTIAL");
my_H5Aclose(hdf5_attribute, "EVALPOTENTIAL");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "NSOFTTYPES", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 4;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "NSOFTTYPES");
my_H5Aclose(hdf5_attribute, "NSOFTTYPES");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "MULTIPLE_NODE_SOFTENING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "MULTIPLE_NODE_SOFTENING");
my_H5Aclose(hdf5_attribute, "MULTIPLE_NODE_SOFTENING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "INDIVIDUAL_GRAVITY_SOFTENING", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 32;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "INDIVIDUAL_GRAVITY_SOFTENING");
my_H5Aclose(hdf5_attribute, "INDIVIDUAL_GRAVITY_SOFTENING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "ADAPTIVE_HYDRO_SOFTENING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "ADAPTIVE_HYDRO_SOFTENING");
my_H5Aclose(hdf5_attribute, "ADAPTIVE_HYDRO_SOFTENING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "RCUT", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 5.0;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "RCUT");
my_H5Aclose(hdf5_attribute, "RCUT");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "FFT_COLUMN_BASED", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "FFT_COLUMN_BASED");
my_H5Aclose(hdf5_attribute, "FFT_COLUMN_BASED");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "CHUNKING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "CHUNKING");
my_H5Aclose(hdf5_attribute, "CHUNKING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "DOUBLEPRECISION", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 1;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "DOUBLEPRECISION");
my_H5Aclose(hdf5_attribute, "DOUBLEPRECISION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "DOUBLEPRECISION_FFTW", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "DOUBLEPRECISION_FFTW");
my_H5Aclose(hdf5_attribute, "DOUBLEPRECISION_FFTW");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "OUTPUT_IN_DOUBLEPRECISION", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "OUTPUT_IN_DOUBLEPRECISION");
my_H5Aclose(hdf5_attribute, "OUTPUT_IN_DOUBLEPRECISION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "OUTPUT_COORDINATES_IN_DOUBLEPRECISION", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "OUTPUT_COORDINATES_IN_DOUBLEPRECISION");
my_H5Aclose(hdf5_attribute, "OUTPUT_COORDINATES_IN_DOUBLEPRECISION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "NGB_TREE_DOUBLEPRECISION", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "NGB_TREE_DOUBLEPRECISION");
my_H5Aclose(hdf5_attribute, "NGB_TREE_DOUBLEPRECISION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "FOF", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "FOF");
my_H5Aclose(hdf5_attribute, "FOF");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "FOF_PRIMARY_LINK_TYPES", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 2;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "FOF_PRIMARY_LINK_TYPES");
my_H5Aclose(hdf5_attribute, "FOF_PRIMARY_LINK_TYPES");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "FOF_SECONDARY_LINK_TYPES", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 1+16+32;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "FOF_SECONDARY_LINK_TYPES");
my_H5Aclose(hdf5_attribute, "FOF_SECONDARY_LINK_TYPES");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SUBFIND", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SUBFIND");
my_H5Aclose(hdf5_attribute, "SUBFIND");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SAVE_HSML_IN_SNAPSHOT", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SAVE_HSML_IN_SNAPSHOT");
my_H5Aclose(hdf5_attribute, "SAVE_HSML_IN_SNAPSHOT");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SUBFIND_CALC_MORE", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SUBFIND_CALC_MORE");
my_H5Aclose(hdf5_attribute, "SUBFIND_CALC_MORE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SUBFIND_EXTENDED_PROPERTIES", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SUBFIND_EXTENDED_PROPERTIES");
my_H5Aclose(hdf5_attribute, "SUBFIND_EXTENDED_PROPERTIES");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SOFTEREQS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SOFTEREQS");
my_H5Aclose(hdf5_attribute, "SOFTEREQS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "PROCESS_TIMES_OF_OUTPUTLIST", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "PROCESS_TIMES_OF_OUTPUTLIST");
my_H5Aclose(hdf5_attribute, "PROCESS_TIMES_OF_OUTPUTLIST");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "VORONOI_DYNAMIC_UPDATE", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "VORONOI_DYNAMIC_UPDATE");
my_H5Aclose(hdf5_attribute, "VORONOI_DYNAMIC_UPDATE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "ENLARGE_DYNAMIC_RANGE_IN_TIME", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "ENLARGE_DYNAMIC_RANGE_IN_TIME");
my_H5Aclose(hdf5_attribute, "ENLARGE_DYNAMIC_RANGE_IN_TIME");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES");
my_H5Aclose(hdf5_attribute, "FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "REDUCE_FLUSH", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "REDUCE_FLUSH");
my_H5Aclose(hdf5_attribute, "REDUCE_FLUSH");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "OUTPUT_CPU_CSV", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "OUTPUT_CPU_CSV");
my_H5Aclose(hdf5_attribute, "OUTPUT_CPU_CSV");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "OUTPUT_CENTER_OF_MASS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "OUTPUT_CENTER_OF_MASS");
my_H5Aclose(hdf5_attribute, "OUTPUT_CENTER_OF_MASS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "OUTPUTPOTENTIAL", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "OUTPUTPOTENTIAL");
my_H5Aclose(hdf5_attribute, "OUTPUTPOTENTIAL");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "HAVE_HDF5", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "HAVE_HDF5");
my_H5Aclose(hdf5_attribute, "HAVE_HDF5");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "HOST_MEMORY_REPORTING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "HOST_MEMORY_REPORTING");
my_H5Aclose(hdf5_attribute, "HOST_MEMORY_REPORTING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM");
my_H5Aclose(hdf5_attribute, "GFM");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_STELLAR_EVOLUTION", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 0;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "GFM_STELLAR_EVOLUTION");
my_H5Aclose(hdf5_attribute, "GFM_STELLAR_EVOLUTION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_CONST_IMF", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 0;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "GFM_CONST_IMF");
my_H5Aclose(hdf5_attribute, "GFM_CONST_IMF");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_PREENRICH", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_PREENRICH");
my_H5Aclose(hdf5_attribute, "GFM_PREENRICH");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_WINDS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_WINDS");
my_H5Aclose(hdf5_attribute, "GFM_WINDS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_COOLING_METAL", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_COOLING_METAL");
my_H5Aclose(hdf5_attribute, "GFM_COOLING_METAL");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_AGN_RADIATION", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_AGN_RADIATION");
my_H5Aclose(hdf5_attribute, "GFM_AGN_RADIATION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_STELLAR_PHOTOMETRICS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_STELLAR_PHOTOMETRICS");
my_H5Aclose(hdf5_attribute, "GFM_STELLAR_PHOTOMETRICS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_OUTPUT_MASK", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 1+2+4+8+16+32+64+256;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "GFM_OUTPUT_MASK");
my_H5Aclose(hdf5_attribute, "GFM_OUTPUT_MASK");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_NORMALIZED_METAL_ADVECTION", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_NORMALIZED_METAL_ADVECTION");
my_H5Aclose(hdf5_attribute, "GFM_NORMALIZED_METAL_ADVECTION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_OUTPUT_BIRTH_POS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_OUTPUT_BIRTH_POS");
my_H5Aclose(hdf5_attribute, "GFM_OUTPUT_BIRTH_POS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_CHEMTAGS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_CHEMTAGS");
my_H5Aclose(hdf5_attribute, "GFM_CHEMTAGS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_SPLITFE", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_SPLITFE");
my_H5Aclose(hdf5_attribute, "GFM_SPLITFE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_RPROCESS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_RPROCESS");
my_H5Aclose(hdf5_attribute, "GFM_RPROCESS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_DISCRETE_ENRICHMENT", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_DISCRETE_ENRICHMENT");
my_H5Aclose(hdf5_attribute, "GFM_DISCRETE_ENRICHMENT");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

my_H5Tclose(atype);
}
#endif
