#include <stdio.h>
#include "arepoconfig.h"
#ifdef HAVE_HDF5
#include <hdf5.h>

hid_t my_H5Acreate(hid_t loc_id, const char *attr_name, hid_t type_id, hid_t space_id, hid_t acpl_id);
hid_t my_H5Screate(H5S_class_t type);
herr_t my_H5Aclose(hid_t attr_id, const char *attr_name);
herr_t my_H5Awrite(hid_t attr_id, hid_t mem_type_id, const void *buf, const char *attr_name);
herr_t my_H5Sclose(hid_t dataspace_id, H5S_class_t type);

herr_t my_H5Tclose(hid_t type_id);

void write_compile_time_options_in_hdf5(hid_t handle)
{
hid_t hdf5_dataspace, hdf5_attribute;
double val;
hid_t atype = H5Tcopy(H5T_C_S1);
H5Tset_size(atype, 1);
hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "PMGRID" , H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 512;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "PMGRID");
my_H5Aclose(hdf5_attribute, "PMGRID");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "INPUT_IN_DOUBLEPRECISION" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "INPUT_IN_DOUBLEPRECISION");
my_H5Aclose(hdf5_attribute, "INPUT_IN_DOUBLEPRECISION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "BECDM" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "BECDM");
my_H5Aclose(hdf5_attribute, "BECDM");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "BECDM_INPUT_PHASE_AS_VX" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "BECDM_INPUT_PHASE_AS_VX");
my_H5Aclose(hdf5_attribute, "BECDM_INPUT_PHASE_AS_VX");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "EVALPOTENTIAL" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "EVALPOTENTIAL");
my_H5Aclose(hdf5_attribute, "EVALPOTENTIAL");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "NTYPES" , H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 6;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "NTYPES");
my_H5Aclose(hdf5_attribute, "NTYPES");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "PERIODIC" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "PERIODIC");
my_H5Aclose(hdf5_attribute, "PERIODIC");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "VORONOI" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "VORONOI");
my_H5Aclose(hdf5_attribute, "VORONOI");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "FORCE_EQUAL_TIMESTEPS" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "FORCE_EQUAL_TIMESTEPS");
my_H5Aclose(hdf5_attribute, "FORCE_EQUAL_TIMESTEPS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SELFGRAVITY" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SELFGRAVITY");
my_H5Aclose(hdf5_attribute, "SELFGRAVITY");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "RCUT" , H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 5.0;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "RCUT");
my_H5Aclose(hdf5_attribute, "RCUT");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "CHUNKING" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "CHUNKING");
my_H5Aclose(hdf5_attribute, "CHUNKING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "DOUBLEPRECISION" , H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 1;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "DOUBLEPRECISION");
my_H5Aclose(hdf5_attribute, "DOUBLEPRECISION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "DOUBLEPRECISION_FFTW" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "DOUBLEPRECISION_FFTW");
my_H5Aclose(hdf5_attribute, "DOUBLEPRECISION_FFTW");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "NGB_TREE_DOUBLEPRECISION" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "NGB_TREE_DOUBLEPRECISION");
my_H5Aclose(hdf5_attribute, "NGB_TREE_DOUBLEPRECISION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "FOF" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "FOF");
my_H5Aclose(hdf5_attribute, "FOF");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "FOF_PRIMARY_LINK_TYPES" , H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 2;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "FOF_PRIMARY_LINK_TYPES");
my_H5Aclose(hdf5_attribute, "FOF_PRIMARY_LINK_TYPES");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "FOF_SECONDARY_LINK_TYPES" , H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 1+16+32;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "FOF_SECONDARY_LINK_TYPES");
my_H5Aclose(hdf5_attribute, "FOF_SECONDARY_LINK_TYPES");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SUBFIND" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SUBFIND");
my_H5Aclose(hdf5_attribute, "SUBFIND");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SAVE_HSML_IN_SNAPSHOT" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SAVE_HSML_IN_SNAPSHOT");
my_H5Aclose(hdf5_attribute, "SAVE_HSML_IN_SNAPSHOT");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SUBFIND_CALC_MORE" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SUBFIND_CALC_MORE");
my_H5Aclose(hdf5_attribute, "SUBFIND_CALC_MORE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SUBFIND_EXTENDED_PROPERTIES" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SUBFIND_EXTENDED_PROPERTIES");
my_H5Aclose(hdf5_attribute, "SUBFIND_EXTENDED_PROPERTIES");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "PROCESS_TIMES_OF_OUTPUTLIST" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "PROCESS_TIMES_OF_OUTPUTLIST");
my_H5Aclose(hdf5_attribute, "PROCESS_TIMES_OF_OUTPUTLIST");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "VORONOI_DYNAMIC_UPDATE" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "VORONOI_DYNAMIC_UPDATE");
my_H5Aclose(hdf5_attribute, "VORONOI_DYNAMIC_UPDATE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "ENLARGE_DYNAMIC_RANGE_IN_TIME" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "ENLARGE_DYNAMIC_RANGE_IN_TIME");
my_H5Aclose(hdf5_attribute, "ENLARGE_DYNAMIC_RANGE_IN_TIME");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "REDUCE_FLUSH" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "REDUCE_FLUSH");
my_H5Aclose(hdf5_attribute, "REDUCE_FLUSH");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "OUTPUT_CPU_CSV" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "OUTPUT_CPU_CSV");
my_H5Aclose(hdf5_attribute, "OUTPUT_CPU_CSV");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "HAVE_HDF5" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "HAVE_HDF5");
my_H5Aclose(hdf5_attribute, "HAVE_HDF5");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "HOST_MEMORY_REPORTING" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "HOST_MEMORY_REPORTING");
my_H5Aclose(hdf5_attribute, "HOST_MEMORY_REPORTING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

my_H5Tclose(atype);
}
#endif
