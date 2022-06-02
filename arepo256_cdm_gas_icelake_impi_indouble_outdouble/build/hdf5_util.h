/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/hdf5_util.h
 * \date        04/2020
 * \author      Simon May
 * \brief
 * \details
 *
 */

#ifndef HDF5_UTIL_H
#define HDF5_UTIL_H

#include <stddef.h>

#ifdef HAVE_HDF5
#include <hdf5.h>

#define my_H5Fcreate(fname, flags, fcpl_id, fapl_id) \
  my_H5Fcreate_fullinfo(fname, flags, fcpl_id, fapl_id, __func__, __FILE__, __LINE__)
#define my_H5Gcreate(loc_id, groupname, size_hint) my_H5Gcreate_fullinfo(loc_id, groupname, size_hint, __func__, __FILE__, __LINE__)
#define my_H5Dcreate(loc_id, datasetname, type_id, space_id, dcpl_id) \
  my_H5Dcreate_fullinfo(loc_id, datasetname, type_id, space_id, dcpl_id, __func__, __FILE__, __LINE__)
#define my_H5Acreate(loc_id, attr_name, type_id, space_id, acpl_id) \
  my_H5Acreate_fullinfo(loc_id, attr_name, type_id, space_id, acpl_id, __func__, __FILE__, __LINE__)
#define my_H5Screate(type) my_H5Screate_fullinfo(type, __func__, __FILE__, __LINE__)
#define my_H5Screate_simple(rank, current_dims, maximum_dims) \
  my_H5Screate_simple_fullinfo(rank, current_dims, maximum_dims, __func__, __FILE__, __LINE__)
#define my_H5Dwrite(dataset_id, mem_type_id, mem_space_id, file_space_id, xfer_plist_id, buf, datasetname)                        \
  my_H5Dwrite_fullinfo(dataset_id, mem_type_id, mem_space_id, file_space_id, xfer_plist_id, buf, datasetname, __func__, __FILE__, \
                       __LINE__)
#define my_H5Awrite(attr_id, mem_type_id, buf, attr_name) \
  my_H5Awrite_fullinfo(attr_id, mem_type_id, buf, attr_name, __func__, __FILE__, __LINE__)
#define my_H5Fopen(fname, flags, fapl_id) my_H5Fopen_fullinfo(fname, flags, fapl_id, __func__, __FILE__, __LINE__)
#define my_H5Dopen(file_id, datasetname) my_H5Dopen_fullinfo(file_id, datasetname, __func__, __FILE__, __LINE__)
#define my_H5Dread(dataset_id, mem_type_id, mem_space_id, file_space_id, xfer_plist_id, buf, datasetname)                        \
  my_H5Dread_fullinfo(dataset_id, mem_type_id, mem_space_id, file_space_id, xfer_plist_id, buf, datasetname, __func__, __FILE__, \
                      __LINE__)
#define my_H5Gopen(loc_id, groupname) my_H5Gopen_fullinfo(loc_id, groupname, __func__, __FILE__, __LINE__)
#define my_H5Aopen_name(loc_id, attr_name) my_H5Aopen_name_fullinfo(loc_id, attr_name, __func__, __FILE__, __LINE__)
#define my_H5Aread(attr_id, mem_type_id, buf, attr_name, size) \
  my_H5Aread_fullinfo(attr_id, mem_type_id, buf, attr_name, size, __func__, __FILE__, __LINE__)

#define my_H5Aclose(attr_id, attr_name) my_H5Aclose_fullinfo(attr_id, attr_name, __func__, __FILE__, __LINE__)
#define my_H5Dclose(dataset_id, datasetname) my_H5Dclose_fullinfo(dataset_id, datasetname, __func__, __FILE__, __LINE__)
#define my_H5Gclose(group_id, groupname) my_H5Gclose_fullinfo(group_id, groupname, __func__, __FILE__, __LINE__)
#define my_H5Fclose(file_id, fname) my_H5Fclose_fullinfo(file_id, fname, __func__, __FILE__, __LINE__)
#define my_H5Sclose(dataspace_id, type) my_H5Sclose_fullinfo(dataspace_id, type, __func__, __FILE__, __LINE__)

#define my_H5Tcopy(type_id) my_H5Tcopy_fullinfo(type_id, __func__, __FILE__, __LINE__)
#define my_H5Tclose(type_id) my_H5Tclose_fullinfo(type_id, __func__, __FILE__, __LINE__)

#define my_H5Sselect_hyperslab(space_id, op, start, stride, count, block) \
  my_H5Sselect_hyperslab_fullinfo(space_id, op, start, stride, count, block, __func__, __FILE__, __LINE__)
#define my_H5Tget_size(datatype_id) my_H5Tget_size_fullinfo(datatype_id, __func__, __FILE__, __LINE__)
#define my_H5Tset_size(datatype_id, size) my_H5Tset_size_fullinfo(datatype_id, size, __func__, __FILE__, __LINE__)

#define my_H5Sset_extent_simple(space_id, rank, current_size, maximum_size, attr_name) \
  my_H5Sset_extent_simple_fullinfo(space_id, rank, current_size, maximum_size, attr_name, __func__, __FILE__, __LINE__)
#define my_H5Dget_space(dataset_id, datasetname) my_H5Dget_space_fullinfo(dataset_id, datasetname, __func__, __FILE__, __LINE__)

hid_t my_H5Dopen_if_existing(hid_t file_id, const char *datasetname);

hid_t my_H5Fcreate_fullinfo(const char *fname, unsigned flags, hid_t fcpl_id, hid_t fapl_id, const char *func, const char *file,
                            int line);
hid_t my_H5Gcreate_fullinfo(hid_t loc_id, const char *groupname, size_t size_hint, const char *func, const char *file, int line);
hid_t my_H5Dcreate_fullinfo(hid_t loc_id, const char *datasetname, hid_t type_id, hid_t space_id, hid_t dcpl_id, const char *func,
                            const char *file, int line);
hid_t my_H5Acreate_fullinfo(hid_t loc_id, const char *attr_name, hid_t type_id, hid_t space_id, hid_t acpl_id, const char *func,
                            const char *file, int line);
hid_t my_H5Screate_fullinfo(H5S_class_t type, const char *func, const char *file, int line);
hid_t my_H5Screate_simple_fullinfo(int rank, const hsize_t *current_dims, const hsize_t *maximum_dims, const char *func,
                                   const char *file, int line);
herr_t my_H5Dwrite_fullinfo(hid_t dataset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t xfer_plist_id,
                            const void *buf, const char *datasetname, const char *func, const char *file, int line);
herr_t my_H5Awrite_fullinfo(hid_t attr_id, hid_t mem_type_id, const void *buf, const char *attr_name, const char *func,
                            const char *file, int line);
hid_t my_H5Fopen_fullinfo(const char *fname, unsigned int flags, hid_t fapl_id, const char *func, const char *file, int line);
hid_t my_H5Dopen_fullinfo(hid_t file_id, const char *datasetname, const char *func, const char *file, int line);
herr_t my_H5Dread_fullinfo(hid_t dataset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t xfer_plist_id,
                           void *buf, const char *datasetname, const char *func, const char *file, int line);
hid_t my_H5Gopen_fullinfo(hid_t loc_id, const char *groupname, const char *func, const char *file, int line);
hid_t my_H5Aopen_name_fullinfo(hid_t loc_id, const char *attr_name, const char *func, const char *file, int line);
herr_t my_H5Aread_fullinfo(hid_t attr_id, hid_t mem_type_id, void *buf, const char *attr_name, hssize_t size, const char *func,
                           const char *file, int line);

herr_t my_H5Aclose_fullinfo(hid_t attr_id, const char *attr_name, const char *func, const char *file, int line);
herr_t my_H5Dclose_fullinfo(hid_t dataset_id, const char *datasetname, const char *func, const char *file, int line);
herr_t my_H5Gclose_fullinfo(hid_t group_id, const char *groupname, const char *func, const char *file, int line);
herr_t my_H5Fclose_fullinfo(hid_t file_id, const char *fname, const char *func, const char *file, int line);
herr_t my_H5Sclose_fullinfo(hid_t dataspace_id, H5S_class_t type, const char *func, const char *file, int line);

hid_t my_H5Tcopy_fullinfo(hid_t type_id, const char *func, const char *file, int line);
herr_t my_H5Tclose_fullinfo(hid_t type_id, const char *func, const char *file, int line);

herr_t my_H5Sselect_hyperslab_fullinfo(hid_t space_id, H5S_seloper_t op, const hsize_t *start, const hsize_t *stride,
                                       const hsize_t *count, const hsize_t *block, const char *func, const char *file, int line);
size_t my_H5Tget_size_fullinfo(hid_t datatype_id, const char *func, const char *file, int line);
herr_t my_H5Tset_size_fullinfo(hid_t datatype_id, size_t size, const char *func, const char *file, int line);

herr_t my_H5Sset_extent_simple_fullinfo(hid_t space_id, int rank, const hsize_t *current_size, const hsize_t *maximum_size,
                                        const char *attr_name, const char *func, const char *file, int line);
hid_t my_H5Dget_space_fullinfo(hid_t dataset_id, const char *datasetname, const char *func, const char *file, int line);

#ifdef HDF5_FILTERS

#define my_H5Pall_filters_avail(plist_id) my_H5Pall_filters_avail_fullinfo(plist_id, __func__, __FILE__, __LINE__)
#define my_H5Pcreate(class_id) my_H5Pcreate_fullinfo(class_id, __func__, __FILE__, __LINE__)
#define my_H5Pclose(plist) my_H5Pclose_fullinfo(plist, __func__, __FILE__, __LINE__)
#define my_H5Pset_chunk(plist, ndims, dim) my_H5Pset_chunk_fullinfo(plist, ndims, dim, __func__, __FILE__, __LINE__)
#define my_H5Pset_shuffle(plist_id) my_H5Pset_shuffle_fullinfo(plist_id, __func__, __FILE__, __LINE__)
#define my_H5Pset_deflate(plist_id, level) my_H5Pset_deflate_fullinfo(plist_id, level, __func__, __FILE__, __LINE__)
#define my_H5Pset_fletcher32(plist_id) my_H5Pset_fletcher32_fullinfo(plist_id, __func__, __FILE__, __LINE__)

htri_t my_H5Pall_filters_avail_fullinfo(hid_t plist_id, const char *func, const char *file, int line);
hid_t my_H5Pcreate_fullinfo(hid_t class_id, const char *func, const char *file, int line);
herr_t my_H5Pclose_fullinfo(hid_t plist, const char *func, const char *file, int line);
herr_t my_H5Pset_chunk_fullinfo(hid_t plist, int ndims, const hsize_t *dim, const char *func, const char *file, int line);
herr_t my_H5Pset_shuffle_fullinfo(hid_t plist_id, const char *func, const char *file, int line);
herr_t my_H5Pset_deflate_fullinfo(hid_t plist_id, unsigned level, const char *func, const char *file, int line);
herr_t my_H5Pset_fletcher32_fullinfo(hid_t plist_id, const char *func, const char *file, int line);
#endif

#endif

#endif
