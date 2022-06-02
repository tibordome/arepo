/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/hdf5_util.c
 * \date        09/2014
 * \author      Federico Marinacci
 * \brief Contains the wrapper functions to the HDF5 library functions.
 * \details The wrapper functions explicitly check for error conditions
 *  and terminate the run if such conditions occur.  The HDF5 error handler
 *  is disabled in case of termination not to repeat the error message of the
 *  handler again at the program exit.
 *
 * \par Major modifications and contributions:
 *
 * - 11.09.2014 first file version
 */

#include "allvars.h"
#include "proto.h"

#ifdef HAVE_HDF5
#ifndef HDF5UTIL_H
#define HDF5UTIL_H

#include <hdf5.h>

/*! \brief Wraps creating a file to give a nice error message.
 *
 *  Calls H5Fcreate.
 *
 *  \param[in] fname File name.
 *  \param[in] flags Flags handed to H5Fcreate.
 *  \param[in] fcpl_id File creation property list identifier, used when
 *             modifying default file meta-data. Use H5P_DEFAULT to specify
 *             default file creation properties.
 *  \param[in] fapl_id File access property list identifier. If parallel file
 *             access is desired, this is a collective call according to the
 *             communicator stored in the fapl_id. Use H5P_DEFAULT for default
 *             file access properties.
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return File identifier.
 */
hid_t my_H5Fcreate_fullinfo(const char *fname, unsigned int flags, hid_t fcpl_id, hid_t fapl_id, const char *func, const char *file,
                            const int line)
{
  hid_t file_id = H5Fcreate(fname, flags, fcpl_id, fapl_id);

#ifndef TOLERATE_WRITE_ERROR
  if(file_id < 0)
    terminate(
        "On Task %d, error detected in HDF5: unable to create file %s\n"
        "  at %s()/%s/line %d.",
        ThisTask, fname, func, file, line);
#endif

  return file_id;
}

/*! \brief Wraps creating a group to give a nice error message.
 *
 *  Calls H5Gcreate.
 *
 *  \param[in] loc_id File or group identifier.
 *  \param[in] groupname Absolute or relative name of the o new group.
 *  \param[in] size_hint Optional parameter indicating the number of bytes to
 *             reserve for the names that will appear in the group. A
 *             conservative estimate could result in multiple system-level
 *             I/O requests to read the group name heap; a liberal estimate
 *             could result in a single large I/O request even when the group
 *             has just a few names. HDF5 stores each name with a null
 *             terminator.
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Group identifier.
 */
hid_t my_H5Gcreate_fullinfo(hid_t loc_id, const char *groupname, size_t size_hint, const char *func, const char *file, const int line)
{
  hid_t group_id = H5Gcreate(loc_id, groupname, size_hint);

#ifndef TOLERATE_WRITE_ERROR
  if(group_id < 0)
    terminate(
        "On Task %d, error detected in HDF5: unable to create group %s\n"
        "  at %s()/%s/line %d.",
        ThisTask, groupname, func, file, line);
#endif

  return group_id;
}

/*! \brief Wraps creating a dataset to give a nice error message.
 *
 *  Calls H5Dcreate.
 *
 *  \param[in] loc_id Identifier of the file or group within which to create
 *             the dataset.
 *  \param[in] datasetname The name of the dataset to create.
 *  \param[in] type_id Identifier of the datatype to use when creating the
 *             dataset.
 *  \param[in] space_id Identifier of the dataspace to use when creating the
 *             dataset.
 *  \param[in] dcpl_id Dataset creation property list identifier.
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Dataset identifier.
 */
hid_t my_H5Dcreate_fullinfo(hid_t loc_id, const char *datasetname, hid_t type_id, hid_t space_id, hid_t dcpl_id, const char *func,
                            const char *file, const int line)
{
  hid_t dataset_id = H5Dcreate(loc_id, datasetname, type_id, space_id, dcpl_id);

#ifndef TOLERATE_WRITE_ERROR
  if(dataset_id < 0)
    terminate(
        "On Task %d, Error detected in HDF5: unable to create dataset %s\n"
        "  at %s()/%s/line %d.",
        ThisTask, datasetname, func, file, line);
#endif

  return dataset_id;
}

/*! \brief Wraps writing a dataset to give a nice error message.
 *
 *  Calls H5Dwrite.
 *
 *  \param[in] dataset_id Identifier of the dataset to write to.
 *  \param[in] mem_type_id Identifier of the memory datatype.
 *  \param[in] mem_space_id Identifier of the memory dataspace.
 *  \param[in] file_space_id Identifier of the dataset's dataspace in the file.
 *  \param[in] xfer_plist_id  Identifier of a transfer property list for this
 *             I/O operation.
 *  \param[in] buf Buffer with data to be written to the file.
 *  \param[in] datasetname Name of dataset (for error message only)
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Status of write operation.
 */
herr_t my_H5Dwrite_fullinfo(hid_t dataset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t xfer_plist_id,
                            const void *buf, const char *datasetname, const char *func, const char *file, const int line)
{
#ifdef TOLERATE_WRITE_ERROR
  if(WriteErrorFlag)
    return 0;
#endif

  herr_t status = H5Dwrite(dataset_id, mem_type_id, mem_space_id, file_space_id, xfer_plist_id, buf);

#ifndef TOLERATE_WRITE_ERROR
  if(status < 0)
    terminate(
        "On Task %d, error detected in HDF5: unable to write dataset %s\n"
        "  at %s()/%s/line %d.",
        ThisTask, datasetname, func, file, line);
#endif

  return status;
}

/*! \brief Wraps creating an attribute to give a nice error message.
 *
 *  \param[in] loc_id Identifier for the object to which the attribute is to be
 *             attached. May be any HDF5 object identifier (group, dataset, or
 *             committed datatype) or an HDF5 file identifier; if loc_id is a
 *             file identifer, the attribute will be attached to that files
 *             root group.
 *  \param[in] attr_name Name of attribute to create.
 *  \param[in] type_id Identifier of datatype for attribute.
 *  \param[in] space_id Identifier of dataspace for attribute.
 *  \param[in] acpl_id Identifier of creation property list (specify
 *             H5P_DEFAULT).
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Attribute identifier.
 */
hid_t my_H5Acreate_fullinfo(hid_t loc_id, const char *attr_name, hid_t type_id, hid_t space_id, hid_t acpl_id, const char *func,
                            const char *file, const int line)
{
  hid_t attribute_id = H5Acreate(loc_id, attr_name, type_id, space_id, acpl_id);

#ifndef TOLERATE_WRITE_ERROR
  if(attribute_id < 0)
    terminate(
        "On Task %d, error detected in HDF5: unable to create attribute %s\n"
        "  at %s()/%s/line %d.",
        ThisTask, attr_name, func, file, line);
#endif

  return attribute_id;
}

/*! \brief Wraps writing an attribute to give a nice error message.
 *
 *  \param[in] attr_id Identifier of an attribute to write.
 *  \param[in] mem_type_id Identifier of the attribute datatype (in memory).
 *  \param[in] buf Data to be written.
 *  \param[in] attr_name Name of attribute (for error message only).
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return status (non-negative if successful).
 */
herr_t my_H5Awrite_fullinfo(hid_t attr_id, hid_t mem_type_id, const void *buf, const char *attr_name, const char *func,
                            const char *file, const int line)
{
#ifdef TOLERATE_WRITE_ERROR
  if(WriteErrorFlag)
    return 0;
#endif

  herr_t status = H5Awrite(attr_id, mem_type_id, buf);

#ifndef TOLERATE_WRITE_ERROR
  if(status < 0)
    terminate(
        "On Task %d, error detected in HDF5: unable to write attribute %s\n"
        "  at %s()/%s/line %d.",
        ThisTask, attr_name, func, file, line);
#endif

  return status;
}

/*! \brief Wraps creating a dataspace to give a nice error message.
 *
 *  \param[in] type Type of dataspace to be created.
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Dataspace identifier if successful.
 */
hid_t my_H5Screate_fullinfo(H5S_class_t type, const char *func, const char *file, const int line)
{
  hid_t dataspace_id = H5Screate(type);

#ifndef TOLERATE_WRITE_ERROR
  if(dataspace_id < 0)
    {
      switch(type)
        {
          case H5S_SCALAR:
            terminate(
                "On Task %d, error detected in HDF5: unable to create a scalar dataspace\n"
                "  at %s()/%s/line %d.",
                ThisTask, func, file, line);
            break;
          case H5S_SIMPLE:
            terminate(
                "On Task %d, error detected in HDF5: unable to create a simple dataspace\n"
                "  at %s()/%s/line %d.",
                ThisTask, func, file, line);
            break;
          default:
            terminate(
                "On Task %d, error detected in HDF5: unknown dataspace type\n"
                "  at %s()/%s/line %d.",
                ThisTask, func, file, line);
            break;
        }
    }
#endif

  return dataspace_id;
}

/*! \brief Wraps creating a simple dataspace to give a nice error message.
 *
 *  \param[in] rank Number of dimensions of dataspace.
 *  \param[in] current_dims Array specifying the size of each dimension.
 *  \param[in] maximum_dims Array specifying the maximum size of each
 *             dimension.
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Dataspace identifier if successful.
 */
hid_t my_H5Screate_simple_fullinfo(int rank, const hsize_t *current_dims, const hsize_t *maximum_dims, const char *func,
                                   const char *file, const int line)
{
  hid_t dataspace_id = H5Screate_simple(rank, current_dims, maximum_dims);

#ifndef TOLERATE_WRITE_ERROR
  if(dataspace_id < 0)
    terminate(
        "On Task %d, error detected in HDF5: unable to create a simple dataspace\n"
        "  at %s()/%s/line %d.",
        ThisTask, func, file, line);
#endif

  return dataspace_id;
}

/*! \brief Wraps opening a file to give a nice error message.
 *
 *  \param[in] fname Name of the file to be opened.
 *  \param[in] flags File access flags. Allowable values are:
 *             H5F_ACC_RDWR -- Allow read and write access to file.
 *             H5F_ACC_RDONLY -- Allow read-only access to file.
 *  \param[in] fapl_id Identifier for the file access properties list. If
 *             parallel file access is desired, this is a collective call
 *             according to the communicator stored in the fapl_id. Use
 *             H5P_DEFAULT for default file access properties.
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return File identifier if successful.
 */
hid_t my_H5Fopen_fullinfo(const char *fname, unsigned int flags, hid_t fapl_id, const char *func, const char *file, const int line)
{
  hid_t file_id = H5Fopen(fname, flags, fapl_id);

  if(file_id < 0)
    terminate(
        "On Task %d, error detected in HDF5: unable to open file %s\n"
        "  at %s()/%s/line %d.",
        ThisTask, fname, func, file, line);

  return file_id;
}

/*! \brief Wraps opening a group to give a nice error message.
 *
 *  \param[in] loc_id File or group identifier within which the group is to be
 *             opened.
 *  \param[in] groupname Name of group.
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Valid group identifier if successful.
 */
hid_t my_H5Gopen_fullinfo(hid_t loc_id, const char *groupname, const char *func, const char *file, const int line)
{
  hid_t group = H5Gopen(loc_id, groupname);

#ifndef TOLERATE_WRITE_ERROR
  if(group < 0)
    terminate(
        "On Task %d, error detected in HDF5: unable to open group %s\n"
        "  at %s()/%s/line %d.",
        ThisTask, groupname, func, file, line);
#endif

  return group;
}

/*! \brief Wraps opening a dataset to give a nice error message.
 *
 *  \param[in] file_id Identifier of the file or group within which the
 *             dataset to be accessed will be found.
 *  \param[in] datasetname Name of the dataset to access.
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Dataset identifier if successful.
 */
hid_t my_H5Dopen_fullinfo(hid_t file_id, const char *datasetname, const char *func, const char *file, const int line)
{
  hid_t dataset = H5Dopen(file_id, datasetname);

#ifndef TOLERATE_WRITE_ERROR
  if(dataset < 0)
    terminate(
        "On Task %d, error detected in HDF5: unable to open dataset %s\n"
        "  at %s()/%s/line %d.",
        ThisTask, datasetname, func, file, line);
#endif

  return dataset;
}

/*! \brief Wraps opening a dataset.
 *
 *  In contrast to #my_H5Dopen(), if the dataset does not exist or opening
 *  fails for other reasons, it does not terminate the run. This is useful
 *  while reading an IC file because in that case a non-exisitng dataset is
 *  put to zero (see also read_ic.c).
 *
 *  \param[in] file_id file_id Identifier of the file or group within which the
 *             dataset to be accessed will be found.
 *  \param[in] datasetname Name of the dataset to access.
 *
 *  \return Dataset identifier if successful; otherwise negative value.
 */
hid_t my_H5Dopen_if_existing(hid_t file_id, const char *datasetname)
{
  /* save error handler and disable it */
  H5E_auto_t errfunc;
  void *client_data;
  H5Eget_auto(&errfunc, &client_data);
  H5Eset_auto(NULL, NULL);

  hid_t dataset = H5Dopen(file_id, datasetname);

  /* reset error handler */
  H5Eset_auto(errfunc, client_data);

  return dataset;
}

/*! \brief Wraps opening an attribute to give a nice error message.
 *
 *  \param[in] loc_id  Identifier of a group, dataset, or named datatype that
 *             attribute is attached to.
 *  \param[in] attr_name Attribute name.
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Returns attribute identifier if successful.
 */
hid_t my_H5Aopen_name_fullinfo(hid_t loc_id, const char *attr_name, const char *func, const char *file, const int line)
{
  hid_t attribute_id = H5Aopen_name(loc_id, attr_name);

#ifndef TOLERATE_WRITE_ERROR
  if(attribute_id < 0)
    terminate(
        "On Task %d, error detected in HDF5: unable to open attribute %s\n"
        "  at %s()/%s/line %d.",
        ThisTask, attr_name, func, file, line);
#endif

  return attribute_id;
}

/*! \brief Wraps reading a dataset to give a nice error message.
 *
 *  \param[in] dataset_id Identifier of the dataset read from.
 *  \param[in] mem_type_id Identifier of the memory datatype.
 *  \param[in] mem_space_id Identifier of the memory dataspace.
 *  \param[in] file_space_id Identifier of the dataset's dataspace in the file.
 *  \param[in] xfer_plist_id Identifier of a transfer property list for this
 *             I/O operation.
 *  \param[out] buf Buffer to receive data read from file.
 *  \param[in] datasetname Name of dataset (only for error message).
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Returns a non-negative value if successful.
 */
herr_t my_H5Dread_fullinfo(hid_t dataset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t xfer_plist_id,
                           void *buf, const char *datasetname, const char *func, const char *file, const int line)
{
  herr_t status = H5Dread(dataset_id, mem_type_id, mem_space_id, file_space_id, xfer_plist_id, buf);
  if(status < 0)
    terminate(
        "On Task %d, error detected in HDF5: unable to read dataset %s\n"
        "  at %s()/%s/line %d.",
        ThisTask, datasetname, func, file, line);
  return status;
}

/*! \brief Wraps makeing a copy of the dataspace to give a nice error message.
 *
 *  \param[in] dataset_id Identifier of the dataset to query.
 *  \param[in] datasetname Name of the dataset (for error message only).
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Dataspace identifier if successful.
 */
hid_t my_H5Dget_space_fullinfo(hid_t dataset_id, const char *datasetname, const char *func, const char *file, const int line)
{
  hid_t status = H5Dget_space(dataset_id);

#ifndef TOLERATE_WRITE_ERROR
  if(status < 0)
    terminate(
        "On Task %d, error detected in HDF5: unable to determine space for dataset %s\n"
        "  at %s()/%s/line %d.",
        ThisTask, datasetname, func, file, line);
#endif

  return status;
}

/*! \brief Wraps reading an attribute to give a nice error message
 *
 *  \param[in] attr_id Identifier of an attribute to read.
 *  \param[in] mem_type_id Identifier of the attribute datatype (in memory).
 *  \param[out] buf Buffer for data to be read.
 *  \param[in] attr_name Name of the attribute.
 *  \param[in] size Size of the attribute.
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Non-negative value if successful.
 */
herr_t my_H5Aread_fullinfo(hid_t attr_id, hid_t mem_type_id, void *buf, const char *attr_name, hssize_t size, const char *func,
                           const char *file, const int line)
{
  hid_t hdf5_space   = H5Aget_space(attr_id);
  hssize_t attr_size = H5Sget_simple_extent_npoints(hdf5_space);
  H5Sclose(hdf5_space);

  if(attr_size != size)
    terminate(
        "On Task %d, error detected in HDF5: mismatch in size for attribute %s, expected size = %lld, actual attribute size = "
        "%lld\n"
        "  at %s()/%s/line %d.",
        ThisTask, attr_name, size, attr_size, func, file, line);

  herr_t status = H5Aread(attr_id, mem_type_id, buf);
  if(status < 0)
    terminate(
        "On Task %d, error detected in HDF5: unable to read attribute %s\n"
        "  at %s()/%s/line %d.",
        ThisTask, attr_name, func, file, line);
  return status;
}

/*! \brief Wraps reseting the size of an existing dataspace to give a nice
 *         error message.
 *
 *  \param[in] space_id Dataspace identifier.
 *  \param[in] rank Rank, or dimensionality, of the dataspace.
 *  \param[in] current_size Array containing current size of dataspace.
 *  \param[in] maximum_size Array containing maximum size of dataspace.
 *  \param[in] attr_name Name of attribute (only for error message).
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Non-negative value if successful.
 */
herr_t my_H5Sset_extent_simple_fullinfo(hid_t space_id, int rank, const hsize_t *current_size, const hsize_t *maximum_size,
                                        const char *attr_name, const char *func, const char *file, const int line)
{
  herr_t status = H5Sset_extent_simple(space_id, rank, current_size, maximum_size);

#ifndef TOLERATE_WRITE_ERROR
  if(status < 0)
    terminate(
        "On Task %d, error detected in HDF5: unable to set extent for attribute %s\n"
        "  at %s()/%s/line %d.",
        ThisTask, attr_name, func, file, line);
#endif

  return status;
}

/*! \brief Wraps closing an attribute to give a nice error message.
 *
 *  \param[in] attr_id Attribute to release access to.
 *  \param[in] attr_name Name of the attribute (for error message only).
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Non-negative value if successful.
 */
herr_t my_H5Aclose_fullinfo(hid_t attr_id, const char *attr_name, const char *func, const char *file, const int line)
{
  herr_t status = H5Aclose(attr_id);

#ifndef TOLERATE_WRITE_ERROR
  if(status < 0)
    terminate(
        "On Task %d, error detected in HDF5: unable to close attribute %s\n"
        "  at %s()/%s/line %d.",
        ThisTask, attr_name, func, file, line);
#endif

  return status;
}

/*! \brief Wraps closing a dataset to give a nice error message.
 *
 *  \param[in] dataset_id Identifier of the dataset to close access to.
 *  \param[in] datasetname Name of the dataset (for error message only).
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Non-negative value if successful.
 */
herr_t my_H5Dclose_fullinfo(hid_t dataset_id, const char *datasetname, const char *func, const char *file, const int line)
{
  herr_t status = H5Dclose(dataset_id);

#ifndef TOLERATE_WRITE_ERROR
  if(status < 0)
    terminate(
        "On Task %d, error detected in HDF5: unable to close dataset %s\n"
        "  at %s()/%s/line %d.",
        ThisTask, datasetname, func, file, line);
#endif

  return status;
}

/*! \brief Wraps closing a group to give a nice error message.
 *
 *  \param[in] group_id Group identifier to release.
 *  \param[in] groupname Name of the group (for error message only).
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Non-negative value if successful.
 */
herr_t my_H5Gclose_fullinfo(hid_t group_id, const char *groupname, const char *func, const char *file, const int line)
{
  herr_t status = H5Gclose(group_id);

#ifndef TOLERATE_WRITE_ERROR
  if(status < 0)
    terminate(
        "On Task %d, error detected in HDF5: unable to close group %s\n"
        "  at %s()/%s/line %d.",
        ThisTask, groupname, func, file, line);
#endif

  return status;
}

/*! \brief Wraps closing a file to give a nice error message.
 *
 *  \param[in] file_id Identifier of a file to terminate access to.
 *  \param[in] fname File  name (for error message only).
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Non-negative value if successful.
 */
herr_t my_H5Fclose_fullinfo(hid_t file_id, const char *fname, const char *func, const char *file, const int line)
{
  herr_t status = H5Fclose(file_id);

#ifndef TOLERATE_WRITE_ERROR
  if(status < 0)
    terminate(
        "On Task %d, error detected in HDF5: unable to close file %s\n"
        "  at %s()/%s/line %d.",
        ThisTask, fname, func, file, line);
#endif
  return status;
}

/*! \brief Wraps releasing and terminating access to a dataspace to give a nice
 *         error message.
 *
 *  \param[in] dataspace_id Identifier of dataspace to release.
 *  \param[in] type type of dataspace (simple, scalar, ...).
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Non-negative value if successful.
 */
herr_t my_H5Sclose_fullinfo(hid_t dataspace_id, H5S_class_t type, const char *func, const char *file, const int line)
{
  herr_t status = H5Sclose(dataspace_id);

#ifndef TOLERATE_WRITE_ERROR
  if(status < 0)
    {
      switch(type)
        {
          case H5S_SCALAR:
            terminate(
                "On Task %d, error detected in HDF5: unable to close a scalar dataspace\n"
                "  at %s()/%s/line %d.",
                ThisTask, func, file, line);
            break;
          case H5S_SIMPLE:
            terminate(
                "On Task %d, error detected in HDF5: unable to close a simple dataspace\n"
                "  at %s()/%s/line %d.",
                ThisTask, func, file, line);
            break;
          default:
            terminate(
                "On Task %d, error detected in HDF5: unknown dataspace type\n"
                "  at %s()/%s/line %d.",
                ThisTask, func, file, line);
            break;
        }
    }
#endif

  return status;
}

/*! \brief Wraps copying an existing datatype to give a nice error message.
 *
 *  \param[in] type_id Identifier of datatype to copy. Can be a datatype
 *             identifier, a predefined datatype (defined in H5Tpublic.h), or
 *             a dataset identifier.
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Datatype identifier if successful.
 */
hid_t my_H5Tcopy_fullinfo(hid_t type_id, const char *func, const char *file, const int line)
{
  hid_t datatype_id = H5Tcopy(type_id);
#ifndef TOLERATE_WRITE_ERROR
  if(datatype_id < 0)
    terminate(
        "On Task %d, error detected in HDF5: could not properly copy datatype\n"
        "  at %s()/%s/line %d.",
        ThisTask, func, file, line);
#endif
  return datatype_id;
}

/*! \brief Wraps closing a datatype to give a nice error message.
 *
 *  \param[in] type_id Identifier of datatype to release.
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Non-negative value if successful.
 */
herr_t my_H5Tclose_fullinfo(hid_t type_id, const char *func, const char *file, const int line)
{
  herr_t status = H5Tclose(type_id);
#ifndef TOLERATE_WRITE_ERROR
  if(status < 0)
    terminate(
        "On Task %d, error detected in HDF5: could not properly close datatype\n"
        "  at %s()/%s/line %d.",
        ThisTask, func, file, line);
#endif
  return status;
}

/*! \brief Wraps selecting a hyperslab to give a nice error message.
 *
 *  \param[in] space_id Identifier of dataspace selection to modify.
 *  \param[in] op Operation to perform on current selection.
 *  \param[in] start Offset of start of hyperslab.
 *  \param[in] stride Hyperslab stride.
 *  \param[in] count Number of blocks included in hyperslab.
 *  \param[in] block Size of block in hyperslab.
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Non-negative value if successful.
 */
herr_t my_H5Sselect_hyperslab_fullinfo(hid_t space_id, H5S_seloper_t op, const hsize_t *start, const hsize_t *stride,
                                       const hsize_t *count, const hsize_t *block, const char *func, const char *file, const int line)
{
  herr_t status = H5Sselect_hyperslab(space_id, op, start, stride, count, block);

#ifndef TOLERATE_WRITE_ERROR
  if(status < 0)
    terminate(
        "On Task %d, error detected in HDF5: could not properly select the chosen hyperslab\n"
        "  at %s()/%s/line %d.",
        ThisTask, func, file, line);
#endif
  return status;
}

/*! \brief Wraps returning the size in bytes of a given datatype to give a nice
 *         error message.
 *
 *  \param[in] datatype_id Identifier of datatype to query.
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return The size of the datatype in bytes.
 */
size_t my_H5Tget_size_fullinfo(hid_t datatype_id, const char *func, const char *file, const int line)
{
  size_t size = H5Tget_size(datatype_id);

#ifndef TOLERATE_WRITE_ERROR
  if(size == 0)
    terminate(
        "On Task %d, error detected in HDF5: unable to determine the size of the given datatype\n"
        "  at %s()/%s/line %d.",
        ThisTask, func, file, line);
#endif
  return size;
}

/*! \brief Wraps setting the size in bytes of a given datatype to give a nice
 *         error message.
 *
 *  \param[in] datatype_id Identifier of datatype for which the size is being
 *             changed.
 *  \param[in] size New datatype size in bytes or H5T_VARIABLE.
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Non-negative value if successful.
 */
herr_t my_H5Tset_size_fullinfo(hid_t datatype_id, size_t size, const char *func, const char *file, const int line)
{
  herr_t status = H5Tset_size(datatype_id, size);

#ifndef TOLERATE_WRITE_ERROR
  if(status < 0)
    terminate(
        "On Task %d, error detected in HDF5: could not properly set the size of the given datatype\n"
        "  at %s()/%s/line %d.",
        ThisTask, func, file, line);
#endif

  return status;
}

#ifdef HDF5_FILTERS
/*! \brief Wraps checking if all hdf5 filters selected for plist_id are
 *         available to give a nice error message.
 *
 *  \param[in] plist_id Dataset or group creation property list identifier.
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Positive value if all filters are available;
 *          0 if at least one filter is not currently available.
 */
htri_t my_H5Pall_filters_avail_fullinfo(hid_t plist_id, const char *func, const char *file, const int line)
{
  htri_t status = H5Pall_filters_avail(plist_id);
  if(status < 0)
    terminate(
        "On Task %d, error detected in HDF5: could not properly verify the availability of all filters\n"
        "  at %s()/%s/line %d.",
        ThisTask, func, file, line);
  return status;
}

/*! \brief Wraps creating the property list of the given property class
 *         identified by class_id to give a nice error message.
 *
 *  \param[in] The class of the property list to create.
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Property list identifier if successful.
 */
hid_t my_H5Pcreate_fullinfo(hid_t class_id, const char *func, const char *file, const int line)
{
  hid_t plist_id = H5Pcreate(class_id);
  if(plist_id < 0)
    terminate(
        "On Task %d, error detected in HDF5: could not create the property list associated to the given property class\n"
        "  at %s()/%s/line %d.",
        ThisTask, func, file, line);
  return plist_id;
}

/*! \brief Wraps closing a property list to give a nice error message.
 *
 *  \param[in] Identifier of the property list to terminate access to.
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 * \return Non-negative value if successful.
 */
herr_t my_H5Pclose_fullinfo(hid_t plist, const char *func, const char *file, const int line)
{
  herr_t status = H5Pclose(plist);
  if(status < 0)
    terminate(
        "On Task %d, error detected in HDF5: could not close the input property list\n"
        "  at %s()/%s/line %d.",
        ThisTask, func, file, line);
  return status;
}

/*! \brief Wraps setting the size of the chunks of a chunked dataset to give a
 *         nice error message.
 *
 *  \param[in] plist Dataset creation property list identifier.
 *  \param[in] ndims The number of dimensions of each chunk.
 *  \param[in] dim An array defining the size, in dataset elements, of each
 *             chunk.
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Non-negative value if successful.
 */
herr_t my_H5Pset_chunk_fullinfo(hid_t plist, int ndims, const hsize_t *dim, const char *func, const char *file, const int line)
{
  herr_t status = H5Pset_chunk(plist, ndims, dim);
  if(status < 0)
    terminate(
        "On Task %d, error detected in HDF5: could not set chunk size for the dataset\n"
        "  at %s()/%s/line %d.",
        ThisTask, func, file, line);
  return status;
}

/*! \brief Wraps setting the use of the shuffle filter to give a nice error
 *         message.
 *
 *  \param[in] plist_id Dataset creation property list identifier.
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Non-negative value if successful.
 */
herr_t my_H5Pset_shuffle_fullinfo(hid_t plist_id, const char *func, const char *file, const int line)
{
  herr_t status = H5Pset_shuffle(plist_id);
  if(status < 0)
    terminate(
        "On Task %d, error detected in HDF5: could not set the shuffle filter in the properties list\n"
        "  at %s()/%s/line %d.",
        ThisTask, func, file, line);
  return status;
}

/*! \brief Wraps setting the use of the deflate compression (gzip) to give a
 *         nice error message.
 *
 *  \param[in] plist_id Dataset or group creation property list identifier.
 *  \param[in] level Compression level.
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Non-negative value if successful.
 */
herr_t my_H5Pset_deflate_fullinfo(hid_t plist_id, unsigned level, const char *func, const char *file, const int line)
{
  herr_t status = H5Pset_deflate(plist_id, level);
  if(status < 0)
    terminate(
        "On Task %d, error detected in HDF5: could not set the deflate compression in the properties list\n"
        "  at %s()/%s/line %d.",
        ThisTask, func, file, line);
  return status;
}

/*! \brief Wraps setting the use of the Fletcher32 checksum to give a nice
 *         error message.
 *
 *  \param[in] plist_id Dataset or group creation property list identifier.
 *  \param[in] func Name of the function that has called the HDF5 routine
 *             (usually given by the __func__ macro).
 *  \param[in] file File where the function that has called the HDF5
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of the file where the HDF5 routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return Non-negative value if successful.
 */
herr_t my_H5Pset_fletcher32_fullinfo(hid_t plist_id, const char *func, const char *file, const int line)
{
  herr_t status = H5Pset_fletcher32(plist_id);
  if(status < 0)
    terminate(
        "On Task %d, error detected in HDF5: could not set the Fletcher32 checksum in the properties list\n"
        "  at %s()/%s/line %d.",
        ThisTask, func, file, line);
  return status;
}
#endif

#endif
#endif
