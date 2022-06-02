/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/fof/fof_io.c
 * \date        MM/YYYY
 * \author
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <gsl/gsl_math.h>
#include <inttypes.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "../allvars.h"
#include "../domain.h"
#include "../gitversion/version.h"
#include "../proto.h"
#include "../subfind/subfind.h"
#include "fof.h"

#ifdef HAVE_HDF5
#include <hdf5.h>
static void fof_subfind_write_header_attributes_in_hdf5(hid_t handle);
#endif

#ifdef FOF

/*! \brief Main routine for group output.
 *
 *  \param[in] num Index of group file (snapshot index for this output).
 *
 *  \return void
 */
void fof_save_groups(const int num)
{
  int filenr, gr, ngrps, masterTask, lastTask;
  double t0, t1;
  char buf[MAXLEN_PATH];

#ifdef FOF_STOREIDS
  fof_subfind_prepare_ID_list();
#endif

  t0 = second();

  CommBuffer = mymalloc("CommBuffer", COMMBUFFERSIZE);

  if(All.SnapFormat < SNAP_FORMAT_GADGET || All.SnapFormat > SNAP_FORMAT_HDF5)
    mpi_printf("Unsupported File-Format. All.SnapFormat=%d\n", All.SnapFormat);

#ifndef HAVE_HDF5
  if(All.SnapFormat == SNAP_FORMAT_HDF5)
    mpi_terminate("Code wasn't compiled with HDF5 support enabled!\n");
#endif

  /* assign processors to output files */
  distribute_file(All.NumFilesPerSnapshot, &filenr, &masterTask, &lastTask);

  if(All.NumFilesPerSnapshot > 1)
    {
      if(ThisTask == 0)
        {
          file_path_sprintf(buf, "%s/groups_%03d", All.OutputDir, num);
          mkdir(buf, MKDIR_MODE);
        }
      MPI_Barrier(MPI_COMM_WORLD);
    }

  if(All.NumFilesPerSnapshot > 1)
    file_path_sprintf(buf, "%s/groups_%03d/fof_tab_%03d.%d", All.OutputDir, num, num, filenr);
  else
    file_path_sprintf(buf, "%s/fof_tab_%03d", All.OutputDir, num);

  ngrps = All.NumFilesPerSnapshot / All.NumFilesWrittenInParallel;
  if((All.NumFilesPerSnapshot % All.NumFilesWrittenInParallel))
    ngrps++;

  for(gr = 0; gr < ngrps; gr++)
    {
      if((filenr / All.NumFilesWrittenInParallel) == gr) /* ok, it's this processor's turn */
        fof_subfind_write_file(buf, masterTask, lastTask);

      MPI_Barrier(MPI_COMM_WORLD);
    }

  myfree(CommBuffer);
  CommBuffer = NULL;

#ifdef FOF_STOREIDS
  myfree(ID_list);
#endif

  t1 = second();

  mpi_printf("FOF: Group catalogues saved. took = %g sec\n", timediff(t0, t1));
}

/*! \brief Prepares ID list for option FOF_STOREIDS.
 *
 *  \return void
 */
void fof_subfind_prepare_ID_list(void)
{
  int i, nids;
  long long totNids;
  double t0, t1;

  t0 = second();

  ID_list = (struct id_list *)mymalloc("ID_list", sizeof(struct id_list) * Nids);

  for(i = 0, nids = 0; i < NumPart; i++)
    {
      if(PS[i].GrNr < TotNgroups)
        {
          if(nids >= Nids)
            terminate("nids >= Nids");

          ID_list[nids].GrNr = PS[i].GrNr;
          ID_list[nids].Type = P[i].Type;
          ID_list[nids].ID   = P[i].ID;
#ifdef SUBFIND
          ID_list[nids].SubNr      = PS[i].SubNr;
          ID_list[nids].BindingEgy = PS[i].BindingEnergy;
#endif
          nids++;
        }
    }

  sumup_large_ints(1, &nids, &totNids);
  if(totNids != TotNids)
    terminate("Task=%d Nids=%d totNids=%lld TotNids=%lld\n", ThisTask, Nids, totNids, TotNids);

    /* sort the particle IDs according to group-number, and optionally subhalo number and binding energy  */
#ifdef SUBFIND
  parallel_sort(ID_list, Nids, sizeof(struct id_list), subfind_compare_ID_list);
#else
  parallel_sort(ID_list, Nids, sizeof(struct id_list), fof_compare_ID_list_GrNrID);
#endif

  t1 = second();
  mpi_printf("FOF/SUBFIND: Particle/cell IDs in groups globally sorted. took = %g sec\n", timediff(t0, t1));
}

/*! \brief Writes a file with name fname containing data from writeTask to
 *         lastTask.
 *
 *  \param[in] fname Filename of the output file.
 *  \param[in] writeTask Task responsible for writing the file.
 *  \param[in] lastTask Last task whose data is still in this file.
 *
 *  \return void
 */
void fof_subfind_write_file(const char *const fname, const int writeTask, const int lastTask)
{
  int bytes_per_blockelement, npart, nextblock;
  int n_for_this_task, n, p, pc, task;
  int blockmaxlen, n_type[3], ntot_type[3], nn[3];
  enum fof_subfind_iofields blocknr;
  char label[IO_LABEL_SIZE];
  int bnr;
  int blksize;
  MPI_Status status;
  FILE *fd = 0;
#ifdef HAVE_HDF5
  hid_t hdf5_file = 0, hdf5_grp[3], hdf5_headergrp = 0, hdf5_dataspace_memory;
  hid_t hdf5_datatype = 0, hdf5_dataspace_in_file = 0, hdf5_dataset = 0;
  hid_t hdf5_paramsgrp = 0, hdf5_configgrp = 0;
  herr_t hdf5_status;
  hsize_t dims[2], count[2], start[2];
  int rank = 0, pcsum = 0;
  char dataset_name[IO_DATASET_NAME_SIZE];
#endif

#define SKIP                                 \
  {                                          \
    my_fwrite(&blksize, sizeof(int), 1, fd); \
  }

  /* determine group/id numbers of each type in file */

  n_type[0] = Ngroups;
  n_type[1] = Nsubgroups;
  n_type[2] = Nids;

  if(ThisTask == writeTask)
    {
      for(n = 0; n < 3; n++)
        ntot_type[n] = n_type[n];

      for(task = writeTask + 1; task <= lastTask; task++)
        {
          MPI_Recv(&nn[0], 3, MPI_INT, task, TAG_LOCALN, MPI_COMM_WORLD, &status);
          for(n = 0; n < 3; n++)
            ntot_type[n] += nn[n];
        }

      for(task = writeTask + 1; task <= lastTask; task++)
        MPI_Send(&ntot_type[0], 3, MPI_INT, task, TAG_N, MPI_COMM_WORLD);
    }
  else
    {
      MPI_Send(&n_type[0], 3, MPI_INT, writeTask, TAG_LOCALN, MPI_COMM_WORLD);
      MPI_Recv(&ntot_type[0], 3, MPI_INT, writeTask, TAG_N, MPI_COMM_WORLD, &status);
    }

  /* fill file header */

  catalogue_header.Ngroups    = ntot_type[0];
  catalogue_header.Nsubgroups = ntot_type[1];
  catalogue_header.Nids       = ntot_type[2];

  catalogue_header.TotNgroups    = TotNgroups;
  catalogue_header.TotNsubgroups = TotNsubgroups;
  catalogue_header.TotNids       = TotNids;

  catalogue_header.num_files = All.NumFilesPerSnapshot;

  catalogue_header.time = All.Time;
  if(All.ComovingIntegrationOn)
    catalogue_header.redshift = 1.0 / All.Time - 1;
  else
    catalogue_header.redshift = 0;
  catalogue_header.HubbleParam = All.HubbleParam;
  catalogue_header.BoxSize     = All.BoxSize;
  catalogue_header.Omega0      = All.Omega0;
  catalogue_header.OmegaLambda = All.OmegaLambda;

#ifdef OUTPUT_IN_DOUBLEPRECISION
  catalogue_header.flag_doubleprecision = 1;
#else
  catalogue_header.flag_doubleprecision = 0;
#endif

  /* open file and write header */

  if(ThisTask == writeTask)
    {
      if(All.SnapFormat == SNAP_FORMAT_HDF5)
        {
#ifdef HAVE_HDF5
          char buf[MAXLEN_PATH];
          file_path_sprintf(buf, "%s.hdf5", fname);
          hdf5_file = my_H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
          mpi_printf("FOF/SUBFIND: writing group catalogue: '%s' (file 1 of %d)\n", fname, All.NumFilesPerSnapshot);
          hdf5_headergrp = my_H5Gcreate(hdf5_file, "/Header", 0);

          hdf5_grp[0] = my_H5Gcreate(hdf5_file, "/Group", 0);
          hdf5_grp[1] = my_H5Gcreate(hdf5_file, "/Subhalo", 0);
          hdf5_grp[2] = my_H5Gcreate(hdf5_file, "/IDs", 0);

          fof_subfind_write_header_attributes_in_hdf5(hdf5_headergrp);

          hdf5_paramsgrp = my_H5Gcreate(hdf5_file, "/Parameters", 0);
          write_parameters_attributes_in_hdf5(hdf5_paramsgrp);

          hdf5_configgrp = my_H5Gcreate(hdf5_file, "/Config", 0);
          write_compile_time_options_in_hdf5(hdf5_configgrp);

#endif
        }
      else
        {
          if(!(fd = fopen(fname, "w")))
            {
              printf("can't open file `%s' for writing snapshot.\n", fname);
              terminate("file open error");
            }

          mpi_printf("FOF/SUBFIND: writing group catalogue: '%s' (file 1 of %d)\n", fname, All.NumFilesPerSnapshot);

          if(All.SnapFormat == SNAP_FORMAT_GADGET_VARIANT)
            {
              blksize = sizeof(int) + 4 * sizeof(char);
              SKIP;
              my_fwrite((void *)"HEAD", sizeof(char), 4, fd);
              nextblock = sizeof(catalogue_header) + 2 * sizeof(int);
              my_fwrite(&nextblock, sizeof(int), 1, fd);
              SKIP;
            }

          blksize = sizeof(catalogue_header);

          SKIP;
          my_fwrite(&catalogue_header, sizeof(catalogue_header), 1, fd);
          SKIP;
        }
    }

  for(bnr = 0; bnr < 1000; bnr++)
    {
      blocknr = (enum fof_subfind_iofields)bnr;

      if(blocknr == IO_FOF_LASTENTRY)
        break;

      if(fof_subfind_blockpresent(blocknr))
        {
          bytes_per_blockelement = fof_subfind_get_bytes_per_blockelement(blocknr);

          blockmaxlen = (int)(COMMBUFFERSIZE / bytes_per_blockelement);

          npart   = fof_subfind_get_particles_in_block(blocknr);
          int grp = fof_subfind_get_dataset_group(blocknr);

          if(npart > 0)
            {
              if(ThisTask == 0)
                {
                  char tmp[IO_DATASET_NAME_SIZE];

                  fof_subfind_get_dataset_name(blocknr, tmp);
                  printf("FOF/SUBFIND: writing block %d (%s)...\n", (int)blocknr, tmp);
                }

              if(ThisTask == writeTask)
                {
                  if(All.SnapFormat == SNAP_FORMAT_GADGET || All.SnapFormat == SNAP_FORMAT_GADGET_VARIANT)
                    {
                      if(All.SnapFormat == SNAP_FORMAT_GADGET_VARIANT)
                        {
                          blksize = sizeof(int) + IO_LABEL_SIZE * sizeof(char);
                          SKIP;
                          fof_subfind_get_Tab_IO_Label(blocknr, label);
                          my_fwrite(label, sizeof(char), IO_LABEL_SIZE, fd);
                          nextblock = npart * bytes_per_blockelement + 2 * sizeof(int);
                          my_fwrite(&nextblock, sizeof(int), 1, fd);
                          SKIP;
                        }

                      blksize = npart * bytes_per_blockelement;
                      SKIP;
                    }
                  else if(All.SnapFormat == SNAP_FORMAT_HDF5)
                    {
#ifdef HAVE_HDF5
                      switch(fof_subfind_get_datatype(blocknr))
                        {
                          case 0:
                            hdf5_datatype = my_H5Tcopy(H5T_NATIVE_INT);
                            break;
                          case 1:
#ifdef OUTPUT_IN_DOUBLEPRECISION
                            hdf5_datatype = my_H5Tcopy(H5T_NATIVE_DOUBLE);
#else
                            hdf5_datatype = my_H5Tcopy(H5T_NATIVE_FLOAT);
#endif
                            break;
                          case 2:
                            hdf5_datatype = my_H5Tcopy(H5T_NATIVE_UINT64);
                            break;
                        }

                      dims[0] = ntot_type[grp];
                      dims[1] = fof_subfind_get_values_per_blockelement(blocknr);
                      if(dims[1] == 1)
                        rank = 1;
                      else
                        rank = 2;

                      fof_subfind_get_dataset_name(blocknr, dataset_name);

                      hdf5_dataspace_in_file = my_H5Screate_simple(rank, dims, NULL);

                      hdf5_dataset = my_H5Dcreate(hdf5_grp[grp], dataset_name, hdf5_datatype, hdf5_dataspace_in_file, H5P_DEFAULT);

                      pcsum = 0;
#endif
                    }
                }

              for(task = writeTask; task <= lastTask; task++)
                {
                  if(task == ThisTask)
                    {
                      n_for_this_task = n_type[grp];

                      for(p = writeTask; p <= lastTask; p++)
                        if(p != ThisTask)
                          MPI_Send(&n_for_this_task, 1, MPI_INT, p, TAG_NFORTHISTASK, MPI_COMM_WORLD);
                    }
                  else
                    MPI_Recv(&n_for_this_task, 1, MPI_INT, task, TAG_NFORTHISTASK, MPI_COMM_WORLD, &status);

                  while(n_for_this_task > 0)
                    {
                      pc = n_for_this_task;

                      if(pc > blockmaxlen)
                        pc = blockmaxlen;

                      if(ThisTask == task)
                        fof_subfind_fill_write_buffer(blocknr, 0, pc);

                      if(ThisTask == writeTask && task != writeTask)
                        MPI_Recv(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, task, TAG_PDATA, MPI_COMM_WORLD, &status);

                      if(ThisTask != writeTask && task == ThisTask)
                        MPI_Ssend(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, writeTask, TAG_PDATA, MPI_COMM_WORLD);

                      if(ThisTask == writeTask)
                        {
                          if(All.SnapFormat == SNAP_FORMAT_HDF5)
                            {
#ifdef HAVE_HDF5
                              start[0] = pcsum;
                              start[1] = 0;

                              count[0] = pc;
                              count[1] = fof_subfind_get_values_per_blockelement(blocknr);
                              pcsum += pc;

                              my_H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET, start, NULL, count, NULL);

                              dims[0]               = pc;
                              dims[1]               = fof_subfind_get_values_per_blockelement(blocknr);
                              hdf5_dataspace_memory = my_H5Screate_simple(rank, dims, NULL);

                              hdf5_status = my_H5Dwrite(hdf5_dataset, hdf5_datatype, hdf5_dataspace_memory, hdf5_dataspace_in_file,
                                                        H5P_DEFAULT, CommBuffer, dataset_name);

                              (void)hdf5_status;

                              my_H5Sclose(hdf5_dataspace_memory, H5S_SIMPLE);
#endif
                            }
                          else
                            {
                              my_fwrite(CommBuffer, bytes_per_blockelement, pc, fd);
                            }
                        }

                      n_for_this_task -= pc;
                    }
                }

              if(ThisTask == writeTask)
                {
                  if(All.SnapFormat == 3)
                    {
#ifdef HAVE_HDF5
                      my_H5Dclose(hdf5_dataset, dataset_name);
                      my_H5Sclose(hdf5_dataspace_in_file, H5S_SIMPLE);
                      my_H5Tclose(hdf5_datatype);
#endif
                    }
                  else
                    SKIP;
                }
            }
        }
    }

  if(ThisTask == writeTask)
    {
      if(All.SnapFormat == SNAP_FORMAT_HDF5)
        {
#ifdef HAVE_HDF5
          my_H5Gclose(hdf5_grp[0], "/Group");
          my_H5Gclose(hdf5_grp[1], "/Subhalo");
          my_H5Gclose(hdf5_grp[2], "/IDs");
          my_H5Gclose(hdf5_headergrp, "/Header");
          my_H5Gclose(hdf5_paramsgrp, "/Parameters");
          my_H5Gclose(hdf5_configgrp, "/Config");

          my_H5Fclose(hdf5_file, fname);
#endif
        }
      else
        fclose(fd);
    }
}

/*! \brief Copies data from global group array to appropriate output buffer.
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *  \param[in] startindex First particle index to be included.
 *  \param[in] pc Particle count; number of particles to be written.
 *
 *  \return void
 */
void fof_subfind_fill_write_buffer(const enum fof_subfind_iofields blocknr, const int startindex, const int pc)
{
  MyOutputFloat *fp = (MyOutputFloat *)CommBuffer;
  int *ip           = (int *)CommBuffer;
#if defined(FOF_STOREIDS) || defined(SUBFIND)
  MyIDType *idp = (MyIDType *)CommBuffer;
#endif

#ifdef FOF_FUZZ_SORT_BY_NEAREST_GROUP
  unsigned long long *llp = (unsigned long long *)CommBuffer;
#endif

  for(int pindex = startindex; pindex < startindex + pc; pindex++)
    {
      switch(blocknr)
        {
          case IO_FOF_LEN:
            *ip++ = Group[pindex].Len;
            break;
          case IO_FOF_MTOT:
            *fp++ = Group[pindex].Mass;
            break;
          case IO_FOF_POS:
            for(int k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = wrap_position_shift(Group[pindex].Pos[k], k);
#else
              *fp++ = wrap_position_shift(Group[pindex].CM[k], k);
#endif
            break;
          case IO_FOF_CM:
            for(int k = 0; k < 3; k++)
              *fp++ = wrap_position_shift(Group[pindex].CM[k], k);
            break;
          case IO_FOF_POSMINPOT:
#if defined(GFM_BIPOLAR_WINDS) && (GFM_BIPOLAR_WINDS == 3)
            for(int k = 0; k < 3; k++)
              *fp++ = Group[pindex].Pos_MinPotential[k];
#endif
            break;
          case IO_FOF_VEL:
            for(int k = 0; k < 3; k++)
              *fp++ = Group[pindex].Vel[k];
            break;
          case IO_FOF_LENTYPE:
            for(int k = 0; k < NTYPES; k++)
              *ip++ = Group[pindex].LenType[k];
            break;
          case IO_FOF_MASSTYPE:
            for(int k = 0; k < NTYPES; k++)
              *fp++ = Group[pindex].MassType[k];
            break;
          case IO_FOF_SFR:
#ifdef USE_SFR
            *fp++ = Group[pindex].Sfr;
#endif
            break;
          case IO_FOF_GASMETAL:
#ifdef GFM_STELLAR_EVOLUTION
            if(Group[pindex].MassType[0] > 0.0)
              *fp++ = Group[pindex].GasMassMetallicity / Group[pindex].MassType[0];
            else
              *fp++ = 0.0;
#endif
            break;
          case IO_FOF_STARMETAL:
#ifdef GFM_STELLAR_EVOLUTION
            if(Group[pindex].MassType[4] > 0.0)
              *fp++ = Group[pindex].StellarMassMetallicity / Group[pindex].MassType[4];
            else
              *fp++ = 0.0;
#endif
            break;
          case IO_FOF_GASMETALELEMENTS:
#ifdef GFM_STELLAR_EVOLUTION
            if(Group[pindex].MassType[0] > 0.0)
              for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                *fp++ = Group[pindex].GasMassMetals[k] / Group[pindex].MassType[0];
            else
              for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                *fp++ = 0.0;
#endif
            break;
          case IO_FOF_STARMETALELEMENTS:
#ifdef GFM_STELLAR_EVOLUTION
            if(Group[pindex].MassType[4] > 0.0)
              for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                *fp++ = Group[pindex].StellarMassMetals[k] / Group[pindex].MassType[4];
            else
              for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                *fp++ = 0.0;
#endif
            break;
          case IO_FOF_GASDUSTMETAL:
#ifdef GFM_DUST
            if(Group[pindex].MassType[0] > 0.0)
              *fp++ = Group[pindex].GasMassDustMetallicity / Group[pindex].MassType[0];
            else
              *fp++ = 0.0;
#endif
            break;
          case IO_FOF_BHMASS:
#ifdef BLACK_HOLES
            *fp++ = Group[pindex].BH_Mass;
#endif
            break;
          case IO_FOF_BHMDOT:
#ifdef BLACK_HOLES
            *fp++ = Group[pindex].BH_Mdot;
#endif
            break;
          case IO_FOF_WINDMASS:
#ifdef GFM_WINDS
            *fp++ = Group[pindex].WindMass;
#endif
            break;
          case IO_FOF_M_MEAN200:
#ifdef SUBFIND
            *fp++ = Group[pindex].M_Mean200;
#endif
            break;
          case IO_FOF_R_MEAN200:
#ifdef SUBFIND
            *fp++ = Group[pindex].R_Mean200;
#endif
            break;

#ifdef SUBFIND_EXTENDED_PROPERTIES
          case IO_FOF_J_MEAN200:
#ifdef SUBFIND
            for(int k = 0; k < 3; k++)
              *fp++ = Group[pindex].J_Mean200[k];
#endif
            break;
          case IO_FOF_JDM_MEAN200:
#ifdef SUBFIND
            for(int k = 0; k < 3; k++)
              *fp++ = Group[pindex].JDM_Mean200[k];
#endif
            break;
          case IO_FOF_JGAS_MEAN200:
#ifdef SUBFIND
            for(int k = 0; k < 3; k++)
              *fp++ = Group[pindex].JGas_Mean200[k];
#endif
            break;
          case IO_FOF_JSTARS_MEAN200:
#ifdef SUBFIND
            for(int k = 0; k < 3; k++)
              *fp++ = Group[pindex].JStars_Mean200[k];
#endif
            break;
          case IO_FOF_MASSTYPE_MEAN200:
#ifdef SUBFIND
            for(int k = 0; k < NTYPES; k++)
              *fp++ = Group[pindex].MassType_Mean200[k];
#endif
            break;
          case IO_FOF_LENTYPE_MEAN200:
#ifdef SUBFIND
            for(int k = 0; k < NTYPES; k++)
              *ip++ = Group[pindex].LenType_Mean200[k];
#endif
            break;
          case IO_FOF_CMFRAC_MEAN200:
#ifdef SUBFIND
            *fp++ = Group[pindex].CMFrac_Mean200;
#endif
            break;
          case IO_FOF_CMFRACTYPE_MEAN200:
#ifdef SUBFIND
            for(int k = 0; k < NTYPES; k++)
              *fp++ = Group[pindex].CMFracType_Mean200[k];
#endif
            break;
          case IO_FOF_J_CRIT200:
#ifdef SUBFIND
            for(int k = 0; k < 3; k++)
              *fp++ = Group[pindex].J_Crit200[k];
#endif
            break;
          case IO_FOF_JDM_CRIT200:
#ifdef SUBFIND
            for(int k = 0; k < 3; k++)
              *fp++ = Group[pindex].JDM_Crit200[k];
#endif
            break;
          case IO_FOF_JGAS_CRIT200:
#ifdef SUBFIND
            for(int k = 0; k < 3; k++)
              *fp++ = Group[pindex].JGas_Crit200[k];
#endif
            break;
          case IO_FOF_JSTARS_CRIT200:
#ifdef SUBFIND
            for(int k = 0; k < 3; k++)
              *fp++ = Group[pindex].JStars_Crit200[k];
#endif
            break;
          case IO_FOF_MASSTYPE_CRIT200:
#ifdef SUBFIND
            for(int k = 0; k < NTYPES; k++)
              *fp++ = Group[pindex].MassType_Crit200[k];
#endif
            break;
          case IO_FOF_LENTYPE_CRIT200:
#ifdef SUBFIND
            for(int k = 0; k < NTYPES; k++)
              *ip++ = Group[pindex].LenType_Crit200[k];
#endif
            break;
          case IO_FOF_CMFRAC_CRIT200:
#ifdef SUBFIND
            *fp++ = Group[pindex].CMFrac_Crit200;
#endif
            break;
          case IO_FOF_CMFRACTYPE_CRIT200:
#ifdef SUBFIND
            for(int k = 0; k < NTYPES; k++)
              *fp++ = Group[pindex].CMFracType_Crit200[k];
#endif
            break;
          case IO_FOF_J_CRIT500:
#ifdef SUBFIND
            for(int k = 0; k < 3; k++)
              *fp++ = Group[pindex].J_Crit500[k];
#endif
            break;
          case IO_FOF_JDM_CRIT500:
#ifdef SUBFIND
            for(int k = 0; k < 3; k++)
              *fp++ = Group[pindex].JDM_Crit500[k];
#endif
            break;
          case IO_FOF_JGAS_CRIT500:
#ifdef SUBFIND
            for(int k = 0; k < 3; k++)
              *fp++ = Group[pindex].JGas_Crit500[k];
#endif
            break;
          case IO_FOF_JSTARS_CRIT500:
#ifdef SUBFIND
            for(int k = 0; k < 3; k++)
              *fp++ = Group[pindex].JStars_Crit500[k];
#endif
            break;
          case IO_FOF_MASSTYPE_CRIT500:
#ifdef SUBFIND
            for(int k = 0; k < NTYPES; k++)
              *fp++ = Group[pindex].MassType_Crit500[k];
#endif
            break;
          case IO_FOF_LENTYPE_CRIT500:
#ifdef SUBFIND
            for(int k = 0; k < NTYPES; k++)
              *ip++ = Group[pindex].LenType_Crit500[k];
#endif
            break;
          case IO_FOF_CMFRAC_CRIT500:
#ifdef SUBFIND
            *fp++ = Group[pindex].CMFrac_Crit500;
#endif
            break;
          case IO_FOF_CMFRACTYPE_CRIT500:
#ifdef SUBFIND
            for(int k = 0; k < NTYPES; k++)
              *fp++ = Group[pindex].CMFracType_Crit500[k];
#endif
            break;
          case IO_FOF_J_TOPHAT200:
#ifdef SUBFIND
            for(int k = 0; k < 3; k++)
              *fp++ = Group[pindex].J_TopHat200[k];
#endif
            break;
          case IO_FOF_JDM_TOPHAT200:
#ifdef SUBFIND
            for(int k = 0; k < 3; k++)
              *fp++ = Group[pindex].JDM_TopHat200[k];
#endif
            break;
          case IO_FOF_JGAS_TOPHAT200:
#ifdef SUBFIND
            for(int k = 0; k < 3; k++)
              *fp++ = Group[pindex].JGas_TopHat200[k];
#endif
            break;
          case IO_FOF_JSTARS_TOPHAT200:
#ifdef SUBFIND
            for(int k = 0; k < 3; k++)
              *fp++ = Group[pindex].JStars_TopHat200[k];
#endif
            break;
          case IO_FOF_MASSTYPE_TOPHAT200:
#ifdef SUBFIND
            for(int k = 0; k < NTYPES; k++)
              *fp++ = Group[pindex].MassType_TopHat200[k];
#endif
            break;
          case IO_FOF_LENTYPE_TOPHAT200:
#ifdef SUBFIND
            for(int k = 0; k < NTYPES; k++)
              *ip++ = Group[pindex].LenType_TopHat200[k];
#endif
            break;
          case IO_FOF_CMFRAC_TOPHAT200:
#ifdef SUBFIND
            *fp++ = Group[pindex].CMFrac_TopHat200;
#endif
            break;
          case IO_FOF_CMFRACTYPE_TOPHAT200:
#ifdef SUBFIND
            for(int k = 0; k < NTYPES; k++)
              *fp++ = Group[pindex].CMFracType_TopHat200[k];
#endif
            break;

          case IO_FOF_EPOT_CRIT200:
#ifdef SUBFIND
            *fp++ = Group[pindex].Epot_Crit200;
#endif
            break;
          case IO_FOF_EKIN_CRIT200:
#ifdef SUBFIND
            *fp++ = Group[pindex].Ekin_Crit200;
#endif
            break;
          case IO_FOF_ETHR_CRIT200:
#ifdef SUBFIND
            *fp++ = Group[pindex].Ethr_Crit200;
#endif
            break;
          case IO_FOF_EPOT_MEAN200:
#ifdef SUBFIND
            *fp++ = Group[pindex].Epot_Mean200;
#endif
            break;
          case IO_FOF_EKIN_MEAN200:
#ifdef SUBFIND
            *fp++ = Group[pindex].Ekin_Mean200;
#endif
            break;
          case IO_FOF_ETHR_MEAN200:
#ifdef SUBFIND
            *fp++ = Group[pindex].Ethr_Mean200;
#endif
            break;
          case IO_FOF_EPOT_TOPHAT200:
#ifdef SUBFIND
            *fp++ = Group[pindex].Epot_TopHat200;
#endif
            break;
          case IO_FOF_EKIN_TOPHAT200:
#ifdef SUBFIND
            *fp++ = Group[pindex].Ekin_TopHat200;
#endif
            break;
          case IO_FOF_ETHR_TOPHAT200:
#ifdef SUBFIND
            *fp++ = Group[pindex].Ethr_TopHat200;
#endif
            break;
          case IO_FOF_EPOT_CRIT500:
#ifdef SUBFIND
            *fp++ = Group[pindex].Epot_Crit500;
#endif
            break;
          case IO_FOF_EKIN_CRIT500:
#ifdef SUBFIND
            *fp++ = Group[pindex].Ekin_Crit500;
#endif
            break;
          case IO_FOF_ETHR_CRIT500:
#ifdef SUBFIND
            *fp++ = Group[pindex].Ethr_Crit500;
#endif
            break;

          case IO_FOF_J:
#ifdef SUBFIND
            for(int k = 0; k < 3; k++)
              *fp++ = Group[pindex].J[k];
#endif
            break;
          case IO_FOF_JDM:
#ifdef SUBFIND
            for(int k = 0; k < 3; k++)
              *fp++ = Group[pindex].JDM[k];
#endif
            break;
          case IO_FOF_JGAS:
#ifdef SUBFIND
            for(int k = 0; k < 3; k++)
              *fp++ = Group[pindex].JGas[k];
#endif
            break;
          case IO_FOF_JSTARS:
#ifdef SUBFIND
            for(int k = 0; k < 3; k++)
              *fp++ = Group[pindex].JStars[k];
#endif
            break;
          case IO_FOF_CMFRAC:
#ifdef SUBFIND
            *fp++ = Group[pindex].CMFrac;
#endif
            break;
          case IO_FOF_CMFRACTYPE:
#ifdef SUBFIND
            for(int k = 0; k < NTYPES; k++)
              *fp++ = Group[pindex].CMFracType[k];
#endif
            break;
          case IO_FOF_EKIN:
#ifdef SUBFIND
            *fp++ = Group[pindex].Ekin;
#endif
            break;
          case IO_FOF_ETHR:
#ifdef SUBFIND
            *fp++ = Group[pindex].Ethr;
#endif
            break;
          case IO_FOF_EPOT:
#ifdef SUBFIND
            *fp++ = Group[pindex].Epot;
#endif
            break;
          case IO_SUB_EKIN:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].Ekin;
#endif
            break;
          case IO_SUB_ETHR:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].Ethr;
#endif
            break;
          case IO_SUB_EPOT:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].Epot;
#endif
            break;
          case IO_SUB_J:
            for(int k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].J[k];
#endif
            break;
          case IO_SUB_JDM:
            for(int k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].Jdm[k];
#endif
            break;
          case IO_SUB_JGAS:
            for(int k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].Jgas[k];
#endif
            break;
          case IO_SUB_JSTARS:
            for(int k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].Jstars[k];
#endif
            break;
          case IO_SUB_JINHALFRAD:
            for(int k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].J_inHalfRad[k];
#endif
            break;
          case IO_SUB_JDMINHALFRAD:
            for(int k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].Jdm_inHalfRad[k];
#endif
            break;
          case IO_SUB_JGASINHALFRAD:
            for(int k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].Jgas_inHalfRad[k];
#endif
            break;
          case IO_SUB_JSTARSINHALFRAD:
            for(int k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].Jstars_inHalfRad[k];
#endif
            break;
          case IO_SUB_JINRAD:
            for(int k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].J_inRad[k];
#endif
            break;
          case IO_SUB_JDMINRAD:
            for(int k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].Jdm_inRad[k];
#endif
            break;
          case IO_SUB_JGASINRAD:
            for(int k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].Jgas_inRad[k];
#endif
            break;
          case IO_SUB_JSTARSINRAD:
            for(int k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].Jstars_inRad[k];
#endif
            break;
          case IO_SUB_CMFRAC:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].CMFrac;
#endif
            break;
          case IO_SUB_CMFRACTYPE:
#ifdef SUBFIND
            for(int k = 0; k < NTYPES; k++)
              *fp++ = SubGroup[pindex].CMFracType[k];
#endif
            break;
          case IO_SUB_CMFRACINHALFRAD:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].CMFrac_inHalfRad;
#endif
            break;
          case IO_SUB_CMFRACTYPEINHALFRAD:
#ifdef SUBFIND
            for(int k = 0; k < NTYPES; k++)
              *fp++ = SubGroup[pindex].CMFracType_inHalfRad[k];
#endif
            break;
          case IO_SUB_CMFRACINRAD:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].CMFrac_inRad;
#endif
            break;
          case IO_SUB_CMFRACTYPEINRAD:
#ifdef SUBFIND
            for(int k = 0; k < NTYPES; k++)
              *fp++ = SubGroup[pindex].CMFracType_inRad[k];
#endif
            break;
#endif

            break;
          case IO_FOF_M_CRIT200:
#ifdef SUBFIND
            *fp++ = Group[pindex].M_Crit200;
#endif
            break;
          case IO_FOF_R_CRIT200:
#ifdef SUBFIND
            *fp++ = Group[pindex].R_Crit200;
#endif
            break;
          case IO_FOF_M_CRIT500:
#ifdef SUBFIND
            *fp++ = Group[pindex].M_Crit500;
#endif
            break;
          case IO_FOF_R_CRIT500:
#ifdef SUBFIND
            *fp++ = Group[pindex].R_Crit500;
#endif
            break;
          case IO_FOF_M_TOPHAT200:
#ifdef SUBFIND
            *fp++ = Group[pindex].M_TopHat200;
#endif
            break;
          case IO_FOF_R_TOPHAT200:
#ifdef SUBFIND
            *fp++ = Group[pindex].R_TopHat200;
#endif
            break;
          case IO_FOF_NSUBS:
#ifdef SUBFIND
            *ip++ = Group[pindex].Nsubs;
#endif
            break;
          case IO_FOF_FIRSTSUB:
#ifdef SUBFIND
            *ip++ = Group[pindex].FirstSub;
#endif
            break;
          case IO_FOF_FUZZOFFTYPE:
#ifdef FOF_FUZZ_SORT_BY_NEAREST_GROUP
            for(int k = 0; k < NTYPES; k++)
              *llp++ = Group[pindex].FuzzOffsetType[k];
#endif
            break;
          case IO_FOF_XRAYLUM:
#ifdef BH_NF_RADIO
            *fp++ = Group[pindex].XrayLum;
#endif
            break;
          case IO_FOF_RADIOLUM:
#ifdef BH_NF_RADIO
            *fp++ = Group[pindex].RadioLum;
#endif
            break;
          case IO_FOF_DENSLVEC:
#if defined(GFM_BIPOLAR_WINDS) && (GFM_BIPOLAR_WINDS == 3)
            for(int k = 0; k < 3; k++)
              *fp++ = Group[pindex].DensGasAngMomentum[k];
#endif
            break;
          case IO_SUB_LEN:
#ifdef SUBFIND
            *ip++ = SubGroup[pindex].Len;
#endif
            break;
          case IO_SUB_MTOT:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].Mass;
#endif
            break;
          case IO_SUB_POS:
#ifdef SUBFIND
            for(int k = 0; k < 3; k++)
              *fp++ = wrap_position_shift(SubGroup[pindex].Pos[k], k);
#endif
            break;
          case IO_SUB_VEL:
#ifdef SUBFIND
            for(int k = 0; k < 3; k++)
              *fp++ = SubGroup[pindex].Vel[k];
#endif
            break;
          case IO_SUB_LENTYPE:
#ifdef SUBFIND
            for(int k = 0; k < NTYPES; k++)
              *ip++ = SubGroup[pindex].LenType[k];
#endif
            break;
          case IO_SUB_MASSTYPE:
#ifdef SUBFIND
            for(int k = 0; k < NTYPES; k++)
              *fp++ = SubGroup[pindex].MassType[k];
#endif
            break;
          case IO_SUB_CM:
#ifdef SUBFIND
            for(int k = 0; k < 3; k++)
              *fp++ = wrap_position_shift(SubGroup[pindex].CM[k], k);
#endif
            break;
          case IO_SUB_SPIN:
            for(int k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].Spin[k];
#endif
            break;
          case IO_SUB_VELDISP:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].SubVelDisp;
#endif
            break;
          case IO_SUB_VMAX:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].SubVmax;
#endif
            break;
          case IO_SUB_VMAXRAD:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].SubVmaxRad;
#endif
            break;
          case IO_SUB_HALFMASSRAD:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].SubHalfMassRad;
#endif
            break;
          case IO_SUB_HALFMASSRADTYPE:
#ifdef SUBFIND
            for(int k = 0; k < NTYPES; k++)
              *fp++ = SubGroup[pindex].SubHalfMassRadType[k];
#endif
            break;
          case IO_SUB_MASSINRAD:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].SubMassInRad;
#endif
            break;
          case IO_SUB_MASSINRADTYPE:
#ifdef SUBFIND
            for(int k = 0; k < NTYPES; k++)
              *fp++ = SubGroup[pindex].SubMassInRadType[k];
#endif
            break;
          case IO_SUB_MASSINHALFRAD:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].SubMassInHalfRad;
#endif
            break;
          case IO_SUB_MASSINHALFRADTYPE:
#ifdef SUBFIND
            for(int k = 0; k < NTYPES; k++)
              *fp++ = SubGroup[pindex].SubMassInHalfRadType[k];
#endif
            break;
          case IO_SUB_MASSINMAXRAD:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].SubMassInMaxRad;
#endif
            break;
          case IO_SUB_MASSINMAXRADTYPE:
#ifdef SUBFIND
            for(int k = 0; k < NTYPES; k++)
              *fp++ = SubGroup[pindex].SubMassInMaxRadType[k];
#endif
            break;
          case IO_SUB_IDMOSTBOUND:
#ifdef SUBFIND
            *idp++ = SubGroup[pindex].SubMostBoundID;
#endif
            break;
          case IO_SUB_GRNR:
#ifdef SUBFIND
            *ip++ = SubGroup[pindex].GrNr;
#endif
            break;
          case IO_SUB_PARENT:
#ifdef SUBFIND
            *ip++ = SubGroup[pindex].SubParent;
#endif
            break;
          case IO_SUB_BFLD_HALO:
#if defined(MHD) && defined(SUBFIND)
            *fp++ = SubGroup[pindex].Bfld_Halo * sqrt(4. * M_PI);
#endif
            break;
          case IO_SUB_BFLD_DISK:
#if defined(MHD) && defined(SUBFIND)
            *fp++ = SubGroup[pindex].Bfld_Disk * sqrt(4. * M_PI);
#endif
            break;
          case IO_SUB_SFR:
#if defined(USE_SFR) && defined(SUBFIND)
            *fp++ = SubGroup[pindex].Sfr;
#endif
            break;
          case IO_SUB_SFRINRAD:
#if defined(USE_SFR) && defined(SUBFIND)
            *fp++ = SubGroup[pindex].SfrInRad;
#endif
            break;
          case IO_SUB_SFRINHALFRAD:
#if defined(USE_SFR) && defined(SUBFIND)
            *fp++ = SubGroup[pindex].SfrInHalfRad;
#endif
            break;
          case IO_SUB_SFRINMAXRAD:
#if defined(USE_SFR) && defined(SUBFIND)
            *fp++ = SubGroup[pindex].SfrInMaxRad;
#endif
            break;
          case IO_SUB_GASMETAL:
#if defined(GFM_STELLAR_EVOLUTION) && defined(SUBFIND)
            if(SubGroup[pindex].SubMassInRadType[0] > 0.0)
              *fp++ = SubGroup[pindex].GasMassMetallicity / SubGroup[pindex].SubMassInRadType[0];
            else
              *fp++ = 0.0;
#endif
            break;
          case IO_SUB_GASMETALHALFRAD:
#if defined(GFM_STELLAR_EVOLUTION) && defined(SUBFIND)
            if(SubGroup[pindex].SubMassInHalfRadType[0] > 0.0)
              *fp++ = SubGroup[pindex].GasMassMetallicityHalfRad / SubGroup[pindex].SubMassInHalfRadType[0];
            else
              *fp++ = 0.0;
#endif
            break;
          case IO_SUB_GASMETALMAXRAD:
#if defined(GFM_STELLAR_EVOLUTION) && defined(SUBFIND)
            if(SubGroup[pindex].SubMassInMaxRadType[0] > 0.0)
              *fp++ = SubGroup[pindex].GasMassMetallicityMaxRad / SubGroup[pindex].SubMassInMaxRadType[0];
            else
              *fp++ = 0.0;
#endif
            break;
          case IO_SUB_GASMETALSFR:
#if defined(GFM_STELLAR_EVOLUTION) && defined(SUBFIND)
            if(SubGroup[pindex].GasMassSfr > 0.0)
              *fp++ = SubGroup[pindex].GasMassMetallicitySfr / SubGroup[pindex].GasMassSfr;
            else
              *fp++ = 0.0;
#endif
            break;
          case IO_SUB_GASMETALSFRWEIGHTED:
#if defined(GFM_STELLAR_EVOLUTION) && defined(SUBFIND)
            if(SubGroup[pindex].Sfr > 0.0)
              *fp++ = SubGroup[pindex].GasMassMetallicitySfrWeighted / SubGroup[pindex].Sfr;
            else
              *fp++ = 0.0;
#endif
            break;
          case IO_SUB_STARMETAL:
#if defined(GFM_STELLAR_EVOLUTION) && defined(SUBFIND)
            if(SubGroup[pindex].SubMassInRadType[4] > 0.0)
              *fp++ = SubGroup[pindex].StellarMassMetallicity / SubGroup[pindex].SubMassInRadType[4];
            else
              *fp++ = 0.0;
#endif
            break;
          case IO_SUB_STARMETALHALFRAD:
#if defined(GFM_STELLAR_EVOLUTION) && defined(SUBFIND)
            if(SubGroup[pindex].SubMassInHalfRadType[4] > 0.0)
              *fp++ = SubGroup[pindex].StellarMassMetallicityHalfRad / SubGroup[pindex].SubMassInHalfRadType[4];
            else
              *fp++ = 0.0;
#endif
            break;
          case IO_SUB_STARMETALMAXRAD:
#if defined(GFM_STELLAR_EVOLUTION) && defined(SUBFIND)
            if(SubGroup[pindex].SubMassInMaxRadType[4] > 0.0)
              *fp++ = SubGroup[pindex].StellarMassMetallicityMaxRad / SubGroup[pindex].SubMassInMaxRadType[4];
            else
              *fp++ = 0.0;
#endif
            break;
          case IO_SUB_GASMETALELEMENTS:
#if defined(GFM_STELLAR_EVOLUTION) && defined(SUBFIND)
            if(SubGroup[pindex].SubMassInRadType[0] > 0.0)
              for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                *fp++ = SubGroup[pindex].GasMassMetals[k] / SubGroup[pindex].SubMassInRadType[0];
            else
              for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                *fp++ = 0.0;
#endif
            break;
          case IO_SUB_GASMETALELEMENTSHALFRAD:
#if defined(GFM_STELLAR_EVOLUTION) && defined(SUBFIND)
            if(SubGroup[pindex].SubMassInHalfRadType[0] > 0.0)
              for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                *fp++ = SubGroup[pindex].GasMassMetalsHalfRad[k] / SubGroup[pindex].SubMassInHalfRadType[0];
            else
              for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                *fp++ = 0.0;
#endif
            break;
          case IO_SUB_GASMETALELEMENTSMAXRAD:
#if defined(GFM_STELLAR_EVOLUTION) && defined(SUBFIND)
            if(SubGroup[pindex].SubMassInMaxRadType[0] > 0.0)
              for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                *fp++ = SubGroup[pindex].GasMassMetalsMaxRad[k] / SubGroup[pindex].SubMassInMaxRadType[0];
            else
              for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                *fp++ = 0.0;
#endif
            break;
          case IO_SUB_GASMETALELEMENTSSFR:
#if defined(GFM_STELLAR_EVOLUTION) && defined(SUBFIND)
            if(SubGroup[pindex].GasMassSfr > 0.0)
              for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                *fp++ = SubGroup[pindex].GasMassMetalsSfr[k] / SubGroup[pindex].GasMassSfr;
            else
              for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                *fp++ = 0.0;
#endif
            break;
          case IO_SUB_GASMETALELEMENTSSFRWEIGHTED:
#if defined(GFM_STELLAR_EVOLUTION) && defined(SUBFIND)
            if(SubGroup[pindex].Sfr > 0.0)
              for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                *fp++ = SubGroup[pindex].GasMassMetalsSfrWeighted[k] / SubGroup[pindex].Sfr;
            else
              for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                *fp++ = 0.0;
#endif
            break;
          case IO_SUB_STARMETALELEMENTS:
#if defined(GFM_STELLAR_EVOLUTION) && defined(SUBFIND)
            if(SubGroup[pindex].SubMassInRadType[4] > 0.0)
              for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                *fp++ = SubGroup[pindex].StellarMassMetals[k] / SubGroup[pindex].SubMassInRadType[4];
            else
              for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                *fp++ = 0.0;
#endif
            break;
          case IO_SUB_STARMETALELEMENTSHALFRAD:
#if defined(GFM_STELLAR_EVOLUTION) && defined(SUBFIND)
            if(SubGroup[pindex].SubMassInHalfRadType[4] > 0.0)
              for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                *fp++ = SubGroup[pindex].StellarMassMetalsHalfRad[k] / SubGroup[pindex].SubMassInHalfRadType[4];
            else
              for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                *fp++ = 0.0;
#endif
            break;
          case IO_SUB_STARMETALELEMENTSMAXRAD:
#if defined(GFM_STELLAR_EVOLUTION) && defined(SUBFIND)
            if(SubGroup[pindex].SubMassInMaxRadType[4] > 0.0)
              for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                *fp++ = SubGroup[pindex].StellarMassMetalsMaxRad[k] / SubGroup[pindex].SubMassInMaxRadType[4];
            else
              for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                *fp++ = 0.0;
#endif
            break;
          case IO_SUB_GASDUSTMETAL:
#if defined(GFM_DUST) && defined(SUBFIND)
            if(SubGroup[pindex].SubMassInRadType[0] > 0.0)
              *fp++ = SubGroup[pindex].GasMassDustMetallicity / SubGroup[pindex].SubMassInRadType[0];
            else
              *fp++ = 0.0;
#endif
            break;
          case IO_SUB_GASDUSTMETALHALFRAD:
#if defined(GFM_DUST) && defined(SUBFIND)
            if(SubGroup[pindex].SubMassInHalfRadType[0] > 0.0)
              *fp++ = SubGroup[pindex].GasMassDustMetallicityHalfRad / SubGroup[pindex].SubMassInHalfRadType[0];
            else
              *fp++ = 0.0;
#endif
            break;
          case IO_SUB_GASDUSTMETALMAXRAD:
#if defined(GFM_DUST) && defined(SUBFIND)
            if(SubGroup[pindex].SubMassInMaxRadType[0] > 0.0)
              *fp++ = SubGroup[pindex].GasMassDustMetallicityMaxRad / SubGroup[pindex].SubMassInMaxRadType[0];
            else
              *fp++ = 0.0;
#endif
            break;
          case IO_SUB_GASDUSTMETALSFR:
#if defined(GFM_DUST) && defined(SUBFIND)
            if(SubGroup[pindex].GasMassSfr > 0.0)
              *fp++ = SubGroup[pindex].GasMassDustMetallicitySfr / SubGroup[pindex].GasMassSfr;
            else
              *fp++ = 0.0;
#endif
            break;
          case IO_SUB_GASDUSTMETALSFRWEIGHTED:
#if defined(GFM_DUST) && defined(SUBFIND)
            if(SubGroup[pindex].Sfr > 0.0)
              *fp++ = SubGroup[pindex].GasMassDustMetallicitySfrWeighted / SubGroup[pindex].Sfr;
            else
              *fp++ = 0.0;
#endif
            break;
          case IO_SUB_BHMASS:
#if defined(BLACK_HOLES) && defined(SUBFIND)
            *fp++ = SubGroup[pindex].BH_Mass;
#endif
            break;
          case IO_SUB_BHMDOT:
#if defined(BLACK_HOLES) && defined(SUBFIND)
            *fp++ = SubGroup[pindex].BH_Mdot;
#endif
            break;
          case IO_SUB_WINDMASS:
#if defined(GFM_WINDS) && defined(SUBFIND)
            *fp++ = SubGroup[pindex].WindMass;
#endif
            break;
          case IO_SUB_H2MASS:
#if defined(SUBFIND_MEASURE_H2MASS) && defined(SUBFIND)
            *fp++ = SubGroup[pindex].H2_Mass;
#endif
            break;
          case IO_SUB_STELLARPHOTOMETRICS:
#if defined(GFM_STELLAR_PHOTOMETRICS) && defined(SUBFIND)
            *fp++ = SubGroup[pindex].Magnitude_U;
            *fp++ = SubGroup[pindex].Magnitude_B;
            *fp++ = SubGroup[pindex].Magnitude_V;
            *fp++ = SubGroup[pindex].Magnitude_K;
            *fp++ = SubGroup[pindex].Magnitude_g;
            *fp++ = SubGroup[pindex].Magnitude_r;
            *fp++ = SubGroup[pindex].Magnitude_i;
            *fp++ = SubGroup[pindex].Magnitude_z;
#endif
            break;
          case IO_SUB_STELLARPHOTOMETRICSRAD:
#if defined(GFM_STELLAR_PHOTOMETRICS) && defined(SUBFIND)
            *fp++ = SubGroup[pindex].SurfaceBrightnessLimitRad;
#endif
            break;
          case IO_SUB_STELLARPHOTOMETRICSMASSINRAD:
#if defined(GFM_STELLAR_PHOTOMETRICS) && defined(SUBFIND)
            *fp++ = SubGroup[pindex].SubMassInPhotRad;
#endif
            break;
          case IO_FOFSUB_IDS:
#ifdef FOF_STOREIDS
            *idp++ = ID_list[pindex].ID;
#endif
            break;

          case IO_FOF_LASTENTRY:
            terminate("should not be reached");
            break;
        }
    }
}

/*! \brief Associates the output variable blocknumber with its name.
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *  \param[out] label Name of field.
 *
 *  \return void
 */
void fof_subfind_get_dataset_name(const enum fof_subfind_iofields blocknr, char *const label)
{
  switch(blocknr)
    {
      case IO_FOF_LEN:
        strcpy(label, "GroupLen");
        break;
      case IO_FOF_MTOT:
        strcpy(label, "GroupMass");
        break;
      case IO_FOF_POS:
        strcpy(label, "GroupPos");
        break;
      case IO_FOF_CM:
        strcpy(label, "GroupCM");
        break;
      case IO_FOF_POSMINPOT:
        strcpy(label, "GroupPosMinPot");
        break;
      case IO_FOF_VEL:
        strcpy(label, "GroupVel");
        break;
      case IO_FOF_LENTYPE:
        strcpy(label, "GroupLenType");
        break;
      case IO_FOF_MASSTYPE:
        strcpy(label, "GroupMassType");
        break;
      case IO_FOF_SFR:
        strcpy(label, "GroupSFR");
        break;
      case IO_FOF_GASMETAL:
        strcpy(label, "GroupGasMetallicity");
        break;
      case IO_FOF_STARMETAL:
        strcpy(label, "GroupStarMetallicity");
        break;
      case IO_FOF_GASMETALELEMENTS:
        strcpy(label, "GroupGasMetalFractions");
        break;
      case IO_FOF_STARMETALELEMENTS:
        strcpy(label, "GroupStarMetalFractions");
        break;
      case IO_FOF_GASDUSTMETAL:
        strcpy(label, "GroupGasDustMetallicity");
        break;
      case IO_FOF_BHMASS:
        strcpy(label, "GroupBHMass");
        break;
      case IO_FOF_BHMDOT:
        strcpy(label, "GroupBHMdot");
        break;
      case IO_FOF_WINDMASS:
        strcpy(label, "GroupWindMass");
        break;
      case IO_FOF_M_MEAN200:
        strcpy(label, "Group_M_Mean200");
        break;
      case IO_FOF_R_MEAN200:
        strcpy(label, "Group_R_Mean200");
        break;

#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_FOF_J_MEAN200:
        strcpy(label, "Group_J_Mean200");
        break;
      case IO_FOF_JDM_MEAN200:
        strcpy(label, "Group_Jdm_Mean200");
        break;
      case IO_FOF_JGAS_MEAN200:
        strcpy(label, "Group_Jgas_Mean200");
        break;
      case IO_FOF_JSTARS_MEAN200:
        strcpy(label, "Group_Jstars_Mean200");
        break;
      case IO_FOF_MASSTYPE_MEAN200:
        strcpy(label, "Group_MassType_Mean200");
        break;
      case IO_FOF_LENTYPE_MEAN200:
        strcpy(label, "Group_LenType_Mean200");
        break;
      case IO_FOF_CMFRAC_MEAN200:
        strcpy(label, "Group_CMFrac_Mean200");
        break;
      case IO_FOF_CMFRACTYPE_MEAN200:
        strcpy(label, "Group_CMFracType_Mean200");
        break;
      case IO_FOF_J_CRIT200:
        strcpy(label, "Group_J_Crit200");
        break;
      case IO_FOF_JDM_CRIT200:
        strcpy(label, "Group_Jdm_Crit200");
        break;
      case IO_FOF_JGAS_CRIT200:
        strcpy(label, "Group_Jgas_Crit200");
        break;
      case IO_FOF_JSTARS_CRIT200:
        strcpy(label, "Group_Jstars_Crit200");
        break;
      case IO_FOF_MASSTYPE_CRIT200:
        strcpy(label, "Group_MassType_Crit200");
        break;
      case IO_FOF_LENTYPE_CRIT200:
        strcpy(label, "Group_LenType_Crit200");
        break;
      case IO_FOF_CMFRAC_CRIT200:
        strcpy(label, "Group_CMFrac_Crit200");
        break;
      case IO_FOF_CMFRACTYPE_CRIT200:
        strcpy(label, "Group_CMFracType_Crit200");
        break;
      case IO_FOF_J_CRIT500:
        strcpy(label, "Group_J_Crit500");
        break;
      case IO_FOF_JDM_CRIT500:
        strcpy(label, "Group_Jdm_Crit500");
        break;
      case IO_FOF_JGAS_CRIT500:
        strcpy(label, "Group_Jgas_Crit500");
        break;
      case IO_FOF_JSTARS_CRIT500:
        strcpy(label, "Group_Jstars_Crit500");
        break;
      case IO_FOF_MASSTYPE_CRIT500:
        strcpy(label, "Group_MassType_Crit500");
        break;
      case IO_FOF_LENTYPE_CRIT500:
        strcpy(label, "Group_LenType_Crit500");
        break;
      case IO_FOF_CMFRAC_CRIT500:
        strcpy(label, "Group_CMFrac_Crit500");
        break;
      case IO_FOF_CMFRACTYPE_CRIT500:
        strcpy(label, "Group_CMFracType_Crit500");
        break;
      case IO_FOF_J_TOPHAT200:
        strcpy(label, "Group_J_TopHat200");
        break;
      case IO_FOF_JDM_TOPHAT200:
        strcpy(label, "Group_Jdm_TopHat200");
        break;
      case IO_FOF_JGAS_TOPHAT200:
        strcpy(label, "Group_Jgas_TopHat200");
        break;
      case IO_FOF_JSTARS_TOPHAT200:
        strcpy(label, "Group_Jstars_TopHat200");
        break;
      case IO_FOF_MASSTYPE_TOPHAT200:
        strcpy(label, "Group_MassType_TopHat200");
        break;
      case IO_FOF_LENTYPE_TOPHAT200:
        strcpy(label, "Group_LenType_TopHat200");
        break;
      case IO_FOF_CMFRAC_TOPHAT200:
        strcpy(label, "Group_CMFrac_TopHat200");
        break;
      case IO_FOF_CMFRACTYPE_TOPHAT200:
        strcpy(label, "Group_CMFracType_TopHat200");
        break;

      case IO_FOF_EPOT_CRIT200:
        strcpy(label, "Group_Epot_Crit200");
        break;
      case IO_FOF_EKIN_CRIT200:
        strcpy(label, "Group_Ekin_Crit200");
        break;
      case IO_FOF_ETHR_CRIT200:
        strcpy(label, "Group_Ethr_Crit200");
        break;
      case IO_FOF_EPOT_MEAN200:
        strcpy(label, "Group_Epot_Mean200");
        break;
      case IO_FOF_EKIN_MEAN200:
        strcpy(label, "Group_Ekin_Mean200");
        break;
      case IO_FOF_ETHR_MEAN200:
        strcpy(label, "Group_Ethr_Mean200");
        break;
      case IO_FOF_EPOT_TOPHAT200:
        strcpy(label, "Group_Epot_TopHat200");
        break;
      case IO_FOF_EKIN_TOPHAT200:
        strcpy(label, "Group_Ekin_TopHat200");
        break;
      case IO_FOF_ETHR_TOPHAT200:
        strcpy(label, "Group_Ethr_TopHat200");
        break;
      case IO_FOF_EPOT_CRIT500:
        strcpy(label, "Group_Epot_Crit500");
        break;
      case IO_FOF_EKIN_CRIT500:
        strcpy(label, "Group_Ekin_Crit500");
        break;
      case IO_FOF_ETHR_CRIT500:
        strcpy(label, "Group_Ethr_Crit500");
        break;

      case IO_FOF_J:
        strcpy(label, "Group_J");
        break;
      case IO_FOF_JDM:
        strcpy(label, "Group_Jdm");
        break;
      case IO_FOF_JGAS:
        strcpy(label, "Group_Jgas");
        break;
      case IO_FOF_JSTARS:
        strcpy(label, "Group_Jstars");
        break;
      case IO_FOF_CMFRAC:
        strcpy(label, "Group_CMFrac");
        break;
      case IO_FOF_CMFRACTYPE:
        strcpy(label, "Group_CMFracType");
        break;
      case IO_FOF_EKIN:
        strcpy(label, "GroupEkin");
        break;
      case IO_FOF_ETHR:
        strcpy(label, "GroupEthr");
        break;
      case IO_FOF_EPOT:
        strcpy(label, "GroupEpot");
        break;
      case IO_SUB_EKIN:
        strcpy(label, "SubhaloEkin");
        break;
      case IO_SUB_ETHR:
        strcpy(label, "SubhaloEthr");
        break;
      case IO_SUB_EPOT:
        strcpy(label, "SubhaloEpot");
        break;
      case IO_SUB_J:
        strcpy(label, "Subhalo_J");
        break;
      case IO_SUB_JDM:
        strcpy(label, "Subhalo_Jdm");
        break;
      case IO_SUB_JGAS:
        strcpy(label, "Subhalo_Jgas");
        break;
      case IO_SUB_JSTARS:
        strcpy(label, "Subhalo_Jstars");
        break;
      case IO_SUB_JINHALFRAD:
        strcpy(label, "Subhalo_JInHalfRad");
        break;
      case IO_SUB_JDMINHALFRAD:
        strcpy(label, "Subhalo_JdmInHalfRad");
        break;
      case IO_SUB_JGASINHALFRAD:
        strcpy(label, "Subhalo_JgasInHalfRad");
        break;
      case IO_SUB_JSTARSINHALFRAD:
        strcpy(label, "Subhalo_JstarsInHalfRad");
        break;
      case IO_SUB_JINRAD:
        strcpy(label, "Subhalo_JInRad");
        break;
      case IO_SUB_JDMINRAD:
        strcpy(label, "Subhalo_JdmInRad");
        break;
      case IO_SUB_JGASINRAD:
        strcpy(label, "Subhalo_JgasInRad");
        break;
      case IO_SUB_JSTARSINRAD:
        strcpy(label, "Subhalo_JstarsInRad");
        break;
      case IO_SUB_CMFRAC:
        strcpy(label, "Subhalo_CMFrac");
        break;
      case IO_SUB_CMFRACTYPE:
        strcpy(label, "Subhalo_CMFracType");
        break;
      case IO_SUB_CMFRACINHALFRAD:
        strcpy(label, "Subhalo_CMFracInHalfRad");
        break;
      case IO_SUB_CMFRACTYPEINHALFRAD:
        strcpy(label, "Subhalo_CMFracTypeInHalfRad");
        break;
      case IO_SUB_CMFRACINRAD:
        strcpy(label, "Subhalo_CMFracInRad");
        break;
      case IO_SUB_CMFRACTYPEINRAD:
        strcpy(label, "Subhalo_CMFracTypeInRad");
        break;
#endif

      case IO_FOF_M_CRIT200:
        strcpy(label, "Group_M_Crit200");
        break;
      case IO_FOF_R_CRIT200:
        strcpy(label, "Group_R_Crit200");
        break;
      case IO_FOF_M_CRIT500:
        strcpy(label, "Group_M_Crit500");
        break;
      case IO_FOF_R_CRIT500:
        strcpy(label, "Group_R_Crit500");
        break;
      case IO_FOF_M_TOPHAT200:
        strcpy(label, "Group_M_TopHat200");
        break;
      case IO_FOF_R_TOPHAT200:
        strcpy(label, "Group_R_TopHat200");
        break;
      case IO_FOF_NSUBS:
        strcpy(label, "GroupNsubs");
        break;
      case IO_FOF_FIRSTSUB:
        strcpy(label, "GroupFirstSub");
        break;
      case IO_FOF_FUZZOFFTYPE:
        strcpy(label, "GroupFuzzOffsetType");
        break;
      case IO_FOF_RADIOLUM:
        strcpy(label, "GroupRadioLuminosity");
        break;
      case IO_FOF_XRAYLUM:
        strcpy(label, "GroupXrayLuminosity");
        break;
      case IO_FOF_DENSLVEC:
        strcpy(label, "GroupDensGasAngMomentum");
        break;
      case IO_SUB_LEN:
        strcpy(label, "SubhaloLen");
        break;
      case IO_SUB_MTOT:
        strcpy(label, "SubhaloMass");
        break;
      case IO_SUB_POS:
        strcpy(label, "SubhaloPos");
        break;
      case IO_SUB_VEL:
        strcpy(label, "SubhaloVel");
        break;
      case IO_SUB_LENTYPE:
        strcpy(label, "SubhaloLenType");
        break;
      case IO_SUB_MASSTYPE:
        strcpy(label, "SubhaloMassType");
        break;
      case IO_SUB_CM:
        strcpy(label, "SubhaloCM");
        break;
      case IO_SUB_SPIN:
        strcpy(label, "SubhaloSpin");
        break;
      case IO_SUB_VELDISP:
        strcpy(label, "SubhaloVelDisp");
        break;
      case IO_SUB_VMAX:
        strcpy(label, "SubhaloVmax");
        break;
      case IO_SUB_VMAXRAD:
        strcpy(label, "SubhaloVmaxRad");
        break;
      case IO_SUB_HALFMASSRAD:
        strcpy(label, "SubhaloHalfmassRad");
        break;
      case IO_SUB_HALFMASSRADTYPE:
        strcpy(label, "SubhaloHalfmassRadType");
        break;
      case IO_SUB_MASSINRAD:
        strcpy(label, "SubhaloMassInRad");
        break;
      case IO_SUB_MASSINHALFRAD:
        strcpy(label, "SubhaloMassInHalfRad");
        break;
      case IO_SUB_MASSINMAXRAD:
        strcpy(label, "SubhaloMassInMaxRad");
        break;
      case IO_SUB_MASSINRADTYPE:
        strcpy(label, "SubhaloMassInRadType");
        break;
      case IO_SUB_MASSINHALFRADTYPE:
        strcpy(label, "SubhaloMassInHalfRadType");
        break;
      case IO_SUB_MASSINMAXRADTYPE:
        strcpy(label, "SubhaloMassInMaxRadType");
        break;
      case IO_SUB_IDMOSTBOUND:
        strcpy(label, "SubhaloIDMostbound");
        break;
      case IO_SUB_GRNR:
        strcpy(label, "SubhaloGrNr");
        break;
      case IO_SUB_PARENT:
        strcpy(label, "SubhaloParent");
        break;
      case IO_SUB_BFLD_HALO:
        strcpy(label, "SubhaloBfldHalo");
        break;
      case IO_SUB_BFLD_DISK:
        strcpy(label, "SubhaloBfldDisk");
        break;
      case IO_SUB_SFR:
        strcpy(label, "SubhaloSFR");
        break;
      case IO_SUB_SFRINRAD:
        strcpy(label, "SubhaloSFRinRad");
        break;
      case IO_SUB_SFRINHALFRAD:
        strcpy(label, "SubhaloSFRinHalfRad");
        break;
      case IO_SUB_SFRINMAXRAD:
        strcpy(label, "SubhaloSFRinMaxRad");
        break;
      case IO_SUB_GASMETAL:
        strcpy(label, "SubhaloGasMetallicity");
        break;
      case IO_SUB_GASMETALHALFRAD:
        strcpy(label, "SubhaloGasMetallicityHalfRad");
        break;
      case IO_SUB_GASMETALMAXRAD:
        strcpy(label, "SubhaloGasMetallicityMaxRad");
        break;
      case IO_SUB_GASMETALSFR:
        strcpy(label, "SubhaloGasMetallicitySfr");
        break;
      case IO_SUB_GASMETALSFRWEIGHTED:
        strcpy(label, "SubhaloGasMetallicitySfrWeighted");
        break;
      case IO_SUB_STARMETAL:
        strcpy(label, "SubhaloStarMetallicity");
        break;
      case IO_SUB_STARMETALHALFRAD:
        strcpy(label, "SubhaloStarMetallicityHalfRad");
        break;
      case IO_SUB_STARMETALMAXRAD:
        strcpy(label, "SubhaloStarMetallicityMaxRad");
        break;
      case IO_SUB_GASMETALELEMENTS:
        strcpy(label, "SubhaloGasMetalFractions");
        break;
      case IO_SUB_GASMETALELEMENTSHALFRAD:
        strcpy(label, "SubhaloGasMetalFractionsHalfRad");
        break;
      case IO_SUB_GASMETALELEMENTSMAXRAD:
        strcpy(label, "SubhaloGasMetalFractionsMaxRad");
        break;
      case IO_SUB_GASMETALELEMENTSSFR:
        strcpy(label, "SubhaloGasMetalFractionsSfr");
        break;
      case IO_SUB_GASMETALELEMENTSSFRWEIGHTED:
        strcpy(label, "SubhaloGasMetalFractionsSfrWeighted");
        break;
      case IO_SUB_STARMETALELEMENTS:
        strcpy(label, "SubhaloStarMetalFractions");
        break;
      case IO_SUB_STARMETALELEMENTSHALFRAD:
        strcpy(label, "SubhaloStarMetalFractionsHalfRad");
        break;
      case IO_SUB_STARMETALELEMENTSMAXRAD:
        strcpy(label, "SubhaloStarMetalFractionsMaxRad");
        break;
      case IO_SUB_GASDUSTMETAL:
        strcpy(label, "SubhaloGasDustMetallicity");
        break;
      case IO_SUB_GASDUSTMETALHALFRAD:
        strcpy(label, "SubhaloGasDustMetallicityHalfRad");
        break;
      case IO_SUB_GASDUSTMETALMAXRAD:
        strcpy(label, "SubhaloGasDustMetallicityMaxRad");
        break;
      case IO_SUB_GASDUSTMETALSFR:
        strcpy(label, "SubhaloGasDustMetallicitySfr");
        break;
      case IO_SUB_GASDUSTMETALSFRWEIGHTED:
        strcpy(label, "SubhaloGasDustMetallicitySfrWeighted");
        break;
      case IO_SUB_BHMASS:
        strcpy(label, "SubhaloBHMass");
        break;
      case IO_SUB_BHMDOT:
        strcpy(label, "SubhaloBHMdot");
        break;
      case IO_SUB_WINDMASS:
        strcpy(label, "SubhaloWindMass");
        break;
      case IO_SUB_H2MASS:
        strcpy(label, "SubhaloH2Mass");
        break;
      case IO_SUB_STELLARPHOTOMETRICS:
        strcpy(label, "SubhaloStellarPhotometrics");
        break;
      case IO_SUB_STELLARPHOTOMETRICSRAD:
        strcpy(label, "SubhaloStellarPhotometricsRad");
        break;
      case IO_SUB_STELLARPHOTOMETRICSMASSINRAD:
        strcpy(label, "SubhaloStellarPhotometricsMassInRad");
        break;
      case IO_FOFSUB_IDS:
        strcpy(label, "ID");
        break;

      case IO_FOF_LASTENTRY:
        terminate("should not be reached");
        break;
    }
}

/*! \brief Is this output field a group or subhalo property?
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *
 *  \return 0: group property; 1 subhalo property; 2: both (unused)
 */
int fof_subfind_get_dataset_group(const enum fof_subfind_iofields blocknr)
{
  switch(blocknr)
    {
      case IO_FOF_LEN:
      case IO_FOF_MTOT:
      case IO_FOF_POS:
      case IO_FOF_CM:
      case IO_FOF_POSMINPOT:
      case IO_FOF_VEL:
      case IO_FOF_LENTYPE:
      case IO_FOF_MASSTYPE:
      case IO_FOF_SFR:
      case IO_FOF_GASMETAL:
      case IO_FOF_STARMETAL:
      case IO_FOF_GASMETALELEMENTS:
      case IO_FOF_STARMETALELEMENTS:
      case IO_FOF_GASDUSTMETAL:
      case IO_FOF_BHMASS:
      case IO_FOF_BHMDOT:
      case IO_FOF_WINDMASS:
      case IO_FOF_M_MEAN200:
      case IO_FOF_R_MEAN200:
      case IO_FOF_M_CRIT200:
      case IO_FOF_R_CRIT200:
      case IO_FOF_M_TOPHAT200:
      case IO_FOF_R_TOPHAT200:
      case IO_FOF_M_CRIT500:
      case IO_FOF_R_CRIT500:
      case IO_FOF_NSUBS:
      case IO_FOF_FIRSTSUB:
      case IO_FOF_FUZZOFFTYPE:
      case IO_FOF_RADIOLUM:
      case IO_FOF_XRAYLUM:
      case IO_FOF_DENSLVEC:
#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_FOF_J_MEAN200:
      case IO_FOF_JDM_MEAN200:
      case IO_FOF_JGAS_MEAN200:
      case IO_FOF_JSTARS_MEAN200:
      case IO_FOF_MASSTYPE_MEAN200:
      case IO_FOF_LENTYPE_MEAN200:
      case IO_FOF_CMFRAC_MEAN200:
      case IO_FOF_CMFRACTYPE_MEAN200:
      case IO_FOF_J_CRIT200:
      case IO_FOF_JDM_CRIT200:
      case IO_FOF_JGAS_CRIT200:
      case IO_FOF_JSTARS_CRIT200:
      case IO_FOF_MASSTYPE_CRIT200:
      case IO_FOF_LENTYPE_CRIT200:
      case IO_FOF_CMFRAC_CRIT200:
      case IO_FOF_CMFRACTYPE_CRIT200:
      case IO_FOF_J_TOPHAT200:
      case IO_FOF_JDM_TOPHAT200:
      case IO_FOF_JGAS_TOPHAT200:
      case IO_FOF_JSTARS_TOPHAT200:
      case IO_FOF_MASSTYPE_TOPHAT200:
      case IO_FOF_LENTYPE_TOPHAT200:
      case IO_FOF_CMFRAC_TOPHAT200:
      case IO_FOF_CMFRACTYPE_TOPHAT200:
      case IO_FOF_J_CRIT500:
      case IO_FOF_JDM_CRIT500:
      case IO_FOF_JGAS_CRIT500:
      case IO_FOF_JSTARS_CRIT500:
      case IO_FOF_MASSTYPE_CRIT500:
      case IO_FOF_LENTYPE_CRIT500:
      case IO_FOF_CMFRAC_CRIT500:
      case IO_FOF_CMFRACTYPE_CRIT500:
      case IO_FOF_J:
      case IO_FOF_JDM:
      case IO_FOF_JGAS:
      case IO_FOF_JSTARS:
      case IO_FOF_CMFRAC:
      case IO_FOF_CMFRACTYPE:
      case IO_FOF_EKIN:
      case IO_FOF_ETHR:
      case IO_FOF_EPOT:
      case IO_FOF_EPOT_CRIT200:
      case IO_FOF_EKIN_CRIT200:
      case IO_FOF_ETHR_CRIT200:
      case IO_FOF_EPOT_MEAN200:
      case IO_FOF_EKIN_MEAN200:
      case IO_FOF_ETHR_MEAN200:
      case IO_FOF_EPOT_TOPHAT200:
      case IO_FOF_EKIN_TOPHAT200:
      case IO_FOF_ETHR_TOPHAT200:
      case IO_FOF_EPOT_CRIT500:
      case IO_FOF_EKIN_CRIT500:
      case IO_FOF_ETHR_CRIT500:
#endif

        return 0;

      case IO_SUB_LEN:
      case IO_SUB_MTOT:
      case IO_SUB_POS:
      case IO_SUB_VEL:
      case IO_SUB_LENTYPE:
      case IO_SUB_MASSTYPE:
      case IO_SUB_CM:
      case IO_SUB_SPIN:
      case IO_SUB_VELDISP:
      case IO_SUB_VMAX:
      case IO_SUB_VMAXRAD:
      case IO_SUB_HALFMASSRAD:
      case IO_SUB_HALFMASSRADTYPE:
      case IO_SUB_MASSINRAD:
      case IO_SUB_MASSINHALFRAD:
      case IO_SUB_MASSINMAXRAD:
      case IO_SUB_MASSINRADTYPE:
      case IO_SUB_MASSINHALFRADTYPE:
      case IO_SUB_MASSINMAXRADTYPE:
      case IO_SUB_IDMOSTBOUND:
      case IO_SUB_GRNR:
      case IO_SUB_PARENT:
      case IO_SUB_BFLD_HALO:
      case IO_SUB_BFLD_DISK:
      case IO_SUB_SFR:
      case IO_SUB_SFRINRAD:
      case IO_SUB_SFRINHALFRAD:
      case IO_SUB_SFRINMAXRAD:
      case IO_SUB_GASMETAL:
      case IO_SUB_GASMETALHALFRAD:
      case IO_SUB_GASMETALMAXRAD:
      case IO_SUB_GASMETALSFR:
      case IO_SUB_GASMETALSFRWEIGHTED:
      case IO_SUB_STARMETAL:
      case IO_SUB_STARMETALHALFRAD:
      case IO_SUB_STARMETALMAXRAD:
      case IO_SUB_GASMETALELEMENTS:
      case IO_SUB_GASMETALELEMENTSHALFRAD:
      case IO_SUB_GASMETALELEMENTSMAXRAD:
      case IO_SUB_GASMETALELEMENTSSFR:
      case IO_SUB_GASMETALELEMENTSSFRWEIGHTED:
      case IO_SUB_STARMETALELEMENTS:
      case IO_SUB_STARMETALELEMENTSHALFRAD:
      case IO_SUB_STARMETALELEMENTSMAXRAD:
      case IO_SUB_GASDUSTMETAL:
      case IO_SUB_GASDUSTMETALHALFRAD:
      case IO_SUB_GASDUSTMETALMAXRAD:
      case IO_SUB_GASDUSTMETALSFR:
      case IO_SUB_GASDUSTMETALSFRWEIGHTED:
      case IO_SUB_BHMASS:
      case IO_SUB_BHMDOT:
      case IO_SUB_WINDMASS:
      case IO_SUB_H2MASS:
      case IO_SUB_STELLARPHOTOMETRICS:
      case IO_SUB_STELLARPHOTOMETRICSRAD:
      case IO_SUB_STELLARPHOTOMETRICSMASSINRAD:
#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_SUB_EKIN:
      case IO_SUB_ETHR:
      case IO_SUB_EPOT:
      case IO_SUB_J:
      case IO_SUB_JDM:
      case IO_SUB_JGAS:
      case IO_SUB_JSTARS:
      case IO_SUB_JINHALFRAD:
      case IO_SUB_JDMINHALFRAD:
      case IO_SUB_JGASINHALFRAD:
      case IO_SUB_JSTARSINHALFRAD:
      case IO_SUB_JINRAD:
      case IO_SUB_JDMINRAD:
      case IO_SUB_JGASINRAD:
      case IO_SUB_JSTARSINRAD:
      case IO_SUB_CMFRAC:
      case IO_SUB_CMFRACTYPE:
      case IO_SUB_CMFRACINHALFRAD:
      case IO_SUB_CMFRACTYPEINHALFRAD:
      case IO_SUB_CMFRACINRAD:
      case IO_SUB_CMFRACTYPEINRAD:
#endif
        return 1;

      case IO_FOFSUB_IDS:
        return 2;

      case IO_FOF_LASTENTRY:
        terminate("reached last entry in switch - strange.");
        break;
    }

  terminate("reached end of function - this should not happen");
  return 0;
}

/*! \brief Returns number of particles of specific field.
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *
 *  \return Number of entries of this property.
 */
int fof_subfind_get_particles_in_block(enum fof_subfind_iofields blocknr)
{
  switch(blocknr)
    {
      case IO_FOF_LEN:
      case IO_FOF_MTOT:
      case IO_FOF_POS:
      case IO_FOF_CM:
      case IO_FOF_POSMINPOT:
      case IO_FOF_VEL:
      case IO_FOF_LENTYPE:
      case IO_FOF_MASSTYPE:
      case IO_FOF_SFR:
      case IO_FOF_GASMETAL:
      case IO_FOF_STARMETAL:
      case IO_FOF_GASMETALELEMENTS:
      case IO_FOF_STARMETALELEMENTS:
      case IO_FOF_GASDUSTMETAL:
      case IO_FOF_BHMASS:
      case IO_FOF_BHMDOT:
      case IO_FOF_WINDMASS:
      case IO_FOF_FUZZOFFTYPE:
      case IO_FOF_RADIOLUM:
      case IO_FOF_XRAYLUM:
      case IO_FOF_DENSLVEC:
        return catalogue_header.Ngroups;

      case IO_FOF_M_MEAN200:
      case IO_FOF_R_MEAN200:
      case IO_FOF_M_CRIT200:
      case IO_FOF_R_CRIT200:
      case IO_FOF_M_TOPHAT200:
      case IO_FOF_R_TOPHAT200:
      case IO_FOF_M_CRIT500:
      case IO_FOF_R_CRIT500:
      case IO_FOF_NSUBS:
      case IO_FOF_FIRSTSUB:

#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_FOF_J_MEAN200:
      case IO_FOF_JDM_MEAN200:
      case IO_FOF_JGAS_MEAN200:
      case IO_FOF_JSTARS_MEAN200:
      case IO_FOF_MASSTYPE_MEAN200:
      case IO_FOF_LENTYPE_MEAN200:
      case IO_FOF_CMFRAC_MEAN200:
      case IO_FOF_CMFRACTYPE_MEAN200:
      case IO_FOF_J_CRIT200:
      case IO_FOF_JDM_CRIT200:
      case IO_FOF_JGAS_CRIT200:
      case IO_FOF_JSTARS_CRIT200:
      case IO_FOF_MASSTYPE_CRIT200:
      case IO_FOF_LENTYPE_CRIT200:
      case IO_FOF_CMFRAC_CRIT200:
      case IO_FOF_CMFRACTYPE_CRIT200:
      case IO_FOF_J_TOPHAT200:
      case IO_FOF_JDM_TOPHAT200:
      case IO_FOF_JGAS_TOPHAT200:
      case IO_FOF_JSTARS_TOPHAT200:
      case IO_FOF_MASSTYPE_TOPHAT200:
      case IO_FOF_LENTYPE_TOPHAT200:
      case IO_FOF_CMFRAC_TOPHAT200:
      case IO_FOF_CMFRACTYPE_TOPHAT200:
      case IO_FOF_J_CRIT500:
      case IO_FOF_JDM_CRIT500:
      case IO_FOF_JGAS_CRIT500:
      case IO_FOF_JSTARS_CRIT500:
      case IO_FOF_MASSTYPE_CRIT500:
      case IO_FOF_LENTYPE_CRIT500:
      case IO_FOF_CMFRAC_CRIT500:
      case IO_FOF_CMFRACTYPE_CRIT500:
      case IO_FOF_J:
      case IO_FOF_JDM:
      case IO_FOF_JGAS:
      case IO_FOF_JSTARS:
      case IO_FOF_CMFRAC:
      case IO_FOF_CMFRACTYPE:
      case IO_FOF_EKIN:
      case IO_FOF_ETHR:
      case IO_FOF_EPOT:
      case IO_FOF_EPOT_CRIT200:
      case IO_FOF_EKIN_CRIT200:
      case IO_FOF_ETHR_CRIT200:
      case IO_FOF_EPOT_MEAN200:
      case IO_FOF_EKIN_MEAN200:
      case IO_FOF_ETHR_MEAN200:
      case IO_FOF_EPOT_TOPHAT200:
      case IO_FOF_EKIN_TOPHAT200:
      case IO_FOF_ETHR_TOPHAT200:
      case IO_FOF_EPOT_CRIT500:
      case IO_FOF_EKIN_CRIT500:
      case IO_FOF_ETHR_CRIT500:
#endif

#ifdef SUBFIND
        return catalogue_header.Ngroups;
#else
        return 0;
#endif

      case IO_SUB_LEN:
      case IO_SUB_MTOT:
      case IO_SUB_POS:
      case IO_SUB_VEL:
      case IO_SUB_LENTYPE:
      case IO_SUB_MASSTYPE:
      case IO_SUB_CM:
      case IO_SUB_SPIN:
      case IO_SUB_VELDISP:
      case IO_SUB_VMAX:
      case IO_SUB_VMAXRAD:
      case IO_SUB_HALFMASSRAD:
      case IO_SUB_HALFMASSRADTYPE:
      case IO_SUB_MASSINRAD:
      case IO_SUB_MASSINHALFRAD:
      case IO_SUB_MASSINMAXRAD:
      case IO_SUB_MASSINRADTYPE:
      case IO_SUB_MASSINHALFRADTYPE:
      case IO_SUB_MASSINMAXRADTYPE:
      case IO_SUB_IDMOSTBOUND:
      case IO_SUB_GRNR:
      case IO_SUB_PARENT:
      case IO_SUB_BFLD_HALO:
      case IO_SUB_BFLD_DISK:
      case IO_SUB_SFR:
      case IO_SUB_SFRINRAD:
      case IO_SUB_SFRINHALFRAD:
      case IO_SUB_SFRINMAXRAD:
      case IO_SUB_GASMETAL:
      case IO_SUB_GASMETALHALFRAD:
      case IO_SUB_GASMETALMAXRAD:
      case IO_SUB_GASMETALSFR:
      case IO_SUB_GASMETALSFRWEIGHTED:
      case IO_SUB_STARMETAL:
      case IO_SUB_STARMETALHALFRAD:
      case IO_SUB_STARMETALMAXRAD:
      case IO_SUB_GASMETALELEMENTS:
      case IO_SUB_GASMETALELEMENTSHALFRAD:
      case IO_SUB_GASMETALELEMENTSMAXRAD:
      case IO_SUB_GASMETALELEMENTSSFR:
      case IO_SUB_GASMETALELEMENTSSFRWEIGHTED:
      case IO_SUB_STARMETALELEMENTS:
      case IO_SUB_STARMETALELEMENTSHALFRAD:
      case IO_SUB_STARMETALELEMENTSMAXRAD:
      case IO_SUB_GASDUSTMETAL:
      case IO_SUB_GASDUSTMETALHALFRAD:
      case IO_SUB_GASDUSTMETALMAXRAD:
      case IO_SUB_GASDUSTMETALSFR:
      case IO_SUB_GASDUSTMETALSFRWEIGHTED:
      case IO_SUB_BHMASS:
      case IO_SUB_BHMDOT:
      case IO_SUB_WINDMASS:
      case IO_SUB_H2MASS:
      case IO_SUB_STELLARPHOTOMETRICS:
      case IO_SUB_STELLARPHOTOMETRICSRAD:
      case IO_SUB_STELLARPHOTOMETRICSMASSINRAD:

#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_SUB_EKIN:
      case IO_SUB_ETHR:
      case IO_SUB_EPOT:
      case IO_SUB_J:
      case IO_SUB_JDM:
      case IO_SUB_JGAS:
      case IO_SUB_JSTARS:
      case IO_SUB_JINHALFRAD:
      case IO_SUB_JDMINHALFRAD:
      case IO_SUB_JGASINHALFRAD:
      case IO_SUB_JSTARSINHALFRAD:
      case IO_SUB_JINRAD:
      case IO_SUB_JDMINRAD:
      case IO_SUB_JGASINRAD:
      case IO_SUB_JSTARSINRAD:
      case IO_SUB_CMFRAC:
      case IO_SUB_CMFRACTYPE:
      case IO_SUB_CMFRACINHALFRAD:
      case IO_SUB_CMFRACTYPEINHALFRAD:
      case IO_SUB_CMFRACINRAD:
      case IO_SUB_CMFRACTYPEINRAD:
#endif

#ifdef SUBFIND
        return catalogue_header.Nsubgroups;
#else
        return 0;
#endif

      case IO_FOFSUB_IDS:
        return catalogue_header.Nids;

      case IO_FOF_LASTENTRY:
        terminate("reached last entry in switch - strange.");
        break;
    }

  terminate("reached end of function - this should not happen");
  return 0;
}

/*! \brief Returns the number of elements per entry of a given property.
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *
 *  \return Number of values per element of the specified property.
 */
int fof_subfind_get_values_per_blockelement(const enum fof_subfind_iofields blocknr)
{
  int values = 0;

  switch(blocknr)
    {
      case IO_FOF_LEN:
      case IO_FOF_NSUBS:
      case IO_FOF_FIRSTSUB:
      case IO_SUB_LEN:
      case IO_SUB_GRNR:
      case IO_SUB_PARENT:
      case IO_FOF_MTOT:
      case IO_FOF_SFR:
      case IO_FOF_GASMETAL:
      case IO_FOF_STARMETAL:
      case IO_FOF_GASDUSTMETAL:
      case IO_FOF_BHMASS:
      case IO_FOF_BHMDOT:
      case IO_FOF_WINDMASS:
      case IO_FOF_M_MEAN200:
      case IO_FOF_R_MEAN200:
      case IO_FOF_M_CRIT200:
      case IO_FOF_R_CRIT200:
      case IO_FOF_M_TOPHAT200:
      case IO_FOF_R_TOPHAT200:
      case IO_FOF_M_CRIT500:
      case IO_FOF_R_CRIT500:
      case IO_FOF_RADIOLUM:
      case IO_FOF_XRAYLUM:
      case IO_SUB_MTOT:
      case IO_SUB_VELDISP:
      case IO_SUB_VMAX:
      case IO_SUB_VMAXRAD:
      case IO_SUB_HALFMASSRAD:
      case IO_SUB_MASSINRAD:
      case IO_SUB_MASSINHALFRAD:
      case IO_SUB_MASSINMAXRAD:
      case IO_SUB_IDMOSTBOUND:
      case IO_SUB_BFLD_HALO:
      case IO_SUB_BFLD_DISK:
      case IO_SUB_SFR:
      case IO_SUB_SFRINRAD:
      case IO_SUB_SFRINHALFRAD:
      case IO_SUB_SFRINMAXRAD:
      case IO_SUB_GASMETAL:
      case IO_SUB_GASMETALHALFRAD:
      case IO_SUB_GASMETALMAXRAD:
      case IO_SUB_GASMETALSFR:
      case IO_SUB_GASMETALSFRWEIGHTED:
      case IO_SUB_STARMETAL:
      case IO_SUB_STARMETALHALFRAD:
      case IO_SUB_STARMETALMAXRAD:
      case IO_SUB_GASDUSTMETAL:
      case IO_SUB_GASDUSTMETALHALFRAD:
      case IO_SUB_GASDUSTMETALMAXRAD:
      case IO_SUB_GASDUSTMETALSFR:
      case IO_SUB_GASDUSTMETALSFRWEIGHTED:
      case IO_SUB_BHMASS:
      case IO_SUB_BHMDOT:
      case IO_SUB_WINDMASS:
      case IO_SUB_H2MASS:
      case IO_SUB_STELLARPHOTOMETRICSRAD:
      case IO_SUB_STELLARPHOTOMETRICSMASSINRAD:
      case IO_FOFSUB_IDS:
#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_FOF_CMFRAC_MEAN200:
      case IO_FOF_CMFRAC_CRIT200:
      case IO_FOF_CMFRAC_TOPHAT200:
      case IO_FOF_CMFRAC_CRIT500:
      case IO_FOF_EPOT_CRIT200:
      case IO_FOF_EKIN_CRIT200:
      case IO_FOF_ETHR_CRIT200:
      case IO_FOF_EPOT_MEAN200:
      case IO_FOF_EKIN_MEAN200:
      case IO_FOF_ETHR_MEAN200:
      case IO_FOF_EPOT_TOPHAT200:
      case IO_FOF_EKIN_TOPHAT200:
      case IO_FOF_ETHR_TOPHAT200:
      case IO_FOF_EPOT_CRIT500:
      case IO_FOF_EKIN_CRIT500:
      case IO_FOF_ETHR_CRIT500:
      case IO_FOF_EKIN:
      case IO_FOF_ETHR:
      case IO_FOF_EPOT:
      case IO_SUB_EKIN:
      case IO_SUB_ETHR:
      case IO_SUB_EPOT:
      case IO_SUB_CMFRAC:
      case IO_SUB_CMFRACINHALFRAD:
      case IO_SUB_CMFRACINRAD:
      case IO_FOF_CMFRAC:
#endif
        values = 1;
        break;

      case IO_FOF_LENTYPE:
      case IO_SUB_LENTYPE:
      case IO_FOF_MASSTYPE:
      case IO_SUB_MASSTYPE:
      case IO_SUB_HALFMASSRADTYPE:
      case IO_SUB_MASSINRADTYPE:
      case IO_SUB_MASSINHALFRADTYPE:
      case IO_SUB_MASSINMAXRADTYPE:
      case IO_FOF_FUZZOFFTYPE:
#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_FOF_CMFRACTYPE:
      case IO_SUB_CMFRACTYPE:
      case IO_SUB_CMFRACTYPEINHALFRAD:
      case IO_SUB_CMFRACTYPEINRAD:
      case IO_FOF_LENTYPE_MEAN200:
      case IO_FOF_LENTYPE_CRIT200:
      case IO_FOF_LENTYPE_CRIT500:
      case IO_FOF_LENTYPE_TOPHAT200:
      case IO_FOF_MASSTYPE_MEAN200:
      case IO_FOF_MASSTYPE_CRIT200:
      case IO_FOF_MASSTYPE_CRIT500:
      case IO_FOF_MASSTYPE_TOPHAT200:
      case IO_FOF_CMFRACTYPE_MEAN200:
      case IO_FOF_CMFRACTYPE_CRIT200:
      case IO_FOF_CMFRACTYPE_CRIT500:
      case IO_FOF_CMFRACTYPE_TOPHAT200:
#endif
        values = NTYPES;
        break;

      case IO_FOF_POS:
      case IO_FOF_CM:
      case IO_FOF_POSMINPOT:
      case IO_FOF_VEL:
      case IO_SUB_POS:
      case IO_SUB_VEL:
      case IO_SUB_CM:
      case IO_SUB_SPIN:
      case IO_FOF_DENSLVEC:
#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_SUB_J:
      case IO_SUB_JDM:
      case IO_SUB_JGAS:
      case IO_SUB_JSTARS:
      case IO_SUB_JINHALFRAD:
      case IO_SUB_JDMINHALFRAD:
      case IO_SUB_JGASINHALFRAD:
      case IO_SUB_JSTARSINHALFRAD:
      case IO_SUB_JINRAD:
      case IO_SUB_JDMINRAD:
      case IO_SUB_JGASINRAD:
      case IO_SUB_JSTARSINRAD:
      case IO_FOF_J_MEAN200:
      case IO_FOF_JDM_MEAN200:
      case IO_FOF_JGAS_MEAN200:
      case IO_FOF_JSTARS_MEAN200:
      case IO_FOF_J_CRIT200:
      case IO_FOF_JDM_CRIT200:
      case IO_FOF_JGAS_CRIT200:
      case IO_FOF_JSTARS_CRIT200:
      case IO_FOF_J_TOPHAT200:
      case IO_FOF_JDM_TOPHAT200:
      case IO_FOF_JGAS_TOPHAT200:
      case IO_FOF_JSTARS_TOPHAT200:
      case IO_FOF_J_CRIT500:
      case IO_FOF_JDM_CRIT500:
      case IO_FOF_JGAS_CRIT500:
      case IO_FOF_JSTARS_CRIT500:
      case IO_FOF_J:
      case IO_FOF_JDM:
      case IO_FOF_JGAS:
      case IO_FOF_JSTARS:
#endif
        values = 3;
        break;

      case IO_FOF_GASMETALELEMENTS:
      case IO_FOF_STARMETALELEMENTS:
      case IO_SUB_GASMETALELEMENTS:
      case IO_SUB_GASMETALELEMENTSHALFRAD:
      case IO_SUB_GASMETALELEMENTSMAXRAD:
      case IO_SUB_GASMETALELEMENTSSFR:
      case IO_SUB_GASMETALELEMENTSSFRWEIGHTED:
      case IO_SUB_STARMETALELEMENTS:
      case IO_SUB_STARMETALELEMENTSHALFRAD:
      case IO_SUB_STARMETALELEMENTSMAXRAD:
#if defined(GFM_STELLAR_EVOLUTION)
        values = GFM_N_CHEM_ELEMENTS;
#endif
        break;

      case IO_SUB_STELLARPHOTOMETRICS:
#ifdef GFM_STELLAR_PHOTOMETRICS
        values = GFM_STELLAR_PHOTOMETRICS_BANDS;
#endif
        break;

      case IO_FOF_LASTENTRY:
        terminate("reached last entry in switch - should not get here");
        break;
    }
  return values;
}

/*! \brief Returns the number of bytes per element of a given property.
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *
 *  \return Number of bytes per element for this property.
 */
int fof_subfind_get_bytes_per_blockelement(const enum fof_subfind_iofields blocknr)
{
  int bytes_per_blockelement = 0;

  switch(blocknr)
    {
      case IO_FOF_LEN:
      case IO_FOF_NSUBS:
      case IO_FOF_FIRSTSUB:
      case IO_SUB_LEN:
      case IO_SUB_GRNR:
      case IO_SUB_PARENT:
        bytes_per_blockelement = sizeof(int);
        break;

      case IO_FOF_LENTYPE:
      case IO_SUB_LENTYPE:
#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_FOF_LENTYPE_MEAN200:
      case IO_FOF_LENTYPE_CRIT200:
      case IO_FOF_LENTYPE_CRIT500:
      case IO_FOF_LENTYPE_TOPHAT200:
#endif
        bytes_per_blockelement = NTYPES * sizeof(int);
        break;

      case IO_FOF_MTOT:
      case IO_FOF_SFR:
      case IO_FOF_GASMETAL:
      case IO_FOF_STARMETAL:
      case IO_FOF_GASDUSTMETAL:
      case IO_FOF_BHMASS:
      case IO_FOF_BHMDOT:
      case IO_FOF_WINDMASS:
      case IO_FOF_M_MEAN200:
      case IO_FOF_R_MEAN200:
      case IO_FOF_M_CRIT200:
      case IO_FOF_R_CRIT200:
      case IO_FOF_M_TOPHAT200:
      case IO_FOF_R_TOPHAT200:
      case IO_FOF_M_CRIT500:
      case IO_FOF_R_CRIT500:
      case IO_FOF_RADIOLUM:
      case IO_FOF_XRAYLUM:
      case IO_SUB_MTOT:
      case IO_SUB_VELDISP:
      case IO_SUB_VMAX:
      case IO_SUB_VMAXRAD:
      case IO_SUB_HALFMASSRAD:
      case IO_SUB_MASSINRAD:
      case IO_SUB_MASSINHALFRAD:
      case IO_SUB_MASSINMAXRAD:
      case IO_SUB_BFLD_HALO:
      case IO_SUB_BFLD_DISK:
      case IO_SUB_SFR:
      case IO_SUB_SFRINRAD:
      case IO_SUB_SFRINHALFRAD:
      case IO_SUB_SFRINMAXRAD:
      case IO_SUB_GASMETAL:
      case IO_SUB_GASMETALHALFRAD:
      case IO_SUB_GASMETALMAXRAD:
      case IO_SUB_GASMETALSFR:
      case IO_SUB_GASMETALSFRWEIGHTED:
      case IO_SUB_STARMETAL:
      case IO_SUB_STARMETALHALFRAD:
      case IO_SUB_STARMETALMAXRAD:
      case IO_SUB_GASDUSTMETAL:
      case IO_SUB_GASDUSTMETALHALFRAD:
      case IO_SUB_GASDUSTMETALMAXRAD:
      case IO_SUB_GASDUSTMETALSFR:
      case IO_SUB_GASDUSTMETALSFRWEIGHTED:
      case IO_SUB_BHMASS:
      case IO_SUB_BHMDOT:
      case IO_SUB_WINDMASS:
      case IO_SUB_H2MASS:
      case IO_SUB_STELLARPHOTOMETRICSRAD:
      case IO_SUB_STELLARPHOTOMETRICSMASSINRAD:
#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_FOF_CMFRAC_MEAN200:
      case IO_FOF_CMFRAC_CRIT200:
      case IO_FOF_CMFRAC_TOPHAT200:
      case IO_FOF_CMFRAC_CRIT500:
      case IO_FOF_CMFRAC:
      case IO_FOF_EKIN:
      case IO_FOF_ETHR:
      case IO_FOF_EPOT:
      case IO_SUB_EKIN:
      case IO_SUB_ETHR:
      case IO_SUB_EPOT:
      case IO_SUB_CMFRAC:
      case IO_SUB_CMFRACINHALFRAD:
      case IO_SUB_CMFRACINRAD:
      case IO_FOF_EPOT_CRIT200:
      case IO_FOF_EKIN_CRIT200:
      case IO_FOF_ETHR_CRIT200:
      case IO_FOF_EPOT_MEAN200:
      case IO_FOF_EKIN_MEAN200:
      case IO_FOF_ETHR_MEAN200:
      case IO_FOF_EPOT_TOPHAT200:
      case IO_FOF_EKIN_TOPHAT200:
      case IO_FOF_ETHR_TOPHAT200:
      case IO_FOF_EPOT_CRIT500:
      case IO_FOF_EKIN_CRIT500:
      case IO_FOF_ETHR_CRIT500:
#endif
        bytes_per_blockelement = sizeof(MyOutputFloat);
        break;

      case IO_FOF_POS:
      case IO_FOF_CM:
      case IO_FOF_POSMINPOT:
      case IO_FOF_VEL:
      case IO_SUB_POS:
      case IO_SUB_VEL:
      case IO_SUB_CM:
      case IO_SUB_SPIN:
      case IO_FOF_DENSLVEC:
#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_SUB_J:
      case IO_SUB_JDM:
      case IO_SUB_JGAS:
      case IO_SUB_JSTARS:
      case IO_SUB_JINHALFRAD:
      case IO_SUB_JDMINHALFRAD:
      case IO_SUB_JGASINHALFRAD:
      case IO_SUB_JSTARSINHALFRAD:
      case IO_SUB_JINRAD:
      case IO_SUB_JDMINRAD:
      case IO_SUB_JGASINRAD:
      case IO_SUB_JSTARSINRAD:
      case IO_FOF_J_MEAN200:
      case IO_FOF_JDM_MEAN200:
      case IO_FOF_JGAS_MEAN200:
      case IO_FOF_JSTARS_MEAN200:
      case IO_FOF_J_CRIT200:
      case IO_FOF_JDM_CRIT200:
      case IO_FOF_JGAS_CRIT200:
      case IO_FOF_JSTARS_CRIT200:
      case IO_FOF_J_TOPHAT200:
      case IO_FOF_JDM_TOPHAT200:
      case IO_FOF_JGAS_TOPHAT200:
      case IO_FOF_JSTARS_TOPHAT200:
      case IO_FOF_J_CRIT500:
      case IO_FOF_JDM_CRIT500:
      case IO_FOF_JGAS_CRIT500:
      case IO_FOF_JSTARS_CRIT500:
      case IO_FOF_J:
      case IO_FOF_JDM:
      case IO_FOF_JGAS:
      case IO_FOF_JSTARS:
#endif
        bytes_per_blockelement = 3 * sizeof(MyOutputFloat);
        break;

      case IO_FOF_MASSTYPE:
      case IO_SUB_MASSTYPE:
      case IO_SUB_HALFMASSRADTYPE:
      case IO_SUB_MASSINRADTYPE:
      case IO_SUB_MASSINHALFRADTYPE:
      case IO_SUB_MASSINMAXRADTYPE:
#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_FOF_MASSTYPE_MEAN200:
      case IO_FOF_MASSTYPE_CRIT200:
      case IO_FOF_MASSTYPE_CRIT500:
      case IO_FOF_MASSTYPE_TOPHAT200:
      case IO_FOF_CMFRACTYPE_MEAN200:
      case IO_FOF_CMFRACTYPE_CRIT200:
      case IO_FOF_CMFRACTYPE_CRIT500:
      case IO_FOF_CMFRACTYPE_TOPHAT200:
      case IO_FOF_CMFRACTYPE:
      case IO_SUB_CMFRACTYPE:
      case IO_SUB_CMFRACTYPEINHALFRAD:
      case IO_SUB_CMFRACTYPEINRAD:
#endif
        bytes_per_blockelement = NTYPES * sizeof(MyOutputFloat);
        break;

      case IO_FOF_GASMETALELEMENTS:
      case IO_FOF_STARMETALELEMENTS:
      case IO_SUB_GASMETALELEMENTS:
      case IO_SUB_GASMETALELEMENTSHALFRAD:
      case IO_SUB_GASMETALELEMENTSMAXRAD:
      case IO_SUB_GASMETALELEMENTSSFR:
      case IO_SUB_GASMETALELEMENTSSFRWEIGHTED:
      case IO_SUB_STARMETALELEMENTS:
      case IO_SUB_STARMETALELEMENTSHALFRAD:
      case IO_SUB_STARMETALELEMENTSMAXRAD:
#if defined(GFM_STELLAR_EVOLUTION)
        bytes_per_blockelement = GFM_N_CHEM_ELEMENTS * sizeof(MyOutputFloat);
#endif
        break;

      case IO_SUB_STELLARPHOTOMETRICS:
#if defined(GFM_STELLAR_PHOTOMETRICS)
        bytes_per_blockelement = GFM_STELLAR_PHOTOMETRICS_BANDS * sizeof(MyOutputFloat);
#endif
        break;

      case IO_SUB_IDMOSTBOUND:
      case IO_FOFSUB_IDS:
        bytes_per_blockelement = sizeof(MyIDType);
        break;

      case IO_FOF_FUZZOFFTYPE:
        bytes_per_blockelement = NTYPES * sizeof(long long);
        break;

      case IO_FOF_LASTENTRY:
        terminate("reached last entry in switch - should not get here");
        break;
    }
  return bytes_per_blockelement;
}

/*! \brief Returns key for datatype of element of a given property.
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *
 *  \return Key for datatype: 0: int, 1: (output)float, 2: long long.
 */
int fof_subfind_get_datatype(const enum fof_subfind_iofields blocknr)
{
  int typekey = 0;

  switch(blocknr)
    {
      case IO_FOF_LEN:
      case IO_FOF_LENTYPE:
      case IO_FOF_NSUBS:
      case IO_FOF_FIRSTSUB:
      case IO_SUB_LEN:
      case IO_SUB_LENTYPE:
      case IO_SUB_GRNR:
      case IO_SUB_PARENT:
#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_FOF_LENTYPE_MEAN200:
      case IO_FOF_LENTYPE_CRIT200:
      case IO_FOF_LENTYPE_CRIT500:
      case IO_FOF_LENTYPE_TOPHAT200:
#endif
        typekey = 0; /* native int */
        break;

      case IO_FOF_MTOT:
      case IO_FOF_POS:
      case IO_FOF_CM:
      case IO_FOF_POSMINPOT:
      case IO_FOF_VEL:
      case IO_FOF_MASSTYPE:
      case IO_FOF_SFR:
      case IO_FOF_GASMETAL:
      case IO_FOF_STARMETAL:
      case IO_FOF_GASMETALELEMENTS:
      case IO_FOF_STARMETALELEMENTS:
      case IO_FOF_GASDUSTMETAL:
      case IO_FOF_BHMASS:
      case IO_FOF_BHMDOT:
      case IO_FOF_WINDMASS:
      case IO_FOF_M_MEAN200:
      case IO_FOF_R_MEAN200:
      case IO_FOF_M_CRIT200:
      case IO_FOF_R_CRIT200:
      case IO_FOF_M_TOPHAT200:
      case IO_FOF_R_TOPHAT200:
      case IO_FOF_M_CRIT500:
      case IO_FOF_R_CRIT500:
      case IO_FOF_RADIOLUM:
      case IO_FOF_XRAYLUM:
      case IO_SUB_MTOT:
      case IO_SUB_POS:
      case IO_SUB_VEL:
      case IO_SUB_MASSTYPE:
      case IO_SUB_CM:
      case IO_SUB_SPIN:
      case IO_SUB_VELDISP:
      case IO_SUB_VMAX:
      case IO_SUB_VMAXRAD:
      case IO_SUB_HALFMASSRAD:
      case IO_SUB_HALFMASSRADTYPE:
      case IO_SUB_MASSINRAD:
      case IO_SUB_MASSINHALFRAD:
      case IO_SUB_MASSINMAXRAD:
      case IO_SUB_MASSINRADTYPE:
      case IO_SUB_MASSINHALFRADTYPE:
      case IO_SUB_MASSINMAXRADTYPE:
      case IO_SUB_BFLD_HALO:
      case IO_SUB_BFLD_DISK:
      case IO_SUB_SFR:
      case IO_SUB_SFRINRAD:
      case IO_SUB_SFRINHALFRAD:
      case IO_SUB_SFRINMAXRAD:
      case IO_SUB_GASMETAL:
      case IO_SUB_GASMETALHALFRAD:
      case IO_SUB_GASMETALMAXRAD:
      case IO_SUB_GASMETALSFR:
      case IO_SUB_GASMETALSFRWEIGHTED:
      case IO_SUB_STARMETAL:
      case IO_SUB_STARMETALHALFRAD:
      case IO_SUB_STARMETALMAXRAD:
      case IO_SUB_GASMETALELEMENTS:
      case IO_SUB_GASMETALELEMENTSHALFRAD:
      case IO_SUB_GASMETALELEMENTSMAXRAD:
      case IO_SUB_GASMETALELEMENTSSFR:
      case IO_SUB_GASMETALELEMENTSSFRWEIGHTED:
      case IO_SUB_STARMETALELEMENTS:
      case IO_SUB_STARMETALELEMENTSHALFRAD:
      case IO_SUB_STARMETALELEMENTSMAXRAD:
      case IO_SUB_GASDUSTMETAL:
      case IO_SUB_GASDUSTMETALHALFRAD:
      case IO_SUB_GASDUSTMETALMAXRAD:
      case IO_SUB_GASDUSTMETALSFR:
      case IO_SUB_GASDUSTMETALSFRWEIGHTED:
      case IO_SUB_BHMASS:
      case IO_SUB_BHMDOT:
      case IO_SUB_WINDMASS:
      case IO_SUB_H2MASS:
      case IO_SUB_STELLARPHOTOMETRICS:
      case IO_SUB_STELLARPHOTOMETRICSRAD:
      case IO_SUB_STELLARPHOTOMETRICSMASSINRAD:
      case IO_FOF_DENSLVEC:
#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_FOF_MASSTYPE_MEAN200:
      case IO_FOF_MASSTYPE_CRIT200:
      case IO_FOF_MASSTYPE_CRIT500:
      case IO_FOF_MASSTYPE_TOPHAT200:
      case IO_FOF_J_MEAN200:
      case IO_FOF_JDM_MEAN200:
      case IO_FOF_JGAS_MEAN200:
      case IO_FOF_JSTARS_MEAN200:
      case IO_FOF_CMFRAC_MEAN200:
      case IO_FOF_CMFRACTYPE_MEAN200:
      case IO_FOF_J_CRIT200:
      case IO_FOF_JDM_CRIT200:
      case IO_FOF_JGAS_CRIT200:
      case IO_FOF_JSTARS_CRIT200:
      case IO_FOF_CMFRAC_CRIT200:
      case IO_FOF_CMFRACTYPE_CRIT200:
      case IO_FOF_J_TOPHAT200:
      case IO_FOF_JDM_TOPHAT200:
      case IO_FOF_JGAS_TOPHAT200:
      case IO_FOF_JSTARS_TOPHAT200:
      case IO_FOF_CMFRAC_TOPHAT200:
      case IO_FOF_CMFRACTYPE_TOPHAT200:
      case IO_FOF_J_CRIT500:
      case IO_FOF_JDM_CRIT500:
      case IO_FOF_JGAS_CRIT500:
      case IO_FOF_JSTARS_CRIT500:
      case IO_FOF_CMFRAC_CRIT500:
      case IO_FOF_CMFRACTYPE_CRIT500:
      case IO_FOF_J:
      case IO_FOF_JDM:
      case IO_FOF_JGAS:
      case IO_FOF_JSTARS:
      case IO_FOF_CMFRAC:
      case IO_FOF_CMFRACTYPE:
      case IO_FOF_EKIN:
      case IO_FOF_ETHR:
      case IO_FOF_EPOT:
      case IO_FOF_EPOT_CRIT200:
      case IO_FOF_EKIN_CRIT200:
      case IO_FOF_ETHR_CRIT200:
      case IO_FOF_EPOT_MEAN200:
      case IO_FOF_EKIN_MEAN200:
      case IO_FOF_ETHR_MEAN200:
      case IO_FOF_EPOT_TOPHAT200:
      case IO_FOF_EKIN_TOPHAT200:
      case IO_FOF_ETHR_TOPHAT200:
      case IO_FOF_EPOT_CRIT500:
      case IO_FOF_EKIN_CRIT500:
      case IO_FOF_ETHR_CRIT500:
      case IO_SUB_EKIN:
      case IO_SUB_ETHR:
      case IO_SUB_EPOT:
      case IO_SUB_J:
      case IO_SUB_JDM:
      case IO_SUB_JGAS:
      case IO_SUB_JSTARS:
      case IO_SUB_JINHALFRAD:
      case IO_SUB_JDMINHALFRAD:
      case IO_SUB_JGASINHALFRAD:
      case IO_SUB_JSTARSINHALFRAD:
      case IO_SUB_JINRAD:
      case IO_SUB_JDMINRAD:
      case IO_SUB_JGASINRAD:
      case IO_SUB_JSTARSINRAD:
      case IO_SUB_CMFRAC:
      case IO_SUB_CMFRACTYPE:
      case IO_SUB_CMFRACINHALFRAD:
      case IO_SUB_CMFRACTYPEINHALFRAD:
      case IO_SUB_CMFRACINRAD:
      case IO_SUB_CMFRACTYPEINRAD:
#endif
        typekey = 1; /* native MyOutputFloat */
        break;

      case IO_SUB_IDMOSTBOUND:
      case IO_FOFSUB_IDS:
#ifdef LONGIDS
        typekey = 2; /* native long long */
#else
        typekey = 0; /* native int */
#endif
        break;

      case IO_FOF_FUZZOFFTYPE:
        typekey = 2; /* native long long */
        break;

      case IO_FOF_LASTENTRY:
        terminate("should not be reached");
        break;
    }

  return typekey;
}

/*! \brief Determines if block is present in the current code configuration.
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *
 *  \return 0: not present; 1: present.
 */
int fof_subfind_blockpresent(const enum fof_subfind_iofields blocknr)
{
  int present = 0;

  switch(blocknr)
    {
      case IO_FOF_LEN:
      case IO_FOF_LENTYPE:
      case IO_FOF_MTOT:
      case IO_FOF_POS:
      case IO_FOF_CM:
      case IO_FOF_VEL:
      case IO_FOF_MASSTYPE:
        present = 1;
        break;

      case IO_FOF_SFR:
      case IO_SUB_SFR:
      case IO_SUB_SFRINRAD:
      case IO_SUB_SFRINHALFRAD:
      case IO_SUB_SFRINMAXRAD:
#ifdef USE_SFR
        present = 1;
#endif
        break;

      case IO_SUB_BFLD_HALO:
      case IO_SUB_BFLD_DISK:
#ifdef MHD
        present = 1;
#endif
        break;

      case IO_FOF_POSMINPOT:
#if defined(GFM_BIPOLAR_WINDS) && (GFM_BIPOLAR_WINDS == 3)
        present = 1;
#endif
        break;

      case IO_FOF_FUZZOFFTYPE:
#ifdef FOF_FUZZ_SORT_BY_NEAREST_GROUP
        present = 1;
#endif
        break;

      case IO_FOF_DENSLVEC:
#if defined(GFM_BIPOLAR_WINDS) && (GFM_BIPOLAR_WINDS == 3)
        present = 1;
#endif
        break;

      case IO_FOF_GASMETAL:
      case IO_FOF_STARMETAL:
      case IO_SUB_GASMETAL:
      case IO_SUB_GASMETALHALFRAD:
      case IO_SUB_GASMETALMAXRAD:
      case IO_SUB_GASMETALSFR:
      case IO_SUB_GASMETALSFRWEIGHTED:
      case IO_SUB_STARMETAL:
      case IO_SUB_STARMETALHALFRAD:
      case IO_SUB_STARMETALMAXRAD:
#ifdef GFM_STELLAR_EVOLUTION
        present = 1;
#endif
        break;

      case IO_FOF_GASMETALELEMENTS:
      case IO_FOF_STARMETALELEMENTS:
      case IO_SUB_GASMETALELEMENTS:
      case IO_SUB_GASMETALELEMENTSHALFRAD:
      case IO_SUB_GASMETALELEMENTSMAXRAD:
      case IO_SUB_GASMETALELEMENTSSFR:
      case IO_SUB_GASMETALELEMENTSSFRWEIGHTED:
      case IO_SUB_STARMETALELEMENTS:
      case IO_SUB_STARMETALELEMENTSHALFRAD:
      case IO_SUB_STARMETALELEMENTSMAXRAD:
#if defined(GFM_STELLAR_EVOLUTION) && !defined(GFM_STELLAR_EVOLUTION_NO_ELEMENTS)
        present = 1;
#endif
        break;

      case IO_FOF_GASDUSTMETAL:
      case IO_SUB_GASDUSTMETAL:
      case IO_SUB_GASDUSTMETALHALFRAD:
      case IO_SUB_GASDUSTMETALMAXRAD:
      case IO_SUB_GASDUSTMETALSFR:
      case IO_SUB_GASDUSTMETALSFRWEIGHTED:
#ifdef GFM_DUST
        present = 1;
#endif
        break;

      case IO_FOF_RADIOLUM:
      case IO_FOF_XRAYLUM:
#ifdef BH_NF_RADIO
        present = 1;
#endif
        break;

      case IO_FOF_BHMASS:
      case IO_FOF_BHMDOT:
      case IO_SUB_BHMASS:
      case IO_SUB_BHMDOT:
#ifdef BLACK_HOLES
        present = 1;
#endif
        break;

      case IO_FOF_WINDMASS:
      case IO_SUB_WINDMASS:
#ifdef GFM_WINDS
        present = 1;
#endif
        break;
      case IO_SUB_H2MASS:
#ifdef SUBFIND_MEASURE_H2MASS
        present = 1;
#endif
        break;

      case IO_SUB_STELLARPHOTOMETRICS:
#ifdef GFM_STELLAR_PHOTOMETRICS
        present = 1;
#endif
        break;
      case IO_SUB_STELLARPHOTOMETRICSRAD:
#ifdef GFM_STELLAR_PHOTOMETRICS
        present = 1;
#endif
        break;
      case IO_SUB_STELLARPHOTOMETRICSMASSINRAD:
#ifdef GFM_STELLAR_PHOTOMETRICS
        present = 1;
#endif
        break;

      case IO_FOF_M_MEAN200:
      case IO_FOF_R_MEAN200:
      case IO_FOF_M_CRIT200:
      case IO_FOF_R_CRIT200:
      case IO_FOF_M_TOPHAT200:
      case IO_FOF_R_TOPHAT200:
      case IO_FOF_M_CRIT500:
      case IO_FOF_R_CRIT500:
      case IO_FOF_NSUBS:
      case IO_FOF_FIRSTSUB:
      case IO_SUB_LEN:
      case IO_SUB_LENTYPE:
      case IO_SUB_MTOT:
      case IO_SUB_POS:
      case IO_SUB_VEL:
      case IO_SUB_MASSTYPE:
      case IO_SUB_CM:
      case IO_SUB_SPIN:
      case IO_SUB_VELDISP:
      case IO_SUB_VMAX:
      case IO_SUB_VMAXRAD:
      case IO_SUB_HALFMASSRAD:
      case IO_SUB_HALFMASSRADTYPE:
      case IO_SUB_MASSINRAD:
      case IO_SUB_MASSINHALFRAD:
      case IO_SUB_MASSINMAXRAD:
      case IO_SUB_MASSINRADTYPE:
      case IO_SUB_MASSINHALFRADTYPE:
      case IO_SUB_MASSINMAXRADTYPE:
      case IO_SUB_IDMOSTBOUND:
      case IO_SUB_GRNR:
      case IO_SUB_PARENT:
#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_FOF_J_MEAN200:
      case IO_FOF_JDM_MEAN200:
      case IO_FOF_JGAS_MEAN200:
      case IO_FOF_JSTARS_MEAN200:
      case IO_FOF_CMFRAC_MEAN200:
      case IO_FOF_CMFRACTYPE_MEAN200:
      case IO_FOF_J_CRIT200:
      case IO_FOF_JDM_CRIT200:
      case IO_FOF_JGAS_CRIT200:
      case IO_FOF_JSTARS_CRIT200:
      case IO_FOF_CMFRAC_CRIT200:
      case IO_FOF_CMFRACTYPE_CRIT200:
      case IO_FOF_J_TOPHAT200:
      case IO_FOF_JDM_TOPHAT200:
      case IO_FOF_JGAS_TOPHAT200:
      case IO_FOF_JSTARS_TOPHAT200:
      case IO_FOF_CMFRAC_TOPHAT200:
      case IO_FOF_CMFRACTYPE_TOPHAT200:
      case IO_FOF_J_CRIT500:
      case IO_FOF_JDM_CRIT500:
      case IO_FOF_JGAS_CRIT500:
      case IO_FOF_JSTARS_CRIT500:
      case IO_FOF_CMFRAC_CRIT500:
      case IO_FOF_CMFRACTYPE_CRIT500:
      case IO_FOF_J:
      case IO_FOF_JDM:
      case IO_FOF_JGAS:
      case IO_FOF_JSTARS:
      case IO_FOF_CMFRAC:
      case IO_FOF_CMFRACTYPE:
      case IO_FOF_EKIN:
      case IO_FOF_ETHR:
      case IO_FOF_EPOT:
      case IO_FOF_MASSTYPE_MEAN200:
      case IO_FOF_MASSTYPE_CRIT200:
      case IO_FOF_MASSTYPE_CRIT500:
      case IO_FOF_MASSTYPE_TOPHAT200:
      case IO_FOF_LENTYPE_MEAN200:
      case IO_FOF_LENTYPE_CRIT200:
      case IO_FOF_LENTYPE_CRIT500:
      case IO_FOF_LENTYPE_TOPHAT200:
      case IO_FOF_EPOT_CRIT200:
      case IO_FOF_EKIN_CRIT200:
      case IO_FOF_ETHR_CRIT200:
      case IO_FOF_EPOT_MEAN200:
      case IO_FOF_EKIN_MEAN200:
      case IO_FOF_ETHR_MEAN200:
      case IO_FOF_EPOT_TOPHAT200:
      case IO_FOF_EKIN_TOPHAT200:
      case IO_FOF_ETHR_TOPHAT200:
      case IO_FOF_EPOT_CRIT500:
      case IO_FOF_EKIN_CRIT500:
      case IO_FOF_ETHR_CRIT500:
      case IO_SUB_EKIN:
      case IO_SUB_ETHR:
      case IO_SUB_EPOT:
      case IO_SUB_J:
      case IO_SUB_JDM:
      case IO_SUB_JGAS:
      case IO_SUB_JSTARS:
      case IO_SUB_JINHALFRAD:
      case IO_SUB_JDMINHALFRAD:
      case IO_SUB_JGASINHALFRAD:
      case IO_SUB_JSTARSINHALFRAD:
      case IO_SUB_JINRAD:
      case IO_SUB_JDMINRAD:
      case IO_SUB_JGASINRAD:
      case IO_SUB_JSTARSINRAD:
      case IO_SUB_CMFRAC:
      case IO_SUB_CMFRACTYPE:
      case IO_SUB_CMFRACINHALFRAD:
      case IO_SUB_CMFRACTYPEINHALFRAD:
      case IO_SUB_CMFRACINRAD:
      case IO_SUB_CMFRACTYPEINRAD:
#endif
#ifdef SUBFIND
        present = 1;
#else
        present = 0;
#endif
        break;

      case IO_FOFSUB_IDS:
#ifdef FOF_STOREIDS
        present = 1;
#else
        present = 0;
#endif
        break;

      case IO_FOF_LASTENTRY:
        terminate("should not be reached");
        break;
    }
  return present;
}

/*! \brief Get the 4 letter IO label for a given output field.
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *  \param[out] label String with the label.
 *
 *  \return void
 */
void fof_subfind_get_Tab_IO_Label(const enum fof_subfind_iofields blocknr, char *const label)
{
  switch(blocknr)
    {
      case IO_FOF_LEN:
        memcpy(label, "FLEN", IO_LABEL_SIZE);
        break;
      case IO_FOF_MTOT:
        memcpy(label, "FMAS", IO_LABEL_SIZE);
        break;
      case IO_FOF_POS:
        memcpy(label, "FPOS", IO_LABEL_SIZE);
        break;
      case IO_FOF_CM:
        memcpy(label, "FGCM", IO_LABEL_SIZE);
        break;
      case IO_FOF_POSMINPOT:
        memcpy(label, "FPMP", IO_LABEL_SIZE);
        break;
      case IO_FOF_VEL:
        memcpy(label, "FVEL", IO_LABEL_SIZE);
        break;
      case IO_FOF_LENTYPE:
        memcpy(label, "FLTY", IO_LABEL_SIZE);
        break;
      case IO_FOF_MASSTYPE:
        memcpy(label, "FMTY", IO_LABEL_SIZE);
        break;
      case IO_FOF_SFR:
        memcpy(label, "FSFR", IO_LABEL_SIZE);
        break;
      case IO_FOF_GASMETAL:
        memcpy(label, "FGZZ", IO_LABEL_SIZE);
        break;
      case IO_FOF_STARMETAL:
        memcpy(label, "FSZZ", IO_LABEL_SIZE);
        break;
      case IO_FOF_GASMETALELEMENTS:
        memcpy(label, "FGZE", IO_LABEL_SIZE);
        break;
      case IO_FOF_STARMETALELEMENTS:
        memcpy(label, "FSZE", IO_LABEL_SIZE);
        break;
      case IO_FOF_GASDUSTMETAL:
        memcpy(label, "FGZD", IO_LABEL_SIZE);
        break;
      case IO_FOF_BHMASS:
        memcpy(label, "FBHM", IO_LABEL_SIZE);
        break;
      case IO_FOF_BHMDOT:
        memcpy(label, "FMDT", IO_LABEL_SIZE);
        break;
      case IO_FOF_WINDMASS:
        memcpy(label, "GWIM", IO_LABEL_SIZE);
        break;
      case IO_FOF_DENSLVEC:
        memcpy(label, "GDGS", IO_LABEL_SIZE);
        break;

      case IO_FOF_M_MEAN200:
        memcpy(label, "FMM2", IO_LABEL_SIZE);
        break;
      case IO_FOF_R_MEAN200:
        memcpy(label, "FRM2", IO_LABEL_SIZE);
        break;
      case IO_FOF_M_CRIT200:
        memcpy(label, "FMC2", IO_LABEL_SIZE);
        break;
      case IO_FOF_R_CRIT200:
        memcpy(label, "FRC2", IO_LABEL_SIZE);
        break;
      case IO_FOF_M_TOPHAT200:
        memcpy(label, "FMT2", IO_LABEL_SIZE);
        break;
      case IO_FOF_R_TOPHAT200:
        memcpy(label, "FRT2", IO_LABEL_SIZE);
        break;
      case IO_FOF_M_CRIT500:
        memcpy(label, "FMC5", IO_LABEL_SIZE);
        break;
      case IO_FOF_R_CRIT500:
        memcpy(label, "FRC5", IO_LABEL_SIZE);
        break;
      case IO_FOF_NSUBS:
        memcpy(label, "FNSH", IO_LABEL_SIZE);
        break;
      case IO_FOF_FIRSTSUB:
        memcpy(label, "FFSH", IO_LABEL_SIZE);
        break;
      case IO_FOF_FUZZOFFTYPE:
        memcpy(label, "FUOF", IO_LABEL_SIZE);
        break;
      case IO_FOF_RADIOLUM:
        memcpy(label, "FRAD", IO_LABEL_SIZE);
        break;
      case IO_FOF_XRAYLUM:
        memcpy(label, "FXRY", IO_LABEL_SIZE);
        break;

      case IO_SUB_LEN:
        memcpy(label, "SLEN", IO_LABEL_SIZE);
        break;
      case IO_SUB_MTOT:
        memcpy(label, "SMAS", IO_LABEL_SIZE);
        break;
      case IO_SUB_POS:
        memcpy(label, "SPOS", IO_LABEL_SIZE);
        break;
      case IO_SUB_VEL:
        memcpy(label, "SVEL", IO_LABEL_SIZE);
        break;
      case IO_SUB_LENTYPE:
        memcpy(label, "SLTY", IO_LABEL_SIZE);
        break;
      case IO_SUB_MASSTYPE:
        memcpy(label, "SMTY", IO_LABEL_SIZE);
        break;
      case IO_SUB_CM:
        memcpy(label, "SCMP", IO_LABEL_SIZE);
        break;
      case IO_SUB_SPIN:
        memcpy(label, "SSPI", IO_LABEL_SIZE);
        break;
      case IO_SUB_VELDISP:
        memcpy(label, "SVDI", IO_LABEL_SIZE);
        break;
      case IO_SUB_VMAX:
        memcpy(label, "SVMX", IO_LABEL_SIZE);
        break;
      case IO_SUB_VMAXRAD:
        memcpy(label, "SVRX", IO_LABEL_SIZE);
        break;
      case IO_SUB_HALFMASSRAD:
        memcpy(label, "SHMR", IO_LABEL_SIZE);
        break;
      case IO_SUB_HALFMASSRADTYPE:
        memcpy(label, "SHMT", IO_LABEL_SIZE);
        break;
      case IO_SUB_MASSINRAD:
        memcpy(label, "SMIR", IO_LABEL_SIZE);
        break;
      case IO_SUB_MASSINHALFRAD:
        memcpy(label, "SMIH", IO_LABEL_SIZE);
        break;
      case IO_SUB_MASSINMAXRAD:
        memcpy(label, "SMIM", IO_LABEL_SIZE);
        break;
      case IO_SUB_MASSINRADTYPE:
        memcpy(label, "SMIT", IO_LABEL_SIZE);
        break;
      case IO_SUB_MASSINHALFRADTYPE:
        memcpy(label, "SMHT", IO_LABEL_SIZE);
        break;
      case IO_SUB_MASSINMAXRADTYPE:
        memcpy(label, "SMMT", IO_LABEL_SIZE);
        break;
      case IO_SUB_IDMOSTBOUND:
        memcpy(label, "SIDM", IO_LABEL_SIZE);
        break;
      case IO_SUB_GRNR:
        memcpy(label, "SGNR", IO_LABEL_SIZE);
        break;
      case IO_SUB_PARENT:
        memcpy(label, "SPRT", IO_LABEL_SIZE);
        break;
      case IO_SUB_BFLD_HALO:
        memcpy(label, "BFDH", IO_LABEL_SIZE);
        break;
      case IO_SUB_BFLD_DISK:
        memcpy(label, "BFDD", IO_LABEL_SIZE);
        break;
      case IO_SUB_SFR:
        memcpy(label, "SSFR", IO_LABEL_SIZE);
        break;
      case IO_SUB_SFRINRAD:
        memcpy(label, "SSFI", IO_LABEL_SIZE);
        break;
      case IO_SUB_SFRINHALFRAD:
        memcpy(label, "SSFH", IO_LABEL_SIZE);
        break;
      case IO_SUB_SFRINMAXRAD:
        memcpy(label, "SSFM", IO_LABEL_SIZE);
        break;
      case IO_SUB_GASMETAL:
        memcpy(label, "SGZZ", IO_LABEL_SIZE);
        break;
      case IO_SUB_GASMETALHALFRAD:
        memcpy(label, "SGZH", IO_LABEL_SIZE);
        break;
      case IO_SUB_GASMETALMAXRAD:
        memcpy(label, "SGZM", IO_LABEL_SIZE);
        break;
      case IO_SUB_GASMETALSFR:
        memcpy(label, "SGSZ", IO_LABEL_SIZE);
        break;
      case IO_SUB_GASMETALSFRWEIGHTED:
        memcpy(label, "SGZW", IO_LABEL_SIZE);
        break;
      case IO_SUB_STARMETAL:
        memcpy(label, "SSZZ", IO_LABEL_SIZE);
        break;
      case IO_SUB_STARMETALHALFRAD:
        memcpy(label, "SSZH", IO_LABEL_SIZE);
        break;
      case IO_SUB_STARMETALMAXRAD:
        memcpy(label, "SSZM", IO_LABEL_SIZE);
        break;
      case IO_SUB_GASMETALELEMENTS:
        memcpy(label, "SGZE", IO_LABEL_SIZE);
        break;
      case IO_SUB_GASMETALELEMENTSHALFRAD:
        memcpy(label, "SGEH", IO_LABEL_SIZE);
        break;
      case IO_SUB_GASMETALELEMENTSMAXRAD:
        memcpy(label, "SGEM", IO_LABEL_SIZE);
        break;
      case IO_SUB_GASMETALELEMENTSSFR:
        memcpy(label, "SGSE", IO_LABEL_SIZE);
        break;
      case IO_SUB_GASMETALELEMENTSSFRWEIGHTED:
        memcpy(label, "SGWE", IO_LABEL_SIZE);
        break;
      case IO_SUB_STARMETALELEMENTS:
        memcpy(label, "SSZE", IO_LABEL_SIZE);
        break;
      case IO_SUB_STARMETALELEMENTSHALFRAD:
        memcpy(label, "SSEH", IO_LABEL_SIZE);
        break;
      case IO_SUB_STARMETALELEMENTSMAXRAD:
        memcpy(label, "SSEM", IO_LABEL_SIZE);
        break;
      case IO_SUB_GASDUSTMETAL:
        memcpy(label, "SGDZ", IO_LABEL_SIZE);
        break;
      case IO_SUB_GASDUSTMETALHALFRAD:
        memcpy(label, "SGDH", IO_LABEL_SIZE);
        break;
      case IO_SUB_GASDUSTMETALMAXRAD:
        memcpy(label, "SGDM", IO_LABEL_SIZE);
        break;
      case IO_SUB_GASDUSTMETALSFR:
        memcpy(label, "SGDS", IO_LABEL_SIZE);
        break;
      case IO_SUB_GASDUSTMETALSFRWEIGHTED:
        memcpy(label, "SGDW", IO_LABEL_SIZE);
        break;
      case IO_SUB_BHMASS:
        memcpy(label, "SBHM", IO_LABEL_SIZE);
        break;
      case IO_SUB_BHMDOT:
        memcpy(label, "SMDT", IO_LABEL_SIZE);
        break;
      case IO_SUB_WINDMASS:
        memcpy(label, "SWIM", IO_LABEL_SIZE);
        break;
      case IO_SUB_H2MASS:
        memcpy(label, "SH2M", IO_LABEL_SIZE);
        break;
      case IO_SUB_STELLARPHOTOMETRICS:
        memcpy(label, "SSPH", IO_LABEL_SIZE);
        break;
      case IO_SUB_STELLARPHOTOMETRICSRAD:
        memcpy(label, "SSPR", IO_LABEL_SIZE);
        break;
      case IO_SUB_STELLARPHOTOMETRICSMASSINRAD:
        memcpy(label, "SSPM", IO_LABEL_SIZE);
        break;
      case IO_FOFSUB_IDS:
        memcpy(label, "PIDS", IO_LABEL_SIZE);
        break;

#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_FOF_J_MEAN200:
        memcpy(label, "FJM2", IO_LABEL_SIZE);
        break;
      case IO_FOF_JDM_MEAN200:
        memcpy(label, "JDM2", IO_LABEL_SIZE);
        break;
      case IO_FOF_JGAS_MEAN200:
        memcpy(label, "JGM2", IO_LABEL_SIZE);
        break;
      case IO_FOF_JSTARS_MEAN200:
        memcpy(label, "JSM2", IO_LABEL_SIZE);
        break;
      case IO_FOF_MASSTYPE_MEAN200:
        memcpy(label, "MTM2", IO_LABEL_SIZE);
        break;
      case IO_FOF_LENTYPE_MEAN200:
        memcpy(label, "LTM2", IO_LABEL_SIZE);
        break;
      case IO_FOF_CMFRAC_MEAN200:
        memcpy(label, "CFM2", IO_LABEL_SIZE);
        break;
      case IO_FOF_CMFRACTYPE_MEAN200:
        memcpy(label, "FTM2", IO_LABEL_SIZE);
        break;
      case IO_FOF_J_CRIT200:
        memcpy(label, "FJC2", IO_LABEL_SIZE);
        break;
      case IO_FOF_JDM_CRIT200:
        memcpy(label, "JDC2", IO_LABEL_SIZE);
        break;
      case IO_FOF_JGAS_CRIT200:
        memcpy(label, "JGC2", IO_LABEL_SIZE);
        break;
      case IO_FOF_JSTARS_CRIT200:
        memcpy(label, "JSC2", IO_LABEL_SIZE);
        break;
      case IO_FOF_MASSTYPE_CRIT200:
        memcpy(label, "MTC2", IO_LABEL_SIZE);
        break;
      case IO_FOF_LENTYPE_CRIT200:
        memcpy(label, "LTC2", IO_LABEL_SIZE);
        break;
      case IO_FOF_CMFRAC_CRIT200:
        memcpy(label, "CFC2", IO_LABEL_SIZE);
        break;
      case IO_FOF_CMFRACTYPE_CRIT200:
        memcpy(label, "FTC2", IO_LABEL_SIZE);
        break;
      case IO_FOF_J_TOPHAT200:
        memcpy(label, "FJT2", IO_LABEL_SIZE);
        break;
      case IO_FOF_JDM_TOPHAT200:
        memcpy(label, "JDT2", IO_LABEL_SIZE);
        break;
      case IO_FOF_JGAS_TOPHAT200:
        memcpy(label, "JGT2", IO_LABEL_SIZE);
        break;
      case IO_FOF_JSTARS_TOPHAT200:
        memcpy(label, "JST2", IO_LABEL_SIZE);
        break;
      case IO_FOF_MASSTYPE_TOPHAT200:
        memcpy(label, "MTT2", IO_LABEL_SIZE);
        break;
      case IO_FOF_LENTYPE_TOPHAT200:
        memcpy(label, "LTT2", IO_LABEL_SIZE);
        break;
      case IO_FOF_CMFRAC_TOPHAT200:
        memcpy(label, "CFT2", IO_LABEL_SIZE);
        break;
      case IO_FOF_CMFRACTYPE_TOPHAT200:
        memcpy(label, "FTT2", IO_LABEL_SIZE);
        break;
      case IO_FOF_J_CRIT500:
        memcpy(label, "FJC5", IO_LABEL_SIZE);
        break;
      case IO_FOF_JDM_CRIT500:
        memcpy(label, "JDC5", IO_LABEL_SIZE);
        break;
      case IO_FOF_JGAS_CRIT500:
        memcpy(label, "JGC5", IO_LABEL_SIZE);
        break;
      case IO_FOF_JSTARS_CRIT500:
        memcpy(label, "JSC5", IO_LABEL_SIZE);
        break;
      case IO_FOF_MASSTYPE_CRIT500:
        memcpy(label, "MTC5", IO_LABEL_SIZE);
        break;
      case IO_FOF_LENTYPE_CRIT500:
        memcpy(label, "LTC5", IO_LABEL_SIZE);
        break;
      case IO_FOF_CMFRAC_CRIT500:
        memcpy(label, "CFC5", IO_LABEL_SIZE);
        break;
      case IO_FOF_CMFRACTYPE_CRIT500:
        memcpy(label, "FTC5", IO_LABEL_SIZE);
        break;
      case IO_FOF_J:
        memcpy(label, "FOFJ", IO_LABEL_SIZE);
        break;
      case IO_FOF_JDM:
        memcpy(label, "FOJD", IO_LABEL_SIZE);
        break;
      case IO_FOF_JGAS:
        memcpy(label, "FOJG", IO_LABEL_SIZE);
        break;
      case IO_FOF_JSTARS:
        memcpy(label, "FOJS", IO_LABEL_SIZE);
        break;
      case IO_FOF_CMFRAC:
        memcpy(label, "FOCF", IO_LABEL_SIZE);
        break;
      case IO_FOF_CMFRACTYPE:
        memcpy(label, "FOFT", IO_LABEL_SIZE);
        break;
      case IO_FOF_EKIN:
        memcpy(label, "EKIN", IO_LABEL_SIZE);
        break;
      case IO_FOF_ETHR:
        memcpy(label, "ETHR", IO_LABEL_SIZE);
        break;
      case IO_FOF_EPOT:
        memcpy(label, "EPOT", IO_LABEL_SIZE);
        break;

      case IO_FOF_EPOT_CRIT200:
        memcpy(label, "EPO1", IO_LABEL_SIZE);
        break;
      case IO_FOF_EKIN_CRIT200:
        memcpy(label, "EKI1", IO_LABEL_SIZE);
        break;
      case IO_FOF_ETHR_CRIT200:
        memcpy(label, "ETH1", IO_LABEL_SIZE);
        break;
      case IO_FOF_EPOT_MEAN200:
        memcpy(label, "EPO2", IO_LABEL_SIZE);
        break;
      case IO_FOF_EKIN_MEAN200:
        memcpy(label, "EKI2", IO_LABEL_SIZE);
        break;
      case IO_FOF_ETHR_MEAN200:
        memcpy(label, "ETH2", IO_LABEL_SIZE);
        break;
      case IO_FOF_EPOT_TOPHAT200:
        memcpy(label, "EPO3", IO_LABEL_SIZE);
        break;
      case IO_FOF_EKIN_TOPHAT200:
        memcpy(label, "EKI3", IO_LABEL_SIZE);
        break;
      case IO_FOF_ETHR_TOPHAT200:
        memcpy(label, "ETH3", IO_LABEL_SIZE);
        break;
      case IO_FOF_EPOT_CRIT500:
        memcpy(label, "EPO4", IO_LABEL_SIZE);
        break;
      case IO_FOF_EKIN_CRIT500:
        memcpy(label, "EKI4", IO_LABEL_SIZE);
        break;
      case IO_FOF_ETHR_CRIT500:
        memcpy(label, "ETH4", IO_LABEL_SIZE);
        break;

      case IO_SUB_EKIN:
        memcpy(label, "SEKN", IO_LABEL_SIZE);
        break;
      case IO_SUB_ETHR:
        memcpy(label, "SETH", IO_LABEL_SIZE);
        break;
      case IO_SUB_EPOT:
        memcpy(label, "SEPT", IO_LABEL_SIZE);
        break;
      case IO_SUB_J:
        memcpy(label, "SUBJ", IO_LABEL_SIZE);
        break;
      case IO_SUB_JDM:
        memcpy(label, "SJDM", IO_LABEL_SIZE);
        break;
      case IO_SUB_JGAS:
        memcpy(label, "SJGS", IO_LABEL_SIZE);
        break;
      case IO_SUB_JSTARS:
        memcpy(label, "SJST", IO_LABEL_SIZE);
        break;
      case IO_SUB_JINHALFRAD:
        memcpy(label, "SJHR", IO_LABEL_SIZE);
        break;
      case IO_SUB_JDMINHALFRAD:
        memcpy(label, "SJDH", IO_LABEL_SIZE);
        break;
      case IO_SUB_JGASINHALFRAD:
        memcpy(label, "SJGH", IO_LABEL_SIZE);
        break;
      case IO_SUB_JSTARSINHALFRAD:
        memcpy(label, "SJSH", IO_LABEL_SIZE);
        break;
      case IO_SUB_JINRAD:
        memcpy(label, "SJMR", IO_LABEL_SIZE);
        break;
      case IO_SUB_JDMINRAD:
        memcpy(label, "SJDR", IO_LABEL_SIZE);
        break;
      case IO_SUB_JGASINRAD:
        memcpy(label, "SJGR", IO_LABEL_SIZE);
        break;
      case IO_SUB_JSTARSINRAD:
        memcpy(label, "SJSR", IO_LABEL_SIZE);
        break;
      case IO_SUB_CMFRAC:
        memcpy(label, "SCMF", IO_LABEL_SIZE);
        break;
      case IO_SUB_CMFRACTYPE:
        memcpy(label, "SCMT", IO_LABEL_SIZE);
        break;
      case IO_SUB_CMFRACINHALFRAD:
        memcpy(label, "SCMH", IO_LABEL_SIZE);
        break;
      case IO_SUB_CMFRACTYPEINHALFRAD:
        memcpy(label, "SCTH", IO_LABEL_SIZE);
        break;
      case IO_SUB_CMFRACINRAD:
        memcpy(label, "SCMR", IO_LABEL_SIZE);
        break;
      case IO_SUB_CMFRACTYPEINRAD:
        memcpy(label, "SCTR", IO_LABEL_SIZE);
        break;
#endif

      case IO_FOF_LASTENTRY:
        terminate("should not be reached");
        break;
    }
}

#ifdef HAVE_HDF5
void fof_subfind_write_header_attributes_in_hdf5(const hid_t handle)
{
  hid_t hdf5_dataspace, hdf5_attribute;

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Ngroups_ThisFile", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &catalogue_header.Ngroups, "Ngroups_ThisFile");
  my_H5Aclose(hdf5_attribute, "Ngroups_ThisFile");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Nsubgroups_ThisFile", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &catalogue_header.Nsubgroups, "Nsubgroups_ThisFile");
  my_H5Aclose(hdf5_attribute, "Nsubgroups_ThisFile");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Nids_ThisFile", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &catalogue_header.Nids, "Nids_ThisFile");
  my_H5Aclose(hdf5_attribute, "Nids_ThisFile");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Ngroups_Total", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &catalogue_header.TotNgroups, "Ngroups_Total");
  my_H5Aclose(hdf5_attribute, "Ngroups_Total");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Nsubgroups_Total", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &catalogue_header.TotNsubgroups, "Nsubgroups_Total");
  my_H5Aclose(hdf5_attribute, "Nsubgroups_Total");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Nids_Total", H5T_NATIVE_INT64, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT64, &catalogue_header.TotNids, "Nids_Total");
  my_H5Aclose(hdf5_attribute, "Nids_Total");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "NumFiles", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &catalogue_header.num_files, "NumFiles");
  my_H5Aclose(hdf5_attribute, "NumFiles");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Time", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &catalogue_header.time, "Time");
  my_H5Aclose(hdf5_attribute, "Time");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Redshift", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &catalogue_header.redshift, "Redshift");
  my_H5Aclose(hdf5_attribute, "Redshift");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "HubbleParam", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &catalogue_header.HubbleParam, "HubbleParam");
  my_H5Aclose(hdf5_attribute, "HubbleParam");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "BoxSize", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &catalogue_header.BoxSize, "BoxSize");
  my_H5Aclose(hdf5_attribute, "BoxSize");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Omega0", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &catalogue_header.Omega0, "Omega0");
  my_H5Aclose(hdf5_attribute, "Omega0");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "OmegaLambda", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &catalogue_header.OmegaLambda, "OmegaLambda");
  my_H5Aclose(hdf5_attribute, "OmegaLambda");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "FlagDoubleprecision", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &catalogue_header.flag_doubleprecision, "FlagDoubleprecision");
  my_H5Aclose(hdf5_attribute, "FlagDoubleprecision");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hid_t atype = my_H5Tcopy(H5T_C_S1);

  my_H5Tset_size(atype, strlen(GIT_COMMIT));
  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Git_commit", atype, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, atype, GIT_COMMIT, "Git_commit");
  my_H5Aclose(hdf5_attribute, "Git_commit");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  my_H5Tset_size(atype, strlen(GIT_DATE));
  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Git_date", atype, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, atype, GIT_DATE, "Git_date");
  my_H5Aclose(hdf5_attribute, "Git_date");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);
}
#endif

#endif
