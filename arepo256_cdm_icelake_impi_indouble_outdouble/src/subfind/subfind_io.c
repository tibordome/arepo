/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/subfind/subfind_io.c
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

#include <gsl/gsl_rng.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "../allvars.h"
#include "../domain.h"
#include "../fof/fof.h"
#include "../proto.h"

#ifdef SUBFIND

#include "subfind.h"

/*! \brief Saves subfind group catalogue to disk.
 *
 *  Note that this routine calls the FoF I/O routines.
 *
 *  \param[in] num Index of this snapshot output.
 *
 *  \return void
 */
void subfind_save_final(int num)
{
  int i, filenr, gr, ngrps, masterTask, lastTask, totsubs;
  char buf[MAXLEN_PATH];
  double t0, t1;

  /* prepare list of ids with assigned group numbers */
#ifdef FOF_STOREIDS
  fof_subfind_prepare_ID_list();
#endif

  t0 = second();

  /* fill in the FirstSub-values */
  for(i = 0, totsubs = 0; i < Ngroups; i++)
    {
      if(i > 0)
        Group[i].FirstSub = Group[i - 1].FirstSub + Group[i - 1].Nsubs;
      else
        Group[i].FirstSub = 0;
      totsubs += Group[i].Nsubs;
    }

  MPI_Allgather(&totsubs, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 1, Send_offset[0] = 0; i < NTask; i++)
    Send_offset[i] = Send_offset[i - 1] + Send_count[i - 1];

  for(i = 0; i < Ngroups; i++)
    {
      if(Group[i].Nsubs > 0)
        Group[i].FirstSub += Send_offset[ThisTask];
      else
        Group[i].FirstSub = -1;
    }

  CommBuffer = (void *)mymalloc("CommBuffer", COMMBUFFERSIZE);

  if(All.SnapFormat < SNAP_FORMAT_GADGET || All.SnapFormat > SNAP_FORMAT_HDF5)
    mpi_printf("Unsupported File-Format All.SnapFormat=%d \n", All.SnapFormat);

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
    file_path_sprintf(buf, "%s/groups_%03d/fof_subhalo_tab_%03d.%d", All.OutputDir, num, num, filenr);
  else
    file_path_sprintf(buf, "%s/fof_subhalo_tab_%03d", All.OutputDir, num);

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

  mpi_printf("SUBFIND: Subgroup catalogues saved. took = %g sec\n", timediff(t0, t1));
}

#endif
