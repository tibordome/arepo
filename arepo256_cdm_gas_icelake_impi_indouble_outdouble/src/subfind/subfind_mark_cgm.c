/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/subfind/subfind_mark_cgm.c
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
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "../allvars.h"
#include "../proto.h"

#if defined(SUBFIND) && defined(REFINEMENT_CGM)

#include "../fof/fof.h"
#include "subfind.h"

static int subfind_mark_cgm_evaluate(int target, int mode, int threadid);

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  double R200;

  int Firstnode;
} data_in;

static data_in *DataGet;

/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in *in, int i, int firstnode)
{
  in->Pos[0] = Group[i].Pos[0];
  in->Pos[1] = Group[i].Pos[1];
  in->Pos[2] = Group[i].Pos[2];
#ifdef REFINEMENT_CGM_USE_R200M
  in->R200 = All.FracRadiusForCGMRefinement * Group[i].R_Mean200;
#else
  in->R200 = All.FracRadiusForCGMRefinement * Group[i].R_Crit200;
#endif

  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
  char dummy;
} data_out;

static data_out *DataResult;

/* routine to store or combine result data */
static void out2particle(data_out *out, int i, int mode) {}

#include "../generic_comm_helpers2.h"

static long long count_marked;
static double volume_marked;

/*! \brief Routine that defines what to do with local particles.
 *
 *  Calls the *_evaluate function in MODE_LOCAL_PARTICLES.
 *
 *  \return void
 */
static void kernel_local(void)
{
  int i;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif

#pragma omp parallel private(i)
  {
    int threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    for(int j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              if(generic_polling_primary(count, Ngroups))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= Ngroups)
          break;

#ifdef REFINEMENT_CGM_USE_R200M
        if(Group[i].M_Mean200 > All.MinMassForCGMRefinement)
#else
        if(Group[i].M_Crit200 > All.MinMassForCGMRefinement)
#endif
          subfind_mark_cgm_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}

/*! \brief Routine that defines what to do with imported particles.
 *
 *  Calls the *_evaluate function in MODE_IMPORTED_PARTICLES.
 *
 *  \return void
 */
static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, cnt = 0;
#pragma omp parallel private(i)
  {
    int threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    while(1)
      {
#pragma omp atomic capture
        i = cnt++;

        if(i >= Nimport)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              generic_polling_secondary();
          }

        count++;
#endif

        subfind_mark_cgm_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void subfind_mark_cgm(void)
{
  for(int p = 0; p < NumGas; p++)
    {
      SphP[p].HighResMassCGM    = 0;
      SphP[p].HighResDensityCGM = 0;
    }

  count_marked  = 0;
  volume_marked = 0;

  /* allocate buffers to arrange communication */
  generic_set_MaxNexport();

  generic_comm_pattern(Ngroups, kernel_local, kernel_imported);

  /* deal with points that may have been exported in the tree construction */
  struct treepoint_data *export_Tree_Points =
      (struct treepoint_data *)mymalloc("export_Tree_Points", Tree_NumPartExported * sizeof(struct treepoint_data));

  /* exchange  data */
  for(int ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Force_Send_count[recvTask] > 0 || Force_Recv_count[recvTask] > 0)
          MPI_Sendrecv(&Tree_Points[Force_Recv_offset[recvTask]], Force_Recv_count[recvTask] * sizeof(struct treepoint_data), MPI_BYTE,
                       recvTask, TAG_DENS_A, &export_Tree_Points[Force_Send_offset[recvTask]],
                       Force_Send_count[recvTask] * sizeof(struct treepoint_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                       MPI_STATUS_IGNORE);
    }

  for(int n = 0; n < Tree_NumPartExported; n++)
    {
      if(export_Tree_Points[n].marked_flag == 1)
        {
          int p = export_Tree_Points[n].index;

          if(P[p].Type == 0)
            {
              if(SphP[p].HighResMassCGM == 0)
                {
                  SphP[p].HighResMassCGM    = P[p].Mass;
                  SphP[p].HighResDensityCGM = SphP[p].Density;

                  count_marked++;
                  volume_marked += SphP[p].Volume;
                }
            }
        }
    }

  myfree(export_Tree_Points);

  MPI_Allreduce(MPI_IN_PLACE, &count_marked, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &volume_marked, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  mpi_printf("REFINEMENT_CGM:  marked %lld cells, total physical volume=%g, expected number of cells in this region at least %g\n",
             (long long)count_marked, volume_marked / All.cf_a3inv, volume_marked / All.TargetGasVolume);
}

/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
static int subfind_mark_cgm_evaluate(int target, int mode, int threadid)
{
  int k, p, no, numnodes, *firstnode;
  double hsml;
  MyDouble *pos;
  struct NODE *current;
  MyDouble dx, dy, dz, dist, r2;

  data_in local, *target_data;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data = &local;

      numnodes  = 1;
      firstnode = NULL;
    }
  else
    {
      target_data = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  pos  = target_data->Pos;
  hsml = target_data->R200;

  for(k = 0; k < numnodes; k++)
    {
      if(mode == MODE_LOCAL_PARTICLES)
        {
          no = Tree_MaxPart; /* root node */
        }
      else
        {
          no = firstnode[k];
          no = Nodes[no].u.d.nextnode; /* open it */
        }

      while(no >= 0)
        {
          if(no < Tree_MaxPart) /* single particle */
            {
              p  = no;
              no = Nextnode[no];

              dist = hsml;
              dx   = FOF_NEAREST_LONG_X(Tree_Pos_list[3 * p + 0] - pos[0]);
              if(dx > dist)
                continue;
              dy = FOF_NEAREST_LONG_Y(Tree_Pos_list[3 * p + 1] - pos[1]);
              if(dy > dist)
                continue;
              dz = FOF_NEAREST_LONG_Z(Tree_Pos_list[3 * p + 2] - pos[2]);
              if(dz > dist)
                continue;
              if(dx * dx + dy * dy + dz * dz > dist * dist)
                continue;

              if(P[p].Type == 0)
                {
                  if(SphP[p].HighResMassCGM == 0)
                    {
                      SphP[p].HighResMassCGM    = P[p].Mass;
                      SphP[p].HighResDensityCGM = SphP[p].Density;

                      count_marked++;
                      volume_marked += SphP[p].Volume;
                    }
                }
            }
          else if(no < Tree_MaxPart + Tree_MaxNodes) /* internal node */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                {
                  if(no < Tree_FirstNonTopLevelNode)
                    /* we reached a top-level node again, which means that we are done with the branch */
                    break;
                }

              current = &Nodes[no];

              no = current->u.d.sibling; /* in case the node can be discarded */

              dist = hsml + 0.5 * current->len;
              dx   = FOF_NEAREST_LONG_X(current->center[0] - pos[0]);
              if(dx > dist)
                continue;
              dy = FOF_NEAREST_LONG_Y(current->center[1] - pos[1]);
              if(dy > dist)
                continue;
              dz = FOF_NEAREST_LONG_Z(current->center[2] - pos[2]);
              if(dz > dist)
                continue;
              /* now test against the minimal sphere enclosing everything */
              dist += FACT1 * current->len;
              if((r2 = (dx * dx + dy * dy + dz * dz)) > dist * dist)
                continue;

              no = current->u.d.nextnode; /* ok, we need to open the node */
            }
          else if(no >= Tree_ImportedNodeOffset) /* point from imported nodelist */
            {
              int n = no - Tree_ImportedNodeOffset;
              no    = Nextnode[no - Tree_MaxNodes];

              dist = hsml;
              dx   = FOF_NEAREST_LONG_X(Tree_Points[n].Pos[0] - pos[0]);
              if(dx > dist)
                continue;
              dy = FOF_NEAREST_LONG_Y(Tree_Points[n].Pos[1] - pos[1]);
              if(dy > dist)
                continue;
              dz = FOF_NEAREST_LONG_Z(Tree_Points[n].Pos[2] - pos[2]);
              if(dz > dist)
                continue;
              if(dx * dx + dy * dy + dz * dz > dist * dist)
                continue;

              if(Tree_Points[n].marked_flag == 0)
                Tree_Points[n].marked_flag = 1;
            }
          else /* pseudo particle */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                terminate("mode == MODE_IMPORTED_PARTICLES");

              if(mode == MODE_LOCAL_PARTICLES)
                tree_treefind_export_node_threads(no, target, threadid);

              no = Nextnode[no - Tree_MaxNodes];
            }
        }
    }

  return 0;
}

#endif
