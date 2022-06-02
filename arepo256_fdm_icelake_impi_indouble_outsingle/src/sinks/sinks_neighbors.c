/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sinks/sinks_neighbors.c
 * \date        07/2017
 * \author      Federico Marinacci
 * \brief       find neighbouring sink particles and detect potential mergers
 * \details     largely based on src/blackhole/blackhole_neighbors.c
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

#include "../allvars.h"
#include "../proto.h"

#ifdef SINKS
#ifdef SINKS_MERGERS

static int sinks_tree_ngb_evaluate(int i, int mode, int threadid);

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyDouble Mass;
  MyFloat Sinks_Hsml;
  MyIDType ID;

  int Firstnode;
} data_in;

static data_in *DataGet;

/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in *in, int i, int firstnode)
{
  int k;

  if(i < NumPart)
    {
      for(k = 0; k < 3; k++)
        in->Pos[k] = P[i].Pos[k];

      in->Mass       = P[i].Mass;
      in->Sinks_Hsml = SKP(i).Sinks_Hsml;
      in->ID         = P[i].ID;
    }
  else
    {
      i -= Tree_ImportedNodeOffset;

      for(k = 0; k < 3; k++)
        in->Pos[k] = Tree_Points[i].Pos[k];

      in->Mass       = Tree_Points[i].Mass;
      in->Sinks_Hsml = TSKP(i).Sinks_Hsml;
      in->ID         = TSKP(i).ID;
    }

  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
  MyFloat dummy;
} data_out;

static data_out *DataResult;

/* routine to store or combine result data */
static void out2particle(data_out *out, int i, int mode) { return; }

#include "../generic_comm_helpers2.h"

static void kernel_local(void)
{
  int i;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif
#pragma omp parallel private(i)
  {
    int j, threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              if(generic_polling_primary(count, Nforces))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= Nforces)
          break;

        int idx = TargetList[i];

        sinks_tree_ngb_evaluate(idx, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}

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

        sinks_tree_ngb_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

/* This routines uses the gravitational tree to search in the Sinks_Hsml neighborhood around each
 * active sinks for other sinks (which become candidates for merging).
 */
void sinks_find_neighboring_sinks(void)
{
  int idx, i, ncount;

  TIMER_START(CPU_SINKS_NGB);

  mpi_printf("SINKS: Begin finding sinks neighbours\n");

  generic_set_MaxNexport();

  /* Create list of targets. We do this here to simplify the treatment of the two possible locations of source points */
  TargetList           = mymalloc("TargetList", (NumPart + Tree_NumPartImported) * sizeof(int));
  Tree_ResultIndexList = mymalloc("Tree_ResultIndexList", Tree_NumPartImported * sizeof(int));

  Nforces = 0;
  for(idx = 0; idx < TimeBinsSinksAccretion.NActiveParticles; idx++)
    {
      i = TimeBinsSinksAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type == 5 && Tree_Task_list[i] == ThisTask)
        TargetList[Nforces++] = i;
    }

  for(i = 0, ncount = 0; i < Tree_NumPartImported; i++)
#ifndef HIERARCHICAL_GRAVITY
    if(Tree_Points[i].ActiveFlag)
#endif
      if(Tree_Points[i].Type == 5)
        {
          Tree_ResultIndexList[i] = ncount++;
          TargetList[Nforces++]   = i + Tree_ImportedNodeOffset;
        }

  generic_comm_pattern(Nforces, kernel_local, kernel_imported);

  myfree(Tree_ResultIndexList);
  myfree(TargetList);

  mpi_printf("SINKS: done with sinks neighbours\n");

  TIMER_STOP(CPU_SINKS_NGB);
}

static int sinks_tree_ngb_evaluate(int target, int mode, int threadid)
{
  int k, numnodes, *firstnode;
  int no;
  double h, h2;

  double dx, dy, dz, r2;
  MyIDType id;
  MyDouble *pos;

  data_in local, *target_data;
  data_out out;

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

  pos = target_data->Pos;
  h   = target_data->Sinks_Hsml * All.HubbleParam / All.cf_atime;
  id  = target_data->ID;

  h2 = h * h;

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
              dx = GRAVITY_NEAREST_X(Tree_Pos_list[3 * no + 0] - pos[0]);
              dy = GRAVITY_NEAREST_Y(Tree_Pos_list[3 * no + 1] - pos[1]);
              dz = GRAVITY_NEAREST_Z(Tree_Pos_list[3 * no + 2] - pos[2]);

              r2 = dx * dx + dy * dy + dz * dz;

              if(r2 < h2)
                {
                  if(P[no].Type == 5) /* we have a potential sink merger */
                    {
                      if(id != P[no].ID)
                        {
                          if(SKP(no).SwallowID < id && P[no].ID < id)
                            SKP(no).SwallowID = id;
                        }
                    }
                }

              no = Nextnode[no];
            }
          else if(no < Tree_MaxPart + Tree_MaxNodes) /* internal node */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                {
                  if(no <
                     Tree_FirstNonTopLevelNode) /* we reached a top-level node again, which means that we are done with the branch */
                    break;
                }

              struct NODE *current = &Nodes[no];

              no = current->u.d.sibling; /* in case the node can be discarded */

              double dist = h + 0.5 * current->len;
              dx          = NGB_PERIODIC_LONG_X(current->center[0] - pos[0]);
              if(dx > dist)
                continue;
              dy = NGB_PERIODIC_LONG_Y(current->center[1] - pos[1]);
              if(dy > dist)
                continue;
              dz = NGB_PERIODIC_LONG_Z(current->center[2] - pos[2]);
              if(dz > dist)
                continue;
              /* now test against the minimal sphere enclosing everything */
              dist += FACT1 * current->len;
              if(dx * dx + dy * dy + dz * dz > dist * dist)
                continue;

              no = current->u.d.nextnode; /* ok, we need to open the node */
            }
          else if(no >= Tree_ImportedNodeOffset) /* point from imported nodelist */
            {
              int n = no - Tree_ImportedNodeOffset;

              dx = GRAVITY_NEAREST_X(Tree_Points[n].Pos[0] - pos[0]);
              dy = GRAVITY_NEAREST_Y(Tree_Points[n].Pos[1] - pos[1]);
              dz = GRAVITY_NEAREST_Z(Tree_Points[n].Pos[2] - pos[2]);

              r2 = dx * dx + dy * dy + dz * dz;

              if(r2 < h2)
                {
                  if(Tree_Points[n].Type == 5) /* we have a potential sink merger */
                    {
                      if(id != TSKP(n).ID)
                        {
                          if(TSKP(n).SwallowID < id && TSKP(n).ID < id)
                            TSKP(n).SwallowID = id;
                        }
                    }
                }

              no = Nextnode[no - Tree_MaxNodes];
            }
          else /* pseudo particle */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                terminate("mode == MODE_IMPORTED_PARTICLES");

              if(target >= 0)
                tree_treefind_export_node_threads(no, target, threadid);

              no = Nextnode[no - Tree_MaxNodes];
              continue;
            }
        }
    }

  out.dummy = 0.0;

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif
#endif
