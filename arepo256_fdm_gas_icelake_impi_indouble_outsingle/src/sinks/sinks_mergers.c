/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sinks/sinks_mergers.c
 * \date        05/2017
 * \author      Federico Marinacci
 * \brief       merger module for sink particles
 * \details     largely based on the BH merger modules in src/blackhole/blackhole_merger.c
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

static int sinks_tree_merger_evaluate(int i, int mode, int threadid);

struct sinks_properties
{
  MyDouble Mass;
  MyFloat Momentum[3];

  int Sinks_CountProgs;
};

static struct sinksresultsimported_data
{
  struct sinks_properties Prop;
  int index;
} * Tree_SinksngbResultsImported;

static struct sinksdeleted_data
{
  int index;
} * Tree_SinksngbDeleted;

static int N_Sinks_swallowed_local, N_Sinks_swallowed_imported;
static char *treeSinksexportflag;

static struct sinks_properties *sinks_accreted;

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
  struct sinks_properties Prop;
} data_out;

static data_out *DataResult;

/* routine to store or combine result data */
static void out2particle(data_out *out, int i, int mode)
{
  int k;

  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
      if(i < NumPart)
        {
          if(P[i].AuxDataID >= NumSinks)
            terminate("SINKS: P[target(=%d)].AuxDataID(=%lld) >= NumSinks=%d", i, (long long)P[i].AuxDataID, NumSinks);

          sinks_accreted[P[i].AuxDataID] = out->Prop;
        }
      else
        {
          int idx = Tree_ResultIndexList[i - Tree_ImportedNodeOffset];

          Tree_SinksngbResultsImported[idx].Prop = out->Prop;
        }
    }
  else /* combine */
    {
      if(i < NumPart)
        {
          if(P[i].AuxDataID >= NumSinks)
            terminate("SINKS: P[i(=%d)].AuxDataID(=%lld) >= NumSinks=%d", i, (long long)P[i].AuxDataID, NumSinks);

          sinks_accreted[P[i].AuxDataID].Mass += out->Prop.Mass;

          for(k = 0; k < 3; k++)
            sinks_accreted[P[i].AuxDataID].Momentum[k] += out->Prop.Momentum[k];

          sinks_accreted[P[i].AuxDataID].Sinks_CountProgs += out->Prop.Sinks_CountProgs;
        }
      else
        {
          int idx = Tree_ResultIndexList[i - Tree_ImportedNodeOffset];

          Tree_SinksngbResultsImported[idx].Prop.Mass += out->Prop.Mass;

          for(k = 0; k < 3; k++)
            Tree_SinksngbResultsImported[idx].Prop.Momentum[k] += out->Prop.Momentum[k];

          Tree_SinksngbResultsImported[idx].Prop.Sinks_CountProgs += out->Prop.Sinks_CountProgs;
        }
    }
}

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

        sinks_tree_merger_evaluate(idx, MODE_LOCAL_PARTICLES, threadid);
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

        sinks_tree_merger_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

/* This routine uses the gravitational tree to search in the Sink_Hsml neighborhood of
 * an active Sink for other Sinks that are to be merged, and carries out the mergers.
 * Only Sinks with SwallowID==0 can swallow other Sinks because they are guaranteed to survive themselves.
 */
void sinks_do_mergers(void)
{
  int idx, i, j, k, n, ncount, nexport, nimport;
  int Ntot_Sinks_swallowed_local, Ntot_Sinks_swallowed_imported;
  int ngrp, recvTask;

  TIMER_START(CPU_SINKS_MERGERS);

  mpi_printf("SINKS: Begin sinks mergers.\n");

  treeSinksexportflag = mymalloc_clear("treeSinksexportflag", Tree_NumSinksImported * sizeof(char));

  N_Sinks_swallowed_local    = 0;
  N_Sinks_swallowed_imported = 0;

  /* allocate temporary variables */
  sinks_accreted = mymalloc_clear("sinks_accreted", NumSinks * sizeof(struct sinks_properties));

  generic_set_MaxNexport();

  /* Create list of is. We do this here to simplify the treatment of the two possible locations of source points */
  TargetList           = mymalloc("TargetList", (NumPart + Tree_NumPartImported) * sizeof(int));
  Tree_ResultIndexList = mymalloc("Tree_ResultIndexList", Tree_NumPartImported * sizeof(int));

  Nforces = 0;

  for(idx = 0; idx < TimeBinsSinksAccretion.NActiveParticles; idx++)
    {
      i = TimeBinsSinksAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(Tree_Task_list[i] == ThisTask)
        if(SKP(i).SwallowID == 0)
          TargetList[Nforces++] = i;
    }

  for(i = 0, ncount = 0; i < Tree_NumPartImported; i++)
#ifndef HIERARCHICAL_GRAVITY
    if(Tree_Points[i].ActiveFlag)
#endif
      if(Tree_Points[i].Type == 5)
        if(TSKP(i).SwallowID == 0)
          {
            Tree_ResultIndexList[i] = ncount++;
            TargetList[Nforces++]   = i + Tree_ImportedNodeOffset;
          }

  Tree_SinksngbResultsImported = mymalloc_clear("Tree_SinksngbResultsImported", ncount * sizeof(struct sinksresultsimported_data));

  generic_comm_pattern(Nforces, kernel_local, kernel_imported);

  /* now communicate the results in Tree_BhngbResultsImported */

  for(j = 0; j < NTask; j++)
    Recv_count[j] = 0;

  for(i = 0, n = 0, k = 0; i < NTask; i++)
    for(j = 0; j < Force_Recv_count[i]; j++, n++)
      {
#ifndef HIERARCHICAL_GRAVITY
        if(Tree_Points[n].ActiveFlag)
#endif
          if(Tree_Points[n].Type == 5)
            if(TSKP(n).SwallowID == 0)
              {
                Tree_SinksngbResultsImported[k].index = Tree_Points[n].index;
                Recv_count[i]++;
                k++;
              }
      }

  MPI_Alltoall(Recv_count, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nexport = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nexport += Send_count[j];
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  struct sinksresultsimported_data *tmp_results = mymalloc("tmp_results", nexport * sizeof(struct sinksresultsimported_data));
  memset(tmp_results, -1, nexport * sizeof(struct sinksresultsimported_data));

  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&Tree_SinksngbResultsImported[Recv_offset[recvTask]],
                       Recv_count[recvTask] * sizeof(struct sinksresultsimported_data), MPI_BYTE, recvTask, TAG_FOF_A,
                       &tmp_results[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct sinksresultsimported_data), MPI_BYTE,
                       recvTask, TAG_FOF_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  for(i = 0; i < nexport; i++)
    {
      int target = tmp_results[i].index;
      if(P[target].AuxDataID >= NumSinks)
        terminate("SINKS: P[target(=%d)].AuxDataID(=%lld) >= NumSinks=%d", target, (long long)P[target].AuxDataID, NumSinks);

      sinks_accreted[P[target].AuxDataID] = tmp_results[i].Prop;
    }

  myfree(tmp_results);
  myfree(Tree_SinksngbResultsImported);

  /* now in case we have swallowed remote sinks out of Tree_Points, need to erase them locally too */
  Tree_SinksngbDeleted = mymalloc("Tree_SinksngbDeleted", sizeof(struct sinksdeleted_data) * N_Sinks_swallowed_imported);

  for(j = 0; j < NTask; j++)
    Recv_count[j] = 0;

  for(i = 0, n = 0, k = 0; i < NTask; i++)
    for(j = 0; j < Force_Recv_count[i]; j++, n++)
      {
        if(Tree_Points[n].Type == 5 && treeSinksexportflag[Tree_Points[n].AuxDataIndex] == 1)
          {
            if(k >= N_Sinks_swallowed_imported)
              terminate("SINKS: k >= N_Sinks_swallowed k=%d N_Sinks_swallowed_imported=%d n=%d i=%d Tree_NumPartImported=%d", k,
                        N_Sinks_swallowed_imported, n, i, Tree_NumPartImported);

            Tree_SinksngbDeleted[k].index = Tree_Points[n].index;

            Recv_count[i]++;
            k++;
          }
      }

  MPI_Alltoall(Recv_count, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nexport = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nexport += Send_count[j];
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  struct sinksdeleted_data *tmp_deleted = mymalloc("tmp_deleted", nexport * sizeof(struct sinksdeleted_data));

  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&Tree_SinksngbDeleted[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct sinksdeleted_data), MPI_BYTE,
                       recvTask, TAG_FOF_A, &tmp_deleted[Send_offset[recvTask]],
                       Send_count[recvTask] * sizeof(struct sinksdeleted_data), MPI_BYTE, recvTask, TAG_FOF_A, MPI_COMM_WORLD,
                       MPI_STATUS_IGNORE);
    }

  for(i = 0; i < nexport; i++)
    {
      int no = tmp_deleted[i].index;
      if(P[no].Type != 5)
        terminate("SINKS: P[no].Type != 5 when deleting sinks");

      P[no].Mass = 0;
      P[no].ID   = 0;
    }

  myfree(tmp_deleted);
  myfree(Tree_SinksngbDeleted);

  myfree(Tree_ResultIndexList);
  myfree(TargetList);

  /* now update momentum of Sinks */
  for(idx = 0; idx < TimeBinsSinksAccretion.NActiveParticles; idx++)
    {
      n = TimeBinsSinksAccretion.ActiveParticleList[idx];
      if(n < 0)
        continue;

      if(P[n].AuxDataID >= NumSinks)
        terminate("SINKS: P[n(=%d)].AuxDataID(=%lld) >= NumSinks=%d", n, (long long)P[n].AuxDataID, NumSinks);

      if(sinks_accreted[P[n].AuxDataID].Mass > 0)
        {
          for(k = 0; k < 3; k++)
            P[n].Vel[k] = (P[n].Vel[k] * P[n].Mass + sinks_accreted[P[n].AuxDataID].Momentum[k]) /
                          (P[n].Mass + sinks_accreted[P[n].AuxDataID].Mass);

          P[n].Mass += sinks_accreted[P[n].AuxDataID].Mass;

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
          if(((1 << P[n].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
            P[n].SofteningType = get_softening_type_from_mass(P[n].Mass);
#endif

          SKP(n).Sinks_CountProgs += sinks_accreted[P[n].AuxDataID].Sinks_CountProgs;
        }
    }

  timebin_cleanup_list_of_active_particles(&TimeBinsSinksAccretion);
  timebin_cleanup_list_of_active_particles(&TimeBinsGravity);

  MPI_Reduce(&N_Sinks_swallowed_imported, &Ntot_Sinks_swallowed_imported, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&N_Sinks_swallowed_local, &Ntot_Sinks_swallowed_local, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  myfree(sinks_accreted);

  myfree(treeSinksexportflag);

  mpi_printf("SINKS: Sinks merging done: %d (%d/%d) sinks particles swallowed\n",
             Ntot_Sinks_swallowed_local + Ntot_Sinks_swallowed_imported, Ntot_Sinks_swallowed_local, Ntot_Sinks_swallowed_imported);

  TIMER_STOP(CPU_SINKS_MERGERS);
}

static int sinks_tree_merger_evaluate(int target, int mode, int threadid)
{
  int k, no, numnodes, *firstnode;
  double h, h2;
  double dx, dy, dz, r2;
  MyIDType id;
  MyDouble *pos;
  MyDouble mass;

  data_in local, *in;
  data_out out;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      in = &local;

      numnodes  = 1;
      firstnode = NULL;
    }
  else
    {
      in = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  pos  = in->Pos;
  mass = in->Mass;
  h    = in->Sinks_Hsml * All.HubbleParam / All.cf_atime;
  id   = in->ID;

  h2 = h * h;

  struct sinks_properties accreted;
  memset(&accreted, 0, sizeof(struct sinks_properties));

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
                  if(P[no].Type == 5)
                    {
                      if((SKP(no).SwallowID == id) && (id != 0) && (P[no].ID != 0)) /* we have a sink merger */
                        {
                          fprintf(FdSinksMergers, "%d %g %llu %g %llu %g\n", ThisTask, All.Time, (long long)id, mass,
                                  (long long)P[no].ID, P[no].Mass);
                          myflush(FdSinksMergers);

                          accreted.Mass += P[no].Mass;

                          for(k = 0; k < 3; k++)
                            accreted.Momentum[k] += P[no].Mass * P[no].Vel[k];

                          accreted.Sinks_CountProgs += SKP(no).Sinks_CountProgs;

                          P[no].Mass = 0;
                          P[no].ID   = 0;

                          N_Sinks_swallowed_local++;
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
                      if((TSKP(n).SwallowID == id) && (id != 0) && (TSKP(n).ID != 0)) /* we have a sink merger */
                        {
                          fprintf(FdSinksMergers, "%d %g %llu %g %llu %g\n", ThisTask, All.Time, (long long)id, mass,
                                  (long long)TSKP(n).ID, Tree_Points[n].Mass);
                          myflush(FdSinksMergers);

                          accreted.Mass += Tree_Points[n].Mass;

                          for(k = 0; k < 3; k++)
                            accreted.Momentum[k] += Tree_Points[n].Mass * TSKP(n).Vel[k];

                          accreted.Sinks_CountProgs += TSKP(n).Sinks_CountProgs;

                          TSKP(n).ID                                       = 0;
                          treeSinksexportflag[Tree_Points[n].AuxDataIndex] = 1;

                          Tree_Points[n].Mass = 0;

                          N_Sinks_swallowed_imported++;
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

  out.Prop = accreted;

  /* store result at the proper place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target].Prop = accreted;

  return 0;
}

#endif
#endif
