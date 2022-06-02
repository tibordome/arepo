/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/dust_derefinement.c
 * \date        MM/YYYY
 * \author      Ryan McKinnon
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

#include "../allvars.h"
#include "../proto.h"

#ifdef DUST_LIVE
#ifdef DL_DEREFINEMENT

#if defined(HIERARCHICAL_GRAVITY) || defined(ALLOW_DIRECT_SUMMATION)
#error \
    "Use of option DL_DEREFINEMENT does not work with HIERARCHICAL_GRAVITY or ALLOW_DIRECT_SUMMATION since we need a full gravity tree for neighbor searches at each time!"
#endif

/* Helper struct used to actually transfer the contents of a particle to be
 * derefined to its destination task containing P. */
static struct derefineimported_data
{
  int ClosestDustIndex;
  int ClosestDustTask;
  MyIDType ClosestDustID;
  int IsDerefinementTarget;

  MyDouble Mass;
  MyDouble Pos[3];
  MyFloat Vel[3];
  MyFloat NumGrains[DL_GRAIN_BINS];
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
  MyFloat BinSlopes[DL_GRAIN_BINS];
#endif
  MyFloat MetalFractions[GFM_N_CHEM_ELEMENTS];
  MyFloat BinMassChgTau;
  MyFloat OldBinMassChgTau;
} * Tree_DerefineImported;

/* Helper struct used to notify the destination of a particle to be derefined,
 * sending info to the task containing the destination's tree entry. */
static struct simple_derefineimported_data
{
  int ClosestDustIndexTree;
  int ClosestDustTaskTree;
  int ClosestDustHasP;
} * Tree_SimpleDerefineImported;

int dust_compare_tasktree(const void *a, const void *b)
{
  if(((struct simple_derefineimported_data *)a)->ClosestDustTaskTree < ((struct simple_derefineimported_data *)b)->ClosestDustTaskTree)
    return -1;

  if(((struct simple_derefineimported_data *)a)->ClosestDustTaskTree > ((struct simple_derefineimported_data *)b)->ClosestDustTaskTree)
    return 1;

  return 0;
}

int dust_compare_task(const void *a, const void *b)
{
  if(((struct derefineimported_data *)a)->ClosestDustTask < ((struct derefineimported_data *)b)->ClosestDustTask)
    return -1;

  if(((struct derefineimported_data *)a)->ClosestDustTask > ((struct derefineimported_data *)b)->ClosestDustTask)
    return 1;

  return 0;
}

void do_derefinement_notification(struct simple_derefineimported_data *in)
{
  int target = in->ClosestDustIndexTree;

  if(in->ClosestDustHasP == 1)
    {
      DTP(target).IsDerefinementTarget = 1;
    }
  else
    {
      Tree_Points[target].IsDerefinementTarget = 1;
    }
}

/* Take the necessary dust particle information from a particle to be
 * derefined, and add it to the particle info of its closest dust neighbor,
 * which lies on this task.
 *
 * This is an entirely local process, because either the two particles both lie
 * on this task, or the remote small particle's info was sent to this task. */
void do_derefinement_transfer(struct derefineimported_data *in)
{
  int target = in->ClosestDustIndex;
  if(P[target].ID != in->ClosestDustID)
    {
      print_particle_info(target);
      terminate("DUST_LIVE: Trying to derefine into the wrong dust particle, P.ID=%llu and ClosestDustID=%llu!\n",
                (unsigned long long)P[target].ID, (unsigned long long)in->ClosestDustID);
    }

  /* Mark that this dust particle was scheduled to be the target of a
   * derefinement. */
  DTP(target).IsDerefinementTarget = 1;
  /* However, do not allow the derefinement to proceed if the incoming particle
   * itself is the target of a derefinement. */
  if(in->IsDerefinementTarget == 1)
    {
      return;
    }

  double m_small = in->Mass;
  double m_big   = P[target].Mass;
  double m_tot   = m_small + m_big;

  P[target].Mass = m_tot;
  for(int j = 0; j < 3; j++)
    {
      double mr_tot    = m_small * in->Pos[j] + m_big * P[target].Pos[j];
      double p_tot     = m_small * in->Vel[j] + m_big * P[target].Vel[j];
      P[target].Pos[j] = mr_tot / m_tot;
      P[target].Vel[j] = p_tot / m_tot;
    }

  for(int j = 0; j < DL_GRAIN_BINS; j++)
    {
      DTP(target).NumGrains[j] += in->NumGrains[j];
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
      DTP(target).BinSlopes[j] += in->BinSlopes[j];
#endif
    }

  for(int j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
    {
      double m_met                  = m_small * in->MetalFractions[j] + m_big * DTP(target).MetalFractions[j];
      DTP(target).MetalFractions[j] = m_met / m_tot;
    }

  /* Needed for calculation of grain size evolution timescales. */
  DTP(target).OrigMass += m_small;
  DTP(target).BinMassChgTau    = fmin(DTP(target).BinMassChgTau, in->BinMassChgTau);
  DTP(target).OldBinMassChgTau = fmin(DTP(target).OldBinMassChgTau, in->OldBinMassChgTau);
  // TODO: timebins? smoothing lengths? etc.
}

/* First, notify the destinations of desired dust particle derefinements that a
 * particle would like to derefine into them.
 *
 * Then, attempt to do each derefinement, skipping ones where the derefining
 * particle is itself the target of a derefinement. */
void exchange_derefine_results(void)
{
  int idx, i, j, k, nexport, nimport;
  int ngrp, recvTask, n;
  int ncount, ncount1, ncount2;
  struct simple_derefineimported_data *tmp_simple_results;
  struct derefineimported_data *tmp_results;

  /* Send messages to closest neighbors */
  for(idx = 0, ncount = 0, ncount1 = 0; idx < TimeBinsDust.NActiveParticles; idx++)
    {
      i = TimeBinsDust.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if((P[i].Type == DUST_LIVE) && (Tree_Task_list[i] == ThisTask) && (P[i].Mass > 0.0))
        {
          if((P[i].Mass < All.DustMinFrac * All.TargetGasMass) && (PShatter[ncount1].HighMassNeighbor == 1) &&
             (PShatter[ncount1].ClosestDustTaskTree != ThisTask))
            {
              ncount++;
            }
          ncount1++;
        }
    }

  for(i = 0; i < Tree_NumPartImported; i++)
#ifndef HIERARCHICAL_GRAVITY
    if(Tree_Points[i].ActiveFlag)
#endif
      if((Tree_Points[i].Type == DUST_LIVE) && (Tree_Points[i].Mass > 0.0))
        {
          if((Tree_Points[i].Mass < All.DustMinFrac * All.TargetGasMass) && (PShatter[ncount1].HighMassNeighbor == 1) &&
             (PShatter[ncount1].ClosestDustTaskTree != ThisTask))
            {
              ncount++;
            }
          ncount1++;
        }

  Tree_SimpleDerefineImported = (struct simple_derefineimported_data *)mymalloc("Tree_SimpleDerefineImported",
                                                                                ncount * sizeof(struct simple_derefineimported_data));

  memset(Tree_SimpleDerefineImported, 0, ncount * sizeof(struct simple_derefineimported_data));

  for(idx = 0, ncount1 = 0, k = 0; idx < TimeBinsDust.NActiveParticles; idx++)
    {
      i = TimeBinsDust.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if((P[i].Type == DUST_LIVE) && (Tree_Task_list[i] == ThisTask) && (P[i].Mass > 0.0))
        {
          if((P[i].Mass < All.DustMinFrac * All.TargetGasMass) && (PShatter[ncount1].HighMassNeighbor == 1))
            {
              struct simple_derefineimported_data tmpin;

              /* These particles in the tree have direct access to P structure */
              tmpin.ClosestDustIndexTree = PShatter[ncount1].ClosestDustIndexTree;
              tmpin.ClosestDustTaskTree  = PShatter[ncount1].ClosestDustTaskTree;
              tmpin.ClosestDustHasP      = PShatter[ncount1].ClosestDustHasP;

              PShatter[ncount1].IsDerefinementOrigin = 1;

              /* Receiving dust particle is on this task, so do the update */
              if(PShatter[ncount1].ClosestDustTaskTree == ThisTask)
                {
                  do_derefinement_notification(&tmpin);
                }
              /* Receiving dust particle is on another task, so prepare for export */
              else
                {
                  Tree_SimpleDerefineImported[k] = tmpin;
                  k++;
                }
            }

          ncount1++;
        }
    }

  for(i = 0, ncount2 = 0; i < Tree_NumPartImported; i++)
#ifndef HIERARCHICAL_GRAVITY
    if(Tree_Points[i].ActiveFlag)
#endif
      if((Tree_Points[i].Type == DUST_LIVE) && (Tree_Points[i].Mass > 0.0))
        {
          if((Tree_Points[i].Mass < All.DustMinFrac * All.TargetGasMass) && (PShatter[ncount1].HighMassNeighbor == 1))
            {
              struct simple_derefineimported_data tmpin;

              /* These particles in the tree do not have direct access to P structure */
              tmpin.ClosestDustIndexTree = PShatter[ncount1].ClosestDustIndexTree;
              tmpin.ClosestDustTaskTree  = PShatter[ncount1].ClosestDustTaskTree;
              tmpin.ClosestDustHasP      = PShatter[ncount1].ClosestDustHasP;

              PShatter[ncount1].IsDerefinementOrigin = 1;

              /* Receiving dust particle is on this task, so do the update */
              if(PShatter[ncount1].ClosestDustTaskTree == ThisTask)
                {
                  do_derefinement_notification(&tmpin);
                }
              /* Receiving dust particle is on another task, so prepare for export */
              else
                {
                  Tree_SimpleDerefineImported[k] = tmpin;
                  k++;
                }
            }

          ncount1++;
          ncount2++;
        }

  if(k != ncount)
    terminate("DUST_LIVE: k != ncount when preparing dust derefinement notification data for export, k=%d, ncount=%d, ThisTask=%d!\n",
              k, ncount, ThisTask);

  /* Need to ensure that entries in Tree_SimpleDerefineImported are sorted in
   * ascending order according to the task where they'll be sent, so that we
   * can properly compute offsets. */
  qsort(Tree_SimpleDerefineImported, ncount, sizeof(struct simple_derefineimported_data), dust_compare_tasktree);

  for(j = 0; j < NTask; j++)
    Recv_count[j] = 0;

  for(k = 0; k < ncount; k++)
    {
      int dest_task = Tree_SimpleDerefineImported[k].ClosestDustTaskTree;
      Recv_count[dest_task]++;
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

  tmp_simple_results =
      (struct simple_derefineimported_data *)mymalloc("tmp_simple_results", nexport * sizeof(struct simple_derefineimported_data));
  memset(tmp_simple_results, -1, nexport * sizeof(struct simple_derefineimported_data));

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&Tree_SimpleDerefineImported[Recv_offset[recvTask]],
                       Recv_count[recvTask] * sizeof(struct simple_derefineimported_data), MPI_BYTE, recvTask, TAG_FOF_A,
                       &tmp_simple_results[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct simple_derefineimported_data),
                       MPI_BYTE, recvTask, TAG_FOF_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  for(i = 0; i < nexport; i++)
    {
      do_derefinement_notification(&tmp_simple_results[i]);
    }

  myfree(tmp_simple_results);
  myfree(Tree_SimpleDerefineImported);

  /* By this point, we have notified the destinations of all desired
   * derefinements.  We can try to proceed with derefinements, skipping
   * derefinements that come from dust particles that are themselves targets of
   * a derefinement. */
  for(idx = 0, ncount = 0, ncount1 = 0; idx < TimeBinsDust.NActiveParticles; idx++)
    {
      i = TimeBinsDust.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if((P[i].Type == DUST_LIVE) && (Tree_Task_list[i] == ThisTask) && (P[i].Mass > 0.0))
        {
          if((PShatter[ncount1].IsDerefinementOrigin == 1) && (PShatter[ncount1].ClosestDustTask != ThisTask))
            {
              ncount++;
            }
          ncount1++;
        }
    }

  for(i = 0; i < Tree_NumPartImported; i++)
#ifndef HIERARCHICAL_GRAVITY
    if(Tree_Points[i].ActiveFlag)
#endif
      if((Tree_Points[i].Type == DUST_LIVE) && (Tree_Points[i].Mass > 0.0))
        {
          if((PShatter[ncount1].IsDerefinementOrigin == 1) && (PShatter[ncount1].ClosestDustTask != ThisTask))
            {
              ncount++;
            }
          ncount1++;
        }

  Tree_DerefineImported =
      (struct derefineimported_data *)mymalloc("Tree_DerefineImported", ncount * sizeof(struct derefineimported_data));

  memset(Tree_DerefineImported, 0, ncount * sizeof(struct derefineimported_data));

  for(idx = 0, ncount1 = 0, k = 0; idx < TimeBinsDust.NActiveParticles; idx++)
    {
      i = TimeBinsDust.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if((P[i].Type == DUST_LIVE) && (Tree_Task_list[i] == ThisTask) && (P[i].Mass > 0.0))
        {
          if(PShatter[ncount1].IsDerefinementOrigin == 1)
            {
              struct derefineimported_data tmpin;

              /* These particles in the tree have direct access to P structure */
              tmpin.ClosestDustIndex     = PShatter[ncount1].ClosestDustIndex;
              tmpin.ClosestDustTask      = PShatter[ncount1].ClosestDustTask;
              tmpin.ClosestDustID        = PShatter[ncount1].ClosestDustID;
              tmpin.IsDerefinementTarget = DTP(i).IsDerefinementTarget;
              tmpin.Mass                 = P[i].Mass;
              tmpin.BinMassChgTau        = DTP(i).BinMassChgTau;
              tmpin.OldBinMassChgTau     = DTP(i).OldBinMassChgTau;
              for(int l = 0; l < 3; l++)
                {
                  tmpin.Pos[l] = P[i].Pos[l];
                  tmpin.Vel[l] = P[i].Vel[l];
                }
              for(int l = 0; l < DL_GRAIN_BINS; l++)
                {
                  tmpin.NumGrains[l] = DTP(i).NumGrains[l];
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
                  tmpin.BinSlopes[l] = DTP(i).BinSlopes[l];
#endif
                }
              for(int l = 0; l < GFM_N_CHEM_ELEMENTS; l++)
                {
                  tmpin.MetalFractions[l] = DTP(i).MetalFractions[l];
                }

              /* Receiving dust particle is on this task, so do the update */
              if(PShatter[ncount1].ClosestDustTask == ThisTask)
                {
                  do_derefinement_transfer(&tmpin);
                }
              /* Receiving dust particle is on another task, so prepare for export */
              else
                {
                  Tree_DerefineImported[k] = tmpin;
                  k++;
                }
            }

          ncount1++;
        }
    }

  for(i = 0, ncount2 = 0; i < Tree_NumPartImported; i++)
#ifndef HIERARCHICAL_GRAVITY
    if(Tree_Points[i].ActiveFlag)
#endif
      if((Tree_Points[i].Type == DUST_LIVE) && (Tree_Points[i].Mass > 0.0))
        {
          if(PShatter[ncount1].IsDerefinementOrigin == 1)
            {
              struct derefineimported_data tmpin;

              /* These particles in the tree do not have direct access to P structure */
              tmpin.ClosestDustIndex     = PShatter[ncount1].ClosestDustIndex;
              tmpin.ClosestDustTask      = PShatter[ncount1].ClosestDustTask;
              tmpin.ClosestDustID        = PShatter[ncount1].ClosestDustID;
              tmpin.IsDerefinementTarget = Tree_Points[i].IsDerefinementTarget;
              tmpin.Mass                 = Tree_Points[i].Mass;
              tmpin.BinMassChgTau        = Tree_Points[i].BinMassChgTau;
              tmpin.OldBinMassChgTau     = Tree_Points[i].OldBinMassChgTau;
              for(int l = 0; l < 3; l++)
                {
                  tmpin.Pos[l] = Tree_Points[i].Pos[l];
                  tmpin.Vel[l] = Tree_Points[i].Vel[l];
                }
              for(int l = 0; l < DL_GRAIN_BINS; l++)
                {
                  tmpin.NumGrains[l] = Tree_Points[i].NumGrains[l];
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
                  tmpin.BinSlopes[l] = Tree_Points[i].BinSlopes[l];
#endif
                }
              for(int l = 0; l < GFM_N_CHEM_ELEMENTS; l++)
                {
                  tmpin.MetalFractions[l] = Tree_Points[i].MetalFractions[l];
                }

              /* Receiving dust particle is on this task, so do the update */
              if(PShatter[ncount1].ClosestDustTask == ThisTask)
                {
                  do_derefinement_transfer(&tmpin);
                }
              /* Receiving dust particle is on another task, so prepare for export */
              else
                {
                  Tree_DerefineImported[k] = tmpin;
                  k++;
                }
            }

          ncount1++;
          ncount2++;
        }

  if(k != ncount)
    terminate("DUST_LIVE: k != ncount when preparing dust derefinement data for export, k=%d, ncount=%d, ThisTask=%d!\n", k, ncount,
              ThisTask);

  /* Need to ensure that entries in Tree_DerefineImported are sorted in
   * ascending order according to the task where they'll be sent, so that we
   * can properly compute offsets. */
  qsort(Tree_DerefineImported, ncount, sizeof(struct derefineimported_data), dust_compare_task);

  for(j = 0; j < NTask; j++)
    Recv_count[j] = 0;

  for(k = 0; k < ncount; k++)
    {
      int dest_task = Tree_DerefineImported[k].ClosestDustTask;
      Recv_count[dest_task]++;
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

  tmp_results = (struct derefineimported_data *)mymalloc("tmp_results", nexport * sizeof(struct derefineimported_data));
  memset(tmp_results, -1, nexport * sizeof(struct derefineimported_data));

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&Tree_DerefineImported[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct derefineimported_data),
                       MPI_BYTE, recvTask, TAG_FOF_A, &tmp_results[Send_offset[recvTask]],
                       Send_count[recvTask] * sizeof(struct derefineimported_data), MPI_BYTE, recvTask, TAG_FOF_A, MPI_COMM_WORLD,
                       MPI_STATUS_IGNORE);
    }

  for(i = 0; i < nexport; i++)
    {
      do_derefinement_transfer(&tmp_results[i]);
    }

  myfree(tmp_results);
  myfree(Tree_DerefineImported);
}

#endif
#endif
