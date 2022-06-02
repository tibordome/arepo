/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/dust_util.c
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../allvars.h"
#include "../proto.h"

#ifdef DUST_LIVE

void start_dust(void)
{
  TIMER_START(CPU_DUST);

  set_cosmo_factors_for_current_time();

  DustParticle = (struct dust_particle *)mymalloc("DustParticle", N_dust * sizeof(struct dust_particle));

  Ndust = 0;
  for(int idx = 0; idx < TimeBinsDust.NActiveParticles; idx++)
    {
      int i = TimeBinsDust.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Ti_Current != All.Ti_Current)
        {
          terminate("how can this be?");
        }

      if((P[i].Type == DUST_LIVE) && (P[i].Mass > 0))
        {
          DustParticle[Ndust].index      = i;
          DustParticle[Ndust].active_idx = idx;
          DustParticle[Ndust].NumNgb     = 0.0;
          DustParticle[Ndust].NormSph    = 0.0;
          DustParticle[Ndust].TotNgbMass = 0.0;
          DustParticle[Ndust].Dhsmlrho   = 0.0;
          for(int k = 0; k < 3; k++)
            {
              DustParticle[Ndust].LocalGasAccel[k] = 0.0;
              DustParticle[Ndust].LocalGradP[k]    = 0.0;
#ifdef DL_DRAG_BACKREACTION
              DustParticle[Ndust].DragMomentum[k] = 0.0;
#endif
            }
#ifdef DL_DRAG_BACKREACTION
          DustParticle[Ndust].DragDustThermal = 0.0;
#endif
          Ndust++;
        }
    }
}

void end_dust(void)
{
  myfree(DustParticle);

  TIMER_STOP(CPU_DUST);
}

#if(defined(DL_GRAIN_BINS) && (defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION))) || \
    defined(DL_DRAG_BACKREACTION)
void begin_dust_search(void)
{
  int idx, i;

  NforcesDust = 0;
  TargetList  = (int *)mymalloc("TargetList", (NumPart + Tree_NumPartImported) * sizeof(int));

  for(idx = 0; idx < TimeBinsDust.NActiveParticles; idx++)
    {
      i = TimeBinsDust.ActiveParticleList[idx];
      if(i < 0)
        continue;

#ifdef DL_DEREFINEMENT
      /* When we find neighbors of active dust particles, this flag will be set
       * to indicate if another dust wants to derefine into the particle. */
      if((P[i].Type == DUST_LIVE) && (P[i].Mass > 0.0))
        {
          DTP(i).IsDerefinementTarget = 0;
        }
#endif

      if((P[i].Type == DUST_LIVE) && (Tree_Task_list[i] == ThisTask) && (P[i].Mass > 0.0))
        {
          TargetList[NforcesDust++] = i;
        }
    }

  for(i = 0; i < Tree_NumPartImported; i++)
#ifndef HIERARCHICAL_GRAVITY
    if(Tree_Points[i].ActiveFlag)
#endif
      if((Tree_Points[i].Type == DUST_LIVE) && (Tree_Points[i].Mass > 0.0))
        {
#ifdef DL_DEREFINEMENT
          Tree_Points[i].IsDerefinementTarget = 0;
#endif
          TargetList[NforcesDust++] = i + Tree_ImportedNodeOffset;
        }

  PShatter = (struct shatter_data_p *)mymalloc("PShatter", NforcesDust * sizeof(struct shatter_data_p));

  memset(&PShatter[0], 0, NforcesDust * sizeof(struct shatter_data_p));
}

void end_dust_search(void)
{
  myfree(PShatter);
  myfree(TargetList);
}

static struct shatterimported_data
{
  int index;
  MyFloat DustDensity;
  MyFloat DustNumNgb;
  MyFloat MinNumNgbDeviationDust;
  MyFloat MaxNumNgbDeviationDust;
  MyFloat DustHsml;
  MyFloat EnclosedMass;
#ifdef DL_DEREFINEMENT
  MyFloat ClosestDustR;
  MyIDType ClosestDustID;
  int HighMassNeighbor;
#endif
} * Tree_ShatterImported;

void exchange_dust_search_results(void)
{
  int idx, i, j, k, nexport, nimport;
  int ngrp, recvTask, n;
  int ncount, ncount1, ncount2;
  struct shatterimported_data *tmp_results;

  for(i = 0, ncount = 0; i < Tree_NumPartImported; i++)
#ifndef HIERARCHICAL_GRAVITY
    if(Tree_Points[i].ActiveFlag)
#endif
      if((Tree_Points[i].Type == DUST_LIVE) && (Tree_Points[i].Mass > 0.0))
        ncount++;

  Tree_ShatterImported = (struct shatterimported_data *)mymalloc("Tree_ShatterImported", ncount * sizeof(struct shatterimported_data));

  memset(Tree_ShatterImported, 0, ncount * sizeof(struct shatterimported_data));

  for(idx = 0, ncount1 = 0; idx < TimeBinsDust.NActiveParticles; idx++)
    {
      i = TimeBinsDust.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if((P[i].Type == DUST_LIVE) && (Tree_Task_list[i] == ThisTask) && (P[i].Mass > 0.0))
        {
          DTP(i).DustNumNgb             = PShatter[ncount1].DustNumNgb;
          DTP(i).DustHsml               = PShatter[ncount1].DustHsml;
          DTP(i).MinNumNgbDeviationDust = PShatter[ncount1].MinNumNgbDeviationDust;
          DTP(i).MaxNumNgbDeviationDust = PShatter[ncount1].MaxNumNgbDeviationDust;
          DTP(i).DustDensity            = PShatter[ncount1].DustDensity;
          DTP(i).EnclosedMass           = PShatter[ncount1].EnclosedMass;
#ifdef DL_DEREFINEMENT
          DTP(i).ClosestDustR     = PShatter[ncount1].ClosestDustR;
          DTP(i).ClosestDustID    = PShatter[ncount1].ClosestDustID;
          DTP(i).HighMassNeighbor = PShatter[ncount1].HighMassNeighbor;
#endif
          ncount1++;
        }
    }

  for(i = 0, ncount2 = 0; i < Tree_NumPartImported; i++)
#ifndef HIERARCHICAL_GRAVITY
    if(Tree_Points[i].ActiveFlag)
#endif
      if((Tree_Points[i].Type == DUST_LIVE) && (Tree_Points[i].Mass > 0.0))
        {
          Tree_ShatterImported[ncount2].DustNumNgb             = PShatter[ncount1].DustNumNgb;
          Tree_ShatterImported[ncount2].DustHsml               = PShatter[ncount1].DustHsml;
          Tree_ShatterImported[ncount2].MinNumNgbDeviationDust = PShatter[ncount1].MinNumNgbDeviationDust;
          Tree_ShatterImported[ncount2].MaxNumNgbDeviationDust = PShatter[ncount1].MaxNumNgbDeviationDust;
          Tree_ShatterImported[ncount2].DustDensity            = PShatter[ncount1].DustDensity;
          Tree_ShatterImported[ncount2].EnclosedMass           = PShatter[ncount1].EnclosedMass;
#ifdef DL_DEREFINEMENT
          Tree_ShatterImported[ncount2].ClosestDustR     = PShatter[ncount1].ClosestDustR;
          Tree_ShatterImported[ncount2].ClosestDustID    = PShatter[ncount1].ClosestDustID;
          Tree_ShatterImported[ncount2].HighMassNeighbor = PShatter[ncount1].HighMassNeighbor;
#endif
          ncount1++;
          ncount2++;
        }

  for(j = 0; j < NTask; j++)
    Recv_count[j] = 0;

  for(i = 0, n = 0, k = 0; i < NTask; i++)
    for(j = 0; j < Force_Recv_count[i]; j++, n++)
      {
#ifndef HIERARCHICAL_GRAVITY
        if(Tree_Points[n].ActiveFlag)
#endif
          if((Tree_Points[n].Type == DUST_LIVE) && (Tree_Points[n].Mass > 0.0))
            {
              Tree_ShatterImported[k].index = Tree_Points[n].index;
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

  tmp_results = (struct shatterimported_data *)mymalloc("tmp_results", nexport * sizeof(struct shatterimported_data));
  memset(tmp_results, -1, nexport * sizeof(struct shatterimported_data));

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&Tree_ShatterImported[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct shatterimported_data),
                       MPI_BYTE, recvTask, TAG_FOF_A, &tmp_results[Send_offset[recvTask]],
                       Send_count[recvTask] * sizeof(struct shatterimported_data), MPI_BYTE, recvTask, TAG_FOF_A, MPI_COMM_WORLD,
                       MPI_STATUS_IGNORE);
    }

  for(i = 0; i < nexport; i++)
    {
      int target                         = tmp_results[i].index;
      DTP(target).DustDensity            = tmp_results[i].DustDensity;
      DTP(target).DustNumNgb             = tmp_results[i].DustNumNgb;
      DTP(target).MinNumNgbDeviationDust = tmp_results[i].MinNumNgbDeviationDust;
      DTP(target).MaxNumNgbDeviationDust = tmp_results[i].MaxNumNgbDeviationDust;
      DTP(target).DustHsml               = tmp_results[i].DustHsml;
      DTP(target).EnclosedMass           = tmp_results[i].EnclosedMass;
#ifdef DL_DEREFINEMENT
      DTP(target).ClosestDustR     = tmp_results[i].ClosestDustR;
      DTP(target).ClosestDustID    = tmp_results[i].ClosestDustID;
      DTP(target).HighMassNeighbor = tmp_results[i].HighMassNeighbor;
#endif
    }

  myfree(tmp_results);
  myfree(Tree_ShatterImported);
}
#endif

#ifdef DL_GRAIN_BINS
int elem_can_be_dust(int k)
{
  if((k == element_index_Carbon) || (k == element_index_Oxygen) || (k == element_index_Magnesium) || (k == element_index_Silicon) ||
     (k == element_index_Iron))
    {
      return 1;
    }
  /* If GFM_NORMALIZED_METAL_ADVECTION is on, OtherMetals are assumed to not
   * deplete onto dust. */
  return 0;
}
#endif

#endif
