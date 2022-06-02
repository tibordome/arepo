/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/domain_rearrange.c
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

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "allvars.h"
#include "domain.h"
#include "proto.h"
#include "voronoi.h"

/*! \brief Gets rid of inactive/eliminated cells and particles.
 *
 *  Cells that were de-refined or turned into star particles are kept in the
 *  SphP array, but flagged as inactive until this point. This routine cleans
 *  up these arrays in order to make sure only active particles/cells are
 *  exported.
 *
 *  \return void
 */
void domain_rearrange_particle_sequence(void)
{
#if defined(USE_SFR) || defined(SINK_PARTICLES)

#ifdef USE_SFR
  if(Stars_converted)
#else
  if(SinksFormedSinceLastDomain)
#endif
    {
      struct particle_data psave;
      peanokey key;

      for(int i = 0; i < NumGas; i++)
        if(P[i].Type != 0) /*If not a gas particle, swap to the end of the list */
          {
            psave = P[i];
            key   = Key[i];

            P[i]   = P[NumGas - 1];
            Key[i] = Key[NumGas - 1];

            SphP[i] = SphP[NumGas - 1];

#ifdef CHIMES
            if(i < NumGas - 1)
              {
                SphP[i].ChimesGasVars.index = i;
                chimes_move_abundances(i, NumGas - 1);
              }
            chimes_set_pointers_to_null(NumGas - 1);
#endif

            P[NumGas - 1]   = psave;
            Key[NumGas - 1] = key;

#if defined(GFM) || defined(SFR_MCS)
            if(P[NumGas - 1].Type == 4)
              StarP[P[NumGas - 1].AuxDataID].PID = NumGas - 1;
#endif
#ifdef BLACK_HOLES
            if(P[NumGas - 1].Type == 5)
              BHP[P[NumGas - 1].AuxDataID].PID = NumGas - 1;
#endif
#ifdef DUST_LIVE
            if(P[NumGas - 1].Type == DUST_LIVE)
              DustP[P[NumGas - 1].AuxDataID].PID = NumGas - 1;
#endif
            NumGas--;
            i--;
          }
          /*Now we have rearranged the particles,
           *we don't need to do it again unless there are more stars*/
#ifdef USE_SFR
      Stars_converted = 0;
#endif
#ifdef SINK_PARTICLES
      SinksFormedSinceLastDomain = 0;
#endif
    }
#endif

#if defined(BLACK_HOLES) || defined(REFINEMENT_MERGE_CELLS) || defined(GFM_WINDS) || defined(GFM_WINDS_LOCAL) || defined(SINKS) || \
    defined(SINK_PARTICLES) || defined(DUST_LIVE) || defined(SFR_MCS)
  int i, count_elim, count_gaselim, count_BHelim, count_windelim, count_dustelim, count_sinkselim;

  count_elim      = 0;
  count_gaselim   = 0;
  count_BHelim    = 0;
  count_windelim  = 0; /* For SFR_MCS, these are star not wind particles */
  count_dustelim  = 0;
  count_sinkselim = 0;

  for(i = 0; i < NumPart; i++)
#ifdef DUST_LIVE
    if((P[i].Mass == 0 && P[i].ID == 0) || (P[i].Type == PTYPE_STARS && P[i].Mass == 0) || (P[i].Type == DUST_LIVE && P[i].Mass == 0))
#elif defined(GFM_RPROCESS_CHANNELS_NS_KICKS)
    if((P[i].Mass == 0 && P[i].ID == 0) || (P[i].Type == PTYPE_STARS && P[i].Mass == 0 && STP(i).NSNS_channel == -1))
#else
    if((P[i].Mass == 0 && P[i].ID == 0) || (P[i].Type == PTYPE_STARS && P[i].Mass == 0))
#endif
      {
#ifdef TRACER_MC
        if(P[i].NumberOfTracers > 0)
          terminate("have found a particle with no mass but tracers\n");
#endif

        if(P[i].Type == PTYPE_GAS)
          {
            P[i]    = P[NumGas - 1];
            SphP[i] = SphP[NumGas - 1];

#ifdef CHIMES
            if(i < NumGas - 1)
              {
                SphP[i].ChimesGasVars.index = i;
                chimes_move_abundances(i, NumGas - 1);
              }
            chimes_set_pointers_to_null(NumGas - 1);
#endif

            Key[i] = Key[NumGas - 1];

            P[NumGas - 1]   = P[NumPart - 1];
            Key[NumGas - 1] = Key[NumPart - 1];

#if defined(GFM) || defined(SFR_MCS)
            if(P[NumGas - 1].Type == 4)
              StarP[P[NumGas - 1].AuxDataID].PID = NumGas - 1;
#endif
#ifdef BLACK_HOLES
            if(P[NumGas - 1].Type == 5)
              BHP[P[NumGas - 1].AuxDataID].PID = NumGas - 1;
#endif
#ifdef SINKS
            if(P[NumGas - 1].Type == 5)
              SinkP[P[NumGas - 1].AuxDataID].PID = NumGas - 1;
#endif
#ifdef DUST_LIVE
            if(P[NumGas - 1].Type == DUST_LIVE)
              DustP[P[NumGas - 1].AuxDataID].PID = NumGas - 1;
#endif
            NumGas--;
            count_gaselim++;
          }
#if defined(GFM) || defined(SFR_MCS)
        else if(P[i].Type == PTYPE_STARS)
          {
            StarP[P[i].AuxDataID]              = StarP[N_star - 1];
            P[StarP[N_star - 1].PID].AuxDataID = P[i].AuxDataID;

            if(i < NumPart - 1)
              {
                P[i]   = P[NumPart - 1];
                Key[i] = Key[NumPart - 1];

                if(P[i].Type == PTYPE_STARS)
                  StarP[P[i].AuxDataID].PID = i;
#ifdef BLACK_HOLES
                if(P[i].Type == 5)
                  BHP[P[i].AuxDataID].PID = i;
#endif
#ifdef SINKS
                if(P[i].Type == 5)
                  SinkP[P[i].AuxDataID].PID = i;
#endif
#ifdef DUST_LIVE
                if(P[i].Type == DUST_LIVE)
                  DustP[P[i].AuxDataID].PID = i;
#endif
              }
            N_star--;
            count_windelim++;
          }
#endif /* GFM || SFR_MCS */

#ifdef BLACK_HOLES
        else if(P[i].Type == 5)
          {
            BHP[P[i].AuxDataID]              = BHP[NumBHs - 1];
            P[BHP[NumBHs - 1].PID].AuxDataID = P[i].AuxDataID;

            if(i < NumPart - 1)
              {
                P[i]   = P[NumPart - 1];
                Key[i] = Key[NumPart - 1];

#if defined(GFM) || defined(SFR_MCS)
                if(P[i].Type == PTYPE_STARS)
                  StarP[P[i].AuxDataID].PID = i;
#endif
                if(P[i].Type == 5)
                  BHP[P[i].AuxDataID].PID = i;
#ifdef SINKS
                if(P[i].Type == 5)
                  SinkP[P[i].AuxDataID].PID = i;
#endif
#ifdef DUST_LIVE
                if(P[i].Type == DUST_LIVE)
                  DustP[P[i].AuxDataID].PID = i;
#endif
              }

            NumBHs--;
            count_BHelim++;
          }
#endif /* BLACK_HOLES */
#ifdef SINKS
        else if(P[i].Type == 5)
          {
            SinkP[P[i].AuxDataID]                = SinkP[NumSinks - 1];
            P[SinkP[NumSinks - 1].PID].AuxDataID = P[i].AuxDataID;

            if(i < NumPart - 1)
              {
                P[i]   = P[NumPart - 1];
                Key[i] = Key[NumPart - 1];

#if defined(GFM) || defined(SFR_MCS)
                if(P[i].Type == PTYPE_STARS)
                  StarP[P[i].AuxDataID].PID = i;
#endif
#ifdef BLACK_HOLES
                if(P[i].Type == 5)
                  BHP[P[i].AuxDataID].PID = i;
#endif
                if(P[i].Type == 5)
                  SinkP[P[i].AuxDataID].PID = i;
#ifdef DUST_LIVE
                if(P[i].Type == DUST_LIVE)
                  DustP[P[i].AuxDataID].PID = i;
#endif
              }

            NumSinks--;
            count_sinkselim++;
          }
#endif /* SINKS */
#ifdef DUST_LIVE
        else if(P[i].Type == DUST_LIVE)
          {
            DustP[P[i].AuxDataID]              = DustP[N_dust - 1];
            P[DustP[N_dust - 1].PID].AuxDataID = P[i].AuxDataID;

            if(i < NumPart - 1)
              {
                P[i]   = P[NumPart - 1];
                Key[i] = Key[NumPart - 1];

#if defined(GFM) || defined(SFR_MCS)
                if(P[i].Type == PTYPE_STARS)
                  StarP[P[i].AuxDataID].PID = i;
#endif
#ifdef BLACK_HOLES
                if(P[i].Type == 5)
                  BHP[P[i].AuxDataID].PID = i;
#endif
#ifdef SINKS
                if(P[i].Type == 5)
                  SinkP[P[i].AuxDataID].PID = i;
#endif
                if(P[i].Type == DUST_LIVE)
                  DustP[P[i].AuxDataID].PID = i;
              }

            N_dust--;
            count_dustelim++;
          }
#endif /* DUST_LIVE */
#ifdef SINK_PARTICLES
        else if(P[i].Type == 5)
          {
            P[i] = P[NumPart - 1];
            // TODO: check if it is ok like this or if we need to update the index of the SinkP structure of the sink eventually
            // associated to P[Numpart - 1]
          }
#endif

        /* TAKE CARE: This code is executed for any of the #ifdef options given
         * at the top. However, without further changes, this eliminates
         * particles from the end of the array *without* swapping the actual
         * particle to be eliminated to the end first, which is *wrong*!
         * Any #ifdef option added at the top must also add behavior to deal
         * with proper elimination of particles before this point! */
        NumPart--;
        i--;

        count_elim++;
      }

#ifdef SINK_PARTICLES
  for(int i = 0; i < NSinksAllTasks; i++)
    if(SinkP[i].ID == 0 && SinkP[i].Mass == 0)
      {
        SinkP[i] = SinkP[NSinksAllTasks - 1];

        int candidate_task = SinkP[i].HomeTask;
        if(ThisTask == candidate_task)
          {
            NSinksThisTask--;

            count_sinkselim++;
            count_elim++;
          }
        NSinksAllTasks--;
        i--;
      }
#endif

  int count[6] = {count_elim, count_gaselim, 0, 0, 0, 0};
  int tot[6] = {0, 0, 0, 0, 0, 0}, nelem = 2;
#ifdef BLACK_HOLES
  count[2] = count_BHelim;
  nelem    = 3;
#endif
#if defined(GFM) || defined(SFR_MCS)
  count[3] = count_windelim;
  nelem    = 4;
#endif
#ifdef DUST_LIVE
  count[4] = count_dustelim;
  nelem    = 5;
#endif
#if defined(SINK_PARTICLES) || defined(SINKS)
  count[5] = count_sinkselim;
  nelem    = 6;
#endif

  MPI_Allreduce(count, tot, nelem, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
#if !defined(DUST_LIVE) && !defined(SFR_MCS)
      printf(
          "DOMAIN: Eliminated %d derefined/swallowed gas cells, merged away %d black holes/sinks, removed %d recoupled wind "
          "particles.\n",
          tot[1], tot[2] + tot[5], tot[3]);
#else
#ifdef SFR_MCS
      printf("DOMAIN: Eliminated %d derefined/swallowed gas cells, merged away %d black holes/sinks, removed %d star particles.\n",
             tot[1], tot[2] + tot[5], tot[3]);
#else
      printf(
          "DOMAIN: Eliminated %d derefined/swallowed gas cells, merged away %d black holes/sinks, removed %d recoupled wind "
          "particles, removed %d dust particles.\n",
          tot[1], tot[2] + tot[5], tot[3], tot[4]);
#endif
#ifdef SINK_PARTICLES
      printf("DOMAIN: Eliminated %d sink particles\n", tot[5]);
#endif
#endif
      myflush(stdout);
    }

  All.TotNumPart -= tot[0];
  All.TotNumGas -= tot[1];
#ifdef BLACK_HOLES
  All.TotNumBHs -= tot[2];
#endif
#if defined(GFM) || defined(SFR_MCS)
  All.TotN_star -= tot[3];
#endif
#ifdef DUST_LIVE
  All.TotN_dust -= tot[4];
#endif
#ifdef SINKS
  All.TotN_sinks -= tot[5];
#endif
#endif
}
