/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sinks/sinks_init.c
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       Sink particles auxiliary routines
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - 04.2017 Federico Marinacci: added more functions to initilize auxiliary
 *   structures for sink accretion
 */

#include "../allvars.h"
#include "../proto.h"

#ifdef SINKS

int NSinks, TotNSinks;
struct SinksAux_struct *SinksAux;

void sinks_begrun(void)
{
  CPU_Step[CPU_MISC] += measure_time();

  MPI_Bcast(&SKD, sizeof(struct SKD_struct), MPI_BYTE, 0, MPI_COMM_WORLD);

  SKD.Flag    = 0;
  SKD.AccRad2 = pow(SKD.AccRad, 2);

  CPU_Step[CPU_SINKS] += measure_time();
}

void sinks_get_num_sinks(void)
{
  int i;
  SKD.Flag     = 0;
  SKD.NumSinks = 0;

  if(!SKD.Flag)
    {
      for(i = 0; i < NumPart; i++)
        if((P[i].Type == 5) && (P[i].Mass != 0) && (P[i].ID != 0))
          SKD.NumSinks++;

      MPI_Allreduce(&SKD.NumSinks, &SKD.TotNumSinks, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }

  SKD.Flag = 1;
}

void sinks_get_active_sinks(void)
{
  int i;
  NSinks = 0;

  for(int idx = 0; idx < TimeBinsSinksAccretion.NActiveParticles; idx++)
    {
      i = TimeBinsSinksAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Ti_Current != All.Ti_Current)
        terminate("how can this be?");

      if((P[i].Type == 5) && (P[i].Mass != 0) && (P[i].ID != 0))
        {
          SinksAux[NSinks].SinksAuxID = i;
          SinksAux[NSinks].MinTimeBin = TIMEBINS;
          SinksAux[NSinks].MaxTimeBin = 0;
          SinksAux[NSinks].dMass      = 0.0;
          SinksAux[NSinks].MassNorm   = 0.0;
          SKP(i).Sinks_Hsml           = get_accretion_radius(P[i].Mass, SKP(i).Temp);
          NSinks++;
          printf("SINKS: Init Id %i Mass %g Temp %g Hsml %g\n", i, P[i].Mass, SKP(i).Temp, SKP(i).Sinks_Hsml);
          printf("SINKS: Testing i %i SKP(i).PID %i P[i].ID %i\n", i, SKP(i).PID, P[i].ID);
        }
      else
        terminate("Strange: sinks is not of the correct type type %d, mass %g", P[i].Type, P[i].Mass);
    }

  if(NSinks > SKD.NumSinks)
    terminate("Number of active sinks %d larger than total sink number %d", NSinks, SKD.NumSinks);

  MPI_Allreduce(&NSinks, &TotNSinks, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}

void sinks_begin(void)
{
  /* reset derefinement criterion for active particles within accretion radius */
  if(SKD.TotNumSinks)
    for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
      {
        int i = TimeBinsHydro.ActiveParticleList[idx];
        if(i < 0)
          continue;

        SphP[i].InAccrRadius = 0;
      }

  SinksAux = mymalloc("SinksAux", SKD.NumSinks * sizeof(struct SinksAux_struct));
}

void sinks_end(void) { myfree(SinksAux); }

void sinks_check_AuxDataID_references(void)
{
  mpi_printf("SINKS: Sinks checks...\n");

  for(int i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 5)
        if(SinkP[P[i].AuxDataID].PID != i)
          {
            printf("SinkP broken: %llu %d %d %d %d %d\n", (long long)P[i].AuxDataID, i, SinkP[P[i].AuxDataID].PID, NumGas, NumSinks,
                   NumPart);
            terminate("SinkP[P[i].AuxDataID].PID!=i\n");
          }
    }

  mpi_printf("SINKS: done.\n");
}

#endif
