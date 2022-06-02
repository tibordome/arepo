/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sinks/sinks_dmass.c
 * \date        04/2017
 * \author      Federico Marinacci
 * \brief       Accretion rate estimate for sink particles
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"

#ifdef SINKS

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyDouble hsml;
  signed char SinkStep;

  int Firstnode;
} data_in;

static data_in *DataGet;

static double MassAccRate, TotMassAccRate;

/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in *in, int i, int firstnode)
{
  for(int k = 0; k < 3; k++)
    in->Pos[k] = P[SinksAux[i].SinksAuxID].Pos[k];

  in->hsml = SKP(SinksAux[i].SinksAuxID).Sinks_Hsml;

  in->SinkStep = P[SinksAux[i].SinksAuxID].TimeBinGrav;

  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
  MyFloat AccretedMass;
  MyFloat TotAccretedMass;
  signed char MinTimeBin;
  signed char MaxTimeBin;
} data_out;

static data_out *DataResult;

/* routine to store or combine result data */
static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
      SinksAux[i].dMass      = out->AccretedMass;
      SinksAux[i].MassNorm   = out->TotAccretedMass;
      SinksAux[i].MinTimeBin = out->MinTimeBin;
      SinksAux[i].MaxTimeBin = out->MaxTimeBin;
    }
  else /* combine */
    {
      SinksAux[i].dMass += out->AccretedMass;
      SinksAux[i].MassNorm += out->TotAccretedMass;

      if(out->MinTimeBin < SinksAux[i].MinTimeBin)
        SinksAux[i].MinTimeBin = out->MinTimeBin;

      if(out->MaxTimeBin > SinksAux[i].MaxTimeBin)
        SinksAux[i].MaxTimeBin = out->MaxTimeBin;
    }
}

static int sinks_dmass_evaluate(int target, int mode, int threadid);

#include "../generic_comm_helpers2.h"

static void kernel_local(void)
{
  int i;
#pragma omp parallel private(i)
  {
    int j, threadid = get_thread_num();

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= NSinks)
          break;

        sinks_dmass_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}

static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, count = 0;
#pragma omp parallel private(i)
  {
    int threadid = get_thread_num();

    while(1)
      {
#pragma omp atomic capture
        i = count++;

        if(i >= Nimport)
          break;

        sinks_dmass_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void sinks_dmass(void)
{
  MassAccRate = 0;

  generic_set_MaxNexport();

  generic_comm_pattern(NSinks, kernel_local, kernel_imported);

  /* Imposing a safety factor on the total accreted mass (each sink cannot grow more than MAX_MASS_FAC of its mass) */
  for(int i = 0; i < NSinks; i++)
    {
      double maxdmass   = MAX_MASS_FAC * P[SinksAux[i].SinksAuxID].Mass;
      SinksAux[i].dMass = fmin(maxdmass, SinksAux[i].dMass);
      MassAccRate += SinksAux[i].dMass;
    }

  MPI_Reduce(&MassAccRate, &TotMassAccRate, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  mpi_printf("SINKS: computed mass accretion rate %g \n", TotMassAccRate);
}

static int sinks_dmass_evaluate(int target, int mode, int threadid)
{
  int j, n, numnodes, *firstnode;
  double r2, rad, rad2, dx, dy, dz;
  signed char min_time_bin = TIMEBINS;
  signed char max_time_bin = 0;
  // signed char sinkstep;
  MyFloat accreted_mass, tot_accreted_mass;
  MyDouble *pos;

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

  pos = in->Pos;
  // sinkstep = in->SinkStep;

  rad  = in->hsml / sqrt(SKD.DistFac);
  rad2 = rad * rad;

  accreted_mass     = 0;
  tot_accreted_mass = 0;

  int nfound = ngb_treefind_variable_threads(pos, rad, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      if(P[j].Type == 0)
        {
          dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
          dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
          dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);

          r2 = dx * dx + dy * dy + dz * dz;

          if((r2 < rad2) && (P[j].Mass != 0) && (P[j].ID != 0))
            {
              SphP[j].InAccrRadius = 1;

              if(P[j].Mass < MIN_TARGET_MASS_FACTOR_FOR_ACC * All.TargetGasMass)
                continue;

              if(TimeBinSynchronized[P[j].TimeBinHydro])
                {
                  accreted_mass += P[j].Mass;

                  if(P[j].TimeBinHydro > max_time_bin)
                    max_time_bin = P[j].TimeBinHydro;

                  if(P[j].TimeBinHydro < min_time_bin)
                    min_time_bin = P[j].TimeBinHydro;
                }

              tot_accreted_mass += P[j].Mass;
            }
        }
    }

  out.AccretedMass    = accreted_mass;
  out.TotAccretedMass = tot_accreted_mass;
  out.MinTimeBin      = min_time_bin;
  out.MaxTimeBin      = max_time_bin;

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif
