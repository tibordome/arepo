/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sinks/sinks_accrete.c
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       Accretion routines for sink particles
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - 04.2017 Federico Marinacci: major rewriting to use the new communication scheme
 *   the new estimate for the accretion rate in src/sinks/sinks_dmass.c
 */

#include "../allvars.h"
#include "../proto.h"

#ifdef SINKS

static double MassAccreted, TotMassAccreted;

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyDouble hsml;
  MyDouble dMass;
  MyDouble MassNorm;
  signed char DtiStep;
  // int TimeBin;

  int Firstnode;
} data_in;

static data_in *DataGet;

/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in *in, int i, int firstnode)
{
  int k;
  // integertime DtMax;
  signed char DtMax;

  for(k = 0; k < 3; k++)
    in->Pos[k] = P[SinksAux[i].SinksAuxID].Pos[k];

  in->hsml = SKP(SinksAux[i].SinksAuxID).Sinks_Hsml;

  DtMax = SinksAux[i].MaxTimeBin;

  // in->TimeBin = P[SinksAux[i].SinksAuxID].TimeBinGrav;
  // in->dMass = SinksAux[i].dMass / SinksAux[i].MaxTimeBin;
  // in->MassNorm = SinksAux[i].MassNorm;

  if(SinksAux[i].MaxTimeBin > 0)
    {
      in->dMass   = SinksAux[i].dMass / DtMax;
      in->DtiStep = imin(P[SinksAux[i].SinksAuxID].TimeBinGrav, DtMax);
    }
  else
    {
      in->dMass   = 0;
      in->DtiStep = P[SinksAux[i].SinksAuxID].TimeBinGrav;
    }

  in->MassNorm = SinksAux[i].MassNorm;

  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
  MyFloat AccretedMass;
  MyFloat AccretedMomentum[3];

} data_out;

static data_out *DataResult;

/* routine to store or combine result data */
static void out2particle(data_out *out, int i, int mode)
{
  int k;

  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
      for(k = 0; k < 3; k++)
        P[SinksAux[i].SinksAuxID].Vel[k] =
            (P[SinksAux[i].SinksAuxID].Vel[k] * P[SinksAux[i].SinksAuxID].Mass + out->AccretedMomentum[k]) /
            (P[SinksAux[i].SinksAuxID].Mass + out->AccretedMass);

      P[SinksAux[i].SinksAuxID].Mass += out->AccretedMass;

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
      if(((1 << P[SinksAux[i].SinksAuxID].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
        P[SinksAux[i].SinksAuxID].SofteningType = get_softening_type_from_mass(P[SinksAux[i].SinksAuxID].Mass);
#endif
    }
  else /* combine */
    {
      for(k = 0; k < 3; k++)
        P[SinksAux[i].SinksAuxID].Vel[k] =
            (P[SinksAux[i].SinksAuxID].Vel[k] * P[SinksAux[i].SinksAuxID].Mass + out->AccretedMomentum[k]) /
            (P[SinksAux[i].SinksAuxID].Mass + out->AccretedMass);

      P[SinksAux[i].SinksAuxID].Mass += out->AccretedMass;

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
      if(((1 << P[SinksAux[i].SinksAuxID].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
        P[SinksAux[i].SinksAuxID].SofteningType = get_softening_type_from_mass(P[SinksAux[i].SinksAuxID].Mass);
#endif
    }
}

static int sinks_accrete_evaluate(int target, int mode, int threadid);

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

        sinks_accrete_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        sinks_accrete_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void sinks_accrete(void)
{
  MassAccreted = 0;

  generic_set_MaxNexport();

  generic_comm_pattern(NSinks, kernel_local, kernel_imported);

  MPI_Reduce(&MassAccreted, &TotMassAccreted, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  // set sink particle accretion time bin
  for(int i = 0; i < NSinks; i++)
    if(SinksAux[i].MinTimeBin < TIMEBINS)
      P[SinksAux[i].SinksAuxID].TimeBinSink = SinksAux[i].MinTimeBin;

  mpi_printf("SINKS: Accreted mass %g.\n", TotMassAccreted);
}

static int sinks_accrete_evaluate(int target, int mode, int threadid)
{
  int j, k, n, numnodes, *firstnode;
  MyFloat accreted_mass, accreted_momentum[3];
  MyDouble *pos;
  double r2, rad, rad2, dx, dy, dz, fac, dmass, norm;

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

  pos   = in->Pos;
  dmass = in->dMass * in->DtiStep;
  // dmass = in->dMass * in->TimeBin;
  norm = in->MassNorm;

  rad  = in->hsml / sqrt(SKD.DistFac);
  rad2 = rad * rad;

  accreted_mass = 0;

  for(j = 0; j < 3; j++)
    accreted_momentum[j] = 0;

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
              if(P[j].Mass < MIN_TARGET_MASS_FACTOR_FOR_ACC * All.TargetGasMass)
                continue;

              double accmass = dmass * P[j].Mass / norm;

              if((dmass > 0) && (norm == 0))
                terminate("Accreted mass is non zero, but normalization is! dmass=%g, norm %g", dmass, norm);

              /* there are no gas particles in the accretion radius, nothing to accrete */
              if(dmass == 0 && norm == 0)
                accmass = 0;

              MassAccreted += accmass;

              fac = (P[j].Mass - accmass) / P[j].Mass;

              // limit accretion
              if(fac < ACCRETION_SAFETY_FACTOR)
                fac = ACCRETION_SAFETY_FACTOR;

              accreted_mass += (1. - fac) * P[j].Mass;

              for(k = 0; k < 3; k++)
                accreted_momentum[k] += (1. - fac) * P[j].Mass * P[j].Vel[k];

              P[j].Mass *= fac;
              SphP[j].Energy *= fac;
              SphP[j].Momentum[0] *= fac;
              SphP[j].Momentum[1] *= fac;
              SphP[j].Momentum[2] *= fac;

#ifdef MAXSCALARS
              for(int s = 0; s < N_Scalar; s++)
                *(MyFloat *)(((char *)(&SphP[j])) + scalar_elements[s].offset_mass) *= fac;
#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
              SphP[j].Entropy *= fac;
#endif
            }
        }
    }

  out.AccretedMass = accreted_mass;
  for(k = 0; k < 3; k++)
    out.AccretedMomentum[k] = accreted_momentum[k];

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif
