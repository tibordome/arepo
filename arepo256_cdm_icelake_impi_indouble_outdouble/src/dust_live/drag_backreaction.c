/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/drag_backreaction.c
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
#ifdef DL_DRAG_BACKREACTION

typedef struct
{
  MyDouble Pos[3];
  MyFloat DragMomentum[3];
  MyFloat DragDustThermal;
  MyFloat Hsml;
  MyFloat TotNgbMass;
  int Firstnode;
} data_in;

static data_in *DataGet;

typedef struct
{
  char dummy;
} data_out;

static data_out *DataResult;

static void particle2in(data_in *in, int i, int firstnode)
{
  for(int k = 0; k < 3; k++)
    {
      in->Pos[k]          = P[DustParticle[i].index].Pos[k];
      in->DragMomentum[k] = DustParticle[i].DragMomentum[k];
    }
  in->DragDustThermal = DustParticle[i].DragDustThermal;
  in->Hsml            = DTP(DustParticle[i].index).Hsml;
  in->TotNgbMass      = DustParticle[i].TotNgbMass;
  in->Firstnode       = firstnode;
}

static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)
    {
    }
  else
    {
    }
}

#include "../generic_comm_helpers2.h"

static void kernel_local(void)
{
  int i;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif

#pragma omp parallel private(idx)
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
              if(generic_polling_primary(count, Ndust))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= Ndust)
          break;

        drag_backreaction_kernel_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        drag_backreaction_kernel_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void drag_backreaction_kernel(void)
{
  long long ntot;

  sumup_large_ints(1, &Ndust, &ntot);
  if(ntot == 0)
    return;

  mpi_printf("DUST_LIVE: Performing drag backreaction on gas cells.\n");

  generic_set_MaxNexport();

  double t0 = second();

  generic_comm_pattern(Ndust, kernel_local, kernel_imported);

  double t1 = second();

  mpi_printf("DUST_LIVE: Done! Calculation took %g sec\n", timediff(t0, t1));
}

int drag_backreaction_kernel_evaluate(int target, int mode, int thread_id)
{
  int numnodes, *firstnode;
  data_in local, *in;
  data_out out;

  double h, h2, h3;
  double dx, dy, dz, r2;
  int k;

  double weight_fac, wk, hinv, hinv3;
#ifndef GFM_TOPHAT_KERNEL
  double u;
#endif

  MyDouble *pos;
  memset(&out, 0, sizeof(data_out));

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

  pos                = in->Pos;
  h                  = in->Hsml;
  MyFloat totngbmass = in->TotNgbMass;

  h2 = h * h;
  h3 = h * h * h;

  hinv = 1.0 / h;
#ifndef TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  for(int n = 0; n < nfound; n++)
    {
      int j = Thread[thread_id].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0) /* skip cells that have been swallowed or dissolved */
        {
          dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
          dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
          dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);

          r2 = dx * dx + dy * dy + dz * dz;

          if(r2 < h2)
            {
#ifndef GFM_TOPHAT_KERNEL
              u = sqrt(r2) * hinv;
              if(u < 0.5)
                wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
              else
                wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

              weight_fac = NORM_COEFF * P[j].Mass * wk * h3 / totngbmass;
#else
              wk         = hinv3 / NORM_COEFF;
              weight_fac = P[j].Mass / totngbmass;
#endif
              /* We'll recompute kinetic energy after updating momentum from drag. */
              double EK_old = 0.5 *
                              (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                               SphP[j].Momentum[2] * SphP[j].Momentum[2]) /
                              P[j].Mass;
              SphP[j].Energy -= EK_old;

              for(k = 0; k < 3; k++)
                {
                  SphP[j].Momentum[k] -= weight_fac * in->DragMomentum[k];
                }

              SphP[j].Energy += weight_fac * in->DragDustThermal;

              double EK_new = 0.5 *
                              (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                               SphP[j].Momentum[2] * SphP[j].Momentum[2]) /
                              P[j].Mass;
              SphP[j].Energy += EK_new;

              /* Negative change in gas kinetic energy from drag update also
               * goes into gas thermal energy.  Effectively, we could have left
               * SphP[j].Energy unchanged when computing kinetic energies, but
               * it's clearer this way. */
              SphP[j].Energy += (EK_old - EK_new);
            }
        }
    }

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif
#endif
