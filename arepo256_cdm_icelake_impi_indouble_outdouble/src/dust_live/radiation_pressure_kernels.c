/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/radiation_pressure_kernels.c
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
#if defined(DL_RADIATION_PRESSURE) || defined(DL_RADIATION_ABSORPTION) || defined(DL_OUTPUT_RT_FLUX)

typedef struct
{
  MyDouble Pos[3];
  MyDouble Mass;
  MyFloat Hsml;
  MyFloat TotNgbMass;
  MyFloat SigmaTot[MRT_BINS];
  MyFloat SigmaTotAbs[MRT_BINS];
  MyFloat DtRP;
#ifdef DL_WINDS
  int IsWind;
#endif
  int Firstnode;
} data_in;

static data_in *DataGet;

typedef struct
{
#ifdef DL_OUTPUT_RT_FLUX
  MyFloat LocalRT_F[MRT_BINS][3];
#endif
#ifdef DL_RADIATION_PRESSURE
  MyFloat DeltaRPMomentum[3];
#endif
} data_out;

static data_out *DataResult;

static void particle2in(data_in *in, int i, int firstnode)
{
  for(int k = 0; k < 3; k++)
    {
      in->Pos[k] = P[DustParticle[i].index].Pos[k];
    }
  in->Mass       = P[DustParticle[i].index].Mass;
  in->Hsml       = DTP(DustParticle[i].index).Hsml;
  in->TotNgbMass = DustParticle[i].TotNgbMass;
  for(int b = 0; b < MRT_BINS; b++)
    {
      in->SigmaTot[b]    = DustParticle[i].SigmaTot[b];
      in->SigmaTotAbs[b] = DustParticle[i].SigmaTotAbs[b];
    }
  int p    = DustParticle[i].index;
  in->DtRP = (0.5 * (P[p].TimeBinHydro ? (((integertime)1) << P[p].TimeBinHydro) : 0) * All.Timebase_interval) / All.cf_hubble_a;
#ifdef DL_WINDS
  in->IsWind = DTP(DustParticle[i].index).IsWind;
#endif
  in->Firstnode = firstnode;
}

static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)
    {
#ifdef DL_OUTPUT_RT_FLUX
      for(int b = 0; b < MRT_BINS; b++)
        {
          for(int k = 0; k < 3; k++)
            {
              DustParticle[i].LocalRT_F[b][k] = out->LocalRT_F[b][k];
            }
        }
#endif
#ifdef DL_RADIATION_PRESSURE
      for(int k = 0; k < 3; k++)
        {
          DustParticle[i].DeltaRPMomentum[k] = out->DeltaRPMomentum[k];
        }
#endif
    }
  else
    {
#ifdef DL_OUTPUT_RT_FLUX
      for(int b = 0; b < MRT_BINS; b++)
        {
          for(int k = 0; k < 3; k++)
            {
              DustParticle[i].LocalRT_F[b][k] += out->LocalRT_F[b][k];
            }
        }
#endif
#ifdef DL_RADIATION_PRESSURE
      for(int k = 0; k < 3; k++)
        {
          DustParticle[i].DeltaRPMomentum[k] += out->DeltaRPMomentum[k];
        }
#endif
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

        radiation_pressure_kernel_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        radiation_pressure_kernel_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void radiation_pressure_kernel(void)
{
  long long ntot;

  sumup_large_ints(1, &Ndust, &ntot);
  if(ntot == 0)
    return;

  mpi_printf("DUST_LIVE: Performing radiation pressure kernels on gas cells.\n");

  generic_set_MaxNexport();

  double t0 = second();

  generic_comm_pattern(Ndust, kernel_local, kernel_imported);

  double t1 = second();

  mpi_printf("DUST_LIVE: Done! Calculation took %g sec\n", timediff(t0, t1));
}

int radiation_pressure_kernel_evaluate(int target, int mode, int thread_id)
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

              if(in->Mass == 0.0)
                continue;

#ifdef DL_WINDS
              if(in->IsWind)
                continue;
#endif
              double c_speed = 2.99792458e10 / All.UnitVelocity_in_cm_per_s;

              for(int b = 0; b < MRT_BINS; b++)
                {
                  double sigma_tot_pr  = in->SigmaTot[b];
                  double sigma_tot_abs = in->SigmaTotAbs[b];
                  double flux_fac      = in->DtRP * c_internal_units * weight_fac * sigma_tot_abs / SphP[j].Volume;

#ifdef DL_RADIATION_ABSORPTION
                  /* IR bin absorption handled later in thermal coupling */
                  if(b < UV_BINS)
                    {
                      SphP[j].Cons_DensPhot[b] *= exp(-flux_fac);
                      SphP[j].DensPhot[b] = SphP[j].Cons_DensPhot[b] / SphP[j].Volume;
                    }
#endif

#if defined(MRT_IR) && defined(DL_RADIATION_ABSORPTION)
                  /* UV photons move to IR bin */
                  if(b < UV_BINS)
                    {
                      /* UV units are 10**63 photons, while we want the IR
                       * units to be internal energy units. */
                      SphP[j].Cons_DensPhot[UV_BINS] += (flux_fac * SphP[j].Cons_DensPhot[b]) * DRT.photons2energy[b];
                      SphP[j].DensPhot[UV_BINS] = SphP[j].Cons_DensPhot[UV_BINS] / SphP[j].Volume;
                    }
#endif

                  for(int k = 0; k < 3; k++)
                    {
                      // TODO: tie to RTNumSubCycles?  comoving factors?  timestep constraint?
#ifdef DL_RADIATION_PRESSURE
                      double dp_k = in->DtRP * (sigma_tot_pr * weight_fac * SphP[j].RT_F[b][k] * DRT.photons2energy[b] / c_speed);
                      out.DeltaRPMomentum[k] += dp_k;
#endif

#ifdef DL_OUTPUT_RT_FLUX
                      out.LocalRT_F[b][k] += weight_fac * SphP[j].RT_F[b][k];
#endif

#ifdef DL_RADIATION_ABSORPTION
                      /* IR bin absorption handled later in thermal coupling */
                      if(b < UV_BINS)
                        {
                          SphP[j].Cons_RT_F[b][k] *= exp(-flux_fac);
                          SphP[j].RT_F[b][k] = SphP[j].Cons_RT_F[b][k] / SphP[j].Volume;
                        }
#endif
                    }
                }
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
