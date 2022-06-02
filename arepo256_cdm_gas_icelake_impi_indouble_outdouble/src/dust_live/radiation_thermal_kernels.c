/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/radiation_thermal_kernels.c
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

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef DUST_LIVE
#ifdef DL_THERMAL_IR

#ifndef MRT_IR
#error "Dust-radiation thermal coupling requires MRT_IR!"
#endif

typedef struct
{
  MyDouble Pos[3];
  MyDouble Mass;
  MyFloat Hsml;
  MyFloat TotNgbMass;
  MyFloat SigmaTot[MRT_BINS];
  MyFloat SigmaTotAbs[MRT_BINS];
  MyFloat SigmaTotGeo;
  MyFloat DtRP;
#ifdef DL_GRAIN_BINS
  MyFloat DustTemp[DL_GRAIN_BINS];
  MyFloat SigmaBinGeo[DL_GRAIN_BINS];
  MyFloat SigmaBinIRAbs[DL_GRAIN_BINS];
#else
  MyFloat DustTemp;
#endif
  MyIDType ID;
#ifdef DL_GRAIN_BINS
  MyFloat NumGrains[DL_GRAIN_BINS];
#endif
#ifdef DL_WINDS
  int IsWind;
#endif
  int Firstnode;
} data_in;

static data_in *DataGet;

typedef struct
{
#ifdef DL_GRAIN_BINS
  MyFloat DustTemp[DL_GRAIN_BINS];
#else
  MyFloat DustTemp;
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
  in->SigmaTotGeo = DustParticle[i].SigmaTotGeo;
  int p           = DustParticle[i].index;
  in->DtRP = (0.5 * (P[p].TimeBinHydro ? (((integertime)1) << P[p].TimeBinHydro) : 0) * All.Timebase_interval) / All.cf_hubble_a;
#ifdef DL_GRAIN_BINS
  for(int k = 0; k < DL_GRAIN_BINS; k++)
    {
      in->DustTemp[k]      = DTP(DustParticle[i].index).DustTemp[k];
      in->SigmaBinGeo[k]   = DustParticle[i].SigmaBinGeo[k];
      in->SigmaBinIRAbs[k] = DustParticle[i].SigmaBinIRAbs[k];
    }
#else
  in->DustTemp = DTP(DustParticle[i].index).DustTemp;
#endif
  in->ID = P[p].ID;
#ifdef DL_GRAIN_BINS for(int k = 0; k < DL_GRAIN_BINS; k++) { in->NumGrains[k] = DTP(DustParticle[i].index).NumGrains[k]; }
#endif
#ifdef DL_WINDS
  in->IsWind = DTP(DustParticle[i].index).IsWind;
#endif
  in->Firstnode = firstnode;
}

static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)
    {
#ifdef DL_GRAIN_BINS
      for(int k = 0; k < DL_GRAIN_BINS; k++)
        DustParticle[i].DustTemp[k] = out->DustTemp[k];
#else
      DustParticle[i].DustTemp = out->DustTemp;
#endif
    }
  else
    {
#ifdef DL_GRAIN_BINS
      for(int k = 0; k < DL_GRAIN_BINS; k++)
        DustParticle[i].DustTemp[k] += out->DustTemp[k];
#else
      DustParticle[i].DustTemp += out->DustTemp;
#endif
    }
}

struct dust_temp_params
{
  double kappa_rho;
  double E_IR;
  double T_g;
  double gas_prefac;
  double c_speed;
};

double dust_temp_func(double T_d, void *params)
{
  struct dust_temp_params *p = (struct dust_temp_params *)params;
  return -p->kappa_rho * (p->c_speed * radiation_constant * pow(T_d, 4) - c_internal_units * p->E_IR) + p->gas_prefac * (p->T_g - T_d);
}

double dust_temp_deriv(double T_d, void *params)
{
  struct dust_temp_params *p = (struct dust_temp_params *)params;
  return -p->kappa_rho * p->c_speed * radiation_constant * 4.0 * pow(T_d, 3) - p->gas_prefac;
}

void dust_temp_fdf(double T_d, void *params, double *y, double *dy)
{
  *y  = dust_temp_func(T_d, params);
  *dy = dust_temp_deriv(T_d, params);
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

        radiation_thermal_kernel_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        radiation_thermal_kernel_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void radiation_thermal_kernel(void)
{
  long long ntot;

  sumup_large_ints(1, &Ndust, &ntot);
  if(ntot == 0)
    return;

  mpi_printf("DUST_LIVE: Performing radiation thermal kernels on gas cells.\n");

  generic_set_MaxNexport();

  double t0 = second();

  generic_comm_pattern(Ndust, kernel_local, kernel_imported);

  double t1 = second();

  mpi_printf("DUST_LIVE: Done! Calculation took %g sec\n", timediff(t0, t1));
}

int radiation_thermal_kernel_evaluate(int target, int mode, int thread_id)
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
              // TODO: comoving factors
              double c_speed = 2.99792458e10 / All.UnitVelocity_in_cm_per_s;

              int b            = UV_BINS;
              double kappa_rho = weight_fac * in->SigmaTotAbs[b] / SphP[j].Volume; /* 1 / internal length */
              double flux_fac  = in->DtRP * c_internal_units * kappa_rho;
#ifdef DL_RADIATION_ABSORPTION
              /* Absorption of IR photons. */
              for(int k = 0; k < 3; k++)
                {
                  SphP[j].Cons_RT_F[b][k] *= exp(-flux_fac);
                  SphP[j].RT_F[b][k] = SphP[j].Cons_RT_F[b][k] / SphP[j].Volume;
                }
#endif

              double alpha_T = 0.5; /* dimensionless */

              double yHe = (1.0 - HYDROGEN_MASSFRAC) / (4.0 * HYDROGEN_MASSFRAC);
#ifdef COOLING
              double mu = (1.0 + 4.0 * yHe) / (1.0 + yHe + SphP[j].Ne);
#else
              double mu  = (1.0 + 4.0 * yHe) / (1.0 + yHe);
#endif
              double T_g = (mu * PROTONMASS) / BOLTZMANN * GAMMA_MINUS1 * SphP[j].Utherm * All.UnitPressure_in_cgs /
                           All.UnitDensity_in_cgs; /* K */
#ifdef USE_SFR
              if(SphP[j].Sfr > 0.0)
                T_g = 1.0e4;
#endif

              double n_H = (SphP[j].Density * HYDROGEN_MASSFRAC * All.cf_a3inv) /
                           (PROTONMASS / All.UnitMass_in_g * All.HubbleParam);  /* 1 / internal length^3 */
              double v_g = sqrt((8.0 * BOLTZMANN * T_g) / (M_PI * PROTONMASS)); /* cm/s */
              v_g /= All.UnitVelocity_in_cm_per_s;                              /* internal velocity */
              /* Distribute some weighted fraction of the dust particle's
               * cross section to this gas cell. */

#ifdef DL_GRAIN_BINS
              for(int g = 0; g < DL_GRAIN_BINS; g++)
                {
                  /* If following grain size bins, just use cross section in the bin. */
                  double n_sigma       = (weight_fac * in->SigmaBinGeo[g]) / SphP[j].Volume; /* 1 / internal length */
                  double kappa_rho_bin = weight_fac * in->SigmaBinIRAbs[g] / SphP[j].Volume; /* 1 / internal length */
                  /* Newton's method won't work if there is no dust in this bin. */
                  if(n_sigma == 0.0 && kappa_rho_bin == 0.0)
                    continue;
                  /* Use the old dust temperature as the starting guess. */
                  double T_d  = in->DustTemp[g];
                  double T_d0 = in->DustTemp[g];
#else
              {
                /* If no grain size bins, use total cross sections. */
                double n_sigma       = (weight_fac * in->SigmaTotGeo) / SphP[j].Volume;  /* 1 / internal length */
                double kappa_rho_bin = weight_fac * in->SigmaTotAbs[b] / SphP[j].Volume; /* 1 / internal length */
                /* Use the old dust temperature as the starting guess. */
                double T_d  = in->DustTemp;
                double T_d0 = in->DustTemp;
#endif

                  double gas_prefac = n_sigma * n_H * v_g * alpha_T * 2.0 *
                                      (BOLTZMANN / All.UnitEnergy_in_cgs * All.HubbleParam); /* internal energy/(volume*time)/K */

                  gsl_function_fdf FDF;
                  struct dust_temp_params params;
                  params.kappa_rho  = kappa_rho_bin;
                  params.E_IR       = SphP[j].DensPhot[b];
                  params.T_g        = T_g;
                  params.gas_prefac = gas_prefac;
                  params.c_speed    = c_speed;

                  FDF.f      = &dust_temp_func;
                  FDF.df     = &dust_temp_deriv;
                  FDF.fdf    = &dust_temp_fdf;
                  FDF.params = &params;

                  const gsl_root_fdfsolver_type *solver_type = gsl_root_fdfsolver_newton;
                  gsl_root_fdfsolver *solver                 = gsl_root_fdfsolver_alloc(solver_type);
                  gsl_root_fdfsolver_set(solver, &FDF, T_d);
                  int status               = GSL_CONTINUE;
                  const int MAX_ITERATIONS = 100;

                  for(int i = 0; i < MAX_ITERATIONS && status == GSL_CONTINUE; i++)
                    {
                      status = gsl_root_fdfsolver_iterate(solver);
                      if(status != GSL_SUCCESS)
                        break;

                      T_d0 = T_d;
                      T_d  = gsl_root_fdfsolver_root(solver);

                      /* Test based on relative error. */
                      status = gsl_root_test_delta(T_d, T_d0, 0, 0.001);
                    }

                  if(status != GSL_SUCCESS)
                    terminate("DUST_LIVE: failed during dust thermal kernel\n");

                  gsl_root_fdfsolver_free(solver);

#ifdef DL_GRAIN_BINS
                  out.DustTemp[g] += weight_fac * T_d;
#else
                out.DustTemp += weight_fac * T_d;
#endif

                  /* Collisional dust-gas exchange */
                  double psi_gd = gas_prefac * (T_d - T_g);           /* internal energy/volume/time */
                  double dE_g   = psi_gd * in->DtRP * SphP[j].Volume; /* internal energy */
                  SphP[j].Energy += dE_g;
                  // TODO: probably don't need line below
                  // SphP[j].Utherm += dE_g / P[j].Mass;

                  /* Dust-radiation exchange */
                  double psi_rd = kappa_rho_bin * (c_speed * radiation_constant * pow(T_d, 4) -
                                                   c_internal_units * SphP[j].DensPhot[b]); /* internal energy/volume/time */
                  double dE_IR  = psi_rd * in->DtRP * SphP[j].Volume;                       /* internal energy */
                  SphP[j].Cons_DensPhot[b] += dE_IR;

#ifdef DL_GRAIN_BINS
                }
#else
              }
#endif
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
