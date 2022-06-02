/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/MRT/RT_source_utils.c
 * \date        09/2017
 * \author      Federico Marinacci, Rahul Kannan, David Barnes, Mark Vogelsberger
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - 10.10.2017 Black hole routines added
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../allvars.h"
#include "../proto.h"
#include "RT.h"

#ifdef MRT_SOURCES

#define PHOT_NORMALIZATION 1e63

#ifdef MRT_STARS

#define LIFETIME 0.01
#define TAUNORM 0.003

double get_photon_released_star(int p, double dt, double lum_tot)
{
  if(P[p].Type != 4)
    terminate("Particle should be a star! Type=%d, index %d", P[p].Type, p);

  double imass = STP(p).InitialMass * All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS);

  double phot = lum_tot * imass * dt * 1e3 * SEC_PER_MEGAYEAR;
  if(phot < 0)
    {
      printf("Warning: injecting negative number of photons?\n");
      return 0;
    }
  return phot;
}

double interp_lum_tab_star(double age_of_star_in_Gyr, double metallicity, int ind)
{
  // I did not set any protection mechanism when N_age or N_metallicity equals 1.
  // fix this in the future, or mind your input. don't let it only contain 1 age or 1 metallicity.

  int indi, indj;
  double deli, delj;
  double lumbin;
  // interp in metallicity
  if(metallicity <= pow(10, LogMetallicity[0]))
    {
      indi = 0;
      deli = 0;
    }
  else
    {
      if(metallicity >= pow(10, LogMetallicity[N_metallicity - 1]))
        {
          indi = N_metallicity - 2;
          deli = 1;
        }
      else
        {
          for(int ind = 0; ind < N_metallicity - 1; ind++)
            {
              deli = (log10(metallicity) - LogMetallicity[ind]) / (LogMetallicity[ind + 1] - LogMetallicity[ind]);
              if(deli >= 0 && deli < 1)
                {
                  indi = ind;
                  break;
                }
            }
        }
    }
  // interp in age
  if(age_of_star_in_Gyr <= pow(10, LogAge[0]))
    {
      indj = 0;
      delj = 0;
    }
  else
    {
      if(age_of_star_in_Gyr >= pow(10, LogAge[N_age - 1]))
        {
          indj = N_age - 2;
          delj = 1;
        }
      else
        {
          for(int ind = 0; ind < N_age - 1; ind++)
            {
              delj = (log10(age_of_star_in_Gyr) - LogAge[ind]) / (LogAge[ind + 1] - LogAge[ind]);
              if(delj >= 0 && delj < 1)
                {
                  indj = ind;
                  break;
                }
            }
        }
    }

  int offset = (indi * N_age + indj) * UV_BINS;
  lumbin     = (1 - deli) * (1 - delj) * lum_tab[offset + ind] + (1 - deli) * delj * lum_tab[offset + UV_BINS + ind] +
           deli * (1 - delj) * lum_tab[offset + UV_BINS * N_age + ind] +
           deli * delj * lum_tab[offset + UV_BINS * N_age + UV_BINS + ind];
  return lumbin;
}

void start_stellar_sources(void)
{
  int idx, i;
  double time_begstep, dt, dtime, dtime_in_Gyr, age_of_star_in_Gyr_endstep, age_of_star_in_Gyr;
  double local_released_photons[UV_BINS], global_released_photons[UV_BINS];
  int global_sources;
  // TIMER_START(CPU_GFM_ENRICH);

  for(int bin = 0; bin < UV_BINS; bin++)
    local_released_photons[bin] = 0;

  //  mpi_printf_rt(0,"GFM_STELLAR_EVOLUTION: GFM_STELLAR_EVOLUTION...\n");

  StarParticle = (struct star_particle *)mymalloc("StarParticle", N_star * sizeof(struct star_particle));

  Nsource = 0;

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Ti_Current != All.Ti_Current)
        {
          terminate("how can this be?");
          drift_particle(i, All.Ti_Current);
        }

      if(P[i].Type == 4 && P[i].Mass > 0 && STP(i).BirthTime > 0)
        {
          if(All.ComovingIntegrationOn)
            time_begstep = All.TimeBegin * exp(All.Ti_begstep[P[i].TimeBinGrav] * All.Timebase_interval);
          else
            time_begstep = All.TimeBegin + All.Ti_begstep[P[i].TimeBinGrav] * All.Timebase_interval;

          dt = (P[i].TimeBinGrav ? (((integertime)1) << P[i].TimeBinGrav) : 0) * All.Timebase_interval;

          dtime = All.cf_atime * dt / All.cf_time_hubble_a;

          dtime *= All.UnitTime_in_s / All.HubbleParam;
          dtime_in_Gyr = 1e-3 * dtime / SEC_PER_MEGAYEAR;

          age_of_star_in_Gyr_endstep = get_time_difference_in_Gyr(
              STP(i).BirthTime, time_begstep); /* (Note: All.Ti_begstep[] has already been advanced for the next step at this point) */
          age_of_star_in_Gyr = age_of_star_in_Gyr_endstep - dtime_in_Gyr;

          if(N_age == 1)  // no age bin, stop emitting photons after LIFETIME
            {
              if(age_of_star_in_Gyr > LIFETIME)
                continue;
              else
                dt = ((age_of_star_in_Gyr_endstep < LIFETIME) ? dtime_in_Gyr : (LIFETIME - age_of_star_in_Gyr));
            }
          else  // will interpolate the number of photons to emit in age
            dt = dtime_in_Gyr;

          StarParticle[Nsource].index   = i;
          StarParticle[Nsource].NumNgb  = 0;
          StarParticle[Nsource].NormSph = 0;

          double lumbin[UV_BINS], lum_tot = 0;
          double metallicity;
          double phot_released_star;
          if(N_age == 1 && N_metallicity == 1)
            {
              for(int ind = 0; ind < UV_BINS; ind++)
                {
                  lumbin[ind] = lum[ind];
                  lum_tot += lumbin[ind];
                }
              /* photon rate is computed for a SFR of 1 Msun/yr -> convert to appropriate units */
              phot_released_star = All.EscapeFraction * get_photon_released_star(i, dt, PhotRate) / PHOT_NORMALIZATION / 1e9 /
                                   TAUNORM;  // Xiaohan: don't divide by 1e9 TAUNORM if use bc03
            }
          else
            {
              metallicity = STP(i).Metallicity;
              for(int ind = 0; ind < UV_BINS; ind++)
                {
                  lumbin[ind] = interp_lum_tab_star(age_of_star_in_Gyr, metallicity, ind);
                  lum_tot += lumbin[ind];
                }
              /* photon rate is computed for a SFR of 1 Msun/yr -> convert to appropriate units */
              phot_released_star =
                  All.EscapeFraction * get_photon_released_star(i, dt, lum_tot) / PHOT_NORMALIZATION;  // use BC03 here
            }

          for(int bin = 0; bin < UV_BINS; bin++)
            {
              StarParticle[Nsource].TotalPhotReleased[bin] = lumbin[bin] / lum_tot * phot_released_star;
              local_released_photons[bin] += StarParticle[Nsource].TotalPhotReleased[bin];
              if(local_released_photons[bin] < 0)
                printf("Warning: injecting a negative number of photons!\n");
            }

          // printf("source %d: i=%d, LogAge=%g, Z=%g, photons=%g\n", Nsource, i, log10(age_of_star_in_Gyr), log10(metallicity),
          // StarParticle[Nsource].TotalPhotReleased[0]);
          Nsource++;
        }
    }

  MPI_Reduce(&local_released_photons, &global_released_photons, UV_BINS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Allreduce(&Nsource, &global_sources, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(All.MRT_On == 0 && global_sources > 0)
    {
      All.MRT_On = 1;
      int i, j, num1;
      for(i = 0; i < NumGas; i++)
        {
          for(num1 = 0; num1 < MRT_BINS; num1++)
            {
              SphP[i].DensPhot[num1]      = MINDENSPHOT;
              SphP[i].Cons_DensPhot[num1] = SphP[i].DensPhot[num1] * SphP[i].Volume;
              if(isnan(SphP[i].Cons_DensPhot[num1]))
                terminate("Something went wrong\n");

              for(j = 0; j < 3; j++)
                {
                  SphP[i].RT_F[num1][j]      = MINDENSPHOT;
                  SphP[i].Cons_RT_F[num1][j] = SphP[i].RT_F[num1][j] * SphP[i].Volume;
                }
            }
        }
    }

  mpi_printf_rt(0, "MRT_STARS: %d stellar radiation sources.\n", global_sources);
  for(int bin = 0; bin < UV_BINS; bin++)
    mpi_printf_rt(0, "MRT_STARS: in total %g photons injected in frequency bin %d\n", global_released_photons[bin], bin);
}

void end_stellar_sources(void)
{
  myfree(StarParticle);

  mpi_printf_rt(0, "RT: stellar sources done.\n");

  // TIMER_STOP(CPU_GFM_ENRICH);
}

#endif /* MRT_STARS */

#ifdef MRT_BH

struct bh_particle *BHParticle;

double get_photon_released_blackholes(int p, double dt)
{
  if(P[p].Type != 5)
    terminate("Particle should be a black hole! Type=%d, index %d", P[p].Type, p);

  double phot = pow(10.0, All.LogAGNLuminosity) * dt * All.UnitTime_in_s / All.UnitEnergy_in_cgs;

  return phot;
}

void start_blackhole_sources(void)
{
  int idx, i;
  double dt, dtime, dtime_in_Gyr;
  double local_released_photons[MRT_BINS], global_released_photons[MRT_BINS];
  int global_sources;

  for(int bin = 0; bin < MRT_BINS; bin++)
    local_released_photons[bin] = 0;

  BHParticle = (struct bh_particle *)mymalloc("BHParticle", NumBHs * sizeof(struct bh_particle));

  Nsource = 0;

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Ti_Current != All.Ti_Current)
        {
          terminate("how can this be?");
          drift_particle(i, All.Ti_Current);
        }

      if(P[i].Type == 5 && P[i].Mass > 0)
        {
          dt = (P[i].TimeBinGrav ? (((integertime)1) << P[i].TimeBinGrav) : 0) * All.Timebase_interval;

          if(All.ComovingIntegrationOn)
            dtime = All.Time * dt / All.cf_time_hubble_a;
          else
            dtime = dt;

          dtime *= All.UnitTime_in_s / All.HubbleParam;
          dtime_in_Gyr = 1e-3 * dtime / SEC_PER_MEGAYEAR;

          BHParticle[Nsource].index  = i;
          BHParticle[Nsource].NumNgb = 0;
#ifndef MRT_BH_OMEGA_WEIGHT
          BHParticle[Nsource].NormSph = 0;
#else
          BHParticle[Nsource].TotalSolidAngle = 0;
#endif
          BHParticle[Nsource].Dhsmlrho = BPP(i).BH_PhotonHsml;
          for(int bin = 0; bin < MRT_BINS; bin++)
            {
              BHParticle[Nsource].TotalPhotReleased[bin] = 0.0;
            }

#ifdef MRT_BH_UV_INJECTION
          for(int bin = 0; bin < UV_BINS; bin++)
            {
              BHParticle[Nsource].TotalPhotReleased[bin] = lum[bin] * get_photon_released_blackholes(i, dt) * All.UnitEnergy_in_cgs /
                                                           (MeanPhotonEnergy[bin] * ELECTRONVOLT_IN_ERGS) / PHOT_NORMALIZATION;
              local_released_photons[bin] = BHParticle[Nsource].TotalPhotReleased[bin];
            }
#endif

#ifdef MRT_BH_IR_INJECTION
          for(int bin = UV_BINS; bin < MRT_BINS; bin++)
            {
              BHParticle[Nsource].TotalPhotReleased[bin] = lum[bin] * get_photon_released_blackholes(i, dt);
              local_released_photons[bin]                = BHParticle[Nsource].TotalPhotReleased[bin];
            }
#endif
          Nsource++;
        }
    }
  MPI_Reduce(&local_released_photons, &global_released_photons, MRT_BINS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&Nsource, &global_sources, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  for(int bin = 0; bin < MRT_BINS; bin++)
    All.AGNPhotonsInjected[bin] += global_released_photons[bin];
  mpi_printf_rt(0, "MRT_BH: %d black hole radiation sources.\n", global_sources);
  for(int bin = 0; bin < MRT_BINS; bin++)
    mpi_printf_rt(0, "MRT_BH: in step/total %g / %g photons injected in frequency bin %d\n", global_released_photons[bin],
                  All.AGNPhotonsInjected[bin], bin);
}

void end_blackhole_sources(void)
{
  myfree(BHParticle);

  mpi_printf_rt(0, "RT: BH sources done.\n");
}

#endif /* MRT_BH */

#endif /* MRT_SOURCES */
