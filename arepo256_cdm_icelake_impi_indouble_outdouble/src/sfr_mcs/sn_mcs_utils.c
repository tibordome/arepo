/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sfr_mcs/sn_mcs.c
 * \date        04/2018
 * \author      Matthew C Smith
 * \brief
 * \details     Originally developed in 2015, ported into main repo 2018.
                Please contact the author before use.
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

#include "../allvars.h"
#include "../proto.h"

#if defined(SFR_MCS) && defined(SN_MCS)

#include <gsl/gsl_randist.h>

#ifdef SN_MCS_NO_ENERGY
#define SN_ENERGY 0.0
#elif defined(SN_MCS_CHANCE)
#define SN_ENERGY (1e51 * SN_MCS_CHANCE)
#else
#define SN_ENERGY 1e51 /* ergs */
#endif

#ifdef SN_MCS_MECHANICAL
#define P_TERMINAL 3e5 /* Msun km/s */
#endif

void init_sne(void)
{
  All.SupernovaEnergy = SN_ENERGY * All.HubbleParam / All.UnitEnergy_in_cgs;
#ifdef SN_MCS_MECHANICAL
  All.SupernovaTerminalMomentum = P_TERMINAL * SOLAR_MASS * 1e5 * All.HubbleParam / (All.UnitMass_in_g * All.UnitVelocity_in_cm_per_s);
#endif

#ifdef IMF_SAMPLING_MCS
  All.SNStarMinMass *= SOLAR_MASS / (All.UnitMass_in_g / All.HubbleParam);
  All.SNStarMaxMass *= SOLAR_MASS / (All.UnitMass_in_g / All.HubbleParam);
#endif
}

#ifdef SN_MCS_LOG
void setup_sn_log(void)
{
  if(RestartFlag > 2)
    return;

  sn_dens_hist = gsl_histogram_alloc(SN_MCS_LOG_N);

  if(RestartFlag != 1)
    {
      gsl_histogram_set_ranges_uniform(sn_dens_hist, SN_MCS_LOG_MIN, SN_MCS_LOG_MAX);
      /* Write bin edges to file */
      if(ThisTask == 0)
        {
          for(size_t i = 0; i <= sn_dens_hist->n; i++)
            fprintf(FdSNdens, "%g ", sn_dens_hist->range[i]);

          fprintf(FdSNdens, "\n");
          fflush(FdSNdens);
        }
    }
}

void sn_add_to_log(int i)
{
  double n;

  n = SphP[i].Density * All.cf_a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
  n *= HYDROGEN_MASSFRAC / PROTONMASS;

  if(gsl_histogram_increment(sn_dens_hist, log10(n)))
    {
      printf("SN_MCS_LOG: OUT OF BOUNDS on Task: %d Time: %g log10(n): %g\n", ThisTask, All.Time, log10(n));
      /* out of bounds, add to top or bottom bin */
      if(log10(n) >= gsl_histogram_max(sn_dens_hist))
        sn_dens_hist->bin[sn_dens_hist->n - 1] += 1;
      else if(log10(n) >= gsl_histogram_min(sn_dens_hist))
        sn_dens_hist->bin[0] += 0;
    }
}

void write_sn_dens_log(void)
{
  double *denshistogram;

  mpi_printf("SN_MCS_LOG: Writing density histogram...\n");

  if(ThisTask == 0)
    denshistogram = (double *)mymalloc("denshistogram", sn_dens_hist->n * sizeof(double));
  else
    denshistogram = NULL;

  /* Sum histograms from all tasks and reset */
  MPI_Reduce(sn_dens_hist->bin, denshistogram, sn_dens_hist->n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  gsl_histogram_reset(sn_dens_hist);

  if(ThisTask == 0)
    {
      fprintf(FdSNdens, "%e", All.Time);
      for(size_t i = 0; i < sn_dens_hist->n; i++)
        fprintf(FdSNdens, " %g", denshistogram[i]);

      fprintf(FdSNdens, "\n");
      fflush(FdSNdens);

      myfree(denshistogram);
    }
}
#endif  // SN_MCS_LOG

#endif
