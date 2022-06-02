/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sfr_mcs/sfr_mcs.c
 * \date        04/2018
 * \author      Matthew C Smith
 * \brief
 * \details     Originally developed in 2015, ported into main repo 2018.
                Please contact the author before use.
 *
 *
 * \par Major modifications and contributions:
 *
 * - 17.01.2022 Re-ported (code had diverged) to reflect new physics described in Smith+2021
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef SFR_MCS

#ifndef USE_SFR
#error "SFR_MCS requires USE_SFR"
#endif

#if defined(JEANS_PRESSURE_LIMIT_MCS) && (defined(JEANS_PRESSURE_LIMIT) || defined(ENFORCE_JEANS_STABILITY_OF_CELLS))
#error \
    "You have defined JEANS_PRESSURE_LIMIT_MCS along with some other Jeans criteria based EOS,\nwhich you probably don't want to do."
#endif

#ifdef GFM
#error "SFR_MCS cannot be compiled with GFM"
#endif

#if defined(SFR_MCS_LOG) && !(defined(SFR_MCS_LOG_N) && defined(SFR_MCS_LOG_MIN) && defined(SFR_MCS_LOG_MAX))
#error "If using SFR_MCS_LOG, you must define SFR_MCS_LOG_N, SFR_MCS_LOG_MIN and SFR_MCS_LOG_MAX"
#endif

#if defined(SN_MCS_LOG) && !(defined(SN_MCS_LOG_N) && defined(SN_MCS_LOG_MIN) && defined(SN_MCS_LOG_MAX))
#error "If using SN_MCS_LOG, you must define SN_MCS_LOG_N, SN_MCS_LOG_MIN and SN_MCS_LOG_MAX"
#endif

void init_star_formation(void)
{
  sfr_criteria_announce();

#ifdef IMF_SAMPLING_MCS
  init_imf();
#endif

  if(RestartFlag != RESTART_RESTART)
    {
#if SFR_MCS_SELECT_CRITERIA == 1
      /* Used (for SMAUG) to only allow SF when local Jeans length < All.SfrCritLength */
      All.SfrCritFactor = GAMMA * M_PI / (All.G * All.SfrCritLength * All.SfrCritLength);
#elif(SFR_MCS_SELECT_CRITERIA == 2) || (SFR_MCS_SELECT_CRITERIA == 3)
      /* Used for only SF when local Jeans mass < All.SfrCritJeansMassN * mcell */
      All.SfrCritFactor = pow(GAMMA, 1.5) * pow(M_PI, 2.5) / (6.0 * pow(All.G, 1.5) * All.SfrCritJeansMassN);
#endif
#if(SFR_MCS_SELECT_CRITERIA == 0) || (SFR_MCS_SELECT_CRITERIA == 3)
      All.DensThreshold *= ((PROTONMASS / HYDROGEN_MASSFRAC) / (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam));
#endif
      if(All.ComovingIntegrationOn)
        All.OverDensThresh = All.CritOverDensity * All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

#ifdef SMAUG_PRESSURE_FLOOR
      All.PolytropeFactor = BOLTZMANN * All.Polytrope_Tstar / (All.UnitEnergy_in_cgs);
      All.PolytropeFactor *= pow((PROTONMASS / All.UnitMass_in_g), (-All.Polytrope_gstar));
      All.PolytropeFactor *=
          pow((All.Polytrope_nstar * All.UnitLength_in_cm * All.UnitLength_in_cm * All.UnitLength_in_cm), (1.0 - All.Polytrope_gstar));
#endif
#if defined(SN_MCS) || defined(HII_MCS) || defined(PE_MCS)
      init_stellar_feedback();
#endif
    }
}

#ifdef SFR_MCS_LOG
void setup_sf_log(void)
{
  if(RestartFlag > 2)
    return;

  sf_dens_hist = gsl_histogram_alloc(SFR_MCS_LOG_N);

  if(RestartFlag != RESTART_RESTART)
    {
      gsl_histogram_set_ranges_uniform(sf_dens_hist, SFR_MCS_LOG_MIN, SFR_MCS_LOG_MAX);
      /* Write bin edges to file */
      if(ThisTask == 0)
        {
          for(size_t i = 0; i <= sf_dens_hist->n; i++)
            fprintf(FdSFdens, "%g ", sf_dens_hist->range[i]);

          fprintf(FdSFdens, "\n");
          fflush(FdSFdens);
        }
    }
}

void sf_add_to_log(double rho)
{
  double n;
  n = rho * All.cf_a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
  n *= HYDROGEN_MASSFRAC / PROTONMASS;
  if(gsl_histogram_increment(sf_dens_hist, log10(n)))
    {
      printf("SFR_MCS_LOG: OUT OF BOUNDS on Task: %d Time: %g log10(n): %g\n", ThisTask, All.Time, log10(n));
      /* out of bounds, add to top or bottom bin */
      if(log10(n) >= gsl_histogram_max(sf_dens_hist))
        sf_dens_hist->bin[sf_dens_hist->n - 1] += 1;
      else if(log10(n) >= gsl_histogram_min(sf_dens_hist))
        sf_dens_hist->bin[0] += 0;
    }
}

void write_sf_dens_log(void)
{
  double *denshistogram;

  mpi_printf("SFR_MCS_LOG: Writing density histogram...\n");

  if(ThisTask == 0)
    denshistogram = (double *)mymalloc("denshistogram", sf_dens_hist->n * sizeof(double));
  else
    denshistogram = NULL;

  /* Sum histograms from all tasks and reset */
  MPI_Reduce(sf_dens_hist->bin, denshistogram, sf_dens_hist->n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  gsl_histogram_reset(sf_dens_hist);

  if(ThisTask == 0)
    {
      fprintf(FdSFdens, "%e", All.Time);
      for(size_t i = 0; i < sf_dens_hist->n; i++)
        fprintf(FdSFdens, " %g", denshistogram[i]);

      fprintf(FdSFdens, "\n");
      fflush(FdSFdens);

      myfree(denshistogram);
    }
}
#endif

void cooling_and_starformation(void)
{
  TIMER_START(CPU_COOLINGSFR);

#ifdef GRACKLE
  mpi_printf("GRACKLE: Cooling active cells\n");
  cool_active_cells();
#endif

  /* clear the SFR stored in the active timebins */
  for(int bin = 0; bin < TIMEBINS; bin++)
    if(TimeBinSynchronized[bin])
      TimeBinSfr[bin] = 0;

  int idx, i;
  double dens, t_ff;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Mass == 0 && P[i].ID == 0)
        continue; /* skip cells that have been swallowed or eliminated */

#ifndef GRACKLE
      cool_cell(i);
#endif

      dens = SphP[i].Density;

      SphP[i].Sfr = 0.0;

#if SFR_MCS_SELECT_CRITERIA == 0
      if(dens * All.cf_a3inv < All.DensThreshold)
        continue;
#elif SFR_MCS_SELECT_CRITERIA == 1
      /* Used (for SMAUG) to only allow SF when local Jeans length < All.SfrCritLength */
      if((dens * dens * All.cf_a3inv / SphP[i].Pressure) < All.SfrCritFactor)
        continue;

#elif SFR_MCS_SELECT_CRITERIA == 2
      if((P[i].Mass * dens * dens * sqrt(All.cf_a3inv) / pow(SphP[i].Pressure, 1.5)) < All.SfrCritFactor)
        continue;
#elif SFR_MCS_SELECT_CRITERIA == 3
      if((dens * All.cf_a3inv < All.DensThreshold) ||
         ((P[i].Mass * dens * dens * sqrt(All.cf_a3inv) / pow(SphP[i].Pressure, 1.5)) < All.SfrCritFactor))
        continue;
#endif

      if(All.ComovingIntegrationOn) /* to protect against SF at too high redshift */
        if(dens < All.OverDensThresh)
          continue;

#ifdef HII_MCS
      /*No SF in ionized region*/
      if((SphP[i].StromgrenSourceID > 0) && (SphP[i].StromgrenSourceID != HII_MCS_IGNORE_FLAG))
        continue;
#endif

#ifdef TURB_APPROX_MCS
      if(SphP[i].TurbSpecEnergy < All.MinTurbSpecEnergy)
        {
          SphP[i].TurbSpecEnergy = All.MinTurbSpecEnergy;
          SphP[i].TurbEnergy     = P[i].Mass * SphP[i].TurbSpecEnergy;
        }
#endif

      t_ff = sqrt(3.0 * M_PI / (32.0 * All.G * dens * All.cf_a3inv));

#if SFR_MCS_RATE_CRITERIA == 0
      SphP[i].Sfr = All.SfrEfficiency * P[i].Mass / t_ff;

#elif SFR_MCS_RATE_CRITERIA == 1
      /*Padoan 2012. eps_ff = eps_w * exp(-1.6 t_ff / t_dyn).
      t_dyn = L/(2*vdisp3D) = R/vdisp3d = sqrt(5/gradv_sq) */
#ifdef TURB_APPROX_MCS
      MyFloat t_dyn  = get_cell_radius(i) / sqrt(2.0 * SphP[i].TurbSpecEnergy);
      MyFloat eps_ff = All.SfrEfficiency * exp(-1.6 * t_ff / t_dyn);
#else
      MyFloat eps_ff = All.SfrEfficiency * exp(-0.716 * t_ff * sqrt(SphP[i].gradv_sq) * All.cf_a2inv);
#endif
      SphP[i].Sfr    = eps_ff * P[i].Mass / t_ff;

#elif SFR_MCS_RATE_CRITERIA == 2
      /*Padoan 2012. eps_ff = eps_w * exp(-1.6 t_ff / t_dyn).
      t_dyn = L/(2*sqrt(vdisp3D^2 + c_s^2) */
      MyFloat r          = get_cell_radius(i);
#ifdef TURB_APPROX_MCS
      MyFloat vdisp3D_sq = 2.0 * SphP[i].TurbSpecEnergy;
#else
      MyFloat vdisp3D_sq = SphP[i].gradv_sq * r * r * All.cf_a2inv / 5.0;
#endif
      MyFloat c_s        = get_sound_speed(i);
      MyFloat t_dyn      = r * All.cf_atime / sqrt(vdisp3D_sq + 3.0 * c_s * c_s);
      MyFloat eps_ff     = All.SfrEfficiency * exp(-1.6 * t_ff / t_dyn);
      SphP[i].Sfr        = eps_ff * P[i].Mass / t_ff;

#elif SFR_MCS_RATE_CRITERIA == 3
      MyFloat r      = get_cell_radius(i);
#ifdef TURB_APPROX_MCS
      MyFloat virial = 10.0 * SphP[i].TurbSpecEnergy * r / (3.0 * All.G * P[i].Mass);
#else
      MyFloat virial = SphP[i].gradv_sq * r * r * r / (3.0 * All.G * P[i].Mass * All.cf_atime);
#endif
      if(virial > 1)
        continue;
      SphP[i].Sfr = All.SfrEfficiency * P[i].Mass / t_ff;

#elif SFR_MCS_RATE_CRITERIA == 4
      MyFloat r      = get_cell_radius(i);
      MyFloat c_s    = get_sound_speed(i);
#ifdef TURB_APPROX_MCS
      MyFloat vdisp2 = 2.0 * SphP[i].TurbSpecEnergy / 3.0 + c_s * c_s;
#else
      MyFloat vdisp2 = SphP[i].gradv_sq * r * r * All.cf_a2inv / 15.0 + c_s * c_s;
#endif
      MyFloat virial = 5.0 * vdisp2 * r * All.cf_atime / (All.G * P[i].Mass);
      if(virial > 1)
        continue;
      SphP[i].Sfr = All.SfrEfficiency * P[i].Mass / t_ff;
#endif

      SphP[i].Sfr *= (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);
      TimeBinSfr[P[i].TimeBinHydro] += SphP[i].Sfr;
    } /* End of loop over active particles */

  TIMER_STOP(CPU_COOLINGSFR);
}

double get_starformation_rate(int i)
{
  double t_ff, dens;
  double rateOfSF = 0.0;

#ifdef HII_MCS
  if((SphP[i].StromgrenSourceID > 0) && (SphP[i].StromgrenSourceID != HII_MCS_IGNORE_FLAG))
    return 0;
#endif

  dens = SphP[i].Density;

#if SFR_MCS_SELECT_CRITERIA == 0
  if(dens * All.cf_a3inv < All.DensThreshold)
    return 0.0;
#elif SFR_MCS_SELECT_CRITERIA == 1
  /* Used (for SMAUG) to only allow SF when local Jeans length < All.SfrCritLength */
  if((dens * dens * All.cf_a3inv / SphP[i].Pressure) < All.SfrCritFactor)
    return 0.0;
#elif SFR_MCS_SELECT_CRITERIA == 2
  if((P[i].Mass * dens * dens * sqrt(All.cf_a3inv) / pow(SphP[i].Pressure, 1.5)) < All.SfrCritFactor)
    return 0.0;
#elif SFR_MCS_SELECT_CRITERIA == 3
  if((dens * All.cf_a3inv < All.DensThreshold) ||
     ((P[i].Mass * dens * dens * sqrt(All.cf_a3inv) / pow(SphP[i].Pressure, 1.5)) < All.SfrCritFactor))
    return 0.0;
#endif

  if(All.ComovingIntegrationOn) /* to protect against SF at too high redshift */
    if(dens < All.OverDensThresh)
      return 0.0;

#ifdef TURB_APPROX_MCS
  if(SphP[i].TurbSpecEnergy < All.MinTurbSpecEnergy)
    {
      SphP[i].TurbSpecEnergy = All.MinTurbSpecEnergy;
      SphP[i].TurbEnergy     = P[i].Mass * SphP[i].TurbSpecEnergy;
    }
#endif

  t_ff = sqrt(3.0 * M_PI / (32.0 * All.G * dens * All.cf_a3inv));

#if SFR_MCS_RATE_CRITERIA == 0
  rateOfSF = All.SfrEfficiency * P[i].Mass / t_ff;

#elif SFR_MCS_RATE_CRITERIA == 1
  /*Padoan 2012. eps_ff = eps_w * exp(-1.6 t_ff / t_dyn).
  t_dyn = L/(2*vdisp3D) = R/vdisp3d = sqrt(5/gradv_sq) */
#ifdef TURB_APPROX_MCS
  MyFloat t_dyn = get_cell_radius(i) / sqrt(2.0 * SphP[i].TurbSpecEnergy);
  MyFloat eps_ff = All.SfrEfficiency * exp(-1.6 * t_ff / t_dyn);
#else
  MyFloat eps_ff = All.SfrEfficiency * exp(-0.716 * t_ff * sqrt(SphP[i].gradv_sq) * All.cf_a2inv);
#endif
  rateOfSF = eps_ff * P[i].Mass / t_ff;

#elif SFR_MCS_RATE_CRITERIA == 2
  /*Padoan 2012. eps_ff = eps_w * exp(-1.6 t_ff / t_dyn).
  t_dyn = L/(2*sqrt(vdisp3D^2 + c_s^2) */
  MyFloat r = get_cell_radius(i);
#ifdef TURB_APPROX_MCS
  MyFloat vdisp3D_sq = 2.0 * SphP[i].TurbSpecEnergy;
#else
  MyFloat vdisp3D_sq = SphP[i].gradv_sq * r * r * All.cf_a2inv / 5.0;
#endif
  MyFloat c_s = get_sound_speed(i);
  MyFloat t_dyn = r * All.cf_atime / sqrt(vdisp3D_sq + c_s * c_s);
  MyFloat eps_ff = All.SfrEfficiency * exp(-1.6 * t_ff / t_dyn);
  rateOfSF = eps_ff * P[i].Mass / t_ff;

#elif SFR_MCS_RATE_CRITERIA == 3
  MyFloat r = get_cell_radius(i);
#ifdef TURB_APPROX_MCS
  MyFloat virial = 10.0 * SphP[i].TurbSpecEnergy * r / (3.0 * All.G * P[i].Mass);
#else
  MyFloat virial = SphP[i].gradv_sq * r * r * r / (3.0 * All.G * P[i].Mass * All.cf_atime);
#endif
  if(virial > 1)
    return 0.0;
  rateOfSF = All.SfrEfficiency * P[i].Mass / t_ff;

#elif SFR_MCS_RATE_CRITERIA == 4
  MyFloat r = get_cell_radius(i);
  MyFloat c_s = get_sound_speed(i);
#ifdef TURB_APPROX_MCS
  MyFloat vdisp2 = 2.0 * SphP[i].TurbSpecEnergy / 3.0 + c_s * c_s;
#else
  MyFloat vdisp2 = SphP[i].gradv_sq * r * r * All.cf_a2inv / 15.0 + c_s * c_s;
#endif
  MyFloat virial = 5.0 * vdisp2 * r * All.cf_atime / (All.G * P[i].Mass);
  if(virial > 1)
    return 0.0;
  rateOfSF = All.SfrEfficiency * P[i].Mass / t_ff;
#endif

  rateOfSF *= (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

  return rateOfSF;
}

#ifndef SN_MCS_INITIAL_DRIVING
#ifdef SFR_MCS_CHECKS
void check_AuxDataID_references_mcs(void)
{
  mpi_printf("SFR_MCS: SFR_MCS Checks...\n");

  for(unsigned int i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 4)
        if(StarP[P[i].AuxDataID].PID != i)
          {
            printf("StarP broken: %llu %u %d %d %d %d\n", (long long)P[i].AuxDataID, i, StarP[P[i].AuxDataID].PID, NumGas, N_star,
                   NumPart);
            terminate("StarP[P[i].AuxDataID].PID!=i\n");
          }
    }

  mpi_printf("SFR_MCS: done.\n");
}
#endif

/* add a star particle to StarP */
void sfr_mcs_add_star(int i, int j, MyDouble mass_of_star, MyFloat birthtime)
{
  if(N_star >= All.MaxPartStar)
    terminate("There is no space left to create new stars. N_star = %d, MaxPartStar = %d", N_star, All.MaxPartStar);

  P[i].AuxDataID = N_star;

  /* zero StarP[] entries */
  memset(&StarP[N_star], 0, sizeof(struct star_particle_data));

  /* set values */
  StarP[N_star].PID       = i;
  StarP[N_star].BirthTime = birthtime;
#if defined(SN_MCS) || defined(HII_MCS) || defined(PE_MCS)
#ifdef SFR_MCS_DELAY
#if SFR_MCS_DELAY == 1
  StarP[N_star].TimeDelay = All.TimeDelayFactor;
#else
  StarP[N_star].TimeDelay = sqrt(3.0 * M_PI / (32.0 * All.G * SphP[j].Density * All.cf_a3inv)) * All.TimeDelayFactor;
  StarP[N_star].TimeDelay *= (All.UnitTime_in_s / SEC_PER_YEAR);
#endif
  StarP[N_star].Age = -StarP[N_star].TimeDelay;
#else
  StarP[N_star].Age = 0.0;
#endif
#endif

#ifdef SFR_MCS_BIRTH_RECORDS
  StarP[N_star].BirthPos[0]  = P[i].Pos[0];
  StarP[N_star].BirthPos[1]  = P[i].Pos[1];
  StarP[N_star].BirthPos[2]  = P[i].Pos[2];
  StarP[N_star].BirthVel[0]  = P[i].Vel[0];
  StarP[N_star].BirthVel[1]  = P[i].Vel[1];
  StarP[N_star].BirthVel[2]  = P[i].Vel[2];
  StarP[N_star].BirthDensity = SphP[j].Density;
#endif
#ifdef SFR_MCS_LOG_DETAILS
#ifdef GRACKLE
  double temp = get_temp_individual_cell_grackle(j);
#else
  double temp       = -1;
#endif
  fprintf(FdSFDetails, "STAR=%llu %g %g %g %g %g %g %g %g %g %g\n", (long long)P[i].ID, All.Time, P[i].Pos[0], P[i].Pos[1],
          P[i].Pos[2], P[i].Vel[0], P[i].Vel[1], P[i].Vel[2], SphP[j].Density, temp, SphP[j].Metallicity);
#endif

#ifdef SN_MCS
  StarP[N_star].N_SN           = 0;
  StarP[N_star].N_SN_cum       = 0;
  StarP[N_star].N_SN_event_cum = 0;
  StarP[N_star].InitialMass    = mass_of_star;

#ifdef SN_MCS_LOCATION_RECORDS
  StarP[N_star].SNPos[0]  = -1;
  StarP[N_star].SNPos[1]  = -1;
  StarP[N_star].SNPos[2]  = -1;
  StarP[N_star].SNTime    = -1;
  StarP[N_star].SNDensity = -1;
#ifdef GRACKLE
  StarP[N_star].SNTemperature = -1;
#endif
#endif
#endif  // SN_MCS

#ifdef HII_MCS
  StarP[N_star].S_Hii = 0.0;
#ifndef IMF_SAMPLING_MCS
  StarP[N_star].photon_it_high = 1;
#endif
#endif

#ifdef PE_MCS
  StarP[N_star].L_FUV = 0.0;
#ifndef IMF_SAMPLING_MCS
  StarP[N_star].fuv_it_high = 1;
#endif
#ifdef PE_MCS_STORE_DUST_COLUMN
  StarP[N_star].DustColumn = 0;
#endif
#endif

#ifdef METALS
  P[i].Metallicity = SphP[j].Metallicity;
#if(defined(SN_MCS) || (defined(HII_MCS) && !defined(HII_MCS_TEST)) || defined(PE_MCS)) && \
    !(defined(SB99_FIXED_Z) || defined(IMF_SAMPLING_MCS))
  StarP[N_star].iz = get_sb99_z_index(SphP[j].Metallicity);
#endif
#endif

  N_star++;
}
#endif

void sfr_criteria_announce(void)
{
#if SFR_MCS_SELECT_CRITERIA == 0
  mpi_printf("SFR_MCS: Using SFR_MCS_SELECT_CRITERIA 0. Fixed density threshold %g\n", All.DensThreshold);
#elif SFR_MCS_SELECT_CRITERIA == 1
  mpi_printf("SFR_MCS: Using SFR_MCS_SELECT_CRITERIA 1. Fixed Jeans length criteria, %g\n", All.SfrCritLength);
#elif SFR_MCS_SELECT_CRITERIA == 2
  mpi_printf("SFR_MCS: Using SFR_MCS_SELECT_CRITERIA 2. Jeans mass criteria, Mcell > %g MJeans for SF\n", All.SfrCritJeansMassN);
#elif SFR_MCS_SELECT_CRITERIA == 3
  mpi_printf(
      "SFR_MCS: Using SFR_MCS_SELECT_CRITERIA 3. Fixed density threshold  %g and Jeans mass criteria, Mcell > %g MJeans for SF\n",
      All.DensThreshold, All.SfrCritJeansMassN);
#endif

#if SFR_MCS_RATE_CRITERIA == 0
  mpi_printf("SFR_MCS: Using SFR_MCS_RATE_CRITERIA 0. Fixed efficiency %g\n", All.SfrEfficiency);
#elif SFR_MCS_RATE_CRITERIA == 1
  mpi_printf("SFR_MCS: Using SFR_MCS_RATE_CRITERIA 1. Padoan2012 with wind efficiency %g\n", All.SfrEfficiency);
#elif SFR_MCS_RATE_CRITERIA == 2
  mpi_printf("SFR_MCS: Using SFR_MCS_RATE_CRITERIA 2. Padoan2012 (including sound speed) with wind efficiency %g\n",
             All.SfrEfficiency);
#elif SFR_MCS_RATE_CRITERIA == 3
  mpi_printf("SFR_MCS: Using SFR_MCS_RATE_CRITERIA 3. Virial criteria with fixed efficiency %g\n", All.SfrEfficiency);
#elif SFR_MCS_RATE_CRITERIA == 4
  mpi_printf("SFR_MCS: Using SFR_MCS_RATE_CRITERIA 4. Virial criteria (including sound speed) with fixed efficiency %g\n",
             All.SfrEfficiency);
#endif
}
#endif /*SFR_MCS*/
