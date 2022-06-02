/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/GFM/stellar_evolution_main.c
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

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef GFM_STELLAR_EVOLUTION

#define WARN_TOLERANCE 0.00001

/*!
 * driver routine for stellar evolution, set up arguments for stellar_evolution(...) call
 */
void do_stellar_evolution(MyFloat age_of_star_in_Gyr, MyFloat dtime_in_Gyr, int iPart, stellar_evolution_data* sed)
{
  int k_elem;
  MyFloat initial_metals[GFM_N_CHEM_ELEMENTS]; /* original element abundances */
  MyFloat metallicity;
  MyDouble initial_mass;
  MyDouble log_min_mass, log_max_mass, log_metallicity;

  if(P[iPart].Type != 4)
    terminate("GFM_STELLAR_EVOLUTION: Not a stellar particle.\n");

#ifdef GFM_RPROCESS_CHANNELS_NS_KICKS
  if(P[iPart].Mass <= 0 && STP(iPart).NSNS_channel == -1)
#else
  if(P[iPart].Mass <= 0)
#endif
    terminate("GFM_STELLAR_EVOLUTION: Bad mass.\n");

  if(STP(iPart).BirthTime <= 0)
    terminate("GFM_STELLAR_EVOLUTION: Wind particle.\n");

  /* initial fields from star particle properties */
  initial_mass = STP(iPart).InitialMass;
  for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
    initial_metals[k_elem] = STP(iPart).MassMetals[k_elem] / P[iPart].Mass;

#ifdef GFM_RPROCESS_CHANNELS_NS_KICKS
  if(STP(iPart).NSNS_channel >= 0)
    for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
      initial_metals[k_elem] = 0;
#endif

  metallicity = STP(iPart).Metallicity;
  if(metallicity > 0)
    log_metallicity = fmax(log10(metallicity), GFM_MIN_METAL);
  else
    log_metallicity = GFM_MIN_METAL;

  /* zero variables and arrays */
  sed->mass_return_flag          = 0;
  sed->total_mass_released       = 0;
  sed->total_metal_mass_released = 0;
  for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
    sed->metal_mass_released[k_elem] = 0;

#if defined(GFM_DUST) || (defined(DUST_LIVE) && defined(DL_PRODUCTION))
  sed->total_dust_mass_released = 0.0;
#endif
#ifdef GFM_DUST
  for(int chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
    {
      for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
        {
          sed->dust_mass_released[chan][k_elem] = 0.0;
        }
    }
#endif
#if defined(DUST_LIVE) && defined(DL_PRODUCTION)
  for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
    {
      sed->dust_mass_released[k_elem] = 0.0;
    }
#endif
#ifdef GFM_CHEMTAGS
  for(k_elem = 0; k_elem < GFM_N_CHEM_TAGS; k_elem++)
    {
      sed->metal_mass_released_chemtags[k_elem] = 0;
    }
#endif
#ifdef GFM_RPROCESS_CHANNELS
  int idx;
  for(idx = 0; idx < GFM_RPROCESS_CHANNELS; idx++)
    {
      sed->rprocess_mass_released[idx] = 0;
      sed->number_of_NSNS[idx]         = 0;
    }
#endif

#ifdef GFM_SNIA_ENERGY_INJECTION
  sed->NumSNIa = 0;
#endif

  sed->number_of_SNIa     = 0;
  sed->number_of_SNII     = 0;
  sed->AGB_mass_released  = 0.0;
  sed->SNIa_mass_released = 0.0;
  sed->SNII_mass_released = 0.0;
#ifdef SMUGGLE_AGB_WINDS
  sed->OB_mass_released = 0.0;
  sed->OB_return_flag   = 0;
#endif

#ifdef SMUGGLE_DISCRETE_SN
  sed->initial_mass        = initial_mass;
  sed->remaining_mass_frac = P[iPart].Mass / sed->initial_mass;
#endif

#ifdef GFM_RPROCESS_CHANNELS_NS_KICKS
  if(STP(iPart).NSNS_channel >= 0)
    {
      if(age_of_star_in_Gyr < STP(iPart).DelayTime && age_of_star_in_Gyr + dtime_in_Gyr >= STP(iPart).DelayTime)
        {
          if(STP(iPart).NSNS_channel >= GFM_RPROCESS_NSNS)
            terminate("Task=%d, part=%d, channel=%d", ThisTask, iPart, STP(iPart).NSNS_channel);

          sed->rprocess_mass_released[STP(iPart).NSNS_channel] =
              All.rp[STP(iPart).NSNS_channel].NSNS_MassPerEvent / All.UnitMass_in_g * SOLAR_MASS * All.HubbleParam;
          sed->number_of_NSNS[STP(iPart).NSNS_channel] = 1;

          if(!gsl_finite(sed->rprocess_mass_released[STP(iPart).NSNS_channel]))
            terminate("Task=%d, iPart=%d, Channel=%d, MassPerEvent=%g\n", ThisTask, iPart, STP(iPart).NSNS_channel,
                      All.rp[STP(iPart).NSNS_channel].NSNS_MassPerEvent);

          // terminate( "NSNS: Task=%d, part=%d, ID=%lld, channel=%d\n", ThisTask, iPart, (long long)P[iPart].ID,
          // STP(iPart).NSNS_channel );

          STP(iPart).NSNS_channel = -1;

          // printf( "NSNS: Task=%d, part=%d, channel=%d, mass=%g\n", ThisTask, iPart, STP(iPart).NSNS_channel,
          // sed->rprocess_mass_released[STP(iPart).NSNS_channel] );
        }
      else
        {
          // printf( "task=%d, ID=%lld, age_of_star_in_Gyr=%g, STP(iPart).DelayTime=%g, age_of_star_in_Gyr + dtime_in_Gyr=%g\n",
          // ThisTask, (long long)P[iPart].ID, age_of_star_in_Gyr, STP(iPart).DelayTime, age_of_star_in_Gyr + dtime_in_Gyr );
          return;
        }
    }
#endif

  /* minimum and maximum mass of stars that will die during this time-step */
  if(P[iPart].Mass > 0)
    {
      if(get_dying_mass_in_Msun(age_of_star_in_Gyr, metallicity) <= 0 ||
         get_dying_mass_in_Msun(age_of_star_in_Gyr + dtime_in_Gyr, metallicity) <= 0)
        terminate("GFM_STELLAR_EVOLUTION: Negative mass for dying stars");

      log_max_mass = log10(get_dying_mass_in_Msun(age_of_star_in_Gyr, metallicity));
      log_min_mass = log10(get_dying_mass_in_Msun(age_of_star_in_Gyr + dtime_in_Gyr, metallicity));

      if(log_min_mass > log_max_mass)
        terminate(
            "GFM_STELLAR_EVOLUTION: Min mass larger than max mass for stellar evolution: i=%d  log_min_mass=%g  log_max_mass=%g  "
            "diff=%g, dtime_in_Gyr=%g\n",
            iPart, log_min_mass, log_max_mass, log_max_mass - log_min_mass, dtime_in_Gyr);

      if(log_min_mass == log_max_mass)
        return;

        /* initialize the IMF */
#ifdef GFM_VARIABLE_IMF
#if(GFM_VARIABLE_IMF == 0)
      define_imf(STP(iPart).DMVelDisp);
#else
#error "GFM_VARIABLE_IMF mode is not ok"
#endif
#endif

      /* perform stellar evolution, returns: total mass, metall mass, element mass, SNIa, SNII rates via sed */
      evolve_AGB(log_min_mass, log_max_mass, log_metallicity, initial_metals, sed);

      sed->AGB_mass_released = sed->total_mass_released;

#if defined(SMUGGLE_DISCRETE_SN) && defined(GFM_DISCRETE_ENRICHMENT)
      double age_SN_in_Gyr   = get_time_difference_in_Gyr(STP(iPart).BirthTime, STP(iPart).SNTime);
      double dtime_SN_in_Gyr = dtime_in_Gyr + (age_of_star_in_Gyr - age_SN_in_Gyr);

      if(dtime_SN_in_Gyr < 0.0)
        terminate("Invalid time step value");

      evolve_SNIa(log_min_mass, log_max_mass, age_SN_in_Gyr, dtime_SN_in_Gyr, metallicity, initial_metals, sed);

      sed->SNIa_mass_released = sed->total_mass_released - sed->AGB_mass_released;

      double log_max_mass_SN = get_dying_mass_in_Msun(age_SN_in_Gyr, metallicity);
      double log_min_mass_SN = get_dying_mass_in_Msun(age_SN_in_Gyr + dtime_SN_in_Gyr, metallicity);

      if(log_max_mass_SN <= 0 || log_min_mass_SN <= 0)
        terminate("GFM_STELLAR_EVOLUTION: Negative mass for dying stars in SNII channel");

      log_max_mass_SN = log10(log_max_mass_SN);
      log_min_mass_SN = log10(log_min_mass_SN);

      if(log_min_mass > log_max_mass)
        terminate(
            "GFM_STELLAR_EVOLUTION: Min mass larger than max mass for SNII evolution: i=%d  log_min_mass=%g  log_max_mass=%g  "
            "diff=%g\n",
            iPart, log_min_mass, log_max_mass, log_max_mass - log_min_mass);

      evolve_SNII(log_min_mass_SN, log_max_mass_SN, log_metallicity, initial_metals, sed);

      sed->SNII_mass_released = sed->total_mass_released - sed->AGB_mass_released - sed->SNIa_mass_released;
#else
      evolve_SNIa(log_min_mass, log_max_mass, age_of_star_in_Gyr, dtime_in_Gyr, metallicity, initial_metals, sed);

      sed->SNIa_mass_released = sed->total_mass_released - sed->AGB_mass_released;

      evolve_SNII(log_min_mass, log_max_mass, log_metallicity, initial_metals, sed);

      sed->SNII_mass_released = sed->total_mass_released - sed->AGB_mass_released - sed->SNIa_mass_released;
#endif

#ifdef SMUGGLE_AGB_WINDS
      double age_OB_in_Gyr   = get_time_difference_in_Gyr(STP(iPart).BirthTime, STP(iPart).OBTime);
      double dtime_OB_in_Gyr = dtime_in_Gyr + (age_of_star_in_Gyr - age_OB_in_Gyr);

      if(dtime_OB_in_Gyr < 0.0)
        terminate("Invalid time step value");

      evolve_OB_winds(log_min_mass, log_max_mass, age_OB_in_Gyr, dtime_OB_in_Gyr, log_metallicity, initial_metals, sed);

      sed->OB_mass_released = sed->total_mass_released - sed->AGB_mass_released - sed->SNIa_mass_released - sed->SNII_mass_released;
#endif

      /* convert mass and metals_released to code units */
      sed->total_mass_released *= initial_mass;
      sed->total_metal_mass_released *= initial_mass;
      for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
        sed->metal_mass_released[k_elem] *= initial_mass;

        // Note there is no:  NSNS_mass_released parameter ... assuming Eu contributes nothing to mass or energy dump
#ifdef GFM_RPROCESS
      evolve_NSNS(log_min_mass, log_max_mass, age_of_star_in_Gyr, dtime_in_Gyr, metallicity, initial_mass, sed);
#endif

#if defined(GFM_RPROCESS_CHANNELS) && !defined(GFM_RPROCESS_CHANNELS_NS_KICKS)
      for(idx = 0; idx < GFM_RPROCESS_CHANNELS; idx++)
        sed->rprocess_mass_released[idx] = sed->total_mass_released * STP(iPart).MassRProcess[idx] / P[iPart].Mass;

      /* only high resolution star particles do r-process injections */
      if(initial_mass < 2. * All.TargetGasMass)
        evolve_NSNS(log_min_mass, log_max_mass, age_of_star_in_Gyr, dtime_in_Gyr, metallicity, initial_mass, sed);
#endif
    }

#if defined(GFM_DUST) || (defined(DUST_LIVE) && defined(DL_PRODUCTION))
  sed->total_dust_mass_released *= initial_mass;
#endif
#ifdef GFM_DUST
  for(int chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
    {
      for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
        {
          sed->dust_mass_released[chan][k_elem] *= initial_mass;
        }
    }
#endif
#if defined(DUST_LIVE) && defined(DL_PRODUCTION)
  for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
    {
      sed->dust_mass_released[k_elem] *= initial_mass;
    }
#endif
  sed->AGB_mass_released *= initial_mass;
  sed->SNIa_mass_released *= initial_mass;
  sed->SNII_mass_released *= initial_mass;

#ifdef SMUGGLE_AGB_WINDS
  sed->OB_mass_released *= initial_mass;
#endif

#ifdef GFM_CHEMTAGS
  sed->metal_mass_released_chemtags[GFM_SNIA_CHEMTAG] *= initial_mass;
  sed->metal_mass_released_chemtags[GFM_SNII_CHEMTAG] *= initial_mass;
  sed->metal_mass_released_chemtags[GFM_AGB_CHEMTAG] *= initial_mass;
#ifdef GFM_RPROCESS
  sed->metal_mass_released_chemtags[GFM_NSNS_CHEMTAG] *= initial_mass;
#endif
#ifdef GFM_SPLITFE
  sed->metal_mass_released_chemtags[GFM_FESNIA_CHEMTAG] *= initial_mass;
  sed->metal_mass_released_chemtags[GFM_FESNII_CHEMTAG] *= initial_mass;
#endif
#endif

  /* do same for number of SN of each type */
  sed->number_of_SNIa *= initial_mass * All.UnitMass_in_g / SOLAR_MASS / All.HubbleParam;
  sed->number_of_SNII *= initial_mass * All.UnitMass_in_g / SOLAR_MASS / All.HubbleParam;

#ifdef GFM_SNIA_ENERGY_INJECTION
  if(initial_mass < 2. * All.TargetGasMass)
    sed->NumSNIa = gsl_ran_poisson(random_generator, sed->number_of_SNIa);
#endif

#ifdef GFM_RPROCESS_CHANNELS
#if GFM_RPROCESS_NSNS < GFM_RPROCESS_CHANNELS
  /* only high resolution star particles do r-process injections */
  if(initial_mass < 2. * All.TargetGasMass)
    for(idx = GFM_RPROCESS_NSNS; idx < GFM_RPROCESS_CHANNELS; idx++)
      {
        double NSNRP   = sed->number_of_SNII * All.rpSN[idx - GFM_RPROCESS_NSNS].RPSN_FractionPerSN;
        int NSNRP_true = gsl_ran_poisson(random_generator, NSNRP);

        sed->number_of_NSNS[idx] += NSNRP_true;
        sed->rprocess_mass_released[idx] +=
            NSNRP_true * All.rpSN[idx - GFM_RPROCESS_NSNS].RPSN_MassPerEvent / All.UnitMass_in_g * SOLAR_MASS * All.HubbleParam;

        if(!gsl_finite(sed->rprocess_mass_released[idx]))
          {
            terminate("SHIT: idx=%d, res=%g, number=%d, mass per event=%g\n", idx, sed->rprocess_mass_released[idx],
                      sed->number_of_NSNS[idx], All.rpSN[idx - GFM_RPROCESS_NSNS].RPSN_MassPerEvent);
          }
      }
#endif
#endif

  if(sed->total_mass_released > P[iPart].Mass)
    {
#ifdef SMUGGLE_AGB_WINDS
      terminate(
          "GFM_STELLAR_EVOLUTION: released mass larger than stellar mass: ID=%lld P[iPart].Mass=%g sed->total_mass_released=%g "
          "initial_mass=%g age_Gyr=%g dt_Gyr=%g SNIImass=%g, SNIaMass=%g, AGBmass=%g OBmass=%g\n",
          (long long)P[iPart].ID, P[iPart].Mass, sed->total_mass_released, initial_mass, age_of_star_in_Gyr, dtime_in_Gyr,
          sed->SNII_mass_released, sed->SNIa_mass_released, sed->AGB_mass_released, sed->OB_mass_released);
#else
      terminate(
          "GFM_STELLAR_EVOLUTION: released mass larger than stellar mass: ID=%lld P[iPart].Mass=%g sed->total_mass_released=%g "
          "initial_mass=%g age_Gyr=%g dt_Gyr=%g",
          (long long)P[iPart].ID, P[iPart].Mass, sed->total_mass_released, initial_mass, age_of_star_in_Gyr, dtime_in_Gyr);
#endif
    }

  if(sed->total_mass_released < 0)
    {
      if(sed->total_mass_released < -WARN_TOLERANCE * initial_mass)
        terminate(
            "GFM_STELLAR_EVOLUTION: negative mass returned: ID=%lld P[iPart].Mass=%g sed->total_mass_released=%g initial_mass=%g "
            "age_Gyr=%g dt_Gyr=%g\n",
            (long long)P[iPart].ID, P[iPart].Mass, sed->total_mass_released, initial_mass, age_of_star_in_Gyr, dtime_in_Gyr);

      sed->mass_return_flag = 1;
    }

  if(sed->total_metal_mass_released > sed->total_mass_released)
    {
      if(sed->total_metal_mass_released > (1.0 + WARN_TOLERANCE) * sed->total_mass_released)
        terminate(
            "GFM_STELLAR_EVOLUTION: released metal mass larger than total released mass: ID=%lld total_metal_mass_released=%g "
            "total_mass_released=%g",
            (long long)P[iPart].ID, sed->total_metal_mass_released, sed->total_mass_released);

      sed->mass_return_flag = 1;
    }

  if(sed->total_metal_mass_released < 0)
    {
      if(sed->total_metal_mass_released < -WARN_TOLERANCE * initial_mass)
        terminate(
            "GFM_STELLAR_EVOLUTION: negative metal mass returned: ID=%lld P[iPart].Mass=%g sed->total_metal_mass_released=%g "
            "initial_mass=%g\n",
            (long long)P[iPart].ID, P[iPart].Mass, sed->total_metal_mass_released, initial_mass);

      sed->mass_return_flag = 1;
    }

  for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++) /* H, He and metals */
    {
      if(sed->metal_mass_released[k_elem] < 0)
        {
#ifndef GFM_NO_NEGATIVE_ELEMENT_MASS_RELEASED
          if(sed->metal_mass_released[k_elem] < -WARN_TOLERANCE * initial_mass)
#endif
            {
              char buf[2000];
              sprintf(buf,
                      "GFM_STELLAR_EVOLUTION: negative metal mass element returned: ID=%lld P[iPart].Mass=%g "
                      "sed->metal_mass_released[k_elem]=%g initial_mass=%g k_elem=%d  age_of_star_in_Gyr=%g  dtime_in_Gyr=%g "
                      "metallicity=%g  initial_metals=",
                      (long long)P[iPart].ID, P[iPart].Mass, sed->metal_mass_released[k_elem], initial_mass, k_elem,
                      age_of_star_in_Gyr, dtime_in_Gyr, metallicity);

              for(int kk = 0; kk < GFM_N_CHEM_ELEMENTS; kk++)
                {
                  char buf1[100];
                  sprintf(buf1, "%g ", initial_metals[kk]);
                  strcat(buf, buf1);
                }
              strcat(buf, " metal_mass_released=");
              for(int kk = 0; kk < GFM_N_CHEM_ELEMENTS; kk++)
                {
                  char buf1[100];
                  sprintf(buf1, "%g ", sed->metal_mass_released[kk]);
                  strcat(buf, buf1);
                }
              strcat(buf, "\n");
              terminate(buf);
            }

          sed->mass_return_flag = 1;
        }
    }

  if(sed->mass_return_flag)
    {
      sed->total_mass_released       = 0;
      sed->total_metal_mass_released = 0;
      sed->AGB_mass_released         = 0;
#ifdef SMUGGLE_AGB_WINDS
      sed->OB_mass_released = 0;
      sed->OB_return_flag   = 0;
#endif
      sed->SNIa_mass_released = 0;
      sed->SNII_mass_released = 0;
      sed->number_of_SNIa     = 0;
      sed->number_of_SNII     = 0;
      for(int k_elem2 = 0; k_elem2 < GFM_N_CHEM_ELEMENTS; k_elem2++)
        sed->metal_mass_released[k_elem2] = 0;
#ifdef GFM_CHEMTAGS
      for(int k_elem3 = 0; k_elem3 < GFM_N_CHEM_TAGS; k_elem3++)
        sed->metal_mass_released_chemtags[k_elem3] = 0;
#endif
    }

#ifdef GFM_DUST
  for(int chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
    {
      for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
        {
          if(sed->dust_mass_released[chan][k_elem] < 0.0)
            {
              sed->dust_mass_released[chan][k_elem] = 0.0;
            }
        }
    }
#endif

#if GFM_STELLAR_EVOLUTION == 1
  sed->total_mass_released = 0;
#endif

#ifdef GFM_DISCRETE_ENRICHMENT
  double delta_M              = sed->total_mass_released / P[iPart].Mass;
  double discrete_return_frac = 0.0001;  // All.GFM_DiscreteEnrichFrac
  double discrete_min_age     = 0.1;     // Gyr, below which all returns are allowed (for fast evolution of early stellar populations)

#ifndef SMUGGLE_DISCRETE_SN
  if(delta_M < discrete_return_frac && age_of_star_in_Gyr > discrete_min_age)
#else
  int SN_events = 0;
  SN_events     = (int)(sed->number_of_SNII + sed->number_of_SNIa);

  /* the additional condition is needed for proper SN distribution sampling;  */
  /* it forces enrichment to be done if one or more SN events occur.          */
  if((delta_M < discrete_return_frac) && (age_of_star_in_Gyr > discrete_min_age) && (SN_events == 0))
#endif
    {
      // zero all enrichment returns and quit early
      sed->mass_return_flag          = 0;
      sed->total_mass_released       = 0;
      sed->total_metal_mass_released = 0;
      for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
        sed->metal_mass_released[k_elem] = 0;

#ifdef GFM_CHEMTAGS
      for(k_elem = 0; k_elem < GFM_N_CHEM_TAGS; k_elem++)
        sed->metal_mass_released_chemtags[k_elem] = 0;
#endif

      sed->number_of_SNIa     = 0;
      sed->number_of_SNII     = 0;
      sed->AGB_mass_released  = 0.0;
      sed->SNIa_mass_released = 0.0;
      sed->SNII_mass_released = 0.0;
#ifdef SMUGGLE_AGB_WINDS
      sed->OB_mass_released = 0.0;
      sed->OB_return_flag   = 0;
#endif
      return;
    }
#endif /* GFM_DISCRETE_ENRICHMENT */

#if defined(DUST_LIVE) && defined(DL_PRODUCTION)
  if(sed->total_dust_mass_released < 0.0)
    sed->total_dust_mass_released = 0.0;

  for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
    {
      if(sed->dust_mass_released[k_elem] < 0.0)
        {
          sed->dust_mass_released[k_elem] = 0.0;
        }
    }
#endif

  /* fractional mass loss */
  double frac = 1.0;
  if(P[iPart].Mass > 0)
    frac -= sed->total_mass_released / P[iPart].Mass;

  /* scale mass and metal content */
  P[iPart].Mass *= frac;

  for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
    STP(iPart).MassMetals[k_elem] *= frac;

#ifdef GFM_CHEMTAGS
  // NOTE: we assume NSNS production mass is negligable - does not factor into the total mass released
  for(int k_elem = 0; k_elem < GFM_N_CHEM_TAGS; k_elem++)
    STP(iPart).MassMetalsChemTags[k_elem] *= frac;
#endif
#ifdef GFM_RPROCESS_CHANNELS
  for(idx = 0; idx < GFM_RPROCESS_CHANNELS; idx++)
    {
      STP(iPart).MassRProcess[idx] *= frac;
    }
#endif
}

void evolve_active_stars(void)
{
  if(All.HighestActiveTimeBin == 0)
    return;

  int i, iel;
  double time_begstep, age_of_star_in_Gyr, dtime_in_Gyr;
  stellar_evolution_data sed;
  double AGBLocalMassReleased  = 0.0;
  double SNIaLocalMassReleased = 0.0;
  double SNIILocalMassReleased = 0.0;

#ifdef SMUGGLE_STAR_FEEDBACK
  MyFloat SNII_Number;
  MyFloat SNIa_Number;
#endif
#ifdef SMUGGLE_RADIATION_FEEDBACK
  double age_of_star_in_Gyr_rad, dtime_in_Gyr_rad;
#endif

#ifdef GFM_DUST
  int chan;
#endif

  /* don't do this on the first timestep, because dtime_in_Gyr=0, but roundoff errors can lead to dtime_in_Gyr<0 */
  if(All.NumCurrentTiStep == 0)
    return;

  /* do local star particles */
  for(i = 0; i < Nstar; i++)
    {
      /* (Note: All.Ti_begstep[] has already been advanced for the next step at this point) */
      if(All.ComovingIntegrationOn)
        time_begstep = All.TimeBegin * exp(All.Ti_begstep[P[StarParticle[i].index].TimeBinGrav] * All.Timebase_interval);
      else
        time_begstep = All.TimeBegin + All.Ti_begstep[P[StarParticle[i].index].TimeBinGrav] * All.Timebase_interval;

      /* we take the age at the beginning of the interval, and dtime is since the last enrichment event */
      age_of_star_in_Gyr = get_time_difference_in_Gyr(STP(StarParticle[i].index).BirthTime, STP(StarParticle[i].index).lastEnrichTime);
      dtime_in_Gyr       = get_time_difference_in_Gyr(STP(StarParticle[i].index).lastEnrichTime, time_begstep);

#ifdef SMUGGLE_RADIATION_FEEDBACK
      dtime_in_Gyr_rad       = dtime_in_Gyr;
      age_of_star_in_Gyr_rad = age_of_star_in_Gyr;
#endif

#ifdef GFM_DISCRETE_ENRICHMENT
      /* we take the age at the beginning of the interval, and dtime is since the last enrichment event */
      age_of_star_in_Gyr = get_time_difference_in_Gyr(STP(StarParticle[i].index).BirthTime, STP(StarParticle[i].index).lastEnrichTime);
      dtime_in_Gyr       = get_time_difference_in_Gyr(STP(StarParticle[i].index).lastEnrichTime, time_begstep);
#endif

      do_stellar_evolution(age_of_star_in_Gyr, dtime_in_Gyr, StarParticle[i].index, &sed);

#ifdef GFM_DISCRETE_ENRICHMENT
      /* if we will do an enrichment event, update the lastEnrichTime to the end of the current timestep */
      if(sed.total_mass_released > 0.0 || sed.mass_return_flag == 1)
        STP(StarParticle[i].index).lastEnrichTime = time_begstep;

#ifdef SMUGGLE_DISCRETE_SN
      /* this ensures that time variable flows normally for discrete SN sampling */
      STP(StarParticle[i].index).SNTime = time_begstep;
#endif

#else
      STP(StarParticle[i].index).lastEnrichTime = time_begstep;
#endif

#ifdef SMUGGLE_AGB_WINDS
      /* if OB winds actually return mass, update the time variable to the end of current time step */
      if(sed.OB_mass_released > 0.0 || sed.mass_return_flag == 1 || sed.OB_return_flag == 1)
        STP(StarParticle[i].index).OBTime = time_begstep;
#endif

      AGBLocalMassReleased += sed.AGB_mass_released;
      SNIaLocalMassReleased += sed.SNIa_mass_released;
      SNIILocalMassReleased += sed.SNII_mass_released;

      if(dtime_in_Gyr > 0)
        {
          STP(StarParticle[i].index).SNIaRate = sed.number_of_SNIa / (dtime_in_Gyr * 1.e9);
          STP(StarParticle[i].index).SNIIRate = sed.number_of_SNII / (dtime_in_Gyr * 1.e9);
        }

      StarParticle[i].TotalMassReleased      = sed.total_mass_released;
      StarParticle[i].TotalMetalMassReleased = sed.total_metal_mass_released;
      for(iel = 0; iel < GFM_N_CHEM_ELEMENTS; iel++)
        {
          StarParticle[i].MetalMassReleased[iel] = sed.metal_mass_released[iel];
        }

#ifdef GFM_DUST
      for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
        {
          for(iel = 0; iel < GFM_N_CHEM_ELEMENTS; iel++)
            {
              StarParticle[i].DustMassReleased[chan][iel] = sed.dust_mass_released[chan][iel];
            }
        }
#endif
#if defined(GFM_DUST) || defined(DUST_LIVE)
      StarParticle[i].NumSNII = sed.number_of_SNII;
#endif
#if defined(DUST_LIVE) && defined(DL_PRODUCTION)
      STP(StarParticle[i].index).DeltaDustMassTot = 0.0;
      for(iel = 0; iel < GFM_N_CHEM_ELEMENTS; iel++)
        {
          STP(StarParticle[i].index).DeltaDustMass[iel] = sed.dust_mass_released[iel];
          STP(StarParticle[i].index).DeltaDustMassTot += sed.dust_mass_released[iel];
        }
      /* Use total mass returned by stars of different types as a proxy for
       * determining dominant contributor to grain size distribution. */
      if(sed.AGB_mass_released > fmax(sed.SNII_mass_released, sed.SNIa_mass_released))
        STP(StarParticle[i].index).DndaType = GSD_DNDA_AGB;
      else if(sed.SNII_mass_released > sed.SNIa_mass_released)
        STP(StarParticle[i].index).DndaType = GSD_DNDA_SNII;
      else
        STP(StarParticle[i].index).DndaType = GSD_DNDA_SNIA;
#endif

#ifdef GFM_CHEMTAGS
      for(iel = 0; iel < GFM_N_CHEM_TAGS; iel++)
        {
          StarParticle[i].MetalMassReleasedChemTags[iel] = sed.metal_mass_released_chemtags[iel];
        }
#endif
#ifdef GFM_RPROCESS_CHANNELS
      for(iel = 0; iel < GFM_RPROCESS_CHANNELS; iel++)
        {
          StarParticle[i].MassReleasedRProcess[iel] = sed.rprocess_mass_released[iel];
          STP(StarParticle[i].index).NRProcessInjections[iel] += sed.number_of_NSNS[iel];
        }
#endif
#ifdef GFM_SNIA_ENERGY_INJECTION
      StarParticle[i].NumSNIa = sed.NumSNIa;
      // STP(StarParticle[i].index).NumSNIa += sed.NumSNIa; // don't add the number here but in enrich and only for the Ia's that
      // actually injected energy
#endif

#ifdef GFM_STELLAR_FEEDBACK
      StarParticle[i].SNIaEnergyReleased  = sed.number_of_SNIa * All.EnergyPerSNIa * (All.HubbleParam / All.UnitEnergy_in_cgs);
      StarParticle[i].AGBMomentumReleased = sed.AGB_mass_released * All.AGBWindVelocity;
#endif
#ifdef GFM_WINDS_LOCAL
      StarParticle[i].WindEnergyReleased =
          sed.number_of_SNII * All.WindEnergyIn1e51erg * GFM_SNII_ENERGY * (All.HubbleParam / All.UnitEnergy_in_cgs);
#endif
#ifdef GFM_INJECT_B_FROM_SN
      double phi               = 2 * M_PI * get_random_number();
      double theta             = acos(2 * get_random_number() - 1);
      double bubbleRadius      = STP(StarParticle[i].index).Hsml;
      double bubbleRadiusCubed = bubbleRadius * bubbleRadius * bubbleRadius;
      double soft              = get_default_softening_of_particletype(0) * All.cf_atime;
      double softCubed         = soft * soft * soft;
      double mstrength = sqrt(3. * All.SupernovaInjectedMagneticEnergyInErgs * All.HubbleParam / All.UnitEnergy_in_cgs * softCubed *
                              (bubbleRadiusCubed + softCubed) / bubbleRadiusCubed) /
                         sqrt(4.0 * M_PI);

      /* to comoving units */
      mstrength *= All.cf_atime * All.cf_atime;

      StarParticle[i].m[0] = (sed.number_of_SNIa + sed.number_of_SNII) * mstrength * sin(theta) * cos(phi);
      StarParticle[i].m[1] = (sed.number_of_SNIa + sed.number_of_SNII) * mstrength * sin(theta) * sin(phi);
      StarParticle[i].m[2] = (sed.number_of_SNIa + sed.number_of_SNII) * mstrength * cos(theta);
#endif

#ifdef SMUGGLE_STAR_FEEDBACK
      SNII_Number = sed.number_of_SNII;
      SNIa_Number = sed.number_of_SNIa;

#ifdef SMUGGLE_SN_COOLING_RADIUS_BOOST
      StarParticle[i].NumSNII               = sed.number_of_SNII;
      StarParticle[i].NumSNIa               = sed.number_of_SNIa;
      StarParticle[i].TotalMassReleasedSNII = sed.SNII_mass_released;
      StarParticle[i].TotalMassReleasedSNIa = sed.SNIa_mass_released;

#ifdef SMUGGLE_AGB_WINDS
      double psi = 5.94e4 / (1. + pow(age_of_star_in_Gyr / 2.5e-3, 1.4) + pow(age_of_star_in_Gyr / 0.01, 5)) +
                   4.83; /*Hopkins+ (2018, Appendix A)*/
      StarParticle[i].AGBWindSpeed         = sqrt(2e12 * psi) / All.UnitVelocity_in_cm_per_s;
      StarParticle[i].TotalMassReleasedAGB = sed.AGB_mass_released + sed.OB_mass_released;
#else
      StarParticle[i].TotalMassReleasedAGB = sed.AGB_mass_released;
#endif

#endif

#ifdef SMUGGLE_VAR_SN_EFF /*PAM*/
      double metallicity_loc = StarParticle[i].AvgMetalNgb;
      double metal_floor     = 0.001 * GFM_SOLAR_METALLICITY;
      if(metallicity_loc < metal_floor)
        metallicity_loc = metal_floor;

      double sn_eff = 0.3 + (3.0 - 0.3) / (1. + metallicity_loc / (0.1 * GFM_SOLAR_METALLICITY));
      SNII_Number *= sn_eff; /* WARNING: now is not the SN number but is weighted by SN efficiency */
#endif
      StarParticle[i].TotalEnergyReleased = compute_SN_energy(SNII_Number + SNIa_Number);
#ifdef SMUGGLE_OUTPUT_STELLAR_FEEDBACK
      STP(StarParticle[i].index).SNII_Num = SNII_Number;
      STP(StarParticle[i].index).Cum_SNII_Num += SNII_Number;
      STP(StarParticle[i].index).SNIa_Num = SNIa_Number;
      STP(StarParticle[i].index).Cum_SNIa_Num += SNIa_Number;
      STP(StarParticle[i].index).FeedbackEnergy    = StarParticle[i].TotalEnergyReleased;
      STP(StarParticle[i].index).TotalMassReleased = StarParticle[i].TotalMassReleased;
      STP(StarParticle[i].index).TotalMassToEnrich = StarParticle[i].NumNgb;
      STP(StarParticle[i].index).MaxFeedRadius     = StarParticle[i].FeedbackRadiusLimiter;
#endif
#endif

#ifdef SMUGGLE_RADIATION_FEEDBACK

      if(age_of_star_in_Gyr_rad >= All.InputTimeHeatRadiationFeedback)
        {
          StarParticle[i].RadiationMomentumReleased = 0.;
          StarParticle[i].RadCoolShutoffTime        = 0.;
          StarParticle[i].GasColumnDensity          = -999.;
          StarParticle[i].RadFeedTau                = -999.;
#ifdef SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION
          StarParticle[i].StromgrenMass   = -999.;
          StarParticle[i].StromgrenRadius = -999.;
          StarParticle[i].Lum             = 0.;
#endif
        }
      else
        {
          // so that feedback is injected exactly InputTimeMomRadiationFeedback
          dtime_in_Gyr_rad = fmin(dtime_in_Gyr_rad, All.InputTimeMomRadiationFeedback - age_of_star_in_Gyr_rad);

          double Lum = All.LumToMassRatioRadiationFeedback * P[StarParticle[i].index].Mass; /* luminosity in code units */
          Lum *=
              (All.UnitMass_in_g / SOLAR_MASS) * SOLAR_LUM / All.HubbleParam; /* in CGS, HubbleParam because of mass scaling, LVS?? */
          double avgGasdens = StarParticle[i].LocISMdens;                     /* comoving density */

          double minGasDist = StarParticle[i].RadFeed_MinGasDist * All.cf_atime; /* physical */
          double hsml       = fmin(STP(StarParticle[i].index).Hsml, StarParticle[i].FeedbackRadiusLimiter);

#ifdef SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION
          /* for stromgren calculation:  feed in physical and return physical */
          StarParticle[i].StromgrenRadius = compute_stromgren_radius(Lum, avgGasdens * All.cf_a3inv, minGasDist, hsml * All.cf_atime);
          StarParticle[i].StromgrenRadius /= All.cf_atime; /*comoving code units */
          /* for strongren calculation: feed in physical and return physical (mass) */
          StarParticle[i].StromgrenMass = compute_stromgren_mass(Lum, avgGasdens * All.cf_a3inv, hsml * All.cf_atime);

          if(age_of_star_in_Gyr_rad <= All.InputTimeMomRadiationFeedback)
            {
              StarParticle[i].RadiationMomentumReleased = Lum / CLIGHT * dtime_in_Gyr_rad * SEC_PER_GIGAYEAR; /* gr (cm/s) h^-1 ... */
              StarParticle[i].Lum                       = Lum;                                                /* gr (cm/s) h^-1 ... */
              StarParticle[i].RadiationMomentumReleased /=
                  (All.UnitVelocity_in_cm_per_s *
                   All.UnitMass_in_g); /* code units, although missing cf_atime, added at use in radiation_stellar_feedback.c */
              StarParticle[i].RadiationMomentumReleased *= All.HubbleParam;

              STP(StarParticle[i].index).RadFeed_Flag = 0; /* inputs all momentum only once */
            }
          else
            {
              StarParticle[i].RadiationMomentumReleased = 0.;
              StarParticle[i].Lum                       = 0.;
            }

#if(SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION == 1)  // Only allow stars to photoionize material once
          if(STP(StarParticle[i].index).PhotoionizationAttempts == 0)
            {
              StarParticle[i].RadCoolShutoffTime =
                  All.InputTimeHeatRadiationFeedback * SEC_PER_GIGAYEAR * All.HubbleParam / All.UnitTime_in_s; /* in code units */
              StarParticle[i].Lum = Lum;
              StarParticle[i].RadiationMomentumReleased =
                  Lum / CLIGHT * All.InputTimeMomRadiationFeedback * SEC_PER_GIGAYEAR; /* gr (cm/s) h^-1 ... */
              StarParticle[i].RadiationMomentumReleased /=
                  (All.UnitVelocity_in_cm_per_s *
                   All.UnitMass_in_g); /* code units, although missing cf_atime, added at use in radiation_stellar_feedback.c */
              StarParticle[i].RadiationMomentumReleased *= All.HubbleParam;
            }
          else
            {
              StarParticle[i].RadiationMomentumReleased = 0.;
              StarParticle[i].Lum                       = 0.;
              StarParticle[i].RadCoolShutoffTime        = 0.0;
              StarParticle[i].StromgrenMass             = -999.;
              StarParticle[i].StromgrenRadius           = -999.;
            }
#else  // Stars photoionize material at each timestep, for the duration of their timestep
          StarParticle[i].RadCoolShutoffTime =
              dtime_in_Gyr_rad * SEC_PER_GIGAYEAR * All.HubbleParam / All.UnitTime_in_s; /* in code units */
#endif

#else
          if(age_of_star_in_Gyr_rad >= All.InputTimeMomRadiationFeedback && STP(StarParticle[i].index).RadFeed_Flag == 0)
            {
              StarParticle[i].RadiationMomentumReleased =
                  Lum / CLIGHT * age_of_star_in_Gyr_rad * SEC_PER_GIGAYEAR; /* gr (cm/s) h^-1 ... */
              StarParticle[i].RadiationMomentumReleased /=
                  (All.UnitVelocity_in_cm_per_s *
                   All.UnitMass_in_g); /* code units, although missing cf_atime, added at use in radiation_stellar_feedback.c */
              StarParticle[i].RadiationMomentumReleased *= All.HubbleParam;

              STP(StarParticle[i].index).RadFeed_Flag = 0;
            }
          else
            StarParticle[i].RadiationMomentumReleased = 0.;
          StarParticle[i].RadCoolShutoffTime = (All.InputTimeHeatRadiationFeedback - age_of_star_in_Gyr_rad) * SEC_PER_GIGAYEAR *
                                               All.HubbleParam / All.UnitTime_in_s; /* in code units */
                                                                                    //#else
#endif
        }

#ifdef SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION
      STP(StarParticle[i].index).StromgrenRadius = StarParticle[i].StromgrenRadius;
      STP(StarParticle[i].index).StromgrenMass   = StarParticle[i].StromgrenMass;
#endif
#ifdef SMUGGLE_RADIATION_FEEDBACK_DEBUG
      STP(StarParticle[i].index).RadiationMomentumReleased =
          StarParticle[i]
              .RadiationMomentumReleased; /*LVS: this is lower limit, because limiter later in radiation_stellar_feedback.c. */
      STP(StarParticle[i].index).Cum_RadiationMomentumReleased +=
          StarParticle[i].RadiationMomentumReleased; /*LVS: this is lower limit, because limiter later in radiation_stellar_feedback.c.
The real injected is stored in StarParticle[i].Cum_RadMomentumRealInjected */
      STP(StarParticle[i].index).RadCoolShutoffTime = StarParticle[i].RadCoolShutoffTime;
      STP(StarParticle[i].index).RadFeedTau         = StarParticle[i].RadFeedTau;
      STP(StarParticle[i].index).GasColumnDensity   = StarParticle[i].GasColumnDensity;
#endif
#endif
    }

#ifdef VERBOSE
  double AGBGlobalMassReleased, SNIaGlobalMassReleased, SNIIGlobalMassReleased;

  MPI_Allreduce(&AGBLocalMassReleased, &AGBGlobalMassReleased, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&SNIaLocalMassReleased, &SNIaGlobalMassReleased, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&SNIILocalMassReleased, &SNIIGlobalMassReleased, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  All.AGBMassReleased += AGBGlobalMassReleased;
  All.SNIaMassReleased += SNIaGlobalMassReleased;
  All.SNIIMassReleased += SNIIGlobalMassReleased;
#endif

#ifdef GFM_RPROCESS_CHANNELS_NS_KICKS
  /* remove NSNS particles that merged */
  int count_removed = 0;
  for(int idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type != 4)
        continue;

      if(P[i].Mass == 0 && STP(i).NSNS_channel == -1)
        {
          timebin_remove_particle(&TimeBinsGravity, idx, P[i].TimeBinGrav);
          count_removed++;
        }
    }

  int count_removed_all;
  MPI_Reduce(&count_removed, &count_removed_all, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  mpi_printf("NSNS: Removed %d NSNS particles that merged and enriched.\n", count_removed_all);
#endif
}

#endif
