/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sfr_mcs/sn_mcs.c
 * \date        08/2018
 * \author     	Matthew C Smith
 * \brief
 * \details     Private to M C Smith, but collaboration encouraged.
                Originally developed in 2015, ported into main repo 2018.
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

#if !defined(METALS) && !defined(SB99_FIXED_Z)
#error "If you are not using METALS, then Starburst99 requires a metallicity, use SB99_FIXED_Z to achieve this."
#endif

#if defined(SN_MCS_INITIAL_DRIVING) && !defined(DO_NOT_CREATE_STAR_PARTICLES)
#error "SN_MCS_INITIAL_DRIVING requires DO_NOT_CREATE_STAR_PARTICLES"
#endif

#if defined(SN_MCS_VARIABLE_EJECTA) && !defined(SB99_FIXED_Z)
#error "SN_MCS_VARIABLE_EJECTA has not yet been implemented with a metallicity dependence. You must use SB99_FIXED_Z for now."
#endif

#ifndef IMF_SAMPLING_MCS
#ifndef SN_MCS_INITIAL_DRIVING
void check_for_supernovae(int *local_n_sn_event, int *global_n_sn_event)
{
  int i, idx;
  int n_sn_event     = 0;
  int n_sn           = 0;
  int n_sn_tot       = 0;
  int n_sn_event_tot = 0;
  int n_sn_i;
  double dt, dtime, t_star, snr_global;
  double m_ej;
#ifndef SN_MCS_SINGLE_INJECTION
  double raw_sn_rate, norm_sn_rate, n_sn_expect;
#else
  double m_sample, prob;
#endif
#ifdef SN_MCS_DEBUG
  int n_sn_init;
#endif
#ifdef SN_MCS_VARIABLE_EJECTA
  double z_ej;
#endif

  mpi_printf("SN_MCS: Checking for SNe\n");

  /* To Do: Loop over star particle struct instead of all active particles */
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type != 4)
        continue;

      if(P[i].ID == 0 && P[i].Mass == 0)
        continue;

      if(All.ComovingIntegrationOn && (P[i].Mass > All.MaxFBStarMass))
        continue;

      t_star = STP(i).Age;

#ifndef SN_MCS_SINGLE_INJECTION
      if((t_star <= sb99.t_min_SN) || (t_star >= sb99.t_max_SN)) /* No SNe go off outside of this time range */
        continue;

      dt    = (P[i].TimeBinGrav ? (((integertime)1) << P[i].TimeBinGrav) : 0) * All.Timebase_interval;
      dtime = dt / All.cf_hubble_a;

#ifndef SB99_FIXED_Z
      raw_sn_rate = get_sn_rate(t_star, STP(i).iz); /* returns log10(SN Msun^-1 yr^-1) */
#else
      raw_sn_rate = get_sn_rate(t_star); /* returns log10(SN Msun^-1 yr^-1) */
#endif
      norm_sn_rate = pow(10.0, raw_sn_rate) *
                     (All.UnitTime_in_s * All.UnitMass_in_g / (All.HubbleParam * All.HubbleParam * SEC_PER_YEAR * SOLAR_MASS)) *
                     STP(i).InitialMass;

      n_sn_expect = norm_sn_rate * dtime;
      n_sn_i      = (int)gsl_ran_poisson(random_generator, n_sn_expect);

#if SN_MCS_CHANCE
      if(n_sn_i > 0)
        {
          if(get_random_number() > (1.0 / SN_MCS_CHANCE))
            n_sn_i = 0;
        }
#endif

#else  // SN_MCS_SINGLE_INJECTION
      if((t_star < All.SNDelay) || (STP(i).N_SN_event_cum > 0)) /* Particle too young or has already been sampled for SNe */
        continue;

      n_sn_i   = 0;
      m_sample = STP(i).InitialMass;

      while(m_sample > All.SNMassUnit)
        {
          n_sn_i++;
          m_sample -= All.SNMassUnit;
        }

      prob = 1.0 - m_sample / All.SNMassUnit;
      if(get_random_number() < prob)
        n_sn_i++;
#endif

#ifdef SN_MCS_DEBUG
      n_sn_init = n_sn_i;
#endif

      if(n_sn_i > 0)
        {
#ifdef SN_MCS_VARIABLE_EJECTA
          get_ejecta_properties(t_star, &m_ej, &z_ej);
          m_ej *= (SOLAR_MASS / All.UnitMass_in_g / All.HubbleParam);
          STP(i).Z_ej = z_ej;
#else
          m_ej = All.SNMassReturn * n_sn_i;
#endif
          if((P[i].Mass - m_ej) < 0.1 * STP(i).InitialMass)
            {
              STP(i).M_ej = P[i].Mass;
              P[i].Mass   = 0.0;
            }
          else
            {
              STP(i).M_ej = m_ej;
              P[i].Mass -= m_ej;
            }
          STP(i).N_SN = n_sn_i;
          STP(i).N_SN_cum += n_sn_i;

          assert(STP(i).M_ej >= 0);
          assert(P[i].Mass >= 0);

#ifdef SN_MCS_DEBUG
#ifndef METALS
          printf(
              "Task %d, Time %g, dtime = %g, ID = %d, t_star %g, raw_sn_rate = %g, InitialMass = %g, mass = %g, norm_sn_rate = %g, "
              "n_sn_expect = %g, n_sn_init = %d, n_sn_i = %d, N_SN = %d, mass returned = %g\n",
              ThisTask, All.Time, dtime, P[i].ID, t_star, raw_sn_rate, STP(i).InitialMass, P[i].Mass, norm_sn_rate, n_sn_expect,
              n_sn_init, n_sn_i, STP(i).N_SN, STP(i).M_ej);
#else
          printf("SN_MCS_DEBUG %d %g %g %d %g %g %g %g %g %g %g %d %d %d %d %g\n", ThisTask, All.Time, dtime, P[i].ID, t_star,
                 raw_sn_rate, STP(i).InitialMass, P[i].Mass, P[i].Metallicity, norm_sn_rate, n_sn_expect, n_sn_init, n_sn_i,
                 STP(i).N_SN, STP(i).N_SN_cum, STP(i).M_ej);
#endif
#endif
          n_sn_event++;
          STP(i).N_SN_event_cum += 1;
          n_sn += n_sn_i;
        }
      else
        {
          STP(i).N_SN = 0;
          STP(i).M_ej = 0;
#ifdef SN_MCS_SINGLE_INJECTION
          STP(i).N_SN_event_cum += 1;  // In this mode, we use this as a flag to check whether SNe sampling has occured yet.
#endif
        }

    } /* Main loop of active particles */

  MPI_Reduce(&n_sn, &n_sn_tot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); /*This is used for log, so only Task 0 needs this */
  MPI_Allreduce(&n_sn_event, &n_sn_event_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); /*Every task needs this to trigger SN routines */

  if(ThisTask == 0)
    {
      if(All.TimeStep > 0)
        snr_global = n_sn_tot / (All.TimeStep / All.cf_time_hubble_a);
      else
        snr_global = 0.0;
      fprintf(FdSnr, "%14e %14i %14i %14e\n", All.Time, n_sn_tot, n_sn_event_tot, snr_global);
      myflush(FdSnr);
    }

  *local_n_sn_event  = n_sn_event;
  *global_n_sn_event = n_sn_event_tot;
}
#else
int check_for_supernovae(void)
{
  int idx, i, n_sn, n_sn_global;
  double t_ff, dt, dtime, p, dx, dy, r;

  n_sn = 0;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type != 0)
        continue;

      if(P[i].ID == 0 && P[i].Mass == 0)
        continue;

      if((SphP[i].Density * All.cf_a3inv) < All.DensThreshold)
        continue;

      dx = P[i].Pos[0] - 0.5 * All.BoxSize;
      dy = P[i].Pos[1] - 0.5 * All.BoxSize;
      r  = sqrt(dx * dx + dy * dy);
      if(r > All.DrivingZoneRadius)
        continue;

      t_ff = sqrt(3.0 * M_PI / (32.0 * All.G * SphP[i].Density * All.cf_a3inv));

      dt    = (P[i].TimeBinGrav ? (((integertime)1) << P[i].TimeBinGrav) : 0) * All.Timebase_interval;
      dtime = dt / All.cf_hubble_a;

      p = All.SfrEfficiency * dtime / t_ff;                       // Like probability of forming star particle
      p *= P[i].Mass / (100.0 * SOLAR_MASS / All.UnitMass_in_g);  // 1 SN per 100 Msun
      if(get_random_number() < p)
        {
          SphP[i].mass_deposited   = All.SNMassReturn;  // n.b. this is only used for setting momentum injection
          SphP[i].energy_deposited = All.SupernovaEnergy * All.cf_atime * All.cf_atime;
          SphP[i].N_SN_hosted += 1;
          n_sn++;
        }
    }

  MPI_Allreduce(&n_sn, &n_sn_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); /*Every task needs this to trigger SN routines */

  return n_sn_global;
}
#endif
#endif  // IMF_SAMPLING

/*Tidy up, reset things for next timestep, largely for safety */
void sn_finish(void)
{
  int i, idx;

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;
      if(P[i].Type != 4)
        continue;
      if(P[i].ID == 0 && P[i].Mass == 0)
        continue;

      /* Prepare 'dead' star particles for deletion */
      if(P[i].Mass == 0)
        {
          P[i].ID     = 0;
          P[i].Vel[0] = 0;
          P[i].Vel[1] = 0;
          P[i].Vel[2] = 0;
          timebin_remove_particle(&TimeBinsGravity, idx, P[i].TimeBinGrav);
        }

      STP(i).N_SN = 0; /* Redundant, these should be taken care of elsewhere, but kept just in case */
      STP(i).M_ej = 0.0;
#ifdef IMF_SAMPLING_MCS
      STP(i).Z_ej = 0.0;
#endif
    }
}

void do_supernovae(void)
{
  TIMER_START(CPU_SN_FEEDBACK);

#ifndef SN_MCS_INITIAL_DRIVING
  int n_sn_event, n_sn_event_global;

#ifndef IMF_SAMPLING_MCS
  check_for_supernovae(&n_sn_event, &n_sn_event_global);
#else
  n_sn_event        = NumSNLocal;
  n_sn_event_global = NumSNGlobal;
#endif

  if(n_sn_event_global > 0) /*i.e. some task somewhere has at least one SN */
    {
      find_sn_host_cells_and_distribute(n_sn_event);
      stellar_feedback_distribute();
    }
#else
  int n_sn_event_global;

  n_sn_event_global = check_for_supernovae();
  if(n_sn_event_global > 0)
    stellar_feedback_distribute();
#endif

  sn_finish(); /*Tidy up */

  TIMER_STOP(CPU_SN_FEEDBACK);
}

#endif
