/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/dust_transfer.c
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
#ifdef DL_GRAIN_BINS

/* Rescale the grain size distribution by a small amount if necessary to
 * ensure integrated mass matches the particle mass. */
void rescale_gsd_mass(int i)
{
  double mass = 0.0;
  int p       = DustParticle[i].index;
  for(int j = 0; j < DL_GRAIN_BINS; j++)
    {
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
      mass += DTP(p).NumGrains[j] * bin_avg_mass(j, DTP(p).NumGrains[j], DTP(p).BinSlopes[j]);
#else
      mass += DTP(p).NumGrains[j] * GSD.AvgMasses[j];
#endif
    }
  double mass_ratio = (mass > 0.0) ? P[p].Mass / mass : 0.0;
  /* Just scale the grain size distribution if possible. */
  for(int j = 0; j < DL_GRAIN_BINS; j++)
    {
      DTP(p).NumGrains[j] *= mass_ratio;
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
      DTP(p).BinSlopes[j] *= mass_ratio;
#endif
    }
}

#if defined(DL_GROWTH) || defined(DL_SPUTTERING)
void do_dust_transfer(void)
{
  int i;
  double dt_base, dt_step;

  start_dust();
  find_drag_cells(Ndust);
  drag_kernel();
  TIMER_STOPSTART(CPU_DUST, CPU_DUST_TRANSFER);

  for(i = 0; i < Ndust; i++)
    {
      int p = DustParticle[i].index;

      dt_base = (P[p].TimeBinHydro ? (((integertime)1) << P[p].TimeBinHydro) : 0) * All.Timebase_interval;

      dt_step = dt_base / All.cf_hubble_a;

      /* Calculate adot from various physical processes.  Keep adot in units
       * of um/Gyr during each update. */
      double adot = 0.0;
#ifdef DL_GROWTH
      /* Hirashita & Voshchinnikov (2014), Equations 5 to 8 */
      double rho_H_ref   = (1.0e3 * PROTONMASS) / (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam);
      double T_ref       = 10.0;
      double Z_ref       = 0.0127;
      double S_ref       = 0.3;
      double S_local     = S_ref;
      double rho_H_local = DustParticle[i].LocalGasDensityH * All.cf_a3inv;
      double T_local     = DustParticle[i].LocalGasTemp;
      double Z_local     = DustParticle[i].LocalGasZ;

      double adot_growth = (Z_local / Z_ref) * (rho_H_local / rho_H_ref) * sqrt(T_local / T_ref) * (S_local / S_ref);
      adot += adot_growth;
#endif
#ifdef DL_SPUTTERING
      /* Tsai & Matthews (1995), Equation 14, converted to um/Gyr */
      double rho_local_cgs = DTP(p).LocalGasDensity * All.cf_a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
      double T_fac         = 1.0 + pow(2.0e6 / DustParticle[i].LocalGasTemp, 2.5);
      double h_fac         = 1009.28; /* cm^3 um / Gyr */
      double adot_sputter  = -h_fac * (rho_local_cgs / PROTONMASS) / T_fac;
      adot += adot_sputter;
#endif

      double t_sub  = 0.0;
      double dt_gsd = DTP(p).OldBinMassChgTau;

      /* Keep track of how much mass is transferred during the possible
       * subcycling.  We need to store this quantity for later in case there
       * aren't enough metals in surrounding gas. */
      DustParticle[i].DeltaMassExpected = 0.0;
      /* Store the initial grain size distribution, so that we can update the
       * grain size distribution in multiple subcycling timesteps if necessary.
       * This way, the final grain size distribution from one timestep can
       * become the new grain size distribution for a subsequent timestep.  We
       * need to know the original grain size distribution for
       * correct_grain_conserved_quantities() after kernel operations. */
      MyFloat old_num[DL_GRAIN_BINS];
      memcpy(old_num, DTP(p).NumGrains, DL_GRAIN_BINS * sizeof(MyFloat));
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
      MyFloat old_slope[DL_GRAIN_BINS];
      memcpy(old_slope, DTP(p).BinSlopes, DL_GRAIN_BINS * sizeof(MyFloat));
#endif

      while(t_sub < dt_step)
        {
          double dt_sub = fmin(All.MaxBinFracMassChg * dt_gsd, dt_step - t_sub);
          /* Since [a] = um and [adot] = um / Gyr, need timestep in Gyr. */
          double dt_sub_Gyr = dt_sub * All.UnitTime_in_s / All.HubbleParam / SEC_PER_GIGAYEAR;

          dt_gsd = update_grain_sizes(i, adot, dt_sub_Gyr, dt_sub);
          t_sub += dt_sub;
        }

      /* Restore DTP variables to their initial state to prepare to finish the
       * transfer in the kernel update. */
      memcpy(DTP(p).NumGrains, old_num, DL_GRAIN_BINS * sizeof(MyFloat));
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
      memcpy(DTP(p).BinSlopes, old_slope, DL_GRAIN_BINS * sizeof(MyFloat));
#endif
    }

  /* Attempt to remove or add metal mass from or to surrounding gas cells. */
  dust_transfer_kernel();

  /* Correct the grain size distribution update in case there weren't enough
   * metals in some gas cells. */
  for(i = 0; i < Ndust; i++)
    {
      update_dust_element_fractions(i);
      correct_grain_conserved_quantities(i);
      check_dust_for_removal(i);
      rescale_gsd_mass(i);
    }

  TIMER_STOPSTART(CPU_DUST_TRANSFER, CPU_DUST);
  end_dust();
}
#endif

#ifdef DL_SNE_DESTRUCTION
void do_dust_transfer_sne(void)
{
  int i;
  double dt_base, dt_step;

  start_dust();
  find_drag_cells(Ndust);
  drag_kernel();
  TIMER_STOPSTART(CPU_DUST, CPU_DUST_TRANSFER);

  for(i = 0; i < Ndust; i++)
    {
      int p = DustParticle[i].index;

      dt_base = (P[p].TimeBinHydro ? (((integertime)1) << P[p].TimeBinHydro) : 0) * All.Timebase_interval;

      dt_step = dt_base / All.cf_hubble_a;

      double t_sub  = 0.0;
      double dt_gsd = DTP(p).OldBinMassChgTau;

      DustParticle[i].DeltaMassExpected = 0.0;
      MyFloat old_num[DL_GRAIN_BINS];
      memcpy(old_num, DTP(p).NumGrains, DL_GRAIN_BINS * sizeof(MyFloat));
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
      MyFloat old_slope[DL_GRAIN_BINS];
      memcpy(old_slope, DTP(p).BinSlopes, DL_GRAIN_BINS * sizeof(MyFloat));
#endif

      while(t_sub < dt_step)
        {
          double dt_sub = fmin(All.MaxBinFracMassChg * dt_gsd, dt_step - t_sub);
          /* Since [a] = um and [adot] = um / Gyr, need timestep in Gyr. */
          double dt_sub_Gyr = dt_sub * All.UnitTime_in_s / All.HubbleParam / SEC_PER_GIGAYEAR;

          /* Supernova destruction happens in a subgrid fashion, not by specifying
           * a form of da/dt as for growth or sputtering. */
          dt_gsd = update_grain_sizes_sne(i, dt_sub_Gyr, dt_sub);
          t_sub += dt_sub;
        }

      /* Restore DTP variables to their initial state to prepare to finish the
       * transfer in the kernel update. */
      memcpy(DTP(p).NumGrains, old_num, DL_GRAIN_BINS * sizeof(MyFloat));
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
      memcpy(DTP(p).BinSlopes, old_slope, DL_GRAIN_BINS * sizeof(MyFloat));
#endif
    }

  /* Attempt to remove or add metal mass from or to surrounding gas cells. */
  dust_transfer_kernel();

  /* Correct the grain size distribution update in case there weren't enough
   * metals in some gas cells. */
  for(i = 0; i < Ndust; i++)
    {
      update_dust_element_fractions(i);
      correct_grain_conserved_quantities(i);
      check_dust_for_removal(i);
      rescale_gsd_mass(i);
    }

  TIMER_STOPSTART(CPU_DUST_TRANSFER, CPU_DUST);
  end_dust();
}
#endif

#if defined(DL_GROWTH) || defined(DL_SPUTTERING) || defined(DL_SNE_DESTRUCTION)
void update_dust_element_fractions(int i)
{
  int p = DustParticle[i].index;

  double elem_masses[GFM_N_CHEM_ELEMENTS];
  for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
    {
      elem_masses[k] = P[p].Mass * DTP(p).MetalFractions[k] + DustParticle[i].DeltaMetalMasses[k];
      /* Reset to zero in the event of a tiny negative mass remaining from all
       * of a species' mass being lost and floating point arithmetic.  This
       * does not affect the particle mass, just the normalized dust metal
       * fractions. */
      if(elem_masses[k] < 0.0)
        elem_masses[k] = 0.0;
    }

  double sum = 0.0;
  for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
    {
      sum += elem_masses[k];
    }
  for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
    {
      DTP(p).MetalFractions[k] = (sum > 0.0 ? elem_masses[k] / sum : 0.0);
    }
}

void check_dust_for_removal(int i)
{
  int p = DustParticle[i].index;

#ifdef DL_REFINEMENT
  /* If refinement defines a natural maximum dust mass, remove dust particles
   * when their mass is sufficiently small relative to that value. */
  if((P[p].Mass <= 0.0) || (P[p].Mass <= 1.0e-15 * All.DustMaxFrac * All.TargetGasMass))
#else
  if(P[p].Mass <= 0.0)
#endif
    {
      P[p].Mass = 0.0;
      /* If we want to remove a dust particle, its idx in TimeBinsDust (stored
       * when filling in DustParticle) is not the same as its idx in
       * TimeBinsGravity.  Manually search TimeBinsGravity for the idx of this
       * dust particle.  This should not be too inefficient because the number
       * of removed dust particles is not huge. */
      int active_grav_idx = -1;
      for(int idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
        {
          int i_grav = TimeBinsGravity.ActiveParticleList[idx];
          if((i_grav >= 0) && (i_grav == p))
            {
              active_grav_idx = idx;
              break;
            }
        }
      if(active_grav_idx == -1)
        terminate("Could not find active dust particle in the gravity time bins!");

      timebin_remove_particle(&TimeBinsDust, DustParticle[i].active_idx, P[p].TimeBinHydro);
      timebin_remove_particle(&TimeBinsGravity, active_grav_idx, P[p].TimeBinGrav);
    }
}
#endif

#if defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION)
void do_shattering_coagulation(void)
{
  start_dust();
  find_drag_cells(Ndust);
  drag_kernel();
  TIMER_STOPSTART(CPU_DUST, CPU_DUST_SHATTER);

  /* Prepare to do dust-dust neighbor searches. */
  begin_dust_search();
  /* Do the dust-dust neighbor searches to get dust densities. */
  dust_findHsml();
  /* May need to communicate some results depending on the tree structure. */
  exchange_dust_search_results();

#ifdef DL_DEREFINEMENT
  // TODO: probably okay to reuse densities in shattering/coagulation calculations
  // TODO: do I need to worry about P.Mass, DTP.MetalFractions changing here but not in DustParticle? I don't think so
  derefine_dust_particles();
#endif

  /* The above calculations obtained dust densities that may be used even when
   * shattering and coagulation are not active (e.g. in supernova destruction).
   * We only need to proceed if shattering and coagulation take place. */

#if defined(DL_SHATTERING) || defined(DL_COAGULATION)
  /* Using dust densities, actually update grain size distributions. */
  /* Entirely local, no mass transfer to gas.  DustP also contains recent
   * CloudFrac estimates. */
  int i;
  double dt_base, dt_step;

  for(i = 0; i < Ndust; i++)
    {
      int p = DustParticle[i].index;

      /* If DL_DEREFINEMENT is on, it's possible that an entry in DustParticle
       * was just derefined and requires us to manually skip it. */
      if(P[p].Mass == 0.0)
        continue;

      dt_base = (P[p].TimeBinHydro ? (((integertime)1) << P[p].TimeBinHydro) : 0) * All.Timebase_interval;

      dt_step = dt_base / All.cf_hubble_a;

      double t_sub  = 0.0;
      double dt_gsd = DTP(p).OldBinMassChgTau;

      while(t_sub < dt_step)
        {
          double dt_sub = fmin(All.MaxBinFracMassChg * dt_gsd, dt_step - t_sub);
          /* Use timestep in Gyr, with mass loss rates calculated as mass per Gyr. */
          double dt_sub_Gyr = dt_sub * All.UnitTime_in_s / All.HubbleParam / SEC_PER_GIGAYEAR;

          double dt_shat = MAX_REAL_NUMBER, dt_coag = MAX_REAL_NUMBER;
#ifdef DL_SHATTERING
          dt_shat = update_grain_sizes_shattering(i, dt_sub_Gyr, dt_sub, GSD_SHATTERING);
#endif
#ifdef DL_COAGULATION
          dt_coag = update_grain_sizes_shattering(i, dt_sub_Gyr, dt_sub, GSD_COAGULATION);
#endif
          dt_gsd = fmin(dt_shat, dt_coag);
          t_sub += dt_sub;
        }
    }
#endif

  /* Free temporary memory. */
  end_dust_search();

  TIMER_STOPSTART(CPU_DUST_SHATTER, CPU_DUST);
  end_dust();
}
#endif

#ifdef DL_SNE_DESTRUCTION
void update_sn_rates(void)
{
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type == 0)
        {
          if(P[i].Mass == 0 && P[i].ID == 0)
            continue;

          double dt_base = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
          double dt_step;

          dt_step = dt_base / All.cf_hubble_a;

          /* Use timestep in Gyr, with mass loss rates calculated as mass per Gyr. */
          double dt_Gyr = dt_step * All.UnitTime_in_s / All.HubbleParam / SEC_PER_GIGAYEAR;
          if(dt_Gyr > 0.0)
            SphP[i].SNRate = SphP[i].Sfr * 1.0e9 * GSD.SNIINumFrac;
          else
            SphP[i].SNRate = 0.0;

          SphP[i].NumSNII = 0.0;
        }
    }
}
#endif

#endif
#endif
