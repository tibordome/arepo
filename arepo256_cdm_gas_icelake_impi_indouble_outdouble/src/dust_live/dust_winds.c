/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/dust_winds.c
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../allvars.h"
#include "../proto.h"

#ifdef DUST_LIVE
#ifdef DL_WINDS

#ifndef GFM_WINDS
#error "Use of option DL_WINDS currently requires GFM_WINDS!"
#endif

void do_dust_winds(void)
{
  start_dust();
  find_drag_cells(Ndust);
  drag_kernel();

  for(int i = 0; i < Ndust; i++)
    {
      int p = DustParticle[i].index;

      double dt_base = (P[p].TimeBinHydro ? (((integertime)1) << P[p].TimeBinHydro) : 0) * All.Timebase_interval;
      double dt_step = dt_base / All.cf_hubble_a;

      check_dust_wind_recouple(p, dt_step);

      /* Don't double-kick dust particles in the wind. */
      if(DTP(p).IsWind)
        continue;

      /* Dust winds assume dM_dust,wind / M_dust = dM_gas,wind / M_gas. */
      double wmr = DustParticle[i].WindMassRate / (All.UnitMass_in_g / All.HubbleParam / SOLAR_MASS) / (SEC_PER_YEAR) *
                   (All.UnitTime_in_s / All.HubbleParam); /* internal mass / time */
      double mass_ratio = (wmr * dt_step) / (DustParticle[i].TotNgbMass);
      double prob       = mass_ratio;

      if(prob == 0.0)
        continue;

      if(prob < 0.0)
        terminate("prob < 0, Task=%d P[p].Mass=%g WindMassRate=%g TotNgbMass=%g dt_step=%g prob=%g\n", ThisTask, P[p].Mass,
                  DustParticle[i].WindMassRate, DustParticle[i].TotNgbMass, dt_step, prob);

      if(prob > 1.0)
        printf(
            "DUST_LIVE: Warning, want to create more dust wind mass than available dust mass.  Task=%d P[p].Mass=%g WindMassRate=%g "
            "TotNgbMass=%g dt_step=%g prob=%g\n",
            ThisTask, P[p].Mass, DustParticle[i].WindMassRate, DustParticle[i].TotNgbMass, dt_step, prob);

      double p_decide = get_random_number();

      if(p_decide < prob)
        {
          double v_wind = DustParticle[i].WindVel;
          dust_wind_kick(p, i, v_wind);

          DTP(p).IsWind       = 1;
          DTP(p).WindTimeLeft = All.WindFreeTravelMaxTimeFactor / All.cf_hubble_a;
        }
    }

  end_dust();
}

void dust_wind_kick(int p, int i, double v)
{
  double norm, dir[3];
  int j;
#ifndef GFM_BIPOLAR_WINDS
  MyFloat theta, phi;
#endif

#ifdef GFM_BIPOLAR_WINDS

#if(GFM_BIPOLAR_WINDS == 0)
  // TODO: Why doesn't this include GravPM?
  dir[0] = P[p].GravAccel[1] * P[p].Vel[2] - P[p].GravAccel[2] * P[p].Vel[1];
  dir[1] = P[p].GravAccel[2] * P[p].Vel[0] - P[p].GravAccel[0] * P[p].Vel[2];
  dir[2] = P[p].GravAccel[0] * P[p].Vel[1] - P[p].GravAccel[1] * P[p].Vel[0];
#else
#if(GFM_BIPOLAR_WINDS == 3)
  dir[0] = DustParticle[i].DensGasAngMomentum[0];
  dir[1] = DustParticle[i].DensGasAngMomentum[1];
  dir[2] = DustParticle[i].DensGasAngMomentum[2];
  if(dir[0] == 0 && dir[1] == 0 && dir[2] == 0)
    {
      double theta = acos(2 * get_random_number() - 1);
      double phi   = 2 * M_PI * get_random_number();
      dir[0]       = sin(theta) * cos(phi);
      dir[1]       = sin(theta) * sin(phi);
      dir[2]       = cos(theta);
    }
#else
  dir[0] = (P[p].GravAccel[1] - DustParticle[i].GroupGravAcc[1]) * (P[p].Vel[2] - DustParticle[i].GroupVel[2]) -
           (P[p].GravAccel[2] - DustParticle[i].GroupGravAcc[2]) * (P[p].Vel[1] - DustParticle[i].GroupVel[1]);
  dir[1] = (P[p].GravAccel[2] - DustParticle[i].GroupGravAcc[2]) * (P[p].Vel[0] - DustParticle[i].GroupVel[0]) -
           (P[p].GravAccel[0] - DustParticle[i].GroupGravAcc[0]) * (P[p].Vel[2] - DustParticle[i].GroupVel[2]);
  dir[2] = (P[p].GravAccel[0] - DustParticle[i].GroupGravAcc[0]) * (P[p].Vel[1] - DustParticle[i].GroupVel[1]) -
           (P[p].GravAccel[1] - DustParticle[i].GroupGravAcc[1]) * (P[p].Vel[0] - DustParticle[i].GroupVel[0]);
#endif
#endif
#else
  theta  = acos(2 * get_random_number() - 1);
  phi    = 2 * M_PI * get_random_number();
  dir[0] = sin(theta) * cos(phi);
  dir[1] = sin(theta) * sin(phi);
  dir[2] = cos(theta);
#endif

  for(j = 0, norm = 0; j < 3; j++)
    norm += dir[j] * dir[j];

  norm = sqrt(norm);
  if(get_random_number() < 0.5)
    norm = -norm;

  if(norm != 0)
    {
      for(j = 0; j < 3; j++)
        dir[j] /= norm;

      for(j = 0; j < 3; j++)
        {
          P[p].Vel[j] += (v * All.cf_atime * dir[j]);
        }
    }
}

void check_dust_wind_recouple(int p, double dt_step)
{
  if(DTP(p).IsWind == 0)
    return;

  DTP(p).WindTimeLeft -= dt_step;

  if((DTP(p).WindTimeLeft <= 0.0) || (DTP(p).LocalGasDensity * All.cf_a3inv < All.WindFreeTravelDensFac * All.PhysDensThresh))
    {
      DTP(p).IsWind       = 0;
      DTP(p).WindTimeLeft = 0.0;
    }
}

#endif
#endif
