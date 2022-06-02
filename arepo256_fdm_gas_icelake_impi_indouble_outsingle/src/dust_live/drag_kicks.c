/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/drag_kicks.c
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
/* For particle index p and bin index i, calculate a quantity (in units of
 * um^2) at grain size a used in determining the effective stopping time-scale.
 * */
double t_s_eff_helper(int p, int i, double a)
{
  double fac1 = DTP(p).NumGrains[i] * a * a * a / 3.0 / GSD.Widths[i];
  double fac2 = 0.0;
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
  fac2 = DTP(p).BinSlopes[i] * (pow(a, 4) / 4.0 - GSD.Midpoints[i] * a * a * a / 3.0);
#endif
  return fac1 + fac2;
}
#endif

void compute_drag_acceleration(void)
{
  for(int i = 0; i < Ndust; i++)
    {
      int p = DustParticle[i].index;
      if(P[p].Mass == 0.0)
        continue;

        /* Calculate quantities and stopping time-scale in physical units.
         * Will convert to comoving when drag acceleration is applied. */
#ifdef DL_GRAIN_BINS
      /* If we're evolving grain size distributions, we need to calculate an
       * effective stopping time-scale. */
      double rho_grain = All.GrainDensity / (All.HubbleParam * All.HubbleParam * All.UnitDensity_in_cgs);
      double a2        = 0.0;
      for(int k = 0; k < DL_GRAIN_BINS; k++)
        {
          a2 += (t_s_eff_helper(p, k, GSD.Edges[k + 1]) - t_s_eff_helper(p, k, GSD.Edges[k]));
        }
      if(a2 <= 0.0)
        {
          /* Assume grains in smallest bin. */
          a2 = (P[p].Mass / GSD.AvgMasses[0]) * (pow(GSD.Edges[1], 3) - pow(GSD.Edges[0], 3)) / (3.0 * GSD.Widths[0]);
        }

      double a_grain = 3.0 * P[p].Mass / (4.0 * M_PI * GSD.InternalDensity * a2); /* um */
      if(a_grain > All.MaxGrainSize)
        a_grain = All.MaxGrainSize;
      if(a_grain < All.MinGrainSize)
        a_grain = All.MinGrainSize;
      a_grain *= 1.0e-4 * All.HubbleParam / All.UnitLength_in_cm; /* internal length, since 1 um = 1.0e-4 cm */
#else
      /* If we're not evolving grain size distributions, just fix the grain
       * density and radius. */
      double rho_grain = 2.4 / (All.HubbleParam * All.HubbleParam * All.UnitDensity_in_cgs);
      double a_grain   = 1.0e-5 * All.HubbleParam / All.UnitLength_in_cm;
#endif
      double rho_gas = DTP(p).LocalGasDensity * All.cf_a3inv;
      double c_s     = DTP(p).LocalSoundSpeed;
      /* Hopkins+ 2016, Equation 2 */
      double t_s = sqrt(M_PI * GAMMA) / sqrt(8.0) * (rho_grain * a_grain) / (rho_gas * c_s);
#ifdef DL_STOPPING_TIME_CORRECTION
      double v2 = 0.0;
      for(int j = 0; j < 3; j++)
        {
          double dvel = (DTP(p).LocalGasVelocity[j] - P[p].Vel[j]) / All.cf_atime;
          v2 += (dvel * dvel);
        }
      double corr = pow(1.0 + 9.0 * M_PI / 128.0 * v2 / (c_s * c_s), -0.5);
      t_s *= corr;
#endif

      DTP(p).StoppingTime = t_s;
      for(int j = 0; j < 3; j++)
        {
          /* Compute physical velocity difference. */
          double dvel         = -(P[p].Vel[j] - DTP(p).LocalGasVelocity[j]) / All.cf_atime;
          DTP(p).DragAccel[j] = dvel / t_s;
        }
    }
}

void do_drag_step(void)
{
  int i, j;
  double dt_base, dt_kick;

  start_dust();
  find_drag_cells(Ndust);
  drag_kernel();
  TIMER_STOPSTART(CPU_DUST, CPU_DUST_DRAG);

  compute_drag_acceleration();

  for(i = 0; i < Ndust; i++)
    {
      int p = DustParticle[i].index;
      if(P[p].Mass == 0.0)
        continue;

#ifdef DL_WINDS
      if(DTP(p).IsWind)
        continue;
#endif

      dt_base = 0.5 * (P[p].TimeBinHydro ? (((integertime)1) << P[p].TimeBinHydro) : 0) * All.Timebase_interval;

      dt_kick = dt_base / All.cf_hubble_a;

      /* Convert drag acceleration back to comoving to apply kicks. */
#ifdef DL_DRAG_SEMI_IMPLICIT
      /* Loren-Aguilar+ 2015, Equation 17 */
#ifdef DL_DRAG_BACKREACTION
      double DGR = DTP(p).DustDensity / DTP(p).LocalGasDensity;
#else
      /* In the test-particle limit, there is no correction to drag from the
       * dust density. */
      double DGR = 0.0;
#endif
      double t_s = DTP(p).StoppingTime / (1.0 + DGR);
      double xi  = (1.0 - exp(-dt_kick / t_s)) / (1.0 + DGR);
      double lam = (dt_kick + t_s) * xi - dt_kick / (1.0 + DGR);
      for(j = 0; j < 3; j++)
        {
          double a_dust = P[p].GravAccel[j];
#ifdef PMGRID
          a_dust += P[p].GravPM[j];
#endif
          /* Physical pressure gradient picks up a^-3 for pressure and a^-1 for
           * gradient. */
          double a_DGj = ((a_dust - DustParticle[i].LocalGasAccel[j]) * (All.cf_a2inv) +
                          (DustParticle[i].LocalGradP[j] * All.cf_a3inv / All.cf_atime) / (DTP(p).LocalGasDensity * All.cf_a3inv));
          double v_DGj = (P[p].Vel[j] - DTP(p).LocalGasVelocity[j]) / All.cf_atime;
          /* Calculate a physical change in velocity. */
          double dvel = lam * a_DGj - xi * v_DGj;
#ifdef DL_DRAG_BACKREACTION
          /* Calculate negative change in dust kinetic energy that is
           * transferred to gas as thermal energy.  See below also. */
          DustParticle[i].DragDustThermal += 0.5 * P[p].Mass * P[p].Vel[j] * P[p].Vel[j];
#endif
          /* Could compute an effective physical acceleration, but we would
           * convert right back to comoving below. */
          P[p].Vel[j] += All.cf_atime * dvel;
#ifdef DL_DRAG_BACKREACTION
          DustParticle[i].DragMomentum[j] += P[p].Mass * All.cf_atime * dvel;
          DustParticle[i].DragDustThermal -= 0.5 * P[p].Mass * P[p].Vel[j] * P[p].Vel[j];
#endif
        }
#else
      /* Simple explicit update, since dust timesteps have been limited
       * to account for stopping time. */
      for(j = 0; j < 3; j++)
        {
          /* DragAccel was computed as a physical acceleration earlier. */
          double dvel = (dt_kick * DTP(p).DragAccel[j]);
#ifdef DL_DRAG_BACKREACTION
          DustParticle[i].DragDustThermal += 0.5 * P[p].Mass * P[p].Vel[j] * P[p].Vel[j];
#endif
          P[p].Vel[j] += All.cf_atime * dvel;
#ifdef DL_DRAG_BACKREACTION
          DustParticle[i].DragMomentum[j] = P[p].Mass * All.cf_atime * dvel;
          DustParticle[i].DragDustThermal -= 0.5 * P[p].Mass * P[p].Vel[j] * P[p].Vel[j];
#endif
        }
#endif
    }

#ifdef DL_DRAG_BACKREACTION
  drag_backreaction_kernel();
  update_primitive_variables();
#endif

  TIMER_STOPSTART(CPU_DUST_DRAG, CPU_DUST);
  end_dust();
}

void do_drag_step_first_half(void) { do_drag_step(); }

void do_drag_step_second_half(void) { do_drag_step(); }

#ifdef DL_DRAG_BACKREACTION
void do_drag_backreaction_search(void)
{
  // TODO: can eliminate redundancy with shattering search
  begin_dust_search();
  dust_findHsml();
  exchange_dust_search_results();
  end_dust_search();
}
#endif

#endif
