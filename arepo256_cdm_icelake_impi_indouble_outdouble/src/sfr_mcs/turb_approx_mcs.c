/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sfr_mcs/turb_approx_mcs.c
 * \date        07/2019
 * \author
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - 17.01.2022 Ported into current codebase
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef TURB_APPROX_MCS

static double vel_gradient_cross_term(int p, int indi, int indj)
{
  MySingle **dvel;
#ifdef TURB_APPROX_MCS_GRAD_UNLIM
  dvel = SphP[p].dvel_unlim;
#else
  dvel = SphP[p].Grad.dvel;
#endif
  return dvel[indi][0] * dvel[indj][0] + dvel[indi][1] * dvel[indj][1] + dvel[indi][2] * dvel[indj][2];
}

void turb_approx_init(void)
{
  for(int i = 0; i < NumPart; i++)
    {
      if(P[i].Type != 0)
        continue;
      if(RestartFlag == 0)
        {
#ifdef TURB_APPROX_MCS_SEED
          MySingle **dvel;
#ifdef TURB_APPROX_MCS_GRAD_UNLIM
          dvel = SphP[i].dvel_unlim;
#else
          dvel = SphP[i].Grad.dvel;
#endif
          double divv = dvel[0][0] + dvel[1][1] + dvel[2][2];
          /* Shear tensor */
          /* Diagonal terms */
          double sfree_xx = dvel[0][0] - divv / 3.0;
          double sfree_yy = dvel[1][1] - divv / 3.0;
          double sfree_zz = dvel[2][2] - divv / 3.0;
          /* Upper triangle */
          double sfree_xy = 0.5 * (dvel[0][1] + dvel[1][0]);
          double sfree_xz = 0.5 * (dvel[0][2] + dvel[2][0]);
          double sfree_yz = 0.5 * (dvel[1][2] + dvel[2][1]);

          double dx = pow(6.0 * SphP[i].Volume / M_PI, 1.0 / 3.0);

          double s2 = sfree_xx * sfree_xx + sfree_yy * sfree_yy + sfree_zz * sfree_zz + 2.0 * sfree_xy * sfree_xy +
                      2.0 * sfree_xz * sfree_xz + 2.0 * sfree_yz * sfree_yz;

          SphP[i].TurbSpecEnergy = 0.5 * dx * dx * s2;
#elif !defined(TURB_APPROX_MCS_KEEP_ICS)
          SphP[i].TurbSpecEnergy = 0.0;
#endif
        }
#ifdef NODEREFINE_BACKGROUND_GRID
      if(SphP[i].Volume > 0.1 * All.MeanVolume)
        SphP[i].TurbSpecEnergy = 0.0;
#endif

      SphP[i].TurbSpecEnergy = fmax(SphP[i].TurbSpecEnergy, All.MinTurbSpecEnergy);
      SphP[i].TurbEnergy     = P[i].Mass * SphP[i].TurbSpecEnergy;

#ifdef TURB_APPROX_MCS_RENORM
      SphP[i].OldVolume = SphP[i].Volume;
#endif
    }
}

/*! \brief Rescale, produce and dissipate subgrid turbulent energy in each cell
 *
 *  This function is called twice each timestep, integrating by half a dt.
 *
 *  \return void
 */
void update_turbulent_energy(void)
{
  TIMER_START(CPU_TURB_APPROX);

  double kturb, dx, rate, dt, dtime;
  double divv, gradnorm2;
  double s_xx, s_yy, s_zz, s_xy, s_xz, s_yz;             /* Shear tensor */
  double sfree_xx, sfree_yy, sfree_zz;                   /* Diagonal components of trace free shear tensor */
  double tau_xx, tau_yy, tau_zz, tau_xy, tau_xz, tau_yz; /* closure tensor */
  MySingle **dvel;

  /* Closure constants from Schmidt+2006 */
  const double c1   = 0.02;
  const double c2   = 0.7;
  const double c3   = ((2.0 / 3.0) * (1.0 - c2));
  const double ceps = 1.58;  // dissipation constant, Schmidt+2014

  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];

      if(i < 0)
        continue;

      if(P[i].Mass == 0 && P[i].ID == 0)
        continue;

#ifdef NODEREFINE_BACKGROUND_GRID
      if(SphP[i].Volume > 0.1 * All.MeanVolume)
        continue;
#endif

#ifdef TURB_APPROX_MCS_RENORM
      /* First, rescale energy to new cell diameter */
      double csnd2 = GAMMA_MINUS1 * SphP[i].Pressure / SphP[i].Density;
      if(csnd2 > (2.0 * SphP[i].TurbSpecEnergy / 3.0))
        SphP[i].TurbSpecEnergy *= pow(SphP[i].Volume / SphP[i].OldVolume, 2. / 9.);  // Kolmogorov volume scaling
      else
        SphP[i].TurbSpecEnergy *= pow(SphP[i].Volume / SphP[i].OldVolume, 1. / 3.);  // Burgers volume scaling

      SphP[i].TurbSpecEnergy = fmax(SphP[i].TurbSpecEnergy, All.MinTurbSpecEnergy);
      SphP[i].OldVolume      = SphP[i].Volume;  // Reset for next pass
#endif

      SphP[i].TurbEnergy = SphP[i].TurbSpecEnergy * P[i].Mass;

      /*Relevant quantities */
      dx    = pow(6.0 * SphP[i].Volume / M_PI, 1.0 / 3.0);
      kturb = SphP[i].TurbEnergy / SphP[i].Volume;  // This is the K that appears in most SGS models

#ifdef TURB_APPROX_MCS_GRAD_UNLIM
      dvel = SphP[i].dvel_unlim;
#else
      dvel = SphP[i].Grad.dvel;
#endif

      /* Get relevant velocity gradient quantities */
      /* Divergence of velocity field */
      divv = dvel[0][0] + dvel[1][1] + dvel[2][2];

      /* Squared norm of the velocity gradient tensor */
      gradnorm2 = 0.0;
      for(int j = 0; j < 3; j++)
        {
          for(int k = 0; k < 3; k++)
            gradnorm2 += dvel[j][k] * dvel[j][k];
        }

      /* Shear tensor */
      /* Diagonal terms */
      s_xx = dvel[0][0];
      s_yy = dvel[1][1];
      s_zz = dvel[2][2];
      /* Upper triangle */
      s_xy = 0.5 * (dvel[0][1] + dvel[1][0]);
      s_xz = 0.5 * (dvel[0][2] + dvel[2][0]);
      s_yz = 0.5 * (dvel[1][2] + dvel[2][1]);

      /* Trace free shear tensor */
      /* Diagonal terms */
      sfree_xx = s_xx - divv / 3.0;
      sfree_yy = s_yy - divv / 3.0;
      sfree_zz = s_zz - divv / 3.0;

      /* Construct the closure */
      /* First the diagonal parts */
      tau_xx = 2.0 * c1 * dx * sqrt(2.0 * SphP[i].Density * kturb) * sfree_xx -
               2.0 * c2 * kturb * vel_gradient_cross_term(i, 0, 0) / gradnorm2 - c3 * kturb;
      tau_yy = 2.0 * c1 * dx * sqrt(2.0 * SphP[i].Density * kturb) * sfree_yy -
               2.0 * c2 * kturb * vel_gradient_cross_term(i, 1, 1) / gradnorm2 - c3 * kturb;
      tau_zz = 2.0 * c1 * dx * sqrt(2.0 * SphP[i].Density * kturb) * sfree_zz -
               2.0 * c2 * kturb * vel_gradient_cross_term(i, 2, 2) / gradnorm2 - c3 * kturb;
      /* Upper triangle */
      tau_xy =
          2.0 * c1 * dx * sqrt(2.0 * SphP[i].Density * kturb) * s_xy - 2.0 * c2 * kturb * vel_gradient_cross_term(i, 0, 1) / gradnorm2;
      tau_xz =
          2.0 * c1 * dx * sqrt(2.0 * SphP[i].Density * kturb) * s_xz - 2.0 * c2 * kturb * vel_gradient_cross_term(i, 0, 2) / gradnorm2;
      tau_yz =
          2.0 * c1 * dx * sqrt(2.0 * SphP[i].Density * kturb) * s_yz - 2.0 * c2 * kturb * vel_gradient_cross_term(i, 1, 2) / gradnorm2;

      /* Now get dkturb/dt */
      rate = 0;
      /* Start with source term */
      rate += tau_xx * s_xx + tau_yy * s_yy + tau_zz * s_zz + 2.0 * tau_xy * s_xy + 2.0 * tau_xz * s_xz + 2.0 * tau_yz * s_yz;

      /* Now adiabatic work term, pressure_turb * divv */
      rate -= 2.0 * kturb * divv / 3.0;

      /* Now dissipation term */
      rate -= ceps * sqrt(kturb / SphP[i].Density) * kturb / dx;

      /* Now turn into dE/dt */
      rate *= SphP[i].Volume;

      /*Note half timestep*/
      dt = 0.5 * (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
      /*  the actual time-step */
      dtime = All.cf_atime * dt / All.cf_time_hubble_a;

      SphP[i].TurbEnergy += rate * dtime;
      SphP[i].TurbSpecEnergy = SphP[i].TurbEnergy / P[i].Mass;
      if(SphP[i].TurbSpecEnergy < All.MinTurbSpecEnergy)
        {
          SphP[i].TurbSpecEnergy = All.MinTurbSpecEnergy;
          SphP[i].TurbEnergy     = SphP[i].TurbSpecEnergy * P[i].Mass;
        }
    }

  TIMER_STOP(CPU_TURB_APPROX);
}

#ifdef TURB_APPROX_MCS_OUTPUT_RATES

double get_turbulent_production_rate_single(int i)
{
  double kturb, dx, rate;
  double divv, gradnorm2;
  double s_xx, s_yy, s_zz, s_xy, s_xz, s_yz;             /* Shear tensor */
  double sfree_xx, sfree_yy, sfree_zz;                   /* Diagonal components of trace free shear tensor */
  double tau_xx, tau_yy, tau_zz, tau_xy, tau_xz, tau_yz; /* closure tensor */
  MySingle **dvel;

  /* Closure constants from Schmidt+2006 */
  const double c1 = 0.02;
  const double c2 = 0.7;
  const double c3 = ((2.0 / 3.0) * (1.0 - c2));

  dx    = pow(6.0 * SphP[i].Volume / M_PI, 1.0 / 3.0);
  kturb = SphP[i].TurbEnergy / SphP[i].Volume;  // This is the K that appears in most SGS models

#ifdef TURB_APPROX_MCS_GRAD_UNLIM
  dvel = SphP[i].dvel_unlim;
#else
  dvel = SphP[i].Grad.dvel;
#endif

  /* Get relevant velocity gradient quantities */
  /* Divergence of velocity field */
  divv = dvel[0][0] + dvel[1][1] + dvel[2][2];

  /* Squared norm of the velocity gradient tensor */
  gradnorm2 = 0.0;
  for(int j = 0; j < 3; j++)
    {
      for(int k = 0; k < 3; k++)
        gradnorm2 += dvel[j][k] * dvel[j][k];
    }

  /* Shear tensor */
  /* Diagonal terms */
  s_xx = dvel[0][0];
  s_yy = dvel[1][1];
  s_zz = dvel[2][2];
  /* Upper triangle */
  s_xy = 0.5 * (dvel[0][1] + dvel[1][0]);
  s_xz = 0.5 * (dvel[0][2] + dvel[2][0]);
  s_yz = 0.5 * (dvel[1][2] + dvel[2][1]);

  /* Trace free shear tensor */
  /* Diagonal terms */
  sfree_xx = s_xx - divv / 3.0;
  sfree_yy = s_yy - divv / 3.0;
  sfree_zz = s_zz - divv / 3.0;

  /* Construct the closure */
  /* First the diagonal parts */
  tau_xx = 2.0 * c1 * dx * sqrt(2.0 * SphP[i].Density * kturb) * sfree_xx -
           2.0 * c2 * kturb * vel_gradient_cross_term(i, 0, 0) / gradnorm2 - c3 * kturb;
  tau_yy = 2.0 * c1 * dx * sqrt(2.0 * SphP[i].Density * kturb) * sfree_yy -
           2.0 * c2 * kturb * vel_gradient_cross_term(i, 1, 1) / gradnorm2 - c3 * kturb;
  tau_zz = 2.0 * c1 * dx * sqrt(2.0 * SphP[i].Density * kturb) * sfree_zz -
           2.0 * c2 * kturb * vel_gradient_cross_term(i, 2, 2) / gradnorm2 - c3 * kturb;
  /* Upper triangle */
  tau_xy =
      2.0 * c1 * dx * sqrt(2.0 * SphP[i].Density * kturb) * s_xy - 2.0 * c2 * kturb * vel_gradient_cross_term(i, 0, 1) / gradnorm2;
  tau_xz =
      2.0 * c1 * dx * sqrt(2.0 * SphP[i].Density * kturb) * s_xz - 2.0 * c2 * kturb * vel_gradient_cross_term(i, 0, 2) / gradnorm2;
  tau_yz =
      2.0 * c1 * dx * sqrt(2.0 * SphP[i].Density * kturb) * s_yz - 2.0 * c2 * kturb * vel_gradient_cross_term(i, 1, 2) / gradnorm2;

  rate = tau_xx * s_xx + tau_yy * s_yy + tau_zz * s_zz + 2.0 * tau_xy * s_xy + 2.0 * tau_xz * s_xz + 2.0 * tau_yz * s_yz;
  rate *= SphP[i].Volume;

  return rate;
}
#endif  // TURB_APPROX_MCS_OUTPUT_RATES
#endif  // TURB_APPROX_MCS
