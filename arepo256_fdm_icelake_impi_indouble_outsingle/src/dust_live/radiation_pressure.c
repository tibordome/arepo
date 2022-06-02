/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/radiation_pressure.c
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
/* 2D interpolation in GSL requires version >= 2.1. */
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef DUST_LIVE
#ifdef DL_RADIATION

#ifndef MRT
#error "Use of option DL_RADIATION requires MRT for radiation calculations!"
#endif

#ifdef MRT_CHEM_SG
#error "Use of option DL_RADIATION not currently compatible with MRT_CHEM_SG!"
#endif

double get_Q(double a, double lam, gsl_spline2d *spline_gra, gsl_spline2d *spline_sil, gsl_interp_accel *accel_a_gra,
             gsl_interp_accel *accel_a_sil, gsl_interp_accel *accel_lam_gra, gsl_interp_accel *accel_lam_sil)
{
  double loga     = log10(a);
  double loglam   = log10(lam);
  double logQ_gra = gsl_spline2d_eval(spline_gra, loga, loglam, accel_a_gra, accel_lam_gra);
  double logQ_sil = gsl_spline2d_eval(spline_sil, loga, loglam, accel_a_sil, accel_lam_sil);
  /* For now, we assume a 50-50 mixture of graphite and silicate grains. */
  return 0.5 * (pow(10.0, logQ_gra) + pow(10.0, logQ_sil));
}

/* Return dimensionless radiation pressure cross section efficiency given grain
 * size a (in microns) and wavelength lam (in microns). */
double get_Q_pr(double a, double lam)
{
  return get_Q(a, lam, DustRP_gra.spline_Q_pr, DustRP_sil.spline_Q_pr, DustRP_gra.a_acc, DustRP_sil.a_acc, DustRP_gra.lam_acc,
               DustRP_sil.lam_acc);
}

/* Return dimensionless absorption cross section efficiency given grain
 * size a (in microns) and wavelength lam (in microns). */
double get_Q_abs(double a, double lam)
{
  return get_Q(a, lam, DustRP_gra.spline_Q_abs, DustRP_sil.spline_Q_abs, DustRP_gra.a_acc_abs, DustRP_sil.a_acc_abs,
               DustRP_gra.lam_acc_abs, DustRP_sil.lam_acc_abs);
}

double get_Q_geo(double a, double lam) { return 1.0; }

#ifdef DL_GRAIN_BINS
/* Calculate cross section (in internal area units) for a dust particle with
 * index p at wavelength lam in grain size bin j. */
double get_sigma_dust_bin(int p, double lam, double (*get_Q)(double a, double lam), int j)
{
  double Q     = get_Q(GSD.Midpoints[j], lam);
  double sigma = DTP(p).NumGrains[j] * M_PI * GSD.Midpoints[j] * GSD.Midpoints[j] * Q;
  /* Conversion from micron^2 to internal area units. */
  double int_in_um = 1.0e4 * All.UnitLength_in_cm / All.HubbleParam;
  return sigma / (int_in_um * int_in_um);
}
#endif

/* Calculate total cross section (in internal area units) for dust particle
 * with index p at wavelength lam. */
double get_sigma_dust(int p, double lam, double (*get_Q)(double a, double lam))
{
#ifdef DL_GRAIN_BINS
  double tot_sigma = 0.0;
  for(int j = 0; j < DL_GRAIN_BINS; j++)
    {
      tot_sigma += get_sigma_dust_bin(p, lam, get_Q, j);
    }
  return tot_sigma;
#else
  /* If we're not evolving grain size distributions, just fix the grain
   * density and radius. */
  double rho_grain = 2.4 / (All.HubbleParam * All.HubbleParam * All.UnitDensity_in_cgs);
  double a_grain   = 1.0e-5 * All.HubbleParam / All.UnitLength_in_cm;
  /* Assume Q(a, lambda) = 1 and just model geometric cross section. */
  return 3.0 * P[p].Mass / (4.0 * a_grain * rho_grain);
#endif
}

double get_sigma_dust_pr(int p, double lam) { return get_sigma_dust(p, lam, get_Q_pr); }

double get_sigma_dust_abs(int p, double lam) { return get_sigma_dust(p, lam, get_Q_abs); }

double get_sigma_dust_geo(int p, double lam) { return get_sigma_dust(p, lam, get_Q_geo); }

// TODO: add radiation pressure timescale / timestep constraint?
void do_dust_radiation_step(void)
{
  start_dust();
  find_drag_cells(Ndust);
  TIMER_STOPSTART(CPU_DUST, CPU_DUST_RP);

  for(int i = 0; i < Ndust; i++)
    {
      int p = DustParticle[i].index;
      if(P[p].Mass == 0.0)
        continue;

#ifdef DL_WINDS
      if(DTP(p).IsWind)
        continue;
#endif

      /* Calculate the radiation pressure cross section for this dust particle
       * in each of the MRT bins.  Loosely speaking, this is a sum of (number
       * of grains) * (pi * a**2) * (Q), for efficiency factor Q. */
      for(int b = 0; b < MRT_BINS; b++)
        {
          DustParticle[i].SigmaTot[b]    = get_sigma_dust_pr(p, DRT.lam[b]);
          DustParticle[i].SigmaTotAbs[b] = get_sigma_dust_abs(p, DRT.lam[b]);
        }

      DustParticle[i].SigmaTotGeo = get_sigma_dust_geo(p, 0.0);

#if defined(DL_THERMAL_IR) && defined(DL_GRAIN_BINS)
      /* If tracking grain size bins, need cross sections in individual bins
       * for thermal coupling. */
      for(int k = 0; k < DL_GRAIN_BINS; k++)
        {
          DustParticle[i].SigmaBinGeo[k]   = get_sigma_dust_bin(p, 0.0, get_Q_geo, k);
          DustParticle[i].SigmaBinIRAbs[k] = get_sigma_dust_bin(p, DRT.lam[UV_BINS], get_Q_abs, k);
        }
#endif
    }

    /* Given the cross sections, actually go calculate interaction between dust
     * and radiation. */
#if defined(DL_RADIATION_PRESSURE) || defined(DL_RADIATION_ABSORPTION) || defined(DL_OUTPUT_RT_FLUX)
  radiation_pressure_kernel();
#endif
#ifdef DL_THERMAL_IR
  radiation_thermal_kernel();
#endif

  /* Update dust particle state after radiation kernel. */
  for(int i = 0; i < Ndust; i++)
    {
      int p = DustParticle[i].index;
      if(P[p].Mass == 0.0)
        continue;

#ifdef DL_WINDS
      if(DTP(p).IsWind)
        continue;
#endif

#ifdef DL_RADIATION_PRESSURE
      for(int j = 0; j < 3; j++)
        {
          double dv_j = DustParticle[i].DeltaRPMomentum[j] / P[p].Mass;
          P[p].Vel[j] += dv_j;
        }
#endif

#ifdef DL_OUTPUT_RT_FLUX
      for(int b = 0; b < MRT_BINS; b++)
        {
          for(int j = 0; j < 3; j++)
            {
              DTP(p).LocalRT_F[b][j] = DustParticle[i].LocalRT_F[b][j];
            }
        }
#endif

#ifdef DL_THERMAL_IR
#ifdef DL_GRAIN_BINS
      for(int j = 0; j < DL_GRAIN_BINS; j++)
        {
          DTP(p).DustTemp[j] = DustParticle[i].DustTemp[j];
        }
#else
      DTP(p).DustTemp = DustParticle[i].DustTemp;
#endif
#endif
    }

  TIMER_STOPSTART(CPU_DUST_RP, CPU_DUST);
  end_dust();
}

void do_dust_radiation_step_first_half(void) { do_dust_radiation_step(); }

void do_dust_radiation_step_second_half(void) { do_dust_radiation_step(); }

#endif
#endif
