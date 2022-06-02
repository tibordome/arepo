/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/MRT/RT_sgchem.c
 * \date        06/2018
 * \author      Rahul Kannan
 * \brief       Routine to couple MRT to SGchem module
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

/* Routine to couple MRT to SGCHEM network - inputs iotnization and heating rates,
absorbs photons and photon flux and adds radiation pressure when option is turned on.
*/
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "RT.h"

#ifdef MRT_CHEM_SG

/* prevent duplicate definitions of constants with src/simplex/sx_def.h */
#ifndef SIMPLEX
static const int RIH   = 0;
static const int HRIH  = 1;
static const int RIH2  = 2;
static const int HRIH2 = 3;
static const int RDH2  = 4;
#endif
static const int RPH = 5;
#ifndef SIMPLEX
static const int RIHE  = 6;
static const int HRIHE = 7;

static const int F056 = 0;
static const int F112 = 1;
static const int F136 = 2;
static const int F152 = 3;
static const int F246 = 4;
#endif

static const int NRATES = 8;

static const double minMassFrac = 1e-20;  // minimum of the mass fraction to prevent NaN
static const double uvPumpFact  = 6.94;
static const double frac_peh    = 0.55;

static double fac;
static double cspeed;

static struct sitestruct
{
  double xH, xH2, xHp, xHe, xHep;
  double initxH, initxH2, initxHp, initxHe, initxHep;
  double density, volume, numdens, nNucleons, numdr;
  double rates[8];
  double nphot[UV_BINS];
  double kappa[UV_BINS];
  double rp[UV_BINS];
  double rp_dust[UV_BINS];
#ifdef MRT_IR
  double reprocessed_EIR;
#endif
} site;

static void initialize_site(int i);
static void solve_for_rates(int i, double dt);
static void calculate_sgchem_rates_individual(int i, double dt);
static void set_out_rates_and_photon_numbers(int i, double dt);

void set_rates_sgchem(void)
{
  int idx, i;
  double dt;

  fac = 1.0 / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitTime_in_s / pow(All.HubbleParam, 3);

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

      if(All.ComovingIntegrationOn)
        dt /= All.cf_hubble_a;

      if(dt <= 0.0)
        continue;

#ifdef MRT_SUBCYCLE
      dt /= ((double)(All.RTNumSubCycles));
#else
      dt *= 0.5;
#endif

      cspeed = 2.99792458e10 / All.UnitVelocity_in_cm_per_s;

      calculate_sgchem_rates_individual(i, dt);
    }

  return;
}

static void calculate_sgchem_rates_individual(int i, double dt)
{
  initialize_site(i);

  double dt_sec = dt * All.UnitTime_in_s;

  solve_for_rates(i, dt_sec);  // pass time in sec

  set_out_rates_and_photon_numbers(i, dt);  // pass time in internal units

  return;
}

static void initialize_site(int i)
{
  memset(site.rates, 0, NRATES * sizeof(double));
  memset(site.rp, 0, UV_BINS * sizeof(double));
  memset(site.rp_dust, 0, UV_BINS * sizeof(double));

  for(int ll = 0; ll < SGCHEM_NUM_ADVECTED_SPECIES; ll++)
    if(SphP[i].TracAbund[ll] < 0.0)
      SphP[i].TracAbund[ll] = 1e-20;

  site.initxH = site.xH = SphP[i].TracAbund[IHATOM];
  site.initxH2 = site.xH2 = SphP[i].TracAbund[IH2];
  site.initxHp = site.xHp = SphP[i].TracAbund[IHP];
  site.initxHe = site.xHe = SphP[i].TracAbund[IHEATOM];
  site.initxHep = site.xHep = SphP[i].TracAbund[IHEP];

  if(isnan(site.xH))
    terminate("xH is nan\n");

  // need something to calculate otherwise enjoy some NaNs
  if(site.xH == 0.0 || site.xHp == 1.0)
    site.xH = minMassFrac;
  if(site.xH2 == 0.0)
    site.xH2 = minMassFrac;
  if(site.xHe == 0.0)
    site.xHe = minMassFrac;

  site.density = SphP[i].Density * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv;

  site.volume = SphP[i].Volume * pow(All.UnitLength_in_cm / All.HubbleParam, 3) / All.cf_a3inv;

  site.numdens = site.density / ((1. + 4. * ABHE) * PROTONMASS);

  site.nNucleons = site.numdens * site.volume;

  for(int kk = 0; kk < UV_BINS; kk++)
    {
      site.nphot[kk] = SphP[i].DensPhot[kk] * SphP[i].Volume * 1e63;
      site.kappa[kk] = SphP[i].KappaIR_R[kk] / All.UnitLength_in_cm;
      // printf("kappa = %g %g %g \n",  site.kappa[kk], SphP[i].KappaIR_P[kk], site.density) ;
    }

#ifdef MRT_IR
  site.reprocessed_EIR = 0.0;
#endif
  return;
}

static void solve_for_rates(int i, double dt)
{
  double sigma, sigmaH, sigmaH2, sigmaHe;
  double nPhot, nPhotH, nPhotH2, nPhotHe;
  double nSpec, nSpecH, nSpecH2, nSpecHe;
  double attRatio;

  double nPhotDis;

  for(int kk = UV_BINS - 1; kk > -1; kk--)
    {
      if(site.nphot[kk] <= 0.0)
        continue;
#if defined(MRT_RADIATION_PRESSURE) || defined(MRT_IR)
      site.rp_dust[kk] = site.nphot[kk] * (1.0 - exp(-dt * CLIGHT * site.kappa[kk])) * MeanPhotonEnergy[kk];
#ifdef MRT_IR
      site.reprocessed_EIR += site.rp_dust[kk];
#endif
#endif
      site.nphot[kk] *= exp(-dt * CLIGHT * site.kappa[kk]);

      if(kk == F056)
        {
          site.rates[RPH] = frac_peh * c_internal_units * SphP[i].DensPhot[F056] * 1e63 * MeanPhotonEnergy[F056] * fac;
          //	  site.nphot[kk] *= exp(- dt * CLIGHT * site.kappa[kk]) ;
          continue;
        }

      // calculate number of available species
      nSpecH  = site.nNucleons * site.xH;
      nSpecH2 = site.nNucleons * site.xH2;
      nSpecHe = site.nNucleons * site.xHe;
      nSpec   = nSpecH + nSpecH2 + nSpecHe;

      //      if(site.xH < 0.0 && All.Time != 0)
      // terminate("Negative densities xH = %g \n", site.xH) ;

      // double fac = ELECTRONVOLT_IN_ERGS / All.UnitEnergy_in_cgs * All.HubbleParam;
      if(nSpec == 0)
        continue;  // return if there are no available species for a reaction

      // some initial constants
      sigmaH = mrt_sigma_HI[kk] * site.xH;

      if(kk == F112)
        sigmaH2 = mrt_sigma_H2[kk] * site.xH2 * (1. + uvPumpFact);
      else
        sigmaH2 = mrt_sigma_H2[kk] * site.xH2;

      sigmaHe = mrt_sigma_HeI[kk] * site.xHe;
      sigma   = sigmaH + sigmaH2 + sigmaHe;

      // printf("%g %g %g %g \n", sigmaH, sigmaH2, sigmaHe, sigma) ;

      if(sigma == 0)
        continue;  // skip this frequency if there is not a cross-section for the reaction

      // caclulate number of attenuated photons

      nPhot   = site.nphot[kk] * (1.0 - exp(-dt * CLIGHT * site.numdens * sigma));
      nPhotH  = nPhot * sigmaH / sigma;
      nPhotH2 = nPhot * sigmaH2 / sigma;
      nPhotHe = nPhot * sigmaHe / sigma;

      if(kk == F112)
        nPhotDis = nPhotH2 / (1. + uvPumpFact);

      if(nPhotH > nSpecH)
        nPhotH = nSpecH;

      if(kk == F112)
        {
          if(nPhotDis > nSpecH2)
            {
              nPhotDis = nSpecH2;
              nPhotH2  = (1. + uvPumpFact) * nPhotDis;
            }
        }
      else if(nPhotH2 > nSpecH2)
        nPhotH2 = nSpecH2;

      if(nPhotHe > nSpecHe)
        nPhotHe = nSpecHe;

      // update rates and mass fractions of H
      if(site.initxH > 0.0)
        {
          attRatio = nPhotH / (site.nNucleons * site.initxH);
          site.rates[RIH] += attRatio / dt;
          site.rates[HRIH] += attRatio / dt * G_HI[kk];

#ifdef MRT_RADIATION_PRESSURE
          site.rp[kk] += nPhotH * P_HI[kk];
#endif
          //      printf("before %g \t", site.xH) ;
          site.xH -= attRatio * site.initxH;
        }
      // printf("after %g \n", site.xH) ;

      // update rates and mass fractions of H2
      if(site.initxH2 > 0.0)
        {
          if(kk == F112)
            attRatio = nPhotDis / (site.nNucleons * site.initxH2);
          else
            attRatio = nPhotH2 / (site.nNucleons * site.initxH2);

          if(kk == F112)
            site.rates[RDH2] += attRatio / dt;
          else
            site.rates[RIH2] += attRatio / dt;

          site.rates[HRIH2] += attRatio / dt * G_H2[kk];

#ifdef MRT_RADIATION_PRESSURE
          site.rp[kk] += nPhotH2 * P_H2[kk];
#endif

          site.xH2 -= attRatio * site.initxH2;
        }

      // update rates and mass fractions of He
      if(site.initxHe > 0.0)
        {
          attRatio = nPhotHe / (site.nNucleons * site.initxHe);
          site.rates[RIHE] += attRatio / dt;
          site.rates[HRIHE] += attRatio / dt * G_HeI[kk];

#ifdef MRT_RADIATION_PRESSURE
          site.rp[kk] += nPhotHe * P_HeI[kk];
#endif

          site.xHe -= attRatio * site.initxHe;
        }

      if(kk == F112)
        nPhotH2 *= 500.0;

      nPhot = nPhotH + nPhotH2 + nPhotHe;

      site.nphot[kk] = (site.nphot[kk] - nPhot < 0) ? 0 : site.nphot[kk] - nPhot;

      if(isnan(site.rates[RIH]))
        terminate("attnRatio = %g xH = %g %g %g\n", attRatio, site.xH, nPhot, site.nNucleons);
    }

  return;
}

static void set_out_rates_and_photon_numbers(int i, double dt)
{
  for(int kk = 0; kk < NRATES; kk++)
    SphP[i].PhotonRates[kk] = site.rates[kk];

#if defined(MRT_RADIATION_PRESSURE)
  double KE_old = (SphP[i].Momentum[0] * SphP[i].Momentum[0] + SphP[i].Momentum[1] * SphP[i].Momentum[1] +
                   SphP[i].Momentum[2] * SphP[i].Momentum[2]) /
                  P[i].Mass / 2.0;
#endif

  for(int ll = 0; ll < UV_BINS; ll++)
    {
      if(site.nphot[ll] / SphP[i].Volume / 1e63 < MINDENSPHOT)
        site.nphot[ll] = MINDENSPHOT * 1e63 * SphP[i].Volume;

      double Nnew  = site.nphot[ll] * 1e-63;
      double ratio = Nnew / SphP[i].Cons_DensPhot[ll];
      SphP[i].Cons_DensPhot[ll] *= ratio;

      double sum = 0.0;
      for(int num1 = 0; num1 < 3; num1++)
        {
          sum += SphP[i].RT_F[ll][num1] * SphP[i].RT_F[ll][num1];
          SphP[i].Cons_RT_F[ll][num1] *= ratio;
        }

#ifdef MRT_RADIATION_PRESSURE
      if(sum == 0.0)
        sum = 1.0;

      for(int num1 = 0; num1 < 3; num1++)
        {
          double frac = SphP[i].RT_F[ll][num1] / sqrt(sum);

          SphP[i].Momentum[num1] += ((site.rp[ll] + site.rp_dust[ll]) / All.UnitEnergy_in_cgs) * frac / cspeed;
        }
#endif
    }

#ifdef MRT_IR
  double A     = dt * SphP[i].KappaIR_R[UV_BINS] * c_internal_units;
  double ratio = 1.0 / (1 + A);
  for(int num1 = 0; num1 < 3; num1++)
    {
#ifdef MRT_RADIATION_PRESSURE
      SphP[i].Momentum[num1] += dt * SphP[i].KappaIR_R[UV_BINS] * SphP[i].RT_F[UV_BINS][num1] * SphP[i].Volume / cspeed;
#endif
      SphP[i].Cons_RT_F[UV_BINS][num1] *= ratio;
    }
  SphP[i].Cons_DensPhot[UV_BINS] += site.reprocessed_EIR / All.UnitEnergy_in_cgs;
#endif

#if defined(MRT_RADIATION_PRESSURE)
  double KE_new = (SphP[i].Momentum[0] * SphP[i].Momentum[0] + SphP[i].Momentum[1] * SphP[i].Momentum[1] +
                   SphP[i].Momentum[2] * SphP[i].Momentum[2]) /
                  P[i].Mass / 2.0;

  SphP[i].Energy += KE_new - KE_old;
#endif

  return;
}

#endif
