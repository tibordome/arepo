/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/MRT/RT_chem_ps2009.c
 * \date        06/2018
 * \author      Rahul Kannan
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 * Solves the chemical network (H and He) following Petkova & Springel 2009.
 * Can be unstable for large variations in the radiation field internsity.
 */

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"
#include "RT.h"

#ifdef MRT_UV_ONLY_DUST
void mrt_update_chemistry_only_dust(void)
{
  int idx, i;
  double dt, dtime;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

      if(All.ComovingIntegrationOn)
        dtime = dt / All.cf_hubble_a;
      else
        dtime = dt;

#ifdef MRT_SUBCYCLE
      double frac = 1.0 / ((double)(All.RTNumSubCycles));
#else
      double frac = 0.5;
#endif

      dtime *= frac;

      for(int j = 0; j < UV_BINS; j++)
        {
          double A     = dtime * SphP[i].KappaIR_P[j] * c_internal_units;
          double ratio = exp(-A);
          // double dd = SphP[i].DensPhot[j] ;

          // dd *= ratio ;

          // double rratio = (1.0+dd) / SphP[i].DensPhot[j] ;

          //	  SphP[i].DensPhot[j] = 1.0 + dd ;
          // SphP[i].Cons_DensPhot[j] = (1.0+dd) * SphP[i].Volume ;

          // SphP[i].Cons_DensPhot[j] *= ratio ;
          SphP[i].Cons_DensPhot[j] *= ratio;
#ifdef MRT_IR
          SphP[i].Cons_DensPhot[UV_BINS] += A * SphP[i].Cons_DensPhot[j];
#endif
          double modF = sqrt(SphP[i].RT_F[j][0] * SphP[i].RT_F[j][0] * SphP[i].RT_F[j][1] * SphP[i].RT_F[j][1] * SphP[i].RT_F[j][2] *
                             SphP[i].RT_F[j][2]);
          if(modF == 0)
            modF = 1.0;
          for(int num1 = 0; num1 < 3; num1++)
            {
              //  SphP[i].RT_F[j][num1] *= ratio ;
              // if(num1==0)
              //	SphP[i].Cons_RT_F[j][num1] = 0.99999999 * c_internal_units * SphP[i].Cons_DensPhot[j] / sqrt(5);
              // else if(num1==1)
              // SphP[i].Cons_RT_F[j][num1] = 0.99999999 * c_internal_units * SphP[i].Cons_DensPhot[j] * 2.0 / sqrt(5) ;
              // else
              // SphP[i].RT_F[j][num1] = MINDENSPHOT ;

              SphP[i].Cons_RT_F[j][num1] *= ratio;
              // SphP[i].Cons_RT_F[j][num1] *= ratio ;
            }
        }
#ifdef MRT_IR
      SphP[i].DensPhot[UV_BINS] = SphP[i].Cons_DensPhot[UV_BINS] / SphP[i].Volume;
#endif
    }
  return;
}
#endif

#if defined(MRT) && defined(MRT_CHEMISTRY_PS2009)

#ifndef MRT_MULTI_FREQUENCY
void mrt_update_chemistry_ps2009(void)
{
  int idx, i;
  double nH, temp, molecular_weight, rho;
  double nHII, nHI, nHI_new, nHII_new;
  double dt, dtime, c_light;
  double A, B, CC;
  double n_gamma;
  double alpha_HII, gamma_HI;
  double fac;
  double all3inv;

#ifdef MRT_INCLUDE_HE
  double alpha_HeII, alpha_HeIII, gamma_HeI, gamma_HeII;
  double nHeII, nHeIII;
  double D, E, F, G, J, L;
  double y_fac;
#endif

#ifndef MRT_MULTI_FREQUENCY
  fac             = 1.0 / All.UnitLength_in_cm / All.UnitLength_in_cm * All.HubbleParam * All.HubbleParam;
  mrt_sigma_HI[0] = 6.3e-18 * fac;
#endif

  if(All.ComovingIntegrationOn)
    all3inv = All.cf_a3inv;
  else
    all3inv = 1.0;

  fac = All.UnitTime_in_s / pow(All.UnitLength_in_cm, 3) * All.HubbleParam * All.HubbleParam;

  c_light = CLIGHT / All.UnitVelocity_in_cm_per_s;

  set_cosmo_factors_for_current_time();

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

      if(All.ComovingIntegrationOn)
        dtime = dt / All.cf_hubble_a;
      else
        dtime = dt;

#ifdef MRT_SUBCYCLE
      double frac = 1.0 / ((double)(All.RTNumSubCycles));
#else
      double frac = 0.5;
#endif

      dtime *= frac;

      rho = SphP[i].Density * All.cf_a3inv;

      nH = HYDROGEN_MASSFRAC * rho / PROTONMASS * All.UnitMass_in_g / All.HubbleParam;

      molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC + 4 * HYDROGEN_MASSFRAC * SphP[i].Ne);

      // molecular_weight = 1.0 ;

      temp = SphP[i].Utherm * GAMMA_MINUS1 * (molecular_weight * PROTONMASS / All.UnitMass_in_g) / (BOLTZMANN / All.UnitEnergy_in_cgs);

      /* collisional ionization rate */
      gamma_HI = 5.85e-11 * sqrt(temp) * exp(-157809.1 / temp) / (1.0 + sqrt(temp / 1e5)) * fac;

      double lHI = 315614.0 / temp;
#ifdef MRT_NO_OTSA
      /* alpha_A recombination coefficient */

      double alpha_HII_A = 1.269e-13 * pow(lHI, 1.503) / pow((1 + pow(lHI / 0.522, 0.47)), 1.923) * fac;
      double alpha_HII_B = 2.753e-14 * pow(lHI, 1.5) / pow((1 + pow(lHI / 2.74, 0.407)), 2.242) * fac;
      alpha_HII          = alpha_HII_A;
#else
      /* alpha_B recombination coefficient */
      alpha_HII = 2.753e-14 * pow(lHI, 1.5) / pow((1 + pow(lHI / 2.74, 0.407)), 2.242) * fac;
#endif

      n_gamma = SphP[i].DensPhot[0] * 1e63 / nH * all3inv;

      if(SphP[i].DensPhot[0] < 10.0 * MINDENSPHOT)
        n_gamma = 0.0;

      double sigma = mrt_sigma_HI[0];

      /* number of photons should be positive */
      if(n_gamma < 0 || isnan(n_gamma))
        {
          printf("NEGATIVE n_gamma: %g %d %d \n", n_gamma, i, ThisTask);
          printf("n_gamma %g mass %g ainv %g a3inv %g \n", SphP[i].DensPhot[0], P[i].Mass, 1. / All.cf_atime, All.cf_a3inv);
          terminate("111");
        }

      A  = dtime * gamma_HI * nH * SphP[i].Ne;
      B  = dtime * c_light * nH * n_gamma * sigma;
      CC = dtime * alpha_HII * nH * SphP[i].Ne;

      /* semi-implicit scheme for ionization */
      nHII = SphP[i].HII + B + A;

      nHII /= 1.0 + B + CC + A;

      if(nHII < 0.0)
        nHII = 0.0;
      if(nHII > 1.0)
        nHII = 1.0;

      if(nHII < 0.0 || nHII > 1.0 || isnan(nHII))
        {
          print_particle_info(i);
          printf("ERROR nHII %14e \n", nHII);
          terminate("333");
        }

      double nH_times_volume = P[i].Mass;
      SphP[i].nHI            = (1.0 - nHII) * nH_times_volume;
      SphP[i].nHII           = nHII * nH_times_volume;
      SphP[i].ne             = nHII * nH_times_volume;

#ifdef MRT_INCLUDE_HE
      /* collisional ionization rate */
      gamma_HeI  = 2.38e-11 * sqrt(temp) * exp(-285335.4 / temp) / (1.0 + sqrt(temp / 1e5)) * fac;
      gamma_HeII = 5.68e-12 * sqrt(temp) * exp(-631515 / temp) / (1.0 + sqrt(temp / 1e5)) * fac;

      double lHeI  = 570670.0 / temp;
      double lHeII = 1263030.0 / temp;

#ifdef MRT_NO_OTSA

      double alpha_HeII_A = 3e-14 * pow(lHeI, 0.654) * fac;
      double alpha_HeII_B = 1.26e-14 * pow(lHeI, 0.75) * fac;

      double alpha_HeIII_A = 2.538e-13 * pow(lHeII, 1.503) / pow((1.0 + pow(lHeII / 0.522, 0.47)), 1.923) * fac;
      double alpha_HeIII_B = 5.506e-14 * pow(lHeII, 1.5) / pow((1.0 + pow(lHeII / 2.74, 0.407)), 2.242) * fac;

      alpha_HeII  = alpha_HeII_A;
      alpha_HeIII = alpha_HeIII_A;
#else
      /* alpha_B recombination coefficient */
      alpha_HeII  = 1.26e-14 * pow(lHeI, 0.75) * fac;
      alpha_HeIII = 5.506e-14 * pow(lHeII, 1.5) / pow((1.0 + pow(lHeII / 2.74, 0.407)), 2.242) * fac;
#endif
      SphP[i].Ne = nHII + SphP[i].HeII + 2.0 * SphP[i].HeIII;

      D = dtime * gamma_HeII * nH * SphP[i].Ne;
      E = dtime * alpha_HeIII * nH * SphP[i].Ne;
      F = dtime * gamma_HeI * nH * SphP[i].Ne;
      J = dtime * alpha_HeII * nH * SphP[i].Ne;
      G = 0.0;
      L = 0.0;

      y_fac = (1.0 - HYDROGEN_MASSFRAC) / 4.0 / HYDROGEN_MASSFRAC;

      nHeII  = SphP[i].HeII / y_fac;
      nHeIII = SphP[i].HeIII / y_fac;

      nHeII = nHeII + G + F - ((G + F - E) / (1.0 + E)) * nHeIII;

      nHeII /= 1.0 + G + F + D + J + L + ((G + F - E) / (1.0 + E)) * (D + L);

      if(isnan(nHeII))
        {
          printf("%g %g %g %g %g \n", dtime, SphP[i].Volume, SphP[i].HeII, SphP[i].Ne, nH);
          terminate("333");
        }
      if(nHeII < 0)
        nHeII = 0.0;
      if(nHeII > 1)
        nHeII = 1.0;

      nHeIII = nHeIII + (D + L) * nHeII;

      nHeIII /= 1.0 + E;

      if(isnan(nHeIII))
        {
          printf("ERROR nHeIII %g %g\n", nHeIII, temp);
          terminate("333");
        }

      if(nHeIII < 0)
        nHeIII = 0.0;
      if(nHeIII > 1)
        nHeIII = 1.0;

      double nHeI = 1.0 - nHeII - nHeIII;
      nHeII *= y_fac;
      nHeIII *= y_fac;
      nHeI *= y_fac;

      double ne = nHII + nHeII + 2.0 * nHeIII;

      if(isnan(ne))
        terminate("%g %g %g %g \n", ne, nHII, nHeII, nHeIII);

      if(nHeI < 0)
        nHeI = 0.0;

      if(nHeI > y_fac)
        nHeI = y_fac;

      SphP[i].nHeI   = nHeI * nH_times_volume;
      SphP[i].ne     = ne * nH_times_volume;
      SphP[i].nHeII  = nHeII * nH_times_volume;
      SphP[i].nHeIII = nHeIII * nH_times_volume;

#endif

      double KK = dtime * c_light * nH;
      double ratio;
      for(int j = 0; j < UV_BINS; j++)
        {
          double sum_KK = mrt_sigma_HI[j] * SphP[i].HI
#ifdef MRT_INCLUDE_HE
                          + mrt_sigma_HeI[j] * SphP[i].HeI + mrt_sigma_HeII[j] * SphP[i].HeII
#endif
              ;

          ratio = 1.0 / (1.0 + KK * sum_KK);
          SphP[i].Cons_DensPhot[j] *= ratio;

#ifdef MRT_NO_OTSA
          SphP[i].Cons_DensPhot[j] +=
              dtime * SphP[i].Volume * (alpha_HII_A - alpha_HII_B) * SphP[i].HII * nH * SphP[i].Ne * nH / 1e63 / all3inv;
#endif
          int num1;
          for(num1 = 0; num1 < 3; num1++)
            SphP[i].Cons_RT_F[j][num1] *= ratio;
        }
    }
}

#else

/*---------------------------------------------------------------------*/
/* if the multi-frequency scheme is used*/
/*---------------------------------------------------------------------*/
void mrt_update_chemistry_ps2009(void)
{
  int idx, i, j;
  double nH, temp, molecular_weight, rho;
  double nHII, c_light, n_gamma;
  double dt, dtime;
  double A, B, CC;
  double alpha_HII, gamma_HI;
  double fac;
  double k_HI;
  double all3inv;

#ifdef MRT_INCLUDE_HE
  double alpha_HeII, alpha_HeIII, gamma_HeI, gamma_HeII;
  double nHeII, nHeIII;
  double D, E, F, G, J, L;
  double k_HeI, k_HeII;
  double y_fac;
#endif

  if(All.ComovingIntegrationOn)
    all3inv = All.cf_a3inv;
  else
    all3inv = 1.0;

  fac = All.UnitTime_in_s / pow(All.UnitLength_in_cm, 3) * All.HubbleParam * All.HubbleParam;

  c_light = CLIGHT / All.UnitVelocity_in_cm_per_s;

  set_cosmo_factors_for_current_time();

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

#if defined(BOUNDARY_REFL_FLUIDSIDE_MINID) && defined(BOUNDARY_REFL_FLUIDSIDE_MAXID) && defined(BOUNDARY_REFL_SOLIDSIDE_MINID) && \
    defined(BOUNDARY_REFL_SOLIDSIDE_MAXID)

      if(P[i].ID >= BOUNDARY_REFL_SOLIDSIDE_MINID && P[i].ID < BOUNDARY_REFL_SOLIDSIDE_MAXID)
        continue;
#endif

      /* get the photo-ionization rates */
      k_HI = 0.0;
#ifdef MRT_INCLUDE_HE
      k_HeI = k_HeII = 0.0;
#endif

      rho = SphP[i].Density * All.cf_a3inv;

      nH = HYDROGEN_MASSFRAC * rho / PROTONMASS * All.UnitMass_in_g / All.HubbleParam;

      molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC + 4 * HYDROGEN_MASSFRAC * SphP[i].Ne);

      temp = SphP[i].Utherm * GAMMA_MINUS1 * (molecular_weight * PROTONMASS / All.UnitMass_in_g * All.HubbleParam) /
             (BOLTZMANN / All.UnitEnergy_in_cgs * All.HubbleParam);

      for(j = 0; j < UV_BINS; j++)
        {
          n_gamma = SphP[i].DensPhot[j] * 1e63 * all3inv;

          k_HI += c_light * mrt_sigma_HI[j] * n_gamma;

#ifdef MRT_INCLUDE_HE

          k_HeI += c_light * mrt_sigma_HeI[j] * n_gamma;
          k_HeII += c_light * mrt_sigma_HeII[j] * n_gamma;
#endif
        }

      dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

      if(All.ComovingIntegrationOn)
        dtime = dt / All.cf_hubble_a;
      else
        dtime = dt;

#ifdef MRT_SUBCYCLE
      double frac = 1.0 / ((double)(All.RTNumSubCycles));
#else
      double frac = 0.5;
#endif

      dtime *= frac;

      /* collisional ionization rate */
      gamma_HI = 5.85e-11 * sqrt(temp) * exp(-157809.1 / temp) / (1.0 + sqrt(temp / 1e5)) * fac;

      double lHI = 315614.0 / temp;
#ifdef MRT_NO_OTSA
      /* alpha_A recombination coefficient */

      double alpha_HII_A = 1.269e-13 * pow(lHI, 1.503) / pow((1 + pow(lHI / 0.522, 0.47)), 1.923) * fac;
      double alpha_HII_B = 2.753e-14 * pow(lHI, 1.5) / pow((1 + pow(lHI / 2.74, 0.407)), 2.242) * fac;
      alpha_HII = alpha_HII_A;
#else
      /* alpha_B recombination coefficient */
      alpha_HII = 2.753e-14 * pow(lHI, 1.5) / pow((1 + pow(lHI / 2.74, 0.407)), 2.242) * fac;
      // alpha_HII = 2.59e-13 * pow(temp / 1e4, -0.7) * fac;
#endif

      A = dtime * gamma_HI * nH * SphP[i].Ne;
      B = dtime * k_HI;
      CC = dtime * alpha_HII * nH * SphP[i].Ne;

      /* semi-implicit scheme for ionization */
      nHII = SphP[i].HII + B + A;

      nHII /= 1.0 + B + CC + A;

      if(nHII < 0)
        nHII = 0.0;

      if(nHII > 1.0)
        nHII = 1.0;

      double nH_times_volume = P[i].Mass;
      SphP[i].nHI = (1.0 - nHII) * nH_times_volume;
      SphP[i].nHII = nHII * nH_times_volume;
      SphP[i].ne = nHII * nH_times_volume;

#ifdef MRT_INCLUDE_HE
      /* collisional ionization rate */
      gamma_HeI = 2.38e-11 * sqrt(temp) * exp(-285335.4 / temp) / (1.0 + sqrt(temp / 1e5)) * fac;
      gamma_HeII = 5.68e-12 * sqrt(temp) * exp(-631515 / temp) / (1.0 + sqrt(temp / 1e5)) * fac;

      double lHeI = 570670.0 / temp;
      double lHeII = 1263030.0 / temp;

#ifdef MRT_NO_OTSA

      double alpha_HeII_A = 3e-14 * pow(lHeI, 0.654) * fac;
      double alpha_HeII_B = 1.26e-14 * pow(lHeI, 0.75) * fac;

      double alpha_HeIII_A = 2.538e-13 * pow(lHeII, 1.503) / pow((1.0 + pow(lHeII / 0.522, 0.47)), 1.923) * fac;
      double alpha_HeIII_B = 5.506e-14 * pow(lHeII, 1.5) / pow((1.0 + pow(lHeII / 2.74, 0.407)), 2.242) * fac;

      alpha_HeII = alpha_HeII_A;
      alpha_HeIII = alpha_HeIII_A;
#else
      /* alpha_B recombination coefficient */
      alpha_HeII  = 1.26e-14 * pow(lHeI, 0.75) * fac;
      alpha_HeIII = 5.506e-14 * pow(lHeII, 1.5) / pow((1.0 + pow(lHeII / 2.74, 0.407)), 2.242) * fac;
#endif

      SphP[i].Ne = nHII + SphP[i].HeII + 2.0 * SphP[i].HeIII;

      D = dtime * gamma_HeII * nH * SphP[i].Ne;
      E = dtime * alpha_HeIII * nH * SphP[i].Ne;
      F = dtime * gamma_HeI * nH * SphP[i].Ne;
      J = dtime * alpha_HeII * nH * SphP[i].Ne;
      G = dtime * k_HeI;
      L = dtime * k_HeII;

      y_fac = (1.0 - HYDROGEN_MASSFRAC) / 4.0 / HYDROGEN_MASSFRAC;

      nHeII = SphP[i].HeII / y_fac;
      nHeIII = SphP[i].HeIII / y_fac;

      nHeII = nHeII + G + F - ((G + F - E) / (1.0 + E)) * nHeIII;

      nHeII /= 1.0 + G + F + D + J + L + ((G + F - E) / (1.0 + E)) * (D + L);

      if(isnan(nHeII))
        {
          printf("%g %g %g %g %g \n", dtime, SphP[i].Volume, SphP[i].HeII, SphP[i].Ne, nH);
          terminate("333");
        }
      if(nHeII < 0)
        nHeII = 0.0;
      if(nHeII > 1)
        nHeII = 1.0;

      nHeIII = nHeIII + (D + L) * nHeII;

      nHeIII /= 1.0 + E;

      if(isnan(nHeIII))
        {
          printf("ERROR nHeIII %g %g %g\n", nHeIII, temp, k_HeII);
          terminate("333");
        }

      if(nHeIII < 0)
        nHeIII = 0.0;
      if(nHeIII > 1)
        nHeIII = 1.0;

      double nHeI = 1.0 - nHeII - nHeIII;
      nHeII *= y_fac;
      nHeIII *= y_fac;
      nHeI *= y_fac;

      double ne = nHII + nHeII + 2.0 * nHeIII;

      if(isnan(ne))
        terminate("%g %g %g %g \n", ne, nHII, nHeII, nHeIII);

      if(nHeI < 0)
        nHeI = 0.0;

      if(nHeI > y_fac)
        nHeI = y_fac;

      SphP[i].nHeI = nHeI * nH_times_volume;
      SphP[i].ne = ne * nH_times_volume;
      SphP[i].nHeII = nHeII * nH_times_volume;
      SphP[i].nHeIII = nHeIII * nH_times_volume;

#endif
      double KK = dtime * c_light * nH;
      double ratio;
      for(j = 0; j < UV_BINS; j++)
        {
          double sum_KK = mrt_sigma_HI[j] * SphP[i].HI
#ifdef MRT_INCLUDE_HE
                          + mrt_sigma_HeI[j] * SphP[i].HeI + mrt_sigma_HeII[j] * SphP[i].HeII
#endif
              ;

#ifdef MRT_IR
          double AA = SphP[i].KappaIR_P[j] / nH;
          sum_KK += AA;

          SphP[i].Cons_DensPhot[UV_BINS] += KK * AA * MeanPhotonEnergy[j] * SphP[i].Cons_DensPhot[j] * 1e63;
#endif
          ratio = 1.0 / (1.0 + KK * sum_KK);

          SphP[i].Cons_DensPhot[j] *= ratio;

#ifdef MRT_NO_OTSA
          if(nu[j] == 13.6)
            SphP[i].Cons_DensPhot[j] +=
                dtime * SphP[i].Volume * ((alpha_HII_A - alpha_HII_B) * nH) * SphP[i].HII * SphP[i].Ne * (nH / 1e63) / all3inv;
          else if(nu[j] == 24.6)
            SphP[i].Cons_DensPhot[j] +=
                dtime * SphP[i].Volume * ((alpha_HeII_A - alpha_HeII_B) * nH) * SphP[i].HeII * SphP[i].Ne * (nH / 1e63) / all3inv;
          else if(nu[j] == 54.4)
            SphP[i].Cons_DensPhot[j] +=
                dtime * SphP[i].Volume * ((alpha_HeIII_A - alpha_HeIII_B) * nH) * SphP[i].HeIII * SphP[i].Ne * (nH / 1e63) / all3inv;

#endif
          int num1;
          for(num1 = 0; num1 < 3; num1++)
            SphP[i].Cons_RT_F[j][num1] *= ratio;
        }
#ifdef MRT_IR
      SphP[i].DensPhot[UV_BINS] = SphP[i].Cons_DensPhot[UV_BINS] / SphP[i].Volume;
#endif
    }
}
#endif

#endif
