/*! * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/MRT/RT_init.c
 * \date        06/2018
 * \author      Rahul Kannan, Federico Marincaci
 * \brief       Initialization routines for cross sections, photoheating rates
 * \details     and radiation pressure terms
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <math.h>
#include <mpi.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../domain.h"
#include "../proto.h"
#include "../voronoi.h"

#include "RT_proto.h"

#ifdef MRT_CHEM_SG

// piecewise Cross-Section and energies for Ionization of H2 ( Liu & Shemansky 2012 )
#define NPW 11  // Number of PieceWise Intervals for H2

// number of elements = NPW
static double SigmaPW[] = {0.09e-18, 1.15e-18, 3.0e-18, 5.0e-18, 6.75e-18, 8.0e-18, 9.0e-18, 9.5e-18, 9.8e-18, 10.1e-18, 9.8e-18};
// number of elements = NPW + 1
static double EnergyPW[] = {15.20, 15.45, 15.70, 15.95, 16.20, 16.40, 16.65, 16.85, 17.00, 17.20, 17.65, 18.10};

static double H2_ionization_crosssection(double energy);

#endif

#ifndef MRT_NO_UV

void mrt_get_sigma(void)
{
  double fac, fac_two;

#ifdef MRT_CHEM_SG
  fac = 1.0;
#else
  fac = 1.0 / All.UnitLength_in_cm / All.UnitLength_in_cm * All.HubbleParam * All.HubbleParam;
#endif

#if !defined(MRT_MULTI_FREQUENCY) && !defined(MRT_NO_UV)
  mrt_sigma_HI[0] = 6.3e-18 * fac;
  P_HI[0]         = 13.6 * ELECTRONVOLT_IN_ERGS / All.UnitEnergy_in_cgs * All.HubbleParam;

#else
  int i, j, integral;
  double e, d_nu, e_start, e_end;
  double sum_HI_sigma, sum_HI_G, sum_H2_sigma, sum_H2_G;
  double I_nu;
  double sig, f;
#ifdef MRT_INCLUDE_HE
  double sum_HeI_sigma, sum_HeII_sigma;
  double sum_HeI_G, sum_HeII_G;
#endif

  integral = 10000;

#ifdef MRT_CHEM_SG
  fac_two  = ELECTRONVOLT_IN_ERGS;
#else
  fac_two    = ELECTRONVOLT_IN_ERGS / All.UnitEnergy_in_cgs * All.HubbleParam;
#endif

#if defined(MRT_MULTI_FREQUENCY) && !defined(MRT_CHEM_SG) && !defined(MRT_PHOTOELECTRIC_HEATING)
  nu[0]    = 13.6;
  nu[1]    = 24.6;
  nu[2]    = 54.4;
  nu[3]    = 100.0;
#endif

#if defined(MRT_MULTI_FREQUENCY) && defined(MRT_PHOTOELECTRIC_HEATING)
  nu[0]    = 1.0;
  nu[1]    = 13.6;
  nu[2]    = 24.6;
  nu[3]    = 54.4;
  nu[4]    = 100.0;
#endif

#if defined(MRT_MULTI_FREQUENCY) && defined(MRT_CHEM_SG)
  nu[0]    = 1.0;
  nu[1]    = 11.2;
  nu[2]    = 13.6;
  nu[3]    = 15.2;
  nu[4]    = 24.6;
  nu[5]    = 100.0;
#endif

  for(i = 0; i < UV_BINS; i++)
    {

#ifdef MRT_CHEM_SG
      sum_H2_sigma = 0.0;
      sum_H2_G     = 0.0;
#endif
      sum_HI_sigma = 0.0;
      sum_HI_G     = 0.0;
#ifdef MRT_INCLUDE_HE
      sum_HeI_G = sum_HeII_G = 0.0;
      sum_HeI_sigma          = 0.0;
      sum_HeII_sigma         = 0.0;
#endif

      e_start = nu[i];

      e_end = nu[i + 1];

      d_nu            = (e_end - e_start) / (float)(integral - 1);

#ifdef MRT_CHEM_SG
      mrt_sigma_H2[i] = 0.0;
      G_H2[i]         = 0.0;
      P_H2[i]         = 0.0;
#endif

      mrt_sigma_HI[i]   = 0.0;
      G_HI[i]           = 0.0;
      P_HI[i]           = 0.0;
#ifdef DURRIVE_BATTERY
      Pelec_HI[i]       = 0.0;
#endif
#ifdef MRT_INCLUDE_HE
      mrt_sigma_HeI[i]  = 0.0;
      mrt_sigma_HeII[i] = 0.0;
      G_HeI[i] = G_HeII[i] = 0.0;
      P_HeI[i] = P_HeII[i] = 0.0;
#ifdef DURRIVE_BATTERY
      Pelec_HeI[i] = Pelec_HeII[i] = 0.0;
#endif /* def DURRIVE_BATTERY */
#endif

      for(j = 0; j < integral; j++)
        {
          /* midpoint integration */
          e = e_start + (j + 0.5) * d_nu;

          I_nu = get_spectrum(e, 0, 0);

          if(nu[i] >= 13.6)
            {
              f = sqrt((e / 13.6) - 1.0);

              if(f < 1.0e-6)
                sig = 6.3e-18 * pow(13.6 / e, 4) * exp(4 * f * f / 3); /* for small enough f use Taylor expansion */
              else
                sig = 6.3e-18 * pow(13.6 / e, 4) * exp(4 - (4 * atan(f) / f)) / (1.0 - exp(-2 * M_PI / f));

              mrt_sigma_HI[i] += d_nu * sig * I_nu / e;

              sum_HI_sigma += d_nu * I_nu / e;

              G_HI[i] += d_nu * sig * (e - 13.6) * I_nu / e;

              P_HI[i] += d_nu * sig * e * I_nu / e;
#ifdef DURRIVE_BATTERY
              Pelec_HI[i] += d_nu * f_mt_HI(e) * sig * e * I_nu / e;
#endif
              sum_HI_G += d_nu * sig * I_nu / e;
            }

#ifdef MRT_CHEM_SG
          if(nu[i] >= 15.2)
            {
              sig = H2_ionization_crosssection(e);

              mrt_sigma_H2[i] += d_nu * sig * I_nu / e;

              sum_H2_sigma += d_nu * I_nu / e;

              G_H2[i] += d_nu * sig * (e - 15.2) * I_nu / e;

              P_H2[i] += d_nu * sig * e * I_nu / e;

              sum_H2_G += d_nu * sig * I_nu / e;
            }
#endif

#ifdef MRT_INCLUDE_HE
          if(nu[i] >= 24.6)
            {
              f = sqrt((e / 24.6) - 1.0);

              if(f < 1.0e-6)
                sig = 7.83e-18 * pow(24.6 / e, 4) * exp(4 * f * f / 3); /* for small enough f use Taylor expansion */
              else
                sig = 7.83e-18 * pow(24.6 / e, 4) * exp(4 - (4 * atan(f) / f)) / (1.0 - exp(-2 * M_PI / f));

              mrt_sigma_HeI[i] += d_nu * sig * I_nu / e;

              sum_HeI_sigma += d_nu * I_nu / e;

              G_HeI[i] += d_nu * sig * (e - 24.6) * I_nu / e;

              P_HeI[i] += d_nu * sig * e * I_nu / e;
#ifdef DURRIVE_BATTERY
              Pelec_HeI[i] += d_nu * f_mt_HeI(e) * sig * e * I_nu / e;
#endif
              sum_HeI_G += d_nu * sig * I_nu / e;
            }

          if(nu[i] >= 54.4)
            {
              f = sqrt((e / 54.4) - 1.0);

              if(f < 1.0e-6)
                sig = 1.58e-18 * pow(54.4 / e, 4) * exp(4 * f * f / 3); /* for small enough f use Taylor expansion */
              else
                sig = 1.58e-18 * pow(54.4 / e, 4) * exp(4 - (4 * atan(f) / f)) / (1.0 - exp(-2 * M_PI / f));

              mrt_sigma_HeII[i] += d_nu * sig * I_nu / e;

              sum_HeII_sigma += d_nu * I_nu / e;

              G_HeII[i] += d_nu * sig * (e - 54.4) * I_nu / e;

              P_HeII[i] += d_nu * sig * e * I_nu / e;
#ifdef DURRIVE_BATTERY
              Pelec_HeII[i] += d_nu * f_mt_HeII(e) * sig * e * I_nu / e;
#endif
              sum_HeII_G += d_nu * sig * I_nu / e;
            }
#endif
        }

      if(nu[i] >= 13.6)
        {
          mrt_sigma_HI[i] *= fac / sum_HI_sigma;
          G_HI[i] *= fac_two / sum_HI_G;
          P_HI[i] *= fac_two / sum_HI_G;
#ifdef DURRIVE_BATTERY
          Pelec_HI[i] *= fac_two / sum_HI_G;
#endif
        }

#ifdef MRT_CHEM_SG
      if(nu[i] >= 11.2 && nu[i] < 13.6)
        {
          mrt_sigma_H2[i] = 2.1e-19 * fac;
          G_HI[i]         = 0.0;
          P_HI[i]         = 0.0;
        }
      if(nu[i] >= 15.2)
        {
          mrt_sigma_H2[i] *= fac / sum_H2_sigma;
          G_H2[i] *= fac_two / sum_H2_G;
          P_H2[i] *= fac_two / sum_H2_G;
        }
#endif

#ifdef MRT_INCLUDE_HE
      if(nu[i] >= 24.6)
        {
          mrt_sigma_HeI[i] *= fac / sum_HeI_sigma;
          G_HeI[i] *= fac_two / sum_HeI_G;
          P_HeI[i] *= fac_two / sum_HeI_G;
#ifdef DURRIVE_BATTERY
          Pelec_HeI[i] *= fac_two / sum_HeI_G;
#endif
        }

      if(nu[i] >= 54.4)
        {
          mrt_sigma_HeII[i] *= fac / sum_HeII_sigma;
          G_HeII[i] *= fac_two / sum_HeII_G;
          P_HeII[i] *= fac_two / sum_HeII_G;
#ifdef DURRIVE_BATTERY
          Pelec_HeII[i] *= fac_two / sum_HeII_G;
#endif
        }
#endif
    }
  /*for(int kk=0;kk<UV_BINS;kk++)
        {
      /*              if(kk==2)
        mrt_sigma_HI[kk] = 5.0e-18 ;
      else
      mrt_sigma_HI[kk] = 0.0 ;

      //if(kk==2)
      //{
      //  mrt_sigma_HI[kk] = 1.63e-18 ;
      //  G_HI[kk] = 6.4*fac_two ;
      //}
      // else
      //{
      //  mrt_sigma_HI[kk] = 0.0 ;
      //  G_HI[kk] = 0.0 ;
      //	}


      /*if(kk==1)
        mrt_sigma_H2[kk] = 2.1e-19 ;
      else if(kk==2)
        mrt_sigma_H2[kk] = 3.6e-18 ;
      else
      mrt_sigma_H2[kk] = 0.0 ;*/
  /*  if(kk == 2)
        mrt_sigma_HI[kk] = 6.3e-18 ;
      else
        mrt_sigma_HI[kk] = 0.0 ;

      mrt_sigma_H2[kk] = 0.0 ;
      mrt_sigma_HeI[kk] = 0.0 ;
      mrt_sigma_HeII[kk] = 0.0 ;

      G_H2[kk] = 0.0 ;
      G_HI[kk] = 0.0 ;
      G_HeI[kk] = 0.0 ;
      G_HeII[kk] = 0.0 ;

      P_H2[kk] = 0.0 ;
      P_HI[kk] = 0.0 ;
      if(kk==2)
        G_HI[kk] = 6.4 * ELECTRONVOLT_IN_ERGS ;
      else
        G_HI[kk] = 0.0 ;

      P_HeI[kk] = 0.0 ;
      P_HeII[kk] = 0.0 ;
      }*/
  if(ThisTask == 0)
    for(i = 0; i < UV_BINS; i++)
      {
        printf("RT SIGMA: %g %g %g | %g %g %g | %g %g %g \n", mrt_sigma_HI[i] / fac, G_HI[i] / fac_two, P_HI[i] / fac_two,
               mrt_sigma_HeI[i] / fac, G_HeI[i] / fac_two, P_HeI[i] / fac_two, mrt_sigma_HeII[i] / fac, G_HeII[i] / fac_two,
               P_HeII[i] / fac_two);
#ifdef MRT_CHEM_SG
        printf("RT SIGMA - CHEM_SG H2 ionization rates: %g %g %g \n", mrt_sigma_H2[i] / fac, G_H2[i] / fac_two, P_H2[i] / fac_two);
#endif
      }

#ifdef DURRIVE_BATTERY
  if(ThisTask == 0)
    {
      printf("DURRIVE BATTERY - Pelec rates (UV_BINS = %i): Pelec_HI, Pelec_HeI, Pelec_HeII\n", UV_BINS);
      for(i = 0; i < UV_BINS; i++)
        printf("   %g %g %g\n", Pelec_HI[i] / fac_two, Pelec_HeI[i] / fac_two, Pelec_HeII[i] / fac_two);
    }
#endif

#endif
}

void mrt_get_luminosities(void)
{
#ifndef MRT_MULTI_FREQUENCY
  lum[0] = 1.;

#ifdef MRT_BH
  for(int bin = 0; bin < MRT_BINS; bin++)
    {
      lum[bin] = 1.0;
      if(bin < UV_BINS)
        MeanPhotonEnergy[bin] = 13.6;
    }
#endif
#else

#ifndef MRT_CHEM_SG
  double fac = ELECTRONVOLT_IN_ERGS / All.UnitEnergy_in_cgs * All.HubbleParam;
#else
  double fac = ELECTRONVOLT_IN_ERGS;
#endif

  int i, j, integral;
  double e, d_nu, e_start, e_end;
  double I_nu;
  double lum_tot = 0.0;
  //  double E_tot = 0.0 ;

  integral = 10000;

  for(i = 0; i < UV_BINS; i++)
    {
      e_start = nu[i];
      e_end   = nu[i + 1];

      d_nu = (e_end - e_start) / (float)(integral - 1);

      lum[i]              = 0.0;
      MeanPhotonEnergy[i] = 0.0;

      for(j = 0; j < integral; j++)
        {
          /* midpoint integration */
          e = e_start + (j + 0.5) * d_nu;

          I_nu = get_spectrum(e, 0, 0);

          lum[i] += d_nu * I_nu / e;
          MeanPhotonEnergy[i] += d_nu * I_nu;

#if(defined(MRT_SOURCES) && defined(MRT_STARS)) || defined(MRT_LOCAL_FEEDBACK)
          for(int i_metallicity = 0; i_metallicity < N_metallicity; i_metallicity++)
            {
              for(int i_age = 0; i_age < N_age; i_age++)
                {
                  I_nu = get_spectrum(e, i_age, i_metallicity);
                  lum_tab[UV_BINS * (i_metallicity * N_age + i_age) + i] += d_nu * I_nu / e;
                }
            }
#endif
        }

      lum_tot += lum[i];
    }

  for(i = 0; i < UV_BINS; i++)
    {
      MeanPhotonEnergy[i] = MeanPhotonEnergy[i] / lum[i] * fac;
      lum[i] /= lum_tot;
    }

  //  MeanPhotonEnergy[2] = 13.6 * fac ;

  if(ThisTask == 0)
    for(i = 0; i < UV_BINS; i++)
      printf("RT: luminosity fractions i = %d \t %g \t %g \t\t\t Mean energies = %g\n", i, lum[i], lum_tot, MeanPhotonEnergy[i] / fac);

#ifdef MRT_SOURCES
  PhotRate = lum_tot;
#endif

#endif
}
#endif

#ifdef MRT_CHEM_SG

static double H2_ionization_crosssection(double energy)
{
  double val;
  if(energy >= 11.2 && energy < 15.2)
    val = 2.47e-18;
  else if(energy >= 15.2 && energy < 18.1)
    {
      for(int i = 0; i < NPW; i++)
        {
          if(energy >= EnergyPW[i] && energy < EnergyPW[i + 1])
            val = SigmaPW[i];
        }
    }
  else
    val = 9.75e-18 * pow(18.1 / energy, 3);

  return val;
}

#endif
