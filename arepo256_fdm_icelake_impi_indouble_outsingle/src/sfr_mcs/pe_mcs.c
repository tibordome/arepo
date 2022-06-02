/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sfr_mcs/pe_mcs.c
 * \date        04/2019
 * \author      Matthew C. Smith
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

#ifdef PE_MCS

#define HABING_UNIT 5.29e-14 /*in erg cm^-3*/

#ifndef PE_MCS_POSTSHIELD
#define PE_MCS_POSTSHIELD 0
#endif

void init_pe(void)
{
  /* Factor to multiply luminosity / r^2 to get energy density in Habing units */
  All.Factor_FUV =
      1.0 / (4.0 * M_PI * CLIGHT * HABING_UNIT * All.UnitLength_in_cm * All.UnitLength_in_cm / All.HubbleParam / All.HubbleParam);
}

double calculate_pe_heating_rate(int i)
{
  double heating_pe;
  double G_eff, D, eps, n;

  G_eff = fmax(SphP[i].G_FUV, All.G_min);

  n = HYDROGEN_MASSFRAC * SphP[i].Density * All.cf_a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam / PROTONMASS;

  assert(SphP[i].Metallicity > 0.0);

  D = get_DGR_relative_solar(SphP[i].Metallicity);

#if(PE_MCS_POSTSHIELD == 1)
  /*Jeans length shielding */
  double N = sqrt(GAMMA * M_PI * SphP[i].Pressure * All.cf_a3inv / All.G);  // proper code M / L^2
  N *= HYDROGEN_MASSFRAC * All.UnitMass_in_g * All.HubbleParam / (All.UnitLength_in_cm * All.UnitLength_in_cm * PROTONMASS);
  G_eff *= exp((-1.33e-21) * D * N);
#endif

#ifndef PE_MCS_FIXED_EPS
  eps = fmin(0.041, 8.71e-3 * pow(n, 0.235));  // Fit to Wolfire+06 fig. 10b
#else
  eps        = All.PhotoelectricHeatingEps;
#endif

  heating_pe = 1.3e-24 * eps * D * G_eff * n;  // erg s^-1 cm^-3

  return heating_pe;
}

/*! \brief Returns DGR/DGR_solar
 *
 *  Returns dust to gas ratio, relative to solar,
 *  from Remy-Ruyer+2014. Values of solar metallicity
 *  etc. are hard coded to the values in that paper,
 *  regardless of what is adopted elsewhere, in order
 *  to get the normalisation correct. The fit adopted
 *  is the X_CO,Z broken power law from table 1. This
 *  has been converted from 12 + log(O/H) as used in
 *  the paper to Z, but gives identical results.
 *
 *
 *  \param[in] met Absolute fractional metal abundance
 *  \return dgr Dust to gas ratio relative to solar
 */
double get_DGR_relative_solar(double met)
{
#ifndef PE_MCS_LINEAR_DGR
  double logZ, y, DGR;

  logZ = log10(met / 0.014);  // logZ in units of Zsun
  if(logZ > -0.59)
    y = 2.21 - logZ;
  else
    y = 0.96 - 3.1 * logZ;

  DGR = pow(10.0, -y);  // Absolute
  DGR *= 162.0;         // Now relative to solar;
#else
  double DGR = met / 0.014;
#endif
  return DGR;
}

#endif  // PE_MCS
