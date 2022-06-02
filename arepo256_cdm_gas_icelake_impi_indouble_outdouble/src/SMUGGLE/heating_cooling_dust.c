/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 * 
 * \file        src/SMUGGLE/heating_cooling_dust.c
 * \date        03/2020
 * \author      Federico Marinacci
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

#ifdef SMUGGLE_DUST_HEATING_COOLING

double get_DustHeatingCoolingRate(double logT, double XH, double logZinSolar)
{
  double cool_rate   = 0.0;
  double Tdust       = 30.;
  double T           = pow(10, logT);
  double metallicity = pow(10., logZinSolar);

  cool_rate = 1.2e-32 * metallicity * sqrt(T) * (T - Tdust) * (1. - 0.8 * exp(-75. / T)) / (XH * XH);

  return cool_rate;
}
#endif
