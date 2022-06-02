/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 * 
 * \file        src/SMUGGLE/cooling_molecules.c
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

#ifdef SMUGGLE_MOLEC_COOLING

double calc_self_shielding_factor_molecular_cooling(double nH, char J_UV)
{
  double factor;

  /* with no UV background molecules can always form; in other words there is no
   * UV radiation and molecular cooling is maximal at any density. Only when UV
   * radiation is present gas must be dense enough to cool via molecules */
  if(J_UV == 0)
    factor = 0.0;
  else
    factor = calc_self_shielding_factor(nH);

  return factor;
}

double get_MolecularCoolingRate(double logT, double nH, double logZinSolar, char J_UV)
{
  double cool_rate = 0.0;
  double T         = pow(10, logT);
  /* cap metallicity to a maximum of 1.5 Z_sun. Usually if lookup table is out
   * of bounds, interpolation is disabled. This mimicks it. */
  double metallicity = fmin(pow(10., logZinSolar), 1.5);

  cool_rate = 2.896e-26 * (1 + metallicity) / (1 + 0.00143 * nH) *
              (0.001 + (0.1 * nH + metallicity * metallicity) / (1. + nH) + 0.09 * nH / (1. + 0.1 * nH)) *
              exp(-(T / 1.58e5) * (T / 1.58e5)) / (pow(T / 125.215, -4.9202) + pow(T / 1349.86, -1.7288) + pow(T / 6450.06, -0.3075));

#ifdef UVB_SELF_SHIELDING
  double sshield = 1.0 - calc_self_shielding_factor_molecular_cooling(nH, J_UV);

  if(sshield < 0.0)
    terminate("Negative self-shielding correction");

  cool_rate *= sshield;
#endif

  return cool_rate;
}
#endif
