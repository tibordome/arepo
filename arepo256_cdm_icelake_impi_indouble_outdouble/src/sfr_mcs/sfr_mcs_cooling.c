/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sfr_mcs/sfr_mcs.c
 * \date        04/2018
 * \author     	Matthew C Smith
 * \brief
 * \details     Allows use of GFM_COOLING_METAL without rest of GFM. This is primarily
                to avoid the issue that currently SFR_MCS only treats metals as one
                species, so we rescale the helium abundance with solar to get
                the hydrogen abundance correct. To avoid messy ifdefs in cooling/cooling.c
                we replicate several functions here. This mainly involves passing static
                variables from cooling/cooling.c through to local copies here. This
                is probably bad practice and should be updated, but it will suffice
                for now.
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"

#if defined(SFR_MCS) && defined(GFM_COOLING_METAL)

#ifdef GFM_AGN_RADIATION
#error "Cannot currently use SFR_MCS, GFM_COOLING_METAL and GFM_AGN_RADIATION"
#endif

#ifdef GFM_UVB_CORRECTIONS
#error "Cannot currently use SFR_MCS, GFM_COOLING_METAL and GFM_UVB_CORRECTIONS"
#endif

#ifdef RADCOOL
#error "Cannot currently use SFR_MCS, GFM_COOLING_METAL and RADCOOL"
#endif

#define HELIUM_PRIM_TO_SOLAR 0.043

static double MetallicityFloor, LogMetallicityFloorInSolar;

/* Passes static variables from cooling.c through */
void metal_cool_mcs_init(double localMetallicityFloor, double localLogMetallicityFloorInSolar)
{
  MetallicityFloor           = localMetallicityFloor;
  LogMetallicityFloorInSolar = localLogMetallicityFloorInSolar;
}

#ifdef UVB_SELF_SHIELDING
/* This replaces the definition in cooling.c for the sake of avoiding ifdefs everywhere
But pc is static to cooling.c, so we pass it by reference */
void update_radiation_state(MyFloat rho, MyFloat metallicity, PhotoCurrent *localpc)
{
  double xh, met, nH, selfshielding_factor;

  /* Get hydrogen mass fraction by scaling helium abundance as solar, up to one solar metallicity */
  xh = HYDROGEN_MASSFRAC - (HELIUM_PRIM_TO_SOLAR / GFM_SOLAR_METALLICITY + 1.0) * fmin(metallicity, 10.0 * GFM_SOLAR_METALLICITY);
  nH = xh * rho / PROTONMASS * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;

  selfshielding_factor = calc_self_shielding_factor(nH);
  localpc->gJH0 *= selfshielding_factor;
  localpc->gJHe0 *= selfshielding_factor;
  localpc->gJHep *= selfshielding_factor;
  localpc->epsH0 *= selfshielding_factor;
  localpc->epsHe0 *= selfshielding_factor;
  localpc->epsHep *= selfshielding_factor;
}
#endif

/* This replaces the definition in cooling.c for the sake of avoiding ifdefs everywhere
But gs is static to cooling.c, so we pass it by reference */
void update_gas_state(MyFloat rho, MyFloat metallicity, GasState *localgs)
{
  if((metallicity > MetallicityFloor) && (metallicity < 10.0 * GFM_SOLAR_METALLICITY))
    localgs->log_MetallicityInSolar = log10(metallicity / GFM_SOLAR_METALLICITY);
  else if(metallicity > 10.0 * GFM_SOLAR_METALLICITY)
    localgs->log_MetallicityInSolar = 10.0;
  else
    localgs->log_MetallicityInSolar = LogMetallicityFloorInSolar;

  /* Get hydrogen mass fraction by scaling helium abundance as solar, up to one solar metallicity */
  localgs->XH =
      HYDROGEN_MASSFRAC - (HELIUM_PRIM_TO_SOLAR / GFM_SOLAR_METALLICITY + 1.0) * fmin(metallicity, 10.0 * GFM_SOLAR_METALLICITY);
  localgs->yhelium = (1 - localgs->XH - metallicity) / (4. * localgs->XH);

  /* convert to physical cgs units, is already physical at this point */
  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
  localgs->log_HydrogenNumberDensity = log10(localgs->XH * rho / PROTONMASS);
#ifdef GFM_LAMBDA
  localgs->partindex = -1;
#endif
}
#endif