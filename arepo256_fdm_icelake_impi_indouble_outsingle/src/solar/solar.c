/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/solar/solar.c
 * \date        11/2013
 * \author      Mark Lykke Winther
 * \brief       Main routine for solar simulations
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "solar.h"
#include "../allvars.h"
#include "../proto.h"
#include "string.h"

#ifdef SOLAR
void do_solar_evolution()
{
#if defined(SOLAR_RADIATIVE_TRANSFER_DIFF) || defined(SOLAR_RADIATIVE_TRANSFER_EDD)
  do_solar_rad_transfer();
#endif
}
#endif
