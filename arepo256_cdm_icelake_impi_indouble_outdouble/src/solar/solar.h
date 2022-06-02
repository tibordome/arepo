/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/solar/solar.h
 * \date        9/2021
 * \author      Mark Lykke Winther
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#ifndef SOLAR_H
#define SOLAR_H

#include "../allvars.h"

#ifdef SOLAR
void do_solar_evolution(void);
#endif

#if defined(SOLAR_RADIATIVE_TRANSFER_DIFF) || defined(SOLAR_RADIATIVE_TRANSFER_EDD)
void do_solar_rad_transfer(void);
#ifdef SOLAR_RADIATIVE_TRANSFER_EDD
void exchange_variables(void);
void point_get_center(int p, double *Center);
double lookup_opacity(double temperature, double density, double X, double Z);
double compute_diffusion_term(int particle, int fullNormalGradients, double *radval);
void radtrans_edd_approx(void);
#endif
#endif

#endif /* SOLAR_H*/
