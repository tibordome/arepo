/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sgs_turbulence/sgs_turbulence_turbulent_production.c
 * \date        07/2017
 * \author      S. Jacob
 * \brief       Functions to compute the production of SGS turbulence energy by the turbulent cascade.
 * \details		This folder contains a wrapper function for the SGS turbulence energy production by the
 * 				turbulent cascade. The function calls the source term for the chosen closure scheme.
 *
 * 				contains functions:
 * 					void do_sgs_energy_turbulent_production(void)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"
#include "string.h"

#ifdef SGS_TURBULENCE_TURBULENT_PRODUCTION

/*! \brief Source function for turbulent energy production
 *
 *  This function calls the source function for turbulent SGS turbulence energy
 *  production for the used closure scheme (for now only eddy viscosity closure).
 *  Then the primitive variables are updated.
 *
 *  \return void
 */
void do_sgs_energy_turbulent_production(void)
{
#ifdef SGS_TURBULENCE_EDDY_VISCOSITY
  do_sgs_energy_turbulent_production_eddy_viscosity();
#endif /* #ifdef SGS_TURBULENCE_EDDY_VISCOSITY */
  update_primitive_variables();
}

#endif /*#ifdef SGS_TURBULENCE_TURBULENT_PRODUCTION*/
