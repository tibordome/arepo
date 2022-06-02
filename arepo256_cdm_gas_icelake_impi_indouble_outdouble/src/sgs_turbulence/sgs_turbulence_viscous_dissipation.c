/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sgs_turbulence/sgs_turbulence_viscous_dissipation.c
 * \date        07/2017
 * \author      S. Jacob
 * \brief       Functions to compute the viscous dissipation of SGS turbulence energy
 * \details		This folder contains a wrapper function for the viscous dissipation of
 * 				SGS turbulence energy. The function calls the source term for the
 * 				chosen closure scheme.
 *
 * 				contains functions:
 * 					void do_sgs_turbulence_viscous_dissipation(void)
 *
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"
#include "string.h"

#ifdef SGS_TURBULENCE_VISCOUS_DISSIPATION

/*! \brief Source function for viscous dissipation
 *
 *  This function calls the source function for viscous dissipation of SGS turbulence
 *  energy for the used closure scheme (for now only eddy viscosity closure).
 *  Then the primitive variables are updated.
 *
 *  \return void
 */
void do_sgs_turbulence_viscous_dissipation(void)
{
#ifdef SGS_TURBULENCE_EDDY_VISCOSITY
  do_sgs_turbulence_viscous_dissipation_eddy_viscosity();
#endif /* #ifdef SGS_TURBULENCE_EDDY_VISCOSITY */
  update_primitive_variables();
}

#endif /*#ifdef SGS_TURBULENCE_VISCOUS_DISSIPATION*/
