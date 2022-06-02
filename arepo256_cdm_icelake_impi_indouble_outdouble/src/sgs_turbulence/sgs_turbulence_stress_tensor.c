
/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sgs_turbulence/sgs_turbulence.c
 * \date        07/2017
 * \author      S. Jacob
 * \brief       Functions to include the SGS turbulence stress tensor
 * \details     This file contains functions that compute the effects of the SGS
 * 				turbulence stress tensor. This includes functions for the gradient
 * 				reconstruction of the velocity at the interface, the computation of
 * 				a viscous kick and the viscous fluxes. Most of these functions
 * 				have to be implemented for each closure scheme individually.
 *
 * 				contains functions:
 *					void face_get_reconstructed_gradients(struct state *st_L, struct state *st_R, struct state_face
 **st_face, struct fluxes *flux) void face_add_sgs_stress_tensor_fluxes(const struct state_face *st, const struct geometry *geom,
 *struct fluxes *flux) void face_extrapolate_sgs_turbulence_viscous_kick(struct state *delta, const struct state *st)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"
#include "string.h"

#ifdef SGS_TURBULENCE_STRESS_TENSOR
/*! \brief	Compute velocity gradients at the interface
 *
 * 	This function determines the velocity gradients at the
 * 	interface. If the gradients are reconstructed from the
 * 	centre of the cell to the interface with second
 * 	derivatives, the average of the reconstructed gradients
 * 	is used. Otherwise, the upwind value is taken.
 *
 * 	\param[in] st_L Left state of the interface
 * 	\param[in] st_R Right state of the interface
 * 	\param[in] st_face State at the interface
 * 	\param[in] flux Flux over the interface
 *
 * 	\return void
 */
void face_get_reconstructed_gradients(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux)
{
#if defined(SECOND_DERIVATIVES) && defined(RECONSTRUCT_GRADIENTS)
  for(int l = 0; l < 3; l++)
    for(int k = 0; k < 3; k++)
      st_face->vel_grad[k][l] = 0.5 * (st_R->gradReconstruct.dvel[k][l] + st_L->gradReconstruct.dvel[k][l]);
#else
  /*If gradients are not being reconstructed and extrapolated, choose the upwind gradients. */

  if(flux->mass > 0)
    {
      for(int l = 0; l < 3; l++)
        for(int k = 0; k < 3; k++)
          st_face->vel_grad[k][l] = st_L->grad->dvel[k][l];
    }
  else
    {
      for(int l = 0; l < 3; l++)
        for(int k = 0; k < 3; k++)
          st_face->vel_grad[k][l] = st_R->grad->dvel[k][l];
    }
#endif
}

/*! \brief Add viscous flux with chosen closure scheme
 *
 * 	This function selects the function for adding the flux from
 * 	the SGS turbulence stress tensor that corresponds to the
 * 	chosen closure scheme. These functions use the time and space extra-
 * 	polated states on the left and right of the interface to compute
 * 	the flux.
 *
 * 	\param[in] st State at the interface
 * 	\param[in] geom Interface geometry
 * 	\param[in] flux Flux over the interface
 *
 *  \return void
 */
void face_add_sgs_stress_tensor_fluxes(const struct state_face *st, const struct geometry *geom, struct fluxes *flux)
{
#ifdef SGS_TURBULENCE_EDDY_VISCOSITY
  face_add_sgs_stress_tensor_eddy_viscosity_fluxes(st, geom, flux);
#endif /* #ifdef SGS_TURBULENCE_EDDY_VISCOSITY */
}

/*! \brief Add viscous kick to time extrapolation
 *
 * 	This function adds a viscous kick due to the turbulence stress tensor
 * 	to the time extrapolation of the primitive variables. This function
 * 	selects the chosen closure scheme.
 *
 * 	\param[in] delta State that contains the changes of the variables due to the time extrapolation
 * 	\param[in] st Left or right state of the interface
 *
 * 	\return void
 */
void face_extrapolate_sgs_turbulence_viscous_kick(struct state *delta, const struct state *st)
{
#if defined(SECOND_DERIVATIVES) && defined(RECONSTRUCT_GRADIENTS)
#ifdef SGS_TURBULENCE_EDDY_VISCOSITY
  face_extrapolate_sgs_turbulence_viscous_kick_eddy_viscosity(delta, st);
#endif /* #ifdef SGS_TURBULENCE_EDDY_VISCOSITY */
#endif /* #if defined(SECOND_DERIVATIVES) && defined(RECONSTRUCT_GRADIENTS) */
}

#endif /* #ifdef SGS_TURBULENCE_STRESS_TENSOR */
