/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/MRT/RT_riemann_rosunov.c
 * \date        06/2018
 * \author      Rahul Kannan
 * \brief       Rosunov Riemann solver that takes the mesh velocity into account while calculating
 * \details     the eigenvalues
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../allvars.h"
#include "../proto.h"
#include "../voronoi.h"

#ifdef MRT_RIEMANN_ROSUNOV

static void rosunov_get_fluxes_from_state_RT(struct state *st, struct fluxes *flux, double vel_face_x)
{
  for(int num1 = 0; num1 < MRT_BINS; num1++)
    {
#ifndef MRT_COMOVING
      flux->DensPhot[num1] = st->RT_F[num1][0] - st->DensPhot[num1] * vel_face_x;

      flux->RT_F[num1][0] = c_internal_units * c_internal_units * st->PT[num1][0][0] - st->RT_F[num1][0] * vel_face_x;
      flux->RT_F[num1][1] = c_internal_units * c_internal_units * st->PT[num1][1][0] - st->RT_F[num1][1] * vel_face_x;
      flux->RT_F[num1][2] = c_internal_units * c_internal_units * st->PT[num1][2][0] - st->RT_F[num1][2] * vel_face_x;
#else
      flux->DensPhot[num1] = st->RT_F[num1][0] + st->velx * st->DensPhot[num1];

      flux->RT_F[num1][0] = c_internal_units * c_internal_units * st->PT[num1][0][0] + st->velx * st->RT_F[num1][0];
      flux->RT_F[num1][1] = c_internal_units * c_internal_units * st->PT[num1][1][0] + st->velx * st->RT_F[num1][1];
      flux->RT_F[num1][2] = c_internal_units * c_internal_units * st->PT[num1][2][0] + st->velx * st->RT_F[num1][2];

#endif
    }
}

static void get_rosunov_fluxes_RT(const struct state *st_L, const struct state *st_R, const struct fluxes *flux_L,
                                  const struct fluxes *flux_R, struct fluxes *rosunov_flux, double S_plus, double S_minus)
{
  if(S_minus >= 0.0)
    {
      for(int num1 = 0; num1 < MRT_BINS; num1++)
        {
          rosunov_flux->DensPhot[num1] = flux_L->DensPhot[num1];
          rosunov_flux->RT_F[num1][0]  = flux_L->RT_F[num1][0];
          rosunov_flux->RT_F[num1][1]  = flux_L->RT_F[num1][1];
          rosunov_flux->RT_F[num1][2]  = flux_L->RT_F[num1][2];
        }
    }
  else if(S_plus <= 0.0)
    {
      for(int num1 = 0; num1 < MRT_BINS; num1++)
        {
          rosunov_flux->DensPhot[num1] = flux_R->DensPhot[num1];
          rosunov_flux->RT_F[num1][0]  = flux_R->RT_F[num1][0];
          rosunov_flux->RT_F[num1][1]  = flux_R->RT_F[num1][1];
          rosunov_flux->RT_F[num1][2]  = flux_R->RT_F[num1][2];
        }
    }
  else
    {
      for(int num1 = 0; num1 < MRT_BINS; num1++)
        {
          rosunov_flux->DensPhot[num1] = (S_plus * flux_L->DensPhot[num1] - S_minus * flux_R->DensPhot[num1] +
                                          S_plus * S_minus * (st_R->DensPhot[num1] - st_L->DensPhot[num1])) /
                                         (S_plus - S_minus);

          rosunov_flux->RT_F[num1][0] = (S_plus * flux_L->RT_F[num1][0] - S_minus * flux_R->RT_F[num1][0] +
                                         S_plus * S_minus * (st_R->RT_F[num1][0] - st_L->RT_F[num1][0])) /
                                        (S_plus - S_minus);
          rosunov_flux->RT_F[num1][1] = (S_plus * flux_L->RT_F[num1][1] - S_minus * flux_R->RT_F[num1][1] +
                                         S_plus * S_minus * (st_R->RT_F[num1][1] - st_L->RT_F[num1][1])) /
                                        (S_plus - S_minus);
          rosunov_flux->RT_F[num1][2] = (S_plus * flux_L->RT_F[num1][2] - S_minus * flux_R->RT_F[num1][2] +
                                         S_plus * S_minus * (st_R->RT_F[num1][2] - st_L->RT_F[num1][2])) /
                                        (S_plus - S_minus);
        }
    }
}

double godunov_flux_3d_rosunov_RT(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux,
                                  double vel_face_x)
{
  int calc_flag = 1;
  for(int num1 = 0; num1 < MRT_BINS; num1++)
    {
      if(st_L->DensPhot[num1] >= 0.0 && st_R->DensPhot[num1] >= 0.0)
        calc_flag *= 1;
      else
        calc_flag *= 0;
    }

  if(calc_flag)
    {
      double S_plus, S_minus;
      struct fluxes flux_L, flux_R;

      /*set minimum and maximum einvelocities of the system*/
      S_plus  = fmax((-vel_face_x + c_internal_units), (-vel_face_x - c_internal_units));
      S_minus = fmin((-vel_face_x + c_internal_units), (-vel_face_x - c_internal_units));

      /* compute fluxes for the left and right states */
      rosunov_get_fluxes_from_state_RT(st_L, &flux_L, vel_face_x);
      rosunov_get_fluxes_from_state_RT(st_R, &flux_R, vel_face_x);

      /* compute Rosunov fluxes */
      get_rosunov_fluxes_RT(st_L, st_R, &flux_L, &flux_R, flux, S_plus, S_minus);
    }
  else
    {
      terminate("Ngamma is negative\n");
    }

  return 0;
}

#endif
