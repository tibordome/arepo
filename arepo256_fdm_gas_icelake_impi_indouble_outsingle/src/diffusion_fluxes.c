/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/diffusion_fluxes.c
 * \date        MM/YYYY
 * \author
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

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"

#define TCOND_MEAN_MOL_WT 18.01528 /* water */

void face_get_gradients(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux)
{
#if defined(VISCOSITY) || defined(THERMAL_CONDUCTION) || defined(TRACER_DIFFUSION)

  int l;
#if defined(SECOND_DERIVATIVES) && defined(RECONSTRUCT_GRADIENTS)

#ifdef VISCOSITY
  int k;
  for(l = 0; l < 3; l++)
    for(k = 0; k < 3; k++)
      st_face->vel_grad[k][l] = 0.5 * (st_R->gradReconstruct.dvel[k][l] + st_L->gradReconstruct.dvel[k][l]);
#endif

#ifdef THERMAL_CONDUCTION
  double prefactor = TCOND_MEAN_MOL_WT * PROTONMASS / BOLTZMANN;
  for(l = 0; l < 3; l++)
    st_face->dTemp[l] =
        0.5 * prefactor *
        (st_R->gradReconstruct.dpress[l] / st_R->rho - st_R->gradReconstruct.drho[l] * st_R->press / (st_R->rho * st_R->rho) +
         st_L->gradReconstruct.dpress[l] / st_L->rho - st_L->gradReconstruct.drho[l] * st_L->press / (st_L->rho * st_L->rho));
#endif

#ifdef TRACER_DIFFUSION
  for(l = 0; l < 3; l++)
    st_face->dConservedTracer[l] = 0.5 * (st_R->gradReconstruct.dtracer[l] * st_R->rho + st_R->gradReconstruct.drho[l] * st_R->tracer +
                                          st_L->gradReconstruct.dtracer[l] * st_L->rho + st_L->gradReconstruct.drho[l] * st_L->tracer);
#endif

#else
  /*If gradients are not being reconstructed and extrapolated, choose the upwind gradients. */

  if(flux->mass > 0)
    {

#ifdef VISCOSITY
      int k;
      for(l = 0; l < 3; l++)
        for(k = 0; k < 3; k++)
          st_face->vel_grad[k][l] = st_L->grad->dvel[k][l];
#endif
#ifdef THERMAL_CONDUCTION
      double prefactor = TCOND_MEAN_MOL_WT * PROTONMASS / BOLTZMANN;
      for(l = 0; l < 3; l++)
        st_face->dTemp[l] =
            prefactor * (st_L->grad->dpress[l] / st_L->rho - st_L->grad->drho[l] * st_L->press / (st_L->rho * st_L->rho));
#endif
#ifdef TRACER_DIFFUSION
      for(l = 0; l < 3; l++)
        st_face->dConservedTracer[l] = st_L->grad->dtracer[l] * st_L->rho + st_L->grad->drho[l] * st_L->tracer;
#endif
    }
  else
    {

#ifdef VISCOSITY
      int k;
      for(l = 0; l < 3; l++)
        for(k = 0; k < 3; k++)
          st_face->vel_grad[k][l] = st_R->grad->dvel[k][l];
#endif
#ifdef THERMAL_CONDUCTION
      double prefactor = TCOND_MEAN_MOL_WT * PROTONMASS / BOLTZMANN;
      for(l = 0; l < 3; l++)
        st_face->dTemp[l] =
            prefactor * (st_R->grad->dpress[l] / st_R->rho - st_R->grad->drho[l] * st_R->press / (st_R->rho * st_R->rho));
#endif
#ifdef TRACER_DIFFUSION
      for(l = 0; l < 3; l++)
        st_face->dConservedTracer[l] = st_R->grad->dtracer[l] * st_R->rho + st_R->grad->drho[l] * st_R->tracer;
#endif
    }
#endif

#endif
}

double local_get_dynvisc_coefficient(double x, double y, double z, double rho, double press)
{
#ifdef VISCOSITY
#ifdef GLOBAL_VISCOSITY
  return All.dyn_visc;
#endif
#ifdef USE_KINEMATIC_VISCOSITY
  return All.KinematicViscosity * rho;
#endif
#ifdef ALPHA_VISCOSITY
  return get_alpha_viscosity(x, y, z, rho, press) * rho;
#endif
#endif
  return 0;
}

double local_get_bulkvisc_coefficient(double x, double y, double z, double rho, double press)
{
#ifdef VISCOSITY
#ifdef GLOBAL_VISCOSITY
  return All.bulk_visc;
#endif
#ifdef USE_KINEMATIC_VISCOSITY
  return 0;
#endif
#ifdef ALPHA_VISCOSITY
  return 0;
#endif
#endif
  return 0;
}

void face_get_viscous_fluxes(struct state_face *st_face, struct fluxes *flux, struct geometry *geo, double dyn_visc, double bulk_visc)
/*Once we have the velocity and velocity gradients (regardless of how we
obtained them, we need to project them onto the normals of each face)*/
/*The effect of moving boundaries is already taken care of in the computation
of the advectve fluxes*/
{
#ifdef VISCOSITY
  double div_vel = st_face->vel_grad[0][0] + st_face->vel_grad[1][1] + st_face->vel_grad[2][2];

  double fac1 = 4.0 / 3.0 * dyn_visc;
  double fac2 = 2.0 / 3.0 * dyn_visc;
  double fac3 = div_vel * bulk_visc;

  flux->momentum[0] +=
      (geo->nx * (fac1 * st_face->vel_grad[0][0] - fac2 * (st_face->vel_grad[1][1] + st_face->vel_grad[2][2]) + fac3) +
       geo->ny * dyn_visc * (st_face->vel_grad[0][1] + st_face->vel_grad[1][0]) +
       geo->nz * dyn_visc * (st_face->vel_grad[0][2] + st_face->vel_grad[2][0]));

  flux->momentum[1] +=
      (geo->ny * (fac1 * st_face->vel_grad[1][1] - fac2 * (st_face->vel_grad[0][0] + st_face->vel_grad[2][2]) + fac3) +
       geo->nx * dyn_visc * (st_face->vel_grad[0][1] + st_face->vel_grad[1][0]) +
       geo->nz * dyn_visc * (st_face->vel_grad[1][2] + st_face->vel_grad[2][1]));

  flux->momentum[2] +=
      (geo->nz * (fac1 * st_face->vel_grad[2][2] - fac2 * (st_face->vel_grad[0][0] + st_face->vel_grad[1][1]) + fac3) +
       geo->nx * dyn_visc * (st_face->vel_grad[0][2] + st_face->vel_grad[2][0]) +
       geo->ny * dyn_visc * (st_face->vel_grad[1][2] + st_face->vel_grad[2][1]));

  flux->energy +=
      (geo->nx *
           (st_face->velx * (fac1 * st_face->vel_grad[0][0] - fac2 * (st_face->vel_grad[1][1] + st_face->vel_grad[2][2]) + fac3) +
            st_face->vely * dyn_visc * (st_face->vel_grad[1][0] + st_face->vel_grad[0][1]) +
            st_face->velz * dyn_visc * (st_face->vel_grad[2][0] + st_face->vel_grad[0][2])) +
       geo->ny *
           (st_face->vely * (fac1 * st_face->vel_grad[1][1] - fac2 * (st_face->vel_grad[0][0] + st_face->vel_grad[2][2]) + fac3) +
            st_face->velx * dyn_visc * (st_face->vel_grad[1][0] + st_face->vel_grad[0][1]) +
            st_face->velz * dyn_visc * (st_face->vel_grad[1][2] + st_face->vel_grad[2][1])) +
       geo->nz *
           (st_face->velz * (fac1 * st_face->vel_grad[2][2] - fac2 * (st_face->vel_grad[0][0] + st_face->vel_grad[1][1]) + fac3) +
            st_face->velx * dyn_visc * (st_face->vel_grad[2][0] + st_face->vel_grad[0][2]) +
            st_face->vely * dyn_visc * (st_face->vel_grad[2][1] + st_face->vel_grad[1][2])));

  CPU_Step[CPU_VISCOUS_FLUXES] += measure_time();
#endif
}

void face_get_conduction_fluxes(struct state_face *st_face, struct fluxes *flux, struct geometry *geo)
{
#ifdef THERMAL_CONDUCTION
  flux->energy += -All.ThermalConductivity * (st_face->dTemp[0] * geo->nx + st_face->dTemp[1] * geo->ny + st_face->dTemp[2] * geo->nz);
#endif
}

void face_get_scalar_diffusion_fluxes(struct state_face *st_face, struct fluxes *flux, struct geometry *geo)
{
#ifdef TRACER_DIFFUSION

  flux->tracer = All.TracerDiffusivity * (st_face->dConservedTracer[0] * geo->nx + st_face->dConservedTracer[1] * geo->ny +
                                          st_face->dConservedTracer[2] * geo->nz);

#endif
}

void face_extrapolate_viscous_kick(struct state *delta, const struct state *st, double dyn_visc, double bulk_visc)
{
#if defined(VISCOSITY) && defined(SECOND_DERIVATIVES)
  /* here we kick the state velocity of the cell by applying the
  viscous forces as they appear in the Navier-Stokes equations i.e.
  as a source term */

  double laplacian_vel[3], grad_div_vel[3];
  struct hessian_data *hessian;

#ifndef MUSCL_HANCOCK /* we use the default Runge Kutta time integration scheme */
  double dt_half = st->dtExtrapolation;
#else
  double dt_half = st->dt_half;
#endif

  /*  dyn_visc = st->dyn_visc;
     bulk_visc = st->bulk_visc; */

  hessian = st->hessian;

  laplacian_vel[0] = hessian->ddvelx[0][0] + hessian->ddvelx[1][1] + hessian->ddvelx[2][2];
  laplacian_vel[1] = hessian->ddvely[0][0] + hessian->ddvely[1][1] + hessian->ddvely[2][2];
  laplacian_vel[2] = hessian->ddvelz[0][0] + hessian->ddvelz[1][1] + hessian->ddvelz[2][2];

  grad_div_vel[0] = hessian->ddvelx[0][0] + hessian->ddvely[1][0] + hessian->ddvelz[2][0];
  grad_div_vel[1] = hessian->ddvelx[0][1] + hessian->ddvely[1][1] + hessian->ddvelz[2][1];
  grad_div_vel[2] = hessian->ddvelx[0][2] + hessian->ddvely[1][2] + hessian->ddvelz[2][2];

  delta->velx += dt_half * (dyn_visc * laplacian_vel[0] + (bulk_visc + dyn_visc / 3.0 * grad_div_vel[0])) / st->rho;

  delta->vely += dt_half * (dyn_visc * laplacian_vel[1] + (bulk_visc + dyn_visc / 3.0 * grad_div_vel[1])) / st->rho;

  delta->velz += dt_half * (dyn_visc * laplacian_vel[2] + (bulk_visc + dyn_visc / 3.0 * grad_div_vel[2])) / st->rho;

  delta->press += dt_half * GAMMA_MINUS1 * get_viscous_dissipation(st, dyn_visc, bulk_visc);
#endif
}

double get_viscous_dissipation(const struct state *st, double dyn_visc, double bulk_visc)
{
#if defined(VISCOSITY) && defined(SECOND_DERIVATIVES)
  /*compute entries of viscous stress tensor*/
  double tau[3][3];
  get_viscous_stress_tensor(tau, st, dyn_visc, bulk_visc);

  return tau[0][0] * st->grad->dvel[0][0] + tau[0][1] * st->grad->dvel[0][1] + tau[0][2] * st->grad->dvel[0][2] +
         tau[1][0] * st->grad->dvel[1][0] + tau[1][1] * st->grad->dvel[1][1] + tau[1][2] * st->grad->dvel[1][2] +
         tau[2][0] * st->grad->dvel[2][0] + tau[2][1] * st->grad->dvel[2][1] + tau[2][2] * st->grad->dvel[2][2];
#endif
  return 0;
}

void get_viscous_stress_tensor(double tau[3][3], const struct state *st, double dyn_visc, double bulk_visc)
{
#if defined(VISCOSITY) && defined(SECOND_DERIVATIVES)
  double div_vel = st->grad->dvel[0][0] + st->grad->dvel[1][1] + st->grad->dvel[2][2];

  double fac1 = 4.0 / 3.0 * dyn_visc;
  double fac2 = 2.0 / 3.0 * dyn_visc;
  double fac3 = div_vel * bulk_visc;

  tau[0][0] = fac1 * st->grad->dvel[0][0] - fac2 * (st->grad->dvel[1][1] + st->grad->dvel[2][2]) + fac3;
  tau[0][1] = dyn_visc * (st->grad->dvel[0][1] + st->grad->dvel[1][0]);
  tau[0][2] = dyn_visc * (st->grad->dvel[2][0] + st->grad->dvel[0][2]);

  tau[1][0] = dyn_visc * (st->grad->dvel[0][1] + st->grad->dvel[1][0]);
  tau[1][1] = fac1 * st->grad->dvel[1][1] - fac2 * (st->grad->dvel[2][2] + st->grad->dvel[0][0]) + fac3;
  tau[1][2] = dyn_visc * (st->grad->dvel[1][2] + st->grad->dvel[2][1]);

  tau[2][0] = dyn_visc * (st->grad->dvel[2][0] + st->grad->dvel[0][2]);
  tau[2][2] = fac1 * st->grad->dvel[2][2] - fac2 * (st->grad->dvel[0][0] + st->grad->dvel[1][1]) + fac3;
  tau[2][1] = dyn_visc * (st->grad->dvel[1][2] + st->grad->dvel[2][1]);
#endif
}

double get_alpha_viscosity(double x, double y, double z, double rho, double press)
{
#if defined(VISCOSITY) && defined(ALPHA_VISCOSITY)
  if(ALPHA_VISCOSITY == 1)  // alpha-viscosity based on pressure
    return All.AlphaCoefficient * press * rho;
  else if(ALPHA_VISCOSITY == 2)  // alpha-viscosity based on position
    {
#ifdef CIRCUMSTELLAR
      return get_circumstellar_alpha_viscosity(x, y, z, rho, press);
#else
      return 0;
#endif
    }
#endif
  return 0;
}
