/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sgs_turbulence/sgs_turbulence.h
 * \date        07/2017
 * \author		S. Jacob
 * \brief		Functions for eddy viscosity closure
 * \details		This file contains all functions that are specific for the eddy
 * 				viscosity closure. This includes functions to compute the
 * 				turbulence stress tensor, turbulent production of SGS turbulence
 * 				energy and viscous dissipation.
 *
 * 				contains functions:
 * 					double calculate_eddy_viscosity(double filter_scale, double rho, double sgstPressure)
 *
 * 					void face_get_eddy_viscosity(const struct state *st_L, const struct state *st_R, struct
 *state_face *st_face, const struct fluxes *flux) void face_get_sgs_stress_tensor_eddy_viscosity(double tau[3][3], const double
 *vel_grad[3][3], double eddy_visc) void state_get_sgs_stress_tensor_eddy_viscosity(double tau[3][3], MySingle vel_grad[3][3], double
 *eddy_visc) void face_add_sgs_stress_tensor_eddy_viscosity_fluxes(const struct state_face *st, const struct geometry *geom, struct
 *fluxes *flux) double get_viscous_dissipation_eddy_viscosity(const struct state *st, double eddy_visc) void
 *face_extrapolate_sgs_turbulence_viscous_kick_eddy_viscosity(struct state *delta, const struct state *st)
 *
 *					void do_sgs_energy_turbulent_production_eddy_viscosity(void)
 *
 *					void do_sgs_turbulence_viscous_dissipation_eddy_viscosity(void)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"
#include "string.h"

#ifdef SGS_TURBULENCE_EDDY_VISCOSITY

/*! \brief Compute eddy viscosity
 *
 * 	This function computes the eddy viscosity from the primitive
 * 	variables.
 *
 * 	\param[in] filter_scale Effective filter scale
 * 	\param[in] rho Density
 * 	\param[in] sgstPressure SGS turbulence pressure
 *
 * 	\return Eddy viscosity, in code units
 */
double calculate_eddy_viscosity(double filter_scale, double rho, double sgstPressure)
{
  return All.SgsTConst.EddyViscosityParameterCnu * filter_scale * sqrt(rho * sgstPressure / (All.SgsTConst.GammaSgsT - 1.));
}

/* functions to include the turbulence stress tensor */
#ifdef SGS_TURBULENCE_STRESS_TENSOR

/*! \brief Determine eddy viscosity at the interface
 *
 * 	This function determines the eddy viscosity at the interface
 * 	from the filter scale and the SGS turbulence pressure
 * 	at the interface.
 *
 * 	\param[in] st_L Left state of the interface
 * 	\param[in] st_R Right state of the interface
 * 	\param[in] st_face State at the interface
 * 	\param[in] flux Flux over the interface
 *
 * 	\return void
 */
void face_get_eddy_viscosity(const struct state *st_L, const struct state *st_R, struct state_face *st_face, const struct fluxes *flux)
{
#ifdef SGS_TURBULENCE_CONSTANT_EDDY_VISCOSITY
  st_face->eddy_visc = All.SgsTConst.ConstantEddyViscosity;
#else  /* #ifdef SGS_TURBULENCE_CONSTANT_EDDY_VISCOSITY */
  double filter_scale = calculate_filter_scale_at_face(st_L, st_R);

  /*take upwind value for sgstPressure at interface*/
  if(flux->mass > 0)
    st_face->sgstPressure = st_L->sgstPressure;
  else
    st_face->sgstPressure = st_R->sgstPressure;

  st_face->eddy_visc = calculate_eddy_viscosity(filter_scale, st_face->rho, st_face->sgstPressure);
#endif /* #ifdef SGS_TURBULENCE_CONSTANT_EDDY_VISCOSITY #else */
}

/*! \brief Determine the stress tensor at the interface
 *
 * 	This function computes the SGS turbulence stress tensor at
 * 	the interface for the eddy viscosity closure.
 *
 * 	\param[in] tau Traceless part of the eddy viscosity stress tensor
 * 				   that is computed in this function
 * 	\param[in] vel_grad Velocity gradients
 * 	\param[in] eddy_visc Eddy viscosity at the interface
 *
 *  \return void
 */
void face_get_sgs_stress_tensor_eddy_viscosity(double tau[3][3], const double vel_grad[3][3], double eddy_visc)
{
  tau[0][0] = 2. / 3. * eddy_visc * (2. * vel_grad[0][0] - vel_grad[1][1] - vel_grad[2][2]);
  tau[0][1] = eddy_visc * (vel_grad[0][1] + vel_grad[1][0]);
  tau[0][2] = eddy_visc * (vel_grad[2][0] + vel_grad[0][2]);

  tau[1][0] = eddy_visc * (vel_grad[0][1] + vel_grad[1][0]);
  tau[1][1] = 2. / 3. * eddy_visc * (2. * vel_grad[1][1] - vel_grad[0][0] - vel_grad[2][2]);
  tau[1][2] = eddy_visc * (vel_grad[1][2] + vel_grad[2][1]);

  tau[2][0] = eddy_visc * (vel_grad[2][0] + vel_grad[0][2]);
  tau[2][1] = eddy_visc * (vel_grad[1][2] + vel_grad[2][1]);
  tau[2][2] = 2. / 3. * eddy_visc * (2. * vel_grad[2][2] - vel_grad[0][0] - vel_grad[1][1]);
}

/*! \brief Determine the stress tensor of the left/right state
 *
 * 	This function computes the stress tensor of the left or right
 * 	state of an interface.
 *
 * 	\param[in] tau Stress tensor that is computed
 * 	\param[in] vel_grad Velocity gradients
 * 	\param[in] eddy_visc Eddy viscosity of a state
 *
 *	\return void
 */
void state_get_sgs_stress_tensor_eddy_viscosity(double tau[3][3], MySingle vel_grad[3][3], double eddy_visc)
{
  tau[0][0] = 2. / 3. * eddy_visc * (2. * vel_grad[0][0] - vel_grad[1][1] - vel_grad[2][2]);
  tau[0][1] = eddy_visc * (vel_grad[0][1] + vel_grad[1][0]);
  tau[0][2] = eddy_visc * (vel_grad[2][0] + vel_grad[0][2]);

  tau[1][0] = eddy_visc * (vel_grad[0][1] + vel_grad[1][0]);
  tau[1][1] = 2. / 3. * eddy_visc * (2. * vel_grad[1][1] - vel_grad[0][0] - vel_grad[2][2]);
  tau[1][2] = eddy_visc * (vel_grad[1][2] + vel_grad[2][1]);

  tau[2][0] = eddy_visc * (vel_grad[2][0] + vel_grad[0][2]);
  tau[2][1] = eddy_visc * (vel_grad[1][2] + vel_grad[2][1]);
  tau[2][2] = 2. / 3. * eddy_visc * (2. * vel_grad[2][2] - vel_grad[0][0] - vel_grad[1][1]);
}

/* \brief Add flux due to the SGS turbulence stress tensor
 *
 * 	This function adds the flux due to the stress tensor
 * 	to the hydro fluxes. It is specific for the
 * 	eddy viscosity closure.
 *
 * 	\param[in] st State at the interface
 * 	\param[in] geom Interface geometry
 * 	\param[in] flux Flux over the interface
 *
 *  \return void
 */
void face_add_sgs_stress_tensor_eddy_viscosity_fluxes(const struct state_face *st, const struct geometry *geom, struct fluxes *flux)
{
  /*compute subgridscale turbulence stress tensor*/
  double tau[3][3];
  face_get_sgs_stress_tensor_eddy_viscosity(tau, st->vel_grad, st->eddy_visc);

  flux->momentum[0] -= geom->nx * tau[0][0] + geom->ny * tau[1][0] + geom->nz * tau[2][0];
  flux->momentum[1] -= geom->nx * tau[0][1] + geom->ny * tau[1][1] + geom->nz * tau[2][1];
  flux->momentum[2] -= geom->nx * tau[0][2] + geom->ny * tau[1][2] + geom->nz * tau[2][2];

  flux->energy -= geom->nx * (tau[0][0] * st->velx + tau[0][1] * st->vely + tau[0][2] * st->velz) +
                  geom->ny * (tau[1][0] * st->velx + tau[1][1] * st->vely + tau[1][2] * st->velz) +
                  geom->nz * (tau[2][0] * st->velx + tau[2][1] * st->vely + tau[2][2] * st->velz);
}

/*! \brief Compute viscous kick for thermal pressure
 *
 * The SGS turbulence stress tensor leads to viscous dissipation
 * that modifies the time extrapolation of the thermal pressure.
 * This contribution to the viscous kick is computed in this
 * function.
 *
 * \param[in] st Left/Right state of the interface
 * \param[in] eddy_visc Eddy viscosity
 *
 * \return Viscous kick for thermal pressure
 */
double get_viscous_dissipation_eddy_viscosity(const struct state *st, double eddy_visc)
{
  double tau[3][3];
  state_get_sgs_stress_tensor_eddy_viscosity(tau, st->grad->dvel, eddy_visc);

  /* tau_ij * partial_j v_i*/
  return tau[0][0] * st->grad->dvel[0][0] + tau[0][1] * st->grad->dvel[0][1] + tau[0][2] * st->grad->dvel[0][2] +
         tau[1][0] * st->grad->dvel[1][0] + tau[1][1] * st->grad->dvel[1][1] + tau[1][2] * st->grad->dvel[1][2] +
         tau[2][0] * st->grad->dvel[2][0] + tau[2][1] * st->grad->dvel[2][1] + tau[2][2] * st->grad->dvel[2][2];
}

/*! \brief Viscous kick for eddy viscosity
 *
 * 	This function computes the viscous kick due to the SGS turbulence
 * 	stress tensor that modifies the time extrapolation of the left and
 * 	right state of an interface.
 *
 * 	\param[in] delta State that contains changes of the variables due to the time extrapolation
 * 	\param[in] st Left/right state of an interface
 *
 * 	\return void
 */
void face_extrapolate_sgs_turbulence_viscous_kick_eddy_viscosity(struct state *delta, const struct state *st)
{
#if defined(SECOND_DERIVATIVES) && defined(RECONSTRUCT_GRADIENTS)
  /*similar to function face_extrapolate_viscouskick() in diffusion_fluxes.c (same if eddy viscosity is constant) */

  /*determine time step*/
#ifndef MUSCL_HANCOCK /* we use the default Runge Kutta time integration scheme */
  double dt_half = st->dtExtrapolation;
#else  /* #ifndef MUSCL_HANCOCK */
  double dt_half      = st->dt_half;
#endif /* #ifndef MUSCL_HANCOCK #else */

  /*calculate eddy viscosity*/
#ifdef SGS_TURBULENCE_CONSTANT_EDDY_VISCOSITY
  double eddy_visc = All.SgsTConst.ConstantEddyViscosity;
#else  /* #ifdef SGS_TURBULENCE_CONSTANT_EDDY_VISCOSITY */
  double filter_scale = calculate_filter_scale(st->volume);
  double eddy_visc    = calculate_eddy_viscosity(filter_scale, st->rho, st->sgstPressure);
#endif /* #ifdef SGS_TURBULENCE_CONSTANT_EDDY_VISCOSITY #else */

  /*compute required derivatives*/
  struct hessian_data *hessian;
  hessian = st->hessian;

  double laplacian_vel[3], grad_div_vel[3];

  laplacian_vel[0] = hessian->ddvelx[0][0] + hessian->ddvelx[1][1] + hessian->ddvelx[2][2];
  laplacian_vel[1] = hessian->ddvely[0][0] + hessian->ddvely[1][1] + hessian->ddvely[2][2];
  laplacian_vel[2] = hessian->ddvelz[0][0] + hessian->ddvelz[1][1] + hessian->ddvelz[2][2];

  grad_div_vel[0] = hessian->ddvelx[0][0] + hessian->ddvely[1][0] + hessian->ddvelz[2][0];
  grad_div_vel[1] = hessian->ddvelx[0][1] + hessian->ddvely[1][1] + hessian->ddvelz[2][1];
  grad_div_vel[2] = hessian->ddvelx[0][2] + hessian->ddvely[1][2] + hessian->ddvelz[2][2];

  /*compute viscous kicks*/
  delta->velx += dt_half * (eddy_visc * laplacian_vel[0] + eddy_visc / 3.0 * grad_div_vel[0]) / st->rho;

  delta->vely += dt_half * (eddy_visc * laplacian_vel[1] + eddy_visc / 3.0 * grad_div_vel[1]) / st->rho;

  delta->velz += dt_half * (eddy_visc * laplacian_vel[2] + eddy_visc / 3.0 * grad_div_vel[2]) / st->rho;

  delta->press += dt_half * GAMMA_MINUS1 * get_viscous_dissipation_eddy_viscosity(st, eddy_visc);

#ifndef SGS_TURBULENCE_CONSTANT_EDDY_VISCOSITY
  /*need some extra terms if viscosity is not constant */
  double fac =
      All.SgsTConst.EddyViscosityParameterCnu * filter_scale * sqrt(st->sgstPressure / (st->rho * (All.SgsTConst.GammaSgsT - 1.0)));
  double div_vel = st->grad->dvel[0][0] + st->grad->dvel[1][1] + st->grad->dvel[2][2];

  delta->velx += dt_half * (fac * (st->grad->drho[0] / st->rho + st->grad->dsgstPressure[0] / st->sgstPressure) *
                                (st->grad->dvel[0][0] - 1. / 3. * div_vel) +
                            fac / 2.0 * (st->grad->drho[1] / st->rho + st->grad->dsgstPressure[1] / st->sgstPressure) *
                                (st->grad->dvel[1][0] + st->grad->dvel[0][1]) +
                            fac / 2.0 * (st->grad->drho[2] / st->rho + st->grad->dsgstPressure[2] / st->sgstPressure) *
                                (st->grad->dvel[2][0] + st->grad->dvel[0][2]));

  delta->vely += dt_half * (fac / 2.0 * (st->grad->drho[0] / st->rho + st->grad->dsgstPressure[0] / st->sgstPressure) *
                                (st->grad->dvel[0][1] + st->grad->dvel[1][0]) +
                            fac * (st->grad->drho[1] / st->rho + st->grad->dsgstPressure[1] / st->sgstPressure) *
                                (st->grad->dvel[1][1] - 1. / 3. * div_vel) +
                            fac / 2.0 * (st->grad->drho[2] / st->rho + st->grad->dsgstPressure[2] / st->sgstPressure) *
                                (st->grad->dvel[2][1] + st->grad->dvel[1][2]));

  delta->velz += dt_half * (fac / 2.0 * (st->grad->drho[0] / st->rho + st->grad->dsgstPressure[0] / st->sgstPressure) *
                                (st->grad->dvel[0][2] + st->grad->dvel[2][0]) +
                            fac / 2.0 * (st->grad->drho[1] / st->rho + st->grad->dsgstPressure[1] / st->sgstPressure) *
                                (st->grad->dvel[1][2] + st->grad->dvel[2][1]) +
                            fac * (st->grad->drho[2] / st->rho + st->grad->dsgstPressure[2] / st->sgstPressure) *
                                (st->grad->dvel[2][2] - 1. / 3. * div_vel));
#endif /* #ifndef SGS_TURBULENCE_CONSTANT_EDDY_VISCOSITY */

#endif /* #if defined(SECOND_DERIVATIVES) && defined(RECONSTRUCT_GRADIENTS) */
}

#endif /* #ifdef SGS_TURBULENCE_STRESS_TENSOR */

/* functions to include turbulent production */
#ifdef SGS_TURBULENCE_TURBULENT_PRODUCTION

/*! \brief Compute turbulent production for the eddy viscosity closure
 *
 * 	This function computes the SGS turbulence energy that is produced
 * 	by the turbulent cascade for the eddy viscosity closure. The
 * 	corresponding energy is added to the SGS turbulence energy and
 * 	removed from the sum of kinetic and internal energy of the cell.
 *
 * 	\return void
 */
void do_sgs_energy_turbulent_production_eddy_viscosity(void)
{
#ifdef SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION
  double prod_sum        = 0;
  double prod_sum_global = 0;
#endif /* #ifdef SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION */

  int err_count = 0, err_global;

  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      /*compute half the time step*/
      double dt_cell = 0.5 * (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

      /*compute eddy viscosity*/
#ifdef SGS_TURBULENCE_CONSTANT_EDDY_VISCOSITY
      double eddy_visc = All.SgsTConst.ConstantEddyViscosity;
#else  /* #ifdef SGS_TURBULENCE_CONSTANT_EDDY_VISCOSITY */
      double filter_scale = calculate_filter_scale(SphP[i].Volume);
      double eddy_visc    = calculate_eddy_viscosity(filter_scale, SphP[i].Density, SphP[i].SgsTData.Pressure);
#endif /* #ifdef SGS_TURBULENCE_CONSTANT_EDDY_VISCOSITY #else */

      /*compute rate of strain tensor S_ij*/
      double S_xx = SphP[i].Grad.dvel[0][0];
      double S_yy = SphP[i].Grad.dvel[1][1];
      double S_zz = SphP[i].Grad.dvel[2][2];

      double S_xy = 0.5 * (SphP[i].Grad.dvel[0][1] + SphP[i].Grad.dvel[1][0]);
      double S_xz = 0.5 * (SphP[i].Grad.dvel[0][2] + SphP[i].Grad.dvel[2][0]);
      double S_yz = 0.5 * (SphP[i].Grad.dvel[1][2] + SphP[i].Grad.dvel[2][1]);

      double trace_S = S_xx + S_yy + S_zz;

      /*compute rate of subgrid-scale turbulence energy production*/
      double rate_of_turbulent_energy_production = dt_cell * SphP[i].Volume * 2 * eddy_visc *
                                                   (S_xx * S_xx + S_yy * S_yy + S_zz * S_zz + 2. * S_xy * S_xy + 2. * S_xz * S_xz +
                                                    2. * S_yz * S_yz - 1. / 3. * trace_S * trace_S);

      /*test whether there is enough resolved kinetic energy in cell*/
      double E_kin = 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);

      if(E_kin >= rate_of_turbulent_energy_production)
        {
          /*apply rate of turbulent energy production*/
          SphP[i].Energy -= rate_of_turbulent_energy_production;
          SphP[i].SgsTData.Energy += rate_of_turbulent_energy_production;

#ifdef SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION
          prod_sum += rate_of_turbulent_energy_production;
#endif /* #ifdef SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION */
#ifdef SGS_TURBULENCE_MAKE_PRODUCTION_DISSIPATION_SPECTRUM
          SphP[i].SgsTData.EgyTurbProd += rate_of_turbulent_energy_production;
#endif /* #ifdef SGS_TURBULENCE_MAKE_PRODUCTION_DISSIPATION_SPECTRUM */
        }
      else
        {
          /*remove all available resolved kinetic energy and set the momentum to zero*/

          err_count++;

          SphP[i].Energy -= E_kin;
          SphP[i].SgsTData.Energy += E_kin;

          SphP[i].Momentum[0] = 0;
          SphP[i].Momentum[1] = 0;
          SphP[i].Momentum[2] = 0;

#ifdef SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION
          prod_sum += E_kin;
#endif /* #ifdef SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION */
#ifdef SGS_TURBULENCE_MAKE_PRODUCTION_DISSIPATION_SPECTRUM
          SphP[i].SgsTData.EgyTurbProd += E_kin;
#endif /* #ifdef SGS_TURBULENCE_MAKE_PRODUCTION_DISSIPATION_SPECTRUM */
        }
    }

#ifdef SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION
  MPI_Allreduce(&prod_sum, &prod_sum_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  All.SgsTConst.SgsTurbulenceProduction += prod_sum_global;
#endif /* #ifdef SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION */

  MPI_Allreduce(&err_count, &err_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(err_global > 0)
    mpi_printf("SGS_TURBULENCE_TURBULENT_PRODUCTION: not enough kinetic energy in %d cells\n", err_global);
}

#endif /* #ifdef SGS_TURBULENCE_TURBULENT_PRODUCTION */

/* functions to include viscous dissipation */
#ifdef SGS_TURBULENCE_VISCOUS_DISSIPATION

/*! \brief Compute viscous dissipation for the eddy viscosity closure
 *
 * 	This function computes the SGS turbulence energy that is dissipated
 * 	if the eddy viscosity closure is used. The corresponding energy is
 * 	removed from the SGS turbulence energy and added to the sum of kinetic
 * 	and internal energy of the cell.
 *
 * 	\return void
 */
void do_sgs_turbulence_viscous_dissipation_eddy_viscosity(void)
{
#ifdef SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION
  double diss_sum        = 0;
  double diss_sum_global = 0;
#endif /* #ifdef SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION */

  int err_count = 0, err_global;

  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      /*compute half the time step*/
      double dt_cell = 0.5 * (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

      /*compute filter scale*/
      double filter_scale = calculate_filter_scale(SphP[i].Volume);

      /*compute dissipated energy*/
      double energy_diss = dt_cell * SphP[i].Volume * All.SgsTConst.EddyViscosityParameterCepsilon * SphP[i].Density *
                           pow(SphP[i].SgsTData.SpecificEnergy, 3. / 2.) / filter_scale;

      /*Test whether enough subgrid-scale energy is available and then apply changes*/
      if(SphP[i].SgsTData.Energy - P[i].Mass * All.SgsTConst.MinimumSgsTSpecificEnergy >= energy_diss)
        {
          SphP[i].Energy += energy_diss;
          SphP[i].SgsTData.Energy -= energy_diss;

#ifdef SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION
          diss_sum += energy_diss;
#endif /* #ifdef SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION */
#ifdef SGS_TURBULENCE_MAKE_PRODUCTION_DISSIPATION_SPECTRUM
          SphP[i].SgsTData.EgyViscDiss += energy_diss;
#endif /* #ifdef SGS_TURBULENCE_MAKE_PRODUCTION_DISSIPATION_SPECTRUM */
        }
      else
        {
          err_count++;

          SphP[i].Energy += SphP[i].SgsTData.Energy - P[i].Mass * All.SgsTConst.MinimumSgsTSpecificEnergy;
          SphP[i].SgsTData.Energy = P[i].Mass * All.SgsTConst.MinimumSgsTSpecificEnergy;

#ifdef SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION
          diss_sum += SphP[i].SgsTData.Energy - P[i].Mass * All.SgsTConst.MinimumSgsTSpecificEnergy;
#endif /* #ifdef SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION */
#ifdef SGS_TURBULENCE_MAKE_PRODUCTION_DISSIPATION_SPECTRUM
          SphP[i].SgsTData.EgyViscDiss += SphP[i].SgsTData.Energy - P[i].Mass * All.SgsTConst.MinimumSgsTSpecificEnergy;
#endif /* #ifdef SGS_TURBULENCE_MAKE_PRODUCTION_DISSIPATION_SPECTRUM */
        }
    }

#ifdef SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION
  MPI_Allreduce(&diss_sum, &diss_sum_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  All.SgsTConst.SgsTurbulenceDissipation += diss_sum_global;
#endif /* #ifdef SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION */

  MPI_Allreduce(&err_count, &err_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(err_global > 0)
    mpi_printf("SGS_TURBULENCE_VISCOUS_DISSIPATION: not enough energy available in %d cells\n", err_global);
}

#endif /* #ifdef SGS_TURBULENCE_VISCOUS_DISSIPATION */

#endif /*#ifdef SGS_TURBULENCE_EDDY_VISCOSITY*/
