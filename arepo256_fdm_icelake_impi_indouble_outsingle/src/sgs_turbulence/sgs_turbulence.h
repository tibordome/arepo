/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sgs_turbulence/sgs_turbulence.h
 * \date        07/2017
 * \author      S. Jacob
 * \brief       Functions of the SGS turbulence model
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#ifndef SGS_TURBULENCE_H
#define SGS_TURBULENCE_H

struct sgs_turbulence_constants
{
  double MinimumSgsTSpecificEnergy;
  double GammaSgsT; /* default is 5./3., depends on dimension */

#ifdef SGS_TURBULENCE_EDDY_VISCOSITY
  double EddyViscosityParameterCnu;
  double EddyViscosityParameterCepsilon;
#ifdef SGS_TURBULENCE_CONSTANT_EDDY_VISCOSITY
  double ConstantEddyViscosity;
#endif /*#ifdef SGS_TURBULENCE_CONSTANT_EDDY_VISCOSITY*/
#endif /*#ifdef SGS_TURBULENCE_EDDY_VISCOSITY*/

#ifdef SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION
  double SgsTurbulenceProduction;
  double SgsTurbulenceDissipation;
#endif /* #ifdef SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION */
};

struct sgs_turbulence_data
{
  MyFloat Energy;         /*energy in subgrid scale turbulence in a cell*/
  MyFloat SpecificEnergy; /*energy in subgrid scale turbulence per unit gas mass in a cell*/
  MyFloat Pressure;       /*Pressure from subgrid scale turbulence*/

#ifdef SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION
  MyFloat dEnergy_Adiab;
#endif /* #ifdef SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION */

#ifdef SGS_TURBULENCE_MAKE_PRODUCTION_DISSIPATION_SPECTRUM
  MyFloat EgyAdiab;
  MyFloat EgyTurbProd;
  MyFloat EgyViscDiss;

  MyFloat DuDt_Adiab;
  MyFloat DuDt_TurbProd;
  MyFloat DuDt_ViscDiss;
#endif /* #ifdef SGS_TURBULENCE_MAKE_PRODUCTION_DISSIPATION_SPECTRUM */
};

/* general functions (in sgs_turbulence.c)*/
void init_sgs_turbulence(void);

double calculate_filter_scale(double vol);
double calculate_filter_scale_at_face(const struct state *st_L, const struct state *st_R);

void do_sgs_turbulence_source_terms_first_half(void);
void do_sgs_turbulence_source_terms_second_half(void);

void log_sgs_turbulence(void);
#if defined(SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION)
void log_sgs_turbulence_production_dissipation(void);
#endif /* #if defined(SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION) */

/*functions to compute the flux due to the stress tensor (in sgs_turbulence_stress_tensor.c)*/
#ifdef SGS_TURBULENCE_STRESS_TENSOR

void face_get_reconstructed_gradients(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux);
void face_add_sgs_stress_tensor_fluxes(const struct state_face *st, const struct geometry *geom, struct fluxes *flux);
void face_extrapolate_sgs_turbulence_viscous_kick(struct state *delta, const struct state *st);

#endif /* #ifdef SGS_TURBULENCE_STRESS_TENSOR */

/* functions to compute turbulent production (in sgs_turbulence_turbulent_production.c) */
#ifdef SGS_TURBULENCE_TURBULENT_PRODUCTION

void do_sgs_energy_turbulent_production(void);

#endif /* #ifdef SGS_TURBULENCE_TURBULENT_PRODUCTION */

/* functions to compute viscous dissipation (in sgs_turbulence_viscous_dissipation.c) */
#ifdef SGS_TURBULENCE_VISCOUS_DISSIPATION

void do_sgs_turbulence_viscous_dissipation(void);

#endif /*#ifdef SGS_TURBULENCE_VISCOUS_DISSIPATION*/

/* functions for the eddy viscosity closure (in sgs_turbulence_eddy_viscosity_closure.c) */
#ifdef SGS_TURBULENCE_EDDY_VISCOSITY

double calculate_eddy_viscosity(double filter_scale, double rho, double sgstPressure);

#ifdef SGS_TURBULENCE_STRESS_TENSOR
void face_get_eddy_viscosity(const struct state *st_L, const struct state *st_R, struct state_face *st_face,
                             const struct fluxes *flux);
void face_add_sgs_stress_tensor_eddy_viscosity_fluxes(const struct state_face *st, const struct geometry *geom, struct fluxes *flux);
void face_extrapolate_sgs_turbulence_viscous_kick_eddy_viscosity(struct state *delta, const struct state *st);
#endif /* #ifdef SGS_TURBULENCE_STRESS_TENSOR */

#ifdef SGS_TURBULENCE_TURBULENT_PRODUCTION
void do_sgs_energy_turbulent_production_eddy_viscosity(void);
#endif /* #ifdef SGS_TURBULENCE_TURBULENT_PRODUCTION */

#ifdef SGS_TURBULENCE_VISCOUS_DISSIPATION
void do_sgs_turbulence_viscous_dissipation_eddy_viscosity(void);
#endif /* #ifdef SGS_TURBULENCE_VISCOUS_DISSIPATION */

#endif /* #ifdef SGS_TURBULENCE_EDDY_VISCOSITY */

#endif /* #ifndef SGS_TURBULENCE_H */
