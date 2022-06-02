/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sgs_turbulence/sgs_turbulence.c
 * \date        07/2017
 * \author      S. Jacob
 * \brief       Basic functions for the SGS turbulence model
 * \details     This file contains the most general functions of the SGS turbulence model. This includes functions
 *  			for the initialisation, calculating effective filter scales, calling individual source functions and writing
 * log files.
 *
 * 				contains functions:
 * 					void init_sgs_turbulence(void)
 *
 * 					double calculate_filter_scale(double vol)
 * 					double calculate_filter_scale_at_face(const struct state *st_L, const struct state *st_R)
 *
 * 					void do_sgs_turbulence_source_terms_first_half(void)
 * 					void do_sgs_turbulence_source_terms_second_half(void)
 *
 * 					void log_sgs_turbulence(void)
 * 					void log_sgs_turbulence_production_dissipation(void) (#if
 * defined(SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION))
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"
#include "string.h"

#ifdef SGS_TURBULENCE

/*! \brief Initialisation of the module SGS_TURBULENCE
 *
 * This function sets the effective adiabatic index for the SGS turbulence energy.
 * If the eddy viscosity closure is used, it sets the two eddy viscosity parameters
 * C_nu and C_epsilon to the default values, if 0 is chosen in the initial conditions.
 * If a constant eddy viscosity is enforced, this constant value is also set here.
 * If the dissipated and produced SGS turbulence energy is written to a log file,
 * the corresponding parameters are initialized to zero.
 *
 * \return void
 */
void init_sgs_turbulence(void)
{
  /*For now: GammaSgsT = 5/3 always, but note that it actually depends on the dimensions*/
  All.SgsTConst.GammaSgsT = 5. / 3.;

#ifdef SGS_TURBULENCE_EDDY_VISCOSITY
  if(All.SgsTConst.EddyViscosityParameterCnu == 0)
    All.SgsTConst.EddyViscosityParameterCnu =
        0.05; /*Schmidt 2006 and References in Living Review: Large Eddy Simulations in Astrophysics*/
  if(All.SgsTConst.EddyViscosityParameterCepsilon == 0)
    All.SgsTConst.EddyViscosityParameterCepsilon = 1.58; /*Schmidt et al. 2014 and references therein*/
  mpi_printf("SGS_TURBULENCE: Cnu: %g, Cepsilon: %g \n", All.SgsTConst.EddyViscosityParameterCnu,
             All.SgsTConst.EddyViscosityParameterCepsilon);
#ifdef SGS_TURBULENCE_CONSTANT_EDDY_VISCOSITY
  All.SgsTConst.ConstantEddyViscosity = 0.002;
#endif /* #ifdef SGS_TURBULENCE_CONSTANT_EDDY_VISCOSITY */
#endif /* #ifdef SGS_TURBULENCE_EDDY_VISCOSITY */

#ifdef SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION
  All.SgsTConst.SgsTurbulenceProduction  = 0;
  All.SgsTConst.SgsTurbulenceDissipation = 0;
#endif /* #ifdef SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION */
}

/*! \brief Calculate the filter scale
 *
 * This function determines the effective filter scale for each cell from its
 * volume in one, two or three dimensions. The cell can be treated as a cube
 * or a sphere.
 *
 *  \param[in] vol Volume of a cell
 *
 *  \return Volume of the cell, in code units
 */
double calculate_filter_scale(double vol)
{
#ifdef SGS_TURBULENCE_VIEW_CELLS_AS_CUBES
/* distinguish 1D, 2D and 3D simulations */
#ifdef ONEDIMS
  return vol;
#elif defined(TWODIMS) /* #ifdef ONEDIMS */
  return pow(vol, 1. / 2.);
#else                  /* #ifdef ONEDIMS #elif */
  return pow(vol, 1. / 3.);
#endif                 /* #ifdef ONEDIMS #elif #else */

#elif defined(SGS_TURBULENCE_VIEW_CELLS_AS_SPHERES) /* SGS_TURBULENCE_VIEW_CELLS_AS_CUBES */
/* distinguish 1D, 2D and 3D simulations */
#ifdef ONEDIMS
  return vol;
#elif defined(TWODIMS) /* #ifdef ONEDIMS */
  return pow(1. / (4. * M_PI) * vol, 1. / 2.);
#else                  /* #ifdef ONEDIMS #elif */
  return pow(3. / (4. * M_PI) * vol, 1. / 3.);
#endif                 /* #ifdef ONEDIMS #elif #else */

#else /* SGS_TURBULENCE_VIEW_CELLS_AS_CUBES #elif */
/*if you end up here, you forgot a flag*/
#error "Choose a way to compute filter scale from volume"
#endif /* SGS_TURBULENCE_VIEW_CELLS_AS_CUBES #elif #else */
}

/*! \brief Compute filter scale at an interface
 *
 * This function computes the filter scale at an interface between two cells
 * with potentially different volumes and thus corresponding filter scales.
 * This function is still preliminary since it assumes that the volumes of
 * all cells are equal.
 *
 * \param[in] st_L Left state of the interface
 * \param[in] st_R Right state of the interface
 *
 * \return Average filter scale on the interface, in code units
 */
double calculate_filter_scale_at_face(const struct state *st_L, const struct state *st_R)
{
  /*for now just assume static, cartesian mesh, so all volumes should be the same*/
  /*print message if this is not true to 1.e-7*/
  if((st_L->volume - st_R->volume) > 1.e-7)
    printf("CALCULATE FILTER SCALE AT FACE: volumes are not the same!");

  /*compute for each state the filter scale and then return average*/
  double filter_scale_L = calculate_filter_scale(st_L->volume);
  double filter_scale_R = calculate_filter_scale(st_R->volume);

  return 0.5 * (filter_scale_L + filter_scale_R);
}

/*! \brief Compute first half of the source terms for SGS turbulence
 *
 * This function calls the individual source functions for the SGS turbulence
 * model for the first half time step of the source term calculation.
 * The possible source terms are turbulent production and viscous dissipation.
 *
 * \return void
 */
void do_sgs_turbulence_source_terms_first_half(void)
{
  mpi_printf("SGS_TURBULENCE: Computing source terms first half\n");

#ifdef SGS_TURBULENCE_TURBULENT_PRODUCTION
  do_sgs_energy_turbulent_production();
#endif /* #ifdef SGS_TURBULENCE_TURBULENT_PRODUCTION */

#ifdef SGS_TURBULENCE_VISCOUS_DISSIPATION
  do_sgs_turbulence_viscous_dissipation();
#endif /* #ifdef SGS_TURBULENCE_VISCOUS_DISSIPATION */
}

/*! \brief Compute second half of the source terms for SGS turbulence
 *
 * This function calls the individual source functions for the SGS turbulence
 * model for the second half time step of the source term calculation. The
 * individual source terms are called in reverse order compared to the first
 * source term computation.
 *
 * \return void
 */
void do_sgs_turbulence_source_terms_second_half(void)
{
  mpi_printf("SGS_TURBULENCE: Computing source terms second half\n");

#ifdef SGS_TURBULENCE_VISCOUS_DISSIPATION
  do_sgs_turbulence_viscous_dissipation();
#endif /* #ifdef SGS_TURBULENCE_VISCOUS_DISSIPATION */

#ifdef SGS_TURBULENCE_TURBULENT_PRODUCTION
  do_sgs_energy_turbulent_production();
#endif /* #ifdef SGS_TURBULENCE_TURBULENT_PRODUCTION */
}

/*! \brief Write a log-file for the SGS turbulence energy
 *
 * This function writes the log-file for SGS turbulence energy.
 * The file contains
 * - the current time,
 * - and the total specific SGS turbulence energy.
 *
 * \return void
 */
void log_sgs_turbulence(void)
{
  /*compute total Ksgs*/
  double Ksgs_sum = 0, Ksgs_global;
  double mass_sum = 0, mass_global;

  for(int i = 0; i < NumGas; i++)
    {
      if(P[i].Mass == 0 && P[i].ID == 0)
        continue;

      Ksgs_sum += P[i].Mass * SphP[i].SgsTData.SpecificEnergy;

      mass_sum += P[i].Mass;
    }

  MPI_Allreduce(&mass_sum, &mass_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&Ksgs_sum, &Ksgs_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fprintf(FdSgsTurbulence, "%g %g\n", All.Time, Ksgs_global / mass_global);
      myflush(FdSgsTurbulence);
    }
}

#if defined(SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION)
/*! \brief Write a log-file for the SGS turbulence energy that is produced and dissipated
 *
 * This function writes a log-file for the different production and dissipation mechanisms
 * of SGS turbulence energy.
 * The file contains
 * - the current time,
 * - the SGS turbulence energy produced by turbulent production,
 * - the SGS turbulence energy lost by viscous dissipation,
 * - and the change of SGS turbulence energy through adiabatic contraction/expansion.
 *
 * \return void
 */
void log_sgs_turbulence_production_dissipation(void)
{
  double dEnergyAdiab = 0, dEnergyAdiab_global;

  for(int i = 0; i < NumGas; i++)
    {
      if(P[i].Mass == 0 && P[i].ID == 0)
        continue;

      dEnergyAdiab += SphP[i].SgsTData.dEnergy_Adiab;
      SphP[i].SgsTData.dEnergy_Adiab = 0;
    }

  MPI_Allreduce(&dEnergyAdiab, &dEnergyAdiab_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fprintf(FdSgsTurbulenceProductionDissipation, "%g %g %g %g \n", All.Time, All.SgsTConst.SgsTurbulenceProduction,
              All.SgsTConst.SgsTurbulenceDissipation, dEnergyAdiab_global);
      myflush(FdSgsTurbulenceProductionDissipation);
    }
}
#endif /* #if defined(SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION) */

#endif /* #ifdef SGS_TURBULENCE */
