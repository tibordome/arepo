/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/driftfac.c
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

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

enum
{
  TABLE_LENGTH = 1000 /**< length of the lookup table used to hold the drift and kick factors */
};

/** table for the cosmological drift factors */
static double DriftTable[TABLE_LENGTH];

/** table for the cosmological kick factor for gravitational forces */
static double GravKickTable[TABLE_LENGTH];

/** table for the cosmological kick factor for hydrodynmical forces */
static double HydroKickTable[TABLE_LENGTH];

static double logTimeBegin;
static double logTimeMax;

/*! \brief Integrand for drift factor calculation.
 *
 *  For cosmological simulations.
 *
 *  \param[in] a Scale factor.
 *  \param[in] param (unused)
 *
 *  \return Integrand for drift factor calculation.
 */
double drift_integ(const double a, void *param)
{
  const double h = hubble_function(a);
  return 1 / (h * pow(a, 3));
}

/*! \brief Integrand for gravitational kick factor calculation.
 *
 *  For cosmological simulations.
 *
 *  \param[in] a Scale factor.
 *  \param[in] param (unused)
 *
 *  \return Integrand for gravitational kick factor calculation.
 */
double gravkick_integ(const double a, void *param)
{
  const double h = hubble_function(a);
  return 1 / (h * pow(a, 2));
}

/*! \brief Integrand for hydrodynamics kick factor calculation.
 *
 *  For cosmological simulations.
 *
 *  \param[in] a Scale factor.
 *  \param[in] param (unused)
 *
 *  \return Integrand for hydrodynamics kick factor calculation.
 */
double hydrokick_integ(const double a, void *param)
{
  const double h = hubble_function(a);
  return 1 / (h * pow(a, 3 * GAMMA_MINUS1) * a);
}

/*! \brief Initializes lookup table for cosmological pre-factors for a drift.
 *
 *  Numerical integrals using the integrand functions defined above.
 *
 *  \return void
 */
void init_drift_table(void)
{
  enum
  {
    WORKSIZE = 100000
  };

#ifdef DARKENERGY_DEBUG
  /*---------- DEBUG ! ------------*/
  FILE *FdDKfac;
  char buf[MAXLEN_PATH];

  file_path_sprintf(buf, "%s/driftkickfac.txt", All.OutputDir);
  FdDKfac = fopen(buf, "w");
  fprintf(FdDKfac, "i a drift GravKick HydroKick\n");
  /*---------- DEBUG ! ------------*/
#endif

  logTimeBegin = log(All.TimeBegin);
  logTimeMax   = log(All.TimeMax);

  gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(WORKSIZE);

  for(int i = 0; i < TABLE_LENGTH; i++)
    {
      double result, abserr;
      gsl_function F;
      F.function = &drift_integ;
      gsl_integration_qag(&F, exp(logTimeBegin), exp(logTimeBegin + ((logTimeMax - logTimeBegin) / TABLE_LENGTH) * (i + 1)), 0, 1.0e-8,
                          WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
      DriftTable[i] = result;

      F.function = &gravkick_integ;
      gsl_integration_qag(&F, exp(logTimeBegin), exp(logTimeBegin + ((logTimeMax - logTimeBegin) / TABLE_LENGTH) * (i + 1)), 0, 1.0e-8,
                          WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
      GravKickTable[i] = result;

      F.function = &hydrokick_integ;
      gsl_integration_qag(&F, exp(logTimeBegin), exp(logTimeBegin + ((logTimeMax - logTimeBegin) / TABLE_LENGTH) * (i + 1)), 0, 1.0e-8,
                          WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
      HydroKickTable[i] = result;

#ifdef DARKENERGY_DEBUG
      /*---------- DEBUG ! ------------*/
      fprintf(FdDKfac, "%d %e %e %e %e \n", i, exp(logTimeBegin + ((logTimeMax - logTimeBegin) / TABLE_LENGTH) * (i + 1)),
              DriftTable[i], GravKickTable[i], HydroKickTable[i]);
      /*---------- DEBUG ! ------------*/
#endif
    }

  gsl_integration_workspace_free(workspace);

#ifdef DARKENERGY_DEBUG
  /*---------- DEBUG ! ------------*/
  fclose(FdDKfac);
  /*---------- DEBUG ! ------------*/
#endif
}

/*! \brief This function linearly interpolates between elements of the lookup
 *         tables for the cosmological drift or kick prefactors.
 *
 *  \param[in] ti    Time.
 *  \param[in] Table Lookup table (DriftTable, GravKickTable, or HydroKickTable).
 */
static double interpolate_lookup_table(const integertime ti, const double *const Table)
{
  const double log_a = logTimeBegin + ti * All.Timebase_interval;
  const double u     = (log_a - logTimeBegin) / (logTimeMax - logTimeBegin) * TABLE_LENGTH;
  size_t i           = u;
  if(i >= TABLE_LENGTH)
    i = TABLE_LENGTH - 1;
  const double table_left  = i == 0 ? 0 : Table[i - 1];
  const double table_right = Table[i];
  const double result      = table_left + (table_right - table_left) * (u - i);
  return result;
}

/*! \brief This function integrates the cosmological drift or kick prefactor
 *         (determined by the given lookup table) for a time step between time0
 *         and time1.
 *
 *  \param[in] time0 Start time.
 *  \param[in] time1 End time.
 *  \param[in] Table Lookup table (DriftTable, GravKickTable, or HydroKickTable).
 */
static double get_cosmo_factor(const integertime time0, const integertime time1, const double *const Table)
{
  /* note: will only be called for cosmological integration */
  myassert(All.ComovingIntegrationOn);
  myassert(time0 >= 0 && time1 >= 0);
  myassert(time1 >= time0);
  const double cf0 = interpolate_lookup_table(time0, Table);
  const double cf1 = interpolate_lookup_table(time1, Table);
  myassert(cf1 >= cf0);
  return cf1 - cf0;
}

/*! \brief This function integrates the cosmological prefactor for a drift
 *         step between time0 and time1. A lookup table is used for reasons
 *         of speed.
 *
 *  \param[in] time0 Start time.
 *  \param[in] time1 End time.
 *
 *  \return \f[ \int_{a_0}^{a_1} \frac{{\rm d}a}{a^3 * H(a)} \f].
 */
double get_drift_factor(const integertime time0, const integertime time1)
{
  myassert(time0 >= 0 && time1 >= 0);
  static integertime last_time0 = -1, last_time1 = -1;
  static double last_value = NAN;

  if(time0 == last_time0 && time1 == last_time1)
    return last_value;

  last_time0        = time0;
  last_time1        = time1;
  return last_value = get_cosmo_factor(time0, time1, DriftTable);
}

/*! \brief This function integrates the cosmological prefactor for a
 *         gravitational kick between time0 and time1. A lookup table is used
 *         for reasons of speed.
 *
 *  \param[in] time0 Start time.
 *  \param[in] time1 End time.
 *
 *   \return Gravkick factor.
 */
double get_gravkick_factor(const integertime time0, const integertime time1)
{
  myassert(time0 >= 0 && time1 >= 0);
  static integertime last_time0 = -1, last_time1 = -1;
  static double last_value = NAN;

  if(time0 == last_time0 && time1 == last_time1)
    return last_value;

  last_time0        = time0;
  last_time1        = time1;
  return last_value = get_cosmo_factor(time0, time1, GravKickTable);
}

/*! \brief This function integrates the cosmological prefactor for a
 *         hydrodynamical kick between time0 and time1. A lookup table is
 *         used for reasons of speed.
 *
 *  \param[in] time0 Start time
 *  \param[in] time1 End time
 *
 *   \return Hydro kick factor.
 */
double get_hydrokick_factor(const integertime time0, const integertime time1)
{
  myassert(time0 >= 0 && time1 >= 0);
  static integertime last_time0 = -1, last_time1 = -1;
  static double last_value = NAN;

  if(time0 == last_time0 && time1 == last_time1)
    return last_value;

  last_time0        = time0;
  last_time1        = time1;
  return last_value = get_cosmo_factor(time0, time1, HydroKickTable);
}
