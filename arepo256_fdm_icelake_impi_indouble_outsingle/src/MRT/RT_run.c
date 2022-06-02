/*! * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/MRT/RT_run.c
 * \date        06/2018
 * \author      Rahul Kannan
 * \brief       Run routine that controls the RT loop
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../domain.h"
#include "../proto.h"
#include "../voronoi.h"
#include "RT.h"
#include "RT_proto.h"

#ifdef MRT

#ifdef MRT_SUBCYCLE
static void mrt_run_subcycle(void);
#else
static void mrt_run_no_subcycle(void);
#endif

void mrt_run(void)
{
#ifndef MRT_SUBCYCLE
  mrt_run_no_subcycle();
#else
  mrt_run_subcycle();
#endif
}

#ifndef MRT_SUBCYCLE

static void mrt_run_no_subcycle(void)
{
#ifdef MRT_SOURCES
  if(All.MRT_On == 0)
    return;
#endif

  exchange_primitive_variables_RT();

  calculate_gradients_RT();

  exchange_primitive_variables_and_gradients_RT();

#ifdef MRT_SETUP_SPECIAL_BOUNDARIES
  mrt_setup(0);
#endif

  compute_interface_fluxes_RT(&Mesh);

  //  update_primitive_variables_RT() ;
}

#else

static void mrt_run_subcycle(void)
{
  TIMER_START(CPU_RT);

  for(int i = 0; i < All.RTNumSubCycles / 2; i++)
    {
      mpi_printf_rt(0, "RT: Sub Cycle %d\n", i + 1);

#ifdef MRT_SETUP_SPECIAL_BOUNDARIES
      mrt_setup(0);
#endif

      update_primitive_variables_RT();

#if defined(MRT_RADIATION_PRESSURE) && !defined(MRT_CHEM_SG)
      do_radiation_pressure_source_terms();
#endif

#ifdef MRT_EQUIL_CHEM_COOL
      cooling_only();
#else

#ifdef MRT_CHEM_SG
      TIMER_START(CPU_RT_CHEM);
      set_rates_sgchem();
      evolve_chemistry();
      TIMER_STOP(CPU_RT_CHEM);
#else

#ifdef MRT_COUPLED_THERMOCHEMISTRY
      mrt_thermochemistry();
#else

#if defined(MRT_CHEMISTRY_PS2011) || defined(MRT_CHEMISTRY_PS2009) || defined(MRT_UV_ONLY_DUST)
      mrt_update_chemistry();
#endif

#ifdef MRT_IR_ONLY_CHEMISTRY
      mrt_IR_chemistry();
#endif

#ifdef MRT_COOLING_HEATING
      cooling_only();
#endif

#endif
#endif
#endif

#ifdef MRT_SOURCES
      if(All.MRT_On == 0)
        continue;
#endif

      exchange_primitive_variables_RT();

      calculate_gradients_RT();

      exchange_primitive_variables_and_gradients_RT();

#ifdef MRT_SETUP_SPECIAL_BOUNDARIES
      mrt_setup(0);
#endif

      compute_interface_fluxes_RT(&Mesh);
    }

  update_primitive_variables_RT();

  TIMER_STOP(CPU_RT);
}

#endif

void mrt_run_sources(void)
{
#if defined(MRT_RADIATION_PRESSURE) && !defined(MRT_CHEM_SG)
  do_radiation_pressure_source_terms();
#endif

#ifdef MRT_EQUIL_CHEM_COOL
  cooling_only();
#else

#ifdef MRT_CHEM_SG
  TIMER_START(CPU_RT_CHEM);
  set_rates_sgchem();
  evolve_chemistry();
  TIMER_STOP(CPU_RT_CHEM);
#else

#ifdef MRT_COUPLED_THERMOCHEMISTRY
  mrt_thermochemistry();
#else

#if defined(MRT_CHEMISTRY_PS2011) || defined(MRT_CHEMISTRY_PS2009) || defined(MRT_UV_ONLY_DUST)
  mrt_update_chemistry();
#endif

#ifdef MRT_IR_ONLY_CHEMISTRY
  mrt_IR_chemistry();
#endif

#ifdef MRT_COOLING_HEATING
  cooling_only();
#endif

#endif
#endif
#endif

  update_primitive_variables_RT();
  update_primitive_variables();
}

#endif
