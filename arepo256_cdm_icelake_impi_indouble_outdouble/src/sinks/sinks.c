/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sinks/sinks.c
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       Sink particle driver
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"

#ifdef SINKS

void sinks(void)
{
  double t0, t1, dt, dt_max;

  if(All.Time == All.TimeBegin || SKD.AccRad == 0)
    return;

  TIMER_START(CPU_SINKS);

  if(SKD.NHThresh)
    {
      t0 = second();

      TIMER_START(CPU_SINKS_ACCRETION);

      sinks_set_constants();

      sinks_get_num_sinks();
      mpi_printf("SINKS: Got sink numbers %d\n", SKD.TotNumSinks);

      sinks_begin();

      sinks_get_active_sinks();
      mpi_printf("SINKS: Got active sink numbers %d\n", TotNSinks);

      sinks_dmass();

      sinks_accrete();  // Debugging, 2016

      sinks_end();

      TIMER_STOP(CPU_SINKS_ACCRETION);

      sinks_create();

      sinks_check_AuxDataID_references();

      t1 = second();

      dt = timediff(t0, t1);

      MPI_Allreduce(&dt, &dt_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      mpi_printf("SINKS: Done! Took %g seconds\n", dt_max);
    }

  TIMER_STOP(CPU_SINKS);

  //  endrun();
}

#ifdef SINKS_MERGERS
void sinks_mergers(void)
{
  TIMER_START(CPU_SINKS);

  set_cosmo_factors_for_current_time();
  sinks_find_neighboring_sinks();
  sinks_do_mergers();

  TIMER_STOP(CPU_SINKS);
}
#endif

#endif
