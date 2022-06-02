/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/longrange.c
 * \date        MM/YYYY
 * \author
 * \brief       Driver routines for computation of long-range gravitational PM force
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

#ifdef PMGRID

/*! \brief Driver routine to call initialization of periodic or/and
 *         non-periodic FFT routines.
 *
 *  \return void
 */
void long_range_init(void)
{
#ifdef GRAVITY_NOT_PERIODIC
  pm_init_nonperiodic();
#else
  pm_init_periodic();
#ifdef TWODIMS
  pm2d_init_periodic();
#endif
#ifdef PLACEHIGHRESREGION
  pm_init_nonperiodic();
#endif
#endif
}

/*! \brief Driver routine to determine the extend of the non-
 *         periodic or high resolution region.
 *
 *  The initialization is done by pm_init_regionsize(). Afterwards,
 *  the convolution kernels are computed by pm_setup_nonperiodic_kernel().
 *
 *  \return void
 */
void long_range_init_regionsize(void)
{
#if defined(GRAVITY_NOT_PERIODIC) || defined(PLACEHIGHRESREGION)
  if(RestartFlag != RESTART_RESTART)
    pm_init_regionsize();
  pm_setup_nonperiodic_kernel();
#endif
}

/*! \brief This function computes the long-range PM force for all particles.
 *
 *  In case of a periodic grid the force is calculated by pmforce_periodic()
 *  otherwise by pmforce_nonperiodic(). If a high resolution region is
 *  specified for the PM force, pmforce_nonperiodic() calculates that force in
 *  both cases.
 *
 *  \return void
 */
void long_range_force(void)
{
  TIMER_START(CPU_PM_GRAVITY);

  for(int i = 0; i < NumPart; i++)
    {
      for(int j = 0; j < 3; j++)
        P[i].GravPM[j] = 0;
#ifdef EVALPOTENTIAL
      P[i].PM_Potential = 0;
#endif
    }

#ifndef SELFGRAVITY
  return;
#endif

#ifdef GRAVITY_NOT_PERIODIC /* non-periodic PM mesh */
  int i = pmforce_nonperiodic(0);
  if(i == 1) /* this is returned if a particle was outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      i = pmforce_nonperiodic(0); /* try again */
    }
  if(i == 1)
    terminate("although we tried to increase the region, somehow we still don't fit all particles in it");
#ifdef PLACEHIGHRESREGION
  i = pmforce_nonperiodic(1);

  if(i == 1) /* this is returned if a particle was outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      /* try again */
      for(int i = 0; i < NumPart; i++)
        for(int j = 0; j < 3; j++)
          P[i].GravPM[j] = 0;
      i = pmforce_nonperiodic(0) + pmforce_nonperiodic(1);
    }
  if(i != 0)
    terminate("although we tried to increase the region, somehow we still don't fit all particles in it");
#endif

  if(All.ComovingIntegrationOn)
    {
      const double fac = 0.5 * pow(All.Hubble, 2) * All.Omega0;
      for(int i = 0; i < NumPart; i++)
        for(int j = 0; j < 3; j++)
          P[i].GravPM[j] += fac * P[i].Pos[j];
    }
  /* Finally, the following factor allows a computation of cosmological
   * simulations with vacuum energy in physical coordinates */
  else
    {
      const double fac = All.OmegaLambda * pow(All.Hubble, 2);
      for(int i = 0; i < NumPart; i++)
        for(int j = 0; j < 3; j++)
          P[i].GravPM[j] += fac * P[i].Pos[j];
    }
#else /* periodic PM mesh */

#ifdef TWODIMS
  pm2d_force_periodic(0);
#else
  pmforce_periodic(0, NULL);
#endif

#ifdef PLACEHIGHRESREGION
  int i = pmforce_nonperiodic(1);

  if(i == 1) /* this is returned if a particle was outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      i = pmforce_nonperiodic(1); /* try again */
    }
  if(i == 1)
    terminate("although we tried to increase the region, we still don't fit all particles in it");
#endif
#endif

  TIMER_STOP(CPU_PM_GRAVITY);

#ifndef LEGACY_DISPLACEMENT_CONSTRAINT
  find_long_range_step_constraint();
#endif
}

#endif
