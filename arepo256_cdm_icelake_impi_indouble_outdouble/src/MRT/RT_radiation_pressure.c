/*! * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/MRT/RT.c
 * \date        MM/YYYY
 * \author      Rahul Kannan, Federico Marincaci, Mark Vogelsberger
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

/*! \file RT.c
 *  \brief  main driver for a moment based RT with the M1 closure.
 *
 *
 */

#include <math.h>
#include <mpi.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../domain.h"
#include "../proto.h"
#include "../voronoi.h"

#include "RT_proto.h"

/*Add the radiation pressure source terms*/
#if defined(MRT_RADIATION_PRESSURE) && !defined(MRT_CHEM_SG)
void do_radiation_pressure_source_terms(void)
{
  TIMER_START(CPU_RT_RP);

  mpi_printf_rt(1, "RT: Adding the radiation pressure source term...\n");
  /*Single frequency d(rho v)/dt = nHI * sigma * F * E ({energy of absorbed photons}) / CLIGHT*/
  int idx, i, num1;
  double a3inv, hubble_a;
  double c_speed = 2.99792458e10 / All.UnitVelocity_in_cm_per_s;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) *
                  All.Timebase_interval; /*factor of 0.5 included because of spilt for RK time integration*/

#ifdef MRT_SUBCYCLE
      dt /= ((double)(All.RTNumSubCycles));
#else
      dt *= 0.5;
#endif

      if(All.ComovingIntegrationOn)
        {
          a3inv    = 1.0 / All.Time / All.Time / All.Time;
          hubble_a = hubble_function(All.Time);
        }
      else
        {
          a3inv = hubble_a = 1.0;
        }

#ifndef MRT_UPDATE_AT_END_OF_STEP
      double KE_old = (SphP[i].Momentum[0] * SphP[i].Momentum[0] + SphP[i].Momentum[1] * SphP[i].Momentum[1] +
                       SphP[i].Momentum[2] * SphP[i].Momentum[2]) /
                      P[i].Mass / 2.0;
#endif

      int j;

#if !defined(MRT_NO_UV) && !defined(MRT_UV_ONLY_DUST)

      //       nH = (HYDROGEN_MASSFRAC * SphP[i].Density * All.cf_a3inv) / (PROTONMASS / All.UnitMass_in_g * All.HubbleParam);

      double molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC + 4 * HYDROGEN_MASSFRAC * SphP[i].Ne);

#ifdef MRT_LEVITATION_TEST
      molecular_weight = 2.33;
#endif

      double nH  = HYDROGEN_MASSFRAC * SphP[i].Density * a3inv / PROTONMASS * All.UnitMass_in_g / All.HubbleParam;
      double nHI = SphP[i].HI * nH;
      // double nHI = nH ;

#ifdef MRT_INCLUDE_HE
      double nHeI   = SphP[i].HeI * nH;
      double nHeII  = SphP[i].HeII * nH;
      double nHeIII = SphP[i].HeIII * nH;
#endif

      for(j = 0; j < UV_BINS; j++)
        {
          double mom_inj = 0.0;
          double n_gamma = SphP[i].DensPhot[j] * 1e63 * SphP[i].Volume;

          double modF = sqrt(SphP[i].RT_F[j][0] * SphP[i].RT_F[j][0] + SphP[i].RT_F[j][1] * SphP[i].RT_F[j][1] +
                             SphP[i].RT_F[j][2] * SphP[i].RT_F[j][2]);
          if(modF == 0.0)
            modF = 1.0;

          mom_inj += nHI * mrt_sigma_HI[j] * P_HI[j];
#ifdef MRT_INCLUDE_HE
          mom_inj += nHeI * mrt_sigma_HeI[j] * P_HeI[j];
          mom_inj += nHeII * mrt_sigma_HeII[j] * P_HeII[j];
#endif

#ifdef MRT_IR
          mom_inj += SphP[i].KappaIR_R[j] * MeanPhotonEnergy[j];
#endif
          if(isnan(mom_inj))
            terminate("rad pressure gone wrong!!! nHI = %g HI = %g \n\n", nHI, SphP[i].HI);

          for(num1 = 0; num1 < 3; num1++)
            {
#ifdef MRT_UPDATE_AT_END_OF_STEP
              SphP[i].RT_mominj[num1] += dt * mom_inj * SphP[i].RT_F[j][num1] * 1e63 * SphP[i].Volume / c_speed;
#else
              SphP[i].Momentum[num1] += dt * mom_inj * SphP[i].RT_F[j][num1] * 1e63 * SphP[i].Volume / c_speed;
#endif
            }
        }
#endif

        /* Only add momentum to gas cells from dust radiation pressure if live
         * dust radiation pressure tracking is not active. */
#if defined(MRT_UV_ONLY_DUST) && !defined(DL_RADIATION_PRESSURE)
      for(j = 0; j < UV_BINS; j++)
        {
          for(num1 = 0; num1 < 3; num1++)
            {
#ifdef MRT_UPDATE_AT_END_OF_STEP
              SphP[i].RT_mominj[num1] += dt * SphP[i].KappaIR_R[j] * SphP[i].RT_F[j][num1] *
                                         (1e63 * MeanPhotonEnergy[j] * ELECTRONVOLT_IN_ERGS / All.UnitEnergy_in_cgs) * SphP[i].Volume /
                                         c_speed;
#else
              SphP[i].Momentum[num1] += dt * SphP[i].KappaIR_R[j] * SphP[i].RT_F[j][num1] *
                                        (1e63 * MeanPhotonEnergy[j] * ELECTRONVOLT_IN_ERGS / All.UnitEnergy_in_cgs) * SphP[i].Volume /
                                        c_speed;
#endif
            }
        }
#endif

#if defined(MRT_IR) && !defined(DL_RADIATION_PRESSURE)
#ifdef BOUNDARY_REFL_SOLIDSIDE_MINID
      if(P[i].ID < BOUNDARY_REFL_SOLIDSIDE_MINID || P[i].ID > BOUNDARY_REFL_SOLIDSIDE_MAXID)
        {
#endif
          for(j = UV_BINS; j < MRT_BINS; j++)
            {
              for(num1 = 0; num1 < 3; num1++)
                {
#ifdef MRT_UPDATE_AT_END_OF_STEP
                  SphP[i].RT_mominj[num1] += dt * SphP[i].KappaIR_R[j] * SphP[i].RT_F[j][num1] * SphP[i].Volume / c_speed;
#else
              SphP[i].Momentum[num1] += dt * SphP[i].KappaIR_R[j] * SphP[i].RT_F[j][num1] * SphP[i].Volume / c_speed;
#endif
                }
            }
#ifdef BOUNDARY_REFL_SOLIDSIDE_MINID
        }
#endif
#endif

#ifndef MRT_UPDATE_AT_END_OF_STEP
      double KE_new = (SphP[i].Momentum[0] * SphP[i].Momentum[0] + SphP[i].Momentum[1] * SphP[i].Momentum[1] +
                       SphP[i].Momentum[2] * SphP[i].Momentum[2]) /
                      P[i].Mass / 2.0;

      SphP[i].Energy += KE_new - KE_old;
#endif
    }
  TIMER_STOP(CPU_RT_RP);
}
#endif
