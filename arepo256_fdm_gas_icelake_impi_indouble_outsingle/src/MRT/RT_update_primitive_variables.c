/*! * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/MRT/RT.c
 * \date        06/2018
 * \author      Rahul Kannan
 * \brief       updates the RT primitive variables.
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

/*! \file RT_semi_HYPRE.c
 *  \brief main driver for an moment based RT with the VET formalism
 *
 *
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
//#include "../proto.h"
#include "../domain.h"
#include "../voronoi.h"
#include "RT.h"
#include "RT_proto.h"

#ifdef MRT

void update_primitive_variables_RT(void)
{
  mpi_printf_rt(0, "RT: Updating primitive variables...\n");

  int idx, i;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

#ifdef MRT_CHEM_SG
      int index = i;

      double carb_abund, oxy_abund;

#ifdef SGCHEM_VARIABLE_Z
      carb_abund = SphP[index].CarbAbund;
      oxy_abund  = SphP[index].OxyAbund;
#else
      carb_abund = All.CarbAbund;
      oxy_abund  = All.OxyAbund;
#endif

      double non_eq_abundances_kk;
      for(int kk = 0; kk < SGCHEM_NUM_SPECIES; kk++)
        {
          non_eq_abundances_kk = SphP[index].TracAbund[kk];

          if(non_eq_abundances_kk < 0.0)
            {
#ifdef DEBUG_SGCHEM
              printf("update_primitive_variables.c negative abundance from advection, species = %d abundance = %g\n", i,
                     non_eq_abundances_kk);
              printf("update_primitive_variables.c: Setting abundance to +1e-20 (this might not help!)\n");
#endif
              non_eq_abundances_kk = 1e-20;
            }

          if(kk == IH2 && non_eq_abundances_kk > 0.5)
            {
#ifdef DEBUG_SGCHEM
              printf("update_primitive_variables.c H2 abundance greater than 0.5; abundance = %g\n", non_eq_abundances_kk);
              non_eq_abundances_kk = 0.5;
#endif
            }

          if(kk == IHP && non_eq_abundances_kk > 1.0)
            {
#ifdef DEBUG_SGCHEM
              printf("update_primitive_variables.c HP abundance greater than 1.0; abundance = %g\n", non_eq_abundances_kk);
#endif
              non_eq_abundances_kk = 1.0;
            }

          if(kk == ICO && non_eq_abundances_kk > fmin(carb_abund, oxy_abund))
            {
#ifdef DEBUG_SGCHEM
              printf("update_primitive_variables.c CO abundance greater than Carbon & Oxygen abundances; abundance = %g\n",
                     non_eq_abundances_kk);
#endif
              non_eq_abundances_kk = fmin(carb_abund, oxy_abund);
            }

          SphP[index].TracAbund[kk]     = non_eq_abundances_kk;
          SphP[index].MassTracAbund[kk] = P[index].Mass * SphP[index].TracAbund[kk];
        }

#endif

#if !defined(MRT_NO_UV) && !defined(MRT_CHEM_SG)
#ifndef MRT_EQUIL_CHEM_COOL

      double y_fac = (1.0 - HYDROGEN_MASSFRAC) / 4.0 / HYDROGEN_MASSFRAC;

      double nH_times_volume = P[i].Mass;

      SphP[i].HII   = SphP[i].nHII / nH_times_volume;
      SphP[i].HeII  = SphP[i].nHeII / nH_times_volume;
      SphP[i].HeIII = SphP[i].nHeIII / nH_times_volume;

      if(SphP[i].HII < 0.0)
        SphP[i].HII = 1e-8;
#ifdef MRT_INCLUDE_HE
      if(SphP[i].HeII < 0.0)
        SphP[i].HeII = 1e-8;

      if(SphP[i].HeIII < 0.0)
        SphP[i].HeIII = 1e-8;
#endif
      if(SphP[i].HII > 1.0)
        SphP[i].HII = 1.0;

#ifdef MRT_INCLUDE_HE
      if(SphP[i].HeII > y_fac)
        {
          SphP[i].HeII  = y_fac - 1e-8;
          SphP[i].HeIII = 1e-8;
        }

      if(SphP[i].HeIII > y_fac)
        {
          SphP[i].HeIII = y_fac - 1e-8;
          SphP[i].HeII  = 1e-8;
        }

      if(SphP[i].HeII + SphP[i].HeIII > y_fac)
        {
          SphP[i].HeIII *= y_fac / (SphP[i].HeII + SphP[i].HeIII);
          SphP[i].HeII *= y_fac / (SphP[i].HeII + SphP[i].HeIII);
        }
#else
      SphP[i].HeII  = 0.0;
      SphP[i].HeIII = 0.0;
      y_fac         = 0.0;
#endif

      SphP[i].nHII   = SphP[i].HII * nH_times_volume;
      SphP[i].nHeII  = SphP[i].HeII * nH_times_volume;
      SphP[i].nHeIII = SphP[i].HeIII * nH_times_volume;

      SphP[i].HI  = 1.0 - SphP[i].HII;
      SphP[i].HeI = y_fac - SphP[i].HeII - SphP[i].HeIII;
      SphP[i].Ne  = SphP[i].HII + SphP[i].HeII + 2.0 * SphP[i].HeIII;

      SphP[i].nHI  = SphP[i].HI * nH_times_volume;
      SphP[i].nHeI = SphP[i].HeI * nH_times_volume;
      SphP[i].ne   = SphP[i].Ne * nH_times_volume;

      if(SphP[i].HII < 0.0 || SphP[i].HeII < 0.0 || SphP[i].HeIII < 0.0 || SphP[i].HII > 1.0 || SphP[i].HeII > y_fac ||
         SphP[i].HeIII > y_fac)
        terminate("%g %g %g %g %d\n", SphP[i].HII, SphP[i].HeII, SphP[i].HeIII, y_fac, P[i].ID);
#endif
#endif

      for(int num1 = 0; num1 < MRT_BINS; num1++)
        {
          SphP[i].OldCons_DensPhot[num1] = SphP[i].Cons_DensPhot[num1];

          SphP[i].DensPhot[num1] = SphP[i].Cons_DensPhot[num1] / SphP[i].Volume;

          if(SphP[i].DensPhot[num1] < 0.0)
            terminate("Photon Energy Density < 0\n");

          if(SphP[i].DensPhot[num1] < MINDENSPHOT)
            {
              SphP[i].DensPhot[num1]      = MINDENSPHOT;
              SphP[i].Cons_DensPhot[num1] = SphP[i].DensPhot[num1] * SphP[i].Volume;
            }

          for(int j = 0; j < 3; j++)
            {
              SphP[i].RT_F[num1][j] = SphP[i].Cons_RT_F[num1][j] / SphP[i].Volume;
              SphP[i].FN[num1][j]   = SphP[i].RT_F[num1][j] / SphP[i].DensPhot[num1];
            }
        }

#if defined(MRT_IR) || defined(MRT_UV_ONLY_DUST) || defined(MRT_CHEM_SG)
      set_kappa_times_rho_IR(i, SphP);
#endif

      set_VET_single(i, SphP);

      for(int num1 = 0; num1 < MRT_BINS; num1++)
        {
          double modF         = sqrt(SphP[i].RT_F[num1][0] * SphP[i].RT_F[num1][0] + SphP[i].RT_F[num1][1] * SphP[i].RT_F[num1][1] +
                                     SphP[i].RT_F[num1][2] * SphP[i].RT_F[num1][2]);
          SphP[i].modFN[num1] = modF / SphP[i].DensPhot[num1];
          if(SphP[i].modFN[num1] > c_internal_units * 1.000000001)
            terminate("%g\t%g\t%g\t%g\t%g\t%g\t ID = %d\n", SphP[i].modFN[num1], c_internal_units, SphP[i].RT_F[num1][0],
                      SphP[i].RT_F[num1][1], SphP[i].RT_F[num1][2], SphP[i].DensPhot[num1], P[i].ID);
          // SphP[i].modFN[num1] = 0.9999999*c_internal_units ;

          if(SphP[i].DensPhot[num1] < 0.0)
            terminate("UPDATE PRIMITIVE VARIABLES Ngamma < 0.0\n");
          // if(isnan(SphP[i].DensPhot[num1]))
          // terminate("%g %g %g \n", SphP[i].Cons_DensPhot[num1], SphP[i].DensPhot[num1], SphP[i].Volume) ;

          if(isnan(SphP[i].RT_F[0][0]))
            terminate("Nan F\n");
        }
    }
}

#ifndef MRT_NO_UV
void initialize_ionic_species(void)
{
  for(int i = 0; i < NumGas; i++)
    {
      double nH_times_volume = P[i].Mass;
      double y_fac           = (1.0 - HYDROGEN_MASSFRAC) / 4.0 / HYDROGEN_MASSFRAC;
      SphP[i].HI             = 0.999999;
      SphP[i].HII            = 0.000001;
      SphP[i].nHI            = SphP[i].HI * nH_times_volume;
      SphP[i].nHII           = SphP[i].HII * nH_times_volume;

      SphP[i].HeI    = 0.999998 * y_fac;
      SphP[i].HeII   = 0.000001 * y_fac;
      SphP[i].HeIII  = 0.000001 * y_fac;
      SphP[i].Ne     = SphP[i].HII + SphP[i].HeII + 2.0 * SphP[i].HeIII;
      SphP[i].nHeI   = SphP[i].HeI * nH_times_volume;
      SphP[i].nHeII  = SphP[i].HeII * nH_times_volume;
      SphP[i].nHeIII = SphP[i].HeIII * nH_times_volume;

      SphP[i].ne = SphP[i].Ne * nH_times_volume;
    }
}
#endif

void set_full_ionization_mrt(int i)
{
#ifndef MRT_NO_UV
  SphP[i].Ne   = (HYDROGEN_MASSFRAC + 1) / 2 / HYDROGEN_MASSFRAC;
  SphP[i].ne   = SphP[i].Ne * P[i].Mass;
  SphP[i].HI   = 0;
  SphP[i].nHI  = SphP[i].HI * P[i].Mass;
  SphP[i].HII  = 1;
  SphP[i].nHII = SphP[i].HII * P[i].Mass;
#ifdef MRT_INCLUDE_HE
  SphP[i].HeI    = 0;
  SphP[i].nHeI   = SphP[i].HeI * P[i].Mass;
  SphP[i].HeII   = 0;
  SphP[i].nHeII  = SphP[i].HeII * P[i].Mass;
  double y_frac  = (1.0 - HYDROGEN_MASSFRAC) / 4.0 / HYDROGEN_MASSFRAC;
  SphP[i].HeIII  = y_frac;
  SphP[i].nHeIII = SphP[i].HeIII * P[i].Mass;
#endif
#endif
}

#endif
