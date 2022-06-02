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

#if defined(MRT_IR) || defined(MRT_UV_ONLY_DUST) || defined(MRT_CHEM_SG)
#ifdef MRT_IR_GRAIN_KAPPA
int IR_N_pts;
double *IR_logT, *IR_logkappaP, *IR_logkappaR;
gsl_interp_accel *accIR_kappaP;
gsl_interp_accel *accIR_kappaR;
gsl_spline *splineIR_kappaP;
gsl_spline *splineIR_kappaR;

/* Read in tabulated grain opacity (both Planck and Rosseland averaged) values
 * at various temperatures.  Data is in log cgs units. */
void read_grain_kappa_data(void)
{
  char setname[MAXLEN_PATH];

  hid_t file_id = my_H5Fopen(All.GrainKappaPath, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t dataset;

  sprintf(setname, "/NData");
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &IR_N_pts, setname);
  my_H5Dclose(dataset, setname);

  IR_logT = (double *)mymalloc("IR_logT", IR_N_pts * sizeof(double));
  sprintf(setname, "/LogTemperature");
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, IR_logT, setname);
  my_H5Dclose(dataset, setname);

  IR_logkappaP = (double *)mymalloc("IR_logkappaP", IR_N_pts * sizeof(double));
  sprintf(setname, "/LogKappaP");
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, IR_logkappaP, setname);
  my_H5Dclose(dataset, setname);

  IR_logkappaR = (double *)mymalloc("IR_logkappaR", IR_N_pts * sizeof(double));
  sprintf(setname, "/LogKappaR");
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, IR_logkappaR, setname);
  my_H5Dclose(dataset, setname);

  accIR_kappaP    = gsl_interp_accel_alloc();
  accIR_kappaR    = gsl_interp_accel_alloc();
  splineIR_kappaP = gsl_spline_alloc(gsl_interp_linear, IR_N_pts);
  splineIR_kappaR = gsl_spline_alloc(gsl_interp_linear, IR_N_pts);

  gsl_spline_init(splineIR_kappaP, IR_logT, IR_logkappaP, IR_N_pts);
  gsl_spline_init(splineIR_kappaR, IR_logT, IR_logkappaR, IR_N_pts);

  my_H5Fclose(file_id, All.GrainKappaPath);
}
#endif

/*Set kappa * rho for every cell*/
void set_kappa_times_rho_IR(
    int i,
    struct sph_particle_data *kappaSphP) /*Expand this function later to get kappa*rho from the properties of the dust in the cell*/
{
#if defined(DUST_LIVE) && defined(DL_RADIATION_PRESSURE)
  /* If modeling dust with live simulation particles, no need to treat dust
   * opacity in gas cells. */
  for(int num1 = 0; num1 < MRT_BINS; num1++)
    {
      kappaSphP[i].KappaIR_R[num1] = kappaSphP[i].KappaIR_P[num1] = 0.0;
    }
  return;
#endif

#ifndef MRT_IR_GRAIN_KAPPA
  for(int num1 = 0; num1 < MRT_BINS; num1++)
    {
#ifdef MRT_CHEM_SG
      double val1 = All.UVKappa;
#ifdef MRT_IR
      if(num1 == UV_BINS)
        val1 = All.IRKappa;
#endif
      kappaSphP[i].KappaIR_R[num1] = kappaSphP[i].KappaIR_P[num1] =
          kappaSphP[i].Density * val1 * All.UnitDensity_in_cgs * All.UnitLength_in_cm;
#else

#ifdef MRT_LEVITATION_TEST
      double meanweight = 2.33;
      double temperature =
          kappaSphP[i].Utherm * All.UnitEnergy_in_cgs / All.UnitMass_in_g * GAMMA_MINUS1 * meanweight * PROTONMASS / BOLTZMANN;

      double val;

      // if(temperature>1000.0)
      // val = 316.0 * All.UnitMass_in_g / pow(All.UnitLength_in_cm, 2.0) ;
      // else
      val = 0.0316 * (temperature / 10.0) * (temperature / 10.0) * All.UnitMass_in_g / pow(All.UnitLength_in_cm, 2.0) *
            kappaSphP[i].Density;

      kappaSphP[i].KappaIR_R[num1] = val;
      kappaSphP[i].KappaIR_P[num1] = 0.1 * val / 0.0316;
#else
      double nH = HYDROGEN_MASSFRAC * kappaSphP[i].Density / PROTONMASS * All.UnitMass_in_g / All.HubbleParam;
      nH *= pow(All.UnitLength_in_cm, -3);

      double meanweight = 1.0;
      double temperature =
          kappaSphP[i].Utherm * All.UnitEnergy_in_cgs / All.UnitMass_in_g * GAMMA_MINUS1 * meanweight * PROTONMASS / BOLTZMANN;

      double val1 = 0.0;
      double tt   = temperature / 1e5;
      // if()
      if(num1 < UV_BINS)
        val1 = 0.1;
      else
        val1 = 0.1;

      kappaSphP[i].KappaIR_R[num1] = kappaSphP[i].KappaIR_P[num1] = val1 * All.UnitLength_in_cm;
#endif
#endif

#ifdef MRT_BH
      kappaSphP[i].KappaIR_R[num1] = kappaSphP[i].KappaIR_P[num1] =
          (1000.0 * All.UnitMass_in_g / All.UnitLength_in_cm / All.UnitLength_in_cm) * kappaSphP[i].Density;
#endif
    }
#else
  double DGR = 0.0;
  for(int chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
    {
      for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
        {
          DGR += SphP[i].MetalsDustFraction[chan][k];
        }
    }

  double meanweight  = 4. / (3 * HYDROGEN_MASSFRAC + 1 + 4 * HYDROGEN_MASSFRAC * SphP[i].Ne);
  double temperature = SphP[i].Utherm * All.UnitEnergy_in_cgs / All.UnitMass_in_g * GAMMA_MINUS1 * meanweight * PROTONMASS / BOLTZMANN;

  double logT = log10(temperature);
  if(logT < IR_logT[0])
    {
      logT = IR_logT[0];
    }
  if(logT > IR_logT[IR_N_pts - 1])
    {
      logT = IR_logT[IR_N_pts - 1];
    }

  for(int num1 = 0; num1 < MRT_BINS; num1++)
    {
      /* Tabulated grain opacities given in cm^2 per unit mass of dust.
       * Multiply by dust-to-gas ratio to get per unit mass of gas. */
      double logkappaP_cgs    = gsl_spline_eval(splineIR_kappaP, logT, accIR_kappaP);
      double kappaP_cgs       = pow(10.0, logkappaP_cgs);
      double kappaP_code      = kappaP_cgs * All.UnitMass_in_g / (All.UnitLength_in_cm * All.UnitLength_in_cm);
      SphP[i].KappaIR_P[num1] = kappaP_code * DGR * SphP[i].Density;

      double logkappaR_cgs    = gsl_spline_eval(splineIR_kappaR, logT, accIR_kappaR);
      double kappaR_cgs       = pow(10.0, logkappaR_cgs);
      double kappaR_code      = kappaR_cgs * All.UnitMass_in_g / (All.UnitLength_in_cm * All.UnitLength_in_cm);
      SphP[i].KappaIR_R[num1] = kappaR_code * DGR * SphP[i].Density;
    }
#endif
}
#endif
