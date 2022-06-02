/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/tgchem/tgchem_step.c
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       Primordial chemistry and cooling network
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#ifdef TGCHEM_TEST
#include "tgchem_0d.h"
#else
#include "../allvars.h"
#include "../proto.h"
#endif

void tgchem_do_cvode_step(step_data* pstep, cell_data* pcell, double dt_step);
double tgchem_bisection_h2(double abh2, double abhii, double energy, double nh, double rad_heat_rate, double dt);
double tgchem_bisection_hii(double abh2, double abhii, double energy, double nh);
double tgchem_eq_h2(double abh2, double abhii, double energy, double nh);
double tgchem_quadratic(int branch, double a, double b, double c);
int tgchem_check_species(cell_data* pcell);
void tgchem_finish_step(step_data* pstep, cell_data* pcell);

void tgchem_step(step_data* pstep, cell_data* pcell)
{
  // Unpack cell data

  double nh       = pcell->nh;
  double dt_hydro = pcell->dt_hydro;

  // Do step

  tgchem_debug(0, dt_hydro, pstep, pcell);

  if(dt_hydro)
    {
      int eq_h2_flag, eq_hii_flag, flag_break = 0;
      double t0, t1, dt_step, dt_old, dt_sum = 0.;
      double abh2, abhii, energy;
      double tot_hii_rate, tot_rad_rate;
      double dt_h2, dt_hii, dt_rad;

      int count = 0;

      while(!flag_break)
        {
          dt_step = fmax(dt_hydro - dt_sum, 0.);

          // Compute rates

          tgchem_rates(1, pstep, pcell);

          tgchem_debug(1, dt_step, pstep, pcell);

          tgchem_rates_aux(pcell, dt_step);

          tgchem_debug(2, dt_step, pstep, pcell);

          // Unpack cell data

          eq_h2_flag  = pcell->eq_flag[1];
          eq_hii_flag = pcell->eq_flag[2];

          abh2   = pcell->species[1];
          abhii  = pcell->species[2];
          energy = pcell->species[TGCHEM_NUM_ABUNDANCES];

          tot_hii_rate = pcell->tot_rate[2];
          tot_rad_rate = pcell->tot_rad_rate;

          dt_h2  = pcell->dt_rate[1];
          dt_hii = pcell->dt_rate[2];
          dt_rad = pcell->dt_rad;

          // Determine step size

          dt_old = dt_step;

          if(eq_h2_flag)
            {
              if(!eq_hii_flag)
                dt_step = fmin(TGCHEM_TOL_EQ_STEPSIZE * dt_hii, dt_step);

              dt_step = fmin(TGCHEM_TOL_EQ_STEPSIZE * dt_rad, dt_step);
            }
          else
            {
              dt_step = fmin(dt_h2, dt_step);
              dt_step = fmin(dt_hii, dt_step);
              dt_step = fmin(dt_rad, dt_step);
            }

          if(dt_step == dt_old)
            flag_break = 1;

          tgchem_debug(3, dt_step, pstep, pcell);

          // Evolve

          if(eq_h2_flag)
            {
              t0 = second();

              if(eq_hii_flag)
                abhii = tgchem_bisection_hii(abh2, abhii, energy, nh);
              else
                abhii += tot_hii_rate * dt_step;

              energy = tgchem_bisection_h2(abh2, abhii, energy, nh, tot_rad_rate, dt_step);

              abh2 = tgchem_eq_h2(abh2, abhii, energy, nh);

              // Update cell data

              pcell->species[1]                     = abh2;
              pcell->species[2]                     = abhii;
              pcell->species[TGCHEM_NUM_ABUNDANCES] = energy;

              t1 = second();

              if(eq_hii_flag)
                {
                  pstep->num_eq_hii++;
                  pstep->dt_eq_hii += timediff(t0, t1);
                }
              else
                {
                  pstep->num_eq_h2++;
                  pstep->dt_eq_h2 += timediff(t0, t1);
                }
            }
          else
            {
              t0 = second();

              tgchem_do_cvode_step(pstep, pcell, dt_step);

              t1 = second();

              pstep->num_neq++;

              pstep->dt_neq += timediff(t0, t1);
            }

          count++;

          dt_sum += dt_step;

          tgchem_debug(4, dt_step, pstep, pcell);
        }

      // Compute H- abundance

      tgchem_rates(1, pstep, pcell);

      pcell->species[0] = pcell->abhm;

      // Fnish step

      tgchem_debug(5, dt_step, pstep, pcell);

      tgchem_finish_step(pstep, pcell);
    }
}

void tgchem_do_cvode_step(step_data* pstep, cell_data* pcell, double dt_step)
{
  int i;
  double tout;

  int id = pcell->cell_id;

  // Initialize user data

  pointer_data pointer;
  void* ppointer = &pointer;

  pointer.pstep = pstep;
  pointer.pcell = pcell;

  int flag = CVodeSetUserData(TGCD.CVODEMem, ppointer);

  if(flag)
    terminate("TGCHEM: Error in CVodeSetUserData! Returned flag: %d, Cell ID: %d\n", flag, id);

  // Initialize species vector

  N_Vector cv_species = N_VNew_Serial(TGCHEM_NUM_SPECIES);

  for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
    NV_Ith_S(cv_species, i) = pcell->species[i];

  // Do CVODE step

  double dt_step_neq = dt_step;
  double dt_sum_neq  = 0.;
  double species_save[TGCHEM_NUM_SPECIES];

  while(dt_sum_neq < dt_step)
    {
      for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
        species_save[i] = NV_Ith_S(cv_species, i);

      flag = CVodeReInit(TGCD.CVODEMem, 0, cv_species);

      if(flag)
        terminate("TGCHEM: Error in CVodeReInit! Returned flag: %d, Cell ID: %d\n", flag, id);

      flag = CVode(TGCD.CVODEMem, dt_step_neq, cv_species, &tout, CV_NORMAL);

      for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
        pcell->species[i] = NV_Ith_S(cv_species, i);

      if(!flag)
        flag = tgchem_check_species(pcell);

      if(flag)
        {
          for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
            NV_Ith_S(cv_species, i) = species_save[i];

          dt_step_neq /= 2.;
        }
      else
        {
          dt_sum_neq += dt_step_neq;
          dt_step_neq = dt_step - dt_sum_neq;
        }
    }

  // Clean up CVODE

  N_VDestroy_Serial(cv_species);
}

int tgchem_dspecies(double time, N_Vector cv_species, N_Vector cv_dspecies, void* ppointer)
{
  int i;
  step_data* pstep = ((pointer_data*)ppointer)->pstep;
  cell_data* pcell = ((pointer_data*)ppointer)->pcell;

  for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
    pcell->species[i] = NV_Ith_S(cv_species, i);

  tgchem_rates(0, pstep, pcell);

  for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
    NV_Ith_S(cv_dspecies, i) = pcell->tot_rate[i];

  NV_Ith_S(cv_dspecies, 0) = 0.;

  return 0;
}

double tgchem_bisection_h2(double abh2, double abhii, double energy, double nh, double rad_heat_rate, double dt)
{
  double abh2_eq, energy_new, energy_old, energy_min, energy_max, func;

  energy_min = TGCHEM_TOL_BISECTION_H2 * energy;
  energy_max = 1. / TGCHEM_TOL_BISECTION_H2 * energy;

  energy_new = energy;

  while(energy_new != energy_min && energy_new != energy_max)
    {
      energy_old = energy_new;

      abh2_eq = tgchem_eq_h2(abh2, abhii, energy_new, nh);

      func = energy - energy_new - TGCHEM_CHI_H2 * (abh2 - abh2_eq) * nh + rad_heat_rate * dt;

      if(func <= 0)
        energy_max = energy_new;
      else
        energy_min = energy_new;

      energy_new = (energy_min + energy_max) / 2.;

      if(fabs(energy_new - energy_old) < TGCHEM_TOL * energy_old)
        break;
    }

  return energy_new;
}

double tgchem_bisection_hii(double abh2, double abhii, double energy, double nh)
{
  double abhii_new, abhii_old, abhii_eq, abhii_min, abhii_max, func;
  double abhi, abe, mu, ntot, gamma, temp, K;

  abhii_min = TGCHEM_TOL_BISECTION_HII * abhii;
  abhii_max = 1. / TGCHEM_TOL_BISECTION_HII * abhii;
  abhii_max = fmin(abhii_max, TGCD.AbMax[2]);

  abhii_new = abhii;

  while(abhii_new != abhii_min && abhii_new != abhii_max)
    {
      abhii_old = abhii_new;

      tgchem_dep_vars(abh2, abhii_new, nh, &abhi, &abe, &mu, &ntot);

      tgchem_gamma_temp(abh2, abe, energy, ntot, &gamma, &temp);

      K = tgchem_chem_coeff(4, temp) / tgchem_chem_coeff(5, temp);

      abhii_eq = K / (1. + K);

      func = abhii_eq - abhii_new;

      if(func <= 0)
        abhii_max = abhii_new;
      else
        abhii_min = abhii_new;

      abhii_new = (abhii_min + abhii_max) / 2.;

      if(fabs(abhii_new - abhii_old) < TGCHEM_TOL * abhii_old)
        break;
    }

  return abhii_new;
}

double tgchem_eq_h2(double abh2, double abhii, double energy, double nh)
{
  double abh = fmax(1. - abhii, 0.);

  double abhi, abe, mu, ntot;

  tgchem_dep_vars(abh2, abhii, nh, &abhi, &abe, &mu, &ntot);

  double gamma, temp;

  tgchem_gamma_temp(abh2, abe, energy, ntot, &gamma, &temp);

  double K0 = tgchem_eq_coeff(0, temp) / nh;

  double abh2_new = tgchem_quadratic(1., 4., -(4. * abh + K0), abh * abh);

  abh2_new = tgchem_limit_abh2(abh2_new);

  return abh2_new;
}

double tgchem_quadratic(int branch, double a, double b, double c)
{
  double b2, d, quad;

  b2 = b * b;

  d = 4. * a * c;

  if(fabs(d) < 1e-10 * b2)
    {
      quad = d / 2. / fabs(b);

      if(!branch)
        quad *= -1.;
    }
  else
    {
      quad = sqrt(b2 - d);

      if(branch)
        quad *= -1.;

      quad -= b;
    }

  quad /= 2. * a;

  return quad;
}

int tgchem_check_species(cell_data* pcell)
{
  int i;

  // Unpack cell data

  int cell_id = pcell->cell_id;

  double species[TGCHEM_NUM_SPECIES];

  for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
    species[i] = pcell->species[i];

  // Do checks

  for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
    {
      if(species[i] != species[i])
        terminate("Species %d is NaN! ID: %d", i, cell_id);

      if(species[i] < 0.)
        {
          if(species[i] < -NV_Ith_S(TGCD.SpeciesTol, i))
            {
              printf("Warning: species %d (y = %g, ID: %d) is below allowed tolerance. Reducing internal timestep.\n", i, species[i],
                     cell_id);

              return 1;
            }

          species[i] = 0.;
        }
    }

  for(i = 0; i < TGCHEM_NUM_ABUNDANCES; i++)
    if(species[i] > TGCD.AbMax[i])
      {
        if((species[i] - TGCD.AbMax[i]) / species[i] > TGCHEM_TOL)
          {
            printf("Warning: species %d (y = %g, ID: %d) is above allowed tolerance. Reducing internal timestep.\n", i, species[i],
                   cell_id);

            return 2;
          }

        species[i] = TGCD.AbMax[i];
      }

  // Update cell data

  for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
    pcell->species[i] = species[i];

  return 0;
}

void tgchem_finish_step(step_data* pstep, cell_data* pcell)
{
#ifndef TGCHEM_TEST
  int i;

  // Unpack cell data

  int idx        = pcell->cell_idx;
  double rho_int = pcell->rho_int;

  double species[TGCHEM_NUM_SPECIES];

  for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
    species[i] = pcell->species[i];

  // Update particle quantities

  var_data var;

  tgchem_var_data(pcell, &var);

  SphP[idx].Gamma = var.gamma;

  for(i = 0; i < TGCHEM_NUM_ABUNDANCES; i++)
    SphP[idx].Abund[i] = species[i];

  double Utherm_new = TGCD.EnergyConvFac * species[TGCHEM_NUM_ABUNDANCES] / rho_int;
  Utherm_new        = fmax(Utherm_new, All.MinEgySpec);

  double du = Utherm_new - SphP[idx].Utherm;

  SphP[idx].Utherm += du;

  SphP[idx].Energy += du * P[idx].Mass * pow(All.cf_atime, 2);
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  SphP[idx].A       = (SphP[idx].Gamma - 1) * SphP[idx].Utherm / pow(SphP[idx].Density * All.cf_a3inv, GAMMA_MINUS1);
  SphP[idx].Entropy = log(SphP[idx].A) * P[idx].Mass;
#endif
  set_pressure_of_cell(idx);
#endif

  pstep->num_cells++;

  return;
}
