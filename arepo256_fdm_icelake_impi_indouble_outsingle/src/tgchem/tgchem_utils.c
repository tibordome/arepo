/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/tgchem/tgchem_utils.c
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

void tgchem_func(double *f, double *x, int n);
void tgchem_jacobian(double *J, double *x, int n);
void tgchem_root(double *A, double *x, double *b, int n);
int tgchem_lu_decomp(double *A, int n);
void tgchem_lower_triangular(double *L, double *b, double *x, int n);
int tgchem_upper_triangular(double *U, double *b, double *x, int n);

void tgchem_step_data(step_data *pstep)
{
#ifdef TGCHEM_TEST
  double temp_cmb = TEMP_CMB * (1. + TGCD.RedShiftInit);
#else
  set_cosmo_factors_for_current_time();

  double temp_cmb   = TEMP_CMB * (1. + All.cf_redshift);
#endif

  pstep->num_neq        = 0;
  pstep->num_eq_h2      = 0;
  pstep->num_eq_hii     = 0;
  pstep->num_rate_calls = 0;
  pstep->num_cells      = 0;
  pstep->dt_neq         = 0.;
  pstep->dt_eq_h2       = 0.;
  pstep->dt_eq_hii      = 0.;
  pstep->temp_cmb       = temp_cmb;
}

void tgchem_cell_data(step_data *pstep, cell_data *pcell, int idx)
{
  int i;
  double species[TGCHEM_NUM_SPECIES];

  // Compute cell data

#ifdef TGCHEM_TEST
  int cell_idx      = idx;
  int cell_id       = pstep->count;
  int cell_task     = 0;
  double rho_int    = 0.;
  double nh         = pstep->nh;
  double hydro_rate = pstep->hydro_rate;
  double dt_hydro   = pstep->dt_hydro;
  double divv       = 0.;

  for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
    species[i] = pcell->species[i];
#else
  int cell_idx      = idx;
  int cell_id       = P[idx].ID;
  int cell_task     = ThisTask;
  double rho_int    = SphP[idx].Density * All.cf_a3inv * All.HubbleParam * All.HubbleParam;
  double rho        = rho_int * All.UnitDensity_in_cgs;
  double nh         = HYDROGEN_MASSFRAC * rho / PROTONMASS;
  double hydro_rate = 0.;
  double dt_hydro   = (P[idx].TimeBinHydro ? (((integertime)1) << P[idx].TimeBinHydro) : 0) * All.Timebase_interval;
  dt_hydro *= All.UnitTime_in_s / All.cf_hubble_a / All.HubbleParam;
  //#ifndef FORCE_EQUAL_TIMESTEPS
  // dt_hydro *= pow(2., 2. * get_random_number() - 1.);
  //#endif
  double divv = SphP[idx].DivVel * (1. + All.cf_redshift) * All.HubbleParam / All.UnitTime_in_s;

  for(i = 0; i < TGCHEM_NUM_ABUNDANCES; i++)
    species[i] = SphP[idx].Abund[i];

  species[TGCHEM_NUM_ABUNDANCES] = SphP[idx].Utherm * rho_int / TGCD.EnergyConvFac;
#endif

  // Update cell data

  pcell->cell_idx   = cell_idx;
  pcell->cell_id    = cell_id;
  pcell->cell_task  = cell_task;
  pcell->rho_int    = rho_int;
  pcell->nh         = nh;
  pcell->nh2        = nh * nh;
  pcell->nh3        = nh * nh * nh;
  pcell->hydro_rate = hydro_rate;
  pcell->dt_hydro   = dt_hydro;
  pcell->divv       = divv;

  for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
    pcell->species[i] = species[i];

  // Compute var data

  var_data var;

  tgchem_var_data(pcell, &var);

  // Unpack var data

  double abh2  = var.abh2;
  double mu    = var.mu;
  double gamma = var.gamma;
  double temp  = var.temp;

  // Compute additional variables

  double t_ff      = tgchem_t_ff(nh);
  double vth       = sqrt(BOLTZMANN * temp / PROTONMASS);
  double csnd      = sqrt(gamma * BOLTZMANN * temp / mu / PROTONMASS);
  double ljeans    = csnd * t_ff;
  double lsob      = ((divv == 0.) ? ljeans : fmin(vth / fabs(divv), ljeans));
  double h2_column = tgchem_limit_h2_column(abh2 * nh * lsob);

  // Update cell data

  pcell->t_ff      = t_ff;
  pcell->vth       = vth;
  pcell->ljeans    = ljeans;
  pcell->h2_column = h2_column;
}

void tgchem_var_data(cell_data *pcell, var_data *pvar)
{
  // Unpack cell data

  double nh     = pcell->nh;
  double abhm   = tgchem_limit_abhm(pcell->species[0]);
  double abh2   = tgchem_limit_abh2(pcell->species[1]);
  double abhii  = tgchem_limit_abhii(pcell->species[2]);
  double energy = tgchem_limit_energy(pcell->species[TGCHEM_NUM_ABUNDANCES]);

  // Compute dependent variables

  double abhi, abe, mu, ntot;

  tgchem_dep_vars(abh2, abhii, nh, &abhi, &abe, &mu, &ntot);

  // Compute gamma and temperature

  double gamma, temp;

  tgchem_gamma_temp(abh2, abe, energy, ntot, &gamma, &temp);

  // Update var data

  pvar->abhm   = abhm;
  pvar->abh2   = abh2;
  pvar->abhii  = abhii;
  pvar->energy = energy;
  pvar->abhi   = abhi;
  pvar->abe    = abe;
  pvar->mu     = mu;
  pvar->ntot   = ntot;
  pvar->gamma  = gamma;
  pvar->temp   = temp;
}

void tgchem_dep_vars(double abh2, double abhii, double nh, double *abhi, double *abe, double *mu, double *ntot)
{
  *abhi = tgchem_limit_abhi(1. - 2. * abh2 - abhii);
  *abe  = abhii;
  *mu   = (1. + 4. * HE_ABUND) / fmax(1. + HE_ABUND - abh2 + abhii, 0.);
  *ntot = fmax(1. + HE_ABUND - abh2 + *abe, 0.) * nh;
}

void tgchem_gamma_temp(double abh2, double abe, double energy, double ntot, double *gamma, double *temp)
{
  // Compute gamma

  double gamma_est = GAMMA;
  double temp_est  = tgchem_limit_temp((gamma_est - 1.) * energy / BOLTZMANN / ntot);

  double x = 6.1e3 / temp_est;

  double abn = fmax(1. + HE_ABUND - abh2 + abe, 0.);
  double ab  = fmax(1. + HE_ABUND - 2 * abh2 + abe, 0.);

  double gvar_h  = 1. / GAMMA_MINUS1;
  double gvar_h2 = 5. / 2. + pow(x, 2) * exp(x) / pow(exp(x) - 1., 2);

  *gamma = 1. + (abn / (ab * gvar_h + abh2 * gvar_h2));

  // Compute temperature

  *temp = tgchem_limit_temp((*gamma - 1.) * energy / BOLTZMANN / ntot);
}

double tgchem_eq_coeff(int i, double temp)
{
  int temp_idx;
  double dtemp, coeff;

  temp_idx = tgchem_compute_temp_idx(temp);

  dtemp = tgchem_compute_dtemp(temp, temp_idx);

  coeff = TGCD.EqTable[i * TGCHEM_NUM_TEMP + temp_idx] + dtemp * TGCD.DEqTable[i * TGCHEM_NUM_TEMP + temp_idx];

  return coeff;
}

double tgchem_chem_coeff(int i, double temp)
{
  int temp_idx;
  double dtemp, coeff;

  temp_idx = tgchem_compute_temp_idx(temp);

  dtemp = tgchem_compute_dtemp(temp, temp_idx);

  coeff = TGCD.ChemTable[i * TGCHEM_NUM_TEMP + temp_idx] + dtemp * TGCD.DChemTable[i * TGCHEM_NUM_TEMP + temp_idx];

  return coeff;
}

double tgchem_cool_coeff(int i, double temp)
{
  int temp_idx;
  double dtemp, coeff;

  temp_idx = tgchem_compute_temp_idx(temp);

  dtemp = tgchem_compute_dtemp(temp, temp_idx);

  coeff = TGCD.CoolTable[i * TGCHEM_NUM_TEMP + temp_idx] + dtemp * TGCD.DCoolTable[i * TGCHEM_NUM_TEMP + temp_idx];

  return coeff;
}

double tgchem_opac_coeff(int i, double temp)
{
  int temp_idx;
  double dtemp, coeff;

  temp_idx = tgchem_compute_temp_idx(temp);

  dtemp = tgchem_compute_dtemp(temp, temp_idx);

  coeff = TGCD.OpacTable[i * TGCHEM_NUM_TEMP + temp_idx] + dtemp * TGCD.DOpacTable[i * TGCHEM_NUM_TEMP + temp_idx];

  return coeff;
}

int tgchem_compute_temp_idx(double temp)
{
  int temp_idx = imin((int)(log10(temp / TGCHEM_TEMP_MIN) / TGCHEM_LOG_DTEMP), TGCHEM_NUM_TEMP - 1);
  return temp_idx;
}

double tgchem_compute_dtemp(double temp, int temp_idx)
{
  double dtemp = temp - TGCD.TempTable[temp_idx];
  return dtemp;
}

double tgchem_limit_abhm(double abhm)
{
  abhm = fmin(fmax(abhm, 0.), TGCD.AbMax[0]);
  return abhm;
}

double tgchem_limit_abh2(double abh2)
{
  abh2 = fmin(fmax(abh2, 0.), TGCD.AbMax[1]);
  return abh2;
}

double tgchem_limit_abhii(double abhii)
{
  abhii = fmin(fmax(abhii, 0.), TGCD.AbMax[2]);
  return abhii;
}

double tgchem_limit_energy(double energy)
{
  energy = fmax(energy, 0.);
  return energy;
}

double tgchem_limit_abhi(double abhi)
{
  abhi = fmax(abhi, 0.);
  return abhi;
}

double tgchem_limit_temp(double temp)
{
  temp = fmin(fmax(temp, TGCHEM_TEMP_MIN), TGCHEM_TEMP_MAX);
  return temp;
}

double tgchem_limit_h2_column(double h2_column)
{
  h2_column = fmin(fmax(h2_column, TGCHEM_H2_COLUMN_MIN), TGCHEM_H2_COLUMN_MAX);
  return h2_column;
}

double tgchem_t_ff(double nh)
{
  double rho  = PROTONMASS * nh / HYDROGEN_MASSFRAC;
  double t_ff = sqrt(3. * M_PI / 32. / GRAVITY / rho);
  return t_ff;
}

void tgchem_lu_solve(void)
{
  int i, j;

  int n  = 2;
  int n2 = n * n;

  /*
  double* A  = (double *) malloc(n2 * sizeof(double));
  double* A_old  = (double *) malloc(n2 * sizeof(double));
  double* A_new  = (double *) malloc(n2 * sizeof(double));
  double* x1  = (double *) malloc(n * sizeof(double));
  double* x2  = (double *) malloc(n * sizeof(double));
  double* B  = (double *) malloc(n * sizeof(double));
  double* B1  = (double *) malloc(n * sizeof(double));
  double* B2  = (double *) malloc(n * sizeof(double));

  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      {
        A[i * n + j] = rand() % 3;
        A_old[i * n + j] = A[i * n + j];

        printf("A[%d][%d] = %g\n", i, j, A[i * n + j]);
      }

  for(i = 0; i < n; i++)
    {
      B[i] = i + 1;

      printf("B[%d] = %g\n", i, B[i]);
    }

  int ret;
  */
  double *x     = (double *)malloc(n * sizeof(double));
  double *x_old = (double *)malloc(n * sizeof(double));
  double *y     = (double *)malloc(n * sizeof(double));
  double *f     = (double *)malloc(n * sizeof(double));
  double *J     = (double *)malloc(n2 * sizeof(double));
  double *J_old = (double *)malloc(n2 * sizeof(double));

  x[0] = -10.;
  x[1] = -5.;

  int iter = 0;
  int flag = 1;

  while(flag)
    {
      for(i = 0; i < n; i++)
        x_old[i] = x[i];

      tgchem_func(f, x, n);

      for(i = 0; i < n; i++)
        printf("f[%d] = %g\n", i, f[i]);

      for(i = 0; i < n; i++)
        f[i] *= -1.;

      tgchem_jacobian(J, x, n);

      memcpy(J_old, J, n2 * sizeof(double));

      // for(i = 0; i < n; i++)
      // for(j = 0; j < n; j++)
      // printf("J[%d][%d] = %g\n", i, j, J_old[i * n + j]);

      tgchem_root(J, x, f, n);

      for(i = 0; i < n; i++)
        {
          y[i] = 0.;

          for(j = 0; j < n; j++)
            y[i] += J_old[i * n + j] * x[j];
        }

      for(i = 0; i < n; i++)
        printf("y[%d] = %g\n", i, y[i]);

      for(i = 0; i < n; i++)
        x[i] = x_old[i] + x[i];

      flag = 0;

      for(i = 0; i < n; i++)
        if(fabs(x[i] - x_old[i]) > 1e-5)
          flag = 1;

      for(i = 0; i < n; i++)
        printf("x[%d] = %g\n", i, x[i]);

      printf("\n\n");

      iter++;
    }
  /*
  ret = tgchem_lu_decomp(A, n);

  if(ret == -1)
    terminate("Matrix A is singular!\n");

  tgchem_lower_triangular(A, B, x1, n);

  for(i = 0; i < n; i++)
    printf("x1[%d] = %g\n", i, x1[i]);

  ret = tgchem_upper_triangular(A, B, x2, n);

  if(ret == -1)
    terminate("Matrix A is singular!\n");

  for(i = 0; i < n; i++)
    printf("x2[%d] = %g\n", i, x2[i]);

  double val;

  double *L  = (double *) malloc(n2 * sizeof(double));
  double *U  = (double *) malloc(n2 * sizeof(double));

  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      {
        if(i > j)
          val = A[i * n + j];
        else if(i < j)
          val = 0.;
        else
          val = 1.;

        L[i * n + j] = val;
      }

  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      {
        if(i <= j)
          val = A[i * n + j];
        else
          val = 0.;

        U[i * n + j] = val;
      }

  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      {
        A_new[i * n + j] = 0.;

        for(k = 0; k < n; k++)
          A_new[i * n + j] += L[i * n + k] * U[k * n + j];
      }

  for(i = 0; i < n; i++)
    {
      B1[i] = 0.;

      for(j = 0; j < n; j++)
        B1[i] += L[i * n + j] * x1[j];
    }

  for(i = 0; i < n; i++)
    {
      B2[i] = 0.;

      for(j = 0; j < n; j++)
        B2[i] += U[i * n + j] * x2[j];
    }





  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      printf("L[%d][%d] = %g\n", i, j, L[i * n + j]);

  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      printf("U[%d][%d] = %g\n", i, j, U[i * n + j]);

  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      printf("A_new[%d][%d] = %g\n", i, j, A_new[i * n + j]);

  for(i = 0; i < n; i++)
    printf("B1[%d] = %g\n", i, B1[i]);

  for(i = 0; i < n; i++)
    printf("B2[%d] = %g\n", i, B2[i]);

  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      {
        if(A_old[i * n + j] - A_new[i * n + j] > 1e-10)
          terminate("error!\n");
      }

  for(i = 0; i < n; i++)
    if(B[i] - B1[i] > 1e-10)
      terminate("error!\n");

  for(i = 0; i < n; i++)
    if(B[i] - B2[i] > 1e-10)
      terminate("error!\n");
  */
}

void tgchem_jacobian(double *J, double *x, int n)
{
  int i, j, idx;

  double eps    = 0.1;
  double *delta = (double *)malloc(n * sizeof(double));
  double *x1    = (double *)malloc(n * sizeof(double));
  double *f     = (double *)malloc(n * sizeof(double));
  double *f1    = (double *)malloc(n * sizeof(double));
  double *f_ij  = (double *)malloc(n * n * sizeof(double));

  tgchem_func(f, x, n);

  for(i = 0; i < n; i++)
    {
      for(j = 0; j < n; j++)
        x1[j] = x[j];

      delta[i] = eps * fabs(x[i]);

      x1[i] = x[i] + delta[i];

      tgchem_func(f1, x1, n);

      for(j = 0; j < n; j++)
        f_ij[j * n + i] = f1[j];
    }

  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      {
        idx = i * n + j;

        J[idx] = (f_ij[idx] - f[i]) / delta[j];
      }

  free(delta);
  free(x1);
  free(f);
  free(f1);
  free(f_ij);
}

void tgchem_func(double *f, double *x, int n)
{
  f[0] = 1. * (1. - x[0]);
  f[1] = 10. * (x[1] - x[0] * x[0]);
  /*
  f[0] = x[0] * x[0] + x[1] * x[1];
  f[1] = x[0] + x[1];
  */
}

void tgchem_root(double *A, double *x, double *b, int n)
{
  int ret = tgchem_lu_decomp(A, n);

  if(ret == -1)
    terminate("Matrix A is singular in LU decomposition!\n");

  double *y = (double *)malloc(n * sizeof(double));

  tgchem_lower_triangular(A, b, y, n);

  ret = tgchem_upper_triangular(A, y, x, n);

  if(ret == -1)
    terminate("Matrix A is singular in upper triangular solver!\n");
}

int tgchem_lu_decomp(double *A, int n)
{
  // Doolittle LU decomposition

  int i, j, k, p;
  double *p_k, *p_row, *p_col;

  for(k = 0, p_k = A; k < n; p_k += n, k++)
    {
      for(j = k; j < n; j++)
        {
          for(p = 0, p_col = A; p < k; p_col += n, p++)
            *(p_k + j) -= *(p_k + p) * *(p_col + j);
        }

      if(*(p_k + k) == 0.)
        return -1;

      for(i = k + 1, p_row = p_k + n; i < n; p_row += n, i++)
        {
          for(p = 0, p_col = A; p < k; p_col += n, p++)
            *(p_row + k) -= *(p_row + p) * *(p_col + k);

          *(p_row + k) /= *(p_k + k);
        }
    }

  return 0;
}

void tgchem_lower_triangular(double *L, double *b, double *x, int n)
{
  // Solves Lx = b

  int i, j;

  x[0] = b[0];

  for(j = 1, L += n; j < n; L += n, j++)
    for(i = 0, x[j] = b[j]; i < j; i++)
      x[j] -= x[i] * *(L + i);
}

int tgchem_upper_triangular(double *U, double *b, double *x, int n)
{
  // Solves Ux = b

  int i, j;

  for(j = n - 1, U += n * (n - 1); j >= 0; U -= n, j--)
    {
      if(*(U + j) == 0.)
        return -1;

      x[j] = b[j];

      for(i = j + 1; i < n; i++)
        x[j] -= x[i] * *(U + i);

      x[j] /= *(U + j);
    }

  return 0;
}

void tgchem_debug(int mode, double dt, step_data *pstep, cell_data *pcell)
{
  int i, j;

  if(TGCD.DebugFlag > -1)
    {
      int id = pcell->cell_id;

      if(id == TGCD.DebugFlag)
        {
          var_data var;

          tgchem_var_data(pcell, &var);

          if(mode == 0)
            {
              printf("\nBeginning of step:\n");

              printf("NH: %g, Gamma: %g, Temp: %g, dt_hydro = %g\n", pcell->nh, var.gamma, var.temp, dt);

              printf("Species: ");

              for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
                printf("%g ", pcell->species[i]);

              printf("\n");
            }
          else if(mode == 1)
            {
              printf("\nBefore substep:\n");

              printf("NH: %g, Gamma: %g, Temp: %g, dt_step = %g\n", pcell->nh, var.gamma, var.temp, dt);

              printf("Species: ");

              for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
                printf("%g ", pcell->species[i]);

              printf("\n");

              printf("Rates:\n");

              for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
                {
                  printf("\n");

                  for(j = 0; j < TGCHEM_NUM_CHEM_RATES; j++)
                    if(pcell->prate[i][j])
                      printf("\t\t Positive rate %d %d: %g\n", i, j, pcell->prate[i][j]);

                  for(j = 0; j < TGCHEM_NUM_CHEM_RATES; j++)
                    if(pcell->nrate[i][j])
                      printf("\t\t Negative rate %d %d: %g\n", i, j, pcell->nrate[i][j]);
                }

              printf("\n");

              for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
                printf("\t\t Total rate %d: %g\n", i, pcell->tot_rate[i]);

              printf("\n");

              printf("\t\t Total chem rate: %g\n", pcell->tot_chem_rate);
              printf("\t\t Total rad rate: %g\n", pcell->tot_rad_rate);
            }
          else if(mode == 2)
            {
              printf("\n");

              for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
                printf("\t\t Timescale %d: %g\n", i, pcell->dt_rate[i]);

              printf("\n");

              printf("\t\t Timescale Chem: %g\n", pcell->dt_chem);
              printf("\t\t Timescale Rad: %g\n", pcell->dt_rad);

              printf("\n");

              for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
                printf("\t\t Eq flags %d: %d\n", i, pcell->eq_flag[i]);

              printf("\n");
            }
          else if(mode == 3)
            {
              printf("\nLimiting time step: dt_step = %g\n", dt);
            }
          else if(mode == 4)
            {
              printf("\nAfter Substep:\n");

              printf("NH: %g, Gamma: %g, Temp: %g\n", pcell->nh, var.gamma, var.temp);

              printf("Species: ");

              for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
                printf("%g ", pcell->species[i]);

              printf("\n");
            }
          else if(mode == 5)
            {
              printf("\nEnd of step:\n");

              printf("NH: %g, Gamma: %g, Temp: %g\n", pcell->nh, var.gamma, var.temp);

              printf("Species: ");

              for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
                printf("%g ", pcell->species[i]);

              printf("\n");
            }
          else
            {
              terminate("TGCHEM: Unknown mode in debug output!\n");
            }
        }
    }
}
