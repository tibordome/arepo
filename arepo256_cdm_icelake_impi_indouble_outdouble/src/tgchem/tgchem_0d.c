/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/tgchem/tgchem_test.c
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

#include "tgchem_0d.h"

double WallClockTime;

struct TGCD_struct TGCD;

void tgchem_0d_init(step_data* pstep, cell_data* pcell);
void tgchem_0d_hydro(step_data* pstep, cell_data* pcell);
void tgchem_0d_out(step_data* pstep, cell_data* pcell, out_data* pout);
void tgchem_0d_write_outfile(step_data* pstep, out_data* pout);

void tgchem_0d_pars(void)
{
  TGCD.ChemMode      = 0;
  TGCD.ChemIOMode    = 0;
  TGCD.ChemH2Mode    = 0;
  TGCD.ChemInitAbH2  = 6.6e-7;
  TGCD.ChemInitAbHII = 2.6e-4;
  TGCD.ChemJ21       = 1e5;

  TGCD.DebugFlag = -1;

  TGCD.CollapseFac  = 1.;
  TGCD.TestStepSize = 1e-1;
  TGCD.NHInit       = 1e-3;
  TGCD.NHFinal      = 1e25;
  TGCD.TempInit     = 2e2;
  TGCD.RedShiftInit = 10.;
}

int main(int argc, char** argv)
{
  step_data step;
  cell_data cell;
  out_data out;

  double t0 = second();

  printf("starting...\n");

  tgchem_0d_pars();

  tgchem_begrun();

  tgchem_init_cvode();

  tgchem_0d_init(&step, &cell);

  tgchem_0d_out(&step, &cell, &out);

  // tgchem_lu_solve();

  // return 0;

  while(step.nh < step.nh_final)
    {
      tgchem_0d_hydro(&step, &cell);

      tgchem_step_data(&step);

      tgchem_cell_data(&step, &cell, 0);

      tgchem_step(&step, &cell);

      tgchem_0d_out(&step, &cell, &out);
    }

  tgchem_0d_write_outfile(&step, &out);

  tgchem_finish_cvode();

  double t1 = second();

  printf("%d iterations, took %g seconds. Done!\n", step.count, timediff(t0, t1));

  return 0;
}

void tgchem_0d_init(step_data* pstep, cell_data* pcell)
{
  double nh       = TGCD.NHInit;
  double nh_final = TGCD.NHFinal;
  double temp     = TGCD.TempInit;
  double abh2     = TGCD.ChemInitAbH2;
  double abhii    = TGCD.ChemInitAbHII;

  double dnh_dt     = TGCD.CollapseFac * nh / tgchem_t_ff(nh);
  double energy     = BOLTZMANN * temp * (1. + HE_ABUND) * nh / GAMMA_MINUS1;
  double hydro_rate = GAMMA * energy * dnh_dt / nh;

  // Initialize step data

  pstep->count      = 0;
  pstep->nh         = nh;
  pstep->nh_final   = nh_final;
  pstep->hydro_rate = hydro_rate;
  pstep->dt_hydro   = 0.;

  // Initialize cell data

  pcell->species[0] = 0.;
  pcell->species[1] = abh2;
  pcell->species[2] = abhii;
  pcell->species[3] = energy;

  // Compute self-consistent H- abundance

  tgchem_step_data(pstep);

  tgchem_cell_data(pstep, pcell, 0);

  tgchem_rates(1, pstep, pcell);

  pcell->species[0] = pcell->abhm;
}

void tgchem_0d_hydro(step_data* pstep, cell_data* pcell)
{
  // Unpack step data

  double nh       = pstep->nh;
  double nh_final = pstep->nh_final;

  // Unpack cell data

  double energy = pcell->species[TGCHEM_NUM_ABUNDANCES];

  // Do hydro

  double dnh_dt     = TGCD.CollapseFac * nh / tgchem_t_ff(nh);
  double hydro_rate = GAMMA * energy * dnh_dt / nh;
  double dt_hydro   = TGCD.TestStepSize * nh / dnh_dt;

  if(nh + dnh_dt * dt_hydro > nh_final)
    {
      dt_hydro = (nh_final - nh) / dnh_dt;

      nh = nh_final;
    }
  else
    nh += dnh_dt * dt_hydro;

  // Update step data

  pstep->nh         = nh;
  pstep->hydro_rate = hydro_rate;
  pstep->dt_hydro   = dt_hydro;

  // Update cell data

  pcell->species[TGCHEM_NUM_ABUNDANCES] += hydro_rate * pstep->dt_hydro;
}

void tgchem_0d_out(step_data* pstep, cell_data* pcell, out_data* pout)
{
  int i, j;

  // Compute data

  tgchem_rates(1, pstep, pcell);

  var_data var;

  tgchem_var_data(pcell, &var);

  // Unpack step data

  int count           = pstep->count;
  long num_sub_steps  = pstep->num_neq + pstep->num_eq_h2 + pstep->num_eq_hii;
  long num_rate_calls = pstep->num_rate_calls;

  // Unpack cell data

  int eq_flag_h2  = pcell->eq_flag[1];
  int eq_flag_hii = pcell->eq_flag[2];
  double nh       = pcell->nh;

  double species[TGCHEM_NUM_SPECIES];

  for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
    species[i] = pcell->species[i];

  double prate[TGCHEM_NUM_SPECIES][TGCHEM_NUM_CHEM_RATES];
  double nrate[TGCHEM_NUM_SPECIES][TGCHEM_NUM_CHEM_RATES];

  for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
    for(j = 0; j < TGCHEM_NUM_CHEM_RATES; j++)
      {
        prate[i][j] = pcell->prate[i][j];
        nrate[i][j] = pcell->nrate[i][j];
      }

  // Unpack var data

  double abh2  = var.abh2;
  double abhii = var.abhii;
  double gamma = var.gamma;
  double temp  = var.temp;

  // Append output variables

  pout->nh[count]    = nh;
  pout->gamma[count] = gamma;
  pout->temp[count]  = temp;

  for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
    pout->species[i * TGCHEM_MAX_NUM_OUT + count] = species[i];

  int k = 0;

  for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
    {
      for(j = 0; j < TGCHEM_NUM_CHEM_RATES; j++)
        {
          pout->rates[k * TGCHEM_MAX_NUM_OUT + count] = prate[i][j];

          k++;
        }

      for(j = 0; j < TGCHEM_NUM_CHEM_RATES; j++)
        {
          pout->rates[k * TGCHEM_MAX_NUM_OUT + count] = nrate[i][j];

          k++;
        }
    }

  // Update step data;

  pstep->count++;

  // Print to screen

  int flag = 0;

  if(TGCD.DebugFlag >= 0)
    {
      if(count <= TGCD.DebugFlag)
        flag = 1;
      else
        flag = 1;
    }
  else
    flag = 1;

  if(flag)
    {
      printf("Step %-4d: (%d|%d|%-3ld|%-3ld)   ", count, eq_flag_h2, eq_flag_hii, num_sub_steps, num_rate_calls);
      printf("NH: %-12g Temp: %-12g Abund: (%g|%g)\n", nh, temp, abh2, abhii);
    }
}

void tgchem_0d_write_outfile(step_data* pstep, out_data* pout)
{
  int i;
  char buf[MAX_STRING_LEN];
  FILE* file;

  // Unpack step data

  int count = pstep->count;

  // Write output

  sprintf(buf, "../../sator/data/tgchem.dat");

  if(!(file = fopen(buf, "w")))
    terminate("Could not open file!\n");

  fwrite(&count, sizeof(int), 1, file);

  int idummy = TGCHEM_NUM_SPECIES;

  fwrite(&idummy, sizeof(int), 1, file);

  idummy = TGCHEM_NUM_CHEM_RATES;

  fwrite(&idummy, sizeof(int), 1, file);

  double dummy;

  for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
    {
      dummy = NV_Ith_S(TGCD.SpeciesTol, i);

      fwrite(&dummy, sizeof(double), 1, file);
    }

  fwrite(pout->nh, sizeof(double), count, file);
  fwrite(pout->gamma, sizeof(double), count, file);
  fwrite(pout->temp, sizeof(double), count, file);

  for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
    fwrite(pout->species + i * TGCHEM_MAX_NUM_OUT, sizeof(double), count, file);

  for(i = 0; i < 2 * TGCHEM_NUM_SPECIES * TGCHEM_NUM_CHEM_RATES; i++)
    fwrite(pout->rates + i * TGCHEM_MAX_NUM_OUT, sizeof(double), count, file);

  fclose(file);
}
