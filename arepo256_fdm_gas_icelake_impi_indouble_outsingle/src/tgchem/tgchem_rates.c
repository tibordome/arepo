/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/tgchem/tgchem_rates.c
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

double tgchem_h2_line_fesc(double nh, double temp, double h2_column);
double tgchem_h2_cie_fesc(double nh);
double tgchem_hi_exc_fesc(double nh);
double tgchem_hm_fesc(double nh, double temp, double ljeans);
double tgchem_h2_fshield(double vth, double h2_column);
double tgchem_planck_opac(double rho, double temp);

void tgchem_rates(int mode, step_data* pstep, cell_data* pcell)
{
  int i, j;

  // Unpack step data

  double temp_cmb = pstep->temp_cmb;

  // Unpack cell data

  double nh        = pcell->nh;
  double nh2       = pcell->nh2;
  double nh3       = pcell->nh3;
  double vth       = pcell->vth;
  double ljeans    = pcell->ljeans;
  double ljeans2   = ljeans * ljeans;
  double h2_column = pcell->h2_column;

  // Compute var data

  var_data var;

  tgchem_var_data(pcell, &var);

  // Unpack var data

  double abh2  = var.abh2;
  double abhii = var.abhii;
  double abhi  = var.abhi;
  double abe   = var.abe;
  double temp  = var.temp;

  // Compute rates

  double prate[TGCHEM_NUM_SPECIES][TGCHEM_NUM_CHEM_RATES];
  double nrate[TGCHEM_NUM_SPECIES][TGCHEM_NUM_CHEM_RATES];

  for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
    for(j = 0; j < TGCHEM_NUM_CHEM_RATES; j++)
      prate[i][j] = nrate[i][j] = 0.;

  // *** H- *** //

  // H + e -> H- + ph
  prate[0][0] = tgchem_chem_coeff(0, temp) * abhi * abe * nh;

  // H + H- -> H2 + e
  nrate[0][0] = tgchem_chem_coeff(1, temp) * abhi * nh;

  double abhm = 0.;

  if(nrate[0][0])
    {
      abhm = prate[0][0] / nrate[0][0];

      abhm = tgchem_limit_abhm(abhm);
    }

  // *** H2 *** //

  // H + H- -> H2 + e
  prate[1][0] = tgchem_chem_coeff(1, temp) * abhi * abhm * nh;

  // 3H -> H + H2
  prate[1][1] = tgchem_chem_coeff(2, temp) * abhi * abhi * abhi * nh2;

  // 2H + H2 -> 2H2
  prate[1][2] = tgchem_chem_coeff(2, temp) / 8. * abhi * abhi * abh2 * nh2;

  // H + H2 -> 3H
  nrate[1][0] = tgchem_chem_coeff(3, temp) * abhi * abh2 * nh;

  // H2 + H2 -> 2H + H2
  nrate[1][1] = tgchem_chem_coeff(3, temp) / 8. * abh2 * abh2 * nh;

  // H2 + ph -> 2H
  if(!TGCD.ChemMode)
    nrate[1][2] = TGCD.LWDissRate * TGCD.ChemJ21 * abh2 * tgchem_h2_fshield(vth, h2_column);
  else
    nrate[1][2] = 0.;

  // *** HII *** //

  // H + e -> H+ + 2e
  prate[2][0] = tgchem_chem_coeff(4, temp) * abhi * abe * nh;

  // H+ + 2e -> H + e + ph
  double k1 = tgchem_chem_coeff(5, temp) * abhii * abe * nh;

  // H+ + 2e -> H + e + ph
  double k2 = tgchem_chem_coeff(6, temp) * abe * abe * abe * nh2;

  double coeff = 1. / (1. + pow(temp / 7e3, 2));

  // nrate[2][0] = pow(k1, coeff) * pow(k2, 1. - coeff);

  nrate[2][0] = k1;

  // *** Energy *** //

  // Chemical heating: H + H- -> H2 + e
  prate[3][0] = TGCHEM_CHI_H2 * tgchem_chem_coeff(1, temp) * abhi * abhm * nh2;

  // Chemical heating: 3H -> H + H2
  prate[3][1] = TGCHEM_CHI_H2 * tgchem_chem_coeff(2, temp) * abhi * abhi * abhi * nh3;

  // Chemical heating: 2H + H2 -> 2H2
  prate[3][2] = TGCHEM_CHI_H2 * tgchem_chem_coeff(2, temp) / 8. * abhi * abhi * abh2 * nh3;

  // Chemical cooling: H + H2 -> 3H
  nrate[3][0] = TGCHEM_CHI_H2 * tgchem_chem_coeff(3, temp) * abhi * abh2 * nh2;

  // Chemical cooling: H2 + H2 -> 2H + H2
  nrate[3][1] = TGCHEM_CHI_H2 * tgchem_chem_coeff(3, temp) / 8. * abh2 * abh2 * nh2;

  // Chemical cooling: H2 + ph -> 2H
  if(!TGCD.ChemMode)
    nrate[3][2] = TGCHEM_CHI_H2 * TGCD.LWDissRate * TGCD.ChemJ21 * abh2 * nh * tgchem_h2_fshield(vth, h2_column);
  else
    nrate[3][2] = 0.;

  // H2 ro-vibrational cooling
  double n0_rate  = tgchem_cool_coeff(0, temp) * abh2 * abhi * nh2;
  double lte_rate = tgchem_cool_coeff(1, temp) * abh2 * nh;

  if(n0_rate)
    nrate[3][3] = lte_rate / (1. + lte_rate / n0_rate) * tgchem_h2_line_fesc(nh, temp, h2_column);
  else
    nrate[3][3] = 0.;

  // H2 collision-induced emission
  nrate[3][4] = tgchem_cool_coeff(2, temp) * abh2 * nh2 * tgchem_h2_cie_fesc(nh);

  // HI electronic excitation cooling
  nrate[3][5] = tgchem_cool_coeff(3, temp) * abhi * abe * nh2 * tgchem_hi_exc_fesc(nh);

  // H- continuum cooling
  nrate[3][6] = TGCHEM_CHI_HM * tgchem_chem_coeff(0, temp) * abhi * abe * nh2 * tgchem_hm_fesc(nh, temp, ljeans);

  // H- continuum cooling (Becerra et al. 2017)

  // Optically thin cooling rate
  double hm_thin_low = tgchem_cool_coeff(4, temp) * abhi * abe * nh2;

  double hm_thin_upp = tgchem_cool_coeff(5, temp) * abhi * abe * nh2;

  // Planck mean opacity

  // Total opacity lower frequency: H- free-free
  double kappa_p_low = tgchem_opac_coeff(0, temp) * abhi * abe * nh2;

  // Total opacity upper frequency: H- bound-free + H- free-free
  double kappa_p_upp = tgchem_opac_coeff(1, temp) * abhm * nh + tgchem_opac_coeff(2, temp) * abhi * abe * nh2;

  // Rosseland mean opacity

  // Total opacity lower frequency: H- free-free
  double kappa_r_low = tgchem_opac_coeff(3, temp) * abhi * abe * nh2;

  // Rayleigh scattering

  // Total opacity upper frequency: Rayleigh scattering
  double kappa_r_upp = tgchem_opac_coeff(4, temp) * abhi * nh;

  // Optically thick cooling rate = Upper frequency + Lower frequency
  nrate[3][6] =
      hm_thin_upp / (1. + 3. * kappa_p_upp * kappa_r_upp * ljeans2) + hm_thin_low / (1. + 3. * kappa_p_low * kappa_r_low * ljeans2);

  // Inverse Compton cooling with the CMB (Peebles 1971)
  nrate[3][7] = 5.65e-36 * pow(temp_cmb, 4) * (temp - temp_cmb) * abe * nh;

  // Compute total rates

  double tot_rate[TGCHEM_NUM_SPECIES];

  for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
    {
      tot_rate[i] = 0.;

      for(j = 0; j < TGCHEM_NUM_CHEM_RATES; j++)
        tot_rate[i] += prate[i][j] - nrate[i][j];
    }

  // Update cell data

  pstep->num_rate_calls++;

  pcell->abhm = abhm;

  for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
    {
      for(j = 0; j < TGCHEM_NUM_CHEM_RATES; j++)
        {
          pcell->prate[i][j] = prate[i][j];
          pcell->nrate[i][j] = nrate[i][j];
        }

      pcell->tot_rate[i] = tot_rate[i];
    }
}

void tgchem_rates_aux(cell_data* pcell, double dt_step)
{
  // Unpack cell data

  int i, j;
  double nh = pcell->nh;
  double species[TGCHEM_NUM_SPECIES];

  for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
    species[i] = pcell->species[i];

  double prate[TGCHEM_NUM_SPECIES][TGCHEM_NUM_CHEM_RATES];
  double nrate[TGCHEM_NUM_SPECIES][TGCHEM_NUM_CHEM_RATES];
  double tot_rate[TGCHEM_NUM_SPECIES];

  for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
    {
      for(j = 0; j < TGCHEM_NUM_CHEM_RATES; j++)
        {
          prate[i][j] = pcell->prate[i][j];
          nrate[i][j] = pcell->nrate[i][j];
        }

      tot_rate[i] = pcell->tot_rate[i];
    }

  // Compute timescales

  double dt_rate[TGCHEM_NUM_SPECIES];

  for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
    dt_rate[i] = species[i] / fabs(tot_rate[i]);

  double tot_chem_rate = 0.;

  for(i = 0; i < 3; i++)
    tot_chem_rate += prate[TGCHEM_NUM_ABUNDANCES][i];

  for(i = 0; i < 3; i++)
    tot_chem_rate += -nrate[TGCHEM_NUM_ABUNDANCES][i];

  double dt_chem = species[TGCHEM_NUM_ABUNDANCES] / fabs(tot_chem_rate);

  double tot_rad_rate = 0.;

  for(i = 3; i < TGCHEM_NUM_CHEM_RATES; i++)
    tot_rad_rate += -nrate[TGCHEM_NUM_ABUNDANCES][i];

  double dt_rad = species[TGCHEM_NUM_ABUNDANCES] / fabs(tot_rad_rate);

  // Compute eq flags

  int eq_flag[TGCHEM_NUM_SPECIES];
  double max_prate, max_nrate, dt_prate, dt_nrate;

  for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
    {
      if(i == 0)
        eq_flag[i] = 1;
      else if(i == TGCHEM_NUM_ABUNDANCES)
        eq_flag[i] = 0;
      else
        {
          if(nh > 1e10)
            {
              max_prate = max_nrate = 0.;

              for(j = 0; j < TGCHEM_NUM_CHEM_RATES; j++)
                {
                  max_prate = fmax(max_prate, prate[i][j]);
                  max_nrate = fmax(max_nrate, nrate[i][j]);
                }

              dt_prate = species[i] / max_prate;
              dt_nrate = species[i] / max_nrate;

              if(dt_prate < TGCHEM_TOL_TRANS * dt_step && dt_nrate < TGCHEM_TOL_TRANS * dt_step)
                eq_flag[i] = 1;
              else
                eq_flag[i] = 0;
            }
          else
            eq_flag[i] = 0;
        }
    }

  // Update cell data

  for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
    {
      pcell->dt_rate[i] = dt_rate[i];
      pcell->eq_flag[i] = eq_flag[i];
    }

  pcell->tot_chem_rate = tot_chem_rate;
  pcell->tot_rad_rate  = tot_rad_rate;

  pcell->dt_chem = dt_chem;
  pcell->dt_rad  = dt_rad;
}

double tgchem_h2_line_fesc(double nh, double temp, double h2_column)
{
  double fesc;

  if(TGCD.ChemH2Mode == 0)
    {
      // Fitting function (based on Ripamonti & Abel 2004)

      double x = nh / TGCD.H2OptThickNHThresh;

      if(x >= 1.)
        fesc = TGCD.H2OptThickConst * x / (pow(x, TGCD.H2OptThickConst) + TGCD.H2OptThickConst - 1.);
      else
        fesc = 1.;
    }
  else if(TGCD.ChemH2Mode == 1)
    {
      // Sobolev method (Yoshida et al. 2006, Clark et al. 2011)

      int temp_idx = tgchem_compute_temp_idx(temp);

      double dtemp = tgchem_compute_dtemp(temp, temp_idx);

      int h2_column_idx = (int)(log10(h2_column / TGCHEM_H2_COLUMN_MIN) / TGCHEM_H2_LOG_DCOLUMN);
      h2_column_idx     = imin(h2_column_idx, TGCHEM_H2_NUM_COLUMN - 1);

      double dh2_column = h2_column - TGCD.H2ColumnTable[h2_column_idx];

      int idx1 = temp_idx * TGCHEM_H2_NUM_COLUMN + h2_column_idx;
      int idx2 = fmin((temp_idx + 1), TGCHEM_NUM_TEMP - 1) * TGCHEM_H2_NUM_COLUMN + h2_column_idx;

      double emiss1 = TGCD.H2SobEmissTable[idx1] + dh2_column * TGCD.H2DSobEmissTable[idx1];
      double emiss2 = TGCD.H2SobEmissTable[idx2] + dh2_column * TGCD.H2DSobEmissTable[idx2];

      double sob_emiss = emiss1 + dtemp * (emiss2 - emiss1);

      fesc = fmin(sob_emiss / tgchem_cool_coeff(1, temp), 1.);
    }
  else if(TGCD.ChemH2Mode == 0)
    {
      // Optically thin

      fesc = 1.;
    }
  else
    terminate("Chosen value for TGCD.ChemMode is not supported!");

  return fesc;
}

double tgchem_h2_cie_fesc(double nh)
{
  // Fitting function (Clark et al. 2011)

  double x = nh / TGCD.CIEOptThickNHThresh;

  double tau = fmax(pow(x, 2.8), 1e-10);

  double fesc = (1. - exp(-tau)) / tau;

  return fesc;
}

double tgchem_hi_exc_fesc(double nh)
{
  // We need something better than this!

  double x = nh / 1e7;

  double fesc = exp(-x);

  return fesc;
}

double tgchem_hm_fesc(double nh, double temp, double ljeans)
{
  // Assumes Planck mean opacity and uniform density cloud

  double rho = PROTONMASS * nh / HYDROGEN_MASSFRAC;

  double kappa = tgchem_planck_opac(rho, temp);

  double tau = fmax(kappa * rho * ljeans, 1e-10);

  double fesc = exp(-tau);

  return fesc;
}

double tgchem_h2_fshield(double vth, double h2_column)
{
  double b5 = vth / 1e5;

  double x = h2_column / 4e14;

  double fshield = 0.965 / pow(1. + x / b5, 1.1) + 0.035 / sqrt(1. + x) * exp(-8.5e-4 * sqrt(1. + x));

  return fshield;
}

double tgchem_planck_opac(double rho, double temp)
{
  int rho_idx, temp_idx, idx1, idx2;
  double log_rho, log_temp, drho, dtemp;
  double opac1, opac2, opac3, opac4;

  log_rho = log10(rho);
  rho_idx = log_rho + 16;

  if(rho_idx < 0)
    {
      rho_idx = 0;
      drho    = 0.;
    }
  else if(rho_idx > TGCHEM_PLANCK_OPAC_NUM_RHO - 1)
    {
      rho_idx = TGCHEM_PLANCK_OPAC_NUM_RHO - 1;
      drho    = 0.;
    }
  else
    drho = log_rho - (rho_idx - 16);

  log_temp = log10(temp);
  temp_idx = (log_temp - 1.8) / 0.1;

  if(temp_idx < 0)
    {
      temp_idx = 0;
      dtemp    = 0.;
    }
  else if(temp_idx > TGCHEM_PLANCK_OPAC_NUM_TEMP - 1)
    {
      temp_idx = TGCHEM_PLANCK_OPAC_NUM_TEMP - 1;
      dtemp    = 0.;
    }
  else
    dtemp = (log_temp - (1.8 + 0.1 * temp_idx)) / 0.1;

  idx1 = TGCHEM_PLANCK_OPAC_NUM_RHO * temp_idx + rho_idx;

  if(rho_idx < TGCHEM_PLANCK_OPAC_NUM_RHO - 1)
    idx2 = idx1 + 1;
  else
    idx2 = idx1;

  opac1 = TGCD.PlanckOpacTable[idx1];
  opac2 = TGCD.PlanckOpacTable[idx2];

  opac3 = opac1 + drho * (opac2 - opac1);

  if(temp_idx < TGCHEM_PLANCK_OPAC_NUM_TEMP - 1)
    {
      idx1 = TGCHEM_PLANCK_OPAC_NUM_RHO * (temp_idx + 1) + rho_idx;

      if(rho_idx < TGCHEM_PLANCK_OPAC_NUM_RHO - 1)
        idx2 = idx1 + 1;
      else
        idx2 = idx1;

      opac1 = TGCD.PlanckOpacTable[idx1];
      opac2 = TGCD.PlanckOpacTable[idx2];

      opac4 = opac1 + drho * (opac2 - opac1);
    }
  else
    opac4 = opac3;

  double log_opac = opac3 + dtemp * (opac4 - opac3);

  double opac = pow(10., log_opac);

  return opac;
}
