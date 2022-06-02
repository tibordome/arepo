#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "proto.h"

double photoheating(double xi, double Gamma, double epsilon, double nH)
{
  return -xi * Gamma * epsilon / nH; /* Note: Heating, so returns negative value */
}

double cosmic_ray_heating(double xi, double zeta_i, double nH)
{
  return -3.2e-11 * xi * zeta_i / nH; /* Note: Heating, so returns negative value */
}

double compton_cooling(double T, double Tcmb, double xe, double nH)
{
  return 1.017e-37 * pow(Tcmb, 4) * (T - Tcmb) * xe / nH; /* Lambda/nH^2 */
}

/* Using Wolfire et al. (1995) */
double photoelectric_grain_heating(double T, double nH, double ne, double n_total, double Z, double dust_G, double extinction)
{
  /* From Wolfire et al. (2003) */
  if((dust_G * exp(-extinction * G0_GAMMA)) == 0 || ne == 0)
    return 0.0;
  else
    {
      double phi_pah = 0.5;
      double psi     = dust_G * exp(-extinction * G0_GAMMA) * pow(T, 0.5) / (ne * phi_pah);
      double epsilon = (4.9e-2 / (1.0 + pow(psi / 1925.0, 0.73))) + (3.7e-2 * pow((T / 1.0e4), 0.7) / (1.0 + (psi / 5.0e3)));
      return -1.3e-24 * epsilon * dust_G * exp(-extinction * G0_GAMMA) * Z / nH;
    }
}

double gas_grain_transfer(double T, double Tgr, double Z)
{
  /* Note: The metallicity Z here is used as a
   * proxy for the dust ratio */
  return 3.8e-33 * pow(T, 0.5) * (T - Tgr) * (1.0 - 0.8 * exp(-75.0 / T)) * Z;
}

double grain_surface_recombination(double T, double Z, double ne, double nH, double dust_G, double extinction)
{
  /* Note: The metallicity Z here is used as a
   * proxy for the dust ratio */
  double psi;
  if(ne != 0.0)
    {
      psi = 2 * dust_G * exp(-extinction * G0_GAMMA) * pow(T, 0.5) / ne;
      return 2.33e-30 * pow(T, 0.94) * pow(psi, (0.74 / pow(T, 0.068))) * Z * ne / nH;
    }
  else
    return 0.0;
}

/* Cooling and heating processes involving molecules */

/* The following returns the collisional de-excitation rates
 * of H2, as given by Draine et al. (1983). */
double H2_collis_deexc_rates(double T, double J, int mode)
{
  double x;
  if(mode == 1) /* H collisions; rot levels */
    {
      if(J == 2.0)
        x = 7.0438e-14 / (BOLTZMANNCGS * T);
      else /* NB: We assume J is either 2 or 3 */
        x = 1.1669e-13 / (BOLTZMANNCGS * T);
      return 4.6e-12 * (2 * J - 1.0) * pow(T, 0.5) * pow(1.0 + x, 0.5) * exp(-(5.01 * x + 0.1187 * (4.0 * J - 2.0)));
    }
  else if(mode == 2) /* H2 collisions; rot levels */
    {
      if(J == 2.0)
        return 4.0e-14 * pow(T, 0.75);
      else if(J == 3.0)
        return 5.0e-15 * T;
      else
        return 0.0;
    }
  else if(mode == 3) /* H collisions; vib levels */
    {
      if(J == 1.0) /* Note: J is actually v here */
        return 1.0e-12 * pow(T, 0.5) * exp(-1000.0 / T);
      else if(J == 2.0)
        return 1.6e-12 * pow(T, 0.5) * exp(-pow(400.0 / T, 2));
      else
        return 0.0;
    }
  else if(mode == 4) /* H2 collisions; vib levels */
    return 6.6e-14 * T * exp(-79.99 / pow(T, 1.0 / 3.0)) * (1.0 - exp(-5860.0 / T));
  else
    return 0.0;
}

/* Using the fits from table 8 of Glover & Abel (2008) */
double H2_rovibrational_cooling(double T, double xH2, double xH, double xHII, double xHe, double xe, double nH_tot)
{
  double LowDens_cool_rate, LTE_cool_rate;
  double T3 = T / 1000.0;
  int T_index;
  double dT;
  if(T <= 100.0)
    LowDens_cool_rate = pow(10.0, -16.818342 + 37.383713 * log10(T3) + 58.145166 * pow(log10(T3), 2) + 48.656103 * pow(log10(T3), 3) +
                                      20.159831 * pow(log10(T3), 4) + 3.847961 * pow(log10(T3), 5)) *
                        xH2 * xH;
  else if(T <= 1000.0)
    LowDens_cool_rate = pow(10.0, -24.311209 + 3.5692468 * log10(T3) - 11.33286 * pow(log10(T3), 2) - 27.850082 * pow(log10(T3), 3) -
                                      21.328264 * pow(log10(T3), 4) - 4.2519023 * pow(log10(T3), 5)) *
                        xH2 * xH;
  else
    LowDens_cool_rate = pow(10.0, -24.311209 + 4.6450521 * log10(T3) - 3.7209846 * pow(log10(T3), 2) + 5.9369081 * pow(log10(T3), 3) -
                                      5.5108047 * pow(log10(T3), 4) + 1.5538288 * pow(log10(T3), 5)) *
                        xH2 * xH;
  LowDens_cool_rate += pow(10.0, -23.962112 + 2.0943374 * log10(T3) - 0.77151436 * pow(log10(T3), 2) + 0.43693353 * pow(log10(T3), 3) -
                                     0.14913216 * pow(log10(T3), 4) - 0.033638326 * pow(log10(T3), 5)) *
                       xH2 * xH2;
  LowDens_cool_rate += pow(10.0, -23.689237 + 2.1892372 * log10(T3) - 0.81520438 * pow(log10(T3), 2) + 0.29036281 * pow(log10(T3), 3) -
                                     0.16596184 * pow(log10(T3), 4) + 0.19191375 * pow(log10(T3), 5)) *
                       xH2 * xHe;
  LowDens_cool_rate += pow(10.0, -21.716699 + 1.3865783 * log10(T3) - 0.37915285 * pow(log10(T3), 2) + 0.11453688 * pow(log10(T3), 3) -
                                     0.23214154 * pow(log10(T3), 4) + 0.058538864 * pow(log10(T3), 5)) *
                       xH2 * xHII;
  if(T <= 200.0)
    LowDens_cool_rate += pow(10.0, -34.286155 - 48.537163 * log10(T3) - 77.121176 * pow(log10(T3), 2) - 51.352459 * pow(log10(T3), 3) -
                                       15.16916 * pow(log10(T3), 4) - 0.98120322 * pow(log10(T3), 5)) *
                         xH2 * xe;
  else
    LowDens_cool_rate += pow(10.0, -22.190316 + 1.5728955 * log10(T3) - 0.213351 * pow(log10(T3), 2) + 0.96149759 * pow(log10(T3), 3) -
                                       0.91023195 * pow(log10(T3), 4) + 0.13749749 * pow(log10(T3), 5)) *
                         xH2 * xe;

  /* Obtain the LTE rate from the molecular
   * cooling tables. */
  get_index_1d_mydbl(chimesRateTables.mol_cooling_table_H2_lte_temperatures, chimesRateTables.mol_cooling_table_dimensions[10],
                     log10(T), &T_index, &dT);
  LTE_cool_rate = pow(10.0, interpol_1d_mydbl(chimesRateTables.mol_cooling_table_H2_lte, T_index, dT)) *
                  (xH2 / nH_tot); /* lambda / nH_tot ^ 2 (in units erg cm^3 s^-1) */

  if(LowDens_cool_rate != 0)
    return LTE_cool_rate / (1.0 + (LTE_cool_rate / LowDens_cool_rate)); /* Lambda/nH^2, in erg cm^3 /s */
  else
    return 0.0;
}

double H2_crit_density(double T, double xH, double xH2)
{
  /* Returns the critical density at which collisional de-excitation
   * of H2 occurs at the same rate as radiative de-excitation. See
   * GJ07 & references therein*/
  double n_cr_H, n_cr_H2;

  if(xH + xH2 == 0)
    return 0.0;
  else
    {
      n_cr_H  = pow(10, (3.0 - 0.416 * log10(T / 1e4) - 0.327 * pow(log10(T / 1e4), 2)));
      n_cr_H2 = pow(10, (4.845 - 1.3 * log10(T / 1e4) + 1.62 * pow(log10(T / 1e4), 2)));
      return (xH + xH2) / ((xH / n_cr_H) + (xH2 / n_cr_H2));
    }
}

/* The following functions calculate the cooling rates
 * from the molecular CO, H2O and OH, using the fitting
 * coefficients in the molecular cooling tables. */

double CO_rotational_cooling(double T, double N, double xm, double xH2, double xHI, double xe, double nH,
                             struct globalVariables *myGlobalVars)
{
  int i, j;
  double dT, dN;
  double Lm, L0, Llte, nhalf, a, neff;

  /* Obtain the indices for T and N. Note that T is irregularly
   * spaced - the Neufeld & Kaufman (1993) suggest you should
   * be using a 2D cubic spline to interpolate on the grid of
   * (T, N). If you want to continue using linear interpolation
   * in log space (which is probably computationally faster), it
   * would be better if you calculate the 2D cubic splines and
   * then retabulate the coefficients on a much finer, and evenly
   * spaced, (T, N) grid.
   */
  get_index_1d_irregular(chimesRateTables.mol_cooling_table_CO_rot_T, chimesRateTables.mol_cooling_table_dimensions[0], log10(T), &i,
                         &dT);
  get_index_1d_mydbl(chimesRateTables.mol_cooling_table_CO_rot_N, chimesRateTables.mol_cooling_table_dimensions[1], log10(N), &j, &dN);

  /* Interpolate the cooling coefficients */
  L0    = pow(10.0, -interpol_1d_mydbl(chimesRateTables.mol_cooling_table_CO_rot_L0, i, dT));
  Llte  = pow(10.0, -interpol_2d_mydbl(chimesRateTables.mol_cooling_table_CO_rot_Llte, i, j, dT, dN));
  nhalf = pow(10.0, interpol_2d_mydbl(chimesRateTables.mol_cooling_table_CO_rot_nhalf, i, j, dT, dN));
  a     = pow(10.0, interpol_2d_mydbl(chimesRateTables.mol_cooling_table_CO_rot_a, i, j, dT, dN));

  neff = nH * (xH2 + 9.857 * pow((T / 1.0e3), 0.25) * xHI + 680.13 * pow(T, -0.25) * xe);

  Lm = 1.0 / ((1.0 / L0) + (neff / Llte) + (1.0 / L0) * pow((neff / nhalf), a) * (1.0 - (nhalf * L0 / Llte)));

  return Lm * xm * xH2; /* Lambda / nH^2, in erg cm^3 s^-1 */
}

double CO_vibrational_cooling(double T, double N, double xm, double xH2, double xHI, double xe, double nH,
                              struct globalVariables *myGlobalVars)
{
  int i, j;
  double dT, dN;
  double Lm, L0, Llte, neff;

  /* Obtain the indices for T and N. */
  get_index_1d_irregular(chimesRateTables.mol_cooling_table_CO_vib_T, chimesRateTables.mol_cooling_table_dimensions[2], log10(T), &i,
                         &dT);
  get_index_1d_mydbl(chimesRateTables.mol_cooling_table_CO_vib_N, chimesRateTables.mol_cooling_table_dimensions[3], log10(N), &j, &dN);

  /* Interpolate the cooling coefficients */
  L0   = 1.83e-26 * T * exp(-68.0 / pow(T, 1.0 / 3.0)) * exp(-3080.0 / T);
  Llte = pow(10.0, -interpol_2d_mydbl(chimesRateTables.mol_cooling_table_CO_vib_Llte, i, j, dT, dN)) * exp(-3080.0 / T);

  neff = nH * (xH2 + 50.0 * xHI + 9035.09 * exp(68.0 / pow(T, 1.0 / 3.0)) * pow(T / 300.0, 0.938) * xe);

  Lm = 1.0 / ((1.0 / L0) + (neff / Llte));

  return Lm * xm * xH2; /* Lambda / nH^2, in erg cm^3 s^-1 */
}

double H2O_rotational_cooling(double T, double N, double xm, double xH2, double xHI, double xe, double nH,
                              struct globalVariables *myGlobalVars)
{
  int i, j;
  double dT, dN;
  double Lm, L0, Llte, nhalf, a, neff;

  if(T >= 100.0)
    {
      /* Obtain the indices for T and N. */
      get_index_1d_irregular(chimesRateTables.mol_cooling_table_H2O_rot_hiT_T, chimesRateTables.mol_cooling_table_dimensions[4],
                             log10(T), &i, &dT);
      get_index_1d_mydbl(chimesRateTables.mol_cooling_table_H2O_rot_hiT_N, chimesRateTables.mol_cooling_table_dimensions[5], log10(N),
                         &j, &dN);

      /* Interpolate the cooling coefficients */
      L0    = pow(10.0, -interpol_1d_mydbl(chimesRateTables.mol_cooling_table_H2O_rot_hiT_L0, i, dT));
      Llte  = pow(10.0, -interpol_2d_mydbl(chimesRateTables.mol_cooling_table_H2O_rot_hiT_Llte, i, j, dT, dN));
      nhalf = pow(10.0, interpol_2d_mydbl(chimesRateTables.mol_cooling_table_H2O_rot_hiT_nhalf, i, j, dT, dN));
      a     = pow(10.0, interpol_2d_mydbl(chimesRateTables.mol_cooling_table_H2O_rot_hiT_a, i, j, dT, dN));

      neff = nH * (xH2 + 10.0 * xHI +
                   pow(10.0, (-8.02 + (15.749 / pow(T, 1.0 / 6.0)) - (47.137 / pow(T, 1.0 / 3.0)) + (76.648 / pow(T, 0.5)) -
                              (60.191 / pow(T, 2.0 / 3.0)))) *
                       xe / (7.4e-12 * pow(T, 0.5)));

      Lm = 1.0 / ((1.0 / L0) + (neff / Llte) + (1.0 / L0) * pow((neff / nhalf), a) * (1.0 - (nhalf * L0 / Llte)));
    }
  else
    {
      /* Obtain the indices for T and N. */
      get_index_1d_irregular(chimesRateTables.mol_cooling_table_H2O_rot_lowT_T, chimesRateTables.mol_cooling_table_dimensions[6],
                             log10(T), &i, &dT);
      get_index_1d_mydbl(chimesRateTables.mol_cooling_table_H2O_rot_lowT_N, chimesRateTables.mol_cooling_table_dimensions[7], log10(N),
                         &j, &dN);

      /* Assume an ortho:para ratio of 3:1 */

      /* Interpolate the cooling coefficients */
      L0    = pow(10.0, -interpol_1d_mydbl(chimesRateTables.mol_cooling_table_H2Oortho_rot_lowT_L0, i, dT));
      Llte  = pow(10.0, -interpol_2d_mydbl(chimesRateTables.mol_cooling_table_H2Oortho_rot_lowT_Llte, i, j, dT, dN));
      nhalf = pow(10.0, interpol_2d_mydbl(chimesRateTables.mol_cooling_table_H2Oortho_rot_lowT_nhalf, i, j, dT, dN));
      a     = pow(10.0, interpol_2d_mydbl(chimesRateTables.mol_cooling_table_H2Oortho_rot_lowT_a, i, j, dT, dN));

      neff = nH * (xH2 + 10.0 * xHI +
                   pow(10.0, (-8.02 + (15.749 / pow(T, 1.0 / 6.0)) - (47.137 / pow(T, 1.0 / 3.0)) + (76.648 / pow(T, 0.5)) -
                              (60.191 / pow(T, 2.0 / 3.0)))) *
                       xe / (7.4e-12 * pow(T, 0.5)));

      Lm = 0.75 / ((1.0 / L0) + (neff / Llte) + (1.0 / L0) * pow((neff / nhalf), a) * (1.0 - (nhalf * L0 / Llte)));

      /* Interpolate the cooling coefficients */
      L0    = pow(10.0, -interpol_1d_mydbl(chimesRateTables.mol_cooling_table_H2Opara_rot_lowT_L0, i, dT));
      Llte  = pow(10.0, -interpol_2d_mydbl(chimesRateTables.mol_cooling_table_H2Opara_rot_lowT_Llte, i, j, dT, dN));
      nhalf = pow(10.0, interpol_2d_mydbl(chimesRateTables.mol_cooling_table_H2Opara_rot_lowT_nhalf, i, j, dT, dN));
      a     = pow(10.0, interpol_2d_mydbl(chimesRateTables.mol_cooling_table_H2Opara_rot_lowT_a, i, j, dT, dN));

      Lm += 0.25 / ((1.0 / L0) + (neff / Llte) + (1.0 / L0) * pow((neff / nhalf), a) * (1.0 - (nhalf * L0 / Llte)));
    }

  return Lm * xm * xH2; /* Lambda / nH^2, in erg cm^3 s^-1 */
}

double H2O_vibrational_cooling(double T, double N, double xm, double xH2, double xHI, double xe, double nH,
                               struct globalVariables *myGlobalVars)
{
  int i, j;
  double dT, dN;
  double Lm, L0, Llte, neff;

  /* Obtain the indices for T and N. */
  get_index_1d_irregular(chimesRateTables.mol_cooling_table_H2O_vib_T, chimesRateTables.mol_cooling_table_dimensions[8], log10(T), &i,
                         &dT);
  get_index_1d_mydbl(chimesRateTables.mol_cooling_table_H2O_vib_N, chimesRateTables.mol_cooling_table_dimensions[9], log10(N), &j,
                     &dN);

  /* Interpolate the cooling coefficients */
  L0   = 1.03e-26 * T * exp(-47.5 / pow(T, 1.0 / 3.0)) * exp(-2325.0 / T);
  Llte = pow(10.0, -interpol_2d_mydbl(chimesRateTables.mol_cooling_table_H2O_vib_Llte, i, j, dT, dN)) * exp(-2325.0 / T);

  neff = nH * (xH2 + 10.0 * xHI + 4.0625e8 * exp(47.5 / pow(T, 1.0 / 3.0)) * pow(T, -0.5) * xe);

  Lm = 1.0 / ((1.0 / L0) + (neff / Llte));

  return Lm * xm * xH2; /* Lambda / nH^2, in erg cm^3 s^-1 */
}

double OH_rotational_cooling(double T, double N, double dv, double xm, double nH, double n, double tau_d)
{
  /* NOTE: the velocity dispersion dv should be in cgs */
  double tau_T, N_tau, c_tau, n_cr, ym;

  N_tau = 1.485e11 * dv;
  tau_T = 4.0 * N / (10.0 * (T / 27.0) * 6.8e-4 * N_tau);
  c_tau = tau_T * pow((2 * PI * log(2.13 + pow(tau_T / exp(1.0), 2.0))), 0.5) /
          ((exp(-tau_d) / (1.0 + pow(tau_d, 2.0))) +
           2.0 * tau_d * pow(log(1 + (tau_T / exp(1.0))), 0.5) * pow(log(tau_T / (tau_d * exp(1.0))), 0.5));
  if(isnan(c_tau) != 0)
    c_tau = 0.0;
  n_cr = 1.5e10 * pow(T / 1.0e3, 0.5);
  ym   = log(1.0 + (c_tau / (1.0 + 10.0 * (n_cr / n))));
  return (n / nH) * xm * (2.0 * pow(BOLTZMANNCGS * T, 2.0) * 2.3e-2 / (n * 27.0 * BOLTZMANNCGS)) *
         ((2.0 + ym + 0.6 * pow(ym, 2.0)) / (1.0 + c_tau + (n_cr / n) + 1.5 * pow(n_cr / n, 0.5))); /* Lambda / nH^2 (erg cm^3 s^-1) */
}

double calculate_mean_molecular_weight(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars)
{
  double denominator = 0.0;
  int i;

  for(i = 0; i < myGlobalVars->totalNumberOfSpecies; i++)
    denominator += myGasVars->abundances[i];

  return (1.0 + myGasVars->element_abundances[0] * 4.0 + myGasVars->element_abundances[1] * 12.0 +
          myGasVars->element_abundances[2] * 14.0 + myGasVars->element_abundances[3] * 16.0 + myGasVars->element_abundances[4] * 20.0 +
          myGasVars->element_abundances[5] * 24.0 + myGasVars->element_abundances[6] * 28.0 + myGasVars->element_abundances[7] * 32.0 +
          myGasVars->element_abundances[8] * 40.0 + myGasVars->element_abundances[9] * 56.0) /
         denominator;
}

double calculate_total_number_density(double *my_abundances, double nH, struct globalVariables *myGlobalVars)
{
  double result = 0.0;
  int i;

  for(i = 0; i < myGlobalVars->totalNumberOfSpecies; i++)
    result += my_abundances[i] * nH;

  return result;
}

double calculate_total_cooling_rate(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars, double HI_column_density,
                                    double HeI_column_density, double HeII_column_density, double H2_column_density,
                                    double CO_column_density, double H2O_column_density, double OH_column_density, double extinction,
                                    struct All_rate_variables_structure *this_all_rates)
{
  int ns, i, T_index, j, k, l, m, p;
  double dT, N_eff, n_cr, n, dT_2, dnHI, dne, dnHII;
  int NHI_index, NHeI_index, NH_eff_index, NHe_eff_index;
  double dNHI, dNHeI, dNH_eff, dNHe_eff, epsilon;
  int chianti_T_index, chianti_ne_index;
  double chianti_dT, chianti_dne;
  double E1, E2, E3, E4, E5, E6;
  double S1, S2, S3, photoion_rate, shieldFactor, dust_G;
  double total_cooling = 0.0;

  get_index_1d_mydbl(chimesRateTables.NonEqIon->Temperatures, chimesRateTables.NonEqIon->N_Temperatures,
                     log10(min(myGasVars->temperature, MAX_TEMPERATURE_FOR_RATES)), &T_index, &dT);
  get_index_1d_mydbl(chimesRateTables.NonEqIon->shieldingColumnDensities, chimesRateTables.NonEqIon->shieldingColumnDimensions[0],
                     log10(max(HI_column_density, 1.0e-50)), &NHI_index, &dNHI);
  get_index_1d_mydbl(chimesRateTables.NonEqIon->shieldingColumnDensities, chimesRateTables.NonEqIon->shieldingColumnDimensions[0],
                     log10(max(HI_column_density + 3.0 * H2_column_density, 1.0e-50)), &NH_eff_index, &dNH_eff);
  get_index_1d_mydbl(chimesRateTables.NonEqIon->shieldingColumnDensities, chimesRateTables.NonEqIon->shieldingColumnDimensions[0],
                     log10(max(HeI_column_density, 1.0e-50)), &NHeI_index, &dNHeI);
  get_index_1d_mydbl(chimesRateTables.NonEqIon->shieldingColumnDensities, chimesRateTables.NonEqIon->shieldingColumnDimensions[0],
                     log10(max(HeI_column_density + 0.75 * HeII_column_density, 1.0e-50)), &NHe_eff_index, &dNHe_eff);

  if(log10(myGasVars->temperature) < chimesRateTables.nei_cooling_temperature[chimesRateTables.nei_cooling_table_dimensions[0] - 1])
    {
      get_index_1d_mydbl(chimesRateTables.nei_cooling_temperature, chimesRateTables.nei_cooling_table_dimensions[0],
                         log10(myGasVars->temperature), &m, &dT_2);
      get_index_1d_mydbl(chimesRateTables.nei_cooling_HIAbundance, chimesRateTables.nei_cooling_table_dimensions[1],
                         log10(max(myGasVars->abundances[myGlobalVars->speciesIndices[HI]] * myGasVars->nH_tot, 1.0e-100)), &j, &dnHI);
      get_index_1d_mydbl(chimesRateTables.nei_cooling_ElectronAbundance, chimesRateTables.nei_cooling_table_dimensions[2],
                         log10(max(myGasVars->abundances[myGlobalVars->speciesIndices[elec]] * myGasVars->nH_tot, 1.0e-100)), &k,
                         &dne);
      get_index_1d_mydbl(chimesRateTables.nei_cooling_HIIAbundance, chimesRateTables.nei_cooling_table_dimensions[3],
                         log10(max(myGasVars->abundances[myGlobalVars->speciesIndices[HII]] * myGasVars->nH_tot, 1.0e-100)), &p,
                         &dnHII);
    }
  if(log10(myGasVars->temperature) <
     chimesRateTables.chianti_cooling_temperature[chimesRateTables.chianti_cooling_table_dimensions[0] - 1])
    {
      get_index_1d_mydbl(chimesRateTables.chianti_cooling_temperature, chimesRateTables.chianti_cooling_table_dimensions[0],
                         log10(myGasVars->temperature), &chianti_T_index, &chianti_dT);
      get_index_1d_mydbl(chimesRateTables.chianti_cooling_ElectronDensity, chimesRateTables.chianti_cooling_table_dimensions[1],
                         log10(max(myGasVars->abundances[myGlobalVars->speciesIndices[elec]] * myGasVars->nH_tot, 1.0e-100)),
                         &chianti_ne_index, &chianti_dne);
    }

  for(ns = 0; ns < chimesRateTables.NonEqIon->N_Elements; ns++)
    {
      for(i = 0; i < chimesRateTables.NonEqIon->N_Ions[ns]; i++)
        {
          if(log10(myGasVars->temperature) >=
                 chimesRateTables.nei_cooling_temperature[chimesRateTables.nei_cooling_table_dimensions[0] - 1] ||
             ((chimesRateTables.NonEqIon->IonIndexBegin[ns] + i) != CI && (chimesRateTables.NonEqIon->IonIndexBegin[ns] + i) != OI))
            {
              if(log10(myGasVars->temperature) >=
                     chimesRateTables.chianti_cooling_temperature[chimesRateTables.chianti_cooling_table_dimensions[0] - 1] ||
                 ((chimesRateTables.NonEqIon->IonIndexBegin[ns] + i) != CII &&
                  (chimesRateTables.NonEqIon->IonIndexBegin[ns] + i) != NII &&
                  (chimesRateTables.NonEqIon->IonIndexBegin[ns] + i) != SiII &&
                  (chimesRateTables.NonEqIon->IonIndexBegin[ns] + i) != FeII))
                total_cooling +=
                    pow(10.0, interpol_1d_mydbl(chimesRateTables.NonEqIon->NonEqRates[ns].cool[i], T_index, dT)) *
                    myGasVars->abundances[myGlobalVars->speciesIndices[chimesRateTables.NonEqIon->IonIndexBegin[ns] + i]] *
                    myGasVars->abundances[myGlobalVars->speciesIndices[elec]];
            }
          if(i != (chimesRateTables.NonEqIon->N_Ions[ns] - 1)) /* The final ion cannot be photoionised further */
            {
              for(l = 0; l < myGlobalVars->N_spectra; l++)
                {
                  /* Calculate photoionisation rate for given UV spectrum. */
                  if(chimesRateTables.NonEqIon->NonEqRates[ns].E_thresh[i] >= 13.6)
                    {
                      if(HI_column_density == 0.0 && H2_column_density == 0.0 && HeI_column_density == 0.0 &&
                         HeII_column_density == 0.0)
                        shieldFactor = 1.0;
                      else
                        {
                          if(chimesRateTables.NonEqIon->NonEqRates[ns].E_thresh[i] < 15.4)
                            {
                              S1 = pow(10.0, interpol_1d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor1D[l][0][i],
                                                                NHI_index, dNHI));
                              S2 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][0][i],
                                                                NH_eff_index, NHeI_index, dNH_eff, dNHeI));
                              S3 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][1][i],
                                                                NH_eff_index, NHe_eff_index, dNH_eff, dNHe_eff));
                            }
                          else if(chimesRateTables.NonEqIon->NonEqRates[ns].E_thresh[i] < 54.42)
                            {
                              S1 = 0.0;
                              S2 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][0][i],
                                                                NH_eff_index, NHeI_index, dNH_eff, dNHeI));
                              S3 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][1][i],
                                                                NH_eff_index, NHe_eff_index, dNH_eff, dNHe_eff));
                            }
                          else
                            {
                              S1 = 0.0;
                              S2 = 0.0;
                              S3 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][1][i],
                                                                NH_eff_index, NHe_eff_index, dNH_eff, dNHe_eff));
                            }
                          shieldFactor = S1 + S2 + S3;
                        }
                      photoion_rate = chimesRateTables.NonEqIon->NonEqRates[ns].sigmaphot[l][i][0] *
                                      myGasVars->isotropic_photon_density[l] * LIGHTSPEED * shieldFactor;
                    }
                  else
                    {
                      /* For these species, UV attenuation is
                       * by dust. The gamma_d values are stored
                       * in the shieldFactor tables. */
                      shieldFactor =
                          exp(-pow(10.0, ((double)chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor1D[l][0][i][0])) * extinction);
                      photoion_rate = chimesRateTables.NonEqIon->NonEqRates[ns].sigmaphot[l][i][0] *
                                      myGasVars->isotropic_photon_density[l] * LIGHTSPEED * shieldFactor;
                    }

                  /* Calculate epsilon, the average excess photon energy per ionisation. */
                  if(chimesRateTables.NonEqIon->NonEqRates[ns].E_thresh[i] < 13.6 ||
                     (HI_column_density == 0.0 && H2_column_density == 0.0 && HeI_column_density == 0.0 && HeII_column_density == 0))
                    epsilon = chimesRateTables.NonEqIon->NonEqRates[ns].epsilon[l][i];
                  else
                    {
                      if(chimesRateTables.NonEqIon->NonEqRates[ns].E_thresh[i] < 15.4)
                        {
                          E1 = pow(10.0, interpol_1d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor1D[l][1][i],
                                                            NHI_index, dNHI));
                          E4 = pow(10.0, interpol_1d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor1D[l][2][i],
                                                            NHI_index, dNHI));
                          E2 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][2][i],
                                                            NH_eff_index, NHeI_index, dNH_eff, dNHeI));
                          E5 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][4][i],
                                                            NH_eff_index, NHeI_index, dNH_eff, dNHeI));
                          E3 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][3][i],
                                                            NH_eff_index, NHe_eff_index, dNH_eff, dNHe_eff));
                          E6 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][5][i],
                                                            NH_eff_index, NHe_eff_index, dNH_eff, dNHe_eff));
                        }
                      else if(chimesRateTables.NonEqIon->NonEqRates[ns].E_thresh[i] < 54.42)
                        {
                          E1 = 0.0;
                          E4 = 0.0;
                          E2 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][2][i],
                                                            NH_eff_index, NHeI_index, dNH_eff, dNHeI));
                          E5 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][4][i],
                                                            NH_eff_index, NHeI_index, dNH_eff, dNHeI));
                          E3 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][3][i],
                                                            NH_eff_index, NHe_eff_index, dNH_eff, dNHe_eff));
                          E6 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][5][i],
                                                            NH_eff_index, NHe_eff_index, dNH_eff, dNHe_eff));
                        }
                      else
                        {
                          E1 = 0.0;
                          E4 = 0.0;
                          E2 = 0.0;
                          E5 = 0.0;
                          E3 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][3][i],
                                                            NH_eff_index, NHe_eff_index, dNH_eff, dNHe_eff));
                          E6 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][5][i],
                                                            NH_eff_index, NHe_eff_index, dNH_eff, dNHe_eff));
                        }
                      epsilon = (E1 + E2 + E3) / (E4 + E5 + E6);
                    }
                  total_cooling += photoheating(
                      myGasVars->abundances[myGlobalVars->speciesIndices[chimesRateTables.NonEqIon->IonIndexBegin[ns] + i]],
                      photoion_rate, epsilon, myGasVars->nH_tot);
                }

              if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI || (chimesRateTables.NonEqIon->IonIndexBegin[ns] + i) == HeI)
                {
                  /* Recall: 20 eV per PRIMARY ionisation */
                  total_cooling +=
                      cosmic_ray_heating(
                          myGasVars->abundances[myGlobalVars->speciesIndices[chimesRateTables.NonEqIon->IonIndexBegin[ns] + i]],
                          this_all_rates->BensRates[ns].cosmicRays[i], myGasVars->nH_tot) /
                      (1.0 + cr_secondary_ionisation(myGasVars->abundances[myGlobalVars->speciesIndices[HII]],
                                                     chimesRateTables.NonEqIon->IonIndexBegin[ns]));
                }
              else
                total_cooling += cosmic_ray_heating(
                    myGasVars->abundances[myGlobalVars->speciesIndices[chimesRateTables.NonEqIon->IonIndexBegin[ns] + i]],
                    this_all_rates->BensRates[ns].cosmicRays[i], myGasVars->nH_tot);
            }
        }

      if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
        {
          /* The Hydrogen photoion arrays also include H- and
           * H2 -> H2+, so we need to now add photoheating
           * from these. */
          for(i = 1; i < 3; i++)
            {
              for(l = 0; l < myGlobalVars->N_spectra; l++)
                {
                  /* Calculate photoionisation rate for given UV spectrum. */
                  if(chimesRateTables.NonEqIon->NonEqRates[ns].E_thresh[i] >= 13.6)
                    {
                      if(HI_column_density == 0.0 && H2_column_density == 0.0 && HeI_column_density == 0.0 &&
                         HeII_column_density == 0.0)
                        shieldFactor = 1.0;
                      else
                        {
                          if(chimesRateTables.NonEqIon->NonEqRates[ns].E_thresh[i] < 15.4)
                            {
                              S1 = pow(10.0, interpol_1d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor1D[l][0][i],
                                                                NHI_index, dNHI));
                              S2 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][0][i],
                                                                NH_eff_index, NHeI_index, dNH_eff, dNHeI));
                              S3 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][1][i],
                                                                NH_eff_index, NHe_eff_index, dNH_eff, dNHe_eff));
                            }
                          else if(chimesRateTables.NonEqIon->NonEqRates[ns].E_thresh[i] < 54.42)
                            {
                              S1 = 0.0;
                              S2 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][0][i],
                                                                NH_eff_index, NHeI_index, dNH_eff, dNHeI));
                              S3 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][1][i],
                                                                NH_eff_index, NHe_eff_index, dNH_eff, dNHe_eff));
                            }
                          else
                            {
                              S1 = 0.0;
                              S2 = 0.0;
                              S3 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][1][i],
                                                                NH_eff_index, NHe_eff_index, dNH_eff, dNHe_eff));
                            }
                          shieldFactor = S1 + S2 + S3;
                        }
                      photoion_rate = chimesRateTables.NonEqIon->NonEqRates[ns].sigmaphot[l][i][0] *
                                      myGasVars->isotropic_photon_density[l] * LIGHTSPEED * shieldFactor;
                    }
                  else
                    {
                      /* For these species, UV attenuation is
                       * by dust. The gamma_d values are stored
                       * in the shieldFactor tables. */
                      shieldFactor =
                          exp(-pow(10.0, ((double)chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor1D[l][0][i][0])) * extinction);
                      photoion_rate = chimesRateTables.NonEqIon->NonEqRates[ns].sigmaphot[l][i][0] *
                                      myGasVars->isotropic_photon_density[l] * LIGHTSPEED * shieldFactor;
                    }

                  /* Calculate epsilon, the average excess photon energy per ionisation. */
                  if(chimesRateTables.NonEqIon->NonEqRates[ns].E_thresh[i] < 13.6 ||
                     (HI_column_density == 0.0 && H2_column_density == 0.0 && HeI_column_density == 0.0 && HeII_column_density == 0))
                    epsilon = chimesRateTables.NonEqIon->NonEqRates[ns].epsilon[l][i];
                  else
                    {
                      if(chimesRateTables.NonEqIon->NonEqRates[ns].E_thresh[i] < 15.4)
                        {
                          E1 = pow(10.0, interpol_1d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor1D[l][1][i],
                                                            NHI_index, dNHI));
                          E4 = pow(10.0, interpol_1d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor1D[l][2][i],
                                                            NHI_index, dNHI));
                          E2 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][2][i],
                                                            NH_eff_index, NHeI_index, dNH_eff, dNHeI));
                          E5 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][4][i],
                                                            NH_eff_index, NHeI_index, dNH_eff, dNHeI));
                          E3 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][3][i],
                                                            NH_eff_index, NHe_eff_index, dNH_eff, dNHe_eff));
                          E6 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][5][i],
                                                            NH_eff_index, NHe_eff_index, dNH_eff, dNHe_eff));
                        }
                      else if(chimesRateTables.NonEqIon->NonEqRates[ns].E_thresh[i] < 54.42)
                        {
                          E1 = 0.0;
                          E4 = 0.0;
                          E2 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][2][i],
                                                            NH_eff_index, NHeI_index, dNH_eff, dNHeI));
                          E5 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][4][i],
                                                            NH_eff_index, NHeI_index, dNH_eff, dNHeI));
                          E3 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][3][i],
                                                            NH_eff_index, NHe_eff_index, dNH_eff, dNHe_eff));
                          E6 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][5][i],
                                                            NH_eff_index, NHe_eff_index, dNH_eff, dNHe_eff));
                        }
                      else
                        {
                          E1 = 0.0;
                          E4 = 0.0;
                          E2 = 0.0;
                          E5 = 0.0;
                          E3 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][3][i],
                                                            NH_eff_index, NHe_eff_index, dNH_eff, dNHe_eff));
                          E6 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[l][5][i],
                                                            NH_eff_index, NHe_eff_index, dNH_eff, dNHe_eff));
                        }
                      epsilon = (E1 + E2 + E3) / (E4 + E5 + E6);
                    }
                  if(i == 1)
                    total_cooling += photoheating(myGasVars->abundances[myGlobalVars->speciesIndices[Hm]], photoion_rate, epsilon,
                                                  myGasVars->nH_tot);
                  else
                    total_cooling += photoheating(myGasVars->abundances[myGlobalVars->speciesIndices[H2]], photoion_rate, epsilon,
                                                  myGasVars->nH_tot);
                }
            }
        }
    }

  /* For CI and OI cooling, we replace Ben's rates with
   * our old per-ion cooling table below 10^4 K */

  if(log10(myGasVars->temperature) < chimesRateTables.nei_cooling_temperature[chimesRateTables.nei_cooling_table_dimensions[0] - 1])
    {
      /* Note that the factors of xe have already been taken into account within the table itself. */
      /* Note also that we have to divide by an EXTRA factor of nH here */
      total_cooling += pow(10, interpol_4d_mydbl(chimesRateTables.nei_cooling_CI, m, j, k, p, dT_2, dnHI, dne, dnHII)) *
                       myGasVars->abundances[myGlobalVars->speciesIndices[CI]] / myGasVars->nH_tot;
      total_cooling += pow(10, interpol_4d_mydbl(chimesRateTables.nei_cooling_OI, m, j, k, p, dT_2, dnHI, dne, dnHII)) *
                       myGasVars->abundances[myGlobalVars->speciesIndices[OI]] / myGasVars->nH_tot;
    }

  /* For CII, NII, SiII & FeII we replace Ben's rates with
   * those calculated using Chianti. */
  if(log10(myGasVars->temperature) <
     chimesRateTables.chianti_cooling_temperature[chimesRateTables.chianti_cooling_table_dimensions[0] - 1])
    {
      total_cooling +=
          pow(10,
              interpol_2d_mydbl(chimesRateTables.chianti_cooling_CII, chianti_T_index, chianti_ne_index, chianti_dT, chianti_dne)) *
          myGasVars->abundances[myGlobalVars->speciesIndices[CII]] * myGasVars->abundances[myGlobalVars->speciesIndices[elec]];
      total_cooling +=
          pow(10,
              interpol_2d_mydbl(chimesRateTables.chianti_cooling_NII, chianti_T_index, chianti_ne_index, chianti_dT, chianti_dne)) *
          myGasVars->abundances[myGlobalVars->speciesIndices[NII]] * myGasVars->abundances[myGlobalVars->speciesIndices[elec]];
      total_cooling +=
          pow(10,
              interpol_2d_mydbl(chimesRateTables.chianti_cooling_SiII, chianti_T_index, chianti_ne_index, chianti_dT, chianti_dne)) *
          myGasVars->abundances[myGlobalVars->speciesIndices[SiII]] * myGasVars->abundances[myGlobalVars->speciesIndices[elec]];
      total_cooling +=
          pow(10,
              interpol_2d_mydbl(chimesRateTables.chianti_cooling_FeII, chianti_T_index, chianti_ne_index, chianti_dT, chianti_dne)) *
          myGasVars->abundances[myGlobalVars->speciesIndices[FeII]] * myGasVars->abundances[myGlobalVars->speciesIndices[elec]];
    }

  /* Now calculate cooling from molecular species (only below T_mol)*/
  if(myGasVars->temperature <= myGlobalVars->T_mol)
    {
      /* Cooling from H2 */
      total_cooling += H2_rovibrational_cooling(
          myGasVars->temperature, myGasVars->abundances[myGlobalVars->speciesIndices[H2]],
          myGasVars->abundances[myGlobalVars->speciesIndices[HI]], myGasVars->abundances[myGlobalVars->speciesIndices[HII]],
          myGasVars->abundances[myGlobalVars->speciesIndices[HeI]], myGasVars->abundances[myGlobalVars->speciesIndices[elec]],
          myGasVars->nH_tot); /*GA08 */
      total_cooling +=
          7.2e-12 *
          (k9(myGasVars->temperature, myGasVars->nH_tot, myGasVars->abundances[myGlobalVars->speciesIndices[HI]],
              myGasVars->abundances[myGlobalVars->speciesIndices[H2]], myGasVars->abundances[myGlobalVars->speciesIndices[HeI]]) *
               myGasVars->abundances[myGlobalVars->speciesIndices[HI]] +
           (this_all_rates->rate_10 / myGasVars->nH_tot) * myGasVars->abundances[myGlobalVars->speciesIndices[H2]]) *
          myGasVars->abundances[myGlobalVars->speciesIndices[H2]]; /* H2 collis diss */
      total_cooling -= 6.4e-13 * this_all_rates->rate_53 * myGasVars->abundances[myGlobalVars->speciesIndices[H2]] /
                       myGasVars->nH_tot; /* H2 photodiss */
      n_cr = H2_crit_density(myGasVars->temperature, myGasVars->abundances[myGlobalVars->speciesIndices[HI]],
                             myGasVars->abundances[myGlobalVars->speciesIndices[H2]]);
      n    = calculate_total_number_density(myGasVars->abundances, myGasVars->nH_tot, myGlobalVars);
      total_cooling -= 2.7e-11 * this_all_rates->rate_53 * myGasVars->abundances[myGlobalVars->speciesIndices[H2]] * (n / (n + n_cr)) /
                       myGasVars->nH_tot; /* H2 UV pumping */
      total_cooling -=
          (2.93e-12 * (this_all_rates->rate_2 / myGasVars->nH_tot) * myGasVars->abundances[myGlobalVars->speciesIndices[Hm]] +
           5.65e-12 * (this_all_rates->rate_4 / myGasVars->nH_tot) * myGasVars->abundances[myGlobalVars->speciesIndices[H2p]]) *
          myGasVars->abundances[myGlobalVars->speciesIndices[HI]] * (n / (n + n_cr)); /* Gas-phase H2 formation */
      total_cooling -= 7.16e-12 * (this_all_rates->rate_60 / myGasVars->nH_tot) *
                       myGasVars->abundances[myGlobalVars->speciesIndices[HI]] * n * (n / (n + n_cr)) /
                       myGasVars->nH_tot; /* Dust H2 formation */
      total_cooling +=
          cosmic_ray_heating(myGasVars->abundances[myGlobalVars->speciesIndices[H2]], 2.0 * myGasVars->cr_rate, myGasVars->nH_tot);

      /* Cooling from CO */
      if(myGlobalVars->element_included[0] == 1 && myGlobalVars->element_included[2] == 1)
        {
          if(myGlobalVars->StaticMolCooling == 1)
            N_eff = 1.0e5 * CO_column_density /
                    (sqrt(3.0 * BOLTZMANNCGS * myGasVars->temperature /
                          (PROTON_MASS * calculate_mean_molecular_weight(myGasVars, myGlobalVars)))); /* units: cm^-2 per km s^-1 */
          else
            N_eff = fabs(1.0e5 * max(myGasVars->abundances[myGlobalVars->speciesIndices[CO]], 0.0) * myGasVars->nH_tot /
                         myGasVars->divVel); /* units: cm^-2 per km s^-1 */
          total_cooling += CO_rotational_cooling(
              myGasVars->temperature, N_eff, myGasVars->abundances[myGlobalVars->speciesIndices[CO]],
              myGasVars->abundances[myGlobalVars->speciesIndices[H2]], myGasVars->abundances[myGlobalVars->speciesIndices[HI]],
              myGasVars->abundances[myGlobalVars->speciesIndices[elec]], myGasVars->nH_tot, myGlobalVars);
          total_cooling += CO_vibrational_cooling(
              myGasVars->temperature, N_eff, myGasVars->abundances[myGlobalVars->speciesIndices[CO]],
              myGasVars->abundances[myGlobalVars->speciesIndices[H2]], myGasVars->abundances[myGlobalVars->speciesIndices[HI]],
              myGasVars->abundances[myGlobalVars->speciesIndices[elec]], myGasVars->nH_tot, myGlobalVars);
        }

      if(myGlobalVars->element_included[2] == 1)
        {
          /* Cooling from H2O */
          if(myGlobalVars->StaticMolCooling == 1)
            N_eff = 1.0e5 * H2O_column_density /
                    (sqrt(3.0 * BOLTZMANNCGS * myGasVars->temperature /
                          (PROTON_MASS * calculate_mean_molecular_weight(myGasVars, myGlobalVars)))); /* units: cm^-2 per km s^-1 */
          else
            N_eff = fabs(1.0e5 * max(myGasVars->abundances[myGlobalVars->speciesIndices[H2O]], 0.0) * myGasVars->nH_tot /
                         myGasVars->divVel); /* units: cm^-2 per km s^-1 */
          total_cooling += H2O_rotational_cooling(
              myGasVars->temperature, N_eff, myGasVars->abundances[myGlobalVars->speciesIndices[H2O]],
              myGasVars->abundances[myGlobalVars->speciesIndices[H2]], myGasVars->abundances[myGlobalVars->speciesIndices[HI]],
              myGasVars->abundances[myGlobalVars->speciesIndices[elec]], myGasVars->nH_tot, myGlobalVars);
          total_cooling += H2O_vibrational_cooling(
              myGasVars->temperature, N_eff, myGasVars->abundances[myGlobalVars->speciesIndices[H2O]],
              myGasVars->abundances[myGlobalVars->speciesIndices[H2]], myGasVars->abundances[myGlobalVars->speciesIndices[HI]],
              myGasVars->abundances[myGlobalVars->speciesIndices[elec]], myGasVars->nH_tot, myGlobalVars);

          /* Cooling from OH */
          total_cooling += OH_rotational_cooling(myGasVars->temperature, OH_column_density,
                                                 sqrt(3.0 * BOLTZMANNCGS * myGasVars->temperature /
                                                      (PROTON_MASS * calculate_mean_molecular_weight(myGasVars, myGlobalVars))),
                                                 max(myGasVars->abundances[myGlobalVars->speciesIndices[OH]], 0.0), myGasVars->nH_tot,
                                                 n, extinction);
        }

      /* Dust grains */
      dust_G = 0.0;
      for(l = 0; l < myGlobalVars->N_spectra; l++)
        dust_G += myGasVars->isotropic_photon_density[l] * LIGHTSPEED * myGasVars->dust_G_parameter[l];

      total_cooling += photoelectric_grain_heating(myGasVars->temperature, myGasVars->nH_tot,
                                                   myGasVars->abundances[myGlobalVars->speciesIndices[elec]] * myGasVars->nH_tot, n,
                                                   myGasVars->metallicity, dust_G, extinction);
      total_cooling += gas_grain_transfer(myGasVars->temperature, myGlobalVars->grain_temperature, myGasVars->metallicity);
      total_cooling += grain_surface_recombination(myGasVars->temperature, myGasVars->metallicity,
                                                   myGasVars->abundances[myGlobalVars->speciesIndices[elec]] * myGasVars->nH_tot,
                                                   myGasVars->nH_tot, dust_G, extinction);
    }

  /* Compton Cooling */
  total_cooling += compton_cooling(myGasVars->temperature, myGlobalVars->cmb_temperature,
                                   myGasVars->abundances[myGlobalVars->speciesIndices[elec]], myGasVars->nH_tot);

  return (total_cooling * pow(myGasVars->nH_tot, 2.0)) - myGasVars->constant_heating_rate;
}

void do_equilibrium_cooling(void *user_data)
{
  int i, maxIter;
  double u, u_old, u_upper, u_lower, du, network_size, LambdaNet, dt;
  struct gasVariables *myGasVars;
  struct globalVariables *myGlobalVars;
  struct Species_Structure *species;
  UserData *data;
  data = (UserData *)user_data;

  myGasVars    = data->myGasVars;
  myGlobalVars = data->myGlobalVars;
  species      = data->species;
  network_size = data->network_size;
  dt           = myGasVars->hydro_timestep;

  myGasVars->temperature = max(myGasVars->temperature, myGasVars->TempFloor);

  if(myGasVars->ForceEqOn == 1)
    set_equilibrium_abundances_from_tables(myGasVars, myGlobalVars);
  else
    set_equilibrium_abundances(data); /* Note: the column densities in 'data' are also updated in this routine. */

  update_rates(myGasVars, myGlobalVars, *(data->HI_column), *(data->H2_column), *(data->HeI_column), *(data->HeII_column),
               *(data->CO_column), *(data->extinction), data->this_all_rates);
  update_T_dependent_rates(myGasVars, myGlobalVars, data->this_all_rates);

  LambdaNet = -calculate_total_cooling_rate(myGasVars, myGlobalVars, *(data->HI_column), *(data->HeI_column), *(data->HeII_column),
                                            *(data->H2_column), *(data->CO_column), *(data->H2O_column), *(data->OH_column),
                                            *(data->extinction), data->this_all_rates);

  if(myGasVars->temperature <= myGasVars->TempFloor && LambdaNet <= 0.0)
    return;

  u = myGasVars->temperature * 1.5 * calculate_total_number_density(myGasVars->abundances, myGasVars->nH_tot, myGlobalVars) *
      BOLTZMANNCGS;
  u_old   = u;
  u_upper = u;
  u_lower = u;

  /* If the cooling rate is small, take explicit solution. */
  if(fabs(LambdaNet * dt) < 0.10 * u_old)
    {
      u = u_old + LambdaNet * dt;
      myGasVars->temperature =
          max(u / (1.5 * calculate_total_number_density(myGasVars->abundances, myGasVars->nH_tot, myGlobalVars) * BOLTZMANNCGS),
              myGasVars->TempFloor);
      if(myGasVars->ForceEqOn == 1)
        set_equilibrium_abundances_from_tables(myGasVars, myGlobalVars);
      else
        set_equilibrium_abundances(data);

      // Check that explicit solution is valid
      update_rates(myGasVars, myGlobalVars, *(data->HI_column), *(data->H2_column), *(data->HeI_column), *(data->HeII_column),
                   *(data->CO_column), *(data->extinction), data->this_all_rates);
      update_T_dependent_rates(myGasVars, myGlobalVars, data->this_all_rates);
      LambdaNet = -calculate_total_cooling_rate(myGasVars, myGlobalVars, *(data->HI_column), *(data->HeI_column), *(data->HeII_column),
                                                *(data->H2_column), *(data->CO_column), *(data->H2O_column), *(data->OH_column),
                                                *(data->extinction), data->this_all_rates);
      if(fabs(LambdaNet * dt) < 0.10 * u_old && fabs(LambdaNet * dt) < 0.10 * u)
        return;
      else
        {
          /* New cooling rate has increased, and explicit
           * solution is no longer valid. Reset and
           * continue with implicit solution. */
          u = u_old;
          myGasVars->temperature =
              max(u / (1.5 * calculate_total_number_density(myGasVars->abundances, myGasVars->nH_tot, myGlobalVars) * BOLTZMANNCGS),
                  myGasVars->TempFloor);
          if(myGasVars->ForceEqOn == 1)
            set_equilibrium_abundances_from_tables(myGasVars, myGlobalVars);
          else
            set_equilibrium_abundances(data);
        }
    }

  i       = 0;
  maxIter = 150;

  /* Bracketing */
  if(u - u_old - LambdaNet * dt < 0.0) /* heating */
    {
      u_upper *= sqrt(1.2);
      u_lower /= sqrt(1.2);

      myGasVars->temperature =
          u_upper / (1.5 * calculate_total_number_density(myGasVars->abundances, myGasVars->nH_tot, myGlobalVars) * BOLTZMANNCGS);
      if(myGasVars->ForceEqOn == 1)
        set_equilibrium_abundances_from_tables(myGasVars, myGlobalVars);
      else
        set_equilibrium_abundances(data);

      update_rates(myGasVars, myGlobalVars, *(data->HI_column), *(data->H2_column), *(data->HeI_column), *(data->HeII_column),
                   *(data->CO_column), *(data->extinction), data->this_all_rates);
      update_T_dependent_rates(myGasVars, myGlobalVars, data->this_all_rates);
      LambdaNet = -calculate_total_cooling_rate(myGasVars, myGlobalVars, *(data->HI_column), *(data->HeI_column), *(data->HeII_column),
                                                *(data->H2_column), *(data->CO_column), *(data->H2O_column), *(data->OH_column),
                                                *(data->extinction), data->this_all_rates);

      while(u_upper - u_old - LambdaNet * dt < 0.0 && i < maxIter)
        {
          u_upper *= 1.2;
          u_lower *= 1.2;

          myGasVars->temperature =
              u_upper / (1.5 * calculate_total_number_density(myGasVars->abundances, myGasVars->nH_tot, myGlobalVars) * BOLTZMANNCGS);
          if(myGasVars->ForceEqOn == 1)
            set_equilibrium_abundances_from_tables(myGasVars, myGlobalVars);
          else
            set_equilibrium_abundances(data);

          update_rates(myGasVars, myGlobalVars, *(data->HI_column), *(data->H2_column), *(data->HeI_column), *(data->HeII_column),
                       *(data->CO_column), *(data->extinction), data->this_all_rates);
          update_T_dependent_rates(myGasVars, myGlobalVars, data->this_all_rates);
          LambdaNet = -calculate_total_cooling_rate(myGasVars, myGlobalVars, *(data->HI_column), *(data->HeI_column),
                                                    *(data->HeII_column), *(data->H2_column), *(data->CO_column), *(data->H2O_column),
                                                    *(data->OH_column), *(data->extinction), data->this_all_rates);

          i++;
        }

      if(i == maxIter)
        printf("WARNING: Problem with eqm cooling finding the upper bound.\n");
    }
  else /* cooling */
    {
      u_upper *= sqrt(1.2);
      u_lower /= sqrt(1.2);

      myGasVars->temperature =
          u_lower / (1.5 * calculate_total_number_density(myGasVars->abundances, myGasVars->nH_tot, myGlobalVars) * BOLTZMANNCGS);
      if(myGasVars->temperature <= myGasVars->TempFloor)
        {
          myGasVars->temperature = myGasVars->TempFloor;
          u_lower                = myGasVars->TempFloor * 1.5 *
                    calculate_total_number_density(myGasVars->abundances, myGasVars->nH_tot, myGlobalVars) * BOLTZMANNCGS;
          if(myGasVars->ForceEqOn == 1)
            set_equilibrium_abundances_from_tables(myGasVars, myGlobalVars);
          else
            set_equilibrium_abundances(data);

          update_rates(myGasVars, myGlobalVars, *(data->HI_column), *(data->H2_column), *(data->HeI_column), *(data->HeII_column),
                       *(data->CO_column), *(data->extinction), data->this_all_rates);
          update_T_dependent_rates(myGasVars, myGlobalVars, data->this_all_rates);
          LambdaNet = -calculate_total_cooling_rate(myGasVars, myGlobalVars, *(data->HI_column), *(data->HeI_column),
                                                    *(data->HeII_column), *(data->H2_column), *(data->CO_column), *(data->H2O_column),
                                                    *(data->OH_column), *(data->extinction), data->this_all_rates);
        }
      else
        {
          if(myGasVars->ForceEqOn == 1)
            set_equilibrium_abundances_from_tables(myGasVars, myGlobalVars);
          else
            set_equilibrium_abundances(data);

          update_rates(myGasVars, myGlobalVars, *(data->HI_column), *(data->H2_column), *(data->HeI_column), *(data->HeII_column),
                       *(data->CO_column), *(data->extinction), data->this_all_rates);
          update_T_dependent_rates(myGasVars, myGlobalVars, data->this_all_rates);
          LambdaNet = -calculate_total_cooling_rate(myGasVars, myGlobalVars, *(data->HI_column), *(data->HeI_column),
                                                    *(data->HeII_column), *(data->H2_column), *(data->CO_column), *(data->H2O_column),
                                                    *(data->OH_column), *(data->extinction), data->this_all_rates);

          while(u_lower - u_old - LambdaNet * dt > 0.0 && i < maxIter)
            {
              u_upper /= 1.2;
              u_lower /= 1.2;

              myGasVars->temperature =
                  u_lower /
                  (1.5 * calculate_total_number_density(myGasVars->abundances, myGasVars->nH_tot, myGlobalVars) * BOLTZMANNCGS);
              if(myGasVars->temperature <= myGasVars->TempFloor)
                {
                  myGasVars->temperature = myGasVars->TempFloor;
                  u_lower                = myGasVars->TempFloor * 1.5 *
                            calculate_total_number_density(myGasVars->abundances, myGasVars->nH_tot, myGlobalVars) * BOLTZMANNCGS;
                  if(myGasVars->ForceEqOn == 1)
                    set_equilibrium_abundances_from_tables(myGasVars, myGlobalVars);
                  else
                    set_equilibrium_abundances(data);

                  update_rates(myGasVars, myGlobalVars, *(data->HI_column), *(data->H2_column), *(data->HeI_column),
                               *(data->HeII_column), *(data->CO_column), *(data->extinction), data->this_all_rates);
                  update_T_dependent_rates(myGasVars, myGlobalVars, data->this_all_rates);
                  LambdaNet = -calculate_total_cooling_rate(
                      myGasVars, myGlobalVars, *(data->HI_column), *(data->HeI_column), *(data->HeII_column), *(data->H2_column),
                      *(data->CO_column), *(data->H2O_column), *(data->OH_column), *(data->extinction), data->this_all_rates);
                  break;
                }

              if(myGasVars->ForceEqOn == 1)
                set_equilibrium_abundances_from_tables(myGasVars, myGlobalVars);
              else
                set_equilibrium_abundances(data);

              update_rates(myGasVars, myGlobalVars, *(data->HI_column), *(data->H2_column), *(data->HeI_column), *(data->HeII_column),
                           *(data->CO_column), *(data->extinction), data->this_all_rates);
              update_T_dependent_rates(myGasVars, myGlobalVars, data->this_all_rates);
              LambdaNet = -calculate_total_cooling_rate(
                  myGasVars, myGlobalVars, *(data->HI_column), *(data->HeI_column), *(data->HeII_column), *(data->H2_column),
                  *(data->CO_column), *(data->H2O_column), *(data->OH_column), *(data->extinction), data->this_all_rates);

              i++;
            }
          if(i == maxIter)
            printf("WARNING: Problem with eqm cooling finding the lower bound.\n");
        }
      if(u_lower - u_old - LambdaNet * dt > 0.0 && i < maxIter)
        return; /* u_lower reached TempFloor, but is still above converged solution. */
    }

  /* Iterate to convergence */
  i = 0;

  do
    {
      u = 0.5 * (u_lower + u_upper);

      myGasVars->temperature =
          u / (1.5 * calculate_total_number_density(myGasVars->abundances, myGasVars->nH_tot, myGlobalVars) * BOLTZMANNCGS);
      if(myGasVars->ForceEqOn == 1)
        set_equilibrium_abundances_from_tables(myGasVars, myGlobalVars);
      else
        set_equilibrium_abundances(data);

      update_rates(myGasVars, myGlobalVars, *(data->HI_column), *(data->H2_column), *(data->HeI_column), *(data->HeII_column),
                   *(data->CO_column), *(data->extinction), data->this_all_rates);
      update_T_dependent_rates(myGasVars, myGlobalVars, data->this_all_rates);
      LambdaNet = -calculate_total_cooling_rate(myGasVars, myGlobalVars, *(data->HI_column), *(data->HeI_column), *(data->HeII_column),
                                                *(data->H2_column), *(data->CO_column), *(data->H2O_column), *(data->OH_column),
                                                *(data->extinction), data->this_all_rates);

      if(u - u_old - LambdaNet * dt > 0.0)
        u_upper = u;
      else
        u_lower = u;

      du = u_upper - u_lower;
      i++;
    }
  while(fabs(du / u) > 1.0e-6 && i < maxIter);

  if(i >= maxIter)
    printf("WARNING: eqm cooling failed to converge.\n");

  return;
}
