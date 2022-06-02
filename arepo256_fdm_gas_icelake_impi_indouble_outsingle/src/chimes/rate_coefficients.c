#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "proto.h"

double k1(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction1, T_index, dT)); }

double k2(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction2, T_index, dT)); }

double k3(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction3, T_index, dT)); }

double k5(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction5, T_index, dT)); }

double k6(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction6, T_index, dT)); }

double k7(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction7, T_index, dT)); }

double k10(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction10, T_index, dT)); }

double k11(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction11, T_index, dT)); }

double k13(int T_index, double dT, double HI_column_density)
{
  /* To determine case A or case B, use HI cross-section at
   * 1 Ryd (i.e. assume recombination radiation to be monochromatic). */
  double sigma_HI = 6.3463e-18; /* From Verner et al. (1996). */

  if((sigma_HI * HI_column_density) < 1.0)
    return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction13A, T_index, dT));
  else
    return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction13B, T_index, dT));
}

double k15(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction15, T_index, dT)); }

double k16(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction16, T_index, dT)); }

double k17(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction17, T_index, dT)); }

double k24(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction24, T_index, dT)); }

double k25rr(int T_index, double dT, double HeI_column_density)
{
  /* HeI cross-section at 24.59 eV (Verner et al. 1996) */
  double sigma_HeI = 7.4347e-18;

  if((sigma_HeI * HeI_column_density) < 1.0)
    return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction25rrA, T_index, dT));
  else
    return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction25rrB, T_index, dT));
}

double k25di(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction25di, T_index, dT)); }

double k26(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction26, T_index, dT)); }

double k27(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction27, T_index, dT)); }

double k83(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction83, T_index, dT)); }

double k87(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction87, T_index, dT)); }

double k97(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction97, T_index, dT)); }

double k101(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction101, T_index, dT)); }

double k106(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction106, T_index, dT)); }

double k107(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction107, T_index, dT)); }

double k112(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction112, T_index, dT)); }

double k113(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction113, T_index, dT)); }

double k114(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction114, T_index, dT)); }

double k116(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction116, T_index, dT)); }

double k117(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction117, T_index, dT)); }

double k120(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction120, T_index, dT)); }

double k121(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction121, T_index, dT)); }

double k122(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction122, T_index, dT)); }

double k124(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction124, T_index, dT)); }

double k128(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction128, T_index, dT)); }

double k129(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction129, T_index, dT)); }

double k130(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction130, T_index, dT)); }

double k131(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction131, T_index, dT)); }

double k133(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction133, T_index, dT)); }

double k134(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction134, T_index, dT)); }

double k135(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction135, T_index, dT)); }

double k136(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction136, T_index, dT)); }

double k137(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction137, T_index, dT)); }

double k138(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction138, T_index, dT)); }

double k139(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction139, T_index, dT)); }

double k142(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction142, T_index, dT)); }

double k147(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction147, T_index, dT)); }

double k150(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction150, T_index, dT)); }

double k186(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction186, T_index, dT)); }

double k187(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction187, T_index, dT)); }

double k188(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction188, T_index, dT)); }

double k189(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction189, T_index, dT)); }

double k190(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction190, T_index, dT)); }

double k191(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction191, T_index, dT)); }

double k192(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction192, T_index, dT)); }

double k193(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction193, T_index, dT)); }

double k194(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction194, T_index, dT)); }

double k195(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction195, T_index, dT)); }

double k196(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction196, T_index, dT)); }

double k197(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction197, T_index, dT)); }

double k198(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction198, T_index, dT)); }

double k199(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction199, T_index, dT)); }

double k200(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction200, T_index, dT)); }

double k201(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction201, T_index, dT)); }

double k202(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction202, T_index, dT)); }

double k203(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction203, T_index, dT)); }

double k204(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction204, T_index, dT)); }

double k205(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction205, T_index, dT)); }

double k206(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction206, T_index, dT)); }

double k207(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction207, T_index, dT)); }

double k220(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction220, T_index, dT)); }

double k221(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction221, T_index, dT)); }

double k222(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction222, T_index, dT)); }

double k223(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction223, T_index, dT)); }

double k225(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction225, T_index, dT)); }

double k226(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction226, T_index, dT)); }

double k227(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction227, T_index, dT)); }

double k228(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction228, T_index, dT)); }

double k229(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction229, T_index, dT)); }

double k230(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction230, T_index, dT)); }

double k231(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction231, T_index, dT)); }

double k232(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction232, T_index, dT)); }

double k233(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction233, T_index, dT)); }

double k234(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction234, T_index, dT)); }

double k235(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction235, T_index, dT)); }

double k236(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction236, T_index, dT)); }

double k237(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction237, T_index, dT)); }

double k238(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction238, T_index, dT)); }

double k283(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction283, T_index, dT)); }

double k44(void) { return 2.1e-9; }

double grain_recomb_jac(double T, double xe, double nH_tot, double dust_G, double extinction, double dust_ratio, double a, double b,
                        double c, double d, double e)
{
  /* Calculates dk/dxe for the grain
   * recombination reactions. */
  if(dust_G == 0.0 || xe == 0.0)
    return 0.0;
  else
    {
      double psi, dpsi_dxe, dk_dpsi;
      psi = dust_G * exp(-extinction * G0_GAMMA) * pow(T, 0.5) / (xe * nH_tot);
      if(psi > 1.0e-90)
        {
          dpsi_dxe = -psi / xe;
          dk_dpsi  = -a * dust_ratio * (b * c * pow(psi, c - 1.0) * (1.0 + d * pow(psi, e)) + b * d * e * pow(psi, c + e - 1.0)) /
                    pow(1.0 + b * pow(psi, c) * (1.0 + d * pow(psi, e)), 2.0);
          return dk_dpsi * dpsi_dxe;
        }
      else
        return 0.0;
    }
}

double k61(double T, double ne, double dust_G, double extinction, double dust_ratio)
{
  if(dust_G == 0.0 || ne == 0.0)
    {
      /* If ne == 0.0 the equation for psi divides by zero,
       * however since elec is a reactant the final rate will
       * be zero anyway, so it doesn't matter what we return */
      return 1.225e-13 * dust_ratio;
    }
  else
    {
      double psi;
      psi = dust_G * exp(-extinction * G0_GAMMA) * pow(T, 0.5) / ne;
      if(psi > 1.0e-90)
        return 1.225e-13 * dust_ratio /
               (1.0 + 8.074e-6 * pow(psi, 1.378) * (1.0 + 5.087e2 * pow(T, 0.01586) * pow(psi, (-0.4723 - 1.102e-5 * log(T)))));
      else
        return 1.225e-13 * dust_ratio;
    }
}

double k63(double T, double ne, double dust_G, double extinction, double dust_ratio)
{
  if(dust_G == 0.0 || ne == 0.0)
    return 5.572e-14 * dust_ratio;
  else
    {
      double psi = dust_G * exp(-extinction * G0_GAMMA) * pow(T, 0.5) / ne;
      if(psi > 1.0e-90)
        return 5.572e-14 * dust_ratio /
               (1.0 + 3.185e-7 * pow(psi, 1.512) * (1.0 + 5.115e3 * pow(T, 3.903e-7) * pow(psi, (-0.4956 - 5.494e-7 * log(T)))));
      else
        return 5.572e-14 * dust_ratio;
    }
}

double k64(double T, double ne, double dust_G, double extinction, double dust_ratio)
{
  if(dust_G == 0.0 || ne == 0.0)
    return 4.558e-13 * dust_ratio;
  else
    {
      double psi = dust_G * exp(-extinction * G0_GAMMA) * pow(T, 0.5) / ne;
      if(psi > 1.0e-90)
        return 4.558e-13 * dust_ratio /
               (1.0 + 6.089e-3 * pow(psi, 1.128) * (1.0 + 4.331e2 * pow(T, 0.04845) * pow(psi, (-0.812 - 1.333e-4 * log(T)))));
      else
        return 4.558e-13 * dust_ratio;
    }
}

double k65(double T, double ne, double dust_G, double extinction, double dust_ratio)
{
  if(dust_G == 0.0 || ne == 0.0)
    return 3.063e-14 * dust_ratio;
  else
    {
      double psi = dust_G * exp(-extinction * G0_GAMMA) * pow(T, 0.5) / ne;
      if(psi > 1.0e-90)
        return 3.063e-14 * dust_ratio /
               (1.0 + 8.074e-6 * pow(psi, 1.378) * (1.0 + 5.087e2 * pow(T, 0.01586) * pow(psi, (-0.4723 - 1.102e-5 * log(T)))));
      else
        return 3.063e-14 * dust_ratio;
    }
}

double k66(double T, double ne, double dust_G, double extinction, double dust_ratio)
{
  if(dust_G == 0.0 || ne == 0.0)
    return 2.166e-14 * dust_ratio;
  else
    {
      double psi = dust_G * exp(-extinction * G0_GAMMA) * pow(T, 0.5) / ne;
      if(psi > 1.0e-90)
        return 2.166e-14 * dust_ratio /
               (1.0 + 5.678e-8 * pow(psi, 1.874) * (1.0 + 4.375e4 * pow(T, 1.635e-6) * pow(psi, (-0.8964 - 7.538e-5 * log(T)))));
      else
        return 2.166e-14 * dust_ratio;
    }
}

double k77(double T, double ne, double dust_G, double extinction, double dust_ratio)
{
  if(dust_G == 0.0 || ne == 0.0)
    return 1.701e-14 * dust_ratio;
  else
    {
      double psi = dust_G * exp(-extinction * G0_GAMMA) * pow(T, 0.5) / ne;
      if(psi > 1.0e-90)
        return 1.701e-14 * dust_ratio /
               (1.0 + 9.554e-8 * pow(psi, 1.851) * (1.0 + 5.763e4 * pow(T, 4.116e-8) * pow(psi, (-0.9456 - 2.198e-5 * log(T)))));
      else
        return 1.701e-14 * dust_ratio;
    }
}

double k79(void) { return 2.6e-9; }

double k80(void)
{
  /* Note that we now use TH85's rate here,
   * as used by Cloudy. */
  return 1.9e-10;
}

double k4(void) { return 6.4e-10; }

double k8(double T, double nHtotal, double xH, double xH2, double xHe)
{
  double ncr_H  = pow(10.0, (3.0 - 0.416 * log10(T / 1.0e4) - 0.327 * pow(log10(T / 1.0e4), 2)));
  double ncr_H2 = pow(10.0, (4.845 - 1.3 * log10(T / 1.0e4) + 1.62 * pow(log10(T / 1.0e4), 2)));
  double ncr_He = pow(10.0, (5.0792 * (1.0 - 1.23e-5 * (T - 2000.0))));

  double ncr_inv = (xH / ncr_H) + (2.0 * xH2 / ncr_H2) + (xHe / ncr_He);
  double ncr;
  double kLTE, k0, intermediate;
  kLTE = 1.91e-9 * pow(T, 0.136) * exp(-53407.1 / T);
  k0   = max(4.49e-9 * pow(T, 0.11) * exp(-101858.0 / T), 1.0e-100);

  if(abs(ncr_inv * nHtotal) < 1.0e-50)
    intermediate = log10(k0);
  else
    {
      ncr          = 1.0 / ncr_inv;
      intermediate = ((nHtotal / ncr) / (1.0 + (nHtotal / ncr))) * log10(kLTE) + (1.0 / (1.0 + (nHtotal / ncr))) * log10(k0);
    }

  return pow(10.0, intermediate);
}

double k9(double T, double nHtotal, double xH, double xH2, double xHe)
{
  double ncr_H   = pow(10.0, (3.0 - 0.416 * log10(T / 1.0e4) - 0.327 * pow(log10(T / 1.0e4), 2)));
  double ncr_H2  = pow(10.0, (4.845 - 1.3 * log10(T / 1.0e4) + 1.62 * pow(log10(T / 1.0e4), 2)));
  double ncr_He  = pow(10.0, (5.0792 * (1.0 - 1.23e-5 * (T - 2000.0))));
  double ncr_inv = (xH / ncr_H) + (2.0 * xH2 / ncr_H2) + (xHe / ncr_He);
  double ncr;
  double kLTE, k0, intermediate;

  kLTE = 3.52e-9 * exp(-43900.0 / T);
  k0   = max(6.67e-12 * pow(T, 0.5) * exp(-(1.0 + (63593.0 / T))), 1.0e-100);

  if(abs(ncr_inv * nHtotal) < 1.0e-50)
    intermediate = log10(k0);
  else
    {
      ncr          = 1.0 / ncr_inv;
      intermediate = ((nHtotal / ncr) / (1.0 + (nHtotal / ncr))) * log10(kLTE) + (1.0 / (1.0 + (nHtotal / ncr))) * log10(k0);
    }

  return pow(10.0, intermediate);
}

double R52(double Go, double extinction) { return 5.7e-10 * Go * exp(-2.37 * extinction) / 1.7; }

double R53(double NH2, double temperature, double *photon_density, double extinction, double b_5, struct globalVariables *myGlobalVars,
           struct gasVariables *myGasVars)
{
  double S_H2, x, Ncrit, omega, alpha, b_therm, b_tot, rate;
  int l;
  if(NH2 < 1.0e7)
    S_H2 = 1.0;
  else
    {
      if(temperature < 3000.0)
        {
          Ncrit = 1.3e14 * (1.0 + pow(temperature / 600.0, 0.8));
          alpha = 1.4;
        }
      else if(temperature < 4000.0)
        {
          Ncrit = pow(temperature / 2.3e7, -3.8);
          alpha = pow(temperature / 4.5e3, -0.8);
        }
      else
        {
          Ncrit = 2.0e14;
          alpha = 1.1;
        }
      omega = 0.013 * pow((1.0 + pow(temperature / 2.7e3, 1.3)), 0.76923) * exp(-pow(temperature / 3.9e3, 14.6));
      x     = NH2 / Ncrit;

      b_therm = pow((2.0 * BOLTZMANNCGS * temperature) / (2.0 * PROTON_MASS), 0.5) / 1.0e5; /* in km s^-1 */
      b_tot   = pow((b_therm * b_therm) + (b_5 * b_5), 0.5);
      S_H2    = (((1.0 - omega) / pow(1.0 + (x / b_tot), alpha)) * exp(-5.0e-7 * (1.0 + x))) +
             ((omega / pow(1.0 + x, 0.5)) * exp(-8.5e-4 * pow(1.0 + x, 0.5)));
    }

  /* The rate is 7.5e-11 * (n / 2.256e-4) * Sd * SH2, where n is the photon number density in the band 12.24 eV to 13.51 eV
   * This rate has been normalised to an average of the rates from Cloudy's big H2 model for densities 1cc to 10^4 cc in
   * the presence of the Black (1987) ISRF. We then assume other spectra scale with n. */
  rate = 0.0;
  for(l = 0; l < myGlobalVars->N_spectra; l++)
    rate += 3.324468e-7 * myGasVars->H2_dissocJ[l] * photon_density[l] * LIGHTSPEED * exp(-3.74 * extinction) * S_H2;
  return rate;
}

/* The following gives the rate coefficient for the formation
 * of H2 on dust grains, as given by CT02. Note that this is
 * the default used by Cloudy */
double k60(double T, double Tgr, double xHI, double dust_ratio)
{
  double intermediate    = 35.399494936611667; /* The (1 + sqrt((Ehc - Es) / (Ehp - Es)))**2 term */
  double beta_over_alpha = 0.25 * intermediate * exp(-200.0 / Tgr);
  double beta            = 3e12 * exp(-320.0 / Tgr);
  double xi              = 1.0 / (1.0 + ((1.3e13 * exp(-1.5e4 / Tgr) * intermediate) / 2e-10));
  double epsilon         = xi / (1.0 + ((0.005 * 1e-10) / (2 * beta)) + beta_over_alpha);
  double SH              = 1.0 / (1.0 + 0.04 * pow((T + Tgr), 0.5) + 0.002 * T + 8e-6 * pow(T, 2));
  double vH              = pow(((8.0 * 1.3806488e-16 * T) / (3.141592653589793 * 1.672621777e-24)), 0.5);

  /* Note that, with the way our reaction structures are set up,
   * k60 will be multiplied by nHI * nHI, when in fact we want it
   * to be multiplied by nHI * nd = nHI * nHtot * DUST_CROSS_SECTION.
   * We therefore need to divide by xHI here to correct for this.
   * If xHI is zero then the reaction rate will be zero anyway - need
   * to explicitly include this to prevent division by zero. */
  if(xHI > 0.0)
    return 0.5 * vH * DUST_CROSS_SECTION * dust_ratio * epsilon * SH / xHI;
  else
    return 0.0;
}

double R70(double zeta_H) { return 2.0 * zeta_H; }

double k84(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction84, T_index, dT)); }

double k85(double T) { return 3.7e-14 * exp(-35.0 / T); }

double k86(double T) { return 7.2e-15; }

double k88(double T, double nHtotal, double xH, double xH2, double xHe)
{
  double ncr_H   = pow(10.0, (3.0 - 0.416 * log10(T / 1.0e4) - 0.327 * pow(log10(T / 1.0e4), 2)));
  double ncr_H2  = pow(10.0, (4.845 - 1.3 * log10(T / 1.0e4) + 1.62 * pow(log10(T / 1.0e4), 2)));
  double ncr_He  = pow(10.0, (5.0792 * (1.0 - 1.23e-5 * (T - 2000.0))));
  double ncr_inv = (xH / ncr_H) + (2.0 * xH2 / ncr_H2) + (xHe / ncr_He);
  double ncr;
  double log_kLTE, log_k0, intermediate;

  log_kLTE = -27.029 + 3.801 * log10(T) - (29487.0 / T);
  log_k0   = -2.729 - 1.75 * log10(T) - (23474.0 / T);

  if(abs(ncr_inv * nHtotal) < 1.0e-50)
    intermediate = log_k0;
  else
    {
      ncr          = 1.0 / ncr_inv;
      intermediate = ((nHtotal / ncr) / (1.0 + (nHtotal / ncr))) * log_kLTE + (1.0 / (1.0 + (nHtotal / ncr))) * log_k0;
    }

  return pow(10.0, intermediate);
}

double R89(double zeta_H) { return 0.22 * zeta_H; }

double k100(void) { return 1.7e-9; }

double k108(void) { return 4.9e-9; }

double R111(double dust_G, double extinction)
{
  /* Rate from NL99 */
  /* gamma_d from van Dishoeck et al. (2006) */
  return dust_G * 1.5e-10 * exp(-3.52 * extinction);
}

double k115(void) { return 1.0e-16; }

double k118(void) { return 3.8e-10; }

double k119(void) { return 4.0e-10; }

double k123(void) { return 6.59e-11; }

double k125(void) { return 6.64e-11; }

double k126(void) { return 1.33e-10; }

double k127(void) { return 8.0e-11; }

double k132(void) { return 1.0e-10; }

double k140(void) { return 2.4e-9; }

double k141(void) { return 2.0e-9; }

double k143(void) { return 7.5e-10; }

double k144(void) { return 1.2e-9; }

double k145(void) { return 3.5e-10; }

double k146(void) { return 1.4e-9; }

double k148(void) { return 1.6e-9; }

double k149(void) { return 7.5e-10; }

double k151(void) { return 4.0e-10; }

double k152(void) { return 4.8e-10; }

double k153(void) { return 1.7e-9; }

double k154(void) { return 1.5e-9; }

double k155(void) { return 8.4e-10; }

double k156(void) { return 1.3e-9; }

double k157(int T_index, double dT) { return pow(10.0, interpol_1d_mydbl(chimesRateTables.RatesTables->reaction157, T_index, dT)); }

double k158(void) { return 1.01e-9; }

double k159(void) { return 6.4e-10; }

double k160(void) { return 5.9e-9; }

double k161(void) { return 9.0e-10; }

double k162(void) { return 1.8e-9; }

double k163(void) { return 1.0e-11; }

double k164(void) { return 3.8e-10; }

double k165(void) { return 6.2e-10; }

double k166(void) { return 9.1e-10; }

double k167(void) { return 5.2e-11; }

double k168(void) { return 2.7e-11; }

double k169(void) { return 1.1e-9; }

double k170(void) { return 2.5e-9; }

double k171(void) { return 1.9e-9; }

double k172(void) { return 1.4e-9; }

double k173(void) { return 7.5e-10; }

double k174(void) { return 1.6e-9; }

double k175(void) { return 2.1e-9; }

double k176(void) { return 1.1e-9; }

double k177(void) { return 6.9e-9; }

double k178(void) { return 2.04e-10; }

double k179(void) { return 2.86e-10; }

double k180(void) { return 6.05e-11; }

double k181(void) { return 2.0e-9; }

double k182(void) { return 3.3e-11; }

double k183(void) { return 1.1e-9; }

double k184(void) { return 5.2e-11; }

double k185(void) { return 7.5e-10; }

double k208(void) { return 1.0e-9; }

double k209(void) { return 1.0e-9; }

double k210(void) { return 1.0e-10; }

double k211(void) { return 5.0e-10; }

double k212(void) { return 1.0e-13; }

double k213(void) { return 5.0e-10; }

double k214(void) { return 5.0e-10; }

double k215(void) { return 7.0e-10; }

double k216(void) { return 5.0e-10; }

double k217(void) { return 2.25e-15; }

double k218(void) { return 1.0e-17; }

double k219(void) { return 1.0e-17; }

double k224(void) { return 1.5e-15; }

double R239(double Go, double extinction) { return 4.9e-13 * Go * exp(-1.8 * extinction) / 1.7; }

double R240(double Go, double extinction) { return 4.9e-13 * Go * exp(-2.3 * extinction) / 1.7; }

double R241(double Go, double extinction) { return 2.4e-7 * Go * exp(-0.9 * extinction) / 1.7; }

double R242(double Go, double extinction) { return 9.2e-10 * Go * exp(-1.72 * extinction) / 1.7; }

double R243(double Go, double extinction) { return 7.7e-10 * Go * exp(-3.2 * extinction) / 1.7; }

double R244(double Go, double extinction) { return 3.3e-10 * Go * exp(-2.94 * extinction) / 1.7; }

double R245(double Go, double extinction) { return 5.8e-10 * Go * exp(-2.21 * extinction) / 1.7; }

double R246(double Go, double extinction) { return 5.9e-10 * Go * exp(-2.3 * extinction) / 1.7; }

double R247(double Go, double extinction) { return 1.4e-10 * Go * exp(-1.7 * extinction) / 1.7; }

double R248(double Go, double extinction) { return 1.0e-9 * Go * exp(-1.7 * extinction) / 1.7; }

double R249(double Go, double extinction) { return 1.0e-9 * Go * exp(-1.7 * extinction) / 1.7; }

double R250(double Go, double extinction) { return 2.4e-10 * Go * exp(-2.57 * extinction) / 1.7; }

double R251(double Go, double extinction) { return 2.4e-7 * Go * exp(-0.5 * extinction) / 1.7; }

double R252(double Go, double extinction) { return 3.7e-10 * Go * exp(-2.24 * extinction) / 1.7; }

double R253(double Go, double extinction) { return 1.6e-12 * Go * exp(-3.1 * extinction) / 1.7; }

double R254(double Go, double extinction) { return 1.1e-11 * Go * exp(-3.5 * extinction) / 1.7; }

double R255(double Go, double extinction) { return 6.0e-10 * Go * exp(-2.2 * extinction) / 1.7; }

double R256(double Go, double extinction) { return 3.2e-11 * Go * exp(-3.9 * extinction) / 1.7; }

double R257(double Go, double extinction)
{
  if(extinction > 15)
    return 5.0e-11 * Go * exp(-2.8 * extinction) / 1.7;
  else
    return 5.0e-11 * Go * exp(-2.55 * extinction + 0.0165 * pow(extinction, 2.0)) / 1.7;
}

double R258(double Go, double extinction)
{
  if(extinction > 15)
    return 5.0e-11 * Go * exp(-2.8 * extinction) / 1.7;
  else
    return 5.0e-11 * Go * exp(-2.55 * extinction + 0.0165 * pow(extinction, 2.0)) / 1.7;
}

double R259(double Go, double extinction)
{
  if(extinction > 15)
    return 5.0e-11 * Go * exp(-2.8 * extinction) / 1.7;
  else
    return 5.0e-11 * Go * exp(-2.55 * extinction + 0.0165 * pow(extinction, 2.0)) / 1.7;
}

double R260(double Go, double extinction)
{
  if(extinction > 15)
    return 1.5e-10 * Go * exp(-2.8 * extinction) / 1.7;
  else
    return 1.5e-10 * Go * exp(-2.55 * extinction + 0.0165 * pow(extinction, 2.0)) / 1.7;
}

double R261(double Go, double extinction)
{
  if(extinction > 15)
    return 2.5e-11 * Go * exp(-2.8 * extinction) / 1.7;
  else
    return 2.5e-11 * Go * exp(-2.55 * extinction + 0.0165 * pow(extinction, 2.0)) / 1.7;
}

double R262(double Go, double extinction)
{
  if(extinction > 15)
    return 2.5e-11 * Go * exp(-2.8 * extinction) / 1.7;
  else
    return 2.5e-11 * Go * exp(-2.55 * extinction + 0.0165 * pow(extinction, 2.0)) / 1.7;
}

double R263(double Go, double extinction)
{
  if(extinction > 15)
    return 7.5e-12 * Go * exp(-2.8 * extinction) / 1.7;
  else
    return 7.5e-12 * Go * exp(-2.55 * extinction + 0.0165 * pow(extinction, 2.0)) / 1.7;
}

double R264(double Go, double extinction)
{
  if(extinction > 15)
    return 2.5e-11 * Go * exp(-2.8 * extinction) / 1.7;
  else
    return 2.5e-11 * Go * exp(-2.55 * extinction + 0.0165 * pow(extinction, 2.0)) / 1.7;
}

double R265(double Go, double extinction) { return 7.6e-10 * Go * exp(-3.87 * extinction) / 1.7; }

double R266(double Go, double extinction) { return 7.9e-10 * Go * exp(-2.13 * extinction) / 1.7; }

double R267(double Go, double NH2, double NCO, double extinction, struct globalVariables *myGlobalVars)
{
  double S_CO, dN_CO, dN_H2;
  int NH2_index, NCO_index;

  if(NCO < pow(10.0, chimesRateTables.NonEqIon->COself_shielding_N[0]))
    {
      NCO_index = 0;
      dN_CO     = 0.0;
    }
  else if(NCO < pow(10.0, chimesRateTables.NonEqIon->COself_shielding_N[chimesRateTables.NonEqIon->shielding_dimensions[0] - 1]))
    {
      get_index_1d_irregular(chimesRateTables.NonEqIon->COself_shielding_N, chimesRateTables.NonEqIon->shielding_dimensions[0],
                             log10(NCO), &NCO_index, &dN_CO);
    }
  else
    {
      NCO_index = chimesRateTables.NonEqIon->shielding_dimensions[0] - 2;
      dN_CO     = 1.0;
    }

  if(NH2 < pow(10.0, chimesRateTables.NonEqIon->H2CO_shielding_N[0]))
    {
      NH2_index = 0;
      dN_H2     = 0.0;
    }
  else if(NH2 < pow(10.0, chimesRateTables.NonEqIon->H2CO_shielding_N[chimesRateTables.NonEqIon->shielding_dimensions[1] - 1]))
    {
      get_index_1d_irregular(chimesRateTables.NonEqIon->H2CO_shielding_N, chimesRateTables.NonEqIon->shielding_dimensions[1],
                             log10(NH2), &NH2_index, &dN_H2);
    }
  else
    {
      NH2_index = chimesRateTables.NonEqIon->shielding_dimensions[1] - 2;
      dN_H2     = 1.0;
    }

  S_CO = pow(10.0, interpol_2d_mydbl(chimesRateTables.NonEqIon->CO_shielding_S, NCO_index, NH2_index, dN_CO, dN_H2));

  return 2.6e-10 * S_CO * Go * exp(-3.53 * extinction) / 1.7;
}

double R268(double zeta_H) { return 0.037 * zeta_H; }

double R269(double zeta_H) { return 6.5e-4 * zeta_H; }

double R270(double zeta_H) { return 6.5 * zeta_H; }

double R272(double zeta_H) { return 4.0e3 * zeta_H; }

double R273(double zeta_H) { return 9.6e2 * zeta_H; }

double R274(double zeta_H) { return 2.7e3 * zeta_H; }

double R275(double zeta_H) { return 2.7e3 * zeta_H; }

double R276(double zeta_H) { return 1.3e3 * zeta_H; }

double R277(double zeta_H) { return 2.8e3 * zeta_H; }

double R278(double zeta_H) { return 5.3e3 * zeta_H; }

double R279(double zeta_H) { return 4.1e3 * zeta_H; }

double R280(double zeta_H) { return 6.4e2 * zeta_H; }

double R281(double zeta_H, double T, double xH2, double xCO)
{
  if(xCO > 0)
    return 0.21 * pow(T, 0.5) * xH2 * pow(xCO, -0.5) * zeta_H;
  else
    return 0.0;
}

double k290(double T, double ne, double dust_G, double extinction, double dust_ratio)
{
  if(dust_G == 0.0 || ne == 0.0)
    return 2.51e-14 * dust_ratio;
  else
    {
      double psi = dust_G * exp(-extinction * G0_GAMMA) * pow(T, 0.5) / ne;
      if(psi > 1.0e-90)
        return 2.51e-14 * dust_ratio /
               (1.0 + 8.116e-8 * pow(psi, 1.864) * (1.0 + 6.17e4 * pow(T, 2.169e-6) * pow(psi, (-0.9605 - 7.232e-5 * log(T)))));
      else
        return 2.51e-14 * dust_ratio;
    }
}

double k291(double T, double ne, double dust_G, double extinction, double dust_ratio)
{
  if(dust_G == 0.0 || ne == 0.0)
    return 3.064e-14 * dust_ratio;
  else
    {
      double psi = dust_G * exp(-extinction * G0_GAMMA) * pow(T, 0.5) / ne;
      if(psi > 1.0e-90)
        return 3.064e-14 * dust_ratio /
               (1.0 + 7.769e-8 * pow(psi, 1.319) * (1.0 + 1.087e2 * pow(T, 3.475e-1) * pow(psi, (-0.479 - 4.689e-2 * log(T)))));
      else
        return 3.064e-14 * dust_ratio;
    }
}

double k292(double T, double ne, double dust_G, double extinction, double dust_ratio)
{
  if(dust_G == 0.0 || ne == 0.0)
    return 1.636e-14 * dust_ratio;
  else
    {
      double psi = dust_G * exp(-extinction * G0_GAMMA) * pow(T, 0.5) / ne;
      if(psi > 1.0e-90)
        return 1.636e-14 * dust_ratio /
               (1.0 + 8.208e-9 * pow(psi, 2.289) * (1.0 + 1.254e5 * pow(T, 1.349e-9) * pow(psi, (-1.1506 - 7.204e-4 * log(T)))));
      else
        return 1.636e-14 * dust_ratio;
    }
}

double k293(double T, double ne, double dust_G, double extinction, double dust_ratio)
{
  if(dust_G == 0.0 || ne == 0.0)
    return 8.27e-14 * dust_ratio;
  else
    {
      double psi = dust_G * exp(-extinction * G0_GAMMA) * pow(T, 0.5) / ne;
      if(psi > 1.0e-90)
        return 8.27e-14 * dust_ratio /
               (1.0 + 2.051e-4 * pow(psi, 1.252) * (1.0 + 1.59e2 * pow(T, 6.072e-2) * pow(psi, (-0.598 - 4.497e-7 * log(T)))));
      else
        return 8.27e-14 * dust_ratio;
    }
}

double cr_secondary_ionisation(double xHII, int ion_index)
{
  /* This function interpolates the tabulated fraction of
   * secondary to primary CR ionisations, as given by
   * Furlanetto & Stoever (2010), as a function of xHII. */
  int x_ion_index;
  double dx;

  if(xHII <= pow(10.0, chimesRateTables.NonEqIon->x_ion_fraction[0]))
    {
      if(ion_index == HI)
        return pow(10.0, chimesRateTables.NonEqIon->n_ion_HI[0]);
      else if(ion_index == HeI)
        return pow(10.0, chimesRateTables.NonEqIon->n_ion_HeI[0]);
      else
        return 0.0;
    }
  else if(xHII >= pow(10.0, chimesRateTables.NonEqIon->x_ion_fraction[chimesRateTables.NonEqIon->secondary_ionisation_dims[0] - 1]))
    {
      if(ion_index == HI)
        return pow(10.0, chimesRateTables.NonEqIon->n_ion_HI[chimesRateTables.NonEqIon->secondary_ionisation_dims[0] - 1]);
      else if(ion_index == HeI)
        return pow(10.0, chimesRateTables.NonEqIon->n_ion_HeI[chimesRateTables.NonEqIon->secondary_ionisation_dims[0] - 1]);
      else
        return 0.0;
    }
  else
    {
      get_index_1d_irregular(chimesRateTables.NonEqIon->x_ion_fraction, chimesRateTables.NonEqIon->secondary_ionisation_dims[0],
                             log10(xHII), &x_ion_index, &dx);
      if(ion_index == HI)
        return pow(10.0, interpol_1d_mydbl(chimesRateTables.NonEqIon->n_ion_HI, x_ion_index, dx));
      else if(ion_index == HeI)
        return pow(10.0, interpol_1d_mydbl(chimesRateTables.NonEqIon->n_ion_HeI, x_ion_index, dx));
      else
        return 0.0;
    }
}

/* Additional Mg CT reactions */

double k294(void) { return 1.1e-9; }

double k295(void) { return 1.2e-9; }

double k296(void) { return 2.9e-9; }

double k297(void) { return 2.8e-10; }

double k303(void) { return 7.5e-10; }

double k304(void) { return 7.5e-10; }
