/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \date        22/10/2018
 * \author      Ondrej Jaura
 * \brief       Direct chemistry for the FiBY project
 * \details     
 */

#include <stdio.h>

#include "../allvars.h"
#include "../proto.h"

#if(SX_CHEMISTRY == 4)

static double sxMeanPhotEnInv[SX_NFREQ];  // inverse of the mean photon energy in the frequency bins

/*!
 * \brief functions that returns Arepo image property value
 */
double sx_chem_image1(int index) { return SphP[index].MassFract[MFHP]; }
double sx_chem_image2(int index) { return SphP[index].Temperature; }
double sx_chem_image_all(int index, int prop) { return SphP[index].MassFract[prop]; }

//======================================================================================================//
//                                    Chemistry initialization                                          //
//======================================================================================================//

/*!
 * \brief cross-section for different transitions
 */
double computeCrossHI(double x, void *p)  // from Verner et al. 1996, [cm^-2]
{
  return (x < 1.0) ? 0
                   : 1.16032e43 * pow(31.6379628 * x - 1.0, 2.0) * pow(TFIH * x, -4.0185) *
                         pow(1.0 + 1.7107e-8 * pow(TFIH * x, 0.5), -2.963);
}
double computeCrossHeI(double x, void *p)  // from Verner et al. 1996, [cm^-2]
{
  return (x < 1.808192) ? 0
                        : 9.492e-16 * (4.15752 + pow(1.4434 - 0.999118 * x, 2.0)) *
                              pow(1.0 + 0.825067 * pow(4.5625 + pow(0.4434 - 0.999118 * x, 2.0), 0.25), -3.188) *
                              pow(4.5625 + pow(0.4434 - 0.999118 * x, 2.0), -1.953);
}
double computeCrossHeII(double x, void *p)  // from Verner et al. 1996, [cm^-2]
{
  return (x < 4.001849)
             ? 0
             : 7.63458e44 * pow(7.90582 * x - 1.0, 2) * pow(TFIH * x, -4.0185) * pow(1.0 + 8.55151e-9 * pow(TFIH * x, 0.5), -2.963);
}

/*!
 * \brief chemistry initialization at the beginning of the simulation
 */
void sx_chem_initialize(void)
{
  int f, s;
  double dFreqInv, energy1, energy2;

  // initialize frequency bins [Hz]
  All.sxFreqBins[0]  = 1.0;  // 13.6 eV - H+
  All.sxFreqBins[1]  = 0.5 * ((TFIHE / TFIH) - All.sxFreqBins[0]) + All.sxFreqBins[0];
  All.sxFreqBins[2]  = TFIHE / TFIH;  // 24.6 eV - He+
  All.sxFreqBins[3]  = 2.0 * All.sxFreqBins[2] - All.sxFreqBins[1];
  All.sxFreqBins[4]  = 2.0 * All.sxFreqBins[3] - All.sxFreqBins[2];
  All.sxFreqBins[5]  = 0.5 * ((TFIIHE / TFIH) - All.sxFreqBins[4]) + All.sxFreqBins[4];
  All.sxFreqBins[6]  = TFIIHE / TFIH;  // 54.4 eV - He++
  All.sxFreqBins[7]  = 2.0 * All.sxFreqBins[6] - All.sxFreqBins[5];
  All.sxFreqBins[8]  = 2.0 * All.sxFreqBins[7] - All.sxFreqBins[6];
  All.sxFreqBins[9]  = 2.0 * All.sxFreqBins[8] - All.sxFreqBins[7];
  All.sxFreqBins[10] = 10.;  // 136  eV - infinity

  // Initialize integrator
  double error, intVal;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
  gsl_function F;

  for(f = 0; f < SX_NFREQ; f++)
    {
      // Calculate mean cross-sections in the bin using integrals
      dFreqInv   = 1 / (All.sxFreqBins[f + 1] - All.sxFreqBins[f]);
      F.function = &computeCrossHI;
      gsl_integration_qags(&F, All.sxFreqBins[f], All.sxFreqBins[f + 1], 0, 1e-7, 1000, w, &intVal, &error);
      All.sxSigma[CSH][f] = intVal * dFreqInv;
      if(f < 2)
        All.sxSigma[CSHE][f] = 0;
      else
        {
          F.function = &computeCrossHeI;
          gsl_integration_qags(&F, All.sxFreqBins[f], All.sxFreqBins[f + 1], 0, 1e-7, 1000, w, &intVal, &error);
          All.sxSigma[CSHE][f] = intVal * dFreqInv;
        }
      if(f < 6)
        All.sxSigma[CSHEP][f] = 0;
      else
        {
          F.function = &computeCrossHeII;
          gsl_integration_qags(&F, All.sxFreqBins[f], All.sxFreqBins[f + 1], 0, 1e-7, 1000, w, &intVal, &error);
          All.sxSigma[CSHEP][f] = intVal * dFreqInv;
        }

      // Calculate excess energies from the ionization (mean of bins)
      // for H
      energy1              = PLANCK * TFIH * (All.sxFreqBins[f] - 1.0);
      energy2              = PLANCK * TFIH * (All.sxFreqBins[f + 1] - 1.0);
      All.sxEnergy[EEH][f] = 0.5 * (energy1 + energy2);
      // for He
      if(All.sxFreqBins[f] >= TFIHE / TFIH)
        {
          energy1               = PLANCK * TFIH * (All.sxFreqBins[f] - TFIHE / TFIH);
          energy2               = PLANCK * TFIH * (All.sxFreqBins[f + 1] - TFIHE / TFIH);
          All.sxEnergy[EEHE][f] = 0.5 * (energy1 + energy2);
        }
      else
        All.sxEnergy[EEHE][f] = 0.0;
      // for He+
      if(All.sxFreqBins[f] >= TFIIHE / TFIH)
        {
          energy1                = PLANCK * TFIH * (All.sxFreqBins[f] - TFIIHE / TFIH);
          energy2                = PLANCK * TFIH * (All.sxFreqBins[f + 1] - TFIIHE / TFIH);
          All.sxEnergy[EEHEP][f] = 0.5 * (energy1 + energy2);
        }
      else
        All.sxEnergy[EEHEP][f] = 0.0;

      // Calculate mean photon energies in frequency bins in [erg]
      sxMeanPhotEnInv[f] = 1. / ((All.sxFreqBins[f] + 0.5 * (All.sxFreqBins[f + 1] - All.sxFreqBins[f])) * TFIH * PLANCK);
    }

  // Close integrator
  gsl_integration_workspace_free(w);

#ifdef SX_DISPLAY_STATS
  mpi_printf("SX-INIT: Frequency bins freqH+ [ ");
  for(f = 0; f < SX_NFREQ + 1; f++)
    mpi_printf("%.2e ", All.sxFreqBins[f]);
  mpi_printf("]\n                       eV [ ");
  for(f = 0; f < SX_NFREQ + 1; f++)
    mpi_printf("%.1f ", All.sxFreqBins[f] * TFIH / ELECTRONVOLT_IN_HZ);
  mpi_printf("] \n");

  mpi_printf("SX-INIT: sigma HP   cm^-2 [ ");
  for(f = 0; f < SX_NFREQ; f++)
    mpi_printf("%.2e ", All.sxSigma[CSH][f]);
  mpi_printf("] \n");
  mpi_printf("SX-INIT: sigma HEP  cm^-2 [ ");
  for(f = 0; f < SX_NFREQ; f++)
    mpi_printf("%.2e ", All.sxSigma[CSHE][f]);
  mpi_printf("] \n");
  mpi_printf("SX-INIT: sigma HEPP cm^-2 [ ");
  for(f = 0; f < SX_NFREQ; f++)
    mpi_printf("%.2e ", All.sxSigma[CSHEP][f]);
  mpi_printf("] \n");

  mpi_printf("SX-INIT: energy HP   erg [ ");
  for(f = 0; f < SX_NFREQ; f++)
    mpi_printf("%.2e ", All.sxEnergy[EEH][f]);
  mpi_printf("] \n");
  mpi_printf("SX-INIT: energy HEP  erg [ ");
  for(f = 0; f < SX_NFREQ; f++)
    mpi_printf("%.2e ", All.sxEnergy[EEHE][f]);
  mpi_printf("] \n");
  mpi_printf("SX-INIT: energy HEPP erg [ ");
  for(f = 0; f < SX_NFREQ; f++)
    mpi_printf("%.2e ", All.sxEnergy[EEHEP][f]);
  mpi_printf("] \n");
#endif
}

/*!
 * \brief initialize chemistry at the beginning of the run
 */
void sx_chem_start_run(void) {}

//======================================================================================================//
//                                     Rate equation                                                    //
//======================================================================================================//

/*
 * \brief change number fraction due to a reaction
 */
void sx_chem_change_fraction(struct sxSite *site, int from, int to, double value)
{
  // DEBUG: when subtracting, we need to prevent small negative values caused by precision rounding
  if(value > site->nFract[from])
    value = site->nFract[from];
  site->nFract[from] -= value;
  site->nFract[to] += value;
}

/*!
 * \brief solve the heating in the site
 */
void sx_chem_change_temperature(struct sxSite *site)
{
  double U, dU, mu;
  double nSpecH, nSpecHP, nSpecHE, nSpecHEP, nSpecHEPP;

  // calculate number of species
  nSpecH    = site->nNucleons * site->nFract[MFH];
  nSpecHP   = site->nNucleons * site->nFract[MFHP];
  nSpecHE   = site->nNucleons * site->nFract[MFHE];
  nSpecHEP  = site->nNucleons * site->nFract[MFHEP];
  nSpecHEPP = site->nNucleons * site->nFract[MFHEPP];

  // calculate the mean molecular weight
  mu = (nSpecH + nSpecHP + 4. * (nSpecHE + nSpecHEP + nSpecHEPP)) / (nSpecH + 2. * nSpecHP + nSpecHE + 2. * nSpecHEP + 3. * nSpecHEPP);

  // calculate the specific energy
  U = site->temp * (1.5 * BOLTZMANN) / (mu * PROTONMASS);

  // calculate the change in the energy
  dU = (site->heating - site->cooling) / (site->density * site->volume);

  // save the new temperature
  site->temp = (U + dU) * (mu * PROTONMASS) / (1.5 * BOLTZMANN);

#ifdef SX_CMB_TEMP_LIMIT
  // constrain the minimum temperature by CMB
  if(site->temp < sxRunData->tempCMB)
    site->temp = sxRunData->tempCMB;
#endif
}

/*
 * \brief recombination factor for different species
 */
double computeRecombHIICaseB(double T)  // from Hui & Gnedin 1997, [cm^3/s]
{
  return 2.753e-14 * pow(2.0 * 1.57807e5 / T, 1.5) / pow((1.0 + pow(0.729927007 * 1.57807e5 / T, 0.407)), 2.242);
}

double computeRecombHeIICaseB(double T)  // fit to Hummer & Storey 1998, [cm^3/s]
{
  return pow(10., -9.79239 - 0.687189 * log10(T));
}

double computeRecombHeIIICaseB(double T)  // from Hui & Gnedin 1997, [cm^3/s]
{
  return 5.506e-14 * pow(2.0 * 6.31515e5 / T, 1.5) / pow((1.0 + pow(0.729927007 * 6.31515e5 / T, 0.407)), 2.242);
}
double computeDielectronicRecombHeII(double T)  // from Cen 1992
{
  return 1.9e-3 * pow(T, -1.5) * exp(-4.7e5 / T) * (1.0 + 0.3 * exp(-9.4e4 / T));
}

/*
 * \brief compute various cooling coefficients
 */
double computeRecombCoolingHIICaseB(double T)  // from Hui & Gnedin 1997, [erg*cm^3/s]
{
  return 3.435e-30 * T * pow(2.0 * 1.57807e5 / T, 1.970) / pow((1.0 + pow(0.88888888889 * 1.57807e5 / T, 0.376)), 3.720);
}
double computeRecombCoolingHeIICaseB(double T)  // fit to Hummer & Storey 1998, [erg*cm^3/s]
{
  double logT = log10(T);
  return pow(10.0, -25.9067 + 0.500511 * logT + -0.0414826 * logT * logT);
}
double computeRecombCoolingHeIIICaseB(double T)  // from Hui & Gnedin 1997, [erg*cm^3/s]
{
  return 8.0 * 3.435e-30 * T * pow(2.0 * 6.31515e5 / T, 1.970) / pow((1.0 + pow(0.888888888889 * 6.31515e5 / T, 0.376)), 3.720);
}
double computeCollIonisationCoolingHI(double T)  // from Theuns et al. 1998
{
  return 2.54e-21 * sqrt(T) * exp(-1.578091e5 / T) * pow(1.0 + sqrt(T / 1e5), -1.0);
}
double computeCollIonisationCoolingHeI(double T)  // from Theuns et al. 1998
{
  return 1.88e-21 * sqrt(T) * exp(-2.853354e5 / T) * pow(1.0 + sqrt(T / 1e5), -1.0);
}
double computeCollIonisationCoolingHeII(double T)  // from Theuns et al. 1998
{
  return 9.90e-22 * sqrt(T) * exp(-6.31515e5 / T) * pow(1.0 + sqrt(T / 1e5), -1.0);
}
double computeCollExcitationCoolingHI(double T)  // from Cen 1992
{
  return 7.5e-19 * pow(1.0 + sqrt(T / 1.e5), -1.0) * exp(-1.18348e5 / T);
}
double computeCollExcitationCoolingHeI(double T)  // from Cen 1992
{
  return 9.1e-27 * pow(T, -0.1687) * pow(1.0 + sqrt(T / 1.e5), -1.0) * exp(-1.3179e4 / T);
}
double computeCollExcitationCoolingHeII(double T)  // from Cen 1992
{
  return 5.43e-17 * pow(T, -0.379) * pow(1.0 + sqrt(T / 1.e5), -1.0) * exp(-4.73638e5 / T);
}
double computeFreeFreeCooling(double T)  // from Theuns et al. 1998
{
  double g_ff = 1.1 + 0.34 * exp(-0.33333333333 * pow(5.5 - log(T) / log(10), 2.0));
  return 1.42e-27 * g_ff * sqrt(T);
}
double computeDielectronicRecombCoolingHeII(double T)  // from Cen 1992
{
  return 1.24e-13 * pow(T, -1.5) * exp(-4.7e5 / T) * (1.0 + 0.3 * exp(-9.4e4 / T));
}

/*
 * \brief Calculate gas recombination
 */
void sx_chem_solve_recombination(struct sxSite *site)
{
  double ne, tmp, dt;
  double nSpecH, nSpecHP, nSpecHE, nSpecHEP, nSpecHEPP, nSpec;
  double nRecombHP, nRecombHEP, nRecombHEPP;

  // calculate time of the subcycle
  dt = site->dt * site->scFactor;

  // calculate number of available species
  nSpecH    = site->nNucleons * site->nFract[MFH];
  nSpecHP   = site->nNucleons * site->nFract[MFHP];
  nSpecHE   = site->nNucleons * site->nFract[MFHE];
  nSpecHEP  = site->nNucleons * site->nFract[MFHEP];
  nSpecHEPP = site->nNucleons * site->nFract[MFHEPP];

  // time * free electron density
  ne  = (nSpecHP + nSpecHEP + 2.0 * nSpecHEPP) / site->volume;
  tmp = dt * ne;

  // sum up all cooling rates
  site->cooling = 0.0;

  site->cooling += nSpecH * (computeCollIonisationCoolingHI(site->temp) + computeCollExcitationCoolingHI(site->temp));

  site->cooling += nSpecHP * (computeRecombCoolingHIICaseB(site->temp) + computeFreeFreeCooling(site->temp));

  site->cooling += nSpecHE * (computeCollIonisationCoolingHeI(site->temp));

  site->cooling += nSpecHEP * (computeCollIonisationCoolingHeII(site->temp) + computeCollExcitationCoolingHeII(site->temp) +
                               computeRecombCoolingHeIICaseB(site->temp) + computeDielectronicRecombCoolingHeII(site->temp) +
				 computeFreeFreeCooling(site->temp));

  site->cooling += nSpecHEPP * (4.0 * computeFreeFreeCooling(site->temp) + computeRecombCoolingHeIIICaseB(site->temp));

  // multiply by electron density and time [s cm^-3]
  site->cooling *= tmp;

  // HeI excitation cooling depends on the square of the electron number density, so need to add it separately
  // NOTE: the use of nSpecHEP here is intentional and correct - see e.g. Cen (1992)
  site->cooling += computeCollExcitationCoolingHeI(site->temp) * nSpecHEP * ne * tmp;

  tmp = nSpecHP + nSpecHEP + nSpecHEPP;
  if(tmp <= 0)
    return;  // return if there are no species for recombination

  // time * free electron density
  tmp = dt * (nSpecHP + nSpecHEP + 2.0 * nSpecHEPP) / site->volume;

  // calculate number of recombinations
  nRecombHP   = nSpecHP * tmp * computeRecombHIICaseB(site->temp);
  nRecombHEP  = nSpecHEP * tmp * (computeRecombHeIICaseB(site->temp) + computeDielectronicRecombHeII(site->temp));
  nRecombHEPP = nSpecHEPP * tmp * computeRecombHeIIICaseB(site->temp);

  // restrict recombinations by the number of available species
  if(nRecombHP > nSpecHP)
    nRecombHP = nSpecHP;
  if(nRecombHEP > nSpecHEP)
    nRecombHEP = nSpecHEP;
  if(nRecombHEPP > nSpecHEPP)
    nRecombHEPP = nSpecHEPP;

  // update mass fractions
  sx_chem_change_fraction(site, MFHP, MFH, nRecombHP / site->nNucleons);
  sx_chem_change_fraction(site, MFHEP, MFHE, nRecombHEP / site->nNucleons);
  sx_chem_change_fraction(site, MFHEPP, MFHEP, nRecombHEPP / site->nNucleons);
}

/*
 * \brief Collisional ionisation coefficients
 */
double computeCollIonisationHI(double T)  // from Theuns et al. 1998
{
  return 1.17e-10 * sqrt(T) * exp(-1.578091e5 / T) * pow(1.0 + sqrt(T / 1.e5), -1.0);
}
double computeCollIonisationHeI(double T)  // from Theuns et al. 1998
{
  return 4.76e-11 * sqrt(T) * exp(-2.853354e5 / T) * pow(1.0 + sqrt(T / 1e5), -1);
}
double computeCollIonisationHeII(double T)  // from Theuns et al. 1998
{
  return 1.14e-11 * sqrt(T) * exp(-6.31515e5 / T) * pow(1.0 + sqrt(T / 1e5), -1);
}

/*
 * \brief Solve frequencies between 13.6+ eV
 */
void sx_chem_solve_ionization(struct sxSite *site)
{
  int f;
  double sigma, sigmaH, sigmaHE, sigmaHEP;
  double nPhot, nPhotH, nPhotHE, nPhotHEP;
  double nSpec, nSpecH, nSpecHP, nSpecHE, nSpecHEP, nSpecHEPP;
  double nColl, nCollH, nCollHE, nCollHEP;
  double nReactH, nReactHE, nReactHEP;
  double tmp;

  // calculate the time of the subcycle
  double dt = site->dt * site->scFactor;

  // loop through all frequencies
  for(f = SX_NFREQ - 1; f > -1; f--)
    // for(f=0; f<SX_NFREQ; f++)
    {
      // calculate number of available species
      nSpecH    = site->nNucleons * site->nFract[MFH];
      nSpecHP   = site->nNucleons * site->nFract[MFHP];
      nSpecHE   = site->nNucleons * site->nFract[MFHE];
      nSpecHEP  = site->nNucleons * site->nFract[MFHEP];
      nSpecHEPP = site->nNucleons * site->nFract[MFHEPP];
      nSpec     = nSpecH + nSpecHE + nSpecHEP;
      if(nSpec == 0)
        continue;  // skip this frequency if there are no available species for a reaction

      // calculate sigma parameters
      sigmaH   = All.sxSigma[CSH][f] * site->nFract[MFH];
      sigmaHE  = All.sxSigma[CSHE][f] * site->nFract[MFHE];
      sigmaHEP = All.sxSigma[CSHEP][f] * site->nFract[MFHEP];
      sigma    = sigmaH + sigmaHE + sigmaHEP;
      if(sigma == 0)
        continue;  // skip this frequency if there is not a cross-section for the reaction

      // caclulate number of photons that ionize
      nPhot    = site->nPhotons[f] * site->scFactor * (1.0 - exp(-1.0 * site->numdr * sigma));
      nPhotH   = nPhot * sigmaH / sigma;
      nPhotHE  = nPhot * sigmaHE / sigma;
      nPhotHEP = nPhot * sigmaHEP / sigma;

      // restrict
      if(nPhotH > nSpecH)
        nPhotH = nSpecH;
      if(nPhotHE > nSpecHE)
        nPhotHE = nSpecHE;
      if(nPhotHEP > nSpecHEP)
        nPhotHEP = nSpecHEP;
      nPhot = nPhotH + nPhotHE + nPhotHEP;

      // add excess energy to the total heating
      site->heating += nPhotH * All.sxEnergy[EEH][f];
      site->heating += nPhotHE * All.sxEnergy[EEHE][f];
      site->heating += nPhotHEP * All.sxEnergy[EEHEP][f];

      // calculate collisional ionisation
      tmp      = dt * (nSpecHP + nSpecHEP + 2.0 * nSpecHEPP) / site->volume;
      nCollH   = tmp * nSpecH * computeCollIonisationHI(site->temp);
      nCollHE  = tmp * nSpecHE * computeCollIonisationHeI(site->temp);
      nCollHEP = tmp * nSpecHEP * computeCollIonisationHeII(site->temp);

      // calculate total number of ionisation reactions
      nReactH   = (nPhotH + nCollH > nSpecH) ? nSpecH : nPhotH + nCollH;
      nReactHE  = (nPhotHE + nCollHE > nSpecHE) ? nSpecHE : nPhotHE + nCollHE;
      nReactHEP = (nPhotHEP + nCollHEP > nSpecHEP) ? nSpecHEP : nPhotHEP + nCollHEP;

      // update chemical abundances
      sx_chem_change_fraction(site, MFH, MFHP, nReactH / site->nNucleons);
      sx_chem_change_fraction(site, MFHE, MFHEP, nReactHE / site->nNucleons);
      sx_chem_change_fraction(site, MFHEP, MFHEPP, nReactHEP / site->nNucleons);

      // update number of photons  (DEBUG: the following line prevents negative values caused by rounding errors)
      site->nPhotons[f] = (site->nPhotons[f] - nPhot < 0) ? 0 : site->nPhotons[f] - nPhot;
#ifdef SX_DISPLAY_STATS
      sxRunData->nAbsPhotons[f] += nPhot;
#endif

      // dump small number of photons
      if(site->nPhotons[f] < site->nNucleons * All.sxMinNumPhotons)
        {
#ifdef SX_DISPLAY_STATS
          sxRunData->nLostPhotons[f] += site->nPhotons[f];
#endif
          site->nPhotons[f] = 0.0;
        }
    }
}

/*!
 * \brief solve the rate equation
 */
void sx_chem_solve_site(struct sxSite *site)
{
  // initialize variables
  site->heating  = 0.0;  // [erg]
  site->cooling  = 0.0;  // [erg]
  site->scFactor = 1.0;

  // Solve recombination
  sx_chem_solve_recombination(site);

  // Solve ionisation
  sx_chem_solve_ionization(site);

  // Change temperature
  sx_chem_change_temperature(site);
}
/*
void sx_chem_solve_site( struct sxSite *site )
{
  double nfInitInv, nfBefore[3], nfRatio=0;
  double nSpec=0,nPhot=0, subTime=0;
  int i;

  // Initial guess for a scFactor
  for(i=0; i<SX_NFREQ; i++)
    nPhot += site->nPhotons[i];
  nSpec = site->nNucleons * ( site->nFract[MFH] + site->nFract[MFHP] + site->nFract[MFHE]
                              + site->nFract[MFHEP] + site->nFract[MFHEPP] );
  //site->scFactor = (nPhot>nSpec*0.05) ? 0.05*nSpec/nPhot : 1.0;
  site->scFactor = 1.0;

  nfInitInv = (site->nFract[MFH] + site->nFract[MFHE] + site->nFract[MFHEP]);
  nfInitInv = (nfInitInv>0.0) ? 1.0/nfInitInv : 0.0;

  // do Sub-Cycling
  do {
    // Constrain scFactor
    if (site->scFactor<1e-4)
      site->scFactor = 1e-4;
    if (subTime+site->scFactor>1.0)
      site->scFactor = 1.0-subTime;

    // Initialize some variables
    site->heating = 0.0;  // [erg]
    site->cooling = 0.0;  // [erg]
    nfBefore[0] = site->nFract[MFH];
    nfBefore[1] = site->nFract[MFHE];
    nfBefore[2] = site->nFract[MFHEP];

    // Solve recombination
    sx_chem_solve_recombination( site );

    // Solve ionisation
    sx_chem_solve_ionization( site );

    // Change temperature
    sx_chem_change_temperature( site );

    // Calculate the change ratio
    nfRatio = (fabs(site->nFract[MFH]-nfBefore[0]) +
               fabs(site->nFract[MFHE]-nfBefore[1]) +
               fabs(site->nFract[MFHEP]-nfBefore[2])) * nfInitInv;

    // Calculate a new scFactor
    subTime += site->scFactor;
    site->scFactor = (nfRatio>0) ? site->scFactor*0.05/nfRatio : 1e-4;
  } while(subTime<1.0);

}
*/

/*!
 * \brief initialize chemistry variables in the site structure
 */
void sx_chem_initialize_site(struct sxSite *site)
{
  int i;

  site->nFract[MFH]    = SphP[site->index].MassFract[MFH];            // fraction of H
  site->nFract[MFHP]   = SphP[site->index].MassFract[MFHP];           // fraction of H+
  site->nFract[MFHE]   = SphP[site->index].MassFract[MFHE] * 0.25;    // fraction of He
  site->nFract[MFHEP]  = SphP[site->index].MassFract[MFHEP] * 0.25;   // fraction of He+
  site->nFract[MFHEPP] = SphP[site->index].MassFract[MFHEPP] * 0.25;  // fraction of He++
  site->temp           = SphP[site->index].Temperature;               // [K]

  site->density   = SphP[site->index].Density * sxRunData->density2cgs;  // density [g/cm^3]
  site->volume    = SphP[site->index].Volume * sxRunData->volume2cgs;    // volume [cm^3]
  site->numdens   = site->density / PROTONMASS;                          // proton number density [1/cm^3]
  site->nNucleons = site->numdens * site->volume;                        // initial number of nucleon species
  site->numdr     = site->numdens * sxCell[site->index].MeanNeighDistance * sxRunData->length2cgs; // [1/cm^2]
  site->heating = 0;  // [erg]
  site->cooling = 0;  // [erg]
}

/*!
 * \brief save changes in the site at the end
 */
void sx_chem_terminate_site(struct sxSite *site)
{
  if(isnan(site->nFract[MFH]) || site->nFract[MFH] < 0)
    {
      printf("SX: ERROR - nFract[MFH] is NaN/negative for the following site:\n");
      sx_debug_particle(site->index, P[site->index].ID);
      sx_debug_site(site, P[site->index].ID);
      terminate("SX: ERROR - nFract[MFH] is NaN/negative");
    }

  // save new values
  SphP[site->index].MassFract[MFH]    = site->nFract[MFH];
  SphP[site->index].MassFract[MFHP]   = site->nFract[MFHP];
  SphP[site->index].MassFract[MFHE]   = site->nFract[MFHE] * 4;
  SphP[site->index].MassFract[MFHEP]  = site->nFract[MFHEP] * 4;
  SphP[site->index].MassFract[MFHEPP] = site->nFract[MFHEPP] * 4;
  SphP[site->index].Temperature       = site->temp;
}

void sx_chem_compute_ionisation_rates(double ionisationRate[SX_NFREQ], int index)
{
#if(SX_SOURCES == 4)

  int f, starIndex;

  starIndex = P[index].AuxDataID;

  for(f = 0; f < SX_NFREQ; f++)
    {
      // convert the spectral emission in [log(erg/s)] to [ph/s]
      ionisationRate[f] = pow(10, StarP[starIndex].SEDs[f]);  // * sxMeanPhotEnInv[f];
    }

#endif
}

void sx_chem_end_run(void)
{
  /* needs to be defined for use in sx_end_run() */
}

void sx_chem_finalize_site(struct sxSite *site)
{
  /* needs to be defined for use in sx_evolve_sites() */
}


#endif
