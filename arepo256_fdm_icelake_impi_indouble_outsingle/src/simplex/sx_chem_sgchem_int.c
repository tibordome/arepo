/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \date        22/10/2018
 * \author      Ondrej Jaura
 * \brief       Calculation of various rates and integrals for the SGChem
 * \details     
 */

#include <gsl/gsl_integration.h>
#include <stdio.h>

#include "../allvars.h"
#include "../proto.h"

#if(SX_CHEMISTRY == 3)

// piecewise Cross-Section and energies for Ionization of H2 ( Liu & Shemansky 2012 )
#define NPW 11  // Number of PieceWise Intervals for H2

// number of elements = NPW
static double sxSigmaPW[] = {0.09e-18, 1.15e-18, 3.0e-18, 5.0e-18, 6.75e-18, 8.0e-18, 9.0e-18, 9.5e-18, 9.8e-18, 10.1e-18, 9.8e-18};
// number of elements = NPW + 1
static double sxEnergyPW[] = {3.675345566117e+15, 3.735795328717e+15, 3.796245091318e+15, 3.856694853919e+15,
                              3.917144616519e+15, 3.9655044266e+15,   4.0259541892e+15,   4.074313999281e+15,
                              4.110583856841e+15, 4.158943666922e+15, 4.267753239603e+15, 4.376562812284e+15};

/*
 Implementation of integrals from "Baczynski et al 2015"

 Some intermediate calculations of values used bellow

 sigmaH_136 = 6.30e-18 (Osterbrock 1989 pg 14.)
 nuH_136 = 3.28798e15 Hz
 sigmaH_136 * (nuH_136)^3 = 2.23938132e29

 sigmaH2_181 = 9.75e-18 (Liu & Shemansky 2012)
 nuH2_181 = 4.3765628e+15 Hz
 sigmaH2_181 * (nuH2_181)^3 = 8.17342552e29
*/

/***********************************************************************************************
 *                                        Cross Sections                                       *
 ***********************************************************************************************/

double sx_topCrossSectionH( double freq, void * p )  // 13.6+ eV -- Osterbrock 1989 pg 14.
{
  double T = *(double *)p;
  return 2.23938132e29 / (freq * (exp(PLANCK * freq / (BOLTZMANN * T)) - 1.0));
}

double sx_topCrossSectionH2( double freq, void * p ) // 18.1+ eV -- Liu & Shemansky 2012
{
  double T = *(double *)p;
  return 8.17342552e29 / (freq * (exp(PLANCK * freq / (BOLTZMANN * T)) - 1.0));
}

double sx_topCrossSectionHe(double freq, void *p)  // 24.6+ eV -- from Verner et al. 1996
{
  double T = *(double *)p;
  return 9.492e-16 * (4.15752 + pow(1.4434 - 3.038699e-16 * freq, 2)) *
         pow(1.0 + 0.825067 * pow(4.5625 + pow(0.4434 - 3.038699e-16 * freq, 2.0), 0.25), -3.188) *
         pow(4.5625 + pow(0.4434 - 3.038699e-16 * freq, 2.0), -1.953) * freq * freq / (exp(PLANCK * freq / (BOLTZMANN * T)) - 1.0);
}

double sx_bottomCrossSection(double freq, void *p)
{
  double T = *(double *)p;
  return freq * freq / (exp(PLANCK * freq / (BOLTZMANN * T)) - 1.0);
}

void sx_chem_compute_cross_sections(double T)
{
  int f;
  double intTop, intBottom, minFreq, maxFreq, error;

  // Initialize integrator
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
  gsl_function F;
  F.params = &T;

  // cross section for UV pumping
  All.sxSigma[S112DH2] = 2.47e-18;
  All.sxSigma[S136DH2] = 2.47e-18;  // DEBUG: shouldn't it be different for different freq?

  // Approximate cross section for absoprtion of LW photons by atomic H
  //
  // Based on the exponential part of the shielding function given in Wolcott-Green & Haiman (2011).
  // We do not attempt to account for the additional power-law dependence on the column density of H
  // and hence underestimate the shielding. However, this effect is only important when the H column
  // density is extremely large, and in these circumstances we're unlikely to be able to resolve the
  // exact location of the edge of the PDR in any event
  All.sxSigma[S112H] = 5.23e-25;

  // cross section for H (13.6 - 15.2 eV)
  minFreq = All.sxFreqBins[F136];
  maxFreq = All.sxFreqBins[F152];

  F.function = &sx_topCrossSectionH;
  gsl_integration_qags(&F, minFreq, maxFreq, 0, 1e-7, 1000, w, &intTop, &error);
  F.function = &sx_bottomCrossSection;
  gsl_integration_qags(&F, minFreq, maxFreq, 0, 1e-7, 1000, w, &intBottom, &error);
  All.sxSigma[S136IH] = intTop / intBottom;

  // cross sections for H and H2 (15.2-24.6 eV)
  minFreq = All.sxFreqBins[F152];
  maxFreq = All.sxFreqBins[F246];

  F.function = &sx_bottomCrossSection;
  gsl_integration_qags(&F, minFreq, maxFreq, 0, 1e-7, 1000, w, &intBottom, &error);

  F.function = &sx_topCrossSectionH;
  gsl_integration_qags(&F, minFreq, maxFreq, 0, 1e-7, 1000, w, &intTop, &error);
  All.sxSigma[S152IH] = intTop / intBottom;

  F.function = &sx_topCrossSectionH2;
  gsl_integration_qags(&F, sxEnergyPW[NPW], maxFreq, 0, 1e-7, 1000, w, &intTop, &error);  // 18.1+ eV
  double intBottomPW = 0.0;
  for(f = 0; f < NPW; f++)  // 15.2 - 18.1 eV
    {
      // here we can simply take advantage of bottomCrossSection multiplied by constant cross section
      F.function = &sx_bottomCrossSection;
      gsl_integration_qags(&F, sxEnergyPW[f], sxEnergyPW[f + 1], 0, 1e-7, 1000, w, &intBottomPW, &error);
      intTop += sxSigmaPW[f] * intBottomPW;
    }
  All.sxSigma[S152IH2] = intTop / intBottom;

  // 24.6+ eV He, H, H2 cross section
  minFreq = All.sxFreqBins[F246];
  maxFreq = All.sxFreqBins[SX_NFREQ];

  F.function = &sx_bottomCrossSection;
  gsl_integration_qags(&F, minFreq, maxFreq, 0, 1e-7, 1000, w, &intBottom, &error);

  F.function = &sx_topCrossSectionH;
  gsl_integration_qags(&F, minFreq, maxFreq, 0, 1e-7, 1000, w, &intTop, &error);
  All.sxSigma[S246IH] = intTop / intBottom;

  F.function = &sx_topCrossSectionH2;
  gsl_integration_qags(&F, minFreq, maxFreq, 0, 1e-7, 1000, w, &intTop, &error);  // 24.6 eV > 18.1 eV
  All.sxSigma[S246IH2] = intTop / intBottom;

  F.function = &sx_topCrossSectionHe;
  gsl_integration_qags(&F, minFreq, maxFreq, 0, 1e-7, 1000, w, &intTop, &error);
  All.sxSigma[S246IHE] = intTop / intBottom;

  // Close integrator
  gsl_integration_workspace_free(w);
}

/***********************************************************************************************
 *                                 Ionization Heating Energies                                 *
 ***********************************************************************************************/

double sx_topMeanEnergyH_136(double freq, void *p)  // 13.6 - 15.2 eV
{
  double T = *(double *)p;
  return 2.23938132e29 * (1.0 - TFIH / freq) / (exp(PLANCK * freq / (BOLTZMANN * T)) - 1.0);
}

double sx_topMeanEnergyH_152(double freq, void *p)  // 15.2+ eV
{
  double T = *(double *)p;
  return 2.23938132e29 * (1.0 - TFIH2 / freq) / (exp(PLANCK * freq / (BOLTZMANN * T)) - 1.0);
}

double sx_bottomMeanEnergyH(double freq, void *p)
{
  double T = *(double *)p;
  return 2.23938132e29 / (PLANCK * freq * (exp(PLANCK * freq / (BOLTZMANN * T)) - 1.0));
}

double sx_topMeanEnergyH2_152(double freq, void *p)  // 15.2 - 18.1 eV
{
  double T = *(double *)p;
  return freq * freq * freq * (1.0 - TFIH2 / freq) / (exp(PLANCK * freq / (BOLTZMANN * T)) - 1.0);
}

double sx_bottomMeanEnergyH2_152(double freq, void *p)  // 15.2 - 18.1 eV
{
  double T = *(double *)p;
  return freq * freq / (PLANCK * (exp(PLANCK * freq / (BOLTZMANN * T)) - 1.0));
}

double sx_topMeanEnergyH2_181(double freq, void *p)  // 18.1+ eV
{
  double T = *(double *)p;
  return 8.17342552e29 * (1.0 - TFIH2 / freq) / (exp(PLANCK * freq / (BOLTZMANN * T)) - 1.0);
}

double sx_bottomMeanEnergyH2_181(double freq, void *p)  // 18.1+ eV
{
  double T = *(double *)p;
  return 8.17342552e29 / (PLANCK * freq * (exp(PLANCK * freq / (BOLTZMANN * T)) - 1.0));
}

double sx_topMeanEnergyHe(double freq, void *p)  // 18.1+ eV
{
  double T = *(double *)p;
  return 9.492e-16 * (4.15752 + pow(1.4434 - 3.038699e-16 * freq, 2)) *
         pow(1.0 + 0.825067 * pow(4.5625 + pow(0.4434 - 3.038699e-16 * freq, 2.0), 0.25), -3.188) *
         pow(4.5625 + pow(0.4434 - 3.038699e-16 * freq, 2.0), -1.953) * freq * freq * (freq - TFIHE) /
         ((exp(PLANCK * freq / (BOLTZMANN * T)) - 1.0));
}

double sx_bottomMeanEnergyHe(double freq, void *p)  // 18.1+ eV
{
  double T = *(double *)p;
  return 9.492e-16 * (4.15752 + pow(1.4434 - 3.038699e-16 * freq, 2)) *
         pow(1.0 + 0.825067 * pow(4.5625 + pow(0.4434 - 3.038699e-16 * freq, 2.0), 0.25), -3.188) *
         pow(4.5625 + pow(0.4434 - 3.038699e-16 * freq, 2.0), -1.953) * freq * freq /
         (PLANCK * (exp(PLANCK * freq / (BOLTZMANN * T)) - 1.0));
}

void sx_chem_compute_mean_photon_energies(double T)
{
  int f;
  double intTop, intBottom, minFreq, maxFreq, error;

  // Initialize integrator
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
  gsl_function F;
  F.params = &T;

  // average photon energies (13.6 - 15.2 eV)
  minFreq = All.sxFreqBins[F136];
  maxFreq = All.sxFreqBins[F152];

  F.function = &sx_topMeanEnergyH_136;
  gsl_integration_qags(&F, minFreq, maxFreq, 0, 1e-7, 1000, w, &intTop, &error);
  F.function = &sx_bottomMeanEnergyH;
  gsl_integration_qags(&F, minFreq, maxFreq, 0, 1e-7, 1000, w, &intBottom, &error);
  All.sxEnergy[E136IH] = intTop / intBottom;

  // average photon energies (15.2-24.6 eV)
  minFreq = All.sxFreqBins[F152];
  maxFreq = All.sxFreqBins[F246];

  F.function = &sx_topMeanEnergyH_152;
  gsl_integration_qags(&F, minFreq, maxFreq, 0, 1e-7, 1000, w, &intTop, &error);
  F.function = &sx_bottomMeanEnergyH;
  gsl_integration_qags(&F, minFreq, maxFreq, 0, 1e-7, 1000, w, &intBottom, &error);
  All.sxEnergy[E152IH] = intTop / intBottom;

  F.function = &sx_topMeanEnergyH2_181;
  gsl_integration_qags(&F, sxEnergyPW[NPW], maxFreq, 0, 1e-7, 1000, w, &intTop, &error);  // 18.1+ eV
  F.function = &sx_bottomMeanEnergyH2_181;
  gsl_integration_qags(&F, sxEnergyPW[NPW], maxFreq, 0, 1e-7, 1000, w, &intBottom, &error);  // 18.1+ eV
  double intBottomPW, intTopPW;
  for(f = 0; f < NPW; f++)  // 15.2 - 18.1 eV
    {
      F.function = &sx_topMeanEnergyH2_152;
      gsl_integration_qags(&F, sxEnergyPW[f], sxEnergyPW[f + 1], 0, 1e-7, 1000, w, &intTopPW, &error);
      F.function = &sx_bottomMeanEnergyH2_152;
      gsl_integration_qags(&F, sxEnergyPW[f], sxEnergyPW[f + 1], 0, 1e-7, 1000, w, &intBottomPW, &error);
      intTop += sxSigmaPW[f] * intTopPW;
      intBottom += sxSigmaPW[f] * intBottomPW;
    }
  All.sxEnergy[E152IH2] = intTop / intBottom;

  // average photon energies (24.6+ eV)
  minFreq = All.sxFreqBins[F246];
  maxFreq = All.sxFreqBins[SX_NFREQ];

  F.function = &sx_topMeanEnergyH_152;
  gsl_integration_qags(&F, minFreq, maxFreq, 0, 1e-7, 1000, w, &intTop, &error);
  F.function = &sx_bottomMeanEnergyH;
  gsl_integration_qags(&F, minFreq, maxFreq, 0, 1e-7, 1000, w, &intBottom, &error);
  All.sxEnergy[E246IH] = intTop / intBottom;

  F.function = &sx_topMeanEnergyH2_181;
  gsl_integration_qags(&F, minFreq, maxFreq, 0, 1e-7, 1000, w, &intTop, &error);  // 18.1+ eV
  F.function = &sx_bottomMeanEnergyH2_181;
  gsl_integration_qags(&F, minFreq, maxFreq, 0, 1e-7, 1000, w, &intBottom, &error);  // 18.1+ eV
  All.sxEnergy[E246IH2] = intTop / intBottom;

  F.function = &sx_topMeanEnergyHe;
  gsl_integration_qags(&F, minFreq, maxFreq, 0, 1e-7, 1000, w, &intTop, &error);  // 18.1+ eV
  F.function = &sx_bottomMeanEnergyHe;
  gsl_integration_qags(&F, minFreq, maxFreq, 0, 1e-7, 1000, w, &intBottom, &error);  // 18.1+ eV
  All.sxEnergy[E246IHE] = intTop / intBottom;

  // Close integrator
  gsl_integration_workspace_free(w);
}

#undef NPW  // Number of PieceWise Intervals for H2

/***********************************************************************************************
 *                              Ionisation rates of a source                                   *
 ***********************************************************************************************/

#if (SX_SOURCES==5)  

double sx_planckDivNu( double freq, void * p )  // B_nu(T)/(h*nu)
{
  double T = *(double *)p;
  // constant factor 2/c^2 is taken out of the integral
  return freq * freq / (exp(PLANCK * freq / (BOLTZMANN * T)) - 1.0);
}

void sx_chem_compute_ionisation_rates(double ionisationRate[SX_NFREQ], int index)
{
  int f;

  double sourcePhotEmition, emissionFactor;
  double accretion, mass, radius, temperature;

  accretion = (SinkP[index].Mass-SinkP[index].MassOld) *
    sxRunData->mass2cgs * sxRunData->runTimeInv;          // mean mass accretion [g/s]
  //printf("SX: (%d) SinkPID %d Mass %.03e MassOld %.03e accretion %.03e ion [",ThisTask,
  // 	 SinkP[index].ID, SinkP[index].Mass, SinkP[index].MassOld, accretion);  
  SinkP[index].MassOld = SinkP[index].Mass;               // reset the MassOld value
  mass = SinkP[index].Mass * sxRunData->mass2cgs;         // [g]

  // testing values
  // accretion = 6.337e22;  // 1e-3 Msol/yr == 6.337e22 g/s
  // mass = 1.4e35;           // 70 Msol == 1.4e35 g

  radius      = 0;                                           // [cm]
  temperature = 0;                                           // [K]
  setup_popiii_tables_();                                    // setting the lookup tables for PopIII stars
  popiii_lookup_(&accretion, &mass, &radius, &temperature);  // calculate radius and temperature

  emissionFactor = 2.2253001121e-21 * 4.0 * M_PI * radius * radius;  // 2/c^2 * star surface [s^2 cm^-2 * cm^2]

  // Initialize integrator
  double error;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
  gsl_function F;
  F.params = &temperature;

  //printf("(%d) SINK: acc %.02e mass %.02e rad %.02e temp %.02e surf %.02e ionRate [ ",
  //     ThisTask, accretion, mass, radius, temperature, 4.0 * M_PI * radius * radius);

  //char buff[50];
  for(f = 0; f < SX_NFREQ; f++)
    {
      if(f == 0)  // for now we will exclude the dust pumping bin
        sourcePhotEmition = 0.0;
      else
        {
          F.function = &sx_planckDivNu;
          gsl_integration_qags(&F, All.sxFreqBins[f], All.sxFreqBins[f + 1], 0, 1e-7, 1000, w, &sourcePhotEmition, &error);
        }
      ionisationRate[f] = sourcePhotEmition * emissionFactor;
      /*
      if (f==2)
	ionisationRate[f] = 1e49;
      else
	ionisationRate[f] = 0.0;
      sprintf(&buff[f*9],"%.02e ",ionisationRate[f]);
      */
    }
  // printf("SX: (%d) ID %d Mass %.02e MassOld %.02e acc %.02e ion [ %s] \n",ThisTask,
  // 	 SinkP[index].ID, SinkP[index].Mass, SinkP[index].MassOld, accretion, buff);  

  //printf("] \n");
  //printf("] ion [ ");
  //for(f=0; f<SX_NFREQ; f++)
  //  printf("%.02e ",ionisationRate[f]);
  //printf("] \n");
  // Close integrator
  gsl_integration_workspace_free(w);
}

#elif (SX_SOURCES==4)

void sx_chem_compute_ionisation_rates( double ionisationRate[SX_NFREQ], int index )
{

  int f;
  
 // DEBUG: For now we are happy with the approximate ionization in one frequency
  for(f = 0; f < SX_NFREQ; f++)
    {
      if(f == 2)
        {
          ionisationRate[f] = 1e49;
        }
      else
        {
          ionisationRate[f] = 0.0;
        }
    }
}

#endif

#ifdef POPIII_SNE
#define N_IMF_BINS 8
void sx_chem_compute_popIII_ionisation_rates( double ionisationRate[SX_NFREQ], int index )
{

  int f, i, j;
  double sourcePhotEmition, emissionFactor;
  double mass, radius, temperature, luminosity;
 
  double m[N_IMF_BINS] = { 9.,15.,25.,40.,60.,80.,120.,200.};
  double log_T[N_IMF_BINS] = {4.622,4.759,4.850,4.900,4.943,4.970,4.981,4.999};
  double log_L[N_IMF_BINS] = {3.709,4.324,4.890,5.420,5.715,5.947,6.243,6.574};

  double sigma_sb =  5.6704e-5;


  
  for(f=0; f<SX_NFREQ; f++)
      ionisationRate[f]=0;
  

  for (j=0; j<SinkP[index].N_sne; j++)
  {
    mass = SinkP[index].stellar_mass[j]* All.UnitMass_in_g / SOLAR_MASS ;
   
    i = 0;
    while (mass > m[i+1])
      i++; 

    /*Linear interpolation between the values*/
    temperature = log_T[i] + ((log_T[i+1]-log_T[i]) / (m[i+1]-m[i])) * (mass - m[i]);
    temperature = pow(10.0,temperature);
    
    luminosity = log_L[i] + ((log_L[i+1]-log_L[i]) / (m[i+1]-m[i])) * (mass - m[i]);
    luminosity = pow(10.0, luminosity)*SOLAR_LUM;

    radius = sqrt(luminosity / (pow(temperature, 4) * RAD_CONST * M_PI * CLIGHT)); // [cm]
    printf("SX POP III STAR: stellar mass %e, radius %e \n", mass, radius/1.0e5);
 
    emissionFactor = 2.2253001121e-21 * 4.0 * M_PI * radius * radius;  // 2/c^2 * star surface [s^2 cm^-2 * cm^2]

    // Initialize integrator
    double error;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_function F;
    F.params = &temperature;

    for(f=0; f<SX_NFREQ; f++)
    {
      if (f==0)  // for now we will exclude the dust pumping bin
        sourcePhotEmition = 0.0;
      else
      {
        F.function = &sx_planckDivNu;
        gsl_integration_qags(&F, All.sxFreqBins[f], All.sxFreqBins[f+1], 0, 1e-7, 1000, w, &sourcePhotEmition, &error);
      }
      ionisationRate[f] += sourcePhotEmition * emissionFactor;
      printf("%.03e ",sourcePhotEmition);
    }
  
  printf("] ion [ ");
  for(f=0; f<SX_NFREQ; f++)
      printf("%.03e ",ionisationRate[f]);
  printf("] \n");

  // Close integrator
  gsl_integration_workspace_free(w);
  }
}
#endif

#endif
