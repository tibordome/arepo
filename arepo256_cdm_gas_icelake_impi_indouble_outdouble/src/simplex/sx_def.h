/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \date        22/10/2018
 * \author      Ondrej Jaura
 * \brief       Pre-processor definitions, settings and constants
 * \details     
 */

#ifndef SX_DEF_H
#define SX_DEF_H

// Specify which image_flags we allow
#if (REFLECTIVE_X==2)&&(REFLECTIVE_Y==2)&&(REFLECTIVE_Z==2)
#define SX_DP_IMAGE_FLAGS 939524097  // !!! adding special flag for outflow
#else
#define SX_DP_IMAGE_FLAGS 1
#endif
#define SX_DC_IMAGE_FLAGS 1

// Default configuration
#ifndef SX_CHEMISTRY
#define SX_CHEMISTRY 3
#endif
#ifndef SX_NUM_ROT
#define SX_NUM_ROT 0
#endif
#ifndef SX_NUM_STEPS
#define SX_NUM_STEPS 0
#endif
#ifndef SX_RUNNING_MODE
#define SX_RUNNING_MODE 2
#endif

// Memory allocation factors
#define SX_ABPP_INIT 0.3
#define SX_ABPP_REALLOC 1.2
#define SX_QPP_INIT 0.2
#define SX_QPP_REALLOC 1.2

// Some other constants
#define SX_NDIM 3                              // number of dimensions
#define SX_MAX_STRAIGHT_ANGLE 90.0             // maximal straight angle during the transfer
#define SX_FACTOR_SA 5.0                       // multiplyier of the solid angle of a directional bin
#define ELECTRONVOLT_IN_HZ 2.417990504024e+14  // [eV] -> [Hz]

// Threshold frequencies [Hz]
#define TFHDUST 1.354074682254e15  // 5.6  eV - dust attenuation of photon electric heating radiation
#define TFDH2 2.708149364507e15    // 11.2 eV - dissociation frequency of H2
#define TFIH 3.28798e15            // 13.6 eV - hydrogen ionisation frequency
#define TFIH2 3.67534e15           // 15.2 eV - ionisation frequency of H2
#define TFIHE 5.9453e15            // 24.6 eV - helium ionisation frequency (from Cen 1992)
#define TFIIHE 1.3158e16           // 54.4 eV - singly ionised helium ionisation frequency (from Cen 1992)

enum sxTimersList {
  TRUN,
#ifdef SX_DISPLAY_TIMERS
  TRUNs, TRUNe, TRUNi,
  TSTEP, TSTEPs, TSTEPe, TSTEPi, TEVOL,
  TEXCHe, TEXCHc,
#endif
#ifdef SX_DISPLAY_TIMERS_SITE
  TEVOLsi, TEVOLss, TEVOLn, TEVOLtd, TEVOLtb, TEVOLi,
#endif
#if (SX_NUM_ROT>1) 
  TROT,
#ifdef SX_DISPLAY_TIMERS
  TROTs, TROTe,
#endif
#endif
  SX_NTIMERS          // Number of timers
};

#if (SX_CHEMISTRY==3) // SGChem 

#if !defined(SINK_SIMPLEX)&&(SX_SOURCES==5)
#define SINK_SIMPLEX       // switch sink particle 'OldMass'
#endif
#ifndef SGCHEM_RT
#define SGCHEM_RT          // switch SGChem
#endif
#ifndef SX_SINGLE_SHOT
#define SX_SINGLE_SHOT 1   // standard setting of photon emmition (0=continous)
#endif

// Photon energy bins
enum sxFreqencyBins
{
  F056,
  F112,
  F136,
  F152,
  F246,
  SX_NFREQ  // number of frequency bins
};

// Rates
enum sxChemicalRates
{
  RIH,
  HRIH,  // Ionisation and heating of H
  RIH2,
  HRIH2,  // Ionisation and heating of H2
  RDH2,
  HRD,  // Dissociation of H2 and heating of dust
  RIHE,
  HRIHE,     // Ionisation and heating of He
  SX_NRATES  // number of chemical rates
};

#define SX_NIMAGES 8  // number of properties to output in a slice image (=SX_NRATES)

// Cross-sections
enum sxCrossSections {
  S112H, S112DH2,           // 11.2 eV - Absorption of LW photons by atomic H
  S136DH2, S136IH,          // 13.6 eV - H2 dissociation and H ionisation
  S152IH, S152IH2,          // 15.2 eV - H, H2 ionisation
  S246IHE, S246IH, S246IH2, // 24.6 eV - He, H, H2 ionisation
  SX_NSIGMA                 // number of cross-sections
};

// Heating energies
enum sxExcessEnergies {
  E136IH,                   // 13.6 eV - H ionisation
  E152IH, E152IH2,          // 15.2 eV - H, H2 ionisation
  E246IHE, E246IH, E246IH2, // 24.6 eV - He, H, H2 ionisation
  SX_NENERGY                // number of energies
};

#elif(SX_CHEMISTRY == 4)  // FiBY chemistry

enum sxFrequencyBins {
  F136, F191, F246,
  F301, F356, F450,
  F544, F638, F733,
  F827, 
  SX_NFREQ              // number of frequency bins
};

enum sxCrossSections
{
  CSH,
  CSHE,
  CSHEP,
  SX_NSIGMA  // number of cross-sections
};

enum sxExcessEnergies
{
  EEH,
  EEHE,
  EEHEP,
  SX_NENERGY  // number of energies
};

enum sxMassFractions
{
  MFH,
  MFHP,
  MFHE,
  MFHEP,
  MFHEPP,
  SX_NMASSFRACT  // number of mass fractions
};

#define SX_NIMAGES 5  // number of properties to output in a slice image (=SX_NMASSFRACT)

#endif

#endif
