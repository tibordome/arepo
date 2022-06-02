/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sfr_mcs/sfr_mcs_proto.h
 * \date        01/2019
 * \author     	Matthew C Smith
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - 17.01.2022 Ported into current codebase
 */

#include <stdio.h>
#include <stdlib.h>
#include "../dtypes.h"
#include "arepoconfig.h"

#ifndef SFR_MCS_VARS_H
#define SFR_MCS_VARS_H

#ifdef SFR_MCS

#ifndef SFR_MCS_SELECT_CRITERIA
#define SFR_MCS_SELECT_CRITERIA 0
#endif

#ifndef SFR_MCS_RATE_CRITERIA
#define SFR_MCS_RATE_CRITERIA 0
#endif

#ifdef IMF_SAMPLING_MCS
#ifndef N_STAR_SLOTS
#define N_STAR_SLOTS 8
#endif
#endif

/* Alternate definition to that in allvars.h as used with GFM */
extern struct star_particle_data
{
  unsigned int PID;
  MyFloat BirthTime;
#ifdef SFR_MCS_BIRTH_RECORDS
  MyFloat BirthPos[3];
  MyFloat BirthVel[3];
  MyFloat BirthDensity;
#endif

#if defined(SN_MCS) || defined(HII_MCS) || defined(PE_MCS)
  MyFloat Age;  // in yrs
#ifdef SFR_MCS_DELAY
  MyFloat TimeDelay;
#endif
  MyFloat InitialMass;
#if !defined(IMF_SAMPLING_MCS) && !defined(SB99_FIXED_Z)
  int iz;
#endif
#endif

#ifdef SN_MCS
  int N_SN;           /**<Current SNe this timestep*/
  int N_SN_cum;       /**<Total SNe experienced*/
  int N_SN_event_cum; /**<Total SNe events experienced*/
  MyFloat M_ej;
#if defined(IMF_SAMPLING_MCS) || defined(SN_MCS_VARIABLE_EJECTA)
  MyFloat Z_ej;
#endif
#ifdef SN_MCS_LOCATION_RECORDS
  MyFloat SNDensity;  // Last SN event experienced by this particle
  MyFloat SNTime;
  MyFloat SNPos[3];
#ifdef GRACKLE
  MyFloat SNTemperature;
#endif
#endif
#endif

#ifdef HII_MCS
  MyFloat S_Hii;
#ifdef HII_MCS_LR
  MyFloat EnergyPerPhoton;
#endif
#ifndef IMF_SAMPLING_MCS
  int photon_it_high;
#else
  MyFloat S_HiiArr[N_STAR_SLOTS];
#ifdef HII_MCS_LR
  MyFloat EnergyPerPhotonArr[N_STAR_SLOTS];
#endif
#endif  // IMF_SAMPLING_MCS
#endif  // HII_MCS

#ifdef PE_MCS
  MyFloat L_FUV;  // This is always in ergs/s
#ifndef IMF_SAMPLING_MCS
  int fuv_it_high;
#else
  MyFloat L_FUVArr[N_STAR_SLOTS];
#endif
#ifdef PE_MCS_STORE_DUST_COLUMN
  MyFloat DustColumn;  // Local dust column density D*N_H, D relative to solar
#endif
#endif

#ifdef IMF_SAMPLING_MCS
  MyFloat MassArr[N_STAR_SLOTS];
  MyFloat LifetimeArr[N_STAR_SLOTS];
#endif

} * StarP, *DomainStarBuf;

#if defined(SN_MCS) || defined(HII_MCS) || defined(PE_MCS)
#ifndef IMF_SAMPLING_MCS
extern struct sb99_data
{
#ifdef SN_MCS
  int N_t_SN;           /* Number of timesteps for SNe*/
  MyFloat *TimestepsSN; /* t bins */
  MyFloat t_min_SN;     /* first timestep */
  MyFloat t_max_SN;     /* last timestep */
#ifndef SB99_FIXED_Z
  MyFloat **RatesSN; /* accessed RatesSN[Z bin index][t bin index] */
#else
  MyFloat *RatesSN;       /* accessed RatesSN[t bin index] */
#ifdef SN_MCS_VARIABLE_EJECTA
  MyFloat *EjectaMass;    /* accessed EjectaMass[t bin index] */
  MyFloat *ZEjecta;       /* accessed ZEjecta[t bin index] */
#endif
#endif
#endif

#ifdef HII_MCS
  int N_t_Photons;           /* Number of timesteps for ionizing photons*/
  MyFloat *TimestepsPhotons; /* t bins */
#ifndef SB99_FIXED_Z
  MyFloat **RatesPhotons; /* accessed RatesPhotons[Z bin index][t bin index] */
#else
  MyFloat *RatesPhotons;  /* accessed RatesPhotons[t bin index] */
#endif
#endif

#ifdef PE_MCS
  int N_t_FUV;           /* Number of timesteps for FUV luminosity*/
  MyFloat *TimestepsFUV; /* t bins */
#ifndef SB99_FIXED_Z
  MyFloat **LuminosityFUV; /* accessed LuminosityFUV[Z bin index][t bin index] */
#else
  MyFloat *LuminosityFUV; /* accessed LuminosityFUV[t bin index] */
#endif
#endif

#ifndef SB99_FIXED_Z
  int N_Z;                /* Number of metallicities */
  MyFloat *Metallicities; /* Z bins */
#endif
} sb99;

#endif //ifndef IMF_SAMPLING_MCS

#ifdef HII_MCS
#define HII_MCS_IGNORE_FLAG ((MyIDType)-1)  // Note that if MyIDType is unsigned, this still works (all bits set true in this case)
#ifdef HII_MCS_ANISO
#ifndef HII_MCS_N_SIDE
#define HII_MCS_N_SIDE 1
#endif
#define HII_MCS_N_PIX (12 * HII_MCS_N_SIDE * HII_MCS_N_SIDE)
#endif
#endif

#ifdef PE_MCS
#ifndef PE_MCS_PRESHIELD
#define PE_MCS_PRESHIELD 0
#endif
#ifndef PE_MCS_POSTSHIELD
#define PE_MCS_POSTSHIELD 0
#endif
#endif

#if(PE_MCS_PRESHIELD == 1)
/** Struct used for passing the parameters during the mesh cell search. */
typedef struct
{
  MyDouble Pos[3];
  int Task;
  MyFloat DustColumn;
} dust_column_search_data;
#endif
#endif  // defined(SN_MCS) || defined(HII_MCS) || defined(PE_MCS)

#ifdef IMF_SAMPLING_MCS
/* In the future a file will be read, for now we hardcode */
extern struct star_properties_table
{
  int N_m; /*=34*/ /*The number of mass points to interpolate from */
  /*The masses in solar masses */
  MyFloat Masses[34];

  /*The lifetimes in years */
  MyFloat Lifetimes[34];

#ifdef SN_MCS
  /*Ratio of progenitor mass to ejecta mass */
  MyFloat EjectaMassRatio[34];

  MyFloat EjectaMetallicity[34];
#endif

#ifdef HII_MCS
  /*Log of the ionizing photon rate in 1e49 s^-1*/
  MyFloat LogIonizingPhotonRate1e49[34];
#ifdef HII_MCS_LR
  /*Mean energy per photon in ionising band, in ergs*/
  MyFloat EnergyPerPhoton[34];
#endif
#endif

#ifdef PE_MCS
  /*Log FUV luminosity in erg s^-1 */
  MyFloat LogLuminosityFUV[34];
#endif
} StarProperties;
#endif

#endif  // SFR_MCS
#endif  // SFR_MCS_VARS_H
