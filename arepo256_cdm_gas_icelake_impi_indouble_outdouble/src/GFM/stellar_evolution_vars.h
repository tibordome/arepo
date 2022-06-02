/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/GFM/stellar_evolution_vars.h
 * \date        MM/YYYY
 * \author
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#ifndef STELLAR_EVOLUTION_VARS_H
#define STELLAR_EVOLUTION_VARS_H

#include "../allvars.h"

/* additional compile time options, that are too specific for general Config.sh file */
#define GFM_TOPHAT_KERNEL /* use top hat for enrichment */
#ifdef GFM_SINGLE_CELL_INJECTION
#define GFM_EXACT_NUMNGB
#endif
// GFM_EXACT_NUMNGB                         /* exact ngb cell number for enrichment/local feedback */
//#define GFM_SNIA_DTD_EFOLDING            /* use e-folding DTD for SNIa rate */
#define GFM_SNIA_DTD_POWERLAW            /* use powerlaw DTD for SNIa rate (default) */
#define GFM_SNIA_DTD_POWERLAW_INDEX 1.12 /* DTD ~ pow(t,GFM_SNIA_DTD_POWERLAW_INDEX) */
#define GFM_STELLAR_EJECTA_TEMPERATURE \
  0 /* stellar mass return is assumed to have internal energy corresponding to this temperature (in Kelvin) */

#ifdef GFM_RPROCESS
#define GFM_NSNS_DTD_POWERLAW_INDEX 1.12 /* DTD ~ pow(t,GFM_NSNS_DTD_POWERLAW_INDEX) */
#endif

/* basic enrichment constants */
#ifdef GFM_STELLAR_EVOLUTION_NO_ELEMENTS
#define GFM_N_CHEM_ELEMENTS 0
#else
#ifdef GFM_NORMALIZED_METAL_ADVECTION
#ifdef GFM_SPROCESS
#define GFM_N_CHEM_ELEMENTS 15 /* number of chemical elements */
#else
#define GFM_N_CHEM_ELEMENTS 10 /* number of chemical elements */
#endif
#else
#ifdef GFM_SPROCESS
#define GFM_N_CHEM_ELEMENTS 14 /* number of chemical elements */
#else
#define GFM_N_CHEM_ELEMENTS 9 /* number of chemical elements */
#endif
#endif
#endif

#define GFM_N_MASS_BINS 200   /* number of IMF mass bins */
#define GFM_EL_NAME_LENGTH 12 /* maximum number of characters in element names */
#define GFM_MIN_METAL -20.0   /* minimum metallicity */

#ifdef GFM_DUST
enum
{
  GFM_DUST_AGB,
  GFM_DUST_SNII,
  GFM_DUST_SNIa,
  GFM_DUST_N_CHANNELS
};
#endif

#ifdef GFM_DUST
/* Will automatically handle any pseudoelement introduced by GFM_NORMALIZED_METAL_ADVECTION. */
#define GFM_DUST_COUNT_SCALARS (GFM_DUST_N_CHANNELS * GFM_N_CHEM_ELEMENTS)
#else
#define GFM_DUST_COUNT_SCALARS 0
#endif

#ifdef GFM_DUST
/* Species masses in amu. */
#define GFM_DUST_AMU_MG 24.305
#define GFM_DUST_AMU_SI 28.0855
#define GFM_DUST_AMU_FE 55.845
#endif

#ifdef GFM_CHEMTAGS
#ifndef GFM_SPLITFE
#ifndef GFM_RPROCESS
#define GFM_N_CHEM_TAGS 3 /* tags: SNIa, SNII, AGB -> note NSNS is assumed proportional to SNIa - don't need direct tagging */
#else
#define GFM_N_CHEM_TAGS 4 /* tags: SNIa, SNII, AGB, NSNS  */
#define GFM_NSNS_CHEMTAG 3
#endif
#define GFM_SNIA_CHEMTAG 0
#define GFM_SNII_CHEMTAG 1
#define GFM_AGB_CHEMTAG 2
#else  // with split fe
#ifndef GFM_RPROCESS
#define GFM_N_CHEM_TAGS \
  5 /* tags: SNIa, SNII, AGB, FESNIA, FESNII -> note NSNS is assumed proportional to SNIa - don't need direct tagging */
#define GFM_FESNIA_CHEMTAG 3
#define GFM_FESNII_CHEMTAG 4
#else
#define GFM_N_CHEM_TAGS 6 /* tags: SNIa, SNII, AGB, NSNS, FESNIA, FESNII  */
#define GFM_NSNS_CHEMTAG 3
#define GFM_FESNIA_CHEMTAG 4
#define GFM_FESNII_CHEMTAG 5
#endif
#define GFM_SNIA_CHEMTAG 0
#define GFM_SNII_CHEMTAG 1
#define GFM_AGB_CHEMTAG 2
#endif
#endif

/* energies */
#define GFM_SNII_ENERGY 1e51 /* SNII energy for local winds */

/* names of yield tables */
#define GFM_AGB_YIELDNAME "AGB.hdf5" /* hard coded yield table names */
#define GFM_SNIa_YIELDNAME "SNIa.hdf5"
#define GFM_SNII_YIELDNAME "SNII.hdf5"
#define GFM_LIFETIMENAME "Lifetimes.hdf5"

/* for mode of IMF integration */
#define INTEGRATE_IMF_NUMBER 0
#define INTEGRATE_IMF_MASS 1
#define INTEGRATE_IMF_YIELD 2

/* initial gas abundances */
#define GFM_INITIAL_ABUNDANCE_HYDROGEN HYDROGEN_MASSFRAC
#define GFM_INITIAL_ABUNDANCE_HELIUM (1. - GFM_INITIAL_ABUNDANCE_HYDROGEN)
#define GFM_INITIAL_ABUNDANCE_CARBON 0
#define GFM_INITIAL_ABUNDANCE_NITROGEN 0
#define GFM_INITIAL_ABUNDANCE_OXYGEN 0
#define GFM_INITIAL_ABUNDANCE_NEON 0
#define GFM_INITIAL_ABUNDANCE_MAGNESIUM 0
#define GFM_INITIAL_ABUNDANCE_SILICON 0
#define GFM_INITIAL_ABUNDANCE_IRON 0
#define GFM_INITIAL_ABUNDANCE_OTHER 0
#define GFM_INITIAL_METALLICITY 0

/* solar gas abundances (from Asplund+ 2009) */
#define GFM_SOLAR_ABUNDANCE_HYDROGEN 0.7388
#define GFM_SOLAR_ABUNDANCE_HELIUM (1. - GFM_SOLAR_ABUNDANCE_HYDROGEN - GFM_SOLAR_METALLICITY)
#define GFM_SOLAR_ABUNDANCE_CARBON 0.0024
#define GFM_SOLAR_ABUNDANCE_NITROGEN 0.0007
#define GFM_SOLAR_ABUNDANCE_OXYGEN 0.0057
#define GFM_SOLAR_ABUNDANCE_NEON 0.0012
#define GFM_SOLAR_ABUNDANCE_MAGNESIUM 0.0007
#define GFM_SOLAR_ABUNDANCE_SILICON 0.0007
#define GFM_SOLAR_ABUNDANCE_IRON 0.0013
#define GFM_SOLAR_ABUNDANCE_OTHER 0

extern int element_index_Hydrogen;
extern int element_index_Helium;
extern int element_index_Carbon;
extern int element_index_Nitrogen;
extern int element_index_Oxygen;
extern int element_index_Neon;
extern int element_index_Magnesium;
extern int element_index_Silicon;
extern int element_index_Iron;

extern char ElementNames[GFM_N_CHEM_ELEMENTS][GFM_EL_NAME_LENGTH];
extern double imf_by_number[GFM_N_MASS_BINS], imf_mass_bin[GFM_N_MASS_BINS], imf_mass_bin_log10[GFM_N_MASS_BINS];
extern double yield_mass_bin[GFM_N_MASS_BINS], imf_integrand_mass[GFM_N_MASS_BINS];
extern double imf_dlog10_Msun;

extern struct YieldTypeIa
{
  int N_ELEMENTS;             /* number of yields */
  char **ElementName;         /* name of species */
  MyFloat *Yield;             /* Yield (in solar masses) */
  MyFloat *spline;            /* Yields for corresponding spline elements (in solar masses) */
  MyFloat TotalMetals_spline; /* Total metal mass */
} yieldsSNIa;

extern struct YieldTypeII_and_AGB
{
  int N_ELEMENTS, N_MASS, N_Z;  /* number of elements, mass, and initial metallicity bins */
  char **ElementName;           /* name of element */
  MyFloat *Mass;                /* mass bins[N_MASS] */
  MyFloat *Metallicity;         /* metallicity bins[N_Z] (total metallicity/total mass) */
  MyFloat **Ejecta;             /* size [N_MASS*N_Z] with index_3d(imass,iz)  CHANGE MADE HERE! */
  MyFloat **Ejecta_spline;      /* size [N_MASS*N_Z] with index_3d(imass,iz)  CHANGE MADE HERE! */
  MyFloat ***Yield;             /* size [N_ELEMENTS*N_MASS*N_Z] with index_3d(iel,imass,iz) */
  MyFloat ***spline;            /* Yields for corresponding spline elements (in solar masses) */
  MyFloat **TotalMetals;        /* Total metal mass */
  MyFloat **TotalMetals_spline; /* Total metal mass */
} yieldsSNII, yieldsAGB;

extern struct Lifetime_Table
{
  int N_MASS, N_Z;      /* number of elements, mass, and initial metallicity bins */
  MyFloat *Mass;        /* mass bins[N_MASS]     */
  MyFloat *Metallicity; /* metallicity bins[N_Z] (total metallicity/total mass) */
  MyFloat **Dyingtime;  /* size [N_MASS*N_Z] with index_2d(imass,iz) */
} Lifetimes;

extern struct star_particle
{
  int index;
  MyDouble NumNgb;
  MyDouble NormSph;
  MyDouble Dhsmlrho;
  MyDouble TotalMassReleased;
  MyFloat TotalMetalMassReleased;
  MyFloat MetalMassReleased[GFM_N_CHEM_ELEMENTS];
#ifdef GFM_RPROCESS_CHANNELS
  MyFloat MassReleasedRProcess[GFM_RPROCESS_CHANNELS];
#endif
#ifdef GFM_DUST
  MyFloat DustMassReleased[GFM_DUST_N_CHANNELS][GFM_N_CHEM_ELEMENTS];
#endif
#if defined(GFM_DUST) || defined(DUST_LIVE) || defined(SMUGGLE_SN_COOLING_RADIUS_BOOST)
  MyFloat NumSNII;
#endif
#ifdef GFM_SNIA_ENERGY_INJECTION
  int NumSNIa;
#endif
#ifdef GFM_CHEMTAGS
  MyFloat MetalMassReleasedChemTags[GFM_N_CHEM_TAGS];
#endif
  MyFloat ClosestNeighbourDistance;
#ifdef GFM_STELLAR_FEEDBACK
  MyDouble SNIaEnergyReleased;
  MyDouble AGBMomentumReleased;
#endif
#ifdef GFM_WINDS_LOCAL
  MyFloat WindEnergyReleased;
#endif
#ifdef GFM_INJECT_B_FROM_SN
  MyDouble m[3];
#endif
#ifdef SMUGGLE_STAR_FEEDBACK
  MyDouble TotalEnergyReleased;
  MyDouble TotalMomentumReleased;
  MyDouble TotalMomentumInjected;
  MyDouble TotalMomentumInjectedAGB;
  MyDouble TotalEnergyInjected;
#ifdef SMUGGLE_SN_COOLING_RADIUS_BOOST
  MyFloat NumSNIa;
  MyDouble TotalMassReleasedSNII;
  MyDouble TotalMassReleasedSNIa;
  MyDouble TotalMassReleasedAGB;
#endif
#ifdef SMUGGLE_AGB_WINDS
  MyDouble AGBWindSpeed;
#endif
#endif
#if defined(SMUGGLE_SN_COOLING_RADIUS_BOOST) || defined(SMUGGLE_RADIATION_FEEDBACK)
  MyDouble LocISMdens;
  MyDouble LocISMdensH;
  MyDouble LocISMmet;
  MyDouble FeedbackRadiusLimiter;
#endif
#ifdef SMUGGLE_RADIATION_FEEDBACK
  MyDouble RadiationMomentumReleased;
  MyFloat StromgrenRadius;
  MyFloat GasColumnDensity;
#ifdef SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION
  MyFloat StromgrenMass;
  MyFloat NormSphRadFeedback_cold;
  MyFloat Lum;
  MyFloat LowestDensity;
  MyFloat LowestDensityDirection_x;
  MyFloat LowestDensityDirection_y;
  MyFloat LowestDensityDirection_z;
#endif
  MyFloat NormSphRadFeedback;
  MyFloat RadCoolShutoffTime;
  MyFloat RadFeedTau;
  MyFloat RadFeed_MinGasDist;
  int RadFeed_NumNgb;
#endif
  MyDouble deltaEKin;
  MyDouble deltaMomentum[3];
#ifdef SMUGGLE_VAR_SN_EFF
  MyDouble AvgMetalNgb;
#endif
#ifdef SMUGGLE_MASS_WEIGHT_SN
  MyDouble TotNgbMass;
#endif
#ifdef SMUGGLE_OMEGA_WEIGHT_SN
  MyDouble TotSolidAngle;
#endif
#ifdef MRT_SOURCES
  MyFloat TotalPhotReleased[MRT_BINS];
#endif
} * StarParticle;

typedef struct
{
  short int mass_return_flag;
  MyFloat total_mass_released;
  MyFloat total_metal_mass_released;
  MyFloat metal_mass_released[GFM_N_CHEM_ELEMENTS];
#if defined(GFM_DUST) || (defined(DUST_LIVE) && defined(DL_PRODUCTION))
  MyFloat total_dust_mass_released;
#endif
#ifdef GFM_RPROCESS_CHANNELS
  int number_of_NSNS[GFM_RPROCESS_CHANNELS];
  MyFloat rprocess_mass_released[GFM_RPROCESS_CHANNELS];
#endif
#ifdef GFM_DUST
  MyFloat dust_mass_released[GFM_DUST_N_CHANNELS][GFM_N_CHEM_ELEMENTS];
#endif
#if defined(DUST_LIVE) && defined(DL_PRODUCTION)
  MyFloat dust_mass_released[GFM_N_CHEM_ELEMENTS];
#endif
#ifdef GFM_CHEMTAGS
  MyFloat metal_mass_released_chemtags[GFM_N_CHEM_TAGS];
#endif
#ifdef GFM_RPROCESS
  MyFloat number_of_NSNS;
#endif
#ifdef GFM_SNIA_ENERGY_INJECTION
  int NumSNIa;
#endif
  MyFloat number_of_SNIa;
  MyFloat number_of_SNII;
  MyFloat AGB_mass_released, SNIa_mass_released, SNII_mass_released;
#ifdef SMUGGLE_DISCRETE_SN
  MyFloat initial_mass;
  MyFloat remaining_mass_frac;
#endif
#ifdef SMUGGLE_AGB_WINDS
  MyFloat OB_mass_released;
  short int OB_return_flag;
#endif
} stellar_evolution_data;

extern int Nstar;

#ifdef GFM_WINDS_STRIPPING
extern int Nwinds_to_strip;
#endif

#endif
