/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/dust_vars.h
 * \date        MM/YYYY
 * \author      Ryan McKinnon
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#ifndef DUST_VARS_H
#define DUST_VARS_H

#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

extern int Ndust;

extern struct dust_particle
{
  int index;
  int active_idx;
  MyFloat NumNgb;
  MyFloat NormSph;
  MyFloat TotNgbMass;
  MyFloat Dhsmlrho;
  MyFloat LocalGasAccel[3];
  MyFloat LocalGradP[3];
#ifdef DL_GRAIN_BINS
  MyFloat LocalGasDensityH;
  MyFloat LocalGasZ;
  MyFloat LocalGasTemp;
  MyFloat NewNumGrains[DL_GRAIN_BINS];
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
  MyFloat NewBinSlopes[DL_GRAIN_BINS];
#endif
  MyFloat DeltaMassExpected;
  MyFloat DeltaMassActual;
  MyFloat DeltaMetalMasses[GFM_N_CHEM_ELEMENTS];
  MyFloat DeltaMomentum[3];
#ifdef DL_SNE_DESTRUCTION
  MyFloat LocalSNPrefactor;
#endif
#endif
#ifdef DL_WINDS
  MyFloat WindMassRate;
  MyFloat WindVel;
#ifdef GFM_BIPOLAR_WINDS
#if(GFM_BIPOLAR_WINDS == 3)
  MyFloat DensGasAngMomentum[3];
#else
  MyFloat GroupVel[3];
  MyFloat GroupGravAcc[3];
#endif
#endif
#endif
#ifdef DL_DRAG_BACKREACTION
  MyFloat DragMomentum[3];
  MyFloat DragDustThermal;
#endif
#ifdef DL_RADIATION
  MyFloat SigmaTot[MRT_BINS];
  MyFloat SigmaTotAbs[MRT_BINS];
  MyFloat SigmaTotGeo;
#ifdef DL_OUTPUT_RT_FLUX
  MyFloat LocalRT_F[MRT_BINS][3];
#endif
#ifdef DL_RADIATION_PRESSURE
  MyFloat DeltaRPMomentum[3];
#endif
#ifdef DL_THERMAL_IR
#ifdef DL_GRAIN_BINS
  MyFloat DustTemp[DL_GRAIN_BINS];
  MyFloat SigmaBinGeo[DL_GRAIN_BINS];
  MyFloat SigmaBinIRAbs[DL_GRAIN_BINS];
#else
  MyFloat DustTemp;
#endif
#endif
#endif
} * DustParticle;

#ifdef DL_GRAIN_BINS
extern struct grain_size_distribution
{
  double InternalDensity; /* Grain density in internal mass / um^3 */
  double *Edges;
  double *Midpoints;
  double *Widths;
  double *AvgMasses;
  double *EdgeMasses;
#if defined(DL_SHATTERING) || defined(DL_COAGULATION)
  double *VelocitiesCNM;                              /* CNM grain velocities in um / s */
  double *VelocitiesWIM;                              /* WIM grain velocities in um / s */
  double I_2_kj_Prefac[DL_GRAIN_BINS][DL_GRAIN_BINS]; /* prefactors for shattering and coagulation integrals */
#ifdef DL_SHATTERING
  double VelShat; /* Shattering threshold velocity in um / s */
#endif
#ifdef DL_COAGULATION
  double VelCoag[DL_GRAIN_BINS][DL_GRAIN_BINS]; /* Coagulation threshold velocity in um / s */
#endif
#endif
#ifdef DL_SNE_DESTRUCTION
  /* Consult equation 11 of Asano+ (2013) for the full definition of this xi
   * fraction.  In short, XiFrac[i][j] describes the differential number
   * fraction of grains shifted from bin j to bin i as a result of one
   * supernova shock. */
  double XiFrac[DL_GRAIN_BINS][DL_GRAIN_BINS];
  double SNIINumFrac; /* Number of SNII per unit solar mass */
#endif
#ifdef DL_PRODUCTION
  /* Initial grain size distributions on a relative scale. */
  double AGB_NumGrains[DL_GRAIN_BINS];
  double SNII_NumGrains[DL_GRAIN_BINS];
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
  double AGB_BinSlopes[DL_GRAIN_BINS];
  double SNII_BinSlopes[DL_GRAIN_BINS];
#endif
  /* Condensation efficiencies for SNII */
  double SNII_CondEff[GFM_N_CHEM_ELEMENTS];
  double ***AGB_CondEff;
  double ***AGB_CondEff_spline;
  int AGB_CondEff_NZ;
  int AGB_CondEff_NM;
  double *AGB_CondEff_Mass;
  double *AGB_CondEff_Metallicity;
#endif
} GSD;
#endif

enum gsd_collision_type
{
  GSD_SHATTERING,
  GSD_COAGULATION
};

#if(defined(DL_GRAIN_BINS) && (defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION))) || \
    defined(DL_DRAG_BACKREACTION)
extern int NforcesDust;

extern struct shatter_data_p
{
  MyFloat DustHsml;
  MyFloat DustNumNgb;
  MyFloat MinNumNgbDeviationDust;
  MyFloat MaxNumNgbDeviationDust;
  MyFloat DustDensity;
  MyFloat EnclosedMass;
#ifdef DL_DEREFINEMENT
  MyFloat ClosestDustR;
  MyIDType ClosestDustID;
  int ClosestDustTask;      /* task containing P info of neighbor */
  int ClosestDustIndex;     /* index into P of neighbor on above task */
  int ClosestDustTaskTree;  /* task containing tree entry of neighbor */
  int ClosestDustIndexTree; /* index into P or Tree_Points of neighbor on above task */
  int ClosestDustHasP;      /* whether we can access P of neighbor on TaskTree */
  int HighMassNeighbor;
  int IsDerefinementOrigin; /* whether this particle would like to derefine */
#endif
} * PShatter;
#endif

#ifdef DL_PRODUCTION
enum gsd_dnda_type
{
  GSD_DNDA_AGB,
  GSD_DNDA_SNII,
  GSD_DNDA_SNIA
};
#endif

#ifdef DL_RADIATION
extern struct dust_radiation_pressure
{
  gsl_spline2d *spline_Q_pr;
  gsl_interp_accel *a_acc;
  gsl_interp_accel *lam_acc;
  gsl_spline2d *spline_Q_abs;
  gsl_interp_accel *a_acc_abs;
  gsl_interp_accel *lam_acc_abs;
} DustRP_gra, DustRP_sil;

extern struct dust_rp_helpers
{
  double meanE[MRT_BINS];
  double photons2energy[MRT_BINS];
  double lam[MRT_BINS];
} DRT;
#endif

#endif
