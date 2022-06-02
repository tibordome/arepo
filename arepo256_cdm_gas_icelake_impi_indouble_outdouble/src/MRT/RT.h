/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/MRT/RT.h
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

#ifndef MRT_RT_H
#define MRT_RT_H

#include "../allvars.h"

#ifdef MRT

#if defined(MRT_RIEMANN_HLLE) || defined(MRT_RIEMANN_HLLE_NEW)
extern double lambda1[101][101];
extern double lambda2[101][101];
extern double lambda3[101][101];
extern double lambda4[101][101];
#endif

#if defined(MRT_IR) && defined(MRT_IR_GRAIN_KAPPA)
extern int IR_N_pts;
extern double *IR_logT, *IR_logkappaP, *IR_logkappaR;
extern gsl_interp_accel *accIR_kappaP;
extern gsl_interp_accel *accIR_kappaR;
extern gsl_spline *splineIR_kappaP;
extern gsl_spline *splineIR_kappaR;
#endif

#if defined(MRT_SOURCES) || defined(MRT_LOCAL_FEEDBACK)
extern int Nsource;
extern double PhotRate;
extern int Nfreq, N_age, N_metallicity;
extern double *PhotEnergy, *LumPhotons;
extern double *LogAge, *LogMetallicity;
extern double *lum_tab;
#endif /* defined(MRT_SOURCES) || defined(MRT_LOCAL_FEEDBACK) */

#ifdef MRT_SINGLE_STAR
extern struct stellarparameters
{
  double *Age;
  double *Nphot;
  double **frac;

  double **sigH2;
  double **EH2;
  double **PH2;

  double **sigH;
  double **EH;
  double **PH;

  double **sigHe;
  double **EHe;
  double **PHe;

} StellarParameters;
#endif

/* required structures */

extern struct rt_face_data
{
  struct geometry geom;
  int face_active;
  int face_responsibility;
  double face_timestep;
  double vel_face[3];
  double vel_face_turned[3];
} * rtfacedata;

#endif /* MRT */

#endif
