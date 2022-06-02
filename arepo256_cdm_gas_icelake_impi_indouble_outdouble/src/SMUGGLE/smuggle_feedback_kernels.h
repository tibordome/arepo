/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/SMUGGLE/smuggle_feedback_kernels.h
 * \date        03/2020
 * \author
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#ifndef SMUGGLE_FEEDBACK_KERNELS_H
#define SMUGGLE_FEEDBACK_KERNELS_H

#ifdef SMUGGLE_STAR_FEEDBACK

typedef struct
{
  MyDouble Pos[3];
  MyDouble Vel[3];
  MyFloat Hsml;
  MyFloat NormSph;
  MyFloat TotalMassReleased;
#ifdef SMUGGLE_MASS_WEIGHT_SN
  MyFloat TotNgbMass; 
#endif
#ifdef SMUGGLE_OMEGA_WEIGHT_SN
  MyFloat TotSolidAngle; 
#endif
#ifdef GFM_WINDS_LOCAL
  MyFloat WindEnergyReleased;
#endif
#ifdef GFM_STELLAR_FEEDBACK
  MyDouble SNIaEnergyReleased;
  MyDouble AGBMomentumReleased;
#endif
#ifdef SMUGGLE_STAR_FEEDBACK
  MyDouble TotalEnergyReleased;
#ifdef SMUGGLE_SN_COOLING_RADIUS_BOOST
  MyFloat n_SNII;
  MyFloat n_SNIa;
  MyFloat TotalMassReleasedSNII;
  MyFloat TotalMassReleasedSNIa;
  MyFloat TotalMassReleasedAGB;
#endif
#ifdef SMUGGLE_AGB_WINDS
  MyFloat AGBWindSpeed;
#endif
  MyFloat NumNgb; 
#endif
#if defined(SMUGGLE_SN_COOLING_RADIUS_BOOST) || defined(SMUGGLE_RADIATION_FEEDBACK)
  MyFloat LocISMdens;
  MyFloat LocISMdensH;
  MyFloat LocISMmet;
  MyFloat FeedbackRadiusLimiter;
#endif
  int Firstnode;
} data_in;

typedef struct
{
  MyDouble TotalEnergyInjected;
  MyDouble TotalMomentumReleased;
  MyDouble TotalMomentumInjected;
#ifdef SMUGGLE_AGB_WINDS
  MyDouble TotalMomentumInjectedAGB;
#endif
} data_out;

int is_doing_stellar_feedback(int i);

#ifdef SMUGGLE_SN_COOLING_RADIUS_BOOST
void cooling_radius_momentum_feedback(int target, int mode, int thread_id, int numnodes, int *firstnode, data_in *in, data_out *out);
#endif

#endif

#endif
