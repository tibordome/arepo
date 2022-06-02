/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/GFM/stellar_feedback_kernels.h
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

#ifndef STELLAR_FEEDBACK_KERNELS_H
#define STELLAR_FEEDBACK_KERNELS_H

#ifndef SMUGGLE_STAR_FEEDBACK
typedef struct
{
  MyDouble Pos[3];
  MyDouble Vel[3];
  MyFloat Hsml;
  MyFloat NormSph;
  MyFloat TotalMassReleased;
#ifdef GFM_WINDS_LOCAL
  MyFloat WindEnergyReleased;
#endif
#ifdef GFM_STELLAR_FEEDBACK
  MyDouble SNIaEnergyReleased;
  MyDouble AGBMomentumReleased;
#endif

  int Firstnode;
} data_in;

typedef struct
{
  char dummy;
} data_out;

int is_doing_stellar_feedback(int i);

#endif

#ifdef GFM_STELLAR_FEEDBACK
void GFM_stellar_feedback(int target, int mode, int thread_id, int numnodes, int *firstnode, data_in *in);
#endif

#ifdef GFM_WINDS_LOCAL
void GFM_winds_local(int target, int mode, int thread_id, int numnodes, int *firstnode, data_in *in);
#endif

#endif
