/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/GFM/stellar_feedback_kernels.c
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

#include "../allvars.h"
#include "../proto.h"

#include "stellar_feedback_kernels.h"

int is_doing_stellar_feedback(int i)
{
#ifdef SMUGGLE_STAR_FEEDBACK
  if(StarParticle[i].TotalEnergyReleased > 0)
    return 1;
  if(StarParticle[i].TotalMassReleasedAGB > 0)
    return 1;
#endif
#ifdef GFM_STELLAR_FEEDBACK
  if(StarParticle[i].SNIaEnergyReleased > 0 || StarParticle[i].AGBMomentumReleased > 0)
    return 1;
#endif
#ifdef GFM_WINDS_LOCAL
  if(StarParticle[i].WindEnergyReleased > 0)
    return 1;
#endif

  return 0;
}

#ifdef GFM_STELLAR_FEEDBACK
void GFM_stellar_feedback(int target, int mode, int thread_id, int numnodes, int *firstnode, data_in *in)
{
  MyDouble inj_mom[3], *pos, SNIaEnergyReleased, AGBMomentumReleased;
  MyFloat normsph;
  double weight_fac, h, h2;
  double dr[3], r2, r;

  pos                 = in->Pos;
  h                   = in->Hsml;
  normsph             = in->NormSph;
  SNIaEnergyReleased  = in->SNIaEnergyReleased * All.cf_atime * All.cf_atime;
  AGBMomentumReleased = in->AGBMomentumReleased * All.cf_atime;

#ifndef GFM_TOPHAT_KERNEL
  double wk;
  MyFloat hinv = 1.0 / h;
#ifndef TWODIMS
  MyFloat hinv3 = hinv * hinv * hinv;
#else
  MyFloat hinv3 = hinv * hinv / boxSize_Z;
#endif
#endif

  h2 = h * h;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  for(int n = 0; n < nfound; n++)
    {
      int j = Thread[thread_id].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0) /* skip cells that have been swallowed or dissolved */
        {
          dr[0] = NEAREST_X(pos[0] - P[j].Pos[0]);
          dr[1] = NEAREST_Y(pos[1] - P[j].Pos[1]);
          dr[2] = NEAREST_Z(pos[2] - P[j].Pos[2]);

          r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];

          if(r2 < h2)
            {
              r = sqrt(r2);

#ifndef GFM_TOPHAT_KERNEL
              double u = r * hinv;

              if(u < 0.5)
                wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
              else
                wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

              weight_fac = SphP[j].Volume * wk / normsph;
#else
              weight_fac = SphP[j].Volume / normsph;
#endif

              /* injected feedback energy */
              MyDouble de_feedback = weight_fac * SNIaEnergyReleased;
              /* injected feedback momentum */
              MyDouble dp_feedback = weight_fac * AGBMomentumReleased;

              inj_mom[0] = -dp_feedback * dr[0] / r;
              inj_mom[1] = -dp_feedback * dr[1] / r;
              inj_mom[2] = -dp_feedback * dr[2] / r;

              /* this accounts both for thermal and kinetic energy from feedback */
              SphP[j].Energy += de_feedback;
              /* momentum due to stellar mass return is injected in radial direction */
              SphP[j].Momentum[0] += inj_mom[0];
              SphP[j].Momentum[1] += inj_mom[1];
              SphP[j].Momentum[2] += inj_mom[2];

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
              SphP[j].Utherm =
                  (SphP[j].Energy - 0.5 *
                                        (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                                         SphP[j].Momentum[2] * SphP[j].Momentum[2]) /
                                        P[j].Mass) /
                  P[j].Mass / (All.cf_atime * All.cf_atime);
              SphP[j].A       = (GAMMA - 1.0) * SphP[j].Utherm / pow(SphP[j].Density * All.cf_a3inv, GAMMA - 1);
              SphP[j].Entropy = log(SphP[j].A) * P[j].Mass;
#endif
            }
        }
    }
}
#endif

#ifdef GFM_WINDS_LOCAL
void GFM_winds_local(int target, int mode, int thread_id, int numnodes, int *firstnode, data_in *in)
{
  MyDouble *pos, WindEnergyReleased;
  MyFloat normsph;
  double weight_fac, h, h2;
  double dr[3], r2, r;

  pos                = in->Pos;
  h                  = in->Hsml;
  normsph            = in->NormSph;
  WindEnergyReleased = in->WindEnergyReleased;

#ifndef GFM_TOPHAT_KERNEL
  double wk;
  MyFloat hinv = 1.0 / h;
#ifndef TWODIMS
  MyFloat hinv3 = hinv * hinv * hinv;
#else
  MyFloat hinv3 = hinv * hinv / boxSize_Z;
#endif
#endif

  h2 = h * h;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  for(int n = 0; n < nfound; n++)
    {
      int j = Thread[thread_id].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0) /* skip cells that have been swallowed or dissolved */
        {
          dr[0] = NEAREST_X(pos[0] - P[j].Pos[0]);
          dr[1] = NEAREST_Y(pos[1] - P[j].Pos[1]);
          dr[2] = NEAREST_Z(pos[2] - P[j].Pos[2]);

          r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];

          if(r2 < h2)
            {
              r = sqrt(r2);

#ifndef GFM_TOPHAT_KERNEL
              double u = r * hinv;

              if(u < 0.5)
                wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
              else
                wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

              weight_fac = SphP[j].Volume * wk / normsph;
#else
              weight_fac = SphP[j].Volume / normsph;
#endif
              SphP[j].WindEnergyReceived += weight_fac * WindEnergyReleased;

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
              SphP[j].Utherm =
                  (SphP[j].Energy - 0.5 *
                                        (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                                         SphP[j].Momentum[2] * SphP[j].Momentum[2]) /
                                        P[j].Mass) /
                  P[j].Mass / (All.cf_atime * All.cf_atime);
              SphP[j].A       = (GAMMA - 1.0) * SphP[j].Utherm / pow(SphP[j].Density * All.cf_a3inv, GAMMA - 1);
              SphP[j].Entropy = log(SphP[j].A) * P[j].Mass;
#endif
            }
        }
    }
}
#endif
