/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/GFM/stellar_feedback.c
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

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

#if defined(GFM_STELLAR_FEEDBACK) || defined(GFM_WINDS_LOCAL) || defined(SMUGGLE_STAR_FEEDBACK)

#ifndef SMUGGLE_STAR_FEEDBACK
#include "stellar_feedback_kernels.h"
#else
#include "../SMUGGLE/smuggle_feedback_kernels.h"
#endif

static data_in *DataGet;
static data_out *DataResult;

static void particle2in(data_in *in, int i, int firstnode)
{
  in->Pos[0] = P[StarParticle[i].index].Pos[0];
  in->Pos[1] = P[StarParticle[i].index].Pos[1];
  in->Pos[2] = P[StarParticle[i].index].Pos[2];

  in->Vel[0] = P[StarParticle[i].index].Vel[0];
  in->Vel[1] = P[StarParticle[i].index].Vel[1];
  in->Vel[2] = P[StarParticle[i].index].Vel[2];

  in->Hsml              = STP(StarParticle[i].index).Hsml;
  in->NormSph           = StarParticle[i].NormSph;
  in->TotalMassReleased = StarParticle[i].TotalMassReleased;

#ifdef SMUGGLE_MASS_WEIGHT_SN
  in->TotNgbMass = StarParticle[i].TotNgbMass;
#endif
#ifdef SMUGGLE_OMEGA_WEIGHT_SN
  in->TotSolidAngle = StarParticle[i].TotSolidAngle;
#endif
#ifdef GFM_WINDS_LOCAL
  in->WindEnergyReleased = StarParticle[i].WindEnergyReleased;
#endif
#ifdef GFM_STELLAR_FEEDBACK
  in->SNIaEnergyReleased  = StarParticle[i].SNIaEnergyReleased;
  in->AGBMomentumReleased = StarParticle[i].AGBMomentumReleased;
#endif

#ifdef SMUGGLE_STAR_FEEDBACK
  in->TotalEnergyReleased = StarParticle[i].TotalEnergyReleased;

#ifdef SMUGGLE_SN_COOLING_RADIUS_BOOST
  in->n_SNII                = StarParticle[i].NumSNII;
  in->n_SNIa                = StarParticle[i].NumSNIa;
  in->TotalMassReleasedSNII = StarParticle[i].TotalMassReleasedSNII;
  in->TotalMassReleasedSNIa = StarParticle[i].TotalMassReleasedSNIa;
  in->TotalMassReleasedAGB  = StarParticle[i].TotalMassReleasedAGB;
#endif

  in->LocISMdens  = StarParticle[i].LocISMdens;  /* local ISM density (code units) */
  in->LocISMdensH = StarParticle[i].LocISMdensH; /* local ISM H density (code units) */
  in->LocISMmet   = StarParticle[i].LocISMmet;   /* local ISM metallicity          */

  in->FeedbackRadiusLimiter = StarParticle[i].FeedbackRadiusLimiter;

  in->NumNgb = StarParticle[i].NumNgb;

#ifdef SMUGGLE_AGB_WINDS
  in->AGBWindSpeed = StarParticle[i].AGBWindSpeed;
#endif

#endif

  in->Firstnode = firstnode;
}

static void out2particle(data_out *out, int i, int mode)
{
#ifdef SMUGGLE_STAR_FEEDBACK
  if(mode == MODE_LOCAL_PARTICLES)
    {
      /* the first 2 entries are not done if not local bacause they are intrinsic of a star particle */
      StarParticle[i].TotalEnergyInjected   = out->TotalEnergyInjected;
      StarParticle[i].TotalMomentumReleased = out->TotalMomentumReleased;
      StarParticle[i].TotalMomentumInjected = out->TotalMomentumInjected;
#ifdef SMUGGLE_AGB_WINDS
      StarParticle[i].TotalMomentumInjectedAGB = out->TotalMomentumInjectedAGB;
#endif
    }
  else
    {
      StarParticle[i].TotalMomentumInjected += out->TotalMomentumInjected;
#ifdef SMUGGLE_AGB_WINDS
      StarParticle[i].TotalMomentumInjectedAGB += out->TotalMomentumInjectedAGB;
#endif
    }
#endif
}

#include "../generic_comm_helpers2.h"

static void kernel_local(void)
{
  int i;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif
#pragma omp parallel private(idx)
  {
    int j, threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              if(generic_polling_primary(count, Nstar))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= Nstar)
          break;

        if(is_doing_stellar_feedback(i))
          stellar_feedback_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}

static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, cnt = 0;
#pragma omp parallel private(i)
  {
    int threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    while(1)
      {
#pragma omp atomic capture
        i = cnt++;

        if(i >= Nimport)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              generic_polling_secondary();
          }

        count++;
#endif

        stellar_feedback_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

/* main routine for doing explicit stellar feedback for ISM model */
void do_stellar_feedback(void)
{
  long long ntot;

  sumup_large_ints(1, &Nstar, &ntot);
  if(ntot == 0)
    return;

  TIMER_STOPSTART(CPU_GFM_ENRICH, CPU_GFM_FEEDBACK);

#ifdef GFM_WINDS_LOCAL
  for(int i = 0; i < NumGas; i++)
    SphP[i].WindEnergyReceived = 0;
#endif

  mpi_printf("GFM_FEEDBACK: Begin stellar feedback calculation.\n");

  generic_set_MaxNexport();

  double t0 = second();

  generic_comm_pattern(Nstar, kernel_local, kernel_imported);

#ifdef SMUGGLE_STAR_FEEDBACK
#ifdef SMUGGLE_OUTPUT_STELLAR_FEEDBACK
  for(int i = 0; i < Nstar; i++)
    {
      STP(StarParticle[i].index).FeedbackEnergy      = StarParticle[i].TotalEnergyInjected / All.cf_atime / All.cf_atime;
      STP(StarParticle[i].index).FeedbackMomentum    = StarParticle[i].TotalMomentumInjected / All.cf_atime;
      STP(StarParticle[i].index).FeedbackMomentumAGB = StarParticle[i].TotalMomentumInjectedAGB / All.cf_atime;
      STP(StarParticle[i].index).Cum_FeedbackEnergy += StarParticle[i].TotalEnergyInjected / All.cf_atime / All.cf_atime;
      STP(StarParticle[i].index).Cum_FeedbackMomentum += StarParticle[i].TotalMomentumReleased / All.cf_atime;
      STP(StarParticle[i].index).Cum_InjFeedbackMomentum += StarParticle[i].TotalMomentumInjected / All.cf_atime;
      STP(StarParticle[i].index).Cum_InjFeedbackMomentumAGB += StarParticle[i].TotalMomentumInjectedAGB / All.cf_atime;
    }
#endif
#endif

  double t1 = second();

  mpi_printf("GFM_FEEDBACK: stellar feedback calculation done took %g sec\n", timediff(t0, t1));

  TIMER_STOPSTART(CPU_GFM_FEEDBACK, CPU_GFM_ENRICH);
}

int stellar_feedback_evaluate(int target, int mode, int thread_id)
{
  int numnodes, *firstnode;
  data_in local, *in;
  data_out out;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      in = &local;

      numnodes  = 1;
      firstnode = NULL;
    }
  else
    {
      in = &DataGet[target];
      generic_get_numnodes(target, &numnodes, &firstnode);
    }

#ifdef SMUGGLE_STAR_FEEDBACK
  out.TotalEnergyInjected      = 0.0;
  out.TotalMomentumReleased    = 0.0;
  out.TotalMomentumInjected    = 0.0;
  out.TotalMomentumInjectedAGB = 0.0;
#endif

#ifdef GFM_STELLAR_FEEDBACK
  GFM_stellar_feedback(target, mode, thread_id, numnodes, firstnode, in);
#endif
#ifdef GFM_WINDS_LOCAL
  GFM_winds_local(target, mode, thread_id, numnodes, firstnode, in);
#endif
#ifdef SMUGGLE_STAR_FEEDBACK
#ifdef SMUGGLE_SN_COOLING_RADIUS_BOOST
  cooling_radius_momentum_feedback(target, mode, thread_id, numnodes, firstnode, in, &out);
#endif
#endif

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif
