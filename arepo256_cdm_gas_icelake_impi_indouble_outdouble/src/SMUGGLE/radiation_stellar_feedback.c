/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 * 
 * \file        src/SMUGGLE/radiation_Stellar_feedback.c
 * \date        03/2020
 * \author      Federico Marinacci
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef SMUGGLE_RADIATION_FEEDBACK

typedef struct
{
  MyDouble Pos[3];
  MyFloat NormSphRadFeedback;
  MyDouble RadiationMomentumReleased;
  MyFloat StromgrenRadius;
  MyFloat RadCoolShutoffTime;
#if defined(SMUGGLE_OMEGA_WEIGHT_SN)
  MyFloat TotSolidAngle;
#endif
#ifdef SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION
  MyFloat NormSphRadFeedback_cold;
  MyFloat StromgrenMass;
  MyFloat Hsml;
  MyFloat Lum;
  MyFloat LowestDensityDirection_x;
  MyFloat LowestDensityDirection_y;
  MyFloat LowestDensityDirection_z;
#endif
  MyFloat FeedbackRadiusLimiter;
  int Firstnode;
} data_in;

static data_in *DataGet;

typedef struct
{
  MyDouble RadTotalMomentumInjected;
  int PhotoionizationEvents;
} data_out;

static data_out *DataResult;

static void particle2in(data_in *in, int i, int firstnode)
{
  in->Pos[0] = P[StarParticle[i].index].Pos[0];
  in->Pos[1] = P[StarParticle[i].index].Pos[1];
  in->Pos[2] = P[StarParticle[i].index].Pos[2];

  in->FeedbackRadiusLimiter = StarParticle[i].FeedbackRadiusLimiter;

#if defined(SMUGGLE_OMEGA_WEIGHT_SN)
  // in->TotSolidAngle = StarParticle[i].TotSolidAngle;
  in->TotSolidAngle = StarParticle[i].NormSph;
#endif
  in->NormSphRadFeedback = StarParticle[i].NormSphRadFeedback;
#ifdef SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION
  in->NormSphRadFeedback_cold = StarParticle[i].NormSphRadFeedback_cold;
#endif
  in->RadiationMomentumReleased = StarParticle[i].RadiationMomentumReleased;
  in->StromgrenRadius           = StarParticle[i].StromgrenRadius;
  in->RadCoolShutoffTime        = StarParticle[i].RadCoolShutoffTime;

#ifdef SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION
  MyFloat strmass   = StarParticle[i].StromgrenMass;
  in->StromgrenMass = strmass;
  in->Lum           = StarParticle[i].Lum;
  StarParticle[i].RadFeedTau =
      StarParticle[i].GasColumnDensity * All.DustOpacityRadiationFeedback * StarParticle[i].LocISMmet / GFM_SOLAR_METALLICITY;
  in->RadiationMomentumReleased *= (1. + StarParticle[i].RadFeedTau);
  in->Hsml = STP(StarParticle[i].index).Hsml;
#endif

  in->Firstnode = firstnode;
}

static void out2particle(data_out *out, int i, int mode)
{
#ifdef SMUGGLE_RADIATION_FEEDBACK_DEBUG
  if(mode == MODE_LOCAL_PARTICLES)
    {
      STP(StarParticle[i].index).Cum_RadMomentumRealInjected += out->RadTotalMomentumInjected;
    }
  else
    {
      STP(StarParticle[i].index).Cum_RadMomentumRealInjected += out->RadTotalMomentumInjected;
    }

  STP(StarParticle[i].index).PhotoionizationEvents += out->PhotoionizationEvents;
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

    for(j                            = 0; j < NTask; j++)
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

#ifdef SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION
        if((StarParticle[i].StromgrenMass > 0.0) && (StarParticle[i].RadiationMomentumReleased > 0.0) &&
           (StarParticle[i].NormSphRadFeedback > 0))
          {
            radiation_feedback_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
            STP(StarParticle[i].index).PhotoionizationAttempts += 1;
          }
#endif
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

        radiation_feedback_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

#ifdef SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION

/* primary routine called from main run loop to do radiation stellar feedback
 * for explicit ISM model */
void do_radiation_stellar_feedback(void)
{
  long long ntot;
  sumup_large_ints(1, &Nstar, &ntot);

  if(ntot == 0)
    return;

  mpi_printf(
      "SMUGGLE_RADIATION_FEEDBACK: Begin radiation stellar feedback "
      "calculation.\n");

  generic_set_MaxNexport();
  double t0 = second();
  generic_comm_pattern(Nstar, kernel_local, kernel_imported);
  double t1 = second();

#ifdef SMUGGLE_RADIATION_FEEDBACK_DEBUG
  for(int i                               = 0; i < Nstar; i++)
    STP(StarParticle[i].index).RadFeedTau = StarParticle[i].RadFeedTau;
#endif

  mpi_printf("SMUGGLE_RADIATION_FEEDBACK: stellar feedback calculation took %g sec\n", timediff(t0, t1));
}

int radiation_feedback_evaluate(int target, int mode, int thread_id)
{
  int numnodes, *firstnode;
  data_in local, *in;
  data_out out;

  double h, h2, hinv, hinv3;
  double dx, dy, dz, r2, r;
  MyDouble *pos;
  MyDouble de_feedback;

  MyDouble RadiationMomentumReleased;
  MyFloat RadCoolShutoffTime;
  MyFloat tot_mass;
  MyDouble dp_feedback, inj_mom[3];
  MyDouble dp_tot, du;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      in        = &local;
      numnodes  = 1;
      firstnode = NULL;
    }
  else
    {
      in = &DataGet[target];
      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  pos = in->Pos;
  h   = in->Hsml;

  out.PhotoionizationEvents    = 0;
  out.RadTotalMomentumInjected = 0;

  MyFloat strom_mass = in->StromgrenMass;                                   /* mass calculated in
                                                                               src/SMUGGLE/radiation_stellar_feedback_util.c */
  MyFloat pi_kernel_mass = in->NormSphRadFeedback_cold;                     /* amount of *cold* gas mass in kernel,
                                                                               available for photoionization */
  RadiationMomentumReleased = in->RadiationMomentumReleased * All.cf_atime; /* now in comoving code units;
                                                                               copied from AGBMomentum in
                                                                               GFM/stellar_feedbac.c */
  RadCoolShutoffTime = in->RadCoolShutoffTime;                              /* time in physical code units */

  MyFloat lum = in->Lum;

#ifdef SMUGGLE_OMEGA_WEIGHT_SN
  MyFloat totsolidangle = in->TotSolidAngle;
#endif

  if(in->NormSphRadFeedback <= 0 && h > 0)
    terminate(
        "tot_mass<=0 should not happen in radiation_stellar_feedback "
        "mode 0: pos %g|%g|%g, tot_mass %g, NumNgb %g, h %g , "
        "RadiationMomentumReleased %g,  ID = %lld\n",
        pos[0], pos[1], pos[2], tot_mass, StarParticle[target].NumNgb, h, 
        RadiationMomentumReleased, (long long) P[StarParticle[target].index].ID);

  h2   = h * h;
  hinv = 1.0 / h;
#ifndef TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif

  dp_tot = 0.0;

  double p_ionize = strom_mass / pi_kernel_mass; /* fraction of *cold* gas mass in the kernel that should
                                                    be ionized */
  double r2lim = in->FeedbackRadiusLimiter * in->FeedbackRadiusLimiter;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  for(int n = 0; n < nfound; n++)
    {
      int j = Thread[thread_id].Ngblist[n];
      if(P[j].Mass > 0 && P[j].ID != 0) /* skip cells that have been swallowed or dissolved */
        {
          dx = NEAREST_X(pos[0] - P[j].Pos[0]);
          dy = NEAREST_Y(pos[1] - P[j].Pos[1]);
          dz = NEAREST_Z(pos[2] - P[j].Pos[2]);

          r2 = dx * dx + dy * dy + dz * dz;

          // if((r2 < h2) && (r2 < r2lim))
          if(r2 < h2)
            {
              r           = sqrt(r2);
              de_feedback = 0.;

#ifdef SMUGGLE_STAR_FEEDBACK
#ifdef BH_BASED_CGM_ZOOM
              if(P[j].Mass < 0.3 * All.TargetGasMass / All.CGM_RefinementFactor)
                continue;
#else
              if(P[j].Mass < 0.3 * All.TargetGasMass)
                continue;
#endif
#endif

#ifdef SMUGGLE_OMEGA_WEIGHT_SN
              double cell_radius = get_cell_radius(j);
              double cell_area   = cell_radius * cell_radius;
              double omega       = 0.5 * (1.0 - sqrt(r2 / (r2 + cell_area)));
              double weight_fac  = omega / totsolidangle;
#endif

              /* First input momentum if applicable */
              if(RadiationMomentumReleased > 0)
                {
                  double Ekin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                                       SphP[j].Momentum[2] * SphP[j].Momentum[2]) /
                                P[j].Mass;
                  dp_feedback = weight_fac * RadiationMomentumReleased;

                  /* injected feedback momentum radially away */
                  inj_mom[0] = -dp_feedback * dx / r;
                  inj_mom[1] = -dp_feedback * dy / r;
                  inj_mom[2] = -dp_feedback * dz / r;

                  /* momentum due to stellar mass return is injected in radial direction
                   */
                  SphP[j].Momentum[0] += inj_mom[0];
                  SphP[j].Momentum[1] += inj_mom[1];
                  SphP[j].Momentum[2] += inj_mom[2];

                  de_feedback = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                                       SphP[j].Momentum[2] * SphP[j].Momentum[2]) /
                                    P[j].Mass -
                                Ekin;

                  dp_tot += dp_feedback; /* LVS cumulative on each stellar part */
                }                        /* if RadMom > 0 */

              double Ekin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                                   SphP[j].Momentum[2] * SphP[j].Momentum[2]) /
                            P[j].Mass;

              double Utherm = (SphP[j].Energy - Ekin) / (All.cf_atime * All.cf_atime * P[j].Mass);

              /* identify gas that could be photoionized, but make sure we can
               * "reselect" gas that has just returned to non-photoionized state */
              if((Utherm < 1.2 * All.PhotoionizationEgySpec) && (SphP[j].GasRadCoolShutoffTime <= 0.0))
                // if((Utherm < 1.2 * All.PhotoionizationEgySpec))
                {
                  double alpha_rec = 2.6e-13;
                  double s = weight_fac * lum / (All.RadiationFeedbackAvgPhotonEnergyineV * ELECTRONVOLT_IN_ERGS); /* phot/sec */

                  double gas_dens = SphP[j].Density * All.HubbleParam * All.HubbleParam * All.UnitDensity_in_cgs * All.cf_a3inv;
                  double nh       = gas_dens / PROTONMASS * SphP[j].MetalsFraction[element_index("Hydrogen")];
                  double rec = SphP[j].MetalsFraction[element_index("Hydrogen")] * nh * alpha_rec * P[j].Mass * All.UnitMass_in_g /
                               All.HubbleParam / PROTONMASS;
                  p_ionize = s / rec;

                  if(get_random_number() < p_ionize)
                    {
                      du = fmax(All.PhotoionizationEgySpec - Utherm, 0.0);         /* physical */
                      de_feedback += All.cf_atime * All.cf_atime * du * P[j].Mass; /* add *comoving* change in temperature to energy */

                      SphP[j].GasRadCoolShutoffTime = RadCoolShutoffTime;
                      out.PhotoionizationEvents += 1;
                    }
                }
              SphP[j].Energy += de_feedback; /* de_feedback has a momentum and heat
                                                component and in comoving code units*/
            }                                /* if r2<rs2) */
        }                                    /* if P[j].ID != 0 */
    }

  out.RadTotalMomentumInjected = dp_tot / All.cf_atime; /* make it in physical units for logging purposes */

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif
#endif
