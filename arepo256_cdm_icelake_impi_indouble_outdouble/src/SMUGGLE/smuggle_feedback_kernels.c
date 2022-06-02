/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 * 
 * \file        src/SMUGGLE/smuggle_feedback_kernels.c
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

#include "../allvars.h"
#include "../proto.h"
#include "smuggle_feedback_kernels.h"

#ifdef SMUGGLE_STAR_FEEDBACK

#ifdef SMUGGLE_SN_COOLING_RADIUS_BOOST
void cooling_radius_momentum_feedback(int target, int mode, int thread_id, int numnodes, int *firstnode, data_in *in, data_out *out)
{
  int n, j, nfound;
  double weight_fac, dr[3], r2, r, h, h2, normsph;
#ifdef SMUGGLE_MASS_WEIGHT_SN
  double wk, u, h3, hinv, hinv3;
#endif
  double n_SNII, n_SNIa;

  // double r2lim = All.FeedbackRadiusLimiter * All.FeedbackRadiusLimiter;
  double rlim = in->FeedbackRadiusLimiter;

  MyDouble *pos;
  MyDouble *vel;
  MyDouble TotalEnergyReleased, TotalMomentumReleased, TotalMassReleased, de, dpSN, dp;
  MyDouble dpSNv[3];
  MyDouble dptot = 0;

  MyDouble TotalMassReleasedAGB, dpAGB;

#ifdef SMUGGLE_AGB_WINDS
  MyDouble TotalEnergyReleasedAGB, TotalMomentumReleasedAGB;
  MyDouble dptotAGB = 0;
  MyDouble dpAGBv[3];
#endif

  pos = in->Pos;
  vel = in->Vel;
  h   = in->Hsml;

  normsph = in->NormSph;
#ifdef SMUGGLE_MASS_WEIGHT_SN
  MyFloat totngbmass = in->TotNgbMass;
#endif
#ifdef SMUGGLE_OMEGA_WEIGHT_SN
  MyFloat totsolidangle = in->TotSolidAngle;
#endif
  n_SNII = in->n_SNII;
  n_SNIa = in->n_SNIa;

  TotalMassReleasedAGB = in->TotalMassReleasedAGB;

  if(TotalMassReleasedAGB < 0.0)
    terminate("Negative mass released by AGB stars %g\n", TotalMassReleasedAGB);

#ifdef SMUGGLE_AGB_WINDS
  TotalEnergyReleasedAGB = 0.5 * TotalMassReleasedAGB * in->AGBWindSpeed * in->AGBWindSpeed * All.cf_atime *
                           All.cf_atime;                                             /* comoving total energy released */
  TotalMomentumReleasedAGB = TotalMassReleasedAGB * in->AGBWindSpeed * All.cf_atime; /* comoving total momentum released */
#endif

  /* nothing to be done for this star particle */
  if((n_SNII <= 0.0) && (n_SNIa <= 0.0) && (TotalMassReleasedAGB <= 0.0))
    return;

  double SNII_velocity = 0.0;
  if(in->TotalMassReleasedSNII > 0.0)
    SNII_velocity = sqrt(2.0 * n_SNII * All.one_SNe_energy / in->TotalMassReleasedSNII); /* physical velocity in code units */

  double SNIa_velocity = 0.0;
  if(in->TotalMassReleasedSNIa > 0.0)
    SNIa_velocity = sqrt(2.0 * n_SNIa * All.one_SNe_energy / in->TotalMassReleasedSNIa); /* physical velocity in code units */

  TotalMassReleased     = in->TotalMassReleasedSNII + in->TotalMassReleasedSNIa;
  TotalEnergyReleased   = (n_SNII + n_SNIa) * All.one_SNe_energy * All.cf_atime * All.cf_atime; /* comoving total energy released */
  TotalMomentumReleased = (in->TotalMassReleasedSNII * SNII_velocity + in->TotalMassReleasedSNIa * SNIa_velocity) *
                          All.cf_atime; /* comoving total momentum released */

  if(TotalMassReleased < 0.0)
    terminate("Negative mass released by SN %g\n", TotalMassReleased);

  /* Expression for terminal momentum in Hopkins+ 2017, this is the maximum
   * boost momentum can have */
  MyFloat terminal_mom =
      4.8e10 * SOLAR_MASS * All.HubbleParam / (All.UnitMass_in_g * All.UnitVelocity_in_cm_per_s); /* physical code units  */
  terminal_mom *= (n_SNII + n_SNIa) * pow(All.FeedbackEfficiency, 13. / 14.) * All.cf_atime;      /* comoving code units  */

  /* properties of the local ISM for cooling radius calc */
  /* to avoid zero values: facZ saturates at 2 below about 0.01 Zsun anyway */
  double z0   = fmax(in->LocISMmet / GFM_SOLAR_METALLICITY, 1e-5);
  double facZ = pow(fmin(2.0, pow(z0, -0.14)), 1.5);
  double n0   = in->LocISMdensH * All.cf_a3inv / PROTONMASS * All.UnitDensity_in_cgs * All.HubbleParam *
              All.HubbleParam; /* physical density in cgs */
  double facn0 = pow(n0, -1. / 7.);
  double rcool = 28.4 * PARSEC * facZ * pow(n0, -3. / 7.) * (n_SNII + n_SNIa) * pow(All.FeedbackEfficiency, 2. / 7.) /
                 All.UnitLength_in_cm * All.HubbleParam / All.cf_atime;

  terminal_mom *= (facZ * facn0);

  h2 = h * h;

  nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[thread_id].Ngblist[n];
      if(P[j].Mass > 0 && P[j].ID != 0) /* skip cells that have been swallowed or dissolved */
        {
          dr[0] = NEAREST_X(pos[0] - P[j].Pos[0]);
          dr[1] = NEAREST_Y(pos[1] - P[j].Pos[1]);
          dr[2] = NEAREST_Z(pos[2] - P[j].Pos[2]);
          r2    = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];

          if(r2 < h2)
            {
              r = sqrt(r2);
#ifdef SMUGGLE_MASS_WEIGHT_SN
              hinv  = 1.0 / h;
              hinv3 = hinv * hinv * hinv;
              h3    = 1.0 / hinv3;
              u     = r * hinv;
              if(u < 0.5)
                wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
              else
                wk       = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
              weight_fac = NORM_COEFF * P[j].Mass * wk * h3 / totngbmass;
#elif defined(SMUGGLE_OMEGA_WEIGHT_SN)
              double cell_radius = get_cell_radius(j);
              double cell_area   = cell_radius * cell_radius;
              double omega       = 0.5 * (1.0 - sqrt(r2 / (r2 + cell_area)));
              weight_fac         = omega / totsolidangle;
              double mweight_fac = omega / normsph;
#else
              weight_fac = SphP[j].Volume / normsph;
#endif

#ifdef BH_BASED_CGM_ZOOM
              if(P[j].Mass < 0.3 * All.TargetGasMass / All.CGM_RefinementFactor)
                continue;
#else
              if(P[j].Mass < 0.3 * All.TargetGasMass)
                continue;
#endif

              double coolfac = 1.0;

              if(rcool > 0.0)
                coolfac = fmin(1., pow(r / rcool, -6.5));

              double boost_fac  = 0.0;
              double bf[2]      = {0.0, 0.0};
              double inj_mom[3] = {0.0, 0.0, 0.0};
              double inj_mass   = 0.0;
              dpSN              = 0.0;
              de                = 0.0;
              dp                = 0.0;

              if(TotalEnergyReleased > 0)
                {
                  if(r < rlim)
                    {
                      // new plan, compte all in SF star then convert to lab frame
                      dpSNv[0] = (-dr[0] * TotalMomentumReleased / r) * weight_fac; /* dp of a single cell in comoving units */
                      dpSNv[1] = (-dr[1] * TotalMomentumReleased / r) * weight_fac; /* dp of a single cell in comoving units */
                      dpSNv[2] = (-dr[2] * TotalMomentumReleased / r) * weight_fac; /* dp of a single cell in comoving units */

                      double momc[3];
                      momc[0] = SphP[j].Momentum[0] - vel[0] * P[j].Mass;
                      momc[1] = SphP[j].Momentum[1] - vel[1] * P[j].Mass;
                      momc[2] = SphP[j].Momentum[2] - vel[2] * P[j].Mass;

                      dpSN = dpSNv[0] * dpSNv[0] + dpSNv[1] * dpSNv[1] + dpSNv[2] * dpSNv[2];

                      double msq   = momc[0] * momc[0] + momc[1] * momc[1] + momc[2] * momc[2];
                      double cdot  = momc[0] * dpSNv[0] + momc[1] * dpSNv[1] + momc[2] * dpSNv[2];
                      double massr = P[j].Mass / (mweight_fac * TotalMassReleased);

                      boost_fac = fmin(sqrt(1. + massr), terminal_mom / TotalMomentumReleased);

                      double sol, boost_fac_max;
                      if(cdot > 0.0)
                        {
                          sol           = cdot / dpSN + sqrt(1. + massr + msq / (dpSN * massr) + cdot * cdot / (dpSN * dpSN));
                          boost_fac_max = (1. + massr + msq / (dpSN * massr)) / sol;
                        }
                      else
                        boost_fac_max = -cdot / dpSN + sqrt(1. + massr + msq / (dpSN * massr) + cdot * cdot / (dpSN * dpSN));

                      boost_fac = fmax(fmin(boost_fac, 0.99 * boost_fac_max), 1.0);
                      bf[0]     = boost_fac;
                      bf[1]     = boost_fac_max;

                      dpSN = sqrt(dpSN);
                      de   = 0.5 * dpSN * dpSN / (TotalMassReleased * mweight_fac);

                      dpSNv[0] *= boost_fac;
                      dpSNv[1] *= boost_fac;
                      dpSNv[2] *= boost_fac;

                      dp = sqrt(dpSNv[0] * dpSNv[0] + dpSNv[1] * dpSNv[1] + dpSNv[2] * dpSNv[2]);

                      inj_mom[0] = dpSNv[0];
                      inj_mom[1] = dpSNv[1];
                      inj_mom[2] = dpSNv[2];
                    }
                }

              inj_mass = TotalMassReleased * mweight_fac;  // mass is always injected but not feedback

              if(inj_mass > 0.0)
                smuggle_inject_snfeed_into_cell(j, inj_mass, de, inj_mom, bf, vel, coolfac);

              dpAGB      = 0.0;
              de         = 0.0;
              inj_mass   = 0.0;
              inj_mom[0] = 0.0;
              inj_mom[1] = 0.0;
              inj_mom[2] = 0.0;
              bf[0]      = weight_fac;
              bf[1]      = mweight_fac;

#ifdef SMUGGLE_AGB_WINDS
              if(TotalEnergyReleasedAGB > 0.0)
                {
                  if(r < rlim)
                    {
                      dpAGBv[0] = (-dr[0] * TotalMomentumReleasedAGB / r) * weight_fac;
                      dpAGBv[1] = (-dr[1] * TotalMomentumReleasedAGB / r) * weight_fac;
                      dpAGBv[2] = (-dr[2] * TotalMomentumReleasedAGB / r) * weight_fac;
                      dpAGB     = sqrt(dpAGBv[0] * dpAGBv[0] + dpAGBv[1] * dpAGBv[1] + dpAGBv[2] * dpAGBv[2]);

                      de = 0.5 * dpAGB * dpAGB / (TotalMassReleasedAGB * mweight_fac);

                      inj_mom[0] = dpAGBv[0];
                      inj_mom[1] = dpAGBv[1];
                      inj_mom[2] = dpAGBv[2];
                    }
                }
#endif

              inj_mass = TotalMassReleasedAGB * mweight_fac;  // mass is always injected, but not feedback

              if(inj_mass > 0.0)
                smuggle_inject_windfeed_into_cell(j, inj_mass, de, inj_mom, bf, vel);

              dptot += dp;
#ifdef SMUGGLE_AGB_WINDS
              dptotAGB += dpAGB;
#endif
            }
        }
    }

  out->TotalEnergyInjected   = TotalEnergyReleased;
  out->TotalMomentumReleased = TotalMomentumReleased;
  out->TotalMomentumInjected = dptot;
#ifdef SMUGGLE_AGB_WINDS
  out->TotalMomentumInjectedAGB = dptotAGB;
#endif
}
#endif

#endif
