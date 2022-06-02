/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 * 
 * \file        src/SMUGGLE/stellar_feedback_util.c
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

#include <gsl/gsl_randist.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef SMUGGLE_STAR_FEEDBACK

#ifndef SMUGGLE_SFR
#warning SMUGGLE_STAR_FEEDBACK requires SMUGGLE_SFR
#endif

#define SN_ENERGY 1.0e51 /* erg */
#define CM_PER_KM 1.0e5
#define SN_MOMENTUM (1.0e10 * SOLAR_MASS) /* gr cm/s */

/*! \file fm_stellar_feedback.c
 *  \brief routines for stellar feedback
 */
MyDouble compute_SN_energy(MyDouble Number_of_SNae)
{
  MyDouble Injected_Energy =
      All.FeedbackEfficiency * Number_of_SNae * SN_ENERGY;      /* a fraction FeedbackEfficiency is given to the gas */
  Injected_Energy *= (All.HubbleParam / All.UnitEnergy_in_cgs); /* to code units */

  return Injected_Energy;
}

/*  \param b = <e_r, P_i>,
 *  \param c = 2[M\DeltaE + \deltaM E_f] = 2J
 *  E_f is the final kinetic energy \DeltaE the feedback energy
 *  \return the magnitude of the injected momentum in radial direction
 *          such that the kinetic energy of the cell is increased by \DeltaE
 */
MyDouble quadratic_equation_solver(MyDouble b, MyDouble c)
{
  MyDouble res;

  if(b <= 0.0)
    res = -(b - sqrt(b * b + c));
  else
    res = c / (b + sqrt(b * b + c));

  return res;
}

double get_feedback_radius_limiter(double gasdens, double metallicity, double initial_mass)
{
#ifdef SMUGGLE_SUPERBUBBLE_LIMITER
  double n0 = fmax(gasdens * All.cf_a3inv / PROTONMASS * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam,
                   0.5);  // physical density in cgs
  double nsn = initial_mass * All.SN_per_Msun * All.UnitMass_in_g / SOLAR_MASS / All.HubbleParam;
  double tOB = get_lifetime_in_Gyr(All.SNII_MinMass_Msun, metallicity) * SEC_PER_GIGAYEAR;

  double max_rad = pow(0.76 * 125 * nsn * GFM_SNII_ENERGY * tOB * tOB / (154 * M_PI * PROTONMASS * n0), 0.2);  // Weaver+ 1977 eq. 21

  max_rad *= All.HubbleParam / (All.cf_atime * All.UnitLength_in_cm);
#else
  double max_rad = All.FeedbackRadiusLimiter / All.cf_atime;
#endif

  return max_rad;
}

#ifdef SMUGGLE_SN_COOLING_RADIUS_BOOST
void smuggle_inject_snfeed_into_cell(int j, double inj_mass, double inj_energy, double inj_mom[3], double boost[2], double vel[3],
                                     double coolfac)
{
  /* compute current kinetic energy */
  double Ekin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                       SphP[j].Momentum[2] * SphP[j].Momentum[2]) /
                P[j].Mass;

  /* computer cell thermal energy */
  double utherm = SphP[j].Energy - Ekin;

  /* transform momentum in star frame  */
  SphP[j].Momentum[0] -= P[j].Mass * vel[0];
  SphP[j].Momentum[1] -= P[j].Mass * vel[1];
  SphP[j].Momentum[2] -= P[j].Mass * vel[2];

  /* compute current kinetic energy */
  Ekin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                SphP[j].Momentum[2] * SphP[j].Momentum[2]) /
         P[j].Mass;

  /* add injected momentum to cell in star frame */
  SphP[j].Momentum[0] += inj_mom[0];
  SphP[j].Momentum[1] += inj_mom[1];
  SphP[j].Momentum[2] += inj_mom[2];

  /* add injected mass to cell */
  P[j].Mass += inj_mass;

  /* initialize variation in kinentic energy */
  double dEkin = -Ekin;

  /* compute current kinetic energy */
  Ekin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                SphP[j].Momentum[2] * SphP[j].Momentum[2]) /
         P[j].Mass;

  /* compute kinetic energy variation */
  dEkin += Ekin;

  /* compute new thermal energy */
  double du;
  if(inj_energy > 0.0)
    du = inj_energy * (1. - dEkin / inj_energy);
  else
    du = -dEkin;

  if(du < 0)
    warn(
        "negative thermal energy variation du=%g, boost %g|%g, inj_mass %g, "
        "inj_energy %g, inj_mom %g|%g|%g\n",
        du, boost[0], boost[1], inj_mass, inj_energy, inj_mom[0], inj_mom[1], inj_mom[2]);

  if(inj_energy > 0.0)
    du *= coolfac;

  utherm += du;

  /* transform momentum in lab frame */
  SphP[j].Momentum[0] += P[j].Mass * vel[0];
  SphP[j].Momentum[1] += P[j].Mass * vel[1];
  SphP[j].Momentum[2] += P[j].Mass * vel[2];

  /* compute new kinetic in lab frame */
  Ekin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                SphP[j].Momentum[2] * SphP[j].Momentum[2]) /
         P[j].Mass;

  /* compute new total energy in lab frame */
  SphP[j].Energy = Ekin + utherm;
}

void smuggle_inject_windfeed_into_cell(int j, double inj_mass, double inj_energy, double inj_mom[3], double boost[2], double vel[3])
{
  /* compute current kinetic energy */
  double Ekin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                       SphP[j].Momentum[2] * SphP[j].Momentum[2]) /
                P[j].Mass;

  /* computer cell thermal energy */
  double utherm = SphP[j].Energy - Ekin;

  /* transform momentum in star frame  */
  SphP[j].Momentum[0] -= P[j].Mass * vel[0];
  SphP[j].Momentum[1] -= P[j].Mass * vel[1];
  SphP[j].Momentum[2] -= P[j].Mass * vel[2];

  /* compute current kinetic energy */
  Ekin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                SphP[j].Momentum[2] * SphP[j].Momentum[2]) /
         P[j].Mass;

  /* add injected momentum to cell in star frame */
  SphP[j].Momentum[0] += inj_mom[0];
  SphP[j].Momentum[1] += inj_mom[1];
  SphP[j].Momentum[2] += inj_mom[2];

  /* add injected mass to cell */
  P[j].Mass += inj_mass;

  /* initialize variation in kinentic energy */
  double dEkin = -Ekin;

  /* compute current kinetic energy */
  Ekin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                SphP[j].Momentum[2] * SphP[j].Momentum[2]) /
         P[j].Mass;

  /* compute kinetic energy variation */
  dEkin += Ekin;

  /* compute new thermal energy */
  double du;
  if(inj_energy > 0.0)
    du = inj_energy * (1. - dEkin / inj_energy);
  else
    du = -dEkin;

  if(du < 0)
    warn(
        "negative thermal energy variation du=%g, weight %g|%g, inj_mass %g, "
        "inj_energy %g, inj_mom %g|%g|%g\n",
        du, boost[0], boost[1], inj_mass, inj_energy, inj_mom[0], inj_mom[1], inj_mom[2]);

  utherm += du;

  /* transform momentum in lab frame */
  SphP[j].Momentum[0] += P[j].Mass * vel[0];
  SphP[j].Momentum[1] += P[j].Mass * vel[1];
  SphP[j].Momentum[2] += P[j].Mass * vel[2];

  /* compute new kinetic in lab frame */
  Ekin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                SphP[j].Momentum[2] * SphP[j].Momentum[2]) /
         P[j].Mass;

  /* compute new total energy in lab frame */
  SphP[j].Energy = Ekin + utherm;
}

#endif

#ifdef SMUGGLE_DISCRETE_SN
double smuggle_sample_poisson_distribution(double lambda)
{
  double k = 0.0;

  if(lambda > 0.0)
    k = (double)gsl_ran_poisson(random_generator_aux, lambda);

  return k;

  /*
  * It uses the invesion of the CDF to sample a poisson distribution
  * for the determination of SN events given a uniform random number
  * generator. This only requires a generation of a single random number
  * per sampling for low values of the parameter lambda
  if(lambda > 10.)
    {
      return gsl_ran_poisson(random_generator_aux, lambda);
    }
  else
    {
      double k = 0;
      double p = exp(-lambda);
      double s = p;
      double u = get_random_number_aux();

      while(u > s)
        {
          k += 1.;
          p *= lambda / k;
          s += p;
        }

      return k;
    }
  */
}
#endif

#endif  // SMUGGLE_STAR_FEEDBACK
