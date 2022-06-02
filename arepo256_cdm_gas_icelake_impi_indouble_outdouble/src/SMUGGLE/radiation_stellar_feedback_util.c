/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 * 
 * \file        src/SMUGGLE/radiation_stellar_feedback_util.c
 * \date        03/2020
 * \author      Federico Marinacci
 * \brief       routines for radiation stellar feedback
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <gsl/gsl_math.h>
#include <math.h>
#include "../allvars.h"
#include "../proto.h"

#ifdef SMUGGLE_RADIATION_FEEDBACK

/* compute the Stromgren Radius using the average density within kernel of star
 * particle */
/* feed though physical units, and return physical units */
MyFloat compute_stromgren_radius(double lum, double avg_gas_dens, double minGasDist, double hsml)
{
  double s, rs, nh;
  double alpha_rec = 2.6e-13; /* recombination coefficient of ionized gas at T=1.e4 --LVS */
  MyFloat rsfinal;

  s = lum / (All.RadiationFeedbackAvgPhotonEnergyineV * ELECTRONVOLT_IN_ERGS); /* phot/sec */

  avg_gas_dens *= All.HubbleParam * All.HubbleParam * All.UnitDensity_in_cgs;
  nh = avg_gas_dens / PROTONMASS * HYDROGEN_MASSFRAC;

  rs = 3. * s / (4. * 3.14159 * alpha_rec * nh * nh);
  rs = pow(rs, 0.333);
  rs *= All.HubbleParam / All.UnitLength_in_cm; /* back to physical code units */

  rsfinal = fmax(rs, minGasDist * 1.5);

  rsfinal = fmin(rsfinal, hsml); /* rsfinal and hsml in physical units */

  return rsfinal;
}

#ifdef SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION
/* compute the mass in the Stromgren Radius using the average density within
 * kernel of star particle */
/* feed through physical units, and return physical (mass) */
MyFloat compute_stromgren_mass(double lum, double avg_gas_dens, double hsml)
{
  double rs = compute_stromgren_radius(lum, avg_gas_dens, 0.0, hsml);
  return 4.18878667 * rs * rs * rs * avg_gas_dens; /* stromgren mass in code units */
}
#endif  /* SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION */

/* This function returns 0 if the cell cannot cool, otherwise returns the time
 * for which cooling is allowed */
MyFloat update_radiation_cooling_shutoff_time(int i, MyFloat dt)
{
  MyFloat dtcool = 0.0;

  SphP[i].GasRadCoolShutoffTime -= dt;  // particle has been evolved for a time step

  if(SphP[i].GasRadCoolShutoffTime < 0.0)  // cooling can happen again
    {
      // dtcool = -SphP[i].GasRadCoolShutoffTime;  // part of the time step for
      // which cooling is allowed
      SphP[i].GasRadCoolShutoffTime = 0;  // reset the value so that cell can undergo normal cooling
    }

  return dtcool;
}
#endif
