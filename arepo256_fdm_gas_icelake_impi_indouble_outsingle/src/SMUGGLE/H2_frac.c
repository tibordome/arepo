/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 * 
 * \file        src/SMUGGLE/H2_frac.c
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

#include <math.h>
#include "../allvars.h"

#ifdef SMUGGLE_COMPUTE_SFR_FROM_H2

//#define SOLAR_METALLICITY 0.0127        /* the same as GFM */

/* computes the molecular Hydrogen mass fraction see Kuhlen et al. 2013 */
double MolecularHFrac(int i)
{
  double fH2 = 1.0;
  double chi, s;
  double density, metallicity;
  double gradD, sigma, tau;

  density = SphP[i].Density * All.cf_a3inv;

#ifdef GFM_STELLAR_EVOLUTION
  metallicity = SphP[i].Metallicity / GFM_SOLAR_METALLICITY;

  if(metallicity < 1.0e-5)
    metallicity = 1e-5;
#else
  metallicity = 1.;
#endif

  /* There is also the possibility of obtaining the cell radius with the
     function get_cell_radius(i) defined in voronoi.c */

  /*  Alternative for the column mass density \Sigma in Kuhlen+2011 */
  /*
  radius = pow(0.75 * SphP[i].Volume / (a3inv * M_PI), 1./3);
  sigma = density * All.HubbleParam * radius;
  tau = 0.067 * ((All.UnitMass_in_g / SOLAR_MASS) * (pow(PARSEC /
  All.UnitLength_in_cm, 2.0))) * metallicity / SOLAR_METALLICITY *
  sigma;
  */

  /*  Alternative for the column density \Sigma in Gnedin, Tassi & Kravtsov 2009
   */

  gradD = sqrt(SphP[i].Grad.drho[0] * SphP[i].Grad.drho[0] + SphP[i].Grad.drho[1] * SphP[i].Grad.drho[1] +
               SphP[i].Grad.drho[2] * SphP[i].Grad.drho[2]);

  /* to physical units */
  gradD /= All.cf_atime * All.cf_atime * All.cf_atime * All.cf_atime;

  if(gradD > 0)
    {
      /* see also McKee & Krumholz (2010, pag. 319)*/
      sigma = density * All.HubbleParam * (density / gradD);
      tau   = 0.067 * ((All.UnitMass_in_g / SOLAR_MASS) * (pow(PARSEC / All.UnitLength_in_cm, 2.0))) * metallicity * sigma;

      /* see also Krumholz, McKee & Tumlinson (2009) */
      chi = 2.3 * (1.0 + 3.1 * pow(metallicity, 0.365)) / 3.0;
      s   = log(1.0 + chi * (0.6 + 0.01 * chi)) / (0.6 * tau);
    }
  else
    s = 0.0;

  fH2 = 1.0 - 0.75 * (s / (1.0 + 0.25 * s));

  if(fH2 < 0.0)
    fH2 = 0.0;

#ifdef SMUGGLE_OUTPUT_OPTICAL_DEPTH
  SphP[i].OpticalDepth = tau;
#endif

  return fH2;
}

#endif /* closes SMUGGLE_COMPUTE_SFR_FROM_H2 */
