/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sfr_mcs/hii_mcs_healpix_utils.c
 * \date        03/2019
 * \author      Matthew Smith
 * \brief       Healpix functions modified from chealpix.c (Gorski & Hivon 2011)
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - 17.01.2022 Ported into current codebase
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef HII_MCS_ANISO

static const double twothird   = 2.0 / 3.0;
//static const double pi         = 3.141592653589793238462643383279502884197;
static const double twopi      = 6.283185307179586476925286766559005768394;
//static const double halfpi     = 1.570796326794896619231321691639751442099;
static const double inv_halfpi = 0.6366197723675813430755350534900574;

/*! Returns the remainder of the division \a v1/v2.
    The result is non-negative.
    \a v1 can be positive or negative; \a v2 must be positive. */
static double fmodulo(double v1, double v2)
{
  if(v1 >= 0)
    return (v1 < v2) ? v1 : fmod(v1, v2);

  double tmp = fmod(v1, v2) + v2;
  return (tmp == v2) ? 0. : tmp;
}

/*! Returns the remainder of the division \a v1/v2.
   The result is non-negative.
   \a v1 can be positive or negative; \a v2 must be positive. */
static int imodulo(int v1, int v2)
{
  int v = v1 % v2;
  return (v >= 0) ? v : v + v2;
}

static int ang2pix_ring_z_phi(double z, double phi)
{
  double za = fabs(z);
  double tt = fmodulo(phi, twopi) * inv_halfpi; /* in [0,4) */

  if(za <= twothird) /* Equatorial region */
    {
      double temp1 = HII_MCS_N_SIDE * (0.5 + tt);
      double temp2 = HII_MCS_N_SIDE * z * 0.75;
      int jp       = (int)(temp1 - temp2); /* index of  ascending edge line */
      int jm       = (int)(temp1 + temp2); /* index of descending edge line */

      /* ring number counted from z=2/3 */
      int ir     = HII_MCS_N_SIDE + 1 + jp - jm; /* in {1,2n+1} */
      int kshift = 1 - (ir & 1);                 /* kshift=1 if ir even, 0 otherwise */

      int ip = (jp + jm - HII_MCS_N_SIDE + kshift + 1) / 2; /* in {0,4n-1} */
      ip     = imodulo(ip, 4 * HII_MCS_N_SIDE);

      return HII_MCS_N_SIDE * (HII_MCS_N_SIDE - 1) * 2 + (ir - 1) * 4 * HII_MCS_N_SIDE + ip;
    }
  else /* North & South polar caps */
    {
      double tp  = tt - (int)(tt);
      double tmp = HII_MCS_N_SIDE * sqrt(3 * (1 - za));

      int jp = (int)(tp * tmp);         /* increasing edge line index */
      int jm = (int)((1.0 - tp) * tmp); /* decreasing edge line index */

      int ir = jp + jm + 1;    /* ring number counted from the closest pole */
      int ip = (int)(tt * ir); /* in {0,4*ir-1} */
      ip     = imodulo(ip, 4 * ir);

      if(z > 0)
        return 2 * ir * (ir - 1) + ip;
      else
        return 12 * HII_MCS_N_SIDE * HII_MCS_N_SIDE - 2 * ir * (ir + 1) + ip;
    }
}

void vec2pix_ring(const double *vec, double *vlen, int *pix)
{
  *vlen = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
  *pix  = ang2pix_ring_z_phi(vec[2] / (*vlen), atan2(vec[1], vec[0]));
}

#endif
