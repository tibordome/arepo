/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/MRT/RT_set_VET.c
 * \date        MM/YYYY
 * \author      Rahul Kannan
 * \brief       Set the Eddington tensor for each cell
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
#include "../domain.h"
#include "../proto.h"
#include "../voronoi.h"

/* Set VET according to M1 closure ; Equation 13 of Rosdahl+2013*/

#if defined(MRT)

#define tiny 1e-8

void set_VET_single(int i, struct sph_particle_data *vetSphP)
{
  //  if(P[i].ID == 40000006)
  //    return ;

  double modF, f, chi, n[3];
  int num1;
  int j, k, l;

  for(num1 = 0; num1 < MRT_BINS; num1++)
    {
      //      assert(MRT_BINS==1) ;

      modF = sqrt(vetSphP[i].RT_F[num1][0] * vetSphP[i].RT_F[num1][0] + vetSphP[i].RT_F[num1][1] * vetSphP[i].RT_F[num1][1] +
                  vetSphP[i].RT_F[num1][2] * vetSphP[i].RT_F[num1][2]);

      if(vetSphP[i].DensPhot[num1] == 0.0)
        f = modF / (c_internal_units);
      else
        f = modF / (c_internal_units * vetSphP[i].DensPhot[num1]);

      if(f > 1.000000001)
        {
          double ratio = 0.999999999 / f;
          f            = 0.9999999999;
          vetSphP[i].RT_F[num1][0] *= ratio;
          vetSphP[i].RT_F[num1][1] *= ratio;
          vetSphP[i].RT_F[num1][2] *= ratio;
        }
        /*
          f=0.0 ;
        */

#ifdef MRT_FLUX_EXTRAPOLATION
      if(f > 1 && f < 1.000000001)
        f = 0.999999999;
      else if(f > 1.000000001)
        terminate("SET VET f>1\n");

      chi = (3.0 + 4.0 * f * f) / (5.0 + 2.0 * sqrt(4.0 - 3.0 * f * f));

      if(f > 1.0)
        terminate("Chi > 1.0 !!!!!, chi = %g, f = %g, DensPhot = %g, modF = %g\n, Fx = %g, Fy = %g, Fz = %g\n", chi, f,
                  vetSphP[i].DensPhot[num1], modF, vetSphP[i].RT_F[num1][0], vetSphP[i].RT_F[num1][1], vetSphP[i].RT_F[num1][2]);

      if(modF == 0.0)
        modF = 1.0;

      for(l = 0; l < 3; l++)
        n[l] = vetSphP[i].RT_F[num1][l] / modF;

      for(j = 0; j < 3; j++)
        {
          for(k = 0; k < 3; k++)
            {
              if(j == k)
                vetSphP[i].PT[num1][j][k] = ((1.0 - chi) / 2.0) + ((3.0 * chi - 1.0) / 2.0) * n[j] * n[k];
              else
                vetSphP[i].PT[num1][j][k] = ((3.0 * chi - 1.0) / 2.0) * n[j] * n[k];

              vetSphP[i].PT[num1][j][k] *= vetSphP[i].DensPhot[num1];
            }
        }
#endif
    }
  return;
}

#endif
