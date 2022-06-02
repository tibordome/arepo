/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \date        22/10/2018
 * \author      Ondrej Jaura
 * \brief       Mathematical functions
 * \details     
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

//======================================================================================================//
//                                  Algebraic functions                                                 //
//======================================================================================================//

/*!
  \brief Take the inner product of two 3D vectors
*/
double sx_math_compute_inner_product(double v1[SX_NDIM], double v2[SX_NDIM])
{
  // length of the first vector
  double length1 = pow(v1[0], 2.) + pow(v1[1], 2.) + pow(v1[2], 2.);
  // length of the second vector
  double length2 = pow(v2[0], 2.) + pow(v2[1], 2.) + pow(v2[2], 2.);

  double one_over_length1 = 1.0 / (sqrt(length1));
  double one_over_length2 = 1.0 / (sqrt(length2));

  double innerProduct = v1[0] * one_over_length1 * v2[0] * one_over_length2 + v1[1] * one_over_length1 * v2[1] * one_over_length2 +
                        v1[2] * one_over_length1 * v2[2] * one_over_length2;

  return innerProduct;
}

/*!
 *  /brief this function generates a random rotational matrix
 *
 *  implementation of the random rotational matrix (Graphics Gems III, Fast Random Rotation Matrices, Jim Arvo, 1991)
 *  M[3][3] is a 3x3 matrix array
 *  Angle "d" ( 0<d<pi ) is restricting the maximum random rotation from the original direction
 */
void sx_math_compute_rotational_matrix(double M[3][3], double d)
{
  d = (1.0 - cos(d)) * 0.5;  // converting angle to the internal parameter

  // NOTE: we have to calculate random values only on one processor and then distribute them to all processors !!!
  double randNums[3];
  if(ThisTask == 0)
    {
      randNums[0] = gsl_rng_uniform(sxRand);
      randNums[1] = gsl_rng_uniform(sxRand);
      randNums[2] = gsl_rng_uniform(sxRand);
    }
  MPI_Bcast(randNums, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  double theta = 2.0 * M_PI * randNums[0] * d;
  double phi   = 2.0 * M_PI * randNums[1];
  double z     = 2.0 * randNums[2] * d;
  double r     = sqrt(z);
  double Vx    = sin(phi) * r;
  double Vy    = cos(phi) * r;
  double Vz    = sqrt(2.0 - z);
  double st    = sin(theta);
  double ct    = cos(theta);
  double Sx    = Vx * ct - Vy * st;
  double Sy    = Vx * st + Vy * ct;
  M[0][0]      = Vx * Sx - ct;
  M[0][1]      = Vx * Sy - st;
  M[0][2]      = Vx * Vz;
  M[1][0]      = Vy * Sx + st;
  M[1][1]      = Vy * Sy - ct;
  M[1][2]      = Vy * Vz;
  M[2][0]      = Vz * Sx;
  M[2][1]      = Vz * Sy;
  M[2][2]      = 1.0 - z;
#ifdef SX_DISPLAY_STATS
  mpi_fprintf(FdSimplex,"Random matrix: theta %.03e phi %.03e z %.03e \n", theta, phi, z);
#endif
}

//======================================================================================================//
//                                  Numerical integration                                               //
//======================================================================================================//

/*!
 * \brief integrate discrete function in array between two values
 */
double sx_math_integrate_discrete_function(double a, double b, double x[], double function[], int nBins)
{
  // value of the integral to be computed
  double total = 0.0;

  unsigned int i = 0;
  _Bool stop     = 0;
  while(!stop)
    {
      if(x[i] >= a)
        {
          // use trapezium rule for unequally segmented data
          total += 0.5 * (x[i + 1] - x[i]) * (function[i] + function[i + 1]);
        }

      i++;

      if(x[i] >= b || i == nBins - 1)  // we need to set maximum number of bins to prevent the byte overflow!!!
        {
          stop = 1;
        }
    }

  return total;
}
