#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "proto.h"

/*
 * ----------------------------------------------------------------------
 * The following routine takes a value, x, and finds its position, i,
 * in a 1D table, along with the discplacement, dx, of x in between
 * the discrete table values. We assume here that the table is
 * evenly spaced.
 * ----------------------------------------------------------------------
 */

void get_index_1d_mydbl(double *table, int ntable, double x, int *i, double *dx)
{
  double denominator;

  denominator = (double)(table[ntable - 1] - table[0]) / (ntable - 1.0);

  if(x <= table[0])
    {
      *i  = 0;
      *dx = 0.0;
    }
  else if(x >= table[ntable - 1])
    {
      *i  = ntable - 2;
      *dx = 1.0;
    }
  else
    {
      *i  = (int)floor((x - table[0]) / denominator);
      *dx = (x - table[*i]) / denominator;
    }
}

/*
 * ----------------------------------------------------------------------
 * This takes a value, x, and finds its position, i, in a 1D table,
 * along with the discplacement, dx, of x in between the discrete table
 * values. However, here we do not assume that the table is evenly spaced.
 * ----------------------------------------------------------------------
 */

void get_index_1d_irregular(double *table, int ntable, double x, int *i, double *dx)
{
  if(x <= table[0])
    {
      *i  = 0;
      *dx = 0;
    }
  else if(x >= table[ntable - 1])
    {
      *i  = ntable - 2;
      *dx = 1;
    }
  else
    {
      *i = 0;
      while(table[*i] < x)
        *i += 1;

      *i -= 1;
      *dx = (x - table[*i]) / (table[(*i) + 1] - table[*i]);
    }
}

/*
 * ----------------------------------------------------------------------
 * Routine to perform a linear interpolation
 * ----------------------------------------------------------------------
 */

double interpol_1d_mydbl(double *table, int i, double dx)
{
  double output;

  if(dx <= 0.0)
    output = table[i];
  else if(dx >= 1.0)
    output = table[i + 1];
  else
    output = (1 - dx) * table[i] + dx * table[i + 1];

  return output;
}

/*
 * ----------------------------------------------------------------------
 * Routine to perform a linear interpolation. It is
 * to be used when the table is in single precision,
 * but you want to cast the result to double precision.
 * ----------------------------------------------------------------------
 */

double interpol_1d_fltdbl(float *table, int i, double dx)
{
  double output;

  if(dx <= 0.0)
    output = (double)table[i];
  else if(dx >= 1.0)
    output = (double)table[i + 1];
  else
    output = (1 - dx) * ((double)table[i]) + dx * ((double)table[i + 1]);

  return output;
}

/*
 * ----------------------------------------------------------------------
 * Routine to perform a bi-linear interpolation
 * ----------------------------------------------------------------------
 */

double interpol_2d_mydbl(double **table, int i, int j, double dx, double dy)
{
  double output, dx_m, dy_m;

  dx_m = 1.0 - dx;
  dy_m = 1.0 - dy;

  output = dx_m * dy_m * table[i][j] + dx_m * dy * table[i][j + 1] + dx * dy_m * table[i + 1][j] + dx * dy * table[i + 1][j + 1];

  return output;
}

/*
 * ----------------------------------------------------------------------
 * Routine to perform a bi-linear interpolation. It is
 * to be used when the table is in single precision,
 * but you want to cast the result to double precision.
 * ----------------------------------------------------------------------
 */

double interpol_2d_fltdbl(float **table, int i, int j, double dx, double dy)
{
  double output, dx_m, dy_m;

  dx_m = 1.0 - dx;
  dy_m = 1.0 - dy;

  output = dx_m * dy_m * ((double)table[i][j]) + dx_m * dy * ((double)table[i][j + 1]) + dx * dy_m * ((double)table[i + 1][j]) +
           dx * dy * ((double)table[i + 1][j + 1]);

  return output;
}

/*
 * ----------------------------------------------------------------------
 * Routine to perform a tri-linear interpolation
 * ----------------------------------------------------------------------
 */

double interpol_3d_mydbl(double ***table, int i, int j, int k, double dx, double dy, double dz)
{
  double output, dx_m, dy_m, dz_m;

  dx_m = 1.0 - dx;
  dy_m = 1.0 - dy;
  dz_m = 1.0 - dz;

  output = dx_m * dy_m * dz_m * table[i][j][k] + dx_m * dy_m * dz * table[i][j][k + 1] + dx_m * dy * dz_m * table[i][j + 1][k] +
           dx_m * dy * dz * table[i][j + 1][k + 1] + dx * dy_m * dz_m * table[i + 1][j][k] + dx * dy_m * dz * table[i + 1][j][k + 1] +
           dx * dy * dz_m * table[i + 1][j + 1][k] + dx * dy * dz * table[i + 1][j + 1][k + 1];

  return output;
}

/*
 * ----------------------------------------------------------------------
 * Routine to perform a tri-linear interpolation, but for a 4D table
 * (Note that the last index is kept fixed). This is used for the
 * equilibrium abundance tables.
 * ----------------------------------------------------------------------
 */

double interpol_3d_special(double ****table, int i, int j, int k, int l, double dx, double dy, double dz)
{
  double output, dx_m, dy_m, dz_m;

  dx_m = 1.0 - dx;
  dy_m = 1.0 - dy;
  dz_m = 1.0 - dz;

  output = dx_m * dy_m * dz_m * table[i][j][k][l] + dx_m * dy_m * dz * table[i][j][k + 1][l] +
           dx_m * dy * dz_m * table[i][j + 1][k][l] + dx_m * dy * dz * table[i][j + 1][k + 1][l] +
           dx * dy_m * dz_m * table[i + 1][j][k][l] + dx * dy_m * dz * table[i + 1][j][k + 1][l] +
           dx * dy * dz_m * table[i + 1][j + 1][k][l] + dx * dy * dz * table[i + 1][j + 1][k + 1][l];

  return output;
}

/*
 * ----------------------------------------------------------------------
 * Routine to perform a quadri-linear interpolation
 * ----------------------------------------------------------------------
 */

double interpol_4d_mydbl(double ****table, int i, int j, int k, int l, double dx, double dy, double dz, double dw)
{
  double output, dx_m, dy_m, dz_m, dw_m;

  dx_m = 1.0 - dx;
  dy_m = 1.0 - dy;
  dz_m = 1.0 - dz;
  dw_m = 1.0 - dw;

  output = dx_m * dy_m * dz_m * dw_m * table[i][j][k][l] + dx_m * dy_m * dz_m * dw * table[i][j][k][l + 1] +
           dx_m * dy_m * dz * dw_m * table[i][j][k + 1][l] + dx_m * dy_m * dz * dw * table[i][j][k + 1][l + 1] +
           dx_m * dy * dz_m * dw_m * table[i][j + 1][k][l] + dx_m * dy * dz_m * dw * table[i][j + 1][k][l + 1] +
           dx_m * dy * dz * dw_m * table[i][j + 1][k + 1][l] + dx_m * dy * dz * dw * table[i][j + 1][k + 1][l + 1] +
           dx * dy_m * dz_m * dw_m * table[i + 1][j][k][l] + dx * dy_m * dz_m * dw * table[i + 1][j][k][l + 1] +
           dx * dy_m * dz * dw_m * table[i + 1][j][k + 1][l] + dx * dy_m * dz * dw * table[i + 1][j][k + 1][l + 1] +
           dx * dy * dz_m * dw_m * table[i + 1][j + 1][k][l] + dx * dy * dz_m * dw * table[i + 1][j + 1][k][l + 1] +
           dx * dy * dz * dw_m * table[i + 1][j + 1][k + 1][l] + dx * dy * dz * dw * table[i + 1][j + 1][k + 1][l + 1];

  return output;
}

/*
 * ----------------------------------------------------------------------
 * Routine to perform a quinti-linear interpolation
 * ----------------------------------------------------------------------
 */

double interpol_5d_mydbl(double *****table, int i, int j, int k, int l, int m, double dx, double dy, double dz, double dw, double dv)
{
  double output, dx_m, dy_m, dz_m, dw_m, dv_m;

  dx_m = 1.0 - dx;
  dy_m = 1.0 - dy;
  dz_m = 1.0 - dz;
  dw_m = 1.0 - dw;
  dv_m = 1.0 - dv;

  output =
      dx_m * dy_m * dz_m * dw_m * dv_m * table[i][j][k][l][m] + dx_m * dy_m * dz_m * dw * dv_m * table[i][j][k][l + 1][m] +
      dx_m * dy_m * dz * dw_m * dv_m * table[i][j][k + 1][l][m] + dx_m * dy_m * dz * dw * dv_m * table[i][j][k + 1][l + 1][m] +
      dx_m * dy * dz_m * dw_m * dv_m * table[i][j + 1][k][l][m] + dx_m * dy * dz_m * dw * dv_m * table[i][j + 1][k][l + 1][m] +
      dx_m * dy * dz * dw_m * dv_m * table[i][j + 1][k + 1][l][m] + dx_m * dy * dz * dw * dv_m * table[i][j + 1][k + 1][l + 1][m] +
      dx * dy_m * dz_m * dw_m * dv_m * table[i + 1][j][k][l][m] + dx * dy_m * dz_m * dw * dv_m * table[i + 1][j][k][l + 1][m] +
      dx * dy_m * dz * dw_m * dv_m * table[i + 1][j][k + 1][l][m] + dx * dy_m * dz * dw * dv_m * table[i + 1][j][k + 1][l + 1][m] +
      dx * dy * dz_m * dw_m * dv_m * table[i + 1][j + 1][k][l][m] + dx * dy * dz_m * dw * dv_m * table[i + 1][j + 1][k][l + 1][m] +
      dx * dy * dz * dw_m * dv_m * table[i + 1][j + 1][k + 1][l][m] + dx * dy * dz * dw * dv_m * table[i + 1][j + 1][k + 1][l + 1][m] +
      dx_m * dy_m * dz_m * dw_m * dv * table[i][j][k][l][m + 1] + dx_m * dy_m * dz_m * dw * dv * table[i][j][k][l + 1][m + 1] +
      dx_m * dy_m * dz * dw_m * dv * table[i][j][k + 1][l][m + 1] + dx_m * dy_m * dz * dw * dv * table[i][j][k + 1][l + 1][m + 1] +
      dx_m * dy * dz_m * dw_m * dv * table[i][j + 1][k][l][m + 1] + dx_m * dy * dz_m * dw * dv * table[i][j + 1][k][l + 1][m + 1] +
      dx_m * dy * dz * dw_m * dv * table[i][j + 1][k + 1][l][m + 1] + dx_m * dy * dz * dw * dv * table[i][j + 1][k + 1][l + 1][m + 1] +
      dx * dy_m * dz_m * dw_m * dv * table[i + 1][j][k][l][m + 1] + dx * dy_m * dz_m * dw * dv * table[i + 1][j][k][l + 1][m + 1] +
      dx * dy_m * dz * dw_m * dv * table[i + 1][j][k + 1][l][m + 1] + dx * dy_m * dz * dw * dv * table[i + 1][j][k + 1][l + 1][m + 1] +
      dx * dy * dz_m * dw_m * dv * table[i + 1][j + 1][k][l][m + 1] + dx * dy * dz_m * dw * dv * table[i + 1][j + 1][k][l + 1][m + 1] +
      dx * dy * dz * dw_m * dv * table[i + 1][j + 1][k + 1][l][m + 1] +
      dx * dy * dz * dw * dv * table[i + 1][j + 1][k + 1][l + 1][m + 1];

  return output;
}
