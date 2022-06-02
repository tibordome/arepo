/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/windtunnel.c
 * \date        MM/YYYY
 * \author
 * \brief       Algorithms needed for windtunnel boundary conditions.
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"

#ifdef WINDTUNNEL

#ifndef WINDTUNNEL_COORD
#error "WINDTUNNEL requires WINDTUNNEL_COORD to be defined"
#endif

#ifdef WINDTUNNEL_EXTERNAL_SOURCE

static int NWindTable;

/*! \brief Structure stores wind properties (at a given time).
 */
static struct wind_table
{
  double time, rho, vel;
} * WindTable;

/*! \brief Interpolates density and velocity.
 *
 *  Interpolation in time between entries of wind table array.
 *
 *  \param[in] t Time.
 *  \param[out] rho Density (as a result of the interpolation).
 *  \param[out] vel Velocity (as a result of the interpolation).
 *
 *  \return void
 */
void interpolate_from_wind_table(double t, double *rho, double *vel)
{
  int binlow, binhigh, binmid;

  if(t < WindTable[0].time || t > WindTable[NWindTable - 1].time)
    {
      terminate("time outside of wind interpolation table");
    }

  binlow  = 0;
  binhigh = NWindTable - 1;

  while(binhigh - binlow > 1)
    {
      binmid = (binhigh + binlow) / 2;
      if(t < WindTable[binmid].time)
        binhigh = binmid;
      else
        binlow = binmid;
    }

  double dt = WindTable[binhigh].time - WindTable[binlow].time;

  if(dt == 0)
    terminate("dt=0");

  double u = (t - WindTable[binlow].time) / dt;

  *rho = (1 - u) * WindTable[binlow].rho + u * WindTable[binhigh].rho;
  *vel = (1 - u) * WindTable[binlow].vel + u * WindTable[binhigh].vel;
}

/*! \brief Read in windtunnel file.
 *
 *  Store the read-in input in WindTable.
 *
 *  \return void
 */
void read_windtunnel_file(void)
{
  FILE *fd;
  double k, p;

  if(!(fd = fopen(All.WindTunnelExternalSourceFile, "r")))
    terminate("can't read file '%s' with windtunnel data", All.WindTunnelExternalSourceFile);

  NWindTable = 0;
  do
    {
      double t, rho, v;
      if(fscanf(fd, " %lg %lg %lg ", &t, &rho, &v) == 3)
        NWindTable++;
      else
        break;
    }
  while(1);

  fclose(fd);

  mpi_printf("found %d rows in input wind table\n", NWindTable);

  WindTable = (struct wind_table *)mymalloc("WindTable", NWindTable * sizeof(struct wind_table));

  if(!(fd = fopen(All.WindTunnelExternalSourceFile, "r")))
    terminate("can't read file '%s' with windtunnel data", All.WindTunnelExternalSourceFile);

  NWindTable = 0;
  do
    {
      double t, rho, v;
      if(fscanf(fd, " %lg %lg %lg ", &t, &rho, &v) == 3)
        {
          WindTable[NWindTable].time = t;
          WindTable[NWindTable].rho  = rho;
          WindTable[NWindTable].vel  = v;
          NWindTable++;
        }
      else
        break;
    }
  while(1);

  fclose(fd);

  /* note: we'll assume that this file is sorted according to time */
}

#endif /* #if defined(WINDTUNNEL) && defined(WINDTUNNEL_EXTERNAL_SOURCE) */

#ifdef WINDTUNNEL_READ_IN_BFIELD

void WindtunnelReadIn_CalculateIndex(float x, float y, float z, int *nx, int *ny, int *nz, float *x_d, float *y_d, float *z_d)
{
  /*For a gas ceel, the indices and offsets (corresponding to the grid representing the read in field) are calculated */

  x /= All.BoxSize;
  y /= All.BoxSize;
  z /= All.BoxSize;

  y = y - All.Time * All.InjectionVelocity / All.BoxSize;  // Subtract the distance the cell has moved since start of simulation (we
                                                           // want the turbulent field to be comoving with fluid)

  x = fmodf(x, 1.0);
  y = fmodf(y, 1.0);
  z = fmodf(z, 1.0);

  if(x < All.WindtunnelReadIn_DX / 2.0)
    x = 1.0 + x;

  if(y < All.WindtunnelReadIn_DX / 2.0)
    y = 1.0 + y;

  if(z < All.WindtunnelReadIn_DX / 2.0)
    z = 1.0 + z;

  float index_x = (x - All.WindtunnelReadIn_DX / 2.0) / All.WindtunnelReadIn_DX;
  float index_y = (y - All.WindtunnelReadIn_DX / 2.0) / All.WindtunnelReadIn_DX;
  float index_z = (z - All.WindtunnelReadIn_DX / 2.0) / All.WindtunnelReadIn_DX;

  /*x_d is the fractional offset from two locations of on the grid. If the gridpoints are spaced as 1.0,1.1,1.2,1.3.. a particle has a
   * coordinate of x=1.239, then x_d = 0.039. It follows the convention from https://en.wikipedia.org/wiki/Trilinear_interpolation */
  *x_d = fmodf(index_x, 1.0);
  *y_d = fmodf(index_y, 1.0);
  *z_d = fmodf(index_z, 1.0);

  /* set the bin numbers */
  *nx = (int)index_x;
  *ny = (int)index_y;
  *nz = (int)index_z;
}

int WindtunnelReadIn_GetFlatIndex(int nx, int ny, int nz, int i, int j, int k)
{
  /* A helper function, used in TrilinearInterpolation, which calculates the flat index from nx,ny,nz and i,j,k. It takes periodic
   * wrapping into account. */
  nx = (nx + i) % NGRID_BFIELD;
  ny = (ny + j) % NGRID_BFIELD;
  nz = (nz + k) % NGRID_BFIELD;
  return nz + NGRID_BFIELD * nx + NGRID_BFIELD * NGRID_BFIELD * ny;
}

float WindtunnelReadIn_TrilinearInterpolation(float *U, int nx, int ny, int nz, float x_d, float y_d, float z_d)
{
  /* U is quantity to be interpolated. We follow conventions and formula from https://en.wikipedia.org/wiki/Trilinear_interpolation 2nd
   * order accuracy in dx is obtained. */

  float c000 = U[WindtunnelReadIn_GetFlatIndex(nx, ny, nz, 0, 0, 0)];
  float c100 = U[WindtunnelReadIn_GetFlatIndex(nx, ny, nz, 1, 0, 0)];
  float c001 = U[WindtunnelReadIn_GetFlatIndex(nx, ny, nz, 0, 0, 1)];
  float c101 = U[WindtunnelReadIn_GetFlatIndex(nx, ny, nz, 1, 0, 1)];
  float c010 = U[WindtunnelReadIn_GetFlatIndex(nx, ny, nz, 0, 1, 0)];
  float c110 = U[WindtunnelReadIn_GetFlatIndex(nx, ny, nz, 1, 1, 0)];
  float c011 = U[WindtunnelReadIn_GetFlatIndex(nx, ny, nz, 0, 1, 1)];
  float c111 = U[WindtunnelReadIn_GetFlatIndex(nx, ny, nz, 1, 1, 1)];

  float c00 = (1.0 - x_d) * c000 + x_d * c100;
  float c01 = (1.0 - x_d) * c001 + x_d * c101;
  float c10 = (1.0 - x_d) * c010 + x_d * c110;
  float c11 = (1.0 - x_d) * c011 + x_d * c111;

  float c0 = c00 * (1.0 - y_d) + c10 * y_d;
  float c1 = c01 * (1.0 - y_d) + c11 * y_d;

  float c = c0 * (1.0 - z_d) + c1 * z_d;

  return c;
}

void WindtunnelReadIn_InitialiseGlobals(void)
{
  int i;

  // Stuff related to read in of turbulent field:
  hid_t file_id, dataset_id, group_id, attr_id;
  herr_t status;
  int NPartType[6];
  float *pos_read, *B_read;

  // Read in stuff
  file_id = H5Fopen(All.WindtunnelReadIn_InputFileName, H5F_ACC_RDONLY, H5P_DEFAULT);

  if(file_id < 0)
    mpi_terminate("File missing e4f23n3.");

  group_id = H5Gopen2(file_id, "/Header", H5P_DEFAULT);

  attr_id = H5Aopen(group_id, "NumPart_ThisFile", H5P_DEFAULT);
  status  = H5Aread(attr_id, H5T_NATIVE_INT, NPartType);

  attr_id = H5Aopen(group_id, "DX", H5P_DEFAULT);

  status = H5Aread(attr_id, H5T_NATIVE_FLOAT, &All.WindtunnelReadIn_DX);

  pos_read = malloc(3 * NPartType[0] * sizeof(float));
  B_read   = malloc(3 * NPartType[0] * sizeof(float));

  dataset_id = H5Dopen2(file_id, "/PartType0/Coordinates", H5P_DEFAULT);  // change to GroupLenType
  status     = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos_read);
  status     = H5Dclose(dataset_id);

  dataset_id = H5Dopen2(file_id, "/PartType0/MagneticFieldRescaled", H5P_DEFAULT);  // change to GroupLenType
  status     = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, B_read);
  status     = H5Dclose(dataset_id);

  for(i = 0; i < NPartType[0]; i++)
    {
      All.WindtunnelReadIn_Bx[i] = B_read[3 * i + 0];
      All.WindtunnelReadIn_By[i] = B_read[3 * i + 1];
      All.WindtunnelReadIn_Bz[i] = B_read[3 * i + 2];
    }

  status = H5Gclose(group_id);
  status = H5Aclose(attr_id);
  status = H5Fclose(file_id);

  if(NPartType[0] != NGRID_BFIELD3)
    {
      mpi_printf("NGRID_BFIELD and NGRID_BFIELD3 are inconsistent with B-field-file!.");
      exit(21389);
    }
  if(All.WindtunnelReadIn_DX != pos_read[3 * NGRID_BFIELD] - pos_read[0])
    {
      mpi_printf("dx is inconsistent, should never happen.");
      exit(23545055);
    }

  free(pos_read);
  free(B_read);
}

#endif

#endif
