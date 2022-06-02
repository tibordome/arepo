/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/radtrandiff.c
 * \date        02/2021
 * \author      E. P. Bellinger, R. Pakmor and M. L. Winther
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
#include "../proto.h"
#include "solar.h"
#include "sys/stat.h"

#if defined(SOLAR_RADIATIVE_TRANSFER_DIFF) || defined(SOLAR_RADIATIVE_TRANSFER_EDD)

#include <gsl/gsl_linalg.h>

#include "HYPRE.h"
#include "HYPRE_krylov.h"
#include "HYPRE_parcsr_ls.h"
#include "_hypre_utilities.h"

/* Solver parameters */
#define RADTRANS_MAX_ITER 100
#define RADTRANS_ACCURACY 1.0e-8
#define RADTRANS_MAX_COLS 1000

/* DATA STRUCTURES */

static struct radtrans_face_data
{
  int cornerFirst;
  int cornerCount;

  int sng;

  int active;
  double dt;
  double area;
  // Variable data missing, find!

  double nx, ny, nz;
  double mx, my, mz;
  double px, py, pz;

  double failWeight;
} * radtrans_face_data;

static struct corner_list
{
  int index;
  double weight;
} * corner_list;

static struct corner_data
{
  int active;
  double matrix[NUMDIMS + 1][NUMDIMS + 1];

  // Parameters needed at the center of the interface
  double kaprho;
  double kaprhoGrad[3];
  // Will be to the fourth power
  double temperature;
  double temperatureGrad[3];

  int tetra;
  int fail;
} * corner_data;

/* DEFINE GLOBAL VALUES */
static int cornerCount;
double *kaprhoExch;
double *temperatureExch;
double ac = 4 * 5.670374 * pow(10, -5);  // erg cm^-2 s^-1 K^-4, from Wikipedia
int CParticle;

/* DEFINE INTERNAL FUNCTIONS */
static void prepare_stuff();
static void add_corner(int tetra, int iface);
static int get_timebin_of_point(int point);
static double get_surface_area_of_cell(int point);
static void compute_least_squares_matrix_at_corners();
static void compute_geometry_of_interface(int iface);
static void compute_quantities_at_corners();
static void free_stuff();
#ifdef SOLAR_RADIATIVE_TRANSFER_DIFF
static void radtrans_diff_approx();
#endif
/* DEFINE SHARED FUNCTIONS */
void exchange_variables();
double lookup_opacity(double temperature, double density, double X, double Z);
double compute_diffusion_term(int particle, int fullNormalGradients, double *radval);
void point_get_center(int p, double *Center);

/* MAIN ROUTINE */

void do_solar_rad_transfer()
{
  TIMER_START(CPU_SOLARRADTRANS);
  double Cradius = 1e12;
  int idx;
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      int i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double dt_cell = 0.5 * (P[i].TimeBinGrav ? (((integertime)1) << P[i].TimeBinGrav) : 0) * All.Timebase_interval;

      if(P[i].Type == 0)
        {
          double radius = sqrt((P[i].Pos[0] - 0.5 * All.BoxSize) * (P[i].Pos[0] - 0.5 * All.BoxSize) +
                               (P[i].Pos[1] - 0.5 * All.BoxSize) * (P[i].Pos[1] - 0.5 * All.BoxSize) +
                               (P[i].Pos[2] - 0.5 * All.BoxSize) * (P[i].Pos[2] - 0.5 * All.BoxSize));
          if(radius < Cradius)
            {
              Cradius   = radius;
              CParticle = i;
            }

          // inject energy into the core
          if(radius < All.CoreRadius)
            {  // cm
              SphP[i].Energy += All.VolumetricHeatingRate * SphP[i].Volume * dt_cell;
            }

          // Newton cooling at the surface
          // todo: change to mass fraction
          // add loop beforehand summing up masses across all cells
          if(radius > All.SurfaceRadiusInner && radius < All.SurfaceRadiusOuter)
            {
              SphP[i].Energy -= All.VolumetricCoolingRate * SphP[i].Volume * dt_cell;
            }
        }
    }
  update_primitive_variables();
  /* Start of Radiative transfer */
  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
    {
      prepare_stuff();

#ifdef SOLAR_RADIATIVE_TRANSFER_DIFF
      radtrans_diff_approx();
      mpi_printf("RADIATIVE_TRANSFER: Solving radiative transfer via diffusion approximation.\n");
#else
#ifdef SOLAR_RADIATIVE_TRANSFER_EDD
      radtrans_edd_approx();
#endif
#endif

      free_stuff();
    }

  update_primitive_variables();
  TIMER_STOP(CPU_SOLARRADTRANS);
}

/* SHARED ROUTINES */

void prepare_stuff()
{
  // Import variables from external cells
  kaprhoExch      = (double *)mymalloc("kaprhoExch", sizeof(double) * Mesh_nimport);
  temperatureExch = (double *)mymalloc("temperatureExch", sizeof(double) * Mesh_nimport);
  exchange_variables();

  // For all interfaces, make a list of all cells sharing a corner with the interface,
  // and save their center positions

  radtrans_face_data = (struct radtrans_face_data *)mymalloc("radtrans_face_data", Mesh.Nvf * sizeof(struct radtrans_face_data));
  corner_list        = (struct corner_list *)mymalloc("corner_list", Mesh.Nvf * 20 * sizeof(struct corner_list));
  corner_data        = (struct corner_data *)mymalloc("corner_data", Mesh.Ndt * sizeof(struct corner_data));

  // Reads "Compute matrices and magnetic fields at corners" but only sets corners inactive?
  int icorner;
  for(icorner = 0; icorner < Mesh.Ndt; icorner++)
    corner_data[icorner].active = 0;

  // Flag interfraces that are needed and compute their timestep
  int iface;
  for(iface = 0; iface < Mesh.Nvf; iface++)
    {
      struct radtrans_face_data *fd = &radtrans_face_data[iface];
      fd->active                    = 0;
      fd->sng                       = 0;

      int p1 = Mesh.VF[iface].p1;
      int p2 = Mesh.VF[iface].p2;

      if(p1 < 0 || p2 < 0)
        continue;

      int timebin1 = get_timebin_of_point(p1);
      int timebin2 = get_timebin_of_point(p2);

      if(!TimeBinSynchronized[timebin1] && !TimeBinSynchronized[timebin2])
        continue;

      if(Mesh.DP[p1].ID < Mesh.DP[p2].ID)
        {
          if(TimeBinSynchronized[timebin1])
            {
              // Lower ID is active, its task is responsible
              if(Mesh.DP[p1].task != ThisTask || Mesh.DP[p1].index >= NumGas)
                continue;
            }
          else
            {
              // Only higher ID is active, its task is responsible
              if(Mesh.DP[p2].task != ThisTask || Mesh.DP[p2].index >= NumGas)
                continue;
            }
        }
      else if(Mesh.DP[p1].ID > Mesh.DP[p2].ID)
        {
          if(TimeBinSynchronized[timebin2])
            {
              // Lower ID is active, its task is responsible
              if(Mesh.DP[p2].task != ThisTask || Mesh.DP[p2].index >= NumGas)
                continue;
            }
          else
            {
              // Only higher ID is active, its task is responsible
              if(Mesh.DP[p1].task != ThisTask || Mesh.DP[p1].index >= NumGas)
                continue;
            }
        }
      else
        {
          // Interface with at least one ghost point

          if(Mesh.DP[p1].task != ThisTask && Mesh.DP[p2].task != ThisTask)
            continue;

          if(Mesh.DP[p1].task != Mesh.DP[p2].task)
            terminate("This should not happen, I think...");

          if(Mesh.DP[p1].index >= NumGas && Mesh.DP[p2].index >= NumGas)
            continue;

          int p;
          if(Mesh.DP[p1].index < NumGas)
            p = p1;
          else
            p = p2;

          if(!TimeBinSynchronized[P[Mesh.DP[p].index].TimeBinHydro])
            continue;
        }

      double surfacearea = fmax(get_surface_area_of_cell(p1), get_surface_area_of_cell(p2));
      if(Mesh.VF[iface].area <= 1e-5 * surfacearea)
        continue;

      // If we make it here, the interface will be included, and this task is responsible for it
      fd->active = 1;

      // Now get its timestep
      int timeBin = imin(timebin1, timebin2);
      if(timeBin == 0)
        fd->dt = 0.;
      else
        fd->dt = (((integertime)1) << timeBin) * All.Timebase_interval / All.cf_hubble_a;
    }

  // Now begin filling informations for each interface
  cornerCount = 0;

  for(iface = 0; iface < Mesh.Nvf; iface++)
    {
      struct radtrans_face_data *fd = &radtrans_face_data[iface];
      int p1                        = Mesh.VF[iface].p1;
      int p2                        = Mesh.VF[iface].p2;

      if(p1 < 0 || p2 < 0 || fd->active == 0)
        {
          fd->area   = 0;
          fd->active = 0;
          continue;
        }

      fd->area        = Mesh.VF[iface].area;
      fd->cornerFirst = -1;
      fd->cornerCount = 0;

      int tt   = Mesh.VF[iface].dt_index;
      tetra *t = &Mesh.DT[tt];

      fd->cornerFirst = cornerCount;

      int nr;
      for(nr = 0; nr < 6; nr++)
        {
          int start_index = t->p[edge_start[nr]];
          int end_index   = t->p[edge_end[nr]];

          if((start_index == p1 && end_index == p2) || (start_index == p2 && end_index == p1))
            break;
        }

      tetra *prev, *next;
      int i, j, k, l, m, ii, jj, kk, ll, nn;

      i = edge_start[nr];
      j = edge_end[nr];
      k = edge_opposite[nr];
      l = edge_nexttetra[nr];

      prev = t;

      do
        {
          nn   = prev->t[l];
          next = &Mesh.DT[nn];

          add_corner(nn, iface);

          for(m = 0, ll = ii = jj = -1; m < 4; m++)
            {
              if(next->p[m] == prev->p[k])
                ll = m;
              if(next->p[m] == prev->p[i])
                ii = m;
              if(next->p[m] == prev->p[j])
                jj = m;
            }

          if(ll < 0 || ii < 0 || jj < 0)
            terminate("Inconsistency!");

          kk = 6 - (ll + ii + jj);

          prev = next;

          i = ii;
          l = ll;
          j = jj;
          k = kk;
        }
      while(next != t);
    }

  compute_least_squares_matrix_at_corners();

  compute_quantities_at_corners();

  for(iface = 0; iface < Mesh.Nvf; iface++)
    compute_geometry_of_interface(iface);

  MPI_Barrier(MPI_COMM_WORLD);
  mpi_printf("RADITAVE_TRANSFER: Preparations done!\n");
}

double lookup_opacity(double temperature, double density, double X, double Z)
{
  // From exercise in stellar models by JCD
  double kape  = 1.6236784e-33 * pow(density, 0.407895) * pow(temperature, 9.28289);
  double kapi  = 7.1548412e13 * pow(density, 0.138316) * pow(temperature, -1.97541);
  double kappa = pow((pow(kape, -1.) + pow(kapi, -1.)), -1.);

  if(kappa > 1.e5)
    mpi_printf("DB WARNING: kappa=%g\n", kappa);

  return kappa;
}

void exchange_variables()
{
  // Exchange Kaprho and Temp^4 to access later
  int listp;
  int i, j, p, task, off;
  int ngrp, recvTask, place;

  struct varDataExch
  {
    double Density;
    double EOSTemperature;
    double comp1;
    double comp2;
  } * tmpDataExch, *tmpDataRecv;

  tmpDataExch = (struct varDataExch *)mymalloc("tmpDataExch", Mesh_nexport * sizeof(struct varDataExch));

  // Prepare the data for export
  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(i = 0; i < NumGasInMesh; i++)
    {
      p = List_InMesh[i];

      listp = List_P[p].firstexport;
      while(listp >= 0)
        {
          if((task = ListExports[listp].origin) != ThisTask)
            {
              place = ListExports[listp].index;
              off   = Mesh_Send_offset[task] + Mesh_Send_count[task]++;

              tmpDataExch[off].Density        = SphP[place].Density;
              tmpDataExch[off].EOSTemperature = SphP[place].EOSTemperature;
              tmpDataExch[off].comp1          = SphP[place].Composition[0];
              tmpDataExch[off].comp2          = SphP[place].Composition[2];
            }
          listp = ListExports[listp].nextexport;
        }
    }

  // Exchange the data
  double kaprho;
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              tmpDataRecv = (struct varDataExch *)mymalloc("tmpDataRecv", Mesh_Recv_count[recvTask] * sizeof(struct varDataExch));

              // Get the values
              MPI_Sendrecv(&tmpDataExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask] * sizeof(struct varDataExch), MPI_BYTE,
                           recvTask, TAG_DENS_A, tmpDataRecv, Mesh_Recv_count[recvTask] * sizeof(struct varDataExch), MPI_BYTE,
                           recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              for(i = 0; i < Mesh_Recv_count[recvTask]; i++)
                {
                  kaprho = tmpDataRecv[i].Density * lookup_opacity(tmpDataRecv[i].EOSTemperature, tmpDataRecv[i].Density,
                                                                   tmpDataRecv[i].comp1, tmpDataRecv[i].comp2);
                  kaprhoExch[Mesh_Recv_offset[recvTask] + i]      = kaprho;
                  temperatureExch[Mesh_Recv_offset[recvTask] + i] = pow(tmpDataRecv[i].EOSTemperature, 4.0);
                }

              myfree(tmpDataRecv);
            }
        }
    }

  myfree(tmpDataExch);
}

int get_timebin_of_point(int point)
{
  int particle = Mesh.DP[point].index;

  if(Mesh.DP[point].task == ThisTask)
    {
      if(particle >= NumGas)
        particle -= NumGas;
      return P[particle].TimeBinHydro;
    }
  else
    return PrimExch[particle].TimeBinHydro;
}

void add_corner(int tt, int iface)
{
  if(cornerCount >= Mesh.Nvf * 20)
    terminate("urg");

  corner_list[cornerCount].index = tt;
  cornerCount++;

  corner_data[tt].active = 1;

  struct radtrans_face_data *fd = &radtrans_face_data[iface];
  fd->cornerCount++;
}

double get_surface_area_of_cell(int point)
{
  int particle = Mesh.DP[point].index;

  if(Mesh.DP[point].task == ThisTask)
    {
      if(particle >= NumGas)
        particle -= NumGas;
      return SphP[particle].ActiveArea;
    }
  else
    return PrimExch[particle].ActiveArea;
}

void compute_least_squares_matrix_at_corners()
{
  int failCount   = 0;
  int activeCount = 0;

  int totMoves = 0;

  int icorner;
  for(icorner = 0; icorner < Mesh.Ndt; icorner++)
    {
      struct corner_data *cd = &corner_data[icorner];

      cd->fail = 0;

      if(!cd->active)
        continue;

      activeCount++;

      cd->tetra = icorner;

      double weights[NUMDIMS + 1];
      int boundary = 0;

      int row, col, k;
      for(row = 0; row < NUMDIMS + 1; row++)
        {
          int point = Mesh.DT[cd->tetra].p[row];

          if(point < 0)
            {
              for(k = 0; k < NUMDIMS + 1; k++)
                for(col = 0; col < NUMDIMS + 1; col++)
                  cd->matrix[k][col] = 0;

              boundary = 1;
              break;
            }

          double Center[3];
          if(!TimeBinSynchronized[Mesh.DP[point].timebin])
            {
              Center[0] = Mesh.DP[point].x;
              Center[1] = Mesh.DP[point].y;
              Center[2] = Mesh.DP[point].z;
            }
          else
            {
              int particle = Mesh.DP[point].index;

              if(particle >= NumGas && Mesh.DP[point].task == ThisTask)
                particle -= NumGas;

              if(Mesh.DP[point].task == ThisTask)
                {
                  for(k = 0; k < 3; k++)
                    Center[k] = SphP[particle].Center[k];
                }
              else
                {
                  for(k = 0; k < 3; k++)
                    Center[k] = PrimExch[particle].Center[k];
                }
            }

          double weight = 1.;
          weights[row]  = weight;

          double dist[NUMDIMS];
          double xtmp, ytmp, ztmp;
          dist[0] = NEAREST_X(Center[0] - Mesh.DTC[icorner].cx);
          dist[1] = NEAREST_Y(Center[1] - Mesh.DTC[icorner].cy);
          dist[2] = NEAREST_Z(Center[2] - Mesh.DTC[icorner].cz);

          cd->matrix[row][0] = weight;
          for(col = 0; col < NUMDIMS; col++)
            cd->matrix[row][col + 1] = dist[col] * weight;
        }

      if(boundary)
        continue;

      double matrixT[(NUMDIMS + 1)][(NUMDIMS + 1)];
      double matrix[(NUMDIMS + 1) * (NUMDIMS + 1)];
      for(row = 0; row < NUMDIMS + 1; row++)
        for(col = 0; col < NUMDIMS + 1; col++)
          {
            matrixT[row][col] = cd->matrix[col][row];
            int idx           = row * (NUMDIMS + 1) + col;
            matrix[idx]       = 0;
            for(k = 0; k < NUMDIMS + 1; k++)
              matrix[idx] += cd->matrix[k][row] * cd->matrix[k][col];
          }

      int s;
      gsl_matrix_view m     = gsl_matrix_view_array(matrix, NUMDIMS + 1, NUMDIMS + 1);
      gsl_permutation *perm = gsl_permutation_alloc(NUMDIMS + 1);
      gsl_linalg_LU_decomp(&m.matrix, perm, &s);

      if(gsl_linalg_LU_det(&m.matrix, s) != 0.)
        {
          double matrix_inverse[(NUMDIMS + 1) * (NUMDIMS + 1)];
          gsl_matrix_view minv = gsl_matrix_view_array(matrix_inverse, NUMDIMS + 1, NUMDIMS + 1);
          gsl_linalg_LU_invert(&m.matrix, perm, &minv.matrix);
          gsl_permutation_free(perm);

          for(row = 0; row < NUMDIMS + 1; row++)
            for(col = 0; col < NUMDIMS + 1; col++)
              {
                cd->matrix[col][row] = 0;
                for(k = 0; k < NUMDIMS + 1; k++)
                  cd->matrix[col][row] += matrix_inverse[row * (NUMDIMS + 1) + k] * matrixT[k][col];

                cd->matrix[col][row] *= weights[row];
              }
        }
      else
        {
          for(row = 0; row < NUMDIMS + 1; row++)
            {
              cd->matrix[row][0] = 1. / (NUMDIMS + 1);
              for(col = 1; col < NUMDIMS + 1; col++)
                cd->matrix[row][col] = 0.;
            }
        }

      for(row = 0; row < NUMDIMS + 1; row++)
        {
          if(cd->matrix[row][0] < -0.01)
            {
              cd->fail = 1;
              break;
            }
        }

      if(cd->fail)
        failCount++;
    }

  int failCountSum, activeCountSum, totMovesSum;
  MPI_Reduce(&failCount, &failCountSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&activeCount, &activeCountSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&totMoves, &totMovesSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  mpi_printf("RADIATIVE_TRANSFER: %d out of %d active corners were bad, totMoves=%d.\n", failCountSum, activeCountSum, totMovesSum);
}

void compute_quantities_at_corners()
{
  int icorner, k;
  for(icorner = 0; icorner < Mesh.Ndt; icorner++)
    {
      struct corner_data *cd = &corner_data[icorner];

      if(!cd->active)
        continue;

      cd->kaprho      = 0;
      cd->temperature = 0;
      for(k = 0; k < NUMDIMS; k++)
        cd->kaprhoGrad[k] = 0;
      cd->temperatureGrad[k] = 0;

      tetra *tt = &Mesh.DT[cd->tetra];
      int p;
      for(p = 0; p < NUMDIMS + 1; p++)
        {
          int point = tt->p[p];

          if(point < 0)
            break;

          int particle = Mesh.DP[point].index;

          if(particle >= NumGas && Mesh.DP[point].task == ThisTask)
            particle -= NumGas;

          double kaprho, temperature;
          if(Mesh.DP[point].task == ThisTask)
            {
              // Opacity approximated by Kramer's approximation
              // https://en.wikipedia.org/wiki/Kramers%27_opacity_law
              kaprho      = SphP[particle].Density * lookup_opacity(SphP[particle].EOSTemperature, SphP[particle].Density,
                                                                    SphP[particle].Composition[0], SphP[particle].Composition[2]);
              temperature = pow(SphP[particle].EOSTemperature, 4.0);
            }
          else
            {
              kaprho      = kaprhoExch[particle];
              temperature = temperatureExch[particle];
            }

          cd->kaprho += kaprho * cd->matrix[p][0];
          cd->temperature += temperature * cd->matrix[p][0];
          for(k = 0; k < NUMDIMS; k++)
            {
              cd->kaprhoGrad[k] += kaprho * cd->matrix[p][k + 1];
              // Gradient of T^4
              cd->temperatureGrad[k] += temperature * cd->matrix[p][k + 1];
            }
        }
    }
}

static double get_vector_distance(double a[3], double b[3])
{
  double dx = a[0] - b[0];
  double dy = a[1] - b[1];
  double dz = a[2] - b[2];
  return sqrt(dx * dx + dy * dy + dz * dz);
}

static double get_area_triangle(double p1[3], double p2[3], double p3[3])
{
  double a = get_vector_distance(p1, p2);
  double b = get_vector_distance(p2, p3);
  double c = get_vector_distance(p3, p1);

  double s    = 0.5 * (a + b + c);
  double prod = s * (s - a) * (s - b) * (s - c);

  if(prod < 0.)
    return 0.;
  else
    return sqrt(prod);
}

void compute_geometry_of_interface(int iface)
{
  struct radtrans_face_data *fd = &radtrans_face_data[iface];

  if(!fd->active)
    return;

  fd->failWeight = 0;

  int p1 = Mesh.VF[iface].p1;
  int p2 = Mesh.VF[iface].p2;

  // ni: Unit vector in direction i, between center of cells
  double nx = Mesh.DP[p2].x - Mesh.DP[p1].x;
  double ny = Mesh.DP[p2].y - Mesh.DP[p1].y;
  double nz = Mesh.DP[p2].z - Mesh.DP[p1].z;
  double nn = sqrt(nx * nx + ny * ny + nz * nz);

  nx /= nn;
  ny /= nn;
  nz /= nn;

  fd->nx = nx;
  fd->ny = ny;
  fd->nz = nz;

  // We need an ortonormal basis
  if(fd->nx != 0 || fd->ny != 0)
    {
      fd->mx = -fd->ny;
      fd->my = fd->nx;
      fd->mz = 0;
    }
  else
    {
      fd->mx = 1;
      fd->my = 0;
      fd->mz = 0;
    }

  double mm = sqrt(fd->mx * fd->mx + fd->my * fd->my + fd->mz * fd->mz);
  fd->mx /= mm;
  fd->my /= mm;
  fd->mz /= mm;

  // Cross product
  fd->px = fd->ny * fd->mz - fd->nz * fd->my;
  fd->py = fd->nz * fd->mx - fd->nx * fd->mz;
  fd->pz = fd->nx * fd->my - fd->ny * fd->mx;

  // Compute weights to get value at center from values at corners
  // Mesh.DTC: Delaunay Triangle Center
  // cold: Location of midway point between corners
  double cold[3];
  cold[0] =
      0.5 * (Mesh.DTC[corner_list[fd->cornerFirst + fd->cornerCount - 1].index].cx + Mesh.DTC[corner_list[fd->cornerFirst].index].cx);
  cold[1] =
      0.5 * (Mesh.DTC[corner_list[fd->cornerFirst + fd->cornerCount - 1].index].cy + Mesh.DTC[corner_list[fd->cornerFirst].index].cy);
  cold[2] =
      0.5 * (Mesh.DTC[corner_list[fd->cornerFirst + fd->cornerCount - 1].index].cz + Mesh.DTC[corner_list[fd->cornerFirst].index].cz);

  // fc: location of face center
  double fc[3];
  fc[0] = Mesh.VF[iface].cx;
  fc[1] = Mesh.VF[iface].cy;
  fc[2] = Mesh.VF[iface].cz;

  int icorner;
  for(icorner = 0; icorner < fd->cornerCount; icorner++)
    {
      struct corner_list *cl = &corner_list[fd->cornerFirst + icorner];

      // p: Corner coordinate
      double p[3];
      p[0] = Mesh.DTC[cl->index].cx;
      p[1] = Mesh.DTC[cl->index].cy;
      p[2] = Mesh.DTC[cl->index].cz;

      // First the area between one halway point, face center and corner
      double area = get_area_triangle(cold, fc, p);

      int icurr = fd->cornerFirst + icorner;
      int inext = fd->cornerFirst + ((icorner + 1) % fd->cornerCount);

      cold[0] = 0.5 * (Mesh.DTC[corner_list[icurr].index].cx + Mesh.DTC[corner_list[inext].index].cx);
      cold[1] = 0.5 * (Mesh.DTC[corner_list[icurr].index].cy + Mesh.DTC[corner_list[inext].index].cy);
      cold[2] = 0.5 * (Mesh.DTC[corner_list[icurr].index].cz + Mesh.DTC[corner_list[inext].index].cz);

      // Then add the area between second halfway point, face center and corner
      area += get_area_triangle(cold, fc, p);
      cl->weight = area / Mesh.VF[iface].area;
    }

  double grady = 0.;
  double gradz = 0.;
  for(icorner = 0; icorner < fd->cornerCount; icorner++)
    {
      struct corner_list *cl = &corner_list[fd->cornerFirst + icorner];
      struct corner_data *cd = &corner_data[cl->index];

      // Align with kaprhoGrad
      int k;
      double grad_corner[3];
      for(k = 0; k < NUMDIMS; k++)
        grad_corner[k] = cd->kaprhoGrad[k];
      for(k = NUMDIMS; k < 3; k++)
        grad_corner[k] = 0.;

      grady += (grad_corner[0] * fd->mx + grad_corner[1] * fd->my + grad_corner[2] * fd->mz) * cl->weight;
      gradz += (grad_corner[0] * fd->px + grad_corner[1] * fd->py + grad_corner[2] * fd->pz) * cl->weight;
    }

  // Orient coordinate system in interface along gradient in kaprho
  // Afterwards, m points towards the largest gradient in the plane formed by m and p
  // n is still as first defined, p will be recalculated os orthogonal on n and m
  double mx = grady * fd->mx + gradz * fd->px;
  double my = grady * fd->my + gradz * fd->py;
  double mz = grady * fd->mz + gradz * fd->pz;

  mm = sqrt(mx * mx + my * my + mz * mz);
  if(mm > 0)
    {
      fd->mx = mx / mm;
      fd->my = my / mm;
      fd->mz = mz / mm;

      fd->px = fd->ny * fd->mz - fd->nz * fd->my;
      fd->py = fd->nz * fd->mx - fd->nx * fd->mz;
      fd->pz = fd->nx * fd->my - fd->ny * fd->mx;
    }

  for(icorner = 0; icorner < fd->cornerCount; icorner++)
    {
      struct corner_list *cl = &corner_list[fd->cornerFirst + icorner];
      struct corner_data *cd = &corner_data[cl->index];

      if(cd->fail)
        fd->failWeight += cl->weight;
    }
}

void point_get_center(int p, double *Center)
{
  if(!TimeBinSynchronized[Mesh.DP[p].timebin])
    {
      Center[0] = Mesh.DP[p].x;
      Center[1] = Mesh.DP[p].y;
      Center[2] = Mesh.DP[p].z;
      return;
    }

  int particle = Mesh.DP[p].index;
  int j;
  if(Mesh.DP[p].task == ThisTask)
    {
      if(particle >= NumGas)
        particle -= NumGas;

      for(j = 0; j < 3; j++)
        Center[j] = SphP[particle].Center[j];
    }
  else
    {
      for(j = 0; j < 3; j++)
        Center[j] = PrimExch[particle].Center[j];
    }
}

double compute_diffusion_term(int particle, int fullNormalGradients, double *radval)
{
  double val = 0.;
  // Simple gradient estimate of T^4
  // mpi_printf("DB: New particle!\n");
  if(!fullNormalGradients)
    {
      // Extract physical values for current particle
      double kaprho, temperature;
      kaprho      = SphP[particle].Density * lookup_opacity(SphP[particle].EOSTemperature, SphP[particle].Density,
                                                            SphP[particle].Composition[0], SphP[particle].Composition[2]);
      temperature = pow(SphP[particle].EOSTemperature, 4.0);

      // Index of first neighbour
      int q = SphP[particle].first_connection;

      if(particle == CParticle && 0)
        {
          // Format: Type, Number, x, y, z, rho, T, X, Y, Z, V, A
          mpi_printf(
              "*************************************************\n"
              "DB1: P, %d, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
              particle, P[particle].Pos[0], P[particle].Pos[1], P[particle].Pos[2], kaprho, SphP[particle].EOSTemperature,
              SphP[particle].Volume, 0.0, 0.0, 0.0, 0.0);
        }

      // Loop over each neighbour
      while(q >= 0)
        {
          int dp   = DC[q].dp_index;
          int vf   = DC[q].vf_index;
          int nbor = Mesh.DP[dp].index;

          struct radtrans_face_data *fd = &radtrans_face_data[vf];

          if(nbor < 0)
            {
              q = DC[q].next;
              continue;
            }

          if(nbor >= NumGas)
            nbor -= NumGas;

          // Compute geometric coefficients
          double dpCenter[3];
          point_get_center(dp, dpCenter);

          double xtmp, ytmp, ztmp;
          double dx    = NEAREST_X(dpCenter[0] - P[particle].Pos[0]);
          double dy    = NEAREST_Y(dpCenter[1] - P[particle].Pos[1]);
          double dz    = NEAREST_Z(dpCenter[2] - P[particle].Pos[2]);
          double dist2 = dx * dx + dy * dy + dz * dz;
          double gradN = pow(dist2, -0.5);
          double grad2 = (dx * fd->nx + dy * fd->ny + dz * fd->nz) / dist2;
          double area  = Mesh.VF[vf].area;

          // Extract physical values for neighbour
          double kaprhonbor, tempnbor, Volume;
          int task;
          double trho, ttemp;
          if(Mesh.DP[dp].task == ThisTask)
            {
              kaprhonbor = SphP[nbor].Density * lookup_opacity(SphP[nbor].EOSTemperature, SphP[nbor].Density,
                                                               SphP[nbor].Composition[0], SphP[nbor].Composition[2]);
              tempnbor   = pow(SphP[nbor].EOSTemperature, 4.0);
              Volume     = SphP[nbor].Volume;
              task       = 1;
              trho       = SphP[nbor].Density;
              ttemp      = SphP[nbor].EOSTemperature;
            }
          else
            {
              kaprhonbor = kaprhoExch[nbor];
              tempnbor   = temperatureExch[nbor];
              Volume     = PrimExch[nbor].Volume;
              task       = 2;
              trho       = -99.99;
              ttemp      = pow(tempnbor, 1. / 4.);
            }

          // Assume simple versions of value and gradient at interface
          double gammaface, tempgradface;
          gammaface    = 1. / kaprho;  // 0.5 *(1./kaprho + 1./kaprhonbor);
          tempgradface = (tempnbor - temperature) * gradN;
          if(particle == CParticle && 0)
            {
              mpi_printf("DB1: N, %d, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", nbor, dpCenter[0], dpCenter[1], dpCenter[2],
                         kaprhonbor, pow(tempnbor, 0.25), Volume, area, fd->nx, fd->ny, fd->nz);
            }

          double radcoord[3] = {(dpCenter[0] - 0.5 * All.BoxSize), (dpCenter[1] - 0.5 * All.BoxSize),
                                (dpCenter[2] - 0.5 * All.BoxSize)};
          double radius      = sqrt(radcoord[0] * radcoord[0] + radcoord[1] * radcoord[1] + radcoord[2] * radcoord[2]);
          double radnorm[3]  = {radcoord[0] / radius, radcoord[1] / radius, radcoord[2] / radius};
          double radgrad     = (dx * radnorm[0] + dy * radnorm[1] + dz * radnorm[2]) * gradN;

          // Add the value to the b vector entry for current particle
          double ival = ac * gammaface * area * tempgradface / (3. * SphP[particle].Volume);
          val -= ival;
          *radval -= radgrad * val;
          // mpi_printf("DB: gammaface=%g, area=%g, tempgradface=%g, temp1=%g, temp2=%g\n", gammaface, area, tempgradface, temperature,
          // tempnbor); mpi_printf("DB: T1=%g, T2=%g, rho=%g, kaprho=%g, val=%g, radius=%g\n", pow(temperature, 1./4.), ttemp,
          // SphP[particle].Density, kaprho, ival, radius/6.9e10); if(abs(val) > 1.e4 && radius > 4e10) terminate("val=%g, Still not
          // there Mark...\n", val);
          // Break for last connection, else continue
          if(q == SphP[particle].last_connection)
            break;

          q = DC[q].next;
        }
    }
  else
    {
      // Full gradient estimate
      // for interfaces of cell
      terminate("Full gradient not implemented for diffusion term!");
    }
  return val;
}

void free_stuff()
{
  // Free temporary arrays
  myfree(corner_data);
  myfree(corner_list);
  myfree(radtrans_face_data);
  myfree(temperatureExch);
  myfree(kaprhoExch);
}

#ifdef SOLAR_RADIATIVE_TRANSFER_DIFF
void radtrans_diff_approx()
{
  // TODO: We want to start with full gradient when implemented
  int fullNormalGradients = 0;
  int idx;
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      int particle = TimeBinsGravity.ActiveParticleList[idx];
      if(particle < 0)
        continue;

      double radval = 0.;

      double Qrad     = compute_diffusion_term(particle, fullNormalGradients, &radval);
      double dt       = 0.5 * (All.HighestActiveTimeBin ? (((integertime)1) << All.HighestActiveTimeBin) : 0) * All.Timebase_interval;
      double dEUpdate = Qrad * SphP[particle].Volume * dt;

      SphP[particle].Qrad = Qrad;
      SphP[particle].Energy += dEUpdate;
      SphP[particle].RadialQrad = radval;
    }
}
#endif
#endif
