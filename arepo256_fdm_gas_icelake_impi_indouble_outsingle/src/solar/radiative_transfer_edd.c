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

#ifdef SOLAR_RADIATIVE_TRANSFER_EDD

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

struct column_data
{
  int index;
  double value;
  int next_column;
  int task;
};

struct matrix_data
{
  int local_row_count;
  int *local_row_index;  /* row index of a particle */
  int *local_part_index; /* particle index of a row */
  int *first_column;
  int *last_column;
  int *offsets;
  int *imported_particle_indizes;

  struct column_data *column_data;
  int cd_count;

  int *external_index; /* stores global index of external particles imported to PrimExch */
};

struct hypre_data
{
  int nRows;
  int *columnCounts;
  int *rows;
  int *cols;
  double *vals;
  double *bval;
  double *xval;
  int *tasks;
};

struct external_element
{
  int task;
  int index;
  int column_index;
  int originaltask;
  double value;
};

struct external_elements
{
  int N_external_elements;
  int Nmax_external_elements;
  struct external_element *elements;
};

struct external_cells
{
  int task;
  int index;
};

/* IMPORT CONSTANTS */
extern double *kaprhoExch;
extern double *temperatureExch;
extern double ac;
extern int CParticle;

/* DEFINE FUNCTIONS */
static void exchange_row_indizes(struct matrix_data *md);
static void allocate_rows_offsets_and_indizes(struct matrix_data *md);
static void compute_rows_offsets_and_indizes(struct matrix_data *md, int fullNormalGradients);
static void set_matrix_coefficients(struct matrix_data *md, struct hypre_data *hd, int fullNormalGradients);
static void apply_eddington_qrad(HYPRE_IJVector *x, struct matrix_data *md, double *dEUpdate);
static int eddington_approx_solve_matrix(struct matrix_data *md, struct hypre_data *hd, int useMultigridPreconditioner,
                                         MPI_Comm radtransComm);
static void add_matrix_element(int point_row, int point_column, double val, struct matrix_data *md, struct external_elements *ee);
static void compute_and_check_residuals(struct matrix_data *md, struct hypre_data *hd, MPI_Comm radtransComm, int *success,
                                        double *QradBackup);
static void free_hypre_data(struct hypre_data *hd);
static void free_rows_offsets_and_indizes(struct matrix_data *md);

/* EDDINGTON ROUTINES */
void radtrans_edd_approx()
{
  // Define data to be collected
  struct matrix_data md;
  struct hypre_data hd;

  // TODO: We want to start with full gradient when implemented
  int fullNormalGradients = 0;
  allocate_rows_offsets_and_indizes(&md);
  compute_rows_offsets_and_indizes(&md, fullNormalGradients);

  set_matrix_coefficients(&md, &hd, fullNormalGradients);

  int needSplit = 0;
  int task;
  for(task = 0; task < NTask; task++)
    if(md.offsets[task] > md.offsets[task + 1] - 1)
      needSplit++;

  if(needSplit == NTask)
    {
      free_hypre_data(&hd);
      free_rows_offsets_and_indizes(&md);
      return;
    }

  mpi_printf("RADIATIVE_TRANSFER: Using %d of %d cores for simple Eddington transfer.\n", NTask - needSplit, NTask);

  MPI_Comm radtransComm;
  if(needSplit)
    MPI_Comm_split(MPI_COMM_WORLD, md.local_row_count > 0 ? 1 : 2, ThisTask, &radtransComm);
  else
    radtransComm = MPI_COMM_WORLD;

  int failed = 0;
  failed     = eddington_approx_solve_matrix(&md, &hd, 0, radtransComm);

  if(failed)
    {
      mpi_printf("RADIATIVE_TRANSFER: Implicit step without multigrid preconditioner failed, trying again with it.\n");
      failed = eddington_approx_solve_matrix(&md, &hd, 1, radtransComm);
      if(failed)
        {
          if(fullNormalGradients)
            {
              mpi_printf("RADIATIVE_TRANSFER: Changing to simple normal gradients and repeating with multigrid preconditioner.\n");

              free_hypre_data(&hd);
              free_rows_offsets_and_indizes(&md);

              fullNormalGradients = 0;
              allocate_rows_offsets_and_indizes(&md);
              compute_rows_offsets_and_indizes(&md, fullNormalGradients);
              set_matrix_coefficients(&md, &hd, fullNormalGradients);

              failed = eddington_approx_solve_matrix(&md, &hd, 1, radtransComm);
            }

          if(failed)
            {
              eddington_approx_solve_matrix(&md, &hd, 2, radtransComm);
              mpi_printf("RADIATIVE_TRANSFER: Everyting failed, stopping. Make Mark do stuff...\n");
              terminate("HYPRE did not converge.");
            }
        }
    }

  free_hypre_data(&hd);
  free_rows_offsets_and_indizes(&md);

  if(needSplit)
    {
      MPI_Barrier(radtransComm);
      MPI_Comm_free(&radtransComm);
    }

  mpi_printf("RADIATIVE_TRANSFER: Successfully applied, now make debugging output!\n");

  MPI_Barrier(MPI_COMM_WORLD);
}

void allocate_rows_offsets_and_indizes(struct matrix_data *md)
{
  // NumGas: Number of gas particles on the LOCAL processor
  md->local_row_index = (int *)mymalloc("local_row_index", sizeof(int) * NumGas);
  int icell;
  for(icell = 0; icell < NumGas; icell++)
    md->local_row_index[icell] = -1;

  md->local_part_index = (int *)mymalloc("local_part_index", sizeof(int) * NumGas);
  md->first_column     = (int *)mymalloc("first_column", sizeof(int) * NumGas);
  md->last_column      = (int *)mymalloc("last_column", sizeof(int) * NumGas);

  // Why 20?
  md->column_data = (struct column_data *)mymalloc("column_data", sizeof(struct column_data) * NumGas * 20 * NUMDIMS);
  md->cd_count    = 0;

  // Mesh_nimport: ?
  md->imported_particle_indizes = (int *)mymalloc("imported_particle_indizes", sizeof(int) * Mesh_nimport);
  for(icell = 0; icell < Mesh_nimport; icell++)
    md->imported_particle_indizes[icell] = -1;

  // ThisTask: number of the local processor
  // NTask: number of processes
  md->offsets = (int *)mymalloc("offsets", sizeof(int) * (NTask + 1));

  md->local_row_count = 0;
}

static int external_cells_compare(const void *a, const void *b)
{
  if(((struct external_cells *)a)->task < (((struct external_cells *)b)->task))
    return -1;
  if(((struct external_cells *)a)->task > (((struct external_cells *)b)->task))
    return +1;

  return 0;
}

static int external_elements_compare(const void *a, const void *b)
{
  if(((struct external_element *)a)->task < (((struct external_element *)b)->task))
    return -1;

  if(((struct external_element *)a)->task > (((struct external_element *)b)->task))
    return +1;

  return 0;
}

void compute_rows_offsets_and_indizes(struct matrix_data *md, int fullNormalGradients)
{
  // Loop thorugh interfaces and mark cells that are needed
  int N_external_cells    = 0;
  int Nmax_external_cells = Mesh.Indi.AllocFacNflux;
  struct external_cells *external_cells;
  external_cells = (struct external_cells *)mymalloc_movable(&external_cells, "external_cells",
                                                             sizeof(struct external_cells) * Nmax_external_cells);

  int iface;
  for(iface = 0; iface < Mesh.Nvf; iface++)
    {
      struct radtrans_face_data *fd = &radtrans_face_data[iface];

      // Check whether the interface is active or it should be skipped
      if(!fd->active)
        continue;

      int k;
      for(k = 0; k < 2; k++)
        {
          int p = (k == 0 ? Mesh.VF[iface].p1 : Mesh.VF[iface].p2);

          if(Mesh.DP[p].task == ThisTask)
            {
              int particle = Mesh.DP[p].index;
              if(particle >= NumGas)
                particle -= NumGas;

              if(md->local_row_index[particle] == -1)
                {
                  // TODO: Why do we need to add this?
                  md->local_row_index[particle]             = md->local_row_count;
                  md->local_part_index[md->local_row_count] = particle;
                  md->first_column[md->local_row_count]     = -1;
                  md->last_column[md->local_row_count]      = -1;
                  md->local_row_count++;
                }
            }
          else
            {
              if(N_external_cells >= Nmax_external_cells)
                {
                  Nmax_external_cells *= ALLOC_INCREASE_FACTOR;
                  external_cells =
                      (struct external_cells *)myrealloc_movable(external_cells, Nmax_external_cells * sizeof(struct external_cells));
                }

              external_cells[N_external_cells].task  = Mesh.DP[p].task;
              external_cells[N_external_cells].index = Mesh.DP[p].originalindex;
              N_external_cells++;
            }
        }

      if(fullNormalGradients)
        {
          int icorner;
          for(icorner = 0; icorner < fd->cornerCount; icorner++)
            {
              struct corner_list *cl = &corner_list[fd->cornerFirst + icorner];
              struct corner_data *cd = &corner_data[cl->index];
              tetra *tt              = &Mesh.DT[cd->tetra];

              int icol;
              for(icol = 0; icol < NUMDIMS + 1; icol++)
                {
                  int p = tt->p[icol];

                  if(Mesh.DP[p].task == ThisTask)
                    {
                      int particle = Mesh.DP[p].index;
                      if(particle >= NumGas)
                        particle -= NumGas;

                      if(md->local_row_index[particle] == -1)
                        {
                          md->local_row_index[particle]             = md->local_row_count;
                          md->local_part_index[md->local_row_count] = particle;
                          md->first_column[md->local_row_count]     = -1;
                          md->last_column[md->local_row_count]      = -1;
                          md->local_row_count++;
                        }
                    }
                  else
                    {
                      if(N_external_cells >= Nmax_external_cells)
                        {
                          Nmax_external_cells *= ALLOC_INCREASE_FACTOR;
                          external_cells = (struct external_cells *)myrealloc_movable(
                              external_cells, Nmax_external_cells * sizeof(struct external_cells));
                        }

                      external_cells[N_external_cells].task  = Mesh.DP[p].task;
                      external_cells[N_external_cells].index = Mesh.DP[p].originalindex;
                      N_external_cells++;
                    }
                }
            }
        }
    }

  // Exchange list of external cells that need to be included
  int i, j, ngrp, recvTask, nimport;
  mysort(external_cells, N_external_cells, sizeof(struct external_cells), external_cells_compare);

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(i = 0; i < N_external_cells; i++)
    Send_count[external_cells[i].task]++;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  struct external_cells *external_cells_get =
      (struct external_cells *)mymalloc("external_cells_get", nimport * sizeof(struct external_cells));

  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              MPI_Sendrecv(&external_cells[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct external_cells), MPI_BYTE,
                           recvTask, TAG_DENS_A, &external_cells_get[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct external_cells), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  // Add cells that are needed from extern but are not already in our local list
  for(i = 0; i < nimport; i++)
    {
      int particle = external_cells_get[i].index;
      if(md->local_row_index[particle] == -1)
        {
          // Still need to add this?
          md->local_row_index[particle]             = md->local_row_count;
          md->local_part_index[md->local_row_count] = particle;
          md->first_column[md->local_row_count]     = -1;
          md->last_column[md->local_row_count]      = -1;
          md->local_row_count++;
        }
    }

  myfree(external_cells_get);
  myfree(external_cells);

  // Compute global offsets
  md->offsets[NTask] = 0;

  MPI_Allgather(&md->local_row_count, 1, MPI_INT, md->offsets, 1, MPI_INT, MPI_COMM_WORLD);

  int sum = 0;
  for(i = 0; i <= NTask; i++)
    {
      int tmp = sum;
      sum += md->offsets[i];
      md->offsets[i] = tmp;
    }

  exchange_row_indizes(md);
}

void exchange_row_indizes(struct matrix_data *md)
{
  int listp;
  int i, j, p, task, off;
  int ngrp, recvTask, place;
  int *tmpIndizesExch, *tmpIndizesRecv;

  tmpIndizesExch = (int *)mymalloc("tmpIndizesExch", Mesh_nexport * sizeof(int));

  // Prepare data for export
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

              if(md->local_row_index[place] == -1)
                tmpIndizesExch[off] = -1;
              else
                tmpIndizesExch[off] = md->local_row_index[place] + md->offsets[ThisTask];
            }
          listp = ListExports[listp].nextexport;
        }
    }

  // Exchange data
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              tmpIndizesRecv = (int *)mymalloc("tmpIndizesRecv", Mesh_Recv_count[recvTask] * sizeof(int));

              // Get the values
              MPI_Sendrecv(&tmpIndizesExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask] * sizeof(int), MPI_BYTE, recvTask,
                           TAG_DENS_A, tmpIndizesRecv, Mesh_Recv_count[recvTask] * sizeof(int), MPI_BYTE, recvTask, TAG_DENS_A,
                           MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              for(i = 0; i < Mesh_Recv_count[recvTask]; i++)
                {
                  if(Mesh_Recv_offset[recvTask] + i >= Mesh_nimport)
                    terminate("blurb");

                  md->imported_particle_indizes[Mesh_Recv_offset[recvTask] + i] = tmpIndizesRecv[i];
                }

              myfree(tmpIndizesRecv);
            }
        }
    }
  myfree(tmpIndizesExch);
}

void set_matrix_coefficients(struct matrix_data *md, struct hypre_data *hd, int fullNormalGradients)
{
  // Set matrix A and vector b for the solver, solving Ax+b=0
  if(md->local_row_count > 0)
    {
      hd->bval = (double *)mymalloc("bval", md->local_row_count * sizeof(double));

      for(int row = 0; row < md->local_row_count; row++)
        {
          double radval = 0;
          int particle  = md->local_part_index[row];
          double bval   = compute_diffusion_term(particle, fullNormalGradients, &radval);

          // Store the computed bval of the row/particle
          hd->bval[row]             = bval;
          SphP[particle].RadialQrad = radval;
        }
    }

  struct external_elements ee;

  ee.N_external_elements    = 0;
  ee.Nmax_external_elements = Mesh.Indi.AllocFacNflux;
  ee.elements               = (struct external_element *)mymalloc_movable(&ee.elements, "external_elements",
                                                                          sizeof(struct external_element) * ee.Nmax_external_elements);

  int row;
  for(row = 0; row < md->local_row_count; row++)
    {
      int particle = md->local_part_index[row];

      struct column_data *cd     = &md->column_data[md->cd_count];
      md->first_column[particle] = md->cd_count;
      md->last_column[particle]  = md->cd_count;
      md->cd_count++;

      if(md->cd_count >= NumGas * 20 * NUMDIMS)
        terminate("baaad");

      cd->next_column = -1;
      cd->value       = 1.0;
      cd->index       = row + md->offsets[ThisTask];
      cd->task        = ThisTask;
    }

  /* Go through all interfaces and fill vector and matrix entries */
  int iface;
  for(iface = 0; iface < Mesh.Nvf; iface++)
    {
      struct radtrans_face_data *fd = &radtrans_face_data[iface];
      int particle1                 = Mesh.DP[Mesh.VF[iface].p1].index;
      int particle2                 = Mesh.DP[Mesh.VF[iface].p2].index;
      if((particle1 == CParticle || particle2 == CParticle) && 0)
        mpi_printf("DB: part1=%d, part2=%d\n", particle1, particle2);

      if(!fd->active)
        continue;
      if(particle1 >= NumGas)
        particle1 -= NumGas;
      if(particle2 >= NumGas)
        particle2 -= NumGas;

      int boundary = 0;

      // TODO: What are the logic behind these statements?
      if(!fullNormalGradients || fd->failWeight > 0 || boundary || fd->sng)
        {
          double weight = 1.;

          if(fullNormalGradients && (!boundary) && (!fd->sng))
            weight = fd->failWeight;

          double Center1[3], Center2[3];
          point_get_center(Mesh.VF[iface].p1, Center1);
          point_get_center(Mesh.VF[iface].p2, Center2);

          double xtmp, ytmp, ztmp;
          double dx = NEAREST_X(Center2[0] - Center1[0]);
          double dy = NEAREST_Y(Center2[1] - Center1[1]);
          double dz = NEAREST_Z(Center2[2] - Center1[2]);

          double dist2 = dx * dx + dy * dy + dz * dz;

          if(dist2 > 0.)
            {
              // Compute gradient and dot product with interface
              double gradN = (dx * fd->nx + dy * fd->ny + dz * fd->nz) / dist2;
              double mval1 = 0;
              double mval2 = 0;

              // Extract physical values for both particles
              double kaprho1, kaprho2;
              if(Mesh.DP[Mesh.VF[iface].p1].task == ThisTask)
                kaprho1 = SphP[particle1].Density * lookup_opacity(SphP[particle1].EOSTemperature, SphP[particle1].Density,
                                                                   SphP[particle1].Composition[0], SphP[particle1].Composition[2]);
              else
                kaprho1 = kaprhoExch[particle1];
              if(Mesh.DP[Mesh.VF[iface].p2].task == ThisTask)
                kaprho2 = SphP[particle2].Density * lookup_opacity(SphP[particle2].EOSTemperature, SphP[particle2].Density,
                                                                   SphP[particle2].Composition[0], SphP[particle2].Composition[2]);
              else
                kaprho2 = kaprhoExch[particle2];

              // Assume simple versions of value and gradient at interface
              double gammaface = 0.5 * (1. / kaprho1 + 1. / kaprho2);

              mval1 += gammaface * fd->area * gradN * kaprho1;
              mval2 += gammaface * fd->area * gradN * kaprho2;
              if((particle1 == CParticle || particle2 == CParticle) && 0)
                mpi_printf("DB: part1=%d, part2=%d, mval1=%g, mval2=%g, gradN=%g, kaprho1=%g, kaprho2=%g, area=%g\n", particle1,
                           particle2, mval1, mval2, gradN, kaprho1, kaprho2, fd->area);
              // TODO: Why Mesh.VF[iface].p1 and not particle1?
              add_matrix_element(Mesh.VF[iface].p1, Mesh.VF[iface].p1, +mval1, md, &ee);
              add_matrix_element(Mesh.VF[iface].p2, Mesh.VF[iface].p2, +mval2, md, &ee);
              add_matrix_element(Mesh.VF[iface].p1, Mesh.VF[iface].p2, -mval2, md, &ee);
              add_matrix_element(Mesh.VF[iface].p2, Mesh.VF[iface].p1, -mval1, md, &ee);
            }
        }

      if(fullNormalGradients && (!boundary) && (!fd->sng))
        {
          // TODO: Full normal gradient of face
          terminate("Full gradient of A not implemented");
          int icorner;
          for(icorner = 0; icorner < fd->cornerCount; icorner++)
            {
              struct corner_list *cl = &corner_list[fd->cornerFirst + icorner];
              struct corner_data *cd = &corner_data[cl->index];
              tetra *tt              = &Mesh.DT[cd->tetra];

              if(cd->fail || cl->weight == 0.)
                continue;

              int icol;
              for(icol = 0; icol < NUMDIMS + 1; icol++)
                {
                  int col = tt->p[icol];

                  // Geometric gradiant of corner contribution
                  double gx =
                      (cd->matrix[icol][1] * fd->nx + cd->matrix[icol][2] * fd->ny + cd->matrix[icol][3] * fd->nz) * cl->weight;

                  double bvali  = 0;
                  double mvalij = 0;
                  double mvalG  = 0;

                  // TODO: 1/3Vi factor in subroutines
                  mvalij -= cd->kaprho * fd->area *
                            (cd->kaprhoGrad[0] * fd->nx + cd->kaprhoGrad[1] * fd->ny + cd->kaprhoGrad[2] * fd->nz) * cl->weight;
                  mvalG += cd->kaprho * cd->kaprho * fd->area * gx;

                  // TODO: I'm not sure this behaviour is correct
                  // add_matrix_element(Mesh.VF[iface].p1, col, +mvalij + mvalG, md, &ee);
                  // add_matrix_element(Mesh.VF[iface].p2, col, -mvalij - mvalG, md, &ee);
                }
            }
        }
    }

  /* Exchange matrix entries for external rows */
  int i, j, nimport, ngrp, recvTask;
  mysort(ee.elements, ee.N_external_elements, sizeof(struct external_element), external_elements_compare);

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(i = 0; i < ee.N_external_elements; i++)
    Send_count[ee.elements[i].task]++;

  if(Send_count[ThisTask] > 0)
    terminate("Send_count[ThisTask]");

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  struct external_element *external_elements_get =
      (struct external_element *)mymalloc("external_elements_get", nimport * sizeof(struct external_element));

  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              MPI_Sendrecv(&ee.elements[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct external_element), MPI_BYTE,
                           recvTask, TAG_DENS_A, &external_elements_get[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct external_element), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  for(i = 0; i < nimport; i++)
    {
      int row          = external_elements_get[i].index;
      int column       = external_elements_get[i].column_index;
      int particle     = md->local_part_index[row - md->offsets[ThisTask]];
      int originaltask = external_elements_get[i].originaltask;

      struct column_data *cd = NULL;

      int entry = md->first_column[particle];
      while(entry >= 0)
        {
          cd = &md->column_data[entry];

          if(cd->index == column)
            break;

          entry = cd->next_column;
        }

      if(!cd || cd->index != column)
        {
          if(md->cd_count >= NumGas * 20 * NUMDIMS)
            terminate("List too long");

          cd = &md->column_data[md->cd_count];
          if(md->first_column[particle] == -1)
            md->first_column[particle] = md->cd_count;
          else
            md->column_data[md->last_column[particle]].next_column = md->cd_count;

          md->last_column[particle] = md->cd_count;
          md->cd_count++;

          cd->index       = column;
          cd->value       = external_elements_get[i].value;
          cd->next_column = -1;
          cd->task        = originaltask;
        }
      else
        {
          cd->value += external_elements_get[i].value;

          if(cd->task != originaltask)
            {
              printf("ThisTask=%d, cd->task=%d, originaltask=%d, entry=%d\n", ThisTask, cd->task, originaltask, entry);
              terminate("Failed!");
            }
        }
    }

  myfree(external_elements_get);
  myfree(ee.elements);

  /* Check if there's a contribution */
  if(md->local_row_count > 0)
    {
      /* Build matrix and send data to HYPRE */
      hd->nRows = md->local_row_count;

      hd->columnCounts = (int *)mymalloc("columnCounts", md->local_row_count * sizeof(int));
      hd->rows         = (int *)mymalloc("rows", md->local_row_count * sizeof(int));
      hd->cols         = (int *)mymalloc("cols", md->cd_count * sizeof(int));

      hd->vals  = (double *)mymalloc("vals", md->cd_count * sizeof(double));
      hd->xval  = (double *)mymalloc("xval", md->local_row_count * sizeof(double));
      hd->tasks = (int *)mymalloc("point", md->cd_count * sizeof(int));

      int colCount = 0;
      for(row = 0; row < md->local_row_count; row++)
        {
          int particle = md->local_part_index[row];

          // Initialization of xval, guessing previous
          hd->xval[row] = SphP[particle].Qrad;
          hd->rows[row] = row + md->offsets[ThisTask];

          hd->columnCounts[row] = 0;

          double fac = 1.;  // Changing this changes the weighting for the global residual!!!

          int column_index = md->first_column[particle];
          while(column_index >= 0)
            {
              struct column_data *cd = &md->column_data[column_index];

              hd->columnCounts[row]++;
              hd->cols[colCount]  = cd->index;
              hd->vals[colCount]  = cd->value * fac;
              hd->tasks[colCount] = cd->task;
              colCount++;

              column_index = cd->next_column;
            }

          hd->bval[row] *= fac;
        }
    }
  else
    {
      hd->nRows        = 0;
      hd->columnCounts = 0;
      hd->rows         = 0;
      hd->cols         = 0;
      hd->vals         = 0;
      hd->bval         = 0;
      hd->xval         = 0;
    }
}

void add_matrix_element(int point_row, int point_column, double val, struct matrix_data *md, struct external_elements *ee)
{
  // Add matrix element contributions to the i,j index (point_row, point_column)
  int column;

  if(Mesh.DP[point_column].task == ThisTask)
    {
      int particle = Mesh.DP[point_column].index;
      if(particle >= NumGas)
        particle -= NumGas;

      column = md->local_row_index[particle] + md->offsets[ThisTask];

      if(column < 0 || column > md->offsets[NTask] - 1)
        {
          printf("Column fail: %d, last offset: %d, particle: %d, local_row_index: %d, offset: %d\n", column, md->offsets[NTask],
                 particle, md->local_row_index[particle], md->offsets[ThisTask]);
          terminate("Failed matrix construction");
        }
    }
  else
    {
      int particle = Mesh.DP[point_column].index;

      if(particle < 0)
        terminate("Particle < 0");
      if(particle >= Mesh_nimport)
        terminate("Particle >= Mesh_nimport");

      column = md->imported_particle_indizes[particle];

      if(column < 0 || column > md->offsets[NTask] - 1)
        {
          printf("Column fail: %d, last offset: %d\n", column, md->offsets[NTask]);
          terminate("Failed matrix construction");
        }
    }

  if(Mesh.DP[point_row].task == ThisTask)
    {
      /* Add to local list */
      int particle = Mesh.DP[point_row].index;
      if(particle >= NumGas)
        particle -= NumGas;

      int entry = md->first_column[particle];
      while(entry >= 0)
        {
          struct column_data *cd = &md->column_data[entry];

          if(cd->index == column)
            break;

          entry = cd->next_column;
        }

      if(entry == -1)
        {
          /* Need new entry */
          if(md->cd_count >= NumGas * 20 * NUMDIMS)
            terminate("Column data list too long");

          struct column_data *cd = &md->column_data[md->cd_count];
          if(md->first_column[particle] == -1)
            md->first_column[particle] = md->cd_count;
          else
            md->column_data[md->last_column[particle]].next_column = md->cd_count;

          md->last_column[particle] = md->cd_count;
          md->cd_count++;

          // TODO: Might reconsider 1/3Vi
          cd->index       = column;
          cd->value       = val / (3. * SphP[particle].Volume);
          cd->next_column = -1;
          cd->task        = Mesh.DP[point_column].task;
        }
      else
        {
          struct column_data *cd = &md->column_data[entry];
          cd->value += val / (3. * SphP[particle].Volume);

          if(cd->task != Mesh.DP[point_column].task)
            terminate("Task mismatch");
        }
    }
  else
    {
      /* Not owned by this task, prepare to send to other task */
      if(ee->N_external_elements >= ee->Nmax_external_elements)
        {
          ee->Nmax_external_elements *= ALLOC_INCREASE_FACTOR;
          ee->elements =
              (struct external_element *)myrealloc_movable(ee->elements, sizeof(struct external_element) * ee->Nmax_external_elements);
        }

      struct external_element *el = &ee->elements[ee->N_external_elements];
      ee->N_external_elements++;

      el->task         = Mesh.DP[point_row].task;
      el->index        = md->imported_particle_indizes[Mesh.DP[point_row].index];
      el->column_index = column;
      el->value        = val / (3. * PrimExch[Mesh.DP[point_row].index].Volume);
      el->originaltask = Mesh.DP[point_column].task;
    }
}

int eddington_approx_solve_matrix(struct matrix_data *md, struct hypre_data *hd, int useMultigridPreconditioner, MPI_Comm radtransComm)
{
  int success = 0;

  // Backup Qrad for each cell, should it fail
  double *QradBackup = (double *)mymalloc("Qrad", hd->nRows * sizeof(double));
  double *dEUpdate   = (double *)mymalloc("dEUpdate", hd->nRows * sizeof(double));
  int row;
  for(row = 0; row < hd->nRows; row++)
    {
      int particle    = md->local_part_index[row];
      QradBackup[row] = SphP[particle].Qrad;
    }

  if(hd->rows > 0)
    {
      // TODO: Timer for solver?

      // Initialize equation elements Ax + b = 0
      HYPRE_IJMatrix A;
      HYPRE_ParCSRMatrix parcsr_A;
      HYPRE_IJVector b;
      HYPRE_ParVector par_b;
      HYPRE_IJVector x;
      HYPRE_ParVector par_x;

      // Initialize the solver, and its preconditions
      HYPRE_Solver solver, precond;

      // Range of row indices in the matrix for the current process
      int ilower = md->offsets[ThisTask];
      int iupper = md->offsets[ThisTask + 1] - 1;

      // Create and initialize elements
      HYPRE_IJMatrixCreate(radtransComm, ilower, iupper, ilower, iupper, &A);
      HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
      HYPRE_IJMatrixInitialize(A);

      HYPRE_IJVectorCreate(radtransComm, ilower, iupper, &b);
      HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize(b);

      HYPRE_IJVectorCreate(radtransComm, ilower, iupper, &x);
      HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize(x);

      // Fill the elements with their values
      HYPRE_IJMatrixSetValues(A, hd->nRows, hd->columnCounts, hd->rows, hd->cols, hd->vals);
      HYPRE_IJVectorSetValues(b, hd->nRows, hd->rows, hd->bval);
      HYPRE_IJVectorSetValues(x, hd->nRows, hd->rows, hd->xval);

      // Finalize assembly of elements and retrieve pointer/"handle"
      HYPRE_IJMatrixAssemble(A);
      HYPRE_IJMatrixGetObject(A, (void **)&parcsr_A);

      HYPRE_IJVectorAssemble(b);
      HYPRE_IJVectorGetObject(b, (void **)&par_b);

      HYPRE_IJVectorAssemble(x);
      HYPRE_IJVectorGetObject(x, (void **)&par_x);

      if(useMultigridPreconditioner == 2)
        {
          // For failed run, print matrix and stop
          mkdir("matrix", MKDIR_MODE);

          HYPRE_IJMatrixPrint(A, "matrix/IJ.out.A");
          HYPRE_IJVectorPrint(b, "matrix/IJ.out.b");
          HYPRE_IJVectorPrint(x, "matrix/IJ.out.x");

          MPI_Barrier(radtransComm);
          terminate("HYPRE did not converge!");
        }

      int nIter;
      double final_res_norm;

      // Initialize the solver and its coefficients
      HYPRE_ParCSRFlexGMRESCreate(radtransComm, &solver);
      HYPRE_ParCSRFlexGMRESSetTol(solver, RADTRANS_ACCURACY);
      HYPRE_ParCSRFlexGMRESSetPrintLevel(solver, 3);
      HYPRE_ParCSRFlexGMRESSetLogging(solver, 1);

      // If the default preconditioner failed, try a more advanced
      if(useMultigridPreconditioner)
        {
          HYPRE_ParCSRFlexGMRESSetKDim(solver, 500);
          HYPRE_ParCSRFlexGMRESSetMaxIter(solver, 500);

          HYPRE_BoomerAMGCreate(&precond);
          HYPRE_BoomerAMGSetCoarsenType(precond, 6);
          HYPRE_BoomerAMGSetInterpType(precond, 6);
          HYPRE_BoomerAMGSetStrongThreshold(precond, 0.25);
          HYPRE_BoomerAMGSetPMaxElmts(precond, 4);
          HYPRE_BoomerAMGSetRelaxType(precond, 18);
          HYPRE_BoomerAMGSetNumSweeps(precond, 1);
          HYPRE_BoomerAMGSetTol(precond, 0.0);
          HYPRE_BoomerAMGSetPrintLevel(precond, 1);
          HYPRE_BoomerAMGSetMaxIter(precond, 2);

          HYPRE_ParCSRFlexGMRESSetPrecond(solver, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, precond);
        }
      else
        {
          HYPRE_ParCSRFlexGMRESSetKDim(solver, 200);
          HYPRE_ParCSRFlexGMRESSetMaxIter(solver, 200);
        }

      // TODO: Insert timers in the following?
      // Setup and solve the solver
      HYPRE_ParCSRFlexGMRESSetup(solver, parcsr_A, par_b, par_x);
      HYPRE_ParCSRFlexGMRESSolve(solver, parcsr_A, par_b, par_x);
      HYPRE_ParCSRFlexGMRESGetNumIterations(solver, &nIter);
      HYPRE_ParCSRFlexGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);

      // Was it successful, in terms of hypre residual? Apply result
      success = final_res_norm < RADTRANS_ACCURACY;
      apply_eddington_qrad(&x, md, dEUpdate);

      // Release assigned elements
      HYPRE_ParCSRFlexGMRESDestroy(solver);
      HYPRE_IJMatrixDestroy(A);
      HYPRE_IJVectorDestroy(b);
      HYPRE_IJVectorDestroy(x);
      if(useMultigridPreconditioner)
        HYPRE_BoomerAMGDestroy(precond);

      int radtransRank;
      MPI_Comm_rank(radtransComm, &radtransRank);
      if(radtransRank == 0)
        printf("RADIATIVE_TRANSFER: HYPRE iter %d, res_norm %g\n", nIter, final_res_norm);
    }

  // Get parameter values from other tasks
  exchange_variables();

  compute_and_check_residuals(md, hd, radtransComm, &success, QradBackup);

  if(0)
    {
      mpi_printf("DB2: P, %d, %g\n", CParticle, SphP[CParticle].Qrad);

      int q = SphP[CParticle].first_connection;
      while(q >= 0)
        {
          int dp   = DC[q].dp_index;
          int nbor = Mesh.DP[dp].index;

          if(nbor < 0)
            {
              q = DC[q].next;
              continue;
            }

          if(nbor >= NumGas)
            nbor -= NumGas;

          mpi_printf("DB2: N, %d, %g\n", nbor, SphP[nbor].Qrad);
          if(q == SphP[CParticle].last_connection)
            break;

          q = DC[q].next;
        }
      mpi_printf("DBEND\n");
      terminate("RADIATIVE_TRANSFER: First time step completed!");
    }
  if(success)
    {
      myfree(dEUpdate);
      myfree(QradBackup);
      return 0;
    }
  else
    {
      // Restore values from backup
      for(row = 0; row < hd->nRows; row++)
        {
          int particle = md->local_part_index[row];

          SphP[particle].Qrad = QradBackup[row];
          SphP[particle].Energy -= dEUpdate[row];
        }
      myfree(dEUpdate);
      myfree(QradBackup);
      return 1;
    }
}

// TODO: This whole routine is just black magic
void compute_and_check_residuals(struct matrix_data *md, struct hypre_data *hd, MPI_Comm radtransComm, int *success,
                                 double *QradBackup)
{
  double residual_max = -1;
  double residual_b   = -1;
  int residual_index  = -1;

  double *residuals = (double *)mymalloc("residuals", hd->nRows * sizeof(double));

  struct residual_exchange_request
  {
    int task;
    int index;
    int row;
    double val;
  } * res_exch_req_send_unordered, *res_exch_req_send, *res_exch_req_get;
  struct residual_exchange_result
  {
    int row;
    double val;
    double x;
  } * res_exch_res_send, *res_exch_res_get;

  res_exch_req_send_unordered = (struct residual_exchange_request *)mymalloc("res_exch_req_send_unordered",
                                                                             md->cd_count * sizeof(struct residual_exchange_request));
  int Nsend                   = 0;

  int row;
  int colCount = 0;
  for(row = 0; row < hd->nRows; row++)
    {
      residuals[row] = -QradBackup[row];  // Old Qrad values

      int col;
      for(col = 0; col < hd->columnCounts[row]; col++)
        {
          int task = hd->tasks[colCount];
          if(task == ThisTask)
            {
              int particle = md->local_part_index[hd->cols[colCount] - md->offsets[task]];
              residuals[row] += hd->vals[colCount] * SphP[particle].Qrad;
            }
          else
            {
              res_exch_req_send_unordered[Nsend].task  = task;
              res_exch_req_send_unordered[Nsend].index = hd->cols[colCount] - md->offsets[task];
              res_exch_req_send_unordered[Nsend].row   = row;
              res_exch_req_send_unordered[Nsend].val   = hd->vals[colCount];
              Nsend++;
            }
          colCount++;
        }
    }

  int sendCount[NTask];
  int s, t;
  for(t = 0; t < NTask; t++)
    sendCount[t] = 0;
  for(s = 0; s < Nsend; s++)
    sendCount[res_exch_req_send_unordered[s].task]++;

  int sendOffset[NTask];
  sendOffset[0] = 0;
  for(t = 1; t < NTask; t++)
    sendOffset[t] = sendOffset[t - 1] + sendCount[t - 1];

  for(t = 0; t < NTask; t++)
    sendCount[t] = 0;

  res_exch_req_send =
      (struct residual_exchange_request *)mymalloc("res_exch_req_send", Nsend * sizeof(struct residual_exchange_request));

  for(s = 0; s < Nsend; s++)
    {
      int task  = res_exch_req_send_unordered[s].task;
      int index = sendOffset[task] + sendCount[task];
      sendCount[task]++;
      res_exch_req_send[index] = res_exch_req_send_unordered[s];
    }

  int recvCount[NTask];
  MPI_Alltoall(sendCount, 1, MPI_INT, recvCount, 1, MPI_INT, MPI_COMM_WORLD);

  int recvOffset[NTask];
  recvOffset[0] = 0;
  for(t = 1; t < NTask; t++)
    recvOffset[t] = recvOffset[t - 1] + recvCount[t - 1];

  int nimport = recvOffset[NTask - 1] + recvCount[NTask - 1];

  res_exch_req_get =
      (struct residual_exchange_request *)mymalloc("res_exch_req_get", nimport * sizeof(struct residual_exchange_request));

  int ngrp;
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(sendCount[recvTask] > 0 || recvCount[recvTask] > 0)
            {
              MPI_Sendrecv(&res_exch_req_send[sendOffset[recvTask]], sendCount[recvTask] * sizeof(struct residual_exchange_request),
                           MPI_BYTE, recvTask, TAG_DENS_A, &res_exch_req_get[recvOffset[recvTask]],
                           recvCount[recvTask] * sizeof(struct residual_exchange_request), MPI_BYTE, recvTask, TAG_DENS_A,
                           MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  res_exch_res_send =
      (struct residual_exchange_result *)mymalloc("res_exch_res_send", nimport * sizeof(struct residual_exchange_result));

  int i;
  for(i = 0; i < nimport; i++)
    {
      res_exch_res_send[i].row = res_exch_req_get[i].row;
      res_exch_res_send[i].val = res_exch_req_get[i].val;

      int particle           = md->local_part_index[res_exch_req_get[i].index];
      res_exch_res_send[i].x = SphP[particle].Qrad;
    }

  res_exch_res_get =
      (struct residual_exchange_result *)mymalloc("res_exch_res_get", nimport * sizeof(struct residual_exchange_result));

  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(sendCount[recvTask] > 0 || recvCount[recvTask] > 0)
            {
              MPI_Sendrecv(&res_exch_res_send[recvOffset[recvTask]], recvCount[recvTask] * sizeof(struct residual_exchange_request),
                           MPI_BYTE, recvTask, TAG_DENS_A, &res_exch_res_get[sendOffset[recvTask]],
                           sendCount[recvTask] * sizeof(struct residual_exchange_request), MPI_BYTE, recvTask, TAG_DENS_A,
                           MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  for(i = 0; i < Nsend; i++)
    {
      residuals[res_exch_res_get[i].row] += res_exch_res_get[i].val * res_exch_res_get[i].x;
    }

  myfree(res_exch_res_get);
  myfree(res_exch_res_send);
  myfree(res_exch_req_get);
  myfree(res_exch_req_send);
  myfree(res_exch_req_send_unordered);

  double ressum = 0;
  double bsum   = 0;

  for(row = 0; row < hd->nRows; row++)
    {
      ressum += residuals[row] * residuals[row];
      bsum += hd->bval[row] * hd->bval[row];

      if(residuals[row] > residual_max)
        {
          residual_max   = residuals[row];
          residual_b     = hd->bval[row];
          residual_index = md->local_part_index[row];
        }
    }

  myfree(residuals);

  double resall, ball;
  MPI_Allreduce(&ressum, &resall, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&bsum, &ball, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  double res = sqrt(resall / ball);
  mpi_printf("RADIATIVE_TRANSFER: Total residual: %g\n", res);

  int gSuccess;
  MPI_Allreduce(success, &gSuccess, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  *success = gSuccess;

  if((res < 1e-6) && !(*success))
    {
      mpi_printf("Residual is larger than %g, but still smaller than %g, so we accept it and continue.\n", RADTRANS_ACCURACY, 1e-6);
      *success = 1;
    }

  struct
  {
    double val;
    int rank;
  } in, out;

  in.val  = residual_max;
  in.rank = ThisTask;

  MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

  if(ThisTask == out.rank && residual_index != -1)
    {
      int p = residual_index;
      printf(
          "RADIATIVE_TRANSFER: Particle %d on task %d has the highest residual res=%g, b=%g. Qrad=%g (before=%g), Volume=%g, "
          "rho=%g, ID=%lld, Pos=%g,%g,%g, Qrad*V/Etot=%g\n",
          p, ThisTask, residual_max, residual_b, SphP[p].Qrad, hd->bval[md->local_row_index[p]], SphP[p].Volume, SphP[p].Density,
          (long long)P[p].ID, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2], SphP[p].Qrad * SphP[p].Volume / SphP[p].Energy);
    }
}

void apply_eddington_qrad(HYPRE_IJVector *x, struct matrix_data *md, double *dEUpdate)
{
  // Extract solution from HYPRE and update the cell energies
  int rows[md->local_row_count];
  double xval[md->local_row_count];

  int row;
  for(row = 0; row < md->local_row_count; row++)
    rows[row] = row + md->offsets[ThisTask];

  HYPRE_IJVectorGetValues(*x, md->local_row_count, rows, xval);

  for(row = 0; row < md->local_row_count; row++)
    {
      int particle  = md->local_part_index[row];
      double dt     = 0.5 * (All.HighestActiveTimeBin ? (((integertime)1) << All.HighestActiveTimeBin) : 0) * All.Timebase_interval;
      dEUpdate[row] = xval[row] * SphP[particle].Volume * dt;

      SphP[particle].Qrad = xval[row];
      SphP[particle].Energy += dEUpdate[row];
    }
}

void free_hypre_data(struct hypre_data *hd)
{
  if(hd->nRows > 0)
    {
      myfree(hd->tasks);
      myfree(hd->xval);
      myfree(hd->vals);
      myfree(hd->cols);
      myfree(hd->rows);
      myfree(hd->columnCounts);
      myfree(hd->bval);
    }
}

void free_rows_offsets_and_indizes(struct matrix_data *md)
{
  // Free arrays of md
  myfree(md->offsets);
  myfree(md->imported_particle_indizes);
  myfree(md->column_data);
  myfree(md->last_column);
  myfree(md->first_column);
  myfree(md->local_part_index);
  myfree(md->local_row_index);
}

#endif
