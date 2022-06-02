/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/pm/pm_non_periodic.c
 * \date        MM/YYYY
 * \author
 * \brief       code for non-periodic FFT to compute long-range PM force
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

#if defined(PMGRID) && (defined(PLACEHIGHRESREGION) || defined(GRAVITY_NOT_PERIODIC))

#if defined(LONG_X) || defined(LONG_Y) || defined(LONG_Z)
#error "LONG_X/Y/Z not supported for the non-periodic FFT gravity code"
#endif
#ifdef GRAVITY_TALLBOX
#error "GRAVITY_TALLBOX not supported for the non-periodic FFT gravity code"
#endif

static fft_plan myplan; /*!< In this structure, various bookkeeping variables for the distributed FFTs are stored */

/*! \var maxfftsize
 *  \brief maximum size of the local fft grid among all tasks
 */
static size_t maxfftsize;

/*! \var rhogrid
 *  \brief This array hold the local part of the density field and
 *  after the FFTs the local part of the potential
 *
 *  \var forcegrid
 *  \brief This array will contain the force field
 *
 *  \var workspace
 *  \brief Workspace array used during the FFTs
 */
static fft_real *rhogrid, *forcegrid, *workspace;

/*! \brief Array containing the FFT of #rhogrid
 *
 *  This pointer points to the same array as #rhogrid,
 *  because in-place FFTs are used.
 */
static fft_complex *fft_of_rhogrid;

static fft_real *kernel[2];
static fft_complex *fft_of_kernel[2];

/* short-cut functions for accessing different 3D arrays */
static inline large_array_offset FI(const int x, const int y, const int z) { return PM_FI(GRID_NP, GRID2_NP, x, y, z); }
#ifdef FFT_COLUMN_BASED
static inline large_array_offset FC(const int c, const int z) { return PM_FC(myplan.base_firstcol, GRID2_NP, c, z); }
#else
static inline large_array_offset TI(const int x, const int y, const int z) { return PM_TI(myplan.nslab_x, GRID_NP, x, y, z); }
#endif

/*! \param Determine particle extent.
 *
 *  This function determines the particle extension of all particles, and for
 *  those types selected with PLACEHIGHRESREGION if this is used, and then
 *  determines the boundaries of the non-periodic FFT-mesh that can be placed
 *  on this region. Note that a sufficient buffer region at the rim of the
 *  occupied part of the mesh needs to be reserved in order to allow a correct
 *  finite differencing using a 4-point formula. In addition, to allow
 *  non-periodic boundaries, the actual FFT mesh used is twice as large in
 *  each dimension compared with GRID_NP.
 *
 *  \return void
 */
void pm_init_regionsize(void)
{
  double meshinner[2], xmin[2][3], xmax[2][3];
  int i, j;

  /* find enclosing rectangle */

  for(j = 0; j < 3; j++)
    {
      xmin[0][j] = xmin[1][j] = 1.0e36;
      xmax[0][j] = xmax[1][j] = -1.0e36;
    }

  for(i = 0; i < NumPart; i++)
    for(j = 0; j < 3; j++)
      {
        if(P[i].Pos[j] > xmax[0][j])
          xmax[0][j] = P[i].Pos[j];
        if(P[i].Pos[j] < xmin[0][j])
          xmin[0][j] = P[i].Pos[j];

#ifdef PLACEHIGHRESREGION
        if(((1 << P[i].Type) & (PLACEHIGHRESREGION)))
          {
            if(P[i].Pos[j] > xmax[1][j])
              xmax[1][j] = P[i].Pos[j];
            if(P[i].Pos[j] < xmin[1][j])
              xmin[1][j] = P[i].Pos[j];
          }
#endif
      }

  MPI_Allreduce(xmin, All.Xmintot, 6, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(xmax, All.Xmaxtot, 6, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  for(j = 0; j < 2; j++)
    {
      All.TotalMeshSize[j] = All.Xmaxtot[j][0] - All.Xmintot[j][0];
      All.TotalMeshSize[j] = fmax(All.TotalMeshSize[j], All.Xmaxtot[j][1] - All.Xmintot[j][1]);
      All.TotalMeshSize[j] = fmax(All.TotalMeshSize[j], All.Xmaxtot[j][2] - All.Xmintot[j][2]);
#ifdef ENLARGEREGION
      All.TotalMeshSize[j] *= ENLARGEREGION;
#endif

      /* symmetrize the box onto the center */
      for(i = 0; i < 3; i++)
        {
          All.Xmintot[j][i] = (All.Xmintot[j][i] + All.Xmaxtot[j][i]) / 2 - All.TotalMeshSize[j] / 2;
          All.Xmaxtot[j][i] = All.Xmintot[j][i] + All.TotalMeshSize[j];
        }
    }

  /* this will produce enough room for zero-padding and buffer region to
     allow finite differencing of the potential  */

  for(j = 0; j < 2; j++)
    {
      meshinner[j] = All.TotalMeshSize[j];
      All.TotalMeshSize[j] *= 2.001 * (GRID_NP) / ((double)(GRID_NP - 2 - 8));
    }

  /* move lower left corner by two cells to allow finite differencing of the potential by a 4-point function */

  for(j = 0; j < 2; j++)
    for(i = 0; i < 3; i++)
      {
        All.Corner[j][i]      = All.Xmintot[j][i] - 2.0005 * All.TotalMeshSize[j] / GRID_NP;
        All.UpperCorner[j][i] = All.Corner[j][i] + (GRID_NP / 2 - 1) * (All.TotalMeshSize[j] / GRID_NP);
      }

#ifdef PLACEHIGHRESREGION
  All.Asmth[1] = ASMTH * All.TotalMeshSize[1] / GRID_NP;
  All.Rcut[1]  = RCUT * All.Asmth[1];

  if(2 * All.TotalMeshSize[1] / GRID_NP < All.Rcut[0])
    {
      All.TotalMeshSize[1] = 2 * (meshinner[1] + 2 * All.Rcut[0]) * (GRID_NP) / ((double)(GRID_NP - 2));

      for(i = 0; i < 3; i++)
        {
          All.Corner[1][i]      = All.Xmintot[1][i] - 1.0001 * All.Rcut[0];
          All.UpperCorner[1][i] = All.Corner[1][i] + (GRID_NP / 2 - 1) * (All.TotalMeshSize[1] / GRID_NP);
        }

      if(2 * All.TotalMeshSize[1] / GRID_NP > All.Rcut[0])
        {
          All.TotalMeshSize[1] = 2 * (meshinner[1] + 2 * All.Rcut[0]) * (GRID_NP) / ((double)(GRID_NP - 10));

          for(i = 0; i < 3; i++)
            {
              All.Corner[1][i]      = All.Xmintot[1][i] - 1.0001 * (All.Rcut[0] + 2 * All.TotalMeshSize[1] / GRID_NP);
              All.UpperCorner[1][i] = All.Corner[1][i] + (GRID_NP / 2 - 1) * (All.TotalMeshSize[1] / GRID_NP);
            }
        }

      All.Asmth[1] = ASMTH * All.TotalMeshSize[1] / GRID_NP;
      All.Rcut[1]  = RCUT * All.Asmth[1];

      mpi_printf("PM-NONPERIODIC: All.Asmth[0]=%g All.Asmth[1]=%g\n", All.Asmth[0], All.Asmth[1]);
    }

  if(All.TotalMeshSize[1] > 0.9 * All.BoxSize)
    {
      if(ThisTask == 0)
        terminate(
            "PM-NONPERIODIC: meshsize=%g, boxsize=%g, does not make any sense for high res FFT. Likely the high res region overlaps "
            "with the boundary of the background periodic box, which is not allowed. To fix this move your high res region to the "
            "center of the box.",
            All.TotalMeshSize[1], All.BoxSize);
    }

  mpi_printf(
      "PM-NONPERIODIC: Allowed region for isolated PM mesh (high-res): (%g|%g|%g)  -> (%g|%g|%g)   ext=%g  totmeshsize=%g  "
      "meshsize=%g\n\n",
      All.Xmintot[1][0], All.Xmintot[1][1], All.Xmintot[1][2], All.Xmaxtot[1][0], All.Xmaxtot[1][1], All.Xmaxtot[1][2], meshinner[1],
      All.TotalMeshSize[1], All.TotalMeshSize[1] / GRID_NP);
#endif
}

/*! \brief Initialization of the non-periodic PM routines.
 *
 *  The plan-files for FFTW are created. Finally, the routine to set-up the
 *  non-periodic Greens function is called.
 *
 *  \return void
 */
void pm_init_nonperiodic(void)
{
  /* set up the FFTW-3 plan files */
  const int ndim[1] = {GRID_NP}; /* dimension of the 1D transforms */

#ifndef FFT_COLUMN_BASED
  const int stride = GRIDz_NP;
#else
  const int stride = 1;
#endif

#ifdef DOUBLEPRECISION_FFTW
  const unsigned int alignflag = 0;
#else
  /* for single precision, the start of our FFT columns is presently only guaranteed to be 8-byte aligned */
  const unsigned int alignflag = FFTW_UNALIGNED;
#endif

  const unsigned int flags = FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag;

  /* temporarily allocate some arrays to make sure that out-of-place plans are created */
  const size_t tmp_len    = 2 * GRID_NP * (size_t)stride;
  fft_real *rhogrid_tmp   = (fft_real *)mymalloc("rhogrid_tmp", tmp_len * sizeof(*rhogrid_tmp));
  fft_real *forcegrid_tmp = (fft_real *)mymalloc("forcegrid_tmp", tmp_len * sizeof(*forcegrid_tmp));

  myplan.forward_plan_zdir =
      FFTW(plan_many_dft_r2c)(1, ndim, 1, rhogrid_tmp, NULL, 1, 0, (fft_complex *)forcegrid_tmp, NULL, 1, 0, flags);

  myplan.forward_plan_xdir = FFTW(plan_many_dft)(1, ndim, 1, (fft_complex *)rhogrid_tmp, 0, stride, 0, (fft_complex *)forcegrid_tmp,
                                                 NULL, stride, 0, FFTW_FORWARD, flags);

  myplan.forward_plan_ydir = FFTW(plan_many_dft)(1, ndim, 1, (fft_complex *)rhogrid_tmp, 0, stride, 0, (fft_complex *)forcegrid_tmp,
                                                 NULL, stride, 0, FFTW_FORWARD, flags);

  myplan.backward_plan_zdir =
      FFTW(plan_many_dft_c2r)(1, ndim, 1, (fft_complex *)rhogrid_tmp, NULL, 1, 0, forcegrid_tmp, NULL, 1, 0, flags);

  myplan.backward_plan_xdir = FFTW(plan_many_dft)(1, ndim, 1, (fft_complex *)rhogrid_tmp, NULL, stride, 0,
                                                  (fft_complex *)forcegrid_tmp, NULL, stride, 0, FFTW_BACKWARD, flags);

  myplan.backward_plan_ydir = FFTW(plan_many_dft)(1, ndim, 1, (fft_complex *)rhogrid_tmp, NULL, stride, 0,
                                                  (fft_complex *)forcegrid_tmp, NULL, stride, 0, FFTW_BACKWARD, flags);

  myfree(forcegrid_tmp);
  myfree(rhogrid_tmp);

#ifndef FFT_COLUMN_BASED
  my_slab_based_fft_init(&myplan, GRID_NP, GRID_NP, GRID_NP);
  maxfftsize = myplan.largest_x_slab * GRID_NP * ((size_t)GRID2_NP);
#else
  my_column_based_fft_init(&myplan, GRID_NP, GRID_NP, GRID_NP);
  maxfftsize = myplan.max_datasize;
#endif

  /* now allocate memory to hold the FFT fields */

  size_t bytes, bytes_tot = 0;

#ifdef GRAVITY_NOT_PERIODIC
  kernel[0] = (fft_real *)mymalloc("kernel[0]", bytes = maxfftsize * sizeof(fft_real));
  bytes_tot += bytes;
  fft_of_kernel[0] = (fft_complex *)kernel[0];
#endif

#ifdef PLACEHIGHRESREGION
  kernel[1] = (fft_real *)mymalloc("kernel[1]", bytes = maxfftsize * sizeof(fft_real));
  bytes_tot += bytes;
  fft_of_kernel[1] = (fft_complex *)kernel[1];
#endif

  mpi_printf("\nPM-NONPERIODIC: Allocated %g MByte for FFT kernel(s).\n\n", bytes_tot / (1024.0 * 1024.0));
}

#ifdef PLACEHIGHRESREGION
/*! \brief Is this a high res particle in high resolution region?
 *
 *  For cosmological zoom simulations.
 *
 *  \param[in] type Parcile type.
 *  \param[in] Pos Position of particle.
 *
 *  \return 0: not high res; 1: high res.
 */
int pmforce_is_particle_high_res(int type, MyDouble *Pos)
{
  int flag = 1;

  if((1 << type) & (PLACEHIGHRESREGION))
    return 1;

#if defined(PLACEHIGHRESREGION) && defined(FORCETEST_TESTFORCELAW) && (FORCETEST_TESTFORCELAW == 1)
  double r2 = 0;
  for(int j = 0; j < 3; j++)
    r2 += pow(Pos[j] - 0.5 * (All.Xmintot[1][j] + All.Xmaxtot[1][j]), 2);

  if(sqrt(r2) > 0.5 * (All.Xmaxtot[1][0] - All.Xmintot[1][0]))
    return 0;
#else

  for(int j = 0; j < 3; j++)
    if(Pos[j] < All.Xmintot[1][j] || Pos[j] > All.Xmaxtot[1][j])
      {
        flag = 0; /* we are outside */
        break;
      }

#endif

  return flag;
}
#endif

#ifdef PM_ZOOM_OPTIMIZED

static void mysort_pmperiodic(void *b, size_t n, size_t s, int (*cmp)(const void *, const void *));
static int pm_periodic_compare_sortindex(const void *a, const void *b);

/*! \brief This structure links the particles to the mesh cells, to which they contribute their mass
 *
 * Each particle will have eight items of this structure in the #part array.
 * For each of the eight mesh cells the CIC assignment will contribute,
 * one item of this struct exists.
 */
static struct part_slab_data
{
  large_array_offset globalindex; /*!< index in the global density mesh */
  large_numpart_type partindex; /*!< contains the local particle index shifted by 2^3, the first three bits encode to which part of the
                                   CIC assignment this item belongs to */
  large_array_offset localindex; /*!< index to a local copy of the corresponding mesh cell of the global density array (used during
                                    local mass and force assignment) */
} * part;                        /*!< array of part_slab_data linking the local particles to their mesh cells */

static size_t *localfield_sendcount, *localfield_first, *localfield_offset, *localfield_recvcount;
static large_array_offset *localfield_globalindex, *import_globalindex;
static fft_real *localfield_data, *import_data;
static large_numpart_type num_on_grid;

/*! \brief Prepares density field for nonperiodic FFTs.
 *
 *  \param[in] grnr (0, 1) 0 if full mesh, 1 if highres grid.
 *
 *  \return void
 */
void pmforce_nonperiodic_zoom_optimized_prepare_density(int grnr)
{
  large_numpart_type i;
  int level, recvTask;
  MPI_Status status;

  double to_slab_fac = GRID_NP / All.TotalMeshSize[grnr];

  part                               = (struct part_slab_data *)mymalloc("part", 8 * (NumPart * sizeof(struct part_slab_data)));
  large_numpart_type *part_sortindex = (large_numpart_type *)mymalloc("part_sortindex", 8 * (NumPart * sizeof(large_numpart_type)));

  int ngrid = 0;

  /* determine the cells each particle accesses */
#pragma omp parallel for
  for(i = 0; i < NumPart; i++)
    {
      MyDouble *pos;

#ifdef CELL_CENTER_GRAVITY
      if(P[i].Type == 0)
        pos = SphP[i].Center;
      else
#endif
        pos = P[i].Pos;

      if(pos[0] < All.Corner[grnr][0] || pos[0] >= All.UpperCorner[grnr][0])
        continue;
      if(pos[1] < All.Corner[grnr][1] || pos[1] >= All.UpperCorner[grnr][1])
        continue;
      if(pos[2] < All.Corner[grnr][2] || pos[2] >= All.UpperCorner[grnr][2])
        continue;

      int slab_x = (int)(to_slab_fac * (pos[0] - All.Corner[grnr][0]));
      int slab_y = (int)(to_slab_fac * (pos[1] - All.Corner[grnr][1]));
      int slab_z = (int)(to_slab_fac * (pos[2] - All.Corner[grnr][2]));
      int myngrid;

#pragma omp atomic capture
      {
        myngrid = ngrid;
        ngrid += 1;
      }

      large_numpart_type index_on_grid = ((large_numpart_type)myngrid) * 8;

      int xx, yy, zz;

      for(xx = 0; xx < 2; xx++)
        for(yy = 0; yy < 2; yy++)
          for(zz = 0; zz < 2; zz++)
            {
              int slab_xx = slab_x + xx;
              int slab_yy = slab_y + yy;
              int slab_zz = slab_z + zz;

              if(slab_xx >= GRID_NP)
                slab_xx -= GRID_NP;
              if(slab_yy >= GRID_NP)
                slab_yy -= GRID_NP;
              if(slab_zz >= GRID_NP)
                slab_zz -= GRID_NP;

              large_array_offset offset = FI(slab_xx, slab_yy, slab_zz);

              part[index_on_grid].partindex   = (i << 3) + (xx << 2) + (yy << 1) + zz;
              part[index_on_grid].globalindex = offset;
              part_sortindex[index_on_grid]   = index_on_grid;
              index_on_grid++;
            }
    }

  /* note: num_on_grid will be 8 times larger than the particle number, but
   * num_field_points will generally be much smaller */
  num_on_grid = ((large_numpart_type)ngrid) * 8;

  /* bring the part-field into the order of the accessed cells. This allows the removal of duplicates */
  mysort_pmperiodic(part_sortindex, num_on_grid, sizeof(large_numpart_type), pm_periodic_compare_sortindex);

  large_array_offset num_field_points;

  if(num_on_grid > 0)
    num_field_points = 1;
  else
    num_field_points = 0;

    /* determine the number of unique field points */
#pragma omp parallel for reduction(+ : num_field_points)
  for(i = 1; i < num_on_grid; i++)
    {
      if(part[part_sortindex[i]].globalindex != part[part_sortindex[i - 1]].globalindex)
        num_field_points++;
    }

  /* allocate the local field */
  localfield_globalindex = (large_array_offset *)mymalloc_movable(&localfield_globalindex, "localfield_globalindex",
                                                                  num_field_points * sizeof(large_array_offset));
  localfield_data        = (fft_real *)mymalloc_movable(&localfield_data, "localfield_data", num_field_points * sizeof(fft_real));
  localfield_first       = (size_t *)mymalloc_movable(&localfield_first, "localfield_first", NTask * sizeof(size_t));
  localfield_sendcount   = (size_t *)mymalloc_movable(&localfield_sendcount, "localfield_sendcount", NTask * sizeof(size_t));
  localfield_offset      = (size_t *)mymalloc_movable(&localfield_offset, "localfield_offset", NTask * sizeof(size_t));
  localfield_recvcount   = (size_t *)mymalloc_movable(&localfield_recvcount, "localfield_recvcount", NTask * sizeof(size_t));

  for(i = 0; i < NTask; i++)
    {
      localfield_first[i]     = 0;
      localfield_sendcount[i] = 0;
    }

  /* establish the cross link between the part[] array and the local list of
   * mesh points. Also, count on which CPU the needed field points are stored. */
  for(i = 0, num_field_points = 0; i < num_on_grid; i++)
    {
      if(i > 0)
        if(part[part_sortindex[i]].globalindex != part[part_sortindex[i - 1]].globalindex)
          num_field_points++;

      part[part_sortindex[i]].localindex = num_field_points;

      if(i > 0)
        if(part[part_sortindex[i]].globalindex == part[part_sortindex[i - 1]].globalindex)
          continue;

      localfield_globalindex[num_field_points] = part[part_sortindex[i]].globalindex;

#ifndef FFT_COLUMN_BASED
      int slab = part[part_sortindex[i]].globalindex / (GRID_NP * GRID2_NP);
      int task = myplan.slab_to_task[slab];
#else
      int task, column = part[part_sortindex[i]].globalindex / (GRID2_NP);

      if(column < myplan.pivotcol)
        task = column / myplan.avg;
      else
        task = (column - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;
#endif

      if(localfield_sendcount[task] == 0)
        localfield_first[task] = num_field_points;

      localfield_sendcount[task]++;
    }
  if(num_on_grid > 0)
    num_field_points++;

  for(i = 1, localfield_offset[0] = 0; i < NTask; i++)
    localfield_offset[i] = localfield_offset[i - 1] + localfield_sendcount[i - 1];

  myfree_movable(part_sortindex);
  part_sortindex = NULL;

  /* now bin the local particle data onto the mesh list */
#pragma omp parallel for
  for(i = 0; i < (large_numpart_type)num_field_points; i++)
    localfield_data[i] = 0;

  for(i = 0; i < num_on_grid; i += 8)
    {
      int pindex = part[i].partindex >> 3;

      MyDouble *pos;

#ifdef CELL_CENTER_GRAVITY
      if(P[pindex].Type == 0)
        pos = SphP[pindex].Center;
      else
#endif
        pos = P[pindex].Pos;

      int slab_x = (int)(to_slab_fac * (pos[0] - All.Corner[grnr][0]));
      int slab_y = (int)(to_slab_fac * (pos[1] - All.Corner[grnr][1]));
      int slab_z = (int)(to_slab_fac * (pos[2] - All.Corner[grnr][2]));

      double dx = to_slab_fac * (pos[0] - All.Corner[grnr][0]) - slab_x;
      double dy = to_slab_fac * (pos[1] - All.Corner[grnr][1]) - slab_y;
      double dz = to_slab_fac * (pos[2] - All.Corner[grnr][2]) - slab_z;

      double weight = P[pindex].Mass;

      localfield_data[part[i + 0].localindex] += weight * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
      localfield_data[part[i + 1].localindex] += weight * (1.0 - dx) * (1.0 - dy) * dz;
      localfield_data[part[i + 2].localindex] += weight * (1.0 - dx) * dy * (1.0 - dz);
      localfield_data[part[i + 3].localindex] += weight * (1.0 - dx) * dy * dz;
      localfield_data[part[i + 4].localindex] += weight * (dx) * (1.0 - dy) * (1.0 - dz);
      localfield_data[part[i + 5].localindex] += weight * (dx) * (1.0 - dy) * dz;
      localfield_data[part[i + 6].localindex] += weight * (dx)*dy * (1.0 - dz);
      localfield_data[part[i + 7].localindex] += weight * (dx)*dy * dz;
    }

  rhogrid = (fft_real *)mymalloc("rhogrid", maxfftsize * sizeof(fft_real));

  /* clear local FFT-mesh density field */
  large_array_offset ii;
#pragma omp parallel for
  for(ii = 0; ii < maxfftsize; ii++)
    rhogrid[ii] = 0;

  /* exchange data and add contributions to the local mesh-path */
  MPI_Alltoall(localfield_sendcount, sizeof(size_t), MPI_BYTE, localfield_recvcount, sizeof(size_t), MPI_BYTE, MPI_COMM_WORLD);

  for(level = 0; level < (1 << PTask); level++) /* note: for level=0, target is the same task */
    {
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
        {
          if(level > 0)
            {
              import_data = (fft_real *)mymalloc("import_data", localfield_recvcount[recvTask] * sizeof(fft_real));
              import_globalindex =
                  (large_array_offset *)mymalloc("import_globalindex", localfield_recvcount[recvTask] * sizeof(large_array_offset));

              if(localfield_sendcount[recvTask] > 0 || localfield_recvcount[recvTask] > 0)
                {
                  myMPI_Sendrecv(localfield_data + localfield_offset[recvTask], localfield_sendcount[recvTask] * sizeof(fft_real),
                                 MPI_BYTE, recvTask, TAG_NONPERIOD_A, import_data, localfield_recvcount[recvTask] * sizeof(fft_real),
                                 MPI_BYTE, recvTask, TAG_NONPERIOD_A, MPI_COMM_WORLD, &status);

                  myMPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
                                 localfield_sendcount[recvTask] * sizeof(large_array_offset), MPI_BYTE, recvTask, TAG_NONPERIOD_B,
                                 import_globalindex, localfield_recvcount[recvTask] * sizeof(large_array_offset), MPI_BYTE, recvTask,
                                 TAG_NONPERIOD_B, MPI_COMM_WORLD, &status);
                }
            }
          else
            {
              import_data        = localfield_data + localfield_offset[ThisTask];
              import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
            }

            /* note: here every element in rhogrid is only accessed once, so there should be no race condition */
#pragma omp parallel for
          for(i = 0; i < (large_numpart_type)localfield_recvcount[recvTask]; i++)
            {
              /* determine offset in local FFT slab */
#ifndef FFT_COLUMN_BASED
              large_array_offset offset =
                  import_globalindex[i] - myplan.first_slab_x_of_task[ThisTask] * GRID_NP * ((large_array_offset)GRID2_NP);
#else
              large_array_offset offset = import_globalindex[i] - myplan.base_firstcol * ((large_array_offset)GRID2_NP);
#endif
              rhogrid[offset] += import_data[i];
            }

          if(level > 0)
            {
              myfree(import_globalindex);
              myfree(import_data);
            }
        }
    }
#ifdef MODGRAV_EFF_MASS
  modgrav_add_amr_nodes_to_nonperiodic_grid(grnr, rhogrid, to_slab_fac, myplan, maxfftsize);
#endif
}

/*! \brief Reads out the force component corresponding to spatial dimension
 *         'dim'.
 *
 *  If dim is negative, potential values are read out and assigned to
 *  particles.
 *
 *  \param[in] grnr Number of grid (0: base, 1 high-res)
 *  \param[in] dim Dimension to be read out
 *             (<0: potential,>=0 force component).
 *
 *  \return void
 */
void pmforce_nonperiodic_zoom_optimized_readout_forces_or_potential(int grnr, int dim)
{
#ifdef EVALPOTENTIAL
  /* factor to get potential */
  double fac = All.G / pow(All.TotalMeshSize[grnr], 4) * pow(All.TotalMeshSize[grnr] / GRID_NP, 3);
#endif

  large_numpart_type i;
  int level, recvTask;
  MPI_Status status;

  fft_real *grid;

  if(dim < 0)
    grid = rhogrid;
  else
    grid = forcegrid;

  double to_slab_fac = GRID_NP / All.TotalMeshSize[grnr];

  for(level = 0; level < (1 << PTask); level++) /* note: for level=0, target is the same task */
    {
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
        {
          if(level > 0)
            {
              import_data = (fft_real *)mymalloc("import_data", localfield_recvcount[recvTask] * sizeof(fft_real));
              import_globalindex =
                  (large_array_offset *)mymalloc("import_globalindex", localfield_recvcount[recvTask] * sizeof(large_array_offset));

              if(localfield_sendcount[recvTask] > 0 || localfield_recvcount[recvTask] > 0)
                {
                  myMPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
                                 localfield_sendcount[recvTask] * sizeof(large_array_offset), MPI_BYTE, recvTask, TAG_NONPERIOD_C,
                                 import_globalindex, localfield_recvcount[recvTask] * sizeof(large_array_offset), MPI_BYTE, recvTask,
                                 TAG_NONPERIOD_C, MPI_COMM_WORLD, &status);
                }
            }
          else
            {
              import_data        = localfield_data + localfield_offset[ThisTask];
              import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
            }

          for(i = 0; i < (large_numpart_type)localfield_recvcount[recvTask]; i++)
            {
#ifndef FFT_COLUMN_BASED
              large_array_offset offset =
                  import_globalindex[i] - myplan.first_slab_x_of_task[ThisTask] * GRID_NP * ((large_array_offset)GRID2_NP);
#else
              large_array_offset offset = import_globalindex[i] - myplan.base_firstcol * ((large_array_offset)GRID2_NP);
#endif
              import_data[i] = grid[offset];
            }

          if(level > 0)
            {
              myMPI_Sendrecv(import_data, localfield_recvcount[recvTask] * sizeof(fft_real), MPI_BYTE, recvTask, TAG_NONPERIOD_A,
                             localfield_data + localfield_offset[recvTask], localfield_sendcount[recvTask] * sizeof(fft_real),
                             MPI_BYTE, recvTask, TAG_NONPERIOD_A, MPI_COMM_WORLD, &status);

              myfree(import_globalindex);
              myfree(import_data);
            }
        }
    }

  /* read out the force/potential values, which all have been assembled in localfield_data */

  int k, ngrid = (num_on_grid >> 3);

#pragma omp parallel for
  for(k = 0; k < ngrid; k++)
    {
      large_numpart_type j = (((large_numpart_type)k) << 3);

      int i = (part[j].partindex >> 3);

      MyDouble *pos;

#ifdef CELL_CENTER_GRAVITY
      if(P[i].Type == 0)
        pos = SphP[i].Center;
      else
#endif
        pos = P[i].Pos;

#ifdef PLACEHIGHRESREGION
      if(grnr == 1)
        if(!(pmforce_is_particle_high_res(P[i].Type, pos)))
          continue;
#endif

      int slab_x = (int)(to_slab_fac * (pos[0] - All.Corner[grnr][0]));
      double dx  = to_slab_fac * (pos[0] - All.Corner[grnr][0]) - slab_x;

      int slab_y = (int)(to_slab_fac * (pos[1] - All.Corner[grnr][1]));
      double dy  = to_slab_fac * (pos[1] - All.Corner[grnr][1]) - slab_y;

      int slab_z = (int)(to_slab_fac * (pos[2] - All.Corner[grnr][2]));
      double dz  = to_slab_fac * (pos[2] - All.Corner[grnr][2]) - slab_z;

      double value = +localfield_data[part[j + 0].localindex] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz) +
                     localfield_data[part[j + 1].localindex] * (1.0 - dx) * (1.0 - dy) * dz +
                     localfield_data[part[j + 2].localindex] * (1.0 - dx) * dy * (1.0 - dz) +
                     localfield_data[part[j + 3].localindex] * (1.0 - dx) * dy * dz +
                     localfield_data[part[j + 4].localindex] * (dx) * (1.0 - dy) * (1.0 - dz) +
                     localfield_data[part[j + 5].localindex] * (dx) * (1.0 - dy) * dz +
                     localfield_data[part[j + 6].localindex] * (dx)*dy * (1.0 - dz) +
                     localfield_data[part[j + 7].localindex] * (dx)*dy * dz;

      if(dim < 0)
        {
#ifdef EVALPOTENTIAL
          P[i].PM_Potential += value * fac;
#endif
        }
      else
        P[i].GravPM[dim] += value;
    }
}

#else

/*
 *  Here come the routines for a different communication algorithm that is better suited for a homogenuously loaded boxes.
 */
static struct partbuf
{
  MyFloat Mass;
  MyFloat Pos[3];
} * partin, *partout;

static size_t nimport, nexport;

static size_t *Sndpm_count, *Sndpm_offset;
static size_t *Rcvpm_count, *Rcvpm_offset;

/*! \brief Prepares density for pm calculation in algorithm optimized for
 *         uniform densities.
 *
 *  \param[in] grnr Number of grid (0: base grid, 1: high res grid).
 *
 *  \return void
 */
void pmforce_nonperiodic_uniform_optimized_prepare_density(int grnr)
{
  int i, j;

  double to_slab_fac = GRID_NP / All.TotalMeshSize[grnr];

  /* We here enlarge NTask such that each thread gets its own cache line for send_count/send_offset.
   * This should hopefully prevent a performance penalty from 'false sharing' for these variables */
  int multiNtask = roundup_to_multiple_of_cacheline_size(NTask * sizeof(size_t)) / sizeof(size_t);

  Sndpm_count  = (size_t *)mymalloc("Sndpm_count", MaxThreads * multiNtask * sizeof(size_t));
  Sndpm_offset = (size_t *)mymalloc("Sndpm_offset", MaxThreads * multiNtask * sizeof(size_t));
  Rcvpm_count  = (size_t *)mymalloc("Rcvpm_count", NTask * sizeof(size_t));
  Rcvpm_offset = (size_t *)mymalloc("Rcvpm_offset", NTask * sizeof(size_t));

  /* determine the slabs/columns each particles accesses */
#pragma omp parallel private(j)
  {
    size_t *send_count = Sndpm_count + get_thread_num() * multiNtask;

    /* each threads needs to do theloop to clear its send_count[] array */
    for(j = 0; j < NTask; j++)
      send_count[j] = 0;

#pragma omp for schedule(static) private(i)
    for(i = 0; i < NumPart; i++)
      {
        MyDouble *pos;

#ifdef CELL_CENTER_GRAVITY
        if(P[i].Type == 0)
          pos = SphP[i].Center;
        else
#endif
          pos = P[i].Pos;

        if(pos[0] < All.Corner[grnr][0] || pos[0] >= All.UpperCorner[grnr][0])
          continue;
        if(pos[1] < All.Corner[grnr][1] || pos[1] >= All.UpperCorner[grnr][1])
          continue;
        if(pos[2] < All.Corner[grnr][2] || pos[2] >= All.UpperCorner[grnr][2])
          continue;

        int slab_x  = (int)(to_slab_fac * (pos[0] - All.Corner[grnr][0]));
        int slab_xx = slab_x + 1;

#ifndef FFT_COLUMN_BASED
        int task0   = myplan.slab_to_task[slab_x];
        int task1   = myplan.slab_to_task[slab_xx];

        send_count[task0]++;
        if(task0 != task1)
          send_count[task1]++;
#else
        int slab_y  = (int)(to_slab_fac * (pos[1] - All.Corner[grnr][1]));
        int slab_yy = slab_y + 1;

        int column0 = slab_x * GRID_NP + slab_y;
        int column1 = slab_x * GRID_NP + slab_yy;
        int column2 = slab_xx * GRID_NP + slab_y;
        int column3 = slab_xx * GRID_NP + slab_yy;

        int task0, task1, task2, task3;

        if(column0 < myplan.pivotcol)
          task0 = column0 / myplan.avg;
        else
          task0 = (column0 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column1 < myplan.pivotcol)
          task1 = column1 / myplan.avg;
        else
          task1 = (column1 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column2 < myplan.pivotcol)
          task2 = column2 / myplan.avg;
        else
          task2 = (column2 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column3 < myplan.pivotcol)
          task3 = column3 / myplan.avg;
        else
          task3 = (column3 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        send_count[task0]++;
        if(task1 != task0)
          send_count[task1]++;
        if(task2 != task1 && task2 != task0)
          send_count[task2]++;
        if(task3 != task0 && task3 != task1 && task3 != task2)
          send_count[task3]++;
#endif
      }
  }

  /* collect thread-specific offset table and collect the results from the other threads */
  for(i = 0, Sndpm_offset[0] = 0; i < NTask; i++)
    for(j = 0; j < MaxThreads; j++)
      {
        int ind_prev, ind = j * multiNtask + i;
        if(ind > 0)
          {
            if(j == 0)
              ind_prev = (MaxThreads - 1) * multiNtask + i - 1;
            else
              ind_prev = ind - multiNtask;

            Sndpm_offset[ind] = Sndpm_offset[ind_prev] + Sndpm_count[ind_prev];
          }
      }

  for(j = 1; j < MaxThreads; j++)
    for(i = 0; i < NTask; i++)
      Sndpm_count[i] += Sndpm_count[i + j * multiNtask];

  MPI_Alltoall(Sndpm_count, sizeof(size_t), MPI_BYTE, Rcvpm_count, sizeof(size_t), MPI_BYTE, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, nexport = 0, Rcvpm_offset[0] = 0, Sndpm_offset[0] = 0; j < NTask; j++)
    {
      nexport += Sndpm_count[j];
      nimport += Rcvpm_count[j];

      if(j > 0)
        {
          Sndpm_offset[j] = Sndpm_offset[j - 1] + Sndpm_count[j - 1];
          Rcvpm_offset[j] = Rcvpm_offset[j - 1] + Rcvpm_count[j - 1];
        }
    }

  /* allocate import and export buffer */
  partin  = (struct partbuf *)mymalloc("partin", nimport * sizeof(struct partbuf));
  partout = (struct partbuf *)mymalloc("partout", nexport * sizeof(struct partbuf));

#pragma omp parallel private(j)
  {
    size_t *send_count  = Sndpm_count + get_thread_num() * multiNtask;
    size_t *send_offset = Sndpm_offset + get_thread_num() * multiNtask;

    for(j = 0; j < NTask; j++)
      send_count[j] = 0;

      /* fill export buffer */
#pragma omp for schedule(static)
    for(i = 0; i < NumPart; i++)
      {
        MyDouble *pos;

#ifdef CELL_CENTER_GRAVITY
        if(P[i].Type == 0)
          pos = SphP[i].Center;
        else
#endif
          pos = P[i].Pos;

        if(pos[0] < All.Corner[grnr][0] || pos[0] >= All.UpperCorner[grnr][0])
          continue;
        if(pos[1] < All.Corner[grnr][1] || pos[1] >= All.UpperCorner[grnr][1])
          continue;
        if(pos[2] < All.Corner[grnr][2] || pos[2] >= All.UpperCorner[grnr][2])
          continue;

        int slab_x  = (int)(to_slab_fac * (pos[0] - All.Corner[grnr][0]));
        int slab_xx = slab_x + 1;

#ifndef FFT_COLUMN_BASED
        int task0   = myplan.slab_to_task[slab_x];
        int task1   = myplan.slab_to_task[slab_xx];

        size_t ind0        = send_offset[task0] + send_count[task0]++;
        partout[ind0].Mass = P[i].Mass;
        for(j = 0; j < 3; j++)
          partout[ind0].Pos[j] = pos[j];

        if(task0 != task1)
          {
            size_t ind1        = send_offset[task1] + send_count[task1]++;
            partout[ind1].Mass = P[i].Mass;
            for(j = 0; j < 3; j++)
              partout[ind1].Pos[j] = pos[j];
          }
#else
        int slab_y  = (int)(to_slab_fac * (pos[1] - All.Corner[grnr][1]));
        int slab_yy = slab_y + 1;

        int column0 = slab_x * GRID_NP + slab_y;
        int column1 = slab_x * GRID_NP + slab_yy;
        int column2 = slab_xx * GRID_NP + slab_y;
        int column3 = slab_xx * GRID_NP + slab_yy;

        int task0, task1, task2, task3;

        if(column0 < myplan.pivotcol)
          task0 = column0 / myplan.avg;
        else
          task0 = (column0 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column1 < myplan.pivotcol)
          task1 = column1 / myplan.avg;
        else
          task1 = (column1 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column2 < myplan.pivotcol)
          task2 = column2 / myplan.avg;
        else
          task2 = (column2 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column3 < myplan.pivotcol)
          task3 = column3 / myplan.avg;
        else
          task3 = (column3 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        size_t ind0        = send_offset[task0] + send_count[task0]++;
        partout[ind0].Mass = P[i].Mass;
        for(j = 0; j < 3; j++)
          partout[ind0].Pos[j] = pos[j];

        if(task1 != task0)
          {
            size_t ind1        = send_offset[task1] + send_count[task1]++;
            partout[ind1].Mass = P[i].Mass;
            for(j = 0; j < 3; j++)
              partout[ind1].Pos[j] = pos[j];
          }
        if(task2 != task1 && task2 != task0)
          {
            size_t ind2        = send_offset[task2] + send_count[task2]++;
            partout[ind2].Mass = P[i].Mass;
            for(j = 0; j < 3; j++)
              partout[ind2].Pos[j] = pos[j];
          }
        if(task3 != task0 && task3 != task1 && task3 != task2)
          {
            size_t ind3        = send_offset[task3] + send_count[task3]++;
            partout[ind3].Mass = P[i].Mass;
            for(j = 0; j < 3; j++)
              partout[ind3].Pos[j] = pos[j];
          }
#endif
      }
  }

  /* collect the send_count[] results from the other threads */
  for(j = 1; j < MaxThreads; j++)
    for(i = 0; i < NTask; i++)
      Sndpm_count[i] += Sndpm_count[i + j * multiNtask];

  int flag_big = 0, flag_big_all;
  for(i = 0; i < NTask; i++)
    if(Sndpm_count[i] * sizeof(struct partbuf) > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
      flag_big = 1;

  /* produce a flag if any of the send sizes is above our transfer limit, in this case we will
   * transfer the data in chunks.
   */
  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  /* exchange particle data */
  myMPI_Alltoallv(partout, Sndpm_count, Sndpm_offset, partin, Rcvpm_count, Rcvpm_offset, sizeof(struct partbuf), flag_big_all,
                  MPI_COMM_WORLD);

  myfree(partout);

  /* allocate density field */
  rhogrid = (fft_real *)mymalloc("rhogrid", maxfftsize * sizeof(fft_real));

  /* clear local FFT-mesh density field */
  large_array_offset ii;
#pragma omp parallel for
  for(ii = 0; ii < maxfftsize; ii++)
    rhogrid[ii] = 0;

#ifndef FFT_COLUMN_BASED
    /* bin particle data onto mesh, in multi-threaded fashion */
#pragma omp parallel private(i)
  {
    int tid = get_thread_num();

    int first_y, count_y;
    subdivide_evenly(GRID_NP, MaxThreads, tid, &first_y, &count_y);
    int last_y = first_y + count_y - 1;

    for(i = 0; i < nimport; i++)
      {
        int slab_y  = (int)(to_slab_fac * (partin[i].Pos[1] - All.Corner[grnr][1]));
        int slab_yy = slab_y + 1;
        double dy   = to_slab_fac * (partin[i].Pos[1] - All.Corner[grnr][1]) - slab_y;
        int flag_slab_y, flag_slab_yy;

        if(slab_y >= first_y && slab_y <= last_y)
          flag_slab_y = 1;
        else
          flag_slab_y = 0;

        if(slab_yy >= first_y && slab_yy <= last_y)
          flag_slab_yy = 1;
        else
          flag_slab_yy = 0;

        if(flag_slab_y || flag_slab_yy)
          {
            double mass = partin[i].Mass;

            int slab_x  = (int)(to_slab_fac * (partin[i].Pos[0] - All.Corner[grnr][0]));
            int slab_z  = (int)(to_slab_fac * (partin[i].Pos[2] - All.Corner[grnr][2]));
            int slab_xx = slab_x + 1;
            int slab_zz = slab_z + 1;

            double dx = to_slab_fac * (partin[i].Pos[0] - All.Corner[grnr][0]) - slab_x;
            double dz = to_slab_fac * (partin[i].Pos[2] - All.Corner[grnr][2]) - slab_z;

            int flag_slab_x, flag_slab_xx;

            if(myplan.slab_to_task[slab_x] == ThisTask)
              {
                slab_x -= myplan.first_slab_x_of_task[ThisTask];
                flag_slab_x = 1;
              }
            else
              flag_slab_x = 0;

            if(myplan.slab_to_task[slab_xx] == ThisTask)
              {
                slab_xx -= myplan.first_slab_x_of_task[ThisTask];
                flag_slab_xx = 1;
              }
            else
              flag_slab_xx = 0;

            if(flag_slab_x)
              {
                if(flag_slab_y)
                  {
                    rhogrid[FI(slab_x, slab_y, slab_z)] += (mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz));
                    rhogrid[FI(slab_x, slab_y, slab_zz)] += (mass * (1.0 - dx) * (1.0 - dy) * (dz));
                  }

                if(flag_slab_yy)
                  {
                    rhogrid[FI(slab_x, slab_yy, slab_z)] += (mass * (1.0 - dx) * (dy) * (1.0 - dz));
                    rhogrid[FI(slab_x, slab_yy, slab_zz)] += (mass * (1.0 - dx) * (dy) * (dz));
                  }
              }

            if(flag_slab_xx)
              {
                if(flag_slab_y)
                  {
                    rhogrid[FI(slab_xx, slab_y, slab_z)] += (mass * (dx) * (1.0 - dy) * (1.0 - dz));
                    rhogrid[FI(slab_xx, slab_y, slab_zz)] += (mass * (dx) * (1.0 - dy) * (dz));
                  }

                if(flag_slab_yy)
                  {
                    rhogrid[FI(slab_xx, slab_yy, slab_z)] += (mass * (dx) * (dy) * (1.0 - dz));
                    rhogrid[FI(slab_xx, slab_yy, slab_zz)] += (mass * (dx) * (dy) * (dz));
                  }
              }
          }
      }
  }

#else

  struct data_cols
  {
    int col0, col1, col2, col3;
    double dx, dy;
  } * aux;

  aux = (struct data_cols *)mymalloc("aux", nimport * sizeof(struct data_cols));

#pragma omp parallel for
  for(i = 0; i < nimport; i++)
    {
      int slab_x = (int)(to_slab_fac * (partin[i].Pos[0] - All.Corner[grnr][0]));
      int slab_xx = slab_x + 1;

      int slab_y = (int)(to_slab_fac * (partin[i].Pos[1] - All.Corner[grnr][1]));
      int slab_yy = slab_y + 1;

      aux[i].dx = to_slab_fac * (partin[i].Pos[0] - All.Corner[grnr][0]) - slab_x;
      aux[i].dy = to_slab_fac * (partin[i].Pos[1] - All.Corner[grnr][1]) - slab_y;

      aux[i].col0 = slab_x * GRID_NP + slab_y;
      aux[i].col1 = slab_x * GRID_NP + slab_yy;
      aux[i].col2 = slab_xx * GRID_NP + slab_y;
      aux[i].col3 = slab_xx * GRID_NP + slab_yy;
    }

#pragma omp parallel private(i)
  {
    int tid = get_thread_num();

    int first_col, last_col, count_col;
    subdivide_evenly(myplan.base_ncol, MaxThreads, tid, &first_col, &count_col);
    last_col = first_col + count_col - 1;
    first_col += myplan.base_firstcol;
    last_col += myplan.base_firstcol;

    for(i = 0; i < nimport; i++)
      {
        int flag0, flag1, flag2, flag3;
        int col0 = aux[i].col0;
        int col1 = aux[i].col1;
        int col2 = aux[i].col2;
        int col3 = aux[i].col3;

        if(col0 >= first_col && col0 <= last_col)
          flag0 = 1;
        else
          flag0 = 0;

        if(col1 >= first_col && col1 <= last_col)
          flag1 = 1;
        else
          flag1 = 0;

        if(col2 >= first_col && col2 <= last_col)
          flag2 = 1;
        else
          flag2 = 0;

        if(col3 >= first_col && col3 <= last_col)
          flag3 = 1;
        else
          flag3 = 0;

        if(flag0 || flag1 || flag2 || flag3)
          {
            double mass = partin[i].Mass;

            double dx = aux[i].dx;
            double dy = aux[i].dy;

            int slab_z = (int)(to_slab_fac * (partin[i].Pos[2] - All.Corner[grnr][2]));
            int slab_zz = slab_z + 1;

            double dz = to_slab_fac * (partin[i].Pos[2] - All.Corner[grnr][2]) - slab_z;

            const fft_plan *const m = &myplan;
            if(flag0)
              {
                rhogrid[FC(col0, slab_z)] += (mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz));
                rhogrid[FC(col0, slab_zz)] += (mass * (1.0 - dx) * (1.0 - dy) * (dz));
              }

            if(flag1)
              {
                rhogrid[FC(col1, slab_z)] += (mass * (1.0 - dx) * (dy) * (1.0 - dz));
                rhogrid[FC(col1, slab_zz)] += (mass * (1.0 - dx) * (dy) * (dz));
              }

            if(flag2)
              {
                rhogrid[FC(col2, slab_z)] += (mass * (dx) * (1.0 - dy) * (1.0 - dz));
                rhogrid[FC(col2, slab_zz)] += (mass * (dx) * (1.0 - dy) * (dz));
              }

            if(flag3)
              {
                rhogrid[FC(col3, slab_z)] += (mass * (dx) * (dy) * (1.0 - dz));
                rhogrid[FC(col3, slab_zz)] += (mass * (dx) * (dy) * (dz));
              }
          }
      }
  }

  myfree(aux);

#endif
}

/*! \brief If dim < 0, this function reads out the potential, otherwise
 *         Cartesian force components.
 *
 *  \param[in] grnr Grid number (0: base grid, 1: high res grid).
 *  \param[in] dim Dimension of component to be read out (< 0: potential).
 *
 *  \return void
 */
void pmforce_nonperiodic_uniform_optimized_readout_forces_or_potential(int grnr, int dim)
{
#ifdef EVALPOTENTIAL
  /* factor to get potential */
  double fac = All.G / pow(All.TotalMeshSize[grnr], 4) * pow(All.TotalMeshSize[grnr] / GRID_NP, 3);
#endif

  double to_slab_fac = GRID_NP / All.TotalMeshSize[grnr];

  double *flistin  = (double *)mymalloc("flistin", nimport * sizeof(double));
  double *flistout = (double *)mymalloc("flistout", nexport * sizeof(double));

  fft_real *grid;

  if(dim < 0)
    grid = rhogrid;
  else
    grid = forcegrid;

  size_t i;
#pragma omp parallel for
  for(i = 0; i < nimport; i++)
    {
      flistin[i] = 0;

      int slab_x = (int)(to_slab_fac * (partin[i].Pos[0] - All.Corner[grnr][0]));
      int slab_y = (int)(to_slab_fac * (partin[i].Pos[1] - All.Corner[grnr][1]));
      int slab_z = (int)(to_slab_fac * (partin[i].Pos[2] - All.Corner[grnr][2]));

      double dx = to_slab_fac * (partin[i].Pos[0] - All.Corner[grnr][0]) - slab_x;
      double dy = to_slab_fac * (partin[i].Pos[1] - All.Corner[grnr][1]) - slab_y;
      double dz = to_slab_fac * (partin[i].Pos[2] - All.Corner[grnr][2]) - slab_z;

      int slab_xx = slab_x + 1;
      int slab_yy = slab_y + 1;
      int slab_zz = slab_z + 1;

#ifndef FFT_COLUMN_BASED
      if(myplan.slab_to_task[slab_x] == ThisTask)
        {
          slab_x -= myplan.first_slab_x_of_task[ThisTask];

          flistin[i] += +grid[FI(slab_x, slab_y, slab_z)] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz) +
                        grid[FI(slab_x, slab_y, slab_zz)] * (1.0 - dx) * (1.0 - dy) * (dz) +
                        grid[FI(slab_x, slab_yy, slab_z)] * (1.0 - dx) * (dy) * (1.0 - dz) +
                        grid[FI(slab_x, slab_yy, slab_zz)] * (1.0 - dx) * (dy) * (dz);
        }

      if(myplan.slab_to_task[slab_xx] == ThisTask)
        {
          slab_xx -= myplan.first_slab_x_of_task[ThisTask];

          flistin[i] += +grid[FI(slab_xx, slab_y, slab_z)] * (dx) * (1.0 - dy) * (1.0 - dz) +
                        grid[FI(slab_xx, slab_y, slab_zz)] * (dx) * (1.0 - dy) * (dz) +
                        grid[FI(slab_xx, slab_yy, slab_z)] * (dx) * (dy) * (1.0 - dz) +
                        grid[FI(slab_xx, slab_yy, slab_zz)] * (dx) * (dy) * (dz);
        }
#else
      int column0 = slab_x * GRID_NP + slab_y;
      int column1 = slab_x * GRID_NP + slab_yy;
      int column2 = slab_xx * GRID_NP + slab_y;
      int column3 = slab_xx * GRID_NP + slab_yy;

      const fft_plan *const m = &myplan;
      if(column0 >= myplan.base_firstcol && column0 <= myplan.base_lastcol)
        {
          flistin[i] += grid[FC(column0, slab_z)] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz) +
                        grid[FC(column0, slab_zz)] * (1.0 - dx) * (1.0 - dy) * (dz);
        }
      if(column1 >= myplan.base_firstcol && column1 <= myplan.base_lastcol)
        {
          flistin[i] +=
              grid[FC(column1, slab_z)] * (1.0 - dx) * (dy) * (1.0 - dz) + grid[FC(column1, slab_zz)] * (1.0 - dx) * (dy) * (dz);
        }

      if(column2 >= myplan.base_firstcol && column2 <= myplan.base_lastcol)
        {
          flistin[i] +=
              grid[FC(column2, slab_z)] * (dx) * (1.0 - dy) * (1.0 - dz) + grid[FC(column2, slab_zz)] * (dx) * (1.0 - dy) * (dz);
        }

      if(column3 >= myplan.base_firstcol && column3 <= myplan.base_lastcol)
        {
          flistin[i] += grid[FC(column3, slab_z)] * (dx) * (dy) * (1.0 - dz) + grid[FC(column3, slab_zz)] * (dx) * (dy) * (dz);
        }
#endif
    }

  /* exchange the potential component data */
  int flag_big = 0, flag_big_all;
  for(i = 0; i < NTask; i++)
    if(Sndpm_count[i] * sizeof(double) > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
      flag_big = 1;

  /* produce a flag if any of the send sizes is above our transfer limit, in this case we will
   * transfer the data in chunks.
   */
  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  /* exchange  data */
  myMPI_Alltoallv(flistin, Rcvpm_count, Rcvpm_offset, flistout, Sndpm_count, Sndpm_offset, sizeof(double), flag_big_all,
                  MPI_COMM_WORLD);

  /* now assign them to the correct particles */
  int multiNtask = roundup_to_multiple_of_cacheline_size(NTask * sizeof(size_t)) / sizeof(size_t);

#pragma omp parallel /* note: this is *not* a parallel for, instead each threads needs to do the whole loop */
  {
    size_t *send_count  = Sndpm_count + get_thread_num() * multiNtask;
    size_t *send_offset = Sndpm_offset + get_thread_num() * multiNtask;

    int j;
    for(j = 0; j < NTask; j++)
      send_count[j] = 0;

    int i;
#pragma omp for schedule(static)
    for(i = 0; i < NumPart; i++)
      {
        MyDouble *pos;

#ifdef CELL_CENTER_GRAVITY
        if(P[i].Type == 0)
          pos = SphP[i].Center;
        else
#endif
          pos = P[i].Pos;

        if(pos[0] < All.Corner[grnr][0] || pos[0] >= All.UpperCorner[grnr][0])
          continue;
        if(pos[1] < All.Corner[grnr][1] || pos[1] >= All.UpperCorner[grnr][1])
          continue;
        if(pos[2] < All.Corner[grnr][2] || pos[2] >= All.UpperCorner[grnr][2])
          continue;

        int slab_x  = (int)(to_slab_fac * (pos[0] - All.Corner[grnr][0]));
        int slab_xx = slab_x + 1;

#ifndef FFT_COLUMN_BASED
        int task0   = myplan.slab_to_task[slab_x];
        int task1   = myplan.slab_to_task[slab_xx];

        double value = flistout[send_offset[task0] + send_count[task0]++];

        if(task0 != task1)
          value += flistout[send_offset[task1] + send_count[task1]++];
#else
        int slab_y = (int)(to_slab_fac * (pos[1] - All.Corner[grnr][1]));
        int slab_yy = slab_y + 1;

        int column0 = slab_x * GRID_NP + slab_y;
        int column1 = slab_x * GRID_NP + slab_yy;
        int column2 = slab_xx * GRID_NP + slab_y;
        int column3 = slab_xx * GRID_NP + slab_yy;

        int task0, task1, task2, task3;

        if(column0 < myplan.pivotcol)
          task0 = column0 / myplan.avg;
        else
          task0 = (column0 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column1 < myplan.pivotcol)
          task1 = column1 / myplan.avg;
        else
          task1 = (column1 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column2 < myplan.pivotcol)
          task2 = column2 / myplan.avg;
        else
          task2 = (column2 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column3 < myplan.pivotcol)
          task3 = column3 / myplan.avg;
        else
          task3 = (column3 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        double value = flistout[send_offset[task0] + send_count[task0]++];

        if(task1 != task0)
          value += flistout[send_offset[task1] + send_count[task1]++];

        if(task2 != task1 && task2 != task0)
          value += flistout[send_offset[task2] + send_count[task2]++];

        if(task3 != task0 && task3 != task1 && task3 != task2)
          value += flistout[send_offset[task3] + send_count[task3]++];
#endif

#ifdef PLACEHIGHRESREGION
        if(grnr == 1)
          if(!(pmforce_is_particle_high_res(P[i].Type, pos)))
            continue;
#endif

        if(dim < 0)
          {
#ifdef EVALPOTENTIAL
            P[i].PM_Potential += value * fac;
#endif
          }
        else
          P[i].GravPM[dim] += value;
      }
  }

  int j;
  /* restore total Sndpm_count */
  for(j = 1; j < MaxThreads; j++)
    for(i = 0; i < NTask; i++)
      Sndpm_count[i] += Sndpm_count[i + j * multiNtask];

  myfree(flistout);
  myfree(flistin);
}
#endif

/*! \brief Calculates the long-range non-periodic forces using the PM method.
 *
 *  The potential is Gaussian filtered with Asmth, given in mesh-cell units.
 *  The potential is finite differenced using a 4-point finite differencing
 *  formula to obtain the force fields, which are then interpolated to the
 *  particle positions. We carry out a CIC charge assignment, and compute the
 *  potenial by Fourier transform methods. The CIC kernel is deconvolved.
 *
 *  \param[in] grnr Grid number (0: base grid, 1 high res grid).
 *
 *  \return 0
 */
int pmforce_nonperiodic(int grnr)
{
  int i, j, flag, flagsum, dim;

  double tstart = second();

  mpi_printf("PM-NONPERIODIC: Starting non-periodic PM calculation (grid=%d)  presently allocated=%g MB).\n", grnr,
             AllocatedBytes / (1024.0 * 1024.0));

#ifndef NUMPART_PER_TASK_LARGE
  if(((long long)NumPart) << 3 >= ((long long)1) << 31)
    terminate("We are dealing with a too large particle number per MPI rank - enabling NUMPART_PER_TASK_LARGE might help.");
#endif

  double fac = All.G / pow(All.TotalMeshSize[grnr], 4) * pow(All.TotalMeshSize[grnr] / GRID_NP, 3); /* to get potential */
  fac *= 1 / (2 * All.TotalMeshSize[grnr] / GRID_NP);                                               /* for finite differencing */

  /* first, check whether all particles lie in the allowed region */
  for(i = 0, flag = 0; i < NumPart; i++)
    {
      MyDouble *pos;

#ifdef CELL_CENTER_GRAVITY
      if(P[i].Type == 0)
        pos = SphP[i].Center;
      else
#endif
        pos = P[i].Pos;

#ifdef PLACEHIGHRESREGION
      if(grnr == 0 || (grnr == 1 && pmforce_is_particle_high_res(P[i].Type, pos)))
#endif
        {
          for(j = 0; j < 3; j++)
            {
              if(pos[j] < All.Xmintot[grnr][j] || pos[j] > All.Xmaxtot[grnr][j])
                {
                  if(flag == 0)
                    {
                      printf("Particle Id=%llu on task=%d with coordinates (%g|%g|%g) lies outside PM mesh.\n",
                             (unsigned long long)P[i].ID, ThisTask, pos[0], pos[1], pos[2]);
                      myflush(stdout);
                    }
                  flag++;
                  break;
                }
            }
        }
    }

  MPI_Allreduce(&flag, &flagsum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(flagsum > 0)
    {
      mpi_printf("PM-NONPERIODIC: In total %d particles were outside allowed range.\n", flagsum);
      return 1; /* error - need to return because particles were outside allowed range */
    }

#ifdef PM_ZOOM_OPTIMIZED
  pmforce_nonperiodic_zoom_optimized_prepare_density(grnr);
#else
  pmforce_nonperiodic_uniform_optimized_prepare_density(grnr);
#endif

  /* allocate the memory to hold the FFT fields */
  forcegrid = (fft_real *)mymalloc("forcegrid", maxfftsize * sizeof(fft_real));

  workspace = forcegrid;

#ifndef FFT_COLUMN_BASED
  fft_of_rhogrid = (fft_complex *)&rhogrid[0];
#else
  fft_of_rhogrid = (fft_complex *)&workspace[0];
#endif

  /* Do the FFT of the density field */
#ifndef FFT_COLUMN_BASED
  my_slab_based_fft(&myplan, &rhogrid[0], &workspace[0], 1);
#else
  my_column_based_fft(&myplan, rhogrid, workspace, 1); /* result is in workspace, not in rhogrid ! */
#endif

  /* multiply with kernel in Fourier space */
  /* multiply with the Fourier transform of the Green's function (kernel) */

  /* multiply with Green's function in order to obtain the potential */

#ifdef FFT_COLUMN_BASED
#pragma omp parallel for private(x, y, z)
  for(large_array_offset ip = 0; ip < myplan.second_transposed_ncells; ip++)
    {
#else
#pragma omp parallel for private(x, y, z)
  for(int x = 0; x < GRID_NP; x++)
    for(int y = myplan.slabstart_y; y < myplan.slabstart_y + myplan.nslab_y; y++)
      for(int z = 0; z < GRIDz_NP; z++)
        {
#endif

#ifndef FFT_COLUMN_BASED
      large_array_offset ip = ((large_array_offset)GRIDz_NP) * (GRID_NP * (y - myplan.slabstart_y) + x) + z;
#endif

      double re = fft_of_rhogrid[ip][0] * fft_of_kernel[grnr][ip][0] - fft_of_rhogrid[ip][1] * fft_of_kernel[grnr][ip][1];
      double im = fft_of_rhogrid[ip][0] * fft_of_kernel[grnr][ip][1] + fft_of_rhogrid[ip][1] * fft_of_kernel[grnr][ip][0];

      fft_of_rhogrid[ip][0] = re;
      fft_of_rhogrid[ip][1] = im;
    }

    /* Do the inverse FFT to get the potential */

#ifndef FFT_COLUMN_BASED
  my_slab_based_fft(&myplan, rhogrid, workspace, -1);
#else
  my_column_based_fft(&myplan, workspace, rhogrid, -1);
#endif

  /* Now rhogrid holds the potential */

#ifdef EVALPOTENTIAL
#ifdef PM_ZOOM_OPTIMIZED
  pmforce_nonperiodic_zoom_optimized_readout_forces_or_potential(grnr, -1);
#else
  pmforce_nonperiodic_uniform_optimized_readout_forces_or_potential(grnr, -1);
#endif
#endif

  /* get the force components by finite differencing of the potential for each dimension,
   * and send the results back to the right CPUs
   */
  for(dim = 2; dim >= 0; dim--) /* Calculate each component of the force. */
    {
      /* we do the x component last, because for differencing the potential in the x-direction, we need to construct the transpose */
#ifndef FFT_COLUMN_BASED
      if(dim == 0)
        my_slab_transposeA(&myplan, rhogrid, forcegrid); /* compute the transpose of the potential field for finite differencing */

#pragma omp parallel for private(x, z)
      for(int y = 2; y < GRID_NP / 2 - 2; y++)
        for(int x = 0; x < myplan.nslab_x; x++)
          if(x + myplan.slabstart_x >= 2 && x + myplan.slabstart_x < GRID_NP / 2 - 2)
            for(int z = 2; z < GRID_NP / 2 - 2; z++)
              {
                int yrr = y, yll = y, yr = y, yl = y;
                int zrr = z, zll = z, zr = z, zl = z;

                switch(dim)
                  {
                    case 0: /* note: for the x-direction, we difference the transposed direction (y) */
                    case 1:
                      yr  = y + 1;
                      yl  = y - 1;
                      yrr = y + 2;
                      yll = y - 2;

                      break;
                    case 2:
                      zr  = z + 1;
                      zl  = z - 1;
                      zrr = z + 2;
                      zll = z - 2;

                      break;
                  }

                if(dim == 0)
                  forcegrid[TI(x, y, z)] = fac * ((4.0 / 3) * (rhogrid[TI(x, yl, zl)] - rhogrid[TI(x, yr, zr)]) -
                                                  (1.0 / 6) * (rhogrid[TI(x, yll, zll)] - rhogrid[TI(x, yrr, zrr)]));
                else
                  forcegrid[FI(x, y, z)] = fac * ((4.0 / 3) * (rhogrid[FI(x, yl, zl)] - rhogrid[FI(x, yr, zr)]) -
                                                  (1.0 / 6) * (rhogrid[FI(x, yll, zll)] - rhogrid[FI(x, yrr, zrr)]));
              }

      if(dim == 0)
        my_slab_transposeB(&myplan, forcegrid, rhogrid); /* reverse the transpose from above */
#else
      fft_real *scratch = NULL, *forcep, *potp;

      if(dim != 2)
        {
          scratch = (fft_real *)mymalloc("scratch", myplan.fftsize * sizeof(fft_real)); /* need a third field as scratch space */
          memcpy(scratch, rhogrid, myplan.fftsize * sizeof(fft_real));

          if(dim == 1)
            my_fft_swap23(&myplan, scratch, forcegrid);
          else
            my_fft_swap13(&myplan, scratch, forcegrid);
        }

      int ncols;
      if(dim == 2)
        ncols = myplan.base_ncol;
      else if(dim == 1)
        ncols = myplan.ncol_XZ;
      else
        ncols = myplan.ncol_YZ;

      large_array_offset i;

#pragma omp parallel for private(forcep, potp)
      for(i = 0; i < ncols; i++)
        {
          if(dim != 2)
            {
              forcep = &scratch[GRID_NP * i];
              potp   = &forcegrid[GRID_NP * i];
            }
          else
            {
              forcep = &forcegrid[GRID2_NP * i];
              potp   = &rhogrid[GRID2_NP * i];
            }

          int z;
          for(z = 2; z < GRID_NP / 2 - 2; z++)
            {
              int zr  = z + 1;
              int zl  = z - 1;
              int zrr = z + 2;
              int zll = z - 2;

              forcep[z] = fac * ((4.0 / 3) * (potp[zl] - potp[zr]) - (1.0 / 6) * (potp[zll] - potp[zrr]));
            }
        }

      if(dim != 2)
        {
          if(dim == 1)
            my_fft_swap23back(&myplan, scratch, forcegrid);
          else
            my_fft_swap13back(&myplan, scratch, forcegrid);

          myfree(scratch);
        }
#endif

#ifdef PM_ZOOM_OPTIMIZED
      pmforce_nonperiodic_zoom_optimized_readout_forces_or_potential(grnr, dim);
#else
      pmforce_nonperiodic_uniform_optimized_readout_forces_or_potential(grnr, dim);
#endif
    }

  /* free stuff */
  myfree(forcegrid);
  myfree(rhogrid);

#ifdef PM_ZOOM_OPTIMIZED
  myfree(localfield_recvcount);
  myfree(localfield_offset);
  myfree(localfield_sendcount);
  myfree(localfield_first);
  myfree(localfield_data);
  myfree(localfield_globalindex);
  myfree(part);
#else
  myfree(partin);
  myfree(Rcvpm_offset);
  myfree(Rcvpm_count);
  myfree(Sndpm_offset);
  myfree(Sndpm_count);
#endif

  double tend = second();

  mpi_printf("PM-NONPERIODIC: done.  (took %g seconds)\n", timediff(tstart, tend));

  return 0;
}

/*! \brief Sets-up the Greens function for the non-periodic potential in real
 *         space, and then converts it to Fourier space by means of an FFT.
 *
 *  \return void
 */
void pm_setup_nonperiodic_kernel(void)
{
  int i, j, k, x, y, z;
  double xx, yy, zz, r, u, fac;

  mpi_printf("PM-NONPERIODIC: Setting up non-periodic PM kernel (GRID_NP=%d)  presently allocated=%g MB).\n", (int)GRID_NP,
             AllocatedBytes / (1024.0 * 1024.0));

  /* now set up kernel and its Fourier transform */

#ifdef GRAVITY_NOT_PERIODIC
  for(i = 0; i < maxfftsize; i++) /* clear local field */
    kernel[0][i] = 0;

#ifndef FFT_COLUMN_BASED
  for(i = myplan.slabstart_x; i < (myplan.slabstart_x + myplan.nslab_x); i++)
    for(j = 0; j < GRID_NP; j++)
      {
#else
  int c;
  for(c = myplan.base_firstcol; c < (myplan.base_firstcol + myplan.base_ncol); c++)
    {
      i = c / GRID_NP;
      j = c % GRID_NP;
#endif
        for(k = 0; k < GRID_NP; k++)
          {
            xx = ((double)i) / GRID_NP;
            yy = ((double)j) / GRID_NP;
            zz = ((double)k) / GRID_NP;

            if(xx >= 0.5)
              xx -= 1.0;
            if(yy >= 0.5)
              yy -= 1.0;
            if(zz >= 0.5)
              zz -= 1.0;

            r = sqrt(xx * xx + yy * yy + zz * zz);

            u = 0.5 * r / (((double)ASMTH) / GRID_NP);

            fac = 1 - erfc(u);

#ifndef FFT_COLUMN_BASED
            size_t ip = FI(i - myplan.slabstart_x, j, k);
#else
          size_t ip = FC(c, k);
#endif
            if(r > 0)
              kernel[0][ip] = -fac / r;
            else
              kernel[0][ip] = -1 / (sqrt(M_PI) * (((double)ASMTH) / GRID_NP));
          }
      }

  {
    fft_real *workspc = (fft_real *)mymalloc("workspc", maxfftsize * sizeof(fft_real));
    /* Do the FFT of the kernel */
#ifndef FFT_COLUMN_BASED
    my_slab_based_fft(&myplan, kernel[0], workspc, 1);
#else
    my_column_based_fft(&myplan, kernel[0], workspc, 1); /* result is in workspace, not in kernel */
    memcpy(kernel[0], workspc, maxfftsize * sizeof(fft_real));
#endif
    myfree(workspc);
  }

#endif

#ifdef PLACEHIGHRESREGION

  for(size_t i = 0; i < maxfftsize; i++) /* clear local field */
    kernel[1][i] = 0;

#ifndef FFT_COLUMN_BASED
  for(i = myplan.slabstart_x; i < (myplan.slabstart_x + myplan.nslab_x); i++)
    for(j = 0; j < GRID_NP; j++)
      {
#else
  int c;
  for(c = myplan.base_firstcol; c < (myplan.base_firstcol + myplan.base_ncol); c++)
    {
      i = c / GRID_NP;
      j = c % GRID_NP;
#endif
        for(k = 0; k < GRID_NP; k++)
          {
            xx = ((double)i) / GRID_NP;
            yy = ((double)j) / GRID_NP;
            zz = ((double)k) / GRID_NP;

            if(xx >= 0.5)
              xx -= 1.0;
            if(yy >= 0.5)
              yy -= 1.0;
            if(zz >= 0.5)
              zz -= 1.0;

            r = sqrt(xx * xx + yy * yy + zz * zz);

            u = 0.5 * r / (((double)ASMTH) / GRID_NP);

            fac = erfc(u * All.Asmth[1] / All.Asmth[0]) - erfc(u);

#ifndef FFT_COLUMN_BASED
            size_t ip = FI(i - myplan.slabstart_x, j, k);
#else
          size_t ip = FC(c, k);
#endif

            if(r > 0)
              kernel[1][ip] = -fac / r;
            else
              {
                fac           = 1 - All.Asmth[1] / All.Asmth[0];
                kernel[1][ip] = -fac / (sqrt(M_PI) * (((double)ASMTH) / GRID_NP));
              }
          }
#ifndef FFT_COLUMN_BASED
      }
#else
    }
#endif

  {
    fft_real *workspc = (fft_real *)mymalloc("workspc", maxfftsize * sizeof(fft_real));
    /* Do the FFT of the kernel */
#ifndef FFT_COLUMN_BASED
    my_slab_based_fft(&myplan, kernel[1], workspc, 1);
#else
    my_column_based_fft(&myplan, kernel[1], workspc, 1); /* result is in workspace, not in kernel */
    memcpy(kernel[1], workspc, maxfftsize * sizeof(fft_real));
#endif
    myfree(workspc);
  }

#endif

  /* deconvolve the Greens function twice with the CIC kernel */
#ifdef FFT_COLUMN_BASED

  large_array_offset ip, ipcell;

#pragma omp parallel for private(x, y, z)
  for(ip = 0; ip < myplan.second_transposed_ncells; ip++)
    {
      ipcell = ip + myplan.transposed_firstcol * GRID_NP;
      y      = ipcell / (GRID_NP * GRIDz_NP);
      int yr = ipcell % (GRID_NP * GRIDz_NP);
      z      = yr / GRID_NP;
      x      = yr % GRID_NP;
#else
#pragma omp parallel for private(y, z)
  for(x = 0; x < GRID_NP; x++)
    for(y = myplan.slabstart_y; y < myplan.slabstart_y + myplan.nslab_y; y++)
      for(z = 0; z < GRIDz_NP; z++)
        {
#endif

      double kx, ky, kz;

      if(x > GRID_NP / 2)
        kx = x - GRID_NP;
      else
        kx = x;
      if(y > GRID_NP / 2)
        ky = y - GRID_NP;
      else
        ky = y;
      if(z > GRID_NP / 2)
        kz = z - GRID_NP;
      else
        kz = z;

      double k2 = kx * kx + ky * ky + kz * kz;

      if(k2 > 0)
        {
          double fx = 1, fy = 1, fz = 1;

          if(kx != 0)
            {
              fx = (M_PI * kx) / GRID_NP;
              fx = sin(fx) / fx;
            }
          if(ky != 0)
            {
              fy = (M_PI * ky) / GRID_NP;
              fy = sin(fy) / fy;
            }
          if(kz != 0)
            {
              fz = (M_PI * kz) / GRID_NP;
              fz = sin(fz) / fz;
            }

          double ff = 1 / (fx * fy * fz);
          ff        = ff * ff * ff * ff;

#ifndef FFT_COLUMN_BASED
          large_array_offset ip = ((large_array_offset)GRIDz_NP) * (GRID_NP * (y - myplan.slabstart_y) + x) + z;
#endif
#if defined(GRAVITY_NOT_PERIODIC)
          fft_of_kernel[0][ip][0] *= ff;
          fft_of_kernel[0][ip][1] *= ff;
#endif
#if defined(PLACEHIGHRESREGION)
          fft_of_kernel[1][ip][0] *= ff;
          fft_of_kernel[1][ip][1] *= ff;
#endif
        }
    }

  /* end deconvolution */
}

#ifdef PM_ZOOM_OPTIMIZED

/*! \brief Sort function for #part array indices.
 *
 * Sorts the indices into the #part array by the global index of the corresponding #part_slab_data struct.
 *
 * \param a index to be compared
 * \param b index to be compared
 * \return sort result
 */
static int pm_periodic_compare_sortindex(const void *a, const void *b)
{
  if(part[*(int *)a].globalindex < part[*(int *)b].globalindex)
    return -1;

  if(part[*(int *)a].globalindex > part[*(int *)b].globalindex)
    return +1;

  return 0;
}

/*! \brief Implements the sorting function for mysort_pmperiodic()
 *
 *  The index array is sorted using a merge sort algorithm.
 *
 *  \param[in, out] b Index array to sort.
 *  \param[in] n Number of elements to sort.
 *  \param[out] t Temporary buffer array.
 *
 *  \return void
 */
static void msort_pmperiodic_with_tmp(large_numpart_type *b, size_t n, large_numpart_type *t)
{
  large_numpart_type *tmp;
  large_numpart_type *b1, *b2;
  size_t n1, n2;

  if(n <= 1)
    return;

  n1 = n / 2;
  n2 = n - n1;
  b1 = b;
  b2 = b + n1;

  msort_pmperiodic_with_tmp(b1, n1, t);
  msort_pmperiodic_with_tmp(b2, n2, t);

  tmp = t;

  while(n1 > 0 && n2 > 0)
    {
      if(part[*b1].globalindex <= part[*b2].globalindex)
        {
          --n1;
          *tmp++ = *b1++;
        }
      else
        {
          --n2;
          *tmp++ = *b2++;
        }
    }

  if(n1 > 0)
    memcpy(tmp, b1, n1 * sizeof(large_numpart_type));

  memcpy(b, t, (n - n2) * sizeof(large_numpart_type));
}

/*! \brief Sort the index array b of n entries using the sort kernel
 * cmp.
 *
 * The parameter s is set to sizeof(int). The index array b
 * is sorted according to the globalindex field of the referenced item in the #part array
 *
 * \param b the index array to sort
 * \param n number of entries in array b
 * \param size of each entry (must be sizeof(int))
 * \param cmp comparator function
 */
static void mysort_pmperiodic(void *b, size_t n, size_t s, int (*cmp)(const void *, const void *))
{
  const size_t size = n * s;

  large_numpart_type *tmp = (large_numpart_type *)mymalloc("tmp", size);

  msort_pmperiodic_with_tmp((large_numpart_type *)b, n, tmp);

  myfree(tmp);
}
#endif

#endif
