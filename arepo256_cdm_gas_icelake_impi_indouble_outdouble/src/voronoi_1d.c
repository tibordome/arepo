/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/voronoi_1d.c
 * \date        MM/YYYY
 * \author
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"

#if defined(VORONOI) && defined(ONEDIMS) && !defined(ONEDIMS_SPHERICAL) /* will only be compiled in 1D case */

/*! \brief Output of Voroioi mesh to file.
 *
 *  Not supported for 1d.
 *
 *  \return void
 */
void write_voronoi_mesh(tessellation *T, char *fname, int writeTask, int lastTask)
{
  terminate("You stupid idiot should not do that!");
}

/*! \brief Initialises 1d tessellation and create all-enclosing segment.
 *
 *  \param[out] T Pointer to tessllation structure which is set and its arrays
 *              are allocated in this routine.
 *
 *  \return void
 */
void initialize_and_create_first_tetra(tessellation *T)
{
  char msg[200];

  if(NTask > 1)
    mpi_terminate("1D code works only for 1 CPU\n");

  T->MaxNdp = NumGas + 4;
  T->MaxNdt = 4 + T->MaxNdp * 2;
  T->MaxNvf = T->MaxNdt;

  if(NumGas == 0)
    {
      sprintf(msg, "NumGas=%d on Task=%d, but need at least one particle!\n", NumGas, ThisTask);
      terminate(msg);
    }

  T->Ndp = 0;
  T->Nvf = 0;
  T->Ndt = 0;

  T->VF = (face *)mymalloc_movable(&T->VF, "VF", T->MaxNvf * sizeof(face));

  T->DP = (point *)mymalloc_movable(&T->DP, "DP", (T->MaxNdp + 5) * sizeof(point));
  T->DP += 5;

  T->DT = (tetra *)mymalloc_movable(&T->DT, "DT", T->MaxNdt * sizeof(tetra));
}

/*! \brief Computes circumcircles in 1d.
 *
 *  Not necessary in 1d. However, this function has to exist for the 1d code
 *  to work.
 *
 *  \param[in] T Pointer to tessllation structure.
 *
 *  \return void
 */
void compute_circumcircles(tessellation *T) {}

/*! \brief Empty funciton in 1d case.
 *
 *  Not necessary in 1d. However, this function has to exist for the 1d code
 *  to work.
 *
 * \return void
 */
void set_integers_for_point(tessellation *T, int pp) {}

int insert_point(tessellation *T, int pp, int ttstart) /* returns a triangle that (currently) contains the point p */ { return 0; }

/*! \brief Wrapper routine to search for ghost cells for boundary cells.
 *
 *  \param[out] T Pointer to tessellation.
 *
 *  \return 0
 */
int voronoi_ghost_search(tessellation *T) { return voronoi_ghost_search_alternative(T); }

/*! \brief Empty funciton in 1d case.
 *
 *  Not necessary in 1d. However, this function has to exist for the 1d code
 *  to work.
 *
 * \return 0
 */
int count_undecided_tetras(tessellation *T) { return 0; }

/*! \brief Searches for ghost cells in 1d Voronoi mesh.
 *
 *  This routine assumes an x ordered cell array.
 *
 *  \param[out] T pointer to tessellation.
 *
 *  \return 0
 */
int voronoi_ghost_search_alternative(tessellation *T)
{
  double xl, xr;
  int index_l, index_r;

#ifdef REFLECTIVE_X
  xl      = -P[0].Pos[0];
  index_l = 0;

  xr      = boxSize_X + (boxSize_X - P[NumGas - 1].Pos[0]);
  index_r = NumGas - 1;
#else
  xl      = P[NumGas - 1].Pos[0] - boxSize_X;
  index_l = NumGas - 1;

  xr      = P[0].Pos[0] + boxSize_X;
  index_r = 0;
#endif

  point *DP = T->DP;

  DP[-1].x     = xl;
  DP[-1].y     = 0;
  DP[-1].z     = 0;
  DP[-1].task  = ThisTask;
  DP[-1].ID    = P[index_l].ID;
  DP[-1].index = index_l + NumGas; /* this is a mirrored local point */
#ifdef REFLECTIVE_X
  DP[-1].image_flags = REFL_X_FLAGS;
#if REFLECTIVE_X == 2
  DP[-1].image_flags |= OUTFLOW_X;
#endif
#endif
  DP[NumGas].x     = xr;
  DP[NumGas].y     = 0;
  DP[NumGas].z     = 0;
  DP[NumGas].task  = ThisTask;
  DP[NumGas].ID    = P[index_r].ID;
  DP[NumGas].index = index_r + NumGas; /* this is a mirrored local point */
#ifdef REFLECTIVE_X
  DP[NumGas].image_flags = REFL_X_FLAGS;
#if REFLECTIVE_X == 2
  DP[NumGas].image_flags |= OUTFLOW_X;
#endif
#endif
  return 0;
}

/*! \brief Computes faces and volume of cells in 1d Voronoi mesh.
 *
 *  Also computes the center of mass.
 *
 *  \return void
 */
void compute_voronoi_faces_and_volumes(void)
{
  int i;

  tessellation *T = &Mesh;

  T->Nvf    = NumGas + 1;
  point *DP = T->DP;
  face *VF  = T->VF;

#pragma omp parallel for private(i)
  for(i = -1; i < NumGas; i++)
    {
      VF[i + 1].p1 = i;
      VF[i + 1].p2 = i + 1;

      VF[i + 1].cx = 0.5 * (DP[i].x + DP[i + 1].x);

      VF[i + 1].cy   = 0;
      VF[i + 1].cz   = 0;
      VF[i + 1].area = 1;
    }

#pragma omp parallel for private(i)
  for(i = 0; i < NumGas; i++)
    {
      SphP[i].Volume    = VF[i + 1].cx - VF[i].cx;
      SphP[i].Center[0] = 0.5 * (VF[i + 1].cx + VF[i].cx);
      SphP[i].Center[1] = 0;
      SphP[i].Center[2] = 0;

      SphP[i].SurfaceArea = 2.;

      /*
         printf("i=%d x=%g volume=%g l=%g r=%g\n", i, P[i].Pos[0], SphP[i].Volume, VF[i].cx, VF[i + 1].cx);
       */
    }
}

static struct voronoi_1D_data
{
  double x;
  int index;
} * mp;

static int *Id;

/*! \brief Sort cells by their position and reorder in P and SphP array.
 *
 *  \return void
 */
void voronoi_1D_order(void)
{
  int i;

  mpi_printf("begin 1D order...\n");

  if(NumGas)
    {
      mp = (struct voronoi_1D_data *)mymalloc("mp", sizeof(struct voronoi_1D_data) * NumGas);
      Id = (int *)mymalloc("Id", sizeof(int) * NumGas);

#pragma omp parallel for private(i)
      for(i = 0; i < NumGas; i++)
        {
          mp[i].index = i;
          mp[i].x     = P[i].Pos[0];
        }

      mysort(mp, NumGas, sizeof(struct voronoi_1D_data), voronoi_1D_compare_key);

#pragma omp parallel for private(i)
      for(i = 0; i < NumGas; i++)
        Id[mp[i].index] = i;

      voronoi_1D_reorder_gas();

      myfree(Id);
      myfree(mp);
    }

  mpi_printf("1D order done.\n");
}

/*! \brief Compare x value of voronoi_1D_data objects.
 *
 *  \param[in] a Pointer to first voronoi_1D_data object.
 *  \param[in] b Pointer to second voronoi_1D_data object.
 *
 *  \return (-1,0,1) -1 if a->x < b->x.
 */
int voronoi_1D_compare_key(const void *a, const void *b)
{
  if(((struct voronoi_1D_data *)a)->x < (((struct voronoi_1D_data *)b)->x))
    return -1;

  if(((struct voronoi_1D_data *)a)->x > (((struct voronoi_1D_data *)b)->x))
    return +1;

  return 0;
}

/*! \brief Order the gas cells according to the index given in the ID array.
 *
 *  \return void
 */
void voronoi_1D_reorder_gas(void)
{
  int i;
  struct particle_data Psave, Psource;
  struct sph_particle_data SphPsave, SphPsource;
  int idsource, idsave, dest;

#ifdef CHIMES
  double *ChimesAbunSave, *ChimesAbunSource;
  ChimesAbunSave =
      (double *)mymalloc_movable(&ChimesAbunSave, "ChimesAbunSave", ChimesGlobalVars.totalNumberOfSpecies * sizeof(double));
  ChimesAbunSource =
      (double *)mymalloc_movable(&ChimesAbunSource, "ChimesAbunSource", ChimesGlobalVars.totalNumberOfSpecies * sizeof(double));

  chimes_update_all_pointers();
#endif

  for(i = 0; i < NumGas; i++)
    {
      if(Id[i] != i)
        {
          Psource    = P[i];
          SphPsource = SphP[i];

#ifdef CHIMES
          memcpy(ChimesAbunSource, SphP[i].ChimesGasVars.abundances, ChimesGlobalVars.totalNumberOfSpecies * sizeof(double));
#endif

          idsource = Id[i];
          dest     = Id[i];

          do
            {
              Psave    = P[dest];
              SphPsave = SphP[dest];
              idsave   = Id[dest];

#ifdef CHIMES
              memcpy(ChimesAbunSave, SphP[dest].ChimesGasVars.abundances, ChimesGlobalVars.totalNumberOfSpecies * sizeof(double));
#endif

              P[dest]    = Psource;
              SphP[dest] = SphPsource;
              Id[dest]   = idsource;

#ifdef CHIMES
              SphP[dest].ChimesGasVars.index = dest;
              chimes_set_pointers(dest);
              memcpy(SphP[dest].ChimesGasVars.abundances, ChimesAbunSource, ChimesGlobalVars.totalNumberOfSpecies * sizeof(double));
#endif

              if(dest == i)
                break;

              Psource    = Psave;
              SphPsource = SphPsave;
              idsource   = idsave;

#ifdef CHIMES
              memcpy(ChimesAbunSource, ChimesAbunSave, ChimesGlobalVars.totalNumberOfSpecies * sizeof(double));
#endif

              dest = idsource;
            }
          while(1);
        }
    }
#ifdef CHIMES
  myfree_movable(ChimesAbunSource);
  myfree_movable(ChimesAbunSave);
#endif
}

#endif
