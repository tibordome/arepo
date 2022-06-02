/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/ngbtree_search.c
 * \date        MM/YYYY
 * \author
 * \brief       searches cells on the neighbor tree
 * \details     This file contains a search routine on the neighbor tree
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
#include <time.h>

#include "allvars.h"
#include "proto.h"

/* temporary particle arrays */
static MyDouble *ngbsearch_nearest_dist;
static MyDouble *ngbsearch_hsml;
static mesh_search_data *searchdata;

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble pos[3];   /* tracer particle position */
  MyDouble hsml;     /* current search radius */
  MyDouble distance; /* nearest neighbor distance */

  int Firstnode;
} data_in;

static data_in *DataGet;

/*! \brief Routine that fills the relevant particle/cell data into the input
 *         structure defined above. Needed by generic_comm_helpers2.
 *
 *  \param[out] in Data structure to fill.
 *  \param[in] i Index of particle in P and SphP arrays.
 *  \param[in] firstnode First note of communication.
 *
 *  \return void
 */
static void particle2in(data_in *in, int i, int firstnode)
{
  in->pos[0] = searchdata[i].Pos[0];
  in->pos[1] = searchdata[i].Pos[1];
  in->pos[2] = searchdata[i].Pos[2];

  in->hsml     = ngbsearch_hsml[i];
  in->distance = ngbsearch_nearest_dist[i];

  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
  MyDouble Distance; /* distance to closest cell on task */
  int Task;
  int Index;
} data_out;

static data_out *DataResult;

/*! \brief Routine to store or combine result data. Needed by
 *         generic_comm_helpers2.
 *
 *  \param[in] out Data to be moved to appropriate variables in global
 *  particle and cell data arrays (P, SphP,...)
 *  \param[in] i Index of particle in P and SphP arrays
 *  \param[in] mode Mode of function: local particles or information that was
 *  communicated from other tasks and has to be added locally?
 *
 *  \return void
 */
static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
      if(out->Index >= 0)
        {
          ngbsearch_nearest_dist[i] = out->Distance;
          searchdata[i].Task        = out->Task;
          searchdata[i].u.Index     = out->Index;
        }
    }
  else /* combine */
    {
      /* closer cell on other task? */
      if(out->Distance < ngbsearch_nearest_dist[i])
        {
          ngbsearch_nearest_dist[i] = out->Distance;
          searchdata[i].Task        = out->Task;
          searchdata[i].u.Index     = out->Index;
        }
    }
}

#define OMIT_GENERIC_COMM_PATTERN_FOR_GIVEN_PARTICLES
#include "generic_comm_helpers2.h"

static int ngbsearch_primary_cell_evaluate(int target, int mode, int threadid);
static int n;

/*! \brief Routine that defines what to do with local particles.
 *
 *  Calls the *_evaluate function in MODE_LOCAL_PARTICLES.
 *
 *  \return void
 */
static void kernel_local(void)
{
  int i;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif

  /* do local particles */
#pragma omp parallel private(i)
  {
    int j, threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              if(generic_polling_primary(count, n))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= n)
          break;

        if(searchdata[i].Task == -1)
          ngbsearch_primary_cell_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}

/*! \brief Routine that defines what to do with imported particles.
 *
 *  Calls the *_evaluate function in MODE_IMPORTED_PARTICLES.
 *
 *  \return void
 */
static void kernel_imported(void)
{
  int i, cnt = 0;
#pragma omp parallel private(i)
  {
    int threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    while(1)
      {
#pragma omp atomic capture
        i = cnt++;

        if(i >= Nimport)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              generic_polling_secondary();
          }

        count++;
#endif

        ngbsearch_primary_cell_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

/*! \brief Searches the cells at the positions in searchdata.
 *
 *  This function searches the cells which are at the positions specified in
 *  searchdata. The Pos field must be set. After the search is performed the
 *  Task and Index field contain the task/index of the cell at position Pos.
 *  If hsmlguess=1 initial search radius is read from Index/Hsml union in
 *  searchdata.
 *
 *  \param[in] searchdata_input Contains the search positions, after function
 *             call the fields Task and Index are set.
 *  \param[in] nn Number of items in searchdata.
 *  \param[in] hsmlguess Guess for initial search radius;
 *             1: from searchdata; else from MeanVolume of cells.
 *  \param[in] verbose More output.
 *
 *  \return void
 */
void find_nearest_meshpoint_global(mesh_search_data *searchdata_input, int nn, int hsmlguess, int verbose)
{
  int i;
  n                      = nn;
  ngbsearch_nearest_dist = (MyDouble *)mymalloc("ngbsearch_nearest_dist", n * sizeof(MyDouble));
  ngbsearch_hsml         = (MyDouble *)mymalloc("ngbsearch_hsml", n * sizeof(MyDouble));
  searchdata             = searchdata_input;

  for(i = 0; i < n; i++)
    {
      ngbsearch_nearest_dist[i] = MAX_REAL_NUMBER;

      if(hsmlguess)
        ngbsearch_hsml[i] = searchdata[i].u.hsmlguess;
      else
        ngbsearch_hsml[i] = 1e-6 * pow(All.MeanVolume, 1.0 / 3);

      searchdata[i].Task = -1;  // None found yet
    }

  generic_set_MaxNexport();

  int ntot, iter = 0;

  /* we will repeat the whole thing for those points where we did not find a nearest neighbor */
  do
    {
      generic_comm_pattern(n, kernel_local, kernel_imported);

      int npleft = 0;

      /* do final operations on results */
      for(i = 0; i < n; i++)
        {
          if(searchdata[i].Task == -1)
            {
              npleft++;
              ngbsearch_hsml[i] *= 2.0;

              if(iter >= MAXITER - 10)
                {
                  printf("i=%d task=%d hsml=%g nearest dist=%g pos=(%g|%g|%g)\n", i, ThisTask, ngbsearch_hsml[i],
                         ngbsearch_nearest_dist[i], searchdata[i].Pos[0], searchdata[i].Pos[1], searchdata[i].Pos[2]);
                  myflush(stdout);
                }
              if(iter > MAXITER)
                terminate("NGBSEARCH: iter > MAXITER");
            }
        }

      /* sum up the left overs */
      MPI_Allreduce(&npleft, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      if(ntot > 0) /* ok, we need to repeat for a few particles */
        {
          iter++;
          if(iter > 0 && ThisTask == 0 && verbose)
            {
              printf("NGBSEARCH: iteration %d: need to repeat for %d points.\n", iter, ntot);
              myflush(stdout);
            }

          if(iter > MAXITER)
            terminate("NGBSEARCH: failed to converge in tracer particles\n");
        }
    }
  while(ntot > 0);

  myfree(ngbsearch_hsml);
  myfree(ngbsearch_nearest_dist);
}

/**
 * performs the search
 *
 * @param target the index of the particle to process(mode 0: in searchdata, mode 1: in NgbSearchDataGet/Result)
 * @param mode either 0 (handle local particles) or 1 (handle particles sent to us)
 * @param nexport number of particles to export to other task (only if mode == 0)
 * @param nsend_local specifies how many particles are send to which task (only if mode == 0)
 * @param searchdata list of positions to search (only if mode == 0)
 */
int ngbsearch_primary_cell_evaluate(int target, int mode, int threadid)
{
  int j;
  int numnodes, *firstnode;
  MyDouble h, distmax;
  MyDouble dx, dy, dz, r;
  MyDouble *pos;
  data_in local, *target_data;
  data_out out;

  int index = -1;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data = &local;

      numnodes  = 1;
      firstnode = NULL;
    }
  else
    {
      target_data = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  pos     = target_data->pos;
  h       = target_data->hsml;
  distmax = target_data->distance;

  int numngb = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

  for(int k = 0; k < numngb; k++)
    {
      j = Thread[threadid].Ngblist[k];

      /* find the closest image in the given box size */
      dx = NEAREST_X(pos[0] - P[j].Pos[0]);
      dy = NEAREST_Y(pos[1] - P[j].Pos[1]);
      dz = NEAREST_Z(pos[2] - P[j].Pos[2]);
      r  = sqrt(dx * dx + dy * dy + dz * dz);
      if(r < distmax && r < h && P[j].ID != 0 && P[j].Mass > 0)
        {
          distmax = r;
          index   = j;
        }
    }

  out.Distance = distmax;
  out.Task     = ThisTask;
  out.Index    = index;

  if(index < 0)
    {
      out.Distance = MAX_REAL_NUMBER;
      out.Task     = -1;
      out.Index    = -1;
    }

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}
