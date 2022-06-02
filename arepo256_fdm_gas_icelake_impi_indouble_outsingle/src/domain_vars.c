/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/domain_vars.c
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

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "allvars.h"
#include "domain.h"
#include "proto.h"
#include "voronoi.h"

struct domain_peano_hilbert_data *mp;

struct local_topnode_data *topNodes, *branchNodes; /*!< points to the root node of the top-level tree */

double totgravcost, totpartcount, gravcost, totsphcost, sphcost;

struct domain_cost_data *DomainLeaveNode;

double fac_work, fac_load, fac_worksph;
double normsum_work, normsum_load, normsum_worksph;

int Nbranch;

/*! toGo[partner] gives the number of particles on the current task that have to go to task 'partner'
 */
int *toGo, *toGoSph;
int *toGet, *toGetSph;
int *list_NumPart;
int *list_NumGas;
int *list_load;
int *list_loadsph;
double *list_work;
double *list_worksph;

#ifdef TRACER_MC
int MaxNumTracer;
int *toGoTracer;
int *toGetTracer;
int *list_N_Tracer;
#endif
#if defined(GFM) || defined(SFR_MCS)
int *toGoStar;
int *toGetStar;
int *list_N_star;
int *list_loadstar;
#endif
#ifdef BLACK_HOLES
int *toGoBHs;
int *toGetBHs;
int *list_NumBHs;
int *list_loadBHs;
#endif
#ifdef SINKS
int *toGoSinks;
int *toGetSinks;
int *list_NumSinks;
int *list_loadSinks;
#endif
#ifdef DUST_LIVE
int *toGoDust;
int *toGetDust;
int *list_N_dust;
int *list_loaddust;
#endif

/*! \brief Allocates lists needed for communication in domain decomposition.
 *
 *  These lists are holding information about other tasks (number of particles,
 *  load, work, etc.).
 *
 *  \return void
 */
void domain_allocate_lists(void)
{
  Key             = (peanokey *)mymalloc_movable(&Key, "domain_key", (sizeof(peanokey) * All.MaxPart));
  toGo            = (int *)mymalloc_movable(&toGo, "toGo", (sizeof(int) * NTask));
  toGoSph         = (int *)mymalloc_movable(&toGoSph, "toGoSph", (sizeof(int) * NTask));
  toGet           = (int *)mymalloc_movable(&toGet, "toGet", (sizeof(int) * NTask));
  toGetSph        = (int *)mymalloc_movable(&toGetSph, "toGetSph", (sizeof(int) * NTask));
  list_NumPart    = (int *)mymalloc_movable(&list_NumPart, "list_NumPart", (sizeof(int) * NTask));
  list_NumGas     = (int *)mymalloc_movable(&list_NumGas, "list_NumGas", (sizeof(int) * NTask));
  list_load       = (int *)mymalloc_movable(&list_load, "list_load", (sizeof(int) * NTask));
  list_loadsph    = (int *)mymalloc_movable(&list_loadsph, "list_loadsph", (sizeof(int) * NTask));
  list_work       = (double *)mymalloc_movable(&list_work, "list_work", (sizeof(double) * NTask));
  list_worksph    = (double *)mymalloc_movable(&list_worksph, "list_worksph", (sizeof(double) * NTask));
  DomainLeaveNode = (struct domain_cost_data *)mymalloc_movable(&DomainLeaveNode, "DomainLeaveNode",
                                                                (MaxTopNodes * sizeof(struct domain_cost_data)));
#ifdef TRACER_MC
  toGoTracer    = (int *)mymalloc_movable(&toGoTracer, "toGoTracer", (sizeof(int) * NTask));
  toGetTracer   = (int *)mymalloc_movable(&toGetTracer, "toGetTracer", (sizeof(int) * NTask));
  list_N_Tracer = (int *)mymalloc_movable(&list_N_Tracer, "list_N_Tracer", (sizeof(int) * NTask));
#endif
#if defined(GFM) || defined(SFR_MCS)
  toGoStar      = (int *)mymalloc_movable(&toGoStar, "toGoStar", (sizeof(int) * NTask));
  toGetStar     = (int *)mymalloc_movable(&toGetStar, "toGetStar", (sizeof(int) * NTask));
  list_N_star   = (int *)mymalloc_movable(&list_N_star, "list_N_star", (sizeof(int) * NTask));
  list_loadstar = (int *)mymalloc_movable(&list_loadstar, "list_loadstar", (sizeof(int) * NTask));
#endif
#ifdef BLACK_HOLES
  toGoBHs      = (int *)mymalloc_movable(&toGoBHs, "toGoBHs", (sizeof(int) * NTask));
  toGetBHs     = (int *)mymalloc_movable(&toGetBHs, "toGetBHs", (sizeof(int) * NTask));
  list_NumBHs  = (int *)mymalloc_movable(&list_NumBHs, "list_NumBHs", (sizeof(int) * NTask));
  list_loadBHs = (int *)mymalloc_movable(&list_loadBHs, "list_loadBHs", (sizeof(int) * NTask));
#endif
#ifdef SINKS
  toGoSinks      = (int *)mymalloc_movable(&toGoSinks, "toGoSinks", (sizeof(int) * NTask));
  toGetSinks     = (int *)mymalloc_movable(&toGetSinks, "toGetSinks", (sizeof(int) * NTask));
  list_NumSinks  = (int *)mymalloc_movable(&list_NumSinks, "list_NumSinks", (sizeof(int) * NTask));
  list_loadSinks = (int *)mymalloc_movable(&list_loadSinks, "list_loadSinks", (sizeof(int) * NTask));
#endif
#ifdef DUST_LIVE
  toGoDust      = (int *)mymalloc_movable(&toGoDust, "toGoDust", (sizeof(int) * NTask));
  toGetDust     = (int *)mymalloc_movable(&toGetDust, "toGetDust", (sizeof(int) * NTask));
  list_N_dust   = (int *)mymalloc_movable(&list_N_dust, "list_N_dust", (sizeof(int) * NTask));
  list_loaddust = (int *)mymalloc_movable(&list_loaddust, "list_loaddust", (sizeof(int) * NTask));
#endif
}

/*! \brief Frees lists needed for communication in domain decomposition.
 *
 *  This routine is the counterpart of domain_allocate_lists(void).
 *  Frees memory of all arrays allocated there, except Key, which is freed
 *  elsewhere (in void domain_Decomposition(void); see domain.c).
 *
 *  \return void
 */
void domain_free_lists(void)
{
#ifdef DUST_LIVE
  myfree(list_loaddust);
  myfree(list_N_dust);
  myfree(toGetDust);
  myfree(toGoDust);
#endif
#ifdef SINKS
  myfree(list_loadSinks);
  myfree(list_NumSinks);
  myfree(toGetSinks);
  myfree(toGoSinks);
#endif
#ifdef BLACK_HOLES
  myfree(list_loadBHs);
  myfree(list_NumBHs);
  myfree(toGetBHs);
  myfree(toGoBHs);
#endif
#if defined(GFM) || defined(SFR_MCS)
  myfree(list_loadstar);
  myfree(list_N_star);
  myfree(toGetStar);
  myfree(toGoStar);
#endif
#ifdef TRACER_MC
  myfree(list_N_Tracer);
  myfree(toGetTracer);
  myfree(toGoTracer);
#endif
  myfree(DomainLeaveNode);
  myfree(list_worksph);
  myfree(list_work);
  myfree(list_loadsph);
  myfree(list_load);
  myfree(list_NumGas);
  myfree(list_NumPart);
  myfree(toGetSph);
  myfree(toGet);
  myfree(toGoSph);
  myfree(toGo);
}
