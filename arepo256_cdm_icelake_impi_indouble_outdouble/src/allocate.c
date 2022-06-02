/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/allocate.c
 * \date        05/2018
 * \author      Volker Springel
 * \brief       Functions to allocate and reallocate global arrays.
 * \details     contains functions
 *                void allocate_memory(void)
 *                void reallocate_memory_maxpart(void)
 *                void reallocate_memory_maxpartsph(void)
 *                void reallocate_memory_maxpartstar(void)
 *                void reallocate_memory_maxpartBHs(void)
 *                void reallocate_memory_maxpartsinks(void)
 *                void reallocate_memory_maxpartdust(void)
 *                void reallocate_memory_maxparttracer(int delta)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 03.05.2018 Function documentation -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

/*! \brief Allocates memory for global arrays.
 *
 *  This routine allocates memory for
 *  - particle storage, both the collisionless and the cells (SPH particles),
 *  - the ordered binary tree of the timeline,
 *  - communication buffers.
 *
 *  \return void
 */
void allocate_memory(void)
{
  int NTaskTimesThreads;

  NTaskTimesThreads = MaxThreads * NTask;

  Exportflag      = (int *)mymalloc("Exportflag", MaxThreads * ((((NTask - 1) / 16) + 1) * 16) * sizeof(int));
  Exportindex     = (int *)mymalloc("Exportindex", NTaskTimesThreads * sizeof(int));
  Exportnodecount = (int *)mymalloc("Exportnodecount", NTaskTimesThreads * sizeof(int));

  Send = (struct send_recv_counts *)mymalloc("Send", sizeof(struct send_recv_counts) * NTask);
  Recv = (struct send_recv_counts *)mymalloc("Recv", sizeof(struct send_recv_counts) * NTask);

  TasksThatSend = (int *)mymalloc("TasksThatSend", sizeof(int) * NTask);
  TasksThatRecv = (int *)mymalloc("TasksThatRecv", sizeof(int) * NTask);

  Send_count  = (int *)mymalloc("Send_count", sizeof(int) * NTaskTimesThreads);
  Send_offset = (int *)mymalloc("Send_offset", sizeof(int) * NTaskTimesThreads);
  Recv_count  = (int *)mymalloc("Recv_count", sizeof(int) * NTask);
  Recv_offset = (int *)mymalloc("Recv_offset", sizeof(int) * NTask);

  Send_count_nodes  = (int *)mymalloc("Send_count_nodes", sizeof(int) * NTask);
  Send_offset_nodes = (int *)mymalloc("Send_offset_nodes", sizeof(int) * NTask);
  Recv_count_nodes  = (int *)mymalloc("Recv_count_nodes", sizeof(int) * NTask);
  Recv_offset_nodes = (int *)mymalloc("Recv_offset_nodes", sizeof(int) * NTask);

  Mesh_Send_count  = (int *)mymalloc("Mesh_Send_count", sizeof(int) * NTask);
  Mesh_Send_offset = (int *)mymalloc("Mesh_Send_offset", sizeof(int) * NTask);
  Mesh_Recv_count  = (int *)mymalloc("Mesh_Recv_count", sizeof(int) * NTask);
  Mesh_Recv_offset = (int *)mymalloc("Mesh_Recv_offset", sizeof(int) * NTask);

  Force_Send_count  = (int *)mymalloc("Force_Send_count", sizeof(int) * NTask);
  Force_Send_offset = (int *)mymalloc("Force_Send_offset", sizeof(int) * NTask);
  Force_Recv_count  = (int *)mymalloc("Force_Recv_count", sizeof(int) * NTask);
  Force_Recv_offset = (int *)mymalloc("Force_Recv_offset", sizeof(int) * NTask);

  mpi_printf("ALLOCATE: initial allocation for MaxPart = %d\n", All.MaxPart);
  P = (struct particle_data *)mymalloc_movable(&P, "P", All.MaxPart * sizeof(struct particle_data));

  mpi_printf("ALLOCATE: initial allocation for MaxPartSph = %d\n", All.MaxPartSph);
  SphP = (struct sph_particle_data *)mymalloc_movable(&SphP, "SphP", All.MaxPartSph * sizeof(struct sph_particle_data));

#ifdef CHIMES
  ChimesAbundances    = (double *)mymalloc_movable(&ChimesAbundances, "ChimesAbundances",
                                                   All.MaxPartSph * ChimesGlobalVars.totalNumberOfSpecies * sizeof(double));
  ChimesPhotonDensity = (double *)mymalloc_movable(&ChimesPhotonDensity, "ChimesPhotonDensity",
                                                   All.MaxPartSph * ChimesGlobalVars.N_spectra * sizeof(double));
  ChimesDustG = (double *)mymalloc_movable(&ChimesDustG, "ChimesDustG", All.MaxPartSph * ChimesGlobalVars.N_spectra * sizeof(double));
  ChimesH2dissocJ =
      (double *)mymalloc_movable(&ChimesH2dissocJ, "ChimesH2dissocJ", All.MaxPartSph * ChimesGlobalVars.N_spectra * sizeof(double));

#ifdef CHIMES_ADVECT_ABUNDANCES
  ChimesIonAdvect = (double *)mymalloc_movable(&ChimesIonAdvect, "ChimesIonAdvect",
                                               All.MaxPartSph * ChimesGlobalVars.totalNumberOfSpecies * sizeof(double));
#endif
#endif

#ifdef TRACER_MC
  mpi_printf("ALLOCATE: initial allocation for MaxPartTracer = %d\n", All.MaxPartTracer);
  TracerLinkedList     = (struct tracer_linked_list *)mymalloc_movable(&TracerLinkedList, "TracerLinkedList",
                                                                       All.MaxPartTracer * sizeof(struct tracer_linked_list));
  TracerLinkedListHeap = (int *)mymalloc_movable(&TracerLinkedListHeap, "TracerLinkedListHeap", All.MaxPartTracer * sizeof(int));
  for(int i = 0; i < All.MaxPartTracer; i++)
    TracerLinkedListHeap[i] = i;
#endif

#if defined(GFM) || defined(SFR_MCS)
  mpi_printf("ALLOCATE: initial allocation for MaxPartStar = %d\n", All.MaxPartStar);
  StarP = (struct star_particle_data *)mymalloc_movable(&StarP, "StarP", All.MaxPartStar * sizeof(struct star_particle_data));
#endif

#ifdef BLACK_HOLES
  mpi_printf("ALLOCATE: initial allocation for MaxPartBHs = %d\n", All.MaxPartBHs);
  BHP = (struct bh_particle_data *)mymalloc_movable(&BHP, "BHP", All.MaxPartBHs * sizeof(struct bh_particle_data));
#endif

#ifdef SINKS
  mpi_printf("ALLOCATE: initial allocation for MaxPartSinks = %d\n", All.MaxPartSinks);
  SinkP = (struct sink_particle_data *)mymalloc_movable(&SinkP, "SinkP", All.MaxPartSinks * sizeof(struct sink_particle_data));
#endif

#ifdef SINK_PARTICLES
  mpi_printf("ALLOCATE: initial allocation for NSinkBufferSize = %d\n", NSinkBufferSize);
  SinkP = (struct global_sink_particle_data *)mymalloc_movable(&SinkP, "SinkP",
                                                               NSinkBufferSize * sizeof(struct global_sink_particle_data));
#endif

#ifdef DUST_LIVE
  mpi_printf("ALLOCATE: initial allocation for MaxPartDust = %d\n", All.MaxPartDust);
  DustP = (struct dust_particle_data *)mymalloc_movable(&DustP, "DustP", All.MaxPartDust * sizeof(struct dust_particle_data));
#endif

#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
  PartSpecialListGlobal = (struct special_particle_data *)mymalloc_movable(&PartSpecialListGlobal, "PartSpecialListGlobal",
                                                                           All.MaxPartSpecial * sizeof(struct special_particle_data));
#endif

#ifdef REFINEMENT_AROUND_DM
  mpi_printf("ALLOCATE: initial allocation for TotPartDM = %d\n", All.TotPartDM);
  DMPartListGlobal =
      (struct refine_dm_data *)mymalloc_movable(&DMPartListGlobal, "DMPartListGlobal", All.TotPartDM * sizeof(struct refine_dm_data));
#endif

#if defined(CIRCUMSTELLAR) &&                                                                                 \
    (defined(CIRCUMSTELLAR_IRRADIATION) || defined(ALPHA_VISCOSITY) || defined(CIRCUMSTELLAR_REFINEMENTS)) && \
    !defined(EXTERNALGRAVITY)
  SourcePartListGlobal = (struct source_particle_data *)mymalloc_movable(&SourcePartListGlobal, "SourcePartListGlobal",
                                                                         All.MaxPartSources * sizeof(struct source_particle_data));
#endif

#if(NUM_THREADS > 1)
  ParticleLocks = (omp_lock_t *)mymalloc_movable(&ParticleLocks, "ParticleLocks", All.MaxPart * sizeof(omp_lock_t));
  int i;
  for(i = 0; i < All.MaxPart; i++)
    omp_init_lock(&ParticleLocks[i]);
#endif

  timebins_allocate(&TimeBinsHydro);
  timebins_allocate(&TimeBinsGravity);

#ifdef TRACER_PARTICLE
  timebins_allocate(&TimeBinsTracer);
#endif

#ifdef BLACK_HOLES
  timebins_allocate(&TimeBinsBHAccretion);
#endif

#ifdef SINKS
  timebins_allocate(&TimeBinsSinksAccretion);
#endif

#ifdef DUST_LIVE
  timebins_allocate(&TimeBinsDust);
#endif

  /* set to zero */
  memset(P, 0, All.MaxPart * sizeof(struct particle_data));
  memset(SphP, 0, All.MaxPartSph * sizeof(struct sph_particle_data));

#ifdef TRACER_MC
  memset(TracerLinkedList, 0, All.MaxPartTracer * sizeof(struct tracer_linked_list));
#endif

#if defined(GFM) || defined(SFR_MCS)
  memset(StarP, 0, All.MaxPartStar * sizeof(struct star_particle_data));
#endif

#ifdef BLACK_HOLES
  memset(BHP, 0, All.MaxPartBHs * sizeof(struct bh_particle_data));
#endif

#ifdef SINKS
  memset(SinkP, 0, All.MaxPartSinks * sizeof(struct sink_particle_data));
#endif

#ifdef SINK_PARTICLES
  memset(SinkP, 0, NSinkBufferSize * sizeof(struct global_sink_particle_data));
#endif

#ifdef DUST_LIVE
  memset(DustP, 0, All.MaxPartDust * sizeof(struct dust_particle_data));
#endif
}

/*! \brief Reallocates memory for particle data.
 *
 *  Reallocates memory for P and TimeBinsGravity arrays.
 *
 *  \return void
 */
void reallocate_memory_maxpart(void)
{
  mpi_printf("ALLOCATE: Changing to MaxPart = %d\n", All.MaxPart);

  P = (struct particle_data *)myrealloc_movable(P, All.MaxPart * sizeof(struct particle_data));
  timebins_reallocate(&TimeBinsGravity);

#ifdef TRACER_PARTICLE
  timebins_reallocate(&TimeBinsTracer);
#endif

#ifdef BLACK_HOLES
  timebins_reallocate(&TimeBinsBHAccretion);
#endif

#ifdef SINKS
  timebins_reallocate(&TimeBinsSinksAccretion);
#endif

#ifdef DUST_LIVE
  timebins_reallocate(&TimeBinsDust);
#endif

#if(NUM_THREADS > 1)
  ParticleLocks = (omp_lock_t *)myrealloc_movable(ParticleLocks, All.MaxPart * sizeof(omp_lock_t));
  int i;
  for(i = 0; i < All.MaxPart; i++)
    omp_init_lock(&ParticleLocks[i]);
#endif
}

/*! \brief Reallocates memory for particle data.
 *
 *  Reallocates memory for P only.
 *
 *  \return void
 */
void reallocate_memory_maxpart_ignore_timebins(void)
{
  mpi_printf("ALLOCATE: Changing to MaxPart = %d\n", All.MaxPart);

  P = (struct particle_data *)myrealloc_movable(P, All.MaxPart * sizeof(struct particle_data));

#if(NUM_THREADS > 1)
  ParticleLocks = (omp_lock_t *)myrealloc_movable(ParticleLocks, All.MaxPart * sizeof(omp_lock_t));
  int i;
  for(i = 0; i < All.MaxPart; i++)
    omp_init_lock(&ParticleLocks[i]);
#endif
}

/*! \brief Reallocate memory for cell data.
 *
 *  Reallocates memory for cells in SphP and TimeBinsHydro arrays.
 *
 *  \return void
 */
void reallocate_memory_maxpartsph(void)
{
  mpi_printf("ALLOCATE: Changing to MaxPartSph = %d\n", All.MaxPartSph);

  SphP = (struct sph_particle_data *)myrealloc_movable(SphP, All.MaxPartSph * sizeof(struct sph_particle_data));
  timebins_reallocate(&TimeBinsHydro);

#ifdef CHIMES
  ChimesAbundances =
      (double *)myrealloc_movable(ChimesAbundances, All.MaxPartSph * ChimesGlobalVars.totalNumberOfSpecies * sizeof(double));
  ChimesPhotonDensity = (double *)myrealloc_movable(ChimesPhotonDensity, All.MaxPartSph * ChimesGlobalVars.N_spectra * sizeof(double));
  ChimesDustG         = (double *)myrealloc_movable(ChimesDustG, All.MaxPartSph * ChimesGlobalVars.N_spectra * sizeof(double));
  ChimesH2dissocJ     = (double *)myrealloc_movable(ChimesH2dissocJ, All.MaxPartSph * ChimesGlobalVars.N_spectra * sizeof(double));

#ifdef CHIMES_ADVECT_ABUNDANCES
  ChimesIonAdvect =
      (double *)myrealloc_movable(ChimesIonAdvect, All.MaxPartSph * ChimesGlobalVars.totalNumberOfSpecies * sizeof(double));
#endif

  chimes_update_all_pointers();
#endif
}

/*! \brief Reallocate memory for cell data.
 *
 *  Reallocates memory for cells in SphP only.
 *
 *  \return void
 */
void reallocate_memory_maxpartsph_ignore_timebins(void)
{
  mpi_printf("ALLOCATE: Changing to MaxPartSph = %d\n", All.MaxPartSph);

  SphP = (struct sph_particle_data *)myrealloc_movable(SphP, All.MaxPartSph * sizeof(struct sph_particle_data));

#ifdef CHIMES
  ChimesAbundances =
      (double *)myrealloc_movable(ChimesAbundances, All.MaxPartSph * ChimesGlobalVars.totalNumberOfSpecies * sizeof(double));
  ChimesPhotonDensity = (double *)myrealloc_movable(ChimesPhotonDensity, All.MaxPartSph * ChimesGlobalVars.N_spectra * sizeof(double));
  ChimesDustG         = (double *)myrealloc_movable(ChimesDustG, All.MaxPartSph * ChimesGlobalVars.N_spectra * sizeof(double));
  ChimesH2dissocJ     = (double *)myrealloc_movable(ChimesH2dissocJ, All.MaxPartSph * ChimesGlobalVars.N_spectra * sizeof(double));

#ifdef CHIMES_ADVECT_ABUNDANCES
  ChimesIonAdvect =
      (double *)myrealloc_movable(ChimesIonAdvect, All.MaxPartSph * ChimesGlobalVars.totalNumberOfSpecies * sizeof(double));
#endif

  chimes_update_all_pointers();
#endif
}

#if defined(GFM) || defined(SFR_MCS)
void reallocate_memory_maxpartstar(void)
{
  mpi_printf("ALLOCATE: Changing to MaxPartStar = %d\n", All.MaxPartStar);

  StarP = (struct star_particle_data *)myrealloc_movable(StarP, All.MaxPartStar * sizeof(struct star_particle_data));
}
#endif

#ifdef BLACK_HOLES
void reallocate_memory_maxpartBHs(void)
{
  mpi_printf("ALLOCATE: Changing to MaxPartBHs = %d\n", All.MaxPartBHs);

  BHP = (struct bh_particle_data *)myrealloc_movable(BHP, All.MaxPartBHs * sizeof(struct bh_particle_data));
}
#endif

#ifdef SINKS
void reallocate_memory_maxpartsinks(void)
{
  mpi_printf("ALLOCATE: Changing to MaxPartSinks = %d\n", All.MaxPartSinks);

  SinkP = (struct sink_particle_data *)myrealloc_movable(SinkP, All.MaxPartSinks * sizeof(struct sink_particle_data));
}
#endif

#ifdef DUST_LIVE
void reallocate_memory_maxpartdust(void)
{
  mpi_printf("ALLOCATE: Changing to MaxPartDust = %d\n", All.MaxPartDust);

  DustP = (struct dust_particle_data *)myrealloc_movable(DustP, All.MaxPartDust * sizeof(struct dust_particle_data));
}
#endif

#ifdef TRACER_MC
/*! \brief Reallocates memory for Monte-Carlo tracers.
 *
 *  Reallocaes the memory in TracerLinkedListHeap and TracerLinkedList.
 *
 *  \return void
 */
void reallocate_memory_maxparttracer(int delta)
{
  if(delta < 0)
    terminate("delta=%d", delta);

  int old = All.MaxPartTracer;
  All.MaxPartTracer += delta;

  mpi_printf("ALLOCATE: Changing to MaxPartTracer = %d\n", All.MaxPartTracer);

  TracerLinkedListHeap = (int *)myrealloc_movable(TracerLinkedListHeap, All.MaxPartTracer * sizeof(int));
  for(int i = old; i < All.MaxPartTracer; i++)
    TracerLinkedListHeap[i] = i;

  TracerLinkedList =
      (struct tracer_linked_list *)myrealloc_movable(TracerLinkedList, All.MaxPartTracer * sizeof(struct tracer_linked_list));
}
#endif
