/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/voronoi_derefinement.c
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

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"

#ifdef VORONOI

#if defined(REFINEMENT_MERGE_CELLS) && !defined(ONEDIMS) && !defined(REFINEMENT_MERGE_PAIRS)

#define DEREF_SA_FAC 1.0e-4

static void derefine_add_ngb(int edge, int i, int j, double area);
static int derefine_compare_seq_DP_ID(const void *a, const void *b);
static int derefine_flux_list_data_compare(const void *a, const void *b);
static void derefine_apply_flux_list(void);
static void derefine_exchange_flag(void);
static void derefine_apply_probe_list(void);
static int derefine_probe_list_data_compare_task(const void *a, const void *b);
static void voronoi_refinement_check_enclosing_box(point *p, double *minpos, double *maxpos);

static struct derefine_particle_data
{
  int Flag;
  int dp_index;
} * deref_SphP;

static struct flagexch
{
  int Flag;
  MyIDType ID;
} * FlagExch;

static struct flag_delaunay_data
{
  int Flag;
} * flag_DP;

static struct seq_delaunay_data
{
  MyFloat rnd;
  int rank, index;
  MyIDType ID;
} * seq_DP;

static struct probe_list_data
{
  int task, index;
  int sendpart;
  int flag;
} * ProbeList;

static struct flux_list_data
{
  int task, index;
  double dM, dP[3];
#ifdef MHD
  double dB[3];
#ifdef MHD_CT
  double dA[3];
#endif
#ifdef MHD_DEDNER
  double dPsi;
#endif
#endif

#ifdef MRT
  double dDensPhot[MRT_BINS];
  double dRT_F[MRT_BINS][3];
#endif

#ifdef OUTPUT_CELL_SPIN
  double dSpin[3];
#endif
#ifndef ISOTHERM_EQS
  double dEnergy;
#endif
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  double dEntropy;
#endif
#ifdef TRACER_FIELD
  double dConservedTracer;
#endif
#if defined(SMUGGLE_STAR_FEEDBACK) && defined(SMUGGLE_OUTPUT_STELLAR_FEEDBACK)
  double dTotEgyFeed;
  double dIntEgyFeed;
  double dKinEgyFeed;
#endif
#ifdef BH_THERMALFEEDBACK
  double dInjected_BH_Energy;
#endif
#ifdef GFM_WINDS_LOCAL
  double dWindEnergyReceived;
#endif
#ifdef MAXSCALARS
  double dConservedScalars[MAXSCALARS];
#endif
#ifdef GFM_CHEMTAGS
  double dMassMetalsChemTags[GFM_N_CHEM_TAGS];
#endif
#ifdef TGCHEM
  double dVolume;
  double dAbund[TGCHEM_NUM_ABUNDANCES];
#endif
#ifdef CHIMES
  double dNHtot;
  double dN_ion[TOTSIZE];
#endif
#ifdef RUNGE_KUTTA_FULL_UPDATE
  struct conservative_variables rk;
#endif
} * FluxList;

static int Nflux, MaxNflux;

static int *first_ngb, *last_ngb, first_free_ngb;
static struct ngb_data
{
#ifdef OPTIMIZE_MEMORY_USAGE
  MyFloat area;
#else
  double area;
#endif
  int index;
  int edge;
  int next_ngb;
} * ngb;

#ifdef REFINEMENT_SPLIT_CELLS
extern char *FlagDoNotRefine;
#endif

/*! \brief Adds cell in list ngb.
 *
 *  \param[in] edge Element 'edge' in ngb.
 *  \param[in] i Index in first_ngb and last_ngb lists.
 *  \param[in] j Element 'index' in ngb.
 *  \param[in] area Element 'area' in ngb.
 *  \param[in] t Element 't' in ngb.
 *  \param[in] nr Element 'nr' in ngb.
 *
 *  \return void
 */
static void derefine_add_ngb(int edge, int i, int j, double area)
{
  if(i >= 0 && j >= 0)
    {
      if(i >= Mesh.Ndp || j >= Mesh.Ndp)
        {
          terminate("i>= Ndp || j>= Ndp");
        }

      if(first_ngb[i] >= 0)
        {
          ngb[last_ngb[i]].next_ngb = first_free_ngb;
          last_ngb[i]               = first_free_ngb;
        }
      else
        {
          first_ngb[i] = last_ngb[i] = first_free_ngb;
        }

      ngb[first_free_ngb].area     = area;
      ngb[first_free_ngb].edge     = edge;
      ngb[first_free_ngb].index    = j;
      ngb[first_free_ngb].next_ngb = -1;
      first_free_ngb++;
    }
}

int do_derefinements(void)
{
  int idx, i, j, k, count, countall;

#ifdef CHIMES
  double nHtot_i, nHtot_p, dNHtot, dN_ion, N_ion_p;
  int abunIndex;
#endif

  TIMER_START(CPU_DEREFINE);

  deref_SphP =
      (struct derefine_particle_data *)mymalloc_movable(&deref_SphP, "deref_SphP", NumGas * sizeof(struct derefine_particle_data));

  FlagExch = (struct flagexch *)mymalloc_movable(&FlagExch, "FlagExch", Mesh_nimport * sizeof(struct flagexch));

  /* first, check whether we have cells to derefine */
  for(idx = 0, count = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
#ifdef REFINEMENT_SPLIT_CELLS
      FlagDoNotRefine[i] = 0;
#endif

      if(i >= NumGas)
        terminate("index of gas cell greater than NumGas");

      deref_SphP[i].Flag     = 0;
      deref_SphP[i].dp_index = -1;

      if(derefine_should_this_cell_be_merged(i, deref_SphP[i].Flag))
        {
          deref_SphP[i].Flag = 1;
          count++;
        }
    }

  MPI_Allreduce(&count, &countall, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  mpi_printf("DEREFINE: Number of cells that want to be de-refined: %d\n", countall);

  if(countall)
    {
      /* tell the ghost cells whether they want to be refined or not */
      derefine_exchange_flag();

      /* let's create an explicit list of the neighbors of each cell */

      first_ngb = (int *)mymalloc("first_ngb", Mesh.Ndp * sizeof(int));
      ngb       = (struct ngb_data *)mymalloc("ngb", 2 * Mesh.Nvf * sizeof(struct ngb_data));

      last_ngb = (int *)mymalloc("last_ngb", Mesh.Ndp * sizeof(int));

      for(i = 0; i < Mesh.Ndp; i++)
        first_ngb[i] = last_ngb[i] = -1;

      for(i = 0, first_free_ngb = 0; i < Mesh.Nvf; i++)
        {
          derefine_add_ngb(i, Mesh.VF[i].p1, Mesh.VF[i].p2, Mesh.VF[i].area);
          derefine_add_ngb(i, Mesh.VF[i].p2, Mesh.VF[i].p1, Mesh.VF[i].area);
        }

      myfree(last_ngb);

      /* we now make a list of the delaunay points that we can sort in a globally unique way */
      flag_DP = (struct flag_delaunay_data *)mymalloc_movable(&flag_DP, "flag_DP", Mesh.Ndp * sizeof(struct flag_delaunay_data));
      seq_DP  = (struct seq_delaunay_data *)mymalloc("seq_DP", Mesh.Ndp * sizeof(struct seq_delaunay_data));

      for(i = 0; i < Mesh.Ndp; i++)
        {
          seq_DP[i].rank  = i;
          seq_DP[i].index = Mesh.DP[i].index;

          if(Mesh.DP[i].task == ThisTask)
            {
              int li = Mesh.DP[i].index;
              if(li < 0)
                {
                  flag_DP[i].Flag = 0;
                  seq_DP[i].ID    = 0;
                  seq_DP[i].rnd   = 0;
                }
              else
                {
                  if(li < NumGas)
                    if(deref_SphP[li].dp_index < 0)
                      deref_SphP[li].dp_index = i; /* only guaranteed to be set for active cells */

                  if(li >= NumGas)
                    li -= NumGas;

                  flag_DP[i].Flag = deref_SphP[li].Flag;
                  seq_DP[i].ID    = P[li].ID;
                  seq_DP[i].rnd   = get_random_number();
                }
            }
          else
            {
              flag_DP[i].Flag = FlagExch[Mesh.DP[i].index].Flag;
              seq_DP[i].ID    = FlagExch[Mesh.DP[i].index].ID;
              seq_DP[i].rnd   = get_random_number();
            }
        }

      /* sort according to ID */
      mysort(seq_DP, Mesh.Ndp, sizeof(struct seq_delaunay_data), derefine_compare_seq_DP_ID);

      /* now let's go through in sorted order. For each cell that is supposed to be refined, check whether any of the
       * neighbors is already refined. If yes, don't allow it to be refined.
       * Also, if there is a neighbour with the same ID, don't refine it, because this must be a mirrored particle
       */

      for(i = 0; i < Mesh.Ndp; i++)
        {
          if(seq_DP[i].ID != 0)
            {
              j = seq_DP[i].rank;

              if(flag_DP[j].Flag == 1) /* this cell is still eligible for derefinement */
                {
                  /* go through its neighbours and check whether one of them is already up for derefinement */

                  int n = 0;
                  k     = first_ngb[j];
                  while(k >= 0)
                    {
                      /* we only need to consider neighboring cells if they are active */
                      int q = ngb[k].index;

                      if(q >= 0)
                        {
                          int timebin;

                          if(Mesh.DP[q].task == ThisTask)
                            {
                              if(Mesh.DP[q].index < NumGas)
                                timebin = P[Mesh.DP[q].index].TimeBinHydro;
                              else
                                timebin = P[Mesh.DP[q].index - NumGas].TimeBinHydro;
                            }
                          else
                            {
#ifndef OPTIMIZE_MESH_MEMORY_FOR_REFINEMENT
                              timebin = PrimExch[Mesh.DP[q].index].TimeBinHydro;
#else
                              timebin = RefExch[Mesh.DP[q].index].TimeBinHydro;
#endif
                            }

                          if(TimeBinSynchronized[timebin])
                            {
                              if(flag_DP[q].Flag == 2 || flag_DP[q].Flag == 3)
                                n++;

                              if(Mesh.DP[q].ID == seq_DP[i].ID) /* same ID, so we have a mirrored particle */
                                n++;
                            }
                        }

                      k = ngb[k].next_ngb;
                    }

                  if(n == 0) /* ok, none have been found. This means this cell is allowed to be refined */
                    flag_DP[j].Flag = 2;
                  else
                    flag_DP[j].Flag = 3;
                }
            }
        }

      myfree(seq_DP);

      /* copy of the refinement flags to the cell structure */
      for(i = 0; i < Mesh.Ndp; i++)
        if(Mesh.DP[i].task == ThisTask && Mesh.DP[i].index >= 0 && Mesh.DP[i].index < NumGas)
          deref_SphP[Mesh.DP[i].index].Flag = flag_DP[i].Flag;

      myfree(flag_DP);

      /* now let's count again how many cells we would like to derefine */

      for(idx = 0, count = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
            continue;

          if(deref_SphP[i].Flag == 2)
            count++;
        }

      int in[2], out[2];
      in[0] = count;

      /* now we carry out an auxiliary check to make sure that we really
         avoid de-refining two neighboring cells.  If such a pair is
         found, both cells will not be derefined. */

      MaxNflux  = Mesh.Indi.AllocFacNflux;
      Nflux     = 0;
      ProbeList = (struct probe_list_data *)mymalloc_movable(&ProbeList, "ProbeList", MaxNflux * sizeof(struct probe_list_data));

      count = 0;

      for(idx = 0, count = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
            continue;

          if(deref_SphP[i].Flag == 2)
            {
              j = deref_SphP[i].dp_index; /* this is the delaunay point of this cell */
              if(j < 0)
                terminate("j < 0");

              k = first_ngb[j];

              int flag = 0;

              while(k >= 0)
                {
                  if(ngb[k].area > DEREF_SA_FAC * SphP[i].SurfaceArea)
                    {
                      int q = ngb[k].index;

                      if(Mesh.DP[q].task == ThisTask)
                        {
                          int p = Mesh.DP[q].index;

                          if(p < 0)
                            terminate("p < 0");

                          if(p >= NumGas) /* this is a local ghost point */
                            p -= NumGas;

                          if(TimeBinSynchronized[P[p].TimeBinHydro])
                            if(deref_SphP[p].Flag == 2)
                              flag++;
                        }
                      else
                        {
                          /* here we have a foreign ghost point */
                          if(Nflux >= MaxNflux)
                            {
                              Mesh.Indi.AllocFacNflux *= ALLOC_INCREASE_FACTOR;
                              MaxNflux = Mesh.Indi.AllocFacNflux;
#ifdef VERBOSE
                              printf("Task=%d: increase memory allocation, MaxNflux=%d Indi.AllocFacNflux=%g\n", ThisTask, MaxNflux,
                                     Mesh.Indi.AllocFacNflux);
#endif
                              ProbeList =
                                  (struct probe_list_data *)myrealloc_movable(ProbeList, MaxNflux * sizeof(struct probe_list_data));

                              if(Nflux >= MaxNflux)
                                terminate("Nflux >= MaxNflux");
                            }

                          ProbeList[Nflux].task     = Mesh.DP[q].task;
                          ProbeList[Nflux].index    = Mesh.DP[q].originalindex;
                          ProbeList[Nflux].sendpart = i;
                          ProbeList[Nflux].flag     = 0;

                          Nflux++;
                        }
                    }
                  k = ngb[k].next_ngb;
                }

              if(flag)
                {
                  /* ups. It looks like a neigboring point is also about to be dissolved. We hence do not
                     dissolve the current point
                   */
                  deref_SphP[i].Flag = 0;
                  count++;
                }
            }
        }

      /* now let's probe on other tasks */

      derefine_apply_probe_list();

      for(i = 0; i < Nflux; i++)
        {
          if(ProbeList[i].flag)
            if(deref_SphP[ProbeList[i].sendpart].Flag == 2)
              {
                deref_SphP[ProbeList[i].sendpart].Flag = 0;
                count++;
              }
        }

      myfree(ProbeList);

      in[1] = count;
      MPI_Reduce(in, out, 2, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      mpi_printf("DEREFINE: Number of cells that we could de-refine: %d, number of cells we exclude from this set:  %d\n", out[0],
                 out[1]);

      /* we now distribute the conserved quantities of the cell among the neighbours */

      MaxNflux = Mesh.Indi.AllocFacNflux;
      Nflux    = 0;
      FluxList = (struct flux_list_data *)mymalloc_movable(&FluxList, "FluxList", MaxNflux * sizeof(struct flux_list_data));
#ifdef TRACER_MC
      start_MC_tracer(N_tracer);
#endif

      for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
            continue;

          if(deref_SphP[i].Flag == 2)
            {
              j = deref_SphP[i].dp_index; /* this is the delaunay point of this cell */
              if(j < 0)
                terminate("j < 0");

              /* construct local mesh without cell that is to be derefined, then with cell, then compute volume difference for
               * surrounding cells */
              initialize_and_create_first_tetra(&DeRefMesh);

              DeRefMesh.DTC = (tetra_center *)mymalloc_movable(&DeRefMesh.DTC, "DeRefDTC", DeRefMesh.MaxNdt * sizeof(tetra_center));
              DeRefMesh.DTF = (char *)mymalloc_movable(&DeRefMesh.DTF, "DeRefDTF", DeRefMesh.MaxNdt * sizeof(char));
              for(k = 0; k < DeRefMesh.Ndt; k++)
                DeRefMesh.DTF[k] = 0;

              int tlast = 0;

              double minpos[3], maxpos[3];
              for(int i = 0; i < 3; i++)
                {
                  minpos[i] = MAX_DOUBLE_NUMBER;
                  maxpos[i] = 0;
                }

              /* add all direct neighbours */
              k = first_ngb[j];
              while(k >= 0)
                {
                  int q = ngb[k].index;

                  voronoi_check_and_increase_dp_array(&DeRefMesh, 1);

                  DeRefMesh.DP[DeRefMesh.Ndp] = Mesh.DP[q];

#ifndef OPTIMIZE_MEMORY_USAGE
                  set_integers_for_point(&DeRefMesh, DeRefMesh.Ndp);
#endif
                  tlast = insert_point(&DeRefMesh, DeRefMesh.Ndp, tlast);

                  voronoi_refinement_check_enclosing_box(&Mesh.DP[q], minpos, maxpos);

                  DeRefMesh.Ndp++;
                  k = ngb[k].next_ngb;
                }

              int NNgb = DeRefMesh.Ndp;

              /* now add a box around them */
              double boxcenter[3];
              double boxsize[3];
              for(int i = 0; i < 3; i++)
                {
                  boxcenter[i] = 0.5 * (maxpos[i] + minpos[i]);
                  boxsize[i]   = 0.75 * (maxpos[i] - minpos[i]);
                }

              voronoi_check_and_increase_dp_array(&DeRefMesh, (1 << NUMDIMS) + (1 + 2 * (NUMDIMS - 1)) * 2 * NUMDIMS);
              for(int k = 0; k < (1 << NUMDIMS); k++)
                {
                  point p;
                  p.x           = boxcenter[0] + (4. * ((k & 1) > 0) - 2.) * boxsize[0];
                  p.y           = boxcenter[1] + (4. * ((k & 2) > 0) - 2.) * boxsize[1];
                  p.z           = boxcenter[2] + (4. * ((k & 4) > 0) - 2.) * boxsize[2];
                  p.task        = ThisTask;
                  p.index       = i;
                  p.image_flags = 0;

                  DeRefMesh.DP[DeRefMesh.Ndp] = p;

#ifndef OPTIMIZE_MEMORY_USAGE
                  set_integers_for_point(&DeRefMesh, DeRefMesh.Ndp);
#endif
                  tlast = insert_point(&DeRefMesh, DeRefMesh.Ndp, tlast);

                  DeRefMesh.Ndp++;
                }

              // for each face
              for(int k = 0; k < 2 * NUMDIMS; k++)
                {
                  int ix = (k / 2 == 0) * (2 * (k % 2) - 1);
                  int iy = (k / 2 == 1) * (2 * (k % 2) - 1);
                  int iz = (k / 2 == 2) * (2 * (k % 2) - 1);

                  for(int l = 0; l < 1 + 2 * (NUMDIMS - 1); l++)
                    {
                      double ivec[3];
                      ivec[0] = ix;
                      ivec[1] = iy;
                      ivec[2] = iz;
#if NUMDIMS == 2
                      if(l > 0)
                        ivec[1 - k / 2] = -0.5 + (l - 1);

#elif NUMDIMS == 3
                      if(l > 0)
                        {
                          ivec[(k / 2 + 1) % 3] = -0.5 + (l - 1) / 2;
                          ivec[(k / 2 + 2) % 3] = -0.5 + (l - 1) % 2;
                        }
#endif
                      point p;
                      p.x           = boxcenter[0] + 2. * ivec[0] * boxsize[0];
                      p.y           = boxcenter[1] + 2. * ivec[1] * boxsize[1];
                      p.z           = boxcenter[2] + 2. * ivec[2] * boxsize[2];
                      p.task        = ThisTask;
                      p.index       = i;
                      p.image_flags = 0;

                      DeRefMesh.DP[DeRefMesh.Ndp] = p;

#ifndef OPTIMIZE_MEMORY_USAGE
                      set_integers_for_point(&DeRefMesh, DeRefMesh.Ndp);
#endif
                      tlast = insert_point(&DeRefMesh, DeRefMesh.Ndp, tlast);

                      DeRefMesh.Ndp++;
                    }
                }

              /* now compute and save the volumes */
              compute_circumcircles(&DeRefMesh);
              double *VolumeNew = (double *)mymalloc_movable(&VolumeNew, "VolumeNew", (DeRefMesh.Ndp + 1) * sizeof(double));
              derefine_refine_compute_volumes(VolumeNew);
              VolumeNew[DeRefMesh.Ndp] = 0;

              /* add the cell that we want to derefine */
              DeRefMesh.DP[DeRefMesh.Ndp] = Mesh.DP[j];
#ifndef OPTIMIZE_MEMORY_USAGE
              set_integers_for_point(&DeRefMesh, DeRefMesh.Ndp);
#endif
              tlast = insert_point(&DeRefMesh, DeRefMesh.Ndp, tlast);
              DeRefMesh.Ndp++;

              for(k = 0; k < DeRefMesh.Ndt; k++)
                DeRefMesh.DTF[k] = 0;
              compute_circumcircles(&DeRefMesh);
              double *VolumeOld = (double *)mymalloc("VolumeOld", DeRefMesh.Ndp * sizeof(double));
              derefine_refine_compute_volumes(VolumeOld);

              /* first, let's establish the surface area sum for this cell */
              double voltot = 0;
              k             = first_ngb[j];
#ifdef TRACER_MC
              int k_last = 0;
#endif
              int abort_derefinement = 0;

              int qq = 0;
              while(k >= 0)
                {
                  if(ngb[k].area > DEREF_SA_FAC * SphP[i].SurfaceArea)
                    {
                      if(VolumeNew[qq] < VolumeOld[qq])
                        {
                          if(VolumeNew[qq] - VolumeOld[qq] < -0.1 * SphP[i].Volume)
                            {
                              printf("DEREFINEMENT ABORT: Task %d, qq=%d, i=%d, SphP[i].Volume=%g\n", ThisTask, qq, i, SphP[i].Volume);

                              double TotVolNew = 0;
                              double TotVolOld = 0;

                              for(int l = 0; l < DeRefMesh.Ndp; l++)
                                {
                                  printf(
                                      "DEREFINEMENT ABORT: Task %d, i=%d, k=%3d, x=%8g, y=%8g, z=%8g, VolNew=%8g, VolOld=%8g, "
                                      "reldiff=%8g, %8g\n",
                                      ThisTask, i, l, DeRefMesh.DP[l].x - boxcenter[0], DeRefMesh.DP[l].y - boxcenter[1],
                                      DeRefMesh.DP[l].z - boxcenter[2], VolumeNew[l], VolumeOld[l],
                                      (VolumeNew[l] - VolumeOld[l]) / VolumeOld[l], (VolumeNew[l] - VolumeOld[l]) / SphP[i].Volume);

                                  if(l < NNgb)
                                    {
                                      TotVolNew += VolumeNew[l];
                                      TotVolOld += VolumeOld[l];
                                    }

                                  if(l == DeRefMesh.Ndp - 1)
                                    TotVolOld += VolumeOld[l];
                                }
                              printf("DEREFINEMENT ABORT: Task %d, TotVolNew=%8g, TotVolOld=%8g, reldiff=%8g\n", ThisTask, TotVolNew,
                                     TotVolOld, (TotVolNew - TotVolOld) / TotVolOld);
                              // terminate("FAILED");
                              abort_derefinement = 1;
                              break;
                            }
                          VolumeNew[qq] = VolumeOld[qq];
                        }
                      voltot += VolumeNew[qq] - VolumeOld[qq];
#ifdef TRACER_MC
                      k_last = k;
#endif
                    }
                  k = ngb[k].next_ngb;
                  qq++;
                }

              if(abort_derefinement)
                {
                  /* abort here, let's try again next timestep, probably that was a cell that has just been refined and the mesh is
                   * highly irregular */
                  myfree(VolumeOld);
                  myfree(VolumeNew);

                  myfree(DeRefMesh.DTF);
                  myfree(DeRefMesh.DTC);
                  DeRefMesh.DTC = NULL;

                  myfree(DeRefMesh.DT);
                  myfree(DeRefMesh.DP - 5);
                  myfree(DeRefMesh.VF);

                  continue;
                }

              /* now, distribute conserved quantities proportional to the gained volume */
              double facsum = 0;
#ifdef TRACER_MC
              MyDouble CellTracerMass = P[i].Mass;
#endif
              k  = first_ngb[j];
              qq = 0;
              while(k >= 0)
                {
                  if(ngb[k].area > DEREF_SA_FAC * SphP[i].SurfaceArea)
                    {
                      int q = ngb[k].index;

                      double fac = (VolumeNew[qq] - VolumeOld[qq]) / voltot;

                      if(fac < 0)
                        terminate("strange");

                      facsum += fac;
#ifdef TRACER_MC
                      int task = Mesh.DP[q].task;
                      int p    = Mesh.DP[q].index;

                      if(task == ThisTask)
                        {
                          if(p >= NumGas)
                            p -= NumGas;
                        }
                      else
                        p = Mesh.DP[q].originalindex;

                      if(k == k_last)
                        {
                          if(fac > 0 && fabs(fac * P[i].Mass / CellTracerMass - 1) > 1.0e-3)
                            warn("fac * P[i].Mass / CellTracermass=%g", fac * P[i].Mass / CellTracerMass);
                          consider_moving_tracers(i, task, p, Mesh.DP[q].ID, 1);
                        }
                      else
                        {
                          if(fac * P[i].Mass / CellTracerMass > 1)
                            {
                              if(fac * P[i].Mass / CellTracerMass - 1 > 1.0e-3)
                                warn("fac * P[i].Mass / CellTracermass=%g", fac * P[i].Mass / CellTracerMass);
                              consider_moving_tracers(i, task, p, Mesh.DP[q].ID, 1);
                            }
                          else
                            consider_moving_tracers(i, task, p, Mesh.DP[q].ID, fac * P[i].Mass / CellTracerMass);
                        }

                      CellTracerMass -= fac * P[i].Mass;
#endif
                      if(Mesh.DP[q].task == ThisTask)
                        {
                          int p = Mesh.DP[q].index;

                          if(p < 0)
                            terminate("p < 0");

                          if(p >= NumGas) /* this is a local ghost point */
                            p -= NumGas;
                          P[p].Mass += fac * P[i].Mass;
                          SphP[p].Momentum[0] += fac * SphP[i].Momentum[0];
                          SphP[p].Momentum[1] += fac * SphP[i].Momentum[1];
                          SphP[p].Momentum[2] += fac * SphP[i].Momentum[2];

#ifdef MHD
                          SphP[p].BConserved[0] += fac * SphP[i].BConserved[0];
                          SphP[p].BConserved[1] += fac * SphP[i].BConserved[1];
                          SphP[p].BConserved[2] += fac * SphP[i].BConserved[2];
#ifdef MHD_CT
                          SphP[p].AConserved[0] += fac * SphP[i].AConserved[0];
                          SphP[p].AConserved[1] += fac * SphP[i].AConserved[1];
                          SphP[p].AConserved[2] += fac * SphP[i].AConserved[2];
#endif
#ifdef MHD_DEDNER
                          SphP[p].PsiConserved += fac * SphP[i].PsiConserved;
#endif
#endif

#ifdef MRT
                          for(int num1 = 0; num1 < MRT_BINS; num1++)
                            {
                              SphP[p].Cons_DensPhot[num1] += fac * SphP[i].Cons_DensPhot[num1];
                              SphP[p].Cons_RT_F[num1][0] += fac * SphP[i].Cons_RT_F[num1][0];
                              SphP[p].Cons_RT_F[num1][1] += fac * SphP[i].Cons_RT_F[num1][1];
                              SphP[p].Cons_RT_F[num1][2] += fac * SphP[i].Cons_RT_F[num1][2];
                            }
#endif

#ifdef OUTPUT_CELL_SPIN
                          SphP[p].Spin[0] += fac * SphP[i].Spin[0];
                          SphP[p].Spin[1] += fac * SphP[i].Spin[1];
                          SphP[p].Spin[2] += fac * SphP[i].Spin[2];
#endif

#ifndef ISOTHERM_EQS
                          SphP[p].Energy += fac * SphP[i].Energy;
#endif
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
                          SphP[p].Entropy += fac * SphP[i].Entropy;
#endif
#ifdef MAXSCALARS
                          for(int s = 0; s < N_Scalar; s++)
                            *(MyFloat *)(((char *)(&SphP[p])) + scalar_elements[s].offset_mass) +=
                                fac * (*(MyFloat *)(((char *)(&SphP[i])) + scalar_elements[s].offset_mass));
#endif
#ifdef GFM_CHEMTAGS
                          for(int k = 0; k < GFM_N_CHEM_TAGS; k++)
                            SphP[p].MassMetalsChemTags[k] += fac * SphP[i].MassMetalsChemTags[k];
#endif
#ifdef TRACER_FIELD
                          SphP[p].ConservedTracer += fac * SphP[i].ConservedTracer;
#endif
#if defined(SMUGGLE_STAR_FEEDBACK) && defined(SMUGGLE_OUTPUT_STELLAR_FEEDBACK)
                          SphP[p].TotEgyFeed += fac * SphP[i].TotEgyFeed;
                          SphP[p].IntEgyFeed += fac * SphP[i].IntEgyFeed;
                          SphP[p].KinEgyFeed += fac * SphP[i].KinEgyFeed;
#endif
#ifdef BH_THERMALFEEDBACK
                          SphP[p].Injected_BH_Energy += fac * SphP[i].Injected_BH_Energy;
#endif
#ifdef GFM_WINDS_LOCAL
                          SphP[p].WindEnergyReceived += fac * SphP[i].WindEnergyReceived;
#endif
#ifdef TGCHEM
                          int m;
                          for(m = 0; m < TGCHEM_NUM_ABUNDANCES; m++)
                            SphP[p].Abund[m] =
                                (SphP[p].Abund[m] * SphP[p].Volume + SphP[i].Abund[m] * volume[q]) / (SphP[p].Volume + volume[q]);
#endif
#ifdef CHIMES
                            /* When cell i is de-refined, we need to distribute its ions and
                             * molecules to its neighbours. We calculate the total number of each
                             * ion/molecule species in cell i (i.e. the ion/molecule density times
                             * the cell volume) and then pass a fraction 'fac' to neighbour p,
                             * according to the gained volume. However, the CHIMES abundances are
                             * recorded as n_ion / n_Htot, so we also need to track the number of
                             * hydrogen nuclei passed from i to p, and then update the abundances
                             * of p accordingly. */
#ifdef GFM_STELLAR_EVOLUTION
                          nHtot_i = SphP[i].MetalsFraction[0] * SphP[i].Density;
                          nHtot_p = SphP[p].MetalsFraction[0] * SphP[p].Density;
#else
                          nHtot_i = HYDROGEN_MASSFRAC * SphP[i].Density;
                          nHtot_p = HYDROGEN_MASSFRAC * SphP[p].Density;
#endif
                          // Convert nH to cgs
                          nHtot_i *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv / PROTONMASS;
                          nHtot_p *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv / PROTONMASS;

                          dNHtot = nHtot_i * SphP[i].Volume * fac;  // number of H nuclei transferred to particle p

                          for(abunIndex = 0; abunIndex < ChimesGlobalVars.totalNumberOfSpecies; abunIndex++)
                            {
                              dN_ion  = SphP[i].ChimesGasVars.abundances[abunIndex] * nHtot_i * SphP[i].Volume * fac;
                              N_ion_p = SphP[p].ChimesGasVars.abundances[abunIndex] * nHtot_p * SphP[p].Volume;
                              SphP[p].ChimesGasVars.abundances[abunIndex] = (N_ion_p + dN_ion) / ((nHtot_p * SphP[p].Volume) + dNHtot);
                            }
#endif  // CHIMES
#ifdef REFINEMENT_SPLIT_CELLS
                          FlagDoNotRefine[p] = 1;
#endif
#ifdef RUNGE_KUTTA_FULL_UPDATE
                          rk_derefinement_add(&SphP[p].rk, &SphP[i].rk, fac);
#endif
                        }
                      else
                        {
                          /* here we have a foreign ghost point */
                          if(Mesh.DP[q].originalindex < 0)
                            terminate("---> task=%d  q=%d j=%d Ndp=%d\n", ThisTask, q, j, Mesh.Ndp);

                          if(Nflux >= MaxNflux)
                            {
                              Mesh.Indi.AllocFacNflux *= ALLOC_INCREASE_FACTOR;
                              MaxNflux = Mesh.Indi.AllocFacNflux;
#ifdef VERBOSE
                              printf("Task=%d: increase memory allocation, MaxNflux=%d Indi.AllocFacNflux=%g\n", ThisTask, MaxNflux,
                                     Mesh.Indi.AllocFacNflux);
#endif
                              FluxList =
                                  (struct flux_list_data *)myrealloc_movable(FluxList, MaxNflux * sizeof(struct flux_list_data));

                              if(Nflux >= MaxNflux)
                                terminate("Nflux >= MaxNflux");
                            }

                          FluxList[Nflux].task  = Mesh.DP[q].task;
                          FluxList[Nflux].index = Mesh.DP[q].originalindex;
                          FluxList[Nflux].dM    = fac * P[i].Mass;
                          FluxList[Nflux].dP[0] = fac * SphP[i].Momentum[0];
                          FluxList[Nflux].dP[1] = fac * SphP[i].Momentum[1];
                          FluxList[Nflux].dP[2] = fac * SphP[i].Momentum[2];
#ifdef MHD
                          FluxList[Nflux].dB[0] = fac * SphP[i].BConserved[0];
                          FluxList[Nflux].dB[1] = fac * SphP[i].BConserved[1];
                          FluxList[Nflux].dB[2] = fac * SphP[i].BConserved[2];
#ifdef MHD_CT
                          FluxList[Nflux].dA[0] = fac * SphP[i].AConserved[0];
                          FluxList[Nflux].dA[1] = fac * SphP[i].AConserved[1];
                          FluxList[Nflux].dA[2] = fac * SphP[i].AConserved[2];
#endif
#ifdef MHD_DEDNER
                          FluxList[Nflux].dPsi = fac * SphP[i].PsiConserved;
#endif
#endif

#ifdef MRT
                          for(int num1 = 0; num1 < MRT_BINS; num1++)
                            {
                              FluxList[Nflux].dDensPhot[num1] = fac * SphP[i].Cons_DensPhot[num1];
                              FluxList[Nflux].dRT_F[num1][0]  = fac * SphP[i].Cons_RT_F[num1][0];
                              FluxList[Nflux].dRT_F[num1][1]  = fac * SphP[i].Cons_RT_F[num1][1];
                              FluxList[Nflux].dRT_F[num1][2]  = fac * SphP[i].Cons_RT_F[num1][2];
                            }
#endif

#ifdef OUTPUT_CELL_SPIN
                          FluxList[Nflux].dSpin[0] = fac * SphP[i].Spin[0];
                          FluxList[Nflux].dSpin[1] = fac * SphP[i].Spin[1];
                          FluxList[Nflux].dSpin[2] = fac * SphP[i].Spin[2];
#endif
#ifndef ISOTHERM_EQS
                          FluxList[Nflux].dEnergy = fac * SphP[i].Energy;
#endif
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
                          FluxList[Nflux].dEntropy = fac * SphP[i].Entropy;
#endif
#ifdef MAXSCALARS
                          for(int s = 0; s < N_Scalar; s++)
                            FluxList[Nflux].dConservedScalars[s] =
                                fac * (*(MyFloat *)(((char *)(&SphP[i])) + scalar_elements[s].offset_mass));
#endif
#ifdef GFM_CHEMTAGS
                          for(int k = 0; k < GFM_N_CHEM_TAGS; k++)
                            FluxList[Nflux].dMassMetalsChemTags[k] = fac * SphP[i].MassMetalsChemTags[k];
#endif

#ifdef TRACER_FIELD
                          FluxList[Nflux].dConservedTracer = fac * SphP[i].ConservedTracer;
#endif
#if defined(SMUGGLE_STAR_FEEDBACK) && defined(SMUGGLE_OUTPUT_STELLAR_FEEDBACK)
                          FluxList[Nflux].dTotEgyFeed = fac * SphP[i].TotEgyFeed;
                          FluxList[Nflux].dIntEgyFeed = fac * SphP[i].IntEgyFeed;
                          FluxList[Nflux].dKinEgyFeed = fac * SphP[i].KinEgyFeed;
#endif
#ifdef BH_THERMALFEEDBACK
                          FluxList[Nflux].dInjected_BH_Energy = fac * SphP[i].Injected_BH_Energy;
#endif
#ifdef GFM_WINDS_LOCAL
                          FluxList[Nflux].dWindEnergyReceived = fac * SphP[i].WindEnergyReceived;
#endif
#ifdef TGCHEM
                          int m;
                          FluxList[Nflux].dVolume = volume[q];
                          for(m = 0; m < TGCHEM_NUM_ABUNDANCES; m++)
                            FluxList[Nflux].dAbund[m] = SphP[i].Abund[m];
#endif
#ifdef CHIMES
#ifdef GFM_STELLAR_EVOLUTION
                          nHtot_i = SphP[i].MetalsFraction[0] * SphP[i].Density;
#else
                          nHtot_i = HYDROGEN_MASSFRAC * SphP[i].Density;
#endif
                          // Convert to nH, in cgs
                          nHtot_i *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv / PROTONMASS;

                          FluxList[Nflux].dNHtot = nHtot_i * SphP[i].Volume * fac;  // number of H nuclei transferred to particle p
                          for(abunIndex = 0; abunIndex < ChimesGlobalVars.totalNumberOfSpecies; abunIndex++)
                            FluxList[Nflux].dN_ion[abunIndex] =
                                SphP[i].ChimesGasVars.abundances[abunIndex] * nHtot_i * SphP[i].Volume * fac;
#endif  // CHIMES
#ifdef RUNGE_KUTTA_FULL_UPDATE
                          rk_derefinement_set(&FluxList[Nflux].rk, &SphP[i].rk, fac);
#endif
                          Nflux++;
                        }
                    }

                  k = ngb[k].next_ngb;
                  qq++;
                }

              myfree(VolumeOld);
              myfree(VolumeNew);

              myfree(DeRefMesh.DTF);
              myfree(DeRefMesh.DTC);
              DeRefMesh.DTC = NULL;

              myfree(DeRefMesh.DT);
              myfree(DeRefMesh.DP - 5);
              myfree(DeRefMesh.VF);

              /* we set the dissolved cell to zero mass and zero ID. It will be eliminated from the list
               * of cells in the next domain decomposition
               */
              P[i].Mass   = 0;
              P[i].ID     = 0;
              P[i].Vel[0] = 0;
              P[i].Vel[1] = 0;
              P[i].Vel[2] = 0;

              SphP[i].VelVertex[0] = 0;
              SphP[i].VelVertex[1] = 0;
              SphP[i].VelVertex[2] = 0;

              timebin_remove_particle(&TimeBinsHydro, idx, P[i].TimeBinHydro);

#ifdef VORONOI_DYNAMIC_UPDATE
              voronoi_remove_connection(i);
#endif
            }
        }

        /* now let's apply the flux-list */
#ifdef TRACER_MC
      finish_MC_tracer();
#endif

      derefine_apply_flux_list();
      myfree(FluxList);

      myfree(ngb);
      myfree(first_ngb);
    }

  myfree(FlagExch);
  myfree(deref_SphP);

  /* remove removed cells from list of active gravity cells */
  timebin_cleanup_list_of_active_particles(&TimeBinsGravity);

  TIMER_STOP(CPU_DEREFINE);

  return countall;
}

void voronoi_refinement_check_enclosing_box(point *p, double *minpos, double *maxpos)
{
  if(p->x > maxpos[0])
    maxpos[0] = p->x;
  if(p->x < minpos[0])
    minpos[0] = p->x;

  if(p->y > maxpos[1])
    maxpos[1] = p->y;
  if(p->y < minpos[1])
    minpos[1] = p->y;

  if(p->z > maxpos[2])
    maxpos[2] = p->z;
  if(p->z < minpos[2])
    minpos[2] = p->z;
}

static void derefine_apply_probe_list(void)
{
  int i, j, p, nimport, ngrp, recvTask;

  /* now exchange the probe-list and apply it where needed */

  mysort(ProbeList, Nflux, sizeof(struct probe_list_data), derefine_probe_list_data_compare_task);

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(i = 0; i < Nflux; i++)
    Send_count[ProbeList[i].task]++;

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

  struct probe_list_data *ProbeListGet = (struct probe_list_data *)mymalloc("ProbeListGet", nimport * sizeof(struct probe_list_data));

  /* exchange particle data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&ProbeList[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct probe_list_data), MPI_BYTE,
                           recvTask, TAG_DENS_A, &ProbeListGet[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct probe_list_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  /* apply the probes */

  for(i = 0; i < nimport; i++)
    {
      p = ProbeListGet[i].index;

      if(TimeBinSynchronized[P[p].TimeBinHydro])
        if(deref_SphP[p].Flag == 2)
          ProbeListGet[i].flag = 1;
    }

  /* send results back */

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&ProbeListGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct probe_list_data), MPI_BYTE,
                           recvTask, TAG_DENS_A, &ProbeList[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct probe_list_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(ProbeListGet);
}

static void derefine_apply_flux_list(void)
{
  int i, j, p, nimport, ngrp, recvTask;

#ifdef CHIMES
  double nHtot_p, N_ion_p;
  int abunIndex;
#endif

  /* now exchange the flux-list and apply it when needed */

  mysort(FluxList, Nflux, sizeof(struct flux_list_data), derefine_flux_list_data_compare);

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(i = 0; i < Nflux; i++)
    Send_count[FluxList[i].task]++;

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

  struct flux_list_data *FluxListGet = (struct flux_list_data *)mymalloc("FluxListGet", nimport * sizeof(struct flux_list_data));

  /* exchange particle data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&FluxList[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct flux_list_data), MPI_BYTE, recvTask,
                           TAG_DENS_A, &FluxListGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct flux_list_data),
                           MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  /* apply the fluxes */

  for(i = 0; i < nimport; i++)
    {
      p = FluxListGet[i].index;

      if(P[p].ID == 0)
        {
#ifndef LONGIDS
          terminate("On task=%d flux to ID=%" MYIDTYPE_PRI ", but this is already deleted (index p=%d)\n", ThisTask, P[p].ID, p);
#else
          terminate("On task=%d flux to ID=%llu, but this is already deleted (index p=%d)\n", ThisTask, P[p].ID, p);
#endif
        }

      P[p].Mass += FluxListGet[i].dM;
      SphP[p].Momentum[0] += FluxListGet[i].dP[0];
      SphP[p].Momentum[1] += FluxListGet[i].dP[1];
      SphP[p].Momentum[2] += FluxListGet[i].dP[2];
#ifdef MHD
      SphP[p].BConserved[0] += FluxListGet[i].dB[0];
      SphP[p].BConserved[1] += FluxListGet[i].dB[1];
      SphP[p].BConserved[2] += FluxListGet[i].dB[2];
#ifdef MHD_CT
      SphP[p].AConserved[0] += FluxListGet[i].dA[0];
      SphP[p].AConserved[1] += FluxListGet[i].dA[1];
      SphP[p].AConserved[2] += FluxListGet[i].dA[2];
#endif
#ifdef MHD_DEDNER
      SphP[p].PsiConserved += FluxListGet[i].dPsi;
#endif
#endif

#ifdef MRT
      for(int num1 = 0; num1 < MRT_BINS; num1++)
        {
          SphP[p].Cons_DensPhot[num1] += FluxListGet[i].dDensPhot[num1];
          SphP[p].Cons_RT_F[num1][0] += FluxListGet[i].dRT_F[num1][0];
          SphP[p].Cons_RT_F[num1][1] += FluxListGet[i].dRT_F[num1][1];
          SphP[p].Cons_RT_F[num1][2] += FluxListGet[i].dRT_F[num1][2];
        }
#endif

#ifdef OUTPUT_CELL_SPIN
      SphP[p].Spin[0] += FluxListGet[i].dSpin[0];
      SphP[p].Spin[1] += FluxListGet[i].dSpin[1];
      SphP[p].Spin[2] += FluxListGet[i].dSpin[2];
#endif

#ifdef MAXSCALARS
      int k;
      for(k = 0; k < N_Scalar; k++)
        *(MyFloat *)(((char *)(&SphP[p])) + scalar_elements[k].offset_mass) += FluxListGet[i].dConservedScalars[k];
#endif

#ifdef GFM_CHEMTAGS
      for(int k = 0; k < GFM_N_CHEM_TAGS; k++)
        SphP[p].MassMetalsChemTags[k] += FluxListGet[i].dMassMetalsChemTags[k];
#endif

#ifndef ISOTHERM_EQS
      SphP[p].Energy += FluxListGet[i].dEnergy;
#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
      SphP[p].Entropy += FluxListGet[i].dEntropy;
#endif

#ifdef TRACER_FIELD
      SphP[p].ConservedTracer += FluxListGet[i].dConservedTracer;
#endif

#if defined(SMUGGLE_STAR_FEEDBACK) && defined(SMUGGLE_OUTPUT_STELLAR_FEEDBACK)
      SphP[p].TotEgyFeed += FluxListGet[i].dTotEgyFeed;
      SphP[p].IntEgyFeed += FluxListGet[i].dIntEgyFeed;
      SphP[p].KinEgyFeed += FluxListGet[i].dKinEgyFeed;
#endif

#ifdef BH_THERMALFEEDBACK
      SphP[p].Injected_BH_Energy += FluxListGet[i].dInjected_BH_Energy;
#endif
#ifdef GFM_WINDS_LOCAL
      SphP[p].WindEnergyReceived += FluxListGet[i].dWindEnergyReceived;
#endif
#ifdef TGCHEM
      int m;
      for(m = 0; m < TGCHEM_NUM_ABUNDANCES; m++)
        SphP[p].Abund[m] = (SphP[p].Abund[m] * SphP[p].Volume + FluxListGet[i].dAbund[m] * FluxListGet[i].dVolume) /
                           (SphP[p].Volume + FluxListGet[i].dVolume);
#endif
#ifdef CHIMES
#ifdef GFM_STELLAR_EVOLUTION
      nHtot_p = SphP[p].MetalsFraction[0] * SphP[p].Density;
#else
      nHtot_p = HYDROGEN_MASSFRAC * SphP[p].Density;
#endif
      // Convert to nH, in cgs
      nHtot_p *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv / PROTONMASS;

      for(abunIndex = 0; abunIndex < ChimesGlobalVars.totalNumberOfSpecies; abunIndex++)
        {
          N_ion_p = SphP[p].ChimesGasVars.abundances[abunIndex] * nHtot_p * SphP[p].Volume;
          SphP[p].ChimesGasVars.abundances[abunIndex] =
              (N_ion_p + FluxListGet[i].dN_ion[abunIndex]) / ((nHtot_p * SphP[p].Volume) + FluxListGet[i].dNHtot);
        }
#endif  // CHIMES
#ifdef REFINEMENT_SPLIT_CELLS
      FlagDoNotRefine[p] = 1;
#endif
#ifdef RUNGE_KUTTA_FULL_UPDATE
      rk_derefinement_add(&SphP[p].rk, &FluxListGet[i].rk, 1.0);
#endif
    }

  myfree(FluxListGet);
}

static int derefine_flux_list_data_compare(const void *a, const void *b)
{
  if(((struct flux_list_data *)a)->task < (((struct flux_list_data *)b)->task))
    return -1;

  if(((struct flux_list_data *)a)->task > (((struct flux_list_data *)b)->task))
    return +1;

  return 0;
}

static int derefine_probe_list_data_compare_task(const void *a, const void *b)
{
  if(((struct probe_list_data *)a)->task < (((struct probe_list_data *)b)->task))
    return -1;

  if(((struct probe_list_data *)a)->task > (((struct probe_list_data *)b)->task))
    return +1;

  return 0;
}

static int derefine_compare_seq_DP_ID(const void *a, const void *b)
{
  if(((struct seq_delaunay_data *)a)->rnd < (((struct seq_delaunay_data *)b)->rnd))
    return -1;

  if(((struct seq_delaunay_data *)a)->rnd > (((struct seq_delaunay_data *)b)->rnd))
    return +1;

  if(((struct seq_delaunay_data *)a)->ID < (((struct seq_delaunay_data *)b)->ID))
    return -1;

  if(((struct seq_delaunay_data *)a)->ID > (((struct seq_delaunay_data *)b)->ID))
    return +1;

  if(((struct seq_delaunay_data *)a)->index < (((struct seq_delaunay_data *)b)->index))
    return -1;

  if(((struct seq_delaunay_data *)a)->index > (((struct seq_delaunay_data *)b)->index))
    return +1;

  if(((struct seq_delaunay_data *)a)->rank < (((struct seq_delaunay_data *)b)->rank))
    return -1;

  if(((struct seq_delaunay_data *)a)->rank > (((struct seq_delaunay_data *)b)->rank))
    return +1;

  return 0;
}

static void derefine_exchange_flag(void)
{
  int listp;
  int i, j, p, task, off;
  int ngrp, recvTask, place;

  struct exchange_data
  {
    int Flag;
    MyIDType ID;
  } * tmpExch, *tmpRecv;

  tmpExch = (struct exchange_data *)mymalloc("tmpExch", Mesh_nexport * sizeof(struct exchange_data));

  /* prepare data for export */
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

              tmpExch[off].Flag = 0;
              tmpExch[off].ID   = P[place].ID;

              if(P[place].Type == 0)
                if(TimeBinSynchronized[P[place].TimeBinHydro])
                  if(!(P[place].Mass == 0 && P[place].ID == 0))
                    tmpExch[off].Flag = deref_SphP[place].Flag;
            }
          listp = ListExports[listp].nextexport;
        }
    }

  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              tmpRecv = (struct exchange_data *)mymalloc("tmpRecv", Mesh_Recv_count[recvTask] * sizeof(struct exchange_data));

              /* get the values */
              MPI_Sendrecv(&tmpExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask] * sizeof(struct exchange_data), MPI_BYTE,
                           recvTask, TAG_DENS_A, tmpRecv, Mesh_Recv_count[recvTask] * sizeof(struct exchange_data), MPI_BYTE, recvTask,
                           TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              for(i = 0; i < Mesh_Recv_count[recvTask]; i++)
                {
                  if(Mesh_Recv_offset[recvTask] + i >= Mesh_nimport)
                    terminate("number of imported mesh points grater than Mesh_nimport");
                  FlagExch[Mesh_Recv_offset[recvTask] + i].Flag = tmpRecv[i].Flag;
                  FlagExch[Mesh_Recv_offset[recvTask] + i].ID   = tmpRecv[i].ID;
                }

              myfree(tmpRecv);
            }
        }
    }

  myfree(tmpExch);
}

#endif

#endif /* VORONOI */
