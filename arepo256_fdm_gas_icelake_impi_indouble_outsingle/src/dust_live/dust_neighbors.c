/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/dust_neighbors.c
 * \date        MM/YYYY
 * \author      Ryan McKinnon
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef DUST_LIVE
#if(defined(DL_GRAIN_BINS) && (defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION))) || \
    (defined(DL_DRAG_BACKREACTION))

#if defined(HIERARCHICAL_GRAVITY) || defined(ALLOW_DIRECT_SUMMATION)
#error \
    "Use of options DL_SNE_DESTRUCTION and DL_SHATTERING and DL_COAGULATION and DL_DRAG_BACKREACTION does not work with HIERARCHICAL_GRAVITY or ALLOW_DIRECT_SUMMATION since we need a full gravity tree for neighbor searches at each time!"
#endif

typedef struct
{
  /* your fields go here */
  MyDouble Pos[3];
  MyFloat Hsml;
  int Firstnode;
} data_in;

static data_in *DataGet;

typedef struct
{
  /* your fields go here */
  MyFloat DustDensity;
  MyFloat DustNumNgb;
  MyFloat EnclosedMass;
#ifdef DL_DEREFINEMENT
  MyFloat ClosestDustR;
  MyIDType ClosestDustID;
  int ClosestDustTask;
  int ClosestDustIndex;
  int ClosestDustTaskTree;
  int ClosestDustIndexTree;
  int ClosestDustHasP;
  int HighMassNeighbor;
#endif
} data_out;

static data_out *DataResult;

static void particle2in(data_in *in, int i, int firstnode)
{
  int k;

  int idx = TargetList[i];

  if(idx < NumPart)
    {
      for(k = 0; k < 3; k++)
        {
          in->Pos[k] = P[idx].Pos[k];
        }
    }
  else
    {
      idx -= Tree_ImportedNodeOffset;
      for(k = 0; k < 3; k++)
        {
          in->Pos[k] = Tree_Points[idx].Pos[k];
        }
    }
  in->Hsml      = PShatter[i].DustHsml;
  in->Firstnode = firstnode;
}

static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
      PShatter[i].DustDensity  = out->DustDensity;
      PShatter[i].DustNumNgb   = out->DustNumNgb;
      PShatter[i].EnclosedMass = out->EnclosedMass;
#ifdef DL_DEREFINEMENT
      PShatter[i].ClosestDustR         = out->ClosestDustR;
      PShatter[i].ClosestDustID        = out->ClosestDustID;
      PShatter[i].ClosestDustTask      = out->ClosestDustTask;
      PShatter[i].ClosestDustIndex     = out->ClosestDustIndex;
      PShatter[i].ClosestDustTaskTree  = out->ClosestDustTaskTree;
      PShatter[i].ClosestDustIndexTree = out->ClosestDustIndexTree;
      PShatter[i].ClosestDustHasP      = out->ClosestDustHasP;
      PShatter[i].HighMassNeighbor     = out->HighMassNeighbor;
#endif
    }
  else /* combine */
    {
      PShatter[i].DustDensity += out->DustDensity;
      PShatter[i].DustNumNgb += out->DustNumNgb;
      PShatter[i].EnclosedMass += out->EnclosedMass;
#ifdef DL_DEREFINEMENT
      /* Only if we found a non-local dust even closer */
      if(out->ClosestDustR < PShatter[i].ClosestDustR)
        {
          PShatter[i].ClosestDustR         = out->ClosestDustR;
          PShatter[i].ClosestDustID        = out->ClosestDustID;
          PShatter[i].ClosestDustTask      = out->ClosestDustTask;
          PShatter[i].ClosestDustIndex     = out->ClosestDustIndex;
          PShatter[i].ClosestDustTaskTree  = out->ClosestDustTaskTree;
          PShatter[i].ClosestDustIndexTree = out->ClosestDustIndexTree;
          PShatter[i].ClosestDustHasP      = out->ClosestDustHasP;
          PShatter[i].HighMassNeighbor     = out->HighMassNeighbor;
        }
#endif
    }
}

#include "../generic_comm_helpers2.h"

static unsigned char *Todo;

static void kernel_local(void)
{
  int idx;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif

#pragma omp parallel private(idx)
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
              if(generic_polling_primary(count, NforcesDust))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        idx = NextParticle++;

        if(idx >= NforcesDust)
          break;

        if(Todo[idx] > 0)
          dust_findHsml_evaluate(idx, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}

static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
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

        dust_findHsml_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void dust_findHsml(void)
{
  int idx, i;
  int iter = 0, npleft;
  MyDouble *Left, *Right;
  double t0, t1;
  long long ntot, npartall;

  sumup_large_ints(1, &NforcesDust, &npartall);
  if(npartall == 0)
    return;

  mpi_printf("DUST_LIVE: Finding dust-dust hsml values\n");

  /* NforcesDust is set in begin_dust_search() */
  Left  = (MyDouble *)mymalloc("Left", NforcesDust * sizeof(MyDouble));
  Right = (MyDouble *)mymalloc("Right", NforcesDust * sizeof(MyDouble));
  Todo  = (unsigned char *)mymalloc("Todo", NforcesDust * sizeof(unsigned char));

  int nforces = 0;
  for(idx = 0; idx < TimeBinsDust.NActiveParticles; idx++)
    {
      i = TimeBinsDust.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if((P[i].Type == DUST_LIVE) && (Tree_Task_list[i] == ThisTask) && (P[i].Mass > 0.0))
        {
          Left[nforces]                            = 0;
          Right[nforces]                           = 0;
          Todo[nforces]                            = 1;
          PShatter[nforces].DustHsml               = DTP(i).DustHsml;
          PShatter[nforces].MinNumNgbDeviationDust = DTP(i).MinNumNgbDeviationDust;
          PShatter[nforces].MaxNumNgbDeviationDust = DTP(i).MaxNumNgbDeviationDust;
          nforces++;
        }
    }

  for(i = 0; i < Tree_NumPartImported; i++)
#ifndef HIERARCHICAL_GRAVITY
    if(Tree_Points[i].ActiveFlag)
#endif
      if((Tree_Points[i].Type == DUST_LIVE) && (Tree_Points[i].Mass > 0.0))
        {
          Left[nforces]                            = 0;
          Right[nforces]                           = 0;
          Todo[nforces]                            = 1;
          PShatter[nforces].DustHsml               = Tree_Points[i].DustHsml;
          PShatter[nforces].MinNumNgbDeviationDust = Tree_Points[i].MinNumNgbDeviationDust;
          PShatter[nforces].MaxNumNgbDeviationDust = Tree_Points[i].MaxNumNgbDeviationDust;
          nforces++;
        }

  for(i = 0; i < nforces; i++)
    {
      /* Reset dust smoothing length if previously there were not enough
       * dust particles. */
      if(PShatter[i].DustHsml >= 0.1 * All.BoxSize)
        PShatter[i].DustHsml = get_default_softening_of_particletype(DUST_LIVE);
      /* If a previous iteration had to increase the maximum allowable
       * enclosed mass to find a high-mass neighbor, relax it back to the
       * default. */
      if(PShatter[i].MaxNumNgbDeviationDust > All.MaxNumNgbDeviationDust)
        PShatter[i].MaxNumNgbDeviationDust = fmax(PShatter[i].MaxNumNgbDeviationDust / 2.0, All.MaxNumNgbDeviationDust);
    }

  generic_set_MaxNexport();

  do
    {
      t0 = second();

      generic_comm_pattern(NforcesDust, kernel_local, kernel_imported);

      /* do final operations on results */
      for(i = 0, npleft = 0; i < NforcesDust; i++)
        {
          if(Todo[i])
            {
              /* Want to find a dust-dust smoothing length that encloses a
               * desired amount of dust mass and, if derefinement is active,
               * encloses at least one high-mass neighbor. */
              int hsml_too_big = PShatter[i].DustHsml >= All.BoxSize / 10.0;
#ifdef DL_DEREFINEMENT
              /* When DL_DEREFINEMENT==0, we want to derefine into particles
               * with mass above the derefinement limit, and it's possible one
               * of these does not exist inside the kernel.  When
               * DL_DEREFINEMENT==1, we derefine into a particle of any mass.
               * Since HighMassNeighbor indicates whether we have a neighbor
               * whose mass is large enough that we can derefine into that
               * particle, if we set
               *
               * int high_mass_neighbor = PShatter[i].HighMassNeighbor;
               *
               * then for DL_DEREFINEMENT==0, we force the kernel to always
               * include a particle that we can derefine into.  (For
               * DL_DEREFINEMENT==1, this value is always 1, since there's
               * always a neighbor in the kernel.)  However, this smoothing
               * length logic is complicated, and the line below skips that
               * constraint. */
              int high_mass_neighbor = 1;
#else
              /* If we're not derefining, we don't need a neighbor with mass
               * above the derefinement limit, so ignore this constraint. */
              int high_mass_neighbor = 1;
#endif

#ifdef DL_PRODUCTION
              /* Unlike gas density estimates, here we use a enclosed mass
               * constraint (instead of number of neighbors) since dust
               * particle masses can vary a bit. */
              double min_enclosed_mass =
                  (All.DesNumNgbDust - PShatter[i].MinNumNgbDeviationDust) * All.DustTargetFrac * All.TargetGasMass;
              double max_enclosed_mass =
                  (All.DesNumNgbDust + PShatter[i].MaxNumNgbDeviationDust) * All.DustTargetFrac * All.TargetGasMass;
              int too_small_enclosed = PShatter[i].EnclosedMass < min_enclosed_mass;
              int too_big_enclosed   = PShatter[i].EnclosedMass > max_enclosed_mass;
              int desired_enclosed =
                  (PShatter[i].EnclosedMass >= min_enclosed_mass) && (PShatter[i].EnclosedMass <= max_enclosed_mass);
#else
              /* If no target dust mass is defined, fall back to an ordinary
               * number of neighbors criterion. */
              double min_enclosed_ngbs = All.DesNumNgbDust - PShatter[i].MinNumNgbDeviationDust;
              double max_enclosed_ngbs = All.DesNumNgbDust + PShatter[i].MaxNumNgbDeviationDust;
              int too_small_enclosed   = PShatter[i].DustNumNgb < min_enclosed_ngbs;
              int too_big_enclosed     = PShatter[i].DustNumNgb > max_enclosed_ngbs;
              int desired_enclosed     = (PShatter[i].DustNumNgb >= min_enclosed_ngbs) && (PShatter[i].DustNumNgb <= max_enclosed_ngbs);
#endif

              int should_hsml_stop = (hsml_too_big) || (desired_enclosed && high_mass_neighbor);
              if(!should_hsml_stop)
                {
                  /* need to redo this particle */
                  npleft++;

                  /* If we don't have enough mass or a high-mass neighbor,
                   * we'll have to increase the smoothing length by setting a
                   * new lower bound. */
                  if(too_small_enclosed)
                    Left[i] = fmax(PShatter[i].DustHsml, Left[i]);
                  else if(desired_enclosed && !high_mass_neighbor)
                    Left[i] = fmax(PShatter[i].DustHsml, Left[i]);
                  else if(too_big_enclosed && !high_mass_neighbor)
                    {
                      Left[i] = fmax(PShatter[i].DustHsml, Left[i]);
#ifdef DL_PRODUCTION
                      /* If we have too much enclosed mass but not yet a
                       * high-mass neighbor, we'll need to simultaneously set a
                       * new lower bound (to increase the smoothing length) and
                       * increase the allowed amount of enclosed mass to at
                       * least the current value, plus some buffer. */
                      PShatter[i].MaxNumNgbDeviationDust =
                          (PShatter[i].EnclosedMass / (All.DustTargetFrac * All.TargetGasMass)) - All.DesNumNgbDust;
#endif
                      PShatter[i].MaxNumNgbDeviationDust *= 1.26;
                    }
                  else /* too_big_enclosed && high_mass_neighbor */
                    {
                      /* Since we have a high-mass neighbor and too much
                       * enclosed mass, try to shrink the smoothing length by
                       * setting a new upper bound. */
                      if(Right[i] != 0)
                        {
                          if(PShatter[i].DustHsml < Right[i])
                            Right[i] = PShatter[i].DustHsml;
                        }
                      else
                        Right[i] = PShatter[i].DustHsml;
                    }

                  if(iter >= MAXITER - 10)
                    {
                      printf("DUST_LIVE: i=%d task=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g EncMass=%g HighMassNeighbor=%d\n",
                             i, ThisTask, PShatter[i].DustHsml, Left[i], Right[i], (double)PShatter[i].DustNumNgb, Right[i] - Left[i],
                             PShatter[i].EnclosedMass, high_mass_neighbor);
                      fflush(stdout);
                    }

                  if(Right[i] > 0 && Left[i] > 0)
                    PShatter[i].DustHsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
                  else
                    {
                      if(Right[i] == 0 && Left[i] == 0)
                        terminate("BAD"); /* can't occur */

                      if(Right[i] == 0 && Left[i] > 0)
                        PShatter[i].DustHsml *= 1.26;

                      if(Right[i] > 0 && Left[i] == 0)
                        PShatter[i].DustHsml /= 1.26;
                    }
                }
              else
                {
                  Todo[i] = 0;
                }
            }
        }

      sumup_large_ints(1, &npleft, &ntot);

      t1 = second();

      if(ntot > 0)
        {
          iter++;

          if(iter > 0 && ThisTask == 0)
            {
              printf("DUST_LIVE: ngb iteration %d: need to repeat for %llu particles. (took %g sec)\n", iter, (unsigned long long)ntot,
                     timediff(t0, t1));
              fflush(stdout);
            }

          if(iter > MAXITER)
            {
              printf("DUST_LIVE: failed to converge in dust-dust neighbour iteration\n");
              fflush(stdout);
              terminate("BAD");
            }
        }
    }
  while(ntot > 0);

  myfree(Todo);
  myfree(Right);
  myfree(Left);

  mpi_printf("DUST_LIVE: Done with dust-dust hsml values\n");
}

void dust_findHsml_evaluate(int target, int mode, int threadid)
{
  int numnodes, *firstnode;
  data_in local, *in;
  data_out out;
  int k;
  int no;
  double wk, h, h2, h3, hinv, hinv3;
  double r, r2, u;
  MyDouble *pos;
  MyDouble dx, dy, dz;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      in = &local;

      numnodes  = 1;
      firstnode = NULL;
    }
  else
    {
      in = &DataGet[target];
      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  memset(&out, 0, sizeof(data_out));
#ifdef DL_DEREFINEMENT
  out.ClosestDustR = MAX_REAL_NUMBER;
#endif

  pos   = in->Pos;
  h     = in->Hsml;
  h2    = h * h;
  h3    = h2 * h;
  hinv  = 1. / h;
  hinv3 = hinv * hinv * hinv;

  for(k = 0; k < numnodes; k++)
    {
      if(mode == MODE_LOCAL_PARTICLES)
        {
          no = Tree_MaxPart; /* root node */
        }
      else
        {
          no = firstnode[k];
          no = Nodes[no].u.d.nextnode; /* open it */
        }

      while(no >= 0)
        {
          if(no < Tree_MaxPart) /* single particle */
            {
              dx = GRAVITY_NEAREST_X(Tree_Pos_list[3 * no + 0] - pos[0]);
              dy = GRAVITY_NEAREST_Y(Tree_Pos_list[3 * no + 1] - pos[1]);
              dz = GRAVITY_NEAREST_Z(Tree_Pos_list[3 * no + 2] - pos[2]);

              r2 = dx * dx + dy * dy + dz * dz;

              if((r2 < h2) && (P[no].Type == DUST_LIVE) && (P[no].Mass > 0.0))
                {
                  // hinv = 1. / h;
                  // hinv3 = hinv * hinv * hinv;

                  r = sqrt(r2);
                  u = r * hinv;

                  if(u < 0.5)
                    {
                      wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                    }
                  else
                    {
                      wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
                    }

                  out.DustDensity += P[no].Mass * wk;
                  out.DustNumNgb += NORM_COEFF * wk * h3;
                  out.EnclosedMass += P[no].Mass;
#ifdef DL_DEREFINEMENT
#if(DL_DEREFINEMENT == 0)
                  /* Only search for dust particles that will not be derefined */
                  if((r < out.ClosestDustR) && (r > 0.0) && (P[no].Mass >= All.DustMinFrac * All.TargetGasMass))
#elif(DL_DEREFINEMENT == 1)
                  /* Search for dust particles with any mass */
                  if((r < out.ClosestDustR) && (r > 0.0))
#endif
                    {
                      out.ClosestDustR         = r;
                      out.ClosestDustID        = P[no].ID;
                      out.ClosestDustTask      = Tree_Task_list[no];
                      out.ClosestDustIndex     = no;
                      out.ClosestDustTaskTree  = ThisTask;
                      out.ClosestDustIndexTree = no;
                      out.ClosestDustHasP      = 1;
                      out.HighMassNeighbor     = 1;
                    }
#endif
                }
              no = Nextnode[no];
            }
          else if(no < Tree_MaxPart + Tree_MaxNodes) /* internal node */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                {
                  if(no <
                     Tree_FirstNonTopLevelNode) /* we reached a top-level node again, which means that we are done with the branch */
                    break;
                }

              struct NODE *current = &Nodes[no];
              no                   = current->u.d.sibling; /* in case the node can be discarded */

              double dist = h + 0.5 * current->len;
              dx          = NGB_PERIODIC_LONG_X(current->center[0] - pos[0]);
              if(dx > dist)
                continue;
              dy = NGB_PERIODIC_LONG_Y(current->center[1] - pos[1]);
              if(dy > dist)
                continue;
              dz = NGB_PERIODIC_LONG_Z(current->center[2] - pos[2]);
              if(dz > dist)
                continue;

              /* now test against the minimal sphere enclosing everything */
              dist += FACT1 * current->len;
              if(dx * dx + dy * dy + dz * dz > dist * dist)
                continue;

              no = current->u.d.nextnode; /* ok, we need to open the node */
            }

          else if(no >= Tree_ImportedNodeOffset) /* point from imported nodelist */
            {
              int n = no - Tree_ImportedNodeOffset;

              dx = GRAVITY_NEAREST_X(Tree_Points[n].Pos[0] - pos[0]);
              dy = GRAVITY_NEAREST_Y(Tree_Points[n].Pos[1] - pos[1]);
              dz = GRAVITY_NEAREST_Z(Tree_Points[n].Pos[2] - pos[2]);

              r2 = dx * dx + dy * dy + dz * dz;

              if((r2 < h2) && (Tree_Points[n].Type == DUST_LIVE) && (Tree_Points[n].Mass > 0.0))
                {
                  // hinv = 1. / h;
                  // hinv3 = hinv * hinv * hinv;

                  r = sqrt(r2);
                  u = r * hinv;

                  if(u < 0.5)
                    {
                      wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                    }
                  else
                    {
                      wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
                    }

                  out.DustDensity += Tree_Points[n].Mass * wk;
                  out.DustNumNgb += NORM_COEFF * wk * h3;
                  out.EnclosedMass += Tree_Points[n].Mass;
#ifdef DL_DEREFINEMENT
#if(DL_DEREFINEMENT == 0)
                  /* Only search for dust particles that will not be derefined */
                  if((r < out.ClosestDustR) && (r > 0.0) && (Tree_Points[n].Mass >= All.DustMinFrac * All.TargetGasMass))
#elif(DL_DEREFINEMENT == 1)
                  /* Search for dust particles with any mass */
                  if((r < out.ClosestDustR) && (r > 0.0))
#endif
                    {
                      out.ClosestDustR         = r;
                      out.ClosestDustID        = Tree_Points[n].DustID;
                      out.ClosestDustTask      = Tree_Points[n].origTask;
                      out.ClosestDustIndex     = Tree_Points[n].index;
                      out.ClosestDustTaskTree  = ThisTask;
                      out.ClosestDustIndexTree = n;
                      out.ClosestDustHasP      = 0;
                      out.HighMassNeighbor     = 1;
                    }
#endif
                }

              no = Nextnode[no - Tree_MaxNodes];
            }
          else /* pseudo particle */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                terminate("mode == MODE_IMPORTED_PARTICLES");

              if(target >= 0)
                tree_treefind_export_node_threads(no, target, threadid);

              no = Nextnode[no - Tree_MaxNodes];
              continue;
            }
        }
    }

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;
}

#endif
#endif
