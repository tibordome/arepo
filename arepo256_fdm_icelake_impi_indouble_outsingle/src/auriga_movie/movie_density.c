/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/movie_auriga/movie_density.c
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

#include <hdf5.h>

#include "../allvars.h"
#include "../proto.h"
#include "movie.h"

#ifdef AURIGA_MOVIE

#define SEARCHFAC 1.26

static void auriga_movie_density_initialize_hsml(void);
static void auriga_movie_density_compute(void);

static int auriga_movie_density_evaluate(int target, int mode, int threadid);

static MyFloat *NumNgb, *Rho, *Dhsml;

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Hsml;
  int Type;

  int Firstnode;
} data_in;

static data_in *DataGet;

/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in *in, int i, int firstnode)
{
#ifdef CELL_CENTER_GRAVITY
  if(P[i].Type == 0)
    {
      in->Pos[0] = SphP[i].Center[0];
      in->Pos[1] = SphP[i].Center[1];
      in->Pos[2] = SphP[i].Center[2];
    }
  else
#endif
    {
      in->Pos[0] = P[i].Pos[0];
      in->Pos[1] = P[i].Pos[1];
      in->Pos[2] = P[i].Pos[2];
    }
  in->Hsml = P[i].Auriga_Movie_Hsml;
  in->Type = P[i].Type;

  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
  MyFloat Ngb;
  MyFloat Rho;
  MyFloat Dhsml;
} data_out;

static data_out *DataResult;

/* routine to store or combine result data */
static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
      NumNgb[i] = out->Ngb;
      Rho[i]    = out->Rho;
      Dhsml[i]  = out->Dhsml;
    }
  else /* combine */
    {
      NumNgb[i] += out->Ngb;
      Rho[i] += out->Rho;
      Dhsml[i] += out->Dhsml;
    }
}

#include "../generic_comm_helpers2.h"

MyFloat *Left, *Right;
int *Active;

static void kernel_local(void)
{
  /* do local particles */
  int i;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif
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
              if(generic_polling_primary(count, Npart))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= NumPart)
          break;

        if(Active[i])
          auriga_movie_density_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        auriga_movie_density_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void auriga_movie_calculate_smoothing_lenghts(void)
{
  mpi_printf("AURIGA MOVIE: computing densities.\n");
  auriga_movie_density_compute();
  mpi_printf("AURIGA MOVIE: densities done.\n");
}

void auriga_movie_density_initialize_hsml(void)
{
  int i;

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type > 0 && P[i].Auriga_Movie_Hsml == 0)
        {
          int no = Father[i];

          if(no >= 0)
            {
              while(Nodes[no].u.d.mass < 48. * P[i].Mass)
                {
                  int p = Nodes[no].u.d.father;
                  if(p < 0)
                    break;
                  no = p;
                }

              double soft            = All.SofteningTable[P[i].SofteningType];
              P[i].Auriga_Movie_Hsml = fmax(fmin(soft, Nodes[no].len), 0.01 * soft);
            }
          else
            {
              P[i].Auriga_Movie_Hsml = All.SofteningTable[P[i].SofteningType];
            }
        }
    }
}

void auriga_movie_density_compute(void)
{
  int i, npleft, iter = 0;
  long long ntot;
  double desnumngb, desnumngbdev, t0, t1;

  CPU_Step[CPU_MISC] += measure_time();

  int Nstar = 0;
  for(i = 0; i < NumPart; i++)
    if(P[i].Type == 4)
      Nstar++;

  long long NstarTot;
  sumup_large_ints(1, &Nstar, &NstarTot);
  int doStars = 1;
  if(NstarTot < 128)
    {
      doStars = 0;
      mpi_printf("AURIGA MOVIE: Not calculating any densities for star particles, because their total number is only %lld.\n",
                 NstarTot);
    }

#ifdef HIERARCHICAL_GRAVITY
  int NActiveParticles_Saved    = TimeBinsGravity.NActiveParticles;
  int *ActiveParticleList_Saved = (int *)mymalloc("ActiveParticleList_Saved", NumPart * sizeof(int));
  for(i = 0; i < NumPart; i++)
    {
      ActiveParticleList_Saved[i] = TimeBinsGravity.ActiveParticleList[i];
    }

  TimeBinsGravity.NActiveParticles = 0;
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type > 0 && P[i].Type != 5)
        {
          TimeBinsGravity.ActiveParticleList[TimeBinsGravity.NActiveParticles] = i;
          TimeBinsGravity.NActiveParticles++;
        }
    }
#endif

  construct_forcetree(0, 1, 0, 0); /* build forcetree for all particles */

  mpi_printf("AURIGA MOVIE: hsml guesses.\n");
  auriga_movie_density_initialize_hsml();
  mpi_printf("AURIGA MOVIE: hsml guesses done.\n");

  NumNgb = (MyFloat *)mymalloc("NumNgb", NumPart * sizeof(MyFloat));
  Rho    = (MyFloat *)mymalloc("Rho", NumPart * sizeof(MyFloat));
  Dhsml  = (MyFloat *)mymalloc("Dhsml", NumPart * sizeof(MyFloat));
  Left   = (MyFloat *)mymalloc("Left", NumPart * sizeof(MyFloat));
  Right  = (MyFloat *)mymalloc("Right", NumPart * sizeof(MyFloat));
  Active = (int *)mymalloc("Active", NumPart * sizeof(int));

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type > 0 && P[i].Type != 5)
        {
          Left[i]   = 0;
          Right[i]  = 0;
          Active[i] = 1;
        }
      else
        Active[i] = 0;

      if(P[i].Type == 4 && !doStars)
        {
          Active[i]              = 0;
          P[i].Auriga_Movie_Hsml = All.BoxSize;
        }
    }

  generic_set_MaxNexport();

  desnumngb    = 48;
  desnumngbdev = 1;

  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {
      t0 = second();

      generic_comm_pattern(NumPart, kernel_local, kernel_imported);

      /* do final operations on results */
      for(i = 0, npleft = 0; i < NumPart; i++)
        {
          if(Active[i])
            {
              if(Rho[i] > 0)
                {
                  Dhsml[i] *= P[i].Auriga_Movie_Hsml / (3. * Rho[i]);
                  if(Dhsml[i] > -0.9)
                    Dhsml[i] = 1. / (1. + Dhsml[i]);
                  else
                    Dhsml[i] = 1.;
                }

              if(NumNgb[i] < (desnumngb - desnumngbdev) || NumNgb[i] > (desnumngb + desnumngbdev))
                {
                  /* need to redo this particle */
                  npleft++;

                  if(Left[i] > 0 && Right[i] > 0)
                    if((Right[i] - Left[i]) < 1.0e-3 * Left[i])
                      {
                        /* this one should be ok */
                        npleft--;
                        Active[i] = 0; /* Mark as inactive */
                        continue;
                      }

                  if(NumNgb[i] < (desnumngb - desnumngbdev))
                    Left[i] = fmax(P[i].Auriga_Movie_Hsml, Left[i]);
                  else
                    {
                      if(Right[i] != 0)
                        {
                          if(P[i].Auriga_Movie_Hsml < Right[i])
                            Right[i] = P[i].Auriga_Movie_Hsml;
                        }
                      else
                        Right[i] = P[i].Auriga_Movie_Hsml;
                    }

                  if(Right[i] > 0 && Left[i] > 0)
                    {
                      P[i].Auriga_Movie_Hsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
                    }
                  else
                    {
                      if(Right[i] == 0 && Left[i] == 0)
                        terminate("should not occur");

                      if(Right[i] == 0 && Left[i] > 0)
                        {
                          if(fabs(NumNgb[i] - desnumngb) < 0.5 * desnumngb)
                            {
                              double fac;
                              if(Dhsml[i] > 0.)
                                fac = 1 - (NumNgb[i] - desnumngb) / (3. * NumNgb[i]) * Dhsml[i];
                              else
                                fac = SEARCHFAC;

                              if(fac < SEARCHFAC)
                                P[i].Auriga_Movie_Hsml *= fac;
                              else
                                P[i].Auriga_Movie_Hsml *= SEARCHFAC;
                            }
                          else
                            P[i].Auriga_Movie_Hsml *= SEARCHFAC;
                        }

                      if(Right[i] > 0 && Left[i] == 0)
                        {
                          if(P[i].Type == 0 && fabs(NumNgb[i] - desnumngb) < 0.5 * desnumngb)
                            {
                              double fac;
                              if(Dhsml[i] > 0.)
                                fac = 1 - (NumNgb[i] - desnumngb) / (3. * NumNgb[i]) * Dhsml[i];
                              else
                                fac = SEARCHFAC;

                              if(fac > 1 / SEARCHFAC)
                                P[i].Auriga_Movie_Hsml *= fac;
                              else
                                P[i].Auriga_Movie_Hsml /= SEARCHFAC;
                            }
                          else
                            P[i].Auriga_Movie_Hsml /= SEARCHFAC;
                        }
                    }

                  if(Right[i] > 0 && Right[i] < 0.01 * All.SofteningTable[P[i].SofteningType])
                    {
                      P[i].Auriga_Movie_Hsml = 0.01 * All.SofteningTable[P[i].SofteningType];
                      Active[i]              = 0;
                      npleft--;
                    }
                }
              else
                {
                  Active[i] = 0; /* Mark as inactive */
                }
            }
        }

      sumup_large_ints(1, &npleft, &ntot);

      t1 = second();

      if(ntot > 0)
        {
          iter++;

          mpi_printf("AURIGA MOVIE DENSITY: ngb iteration %3d: need to repeat for %12lld particles. (took %g sec)\n", iter, ntot,
                     timediff(t0, t1));

          if(iter > 100)
            {
              for(i = 0; i < NumPart; i++)
                {
                  if(Active[i])
                    {
                      printf(
                          "AURIGA MOVIE DENSITY: After iter=%d still not converged for particle ID=%u, type=%d, hsml=%g, Left=%g, "
                          "Right=%g, NumNgb=%g, Dhsml=%g\n",
                          iter, P[i].ID, P[i].Type, P[i].Auriga_Movie_Hsml, Left[i], Right[i], NumNgb[i], Dhsml[i]);
                    }
                }
            }

          if(iter > MAXITER)
            terminate("failed to converge in neighbour iteration in density()\n");
        }
    }
  while(ntot > 0);

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type > 0)
        {
          P[i].Auriga_Movie_Density = Rho[i];
        }
    }

  myfree(Active);
  myfree(Right);
  myfree(Left);
  myfree(Dhsml);
  myfree(Rho);
  myfree(NumNgb);

  myfree(Father);
  myfree(Nextnode);

#ifdef BLACK_HOLES
  myfree(Tree_AuxBH_Points);
#endif
  myfree(Tree_Points);
  force_treefree();

#ifdef HIERARCHICAL_GRAVITY
  for(i = 0; i < NumPart; i++)
    {
      TimeBinsGravity.ActiveParticleList[i] = ActiveParticleList_Saved[i];
    }
  myfree(ActiveParticleList_Saved);
  TimeBinsGravity.NActiveParticles = NActiveParticles_Saved;
#endif
}

/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */

static int auriga_movie_density_evaluate(int target, int mode, int threadid)
{
  int k, numnodes, *firstnode, type, ptype;
  double hsml;
  double numngb, rho, dhsml;
  MyDouble *pos;
  int no, p;
  struct NODE *current;
  double dx, dy, dz, r2, mass;
  double h2, h3, hinv, hinv3, hinv4, r, u, wk, dwk;

  data_in local, *target_data;
  data_out out;

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

  pos  = target_data->Pos;
  hsml = target_data->Hsml;
  type = target_data->Type;

  h2    = hsml * hsml;
  h3    = h2 * hsml;
  hinv  = 1.0 / hsml;
  hinv3 = hinv * hinv * hinv;
  hinv4 = hinv3 * hinv;

  numngb = 0;
  rho    = 0;
  dhsml  = 0;

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
              p  = no;
              no = Nextnode[no];

              dx = NEAREST_X(Tree_Pos_list[3 * p + 0] - pos[0]);
              if(dx > hsml)
                continue;
              dy = NEAREST_Y(Tree_Pos_list[3 * p + 1] - pos[1]);
              if(dy > hsml)
                continue;
              dz = NEAREST_Z(Tree_Pos_list[3 * p + 2] - pos[2]);
              if(dz > hsml)
                continue;

              if((r2 = (dx * dx + dy * dy + dz * dz)) > hsml * hsml)
                continue;

              mass  = P[p].Mass;
              ptype = P[p].Type;
            }
          else if(no < Tree_MaxPart + Tree_MaxNodes) /* internal node */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                {
                  if(no <
                     Tree_FirstNonTopLevelNode) /* we reached a top-level node again, which means that we are done with the branch */
                    break;
                }

              current = &Nodes[no];

              no = current->u.d.sibling; /* in case the node can be discarded */

              double dist = hsml + 0.5 * current->len;

              dx = (MyFloat)NEAREST_X(current->center[0] - pos[0]);
              if(dx > dist)
                continue;
              dy = (MyFloat)NEAREST_Y(current->center[1] - pos[1]);
              if(dy > dist)
                continue;
              dz = (MyFloat)NEAREST_Z(current->center[2] - pos[2]);
              if(dz > dist)
                continue;
              /* now test against the minimal sphere enclosing everything */
              dist += FACT1 * current->len;
              if(dx * dx + dy * dy + dz * dz > dist * dist)
                continue;

              no = current->u.d.nextnode; /* ok, we need to open the node */
              continue;
            }
          else if(no >= Tree_ImportedNodeOffset) /* point from imported nodelist */
            {
              int n = no - Tree_ImportedNodeOffset;
              no    = Nextnode[no - Tree_MaxNodes];

              dx = NEAREST_X(Tree_Points[n].Pos[0] - pos[0]);
              if(dx > hsml)
                continue;
              dy = NEAREST_Y(Tree_Points[n].Pos[1] - pos[1]);
              if(dy > hsml)
                continue;
              dz = NEAREST_Z(Tree_Points[n].Pos[2] - pos[2]);
              if(dz > hsml)
                continue;

              if((r2 = (dx * dx + dy * dy + dz * dz)) > hsml * hsml)
                continue;

              mass  = Tree_Points[n].Mass;
              ptype = Tree_Points[n].Type;

              p = -1;
            }
          else /* pseudo particle */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                terminate("can't be");

              if(target >= 0) /* if no target is given, export will not occur */
                tree_treefind_export_node_threads(no, target, threadid);

              no = Nextnode[no - Tree_MaxNodes];
              continue;
            }

          if(ptype == type)
            {
              if(r2 < h2)
                {
                  r = sqrt(r2);
                  u = r * hinv;

                  if(u < 0.5)
                    {
                      wk  = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                      dwk = hinv4 * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
                    }
                  else
                    {
                      wk  = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
                      dwk = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
                    }

                  rho += mass * wk;
                  numngb += 4. / 3. * M_PI * wk * h3;
                  dhsml -= mass * hinv * (3. * wk + r * dwk);
                }
            }
        }
    }

  out.Ngb   = numngb;
  out.Rho   = rho;
  out.Dhsml = dhsml;

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif
