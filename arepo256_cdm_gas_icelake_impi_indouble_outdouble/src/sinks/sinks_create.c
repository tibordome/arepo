/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sinks/sinks_create.c
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       Sink particle creation
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - 05.2017 Federico Marinacci extensive rewtriting to exploit the new
 *   communication structure
 */

#include "../allvars.h"
#include "../proto.h"

#ifdef SINKS

void sinks_find_candidate(void);
void sinks_get_center_of_mass(void);
void sinks_add_sink(void);

struct ACC_struct *ACC;
static int NewSinks;

void sinks_find_candidate(void)
{
  int idx, i, j, flag, totflag;
  double nh, mu;

  struct
  {
    double nh_max;
    int task;

  } local, global;

  flag = 0;

  local.nh_max = 0.0;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      nh = SKD.NHFac * SphP[i].Density;

      if(flag == 0 || nh > local.nh_max)
        {
          flag = 1;

          local.nh_max = nh;

          SKD.Index   = i;
          SKD.ID      = P[i].ID;
          SKD.TimeBin = P[i].TimeBinHydro;
          mu          = (1. + 4. * HE_ABUND) / fmax(1. + HE_ABUND - SphP[i].Abund[1] + SphP[i].Abund[2], 0.);
          SKD.Temp    = (SphP[i].Gamma - 1) * SphP[i].Utherm / BOLTZMANN * All.UnitEnergy_in_cgs / All.UnitMass_in_g * mu * PROTONMASS;
          // SKD.Temp = 6000;
          SKD.AccRad = get_accretion_radius(P[i].Mass, SKD.Temp);

          for(j = 0; j < 3; j++)
            SKD.Pos[j] = P[i].Pos[j];
        }
    }

  MPI_Allreduce(&flag, &totflag, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(totflag)
    {
      local.task = ThisTask;

      MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

      SKD.Task  = global.task;
      SKD.NHMax = global.nh_max;

      MPI_Bcast(&SKD.Index, 1, MPI_INT, SKD.Task, MPI_COMM_WORLD);
      MPI_Bcast(&SKD.ID, 1, MPI_INT, SKD.Task, MPI_COMM_WORLD);
      MPI_Bcast(SKD.Pos, 3, MPI_DOUBLE, SKD.Task, MPI_COMM_WORLD);
      MPI_Bcast(&SKD.TimeBin, 1, MPI_INT, SKD.Task, MPI_COMM_WORLD);
      MPI_Bcast(&SKD.AccRad, 1, MPI_DOUBLE, SKD.Task, MPI_COMM_WORLD);
      MPI_Bcast(&SKD.Temp, 1, MPI_DOUBLE, SKD.Task, MPI_COMM_WORLD);
    }
  else
    SKD.NHMax = 0;
}

void sinks_get_center_of_mass(void)
{
  int i, j;

  if(ThisTask == SKD.Task)
    {
      for(i = 0; i < 3; i++)
        SKD.CoMPos[i] = SKD.CoMVel[i] = 0;

      SKD.Mass = 0;

      for(i = 0; i < NewSinks; i++)
        {
          for(j = 0; j < 3; j++)
            {
              SKD.CoMPos[j] += ACC[i].Pos[j];
              SKD.CoMVel[j] += ACC[i].Vel[j];
            }

          SKD.Mass += ACC[i].Mass;
        }

      for(i = 0; i < 3; i++)
        {
          SKD.CoMPos[i] /= SKD.Mass;
          SKD.CoMVel[i] /= SKD.Mass;
        }
    }
}

void sinks_add_sink(void)
{
  int i, igas, isink;

  if(ThisTask == SKD.Task)
    {
      igas  = SKD.Index;
      isink = NumPart;

      if(igas >= isink)
        terminate("This is not possible! isink %d < igas %d", isink, igas);

      if(isink >= All.MaxPart)
        terminate("There is no space to add the sink: isink %d, All.MaxPart %d", isink, All.MaxPart);

      if(NumSinks + 1 > All.MaxPartSinks)
        terminate("There is no space to add the SinkP struct: NumSink + 1=%d, All.MaxPartSinks=%d", NumSinks + 1, All.MaxPartSinks);

      P[isink] = P[igas];

      for(i = 0; i < 3; i++)
        {
          P[isink].Pos[i] = SKD.CoMPos[i];
          P[isink].Vel[i] = SKD.CoMVel[i];
        }

      P[isink].Mass          = SKD.Mass;
      P[isink].Type          = 5;
      P[isink].SofteningType = All.SofteningTypeOfPartType[P[isink].Type];

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
      if(((1 << P[isink].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
        P[isink].SofteningType = get_softening_type_from_mass(P[isink].Mass);
#endif

      P[isink].ID = SKD.ID;

      P[isink].TimeBinGrav  = SKD.TimeBin;
      P[isink].TimeBinHydro = P[isink].TimeBinGrav;
      P[isink].TimeBinSink  = P[isink].TimeBinHydro;
      P[isink].AuxDataID    = NumSinks;

      SKP(isink).PID        = isink;
      SKP(isink).Sinks_Hsml = SKD.AccRad;
      SKP(isink).Temp       = SKD.Temp;
      SKP(isink).SwallowID  = 0;

      timebin_add_particle(&TimeBinsGravity, isink, -1, P[isink].TimeBinGrav,
                           TimeBinSynchronized[P[isink].TimeBinGrav]);  // Marinacci, 2016
      timebin_add_particle(&TimeBinsSinksAccretion, isink, -1, P[isink].TimeBinHydro, TimeBinSynchronized[P[isink].TimeBinHydro]);

      NumPart++;
      NumSinks++;

      printf("SINKS: Mass %g ID %d Pos %g|%g|%g Vel %g|%g|%g\n", P[isink].Mass, P[isink].ID, P[isink].Pos[0], P[isink].Pos[1],
             P[isink].Pos[2], P[isink].Vel[0], P[isink].Vel[1], P[isink].Vel[2]);
      printf("SINKS: isink %i AccRad %g Temp %g\n", isink, SKD.AccRad, SKD.Temp);
      printf("SINKS: removed cells %d\n", ACC[0].rcells);
    }

  SKD.TotNumSinks++;
  All.TotNumPart++;
  All.TotN_sinks++;
}

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  double Pos[3];
  double hsml;

  int Firstnode;
} data_in;

static data_in *DataGet;

/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in *in, int i, int firstnode)
{
  for(int k = 0; k < 3; k++)
    in->Pos[k] = SKD.Pos[k];

  in->hsml = SKD.AccRad;

  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
  double Pos[3];
  double Vel[3];
  double Mass;
  int rcells;
} data_out;

static data_out *DataResult;

/* routine to store or combine result data */
static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
      for(int k = 0; k < 3; k++)
        {
          ACC[i].Pos[k] = out->Pos[k];
          ACC[i].Vel[k] = out->Vel[k];
        }

      ACC[i].Mass   = out->Mass;
      ACC[i].rcells = out->rcells;
    }
  else /* combine */
    {
      for(int k = 0; k < 3; k++)
        {
          ACC[i].Pos[k] += out->Pos[k];
          ACC[i].Vel[k] += out->Vel[k];
        }

      ACC[i].Mass += out->Mass;
      ACC[i].rcells += out->rcells;
    }
}

static int sinks_create_evaluate(int target, int mode, int threadid);

#include "../generic_comm_helpers2.h"

static void kernel_local(void)
{
  int i;
#pragma omp parallel private(i)
  {
    int j, threadid = get_thread_num();

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= NewSinks)
          break;

        sinks_create_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}

static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, count = 0;
#pragma omp parallel private(i)
  {
    int threadid = get_thread_num();

    while(1)
      {
#pragma omp atomic capture
        i = count++;

        if(i >= Nimport)
          break;

        sinks_create_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void sinks_create(void)
{
  TIMER_START(CPU_SINKS_CREATE);

  NewSinks = 0;

  sinks_find_candidate();

  mpi_printf(
      "SINKS: task %d SKD task %d SKD index %d SKD ID %d SKD NHMax %g SKD NThresh %g SKD pos %g|%g|%g SKD TimeBin %d SKD AccRad %g "
      "SKD Temp %g\n",
      ThisTask, SKD.Task, SKD.Index, SKD.ID, SKD.NHMax, SKD.NHThresh, SKD.Pos[0], SKD.Pos[1], SKD.Pos[2], SKD.TimeBin, SKD.AccRad,
      SKD.Temp);

  if(SKD.NHMax > SKD.NHThresh)
    {
      if(ThisTask == SKD.Task)
        NewSinks = 1;

      ACC = mymalloc("ACC", NewSinks * sizeof(struct ACC_struct));

      generic_set_MaxNexport();

      generic_comm_pattern(NewSinks, kernel_local, kernel_imported);

#ifndef SGCHEM
      timebin_cleanup_list_of_active_particles(&TimeBinsHydro);
      timebin_cleanup_list_of_active_particles(&TimeBinsGravity);
#endif

      sinks_get_center_of_mass();

      sinks_add_sink();

      myfree(ACC);

      mpi_printf("SINKS: Created a new sink.\n");
    }

  TIMER_STOP(CPU_SINKS_CREATE);
}

static int sinks_create_evaluate(int target, int mode, int threadid)
{
  int j, n, numnodes, *firstnode;
  int rcells = 0;
  double r2, rad, rad2, dx, dy, dz;
  double pos[3] = {0.0, 0.0, 0.0};
  double vel[3] = {0.0, 0.0, 0.0};
  double acc[3] = {0.0, 0.0, 0.0};
  double mass   = 0.0;

  data_in local, *in;
  data_out out;

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

  rad  = in->hsml / sqrt(SKD.DistFac);
  rad2 = rad * rad;

  int nfound = ngb_treefind_variable_threads(in->Pos, rad, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      if(P[j].Type == 0)
        {
          dx = NGB_PERIODIC_LONG_X(in->Pos[0] - P[j].Pos[0]);
          dy = NGB_PERIODIC_LONG_Y(in->Pos[1] - P[j].Pos[1]);
          dz = NGB_PERIODIC_LONG_Z(in->Pos[2] - P[j].Pos[2]);

          r2 = dx * dx + dy * dy + dz * dz;

          if((r2 < rad2) && (P[j].Mass != 0) && (P[j].ID != 0) && (TimeBinSynchronized[P[j].TimeBinHydro] > 0) &&
             (P[j].TimeBinHydro <= SKD.TimeBin))
            {
              for(int k = 0; k < 3; k++)
                {
                  pos[k] += P[j].Mass * P[j].Pos[k];
                  vel[k] += P[j].Mass * P[j].Vel[k];
                }

              mass += P[j].Mass;

              P[j].Mass = 0;
              P[j].ID   = 0;

#ifdef VORONOI_DYNAMIC_UPDATE
              voronoi_remove_connection(j);
#endif
              rcells++;
            }
        }
    }

  for(int k = 0; k < 3; k++)
    {
      out.Pos[k] = pos[k];
      out.Vel[k] = vel[k];
    }

  out.Mass   = mass;
  out.rcells = rcells;

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif
