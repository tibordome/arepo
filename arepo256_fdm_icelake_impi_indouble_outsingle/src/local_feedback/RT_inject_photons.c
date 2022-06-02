/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file
 * \date        10/2014
 * \author      Rahul Kannan, Christine Simpson
 * \brief
 * \details     Injects thermal & kinetic feedback in sphere
 * \notes       Please contact Christine Simpson before using
 *              this routine christine.simpson@gmail.com
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <math.h>
#include "../MRT/RT.h"
#include "../allvars.h"
#include "../proto.h"

#if defined(LOCAL_FEEDBACK) && defined(LOCAL_FEEDBACK_PARTICLES) && defined(MRT_LOCAL_FEEDBACK)

#include "../parallel_logs.h"

static int calc; /* determines phase of calculation */
static int find_feedback_cells_evaluate(int target, int mode, int thread_id);

static struct plog_data plog;

/* communication structures */
typedef struct
{
  MyDouble Pos[3];
  MyDouble Vel[3];
  MyDouble Volume;

  MyDouble total_vol;
  MyDouble total_solid_angle;
  MyDouble h;

  MyDouble injected_photons[MRT_BINS];

  int Firstnode;
} data_in;

static data_in *DataGet;

/* temporary result structures */
static struct temp_celldata
{
  char flag_feedback;

  MyDouble total_vol;

  MyDouble total_solid_angle;

  MyDouble h, h_h, h_l, dh;

  int Ncells;

  MyDouble injected_photons[MRT_BINS];

  /* If feedback is centered on particles, record info about closest hydro cell */

  MyDouble NgbDistance;
  MyDouble NgbVolume;
  MyDouble NgbVel[3];

} * TempCellData;

/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in *in, int i, int firstnode)
{
  int k;

  for(k = 0; k < 3; k++)
    in->Pos[k] = P[i].Pos[k];

  if(TempCellData[i].Ncells == 0)
    in->Volume = 0.0;
  else
    in->Volume = TempCellData[i].NgbVolume;

#if defined(EXTERNALSHEARBOX_KSRATE_RANDOM) || defined(EXTERNALSHEARBOX_MIXED_INJECTION)
  for(k = 0; k < 3; k++)
    in->Vel[k] = TempCellData[i].NgbVel[k];
#else
  for(k = 0; k < 3; k++)
    in->Vel[k] = P[i].Vel[k];
#endif

  in->total_vol         = TempCellData[i].total_vol;
  in->total_solid_angle = TempCellData[i].total_solid_angle;
  in->h                 = TempCellData[i].h;

  if(calc > 0)
    {
      for(int num = 0; num < MRT_BINS; num++)
        in->injected_photons[num] = TempCellData[i].injected_photons[num];
    }

  in->Firstnode = firstnode;
}

typedef struct
{
  MyDouble total_vol;

  MyDouble total_solid_angle;
  MyDouble aterm1, b, c, tot_ekin_b, h;
  MyDouble resid[3];
  int Ncells;

  MyDouble NgbDistance;
  MyDouble NgbVolume;
  MyDouble NgbRadius;
  MyDouble NgbVel[3];

  MyDouble injected_photons[MRT_BINS];

} data_out;

static data_out *DataResult;

/* routine to store or combine result data */
static void out2particle(data_out *out, int i, int mode)
{
  int k;
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
      if(calc == 0)
        {
          TempCellData[i].total_vol         = out->total_vol;
          TempCellData[i].total_solid_angle = out->total_solid_angle;
          TempCellData[i].Ncells            = out->Ncells;
#if defined(LOCAL_FEEDBACK_PARTICLES) || defined(EXTERNALSHEARBOX_KSRATE_RANDOM) || defined(EXTERNALSHEARBOX_MIXED_INJECTION)
          TempCellData[i].NgbDistance = out->NgbDistance;
          TempCellData[i].NgbVolume   = out->NgbVolume;
          for(k = 0; k < 3; k++)
            TempCellData[i].NgbVel[k] = out->NgbVel[k];
#endif
        }
      else
        {
          for(int num = 0; num < MRT_BINS; num++)
            TempCellData[i].injected_photons[num] = out->injected_photons[num];
        }
    }
  else /* combine */
    {
      if(calc == 0)
        {
          TempCellData[i].total_vol += out->total_vol;
          TempCellData[i].total_solid_angle += out->total_solid_angle;
          TempCellData[i].Ncells += out->Ncells;
#if defined(LOCAL_FEEDBACK_PARTICLES) || defined(EXTERNALSHEARBOX_KSRATE_RANDOM) || defined(EXTERNALSHEARBOX_MIXED_INJECTION)
          if(out->NgbDistance < TempCellData[i].NgbDistance)
            {
              TempCellData[i].NgbDistance = out->NgbDistance;
              TempCellData[i].NgbVolume   = out->NgbVolume;
              for(k = 0; k < 3; k++)
                TempCellData[i].NgbVel[k] = out->NgbVel[k];
            }
#endif
        }
      else
        {
          for(int num = 0; num < MRT_BINS; num++)
            TempCellData[i].injected_photons[num] = out->injected_photons[num];
        }
    }
}

#include "../generic_comm_helpers2.h"

/* Routine that injects feedback either around previously created star particles, or if star particles
   are not used, instantaneously injected based on local gas properties.
   If instantaneous injection is used, the injection rate is computed either based on the KS rate or
   based on the local dynamical time.
   The total amount of energy injected per SN event, and the fraction of that energy that takes a kinetic
   form, is set by run time parameters. */

static void kernel_local(void)
{
  int i, idx;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif

#pragma omp parallel private(i, idx)
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
              if(generic_polling_primary(count, TimeBinsGravity.NActiveParticles))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        idx = NextParticle++;
        if(idx >= TimeBinsGravity.NActiveParticles)
          break;

        i = TimeBinsGravity.ActiveParticleList[idx];
        if(i < 0)
          continue;

        if(TempCellData[i].flag_feedback == 1)
          find_feedback_cells_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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
        find_feedback_cells_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void inject_photons(void)
{
  CPU_Step[CPU_MISC] += measure_time();
  mpi_printf("CMS_INJECT_PHOTONS: Starting to inject photons\n");

  int i;
  double rho, sfr, dz, dt, prob, sech, deltaM;
  double rand;

#ifdef LOCAL_FEEDBACK_PARTICLES
  double tform, Age, N_SN;
#endif

  TempCellData = mymalloc("TemCellData", NumPart * sizeof(struct temp_celldata));
  memset(TempCellData, 0, NumPart * sizeof(struct temp_celldata));

  parallel_log_allocate(&plog, 0);

  int count_spawned = 0, count_converted = 0, count_feedback = 0;

  int active, Active;

  active = Active = 0;
  /*** Let's first flag the cells or paritcles that inject a feedback event ***/
  /* loop over all active cells and particles */
  int idx;
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      TempCellData[i].flag_feedback = 0;
      TempCellData[i].Ncells        = 0;
#ifdef LOCAL_FEEDBACK_PARTICLES /* Feedback events created from star particles that were previously created */
      if(P[i].Type == 4)
        {
          if(P[i].TagSNe == 0 || P[i].FeedbackDone > 0)
            continue;

          TempCellData[i].flag_feedback = 1;
          active++;
        }
#endif
    }

  MPI_Allreduce(&active, &Active, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  mpi_printf("CMS_INJECT_PHOTONS: Finished selecting %d high mass stars\n", Active);
  /*** Now, let's go for each of the selected cells three times over all neighbours in a search sphere.
       on iteration -1: Find closest gas cell
       on iteration 0: Just determine total volume of the neighbours cells in the aperture
       on iteration 1: Compute momentum we will add to the gas
       on iteration 2: Inject energy, mass and momentum
  ***/

  generic_set_MaxNexport();

  int j, nexport, nimport, place, ngrp, recvTask, ndone_flag, ndone, dummy;
  int repeat_calc0 = 0;
  double cell_radius, h, h2;

  for(calc = 0; calc < 2; calc++)
    {
      int idx;
      int enlarge_radius = 1;
      /* loop over all active cells and particles */

      for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
        {
          i = TimeBinsGravity.ActiveParticleList[idx];
          if(i < 0)
            continue;

          double dt;
#ifndef MRT_INJECT_PHOTONS_EVERY_STEP
          dt = (P[i].TimeBinGrav ? (((integertime)1) << P[i].TimeBinGrav) : 0) * All.Timebase_interval;
#else
          dt = All.TimeStep;
#endif

          if(TempCellData[i].flag_feedback)
            {
              /* on first pass, tie radius to softening length */
              if(calc == 0)
                {
#if defined(LOCAL_FEEDBACK_PARTICLES) || defined(EXTERNALSHEARBOX_KSRATE_RANDOM) || defined(EXTERNALSHEARBOX_MIXED_INJECTION)
                  if(P[i].Type == 4)
                    {
#ifndef MRT_INJECT_PHOTONS_EVERY_STEP
                      cell_radius = 25.0 * All.MinimumComovingHydroSoftening;
#else
                      cell_radius = All.BoxSize / 20.0;
#endif
                    }
                  else
                    cell_radius = pow(SphP[i].Volume * 3.0 / 4.0 / M_PI, 1.0 / 3.0);
#else
                  cell_radius = pow(SphP[i].Volume * 3.0 / 4.0 / M_PI, 1.0 / 3.0);
#endif
                  TempCellData[i].h_l = 0.0;
                  TempCellData[i].h_h = 6.0 * cell_radius;
                  TempCellData[i].dh  = 0.5 * (TempCellData[i].h_h - TempCellData[i].h_l);
                  TempCellData[i].h   = 0.5 * (TempCellData[i].h_h + TempCellData[i].h_l);
                }
              h2 = TempCellData[i].h * TempCellData[i].h;

#ifdef LOCAL_FEEDBACK_PARTICLES
              if(calc >= 1)
                if(TempCellData[i].h == 0.0)
                  terminate("CMS_FEEDBACK[%d]: zero radius feedback sphere: id = %d; NgbDistance = %g; NgbVolume = %g\n", ThisTask,
                            P[i].ID, TempCellData[i].NgbDistance, TempCellData[i].NgbVolume);
#endif
              if(calc == 1)
                {
                  tform = P[i].StellarAge;
                  Age   = All.Time - tform;
                  for(int num = 0; num < MRT_BINS; num++)
                    TempCellData[i].injected_photons[num] = inject_photons_from_star(num, dt, Age);

                  count_feedback++;
                }
            }
        }  // end loop over particles

      while(enlarge_radius)
        {
          generic_comm_pattern(TimeBinsGravity.NActiveParticles, kernel_local, kernel_imported);

          enlarge_radius = 0;
          if(calc == 0)
            {
              for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
                {
                  i = TimeBinsGravity.ActiveParticleList[idx];
                  if(i < 0)
                    continue;

                  if(TempCellData[i].flag_feedback && (TempCellData[i].Ncells < 14 || TempCellData[i].Ncells > 16))
                    {
                      enlarge_radius = 1;

                      if(TempCellData[i].Ncells < 14)
                        {
                          TempCellData[i].h_l += TempCellData[i].dh;
                          TempCellData[i].h_h += TempCellData[i].dh;
                          TempCellData[i].h = 0.5 * (TempCellData[i].h_h + TempCellData[i].h_l);
                        }
                      else if(TempCellData[i].Ncells > 16)
                        {
                          TempCellData[i].h_h -= TempCellData[i].dh;
                          TempCellData[i].dh *= 0.5;
                          TempCellData[i].h = 0.5 * (TempCellData[i].h_h + TempCellData[i].h_l);
                        }
                    }
                }
            }
          MPI_Allreduce(MPI_IN_PLACE, &enlarge_radius, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
    }

  /* loop over all active cells and particles */
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(TempCellData[i].flag_feedback)
        {
          TempCellData[i].flag_feedback = 0;
          // double resid = sqrt(TempCellData[i].resid[0] * TempCellData[i].resid[0] + TempCellData[i].resid[1] *
          // TempCellData[i].resid[1] + TempCellData[i].resid[2] * TempCellData[i].resid[2]); pl_fprintf( &plog, "CMS_FEEDBACK:    Cell
          // %d: resid/deltap %f \n", P[i].ID, resid / TempCellData[i].deltap);
#ifdef EXTERNALSHEARBOX_KSRATE_RANDOM
          P[i].ID   = 0;
          P[i].Mass = 0;
          timebin_remove_particle(&TimeBinsGravity, idx, P[i].TimeBinGrav);
#endif
        }
    }

  /* Print statement */
  int countall_feedback;
  MPI_Allreduce(&count_feedback, &countall_feedback, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  // if(countall_feedback > 0)
  mpi_printf("CMS_INJECT_PHOTONS: Injected %d events\n", countall_feedback);

  parallel_log_clear(&plog, PLOG_LOCAL_FEEDBACK);
  parallel_log_free(&plog);

  myfree(TempCellData);
}

/*! This function represents the core of the star density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
static int find_feedback_cells_evaluate(int target, int mode, int thread_id)
{
  int j, n, k;
  double h2;
  double dx, dy, dz, r2;
  double pi = 3.14159265358979323846;
  /* mode is */
#ifdef LOCAL_FEEDBACK_PARTICLES
  double density_threshold = 10.0;
  density_threshold *= PROTONMASS * All.UnitLength_in_cm * All.UnitLength_in_cm * All.UnitLength_in_cm / All.UnitMass_in_g;
#endif

  MyDouble *pos, *frame_vel, target_vol, vol;
  double local_resid[3];
  double cell_radius;

  double total_vol, fkin, fcr, SNEnergy, SNMassReturn, deltap, SNMetalReturn, injected_photons[MRT_BINS], total_solid_angle;
  double aterm1, a, b, c, tot_ekin_b;
  /* determine input parameters for particle */

  int numnodes, *firstnode;
  data_in local, *in;
  data_out out;

  if(mode == MODE_LOCAL_PARTICLES) /* local particle */
    {
      particle2in(&local, target, 0);
      in = &local;

      numnodes  = 1;
      firstnode = NULL;
    }
  else /* imported particle */
    {
      in = &DataGet[target];
      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  pos       = in->Pos;
  frame_vel = in->Vel;
  double h  = in->h;
  if(calc > 0)
    {
      total_vol         = in->total_vol;
      total_solid_angle = in->total_solid_angle;
      for(int num = 0; num < MRT_BINS; num++)
        injected_photons[num] = in->injected_photons[num];
    }

#if !defined(LOCAL_FEEDBACK_PARTICLES) && !defined(EXTERNALSHEARBOX_KSRATE_RANDOM) && !defined(EXTERNALSHEARBOX_MIXED_INJECTION)
  target_vol = in->Volume;
#endif

  double cell_aterm1, cell_b, cell_c, cell_tot_ekin_b;
  double local_vol_sum, local_aterm1, local_b, local_c, local_tot_ekin_b, local_solid_angle_sum;
  int local_Ncells, ingb;
  double ngb_distance2;
  double local_photons[MRT_BINS];
#if defined(LOCAL_FEEDBACK_PARTICLES) || defined(EXTERNALSHEARBOX_KSRATE_RANDOM) || defined(EXTERNALSHEARBOX_MIXED_INJECTION)
  ngb_distance2 = boxSize_X * boxSize_Z;
  ingb          = -1;
#endif

  if(calc == 0)
    {
      local_vol_sum         = 0.0;
      local_solid_angle_sum = 0.0;
      local_Ncells          = 0;
    }

  if(calc == 1)
    {
      for(int num = 0; num < MRT_BINS; num++)
        local_photons[num] = 0.0;
    }
  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  h2 = h * h;
  for(n = 0; n < nfound; n++)
    {
      j = Thread[thread_id].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0 && P[j].Type == PTYPE_GAS)
        {
          double cellrad  = get_cell_radius(j);
          double cellarea = M_PI * cellrad * cellrad;
          dx              = nearest_x(pos[0] - SphP[j].Center[0]);
          dy              = nearest_y(pos[1] - SphP[j].Center[1]);
          dz              = nearest_z(pos[2] - SphP[j].Center[2]);

          double dd[3] = {dx, dy, dz};

          r2 = dx * dx + dy * dy + dz * dz;
#if defined(LOCAL_FEEDBACK_PARTICLES) || defined(EXTERNALSHEARBOX_KSRATE_RANDOM) || defined(EXTERNALSHEARBOX_MIXED_INJECTION)
          if(r2 < ngb_distance2 && calc == 0)
            {
              ngb_distance2 = r2;
              ingb          = j;
            }
#endif
          if(r2 < h2)
            {
              /* compute total volume for normalization purposes */

              if(calc == 0)
                {
                  local_solid_angle_sum += 0.5 * (1.0 - 1.0 / sqrt(1.0 + cellarea / M_PI / r2));
                  local_vol_sum += SphP[j].Volume;
                  local_Ncells += 1;
                  // local_photons = injected_photons ;
                }
              /* compute kinetic energy injection constants */

              if(calc == 1)
                {
                  for(int num = 0; num < MRT_BINS; num++)
                    local_photons[num] = injected_photons[num];

                  vol = SphP[j].Volume;

                  double fvol = vol / total_vol;

                  double frac = fvol;  // 0.5*(1.0 - 1.0/sqrt(1.0+cellarea/M_PI/r2)) / total_solid_angle;
                  // printf("%g %g %g\n", injected_photons[0], frac, frac*injected_photons[0]) ;
                  for(int num = 0; num < MRT_BINS; num++)
                    {
                      SphP[j].Cons_DensPhot[num] += frac * injected_photons[num];
                      //		    for(int dims=0;dims<3;dims++)
                      // SphP[j].Cons_RT_F[num][dims] += 0.9999999*c_internal_units*frac*injected_photons[num] * dd[dims]/sqrt(r2) ;
                    }
                }
            }
        }
    }

  /* now store the results at the appropriate place */
  if(calc == 0)
    {
#if defined(LOCAL_FEEDBACK_PARTICLES) || defined(EXTERNALSHEARBOX_KSRATE_RANDOM) || defined(EXTERNALSHEARBOX_MIXED_INJECTION)
      out.NgbDistance = sqrt(ngb_distance2);
      if(ingb >= 0)
        {
          out.NgbVolume = SphP[ingb].Volume;
          for(k = 0; k < 3; k++)
            out.NgbVel[k] = P[ingb].Vel[k];
        }
      else
        {
          out.NgbVolume = 0.0;
          for(k = 0; k < 3; k++)
            out.NgbVel[k] = 0.0;
        }
#else
      out.NgbDistance = 0.0;
      out.NgbVolume   = target_vol;
      for(k = 0; k < 3; k++)
        out.NgbVel[k] = frame_vel[k];
#endif
      out.Ncells            = local_Ncells;
      out.total_vol         = local_vol_sum;
      out.total_solid_angle = local_solid_angle_sum;
    }
  else
    {
      for(int num = 0; num < MRT_BINS; num++)
        out.injected_photons[num] = local_photons[num];
    }

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif
