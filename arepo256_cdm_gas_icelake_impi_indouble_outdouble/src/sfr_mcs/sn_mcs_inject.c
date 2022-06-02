/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sfr_mcs/sn_mcs_inject.c
 * \date        08/2018
 * \author      Matthew C Smith
 * \brief
 * \details     Private to M C Smith, but collaboration encouraged.
                Originally developed in 2015, ported into main repo 2018.
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

#if defined(SFR_MCS) && defined(SN_MCS)

#if !defined(VORONOI_DYNAMIC_UPDATE)
#error "SN_MCS requires VORONOI_DYNAMIC_UPDATE."
#endif

#ifdef TRACER_MC_CHECKS
#error "TRACER_MC_CHECKS does not work with SN_MCS"
#endif

#define SN_MCS_NNGB_MAX_DEFAULT 64 /* This should never really be reached */
#if !defined(SN_MCS_NNGB_MAX) || (SN_MCS_NNGB_MAX < SN_MCS_NNGB_MAX_DEFAULT)
#define SN_MCS_NNGB_MAX SN_MCS_NNGB_MAX_DEFAULT
#else
#define SN_MCS_NNGB_MAX_ADVISORY
#endif

#ifndef SN_MCS_INITIAL_DRIVING
static struct snhost_in
{
  int HostTask;
  int HostIndex;
#ifdef SN_MCS_LOCATION_RECORDS
  int OriginTask;
  int OriginStarIndex;
#endif
  double mass_deposited;
#if defined(IMF_SAMPLING_MCS) || defined(SN_MCS_VARIABLE_EJECTA)
  double metal_deposited;
#endif
  double energy_deposited;
  double starvel[3];
#ifdef SN_MCS_LOG_DETAILS
  MyIDType StarID;
  double StarPos[3];
#endif
} * SnhostIn, *SnhostGet;

#ifdef SN_MCS_LOCATION_RECORDS
static struct snhostinfo_in
{
#ifdef SN_MCS_LOCATION_RECORDS_CHECKS
  int OriginTask;
#endif
  int OriginStarIndex;
  MyFloat host_density;
  MyFloat host_temperature;
} * SnhostinfoIn, *SnhostinfoGet;
#endif
#endif

static struct iso_weights
{
  double r;
  double x_plus[3];
  double x_minus[3];
  double omega;
  double w[3];
  int Task;
  int Index;
} * Neighbour_weights;

static struct snvar_in
{
  int Task;
  int Index;
  double dm;
#if defined(IMF_SAMPLING_MCS) || defined(SN_MCS_VARIABLE_EJECTA)
  double dmZ;
#endif
  double dp[3];
  double dE;
#ifdef SN_MCS_MECHANICAL
  double E_51;
  double p_ej_tot;
#endif
} * SnvarIn, *SnvarGet;

/* Functions used with qsort */
static int snhost_compare(const void *a, const void *b);
static int snvar_compare(const void *a, const void *b);

#ifndef SN_MCS_INITIAL_DRIVING
/* This function finds the host cells of SNe. Distributes SNe variables to hosts,
whether local or on other tasks.
int n_sn_on_task: the number of star particles with a SNe event this timestep on
this task. */
void find_sn_host_cells_and_distribute(int n_sn_on_task)
{
  int i, idx, j, stid, searchid;
  int n_sn_processed, nexport, nimport, ngrp, recvTask;
  mesh_search_data *searchdata;
  int *starindex;
  MPI_Status status;

  mpi_printf("SN_MCS: Finding SN host cells...\n");

  searchdata = (mesh_search_data *)mymalloc("searchdata", n_sn_on_task * sizeof(mesh_search_data));
  starindex  = (int *)mymalloc("starindex", n_sn_on_task * sizeof(int));

#ifdef TRACER_MC
  start_MC_tracer(N_tracer); /* allocate buffer for tracer exchange */
#endif

  n_sn_processed = 0;

  /* Load search data for all star particles with a SNe event this timestep,
  then search. */
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if((P[i].Type == 4) && (STP(i).N_SN > 0))
        {
          searchdata[n_sn_processed].Pos[0] = P[i].Pos[0];
          searchdata[n_sn_processed].Pos[1] = P[i].Pos[1];
          searchdata[n_sn_processed].Pos[2] = P[i].Pos[2];
          starindex[n_sn_processed]         = i;
          n_sn_processed++;
        }
    }

  assert(n_sn_processed == n_sn_on_task); /* Check the number is correct */

  find_nearest_meshpoint_global(searchdata, n_sn_processed, 0, 0);

  /* Check how many particles need to be exported */
  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(i = 0, nexport = 0; i < n_sn_processed; i++)
    {
      if(searchdata[i].Task != ThisTask)
        {
          Send_count[searchdata[i].Task] += 1;
          nexport++;
        }
    }

  /* Send the export counts to other tasks and get receive counts */
  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(i = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; i < NTask; i++)
    {
      nimport += Recv_count[i];

      if(i > 0)
        {
          Send_offset[i] = Send_offset[i - 1] + Send_count[i - 1];
          Recv_offset[i] = Recv_offset[i - 1] + Recv_count[i - 1];
        }
    }

    /* Allocate export/import buffers */
#ifdef SN_MCS_LOCATION_RECORDS
  SnhostinfoIn = (struct snhostinfo_in *)mymalloc("SnhostinfoGet", nimport * sizeof(struct snhostinfo_in));
#endif
  SnhostGet = (struct snhost_in *)mymalloc("SnhostGet", nimport * sizeof(struct snhost_in));
  SnhostIn  = (struct snhost_in *)mymalloc("SnhostIn", nexport * sizeof(struct snhost_in));

  /* Distribute SN variables either locally or to the export buffer */
  for(i = 0, j = 0; i < n_sn_processed; i++)
    {
      stid = starindex[i];
      if(searchdata[i].Task == ThisTask)
        {
          if(P[searchdata[i].u.Index].ID != 0)
            {
              searchid = searchdata[i].u.Index;
              SphP[searchid].mass_deposited += STP(stid).M_ej;
#if defined(IMF_SAMPLING_MCS) || defined(SN_MCS_VARIABLE_EJECTA)
              SphP[searchid].metal_deposited += STP(stid).M_ej * STP(stid).Z_ej;
#endif
              SphP[searchid].energy_deposited += STP(stid).N_SN * All.SupernovaEnergy * All.cf_atime * All.cf_atime;
              SphP[searchid].starvel[0] += P[stid].Vel[0];
              SphP[searchid].starvel[1] += P[stid].Vel[1];
              SphP[searchid].starvel[2] += P[stid].Vel[2];
              SphP[searchid].N_SN_hosted += 1;
#ifdef TRACER_MC
              consider_moving_tracers_local(stid, searchid, fmax(STP(stid).M_ej / (P[stid].Mass + STP(stid).M_ej), 1.0));
#endif
#ifdef SN_MCS_LOG
              sn_add_to_log(searchid);
#endif
#ifdef SN_MCS_LOCATION_RECORDS
              STP(stid).SNPos[0]  = P[stid].Pos[0];
              STP(stid).SNPos[1]  = P[stid].Pos[1];
              STP(stid).SNPos[2]  = P[stid].Pos[2];
              STP(stid).SNTime    = All.Time;
              STP(stid).SNDensity = SphP[searchid].Density;
#ifdef GRACKLE
              STP(stid).SNTemperature = get_temp_individual_cell_grackle(searchid);
#endif
#endif
#ifdef SN_MCS_LOG_DETAILS
#ifdef GRACKLE
              double temp = get_temp_individual_cell_grackle(searchid);
#else
              double temp = -1;
#endif
              fprintf(FdSNDetails, "SN=%llu %g %g %g %g %g %g %g %g %g %g %g  %g %d\n", (long long)P[stid].ID, All.Time,
                      P[stid].Pos[0], P[stid].Pos[1], P[stid].Pos[2], P[stid].Vel[0], P[stid].Vel[1], P[stid].Vel[2],
                      SphP[searchid].Density, temp, SphP[searchid].Metallicity, STP(stid).M_ej,
                      STP(stid).N_SN * All.SupernovaEnergy * All.cf_atime * All.cf_atime, 0);
#endif
            }
        }
      else
        {
          SnhostIn[j].HostTask       = searchdata[i].Task;
          SnhostIn[j].HostIndex      = searchdata[i].u.Index;
          SnhostIn[j].mass_deposited = STP(stid).M_ej;
#if defined(IMF_SAMPLING_MCS) || defined(SN_MCS_VARIABLE_EJECTA)
          SnhostIn[j].metal_deposited = STP(stid).M_ej * STP(stid).Z_ej;
#endif
          SnhostIn[j].energy_deposited = STP(stid).N_SN * All.SupernovaEnergy * All.cf_atime * All.cf_atime;
          SnhostIn[j].starvel[0]       = P[stid].Vel[0];
          SnhostIn[j].starvel[1]       = P[stid].Vel[1];
          SnhostIn[j].starvel[2]       = P[stid].Vel[2];
#ifdef TRACER_MC
          consider_moving_tracers(stid, searchdata[i].Task, searchdata[i].u.Index, 0,
                                  fmax(STP(stid).M_ej / (P[stid].Mass + STP(stid).M_ej), 1.0));
#endif
#ifdef SN_MCS_LOG_DETAILS
          SnhostIn[j].StarID     = P[stid].ID;
          SnhostIn[j].StarPos[0] = P[stid].Pos[0];
          SnhostIn[j].StarPos[1] = P[stid].Pos[1];
          SnhostIn[j].StarPos[2] = P[stid].Pos[2];
#endif
#ifdef SN_MCS_LOCATION_RECORDS
          STP(stid).SNPos[0] = P[stid].Pos[0];
          STP(stid).SNPos[1] = P[stid].Pos[1];
          STP(stid).SNPos[2] = P[stid].Pos[2];
          STP(stid).SNTime   = All.Time;
#ifdef SN_MCS_LOCATION_RECORDS_CHECKS
          SnhostIn[j].OriginTask = ThisTask;
#endif
          SnhostIn[j].OriginStarIndex = P[stid].AuxDataID;
#endif
          j++;
        }

      /* Reset star properties */
      STP(stid).M_ej = 0;
      STP(stid).N_SN = 0;
#if defined(IMF_SAMPLING_MCS) || defined(SN_MCS_VARIABLE_EJECTA)
      STP(stid).Z_ej = 0;
#endif
    }

  assert(j == nexport); /* Another check that buffers have been filled correctly */

  /* Get export buffer into correct order */
  qsort(SnhostIn, nexport, sizeof(struct snhost_in), snhost_compare);

  /* Exchange particle data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&SnhostIn[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct snhost_in), MPI_BYTE, recvTask,
                           TAG_DENS_A, &SnhostGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct snhost_in), MPI_BYTE,
                           recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
            }
        }
    }

  myfree(SnhostIn);

  /* Now distribute imported SN variables */
  for(i = 0; i < nimport; i++)
    {
      if(P[SnhostGet[i].HostIndex].ID != 0)
        {
          assert(SnhostGet[i].HostTask == ThisTask);   /* Check we have right task */
          assert(P[SnhostGet[i].HostIndex].Type == 0); /* Check we have a gas cell */
          SphP[SnhostGet[i].HostIndex].mass_deposited += SnhostGet[i].mass_deposited;
          SphP[SnhostGet[i].HostIndex].energy_deposited += SnhostGet[i].energy_deposited;
#if defined(IMF_SAMPLING_MCS) || defined(SN_MCS_VARIABLE_EJECTA)
          SphP[SnhostGet[i].HostIndex].metal_deposited += SnhostGet[i].metal_deposited;
#endif
          SphP[SnhostGet[i].HostIndex].starvel[0] += SnhostGet[i].starvel[0];
          SphP[SnhostGet[i].HostIndex].starvel[1] += SnhostGet[i].starvel[1];
          SphP[SnhostGet[i].HostIndex].starvel[2] += SnhostGet[i].starvel[2];
          SphP[SnhostGet[i].HostIndex].N_SN_hosted += 1;
#ifdef SN_MCS_LOG
          sn_add_to_log(SnhostGet[i].HostIndex);
#endif
#ifdef SN_MCS_LOCATION_RECORDS
#ifdef SN_MCS_LOCATION_RECORDS_CHECKS
          SnhostinfoIn[i].OriginTask = SnhostGet[i].OriginTask;
#endif
          SnhostinfoIn[i].OriginStarIndex = SnhostGet[i].OriginStarIndex;
          SnhostinfoIn[i].host_density    = SphP[SnhostGet[i].HostIndex].Density;
#ifdef GRACKLE
          SnhostinfoIn[i].host_temperature = get_temp_individual_cell_grackle(SnhostGet[i].HostIndex);
#endif
#endif
#ifdef SN_MCS_LOG_DETAILS
#ifdef GRACKLE
          double temp = get_temp_individual_cell_grackle(SnhostGet[i].HostIndex);
#else
          double temp = -1;
#endif
          fprintf(FdSNDetails, "SN=%llu %g %g %g %g %g %g %g %g %g %g %g %g %d\n", (long long)SnhostGet[i].StarID, All.Time,
                  SnhostGet[i].StarPos[0], SnhostGet[i].StarPos[1], SnhostGet[i].StarPos[2], SnhostGet[i].starvel[0],
                  SnhostGet[i].starvel[1], SnhostGet[i].starvel[2], SphP[SnhostGet[i].HostIndex].Density, temp,
                  SphP[SnhostGet[i].HostIndex].Metallicity, SnhostGet[i].mass_deposited, SnhostGet[i].energy_deposited, 1);
#endif
        }
    }

  myfree(SnhostGet);

#ifdef SN_MCS_LOCATION_RECORDS
  SnhostinfoGet = (struct snhostinfo_in *)mymalloc("SnhostinfoIn", nexport * sizeof(struct snhostinfo_in));

  /* Exchange particle data, with sending and receiving flipped from previous,
  since all communication will be back to the tasks from which the star particles
  originated i.e. using Recv_count as the send count etc.*/
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&SnhostinfoIn[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct snhostinfo_in), MPI_BYTE,
                           recvTask, TAG_DENS_A, &SnhostinfoGet[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct snhostinfo_in), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                           &status);
            }
        }
    }

  for(i = 0; i < nexport; i++)
    {
#ifdef SN_MCS_LOCATION_RECORDS_CHECKS
      assert(SnhostinfoGet[i].OriginTask == ThisTask);
#endif
      StarP[SnhostinfoGet[i].OriginStarIndex].SNDensity = SnhostinfoGet[i].host_density;
#ifdef GRACKLE
      StarP[SnhostinfoGet[i].OriginStarIndex].SNTemperature = SnhostinfoGet[i].host_temperature;
#endif
    }

  myfree(SnhostinfoGet);
  myfree(SnhostinfoIn);

#endif  // SN_MCS_LOCATION_RECORDS

#ifdef TRACER_MC
  finish_MC_tracer();
#endif

  myfree(starindex);
  myfree(searchdata);
}

/* Sorts by task, then by index */
static int snhost_compare(const void *a, const void *b)
{
  if(((struct snhost_in *)a)->HostTask < (((struct snhost_in *)b)->HostTask))
    return -1;

  if(((struct snhost_in *)a)->HostTask > (((struct snhost_in *)b)->HostTask))
    return +1;

  if(((struct snhost_in *)a)->HostIndex < (((struct snhost_in *)b)->HostIndex))
    return -1;

  if(((struct snhost_in *)a)->HostIndex > (((struct snhost_in *)b)->HostIndex))
    return +1;

  return 0;
}
#endif

void stellar_feedback_distribute(void)
{
  int i, idx, j, k;
  int n_sn_host_on_task, nexport, nimport, ngrp, recvTask;
  int q, dpi, vf, particle, n_neighbour, index;
  double dx, dy, dz;
  double sum_omega_x_plus[3], sum_omega_x_minus[3], factor_plus[3], factor_minus[3], sum_weights;
  double p_ej_tot, mod_weight, dm, u0, u1;
  MPI_Status status;
#ifndef SN_MCS_NO_ENERGY
  double dv[3], kin_energy, dp_sq0, dp[3];
#endif
#ifdef SN_MCS_MECHANICAL
  double n, p_terminal, E_51, Zsun, boost_factor;
#endif
#if defined(GRACKLE) && !defined(GRACKLE_TAB)
  double spec_norm;
#endif

  mpi_printf("SN_MCS: Distributing feedback...\n");

  /* Find number of SNe hosting cells on this task, also compute average of velocity vector */
  n_sn_host_on_task = 0;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
      if(SphP[i].N_SN_hosted > 0)
        {
          SphP[i].starvel[0] /= SphP[i].N_SN_hosted;
          SphP[i].starvel[1] /= SphP[i].N_SN_hosted;
          SphP[i].starvel[2] /= SphP[i].N_SN_hosted;
          n_sn_host_on_task++;
        }
    }

  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

#ifdef TRACER_MC
  start_MC_tracer(N_tracer); /* allocate buffer for tracer exchange */
#endif

  SnvarIn = (struct snvar_in *)mymalloc("SnvarIn", SN_MCS_NNGB_MAX * n_sn_host_on_task * sizeof(struct snvar_in));

  Neighbour_weights = (struct iso_weights *)mymalloc("Neighbour_weights", SN_MCS_NNGB_MAX * sizeof(struct iso_weights));

  nexport = 0;

  /* Find cells hosting SNe and compute weights */
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
      if(SphP[i].N_SN_hosted > 0)
        {
          q           = SphP[i].first_connection;
          n_neighbour = 0;

          while(q >= 0)
            {
              if(n_neighbour == SN_MCS_NNGB_MAX)
                terminate("Number of neighbouring cells too large, increase SN_MCS_NNGB_MAX. n_neighbour %d\n", n_neighbour);

              dpi      = DC[q].dp_index;
              vf       = DC[q].vf_index;
              particle = Mesh.DP[dpi].index;

              if(particle < 0)
                {
                  q = DC[q].next;
                  continue;
                }

              if(Mesh.DP[dpi].task == ThisTask)
                {
                  if(particle >= NumGas && Mesh.DP[dpi].task == ThisTask)
                    particle -= NumGas;

                  Neighbour_weights[n_neighbour].Index = particle;
                }
              else
                Neighbour_weights[n_neighbour].Index = Mesh.DP[dpi].originalindex;

              Neighbour_weights[n_neighbour].Task = Mesh.DP[dpi].task;

              dx = Mesh.DP[dpi].x - P[i].Pos[0];
              dy = Mesh.DP[dpi].y - P[i].Pos[1];
              dz = Mesh.DP[dpi].z - P[i].Pos[2];

              Neighbour_weights[n_neighbour].r = sqrt(dx * dx + dy * dy + dz * dz);

              Neighbour_weights[n_neighbour].x_plus[0] = fmax(dx, 0) / Neighbour_weights[n_neighbour].r;
              Neighbour_weights[n_neighbour].x_plus[1] = fmax(dy, 0) / Neighbour_weights[n_neighbour].r;
              Neighbour_weights[n_neighbour].x_plus[2] = fmax(dz, 0) / Neighbour_weights[n_neighbour].r;

              Neighbour_weights[n_neighbour].x_minus[0] = fmin(dx, 0) / Neighbour_weights[n_neighbour].r;
              Neighbour_weights[n_neighbour].x_minus[1] = fmin(dy, 0) / Neighbour_weights[n_neighbour].r;
              Neighbour_weights[n_neighbour].x_minus[2] = fmin(dz, 0) / Neighbour_weights[n_neighbour].r;

              Neighbour_weights[n_neighbour].omega =
                  0.5 * (1.0 - 1.0 / (sqrt(1.0 + (4.0 * Mesh.VF[vf].area /
                                                  (M_PI * Neighbour_weights[n_neighbour].r * Neighbour_weights[n_neighbour].r)))));

              n_neighbour++;

              if(q == SphP[i].last_connection)
                break;

              q = DC[q].next;
            }

#ifdef SN_MCS_NNGB_MAX_ADVISORY
          if(n_neighbour > SN_MCS_NNGB_MAX_DEFAULT)
            printf(
                "SN_MCS: Warning, SN_MCS_NNGB_MAX_DEFAULT (%d) exceeded on Task %d at Sync-point %d. n_neighbour = %d, ID = %d, pos "
                "%g %g %g rho %g\n",
                SN_MCS_NNGB_MAX_DEFAULT, ThisTask, All.NumCurrentTiStep, n_neighbour, P[i].ID, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2],
                SphP[i].Density);
#endif

          for(k = 0; k < 3; k++)
            {
              sum_omega_x_plus[k]  = 0;
              sum_omega_x_minus[k] = 0;
            }

          for(j = 0; j < n_neighbour; j++)
            {
              for(k = 0; k < 3; k++)
                {
                  sum_omega_x_plus[k] += Neighbour_weights[j].omega * fabs(Neighbour_weights[j].x_plus[k]);
                  sum_omega_x_minus[k] += Neighbour_weights[j].omega * fabs(Neighbour_weights[j].x_minus[k]);
                }
            }

          for(k = 0; k < 3; k++)
            {
              factor_plus[k] =
                  sqrt(0.5 * (1.0 + sum_omega_x_minus[k] * sum_omega_x_minus[k] / (sum_omega_x_plus[k] * sum_omega_x_plus[k])));
              factor_minus[k] =
                  sqrt(0.5 * (1.0 + sum_omega_x_plus[k] * sum_omega_x_plus[k] / (sum_omega_x_minus[k] * sum_omega_x_minus[k])));
            }

          for(j = 0; j < n_neighbour; j++)
            {
              for(k = 0; k < 3; k++)
                {
                  Neighbour_weights[j].w[k] =
                      Neighbour_weights[j].x_plus[k] * factor_plus[k] + Neighbour_weights[j].x_minus[k] * factor_minus[k];
                  Neighbour_weights[j].w[k] *= Neighbour_weights[j].omega;
                }
            }

          sum_weights = 0;

          for(j = 0; j < n_neighbour; j++)
            {
              sum_weights +=
                  sqrt(Neighbour_weights[j].w[0] * Neighbour_weights[j].w[0] + Neighbour_weights[j].w[1] * Neighbour_weights[j].w[1] +
                       Neighbour_weights[j].w[2] * Neighbour_weights[j].w[2]);
            }

          for(j = 0; j < n_neighbour; j++)
            {
              Neighbour_weights[j].w[0] *= (1.0 - All.HostCellFeedbackFraction) / sum_weights;
              Neighbour_weights[j].w[1] *= (1.0 - All.HostCellFeedbackFraction) / sum_weights;
              Neighbour_weights[j].w[2] *= (1.0 - All.HostCellFeedbackFraction) / sum_weights;
#ifdef SN_MCS_WEIGHTS_VERBOSE
              printf("SN_MCS_WEIGHTS_VERBOSE: %d %d %g %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", ThisTask,
                     All.NumCurrentTiStep, All.Time, P[i].ID, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], Neighbour_weights[j].r,
                     factor_plus[0], factor_plus[1], factor_plus[2], factor_minus[0], factor_minus[1], factor_minus[2],
                     Neighbour_weights[j].omega, Neighbour_weights[j].x_plus[0], Neighbour_weights[j].x_plus[1],
                     Neighbour_weights[j].x_plus[2], Neighbour_weights[j].x_minus[0], Neighbour_weights[j].x_minus[1],
                     Neighbour_weights[j].x_minus[2], sum_weights);
#endif
            }

          p_ej_tot = sqrt(2.0 * All.SNKineticRatio * SphP[i].energy_deposited * SphP[i].mass_deposited);
#ifdef SN_MCS_MECHANICAL
          E_51 = (SphP[i].energy_deposited * All.SNKineticRatio) / (All.SupernovaEnergy * All.cf_atime * All.cf_atime);
#endif

          for(j = 0; j < n_neighbour; j++)
            {
              if(Neighbour_weights[j].Task == ThisTask)
                {
                  mod_weight = sqrt(Neighbour_weights[j].w[0] * Neighbour_weights[j].w[0] +
                                    Neighbour_weights[j].w[1] * Neighbour_weights[j].w[1] +
                                    Neighbour_weights[j].w[2] * Neighbour_weights[j].w[2]);

                  index = Neighbour_weights[j].Index;
                  if(P[index].ID == 0)
                    {
                      printf("ID0 NGB INJECT LOCAL: Task %d Sync-point %d index %d ID %d NumGas %d Mass %g\n", ThisTask,
                             All.NumCurrentTiStep, index, P[index].ID, NumGas, P[index].Mass);
                      continue;
                    }

                  u0 = (SphP[index].Energy -
                        0.5 *
                            (SphP[index].Momentum[0] * SphP[index].Momentum[0] + SphP[index].Momentum[1] * SphP[index].Momentum[1] +
                             SphP[index].Momentum[2] * SphP[index].Momentum[2]) /
                            P[index].Mass) /
                       P[index].Mass;
#ifdef MCS_UTHERM_CATCH
                  /*Debug code */
                  if(u0 < 0)
                    {
                      printf(
                          "UTHERM_CATCH_1: u0 < 0: Sync-point %d Time %g Task %d Energy %g Momentum %g %g %g Mass %g u0 %g index %d "
                          "ID %d\n",
                          All.NumCurrentTiStep, All.Time, ThisTask, SphP[index].Energy, SphP[index].Momentum[0],
                          SphP[index].Momentum[1], SphP[index].Momentum[2], P[index].Mass, u0, index, P[index].ID);
                      SphP[index].Energy += (-u0 + All.MinEgySpec * All.cf_atime * All.cf_atime) * P[index].Mass;
                    }
#endif

                  dm = SphP[i].mass_deposited * mod_weight;

#ifndef SN_MCS_NO_ENERGY

                  SphP[index].Energy += mod_weight * SphP[i].energy_deposited;

                  dp[0] = Neighbour_weights[j].w[0] * p_ej_tot;
                  dp[1] = Neighbour_weights[j].w[1] * p_ej_tot;
                  dp[2] = Neighbour_weights[j].w[2] * p_ej_tot;

                  if(dm > 0) /*Guard against machine precision errors, occasionally occurs if omega ~ 0 becuase of distorted cell*/
                    {
                      dp_sq0 = dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2];

                      /* Now move from star rest frame to coordinate frame */
                      dp[0] += dm * SphP[i].starvel[0];
                      dp[1] += dm * SphP[i].starvel[1];
                      dp[2] += dm * SphP[i].starvel[2];

                      SphP[index].Energy += (dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2] - dp_sq0) / (2.0 * dm);

#ifdef SN_MCS_MECHANICAL
                      n = SphP[index].Density * All.cf_a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
                      n *= HYDROGEN_MASSFRAC / PROTONMASS;
                      Zsun       = fmax(SphP[index].Metallicity / 0.0127, 0.01);
                      p_terminal = All.SupernovaTerminalMomentum * pow(E_51, 16.0 / 17.0) * pow(n, -2.0 / 17.0) * pow(Zsun, -0.14) *
                                   All.cf_atime;
                      boost_factor = fmin(sqrt(1.0 + P[index].Mass / dm), (p_terminal / p_ej_tot));
                      if(boost_factor < 1.0)
                        {
                          terminate(
                              "SN_MCS SN_MCS_MECHANICAL: boost_factor < 1, this shouldn't happen! p_terminal %g p_ej_tot %g "
                              "P[index].Mass %g dm %g n %g Zsun %g E_51 %g\n",
                              p_terminal, p_ej_tot / All.cf_atime, P[index].Mass, dm, n, Zsun, E_51);
                        }

                      dp[0] *= boost_factor;
                      dp[1] *= boost_factor;
                      dp[2] *= boost_factor;
#endif
                      SphP[index].Momentum[0] += dp[0];
                      SphP[index].Momentum[1] += dp[1];
                      SphP[index].Momentum[2] += dp[2];
                    }
#endif
#ifndef SN_MCS_INITIAL_DRIVING
                  P[index].Mass += dm;

#ifdef TRACER_MC
                  consider_moving_tracers_ejecta_local(i, index, dm / SphP[i].mass_deposited);
#endif

#ifdef METALS
#ifdef MCS_METAL_ERROR_CATCH
                  /*Debug code*/
                  if((P[index].Metallicity < 0.0) || (SphP[index].Metallicity < 0.0) || (SphP[index].MassMetallicity < 0.0))
                    {
                      mpi_printf(
                          "MCS_METAL_ERROR_CATCH: sn_evaluate Task %d Time %g ID %d P[index].Mass %g P[index].Metallicity %g "
                          "SphP[index].Metallicity %g SphP[index].MassMetallicity %g Density %g Pos %g %g %g\n",
                          ThisTask, All.Time, P[index].ID, P[index].Mass, P[index].Metallicity, SphP[index].Metallicity,
                          SphP[index].MassMetallicity, SphP[index].Density, P[index].Pos[0], P[index].Pos[1], P[index].Pos[2]);
#ifdef MIN_METALLICITY_ON_STARTUP
                      P[index].Metallicity        = All.MinimumMetallicityOnStartUp;
                      SphP[index].Metallicity     = All.MinimumMetallicityOnStartUp;
                      SphP[index].MassMetallicity = P[index].Mass * All.MinimumMetallicityOnStartUp;
#else
                      P[index].Metallicity        = 0.0;
                      SphP[index].Metallicity     = 0.0;
                      SphP[index].MassMetallicity = 0.0;
#endif
                    }
#endif
#if !defined(IMF_SAMPLING_MCS) && !defined(SN_MCS_VARIABLE_EJECTA)
                  SphP[index].MassMetallicity += dm * All.SNEjectaMetallicity;
#else
                  SphP[index].MassMetallicity += SphP[i].metal_deposited * mod_weight;
#endif
                  SphP[index].Metallicity = SphP[index].MassMetallicity / P[index].Mass;

                  assert((SphP[index].Metallicity >= 0) && (SphP[index].Metallicity <= 1.0));
                  P[index].Metallicity = SphP[index].Metallicity;

#if defined(GRACKLE) && !defined(GRACKLE_TAB)
                  spec_norm = P[index].Metallicity;
                  for(int spec = 0; spec < GRACKLE_SPECIES_NUMBER; spec++)
                    spec_norm += SphP[index].GrackleSpeciesFraction[spec];
                  spec_norm += SphP[index].e_frac;

                  for(int spec = 0; spec < GRACKLE_SPECIES_NUMBER; spec++)
                    {
                      SphP[index].GrackleSpeciesFraction[spec] /= spec_norm;
                      SphP[index].GrackleSpeciesMass[spec] = SphP[index].GrackleSpeciesFraction[spec] * P[index].Mass;
                    }
                  SphP[index].e_frac /= spec_norm;
                  SphP[index].e_mass = SphP[index].e_frac * P[i].Mass;
#endif
#endif
#endif

#ifdef SN_MCS_CELL_VERBOSE
                  printf("SN_MCS_CELL_VERBOSE: %d %d %g %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g\n", ThisTask,
                         All.NumCurrentTiStep, All.Time, P[i].ID, SphP[i].N_SN_hosted, n_neighbour, j + 1, Neighbour_weights[j].Task,
                         P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], SphP[i].starvel[0], SphP[i].starvel[1], SphP[i].starvel[2], dm,
                         mod_weight * SphP[i].energy_deposited * (1.0 - All.SNKineticRatio), dp[0], dp[1], dp[2],
                         Neighbour_weights[j].r);
#endif
                  u1 = (SphP[index].Energy -
                        0.5 *
                            (SphP[index].Momentum[0] * SphP[index].Momentum[0] + SphP[index].Momentum[1] * SphP[index].Momentum[1] +
                             SphP[index].Momentum[2] * SphP[index].Momentum[2]) /
                            P[index].Mass) /
                       P[index].Mass;
                  if(u1 < u0)
                    SphP[index].Energy = u0 * P[index].Mass + 0.5 *
                                                                  (SphP[index].Momentum[0] * SphP[index].Momentum[0] +
                                                                   SphP[index].Momentum[1] * SphP[index].Momentum[1] +
                                                                   SphP[index].Momentum[2] * SphP[index].Momentum[2]) /
                                                                  P[index].Mass;
                }
              else
                {
                  mod_weight = sqrt(Neighbour_weights[j].w[0] * Neighbour_weights[j].w[0] +
                                    Neighbour_weights[j].w[1] * Neighbour_weights[j].w[1] +
                                    Neighbour_weights[j].w[2] * Neighbour_weights[j].w[2]);

                  dm = SphP[i].mass_deposited * mod_weight;

#ifndef SN_MCS_NO_ENERGY

                  SnvarIn[nexport].dE = mod_weight * SphP[i].energy_deposited;

                  dp[0] = Neighbour_weights[j].w[0] * p_ej_tot;
                  dp[1] = Neighbour_weights[j].w[1] * p_ej_tot;
                  dp[2] = Neighbour_weights[j].w[2] * p_ej_tot;

                  if(dm > 0) /*Guard against machine precision errors, occasionally occurs if omega ~ 0 becuase of distorted cell*/
                    {
                      dp_sq0 = dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2];

                      /* Now move from star rest frame to coordinate frame */
                      dp[0] += dm * SphP[i].starvel[0];
                      dp[1] += dm * SphP[i].starvel[1];
                      dp[2] += dm * SphP[i].starvel[2];

                      SnvarIn[nexport].dE += (dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2] - dp_sq0) / (2.0 * dm);
                    }

                  SnvarIn[nexport].dp[0] = dp[0];
                  SnvarIn[nexport].dp[1] = dp[1];
                  SnvarIn[nexport].dp[2] = dp[2];

#ifdef SN_MCS_MECHANICAL
                  SnvarIn[nexport].p_ej_tot = p_ej_tot;
                  SnvarIn[nexport].E_51     = E_51;
#endif
#endif

                  SnvarIn[nexport].dm = dm;
#if defined(IMF_SAMPLING_MCS) || defined(SN_MCS_VARIABLE_EJECTA)
                  SnvarIn[nexport].dmZ = SphP[i].metal_deposited * mod_weight;
#endif
                  SnvarIn[nexport].Task = Neighbour_weights[j].Task;
                  Send_count[Neighbour_weights[j].Task] += 1;
                  SnvarIn[nexport].Index = Neighbour_weights[j].Index;
#ifdef TRACER_MC
                  consider_moving_tracers_ejecta(i, Neighbour_weights[j].Task, Neighbour_weights[j].Index,
                                                 dm / SphP[i].mass_deposited);
#endif
                  nexport++;
#ifdef SN_MCS_CELL_VERBOSE
                  printf("SN_MCS_CELL_VERBOSE: %d %d %g %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g\n", ThisTask,
                         All.NumCurrentTiStep, All.Time, P[i].ID, SphP[i].N_SN_hosted, n_neighbour, j + 1, Neighbour_weights[j].Task,
                         P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], SphP[i].starvel[0], SphP[i].starvel[1], SphP[i].starvel[2], dm,
                         mod_weight * SphP[i].energy_deposited * (1.0 - All.SNKineticRatio), dp[0], dp[1], dp[2],
                         Neighbour_weights[j].r);
#endif
                }
            }

          /* Give host cell energy, mass and metals */
          if(P[i].ID != 0)
            {
#ifndef SN_MCS_INITIAL_DRIVING
              P[i].Mass += SphP[i].mass_deposited * All.HostCellFeedbackFraction;
#endif
#ifndef SN_MCS_NO_ENERGY
              SphP[i].Energy += SphP[i].energy_deposited * All.HostCellFeedbackFraction;
              dv[0] = SphP[i].starvel[0] - P[i].Vel[0];
              dv[1] = SphP[i].starvel[1] - P[i].Vel[1];
              dv[2] = SphP[i].starvel[2] - P[i].Vel[2];
              kin_energy =
                  All.HostCellFeedbackFraction * 0.5 * SphP[i].mass_deposited * (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
              SphP[i].Energy += kin_energy;
#endif

#ifndef SN_MCS_INITIAL_DRIVING
#ifdef METALS
#ifdef MCS_METAL_ERROR_CATCH
              /*Debug code*/
              if((P[i].Metallicity < 0.0) || (SphP[i].Metallicity < 0.0) || (SphP[i].MassMetallicity < 0.0))
                {
                  mpi_printf(
                      "MCS_METAL_ERROR_CATCH: sn_evaluate Task %d Time %g ID %d P[i].Mass %g P[i].Metallicity %g SphP[i].Metallicity "
                      "%g "
                      "SphP[i].MassMetallicity %g Density %g Pos %g %g %g\n",
                      ThisTask, All.Time, P[i].ID, P[i].Mass, P[i].Metallicity, SphP[i].Metallicity, SphP[i].MassMetallicity,
                      SphP[i].Density, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
#ifdef MIN_METALLICITY_ON_STARTUP
                  P[i].Metallicity        = All.MinimumMetallicityOnStartUp;
                  SphP[i].Metallicity     = All.MinimumMetallicityOnStartUp;
                  SphP[i].MassMetallicity = P[i].Mass * All.MinimumMetallicityOnStartUp;
#else
                  P[i].Metallicity = 0.0;
                  SphP[i].Metallicity = 0.0;
                  SphP[i].MassMetallicity = 0.0;
#endif
                }
#endif
#if !defined(IMF_SAMPLING_MCS) && !defined(SN_MCS_VARIABLE_EJECTA)
              SphP[i].MassMetallicity += SphP[i].mass_deposited * All.HostCellFeedbackFraction * All.SNEjectaMetallicity;
#else
              SphP[i].MassMetallicity += SphP[i].metal_deposited * All.HostCellFeedbackFraction;
#endif
              SphP[i].Metallicity = SphP[i].MassMetallicity / P[i].Mass;
              assert((SphP[i].Metallicity >= 0) && (SphP[i].Metallicity <= 1.0));
              P[i].Metallicity = SphP[i].Metallicity;
#endif
#if defined(GRACKLE) && !defined(GRACKLE_TAB)
              spec_norm = P[i].Metallicity;
              for(int spec = 0; spec < GRACKLE_SPECIES_NUMBER; spec++)
                spec_norm += SphP[i].GrackleSpeciesFraction[spec];
              spec_norm += SphP[i].e_frac;

              for(int spec = 0; spec < GRACKLE_SPECIES_NUMBER; spec++)
                {
                  SphP[i].GrackleSpeciesFraction[spec] /= spec_norm;
                  SphP[i].GrackleSpeciesMass[spec] = SphP[i].GrackleSpeciesFraction[spec] * P[i].Mass;
                }
              SphP[i].e_frac /= spec_norm;
              SphP[i].e_mass = SphP[i].e_frac * P[i].Mass;
#endif
#endif

#ifdef SN_MCS_CELL_VERBOSE
              printf("SN_MCS_CELL_VERBOSE: %d %d %g %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g\n", ThisTask,
                     All.NumCurrentTiStep, All.Time, P[i].ID, SphP[i].N_SN_hosted, n_neighbour, 0, ThisTask, P[i].Pos[0], P[i].Pos[1],
                     P[i].Pos[2], SphP[i].starvel[0], SphP[i].starvel[1], SphP[i].starvel[2],
                     SphP[i].mass_deposited * All.HostCellFeedbackFraction, SphP[i].energy_deposited * All.HostCellFeedbackFraction,
                     0.0, 0.0, 0.0, 0.0);
#endif
            }
          else
            {
              /*Debug code*/
              printf("ID0 NGB HOSTINJECT LOCAL: Task %d Sync-point %d index %d ID %d NumGas %d Mass %g\n", ThisTask,
                     All.NumCurrentTiStep, i, P[i].ID, NumGas, P[i].Mass);
            }

          /* Reset host quantities */
          SphP[i].N_SN_hosted      = 0;
          SphP[i].mass_deposited   = 0;
          SphP[i].energy_deposited = 0;
#if defined(IMF_SAMPLING_MCS) || defined(SN_MCS_VARIABLE_EJECTA)
          SphP[i].metal_deposited = 0;
#endif
          SphP[i].starvel[0] = 0;
          SphP[i].starvel[1] = 0;
          SphP[i].starvel[2] = 0;
#ifdef TRACER_MC
          finalise_host_tracers(i);
#endif
        }
    }

  myfree(Neighbour_weights);

  /* Get export buffer into correct order */
  qsort(SnvarIn, nexport, sizeof(struct snvar_in), snvar_compare);

  /* Send the export counts to other tasks and get receive counts */
  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(i = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; i < NTask; i++)
    {
      nimport += Recv_count[i];

      if(i > 0)
        {
          Send_offset[i] = Send_offset[i - 1] + Send_count[i - 1];
          Recv_offset[i] = Recv_offset[i - 1] + Recv_count[i - 1];
        }
    }

  SnvarGet = (struct snvar_in *)mymalloc("SnvarGet", nimport * sizeof(struct snvar_in));

  /* Exchange particle data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&SnvarIn[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct snvar_in), MPI_BYTE, recvTask,
                           TAG_DENS_A, &SnvarGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct snvar_in), MPI_BYTE,
                           recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
            }
        }
    }

  /* Distribute imported quantitites */
  for(i = 0; i < nimport; i++)
    {
      assert(SnvarGet[i].Task == ThisTask); /* Check we have right task */
      index = SnvarGet[i].Index;
      if(index > NumGas)
        index -= NumGas;
      assert(P[index].Type == 0); /* Check we have a gas cell */
      if(P[index].ID == 0)
        {
          /*Debug code*/
          printf("ID0 NGB INJECT LOCAL: Task %d Sync-point %d index %d ID %d NumGas %d Mass %g\n", ThisTask, All.NumCurrentTiStep,
                 index, P[index].ID, NumGas, P[index].Mass);
          continue;
        }

      u0 = (SphP[index].Energy -
            0.5 *
                (SphP[index].Momentum[0] * SphP[index].Momentum[0] + SphP[index].Momentum[1] * SphP[index].Momentum[1] +
                 SphP[index].Momentum[2] * SphP[index].Momentum[2]) /
                P[index].Mass) /
           P[index].Mass;
#ifdef MCS_UTHERM_CATCH
      /*Debug code */
      if(u0 < 0)
        {
          printf("UTHERM_CATCH_2: u0 < 0: Sync-point %d Time %g Task %d Energy %g Momentum %g %g %g Mass %g u0 %g index %d ID %d\n",
                 All.NumCurrentTiStep, All.Time, ThisTask, SphP[index].Energy, SphP[index].Momentum[0], SphP[index].Momentum[1],
                 SphP[index].Momentum[2], P[index].Mass, u0, index, P[index].ID);
          SphP[index].Energy += (-u0 + All.MinEgySpec * All.cf_atime * All.cf_atime) * P[index].Mass;
        }
#endif

#ifndef SN_MCS_NO_ENERGY
      SphP[index].Energy += SnvarGet[i].dE;

#ifdef SN_MCS_MECHANICAL
      if(SnvarGet[i].dm > 0) /*Guard against machine precision errors, occasionally occurs if omega ~ 0 becuase of distorted cell*/
        {
          n = SphP[index].Density * All.cf_a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
          n *= HYDROGEN_MASSFRAC / PROTONMASS;
          Zsun       = fmax(SphP[index].Metallicity / 0.0127, 0.01);
          p_terminal = All.SupernovaTerminalMomentum * pow(SnvarGet[i].E_51, 16.0 / 17.0) * pow(n, -2.0 / 17.0) * pow(Zsun, -0.14) *
                       All.cf_atime;
          boost_factor = fmin(sqrt(1.0 + P[index].Mass / SnvarGet[i].dm), (p_terminal / SnvarGet[i].p_ej_tot));
          if(boost_factor < 1.0)
            {
              terminate(
                  "SN_MCS SN_MCS_MECHANICAL: boost_factor < 1, this shouldn't happen! p_terminal %g p_ej_tot %g P[index].Mass %g dm "
                  "%g n %g Zsun %g E_51 %g\n",
                  p_terminal, SnvarGet[i].p_ej_tot / All.cf_atime, P[index].Mass, SnvarGet[i].dm, n, Zsun, SnvarGet[i].E_51);
            }
          SnvarGet[i].dp[0] *= boost_factor;
          SnvarGet[i].dp[1] *= boost_factor;
          SnvarGet[i].dp[2] *= boost_factor;
        }
#endif
      SphP[index].Momentum[0] += SnvarGet[i].dp[0];
      SphP[index].Momentum[1] += SnvarGet[i].dp[1];
      SphP[index].Momentum[2] += SnvarGet[i].dp[2];

#endif

#ifndef SN_MCS_INITIAL_DRIVING
      P[index].Mass += SnvarGet[i].dm;

#ifdef METALS
#ifdef MCS_METAL_ERROR_CATCH
      /*Debug code*/
      if((P[index].Metallicity < 0.0) || (SphP[index].Metallicity < 0.0) || (SphP[index].MassMetallicity < 0.0))
        {
          mpi_printf(
              "MCS_METAL_ERROR_CATCH: sn_evaluate Task %d Time %g ID %d P[i].Mass %g P[i].Metallicity %g SphP[i].Metallicity %g "
              "SphP[i].MassMetallicity %g Density %g Pos %g %g %g\n",
              ThisTask, All.Time, P[i].ID, P[index].Mass, P[index].Metallicity, SphP[index].Metallicity, SphP[index].MassMetallicity,
              SphP[index].Density, P[index].Pos[0], P[index].Pos[1], P[index].Pos[2]);
#ifdef MIN_METALLICITY_ON_STARTUP
          P[index].Metallicity        = All.MinimumMetallicityOnStartUp;
          SphP[index].Metallicity     = All.MinimumMetallicityOnStartUp;
          SphP[index].MassMetallicity = P[index].Mass * All.MinimumMetallicityOnStartUp;
#else
          P[index].Metallicity = 0.0;
          SphP[index].Metallicity = 0.0;
          SphP[index].MassMetallicity = 0.0;
#endif
        }
#endif
#if !defined(IMF_SAMPLING_MCS) && !defined(SN_MCS_VARIABLE_EJECTA)
      SphP[index].MassMetallicity += SnvarGet[i].dm * All.SNEjectaMetallicity;
#else
      SphP[index].MassMetallicity += SnvarGet[i].dmZ;
#endif
      SphP[index].Metallicity = SphP[index].MassMetallicity / P[index].Mass;
      assert((SphP[index].Metallicity >= 0) && (SphP[index].Metallicity <= 1.0));
      P[index].Metallicity = SphP[index].Metallicity;

#if defined(GRACKLE) && !defined(GRACKLE_TAB)
      spec_norm = P[index].Metallicity;
      for(int spec = 0; spec < GRACKLE_SPECIES_NUMBER; spec++)
        spec_norm += SphP[index].GrackleSpeciesFraction[spec];
      spec_norm += SphP[index].e_frac;

      for(int spec = 0; spec < GRACKLE_SPECIES_NUMBER; spec++)
        {
          SphP[index].GrackleSpeciesFraction[spec] /= spec_norm;
          SphP[index].GrackleSpeciesMass[spec] = SphP[index].GrackleSpeciesFraction[spec] * P[index].Mass;
        }
      SphP[index].e_frac /= spec_norm;
      SphP[index].e_mass = SphP[index].e_frac * P[i].Mass;
#endif
#endif
#endif

      u1 = (SphP[index].Energy -
            0.5 *
                (SphP[index].Momentum[0] * SphP[index].Momentum[0] + SphP[index].Momentum[1] * SphP[index].Momentum[1] +
                 SphP[index].Momentum[2] * SphP[index].Momentum[2]) /
                P[index].Mass) /
           P[index].Mass;
      if(u1 < u0)
        SphP[index].Energy = u0 * P[index].Mass + 0.5 *
                                                      (SphP[index].Momentum[0] * SphP[index].Momentum[0] +
                                                       SphP[index].Momentum[1] * SphP[index].Momentum[1] +
                                                       SphP[index].Momentum[2] * SphP[index].Momentum[2]) /
                                                      P[index].Mass;
    }

  myfree(SnvarGet);
  myfree(SnvarIn);

#ifdef TRACER_MC
  finish_MC_tracer();
#endif

#ifdef SN_MCS_LOG_DETAILS
  myflush(FdSNDetails);
#endif

  update_primitive_variables();
}

/* Sorts by task, then by index */
static int snvar_compare(const void *a, const void *b)
{
  if(((struct snvar_in *)a)->Task < (((struct snvar_in *)b)->Task))
    return -1;

  if(((struct snvar_in *)a)->Task > (((struct snvar_in *)b)->Task))
    return +1;

  if(((struct snvar_in *)a)->Index < (((struct snvar_in *)b)->Index))
    return -1;

  if(((struct snvar_in *)a)->Index > (((struct snvar_in *)b)->Index))
    return +1;

  return 0;
}

#ifdef TRACER_MC
/*! \brief Monte Carlo exchange of ejecta tracers from a local cell/particle to a particle/cell that may be local or remote.
 *
 * \param p P_index of original parent.
 * \param pother_task Task number where pother_ID potential parent exists.
 * \param pother_index Pindex of the potential parent on the task given by pother_task.
 * \param prob probability whether or not to move each child tracer.
 */
int consider_moving_tracers_ejecta(int p, int pother_task, int pother_index, double prob)
{
  int nmoved = 0;
  int next   = P[p].TracerHead;

  if(prob < 0.0 || prob > 1.0)
    terminate("TRACER_MC: bad prob=%g\n", prob);

  while(next >= 0)
    {
      int move = next;
      next     = TracerLinkedList[next].Next;

      /* We are only considering tracers that just arrived as ejecta */
      if(TracerLinkedList[move].EjectaFlag == 0)
        continue;

      double rndnum = get_random_number();

      if(rndnum < prob)
        {
          TracerLinkedList[move].EjectaFlag = 0;
#if defined(TRACER_MC_NUM_FLUID_QUANTITIES) && (TRACER_MC_EJECTA_COUNTER)
          TracerLinkedList[move].fluid_quantities[TracerMCEjectaCounterIndex] += 1;
#endif
          add_tracer_to_TFL(p, move, pother_task, pother_index, 0);
          nmoved++;
        }
    }

  return nmoved;
}

/*! \brief Monte Carlo exchange of ejecta tracers between two cells/particles that are certainly local.
 *
 *  \param p_from P_index to move tracers from
 *  \param p_to P_index to move tracers to
 *  \param prob probability whether or not to move each child tracer
 */
int consider_moving_tracers_ejecta_local(int p_from, int p_to, double prob)
{
  int Nmoved = 0;
  int next   = P[p_from].TracerHead;

  if(prob < 0.0 || prob > 1.0)
    terminate("TRACER_MC: bad prob=%g\n", prob);

  while(next >= 0)
    {
      int move = next;
      next     = TracerLinkedList[next].Next;

      /* We are only considering tracers that just arrived as ejecta */
      if(TracerLinkedList[move].EjectaFlag == 0)
        continue;

      double rndnum = get_random_number();

      if(rndnum < prob) /* move tracer */
        {
          TracerLinkedList[move].EjectaFlag = 0;
#if defined(TRACER_MC_NUM_FLUID_QUANTITIES) && (TRACER_MC_EJECTA_COUNTER)
          TracerLinkedList[move].fluid_quantities[TracerMCEjectaCounterIndex] += 1;
#endif
          move_tracer_between_parents(p_from, p_to, move);
          Nmoved++;
        }
    }

  return Nmoved;
}

void finalise_host_tracers(int p)
{
  int next = P[p].TracerHead;

  while(next >= 0)
    {
      int move = next;
      next     = TracerLinkedList[next].Next;

      if(TracerLinkedList[move].EjectaFlag == 1)
        {
          TracerLinkedList[move].EjectaFlag = 0;
#if defined(TRACER_MC_NUM_FLUID_QUANTITIES) && (TRACER_MC_EJECTA_COUNTER)
          TracerLinkedList[move].fluid_quantities[TracerMCEjectaCounterIndex] += 1;
#endif
        }
    }
}
#endif  // TRACER_MC

#endif
