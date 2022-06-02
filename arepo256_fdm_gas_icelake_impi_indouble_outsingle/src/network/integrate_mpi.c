/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/network/integrate_mpi.c
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

#include <stdlib.h>
#include "string.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "../allvars.h"
#include "../proto.h"
#include "integrate.h"
#include "network_solver.h"

struct network_cell_data
{
  int origTask;
  int origIndex;
  double Density;
  double Utherm;
  double Temperature;
  double Composition[EOS_NSPECIES];
  double dt;
  double dUtherm;
  double RunTime;
  double Sign;
};

struct network_work_data
{
  int task;
  int index;
  double RunTime;
};

struct network_timing_data
{
  double RunTimeMin;
  double RunTimeMax;
  double RunTimeTot;
  int NCells;
};

static int work_data_local_sort(const void *a, const void *b)
{
  if(((struct network_work_data *)a)->RunTime > ((struct network_work_data *)b)->RunTime)
    return -1;

  if(((struct network_work_data *)a)->RunTime < ((struct network_work_data *)b)->RunTime)
    return +1;

  return 0;
}

static int cell_data_local_sort(const void *a, const void *b)
{
  if(((struct network_cell_data *)a)->origTask > ((struct network_cell_data *)b)->origTask)
    return +1;

  if(((struct network_cell_data *)a)->origTask < ((struct network_cell_data *)b)->origTask)
    return -1;

  return 0;
}

static void network_timing_init(struct network_timing_data *td)
{
  td->RunTimeMin = MAX_DOUBLE_NUMBER;
  td->RunTimeMax = 0;
  td->RunTimeTot = 0;
  td->NCells     = 0;
};

static void network_timing_add_cell(struct network_timing_data *td, double RunTime)
{
  if(RunTime < td->RunTimeMin)
    td->RunTimeMin = RunTime;
  if(RunTime > td->RunTimeMax)
    td->RunTimeMax = RunTime;
  td->RunTimeTot += RunTime;
  td->NCells++;
}

static void network_timing_summarise(struct network_timing_data *td)
{
  double RunTimeMin;
  double RunTimeMax;
  double RunTimeTot;
  long long NCells;

  MPI_Reduce(&td->RunTimeMin, &RunTimeMin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&td->RunTimeMax, &RunTimeMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&td->RunTimeTot, &RunTimeTot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  sumup_large_ints(1, &td->NCells, &NCells);

  mpi_printf("NUCLEAR NETWORK: slowest cell=%gs, average=%g\n", RunTimeMax, RunTimeTot / NCells);
}

static void network_do_single_cell_mpi(struct network_cell_data *cell, double dt_cell, double *NewComposition);

static double nuclear_network_root_mpi(double dt_cell, void *params)
{
  struct network_cell_data *cell = (struct network_cell_data *)params;

  double NewComposition[EOS_NSPECIES];
  network_do_single_cell_mpi(cell, dt_cell, NewComposition);

  struct eos_result res;
  double Temperature = cell->Temperature;
  eos_calc_egiven(cell->Density, NewComposition, cell->Utherm + cell->dUtherm, &Temperature, &res);

  double dlnT = log(Temperature) - log(cell->Temperature);

  return dlnT - cell->Sign * All.NuclearNetworkMaxTempChange;
}

void network_do_cells_mpi(int *DoNetwork)
{
  CPU_Step[CPU_MISC] += measure_time();

  double integrate_start = second();

  int CountLocal = 0;
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      if(DoNetwork[idx])
        {
          int i = TimeBinsHydro.ActiveParticleList[idx];
          if(P[i].TimeBinHydro == 0)
            {
              DoNetwork[idx] = 0;
            }
          else
            {
              CountLocal++;
            }
        }
    }

  long long CountTot;
  sumup_large_ints(1, &CountLocal, &CountTot);
  if(CountTot == 0)
    {
      mpi_printf("NUCLEAR NETWORK: Nothing do do!\n");
      return;
    }

  /* collect stuff on the first core of each node, they will only do communication, the other cores will do the work */
  MPI_Comm COMM_NODE;
  char node_name[MPI_MAX_PROCESSOR_NAME];
  int node_name_length, i;
  int node_colour = 0;
  MPI_Get_processor_name(node_name, &node_name_length);
  for(i = 0; i < node_name_length; i++)
    node_colour += (i + 1) * node_name[i];

  int ThisTask_node;
  int NTask_node;
  MPI_Comm_split(MPI_COMM_WORLD, node_colour, 0, &COMM_NODE);
  MPI_Comm_rank(COMM_NODE, &ThisTask_node);
  MPI_Comm_size(COMM_NODE, &NTask_node);
  MPI_Barrier(MPI_COMM_WORLD);

  struct network_cell_data *data_local =
      (struct network_cell_data *)mymalloc("data_local", CountLocal * sizeof(struct network_cell_data));

  int Count = 0;
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      if(DoNetwork[idx])
        {
          int i = TimeBinsHydro.ActiveParticleList[idx];

          data_local[Count].Density     = SphP[i].Density;
          data_local[Count].Utherm      = SphP[i].Utherm;
          data_local[Count].Temperature = SphP[i].EOSTemperature;
          for(int k = 0; k < EOS_NSPECIES; k++)
            data_local[Count].Composition[k] = SphP[i].Composition[k];

          data_local[Count].dt      = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
          data_local[Count].RunTime = SphP[i].NetworkRunTime;

          data_local[Count].origTask  = ThisTask;
          data_local[Count].origIndex = i;
          Count++;
        }
    }

  assert(Count == CountLocal);

  int CountNode;
  MPI_Reduce(&CountLocal, &CountNode, 1, MPI_INT, MPI_SUM, 0, COMM_NODE);

  int DataInCount[NTask_node];
  MPI_Gather(&CountLocal, 1, MPI_INT, DataInCount, 1, MPI_INT, 0, COMM_NODE);

  int DataInOffset[NTask_node];
  struct network_cell_data *data_node;
  if(ThisTask_node == 0)
    {
      data_node = (struct network_cell_data *)mymalloc("data_node", CountNode * sizeof(struct network_cell_data));

      DataInOffset[0] = 0;
      for(int i = 0; i < NTask_node; i++)
        {
          DataInCount[i] *= sizeof(struct network_cell_data);
          if(i > 0)
            DataInOffset[i] = DataInOffset[i - 1] + DataInCount[i - 1];
        }
    }

  MPI_Gatherv(data_local, CountLocal * sizeof(struct network_cell_data), MPI_BYTE, data_node, DataInCount, DataInOffset, MPI_BYTE, 0,
              COMM_NODE);

  int TasksDone = 1;
  int CountTask[NTask];
  int CountTaskAll                          = 0;
  struct network_work_data *work_data_local = NULL;
  struct network_work_data *work_data_all   = NULL;

  int IsWorker = 0;
  if((ThisTask_node > 0 || CountNode == 0) && (ThisTask > 0))
    IsWorker = 1;

  int WorkerCount;
  MPI_Reduce(&IsWorker, &WorkerCount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  mpi_printf("NUCLEAR NETWORK: Using %d tasks for work and %d for data management.\n", WorkerCount, NTask - WorkerCount);

  if(!IsWorker)
    {
      work_data_local = (struct network_work_data *)mymalloc("work_data_local", CountNode * sizeof(struct network_work_data));
      for(int i = 0; i < CountNode; i++)
        {
          work_data_local[i].task    = ThisTask;
          work_data_local[i].index   = i;
          work_data_local[i].RunTime = data_node[i].RunTime;
        }

      if(ThisTask == 0)
        {
          CountTask[0] = CountNode;
          CountTaskAll = CountNode;
          for(int task = 1; task < NTask; task++)
            {
              MPI_Status status;
              MPI_Recv(&CountTask[task], 1, MPI_INT, task, 666000, MPI_COMM_WORLD, &status);

              CountTaskAll += CountTask[task];
            }

          printf("NUCLEAR NETWORK: Collecting %d cells for sorting.\n", CountTaskAll);
          work_data_all = (struct network_work_data *)mymalloc("work_data_all", CountTaskAll * sizeof(struct network_work_data));

          int Offset = 0;
          for(int task = 0; task < NTask; task++)
            {
              if(CountTask[task] == 0)
                continue;

              if(task > 0)
                TasksDone++;

              if(task == 0)
                {
                  memcpy(work_data_all, work_data_local, CountTask[task] * sizeof(struct network_work_data));
                  Offset += CountTask[task];
                  continue;
                }

              MPI_Recv(&work_data_all[Offset], CountTask[task] * sizeof(struct network_work_data), MPI_BYTE, task, 666001,
                       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
              Offset += CountTask[task];
            }

          mysort(work_data_all, CountTaskAll, sizeof(struct network_work_data), work_data_local_sort);
        }
      else
        {
          MPI_Send(&CountNode, 1, MPI_INT, 0, 666000, MPI_COMM_WORLD);
          MPI_Send(work_data_local, CountNode * sizeof(struct network_work_data), MPI_BYTE, 0, 666001, MPI_COMM_WORLD);
        }
    }
  else
    {
      CountNode = 0;
      MPI_Send(&CountNode, 1, MPI_INT, 0, 666000, MPI_COMM_WORLD);
    }

  int CellWorkIndex = 0;  // only relevant on Task 0

  struct network_cell_data *cell_data = (struct network_cell_data *)mymalloc_movable(&cell_data, "cell_data", 0);
  int cell_count                      = 0;

  struct network_timing_data td;
  network_timing_init(&td);

  mpi_printf("NUCLEAR NETWORK: starting actual computation...\n");

  while(1)
    {
      if(!IsWorker)
        {
          /* check for any new requests */
          MPI_Status status;
          MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

          if(status.MPI_TAG == 666002)
            {
              if(ThisTask != 0)
                terminate("This should not happen!");

              int ping;
              MPI_Recv(&ping, 1, MPI_INT, status.MPI_SOURCE, 666002, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
              if(ping != -1)
                terminate("This should not happen either!");

              if(CellWorkIndex < CountTaskAll)
                {
                  // Request for new cell to work on, just send the top of the list
                  MPI_Send(&work_data_all[CellWorkIndex], sizeof(struct network_work_data), MPI_BYTE, status.MPI_SOURCE, 666003,
                           MPI_COMM_WORLD);
                  CellWorkIndex++;
                }
              else
                {
                  // tell it to stop working
                  struct network_work_data data_done;
                  data_done.task    = -1;
                  data_done.index   = -1;
                  data_done.RunTime = 0;

                  MPI_Send(&data_done, sizeof(struct network_work_data), MPI_BYTE, status.MPI_SOURCE, 666003, MPI_COMM_WORLD);

                  TasksDone++;

                  if(TasksDone == NTask)
                    {
                      /* send all other data tasks a message that we are done */
                      for(int task = 1; task < NTask; task++)
                        {
                          if(CountTask[task] > 0)
                            {
                              int ping = -1;
                              MPI_Send(&ping, 1, MPI_INT, task, 666006, MPI_COMM_WORLD);
                            }
                        }

                      // we are done as well!
                      break;
                    }
                }
            }
          else if(status.MPI_TAG == 666004)
            {
              /* a request for cell data */
              int idx;
              MPI_Recv(&idx, 1, MPI_INT, status.MPI_SOURCE, 666004, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
              MPI_Send(&data_node[idx], sizeof(struct network_cell_data), MPI_BYTE, status.MPI_SOURCE, 666005, MPI_COMM_WORLD);
            }
          else if(status.MPI_TAG == 666006)
            {
              if(status.MPI_SOURCE != 0)
                terminate("This should not happen!");

              int ping;
              MPI_Recv(&ping, 1, MPI_INT, status.MPI_SOURCE, 666006, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
              break;
            }
        }
      else
        {
          /* ping master task to figure out what we should do */
          int ping = -1;
          MPI_Send(&ping, 1, MPI_INT, 0, 666002, MPI_COMM_WORLD);

          struct network_work_data work_data;
          MPI_Recv(&work_data, sizeof(struct network_work_data), MPI_BYTE, 0, 666003, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

          /* check if there is any work left for us */
          if(work_data.task < 0)
            break;

          cell_data = (struct network_cell_data *)myrealloc_movable(cell_data, (cell_count + 1) * sizeof(struct network_cell_data));

          /* get cell data that we should work on */
          MPI_Send(&work_data.index, 1, MPI_INT, work_data.task, 666004, MPI_COMM_WORLD);
          MPI_Recv(&cell_data[cell_count], sizeof(struct network_cell_data), MPI_BYTE, work_data.task, 666005, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);

          /* run nuclear network */
          double network_start = second();

#ifdef NUCLEAR_NETWORK_LIMIT_COMPOSITION_CHANGE
          double NewComposition[EOS_NSPECIES];
          network_do_single_cell_mpi(&cell_data[cell_count], cell_data[cell_count].dt_cell, NewComposition);

          double dComposition[EOS_NSPECIES];
          double fac = 1.0;

          for(int k = 0; k < EOS_NSPECIES; k++)
            {
              dComposition[k] = NewComposition[k] - cell_data[cell_count].Composition[k];
              if(fabs(dComposition[k]) > All.NuclearNetworkMaxCompositionChange)
                fac = fmin(fac, All.NuclearNetworkMaxCompositionChange / fabs(dComposition[k]));
            }

          for(int k = 0; k < EOS_NSPECIES; k++)
            cell_data[cell_count].Composition[k] += fac * dComposition[k];
          cell_data[cell_count].dUtherm *= fac;
#else
          double NewComposition[EOS_NSPECIES];
          network_do_single_cell_mpi(&cell_data[cell_count], cell_data[cell_count].dt, NewComposition);

          struct eos_result res;
          double Temperature = cell_data[cell_count].Temperature;
          eos_calc_egiven(cell_data[cell_count].Density, NewComposition, cell_data[cell_count].Utherm + cell_data[cell_count].dUtherm,
                          &Temperature, &res);

          double dlnT = log(Temperature) - log(cell_data[cell_count].Temperature);

          if(fabs(dlnT) > All.NuclearNetworkMaxTempChange)
            {
              gsl_function F;
              F.function = nuclear_network_root_mpi;

              if(dlnT < 0)
                cell_data[cell_count].Sign = -1;
              else
                cell_data[cell_count].Sign = +1;
              F.params = &cell_data[cell_count];

              const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
              gsl_root_fsolver *s            = gsl_root_fsolver_alloc(T);
              gsl_root_fsolver_set(s, &F, cell_data[cell_count].dt * 1e-10, cell_data[cell_count].dt);
              // gsl_root_fsolver_set(s, &F, 0., cell_data[cell_count].dt);

              int status;
              int iter = 0;
              do
                {
                  iter++;
                  status                   = gsl_root_fsolver_iterate(s);
                  cell_data[cell_count].dt = gsl_root_fsolver_root(s);
                  double x_lo              = gsl_root_fsolver_x_lower(s);
                  double x_hi              = gsl_root_fsolver_x_upper(s);
                  status                   = gsl_root_test_interval(x_lo, x_hi, 0, 0.01);
                }
              while(status == GSL_CONTINUE && iter < 1000);

              if(iter == 1000)
                terminate("Network not converged for cell count=%d", cell_count);

              gsl_root_fsolver_free(s);

              network_do_single_cell_mpi(&cell_data[cell_count], cell_data[cell_count].dt, NewComposition);
            }

          for(int k = 0; k < EOS_NSPECIES; k++)
            cell_data[cell_count].Composition[k] = NewComposition[k];
#endif

          cell_data[cell_count].RunTime = timediff(network_start, second());
          network_timing_add_cell(&td, cell_data[cell_count].RunTime);
          cell_count++;
        }
    }

  CPU_Step[CPU_NETWORK_INTEGRATION] += measure_time();
  MPI_Barrier(MPI_COMM_WORLD);
  CPU_Step[CPU_NETWORK_IMBALANCE] += measure_time();

  /* all cells done, let's redistribute the results */
  mysort(cell_data, cell_count, sizeof(struct network_cell_data), cell_data_local_sort);

  int Send_Count[NTask];
  for(int j = 0; j < NTask; j++)
    Send_Count[j] = 0;

  for(int i = 0; i < cell_count; i++)
    Send_Count[cell_data[i].origTask]++;

  int Recv_Count[NTask];
  MPI_Alltoall(Send_Count, 1, MPI_INT, Recv_Count, 1, MPI_INT, MPI_COMM_WORLD);

  int Send_Offset[NTask], Recv_Offset[NTask];
  Recv_Offset[0]  = 0;
  Send_Offset[0]  = 0;
  int Local_Count = 0;
  for(int j = 0; j < NTask; j++)
    {
      Local_Count += Recv_Count[j];

      if(j > 0)
        {
          Send_Offset[j] = Send_Offset[j - 1] + Send_Count[j - 1];
          Recv_Offset[j] = Recv_Offset[j - 1] + Recv_Count[j - 1];
        }
    }

  struct network_cell_data *cell_results =
      (struct network_cell_data *)mymalloc("cell_results", Local_Count * sizeof(struct network_cell_data));

  memcpy(&cell_results[Recv_Offset[ThisTask]], &cell_data[Send_Offset[ThisTask]],
         Send_Count[ThisTask] * sizeof(struct network_cell_data));

  for(int ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_Count[recvTask] > 0 || Recv_Count[recvTask] > 0)
            {
              MPI_Sendrecv(&cell_data[Send_Offset[recvTask]], Send_Count[recvTask] * sizeof(struct network_cell_data), MPI_BYTE,
                           recvTask, TAG_DENS_B, &cell_results[Recv_Offset[recvTask]],
                           Recv_Count[recvTask] * sizeof(struct network_cell_data), MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  for(int idx = 0; idx < Local_Count; idx++)
    {
      int i = cell_results[idx].origIndex;

      for(int k = 0; k < EOS_NSPECIES; k++)
        {
          SphP[i].Composition[k]     = cell_results[idx].Composition[k];
          SphP[i].MassComposition[k] = SphP[i].Composition[k] * P[i].Mass;
        }

      SphP[i].Energy += cell_results[idx].dUtherm * P[i].Mass;
      SphP[i].dedt           = cell_results[idx].dUtherm / cell_results[idx].dt;
      SphP[i].NetworkRunTime = cell_results[idx].RunTime;
    }

  network_timing_summarise(&td);

  myfree(cell_results);
  myfree(cell_data);

  if(ThisTask == 0)
    myfree(work_data_all);

  if(!IsWorker)
    myfree(work_data_local);

  if(ThisTask_node == 0)
    myfree(data_node);

  myfree(data_local);

  mpi_printf("NUCLEAR NETWORK: done in %gs.\n", timediff(integrate_start, second()));

  CPU_Step[CPU_MISC] += measure_time();
}

void network_do_single_cell_mpi(struct network_cell_data *cell, double dt_cell, double *NewComposition)
{
  int threadid   = get_thread_num();
  int do_network = 1;

#ifdef NETWORK_NSE
  if(cell->EOSTemperature >= All.NetworkNSEThreshold)
    {
      for(int k = 0; k < EOS_NSPECIES; k++)
        NewComposition[k] = cell->Composition[k];

      double Temperature = cell->Temperature;
      network_nse_integrate_ye(cell->Utherm, cell->Density * All.UnitDensity_in_cgs, NewComposition, dt_cell * All.UnitTime_in_s,
                               &cell->dUtherm, &Temperature, &All.nd_nse, &All.nw_nse[threadid]);

      if(cell->Temperature > 0.9 * All.NetworkNSEThreshold)
        {
          // no need to do network
          for(int k = 0; k < EOS_NSPECIES; k++)
            cell->Composition[k] = NewComposition[k];
          do_network = 0;
        }
    }
#endif
  if(do_network)
    {
      for(int k = 0; k < EOS_NSPECIES; k++)
        NewComposition[k] = cell->Composition[k];

      network_integrate(cell->Temperature, cell->Density * All.UnitDensity_in_cgs, NewComposition, dt_cell * All.UnitTime_in_s,
                        &cell->dUtherm, &All.nd, &All.nw[threadid]);
    }

  cell->dUtherm *= All.UnitEnergy_in_cgs / All.UnitTime_in_s * dt_cell;
}
