/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/mpi_utils/healthtest.c
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

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

#define TEST_PACKET_SIZE_IN_MB 5
#define WORK_LOOP_COUNTER 30000000

#ifndef MAX_VARIATION_TOLERANCE
#define MAX_VARIATION_TOLERANCE 0.1
#endif

#ifndef DISABLE_HEALTHTEST
static const int bytecount = TEST_PACKET_SIZE_IN_MB * 1024L * 1024L;

static double measure_cpu_performance(MPI_Comm Communicator);
static double measure_hyper_cube_speed(const char *tag, MPI_Comm Communicator);
#endif

void healthtest(void)
{
  mpi_printf("\n");
#ifdef DISABLE_HEALTHTEST
  mpi_printf("HEALTHTEST: Disabled by config option DISABLE_HEALTHTEST\n");
#else
  measure_cpu_performance(MPI_COMM_WORLD);

#ifdef DISABLE_HEALTHTEST_FULL_HYPERCUBE
  mpi_printf("HEALTHTEST: %25s Disabled by config option DISABLE_HEALTHTEST_FULL_HYPERCUBE\n", "Full hypercube:");
#else
  /* Let's take a look at the communication speed in a global all-to-all data
   * exchange realized through pairwise exchanges along a hypercube */
  if(NTask > 1)
    measure_hyper_cube_speed("Full hypercube:", MPI_COMM_WORLD);
#endif

  /* Let's take a look at inter-node communication speed */
  if(NumNodes > 1)
    {
      const int CommSplitColor = RankInThisNode != 0;
      MPI_Comm comm;
      MPI_Comm_split(MPI_COMM_WORLD, CommSplitColor, ThisTask, &comm);

      if(RankInThisNode == 0)
        measure_hyper_cube_speed("Internode cube:", comm);

      MPI_Comm_free(&comm);
    }

  /* Now look at intra-node communication speed */
  if(NumNodes < NTask)
    {
      const int CommSplitColor = ThisNode;
      MPI_Comm comm;
      MPI_Comm_split(MPI_COMM_WORLD, CommSplitColor, ThisTask, &comm);

      measure_hyper_cube_speed("Intranode cube, 1st node:", comm);

      MPI_Comm_free(&comm);
    }
#endif
  mpi_printf("\n");
}

#ifndef DISABLE_HEALTHTEST
static double measure_cpu_performance(MPI_Comm Communicator)
{
  int loc_ntask, loc_thistask, loc_ptask;

  double ta = second();

  MPI_Comm_rank(Communicator, &loc_thistask);
  MPI_Comm_size(Communicator, &loc_ntask);

  for(loc_ptask = 0; loc_ntask > (1 << loc_ptask); loc_ptask++)
    ;

  double sum = 0;

  MPI_Barrier(Communicator);

  double t0 = second();

  /* do some computationally intense (but useless) work for a while */
  for(int i = 0; i < WORK_LOOP_COUNTER; i++)
    sum += sin((i + 0.1) / WORK_LOOP_COUNTER) / (2.0 + cos(i - 0.1) / WORK_LOOP_COUNTER);

  double t1 = second();

  double tperf = timediff(t0, t1), tperfsum;

  MPI_Allreduce(&tperf, &tperfsum, 1, MPI_DOUBLE, MPI_SUM, Communicator);
  double tavg = tperfsum / loc_ntask;

  struct
  {
    double t;
    int rank;
  } local = {tperf, ThisTask}, localnode = {tperf, ThisNode}, min_time, max_time, min_timenode, max_timenode;

  MPI_Allreduce(&local, &min_time, 1, MPI_DOUBLE_INT, MPI_MINLOC, Communicator);
  MPI_Allreduce(&local, &max_time, 1, MPI_DOUBLE_INT, MPI_MAXLOC, Communicator);

  MPI_Allreduce(&localnode, &min_timenode, 1, MPI_DOUBLE_INT, MPI_MINLOC, Communicator);
  MPI_Allreduce(&localnode, &max_timenode, 1, MPI_DOUBLE_INT, MPI_MAXLOC, Communicator);

  double variation = (max_time.t - min_time.t) / tavg;

  double tb = second();

  mpi_printf(
      "HEALTHTEST: %25s  %8.3f sec            %7.3f%%  variation   | Best=%g on Task=%d/Node=%d, Worst=%g on Task=%d/Node=%d, test "
      "took %g sec (sum=%g)\n",
      "CPU performance:", tavg, 100.0 * variation, min_time.t, min_time.rank, min_timenode.rank, max_time.t, max_time.rank,
      max_timenode.rank, timediff(ta, tb), sum);

  if(variation > MAX_VARIATION_TOLERANCE && ThisTask == 0)
    warn(
        "Consider stopping: the performance variation=%g of the CPUs lies above the prescribed tolerance "
        "MAX_VARIATION_TOLERANCE=%g, possibly indicating a machine problem.\n",
        variation, MAX_VARIATION_TOLERANCE);

  return sum;
}

static double measure_hyper_cube_speed(const char *const tag, MPI_Comm Communicator)
{
  double ta = second();

  int loc_ntask, loc_thistask, loc_ptask;

  MPI_Comm_rank(Communicator, &loc_thistask);
  MPI_Comm_size(Communicator, &loc_ntask);

  for(loc_ptask = 0; loc_ntask > (1 << loc_ptask); loc_ptask++)
    ;

  double tall = 0;
  int count   = 0;

  char *sendbuf = (char *)mymalloc_clear("send", bytecount * sizeof(char));
  char *recvbuf = (char *)mymalloc_clear("recv", bytecount * sizeof(char));

  /* exchange the test data */
  for(int ngrp = 1; ngrp < (1 << loc_ptask); ngrp++)
    {
      int recvTask = loc_thistask ^ ngrp;

      MPI_Barrier(Communicator);

      if(recvTask < loc_ntask)
        {
          double t0 = second();
          MPI_Sendrecv(sendbuf, bytecount, MPI_BYTE, recvTask, TAG_DENS_A, recvbuf, bytecount, MPI_BYTE, recvTask, TAG_DENS_A,
                       Communicator, MPI_STATUS_IGNORE);
          double t1 = second();

          tall += timediff(t0, t1);
          count++;
        }
    }

  myfree(recvbuf);
  myfree(sendbuf);

  double tperf = 0.5 * tall / count, tperfsum;

  MPI_Allreduce(&tperf, &tperfsum, 1, MPI_DOUBLE, MPI_SUM, Communicator);
  double tavg = tperfsum / loc_ntask;

  struct
  {
    double t;
    int rank;
  } local = {tperf, ThisTask}, localnode = {tperf, ThisNode}, min_time, max_time, min_timenode, max_timenode;

  MPI_Allreduce(&local, &min_time, 1, MPI_DOUBLE_INT, MPI_MINLOC, Communicator);
  MPI_Allreduce(&local, &max_time, 1, MPI_DOUBLE_INT, MPI_MAXLOC, Communicator);

  MPI_Allreduce(&localnode, &min_timenode, 1, MPI_DOUBLE_INT, MPI_MINLOC, Communicator);
  MPI_Allreduce(&localnode, &max_timenode, 1, MPI_DOUBLE_INT, MPI_MAXLOC, Communicator);

  double tb = second();

  double variation = (bytecount / min_time.t - bytecount / max_time.t) / (bytecount / tavg);

  mpi_printf(
      "HEALTHTEST: %25s  %8.1f MB/s per pair  %7.3f%%  variation   | Best=%g on Task=%d/Node=%d, Worst=%g on Task=%d/Node=%d, test "
      "took %g sec\n",
      tag, bytecount / tavg * TO_MBYTE_FAC, 100.0 * variation, bytecount / min_time.t * TO_MBYTE_FAC, min_time.rank, min_timenode.rank,
      bytecount / max_time.t * TO_MBYTE_FAC, max_time.rank, max_timenode.rank, timediff(ta, tb));

  if(variation > MAX_VARIATION_TOLERANCE && ThisTask == 0)
    warn(
        "\nThe performance variation=%g of the communication speed lies above the prescribed tolerance MAX_VARIATION_TOLERANCE=%g, "
        "possibly indicating a machine problem.\n",
        variation, MAX_VARIATION_TOLERANCE);

  return tavg;
}
#endif
