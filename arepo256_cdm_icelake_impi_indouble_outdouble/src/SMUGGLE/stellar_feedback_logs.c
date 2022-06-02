/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 * 
 * \file        src/SMUGGLE/stellar_feedback_logs.c
 * \date        03/2020
 * \author      Federico Marinacci
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "../allvars.h"
#include "../proto.h"

#if defined(SMUGGLE_STAR_FEEDBACK) && defined(SMUGGLE_OUTPUT_STELLAR_FEEDBACK)

void output_stellar_feedback_statistics(void)
{
  mpi_printf("SMUGGLE_STAR_FEEDBACK: Writing feedback statistics.\n");

  /* star variables */
  MyFloat tot_feed_egy_time_step = 0.0, global_tot_feed_egy_time_step = 0.0;
  MyFloat SNII_Num_time_step = 0.0, global_SNII_Num_time_step = 0.0;
  MyFloat SNIa_Num_time_step = 0.0, global_SNIa_Num_time_step = 0.0;

  /* stellar feedback data (only active star particles are counted) */
  for(int i = 0; i < Nstar; i++)
    {
      tot_feed_egy_time_step += STP(StarParticle[i].index).FeedbackEnergy;
      SNII_Num_time_step += STP(StarParticle[i].index).SNII_Num;
      SNIa_Num_time_step += STP(StarParticle[i].index).SNIa_Num;
    }

  MPI_Reduce(&tot_feed_egy_time_step, &global_tot_feed_egy_time_step, 1, MPI_MYFLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&SNII_Num_time_step, &global_SNII_Num_time_step, 1, MPI_MYFLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&SNIa_Num_time_step, &global_SNIa_Num_time_step, 1, MPI_MYFLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fprintf(FdFeedback, "%e %e %e %e\n", All.Time, global_SNII_Num_time_step, global_SNIa_Num_time_step,
              global_tot_feed_egy_time_step);
      myflush(FdFeedback);
    }

  mpi_printf("SMUGGLE_STAR_FEEDBACK: Done.\n");
}

#endif
