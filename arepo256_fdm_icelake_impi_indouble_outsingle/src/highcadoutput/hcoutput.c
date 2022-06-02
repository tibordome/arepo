/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/highcadoutput/hcoutput.c
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

#include "../allvars.h"
#include "../proto.h"

#include "hcoutput.h"

#ifdef HCOUTPUT

void hcoutput_init(void)
{
  /* read list of render times */
  if(RestartFlag != RESTART_AURIGA_MOVIE)
    {
      if(ThisTask == 0)
        {
          FILE *fd;
          int count;
          char buf[512], msg[512];

          if(!(fd = fopen(All.HCOutput_OutputListFilename, "r")))
            {
              printf("HCOUTPUT: can't read output list in file '%s'\n", All.HCOutput_OutputListFilename);
              return;
            }

          All.HCOutput_OutputListLength = 0;

          while(1)
            {
              if(fgets(buf, 500, fd) != buf)
                break;

              count = sscanf(buf, " %lg ", &All.HCOutput_OutputListTimes[All.HCOutput_OutputListLength]);

              if(count == 1)
                {
                  if(All.HCOutput_OutputListLength >= MAX_HCSNIPS_OUTPUT_TIMES)
                    {
                      sprintf(msg, "\nHCOUTPUT: too many entries in output-list. You should increase MAX_HCSNIPS_OUTPUT_TIMEST=%d.\n",
                              (int)MAX_HCSNIPS_OUTPUT_TIMES);
                      terminate(msg);
                    }

                  All.HCOutput_OutputListLength++;
                }
            }

          fclose(fd);

          printf("\nHCOUTPUT: found %d times in output-list.\n", All.HCOutput_OutputListLength);
        }
      MPI_Bcast(&All.HCOutput_OutputListLength, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(All.HCOutput_OutputListTimes, All.HCOutput_OutputListLength, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

  All.HCOutput_RadialCut = All.HCOutput_RadialCut * All.HubbleParam / All.Time;
}

void hcoutput_check_output(void)
{
  if(All.HighestActiveTimeBin != All.HighestOccupiedTimeBin)
    return;

  mpi_printf("All.Time=%lg,All.HCOutput_NextOutputTime=%lg\n", All.Time, All.HCOutput_NextOutputTime);
  if(All.Time >= All.HCOutput_NextOutputTime)
    {
      mpi_printf("HCOUTPUT: Writing frame %d at t=%g, target time was %g.\n", All.HCSnipShotFileCount, All.Time,
                 All.HCOutput_NextOutputTime);
      hcoutput_dump();

      if(All.HCSnipShotFileCount < All.HCOutput_OutputListLength)
        {
          All.HCOutput_NextOutputTime = All.HCOutput_OutputListTimes[All.HCSnipShotFileCount];

          if(All.HCSnipShotFileCount + 1 < All.HCOutput_OutputListLength)
            All.HCOutput_NextOutputDelta =
                0.1 * (All.HCOutput_OutputListTimes[All.HCSnipShotFileCount + 1] - All.HCOutput_NextOutputTime);
          else
            All.HCOutput_NextOutputDelta = All.TimeMax;
        }
      else
        {
          All.HCOutput_NextOutputTime  = All.TimeMax;
          All.HCOutput_NextOutputDelta = All.TimeMax;
        }
    }
}

/* computes projections at current time (on the fly) or for given snapshot */
void hcoutput_dump(void)
{
  double center[3], cmvel[3];
  int outNum;

  hcoutput_get_center_guess_onthefly(center);

  outNum = All.HCSnipShotFileCount;
  All.HCSnipShotFileCount++;

  hcoutput_compute_center(center, cmvel, 5);
  DumpFlag = DUMP_HCOUTPUT;
  savepositions(outNum, 1);
}

void hcoutput_get_center_guess_onthefly(double center[3])
{
  All.HCOutput_Halo_Initialized = 1;

  double cm[3];
  double mass;

  int k;
  for(k = 0; k < 3; k++)
    cm[k] = 0;
  mass = 0;

  int i;
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].HCOutput_Center_Marker)
        {
          double dx, dy, dz;

          if(P[i].Type == 0)
            {
              dx = NEAREST_X(SphP[i].Center[0] - 0.5 * All.BoxSize);
              dy = NEAREST_Y(SphP[i].Center[1] - 0.5 * All.BoxSize);
              dz = NEAREST_Z(SphP[i].Center[2] - 0.5 * All.BoxSize);
            }
          else
            {
              dx = NEAREST_X(P[i].Pos[0] - 0.5 * All.BoxSize);
              dy = NEAREST_Y(P[i].Pos[1] - 0.5 * All.BoxSize);
              dz = NEAREST_Z(P[i].Pos[2] - 0.5 * All.BoxSize);
            }

          cm[0] += dx * P[i].Mass;
          cm[1] += dy * P[i].Mass;
          cm[2] += dz * P[i].Mass;
          mass += P[i].Mass;
        }
    }

  double allcm[3], allmass;
  MPI_Allreduce(&mass, &allmass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(cm, allcm, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(allmass == 0)
    terminate("Fail.");

  for(k = 0; k < 3; k++)
    {
      center[k] = 0.5 * All.BoxSize + allcm[k] / allmass;
    }
}

void hcoutput_compute_center(double center[3], double cmvel[3], int niter)
{
  double radius  = All.HCOutput_CenterRadius / All.HubbleParam * All.Time; /* convert from physical to internal code units */
  double radius2 = radius * radius;

  /* recenter, radius 20kpc physical */
  int iter, k;
  for(iter = 0; iter < niter; iter++)
    {
      double cm[3], mass;
      int j;

      mass = 0;
      for(j = 0; j < 3; j++)
        cm[j] = 0;

      int i;
      for(i = 0; i < NumPart; i++)
        {
          double dx, dy, dz;

          if(P[i].Type == 0)
            {
              dx = NEAREST_X(SphP[i].Center[0] - center[0]);
              dy = NEAREST_Y(SphP[i].Center[1] - center[1]);
              dz = NEAREST_Z(SphP[i].Center[2] - center[2]);
            }
          else
            {
              dx = NEAREST_X(P[i].Pos[0] - center[0]);
              dy = NEAREST_Y(P[i].Pos[1] - center[1]);
              dz = NEAREST_Z(P[i].Pos[2] - center[2]);
            }
          double r2 = dx * dx + dy * dy + dz * dz;

          if(r2 < radius2)
            {
              mass += P[i].Mass;
              cm[0] += dx * P[i].Mass;
              cm[1] += dy * P[i].Mass;
              cm[2] += dz * P[i].Mass;
            }
        }

      double allcm[3], allmass;
      MPI_Allreduce(&mass, &allmass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      if(allmass == 0)
        {
          /* radius is too small, stop with current center */
          mpi_printf("HCOUTPUT: iter=%d, no mass in sphere of radius %g, using current center estimate of %g,%g,%g.\n", iter,
                     sqrt(radius2), center[0], center[1], center[2]);
          MPI_Bcast(center, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
          break;
        }

      MPI_Reduce(cm, allcm, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

      if(ThisTask == 0)
        {
          for(j = 0; j < 3; j++)
            center[j] += allcm[j] / allmass;

          mpi_printf("HCOUTPUT: iter=%d, Mass in CenterRadius: %g, delta=%g,%g,%g, new center: %g,%g,%g\n", iter, allmass,
                     allcm[0] / allmass, allcm[1] / allmass, allcm[2] / allmass, center[0], center[1], center[2]);
        }

      MPI_Bcast(center, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      for(k = 0; k < 3; k++)
        All.HCOutput_Halo_Center[k] = center[k];

      mpi_printf("All.HCOutput_RadialCut=%lg\n", All.HCOutput_RadialCut);
      mpi_printf("All.HCOutput_Halo_Center=%lg|%lg|%lg\n", All.HCOutput_Halo_Center[0], All.HCOutput_Halo_Center[1],
                 All.HCOutput_Halo_Center[2]);
      radius2 *= 0.7 * 0.7;
    }

  radius2 = radius * radius;
}

#endif
