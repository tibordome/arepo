/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/movie_auriga/movie.c
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

#include "movie.h"

#ifdef AURIGA_MOVIE

char *aum_projection_names[aum_count] = {
    "aum_gas_rho",   "aum_gas_vel", "aum_gas_velx", "aum_gas_vely", "aum_gas_velz",         "aum_gas_bfld", "aum_gas_bx",
    "aum_gas_by",    "aum_gas_bz",  "aum_gas_temp", "aum_gas_pres", "aum_gas_metallicity",  "aum_dm_rho",   "aum_dm_rho2",
    "aum_stars_rho", "aum_stars_u", "aum_stars_g",  "aum_stars_r",  "aum_stars_metallicity"};

static void auriga_movie_compute_thin_global_projections(int outNum);
static void auriga_movie_compute_halo_projections(int outNum, double center[3], double cmvel[3]);
static void auriga_movie_compute_galaxy_projections(int outNum, double center[3], double cmvel[3]);
static void auriga_movie_compute_sky_projections(int outNum, double center[3], double rotation[3][3], double cmvel[3]);

void auriga_movie_init(void)
{
  /* read list of render times */
  if(RestartFlag != RESTART_AURIGA_MOVIE)
    {
      if(ThisTask == 0)
        {
          FILE *fd;
          int count;
          char buf[512], msg[512];

          if(!(fd = fopen(All.Auriga_Movie_OutputListFilename, "r")))
            {
              printf("AURIGA MOVIE: can't read output list in file '%s'\n", All.Auriga_Movie_OutputListFilename);
              return;
            }

          All.Auriga_Movie_OutputListLength = 0;

          while(1)
            {
              if(fgets(buf, 500, fd) != buf)
                break;

              count = sscanf(buf, " %lg ", &All.Auriga_Movie_OutputListTimes[All.Auriga_Movie_OutputListLength]);

              if(count == 1)
                {
                  if(All.Auriga_Movie_OutputListLength >= MAX_MOVIE_OUTPUT_TIMES)
                    {
                      sprintf(msg,
                              "\nAURIGA MOVIE: too many entries in output-list. You should increase MAX_MOVIE_OUTPUT_TIMEST=%d.\n",
                              (int)MAX_MOVIE_OUTPUT_TIMES);
                      terminate(msg);
                    }

                  All.Auriga_Movie_OutputListLength++;
                }
            }

          fclose(fd);

          printf("\nAURIGA MOVIE: found %d times in output-list.\n", All.Auriga_Movie_OutputListLength);
        }
      MPI_Bcast(&All.Auriga_Movie_OutputListLength, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(All.Auriga_Movie_OutputListTimes, All.Auriga_Movie_OutputListLength, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
}

void auriga_movie_check_output(void)
{
  if(All.HighestActiveTimeBin != All.HighestOccupiedTimeBin)
    return;

  if(All.Time >= All.Auriga_Movie_NextOutputTime)
    {
      mpi_printf("AURIGA MOVIE: Writing frame %d at t=%g, target time was %g.\n", All.Auriga_Movie_Num, All.Time,
                 All.Auriga_Movie_NextOutputTime);
      auriga_movie_make();

      if(All.Auriga_Movie_Num < All.Auriga_Movie_OutputListLength)
        {
          All.Auriga_Movie_NextOutputTime = All.Auriga_Movie_OutputListTimes[All.Auriga_Movie_Num];

          if(All.Auriga_Movie_Num + 1 < All.Auriga_Movie_OutputListLength)
            All.Auriga_Movie_NextOutputDelta =
                0.1 * (All.Auriga_Movie_OutputListTimes[All.Auriga_Movie_Num + 1] - All.Auriga_Movie_NextOutputTime);
          else
            All.Auriga_Movie_NextOutputDelta = All.TimeMax;
        }
      else
        {
          All.Auriga_Movie_NextOutputTime  = All.TimeMax;
          All.Auriga_Movie_NextOutputDelta = All.TimeMax;
        }
    }
}

/* computes projections at current time (on the fly) or for given snapshot */
void auriga_movie_make(void)
{
  TIMER_START(CPU_AURIGA_MOVIE);

  double center[3], rotation[3][3], cmvel[3];
  int outNum, halo;

  mpi_printf("AURIGA MOVIE Rotation matrix: %8g %8g %8g\n", All.Auriga_Movie_Galaxy_Rotation[0][0],
             All.Auriga_Movie_Galaxy_Rotation[1][0], All.Auriga_Movie_Galaxy_Rotation[2][0]);
  mpi_printf("AURIGA MOVIE Rotation matrix: %8g %8g %8g\n", All.Auriga_Movie_Galaxy_Rotation[0][1],
             All.Auriga_Movie_Galaxy_Rotation[1][1], All.Auriga_Movie_Galaxy_Rotation[2][1]);
  mpi_printf("AURIGA MOVIE Rotation matrix: %8g %8g %8g\n", All.Auriga_Movie_Galaxy_Rotation[0][2],
             All.Auriga_Movie_Galaxy_Rotation[1][2], All.Auriga_Movie_Galaxy_Rotation[2][2]);

  if(RestartFlag == RESTART_AURIGA_MOVIE)
    {
      halo   = auriga_movie_get_center_guess_groupcatalogue(center);
      outNum = RestartSnapNum;
    }
  else
    {
      /* on the fly mode */
      double z = 1. / All.Time - 1.;
      if(z <= 10.)
        {
          auriga_movie_get_center_guess_onthefly(center);
          halo = 1;
        }
      else
        {
          halo = 0;
        }

      outNum = All.Auriga_Movie_Num;
      All.Auriga_Movie_Num++;
    }

  auriga_movie_calculate_smoothing_lenghts();

  auriga_movie_compute_thin_global_projections(outNum);

  if(halo)
    {
      int res = auriga_movie_compute_center_and_rotation(center, rotation, cmvel, 5);
      if(res > 0)
        {
          /* boxsize? 1Mpc? */
          auriga_movie_compute_halo_projections(outNum, center, cmvel);

          /* boxsize? 100kpc? */
          auriga_movie_compute_galaxy_projections(outNum, center, cmvel);

          /* at 8kpc on intermediate axis on side with smaller x coordinate? */
          if(res == 1)
            auriga_movie_compute_sky_projections(outNum, center, rotation, cmvel);
        }
    }

  TIMER_STOP(CPU_AURIGA_MOVIE);
}

void auriga_movie_compute_thin_global_projections(int outNum)
{
  double center[3], cmvel[3];
  auriga_movie_get_center_highres(center, cmvel);

  struct auriga_movie_projection proj;
  auriga_movie_projection_init(&proj, center, All.Auriga_Movie_Galaxy_Rotation, cmvel, 2.0, 2.0, 3840, 2160);

  auriga_movie_projection_project_particles(&proj);

  auriga_movie_projection_save(&proj, outNum, All.Auriga_Movie_Directory, "highres_region");

  auriga_movie_projection_free(&proj);
}

void auriga_movie_compute_halo_projections(int outNum, double center[3], double cmvel[3])
{
  struct auriga_movie_projection proj;
  auriga_movie_projection_init(&proj, center, All.Auriga_Movie_Galaxy_Rotation, cmvel, 0.4 / All.Time * All.HubbleParam,
                               0.4 / All.Time * All.HubbleParam, 3840, 2160);

  auriga_movie_projection_project_particles(&proj);

  auriga_movie_projection_save(&proj, outNum, All.Auriga_Movie_Directory, "halo");

  auriga_movie_projection_free(&proj);
}

void auriga_movie_compute_galaxy_projections(int outNum, double center[3], double cmvel[3])
{
  struct auriga_movie_projection proj;
  auriga_movie_projection_init(&proj, center, All.Auriga_Movie_Galaxy_Rotation, cmvel, 0.1 / All.Time * All.HubbleParam,
                               0.1 / All.Time * All.HubbleParam, 3840, 2160);

  auriga_movie_projection_project_particles(&proj);

  auriga_movie_projection_save(&proj, outNum, All.Auriga_Movie_Directory, "galaxy");

  auriga_movie_projection_free(&proj);
}

void auriga_movie_compute_sky_projections(int outNum, double center[3], double rotation[3][3], double cmvel[3])
{
  struct auriga_movie_projection proj;
  auriga_movie_projection_init(&proj, center, rotation, cmvel, 0, 0.2 / All.Time * All.HubbleParam, 3840, 0);
  // auriga_movie_projection_init( &proj, center, rotation, cmvel, 0, 0.2 / All.Time * All.HubbleParam, 64, 0 );

  auriga_movie_projection_project_particles(&proj);

  auriga_movie_projection_save(&proj, outNum, All.Auriga_Movie_Directory, "allsky");

  auriga_movie_projection_free(&proj);
}
#endif
