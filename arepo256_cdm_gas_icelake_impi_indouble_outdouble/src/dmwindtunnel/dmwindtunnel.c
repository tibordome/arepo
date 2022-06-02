/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dmwindtunnel/dmwindtunnel.c
 * \date        12/2020
 * \author      Philip Mocz, Martin Sparre, Hayden Foote
 * \brief       includes functions for DM windtunnel
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - 08.02.2021 HF added functions for reading and interpolating the wind table from a file
 * - 01.11.2022 HF added star particles to the windtunnel
 *
 */

#include <gmp.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

// If we need to read the wind params from a file, here are the required functions
#if(defined(DM_WINDTUNNEL) && defined(DM_WINDTUNNEL_EXTERNAL_SOURCE))

static int NWindTable;

/*! \brief Structure stores wind properties (at a given time).
 */
static struct wind_table
{
  double time, rho, vel, sigma;
} * WindTable;

/*! \brief Interpolates density, velocity, and velocity dispersion.
 *
 *  Interpolation in time between entries of wind table array.
 *
 *  \param[in] t Time.
 *  \param[out] rho Density (as a result of the interpolation).
 *  \param[out] vel Velocity (as a result of the interpolation).
 *  \param[out] sigma Velocity dispersion (as a result of the interpolation).
 *
 *  \return void
 */
void interpolate_from_dm_wind_table(double t, double *rho, double *vel, double *sigma)
{
  int binlow, binhigh, binmid;
  // Check to make sure the current sim time is within the wind table
  if(t < WindTable[0].time || t > WindTable[NWindTable - 1].time)
    {
      terminate("time outside of wind interpolation table");
    }
  // Indices of lowest and highest time bin
  binlow  = 0;
  binhigh = NWindTable - 1;
  // Find the bin where the sim time lies with bisection
  while(binhigh - binlow > 1)
    {
      binmid = (binhigh + binlow) / 2;  // find the middle index
      if(t < WindTable[binmid].time)    // if the time is lower than the midpoint
        binhigh = binmid;               // make the midpoint the new high bin
      else
        binlow = binmid;  // make the midpoint the new low bin
    }
  // length of the time bin in sim units
  double dt = WindTable[binhigh].time - WindTable[binlow].time;
  // Throw an error and stop the sim if dt is zero
  if(dt == 0)
    terminate("dt=0");
  // Fraction of bin (from bottom) where current sim time is
  double u = (t - WindTable[binlow].time) / dt;
  // Linear interpolation of the wind table to set the wind parameters for this timestep
  double rho_step   = (1 - u) * WindTable[binlow].rho + u * WindTable[binhigh].rho;
  double vel_step   = (1 - u) * WindTable[binlow].vel + u * WindTable[binhigh].vel;
  double sigma_step = (1 - u) * WindTable[binlow].sigma + u * WindTable[binhigh].sigma;
  // Set the windtunnel parameters
  *rho   = rho_step;
  *vel   = vel_step;
  *sigma = sigma_step;
}

/*! \brief Read in dark matter windtunnel file.
 *
 *  Store the read-in input in WindTable.
 *
 *  \return void
 */
void read_dmwindtunnel_file(void)
{
  FILE *fd;  // file name
  double k, p;
  // Make sure we can open the wind file. If we can't, terminate.
  if(!(fd = fopen(All.DMWindtunnelExternalSourceFile, "r")))
    {
      char buf[1000];
      sprintf(buf, "can't read file '%s' with dark matter windtunnel data on task %d\n", buf, ThisTask);
      terminate(buf);
    }

  NWindTable = 0;  // Row counter
  do
    {
      double t, rho, v, sigma;  // parameters that should be in each row
      if(fscanf(fd, " %lg %lg %lg %lg ", &t, &rho, &v, &sigma) == 4)
        NWindTable++;  // Count rows
      else
        break;
    }
  while(1);

  fclose(fd);
  // Let the user know how many rows were found
  mpi_printf("found %d rows in input dark matter wind table\n", NWindTable);
  // Allocate memory for the wind table based on how many rows were found
  WindTable = (struct wind_table *)mymalloc("WindTable", NWindTable * sizeof(struct wind_table));
  // Make sure we can still read the wind file
  if(!(fd = fopen(All.DMWindtunnelExternalSourceFile, "r")))
    {
      char buf[1000];
      sprintf(buf, "can't read file '%s' with dark matter windtunnel data on task %d\n", buf, ThisTask);
      terminate(buf);
    }
  // reset the row counter
  NWindTable = 0;
  do  // read the wind file into the wind table
    {
      double t, rho, v, sigma;
      if(fscanf(fd, " %lg %lg %lg %lg ", &t, &rho, &v, &sigma) == 4)
        {
          WindTable[NWindTable].time  = t;
          WindTable[NWindTable].rho   = rho;
          WindTable[NWindTable].vel   = v;
          WindTable[NWindTable].sigma = sigma;
          NWindTable++;
        }
      else
        break;
    }
  while(1);

  fclose(fd);

  /* note: we'll assume that this file is sorted according to time (in simulation units)*/
}
#endif  // (defined(DM_WINDTUNNEL) && defined(DM_WINDTUNNEL_EXTERNAL_SOURCE))

// If DM_WINDTUNNEL is defined, function for updating the BCs
#ifdef DM_WINDTUNNEL

void apply_windtunnel_bcs(void)
{
  long long NInjectionRegionAll;
  long long NInjectionRegion = 0;

  for(int p = 0; p < NumPart; p++)
    if(P[p].Type == PTYPE_HALO || P[p].Type == PTYPE_STARS)
      if(WRAP_Y(P[p].Pos[1] - All.GlobalDisplacementVector[1]) < All.DMWindtunnelInjectionRegion &&
         WRAP_Y(P[p].Pos[1] - All.GlobalDisplacementVector[1]) > 0.0)
        NInjectionRegion++;

  MPI_Allreduce(&NInjectionRegion, &NInjectionRegionAll, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
  // Print out windtunnel parameters at this step
  mpi_printf("DM_WINDTUNNEL: Density = %e, Velocity = %e, Sigma = %e \n", All.DMWindtunnelInjectionDensity, All.DMWindtunnelVY,
             All.DMWindtunnelSigmaVY);
#ifdef DM_WINDTUNNEL_STARS
  mpi_printf("DM_WINDTUNNEL_STARS: Density = %e, Velocity = %e, Sigma = %e \n", All.StarWindtunnelInjectionDensity, All.DMWindtunnelVY,
             All.StarWindtunnelSigmaVY);
#endif  // DM_WINDTUNNEL_STARS
  for(int p = 0; p < NumPart; p++)
    {
      if(P[p].Type == PTYPE_HALO || P[p].Type == PTYPE_STARS)
        {
          if(WRAP_Y(P[p].Pos[1] - All.GlobalDisplacementVector[1]) > All.DMWindtunnelInjectionRegion ||
             WRAP_Y(P[p].Pos[1] - All.GlobalDisplacementVector[1]) < 0.0)
            P[p].DMWindtunnel_RecentlyUpdated = 0;
          // Modify all particles in the injection region if they have not yet been updated
          // else if( WRAP_Y(P[p].Pos[1] - All.GlobalDisplacementVector[1]) < All.DMWindtunnelInjectionRegion && WRAP_Y(P[p].Pos[1] -
          // All.GlobalDisplacementVector[1]) > 0.0 && P[p].DMWindtunnel_RecentlyUpdated == 0 ){ Modify all particles in the injection
          // region
          else if(WRAP_Y(P[p].Pos[1] - All.GlobalDisplacementVector[1]) < All.DMWindtunnelInjectionRegion &&
                  WRAP_Y(P[p].Pos[1] - All.GlobalDisplacementVector[1]) > 0.0)
            {
              /* All velocity dispersions are set to SigmaVY for now */
              /* Normal velocity distribution sampling */
              if(P[p].Type == 1)
                {  // FOR DARK MATTER
                  P[p].Vel[0] = gsl_ran_gaussian(random_generator, All.DMWindtunnelSigmaVY);
                  P[p].Vel[1] = All.DMWindtunnelVY + gsl_ran_gaussian(random_generator, All.DMWindtunnelSigmaVY);
                  P[p].Vel[2] = gsl_ran_gaussian(random_generator, All.DMWindtunnelSigmaVY);

                  /* Tophat velocity distribution sampling */
                  /* P[p].Vel[0] = All.DMWindtunnelSigmaVY*get_random_number() - 0.5*All.DMWindtunnelSigmaVY;
                  P[p].Vel[1] = All.DMWindtunnelVY + All.DMWindtunnelSigmaVY*get_random_number() - 0.5*All.DMWindtunnelSigmaVY;
                  P[p].Vel[2] = All.DMWindtunnelSigmaVY*get_random_number() - 0.5*All.DMWindtunnelSigmaVY; */

                  P[p].Pos[0] = All.BoxSize * get_random_number();  // random number between 0 and BoxSize
                  P[p].Pos[2] = All.BoxSize * get_random_number();  // Random number between 0 and BoxSize

                  P[p].Mass = All.DMWindtunnelInjectionDensity * All.DMWindtunnelInjectionRegion * All.BoxSize * All.BoxSize /
                              NInjectionRegionAll;

                  P[p].DMWindtunnel_RecentlyUpdated = 1;
                }
#ifdef DM_WINDTUNNEL_STARS
              if(P[p].Type == PTYPE_STARS)
                {  // FOR STARS
                  P[p].Vel[0] = gsl_ran_gaussian(random_generator, All.StarWindtunnelSigmaVY);
                  P[p].Vel[1] = All.DMWindtunnelVY + gsl_ran_gaussian(random_generator, All.StarWindtunnelSigmaVY);
                  P[p].Vel[2] = gsl_ran_gaussian(random_generator, All.StarWindtunnelSigmaVY);

                  /* Tophat velocity distribution sampling */
                  /* P[p].Vel[0] = All.DMWindtunnelSigmaVY*get_random_number() - 0.5*All.DMWindtunnelSigmaVY;
                  P[p].Vel[1] = All.DMWindtunnelVY + All.DMWindtunnelSigmaVY*get_random_number() - 0.5*All.DMWindtunnelSigmaVY;
                  P[p].Vel[2] = All.DMWindtunnelSigmaVY*get_random_number() - 0.5*All.DMWindtunnelSigmaVY; */

                  P[p].Pos[0] = All.BoxSize * get_random_number();  // random number between 0 and BoxSize
                  P[p].Pos[2] = All.BoxSize * get_random_number();  // Random number between 0 and BoxSize

                  P[p].Mass = All.StarWindtunnelInjectionDensity * All.DMWindtunnelInjectionRegion * All.BoxSize * All.BoxSize /
                              NInjectionRegionAll;

                  P[p].DMWindtunnel_RecentlyUpdated = 1;
                }
#endif  // DM_WINDTUNNEL_STARS
            }
        }
    }
}

#endif
