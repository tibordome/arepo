/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sinks/sinks_aux.c
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       Sink particles
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"

#ifdef SINKS

void sinks_set_constants(void)
{
  set_cosmo_factors_for_current_time();

  SKD.DistFac = pow(All.cf_atime / All.HubbleParam, 2);

  SKD.NHFac = HYDROGEN_MASSFRAC * All.UnitDensity_in_cgs * All.cf_a3inv * pow(All.HubbleParam, 2) / PROTONMASS;

#ifdef SGCHEM
  SKD.DistFac = 1.0;
  SKD.NHFac   = HYDROGEN_MASSFRAC * All.UnitDensity_in_cgs / PROTONMASS;

  mpi_printf("Set sink constants %d %g %g\n", ThisTask, SKD.NHFac, All.UnitDensity_in_cgs / PROTONMASS);

#endif
}

void write_sink_data(int xaxis, int yaxis, int zaxis, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{
  int i;
  char buf[MAXLEN_PATH];
  FILE *fd = 0;
  file_path_sprintf(buf, "%s/sinks_%03d", All.OutputDir, RestartSnapNum);

  fd = open_file(buf);

  int num_sinks = 0;
  double sink_pos[2];

  for(i = 0; i < SKD.TotNumSinks; i++)
    {
      if(SKD.SINK[i].Pos[xaxis] >= xmin && SKD.SINK[i].Pos[xaxis] < xmax && SKD.SINK[i].Pos[yaxis] >= ymin &&
         SKD.SINK[i].Pos[yaxis] < ymax && SKD.SINK[i].Pos[zaxis] >= zmin && SKD.SINK[i].Pos[zaxis] < zmax)
        num_sinks++;
    }

  my_fwrite(&num_sinks, sizeof(int), 1, fd);

  for(i = 0; i < SKD.TotNumSinks; i++)
    {
      if(SKD.SINK[i].Pos[xaxis] >= xmin && SKD.SINK[i].Pos[xaxis] < xmax && SKD.SINK[i].Pos[yaxis] >= ymin &&
         SKD.SINK[i].Pos[yaxis] < ymax && SKD.SINK[i].Pos[zaxis] >= zmin && SKD.SINK[i].Pos[zaxis] < zmax)
        {
          sink_pos[0] = (SKD.SINK[i].Pos[xaxis] - xmin) / (xmax - xmin);
          sink_pos[1] = (SKD.SINK[i].Pos[yaxis] - ymin) / (ymax - ymin);

          my_fwrite(&sink_pos, sizeof(double), 2, fd);
          my_fwrite(&SKD.SINK[i].Mass, sizeof(double), 1, fd);
        }
    }

  fclose(fd);
}

double get_accretion_radius(double mass, double temp)
{
  double n_bondi, acc_radius;
  double min_acc_radius = 1e-4;

  n_bondi = 4.66e13 * (temp / 6000) * (temp / 6000) * (temp / 6000) / (mass * All.UnitMass_in_g / SOLAR_MASS) /
            (mass * All.UnitMass_in_g / SOLAR_MASS);

  if(SKD.NHThresh >= n_bondi)
    acc_radius = 17 * (mass * All.UnitMass_in_g / SOLAR_MASS) * (6000 / temp) * (ASTRONOMICAL_UNIT / All.UnitLength_in_cm);
  else
    acc_radius = sqrt(1.49e17 / SKD.NHThresh) * sqrt(temp / 6000) * (ASTRONOMICAL_UNIT / All.UnitLength_in_cm);

  /*if (acc_radius < min_acc_radius)
    return min_acc_radius;
  else
    return acc_radius;*/
  return acc_radius;
}

#endif
