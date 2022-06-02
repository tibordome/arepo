/*! * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/MRT/RT_spectra.c
 * \date        11/2017
 * \author      Federico Marinacci
 * \brief       Collection of routines for reading and initializing spectra for the RT module
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"
#include "RT.h"

#ifdef MRT

static double blackbody_spectrum(double e)
{
  double T_eff = 1e5;
  double hc    = CLIGHT * PLANCK;
  double I_nu  = 2.0 * pow(e * ELECTRONVOLT_IN_ERGS, 3) / (hc * hc) /
                (exp(e * ELECTRONVOLT_IN_ERGS / (BOLTZMANN * T_eff)) - 1.0);  // Blackbody with T_eff

  return I_nu;
}

#if defined(MRT_SOURCES) || defined(MRT_LOCAL_FEEDBACK)

#ifndef HAVE_HDF5
#error "Need HAVE_HDF5 enabled to load stellar spectrum tables"
#endif

double PhotRate;
int Nfreq, N_age, N_metallicity;
double *PhotEnergy, *LumPhotons;
double *LogAge, *LogMetallicity;
double *lum_tab;  // remember to free these variables when the program finishes

void mrt_load_spectrum_table(void)
{
  hid_t file_id, dataset;
  char fname[MAXLEN_PATH];

  /* read in table dimensions for source spectrum */
  /* very important: Energy in eV, Lnu in erg/s/Hz */
  sprintf(fname, "%s", All.SpectrumTablePath);

  file_id = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  /* read in number of frequencies */
  dataset = my_H5Dopen(file_id, "Number_of_freq");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Nfreq, "Number_of_freq");
  my_H5Dclose(dataset, "Number_of_freq");

  /* read in number of ages to interpolate from */
  dataset = my_H5Dopen_if_existing(file_id, "N_age");
  if(dataset < 0)
    {
      // printf("%s\n", );
      N_age = 1;
    }
  else
    {
      my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &N_age, "N_age");
      my_H5Dclose(dataset, "N_age");

      /* read in the age array */
      LogAge  = (double *)mymalloc("LogAge", N_age * sizeof(double));
      dataset = my_H5Dopen(file_id, "age");
      my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, LogAge, "age");
      my_H5Dclose(dataset, "age");
      for(int i = 0; i < N_age; i++)
        {
          if(LogAge[i] < 1e-6)
            LogAge[i] = 1e-6;
          LogAge[i] = log10(LogAge[i]);
        }
    }

  /* read in number of metallicities to interpolate from */
  dataset = my_H5Dopen_if_existing(file_id, "N_metallicity");
  if(dataset < 0)
    {
      // printf("%s\n", );
      N_metallicity = 1;
    }
  else
    {
      my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &N_metallicity, "N_metallicity");
      my_H5Dclose(dataset, "N_metallicity");

      /* read in the metallicity array */
      LogMetallicity = (double *)mymalloc("LogMetallicity", N_metallicity * sizeof(double));
      dataset        = my_H5Dopen(file_id, "metallicity");
      my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, LogMetallicity, "metallicity");
      my_H5Dclose(dataset, "metallicity");
      for(int i = 0; i < N_metallicity; i++)
        LogMetallicity[i] = log10(LogMetallicity[i]);
    }

  /* allocate luminosity table */
  lum_tab = (double *)mymalloc("lum_tab", UV_BINS * N_age * N_metallicity * sizeof(double));
  for(int i = 0; i < UV_BINS * N_age * N_metallicity; i++)
    lum_tab[i] = 0;

  /* allocate and read energy array */
  PhotEnergy = (double *)mymalloc("PhotEnergy", Nfreq * sizeof(double));
  dataset    = my_H5Dopen(file_id, "Energy");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, PhotEnergy, "Energy");
  my_H5Dclose(dataset, "Energy");

  /* allocate and read photon array */
  LumPhotons = (double *)mymalloc("LumPhotons", Nfreq * N_age * N_metallicity * sizeof(double));
  dataset    = my_H5Dopen(file_id, "Lnu");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, LumPhotons, "Lnu");
  my_H5Dclose(dataset, "Lnu");

  my_H5Fclose(file_id, fname);

  /* from erg/s/Hz to phot/s */
  for(int i = 0; i < Nfreq * N_age * N_metallicity; i++)
    LumPhotons[i] /= PLANCK;

  mpi_printf_rt(0, "RT: Source spectrum loaded...\n");
}

void mrt_free_spectrum_table(void)
{
  myfree(LumPhotons);
  myfree(PhotEnergy);

  mpi_printf_rt(0, "RT: Source spectrum initialized...\n");
}

#if defined(MRT_STARS) || defined(MRT_LOCAL_FEEDBACK)
double star_spectrum(double e, int i_age, int i_metallicity)
{
  int idx    = 0;
  double res = 0.0;
  double de;
  // PhotRate = 0.0;

  // find index for linear spectrum interpolation
  while(PhotEnergy[idx] <= e)
    idx++;

  // above/below the upper/lower boundary -> no photons
  if((idx >= Nfreq) || (idx <= 0))
    return res;

  de         = (e - PhotEnergy[idx - 1]) / (PhotEnergy[idx] - PhotEnergy[idx - 1]);
  int offset = i_metallicity * N_age * Nfreq + i_age * Nfreq;
  res        = (1. - de) * LumPhotons[offset + idx - 1] + de * LumPhotons[offset + idx];

  return res;
}
#endif

#endif

double get_spectrum(double e, int i_age, int i_metallicity)
{
#if(defined(MRT_SOURCES) && defined(MRT_STARS))
  return star_spectrum(e, i_age, i_metallicity);
#else
  return blackbody_spectrum(e);
#endif
}

#endif
