/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/timestep.h
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

#ifndef TIMESTEP_H
#define TIMESTEP_H

#include "allvars.h"

#ifdef ENLARGE_DYNAMIC_RANGE_IN_TIME
typedef long long integertime;
#define MPI_INTEGERTIME MPI_LONG_LONG
#define INTEGERTIME_PRI "lld"
#define TIMEBINS 60
#else
typedef int integertime;
#define MPI_INTEGERTIME MPI_INT
#define INTEGERTIME_PRI "d"
#define TIMEBINS 29
#endif
/* The simulated timespan is mapped onto the integer interval [0, TIMESPAN],
 * where TIMESPAN needs to be a power of 2. */
#define TIMEBASE (((integertime)1) << TIMEBINS)

enum timstep_types
{
  TSTP_UNKNOWN = 0,
  TSTP_GRAVITY,
  TSTP_GRAVITY_EXTERNAL,
  TSTP_SIDM,
  TSTP_OTVET,
  TSTP_BECDM,
  TSTP_MAXTIMESTEP,
  TSTP_MOVIE,
  TSTP_HCOUTPUT,
  TSTP_DISPLACEMENT,
  TSTP_BLACKHOLE,
  TSTP_SINKS,
  TSTP_COURANT,
  TSTP_COURANT_TREE,
  TSTP_DEDNER,
  TSTP_SFR,
  TSTP_MRT,
  TSTP_NONIDEALMHD,
  TSTP_POWELL,
  TSTP_SINK_INTERACTION,
  TSTP_OUTPUTLIMIT,
  TSTP_CR_DIFFUSION,
  TSTP_BRAGINSKII_VISCOSITY,
  TSTP_METAL_DIFFUSION,
  TSTP_CONDUCTION,
  TSTP_OHMIC_DIFFUSION,
  TSTP_DUST_LIVE,
  TSTP_STELLAR_EVO
};

struct TimeBinData
{
  int NActiveParticles;
  long long GlobalNActiveParticles;
  int *ActiveParticleList;
  int TimeBinCount[TIMEBINS];

  int FirstInTimeBin[TIMEBINS];
  int LastInTimeBin[TIMEBINS];
  int *NextInTimeBin;
  int *PrevInTimeBin;
  char Name[100];
  int *MaxPart;
};

void find_timesteps_without_gravity(void);
void update_timesteps_from_gravity(void);
integertime get_timestep_gravity(int p);
integertime get_timestep_hydro(int p);
integertime get_timestep_bhaccretion(int p);
integertime get_timestep_stellarevolution(int p);
integertime get_timestep_sinksaccretion(int p);
integertime get_timestep_pm(void);
int test_if_grav_timestep_is_too_large(int p, int bin);
void validate_timestep(double dt, integertime ti_step, int p, int type_tstp);
void find_dt_displacement_constraint(double hfac);
int get_timestep_bin(integertime ti_step);
double get_time_difference_in_Gyr(double a0, double a1);

/* TimeBinData stuff */
void timebins_init(struct TimeBinData *tbData, const char *name, int *MaxPart);
void timebins_allocate(struct TimeBinData *tbData);
void timebins_reallocate(struct TimeBinData *tbData);
int timebins_get_bin_and_do_validity_checks(integertime ti_step, int bin_old);
void timebin_move_particle(struct TimeBinData *tbData, int p, int timeBin_old, int timeBin_new);
void timebin_add_particle(struct TimeBinData *tbData, int i_new, int i_old, int timeBin, int addToListOfActiveParticles);
void timebin_remove_particle(struct TimeBinData *tbData, int idx, int bin);
void timebin_cleanup_list_of_active_particles(struct TimeBinData *tbData);
void timebin_move_sfr(int p, int timeBin_old, int timeBin_new);
void timebin_move_bh(int p, int timeBin_old, int timeBin_new);
void timebin_make_list_of_active_particles_up_to_timebin(struct TimeBinData *tbData, int timebin);
void timebin_add_particles_of_timebin_to_list_of_active_particles(struct TimeBinData *tbData, int timebin);

#endif /* TIMESTEP_H */
