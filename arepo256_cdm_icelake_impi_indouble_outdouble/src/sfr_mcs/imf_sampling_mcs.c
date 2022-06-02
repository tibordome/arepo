/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sfr_mcs/imf_sampling_mcs.c
 * \date        08/2019
 * \author     	Matthew C Smith
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - 17.01.2022 Ported into current codebase
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

#if defined(SFR_MCS) && defined(IMF_SAMPLING_MCS)

static MyFloat EnvelopeNorm;
static MyFloat EnvelopeScale;
static MyFloat IMFNorm;

static MyFloat kroupa_single(MyFloat m);
static MyFloat envelope_single(MyFloat m);
static MyFloat inverse_envelope(MyFloat e);
static MyFloat imf_draw();

/*Setup Kroupa IMF sampler.
For simplicity we only model Kroupa above 0.08 Msun.
This does rejection sampling using a power law envelope.
The draws from the envelope uses inverse transform sampling*/
void init_imf(void)
{
  if(All.MinimumIMFStarMass < 0.08)
    terminate("IMF_SAMPLING_MCS: MinimumIMFStarMass is below the lower knee of 0.08 Msun which we do not model currently\n");

  EnvelopeNorm = (-10.0 / 13.0) * (pow(All.MaximumIMFStarMass, -1.3) - pow(All.MinimumIMFStarMass, -1.3));
  IMFNorm      = (-10.0 / 3.0) * (pow(0.5, -0.3) - pow(All.MinimumIMFStarMass, -0.3)) -
            (5.0 / 13.0) * (pow(All.MaximumIMFStarMass, -1.3) - pow(0.5, -1.3));

  EnvelopeScale = 1.05 * kroupa_single(0.5) / envelope_single(0.5);

  if(RestartFlag != 1)
    {
      All.StarMassReservoir = 0.0;

      All.MinimumImportantStellarMass *= SOLAR_MASS / (All.UnitMass_in_g / All.HubbleParam);
    }
}

static MyFloat inverse_envelope(MyFloat e)
{
  MyFloat m = -1.3 * EnvelopeNorm * e + pow(All.MinimumIMFStarMass, -1.3);
  m         = pow(m, (-10.0 / 13.0));

  return m;
}

static MyFloat envelope_single(MyFloat m)
{
  if(m < All.MinimumIMFStarMass)
    return 0;
  if(m > All.MaximumIMFStarMass)
    return 0;

  MyFloat e = pow(m, -2.3) / EnvelopeNorm;

  return e;
}

static MyFloat kroupa_single(MyFloat m)
{
  if(m < All.MinimumIMFStarMass)
    return 0;
  if(m > All.MaximumIMFStarMass)
    return 0;

  if(m < 0.5)
    return pow(m, -1.3) / IMFNorm;
  else
    return pow(m, -2.3) / (2.0 * IMFNorm);
}

static MyFloat imf_draw(void)
{
  MyFloat m, e, f, u;

  int counter = 0;
  while(1)
    {
      /*Sample from envelope function using inverse transform*/
      e = get_random_number();
      m = inverse_envelope(e);

      /*Now accept/reject m*/
      f = kroupa_single(m);
      u = get_random_number();
      u *= envelope_single(m);
      u *= EnvelopeScale;

      if(u < f)
        break;
      counter++;
    }

  m *= SOLAR_MASS / (All.UnitMass_in_g / All.HubbleParam);
  return m;
}

#ifdef IMF_SAMPLING_MCS_VERBOSE
MyFloat sample_imf(MyFloat *results, MyFloat target)
#else
void sample_imf(MyFloat *results, MyFloat target)
#endif
{
  /*Zero the results*/
  memset(results, 0, N_STAR_SLOTS * sizeof(MyFloat));

  /* This is a low res particle that shouldn't be doing any feedback */
  if(target > All.MaxFBStarMass)
#ifndef IMF_SAMPLING_MCS_VERBOSE
    return;
#else
    return 0;
#endif

  MyFloat adjusted_target = target + All.StarMassReservoir;

  if(adjusted_target <= 0)
    {
      All.StarMassReservoir += target;
#ifndef IMF_SAMPLING_MCS_VERBOSE
      return;
#else
      return 0;
#endif
    }
  else
    {
      All.StarMassReservoir = 0;
    }

  MyFloat current_mass = 0;
  int slot             = 0;
  MyFloat draw;
#ifdef IMF_SAMPLING_MCS_VERBOSE
  MyFloat unimportant = 0;
#endif
  while(current_mass < adjusted_target)
    {
      draw = imf_draw();
      current_mass += draw;
      if(draw > All.MinimumImportantStellarMass)
        {
          if(slot == N_STAR_SLOTS)
          {
            printf("IMF_SAMPLING_MCS: target %g adjusted_target %g current_mass %g \n",
              target,adjusted_target,current_mass);
            for(int islot = 0; islot < slot; islot++)
              printf("Slot %d of %d: %g\n",islot,slot,results[islot]);
            terminate("IMF_SAMPLING_MCS: We ran out of slots to store a stellar mass.\n");
          }

          results[slot] = draw;
          slot++;
        }
#ifdef IMF_SAMPLING_MCS_VERBOSE
      else
        unimportant += draw;
#endif
    }

  All.StarMassReservoir = adjusted_target - current_mass;
#ifdef IMF_SAMPLING_MCS_VERBOSE
  return unimportant;
#endif
}

void do_imf_sampling(int old_N_star)
{
  MyFloat *mpart_arr; /*The particle masses*/
  MyFloat *mstar_arr; /*The stellar masses*/
  int n_sample = N_star - old_N_star;
  int counts[NTask];

  MPI_Status status;

  TIMER_STOPSTART(CPU_COOLINGSFR, CPU_IMF_SAMPLE);

  /*Task 0 receives the number of particles that each task requires to be sampled*/
  MPI_Gather(&n_sample, 1, MPI_INT, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /*Task 0 does all the sampling, receiving requests from other tasks */
  if(ThisTask == 0)
    {
      /*Task 0 will do its own sampling first before doing other tasks*/
      if(n_sample > 0)
        {
          for(int i = old_N_star; i < N_star; i++)
            {
#ifndef IMF_SAMPLING_MCS_VERBOSE
              sample_imf(&StarP[i].MassArr[0], P[StarP[i].PID].Mass);
#else
              MyFloat old_res     = All.StarMassReservoir;
              MyFloat unimportant = sample_imf(&StarP[i].MassArr[0], P[StarP[i].PID].Mass);
              mpi_printf(
                  "IMF_SAMPLING_MCS: IMF draw for task %d. t = %g syncpoint = %d oldStarMassReservoir %g target %g adjusted_target %g "
                  "unimportant %g",
                  0, All.Time, All.NumCurrentTiStep, old_res, P[StarP[i].PID].Mass, P[StarP[i].PID].Mass + old_res, unimportant);
              MyFloat total = 0;
              for(int n = 0; n < N_STAR_SLOTS; n++)
                {
                  total += StarP[i].MassArr[n];
                  mpi_printf("slot n %g ", StarP[i].MassArr[n]);
                }
              mpi_printf("total %g overshoot %g newStarMassReservoir %g\n", total, total - P[StarP[i].PID].Mass - old_res,
                         All.StarMassReservoir);
#endif
            }
        }

      /*Now we turn to other tasks*/

      /*mpart_arr must fit the largest task*/
      int n_sample_max = 0;
      for(int task = 1; task < NTask; task++)
        {
          if(counts[task] > n_sample_max)
            n_sample_max = counts[task];
        }

      mpart_arr = (MyFloat *)mymalloc("mpart_arr", n_sample_max * sizeof(MyFloat));
      mstar_arr = (MyFloat *)mymalloc("mstar_arr", N_STAR_SLOTS * n_sample_max * sizeof(MyFloat));

      /* Now do each task in turn */
      for(int task = 1; task < NTask; task++)
        {
          /*Is this task requesting a sample?*/
          if(counts[task] == 0)
            continue;

          /*Receive the masses of the stars the task would like sampling*/
          MPI_Recv(&mpart_arr[0], counts[task], MPI_MYFLOAT, task, TAG_IMF_PART, MPI_COMM_WORLD, &status);

          /*Do the sampling*/
          for(int i = 0; i < counts[task]; i++)
            {
#ifndef IMF_SAMPLING_MCS_VERBOSE
              sample_imf(&mstar_arr[N_STAR_SLOTS * i], mpart_arr[i]);
#else
              MyFloat old_res     = All.StarMassReservoir;
              MyFloat unimportant = sample_imf(&mstar_arr[N_STAR_SLOTS * i], mpart_arr[i]);
              mpi_printf(
                  "IMF_SAMPLING_MCS: IMF draw for task %d. t = %g syncpoint = %d oldStarMassReservoir %g target %g adjusted_target %g "
                  "unimportant %g ",
                  task, All.Time, All.NumCurrentTiStep, old_res, mpart_arr[i], mpart_arr[i] + old_res, unimportant);
              MyFloat total = 0;
              for(int n = 0; n < N_STAR_SLOTS; n++)
                {
                  total += mstar_arr[N_STAR_SLOTS * i + n];
                  mpi_printf("slot n %g ", mstar_arr[N_STAR_SLOTS * i + n]);
                }
              mpi_printf("total %g overshoot %g newStarMassReservoir %g\n", total, -total + mpart_arr[i] + old_res,
                         All.StarMassReservoir);
#endif
            }

          /*Return the results*/
          MPI_Send(&mstar_arr[0], N_STAR_SLOTS * counts[task], MPI_MYFLOAT, task, TAG_IMF_STAR, MPI_COMM_WORLD);
        }

      myfree(mstar_arr);
      myfree(mpart_arr);
    }
  else /*Other tasks*/
    {
      if(n_sample > 0)
        {
          mpart_arr = (MyFloat *)mymalloc("mpart_arr", n_sample * sizeof(MyFloat));
          mstar_arr = (MyFloat *)mymalloc("mstar_arr", N_STAR_SLOTS * n_sample * sizeof(MyFloat));

          for(int i = old_N_star; i < N_star; i++)
            mpart_arr[i - old_N_star] = P[StarP[i].PID].Mass;

          /*Send the particle masses to task 0*/
          MPI_Send(&mpart_arr[0], n_sample, MPI_MYFLOAT, 0, TAG_IMF_PART, MPI_COMM_WORLD);

          /*Receive the individual star particle masses back*/
          MPI_Recv(&mstar_arr[0], N_STAR_SLOTS * n_sample, MPI_MYFLOAT, 0, TAG_IMF_STAR, MPI_COMM_WORLD, &status);

          /*Put them in the right place on the StarP struct*/
          int j = 0;
          for(int i = old_N_star; i < N_star; i++)
            {
              for(int k = 0; k < N_STAR_SLOTS; k++)
                {
                  StarP[i].MassArr[k] = mstar_arr[j];
                  j++;
                }
            }

          myfree(mstar_arr);
          myfree(mpart_arr);
        }
    }

  set_star_properties(old_N_star);

  TIMER_STOPSTART(CPU_IMF_SAMPLE, CPU_COOLINGSFR);
}

void set_star_properties(int old_N_star)
{
  /*Loop over new star particles*/
  for(int i = old_N_star; i < N_star; i++)
    {
      /*Loop over the star slots*/
      for(int j = 0; j < N_STAR_SLOTS; j++)
        {
          /*There is a star in this slot*/
          if(StarP[i].MassArr[j] > 0)
            {
              set_individual_star_properties(i, j);
            }
          else
            {
              StarP[i].LifetimeArr[j] = MAX_REAL_NUMBER;
#ifdef HII_MCS
              StarP[i].S_HiiArr[j] = 0;
#ifdef HII_MCS_LR
              StarP[i].EnergyPerPhotonArr[j] = 0;
#endif
#endif
#ifdef PE_MCS
              StarP[i].L_FUVArr[j] = 0;
#endif
            }
        }
    }
}

void set_individual_star_properties(int auxid, int nslot)
{
  MyFloat mstar = StarP[auxid].MassArr[nslot];

  mstar /= SOLAR_MASS / (All.UnitMass_in_g / All.HubbleParam);  // Now in solar mass

  int i_low;
  MyFloat delta;
  get_star_mass_interp_index(mstar, StarProperties.Masses, StarProperties.N_m, &i_low, &delta);

#if defined(SN_MCS) && defined(SN_MCS_PROMPT)
  /* If star is SN progenitor, make its lifetime ~0. */
  if((StarP[auxid].MassArr[nslot] >= All.SNStarMinMass) && (StarP[auxid].MassArr[nslot] <= All.SNStarMaxMass))
    StarP[auxid].LifetimeArr[nslot] = All.SNLifetime;
  else
#endif
    {
      StarP[auxid].LifetimeArr[nslot] =
          StarProperties.Lifetimes[i_low] + (StarProperties.Lifetimes[i_low + 1] - StarProperties.Lifetimes[i_low]) * delta;
    }

#ifdef HII_MCS
  MyFloat logS_Hii = StarProperties.LogIonizingPhotonRate1e49[i_low] +
                     (StarProperties.LogIonizingPhotonRate1e49[i_low + 1] - StarProperties.LogIonizingPhotonRate1e49[i_low]) * delta;
  if(logS_Hii > -10)
    {
      StarP[auxid].S_HiiArr[nslot] = pow(10.0, logS_Hii);
#ifdef HII_MCS_LR
      StarP[auxid].EnergyPerPhotonArr[nslot] =
          StarProperties.EnergyPerPhoton[i_low] +
          (StarProperties.EnergyPerPhoton[i_low + 1] - StarProperties.EnergyPerPhoton[i_low]) * delta;
#endif
    }
  else
    {
      StarP[auxid].S_HiiArr[nslot] = 0.0;
#ifdef HII_MCS_LR
      StarP[auxid].EnergyPerPhotonArr[nslot] = 0.0;
#endif
    }
#endif

#ifdef PE_MCS
  MyFloat logL_FUV = StarProperties.LogLuminosityFUV[i_low] +
                     (StarProperties.LogLuminosityFUV[i_low + 1] - StarProperties.LogLuminosityFUV[i_low]) * delta;
  if(logL_FUV > 30)
    StarP[auxid].L_FUVArr[nslot] = pow(10.0, logL_FUV);
  else
    StarP[auxid].L_FUVArr[nslot] = 0.0;

#endif
}

void get_star_mass_interp_index(MyFloat m, MyFloat *masses, int N_m, int *i_low, MyFloat *delta)
{
  MyFloat deltam, dm_local;
  int i1, i2;

  if(m > masses[0])
    {
      for(i1 = 0; i1 < N_m - 1 && m > masses[i1 + 1]; i1++)
        ;

      i2 = i1 + 1;

      if(i2 >= N_m)
        i2 = N_m - 1;

      if(m >= masses[0] && m <= masses[N_m - 1])
        dm_local = m - masses[i1];
      else
        dm_local = 0;

      deltam = masses[i2] - masses[i1];

      if(deltam > 0)
        dm_local = dm_local / deltam;
      else
        dm_local = 0;
    }
  else
    {
      i1       = 0;
      i2       = 0;
      dm_local = 0.0;
    }

  *i_low = i1;
  *delta = dm_local;
}

#ifdef SN_MCS
void snII_yields(MyFloat mstar, MyFloat *m_ej, MyFloat *m_Z)
{
  MyFloat mstar_sol = mstar / SOLAR_MASS / (All.UnitMass_in_g / All.HubbleParam);

  int i_low;
  MyFloat delta;
  get_star_mass_interp_index(mstar_sol, StarProperties.Masses, StarProperties.N_m, &i_low, &delta);

  MyFloat ejecta_ratio = StarProperties.EjectaMassRatio[i_low] +
                         (StarProperties.EjectaMassRatio[i_low + 1] - StarProperties.EjectaMassRatio[i_low]) * delta;
  MyFloat m_ej_temp = mstar * ejecta_ratio;

  MyFloat ejecta_Z = StarProperties.EjectaMetallicity[i_low] +
                     (StarProperties.EjectaMetallicity[i_low + 1] - StarProperties.EjectaMetallicity[i_low]) * delta;
  MyFloat m_Z_temp = m_ej_temp * ejecta_Z;

  *m_ej = m_ej_temp;
  *m_Z  = m_Z_temp;
}
#endif

void init_star_properties_tables(void)
{
  /* In the future a file will be read, for now we hardcode */
  StarProperties.N_m = 33; /*The number of mass points to interpolate from */
  /*The masses in solar masses */

  MyFloat tempMasses[33] = {5.0,  5.2,  5.6,  5.8,  6.0,  6.2,  6.4,  7.0,  9.0,  10.0, 11.0, 12.0, 14.0, 16.0, 18.0,  20.0, 24.0,
                            28.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 90.0, 95.0, 100.0, 120.0};

  /*We have to initialize element by element because initialization is happening after declaration*/
  for(int i = 0; i < 33; i++)
    StarProperties.Masses[i] = tempMasses[i];

  /*The lifetimes in years */
  MyFloat tempLifetimes[33] = {99423290.7876, 91315332.4588, 78031681.2869, 72539819.7104, 67728937.6483, 63366685.0301, 59449341.283,
                               49811836.6924, 31539957.0015, 26318102.2777, 22434520.8189, 19609068.5138, 15625200.1076, 13080112.4102,
                               11314328.2251, 10026543.7011, 8287687.77356, 7168695.96885, 6749434.76111, 5957569.89431, 5403636.29958,
                               4992858.91544, 4671951.79729, 4419904.92876, 4209794.03386, 4033260.74174, 3890149.38259, 3765845.19626,
                               3663929.11229, 3476917.64086, 3404682.30599, 3334621.73399, 3126614.953};

  for(int i = 0; i < 33; i++)
    StarProperties.Lifetimes[i] = tempLifetimes[i];

#ifdef SN_MCS
  /*Ratio of progenitor mass to ejecta mass */
  MyFloat tempEjectaMassRatio[33] = {0.9069, 0.9069, 0.9069, 0.9069, 0.9069, 0.9069, 0.9069, 0.9069, 0.9069, 0.9069, 0.9069,
                                     0.9069, 0.9045, 0.9049, 0.9107, 0.9165, 0.9321, 0.9330, 0.9310, 0.9420, 0.9420, 0.9420,
                                     0.9420, 0.9420, 0.9420, 0.9420, 0.9420, 0.9420, 0.9420, 0.9420, 0.9420, 0.9420, 0.9420};

  for(int i = 0; i < 33; i++)
    StarProperties.EjectaMassRatio[i] = tempEjectaMassRatio[i];

  MyFloat tempEjectaMetallicity[33] = {0.0746, 0.0746, 0.0746, 0.0746, 0.0746, 0.0746, 0.0746, 0.0746, 0.0746, 0.0746, 0.0746,
                                       0.0746, 0.0783, 0.0911, 0.1091, 0.1271, 0.1629, 0.1847, 0.1933, 0.2296, 0.2296, 0.2296,
                                       0.2296, 0.2296, 0.2296, 0.2296, 0.2296, 0.2296, 0.2296, 0.2296, 0.2296, 0.2296, 0.2296};

  for(int i = 0; i < 33; i++)
    StarProperties.EjectaMetallicity[i] = tempEjectaMetallicity[i];
#endif

#ifdef HII_MCS
  /*Log of the ionizing photon rate in 1e49 s^-1*/
  MyFloat tempLogIonizingPhotonRate1e49[33] = {-30.0000, -30.0000, -30.0000, -30.0000, -30.0000, -30.0000, -30.0000, -3.8882, -3.0758,
                                               -2.5492,  -2.1753,  -1.8761,  -1.3779,  -1.0848,  -0.8734,  -0.6975,  -0.3855, -0.1869,
                                               -0.1006,  0.0589,   0.1981,   0.3107,   0.4061,   0.4663,   0.5358,   0.6419,  0.7074,
                                               0.7694,   0.8085,   0.9093,   0.9460,   0.9891,   0.9891};

  for(int i = 0; i < 33; i++)
    StarProperties.LogIonizingPhotonRate1e49[i] = tempLogIonizingPhotonRate1e49[i];

#ifdef HII_MCS_LR
  /*Mean energy for photon in ionising bank (< 912 Angstrom) in eV*/
  MyFloat tempEnergyPerPhoton[33] = {16.009,  16.0915, 16.3105, 16.4351, 16.5372, 16.5846, 16.6321, 16.7637, 17.2051, 17.3889, 17.5668,
                                     17.7448, 18.0699, 18.3567, 18.6096, 18.8387, 19.2233, 19.5393, 19.6974, 20.0141, 20.26,   20.5101,
                                     20.7313, 20.9091, 21.0777, 21.2201, 21.3649, 21.4688, 21.5704, 21.7517, 21.8583, 21.9355};

  /*Load and convert to ergs*/
  for(int i = 0; i < 33; i++)
    StarProperties.EnergyPerPhoton[i] = tempEnergyPerPhoton[i] * ELECTRONVOLT_IN_ERGS;
#endif
#endif

#ifdef PE_MCS
  /*Log FUV luminosity in erg s^-1 */
  MyFloat tempLogLuminosityFUV[33] = {35.3972, 35.4746, 35.6293, 35.7067, 35.8413, 35.9828, 36.1242, 36.3111, 36.6763,
                                      36.8011, 36.9071, 36.9954, 37.1444, 37.2823, 37.3951, 37.4809, 37.6726, 37.7798,
                                      37.8213, 37.9219, 38.0135, 38.0821, 38.1441, 38.1732, 38.2209, 38.3457, 38.3892,
                                      38.4363, 38.4806, 38.5536, 38.5915, 38.6242, 38.6242};

  for(int i = 0; i < 33; i++)
    StarProperties.LogLuminosityFUV[i] = tempLogLuminosityFUV[i];
#endif
}

#endif
