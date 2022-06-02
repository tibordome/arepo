/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sfr_mcs/feedback_mcs_utils.c
 * \date        02/2019
 * \author      Matthew C. Smith
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - 17.01.2022 Ported into current codebase
 */

#include <math.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef PE_MCS
#ifndef PE_MCS_PRESHIELD
#define PE_MCS_PRESHIELD 0
#endif
#endif

#if defined SN_MCS_PROMPT && !(defined(SN_MCS) && defined(IMF_SAMPLING_MCS))
#error "SN_MCS_PROMPT requires SN_MCS and IMF_SAMPLING_MCS"
#endif

#if defined(SFR_MCS) && \
    ((defined(SN_MCS) && !defined(SN_MCS_INITIAL_DRIVING)) || (defined(HII_MCS) && !defined(HII_MCS_TEST)) || defined(PE_MCS))

void update_stellar_ages(void)
{
  MyFloat t;
  for(int auxid = 0; auxid < N_star; auxid++)
    {
      t = get_time_difference_in_Gyr(StarP[auxid].BirthTime, All.Time) * 1.0e09;
#ifndef SFR_MCS_DELAY
      StarP[auxid].Age = t;
#else
      StarP[auxid].Age = t - StarP[auxid].TimeDelay;
#endif
    }
}

#ifndef IMF_SAMPLING_MCS
void read_sb99_tables(void)
{
  int i;
  double *tmpdbl;
  hid_t file_id, dataset;
#ifndef SB99_FIXED_Z
  hid_t datatype;
#endif
  char fname[MAXLEN_PATH], setname[MAXLEN_PATH];

#ifndef SB99_FIXED_Z
  int j;
#endif
  file_id = my_H5Fopen(All.SB99TablesPath, H5F_ACC_RDONLY, H5P_DEFAULT);

  /* read number of elements */

#ifdef SN_MCS
  dataset = my_H5Dopen(file_id, "Number_of_timesteps_supernova");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &sb99.N_t_SN, "Number_of_timesteps_supernova");
  my_H5Dclose(dataset, "Number_of_timesteps_supernova");
#endif

#ifdef HII_MCS
  dataset = my_H5Dopen(file_id, "Number_of_timesteps_photons");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &sb99.N_t_Photons, "Number_of_timesteps_photons");
  my_H5Dclose(dataset, "Number_of_timesteps_photons");
#endif

#ifdef PE_MCS
  dataset = my_H5Dopen(file_id, "Number_of_timesteps_fuv");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &sb99.N_t_FUV, "Number_of_timesteps_fuv");
  my_H5Dclose(dataset, "Number_of_timesteps_fuv");
#endif

#ifndef SB99_FIXED_Z
  dataset = my_H5Dopen(file_id, "Number_of_metallicities");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &sb99.N_Z, "Number_of_metallicities");
  my_H5Dclose(dataset, "Number_of_metallicities");
#endif
  mpi_printf("SFR_MCS: Starburst99 header complete...\n");

  /* Allocate arrays */
#ifdef SN_MCS
  sb99.TimestepsSN = (MyFloat *)mymalloc("sb99.TimestepsSN", sb99.N_t_SN * sizeof(MyFloat));
#endif
#ifdef HII_MCS
  sb99.TimestepsPhotons = (MyFloat *)mymalloc("sb99.TimestepsPhotons", sb99.N_t_Photons * sizeof(MyFloat));
#endif
#ifdef PE_MCS
  sb99.TimestepsFUV = (MyFloat *)mymalloc("sb99.TimestepsFUV", sb99.N_t_FUV * sizeof(MyFloat));
#endif

#ifndef SB99_FIXED_Z
  sb99.Metallicities = (MyFloat *)mymalloc("sb99.Metallicities", sb99.N_Z * sizeof(MyFloat));
#ifdef SN_MCS
  sb99.RatesSN = (MyFloat **)mymalloc("sb99.RatesSN", sb99.N_Z * sizeof(MyFloat *));
  for(i = 0; i < sb99.N_Z; i++)
    sb99.RatesSN[i] = (MyFloat *)mymalloc("sb99.RatesSN", sb99.N_t_SN * sizeof(MyFloat));
#endif
#ifdef HII_MCS
  sb99.RatesPhotons = (MyFloat **)mymalloc("sb99.RatesPhotons", sb99.N_Z * sizeof(MyFloat *));
  for(i = 0; i < sb99.N_Z; i++)
    sb99.RatesPhotons[i] = (MyFloat *)mymalloc("sb99.RatesPhotons", sb99.N_t_Photons * sizeof(MyFloat));
#endif
#ifdef PE_MCS
  sb99.LuminosityFUV = (MyFloat **)mymalloc("sb99.LuminosityFUV", sb99.N_Z * sizeof(MyFloat *));
  for(i = 0; i < sb99.N_Z; i++)
    sb99.LuminosityFUV[i] = (MyFloat *)mymalloc("sb99.LuminosityFUV", sb99.N_t_FUV * sizeof(MyFloat));
#endif
#else
#ifdef SN_MCS
  sb99.RatesSN       = (MyFloat *)mymalloc("sb99.RatesSN", sb99.N_t_SN * sizeof(MyFloat));
#ifdef SN_MCS_VARIABLE_EJECTA
  sb99.EjectaMass    = (MyFloat *)mymalloc("sb99.EjectaMass", sb99.N_t_SN * sizeof(MyFloat));
  sb99.ZEjecta       = (MyFloat *)mymalloc("sb99.ZEjecta", sb99.N_t_SN * sizeof(MyFloat));
#endif
#endif
#ifdef HII_MCS
  sb99.RatesPhotons  = (MyFloat *)mymalloc("sb99.RatesPhotons", sb99.N_t_Photons * sizeof(MyFloat));
#endif
#ifdef PE_MCS
  sb99.LuminosityFUV = (MyFloat *)mymalloc("sb99.LuminosityFUV", sb99.N_t_FUV * sizeof(MyFloat));
#endif
#endif

  mpi_printf("SFR_MCS: Starburst99 arrays allocated...\n");

#ifndef SB99_FIXED_Z
  tmpdbl  = (double *)mymalloc("Metallicities", sb99.N_Z * sizeof(double));
  dataset = my_H5Dopen(file_id, "Metallicities");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Metallicities");
  my_H5Dclose(dataset, "Metallicities");
  for(i = 0; i < sb99.N_Z; i++)
    sb99.Metallicities[i] = tmpdbl[i];
  myfree(tmpdbl);

  char *tempname[sb99.N_Z];
  datatype = my_H5Tcopy(H5T_C_S1);
  my_H5Tset_size(datatype, H5T_VARIABLE);
  dataset = my_H5Dopen(file_id, "Metallicity_Names");
  my_H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempname, "Metallicity_Names");
  my_H5Dclose(dataset, "Metallicity_Names");
  my_H5Tclose(datatype);

#endif

#ifdef SN_MCS
  tmpdbl  = (double *)mymalloc("TimestepsSN", sb99.N_t_SN * sizeof(double));
  dataset = my_H5Dopen(file_id, "TimestepsSN");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "TimestepsSN");
  my_H5Dclose(dataset, "TimestepsSN");
  for(i = 0; i < sb99.N_t_SN; i++)
    sb99.TimestepsSN[i] = tmpdbl[i];
  myfree(tmpdbl);

  sb99.t_min_SN = sb99.TimestepsSN[0];
  sb99.t_max_SN = sb99.TimestepsSN[sb99.N_t_SN - 1];

#ifndef SB99_FIXED_Z

  tmpdbl = (double *)mymalloc("RatesSN", sb99.N_t_SN * sizeof(double));
  for(i = 0; i < sb99.N_Z; i++)
    {
      sprintf(setname, "/RatesSN/%s", tempname[i]);
      dataset = my_H5Dopen(file_id, setname);
      my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, setname);
      my_H5Dclose(dataset, setname);

      for(j = 0; j < sb99.N_t_SN; j++)
        sb99.RatesSN[i][j] = tmpdbl[j];
    }
  myfree(tmpdbl);

#else
  tmpdbl = (double *)mymalloc("RatesSN", sb99.N_t_SN * sizeof(double));
  sprintf(setname, "/RatesSN/%s", All.SB99_metallicity_name);
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, setname);
  my_H5Dclose(dataset, setname);
  for(i = 0; i < sb99.N_t_SN; i++)
    sb99.RatesSN[i] = tmpdbl[i];
  myfree(tmpdbl);

#ifdef SN_MCS_VARIABLE_EJECTA
  tmpdbl = (double *)mymalloc("EjectaMass", sb99.N_t_SN * sizeof(double));
  sprintf(setname, "/EjectaMass/%s", All.SB99_metallicity_name);
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, setname);
  my_H5Dclose(dataset, setname);
  for(i = 0; i < sb99.N_t_SN; i++)
    sb99.EjectaMass[i] = tmpdbl[i];
  myfree(tmpdbl);

  tmpdbl = (double *)mymalloc("ZEjecta", sb99.N_t_SN * sizeof(double));
  sprintf(setname, "/ZEjecta/%s", All.SB99_metallicity_name);
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, setname);
  my_H5Dclose(dataset, setname);
  for(i = 0; i < sb99.N_t_SN; i++)
    sb99.ZEjecta[i] = tmpdbl[i];
  myfree(tmpdbl);
#endif

#endif
#endif  // SN_MCS

#ifdef HII_MCS
  tmpdbl  = (double *)mymalloc("TimestepsPhotons", sb99.N_t_Photons * sizeof(double));
  dataset = my_H5Dopen(file_id, "TimestepsPhotons");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "TimestepsPhotons");
  my_H5Dclose(dataset, "TimestepsPhotons");
  for(i = 0; i < sb99.N_t_Photons; i++)
    sb99.TimestepsPhotons[i] = tmpdbl[i];
  myfree(tmpdbl);

#ifndef SB99_FIXED_Z

  tmpdbl = (double *)mymalloc("RatesPhotons", sb99.N_t_Photons * sizeof(double));
  for(i = 0; i < sb99.N_Z; i++)
    {
      sprintf(setname, "/RatesPhotons/%s", tempname[i]);
      dataset = my_H5Dopen(file_id, setname);
      my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, setname);
      my_H5Dclose(dataset, setname);

      for(j = 0; j < sb99.N_t_Photons; j++)
        sb99.RatesPhotons[i][j] = tmpdbl[j];
    }
  myfree(tmpdbl);
#else
  tmpdbl = (double *)mymalloc("RatesPhotons", sb99.N_t_Photons * sizeof(double));
  sprintf(setname, "/RatesPhotons/%s", All.SB99_metallicity_name);
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, setname);
  my_H5Dclose(dataset, setname);
  for(i = 0; i < sb99.N_t_Photons; i++)
    sb99.RatesPhotons[i] = tmpdbl[i];
  myfree(tmpdbl);
#endif
#endif  // HII_MCS

#ifdef PE_MCS
  tmpdbl  = (double *)mymalloc("TimestepsFUV", sb99.N_t_FUV * sizeof(double));
  dataset = my_H5Dopen(file_id, "TimestepsFUV");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "TimestepsFUV");
  my_H5Dclose(dataset, "TimestepsFUV");
  for(i = 0; i < sb99.N_t_FUV; i++)
    sb99.TimestepsFUV[i] = tmpdbl[i];
  myfree(tmpdbl);

#ifndef SB99_FIXED_Z

  tmpdbl = (double *)mymalloc("LuminosityFUV", sb99.N_t_FUV * sizeof(double));
  for(i = 0; i < sb99.N_Z; i++)
    {
      sprintf(setname, "/LuminosityFUV/%s", tempname[i]);
      dataset = my_H5Dopen(file_id, setname);
      my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, setname);
      my_H5Dclose(dataset, setname);

      for(j = 0; j < sb99.N_t_FUV; j++)
        sb99.LuminosityFUV[i][j] = tmpdbl[j];
    }
  myfree(tmpdbl);
#else
  tmpdbl = (double *)mymalloc("LuminosityFUV", sb99.N_t_FUV * sizeof(double));
  sprintf(setname, "/LuminosityFUV/%s", All.SB99_metallicity_name);
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, setname);
  my_H5Dclose(dataset, setname);
  for(i = 0; i < sb99.N_t_FUV; i++)
    sb99.LuminosityFUV[i] = tmpdbl[i];
  myfree(tmpdbl);
#endif
#endif  // PE_MCS

  my_H5Fclose(file_id, fname);

  mpi_printf("SFR_MCS: Starburst99 setup complete...\n");
}

#ifndef SB99_FIXED_Z
int get_sb99_z_index(MyFloat metallicity)
{
  int left;
  MyFloat middle_val;

  left = 0;

  while(left < (sb99.N_Z - 1))
    {
      middle_val = (sb99.Metallicities[left] + sb99.Metallicities[left + 1]) / 2.0;
      if(metallicity < middle_val)
        return left;

      left++;
    }
  return left;
}
#endif

void get_sb99_t_indicies(MyFloat t, MyFloat *timesteps, int N_t, int *it_low, int *it_high, MyFloat *delta)
{
  MyFloat deltat, dt_local;
  int i1, i2;

  if(t > timesteps[0])
    {
      for(i1 = 0; i1 < N_t - 1 && t > timesteps[i1 + 1]; i1++)
        ;

      i2 = i1 + 1;

      if(i2 >= N_t)
        i2 = N_t - 1;

      if(t >= timesteps[0] && t <= timesteps[N_t - 1])
        dt_local = t - timesteps[i1];
      else
        dt_local = 0;

      deltat = timesteps[i2] - timesteps[i1];

      if(delta > 0)
        dt_local = dt_local / deltat;
      else
        dt_local = 0;
    }
  else
    {
      i1       = 0;
      i2       = 0;
      dt_local = 0.0;
    }

  *it_low  = i1;
  *it_high = i2;
  *delta   = dt_local;
}

#ifdef SN_MCS
#ifndef SN_MCS_SINGLE_INJECTION
#ifndef SB99_FIXED_Z
double get_sn_rate(MyFloat t, int iz)
{
  int it_low, it_high;
  double delta, rate;

  get_sb99_t_indicies(t, sb99.TimestepsSN, sb99.N_t_SN, &it_low, &it_high, &delta);
  rate = sb99.RatesSN[iz][it_low] + (sb99.RatesSN[iz][it_high] - sb99.RatesSN[iz][it_low]) * delta;
  return rate;
}
#else
double get_sn_rate(MyFloat t)
{
  int it_low, it_high;
  double delta, rate;

  get_sb99_t_indicies(t, sb99.TimestepsSN, sb99.N_t_SN, &it_low, &it_high, &delta);
  rate = sb99.RatesSN[it_low] + (sb99.RatesSN[it_high] - sb99.RatesSN[it_low]) * delta;
  return rate;
}

#ifdef SN_MCS_VARIABLE_EJECTA
/*Returns ejecta mass in solar masses */
void get_ejecta_properties(MyFloat t, double *m_ej, double *z_ej)
{
  int it_low, it_high;
  double delta;

  get_sb99_t_indicies(t, sb99.TimestepsSN, sb99.N_t_SN, &it_low, &it_high, &delta);

  *m_ej = sb99.EjectaMass[it_low] + (sb99.EjectaMass[it_high] - sb99.EjectaMass[it_low]) * delta;
  *z_ej = sb99.ZEjecta[it_low] + (sb99.ZEjecta[it_high] - sb99.ZEjecta[it_low]) * delta;
}
#endif  // SN_MCS_VARIABLE_EJECTA
#endif
#endif  /*ifndef SN_MCS_SINGLE_INJECTION */
#endif  // SN_MCS

#if defined(HII_MCS) && !defined(HII_MCS_TEST)
/* Sets StarP[].S_Hii and returns number of star particles with S_Hii > 0 */
int update_photon_rates(void)
{
  int auxid, it_high, nactive;
#ifndef SB99_FIXED_Z
  int iz;
#endif
  MyFloat t, deltat, dt_local, lograte;

  for(auxid = 0, nactive = 0; auxid < N_star; auxid++)
    {
      if(All.ComovingIntegrationOn && (P[StarP[auxid].PID].Mass > All.MaxFBStarMass))
        {
          StarP[auxid].S_Hii = 0.0;
          continue;
        }

      it_high = StarP[auxid].photon_it_high;

      if(it_high >= sb99.N_t_Photons)
        {
          StarP[auxid].S_Hii = 0.0;
          continue;
        }

      t = StarP[auxid].Age;

      if(t > sb99.TimestepsPhotons[it_high])
        {
          for(; (it_high < sb99.N_t_Photons) && (t > sb99.TimestepsPhotons[it_high]); it_high++)
            ;

          StarP[auxid].photon_it_high = it_high;

          if(it_high >= sb99.N_t_Photons)
            {
              StarP[auxid].S_Hii = 0.0;
              continue;
            }
        }

      deltat   = sb99.TimestepsPhotons[it_high] - sb99.TimestepsPhotons[it_high - 1];
      dt_local = t - sb99.TimestepsPhotons[it_high - 1];
      dt_local = dt_local / deltat;
#ifndef SB99_FIXED_Z
      iz      = StarP[auxid].iz;
      lograte = sb99.RatesPhotons[iz][it_high - 1] + (sb99.RatesPhotons[iz][it_high] - sb99.RatesPhotons[iz][it_high - 1]) * dt_local;
#else
      lograte = sb99.RatesPhotons[it_high - 1] + (sb99.RatesPhotons[it_high] - sb99.RatesPhotons[it_high - 1]) * dt_local;
#endif
      if(lograte > 0)
        {
          /* Note rescaling by 1e49 to avoid issues when MyFloat = float */
          StarP[auxid].S_Hii =
              pow(10.0, lograte - 49.0) * StarP[auxid].InitialMass * (All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS));
          StarP[auxid].S_Hii *= All.HiiAttenFac;
          nactive++;
        }
      else
        StarP[auxid].S_Hii = 0.0;
    }

  return nactive;
}
#endif  // HII_MCS

#ifdef PE_MCS
/* Sets StarP[].L_FUV */
void update_FUV_luminosities(void)
{
  TIMER_START(CPU_PE_FEEDBACK);

  int it_high, i;
#ifndef SB99_FIXED_Z
  int iz;
#endif
  MyFloat t, deltat, dt_local, lograte;

#if(PE_MCS_PRESHIELD == 1)
  dust_column_search_data *searchdata = (dust_column_search_data *)mymalloc("searchdata", N_star * sizeof(dust_column_search_data));
  int nsearch                         = 0;
  MyFloat dustcol;
#endif

  for(int idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type != 4)
        continue;

      if(P[i].ID == 0 && P[i].Mass == 0)
        {
          STP(i).L_FUV = 0.0;
          continue;
        }

      if(All.ComovingIntegrationOn && (P[i].Mass > All.MaxFBStarMass))
        {
          STP(i).L_FUV = 0.0;
          continue;
        }

      it_high = STP(i).fuv_it_high;

      if(it_high >= sb99.N_t_FUV)
        {
          STP(i).L_FUV = 0.0;
          continue;
        }

      t = STP(i).Age;

      if(t > sb99.TimestepsFUV[it_high])
        {
          for(; (it_high < sb99.N_t_FUV) && (t > sb99.TimestepsFUV[it_high]); it_high++)
            ;

          STP(i).fuv_it_high = it_high;

          if(it_high >= sb99.N_t_FUV)
            {
              STP(i).L_FUV = 0.0;
              continue;
            }
        }

      deltat   = sb99.TimestepsFUV[it_high] - sb99.TimestepsFUV[it_high - 1];
      dt_local = t - sb99.TimestepsFUV[it_high - 1];
      dt_local = dt_local / deltat;
#ifndef SB99_FIXED_Z
      iz = STP(i).iz;
      lograte =
          sb99.LuminosityFUV[iz][it_high - 1] + (sb99.LuminosityFUV[iz][it_high] - sb99.LuminosityFUV[iz][it_high - 1]) * dt_local;
#else
      lograte = sb99.LuminosityFUV[it_high - 1] + (sb99.LuminosityFUV[it_high] - sb99.LuminosityFUV[it_high - 1]) * dt_local;
#endif
      if(lograte > -15)
        {
          STP(i).L_FUV = pow(10.0, lograte) * STP(i).InitialMass * (All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS));
#if(PE_MCS_PRESHIELD == 1)
          searchdata[nsearch].Pos[0] = P[i].Pos[0];
          searchdata[nsearch].Pos[1] = P[i].Pos[1];
          searchdata[nsearch].Pos[2] = P[i].Pos[2];
          nsearch++;
#endif
        }
      else
        STP(i).L_FUV = 0.0;
    }

#if(PE_MCS_PRESHIELD == 1)
  estimate_dust_column(searchdata, nsearch, 0);

  for(int idx = 0, nsearch = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type != 4)
        continue;

      if(STP(i).L_FUV > 0.0)
        {
          dustcol = searchdata[nsearch].DustColumn;
          if(!isfinite(dustcol))
            {
              printf("PE_MCS DEBUG!!! Sync-point %d Task %d ID %d i %d dustcol %g nsearch %d N_star %d L_FUV %g\n",
                     All.NumCurrentTiStep, ThisTask, P[i].ID, i, dustcol, nsearch, N_star, STP(i).L_FUV);
              terminate("PE_MCS Encountered non-finited dustcol, terminating\n");
            }
          STP(i).L_FUV *= exp((-1.33e-21) * searchdata[nsearch].DustColumn);
#ifdef PE_MCS_STORE_DUST_COLUMN
          STP(i).DustColumn = searchdata[nsearch].DustColumn;
#endif
          nsearch++;
        }
    }

  myfree(searchdata);
#endif

  TIMER_STOP(CPU_PE_FEEDBACK);
}
#endif  // PE_MCS

#else  // IMF_SAMPLING_MCS

void check_for_dead_stars(void)
{
  TIMER_START(CPU_IMF_FEEDBACK);

  MyFloat age;
  int i, auxid;
#ifdef SN_MCS
  MyFloat m_ej, m_ej_temp;
  MyFloat m_Z, m_Z_temp;
  int n_sn_i;
  int n_sn = 0;
  int n_sn_tot;
  NumSNLocal = NumSNGlobal = 0;
#endif

  for(int idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type != 4)
        continue;

      if((P[i].Mass == 0) && (P[i].ID == 0))
        continue;

      auxid = P[i].AuxDataID;

      age = StarP[auxid].Age;
#ifdef SN_MCS
      m_ej = 0;
      m_Z = 0;
      n_sn_i = 0;
#endif
      for(int n = 0; n < N_STAR_SLOTS; n++)
        {
          if(StarP[auxid].LifetimeArr[n] < age)
            {
#ifdef HII_MCS
              StarP[auxid].S_HiiArr[n] = 0;
#ifdef HII_MCS_LR
              StarP[auxid].EnergyPerPhotonArr[n] = 0;
#endif
#endif
#ifdef PE_MCS
              StarP[auxid].L_FUVArr[n] = 0;
#endif
#ifdef SN_MCS
              if((StarP[auxid].MassArr[n] >= All.SNStarMinMass) && (StarP[auxid].MassArr[n] <= All.SNStarMaxMass))
                {
#ifdef SN_MCS_CHANCE
                  if(get_random_number() < (1.0 / SN_MCS_CHANCE))
                    {
#endif
                      snII_yields(StarP[auxid].MassArr[n], &m_ej_temp, &m_Z_temp);
                      m_ej += m_ej_temp;
                      m_Z += m_Z_temp;
                      n_sn_i++;
                      n_sn++;
#ifdef SN_MCS_CHANCE
                    }
#endif
                }
#endif
              StarP[auxid].LifetimeArr[n] = MAX_REAL_NUMBER;
            }
        }
#ifdef SN_MCS
      if(n_sn_i > 0)
        {
          NumSNLocal++;
          StarP[auxid].N_SN = n_sn_i;
          StarP[auxid].N_SN_cum += n_sn_i;
          StarP[auxid].Z_ej = m_Z / m_ej;
          m_ej = fmin(m_ej, P[i].Mass);
          StarP[auxid].M_ej = m_ej;
          P[i].Mass -= StarP[auxid].M_ej;
          if(P[i].Mass < 0.01 * StarP[auxid].InitialMass)
            {
              StarP[auxid].M_ej += P[i].Mass;
              P[i].Mass = 0;
            }
        }
      else
        {
          StarP[auxid].N_SN = 0;
          StarP[auxid].M_ej = 0;
        }
#endif
    }
#ifdef SN_MCS
  MPI_Reduce(&n_sn, &n_sn_tot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);          /*This is used for log, so only Task 0 needs this */
  MPI_Allreduce(&NumSNLocal, &NumSNGlobal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); /*Every task needs this to trigger SN routines */

  if(ThisTask == 0)
    {
      MyFloat snr_global;
      if(All.TimeStep > 0)
        snr_global = n_sn_tot / (All.TimeStep / All.cf_time_hubble_a);
      else
        snr_global = 0.0;
      fprintf(FdSnr, "%14e %14i %14i %14e\n", All.Time, n_sn_tot, NumSNGlobal, snr_global);
      myflush(FdSnr);
    }

#endif

  TIMER_STOP(CPU_IMF_FEEDBACK);
}

#if defined(HII_MCS) && !defined(HII_MCS_TEST)
/* Sets StarP[].S_Hii and returns number of star particles with S_Hii > 0 */
int update_photon_rates(void)
{
  int auxid, nactive;

  for(auxid = 0, nactive = 0; auxid < N_star; auxid++)
    {
      StarP[auxid].S_Hii = 0;
#ifdef HII_MCS_LR
      StarP[auxid].EnergyPerPhoton = 0;
#endif

#ifdef SFR_MCS_DELAY
      if(StarP[auxid].Age < 0)  // star not active yet
        continue;
#endif

      for(int n = 0; n < N_STAR_SLOTS; n++)
        {
          StarP[auxid].S_Hii += StarP[auxid].S_HiiArr[n] * All.HiiAttenFac;
#ifdef HII_MCS_LR
          StarP[auxid].EnergyPerPhoton += StarP[auxid].S_HiiArr[n] * All.HiiAttenFac * StarP[auxid].EnergyPerPhotonArr[n];
#endif
        }

      if(StarP[auxid].S_Hii > 0)
        nactive++;

#ifdef HII_MCS_LR
      StarP[auxid].EnergyPerPhoton /= StarP[auxid].S_Hii;
#endif
    }

  return nactive;
}
#endif

#ifdef PE_MCS
/* Sets StarP[].L_FUV */
void update_FUV_luminosities(void)
{
  TIMER_START(CPU_PE_FEEDBACK);

  int i;

#if(PE_MCS_PRESHIELD == 1)
  dust_column_search_data *searchdata = (dust_column_search_data *)mymalloc("searchdata", N_star * sizeof(dust_column_search_data));
  int nsearch = 0;
  MyFloat dustcol;
#endif

  for(int idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type != 4)
        continue;

      STP(i).L_FUV = 0;

#ifdef SFR_MCS_DELAY
      if(STP(i).Age < 0)  // This star is not active yet.
        continue;
#endif

      for(int n = 0; n < N_STAR_SLOTS; n++)
        STP(i).L_FUV += STP(i).L_FUVArr[n];

#if(PE_MCS_PRESHIELD == 1)
      if(STP(i).L_FUV > 0)
        {
          searchdata[nsearch].Pos[0] = P[i].Pos[0];
          searchdata[nsearch].Pos[1] = P[i].Pos[1];
          searchdata[nsearch].Pos[2] = P[i].Pos[2];
          nsearch++;
        }
#endif
    }

#if(PE_MCS_PRESHIELD == 1)
  estimate_dust_column(searchdata, nsearch, 0);

  for(int idx = 0, nsearch = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type != 4)
        continue;

      if(STP(i).L_FUV > 0.0)
        {
          dustcol = searchdata[nsearch].DustColumn;
          if(!isfinite(dustcol))
            {
              printf("PE_MCS DEBUG!!! Sync-point %d Task %d ID %d i %d dustcol %g nsearch %d N_star %d L_FUV %g\n",
                     All.NumCurrentTiStep, ThisTask, P[i].ID, i, dustcol, nsearch, N_star, STP(i).L_FUV);
              terminate("PE_MCS Encountered non-finited dustcol, terminating\n");
            }
          STP(i).L_FUV *= exp((-1.33e-21) * searchdata[nsearch].DustColumn);
#ifdef PE_MCS_STORE_DUST_COLUMN
          STP(i).DustColumn = searchdata[nsearch].DustColumn;
#endif
          nsearch++;
        }
    }

  myfree(searchdata);
#endif

  TIMER_STOP(CPU_PE_FEEDBACK);
}
#endif  // PE_MCS
#endif  // IMF_SAMPLING_MCS
#endif
