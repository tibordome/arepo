/*! * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/MRT/RT.c
 * \date        MM/YYYY
 * \author      Rahul Kannan, Federico Marincaci, Mark Vogelsberger
 * \brief       main driver for a moment based RT with the M1 closure.
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <math.h>
#include <mpi.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../domain.h"
#include "../proto.h"
#include "../voronoi.h"

#include "RT_proto.h"

#ifdef MRT

#if defined(MRT_RIEMANN_HLLE) || defined(MRT_RIEMANN_HLLE_NEW)
double lambda1[101][101];
double lambda2[101][101];
double lambda3[101][101];
double lambda4[101][101];
#endif

#ifdef MRT_SINGLE_STAR

static int Nlines;
static struct stellarparameters starpar;

static double interpolate_stellar_params(int bin, double age);

#endif

/* function to select the type of chemistry that is used for the simulation */

void mrt_update_chemistry(void)
{
  TIMER_START(CPU_RT_CHEM);
#if defined(MRT_CHEMISTRY_PS2009) && defined(MRT_CHEMISTRY_PS2011)
  terminate("Only one of MRT_CHEMISTRY_2009 or MRT_CHEMISTRY_2011 can be defined\n");
#endif

#if defined(MRT_CHEMISTRY_PS2009) && defined(MRT_UV_ONLY_DUST)
  terminate("Only one of MRT_CHEMISTRY_2009/2011 or MRT_UV_DUST_ONLY can be defined\n");
#endif

#ifndef MRT_NO_UV
#ifdef MRT_CHEMISTRY_PS2011
  mrt_update_chemistry_ps2011();
  mpi_printf_rt(1, "RT: Updated Chemistry with PS2011\n");
#endif
#ifdef MRT_CHEMISTRY_PS2009
  mrt_update_chemistry_ps2009();
  mpi_printf_rt(1, "RT: Updated Chemistry with PS2009\n");
#endif
#endif

#ifdef MRT_UV_ONLY_DUST
  mrt_update_chemistry_only_dust();
  mpi_printf_rt(1, "RT: Updated UV Dust only chemistry\n");
#endif
  TIMER_STOP(CPU_RT_CHEM);
}

/* initialization function called at the start of simulation - sets the values of the necessary global variables */

void init_RT(void)
{
  mpi_printf_rt(1, "BEGRUN RT: UV bins = %d \t IR bins = %d \t Total bins = %d \t Number of gradients = %d\n", UV_BINS, IR_BINS,
                MRT_BINS, COUNT_GRAD_MRT);
#ifdef MRT_SUBCYCLE
  if(All.RTNumSubCycles < 2)
    terminate("MRT_SUBCYCE requires atleast two subcyles per hydro step \n");
  if(All.RTNumSubCycles % 2 != 0)
    terminate("Please specify the number of subcycle as a power of 2 to maximise efficiency\n");
#endif
#ifdef MRT_CONSTANT_KAPPA
  c_internal_units = 1.0;
#else
  c_internal_units = CLIGHT / All.UnitVelocity_in_cm_per_s;
#endif

#ifdef MRT_IR
  radiation_constant = RADIATION_CONSTANT * pow(All.UnitLength_in_cm, 3) / All.UnitEnergy_in_cgs;
  mpi_printf_rt(1, "RT: INIT : Radiation Constant = %le\n", radiation_constant);
#endif

#ifdef MRT_RIEMANN_HLLE
  readin_hlle_eingenvalues();
#endif

#ifdef MRT_SINGLE_STAR
  read_stellar_table();
#endif

#ifdef MRT_IR_GRAIN_KAPPA
  mpi_printf_rt(1, "MRT: Reading grain kappa data.\n");
  read_grain_kappa_data();
#endif

  mpi_printf_rt(1, "MRT: CLIGHT = %le\n", c_internal_units);
}

/* function call to inject photons if MRT_LOCAL_FEEDBACK option is used */
#ifdef MRT_LOCAL_FEEDBACK
double inject_photons_from_star(int num, double dt, double Age)
{
#ifdef MRT_SINGLE_STAR
  return dt * All.UnitTime_in_s * interpolate_stellar_params(num, Age);
#else
  double Nphotons, IREnergy;
  if(Age <= 3.0)
    {
      Nphotons = 3.0e-14;
      IREnergy = 4.14e35;
      lum[0]   = 0.671394060961;
      lum[1]   = 0.136192500759;
      lum[2]   = 0.0294087494618;
      lum[3]   = 0.119223360778;
      lum[4]   = 0.0437813280392;
    }
  else
    {
      Nphotons = 2.2e-14;
      IREnergy = 1.75e36;
      lum[0]   = 0.892053121483;
      lum[1]   = 0.0824284042731;
      lum[2]   = 0.00424745279675;
      lum[3]   = 0.0123438314721;
      lum[4]   = 0.0089271899748;
    }

  if(num < UV_BINS)
    return lum[num] * dt * All.UnitTime_in_s * Nphotons;
  else
    return dt * All.UnitTime_in_s * IREnergy / All.UnitEnergy_in_cgs;
#endif
}
#endif

/* Rudimentary implementation for adding fluxes from sources (stars/BHs) */
void add_source_fluxes(void) /* Correct this */
{
#ifdef MRT_SOURCES

#if defined(MRT_STARS) || defined(MRT_BH)
  mpi_printf_rt(1, "RT: Adding source fluxes....\n");

#ifdef MRT_STARS
  do_ionizing_stellar_sources();
#endif

#ifdef MRT_BH
#ifdef MRT_BH_PULSED
  if(((int)(All.Time * 1000 / All.AGNPulseTime) % 2) == 0)
    {
      mpi_printf_rt(1, "MRT_BH_PULSED: BH radiation ON time=%g\n", All.Time);
      do_ionizing_blackhole_sources();
    }
  else
    mpi_printf_rt(1, "MRT_BH_PULSED: BH radiation OFF time=%g\n", All.Time);

#else
  do_ionizing_blackhole_sources();
#endif
#endif

#endif
#else

  //#ifndef MRT_NO_UV

  if(!TimeBinSynchronized[All.HighestOccupiedTimeBin])
    return;

  mpi_printf_rt(1, "RT: Adding source fluxes....\n");

  double a4 = 1;
  if(All.ComovingIntegrationOn)
    a4 = All.Time * All.Time * All.Time * All.Time;

  double dt = 0.0;

  int i, j, idx;

  double X, Y, Z, radius, xcm, ycm, zcm, Xcm, Ycm, Zcm;

  int nincells, Nincells;

  xcm = ycm = zcm = 0.0;
  nincells        = 0;

  dt = (All.HighestActiveTimeBin ? (((integertime)1) << All.HighestActiveTimeBin) : 0) * All.Timebase_interval;
  if(All.ComovingIntegrationOn)
    dt /= All.cf_hubble_a;

  for(i = 0; i < NumGas; i++)
    {
      X      = P[i].Pos[0] - All.BoxSize / 2.0;
      Y      = P[i].Pos[1] - All.BoxSize / 2.0;
#ifndef TWODIMS
      Z      = P[i].Pos[2] - All.BoxSize / 2.0;
      radius = sqrt(X * X + Y * Y + Z * Z);
#else
      radius = sqrt(X * X + Y * Y);
#endif

      if(radius < 2.0)
        {
          nincells++;
          xcm += P[i].Pos[0];
          ycm += P[i].Pos[1];
          zcm += P[i].Pos[2];
        }
      //}
    }

  MPI_Allreduce(&nincells, &Nincells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&xcm, &Xcm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&ycm, &Ycm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&zcm, &Zcm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  mpi_printf_rt(0, "RT: %d number of cells injected with \t", Nincells);
  mpi_printf_rt(0, "Center of mass of cells x = %g, y = %g, z = %g\n", Xcm / ((double)Nincells), Ycm / ((double)Nincells),
                Zcm / ((double)Nincells));

  for(i = 0; i < NumGas; i++)
    {
      X      = P[i].Pos[0] - All.BoxSize / 2.0;
      Y      = P[i].Pos[1] - All.BoxSize / 2.0;
#ifndef TWODIMS
      Z      = P[i].Pos[2] - All.BoxSize / 2.0;
#endif

#ifndef TWODIMS
      radius = sqrt(X * X + Y * Y + Z * Z);
#else
      radius = sqrt(X * X + Y * Y);
#endif

      if(radius < 2.0)
        {
          for(int num1 = 0; num1 < MRT_BINS; num1++)
            {
              //      if(num1==2)
              SphP[i].Cons_DensPhot[num1] += lum[num1] * dt * All.UnitTime_in_s / All.HubbleParam * 5e-15 / ((double)(Nincells));
            }
          // else if(num1==1)
          //  SphP[i].Cons_DensPhot[num1] +=  dt * All.UnitTime_in_s * 3e-15 / ((double) (Nincells)) ;
        }
    }

#endif
}

/* calculate the divergence of the photn flux - needed for the prediction step in the finite volume solver */
#ifdef MRT_TIME_EXTRAPOLATION
void calculate_div_F(struct state *dl, struct state *st)
{
  struct rt_grad_data *rtgrad = st->rtgrad;

  int num1, j;
  for(num1 = 0; num1 < MRT_BINS; num1++)
    {
      double trace  = 0.0;
      double top    = st->DensPhot[num1];
      double bottom = 1.0 / st->DensPhot[num1];
      {
        for(j = 0; j < 3; j++)
          trace += top * rtgrad->dFN[num1][j][j] + bottom * st->RT_F[num1][j] * rtgrad->dDensPhot[num1][j];
      }
      dl->divF[num1] = trace;
    }
  return;
}
#endif

void RT_initialize_cell(int i)
{
#ifdef MRT_UPDATE_AT_END_OF_STEP
#if defined(MRT_COOLING_HEATING) || defined(MRT_CHEM_SG)
  SphP[i].RT_utherm  = SphP[i].Utherm;
  SphP[i].RT_dutherm = 0.0;
#endif

#ifdef MRT_RADIATION_PRESSURE
  SphP[i].RT_mominj[0] = SphP[i].RT_mominj[1] = SphP[i].RT_mominj[2] = 0.0;
#endif
#endif

  double nH_times_volume = P[i].Mass;

  double y_fac  = (1.0 - HYDROGEN_MASSFRAC) / 4.0 / HYDROGEN_MASSFRAC;
  double nh_fac = 1.0 / (1.0 + y_fac);

  if(RestartFlag == RESTART_IC)
    {
      SphP[i].HI    = 0.999999;
      SphP[i].HII   = 0.000001;
      SphP[i].HeI   = 0.999998 * y_fac;
      SphP[i].HeII  = 0.000001 * y_fac;
      SphP[i].HeIII = 0.000001 * y_fac;
    }

  SphP[i].nHI    = SphP[i].HI * nH_times_volume;
  SphP[i].nHII   = SphP[i].HII * nH_times_volume;
  SphP[i].nHeI   = SphP[i].HeI * nH_times_volume;
  SphP[i].nHeII  = SphP[i].HeII * nH_times_volume;
  SphP[i].nHeIII = SphP[i].HeIII * nH_times_volume;
  SphP[i].Ne     = SphP[i].HII + SphP[i].HeII + 2.0 * SphP[i].HeIII;
  SphP[i].ne     = SphP[i].Ne * nH_times_volume;

#ifdef MRT_INIT_IONIZATION
  SetInitGasState(i);
#endif

#ifdef MRT_COMOVING
  SphP[i].Old_Vel[0] = P[i].Vel[0];
  SphP[i].Old_Vel[1] = P[i].Vel[1];
  SphP[i].Old_Vel[2] = P[i].Vel[2];
#endif

#ifndef MRT_LEVITATION_TEST
  for(int num1 = 0; num1 < MRT_BINS; num1++)
    {
      if(RestartFlag == RESTART_IC)
        SphP[i].DensPhot[num1] = MINDENSPHOT;
      SphP[i].Cons_DensPhot[num1] = SphP[i].DensPhot[num1] * SphP[i].Volume;

      if(isnan(SphP[i].Cons_DensPhot[num1]))
        terminate("Something went wrong\n");

      for(int j = 0; j < 3; j++)
        {
          if(RestartFlag == RESTART_IC)
            SphP[i].RT_F[num1][j] = MINDENSPHOT;
          SphP[i].Cons_RT_F[num1][j] = SphP[i].RT_F[num1][j] * SphP[i].Volume;
        }
    }
#else
  for(int num1 = 0; num1 < MRT_BINS; num1++)
    {
      SphP[i].DensPhot[num1] = 1.03e4 * pow(All.UnitLength_in_cm, 2) * All.UnitTime_in_s / All.UnitEnergy_in_cgs / c_internal_units;
      SphP[i].Cons_DensPhot[num1] = SphP[i].DensPhot[num1] * SphP[i].Volume;
      if(isnan(SphP[i].Cons_DensPhot[num1]))
        terminate("Something went wrong\n");

      for(int j = 0; j < 3; j++)
        {
          if(j == 1)
            SphP[i].RT_F[num1][j] = 0.9999999 * c_internal_units * SphP[i].DensPhot[num1];
          else
            SphP[i].RT_F[num1][j] = MINDENSPHOT;
        }
    }
#endif

#ifdef MRT_SETUP_SPECIAL_BOUNDARIES
  for(int num1 = 0; num1 < MRT_BINS; num1++)
    {
      if(RestartFlag == RESTART_IC)
        SphP[i].DensPhot[num1] = MINDENSPHOT;

      SphP[i].Cons_DensPhot[num1] = SphP[i].DensPhot[num1] * SphP[i].Volume;
      if(isnan(SphP[i].Cons_DensPhot[num1]))
        terminate("Something went wrong\n");

      for(int j = 0; j < 3; j++)
        {
          if(RestartFlag == RESTART_IC)
            SphP[i].RT_F[num1][j] = MINDENSPHOT;
          SphP[i].Cons_RT_F[num1][j] = SphP[i].RT_F[num1][j] * SphP[i].Volume;
        }
    }
#endif
  if(i == 0)
    mpi_printf("RT: INIT: Initialized conservative quantities\n");
}

void mpi_printf_rt(const int flag, const char *fmt, ...)
{
  if(ThisTask == 0)
    {
      va_list l;
      va_start(l, fmt);
#ifdef MRT_REDUCE_OUTPUT
      if(flag == 1)
        {
          vprintf(fmt, l);
          myflush(stdout);
        }
#else
      vprintf(fmt, l);
      myflush(stdout);
#endif
      va_end(l);
    }
}

#ifdef MRT_SINGLE_STAR

static double interpolate_stellar_params(int bin, double age)
{
  int imin, imax;

  double dx;

  if(age < starpar.Age[0])
    {
      imin = 0;
      imax = 0;
      dx   = 1.0;
    }
  else if(age > starpar.Age[Nlines - 1])
    {
      imin = Nlines - 1;
      imax = Nlines - 1;
      dx   = 0.0;
    }
  else
    {
      for(int i = 0; i < Nlines; i++)
        {
          if(starpar.Age[i] > age)
            {
              imax = i;
              imin = i - 1;
              dx   = (starpar.Age[imax] - age) / (starpar.Age[imax] - starpar.Age[imin]);
              break;
            }
        }
    }

  double Nphot = starpar.Nphot[imin] * dx + starpar.Nphot[imax] * (1.0 - dx);
  double frac  = starpar.frac[bin][imin] * dx + starpar.frac[bin][imax] * (1.0 - dx);

  mrt_sigma_H2[bin] = starpar.sigH2[bin][imin] * dx + starpar.sigH2[bin][imax] * (1.0 - dx);
  G_H2[bin]         = starpar.EH2[bin][imin] * dx + starpar.EH2[bin][imax] * (1.0 - dx);
  P_H2[bin]         = starpar.PH2[bin][imin] * dx + starpar.PH2[bin][imax] * (1.0 - dx);

  mrt_sigma_HI[bin] = starpar.sigH[bin][imin] * dx + starpar.sigH[bin][imax] * (1.0 - dx);
  G_HI[bin]         = starpar.EH[bin][imin] * dx + starpar.EH[bin][imax] * (1.0 - dx);
  P_HI[bin]         = starpar.PH[bin][imin] * dx + starpar.PH[bin][imax] * (1.0 - dx);

  mrt_sigma_HeI[bin] = starpar.sigHe[bin][imin] * dx + starpar.sigHe[bin][imax] * (1.0 - dx);
  G_HeI[bin]         = starpar.EHe[bin][imin] * dx + starpar.EHe[bin][imax] * (1.0 - dx);
  P_HeI[bin]         = starpar.PHe[bin][imin] * dx + starpar.PHe[bin][imax] * (1.0 - dx);

  return Nphot * frac;
}

void read_stellar_table(void)
{
  FILE *fp;
  fp = fopen(All.StellarParamFile, "r");

  double fac, fac_two;

#ifdef MRT_CHEM_SG
  fac = 1.0;
#else
  fac     = 1.0 / All.UnitLength_in_cm / All.UnitLength_in_cm * All.HubbleParam * All.HubbleParam;
#endif

#ifdef MRT_CHEM_SG
  fac_two = ELECTRONVOLT_IN_ERGS;
#else
  fac_two = ELECTRONVOLT_IN_ERGS / All.UnitEnergy_in_cgs * All.HubbleParam;
#endif

  double Tmax;

  fscanf(fp, "%d %le", &Nlines, &Tmax);

  All.LocalFeedbackSNTimeDelay = Tmax * 31556926 / All.UnitTime_in_s;

  starpar.frac = (double **)mymalloc("FractionOfPhotons", UV_BINS * sizeof(double *));

  starpar.sigH2 = (double **)mymalloc("Sigma_H2", UV_BINS * sizeof(double *));
  starpar.EH2   = (double **)mymalloc("Energy_H2", UV_BINS * sizeof(double *));
  starpar.PH2   = (double **)mymalloc("Momentum_H2", UV_BINS * sizeof(double *));

  starpar.sigH = (double **)mymalloc("Sigma_H", UV_BINS * sizeof(double *));
  starpar.EH   = (double **)mymalloc("Energy_H", UV_BINS * sizeof(double *));
  starpar.PH   = (double **)mymalloc("Momentum_H", UV_BINS * sizeof(double *));

  starpar.sigHe = (double **)mymalloc("Sigma_He", UV_BINS * sizeof(double *));
  starpar.EHe   = (double **)mymalloc("Energy_He", UV_BINS * sizeof(double *));
  starpar.PHe   = (double **)mymalloc("Momentum_He", UV_BINS * sizeof(double *));

  for(int i = 0; i < UV_BINS; i++)
    {
      starpar.Age     = (double *)mymalloc("AgeOfStars", Nlines * sizeof(double));
      starpar.Nphot   = (double *)mymalloc("NumberOfPhotons", Nlines * sizeof(double));
      starpar.frac[i] = (double *)mymalloc("FractionOfPhotons", Nlines * sizeof(double));

      starpar.sigH2[i] = (double *)mymalloc("Sigma_H2", Nlines * sizeof(double));
      starpar.EH2[i]   = (double *)mymalloc("Energy_H2", Nlines * sizeof(double));
      starpar.PH2[i]   = (double *)mymalloc("Momentum_H2", Nlines * sizeof(double));

      starpar.sigH[i] = (double *)mymalloc("Sigma_H", Nlines * sizeof(double));
      starpar.EH[i]   = (double *)mymalloc("Energy_H", Nlines * sizeof(double));
      starpar.PH[i]   = (double *)mymalloc("Momentum_H", Nlines * sizeof(double));

      starpar.sigHe[i] = (double *)mymalloc("Sigma_He", Nlines * sizeof(double));
      starpar.EHe[i]   = (double *)mymalloc("Energy_He", Nlines * sizeof(double));
      starpar.PHe[i]   = (double *)mymalloc("Momentum_He", Nlines * sizeof(double));
    }

  for(int num2 = 0; num2 < Nlines; num2++)
    {
      fscanf(fp, "%le %le", &starpar.Age[num2], &starpar.Nphot[num2]);
      starpar.Age[num2] *= 31556926 / All.UnitTime_in_s;
      for(int num1 = 0; num1 < UV_BINS; num1++)
        {
          fscanf(fp, "%le %le %le %le %le %le %le %le %le %le ", &starpar.frac[num1][num2], &starpar.sigH2[num1][num2],
                 &starpar.EH2[num1][num2], &starpar.PH2[num1][num2], &starpar.sigH[num1][num2], &starpar.EH[num1][num2],
                 &starpar.PH[num1][num2], &starpar.sigHe[num1][num2], &starpar.EHe[num1][num2], &starpar.PHe[num1][num2]);

          starpar.sigH2[num1][num2] *= fac;
          starpar.sigH[num1][num2] *= fac;
          starpar.sigHe[num1][num2] *= fac;

          starpar.EH2[num1][num2] *= fac_two;
          starpar.EH[num1][num2] *= fac_two;
          starpar.EHe[num1][num2] *= fac_two;

          starpar.PH2[num1][num2] *= fac_two;
          starpar.PH[num1][num2] *= fac_two;
          starpar.PHe[num1][num2] *= fac_two;
        }
    }

  fclose(fp);

  mpi_printf_rt(0, "RT: Read in the stellar parameter file \n");
}
#endif

#endif /* MRT */
