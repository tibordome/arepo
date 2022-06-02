/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/network/integrate.c
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

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <stdlib.h>

#include "../allvars.h"
#include "../proto.h"
#include "network_solver.h"
#include "string.h"

#include "integrate.h"

static const double conv = 1.602177e-12 * 1.0e3 * 6.0221367e23; /* eV2erg * 1.0e3 [keV] * avogadro */

static void network_do_single_particle(int i, double dt, double *dUtherm, double *composition);
static void network_fix_abundances(double *x, double *denergy, const struct network_data *nd);

static double nuclear_network_root(double dt_cell, void *params)
{
  int i       = *(int *)params;
  double sign = 1;
  if(i >= NumGas)
    {
      i -= NumGas;
      sign = -1;
    }

  double Composition[EOS_NSPECIES], dUtherm;
  network_do_single_particle(i, dt_cell, &dUtherm, Composition);

  struct eos_result res;
  double Temperature = SphP[i].EOSTemperature;
  eos_calc_egiven(SphP[i].Density, Composition, SphP[i].Utherm + dUtherm, &Temperature, &res);

  double dlnT = log(Temperature) - log(SphP[i].EOSTemperature);

  // printf( "cell=%d, sign=%g, dt_cell=%g, T=%g, Told=%g, dlnT=%g, res=%g\n", i, sign, dt_cell, Temperature, SphP[i].EOSTemperature,
  // dlnT, dlnT - sign * All.NuclearNetworkMaxTempChange );

  return dlnT - sign * All.NuclearNetworkMaxTempChange;
}

void network_main(int *timeBin)
{
  CPU_Step[CPU_MISC] += measure_time();

  mpi_printf("Doing nuclear network.\n");
  MPI_Barrier(MPI_COMM_WORLD);

  double tstart = second();
  int idx, i, *DoNetwork;

  DoNetwork = (int *)mymalloc("DoNetwork", sizeof(int) * TimeBinsHydro.NActiveParticles);
  memset(DoNetwork, 0, sizeof(int) * TimeBinsHydro.NActiveParticles);

  /* update DivVel and pressure gradient */
  exchange_primitive_variables();
  calculate_gradients();

  int networkCount = 0;
  int shockCount   = 0;

  /* select the particles for which we have to do nuclear burning */
#pragma omp parallel for private(idx, i) reduction(+ : networkCount, shockCount)
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      SphP[i].dedt = 0;

      if(SphP[i].EOSTemperature > All.NetworkTempThreshold && SphP[i].Density > 1e4)
        {
#ifdef NUCLEAR_NETWORK_DISABLE_BURNING_IN_SHOCK
#ifdef NUCLEAR_NETWORK_USE_SHOCKFINDER
          if(SphP[i].Machnumber < 1.)
#else
          double gradp = sqrt(SphP[i].Grad.dpressU[0] * SphP[i].Grad.dpressU[0] + SphP[i].Grad.dpressU[1] * SphP[i].Grad.dpressU[1] +
                              SphP[i].Grad.dpressU[2] * SphP[i].Grad.dpressU[2]);
          if(!(SphP[i].DivVel < 0 && gradp * get_cell_radius(i) / SphP[i].Pressure > 0.66))
#endif
#endif
            {
              DoNetwork[idx] = 1;
              networkCount++;

              double dUtherm;
              network_fix_abundances(SphP[i].Composition, &dUtherm, &All.nd);

              if(fabs(dUtherm) > 0.)
                {
                  dUtherm *= All.UnitEnergy_in_cgs;
                  SphP[i].Energy += dUtherm * P[i].Mass;
                  SphP[i].Utherm += dUtherm;

                  for(int k = 0; k < EOS_NSPECIES; k++)
                    SphP[i].MassComposition[k] = SphP[i].Composition[k] * P[i].Mass;

                  struct eos_result res;
                  eos_calc_egiven(SphP[i].Density, SphP[i].Composition, SphP[i].Utherm, &SphP[i].EOSTemperature, &res);
                }
            }
#ifdef NUCLEAR_NETWORK_DISABLE_BURNING_IN_SHOCK
          else
            {
              shockCount++;
            }
#endif
        }
    }

  int networkCountTot, NActiveParticlesTot, shockCountTot;
  MPI_Reduce(&networkCount, &networkCountTot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&TimeBinsHydro.NActiveParticles, &NActiveParticlesTot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&shockCount, &shockCountTot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  mpi_printf("NUCLEAR NETWORK: %d out of %d cells do nuclear network, %d disabled by shock detection.\n", networkCountTot,
             NActiveParticlesTot, shockCountTot);

  double tcost[NUM_THREADS];
  int tcount[NUM_THREADS], tmaxiter[NUM_THREADS], tsumiter[NUM_THREADS];

#ifdef NUCLEAR_NETWORK_TIMESTEP_LIMITER
  if(timeBin && *timeBin > 0)
    {
      int network_done    = 0;
      double *Composition = (double *)mymalloc("Composition", TimeBinsHydro.NActiveParticles * EOS_NSPECIES * sizeof(double));
      double *dUtherm     = (double *)mymalloc("dUtherm", TimeBinsHydro.NActiveParticles * sizeof(double));

      do
        {
          double dt = (*timeBin ? (((integertime)1) << *timeBin) : 0) * All.Timebase_interval;

          int stop = 0;
#pragma omp parallel for private(idx, i) schedule(dynamic)
          for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
            {
              i = TimeBinsHydro.ActiveParticleList[idx];
#pragma omp flush(stop)
              if(i < 0 || DoNetwork[idx] == 0 || stop)
                continue;

              network_do_single_particle(i, dt, &dUtherm[idx], &Composition[idx * EOS_NSPECIES]);

              int k;
              for(k = 0; k < EOS_NSPECIES; k++)
                {
                  if(fabs(dUtherm[idx]) / SphP[i].Utherm > All.NuclearNetworkMaxEnergyDiff)
                    stop = 1;
                }
            }

          if(!stop)
            {
#pragma omp parallel for private(idx, i) schedule(dynamic)
              for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
                {
                  i = TimeBinsHydro.ActiveParticleList[idx];
                  if(i < 0 || DoNetwork[idx] == 0)
                    continue;

                  int k;
                  for(k = 0; k < EOS_NSPECIES; k++)
                    {
                      SphP[i].Composition[k]     = Composition[idx * EOS_NSPECIES + k];
                      SphP[i].MassComposition[k] = SphP[i].Composition[k] * P[i].Mass;
                    }

                  SphP[i].Energy += dUtherm[idx] * P[i].Mass;
                  SphP[i].dedt = dUtherm[idx] / dt;
                }
              network_done = 1;
            }
          else
            {
              (*timeBin)--;
            }
        }
      while(!network_done);

      myfree(dUtherm);
      myfree(Composition);
    }
  else
#endif
    {
#ifdef NUCLEAR_NETWORK_ALT_MPI_PARALLELISATION
      if(NTask == 1)
#endif
        {
          for(int th = 0; th < NUM_THREADS; th++)
            {
              tcost[th]    = 0.;
              tcount[th]   = 0;
              tmaxiter[th] = 0;
              tsumiter[th] = 0;
            }

#pragma omp parallel for private(idx, i) schedule(dynamic)
          for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
            {
              i = TimeBinsHydro.ActiveParticleList[idx];
              if(i < 0 || DoNetwork[idx] == 0)
                continue;

              if(P[i].TimeBinHydro == 0)
                continue;

              int threadid        = get_thread_num();
              double tracer_start = second();

              double dt_cell = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

#ifdef NUCLEAR_NETWORK_LIMIT_COMPOSITION_CHANGE
              double Composition[EOS_NSPECIES], dUtherm;
              network_do_single_particle(i, dt_cell, &dUtherm, Composition);

              double dComposition[EOS_NSPECIES];
              double fac = 1.0;

              int k;
              for(k = 0; k < EOS_NSPECIES; k++)
                {
                  dComposition[k] = Composition[k] - SphP[i].Composition[k];
                  if(fabs(dComposition[k]) > All.NuclearNetworkMaxCompositionChange)
                    fac = fmin(fac, All.NuclearNetworkMaxCompositionChange / fabs(dComposition[k]));
                }

              for(k = 0; k < EOS_NSPECIES; k++)
                {
                  SphP[i].Composition[k] += fac * dComposition[k];
                  SphP[i].MassComposition[k] = SphP[i].Composition[k] * P[i].Mass;
                }

              SphP[i].Energy += fac * dUtherm * P[i].Mass;
              SphP[i].dedt = dUtherm / dt_cell;
#else
          double Composition[EOS_NSPECIES], dUtherm;
          network_do_single_particle(i, dt_cell, &dUtherm, Composition);

          struct eos_result res;
          double Temperature = SphP[i].EOSTemperature;
          eos_calc_egiven(SphP[i].Density, Composition, SphP[i].Utherm + dUtherm, &Temperature, &res);

          double dlnT = log(Temperature) - log(SphP[i].EOSTemperature);

          if(fabs(dlnT) > All.NuclearNetworkMaxTempChange)
            {
              gsl_function F;
              F.function = nuclear_network_root;

              int p    = dlnT < 0 ? NumGas + i : i;
              F.params = &p;

              const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
              gsl_root_fsolver *s            = gsl_root_fsolver_alloc(T);
              gsl_root_fsolver_set(s, &F, dt_cell * 1e-10, dt_cell);

              int status;
              int iter = 0;
              do
                {
                  iter++;
                  status      = gsl_root_fsolver_iterate(s);
                  dt_cell     = gsl_root_fsolver_root(s);
                  double x_lo = gsl_root_fsolver_x_lower(s);
                  double x_hi = gsl_root_fsolver_x_upper(s);
                  status      = gsl_root_test_interval(x_lo, x_hi, 0, 0.01);
                }
              while(status == GSL_CONTINUE && iter < 1000);

              if(iter == 1000)
                terminate("Network not converged for cell ID=%d", P[i].ID);

              tmaxiter[threadid] = imax(tmaxiter[threadid], iter);
              tsumiter[threadid] += iter;

              gsl_root_fsolver_free(s);
              network_do_single_particle(i, dt_cell, &dUtherm, Composition);
            }
          else
            {
              tmaxiter[threadid] = imax(tmaxiter[threadid], 1);
              tsumiter[threadid] += 1;
            }

          for(int k = 0; k < EOS_NSPECIES; k++)
            {
              SphP[i].Composition[k]     = Composition[k];
              SphP[i].MassComposition[k] = SphP[i].Composition[k] * P[i].Mass;
            }

          SphP[i].Energy += dUtherm * P[i].Mass;
          SphP[i].dedt = dUtherm / dt_cell;

          tcount[threadid]++;
          tcost[threadid] += timediff(tracer_start, second());
#endif
            }
        }
#ifdef NUCLEAR_NETWORK_ALT_MPI_PARALLELISATION
      else
        {
          network_do_cells_mpi(DoNetwork);
        }
#endif
    }

  myfree(DoNetwork);

  CPU_Step[CPU_NETWORK_INTEGRATION] += measure_time();
  update_primitive_variables();

  double maxcost = 0.;
  int maxiter    = 0;
  int totiter    = 0;
  int totcount   = 0;
  for(int th = 0; th < NUM_THREADS; th++)
    {
      if(tcost[th] > maxcost)
        maxcost = tcost[th];
      if(tmaxiter[th] > maxiter)
        maxiter = tmaxiter[th];
      totiter += tsumiter[th];
      totcount += tcount[th];
    }

  double losttime = 0.;
  for(int th = 0; th < NUM_THREADS; th++)
    {
      losttime += maxcost - tcost[th];
    }

  double tend = second();
  if(NTask == 1)
    mpi_printf("NUCLEAR NETWORK: done in %gs, slowest thread=%gs, wait fraction=%g, maxiter=%d, average iterations=%g.\n",
               timediff(tstart, tend), maxcost, losttime / (NUM_THREADS * maxcost), maxiter, (double)totiter / (double)totcount);
  MPI_Barrier(MPI_COMM_WORLD);

  CPU_Step[CPU_NETWORK_IMBALANCE] += measure_time();
}

void network_do_single_particle(int i, double dt, double *dUtherm, double *composition)
{
  double dedt = 0.;
  int k;

  int threadid = get_thread_num();

  int do_network = 1;
#ifdef NETWORK_NSE
  if(SphP[i].EOSTemperature >= All.NetworkNSEThreshold)
    {
      int k;
      double Composition[EOS_NSPECIES];
      for(k = 0; k < EOS_NSPECIES; k++)
        Composition[k] = SphP[i].Composition[k];
      double Temperature = SphP[i].EOSTemperature;

      network_nse_integrate_ye(SphP[i].Utherm, SphP[i].Density * All.UnitDensity_in_cgs, Composition, dt * All.UnitTime_in_s, &dedt,
                               &Temperature, &All.nd_nse, &All.nw_nse[threadid]);

      if(Temperature > 0.9 * All.NetworkNSEThreshold)
        {
          for(k = 0; k < EOS_NSPECIES; k++)
            composition[k] = Composition[k];
          do_network = 0;
        }
    }
#endif
  if(do_network)
    {
      double Composition[EOS_NSPECIES];
      for(k = 0; k < EOS_NSPECIES; k++)
        Composition[k] = SphP[i].Composition[k];

      network_integrate(SphP[i].EOSTemperature, SphP[i].Density * All.UnitDensity_in_cgs, Composition, dt * All.UnitTime_in_s, &dedt,
                        &All.nd, &All.nw[threadid]);

      for(k = 0; k < EOS_NSPECIES; k++)
        composition[k] = Composition[k];
    }

  *dUtherm = dedt * All.UnitEnergy_in_cgs / All.UnitTime_in_s * dt;
}

void network_normalize(double *x, double *e, const struct network_data *nd, struct network_workspace *nw)
{
  double sum, xnew;
  int i;

  sum = 0;
  for(i = 0; i < nd->nuc_count; i++)
    {
      sum += x[i];
    }

  if(e)
    {
      for(i = 0; i < nd->nuc_count; i++)
        {
          xnew = x[i] / sum;
          *e -= (xnew - x[i]) * nd->nucdata[i].exm * conv;
          x[i] = xnew;
        }
    }
  else
    {
      for(i = 0; i < nd->nuc_count; i++)
        x[i] /= sum;
    }
}

int network_integrate(double temp, double rho, double *x, double dt, double *dedt, const struct network_data *nd,
                      struct network_workspace *nw)
{
  double *y;
  double sum;
  int i;

  if(dt == 0 || temp < 1e7)
    {
      *dedt = 0;
      return 0;
    }

  /* calculate number densities */
  y = nw->y;
#ifdef NETWORK_SEP_YZ
  y[nd->iYz] = 0.0;
#endif
  for(i = 0; i < nd->nuc_count; i++)
    {
      y[i] = x[i] / nd->nucdata[i].na;
#ifdef NETWORK_SEP_YZ
      y[nd->iYz] += y[i] * gsl_pow_2(nd->nucdata[i].nz);
#endif
    }

#if NETWORK_VARIABLE
  y[nd->iTemp] = temp;
#endif
#if NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP
  y[nd->iRho] = rho;
#endif

  /* run network */
  network_solver_integrate(temp, rho, y, dt, nd, nw);

  /* normalise */
  sum = 0.0;
  for(i = 0; i < nd->nuc_count; i++)
    {
      if(y[i] > 1.0 / nd->nucdata[i].na)
        y[i] = 1.0 / nd->nucdata[i].na;
      if(y[i] < 1e-30)
        y[i] = 1e-30;
      sum += y[i] * nd->nucdata[i].na;
    }
  for(i = 0; i < nd->nuc_count; i++)
    {
      y[i] /= sum;
    }
  /* calculate change of mass fractions and energy release */
  *dedt = 0;
  for(i = 0; i < nd->nuc_count; i++)
    {
      double dx = (y[i] * nd->nucdata[i].na - x[i]) / dt;
      *dedt -= dx / nd->nucdata[i].na * nd->nucdata[i].exm;
      x[i] = y[i] * nd->nucdata[i].na;
    }
  *dedt *= conv;

  return 0;
}

void network_fix_abundances(double *x, double *denergy, const struct network_data *nd)
{
  double xold[nd->nuc_count];
  double sum = 0.0;
  for(int i = 0; i < nd->nuc_count; i++)
    {
      xold[i] = x[i];
      if(x[i] > 1.0)
        x[i] = 1.0;
      if(x[i] < 1e-30)
        x[i] = 1e-30;
      sum += x[i];
    }
  for(int i = 0; i < nd->nuc_count; i++)
    x[i] /= sum;

  /* calculate energy release */
  *denergy = 0;
  for(int i = 0; i < nd->nuc_count; i++)
    *denergy -= (xold[i] - x[i]) * nd->nucdata[i].exm;
  *denergy *= conv;
}

void network_composition_statistics(void)
{
  int i, k;

  double mass;
  double masscomp[EOS_NSPECIES];

  mass = 0;
  for(k = 0; k < EOS_NSPECIES; k++)
    masscomp[k] = 0;

  for(i = 0; i < NumGas; i++)
    if(P[i].Mass > 0 && P[i].ID > 0)
      {
        mass += P[i].Mass;
        for(k = 0; k < EOS_NSPECIES; k++)
          masscomp[k] += SphP[i].MassComposition[k];
      }

  double masstot;
  double masscomptot[EOS_NSPECIES];

  MPI_Reduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(masscomp, masscomptot, EOS_NSPECIES, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fprintf(FdNetwork, "%8g ", All.Time);
      for(k = 0; k < EOS_NSPECIES; k++)
        fprintf(FdNetwork, "%8g ", masscomptot[k] / masstot);
      fprintf(FdNetwork, "\n");
      myflush(FdNetwork);
    }
}

struct massrad
{
  double mass;
  double radius;
};

int compare_radii(const void *a, const void *b)
{
  if(((struct massrad *)a)->radius < ((struct massrad *)b)->radius)
    return -1;

  if(((struct massrad *)a)->radius > ((struct massrad *)b)->radius)
    return +1;

  return 0;
}

void detonate(void)
{
  int i;
  int count;
  /*
  int counts[NTask];
  MPI_Gather( &NumGas, 1, MPI_INT, counts, 1, MPI_INT, 0, MPI_COMM_WORLD );

  struct massrad *massrad;

  int totcount = 0;
  if(ThisTask == 0)
    {
      for(i = 0; i < NTask; i++)
        totcount += counts[i];
      massrad = mymalloc( "massrad", totcount * sizeof(struct massrad) );
    }
  else
    {
      massrad = mymalloc( "massrad", NumGas * sizeof(struct massrad) );
    }

  if(ThisTask == 0)
    {
      for(i = 0; i < NumGas; i++)
        {
          massrad[i].mass = P[i].Mass;

          double dx = SphP[i].Center[0] - 0.5*All.BoxSize;
          double dy = SphP[i].Center[1] - 0.5*All.BoxSize;
          double dz = SphP[i].Center[2] - 0.5*All.BoxSize;

          massrad[i].radius = sqrt( dx*dx + dy*dy + dz*dz );
        }

      int count = NumGas;

      int task;
      MPI_Status status;
      for(task = 1; task < NTask; task++)
        {
          MPI_Recv( &massrad[count], counts[task] * sizeof(struct massrad), MPI_BYTE, task, 666, MPI_COMM_WORLD, &status );
          count += counts[task];
        }
    }
  else
    {
      for(i = 0; i < NumGas; i++)
        {
          massrad[i].mass = P[i].Mass;

          double dx = SphP[i].Center[0] - 0.5*All.BoxSize;
          double dy = SphP[i].Center[1] - 0.5*All.BoxSize;
          double dz = SphP[i].Center[2] - 0.5*All.BoxSize;

          massrad[i].radius = sqrt( dx*dx + dy*dy + dz*dz );
        }

      MPI_Send( massrad, NumGas * sizeof(struct massrad), MPI_BYTE, 0, 666, MPI_COMM_WORLD );
    }

  double detradius = -1.;
  if(ThisTask == 0)
    {
      qsort( massrad, totcount, sizeof(struct massrad), compare_radii );

      double msum = 0;
      for(i = 0; i < totcount; i++)
        {
          msum += massrad[i].mass;

          if(msum/SOLAR_MASS >= 0.005)
            {
              detradius = massrad[i].radius;
              break;
            }
        }
    }

  myfree( massrad );

  MPI_Bcast( &detradius, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

  for(i = 0; i < NumGas; i++)
    {
      double dx = SphP[i].Center[0] - 0.5*All.BoxSize;
      double dy = SphP[i].Center[1] - 0.5*All.BoxSize;
      double dz = SphP[i].Center[2] - 0.5*All.BoxSize;

      double rad = sqrt( dx*dx + dy*dy + dz*dz );

      if(rad < detradius)
        {
          SphP[i].Utherm += 5e16;
          SphP[i].Energy = P[i].Mass * SphP[i].Utherm + 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] +
                                                                           P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
          count++;
        }
    }
  */

#ifdef NUCLEAR_NETWORK_DETONATE_POSITION
  // get radius for detonation from file
#define LENGTH 30  // length of a row in Core.txt
#define NUMBER 2   // number of rows in Core.txt
  FILE *fp;
  char line[NUMBER][LENGTH];
  char *Mass;
  int d;
  fp = fopen("Core_Radius.txt", "r");
  if(fp == NULL)
    {
      printf("File can not be opened.\n");
      if(!fp)
        perror("fopen");
    }

  for(d = 0; d < NUMBER; d++)
    {  // reads file line for line
      fgets(line[d], LENGTH, fp);
      //    printf("%i) %s\n", d, line[d]);
    }

  Mass           = line[1];  // value of CoreMass is put into variable
  float CoreMass = atof(Mass);
  fclose(fp);

#ifdef ONEDIMS
  count = 0;
  for(i = 0; i < NumGas; i++)
    {
      msum += P[i].Mass;
      float MassMax, MassMin;
      MassMax = CoreMass / SOLAR_MASS + 0.002;
      MassMin = CoreMass / SOLAR_MASS - 0.002;
      if((msum / SOLAR_MASS < MassMax) && (msum / SOLAR_MASS > MassMin))
        {
          SphP[i].Utherm += 5e16;
          SphP[i].Energy = P[i].Mass * SphP[i].Utherm +
                           0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
          count++;
        }
    }

  int countall = 0;
  MPI_Reduce(&count, &countall, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  mpi_printf("Ignited detonation in %d cells.\n", countall);
#endif

#ifndef ONEDIMS
  detonate_shell_point(CoreMass);
  mpi_printf("Core mass: %g\n", CoreMass);
  mpi_printf("\nShell detonation\n");
#endif

#endif

#ifndef NUCLEAR_NETWORK_DETONATE_POSITION
  count = 0;
  for(i = 0; i < NumGas; i++)
    {
      double rad    = P[i].Pos[0];
      double maxrad = 4e7;
      if(rad < maxrad)
        {
          SphP[i].EOSTemperature = 2e9 - rad / maxrad * 0.8e9;
          struct eos_result res;
          eos_calc_tgiven(SphP[i].Density, SphP[i].Composition, SphP[i].EOSTemperature, &res);
          SphP[i].Utherm = res.e.v;

          // P[i].Vel[0] += P[i].Pos[0] / rad * 5e8 * rad / maxrad;
          // SphP[i].Momentum[0] = P[i].Mass * P[i].Vel[0];

          SphP[i].Energy = P[i].Mass * SphP[i].Utherm +
                           0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
          count++;
        }
    }

  int countall = 0;
  MPI_Reduce(&count, &countall, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  mpi_printf("Ignited detonation in %d cells.\n", countall);
#endif
}

#ifdef NUCLEAR_NETWORK_DETONATE_CORE
void detonate_core(void)
{
  double massdet = 0;
  int j          = 0;
  double Number  = TimeBinsGravity.GlobalNActiveParticles;
  mpi_printf("\nNumber %g\n", Number);
  int ii = 0;

  for(ii = 0; ii < Number; ii++)
    {
      //    printf("\ni %d", i);
      //    mpi_printf("\nchecking values");
      if(SphP[ii].Density > 6.5e6 && SphP[ii].Temperature > 1.10e8 && kk == 0 && P[ii].Pos[0] < 5.2e7)
        {
          for(j = 0; j < ii; j++)
            {
              massdet += P[j].Mass;
            }
          int count = 0;
          int counts[NTask];
          double msum = 0;
          int i;

          MPI_Gather(&NumGas, 1, MPI_INT, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);
          double CoreMass = DetMass / SOLAR_MASS;  // /SOLAR_MASS;
          double MassMax, MassMin;
          MassMax = CoreMass + 0.0012;  // asymmetric because of transition area at masses < Coremass
          MassMin = CoreMass - 0.0012;
          if(MassMin < 0)
            {
              MassMin = 0;
            }
          //          mpi_printf( "DetMass %e\n" , DetMass );
          //          mpi_printf( "CoreMass %e\n" , CoreMass );
          //          mpi_printf( "MassMax %e\n" , MassMax );
          //          mpi_printf( "MassMin %e\n" , MassMin );

          for(i = 0; i < NumGas; i++)
            {
              msum += P[i].Mass;
              if((msum / SOLAR_MASS >= MassMin) && (msum / SOLAR_MASS < MassMax))
                {
                  SphP[i].Utherm += 5e16;
                  SphP[i].Energy =
                      P[i].Mass * SphP[i].Utherm +
                      0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
                  count++;
                }
            }

          MPI_Bcast(&CoreMass, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

          int countall = 0;
          MPI_Reduce(&count, &countall, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
          mpi_printf("Ignited detonation in %d cells.\n", countall);
          produce_dump();
          mpi_printf("\nCore detonation\n");
          kk++;
        }
    }
}
#endif
#ifdef NUCLEAR_NETWORK_DETONATE_POSITION
void detonate_shell_point(double CoreMass)
{
  int counts[NTask];
  int i;

  MPI_Gather(&NumGas, 1, MPI_INT, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

  struct massrad *massrad;

  int totcount = 0;
  if(ThisTask == 0)
    {
      for(i = 0; i < NTask; i++)
        totcount += counts[i];
      massrad = mymalloc("massrad", totcount * sizeof(struct massrad));
    }
  else
    {
      massrad = mymalloc("massrad", NumGas * sizeof(struct massrad));
    }

  if(ThisTask == 0)
    {
      for(i = 0; i < NumGas; i++)
        {
          massrad[i].mass = P[i].Mass;

          double dx = SphP[i].Center[0] - 0.5 * All.BoxSize;
          double dy = SphP[i].Center[1] - 0.5 * All.BoxSize;
          double dz = SphP[i].Center[2] - 0.5 * All.BoxSize;

          massrad[i].radius = sqrt(dx * dx + dy * dy + dz * dz);
        }

      int count = NumGas;

      int task;
      MPI_Status status;
      for(task = 1; task < NTask; task++)
        {
          MPI_Recv(&massrad[count], counts[task] * sizeof(struct massrad), MPI_BYTE, task, 666, MPI_COMM_WORLD, &status);
          count += counts[task];
        }
    }
  else
    {
      for(i = 0; i < NumGas; i++)
        {
          massrad[i].mass = P[i].Mass;

          double dx = SphP[i].Center[0] - 0.5 * All.BoxSize;
          double dy = SphP[i].Center[1] - 0.5 * All.BoxSize;
          double dz = SphP[i].Center[2] - 0.5 * All.BoxSize;

          massrad[i].radius = sqrt(dx * dx + dy * dy + dz * dz);
        }

      MPI_Send(massrad, NumGas * sizeof(struct massrad), MPI_BYTE, 0, 666, MPI_COMM_WORLD);
    }

  double detradius = -1.;
  if(ThisTask == 0)
    {
      qsort(massrad, totcount, sizeof(struct massrad), compare_radii);

      double msum = 0;
      for(i = 0; i < totcount; i++)
        {
          msum += massrad[i].mass;

          //          mpi_printf("Core mass2: %g\n", CoreMass);
          float MassMax, MassMin;
          MassMax = CoreMass / SOLAR_MASS + 0.005;
          MassMin = CoreMass / SOLAR_MASS - 0.005;
          if((msum / SOLAR_MASS < MassMax) && (msum / SOLAR_MASS > MassMin))
            {
              detradius = massrad[i].radius;
              break;
            }
        }
    }
  myfree(massrad);

  MPI_Bcast(&detradius, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  int count = 0;
  for(i = 0; i < NumGas; i++)
    {
      double dx = SphP[i].Center[0] - 0.5 * All.BoxSize;
      double dy = SphP[i].Center[1] - 0.5 * All.BoxSize;
      double dz = SphP[i].Center[2] - 0.5 * All.BoxSize;

      double rad    = sqrt(dx * dx + dy * dy + dz * dz);
      double maxrad = detradius + 0.04 * detradius;
      double minrad = detradius - 0.04 * detradius;
      double maxdx  = 3e7;
      double mindx  = -3e7;
      double maxdy  = 3e7;
      double mindy  = -3e7;
      //      mpi_printf("detradius %g", detradius);

      if((rad <= maxrad) && (rad >= minrad) && (dx <= maxdx) && (dx >= mindx) && (dy <= maxdy) && (dy >= mindy) && (dz > 0))
        {
          SphP[i].Utherm += 5e16;
          SphP[i].Energy = P[i].Mass * SphP[i].Utherm +
                           0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
          count++;
        }
      //      mpi_printf( "Detonation radius %g, max %g, min %g\n", detradius, maxrad, minrad );
    }

  mpi_printf("Detonation radius: %g\n", detradius);
  int countall = 0;
  MPI_Reduce(&count, &countall, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  mpi_printf("Ignited detonation in %d cells in shell.\n", countall);
}
#endif
