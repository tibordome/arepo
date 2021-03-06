/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/global.c
 * \date        MM/YYYY
 * \author
 * \brief       Routines to compute statistics of the global state of the system
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

#ifdef COSMIC_RAYS
static void cr_energy_statistics(void);
#endif

#ifdef GW_SIGNAL
static void compute_gw_signal(void);
#endif

/*! \brief Computes new global statistics if needed (call of
 *         energy_statistics()).
 *
 *  \return void
 */
void compute_statistics(void)
{
  /* check whether we want a full energy statistics */
  if((All.Time - All.TimeLastStatistics) >= All.TimeBetStatistics &&
     All.HighestActiveTimeBin == All.HighestOccupiedTimeBin) /* allow only top-level synchronization points */
    {
      TIMER_START(CPU_LOGS);

      energy_statistics(); /* compute and output energy statistics */
#ifdef BINARYLOG
      binary_statistics();
#endif

#ifdef NUCLEAR_NETWORK
      network_composition_statistics();
#endif

#ifdef GENERAL_RELATIVITY
      general_relativity_statistics();
#endif

#ifdef COSMIC_RAYS
      cr_energy_statistics();
#endif

#ifdef GW_SIGNAL
      compute_gw_signal();
#endif

      All.TimeLastStatistics += All.TimeBetStatistics;

      TIMER_STOP(CPU_LOGS);
    }
}

/*! \brief Compute global statistics of the system.
 *
 *  This function first calls a computation of various global
 *  quantities of the particle distribution
 *  (compute_global_quantities_of_system() ), and then writes some statistics
 *  about the energies of the various particle types to the file FdEnergy
 *  (energy.txt).
 *
 *  \return void
 */
void energy_statistics(void)
{
  double egyinj_tot;

  compute_global_quantities_of_system();

  MPI_Reduce(&EgyInjection, &egyinj_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fprintf(FdEnergy, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", All.Time,
              SysState.EnergyInt, SysState.EnergyPot, SysState.EnergyKin, SysState.EnergyIntComp[0], SysState.EnergyPotComp[0],
              SysState.EnergyKinComp[0], SysState.EnergyIntComp[1], SysState.EnergyPotComp[1], SysState.EnergyKinComp[1],
              SysState.EnergyIntComp[2], SysState.EnergyPotComp[2], SysState.EnergyKinComp[2], SysState.EnergyIntComp[3],
              SysState.EnergyPotComp[3], SysState.EnergyKinComp[3], SysState.EnergyIntComp[4], SysState.EnergyPotComp[4],
              SysState.EnergyKinComp[4], SysState.EnergyIntComp[5], SysState.EnergyPotComp[5], SysState.EnergyKinComp[5],
              SysState.MassComp[0], SysState.MassComp[1], SysState.MassComp[2], SysState.MassComp[3], SysState.MassComp[4],
              SysState.MassComp[5], egyinj_tot);

      myflush(FdEnergy);
    }
}

#ifdef COSMIC_RAYS
void cr_energy_statistics(void)
{
  double localEnergies[10], globalEnergies[10];

  for(int k = 0; k < 10; k++)
    localEnergies[k] = 0;

  for(int i = 0; i < NumGas; i++)
    localEnergies[0] += SphP[i].CR_Energy;

  localEnergies[1] = All.TotalCREnergyInjected;
  localEnergies[2] = All.TotalCREnergyCooledHadronic;
  localEnergies[3] = All.TotalCREnergyCooledCoulomb;
  localEnergies[4] = All.TotalCREnergyCooled;

#ifdef COSMIC_RAYS_EXTRA_DIAGNOSTICS
  localEnergies[5] = All.TotalCREnergyChangeAdiabatic;
  localEnergies[6] = All.TotalCREnergyErrorDiffusion;
  localEnergies[7] = All.TotalCREnergyLossSfr;
  localEnergies[8] = All.TotalCREnergyUpdatePrims;
  localEnergies[9] = All.TotalCREnergyErrorStreaming;
#endif

  MPI_Reduce(&localEnergies[0], &globalEnergies[0], 9, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
#ifdef COSMIC_RAYS_EXTRA_DIAGNOSTICS
      fprintf(FdCREnergy, "%g %g %g %g %g %g %g %g %g %g %g\n", All.Time, globalEnergies[0], globalEnergies[1], globalEnergies[2],
              globalEnergies[3], globalEnergies[4], globalEnergies[5], globalEnergies[6], globalEnergies[7], globalEnergies[8],
              globalEnergies[9]);
#else
      fprintf(FdCREnergy, "%g %g %g %g %g %g\n", All.Time, globalEnergies[0], globalEnergies[1], globalEnergies[2], globalEnergies[3],
              globalEnergies[4]);
#endif
      myflush(FdCREnergy);
    }
}
#endif

#ifdef GW_SIGNAL
void compute_gw_signal(void)
{
  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
    {
      double Mass = 0;
      double CenterOfMass[3];
      for(int i = 0; i < 3; i++)
        CenterOfMass[i] = 0;

      for(int i = 0; i < NumPart; i++)
        {
          if(P[i].Mass == 0 && P[i].ID == 0)
            continue;

#ifdef TRACER_PARTICLE
          if(P[i].Type == TRACER_PARTICLE)
            continue;
#endif

          double *Center;
          if(P[i].Type == 0)
            Center = SphP[i].Center;
          else
            Center = P[i].Pos;

          Mass += P[i].Mass;
          for(int k = 0; k < 3; k++)
            CenterOfMass[k] += P[i].Mass * Center[k];
        }

      double MassAll, CenterOfMassAll[3];
      MPI_Allreduce(&Mass, &MassAll, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(CenterOfMass, CenterOfMassAll, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      for(int i = 0; i < 3; i++)
        CenterOfMassAll[i] /= MassAll;

      double A[6];
      for(int i = 0; i < 6; i++)
        A[i] = 0;

      for(int i = 0; i < NumPart; i++)
        {
          if(P[i].Mass == 0 && P[i].ID == 0)
            continue;

#ifdef TRACER_PARTICLE
          if(P[i].Type == TRACER_PARTICLE)
            continue;
#endif
          double acc[3];
          for(int k = 0; k < 3; k++)
            acc[k] = P[i].GravAccel[k];
#ifdef PMGRID
          for(int k = 0; k < 3; k++)
            acc[k] += P[i].GravPM[k];
#endif

          double *Center;
          if(P[i].Type == 0)
            Center = SphP[i].Center;
          else
            Center = P[i].Pos;

          A[0] += P[i].Mass * (2. * P[i].Vel[0] * P[i].Vel[0] + (Center[0] - CenterOfMassAll[0]) * acc[0] +
                               (Center[0] - CenterOfMassAll[0]) * acc[0]);
          A[1] += P[i].Mass * (2. * P[i].Vel[1] * P[i].Vel[1] + (Center[1] - CenterOfMassAll[1]) * acc[1] +
                               (Center[1] - CenterOfMassAll[1]) * acc[1]);
          A[2] += P[i].Mass * (2. * P[i].Vel[2] * P[i].Vel[2] + (Center[2] - CenterOfMassAll[2]) * acc[2] +
                               (Center[2] - CenterOfMassAll[2]) * acc[2]);
          A[3] += P[i].Mass * (2. * P[i].Vel[0] * P[i].Vel[1] + (Center[0] - CenterOfMassAll[0]) * acc[1] +
                               (Center[1] - CenterOfMassAll[1]) * acc[0]);
          A[4] += P[i].Mass * (2. * P[i].Vel[0] * P[i].Vel[2] + (Center[0] - CenterOfMassAll[0]) * acc[2] +
                               (Center[2] - CenterOfMassAll[2]) * acc[0]);
          A[5] += P[i].Mass * (2. * P[i].Vel[1] * P[i].Vel[2] + (Center[1] - CenterOfMassAll[1]) * acc[2] +
                               (Center[2] - CenterOfMassAll[2]) * acc[1]);
        }

      double AAll[6];
      MPI_Reduce(A, AAll, 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

      if(ThisTask == 0)
        {
          for(int i = 0; i < 6; i++)
            AAll[i] *= GRAVITY / pow(CLIGHT, 4) * All.UnitMass_in_g * pow(All.UnitVelocity_in_cm_per_s, 2);

          fprintf(FdGW, "%24.14e %24.14e %24.14e %24.14e %24.14e %24.14e %24.14e %24.14e %24.14e %24.14e %24.14e\n", All.Time, MassAll,
                  CenterOfMassAll[0], CenterOfMassAll[1], CenterOfMassAll[2], AAll[0], AAll[1], AAll[2], AAll[3], AAll[4], AAll[5]);
          fflush(FdGW);
        }
    }
}
#endif

/*! \brief This routine computes various global properties of the particle
 *         distribution and stores the result in the struct `SysState'.
 *
 *  Currently, not all the information that's computed here is
 *  actually used (e.g. momentum is not really used anywhere),
 *  just the energies are written to a log-file every once in a while.
 *
 *  \return void
 */
void compute_global_quantities_of_system(void)
{
  int i, j, n;
  struct state_of_system sys;
  double egyspec, vel[3];

  for(n = 0; n < NTYPES; n++)
    {
      sys.MassComp[n] = sys.EnergyKinComp[n] = sys.EnergyPotComp[n] = sys.EnergyIntComp[n] = 0;

      for(j = 0; j < 4; j++)
        sys.CenterOfMassComp[n][j] = sys.MomentumComp[n][j] = sys.AngMomentumComp[n][j] = 0;
    }

  for(i = 0; i < NumPart; i++)
    {
      sys.MassComp[P[i].Type] += P[i].Mass;

#if defined(SELFGRAVITY)
#ifdef EVALPOTENTIAL
#ifndef EXACT_GRAVITY_FOR_PARTICLE_TYPE
      sys.EnergyPotComp[P[i].Type] +=
          0.5 * P[i].Mass * (P[i].Potential + All.G * P[i].Mass / (All.ForceSoftening[P[i].SofteningType] / 2.8)) / All.cf_atime;
#else
      /* ignore self-contribution from gravity if exact gravity is used */
      if(P[i].Type == EXACT_GRAVITY_FOR_PARTICLE_TYPE)
        sys.EnergyPotComp[P[i].Type] += 0.5 * P[i].Mass * P[i].Potential / All.cf_atime;
      else
        sys.EnergyPotComp[P[i].Type] +=
            0.5 * P[i].Mass * (P[i].Potential + All.G * P[i].Mass / (All.ForceSoftening[P[i].SofteningType] / 2.8)) / All.cf_atime;
#endif
#endif
#endif

#if defined(EXTERNALGRAVITY)
#if defined(SELFGRAVITY)
      sys.EnergyPotComp[P[i].Type] += 0.5 * P[i].Mass * P[i].ExtPotential; /* note: ExtPotential already included on P[].p.Potential,
                                                                              that's why only 0.5 is needed here to recover the rest */
#else
      sys.EnergyPotComp[P[i].Type] += 1.0 * P[i].Mass * P[i].ExtPotential;
#endif
#endif

      if(P[i].Type == 0)
        {
          for(j = 0; j < 3; j++)
            {
              vel[j] = P[i].Vel[j];
            }

          sys.EnergyKinComp[0] += 0.5 * P[i].Mass * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);

          egyspec = SphP[i].Utherm;

          sys.EnergyIntComp[0] += P[i].Mass * egyspec;
        }
      else
        {
          for(j = 0; j < 3; j++)
            {
              vel[j] = P[i].Vel[j];
            }
          sys.EnergyKinComp[P[i].Type] += 0.5 * P[i].Mass * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]) * All.cf_a2inv;
        }

      for(j = 0; j < 3; j++)
        {
          sys.MomentumComp[P[i].Type][j] += P[i].Mass * vel[j];
          sys.CenterOfMassComp[P[i].Type][j] += P[i].Mass * P[i].Pos[j];
        }

      sys.AngMomentumComp[P[i].Type][0] += P[i].Mass * (P[i].Pos[1] * vel[2] - P[i].Pos[2] * vel[1]);
      sys.AngMomentumComp[P[i].Type][1] += P[i].Mass * (P[i].Pos[2] * vel[0] - P[i].Pos[0] * vel[2]);
      sys.AngMomentumComp[P[i].Type][2] += P[i].Mass * (P[i].Pos[0] * vel[1] - P[i].Pos[1] * vel[0]);
    }

  /* some the stuff over all processors */
  MPI_Reduce(&sys.MassComp[0], &SysState.MassComp[0], NTYPES, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sys.EnergyPotComp[0], &SysState.EnergyPotComp[0], NTYPES, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sys.EnergyIntComp[0], &SysState.EnergyIntComp[0], NTYPES, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sys.EnergyKinComp[0], &SysState.EnergyKinComp[0], NTYPES, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sys.MomentumComp[0][0], &SysState.MomentumComp[0][0], NTYPES * 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sys.AngMomentumComp[0][0], &SysState.AngMomentumComp[0][0], NTYPES * 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sys.CenterOfMassComp[0][0], &SysState.CenterOfMassComp[0][0], NTYPES * 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      for(i = 0; i < NTYPES; i++)
        SysState.EnergyTotComp[i] = SysState.EnergyKinComp[i] + SysState.EnergyPotComp[i] + SysState.EnergyIntComp[i];

      SysState.Mass = SysState.EnergyKin = SysState.EnergyPot = SysState.EnergyInt = SysState.EnergyTot = 0;

      for(j = 0; j < 3; j++)
        SysState.Momentum[j] = SysState.AngMomentum[j] = SysState.CenterOfMass[j] = 0;

      for(i = 0; i < NTYPES; i++)
        {
          SysState.Mass += SysState.MassComp[i];
          SysState.EnergyKin += SysState.EnergyKinComp[i];
          SysState.EnergyPot += SysState.EnergyPotComp[i];
          SysState.EnergyInt += SysState.EnergyIntComp[i];
          SysState.EnergyTot += SysState.EnergyTotComp[i];

          for(j = 0; j < 3; j++)
            {
              SysState.Momentum[j] += SysState.MomentumComp[i][j];
              SysState.AngMomentum[j] += SysState.AngMomentumComp[i][j];
              SysState.CenterOfMass[j] += SysState.CenterOfMassComp[i][j];
            }
        }

      for(i = 0; i < NTYPES; i++)
        for(j = 0; j < 3; j++)
          if(SysState.MassComp[i] > 0)
            SysState.CenterOfMassComp[i][j] /= SysState.MassComp[i];

      for(j = 0; j < 3; j++)
        if(SysState.Mass > 0)
          SysState.CenterOfMass[j] /= SysState.Mass;

      for(i = 0; i < NTYPES; i++)
        {
          SysState.CenterOfMassComp[i][3] = SysState.MomentumComp[i][3] = SysState.AngMomentumComp[i][3] = 0;
          for(j = 0; j < 3; j++)
            {
              SysState.CenterOfMassComp[i][3] += SysState.CenterOfMassComp[i][j] * SysState.CenterOfMassComp[i][j];
              SysState.MomentumComp[i][3] += SysState.MomentumComp[i][j] * SysState.MomentumComp[i][j];
              SysState.AngMomentumComp[i][3] += SysState.AngMomentumComp[i][j] * SysState.AngMomentumComp[i][j];
            }
          SysState.CenterOfMassComp[i][3] = sqrt(SysState.CenterOfMassComp[i][3]);
          SysState.MomentumComp[i][3]     = sqrt(SysState.MomentumComp[i][3]);
          SysState.AngMomentumComp[i][3]  = sqrt(SysState.AngMomentumComp[i][3]);
        }

      SysState.CenterOfMass[3] = SysState.Momentum[3] = SysState.AngMomentum[3] = 0;

      for(j = 0; j < 3; j++)
        {
          SysState.CenterOfMass[3] += SysState.CenterOfMass[j] * SysState.CenterOfMass[j];
          SysState.Momentum[3] += SysState.Momentum[j] * SysState.Momentum[j];
          SysState.AngMomentum[3] += SysState.AngMomentum[j] * SysState.AngMomentum[j];
        }

      SysState.CenterOfMass[3] = sqrt(SysState.CenterOfMass[3]);
      SysState.Momentum[3]     = sqrt(SysState.Momentum[3]);
      SysState.AngMomentum[3]  = sqrt(SysState.AngMomentum[3]);
    }

  /* give everyone the result, maybe the want to do something with it */
  MPI_Bcast(&SysState, sizeof(struct state_of_system), MPI_BYTE, 0, MPI_COMM_WORLD);
}

/*! \brief Compute binary statistics of the system.
 *
 * This function computes positions and velocities for binary systems
 */
#ifdef BINARYLOG_FOR_MERGERS
struct star_data
{
  double center[3];
  double vel[3];
  double mass;
};

static void get_star_properties(struct star_data *star, int ipass)
{
  double mass = 0;
  double center[3], mom[3];
  for(int k = 0; k < 3; k++)
    {
      center[k] = 0;
      mom[k]    = 0;
    }

  for(int i = 0; i < NumGas; i++)
    {
      if(P[i].Type != 0 || (P[i].Mass == 0 && P[i].ID == 0))
        continue;

      double cmass = P[i].Mass * SphP[i].PScalars[ipass];
      mass += cmass;
      for(int k = 0; k < 3; k++)
        {
          center[k] += SphP[i].Center[k] * cmass;
          mom[k] += P[i].Vel[k] * cmass;
        }
    }
  double mtot, centertot[3], momtot[3];
  MPI_Allreduce(&mass, &mtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(center, centertot, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(mom, momtot, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  star->mass = mtot;
  if(star->mass > 0)
    {
      for(int k = 0; k < 3; k++)
        {
          star->center[k] = centertot[k] / star->mass;
          star->vel[k]    = momtot[k] / star->mass;
        }
    }
  else
    {
      for(int k = 0; k < 3; k++)
        {
          star->center[k] = 0;
          star->vel[k]    = 0;
        }
    }
}
#endif
#ifdef BINARYLOG
void binary_statistics(void)
{
  int i, j;
  double pi = 3.141592653589;

  /* get position and velocity of RG core and companion */
  double pos_rgcore[3], pos_companion[3];
  double vel_rgcore[3], vel_companion[3];
  double m_rgcore, m_companion;
  for(j = 0; j < 3; j++)
    pos_rgcore[j] = pos_companion[j] = vel_rgcore[j] = vel_companion[j] = 0.0;
  m_rgcore    = 0;
  m_companion = 0;
#ifndef BINARYLOG_FOR_MERGERS
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0)
        continue;
      if(P[i].ID == ID_RGCORE)
        {
          for(j = 0; j < 3; j++)
            {
              pos_rgcore[j] = P[i].Pos[j];
              vel_rgcore[j] = P[i].Vel[j];
            }
          m_rgcore = P[i].Mass;
        }
      if(P[i].ID == ID_RGCORE + 1)
        {
          for(j = 0; j < 3; j++)
            {
              pos_companion[j] = P[i].Pos[j];
              vel_companion[j] = P[i].Vel[j];
            }
          m_companion = P[i].Mass;
        }
    }

  /* communicate everything */
  double buf[3];
  MPI_Reduce(pos_rgcore, buf, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  for(j = 0; j < 3; j++)
    pos_rgcore[j] = buf[j];
  MPI_Reduce(vel_rgcore, buf, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  for(j = 0; j < 3; j++)
    vel_rgcore[j] = buf[j];
  MPI_Reduce(pos_companion, buf, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  for(j = 0; j < 3; j++)
    pos_companion[j] = buf[j];
  MPI_Reduce(vel_companion, buf, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  for(j = 0; j < 3; j++)
    vel_companion[j] = buf[j];
  MPI_Reduce(&m_rgcore, buf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  m_rgcore = buf[0];
  MPI_Reduce(&m_companion, buf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  m_companion = buf[0];

#else
  struct star_data star1, star2;
  get_star_properties(&star1, 0);  // Modify passive scalar to the ones that fit your stars
  get_star_properties(&star2, 1);
  for(j = 0; j < 3; j++)
    {
      pos_rgcore[j] = star1.center[j];
      vel_rgcore[j] = star1.vel[j];
      pos_companion[j] = star2.center[j];
      vel_companion[j] = star2.vel[j];
    }
  m_rgcore = star1.mass;
  m_companion = star2.mass;
#endif

  if(ThisTask == 0)
    {
      double distance;
      distance = sqrt((pos_rgcore[0] - pos_companion[0]) * (pos_rgcore[0] - pos_companion[0]) +
                      (pos_rgcore[1] - pos_companion[1]) * (pos_rgcore[1] - pos_companion[1]) +
                      (pos_rgcore[2] - pos_companion[2]) * (pos_rgcore[2] - pos_companion[2]));

      double period;
      period = sqrt(GRAVITY * (m_rgcore + m_companion) / (distance * distance * distance));
      period = 2.0 * pi / period;

      fprintf(FdBinary,
              "%24.14e %24.14e %24.14e %24.14e %24.14e %24.14e %24.14e %24.14e %24.14e %24.14e %24.14e %24.14e %24.14e %24.14e "
              "%24.14e %24.14e %24.14e %24.14e %24.14e %24.14e %24.14e %24.14e %24.14e %24.14e %24.14e %24.14e\n",
              All.Time, pos_rgcore[0], pos_rgcore[1], pos_rgcore[2], pos_companion[0], pos_companion[1], pos_companion[2],
              vel_rgcore[0], vel_rgcore[1], vel_rgcore[2], vel_companion[0], vel_companion[1], vel_companion[2], m_rgcore, m_companion,
              distance, period, SysState.CenterOfMass[0], SysState.CenterOfMass[1], SysState.CenterOfMass[2], SysState.Momentum[0],
              SysState.Momentum[1], SysState.Momentum[2], SysState.AngMomentum[0], SysState.AngMomentum[1], SysState.AngMomentum[2]);

      myflush(FdBinary);
    }
}
#endif
