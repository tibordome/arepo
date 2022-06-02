/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/io_fields.c
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

#include <errno.h>
#include <math.h>
#include <mpi.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "allvars.h"
#include "proto.h"

/* -- user defined functions: general -- */

#ifdef OUTPUT_TASK
/*! \brief Output of the task the particles are at.
 *
 *  \param[in] particle (unused)
 *  \param[in] components (unused)
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode (unused)
 *
 *  \return void
 */
static void io_func_task(int particle, int components, void *out_buffer, int mode) { ((int *)out_buffer)[0] = ThisTask; }
#endif

#ifdef OUTPUT_TIMEBIN_HYDRO
/*! \brief Output function of the timebin corresponding to the hydrodynamic
 *         timestep.
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components (unused)
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode (unused)
 *
 *  \return void
 */
static void io_func_timebin_hydro(int particle, int components, void *out_buffer, int mode)
{
  ((int *)out_buffer)[0] = P[particle].TimeBinHydro;
}
#endif

#ifdef OUTPUTTIMESTEP
/*! \brief Output function of the hydrodynamic timestep.
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components (unused)
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode (unused)
 *
 *  \return void
 */
static void io_func_timestep(int particle, int components, void *out_buffer, int mode)
{
  ((MyOutputFloat *)out_buffer)[0] =
      (P[particle].TimeBinHydro ? (((integertime)1) << P[particle].TimeBinHydro) : 0) * All.Timebase_interval;
}
#endif

#ifdef OUTPUT_SOFTENINGS
/*! \brief Output function of the force softening.
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components (unused)
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode Mode (output)
 *
 *  \return void
 */
static void io_func_softenings(int particle, int components, void *out_buffer, int mode)
{
  ((MyOutputFloat *)out_buffer)[0] = All.ForceSoftening[P[particle].SofteningType];
}
#endif

/*! \brief IO function of the particle positions.
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components Number of entries in array.
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode Mode (0: output, 1: input).
 *
 *  \return void
 */
void io_func_pos(const int particle, const int components, void *const buffer, const int mode)
{
  if(mode == 0)
    {
      double pp[3];
      myassert(0 <= components && components <= 3);
      for(int k = 0; k < components; k++)
        pp[k] = P[particle].Pos[k];
#ifdef GRAVITY_NOT_PERIODIC
      const bool gravity_periodic = false;
#else
      const bool gravity_periodic         = true;
#endif
      if(gravity_periodic || P[particle].Type == 0)
        domain_displacePosition(pp, DISPLACE_POSITION_BACKWARD);
      else
        domain_shiftPosition(pp, DISPLACE_POSITION_BACKWARD);
#ifdef OUTPUT_COORDINATES_IN_DOUBLEPRECISION
      const bool output_coordinates_in_dp = true;
#else
      const bool output_coordinates_in_dp = false;
#endif
      MyOutputFloat *ppo = (MyOutputFloat *)buffer;
      double *ppd        = (double *)buffer;
      for(int k = 0; k < components; k++)
        {
          if(output_coordinates_in_dp && DumpFlag != DUMP_BOTH_MINI)
            ppd[k] = pp[k];
          else
            ppo[k] = pp[k];
        }
    }
  else
    {
#ifdef READ_COORDINATES_IN_DOUBLE
      double *in_buffer = (double *)buffer;
#else
      MyInputFloat *in_buffer             = (MyInputFloat *)buffer;
#endif
      for(int k = 0; k < components; k++)
        P[particle].Pos[k] = in_buffer[k];
    }
}

#ifdef OUTPUT_CENTER_OF_MASS
/*! \brief IO function of the particle positions.
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components Number of entries in array.
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode Mode (0: output, 1: input).
 *  \return void
 */

static void io_func_center_of_mass(const int particle, const int components, void *const buffer, const int mode)
{
  if(mode == 0)
    {
      MyOutputFloat *pp = (MyOutputFloat *)buffer;
      for(int k = 0; k < components; k++)
        {
          pp[k] = domain_shiftCoordinate(SphP[particle].Center[k], k, DISPLACE_POSITION_BACKWARD);

#ifdef GRAVITY_NOT_PERIODIC
          if(P[particle].Type != 0)
            continue;
#endif
          if(k == 0)
            pp[k] = WRAP_X(pp[k]);
          if(k == 1)
            pp[k] = WRAP_Y(pp[k]);
          if(k == 2)
            pp[k] = WRAP_Z(pp[k]);
        }
    }
  else
    {
      for(int k = 0; k < components; k++)
        SphP[particle].Center[k] = ((MyInputFloat *)buffer)[k];
    }
}
#endif

/*! \brief IO function for velocities.
 *
 *  Note the different factors of scalefactor in the output than in the code!
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components Number of entries in array.
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode Mode 0: output, 1: input.
 *
 *  \return void
 */
static void io_func_vel(const int particle, const int components, void *const buffer, const int mode)
{
  if(mode == 0)
    {
      for(int k = 0; k < components; k++)
        {
          ((MyOutputFloat *)buffer)[k] = P[particle].Vel[k];
          ((MyOutputFloat *)buffer)[k] *= sqrt(All.cf_a3inv); /* we are dealing with p = a^2 * xdot */
        }
    }
  else
    {
      for(int k = 0; k < components; k++)
        {
          P[particle].Vel[k] = ((MyInputFloat *)buffer)[k];
        }
    }
}

#ifdef OUTPUTACCELERATION
/*! \brief IO function for gravitational accelerations.
 *
 *  Note different a factors in output than in code.
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components Number of entries in array.
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode Mode 0: output, 1: input.
 *
 *  \return void
 */
static void io_func_accel(int particle, int components, void *out_buffer, int mode)
{
  int k;

  if(mode == 0)
    {
      if(RestartFlag != RESTART_SNAP_CONVERSION)
        for(k = 0; k < 3; k++)
          ((MyOutputFloat *)out_buffer)[k] = All.cf_a2inv * P[particle].GravAccel[k];
      else
        for(k = 0; k < 3; k++)
          ((MyOutputFloat *)out_buffer)[k] = P[particle].GravAccel[k];
#ifdef PMGRID
      if(RestartFlag != RESTART_SNAP_CONVERSION)
        for(k = 0; k < 3; k++)
          ((MyOutputFloat *)out_buffer)[k] += All.cf_a2inv * P[particle].GravPM[k];
      else
        for(k = 0; k < 3; k++)
          ((MyOutputFloat *)out_buffer)[k] += P[particle].GravPM[k];
#endif
    }
  else
    {
      for(k = 0; k < 3; k++)
        P[particle].GravAccel[k] = ((MyOutputFloat *)out_buffer)[k];
    }
}
#endif

/* -- user defined functions: shock finder -- */

#if defined(SHOCK_FINDER_BEFORE_OUTPUT_MORE) && !defined(COSMIC_RAYS)
static void io_func_pre_shock_velocity(int particle, int components, void *out_buffer, int mode)
{
  int k;

  for(k = 0; k < components; k++)
    {
      ((MyOutputFloat *)out_buffer)[k] = SphP[particle].VpreShock[k] * sqrt(All.cf_a3inv);
    }
}

static void io_func_post_shock_velocity(int particle, int components, void *out_buffer, int mode)
{
  int k;

  for(k = 0; k < components; k++)
    {
      ((MyOutputFloat *)out_buffer)[k] = SphP[particle].VpostShock[k] * sqrt(All.cf_a3inv);
    }
}
#endif

/* -- user defined functions: GFM -- */

#ifdef OUTPUT_BLACK_HOLE_TIMESTEP
static void io_func_bh_timestep(int particle, int components, void *buffer, int mode)
{
  int k;
  MyOutputFloat timesteps[3], acc[3];

  for(k = 0; k < 3; k++)
    {
      acc[k] = All.cf_a2inv * P[particle].GravAccel[k];
#ifdef PMGRID
      acc[k] += All.cf_a2inv * P[particle].GravPM[k];
#endif
    }

  acc[0] = sqrt(acc[0] * acc[0] + acc[1] * acc[1] + acc[2] * acc[2]);

  timesteps[0] = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * All.ForceSoftening[P[particle].SofteningType] / 2.8 / acc[0]);
  timesteps[1] = 0.25 * BPP(particle).BH_Mass / BPP(particle).BH_Mdot;
  timesteps[2] = All.CourantFac * BPP(particle).BH_DtGasNeighbor;

  for(k = 0; k < 3; k++)
    ((MyOutputFloat *)buffer)[k] = timesteps[k];
}
#endif

#ifdef GFM_DUST
static void io_func_gfm_dust_sputter_tau(int particle, int components, void *buffer, int mode)
{
  ((MyOutputFloat *)buffer)[0] = gfm_get_dust_thermal_sput_tau(particle);
}
#endif

#ifdef GFM_COOLING_METAL
static void io_func_gfm_coolrate(int particle, int components, void *buffer, int mode)
{
  double coolrate, ne, nh0;

  ne = SphP[particle].Ne;
  SetOutputGasState(particle, &ne, &nh0, &coolrate);

  ((MyOutputFloat *)buffer)[0] = coolrate;
}
#endif

#ifdef OUTPUTCOOLRATE
/*! \brief Output function of cooling rate.
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] (unused)
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode (unused)
 *
 *  \return void
 */
static void io_func_coolrate(int particle, int components, void *buffer, int mode)
{
  double tcool, ne, nh0, coolrate;

  ne = SphP[particle].Ne;
  SetOutputGasState(particle, &ne, &nh0, &coolrate);

  /* get cooling time */
  tcool = GetCoolingTime(SphP[particle].Utherm, SphP[particle].Density * All.cf_a3inv, &ne);

  /* convert cooling time with current thermal energy to du/dt */
  if(tcool != 0)
    ((MyOutputFloat *)buffer)[0] = SphP[particle].Utherm / tcool;
  else
    ((MyOutputFloat *)buffer)[0] = 0;
}
#endif

#ifdef GFM_STELLAR_PHOTOMETRICS
static void io_func_gfm_stellar_photometrics(int particle, int components, void *buffer, int mode)
{
  stellar_photometrics st_photo;

  assign_stellar_photometrics(particle, &st_photo);

  ((MyOutputFloat *)buffer)[0] = st_photo.Magnitude_U;
  ((MyOutputFloat *)buffer)[1] = st_photo.Magnitude_B;
  ((MyOutputFloat *)buffer)[2] = st_photo.Magnitude_V;
  ((MyOutputFloat *)buffer)[3] = st_photo.Magnitude_K;
  ((MyOutputFloat *)buffer)[4] = st_photo.Magnitude_g;
  ((MyOutputFloat *)buffer)[5] = st_photo.Magnitude_r;
  ((MyOutputFloat *)buffer)[6] = st_photo.Magnitude_i;
  ((MyOutputFloat *)buffer)[7] = st_photo.Magnitude_z;
}
#endif

#ifdef GFM_STELLAR_EVOLUTION
static void io_func_gfm_metals(int particle, int components, void *buffer, int mode)
{
  int k;

  if(mode == 0)
    {
      if(P[StarP[particle].PID].Mass > 0)
        {
          for(k = 0; k < components; k++)
            ((MyOutputFloat *)buffer)[k] = StarP[particle].MassMetals[k] / P[StarP[particle].PID].Mass;
        }
      else
        {
          for(k = 0; k < components; k++)
            ((MyOutputFloat *)buffer)[k] = StarP[particle].MassMetals[k];
        }
    }
  else
    {
      for(k = 0; k < components; k++)
        StarP[particle].MassMetals[k] = ((MyInputFloat *)buffer)[k];
    }
}
#endif

#ifdef GFM_RPROCESS_CHANNELS
static void io_func_gfm_rprocess(int particle, int components, void *buffer, int mode)
{
  int k;

  if(mode == 0)
    {
      if(P[particle].Type == 4)
        {
          if(P[particle].Mass > 0)
            {
              for(k = 0; k < components; k++)
                ((MyOutputFloat *)buffer)[k] = STP(particle).MassRProcess[k] / P[particle].Mass;
            }
          else
            {
              for(k = 0; k < components; k++)
                ((MyOutputFloat *)buffer)[k] = STP(particle).MassRProcess[k];
            }
        }
      else
        {
          /* type 0 */
          for(k = 0; k < components; k++)
            ((MyOutputFloat *)buffer)[k] = SphP[particle].MassFractionRProcess[k];
        }
    }
  /* we don't read this field at the moment */
}
#endif

/* -- user defined functions: tracers -- */

#ifdef TRACER_MC
/*! \brief Output function of number of tracers in cell
 *
 *  \param[in] particle Index of particle/cell
 *  \param[in] (unused)
 *  \param[out] out_buffer File output buffer
 *  \param[in] mode (unused)
 *
 *  \return void
 */
static void io_func_tracermc_numtracers(int particle, int components, void *buffer, int mode)
{
  if(mode == 0)
    ((int *)buffer)[0] = get_number_of_tracers(particle);
  else
    return;
}

/*! \brief IO function for tracer id
 *
 *  \param[in] particle Index of particle/cell
 *  \param[in] components (unused)
 *  \param[out] out_buffer File output buffer
 *  \param[in] mode Mode (0: output, 1: input)
 *
 *  \return void
 */
static void io_func_tracer_id(int particle, int components, void *buffer, int mode)
{
  if(mode == 0)
    {
      ((MyIDType *)buffer)[0] = TracerLinkedList[particle].ID;
    }
  else
    {
      TracerLinkedList[particle].ID = ((MyIDType *)buffer)[0];
    }
}

/*! \brief IO function for the ID of the parent cell of a tracer
 *
 *  \param[in] particle Index of particle/cell
 *  \param[in] components (unused)
 *  \param[out] buffer File IO buffer
 *  \param[in] mode Mode (0: output, 1: input)
 *
 *  \return void
 */
static void io_func_tracer_parent_id(int particle, int components, void *buffer, int mode)
{
  if(mode == 0)
    {
      ((MyIDType *)buffer)[0] = tracer_cellids[particle];
    }
  else
    {
      tracer_cellids[particle] = ((MyIDType *)buffer)[0];
#ifdef TRACER_MC_CHECKS
      TracerLinkedList[particle].ParentID = tracer_cellids[particle];
#endif
    }
}

#ifdef TRACER_MC_NUM_FLUID_QUANTITIES
/*! \brief IO function for fluid quantities stored in Monte Carlo tracers
 *
 *  \param[in] particle Index of particle/cell
 *  \param[in] components Number of entries in array
 *  \param[out] buffer File IO buffer
 *  \param[in] mode Mode (0: output, 1: input)
 *
 *  \return void
 */
static void io_func_tracer_fluid(int particle, int components, void *buffer, int mode)
{
  int k;

  if(mode == 0)
    {
      for(k = 0; k < components; k++)
        ((MyOutputFloat *)buffer)[k] = TracerLinkedList[particle].fluid_quantities[k];
    }
  else
    {
      for(k = 0; k < components; k++)
        TracerLinkedList[particle].fluid_quantities[k] = ((MyInputFloat *)buffer)[k];
    }
}
#endif
#endif

/* -- user defined functions: gas properties -- */

#if defined(COOLING) || defined(OTVET)
#if !defined(CHIMES) && !defined(GRACKLE)
/*! \brief IO function of the electron number density.
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components (unused)
 *  \param[out] out_buffer File IO buffer.
 *  \param[in] mode Mode 0: output, 1: input.
 *
 *  \return void
 */
static void io_func_ne(int particle, int components, void *buffer, int mode)
{
  if(mode == 0)
    {
#ifdef RT_COOLING_PHOTOHEATING
      ((MyOutputFloat *)buffer)[0] = SphP[particle].n_elec;
      return;
#endif

#if !defined(COOLING) && !defined(SIMPLE_COOLING) && defined(OTVET)
      ((MyOutputFloat *)buffer)[0] = SphP[particle].Ne;
      return;
#endif

      /* normal code path: calculate Ne accounting for GFM options and USE_SFR */
      double ne = SphP[particle].Ne;

#if defined(USE_SFR) && (!defined(LOCAL_FEEDBACK) && !defined(SMUGGLE_STAR_FEEDBACK) && !defined(SFR_MCS))
      /* reproduces previous behavior that Ne is updated prior to output only for Sfr>0 cells */
      /* if this is unwanted (or redundant) this if() condition should be removed */
      double nh0, coolrate;
      if(get_starformation_rate(particle) > 0)
        SetOutputGasState(particle, &ne, &nh0, &coolrate);
#endif

      ((MyOutputFloat *)buffer)[0] = ne;
    }
  else
    {
      SphP[particle].Ne = ((MyInputFloat *)buffer)[0];
    }
}
#endif /* !(CHIMES) && !(GRACKLE) */
#endif

#if defined(COOLING) && !defined(OTVET) && !defined(CHIMES) & !defined(GRACKLE)
/*! \brief Output function for neutral hydrogen fraction.
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components (unused)
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode (unused)
 *
 *  \return void
 */
static void io_func_nh(int particle, int components, void *buffer, int mode)
{
  double ne, nh0, coolrate;

#if defined(RT_COOLING_PHOTOHEATING)
  ((MyOutputFloat *)buffer)[0] = SphP[particle].nHI;
  return;
#endif

  ne = SphP[particle].Ne;
  SetOutputGasState(particle, &ne, &nh0, &coolrate);

  ((MyOutputFloat *)buffer)[0] = nh0;
}
#endif

#if defined(COOLING) && defined(OUTPUT_HE_IONIZATION_STATE) && !defined(OTVET) && !defined(CHIMES) && !defined(GRACKLE)
/*! \brief Output function for neutral helium fraction.
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components (unused)
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode (unused)
 *
 *  \return void
 */
static void io_func_nHe0(int particle, int components, void *buffer, int mode)
{
  double ne, nh0, nHe0, nHep, nHepp, coolrate;

  ne = SphP[particle].Ne;
  SetOutputGasState(particle, &ne, &nh0, &coolrate);

  getHeState(&nHe0, &nHep, &nHepp);

  ((MyOutputFloat *)buffer)[0] = nHe0;
}

/*! \brief Output function for singly-ionized helium fraction.
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components (unused)
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode (unused)
 *
 *  \return void
 */
static void io_func_nHep(int particle, int components, void *buffer, int mode)
{
  double ne, nh0, nHe0, nHep, nHepp, coolrate;

  ne = SphP[particle].Ne;
  SetOutputGasState(particle, &ne, &nh0, &coolrate);

  getHeState(&nHe0, &nHep, &nHepp);

  ((MyOutputFloat *)buffer)[0] = nHep;
}

/*! \brief Output function for doubly-ionized helium fraction.
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components (unused)
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode (unused)
 *
 *  \return void
 */
static void io_func_nHepp(int particle, int components, void *buffer, int mode)
{
  double ne, nh0, nHe0, nHep, nHepp, coolrate;

  ne = SphP[particle].Ne;
  SetOutputGasState(particle, &ne, &nh0, &coolrate);

  getHeState(&nHe0, &nHep, &nHepp);

  ((MyOutputFloat *)buffer)[0] = nHepp;
}
#endif

#ifdef USE_SFR
/*! \brief IO function for star formation rate.
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components (unused)
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode Mode 0: output, 1: input.
 *
 *  \return void
 */
static void io_func_sfr(int particle, int components, void *buffer, int mode)
{
  if(mode == 0)
    {
      if(RestartFlag > RESTART_SNAPSHOT) /* post-processing, don't change the value of Sfr */
        {
          ((MyOutputFloat *)buffer)[0] = SphP[particle].Sfr;
          return;
        }

#ifdef LOCAL_FEEDBACK
      ((MyOutputFloat *)buffer)[0] = SphP[particle].Sfr;
#else
      ((MyOutputFloat *)buffer)[0] = get_starformation_rate(particle);
#endif
    }
  else
    {
      SphP[particle].Sfr = ((MyInputFloat *)buffer)[0];
    }
}
#endif

#ifdef SMUGGLE_OUTPUT_SF_PROBABILITY
static void io_func_sf_prob(int particle, int components, void *buffer, int mode)
{
  ((MyOutputFloat *)buffer)[0] = compute_sf_probability(particle);
}
#endif

#ifdef GRACKLE
static void io_func_grackle_temp(int particle, int components, void *buffer, int mode)
{
  ((MyOutputFloat *)buffer)[0] = get_temp_individual_cell_grackle(particle);
}

static void io_func_grackle_cooltime(int particle, int components, void *buffer, int mode)
{
  ((MyOutputFloat *)buffer)[0] = get_cooling_time_individual_cell_grackle(particle);
}
#endif

/* -- user defined functions: other -- */

#ifdef SGCHEM
static void io_func_sgchem(int particle, int components, void *buffer, int mode)
{
  int k;
  double carb_abund, oxy_abund, m_abund;

  if(mode == 0)
    {
      for(k = 0; k < components; k++)
        ((MyOutputFloat *)buffer)[k] = SphP[particle].TracAbund[k];
    }
  else
    {
      for(k = 0; k < components; k++)
        SphP[particle].TracAbund[k] = ((MyInputFloat *)buffer)[k];

        /* Get local elemental abundances of C, O, M -- already set in IO_SGCHEM_METALS or in the parameter file */
#ifdef SGCHEM_VARIABLE_Z
      carb_abund = SphP[particle].CarbAbund;
      oxy_abund  = SphP[particle].OxyAbund;
      m_abund    = SphP[particle].MAbund;
#else
      carb_abund                   = All.CarbAbund;
      oxy_abund                    = All.OxyAbund;
      m_abund                      = All.MAbund;
#endif

      /* Infer remaining abundances from conservation laws */
#if CHEMISTRYNETWORK == 1
      SphP[particle].TracAbund[IHATOM] =
          1.0 - 2.0 * SphP[particle].TracAbund[IH2] - SphP[particle].TracAbund[IHP] - SphP[particle].TracAbund[IHD];
      SphP[particle].TracAbund[IHEATOM] = HE_ABUND - SphP[particle].TracAbund[IHEP] - SphP[particle].TracAbund[IHEPP];
      SphP[particle].TracAbund[IDATOM]  = All.DeutAbund - SphP[particle].TracAbund[IDP] - SphP[particle].TracAbund[IHD];
#endif
#if CHEMISTRYNETWORK == 4
      SphP[particle].TracAbund[IHATOM] = 1.0 - 2.0 * SphP[particle].TracAbund[IH2] - SphP[particle].TracAbund[IHP];
#endif
#if CHEMISTRYNETWORK == 5
      SphP[particle].TracAbund[IHATOM] = 1.0 - 2.0 * SphP[particle].TracAbund[IH2] - SphP[particle].TracAbund[IHP];
      SphP[particle].TracAbund[ICP]    = carb_abund - SphP[particle].TracAbund[ICO];
#endif
#if CHEMISTRYNETWORK == 7
      SphP[particle].TracAbund[IHATOM]  = 1.0 - 2.0 * SphP[particle].TracAbund[IH2] - SphP[particle].TracAbund[IHP];
      SphP[particle].TracAbund[ICP]     = carb_abund - SphP[particle].TracAbund[ICO];
      SphP[particle].TracAbund[IHEATOM] = HE_ABUND - SphP[particle].TracAbund[IHEP];
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
      SphP[particle].TracAbund[IHATOM] = 1.0 - 2.0 * SphP[particle].TracAbund[IH2] - SphP[particle].TracAbund[IHP] -
                                         SphP[particle].TracAbund[ICHX] - SphP[particle].TracAbund[IOHX] -
                                         SphP[particle].TracAbund[IHCOP];
      SphP[particle].TracAbund[IHEATOM] = HE_ABUND - SphP[particle].TracAbund[IHEP];
      SphP[particle].TracAbund[ICATOM]  = carb_abund - SphP[particle].TracAbund[ICP] - SphP[particle].TracAbund[ICHX] -
                                         SphP[particle].TracAbund[ICO] - SphP[particle].TracAbund[IHCOP];
      SphP[particle].TracAbund[IOATOM] =
          oxy_abund - SphP[particle].TracAbund[IOHX] - SphP[particle].TracAbund[ICO] - SphP[particle].TracAbund[IHCOP];
      SphP[particle].TracAbund[IMATOM] = m_abund - SphP[particle].TracAbund[IMP];
#endif
      for(k = 0; k < SGCHEM_NUM_ADVECTED_SPECIES; k++)
        SphP[particle].MassTracAbund[k] = SphP[particle].TracAbund[k] * P[particle].Mass;
    }
}
#endif

#ifdef SGCHEM_VARIABLE_Z
static void io_func_sgchem_metals(int particle, int components, void *buffer, int mode)
{
  if(mode == 0)
    {
      ((MyOutputFloat *)buffer)[0] = SphP[particle].CarbAbund;
      ((MyOutputFloat *)buffer)[1] = SphP[particle].OxyAbund;
      ((MyOutputFloat *)buffer)[2] = SphP[particle].MAbund;
      ((MyOutputFloat *)buffer)[3] = SphP[particle].ZAtom;
    }
  else
    {
      SphP[particle].CarbAbund = ((MyInputFloat *)buffer)[0];
      SphP[particle].OxyAbund  = ((MyInputFloat *)buffer)[1];
      SphP[particle].MAbund    = ((MyInputFloat *)buffer)[2];
      SphP[particle].ZAtom     = ((MyInputFloat *)buffer)[3];
    }
}
#endif

#if defined(OUTPUT_CURLVEL) || defined(ADJ_BOX_POWERSPEC)
/*! \brief Output function for curl of velocity field.
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components (unused)
 *  \param[out] out_buffer File IO buffer.
 *  \param[in] mode Mode 0: output.
 *
 *  \return void
 */
static void io_func_curlvel(int particle, int components, void *out_buffer, int mode)
{
  if(mode == 0)
    {
      ((MyOutputFloat *)out_buffer)[0] = SphP[particle].CurlVel;
    }
  else
    {
#ifdef ADJ_BOX_POWERSPEC
      if(RestartFlag == RESTART_GAS_VELOCITY_POWER_SPECTRUM)
        SphP[particle].CurlVel = ((MyOutputFloat *)out_buffer)[0];
#endif
    }
}
#endif

#ifdef METALS
static void io_func_metallicity(int particle, int components, void *out_buffer, int mode)
{
  if(mode == 0)
    {
      if(P[particle].Type == 0)
        ((MyOutputFloat *)out_buffer)[0] = SphP[particle].Metallicity;
      else
        ((MyOutputFloat *)out_buffer)[0] = P[particle].Metallicity;
    }
  else
    {
      P[particle].Metallicity = ((MyOutputFloat *)out_buffer)[0];
#if(defined(SN_MCS) || (defined(HII_MCS) && !defined(HII_MCS_TEST)) || defined(PE_MCS)) && \
    !(defined(SN_MCS_INITIAL_DRIVING) || defined(IMF_SAMPLING_MCS) || defined(SB99_FIXED_Z))
      if(P[particle].Type == 4)
        STP(particle).iz = get_sb99_z_index(P[particle].Metallicity);
#endif
    }
}
#endif

#ifdef MRT_LOCAL_FEEDBACK
static void io_func_meandensity(int particle, int components, void *out_buffer, int mode)
{
  if(mode == 0)
    {
      if(P[particle].Type == 4)
        ((MyOutputFloat *)out_buffer)[0] = P[particle].MeanDensity;
    }
}
#endif

#ifdef RT_ADVECT
static void io_func_rt_photons(int particle, int components, void *out_buffer, int mode)
{
  int k;

#ifdef RT_COMBINE_N_DIR_IN_OUTPUT
  for(k = 0; k < RT_N_DIR; k++)
    ((MyOutputFloat *)out_buffer)[0] = ((MyOutputFloat *)out_buffer)[0] + SphP[particle].Photons[k] / 1e53;
#else
  for(k = 0; k < RT_N_DIR; k++)
    ((MyOutputFloat *)out_buffer)[k] = SphP[particle].Photons[k] / 1e53;
#endif
}

static void io_func_rt_photdensity(int particle, int components, void *out_buffer, int mode)
{
  int k;
#ifdef RT_COMBINE_N_DIR_IN_OUTPUT
  for(k = 0; k < RT_N_DIR; k++)
    ((MyOutputFloat *)out_buffer)[0] = ((MyOutputFloat *)out_buffer)[0] + SphP[particle].DensPhot[k] / 1e53;
#else
  for(k = 0; k < RT_N_DIR; k++)
    ((MyOutputFloat *)out_buffer)[k] = SphP[particle].DensPhot[k] / 1e53;
#endif
}
#endif

#ifdef MRT
static void io_func_mrt_photondensity(int particle, int components, void *out_buffer, int mode)
{
  int k;
  if(mode == 0)
    {
      for(k = 0; k < MRT_BINS; k++)
        ((MyOutputFloat *)out_buffer)[k] = SphP[particle].DensPhot[k];
    }
  else
    {
      for(k = 0; k < MRT_BINS; k++)
        SphP[particle].DensPhot[k] = ((MyOutputFloat *)out_buffer)[k];
    }
}

#ifdef MRT_OUTPUT_FLUX
static void io_func_mrt_rtfx(int particle, int components, void *out_buffer, int mode)
{
  int k;
  if(mode == 0)
    {
      for(k = 0; k < MRT_BINS; k++)
        ((MyOutputFloat *)out_buffer)[k] = SphP[particle].RT_F[k][0];
    }
  else
    {
      for(k = 0; k < MRT_BINS; k++)
        SphP[particle].RT_F[k][0] = ((MyOutputFloat *)out_buffer)[k];
    }
}

static void io_func_mrt_rtfy(int particle, int components, void *out_buffer, int mode)
{
  int k;
  if(mode == 0)
    {
      for(k = 0; k < MRT_BINS; k++)
        ((MyOutputFloat *)out_buffer)[k] = SphP[particle].RT_F[k][1];
    }
  else
    {
      for(k = 0; k < MRT_BINS; k++)
        SphP[particle].RT_F[k][1] = ((MyOutputFloat *)out_buffer)[k];
    }
}

static void io_func_mrt_rtfz(int particle, int components, void *out_buffer, int mode)
{
  int k;
  if(mode == 0)
    {
      for(k = 0; k < MRT_BINS; k++)
        ((MyOutputFloat *)out_buffer)[k] = SphP[particle].RT_F[k][2];
    }
  else
    {
      for(k = 0; k < MRT_BINS; k++)
        SphP[particle].RT_F[k][2] = ((MyOutputFloat *)out_buffer)[k];
    }
}
#endif

#ifdef MRT_IR_PHOTON_TRAPPING
static void io_func_mrt_photons(int particle, int components, void *out_buffer, int mode)
{
  int k;
  for(k = UV_BINS; k < (UV_BINS + IR_BINS); k++)
    ((MyOutputFloat *)out_buffer)[k] = SphP[particle].Trapped_DensPhot[k];
}
#endif
#endif

#if defined(DL_RADIATION) && defined(DL_OUTPUT_RT_FLUX)
static void io_func_dust_mrt_rtfx(int particle, int components, void *out_buffer, int mode)
{
  int k;
  for(k = 0; k < MRT_BINS; k++)
    ((MyOutputFloat *)out_buffer)[k] = DustP[particle].LocalRT_F[k][0];
}

static void io_func_dust_mrt_rtfy(int particle, int components, void *out_buffer, int mode)
{
  int k;
  for(k = 0; k < MRT_BINS; k++)
    ((MyOutputFloat *)out_buffer)[k] = DustP[particle].LocalRT_F[k][1];
}

static void io_func_dust_mrt_rtfz(int particle, int components, void *out_buffer, int mode)
{
  int k;
  for(k = 0; k < MRT_BINS; k++)
    ((MyOutputFloat *)out_buffer)[k] = DustP[particle].LocalRT_F[k][2];
}
#endif

#ifdef OUTPUT_VORTICITY
/*! \brief Output function of vorticity (calculated from velocity spatial
 *         derivatives).
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components (unused)
 *  \param[out] out_buffer File IO buffer.
 *  \param[in] mode Mode 0: output
 *
 *  \return void
 */
static void io_func_vorticity(int particle, int components, void *out_buffer, int mode)
{
  if(mode == 0)
    {
      ((MyOutputFloat *)out_buffer)[0] = SphP[particle].Grad.dvel[2][1] - SphP[particle].Grad.dvel[1][2];
      ((MyOutputFloat *)out_buffer)[1] = SphP[particle].Grad.dvel[0][2] - SphP[particle].Grad.dvel[2][0];
      ((MyOutputFloat *)out_buffer)[2] = SphP[particle].Grad.dvel[1][0] - SphP[particle].Grad.dvel[0][1];
    }
  else
    {
#ifdef ADJ_BOX_POWERSPEC
      if(RestartFlag == RESTART_GAS_VELOCITY_POWER_SPECTRUM)
        for(int k = 0; k < 3; k++)
          SphP[particle].Vorticity[k] = ((MyOutputFloat *)out_buffer)[k];
#endif
    }
}
#endif

#ifdef BRAGINSKII_OUTPUT_PRESSURE_ANISOTROPY
/*! \brief Output function for pressure anisotropy (calculated from velocity
 *       spatial derivatives, density, pressure and magnetic field direction).
 *
 *    Δp = ρ ν_∥ (3bb:∇v - ∇⋅v)
 *
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components (unused)
 *  \param[out] out_buffer File IO buffer.
 *  \param[in] mode Mode 0: output (unused)
 *
 *  \return void
 */
static void io_func_pressure_anisotropy(int particle, int components, void *out_buffer, int mode)
{
  double B2  = 0.0;
  double res = 0;

  // B² Magnetic field strength squared
  for(int i = 0; i < 3; i++)
    B2 += SphP[particle].B[i] * SphP[particle].B[i];

  // b magnetic field direction
  double b[3];
  for(int i = 0; i < 3; i++)
    b[i] = SphP[particle].B[i] / sqrt(B2);

  // 3 bb:∇v = 3 ∑ᵢⱼ bᵢbⱼ ∂ᵢvⱼ
  for(int i = 0; i < 3; i++)
    {
      for(int j = 0; j < 3; j++)
        {
          res += 3.0 * b[i] * b[j] * SphP[particle].Grad.dvel[i][j];
        }
    }

  // Subtract velocity divergence (∇⋅v)
  for(int i = 0; i < 3; i++)
    {
      res -= SphP[particle].Grad.dvel[i][i];
    }

// Calculate viscosity coefficieint ν_∥
#ifdef BRAGINSKII_SPITZER
  double T  = SphP[particle].Pressure / SphP[particle].Density;
  double nu = All.BragViscosityCoefficient * pow(T, 2.5) / SphP[particle].Density;
#else
  double nu = All.BragViscosityCoefficient;
#endif

  // Return Δp = ρ ν_∥ (3bb:∇v - ∇⋅v)
  ((MyOutputFloat *)out_buffer)[0] = nu * SphP[particle].Density * res;
}
#endif

#ifdef BRAGINSKII_OUTPUT_NUM_SUBSTEPS
/*! \brief Output function of the number of substeps used
 *         for Braginskii visosity (when using RKL2 or subcycling)
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components (unused)
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode (unused)
 *
 *  \return void
 */
static void io_func_brag_substeps(int particle, int components, void *out_buffer, int mode)
{
  ((int *)out_buffer)[0] = SphP[particle].BragViscositySubsteps;
}
#endif

#ifdef OUTPUT_CELL_SPIN
/*! \brief Output function of cell spin.
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components (unused)
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode (unused)
 *
 *  \return void
 */
static void io_func_cell_spin(int particle, int components, void *out_buffer, int mode)
{
  ((MyOutputFloat *)out_buffer)[0] = SphP[particle].Spin[0] / P[particle].Mass;
  ((MyOutputFloat *)out_buffer)[1] = SphP[particle].Spin[1] / P[particle].Mass;
  ((MyOutputFloat *)out_buffer)[2] = SphP[particle].Spin[2] / P[particle].Mass;
}
#endif

#ifdef MHD
/*! \brief IO function for magnetic field.
 *
 *  Note that the output is in Gauss unit system (in code units) while the
 *  internal B-field is in Heaviside-Lorentz system (FACTOR of sqrt(4 PI)!).
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components (unused)
 *  \param[out] out_buffer File IO buffer.
 *  \param[in] mode Mode 0: output, 1: input.
 *
 *  \return void
 */
static void io_func_bfield(int particle, int components, void *out_buffer, int mode)
{
  int k;

  if(mode == 0)
    {
      /* writing: convert from Heavyside-Lorentz to Gauss */
      for(k = 0; k < 3; k++)
        ((MyOutputFloat *)out_buffer)[k] = SphP[particle].B[k] * sqrt(4. * M_PI);
    }
  else
    {
      /* reading: convert from Gauss to Heavyside-Lorentz */
      for(k = 0; k < 3; k++)
        SphP[particle].B[k] = ((MyInputFloat *)out_buffer)[k] / sqrt(4. * M_PI);
    }
}
#endif

#ifdef GFM_DUST
static void io_func_gfm_dust_agb(int particle, int components, void *buffer, int mode)
{
  int k;

  if(mode == 0)
    {
      for(k = 0; k < components; k++)
        ((MyOutputFloat *)buffer)[k] = SphP[particle].MetalsDustFraction[GFM_DUST_AGB][k];
    }
  else
    {
      for(k = 0; k < components; k++)
        SphP[particle].MetalsDustFraction[GFM_DUST_AGB][k] = ((MyOutputFloat *)buffer)[k];
    }
}

static void io_func_gfm_dust_agb2(int particle, int components, void *buffer, int mode)
{
  int k;

  if(mode == 0)
    {
      for(k = 0; k < components; k++)
        ((MyOutputFloat *)buffer)[k] =
            StarP[particle].MassMetals[k] / P[StarP[particle].PID].Mass * StarP[particle].InitialDustFractions[GFM_DUST_AGB][k];
    }
}

static void io_func_gfm_dust_snii(int particle, int components, void *buffer, int mode)
{
  int k;

  if(mode == 0)
    {
      for(k = 0; k < components; k++)
        ((MyOutputFloat *)buffer)[k] = SphP[particle].MetalsDustFraction[GFM_DUST_SNII][k];
    }
  else
    {
      for(k = 0; k < components; k++)
        SphP[particle].MetalsDustFraction[GFM_DUST_SNII][k] = ((MyOutputFloat *)buffer)[k];
    }
}

static void io_func_gfm_dust_snii2(int particle, int components, void *buffer, int mode)
{
  int k;

  if(mode == 0)
    {
      for(k = 0; k < components; k++)
        ((MyOutputFloat *)buffer)[k] =
            StarP[particle].MassMetals[k] / P[StarP[particle].PID].Mass * StarP[particle].InitialDustFractions[GFM_DUST_SNII][k];
    }
}

static void io_func_gfm_dust_snia(int particle, int components, void *buffer, int mode)
{
  int k;

  if(mode == 0)
    {
      for(k = 0; k < components; k++)
        ((MyOutputFloat *)buffer)[k] = SphP[particle].MetalsDustFraction[GFM_DUST_SNIa][k];
    }
  else
    {
      for(k = 0; k < components; k++)
        SphP[particle].MetalsDustFraction[GFM_DUST_SNIa][k] = ((MyOutputFloat *)buffer)[k];
    }
}

static void io_func_gfm_dust_snia2(int particle, int components, void *buffer, int mode)
{
  int k;

  if(mode == 0)
    {
      for(k = 0; k < components; k++)
        ((MyOutputFloat *)buffer)[k] =
            StarP[particle].MassMetals[k] / P[StarP[particle].PID].Mass * StarP[particle].InitialDustFractions[GFM_DUST_SNIa][k];
    }
}

static void io_func_gfm_dust_metal(int particle, int components, void *buffer, int mode)
{
  int k;

  if(mode == 0)
    {
      ((MyOutputFloat *)buffer)[0] = 0.0;
      for(k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
        {
          ((MyOutputFloat *)buffer)[0] += SphP[particle].MetalsDustFraction[GFM_DUST_AGB][k];
          ((MyOutputFloat *)buffer)[0] += SphP[particle].MetalsDustFraction[GFM_DUST_SNII][k];
          ((MyOutputFloat *)buffer)[0] += SphP[particle].MetalsDustFraction[GFM_DUST_SNIa][k];
        }
    }
}
static void io_func_gfm_dust_metal2(int particle, int components, void *buffer, int mode)
{
  int k;

  if(mode == 0)
    {
      ((MyOutputFloat *)buffer)[0] = 0.0;
      for(k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
        {
          /* Includes gas-phase and dust contributions. */
          double metal_frac_k = StarP[particle].MassMetals[k] / P[StarP[particle].PID].Mass;

          ((MyOutputFloat *)buffer)[0] += (metal_frac_k * StarP[particle].InitialDustFractions[GFM_DUST_AGB][k]);
          ((MyOutputFloat *)buffer)[0] += (metal_frac_k * StarP[particle].InitialDustFractions[GFM_DUST_SNII][k]);
          ((MyOutputFloat *)buffer)[0] += (metal_frac_k * StarP[particle].InitialDustFractions[GFM_DUST_SNIa][k]);
        }
    }
}
#endif

#ifdef GFM_CHEMTAGS
static void io_func_gfm_chemtags_gas(int particle, int components, void *out_buffer, int mode)
{
  if(mode == 0)
    {
      for(int k = 0; k < components; k++)
        ((MyOutputFloat *)out_buffer)[k] = SphP[particle].MassMetalsChemTags[k] / P[particle].Mass;
    }
  else
    {
      for(int k = 0; k < GFM_N_CHEM_TAGS; k++)
        SphP[particle].MassMetalsChemTags[k] =
            ((MyOutputFloat *)out_buffer)[k]; /* we will multiple this with the mass in the read-routine  */
    }
}

static void io_func_gfm_chemtags(int particle, int components, void *out_buffer, int mode)
{
  if(mode == 0)
    {
      for(int k = 0; k < components; k++)
        ((MyOutputFloat *)out_buffer)[k] = StarP[particle].MassMetalsChemTags[k] / P[StarP[particle].PID].Mass;
    }
  else
    {
      for(int k = 0; k < GFM_N_CHEM_TAGS; k++)
        StarP[particle].MassMetalsChemTags[k] =
            ((MyOutputFloat *)out_buffer)[k]; /* we will multiple this with the mass in the read-routine once StarP[].PID is set */
    }
}
#endif

#if defined(TREECOLV2) && defined(OUTPUTCOL)
static void io_func_treecolumn(int particle, int components, void *out_buffer, int mode)
{
  int k;
  int offset = 0;

  for(k = 0; k < NPIX; k++)
    ((MyOutputFloat *)out_buffer)[k] = SphP[particle].Projection[k];
  offset += NPIX;

#ifdef TREECOLV2_H2
  for(k = 0; k < NPIX; k++)
    ((MyOutputFloat *)out_buffer)[offset + k] = SphP[particle].ProjectionH2[k];
  offset += NPIX;
#endif

#ifdef TREECOLV2_CO
  for(k = 0; k < NPIX; k++)
    ((MyOutputFloat *)out_buffer)[offset + k] = SphP[particle].ProjectionCO[k];
  offset += NPIX;
#endif

#ifdef TREECOLV2_C
  for(k = 0; k < NPIX; k++)
    ((MyOutputFloat *)out_buffer)[offset + k] = SphP[particle].ProjectionC[k];
#endif
}
#endif

#ifdef OUTPUT_DENSTROPHY
static void io_func_denstrophy(int particle, int components, void *out_buffer, int mode)
{
  double curl_rhovel_x = sqrt(SphP[particle].Density) * (SphP[particle].Grad.dvel[2][1] - SphP[particle].Grad.dvel[1][2]) -
                         1. / (2. * sqrt(SphP[particle].Density)) *
                             (P[particle].Vel[1] * SphP[particle].Grad.drho[2] - P[particle].Vel[2] * SphP[particle].Grad.drho[1]);
  double curl_rhovel_y = sqrt(SphP[particle].Density) * (SphP[particle].Grad.dvel[0][2] - SphP[particle].Grad.dvel[2][0]) -
                         1. / (2. * sqrt(SphP[particle].Density)) *
                             (P[particle].Vel[2] * SphP[particle].Grad.drho[0] - P[particle].Vel[0] * SphP[particle].Grad.drho[2]);
  double curl_rhovel_z = sqrt(SphP[particle].Density) * (SphP[particle].Grad.dvel[1][0] - SphP[particle].Grad.dvel[0][1]) -
                         1. / (2. * sqrt(SphP[particle].Density)) *
                             (P[particle].Vel[0] * SphP[particle].Grad.drho[1] - P[particle].Vel[1] * SphP[particle].Grad.drho[0]);

  ((MyOutputFloat *)out_buffer)[0] =
      0.5 * (curl_rhovel_x * curl_rhovel_x + curl_rhovel_y * curl_rhovel_y + curl_rhovel_z * curl_rhovel_z);
}
#endif

#ifdef CHIMES
static void io_func_chimes_abundances(int particle, int components, void *buffer, int mode)
{
  int k;

  if(mode == 0)
    {
      chimes_set_pointers(particle);
      for(k = 0; k < components; k++)
        ((MyOutputFloat *)buffer)[k] = SphP[particle].ChimesGasVars.abundances[k];
    }
#ifndef CHIMES_INITIALISE_IN_EQM
  else
    {
      SphP[particle].ChimesGasVars.index = particle;
      chimes_set_pointers(particle);
      for(k = 0; k < components; k++)
        SphP[particle].ChimesGasVars.abundances[k] = ((MyInputFloat *)buffer)[k];
    }
#endif
}

static void io_func_chimes_mu(int particle, int components, void *buffer, int mode)
{
  if(mode == 0)
    ((MyOutputFloat *)buffer)[0] = calculate_mean_molecular_weight(&(SphP[particle].ChimesGasVars), &ChimesGlobalVars);
}
#endif /* CHIMES */

#ifdef SFR_MCS
static void io_func_stellar_birthtime(int particle, int components, void *buffer, int mode)
{
  if(mode == 0)
    ((MyOutputFloat *)buffer)[0] = StarP[particle].BirthTime;
  else
    {
      StarP[particle].BirthTime = ((MyInputFloat *)buffer)[0];
#if defined(SN_MCS) || (defined(HII_MCS) && !defined(HII_MCS_TEST)) || defined(PE_MCS)
      StarP[particle].Age = get_time_difference_in_Gyr(StarP[particle].BirthTime, All.Time) * 1e9;
#if defined(HII_MCS) && !(defined(HII_MCS_TEST) || defined(IMF_SAMPLING_MCS))
      int it_high;
      for(it_high = 1; it_high < sb99.N_t_Photons && StarP[particle].Age > sb99.TimestepsPhotons[it_high]; it_high++)
        ;

      StarP[particle].photon_it_high = it_high;
#endif
#if defined(PE_MCS) && !defined(IMF_SAMPLING_MCS)
      int it_fuv_high;
      for(it_fuv_high = 1; it_fuv_high < sb99.N_t_FUV && StarP[particle].Age > sb99.TimestepsFUV[it_fuv_high]; it_fuv_high++)
        ;

      StarP[particle].fuv_it_high = it_fuv_high;
#endif
#endif
    }
}
static void io_func_stellar_age_in_Gyr(int particle, int components, void *buffer, int mode)
{
  if(mode == 0)
    ((MyOutputFloat *)buffer)[0] = (MyOutputFloat)get_time_difference_in_Gyr(StarP[particle].BirthTime, All.Time);
  /* This field is written but not read */
}
#ifdef PE_MCS
static void io_func_pe_heating(int particle, int components, void *buffer, int mode)
{
  if(mode == 0)
    {
      ((MyOutputFloat *)buffer)[0] = (MyOutputFloat)calculate_pe_heating_rate(particle);
    }
}
#endif
#endif /* SFR_MCS */

#ifdef TURB_APPROX_MCS_OUTPUT_RATES
static void io_func_turbulent_production_rate(int particle, int components, void *buffer, int mode)
{
  if(mode == 0)
    ((MyOutputFloat *)buffer)[0] = (MyOutputFloat)get_turbulent_production_rate_single(particle);
}

static void io_func_turbulent_dissipation_rate(int particle, int components, void *buffer, int mode)
{
  if(mode == 0)
    {
      double kturb                 = SphP[particle].TurbEnergy / SphP[particle].Volume;
      double dx                    = pow(6.0 * SphP[particle].Volume / M_PI, 1.0 / 3.0);
      ((MyOutputFloat *)buffer)[0] = -1.58 * SphP[particle].Volume * sqrt(kturb / SphP[particle].Density) * kturb / dx;
      ;
    }
}

static void io_func_turbulent_adiabatic_rate(int particle, int components, void *buffer, int mode)
{
  if(mode == 0)
    {
      double divv = SphP[particle].dvel_unlim[0][0] + SphP[particle].dvel_unlim[1][1] + SphP[particle].dvel_unlim[2][2];
      ((MyOutputFloat *)buffer)[0] = -2.0 * SphP[particle].TurbEnergy * divv / 3.0;
    }
}
#endif /* TURB_APPROX_MCS_OUTPUT_RATES */

/*! \brief Function for field registering.
 *
 *  For init_field arguments read the description of init_field.
 *  Don't forget to add the new IO_FLAG to allvars.h.
 *
 *  \return void
 */
void init_io_fields(void)
{
  /* ALL TYPES */

#ifdef OUTPUT_COORDINATES_IN_DOUBLEPRECISION
  enum types_in_file pos_out = FILE_DOUBLE;
#else
  enum types_in_file pos_out = FILE_MY_IO_FLOAT;
#endif
#ifdef READ_COORDINATES_IN_DOUBLE
  enum types_in_file pos_in = FILE_DOUBLE;
#else
  enum types_in_file pos_in = FILE_MY_IO_FLOAT;
#endif
  init_field(IO_POS, "POS ", "Coordinates", MEM_DOUBLE, pos_out, pos_in, 3, A_NONE, 0, io_func_pos, ALL_TYPES);
  init_units(IO_POS, 1., -1., 1., 0., 0., All.UnitLength_in_cm);

  init_field(IO_POS_MINI, "POS ", "Coordinates", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 3, A_NONE, 0, io_func_pos, ALL_TYPES);
  init_units(IO_POS_MINI, 1., -1., 1., 0., 0., All.UnitLength_in_cm);
  /* second IO tag output to mini-snaps always in single precision */
  init_snapshot_type(IO_POS_MINI, SN_MINI_ONLY);

  /* particle velocities */
  init_field(IO_VEL, "VEL ", "Velocities", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_NONE, 0, io_func_vel, ALL_TYPES);
  /* sqrt(a)*km/s */
  init_units(IO_VEL, 0.5, 0., 0., 0., 1., All.UnitVelocity_in_cm_per_s);
  init_snapshot_type(IO_VEL, SN_MINI);

  init_field(IO_ID, "ID  ", "ParticleIDs", MEM_MY_ID_TYPE, FILE_MY_ID_TYPE, FILE_MY_ID_TYPE, 1, A_P, &P[0].ID, 0, ALL_TYPES);
  init_units(IO_ID, 0, 0, 0, 0, 0, 0);
  init_snapshot_type(IO_ID, SN_MINI);

  /* particle mass */
  init_field(IO_MASS, "MASS", "Masses", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_P, &P[0].Mass, 0,
             SET_IN_GET_PARTICLES_IN_BLOCK);
  init_units(IO_MASS, 0., -1., 0., 1., 0., All.UnitMass_in_g);
  init_snapshot_type(IO_MASS, SN_MINI);

#ifdef OUTPUTPOTENTIAL
  /* gravitational potential */
  init_field(IO_POT, "POT ", "Potential", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_P, &P[0].Potential, 0, ALL_TYPES);
  init_units(IO_POT, -1.0, 0.0, 0.0, 0.0, 2.0, All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s); /* (km/s)^2/a */

  init_field(IO_POT_MINI, "POT ", "Potential", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_P, &P[0].Potential, 0,
             STARS_ONLY | BHS_ONLY);
  init_units(IO_POT_MINI, -1.0, 0.0, 0.0, 0.0, 2.0, All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s);
  /* second IO tag output to mini-snaps for stars/BHs only */
  init_snapshot_type(IO_POT_MINI, SN_MINI_ONLY);
#endif

#ifdef COMPUTE_VORONOI_DM_DENSITY_IN_POSTPROC
  init_field(IO_VDM, "DMVD", "DM_VoronoiDensity", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_P, &P[0].DM_VoronoiDensity, 0,
             ALL_TYPES);
#endif

  /* GAS CELLS */

  /* internal energy */
  init_field(IO_U, "U   ", "InternalEnergy", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].Utherm, 0,
             GAS_ONLY);
  init_units(IO_U, 0., 0., 0., 0., 2., All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s);
  init_snapshot_type(IO_U, SN_MINI);

  /* particle density */
  init_field(IO_RHO, "RHO ", "Density", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].Density, 0, GAS_ONLY);
  init_units(IO_RHO, -3., 2., -3., 1., 0., All.UnitDensity_in_cgs);
  init_snapshot_type(IO_RHO, SN_MINI);

#ifdef OUTPUT_PRESSURE
  init_field(IO_PRESSURE, "PRES", "Pressure", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].Pressure, 0, GAS_ONLY);
  init_units(IO_PRESSURE, -3.0, 2.0, -3.0, 1.0, 2.0,
             All.UnitDensity_in_cgs * All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s);
#endif

#ifdef BRAGINSKII_OUTPUT_PRESSURE_ANISOTROPY
  init_field(IO_PRESSURE_ANISOTROPY, "PRESANIS", "PressureAnisotropy", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_NONE, 0,
             io_func_pressure_anisotropy, GAS_ONLY);
  init_units(IO_PRESSURE_ANISOTROPY, -3.0, 2.0, -3.0, 1.0, 2.0,
             All.UnitDensity_in_cgs * All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s);
#endif
#ifdef BRAGINSKII_OUTPUT_NUM_SUBSTEPS
  init_field(IO_BRAG_SUBSTEPS, "BragSubsteps", "BragViscositySubsteps", MEM_NONE, FILE_INT, FILE_NONE, 1, A_NONE, 0,
             io_func_brag_substeps, GAS_ONLY);
#endif

#ifdef OUTPUT_ENTROPY
  init_field(IO_ENTROPY, "S   ", "Entropy", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].Entropy, 0, GAS_ONLY);
#endif

#ifdef OUTPUT_CSND
  init_field(IO_CSND, "CSND", "SoundSpeed", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].Csnd, 0, GAS_ONLY);
#endif

#ifdef OUTPUT_MACHNUM
  init_field(IO_MACH, "MACH", "MachNumber", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].MaxMachNumber, 0, GAS_ONLY);
  /* dimensionless ratio of velocities */
  init_units(IO_MACH, 0, 0, 0, 0, 0, 0);
#endif

#if defined(COOLING) && !defined(SIMPLE_COOLING) && !defined(OTVET) && !defined(GRACKLE) && !defined(CHIMES)
  /* electron abundance */
  init_field(IO_NE, "NE  ", "ElectronAbundance", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_NONE, 0, io_func_ne, GAS_ONLY);
  init_units(IO_NE, 0, 0, 0, 0, 0, 0); /* dimensionless fraction */
  init_snapshot_type(IO_NE, SN_MINI);

  /* neutral hydrogen fraction */
  init_field(IO_NH, "NH  ", "NeutralHydrogenAbundance", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_NONE, 0, io_func_nh, GAS_ONLY);
  /* dimensionless fraction */
  init_units(IO_NH, 0, 0, 0, 0, 0, 0);

#ifdef OUTPUT_HE_IONIZATION_STATE
  /* neutral helium fraction */
  init_field(IO_NHE0, "NHeO", "NeutralHeliumAbundance", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_NONE, 0, io_func_nHe0, GAS_ONLY);
  /* dimensionless fraction */
  init_units(IO_NHE0, 0, 0, 0, 0, 0, 0);

  /* singly-ionized helium fraction */
  init_field(IO_NHEP, "NHep", "SinglyIonizedHeliumAbundance", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_NONE, 0, io_func_nHep,
             GAS_ONLY);
  /* dimensionless fraction */
  init_units(IO_NHEP, 0, 0, 0, 0, 0, 0);

  /* doubly-ionized helium fraction */
  init_field(IO_NHEPP, "Hepp", "DoublyIonizedHeliumAbundance", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_NONE, 0, io_func_nHepp,
             GAS_ONLY);
  /* dimensionless fraction */
  init_units(IO_NHEPP, 0, 0, 0, 0, 0, 0);
#endif

#endif

#ifdef USE_SFR
  /* star formation rate */
  init_field(IO_SFR, "SFR ", "StarFormationRate", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_NONE, 0, io_func_sfr, GAS_ONLY);
  init_units(IO_SFR, 0.0, 0.0, -1.0, 1.0, 1.0, SOLAR_MASS / SEC_PER_YEAR); /* Msun/yr */
  init_snapshot_type(IO_SFR, SN_MINI);
#endif

#ifdef OUTPUT_DIVVEL
  init_field(IO_DIVVEL, "DIVV", "VelocityDivergence", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].DivVel, 0,
             GAS_ONLY);
#endif

#if defined(OUTPUT_CURLVEL) || defined(ADJ_BOX_POWERSPEC)
  init_field(IO_CURLVEL, "ROTV", "VelocityCurl", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_NONE, 0, io_func_curlvel,
             GAS_ONLY);
#endif

#ifdef OUTPUT_COOLHEAT
  init_field(IO_COOLHEAT, "COHE", "CoolingHeatingEnergy", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].CoolHeat, 0,
             GAS_ONLY);
#endif

#ifdef OUTPUT_SURFACE_AREA
  init_field(IO_SAREA, "AREA", "SurfaceArea", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].SurfaceArea, 0,
             GAS_ONLY);
  init_field(IO_NFACES, "NFAC", "NumFacesCell", MEM_INT, FILE_INT, FILE_INT, 1, A_SPHP, &SphP[0].CountFaces, 0, GAS_ONLY);
#endif

#ifdef SMUGGLE_OUTPUT_SF_PROBABILITY
  init_field(IO_SF_PROBABILITY, "SFPR", "SFProbability", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_NONE, 0, io_func_sf_prob,
             GAS_ONLY);
#endif

#ifdef OUTPUTCOOLRATE
  init_field(IO_COOLRATE, "COOR", "CoolingRate", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_NONE, 0, io_func_coolrate, GAS_ONLY);
#endif

#ifdef MEASURE_DISSIPATION_RATE
  init_field(IO_DUDT, "DUDT", "DissipationRate_DuDt", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].DuDt, 0,
             GAS_ONLY);
#endif

#ifdef SMUGGLE_OUTPUT_MOLECULAR_FRACTION
  init_field(IO_MOLECULAR_FRAC, "HMOL", "MolecularHFrac", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].MolecularFrac,
             0, GAS_ONLY);
#endif

#ifdef OUTPUT_VORTICITY
  init_field(IO_VORT, "VORT", "Vorticity", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_NONE, 0, io_func_vorticity, GAS_ONLY);
#endif

#ifdef OUTPUT_CELL_SPIN
  init_field(IO_SPIN, "SPIN", "CellSpin", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, 3, A_NONE, 0, io_func_cell_spin, GAS_ONLY);
  init_field(IO_LPOS, "LPOS", "CellSpinCenter", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_SPHP, &SphP[0].CenterOffset[0],
             0, GAS_ONLY);
#endif

#ifdef OUTPUT_QRAD
  init_field(IO_QRAD, "QRAD", "RadiativeEnergyTransport", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].Qrad, 0,
             GAS_ONLY);
  init_field(IO_RADIALQRAD, "RADIALQRAD", "RadialRadTrans", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].RadialQrad,
             0, GAS_ONLY);
#endif

  /* GAS CELLS (GRADIENTS) */

#ifdef OUTPUT_PRESSURE_GRADIENT
#ifdef TGSET
  init_field(IO_GRADP, "GRAP", "PressureGradient", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_SPHP,
             &SphP[0].Grad.dpress[0], 0, GAS_ONLY);
#else
  init_field(IO_GRADP, "GRAP", "PressureGradient", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_NONE, 3, A_SPHP, &SphP[0].Grad.dpress[0], 0,
             GAS_ONLY);
#endif
#endif

#ifdef OUTPUT_DENSITY_GRADIENT
  init_field(IO_GRADR, "GRAR", "DensityGradient", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_NONE, 3, A_SPHP, &SphP[0].Grad.drho[0], 0,
             GAS_ONLY);
#endif

#ifdef OUTPUT_VELOCITY_GRADIENT
  init_field(IO_GRADV, "GRAV", "VelocityGradient", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_NONE, 9, A_SPHP, &SphP[0].Grad.dvel[0][0], 0,
             GAS_ONLY);
#endif

#ifdef OUTPUT_BFIELD_GRADIENT
  init_field(IO_GRADB, "GRAB", "BfieldGradient", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_NONE, 9, A_SPHP, &SphP[0].Grad.dB[0][0], 0,
             GAS_ONLY);
#endif

  /* GAS CELLS (MESH PROPERTIES) */

#ifdef OUTPUT_VOLUME
  init_field(IO_VOL, "VOL ", "Volume", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].Volume, 0, GAS_ONLY);
  init_units(IO_VOL, 3., -3., 3., 0., 0., All.UnitLength_in_cm * All.UnitLength_in_cm * All.UnitLength_in_cm);
#endif

#ifdef OUTPUT_VERTEX_VELOCITY
  init_field(IO_VERTEXVEL, "VEVE", "VertexVelocity", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_SPHP,
             &SphP[0].VelVertex[0], 0, GAS_ONLY);
  init_units(IO_VERTEXVEL, 0.5, 0., 0., 0., 1., All.UnitVelocity_in_cm_per_s);
#endif

#ifdef OUTPUT_VERTEX_VELOCITY_DIVERGENCE
  init_field(IO_DIVVERTEXVEL, "DIVV", "VertexVelocityDiv", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP,
             &SphP[0].DivVelVertex, 0, GAS_ONLY);
#endif

#ifdef OUTPUT_MESH_FACE_ANGLE
  init_field(IO_FACEANGLE, "FACA", "MaxFaceAngle", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].MaxFaceAngle, 0,
             GAS_ONLY);
#endif

#ifdef OUTPUT_CENTER_OF_MASS
  init_field(IO_CM, "CMCE", "CenterOfMass", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_NONE, 0, io_func_center_of_mass,
             GAS_ONLY);
  init_units(IO_CM, 1., -1., 1., 0., 0., All.UnitLength_in_cm);
#endif

  /* GAS CELLS (GFM) */

#if defined(GFM_OUTPUT_MASK) && ((GFM_OUTPUT_MASK)&1) /* is the condition on GFM_COOLING_METAL really needed? TBD */
#if defined(GFM_COOLING_METAL) && (defined(GFM_WINDS_VARIABLE) || defined(GFM_WINDS_LOCAL))
  init_field(IO_GFM_WINDHOSTMASS, "GWHM", "GFM_WindHostHaloMass", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP,
             &SphP[0].w.HostHaloMass, 0, GAS_ONLY);
  init_units(IO_GFM_WINDHOSTMASS, 0., -1., 0., 1., 0., All.UnitMass_in_g);
#endif
#if defined(GFM_COOLING_METAL) && defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE == 1)
  init_field(IO_GFM_WINDHOSTDISP, "GWDV", "GFM_WindDMVelDisp", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP,
             &SphP[0].w.DMVelDisp, 0, GAS_ONLY);
  init_units(IO_GFM_WINDHOSTDISP, 0.0, 0.0, 0.0, 0.0, 1.0, All.UnitVelocity_in_cm_per_s);
#endif
#endif

#if defined(GFM_COOLING_METAL) && ((GFM_OUTPUT_MASK)&2)
  init_field(IO_GFM_COOLRATE, "GCOL", "GFM_CoolingRate", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_NONE, 0, io_func_gfm_coolrate,
             GAS_ONLY);
  init_units(IO_GFM_COOLRATE, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0); /* rate of energy change (Lambda_net/n_H^2) in cgs units of cm^3*erg/s */
#endif

#if defined(GFM_AGN_RADIATION) && ((GFM_OUTPUT_MASK)&64)
  init_field(IO_GFM_AGN_RADIATION, "AGNR", "GFM_AGNRadiation", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP,
             &SphP[0].AGNBolIntensity, 0, GAS_ONLY);
  init_units(IO_GFM_AGN_RADIATION, 0.0, 0.0, -3.0, 1.0, 3.0, 1.0); /* bolomeric intensity in physical cgs units of erg/s/cm^2 */
#endif

#ifdef TURBULENT_METALDIFFUSION
#ifndef GFM_STELLAR_EVOLUTION
  init_field(IO_GFM_METALS, "GMET", "GFM_Metals", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, GFM_N_CHEM_ELEMENTS, A_SPHP,
             &SphP[0].MetalsFraction, 0, GAS_ONLY);
  init_units(IO_GFM_METALS, 0, 0, 0, 0, 0, 0); /* dimensionless abundance ratios */
#endif
#endif

  /* STARS */

  if(All.StarformationOn)
    {
#ifdef GFM
      init_field(IO_GFM_AGE, "GAGE", "GFM_StellarFormationTime", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_STARP,
                 &StarP[0].BirthTime, 0, STARS_ONLY);
      init_units(IO_GFM_AGE, 0, 0, 0, 0, 0, 0); /* scale factor */
      init_snapshot_type(IO_GFM_AGE, SN_MINI);
#endif

#if defined(GFM_STELLAR_EVOLUTION) && ((GFM_OUTPUT_MASK)&4)
      init_field(IO_GFM_INITIAL_MASS, "GIMA", "GFM_InitialMass", MEM_MY_DOUBLE, FILE_DOUBLE, FILE_DOUBLE, 1, A_STARP,
                 &StarP[0].InitialMass, 0, STARS_ONLY);
      init_units(IO_GFM_INITIAL_MASS, 0., -1., 0., 1., 0., All.UnitMass_in_g);
      init_snapshot_type(IO_GFM_INITIAL_MASS, SN_MINI);
#endif

#ifdef GFM_STELLAR_EVOLUTION
      int stars_wind_bitmask = STARS_ONLY;
#ifdef GFM_WINDS_SAVE_PARTTYPE
      /* since wind particles have their own particle type, we can selectively include PT4 (star) fields */
      stars_wind_bitmask |= pow(2, GFM_WINDS_SAVE_PARTTYPE);
#endif

#if((GFM_OUTPUT_MASK)&8)
      init_field(IO_GFM_METALLICITY, "GZ  ", "GFM_Metallicity", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP,
                 &SphP[0].Metallicity, 0, GAS_ONLY);
      init_units(IO_GFM_METALLICITY, 0, 0, 0, 0, 0, 0); /* dimensionless ratio of masses */
      init_snapshot_type(IO_GFM_METALLICITY, SN_MINI);

      init_field(IO_GFM_METALLICITY2, "GZ2 ", "GFM_Metallicity", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_STARP,
                 &StarP[0].Metallicity, 0, stars_wind_bitmask);
      init_units(IO_GFM_METALLICITY2, 0, 0, 0, 0, 0, 0); /* dimensionless ratio of masses */
      init_snapshot_type(IO_GFM_METALLICITY2, SN_MINI);
#endif

#if((GFM_OUTPUT_MASK)&16) && !defined(GFM_STELLAR_EVOLUTION_NO_ELEMENTS)
      init_field(IO_GFM_METALS, "GMET", "GFM_Metals", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, GFM_N_CHEM_ELEMENTS, A_SPHP,
                 &SphP[0].MetalsFraction, 0, GAS_ONLY);
      init_units(IO_GFM_METALS, 0, 0, 0, 0, 0, 0); /* dimensionless abundance ratios */

      init_field(IO_GFM_METALS2, "GME2", "GFM_Metals", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, GFM_N_CHEM_ELEMENTS, A_STARP, 0,
                 io_func_gfm_metals, stars_wind_bitmask);
      init_units(IO_GFM_METALS2, 0, 0, 0, 0, 0, 0); /* dimensionless abundance ratios */
#endif

#if defined(GFM_CHEMTAGS) && ((GFM_OUTPUT_MASK)&256)
      init_field(IO_GFM_CHEM_TAGS, "GCTG", "GFM_MetalsTagged", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, GFM_N_CHEM_TAGS, A_SPHP,
                 0, io_func_gfm_chemtags_gas, GAS_ONLY);
      init_field(IO_GFM_CHEM_TAGS2, "GCT2", "GFM_MetalsTagged", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, GFM_N_CHEM_TAGS, A_STARP,
                 0, io_func_gfm_chemtags, STARS_ONLY);
#endif
#endif /* GFM_STELLAR_EVOLUTION */

#ifdef GFM_RPROCESS_CHANNELS
      init_field(IO_GFM_RPROCESS, "XRPR", "GFM_RProcess", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, GFM_RPROCESS_CHANNELS, A_NONE,
                 0, io_func_gfm_rprocess, GAS_ONLY | STARS_ONLY);
      init_field(IO_GFM_RPROCESS_COUNT, "NRPR", "GFM_NSNS_Count", MEM_INT, FILE_INT, FILE_NONE, GFM_RPROCESS_CHANNELS, A_STARP,
                 &StarP[0].NRProcessInjections[0], 0, STARS_ONLY);
#endif
#ifdef GFM_SNIA_ENERGY_INJECTION
      init_field(IO_GFM_RPROCESS_COUNT, "NRIA", "GFM_SNIa_Count", MEM_INT, FILE_INT, FILE_NONE, 1, A_STARP, &StarP[0].NumSNIa, 0,
                 STARS_ONLY);
#endif

#if defined(GFM_STELLAR_PHOTOMETRICS) && ((GFM_OUTPUT_MASK)&32)
      init_field(IO_GFM_STELLAR_PHOTOMETRICS, "GSPH", "GFM_StellarPhotometrics", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE,
                 GFM_STELLAR_PHOTOMETRICS_BANDS, A_NONE, 0, io_func_gfm_stellar_photometrics, STARS_ONLY);
      init_units(IO_GFM_STELLAR_PHOTOMETRICS, 0, 0, 0, 0, 0, 0); /* magnitudes */
#endif

#ifdef GFM_LAMBDA
      init_field(IO_GFM_LAMBDA, "GLAM", "GFM_Lambda", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 6, A_SPHP, &SphP[0].Lambda[0], 0,
                 GAS_ONLY);
      init_units(IO_GFM_LAMBDA, 0, 0, 0, 0, 0, 0);
#endif

#ifdef GFM_DUST
#ifdef GFM_DUST_CAP
      init_field(IO_GFM_DUST_CAPMASS, "GDCM", "GFM_DustCapMass", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP,
                 &SphP[0].DustMassCap, 0, GAS_ONLY);
      init_units(IO_GFM_DUST_CAPMASS, 0, 0, 0, 0, 0, 0);
#endif

      init_field(IO_GFM_DUST_TAUGROWTH, "TAGR", "GFM_DustTauGrowth", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP,
                 &SphP[0].DustTauGrowth, 0, GAS_ONLY);
      init_units(IO_GFM_DUST_TAUGROWTH, 0, 0, 0, 0, 0, 0);

#ifdef GFM_DUST_SPUTTERING
#if(GFM_DUST_SPUTTERING == 1)
      init_field(IO_GFM_DUST_TAUSPUTTER, "TASP", "GFM_DustTauSputter", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, 0,
                 io_func_gfm_dust_sputter_tau, GAS_ONLY);
      init_units(IO_GFM_DUST_TAUSPUTTER, 0, 0, 0, 0, 0, 0);
#endif
#endif
      init_field(IO_GFM_DUST_METALLICITY, "GDZ ", "GFM_DustMetallicity", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, 0,
                 io_func_gfm_dust_metal, GAS_ONLY);
      init_units(IO_GFM_DUST_METALLICITY, 0, 0, 0, 0, 0, 0); /* dimensionless abundance ratios */

      init_field(IO_GFM_DUST_METALLICITY2, "GDZ2", "GFM_DustMetallicity", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_STARP, 0,
                 io_func_gfm_dust_metal2, stars_wind_bitmask);
      init_units(IO_GFM_DUST_METALLICITY2, 0, 0, 0, 0, 0, 0); /* dimensionless abundance ratios */
#endif

#if defined(GFM_DUST) && ((GFM_OUTPUT_MASK)&128)
      init_field(IO_GFM_DUST_AGB, "GDAG", "GFM_DustAGB", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, GFM_N_CHEM_ELEMENTS, A_SPHP, 0,
                 io_func_gfm_dust_agb, GAS_ONLY);
      init_units(IO_GFM_DUST_AGB, 0, 0, 0, 0, 0, 0); /* dimensionless abundance ratios */
      init_field(IO_GFM_DUST_SNII, "GDII", "GFM_DustSNII", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, GFM_N_CHEM_ELEMENTS, A_SPHP,
                 0, io_func_gfm_dust_snii, GAS_ONLY);
      init_units(IO_GFM_DUST_SNII, 0, 0, 0, 0, 0, 0); /* dimensionless abundance ratios */
      init_field(IO_GFM_DUST_SNIa, "GDIa", "GFM_DustSNIa", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, GFM_N_CHEM_ELEMENTS, A_SPHP,
                 0, io_func_gfm_dust_snia, GAS_ONLY);
      init_units(IO_GFM_DUST_SNIa, 0, 0, 0, 0, 0, 0); /* dimensionless abundance ratios */

      init_field(IO_GFM_DUST_AGB2, "GDA2", "GFM_DustAGB", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, GFM_N_CHEM_ELEMENTS, A_STARP,
                 0, io_func_gfm_dust_agb2, stars_wind_bitmask);
      init_units(IO_GFM_DUST_AGB2, 0, 0, 0, 0, 0, 0); /* dimensionless abundance ratios */
      init_field(IO_GFM_DUST_SNII2, "GDI2", "GFM_DustSNII", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, GFM_N_CHEM_ELEMENTS, A_STARP,
                 0, io_func_gfm_dust_snii2, stars_wind_bitmask);
      init_units(IO_GFM_DUST_SNII2, 0, 0, 0, 0, 0, 0); /* dimensionless abundance ratios */
      init_field(IO_GFM_DUST_SNIa2, "GDJ2", "GFM_DustSNIa", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, GFM_N_CHEM_ELEMENTS, A_STARP,
                 0, io_func_gfm_dust_snia2, stars_wind_bitmask);
      init_units(IO_GFM_DUST_SNIa2, 0, 0, 0, 0, 0, 0); /* dimensionless abundance ratios */
#endif
    }

#if defined(SMUGGLE_STAR_FEEDBACK) || defined(GFM_OUTPUT_BIRTH_POS)
  init_field(IO_GFM_BIRTH_POS, "BIRP", "BirthPos", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_NONE, 3, A_STARP, &StarP[0].BirthPos[0], 0,
             STARS_ONLY);
  init_units(IO_GFM_BIRTH_POS, 1., -1., 1., 0., 0., All.UnitLength_in_cm);

  init_field(IO_GFM_BIRTH_VEL, "BIRV", "BirthVel", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_NONE, 3, A_STARP, &StarP[0].BirthVel[0], 0,
             STARS_ONLY);
  init_units(IO_GFM_BIRTH_VEL, 0.5, 0., 0., 0., 1., All.UnitVelocity_in_cm_per_s);

  init_field(IO_GFM_BIRTH_RHO, "BIRD", "BirthDensity", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_STARP, &StarP[0].BirthDensity,
             0, STARS_ONLY);
  init_units(IO_GFM_BIRTH_RHO, -3., 2., -3., 1., 0., All.UnitDensity_in_cgs);
#endif

  /* GAS/STARS (NON-GFM MODELS) */

#ifdef SMUGGLE_OUTPUT_STELLAR_FEEDBACK
  init_field(IO_GFM_SNII_NUM, "SNII", "SNIINumber", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_STARP, &StarP[0].SNII_Num, 0,
             STARS_ONLY);
  init_field(IO_GFM_CUM_SNII_NUM, "CSII", "CumSNIINumber", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_STARP,
             &StarP[0].Cum_SNII_Num, 0, STARS_ONLY);
  init_field(IO_GFM_SNIa_NUM, "SNIa", "SNIaNumber", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_STARP, &StarP[0].SNIa_Num, 0,
             STARS_ONLY);
  init_field(IO_GFM_CUM_SNIa_NUM, "CSIa", "CumSNIaNumber", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_STARP,
             &StarP[0].Cum_SNIa_Num, 0, STARS_ONLY);
  init_field(IO_GFM_FEED_ENERGY, "FERG", "FeedbackEnergy", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_STARP,
             &StarP[0].FeedbackEnergy, 0, STARS_ONLY);
  init_field(IO_GFM_FEED_MOMENTUM, "FMOM", "FeedbackMomentum", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_STARP,
             &StarP[0].FeedbackMomentum, 0, STARS_ONLY);
  init_field(IO_GFM_FEED_MOMENTUM_AGB, "FAGB", "FeedbackMomentumAGB", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_STARP,
             &StarP[0].FeedbackMomentumAGB, 0, STARS_ONLY);
  init_field(IO_GFM_CUM_FEED_MOMENTUM, "CMOM", "CumFeedbackMomentum", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_STARP,
             &StarP[0].Cum_FeedbackMomentum, 0, STARS_ONLY);
  init_field(IO_GFM_CUM_INJ_MOMENTUM, "CMIN", "CumInjFeedbackMomentum", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_STARP,
             &StarP[0].Cum_InjFeedbackMomentum, 0, STARS_ONLY);
  init_field(IO_GFM_CUM_INJ_MOMENTUM_AGB, "CAGB", "CumInjFeedbackMomentumAGB", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_STARP,
             &StarP[0].Cum_InjFeedbackMomentumAGB, 0, STARS_ONLY);
  init_field(IO_GFM_MASS_RELEASED, "MREL", "MassReleased", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_STARP,
             &StarP[0].TotalMassReleased, 0, STARS_ONLY);
  init_field(IO_HSML_STARS, "SSML", "StarsHsml", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_STARP, &StarP[0].Hsml, 0, STARS_ONLY);
  init_field(IO_MAX_RAD_STARS, "SMRA", "StarsMaxRadius", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_STARP,
             &StarP[0].MaxFeedRadius, 0, STARS_ONLY);
  init_field(IO_LOC_DENS_STARS, "LISM", "LocISMDens", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_STARP, &StarP[0].LocISMdens, 0,
             STARS_ONLY);
#endif

#ifdef SMUGGLE_RADIATION_FEEDBACK
  init_field(IO_STROMGREN_RADIUS, "STRM", "StromgrenRadius", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_STARP,
             &StarP[0].StromgrenRadius, 0, STARS_ONLY);
  init_field(IO_TAU_RADFEED, "TAUD", "RadFeedTau", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_STARP, &StarP[0].RadFeedTau, 0,
             STARS_ONLY);
  init_field(IO_NUMNGB_RADFEED, "RFNN", "RadFeed_NumNgb", MEM_INT, FILE_INT, FILE_NONE, 1, A_STARP, &StarP[0].RadFeed_NumNgb, 0,
             STARS_ONLY);
#endif

#ifdef SMUGGLE_RADIATION_FEEDBACK_DEBUG
  init_field(IO_COOL_DELAY_RADFEED, "RFCD", "RadCoolShutoffTime", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_STARP,
             &StarP[0].RadCoolShutoffTime, 0, STARS_ONLY);
  init_field(IO_COOL_DELAY_GAS_RADFEED, "RFCG", "GasRadCoolShutoffTime", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP,
             &SphP[0].GasRadCoolShutoffTime, 0, GAS_ONLY);
  init_field(IO_RADIATION_MOMENTUM, "RADE", "RadiationMomentumReleased", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_STARP,
             &StarP[0].RadiationMomentumReleased, 0, STARS_ONLY);
  init_field(IO_CUM_RADIATION_MOMENTUM, "CRMO", "Cum_RadiationMomentumReleased", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 1,
             A_STARP, &StarP[0].Cum_RadiationMomentumReleased, 0, STARS_ONLY);
  init_field(IO_CUM_RADIATION_MOMENTUM_REAL, "CRMI", "Cum_RadMomentumRealInjected", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 1,
             A_STARP, &StarP[0].Cum_RadMomentumRealInjected, 0, STARS_ONLY);
  init_field(IO_NORMSPH_RADFEED, "RFNS", "NormSphRadFeedback", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_STARP,
             &StarP[0].NormSphRadFeedback, 0, STARS_ONLY);
  init_field(IO_COLUMN_DENSITY, "DCOL", "ColumnDensity", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_STARP,
             &StarP[0].GasColumnDensity, 0, STARS_ONLY);
#endif

#ifdef SMUGGLE_VAR_SN_EFF
  init_field(IO_VAREFF_SN_METAL, "VSNZ", "AvgMetalNgb", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_STARP, &StarP[0].AvgMetalNgb,
             0, STARS_ONLY);
#endif

#ifdef SMUGGLE_MASS_WEIGHT_SN
  init_field(IO_MASS_WEIGHT_SN, "NGBM", "TotNgbMass", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_STARP, &StarP[0].TotNgbMass, 0,
             STARS_ONLY);
#endif

#ifdef SMUGGLE_OUTPUT_OPTICAL_DEPTH /* SMUGGLE_SFR */
  init_field(IO_COOL_DELAY_GAS_RADFEED, "TAU ", "OpticalDepth", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP,
             &SphP[0].OpticalDepth, 0, GAS_ONLY);
#endif
#ifdef SMUGGLE_OUTPUT_VIRIAL_PARAM /* SMUGGLE_SFR */
  init_field(IO_VIRIAL_PARAM, "ALPH", "VirialParameter", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].VirialParam, 0,
             GAS_ONLY);
#endif

#if defined(SMUGGLE_RADPRESS_OPT_THIN) || defined(SMUGGLE_RADPRESS_OPT_THICK)
  init_field(IO_RP_MOM, "RADP", "RadiationPressureMoment", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_SPHP,
             &SphP[0].RadPress[0], 0, GAS_ONLY);
#endif

#ifdef SFR_MCS
  init_field(IO_MCS_BIRTH_TIME, "BIRT", "StellarFormationTime", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_STARP, 0,
             io_func_stellar_birthtime, STARS_ONLY);
  init_units(IO_MCS_BIRTH_TIME, 0, 0, 0, 0, 0, 0);
  init_field(IO_MCS_AGE, "SAGE", "StellarAgeGyr", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_STARP, 0,
             io_func_stellar_age_in_Gyr, STARS_ONLY);
  init_units(IO_MCS_AGE, 0, 0, 1, 0, -1, SEC_PER_GIGAYEAR);

#if SFR_MCS_RATE_CRITERIA > 0
  init_field(IO_MCS_GRADV_SQ, "GRV2", "GradVelSquared", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].gradv_sq,
             0, GAS_ONLY);
#endif

#ifdef SFR_MCS_BIRTH_RECORDS
  init_field(IO_MCS_BIRTH_POS, "BIRP", "BirthPos", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_STARP,
             &StarP[0].BirthPos[0], 0, STARS_ONLY);
  init_units(IO_MCS_BIRTH_POS, 1., -1., 1., 0., 0., All.UnitLength_in_cm);
  init_field(IO_MCS_BIRTH_VEL, "BIRV", "BirthVel", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_STARP,
             &StarP[0].BirthVel[0], 0, STARS_ONLY);
  init_units(IO_MCS_BIRTH_VEL, 0.5, 0., 0., 0., 1., All.UnitVelocity_in_cm_per_s);
  init_field(IO_MCS_BIRTH_RHO, "BIRD", "BirthDensity", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_STARP,
             &StarP[0].BirthDensity, 0, STARS_ONLY);
  init_units(IO_MCS_BIRTH_RHO, -3., 2., -3., 1., 0., All.UnitDensity_in_cgs);
#endif

#ifdef IMF_SAMPLING_MCS
  init_field(IO_IMF_MCS_MASSES, "IMFM", "StellarMassArray", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, N_STAR_SLOTS, A_STARP,
             &StarP[0].MassArr[0], 0, STARS_ONLY);
  init_units(IO_IMF_MCS_MASSES, 0., -1., 0., 1., 0., All.UnitMass_in_g);
  init_field(IO_IMF_MCS_LIFETIMES, "IMFL", "StellarLifetimeArray", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, N_STAR_SLOTS,
             A_STARP, &StarP[0].LifetimeArr[0], 0, STARS_ONLY);
  init_units(IO_IMF_MCS_LIFETIMES, 0., 0., 1., 0., -1., SEC_PER_YEAR);
#ifdef HII_MCS
  init_field(IO_IMF_MCS_PHOTON_RATES, "IMFI", "StellarIonisingPhotonRate1e49Array", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT,
             N_STAR_SLOTS, A_STARP, &StarP[0].S_HiiArr[0], 0, STARS_ONLY);
  init_units(IO_IMF_MCS_PHOTON_RATES, 0., 0., -1, 0, 1., 1e49);
#ifdef HII_MCS_LR
  init_field(IO_IMF_MCS_PHOTON_ENERGY, "IMFE", "StellarIonisingEnergyPerPhotonArray", MEM_MY_SINGLE, FILE_MY_IO_FLOAT,
             FILE_MY_IO_FLOAT, N_STAR_SLOTS, A_STARP, &StarP[0].EnergyPerPhotonArr[0], 0, STARS_ONLY);
#endif
#endif
#ifdef PE_MCS
  init_field(IO_IMF_MCS_LFUVS, "IMFF", "StellarLuminosityFUVArray", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, N_STAR_SLOTS,
             A_STARP, &StarP[0].L_FUVArr[0], 0, STARS_ONLY);
#endif
#endif  // IMF_SAMPLING_MCS

#ifdef SN_MCS
  init_field(IO_MCS_N_SN, "NSN ", "NumberOfSupernovae", MEM_INT, FILE_INT, FILE_INT, 1, A_STARP, &StarP[0].N_SN_cum, 0, STARS_ONLY);
  init_units(IO_MCS_N_SN, 0, 0, 0, 0, 0, 0);
  init_field(IO_MCS_N_SN_EVENT, "NSNE", "NumberOfSupernovaEvents", MEM_INT, FILE_INT, FILE_INT, 1, A_STARP, &StarP[0].N_SN_event_cum,
             0, STARS_ONLY);
  init_units(IO_MCS_N_SN_EVENT, 0, 0, 0, 0, 0, 0);
  init_field(IO_MCS_INITIAL_MASS, "IMAS", "InitialMass", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_STARP,
             &StarP[0].InitialMass, 0, STARS_ONLY);
  init_units(IO_MCS_INITIAL_MASS, 0., -1., 0., 1., 0., All.UnitMass_in_g);
#ifdef SN_MCS_LOCATION_RECORDS
  init_field(IO_MCS_SN_TIME, "SNT ", "SupernovaTime", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_STARP, &StarP[0].SNTime,
             0, STARS_ONLY);
  init_units(IO_MCS_SN_TIME, 0, 0, 0, 0, 0, 0);
  init_field(IO_MCS_SN_POS, "SNP ", "SupernovaPos", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_STARP, &StarP[0].SNPos[0],
             0, STARS_ONLY);
  init_units(IO_MCS_SN_POS, 1., -1., 1., 0., 0., All.UnitLength_in_cm);
  init_field(IO_MCS_SN_RHO, "SNLD", "SupernovaLocalDensity", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_STARP,
             &StarP[0].SNDensity, 0, STARS_ONLY);
  init_units(IO_MCS_SN_RHO, -3., 2., -3., 1., 0., All.UnitDensity_in_cgs);
#ifdef GRACKLE
  init_field(IO_MCS_SN_TEMP, "SNLT", "SupernovaLocalTemperature", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_STARP,
             &StarP[0].SNTemperature, 0, STARS_ONLY);
  init_units(IO_MCS_SN_TEMP, 0, 0, 0, 0, 0, 0);
#endif  // GRACKLE
#endif  // SN_MCS_LOCATION_RECORDS
#endif  // SN_MCS

#ifdef HII_MCS
  init_field(IO_MCS_STROMGREN_SOURCE_ID, "STID", "StromgrenSourceID", MEM_MY_ID_TYPE, FILE_MY_ID_TYPE, FILE_NONE, 1, A_SPHP,
             &SphP[0].StromgrenSourceID, 0, GAS_ONLY);
  init_field(IO_MCS_R_STROMGREN, "STRG", "StromgrenRadius", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP,
             &SphP[0].R_Stromgren, 0, GAS_ONLY);
  init_field(IO_MCS_PHOTON_RATE, "STRR", "IonisingPhotonRate1e49", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_STARP,
             &StarP[0].S_Hii, 0, STARS_ONLY);
#ifdef HII_MCS_LR
  init_field(IO_MCS_EDENS_HII, "EDHI", "EnergyDensityHii", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP,
             &SphP[0].EnergyDensHii, 0, GAS_ONLY);
#endif
#ifdef HII_MCS_RECORDS
#ifdef HII_MCS_ANISO
  init_field(IO_MCS_R_STROMGREN_ARR, "STRA", "StromgrenRadiusArray", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, HII_MCS_N_PIX,
             A_SPHP, &SphP[0].R_StromgrenArr[0], 0, GAS_ONLY);
#endif
  init_field(IO_MCS_N_HII_SOURCE, "NHIS", "NumberOfIonisingStars", MEM_INT, FILE_INT, FILE_INT, 1, A_SPHP, &SphP[0].N_Sources, 0,
             GAS_ONLY);
  init_field(IO_MCS_HII_SOURCE_IDS, "HSID", "IonisingStarIDs", MEM_MY_ID_TYPE, FILE_MY_ID_TYPE, FILE_MY_ID_TYPE, 3, A_SPHP,
             &SphP[0].StarIDArr, 0, GAS_ONLY);
  init_field(IO_MCS_HII_RECOM_RATE, "HIRR", "HiiRecombinationRate1e49", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP,
             &SphP[0].HiiRecombinationRate, 0, GAS_ONLY);
#endif  // HII_MCS_RECORDS
#endif  // HII_MCS

#ifdef PE_MCS
  init_field(IO_MCS_GFUV, "GFUV", "G_Habing", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].G_FUV, 0,
             GAS_ONLY);
  init_field(IO_MCS_PEHR, "PEHR", "PhotoelectricHeatingRate", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, 0,
             io_func_pe_heating, GAS_ONLY);
  init_field(IO_MCS_LUFV, "LUFV", "LuminosityFUV", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_STARP, &StarP[0].L_FUV, 0,
             STARS_ONLY);
#ifdef PE_MCS_STORE_DUST_COLUMN
  init_field(IO_MCS_DCOL, "DCOL", "DustColumnDensity", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_STARP,
             &StarP[0].DustColumn, 0, STARS_ONLY);
#endif
#endif  // PE_MCS
#endif  // SFR_MCS

#ifdef TURB_APPROX_MCS
  init_field(IO_MCS_TURBSPEC, "TSPC", "TurbulentSpecEnergy", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP,
             &SphP[0].TurbSpecEnergy, 0, GAS_ONLY);
#ifdef TURB_APPROX_MCS_OUTPUT_RATES
  init_field(IO_MCS_DVEL_UNLIM, "DVUN", "VelocityGradientUnlimited", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 9, A_SPHP,
             &SphP[0].dvel_unlim[0][0], 0, GAS_ONLY);
  init_field(IO_MCS_TURB_PROD, "TPRD", "TurbulentProductionRate", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, 0,
             io_func_turbulent_production_rate, GAS_ONLY);
  init_field(IO_MCS_TURB_DISS, "TDIS", "TurbulentDissipationRate", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, 0,
             io_func_turbulent_dissipation_rate, GAS_ONLY);
  init_field(IO_MCS_TURB_ADIA, "TADI", "TurbulentAdiabaticRate", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, 0,
             io_func_turbulent_adiabatic_rate, GAS_ONLY);
#endif
#endif

  /* BLACK HOLES */

#if defined(BLACK_HOLES)
  init_field(IO_BHMASS, "BHMA", "BH_Mass", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP, &BHP[0].BH_Mass, 0, BHS_ONLY);
  init_units(IO_BHMASS, 0., -1., 0., 1., 0., All.UnitMass_in_g);

  init_field(IO_BHMDOT, "BHMD", "BH_Mdot", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP, &BHP[0].BH_Mdot, 0, BHS_ONLY);
  init_units(IO_BHMDOT, 0.0, 0.0, 1.0, -1.0, 1.0,
             All.UnitMass_in_g / All.UnitLength_in_cm * All.UnitVelocity_in_cm_per_s); /* mass/time */

  init_field(IO_BHHSML, "BHHS", "BH_Hsml", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP, &BHP[0].BH_Hsml, 0, BHS_ONLY);
  init_units(IO_BHHSML, 1., -1., 1., 0., 0., All.UnitLength_in_cm);

#ifdef MRT_BH
  init_field(IO_BHPHSML, "BPHS", "BH_PhotonHsml", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP, &BHP[0].BH_PhotonHsml, 0,
             BHS_ONLY);
  init_units(IO_BHPHSML, 1., -1., 1., 0., 0., All.UnitLength_in_cm);
#endif

  init_field(IO_BHPRESS, "BHPR", "BH_Pressure", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP, &BHP[0].BH_Pressure, 0,
             BHS_ONLY);
  init_units(
      IO_BHPRESS, 1.0, 2.0, 1.0, 1.0, 2.0,
      All.UnitMass_in_g * All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s * All.UnitLength_in_cm); /* mass/length/time^2 */

  init_field(IO_BHU, "BHU ", "BH_U", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP, &BHP[0].BH_U, 0, BHS_ONLY);
  init_units(IO_BHU, 0., 0., 0., 0., 2., All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s);

  init_field(IO_BHRHO, "BHRO", "BH_Density", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP, &BHP[0].BH_Density, 0,
             BHS_ONLY);
  init_units(IO_BHRHO, -3., 2., -3., 1., 0., All.UnitDensity_in_cgs);

  init_field(IO_BHPROGS, "BHPR", "BH_Progs", MEM_INT, FILE_INT, FILE_INT, 1, A_BHP, &BHP[0].BH_CountProgs, 0, BHS_ONLY);
  init_units(IO_BHPROGS, 0, 0, 0, 0, 0, 0); /* dimensionless counter */

  init_field(IO_BHCMQM, "BCMQ", "BH_CumMassGrowth_QM", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP,
             &BHP[0].BH_CumMass_QM, 0, BHS_ONLY);
  init_field(IO_BHCMRM, "BCMR", "BH_CumMassGrowth_RM", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP,
             &BHP[0].BH_CumMass_RM, 0, BHS_ONLY);
  init_units(IO_BHCMQM, 0., -1., 0., 1., 0., All.UnitMass_in_g);
  init_units(IO_BHCMRM, 0., -1., 0., 1., 0., All.UnitMass_in_g);

  init_field(IO_BHCERM, "BCER", "BH_CumEgyInjection_RM", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP,
             &BHP[0].BH_CumEgy_RM, 0, BHS_ONLY);
  init_field(IO_BHCEQM, "BCEQ", "BH_CumEgyInjection_QM", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP,
             &BHP[0].BH_CumEgy_QM, 0, BHS_ONLY);
  init_field(IO_BHEGYL, "BCCL", "BH_MPB_CumEgyLow", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP,
             &BHP[0].BH_MPB_CumEgyLow, 0, BHS_ONLY);
  init_field(IO_BHEGYH, "BCCH", "BH_MPB_CumEgyHigh", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP,
             &BHP[0].BH_MPB_CumEgyHigh, 0, BHS_ONLY);
  init_units(IO_BHCERM, 2.0, -1.0, 0.0, 1.0, 2.0,
             All.UnitMass_in_g * All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s); /* mass*length^2/time^2 */
  init_units(IO_BHCEQM, 2.0, -1.0, 0.0, 1.0, 2.0, All.UnitMass_in_g * All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s);
  init_units(IO_BHEGYL, 2.0, -1.0, 0.0, 1.0, 2.0, All.UnitMass_in_g * All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s);
  init_units(IO_BHEGYH, 2.0, -1.0, 0.0, 1.0, 2.0, All.UnitMass_in_g * All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s);

  init_field(IO_BHMDOTBONDI, "BHBO", "BH_MdotBondi", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP, &BHP[0].BH_MdotBondi,
             0, BHS_ONLY);
  init_units(IO_BHMDOTBONDI, 0.0, 0.0, 1.0, -1.0, 1.0,
             All.UnitMass_in_g / All.UnitLength_in_cm * All.UnitVelocity_in_cm_per_s); /* mass/time */

  init_field(IO_BHMDOTEDDIN, "BHED", "BH_MdotEddington", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP,
             &BHP[0].BH_MdotEddington, 0, BHS_ONLY);
  init_units(IO_BHMDOTEDDIN, 0.0, 0.0, 1.0, -1.0, 1.0,
             All.UnitMass_in_g / All.UnitLength_in_cm * All.UnitVelocity_in_cm_per_s); /* mass/time */
#endif

#ifdef BH_BUBBLES
  init_field(IO_BHMBUB, "BHMB", "BH_Mass_bubbles", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP, &BHP[0].BH_Mass_bubbles,
             0, BHS_ONLY);
  init_units(IO_BHMBUB, 0., -1., 0., 1., 0., All.UnitMass_in_g);

  init_field(IO_BHMINI, "BHMI", "BH_Mass_ini", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP, &BHP[0].BH_Mass_ini, 0,
             BHS_ONLY);
  init_units(IO_BHMINI, 0., -1., 0., 1., 0., All.UnitMass_in_g);
#endif

#ifdef BH_NF_RADIO
  init_field(IO_BHMDOTQUASAR, "BHMQ", "BH_Mdot_Quasar", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP,
             &BHP[0].BH_Mdot_quasar, 0, BHS_ONLY);
  init_field(IO_BHMDOTRADIO, "BHMR", "BH_Mdot_Radio", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP,
             &BHP[0].BH_Mdot_radio, 0, BHS_ONLY);
  init_field(IO_BHVVIR, "BHVV", "BH_HaloVvir", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP, &BHP[0].BH_HaloVvir, 0,
             BHS_ONLY);
  init_field(IO_BHXRAYLUM, "BXRY", "BH_XrayLuminosity", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP, &BHP[0].BH_XrayLum,
             0, BHS_ONLY);
  init_field(IO_BHRADIOLUM, "BRAL", "BH_RadioLuminosity", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP,
             &BHP[0].BH_RadioLum, 0, BHS_ONLY);
#endif

#if defined(BLACK_HOLES) && defined(GFM_AGN_RADIATION)
  init_field(IO_BHHOSTHALOMASS, "HHA2", "BH_HostHaloMass", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP,
             &BHP[0].HostHaloMass, 0, BHS_ONLY);
#endif

#if defined(BLACK_HOLES) && defined(BH_FRICTION)
  init_field(IO_BH_FRC_MINPOT, "BHMP", "BH_MinPot", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP, &BHP[0].BH_MinPot, 0,
             BHS_ONLY);
  init_field(IO_BH_FRC_POTPOS, "BHPP", "BH_MinPotPos", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_BHP,
             &BHP[0].BH_MinPotPos[0], 0, BHS_ONLY);
  init_field(IO_BH_FRC_MINPOT_EXT, "BHME", "BH_MinPotExt", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP,
             &BHP[0].BH_MinPot_Extended, 0, BHS_ONLY);
  init_field(IO_BH_FRC_POTPOS_EXT, "BHPE", "BH_MinPotPosExt", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_BHP,
             &BHP[0].BH_MinPotPos_Extended[0], 0, BHS_ONLY);
  init_field(IO_BH_FRC_POTVEL, "BHPV", "BH_MinPotVel", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_BHP,
             &BHP[0].BH_MinPotVel[0], 0, BHS_ONLY);
  init_field(IO_BH_FRC_RHOTOT, "BHRT", "BH_RhoTot", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP, &BHP[0].BH_RhoTot, 0,
             BHS_ONLY);
#endif

#if defined(BLACK_HOLES) && defined(BH_SPIN_EVOLUTION)
  init_field(IO_BHANGMOMGASCELLS, "BHAMG", "BH_AngMomGasCells", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_BHP,
             &BHP[0].BH_AngMomGasCells[0], 0, BHS_ONLY);
  init_units(IO_BHANGMOMGASCELLS, 0., 0., 0., 0., 0., 1.);
  init_field(IO_BHSPINMODEL, "BHSM", "BH_SpinModel", MEM_INT, FILE_INT, FILE_INT, 1, A_BHP, &BHP[0].BH_SpinModel, 0, BHS_ONLY);
  init_units(IO_BHSPINMODEL, 0., 0., 0., 0., 0., 1.);
  init_field(IO_BHSPINPARAMETER, "BHSP", "BH_SpinParameter", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP,
             &BHP[0].BH_SpinParameter, 0, BHS_ONLY);
  init_units(IO_BHSPINPARAMETER, 0., 0., 0., 0., 0., 1.);
  init_field(IO_BHSPINORIENTATION, "BHSO", "BH_SpinOrientation", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_BHP,
             &BHP[0].BH_SpinOrientation[0], 0, BHS_ONLY);
  init_units(IO_BHSPINORIENTATION, 0., 0., 0., 0., 0., 1.);
#endif

#ifdef OUTPUT_BLACK_HOLE_TIMESTEP
  init_field(IO_BHTIMESTEP, "BHTS", "BH_TimeStep", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, 3, A_NONE, 0, io_func_bh_timestep, BHS_ONLY);
#endif

#ifdef BH_USE_ALFVEN_SPEED_IN_BONDI
  init_field(IO_BHBPRESS, "BHBP", "BH_BPressure", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_BHP, &BHP[0].BH_Bpress, 0,
             BHS_ONLY);
#endif

  /* LIVE DUST */

#ifdef DUST_LIVE
  init_field(IO_DUST_GASHSML, "DLHS", "Dust_GasHsml", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_DUSTP, &DustP[0].Hsml, 0,
             DUST_ONLY);
  init_units(IO_DUST_GASHSML, 1., -1., 1., 0., 0., All.UnitLength_in_cm);

  init_field(IO_DUST_GASDENSITY, "DGDE", "Dust_GasDensity", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_DUSTP,
             &DustP[0].LocalGasDensity, 0, DUST_ONLY);
  init_units(IO_DUST_GASDENSITY, -3., 2., -3., 1., 0., All.UnitDensity_in_cgs);

#ifdef DL_GRAIN_BINS
  init_field(IO_DUST_GSD, "DGSD", "Dust_NumGrains", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, DL_GRAIN_BINS, A_DUSTP,
             &DustP[0].NumGrains, 0, DUST_ONLY);
  init_units(IO_DUST_GSD, 0, 0, 0, 0, 0, 0);

#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
  init_field(IO_DUST_GSD_SLOPES, "DGBS", "Dust_BinSlopes", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, DL_GRAIN_BINS, A_DUSTP,
             &DustP[0].BinSlopes, 0, DUST_ONLY);
  init_units(IO_DUST_GSD_SLOPES, 0, 0, 0, 0, 0, 0);
#endif

  init_field(IO_DUST_METALFRACS, "DMFR", "Dust_MetalFractions", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, GFM_N_CHEM_ELEMENTS,
             A_DUSTP, &DustP[0].MetalFractions, 0, DUST_ONLY);
  init_units(IO_DUST_METALFRACS, 0, 0, 0, 0, 0, 0);
#endif

#if(defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION)) || (defined(DL_DRAG_BACKREACTION))
  init_field(IO_DUST_DENSITY, "DDEN", "Dust_DustDensity", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_DUSTP,
             &DustP[0].DustDensity, 0, DUST_ONLY);
  init_units(IO_DUST_DENSITY, -3., 2., -3., 1., 0., All.UnitDensity_in_cgs);

  init_field(IO_DUST_DUSTHSML, "DLHD", "Dust_DustHsml", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_DUSTP,
             &DustP[0].DustHsml, 0, DUST_ONLY);
  init_units(IO_DUST_DUSTHSML, 1., -1., 1., 0., 0., All.UnitLength_in_cm);

  init_field(IO_CLOUD_FRACTION, "CDFR", "Cloud_Fraction", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_DUSTP,
             &DustP[0].LocalCloudFrac, 0, DUST_ONLY);
#endif

#ifdef DL_DEREFINEMENT
  init_field(IO_DUST_NGBID, "DNGB", "Dust_NgbID", MEM_MY_ID_TYPE, FILE_MY_ID_TYPE, FILE_MY_ID_TYPE, 1, A_DUSTP,
             &DustP[0].ClosestDustID, 0, DUST_ONLY);
#endif

#ifdef DL_WINDS
  init_field(IO_DUST_ISWIND, "DISW", "Dust_IsWind", MEM_INT, FILE_INT, FILE_INT, 1, A_DUSTP, &DustP[0].IsWind, 0, DUST_ONLY);
  init_field(IO_DUST_WINDTIME, "DWTM", "Dust_WindTimeLeft", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_DUSTP,
             &DustP[0].WindTimeLeft, 0, DUST_ONLY);
#endif

#if defined(DL_RADIATION) && defined(DL_OUTPUT_RT_FLUX)
  init_field(IO_DUST_RTFX, "DRFX", "Dust_RTFluxX", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, MRT_BINS, A_DUSTP, 0, io_func_dust_mrt_rtfx,
             DUST_ONLY);
  init_field(IO_DUST_RTFY, "DRFY", "Dust_RTFluxY", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, MRT_BINS, A_DUSTP, 0, io_func_dust_mrt_rtfy,
             DUST_ONLY);
  init_field(IO_DUST_RTFZ, "DRFZ", "Dust_RTFluxZ", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, MRT_BINS, A_DUSTP, 0, io_func_dust_mrt_rtfz,
             DUST_ONLY);
#endif

#if defined(DL_RADIATION) && defined(DL_THERMAL_IR)
#ifdef DL_GRAIN_BINS
  init_field(IO_DUST_TEMPERATURE, "DTEM", "Dust_Temperature", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, DL_GRAIN_BINS, A_DUSTP,
             &DustP[0].DustTemp, 0, DUST_ONLY);
#else
  init_field(IO_DUST_TEMPERATURE, "DTEM", "Dust_Temperature", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_DUSTP,
             &DustP[0].DustTemp, 0, DUST_ONLY);
#endif
#endif
#endif

  /* DIAGNOSTIC */

#ifdef OUTPUT_TASK
  init_field(IO_TASK, "TASK", "task", MEM_INT, FILE_INT, FILE_NONE, 1, A_NONE, 0, io_func_task, GAS_ONLY);
#endif

#ifdef OUTPUT_TIMEBIN_HYDRO
  init_field(IO_TIMEBIN_HYDRO, "TBH", "TimebinHydro", MEM_NONE, FILE_INT, FILE_NONE, 1, A_NONE, 0, io_func_timebin_hydro, GAS_ONLY);
#endif

#ifdef OUTPUTTIMESTEP
  init_field(IO_TSTP, "TSTP", "TimeStep", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_NONE, 0, io_func_timestep, ALL_TYPES);
#endif

#ifdef OUTPUTACCELERATION
  init_field(IO_ACCEL, "ACCE", "Acceleration", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_NONE, 0, io_func_accel, ALL_TYPES);
#endif

#ifdef OUTPUT_SOFTENINGS
  init_field(IO_SOFTENING, "SOFT", "Softenings", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_NONE, 0, io_func_softenings, ALL_TYPES);
#endif

#ifdef OUTPUTGRAVINTERACTIONS
  init_field(IO_GRAVITERACTIONS, "GINT", "GravityInteractions", MEM_INT, FILE_INT, FILE_NONE, 1, A_SPHP, &SphP[0].GravInteractions, 0,
             ALL_TYPES);
#endif

  /* DG */

#ifdef DG
  /* FIXME add cosmological factors to DG output (and everywhere else, Kevin) */
  init_field(IO_DG_W0, "dgw0", "dgw0", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, NOF_BASE_FUNCTIONS, A_NONE, 0, io_func_dgw0,
             1);
  init_units(IO_DG_W0, 0., 0., -3, 1., 0., All.UnitDensity_in_cgs);
  init_field(IO_DG_W1, "dgw1", "dgw1", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, NOF_BASE_FUNCTIONS, A_NONE, 0, io_func_dgw1,
             1);
  init_units(IO_DG_W1, 0., 0., -3, 1., 1., All.UnitDensity_in_cgs * All.UnitVelocity_in_cm_per_s);
  init_field(IO_DG_W2, "dgw2", "dgw2", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, NOF_BASE_FUNCTIONS, A_NONE, 0, io_func_dgw2,
             1);
  init_units(IO_DG_W2, 0., 0., -3, 1., 1., All.UnitDensity_in_cgs * All.UnitVelocity_in_cm_per_s);
  init_field(IO_DG_W3, "dgw3", "dgw3", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, NOF_BASE_FUNCTIONS, A_NONE, 0, io_func_dgw3,
             1);
  init_units(IO_DG_W3, 0., 0., -3, 1., 1., All.UnitDensity_in_cgs * All.UnitVelocity_in_cm_per_s);
  init_field(IO_DG_W4, "dgw4", "dgw4", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, NOF_BASE_FUNCTIONS, A_NONE, 0, io_func_dgw4,
             1);
  init_units(IO_DG_W4, 0., 0., -3, 1., 2., All.UnitDensity_in_cgs * All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s);
#endif

#ifdef OUTPUT_DG_ACCELERATION
  init_field(IO_DG_ACCEL, "dgac", "dg_accel", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 3, A_NONE, 0, io_func_dg_accel, 1);
#endif

#ifdef OUTPUT_DG_ANGULAR_MOMENTUM
  init_field(IO_DG_AM, "dgam", "dg_amomentum", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_NONE, 0, io_func_dg_angular_momentum, 1);
#endif

#ifdef OUTPUT_DG_SPIN
  init_field(IO_DG_SPIN, "dgsp", "dg_spin", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_NONE, 0, io_func_dg_spin, 1);
#endif

#ifdef OUTPUT_DG_TIMESTEP
  init_field(IO_DG_DT, "dgdt", "dg_timestep", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_NONE, 0, io_func_dg_timestep, 1);
#endif

#ifdef OUTPUT_DG_L1_NORM
  init_field(IO_DG_L1, "dgno", "dg_norm", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_NONE, 0, io_func_dg_norm, 1);
#endif

#ifdef OUTPUT_DG_TEMPERATURE
  init_field(IO_DG_T, "dgte", "dgw_temperature", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, NOF_BASE_FUNCTIONS, A_NONE, 0,
             io_func_dgw_temperature, 1);
#endif

#ifdef OUTPUT_DG_U
  init_field(IO_DG_U, "dgut", "dgw_u", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, NOF_BASE_FUNCTIONS, A_NONE, 0, io_func_dgw_u, 1);
#endif

#ifdef OUTPUT_DG_INFLOW_BOUNDARIES
  init_field(IO_DG_IB, "dgib", "dg_inflow", MEM_INT, FILE_INT, FILE_INT, 1, 0, A_SPHP, &SphP[0].Inflow_boundaries, 1);
#endif

#ifdef OUTPUT_DG_DISCONTINUITIES
  init_field(IO_DG_DC, "dgdc", "dg_discontinuity", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].Discontinuity,
             0, 1);
#endif

#ifdef OUTPUT_DG_MACHNUMS
  init_field(IO_DG_MN, "dgmn", "dg_machnumber", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].Machnum, 0, 1);
#endif

#ifdef OUTPUT_MIN_ANGLES
  init_field(IO_DG_MA_X, "dgmx", "dg_min_angles_x", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 5, A_SPHP,
             &SphP[0].min_angles_x[0], 0, 1);
  init_field(IO_DG_MA_Y, "dgmy", "dg_min_angles_y", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 5, A_SPHP,
             &SphP[0].min_angles_y[0], 0, 1);
  init_field(IO_DG_MA_Z, "dgmz", "dg_min_angles_z", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 5, A_SPHP,
             &SphP[0].min_angles_z[0], 0, 1);
#endif

  /* AMR */

#ifdef AMR
  init_field(IO_AMR_LEVEL, "lvl", "amrlevel", MEM_INT, FILE_INT, FILE_INT, 1, A_SPHP, &SphP[0].Level, 0, GAS_ONLY);
#endif

#ifdef OUTPUT_AMR_REFFLAG
  init_field(IO_AMR_REFFLAG, "ref", "refflag", MEM_INT, FILE_INT, FILE_INT, 1, A_SPHP, &SphP[0].refflag, 0, 1);
#endif

  /* SHOCK FINDER */

#if defined(SHOCK_FINDER_BEFORE_OUTPUT) || defined(SHOCK_FINDER_ON_THE_FLY)
  init_field(IO_SHOCK_MACHNUM, "MACH", "Machnumber", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].Machnumber,
             0, 1);
  init_field(IO_SHOCK_EDISS, "EDIS", "EnergyDissipation", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP,
             &SphP[0].EnergyDissipation, 0, 1);
#endif

#if defined(SHOCK_FINDER_BEFORE_OUTPUT_MORE) || defined(SHOCK_FINDER_ON_THE_FLY_MORE)
  init_field(IO_SHOCK_AREA, "AREA", "ShockSurfacArea", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP,
             &SphP[0].ShockSurfaceArea, 0, 1);
  init_field(IO_SHOCK_DIVVEL, "DIVV", "Divvel", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].Divvel, 0, 1);
  init_field(IO_SHOCK_DIR, "SDIR", "ShockDirection", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_SPHP, &SphP[0].ShockDir[0],
             0, 1);
#ifdef SHOCK_FINDER_BEFORE_OUTPUT_MORE
  init_field(IO_SHOCK_ZONE_FLAG, "ZONE", "ZoneFlag", MEM_INT, FILE_INT, FILE_INT, 1, A_SPHP, &SphP[0].ZoneFlag, 0, 1);
#else
  init_field(IO_SHOCK_ZONE, "SZ", "ShockZone", MEM_INT, FILE_INT, FILE_INT, 1, A_SPHP, &SphP[0].ShockZone, 0, 1);
  init_field(IO_SHOCK_SURFACE, "SS", "ShockSurface", MEM_INT, FILE_INT, FILE_INT, 1, A_SPHP, &SphP[0].ShockSurface, 0, 1);
#endif
#endif

#if defined(SHOCK_FINDER_BEFORE_OUTPUT_MORE) && !defined(COSMIC_RAYS)
  init_field(IO_SHOCK_C_PRE, "CPRE", "PreShockSoundSpeed", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP,
             &SphP[0].CpreShock, 0, 1);
  init_field(IO_SHOCK_RHO_PRE, "RPRE", "PreShockDensity", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP,
             &SphP[0].RhopreShock, 0, 1);
  init_field(IO_SHOCK_P_PRE, "PPRE", "PreShockPressure", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP,
             &SphP[0].PpreShock, 0, 1);
  init_field(IO_SHOCK_RHO_POST, "RPOS", "PostShockDensity", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP,
             &SphP[0].RhopostShock, 0, 1);
  init_field(IO_SHOCK_P_POST, "PPOS", "PostShockPressure", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP,
             &SphP[0].PpostShock, 0, 1);
  init_field(IO_SHOCK_T_PRE, "TPRE", "PreShockTemperature", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP,
             &SphP[0].TpreShock, 0, 1);
  init_field(IO_SHOCK_V_POST, "VPOS", "PostShockVelocity", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, 3, A_NONE, 0,
             io_func_post_shock_velocity, 1);
  init_field(IO_SHOCK_V_PRE, "VPRE", "PreShockVelocity", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, 3, A_NONE, 0,
             io_func_pre_shock_velocity, 1);
#endif

  /* MHD */

#ifdef MHD
  enum types_in_file mhd_read = FILE_MY_IO_FLOAT;
#if defined(MHD_SEEDFIELD) || defined(MHD_SEEDPSPEC) || defined(BINISET)
  if(RestartFlag == RESTART_IC)
    mhd_read = FILE_NONE; /* magnetic field not expected in ICs */
#endif

  /* magnetic field  */
  init_field(IO_BFLD, "BFLD", "MagneticField", MEM_NONE, FILE_MY_IO_FLOAT, mhd_read, 3, A_NONE, 0, io_func_bfield, GAS_ONLY);
  /* divergence of magnetic field  */
  init_field(IO_DIVB, "DIVB", "MagneticFieldDivergence", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].DivB, 0,
             GAS_ONLY);
  /* alternative version of divergence of magnetic field  */
  init_field(IO_DIVBALT, "DVBA", "MagneticFieldDivergenceAlternative", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP,
             &SphP[0].DivBalt, 0, GAS_ONLY);

#ifdef MHD_DEDNER
  init_field(IO_BPSI, "BPSI", "MagneticFieldPsi", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].Psi, 0,
             GAS_ONLY);
#endif
#ifdef MHD_DEDNER_VARIABLE_SPEED
  init_field(IO_VDED, "VDED", "DednerSpeed", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].DednerSpeed, 0,
             GAS_ONLY);
#endif
#endif

#ifdef MHD_CT
  enum types_in_file mhd_ct_read = FILE_MY_IO_FLOAT;
#if defined(MHD_CT_IC) || defined(MHD_SEEDFIELD)
  if(RestartFlag == RESTART_IC)
    /* magnetic vector potential not expected in ICs */
    mhd_ct_read = FILE_NONE;
#endif
  /* Magnetic Vector Potential */
  init_field(IO_AFLD, "AFLD", "MagneticVectorPotential", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, mhd_ct_read, 3, A_SPHP, &SphP[0].A[0], 0,
             GAS_ONLY);
#endif

#ifdef BECDM
  /* output type 1 particles only */
  init_field(IO_PSIR, "PSIR", "PsiRe", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_P, &P[0].PsiRe, 0, 2);
  init_field(IO_PSII, "PSII", "PsiIm", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_P, &P[0].PsiIm, 0, 2);
#endif

#ifdef NON_IDEAL_MHD
#if defined(OHMIC_DIFFUSION) || defined(IMPLICIT_OHMIC_DIFFUSION)
  /* Magnetic Vector curl */
  init_field(IO_CURLB, "CURB", "CurlB", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_SPHP, &SphP[0].CurlB[0], 0, GAS_ONLY);
#endif
#endif

  /* TRACERS */

#ifdef PASSIVE_SCALARS
  init_field(IO_PASS, "PASS", "PassiveScalars", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, PASSIVE_SCALARS, A_SPHP,
             &SphP[0].PScalars[0], 0, GAS_ONLY);
#endif

#ifdef TRACER_FIELD
  init_field(IO_TRACER, "TRCE", "TracerField", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].ConservedTracer,
             0, GAS_ONLY);
#endif

#if defined(TRACER_PARTICLE) && defined(TRACER_PART_NUM_FLUID_QUANTITIES)
  init_field(IO_TRACER_PARTICLE, "TRFQ", "FluidQuantities", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE,
             TRACER_PART_NUM_FLUID_QUANTITIES, A_P, &P[0].fluid_quantities[0], 0, pow(2, TRACER_PARTICLE));
#endif

#ifdef TRACER_MC

#ifdef OUTPUT_MCTRNUM
  init_field(IO_TRACER_MC_NumTracers, "TRNT", "NumTracers", MEM_INT, FILE_INT, FILE_INT, 1, A_TLL, 0, io_func_tracermc_numtracers,
             TRACER_MC_PARENTS);
  /*  init_field(IO_TRACER_MC_NumTracers, "TRNT", "NumTracers", MEM_INT, FILE_INT, FILE_NONE, 1, A_NONE, 0,
   * io_func_tracermc_numtracers, TRACER_MC_PARENTS);*/
  /*      init_field(IO_TRACER_MC_NumTracers, "TRNT", "NumTracers", MEM_INT, FILE_INT, FILE_INT, 1, A_TLL, 0,
   * io_func_tracermc_numtracers, GAS_ONLY); */
#endif

  init_field(IO_TRACER_MC_ID, "TRID", "TracerID", MEM_NONE, FILE_MY_ID_TYPE, FILE_MY_ID_TYPE, 1, A_TLL, 0, io_func_tracer_id,
             TRACER_MC_ONLY);
  init_snapshot_type(IO_TRACER_MC_ID, SN_MINI);

  init_field(IO_TRACER_MC_ParentID, "TRPI", "ParentID", MEM_NONE, FILE_MY_ID_TYPE, FILE_MY_ID_TYPE, 1, A_TLL, 0,
             io_func_tracer_parent_id, TRACER_MC_ONLY);
  init_snapshot_type(IO_TRACER_MC_ParentID, SN_MINI);

#ifdef TRACER_MC_NUM_FLUID_QUANTITIES
  init_field(IO_TRACER_MC_FluidQuantities, "TMFQ", "FluidQuantities", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT,
             TRACER_MC_NUM_FLUID_QUANTITIES, A_TLL, 0, io_func_tracer_fluid, TRACER_MC_ONLY);
#endif
#endif

  /* SIDM */

#ifdef SIDM
  init_field(IO_SIDM_PSUM, "SIPS", "SIDM_Psum", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, SIDM_REACTIONS, A_P, &P[0].sidm_PSum[0], 0,
             SIDM);
  init_field(IO_SIDM_NUMTOTALSCATTER, "SINT", "SIDM_NumTotalScatter", MEM_INT, FILE_INT, FILE_NONE, SIDM_REACTIONS, A_P,
             &P[0].sidm_NumTotalScatter[0], 0, SIDM);
  init_field(IO_SIDM_NUMNGB, "SINN", "SIDM_NumNgb", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_P, &P[0].sidm_NumNgb, 0, SIDM);
  init_field(IO_SIDM_HSML, "SIHS", "SIDM_Hsml", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_P, &P[0].sidm_Hsml, 0, SIDM);
  init_field(IO_SIDM_DENSITY, "SIDE", "SIDM_Density", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, SIDM_STATES, A_P,
             &P[0].sidm_Density[0], 0, SIDM);
  init_field(IO_SIDM_VELDISP, "SIVD", "SIDM_VelDisp", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, SIDM_STATES, A_P,
             &P[0].sidm_VelDisp[0], 0, SIDM);

  enum types_in_file sidm_read_state = FILE_NONE;
  init_field(IO_SIDM_STATE, "SISA", "SIDM_State", MEM_INT, FILE_INT, sidm_read_state, 1, A_P, &P[0].sidm_State, 0, SIDM);
#endif

  /* SGCHEM */

#ifdef SGCHEM
#if CHEMISTRYNETWORK > 1
  init_field(IO_DUSTTEMP, "DUST", "DustTemperature", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].DustTemp, 0,
             GAS_ONLY);
#endif
#ifdef SGCHEM_VARIABLE_Z
  init_field(IO_SGCHEM_DUSTTOGAS, "DTOG", "DusttoGasRatio", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP,
             &SphP[0].DustToGasRatio, 0, GAS_ONLY);
  init_field(IO_METL, "METL", "ElementAbundances", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 4, A_NONE, 0, io_func_sgchem_metals,
             GAS_ONLY);
#endif
  init_field(IO_CHEM, "CHEM", "ChemicalAbundances", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, SGCHEM_NUM_SPECIES, A_NONE, 0,
             io_func_sgchem, GAS_ONLY);
#ifdef VARIABLE_GAMMA
  init_field(IO_GAMMA, "GAMM", "Gamma", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].GammaE, 0, GAS_ONLY);
#endif
#endif

#ifdef SGCHEM_DUMP_THERMAL_RATES
  init_field(IO_SGCHEM_THERMAL_RATES, "SGTR", "SGCHEM_HeatCoolRates", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE,
             SGCHEM_NUM_THERMAL_RATES, A_SPHP, &SphP[0].HeatCoolRates[0], 0, GAS_ONLY);
#endif

#ifdef SGCHEM_OUTPUT_COOLTIME
  init_field(IO_COOLTIME, "CLTM", "CoolTime", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].CoolTime, 0, GAS_ONLY);
#endif

#if defined(OUTPUTCOL) && defined(TREECOLV2)
  int nmaps = 1;
#ifdef TREECOLV2_H2
  nmaps += 1;
#endif
#ifdef TREECOLV2_CO
  nmaps += 1;
#endif
#ifdef TREECOLV2_C
  nmaps += 1;
#endif

  init_field(IO_TREECOLUMN, "COLN", "TREECOL_ColDensMap", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, nmaps * NPIX, A_NONE, 0,
             io_func_treecolumn, GAS_ONLY);
#endif

#if defined(EVALPOTENTIAL) && defined(SINK_PARTICLES)
  init_field(IO_POTPEAK, "PEAK", "PotentialPeak", MEM_INT, FILE_INT, FILE_NONE, 1, A_SPHP, &SphP[0].PotentialPeak, 0, GAS_ONLY);
#endif

  /* TG */

#ifdef TGCHEM
  if(TGCD.ChemMode >= 0)
    {
      init_field(IO_CHEM, "CHEM", "ChemicalAbundances", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, TGCHEM_NUM_ABUNDANCES, A_SPHP,
                 &SphP[0].Abund[0], 0, GAS_ONLY);
      init_field(IO_GAMMA, "GAMM", "Gamma", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].Gamma, 0, GAS_ONLY);
    }
#endif

#if defined(HEALRAY) || defined(TGCHEM)
  enum types_in_file healray_read  = FILE_NONE;
  enum types_in_file healray_write = FILE_NONE;

#ifdef HEALRAY
  if(HRD.SourceFlag >= 0 && (HRD.RayIOMode == 1 || HRD.RayIOMode == 3))
    healray_write = FILE_MY_IO_FLOAT;
  if(HRD.SourceFlag >= 0 && (HRD.RayIOMode == 1 || HRD.RayIOMode == 2))
    healray_read = FILE_MY_IO_FLOAT;
#endif
#ifdef TGCHEM
  if(TGCD.ChemMode >= 0 && (TGCD.ChemIOMode == 1 || TGCD.ChemIOMode == 3))
    healray_write = FILE_MY_IO_FLOAT;
  if(TGCD.ChemMode >= 0 && (TGCD.ChemIOMode == 1 || TGCD.ChemIOMode == 2))
    healray_read = FILE_MY_IO_FLOAT;
#endif

  init_field(IO_HEATRATE, "FESC", "EscFrac", MEM_DOUBLE, healray_write, healray_read, 1, A_SPHP, &SphP[0].EscFrac, 0, GAS_ONLY);
#endif

#ifdef HEALRAY
  init_field(IO_HEATRATE, "HEAT", "HeatRate", MEM_DOUBLE, healray_write, healray_read, 1, A_SPHP, &SphP[0].HeatRate, 0, GAS_ONLY);
#endif

  /* COSMIC_RAYS */

#ifdef COSMIC_RAYS
  enum types_in_file cr_read = FILE_NONE;
#ifdef COSMIC_RAYS_IN_ICS
  if(RestartFlag == RESTART_IC)
    cr_read = FILE_MY_IO_FLOAT;
#endif
  if(RestartFlag >= 2)
    cr_read = FILE_MY_IO_FLOAT;

  init_field(IO_CRENERGY, "CREN", "CosmicRaySpecificEnergy", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, cr_read, 1, A_SPHP,
             &SphP[0].CR_SpecificEnergy, 0, GAS_ONLY);

#endif

#ifdef COSMIC_RAYS
  /* TODO: this group was not fully implemented */
  // init_field(IO_CR_C0, "CRC0", "CR_C0", MEM_NONE, FILE_MY_OUTPUT_FLOAT, FILE_NONE, 1, A_NONE, 0, 0, GAS_ONLY);
  // init_field(IO_CR_Q0, "CRQ0", "CR_Q0", MEM_NONE, FILE_MY_OUTPUT_FLOAT, FILE_NONE, 1, A_NONE, 0, 0, GAS_ONLY);
  // init_field(IO_CR_P0, "CRP0", "CR_P0", MEM_NONE, FILE_MY_OUTPUT_FLOAT, FILE_NONE, 1, A_NONE, 0, 0, GAS_ONLY);
  // init_field(IO_CR_E0, "CRE0", "CR_E0", MEM_NONE, FILE_MY_OUTPUT_FLOAT, FILE_NONE, 1, A_NONE, 0, 0, GAS_ONLY);
  // init_field(IO_CR_n0, "CRn0", "CR_n0", MEM_NONE, FILE_MY_OUTPUT_FLOAT, FILE_NONE, 1, A_NONE, 0, 0, GAS_ONLY);
  // init_field(IO_CR_ThermalizationTime, "CRco", "CR_ThermalizationTime", MEM_NONE, FILE_MY_OUTPUT_FLOAT, FILE_NONE, 1, A_NONE, 0, 0,
  // GAS_ONLY); init_field(IO_CR_DissipationTime, "CRdi", "CR_DissipationTime", MEM_NONE, FILE_MY_OUTPUT_FLOAT, FILE_NONE, 1, A_NONE,
  // 0, 0, GAS_ONLY);
#endif

#ifdef OUTPUT_CR_PRESSURE_GRADIENT
  init_field(IO_GRADPCR, "GRAC", "CRPressureGradient", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_SPHP,
             &SphP[0].Grad.dcrPressure[0], 0, GAS_ONLY);
#endif

  /* OTVET, RT, FLD */

#if defined(OTVET) && !defined(COOLING)
  /* electron abundance */
  init_field(IO_NE, "NE  ", "ElectronAbundance", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_NONE, 0, io_func_ne, GAS_ONLY);
  /* dimensionless fraction */
  init_units(IO_NE, 0, 0, 0, 0, 0, 0);
#endif

#if defined(OTVET) || (defined(MRT) && !defined(MRT_CHEM_SG))
#ifdef OTVET
  init_field(IO_OTVET_GAMMA, "GAMM", "Gamma", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, OT_N_BINS, A_SPHP, &SphP[0].DensPhot,
             0, GAS_ONLY);
#endif
  /* neutral hydrogen abundance */
  init_field(IO_OTVET_HI, "HI  ", "HI_Fraction", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].HI, 0,
             GAS_ONLY);
  /* ionized hydrogen abundance */
  init_field(IO_OTVET_HII, "HII ", "HII_Fraction", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].HII, 0,
             GAS_ONLY);
#endif

#if(defined(OTVET) && defined(OTVET_INCLUDE_HE)) || (defined(MRT) && defined(MRT_INCLUDE_HE) && !defined(MRT_CHEM_SG))
  /* neutral Helium */
  init_field(IO_OTVET_HeI, "HeI ", "HeI_Fraction", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].HeI, 0,
             GAS_ONLY);
  /* ionized Helium */
  init_field(IO_OTVET_HeII, "HeII", "HeII_Fraction", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].HeII, 0,
             GAS_ONLY);
  /* doubly ionized Helium */
  init_field(IO_OTVET_HeIII, "He3I", "HeIII_Fraction", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].HeIII, 0,
             GAS_ONLY);
#endif

#if defined(OTVET) && defined(OTVET_SCATTER_SOURCE) && defined(OTVET_OUTPUT_SOURCEHSML)
  init_field(IO_OTVET_SRC_HSML, "SHSM", "SourceHsml", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_P, &P[0].OtvetHsml, 0,
             ALL_TYPES);
  init_field(IO_OTVET_SRC_HSML, "SRHO", "SourceGasDensity", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_P,
             &P[0].OtvetGasDensity, 0, ALL_TYPES);
#endif

#if defined(OTVET) && defined(OTVET_OUTPUT_ET)
  init_field(IO_OTVET_ET, "EDDT", "EddingtonTensor", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 6, A_SPHP, &SphP[0].ET[0], 0,
             GAS_ONLY);
#endif

#if defined(RT_SHORT_CHARACTERISTICS) || defined(RT_ADVECT) || defined(MRT)
  /* neutral fraction */
  init_field(IO_nHI, "nHI ", "NeutralFraction", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].nHI, 0, GAS_ONLY);
#endif

#ifdef RT_ADVECT
  int rt_components = RT_N_DIR;

#ifdef RT_COMBINE_N_DIR_IN_OUTPUT
  rt_components = 1;
#endif
  init_field(IO_PHOTONS, "PHOT", "Photons", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, rt_components, A_NONE, 0, io_func_rt_photons,
             GAS_ONLY);
  init_field(IO_PHOTDENSITY, "PHOD", "PhotonDensity", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, rt_components, A_NONE, 0,
             io_func_rt_photdensity, GAS_ONLY);
#endif

#ifdef MRT
  enum types_in_file photon_read = FILE_MY_IO_FLOAT;
  if(RestartFlag == RESTART_IC)
    photon_read = FILE_NONE; /* photons not present in ICs */
  init_field(IO_PHOTDENSITY, "PHOD", "PhotonDensity", MEM_NONE, FILE_MY_IO_FLOAT, photon_read, MRT_BINS, A_NONE, 0,
             io_func_mrt_photondensity, GAS_ONLY);
#ifdef MRT_OUTPUT_FLUX
  init_field(IO_RTFX, "RTFX", "RTFluxX", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, MRT_BINS, A_NONE, 0, io_func_mrt_rtfx,
             GAS_ONLY);
  init_field(IO_RTFY, "RTFY", "RTFluxY", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, MRT_BINS, A_NONE, 0, io_func_mrt_rtfy,
             GAS_ONLY);
  init_field(IO_RTFZ, "RTFZ", "RTFluxZ", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, MRT_BINS, A_NONE, 0, io_func_mrt_rtfz,
             GAS_ONLY);
#endif
#ifdef MRT_IR_PHOTON_TRAPPING
  init_field(IO_PHOTONS, "TRPH", "TrappedPhotons", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, MRT_BINS, A_NONE, 0, io_func_mrt_photons,
             GAS_ONLY);
#endif
#if defined(MRT_IR) && !defined(MRT_CHEM_SG)
  init_field(IO_FLD_KAPPA_P, "MIRR", "KappaIR_R", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, MRT_BINS, A_SPHP, &SphP[0].KappaIR_R,
             0, 1);
  init_field(IO_FLD_KAPPA_R, "MIRP", "KappaIR_P", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, MRT_BINS, A_SPHP, &SphP[0].KappaIR_P,
             0, 1);
#endif
#endif

#if defined(RT_ADVECT) && defined(RT_INCLUDE_HE)
  /* He neutral fraction */
  init_field(IO_nHeI, "nHe2", "HeNeutralFraction", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].nHeI, 0, GAS_ONLY);
  /* HeII fraction */
  init_field(IO_nHeII, "nHe3", "HeIIFraction", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].nHeII, 0, GAS_ONLY);
#endif

#if defined(RT_SHORT_CHARACTERISTICS) && defined(RT_OUTPUT_COL_DENS)
  init_field(IO_nHeI, "ColD", "ColumnDensity", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].ColumnDensity, 0,
             GAS_ONLY);
#endif

#ifdef FLD
  init_field(IO_FLD, "GAMM", "Gamma", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].n_gamma, 0, 1);
  init_field(IO_FLD_LAMBDA, "LAMD", "Lambda", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].Lambda, 0, 1);
  init_field(IO_FLD_KAPPA_P, "KA_P", "Kappa_P", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].Kappa_P, 0, 1);
  init_field(IO_FLD_KAPPA_R, "KA_R", "Kappa_R", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].Kappa_R, 0, 1);
  init_field(IO_FLD_GRAD, "GRAG", "gradgamma", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_SPHP, &SphP[0].Grad.dngamma, 0, 1);
#ifdef FLD_CONES
  init_field(IO_FLD_GAMMAS, "GAMS", "gammas", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, FLD_NCONES, A_SPHP, &SphP[0].gammas, 0,
             1);
#endif
#endif

  /* RADCOOL */

#ifdef RADCOOL
  init_field(IO_PHIOS, "OLDS", "RadiationFieldIntensity-OldStars", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].Phios,
             0, GAS_ONLY);
  init_field(IO_PHINS, "NEWS", "RadiationFieldIntensity-NewStars", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].Phins,
             0, GAS_ONLY);
#endif

#ifdef RADCOOL_HOTHALO
  init_field(IO_PHIT6, "GAT6", "RadiationFieldIntensity-GasLogTempMean6", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP,
             &SphP[0].PhiT6, 0, GAS_ONLY);
  init_field(IO_PHIT6, "GAT7", "RadiationFieldIntensity-GasLogTempMean7", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP,
             &SphP[0].PhiT7, 0, GAS_ONLY);
  init_field(IO_PHIT6, "GAT8", "RadiationFieldIntensity-GasLogTempMean8", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP,
             &SphP[0].PhiT8, 0, GAS_ONLY);
  init_field(IO_TEMP, "TEMP", "Temperature", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].Temperature, 0, GAS_ONLY);
#endif

  /* SIMPLEX */

#ifdef SIMPLEX
#if SX_CHEMISTRY == 3
  init_field(IO_SX_PHOTON_RATES, "PHR", "PhotonRates", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, SX_NRATES, A_SPHP,
             &SphP[0].sxPhotonRates[0], 0, GAS_ONLY);
#ifdef SX_OUTPUT_FLUX
  init_field(IO_SX_PHOTON_FLUX, "PHF", "PhotonFlux", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, SX_NFREQ, A_SPHP,
             &SphP[0].sxPhotonFlux[0], 0, GAS_ONLY);
#endif
#elif SX_CHEMISTRY == 4
  init_field(IO_SX_MASSFRACT, "MFRC", "MassFractions", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, SX_NMASSFRACT, A_SPHP,
             &SphP[0].MassFract, 0, GAS_ONLY);
  init_field(IO_SX_TEMPERATURE, "TEMP", "Temperature", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].Temperature,
             0, GAS_ONLY);
#if(SX_SOURCES == 4)
  init_field(IO_SX_SEDS, "SEDs", "SEDs", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, SX_NFREQ, A_STARP, &StarP[0].SEDs, 0,
             STARS_ONLY);
#endif
#endif
#endif

  /* GRACKLE */

#ifdef GRACKLE
  init_field(IO_GRACKLE_TEMP, "TEMP", "GrackleTemperature", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_NONE, 0, io_func_grackle_temp,
             GAS_ONLY);
  init_field(IO_GRACKLE_COOL_TIME, "GCLT", "GrackleCoolTime", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_NONE, 0,
             io_func_grackle_cooltime, GAS_ONLY);
#endif

#if defined(GRACKLE) && !defined(GRACKLE_TAB)
  enum types_in_file grackle_read_in = FILE_NONE;
#ifdef GRACKLE_ABUNDANCE_IN_ICS
  grackle_read_in = FILE_MY_IO_FLOAT;
#endif
  init_field(IO_GRACKLE_E, "ELEC", "ElectronMassFraction", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, grackle_read_in, 1, A_SPHP, &SphP[0].e_frac,
             0, GAS_ONLY);
  init_field(IO_GRACKLE_HI, "HI  ", "HIMassFraction", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, grackle_read_in, 1, A_SPHP,
             &SphP[0].GrackleSpeciesFraction[0], 0, GAS_ONLY);
  init_field(IO_GRACKLE_HII, "HII ", "HIIMassFraction", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, grackle_read_in, 1, A_SPHP,
             &SphP[0].GrackleSpeciesFraction[1], 0, GAS_ONLY);
  init_field(IO_GRACKLE_HeI, "HEI ", "HeIMassFraction", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, grackle_read_in, 1, A_SPHP,
             &SphP[0].GrackleSpeciesFraction[2], 0, GAS_ONLY);
  init_field(IO_GRACKLE_HeII, "HEII", "HeIIMassFraction", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, grackle_read_in, 1, A_SPHP,
             &SphP[0].GrackleSpeciesFraction[3], 0, GAS_ONLY);
  init_field(IO_GRACKLE_HeIII, "HE3I", "HeIIIMassFraction", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, grackle_read_in, 1, A_SPHP,
             &SphP[0].GrackleSpeciesFraction[4], 0, GAS_ONLY);
#endif

#if defined(GRACKLE) && defined(GRACKLE_H2)
  init_field(IO_GRACKLE_HM, "HM  ", "HMMassFraction", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, grackle_read_in, 1, A_SPHP,
             &SphP[0].GrackleSpeciesFraction[5], 0, GAS_ONLY);
  init_field(IO_GRACKLE_H2I, "H2I ", "H2IMassFraction", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, grackle_read_in, 1, A_SPHP,
             &SphP[0].GrackleSpeciesFraction[6], 0, GAS_ONLY);
  init_field(IO_GRACKLE_H2II, "H2II", "H2IIMassFraction", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, grackle_read_in, 1, A_SPHP,
             &SphP[0].GrackleSpeciesFraction[7], 0, GAS_ONLY);
#endif

#if defined(GRACKLE) && defined(GRACKLE_D)
  init_field(IO_GRACKLE_DI, "DI  ", "DIMassFraction", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, grackle_read_in, 1, A_SPHP,
             &SphP[0].GrackleSpeciesFraction[8], 0, GAS_ONLY);
  init_field(IO_GRACKLE_DII, "DII ", "DIIMassFraction", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, grackle_read_in, 1, A_SPHP,
             &SphP[0].GrackleSpeciesFraction[9], 0, GAS_ONLY);
  init_field(IO_GRACKLE_HDI, "HDI ", "HDIMassFraction", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, grackle_read_in, 1, A_SPHP,
             &SphP[0].GrackleSpeciesFraction[10], 0, GAS_ONLY);
#endif

  /* CHIMES */
#ifdef CHIMES
  init_field(IO_CHIMES_ABUNDANCES, "CHIM", "ChimesAbundances", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT,
             ChimesGlobalVars.totalNumberOfSpecies, A_NONE, 0, io_func_chimes_abundances, GAS_ONLY);
  init_units(IO_CHIMES_ABUNDANCES, 0., 0., 0., 0., 0., 0.);
  init_snapshot_type(IO_CHIMES_ABUNDANCES, SN_FULL);

  init_field(IO_CHIMES_MU, "CHMU", "ChimesMu", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_NONE, 0, io_func_chimes_mu, GAS_ONLY);
  init_units(IO_CHIMES_MU, 0., 0., 0., 0., 0., 0.);
  init_snapshot_type(IO_CHIMES_MU, SN_FULL);
#endif /* CHIMES */

  /* OTHER */

#ifdef SAVE_HSML_IN_SNAPSHOT
  init_field(IO_SUBFINDDENSITY, "SFDE", "SubfindDensity", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_PS, &PS[0].SubfindDensity, 0,
             ALL_TYPES);
  init_units(IO_SUBFINDDENSITY, -3., 2., -3., 1., 0., All.UnitDensity_in_cgs);
  init_snapshot_type(IO_SUBFINDDENSITY, SN_NO_SUBBOX);

  init_field(IO_SUBFINDDMDENSITY, "SFDD", "SubfindDMDensity", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_PS,
             &PS[0].SubfindDMDensity, 0, ALL_TYPES);
  init_units(IO_SUBFINDDMDENSITY, -3., 2., -3., 1., 0., All.UnitDensity_in_cgs);
  init_snapshot_type(IO_SUBFINDDMDENSITY, SN_NO_SUBBOX);

  init_field(IO_SUBFINDHSML, "SFHS", "SubfindHsml", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_PS, &PS[0].SubfindHsml, 0,
             ALL_TYPES);
  init_units(IO_SUBFINDHSML, 1., -1., 1., 0., 0., All.UnitLength_in_cm);
  init_snapshot_type(IO_SUBFINDHSML, SN_NO_SUBBOX);

  init_field(IO_SUBFINDVELDISP, "SFVD", "SubfindVelDisp", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_PS, &PS[0].SubfindVelDisp, 0,
             ALL_TYPES);
  init_units(IO_SUBFINDVELDISP, 0.0, 0.0, 0.0, 0.0, 1.0, All.UnitVelocity_in_cm_per_s);
  init_snapshot_type(IO_SUBFINDVELDISP, SN_NO_SUBBOX);
#endif

#if defined(FOF_FUZZ_SORT_BY_NEAREST_GROUP) && (FOF_FUZZ_SORT_BY_NEAREST_GROUP == 1)
  init_field(IO_GROUPNR, "GROU", "GroupNr", MEM_INT, FILE_INT, FILE_NONE, 1, A_PS, &PS[0].GroupNr, 0, ALL_TYPES);
  init_units(IO_GROUPNR, 0, 0, 0, 0, 0, 0);
#endif

#if defined(REFINEMENT_HIGH_RES_GAS) && !defined(TGSET)
  init_field(IO_HIGHRESMASS, "HRGM", "HighResGasMass", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP,
             &SphP[0].HighResMass, 0, GAS_ONLY);
#endif
#if defined(REFINEMENT_CGM)
  init_field(IO_HIGHRESMASSCGM, "HRCM", "HighResGasMassCGM", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP,
             &SphP[0].HighResMassCGM, 0, GAS_ONLY);
#endif
#if defined(REFINEMENT_HIGH_RES_GAS)
  init_field(IO_ALLOWREFINEMENT, "REF ", "AllowRefinement", MEM_INT, FILE_INT, FILE_INT, 1, A_SPHP, &SphP[0].AllowRefinement, 0,
             GAS_ONLY);
#endif

#if defined(EOS_DEGENERATE) || defined(EOS_OPAL) || defined(EOS_PASSIVE)
  init_field(IO_EOSXNUC, "XNUC", "NuclearComposition", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, EOS_NSPECIES, A_SPHP,
             &SphP[0].Composition[0], 0, GAS_ONLY);
#endif

#if defined(EOS_DEGENERATE) || defined(EOS_OPAL)
  init_field(IO_TEMP, "TEMP", "Temperature", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].EOSTemperature, 0, GAS_ONLY);
#endif

#ifdef NUCLEAR_NETWORK
  init_field(IO_NUC_DEDT, "DEDT", "NuclearEnergyGenerationRate", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP,
             &SphP[0].dedt, 0, GAS_ONLY);
#endif

#if defined(VS_TURB) || defined(AB_TURB)
  init_field(IO_VSTURB_DISS, "VSDI", "TurbulenceDissipation", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP,
             &SphP[0].DuDt_diss, 0, GAS_ONLY);
  init_field(IO_VSTURB_DISS, "VSDR", "TurbulenceDriving", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].DuDt_drive,
             0, GAS_ONLY);
#endif

#ifdef REFINEMENT_RPS
  init_field(IO_RPS, "RPS ", "RPSGalaxyMass", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].RPSGalaxyMass, 0,
             GAS_ONLY);
#endif

#ifdef STELLARAGE
  /* TODO: STARS_ONLY not used (old logic still used in get_particles_in_block()) */
  init_field(IO_AGE, "AGE ", "StellarFormationTime", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_P, &P[0].StellarAge, 0,
             STARS_ONLY);
#endif

#ifdef LOCAL_FEEDBACK_PARTICLES
  /* TODO: STARS_ONLY not used (old logic still used in get_particles_in_block()) */
  /* number of feedback events done */
  init_field(IO_LOCFBEVENT, "FBDO", "FeedbackDone", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_P, &P[0].FeedbackDone, 0,
             STARS_ONLY);
#endif

#ifdef OUTPUT_STICKYFLAGS
  /* Flag for sticky cells */
  init_field(IO_STKY, "STKY", "StickyFlag", MEM_INT, FILE_INT, FILE_INT, 1, A_SPHP, &SphP[0].StickyFlag, 0, GAS_ONLY);
#endif

#ifdef METALS
  /* TODO: GAS_AND_STARS not used (old logic still used in get_particles_in_block()) */
  /* Gas and star metallicity */
  init_field(IO_Z, "Z   ", "Metallicity", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_NONE, 0, io_func_metallicity,
             GAS_AND_STARS);
#endif

#ifdef SGS_TURBULENCE
  enum types_in_file sgs_t_read = FILE_NONE;
#ifdef SGS_TURBULENCE_IN_ICS
  if(RestartFlag == RESTART_IC)
    sgs_t_read = FILE_MY_IO_FLOAT;
#endif
  init_field(IO_SGS_T_SPECIFICENERGY, "KSGS", "SgsTurbulenceSpecificEnergy", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1,
             A_SPHP, &SphP[0].SgsTData.SpecificEnergy, 0, GAS_ONLY);
#ifdef SGS_TURBULENCE_OUTPUT_SGS_PRESSURE
  init_field(IO_SGS_T_PRESSURE, "PSGS", "SgsTurbulencePressure", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP,
             &SphP[0].SgsTData.Pressure, 0, GAS_ONLY);
#endif
#endif

#if defined(SECOND_DERIVATIVES) && defined(OUTPUT_HESSIAN)
  init_field(IO_HESSIAN_VELX, "HVELX", "HessianVelocityX", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_NONE, 9, A_SPHP,
             &SphP[0].Hessian.ddvelx[0][0], 0, GAS_ONLY);
  init_field(IO_HESSIAN_VELY, "HVELY", "HessianVelocityY", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_NONE, 9, A_SPHP,
             &SphP[0].Hessian.ddvely[0][0], 0, GAS_ONLY);
  init_field(IO_HESSIAN_VELZ, "HVELZ", "HessianVelocityZ", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_NONE, 9, A_SPHP,
             &SphP[0].Hessian.ddvelz[0][0], 0, GAS_ONLY);
#endif

#ifdef OUTPUT_DENSTROPHY
  init_field(IO_DENSTROPHY, "DNST", "Denstrophy", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_NONE, 0, io_func_denstrophy,
             GAS_ONLY);
#endif

#ifdef OUTPUT_SGS_T_PRESSURE_GRADIENT
  init_field(IO_GRADPSGS, "GPSGS", "SgsTPressureGradient", MEM_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_SPHP,
             &SphP[0].Grad.dsgstPressure[0], 0, GAS_ONLY);
#endif

#ifdef BOUNDARY_FLAG
  init_field(IO_BFLAG, "BFLG", "BoundaryFlag", MEM_INT, FILE_INT, FILE_INT, 1, A_SPHP, &SphP[0].BoundaryFlag, 0, GAS_ONLY);
#endif

#if(defined(BIERMANN_BATTERY) || defined(DURRIVE_BATTERY))
  init_field(IO_NELEC, "NEL ", "ElectronNumberDensity", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].n_elec,
             0, GAS_ONLY);
#endif
#ifdef BIERMANN_BATTERY
  init_field(IO_PELEC, "PEL ", "ElectronPressure", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].p_elec, 0,
             GAS_ONLY);
#endif
#ifdef DURRIVE_BATTERY
  init_field(IO_PDOTELEC, "PDEL", "ElectronMomentumTransferRate", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_SPHP,
             &SphP[0].pdot_elec, 0, GAS_ONLY);
#endif
#ifdef MAGNETIC_BATTERIES_OUTPUT_GRADIENTS
  init_field(IO_DNELEC, "DNEL", "GradElectronNumberDensity", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_SPHP,
             &SphP[0].Grad_n_elec, 0, GAS_ONLY);
  init_field(IO_DPELEC, "DPEL", "GradElectronPressure", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_SPHP,
             &SphP[0].Grad_p_elec, 0, GAS_ONLY);
  init_field(IO_DPDOTELEC, "DPDE", "GradElectronMomentumTransferRate", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 9, A_SPHP,
             &SphP[0].Grad_pdot_elec, 0, GAS_ONLY);
#endif
}
