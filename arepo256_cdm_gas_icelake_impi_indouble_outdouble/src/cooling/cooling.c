/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/cooling/cooling.c
 * \date        MM/YYYY
 * \author
 * \brief       Module for gas radiative cooling
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

#include "../allvars.h"
#include "../proto.h"

#if defined(COOLING) && !defined(GRACKLE) && !defined(CHIMES)

#ifdef RADCOOL
static double Tmin = 2.0; /*Cooling table only tabulated from log10(T)=2.0 */
#else
static double Tmin = 0.0; /**< min temperature in log10 */
#endif
static double Tmax = 9.0; /**< max temperature in log10 */
static double deltaT;     /**< log10 of temperature spacing in the interpolation tables */

static GasState gs; /**< gas state */

static RateTable *RateT; /**< tabulated rates */

static PhotoTable *PhotoTUVB; /**< photo-ionization/heating rate table for UV background */
#ifdef GFM_AGN_RADIATION
static PhotoTable *PhotoTAGN; /**< photo-ionization/heating rate table for AGN background */
#endif

#ifdef RADCOOL
static PhotoTable *PhotoTRAD_OS;
static PhotoTable *PhotoTRAD_NS;
static double fac_rad_mass, fac_rad_length, fac_rad_length_init, fac_rad, fac_rad_ns, fac_rad_timeon;
#ifdef RADCOOL_HOTHALO
static PhotoTable *PhotoTRAD_T6;
static PhotoTable *PhotoTRAD_T7;
static PhotoTable *PhotoTRAD_T8;
static double fac_rad_dens, fac_rad_massdens_mp2_cm5;
#endif
#endif

#ifdef UVB_SELF_SHIELDING
static SelfShieldingTable *SelfShieldingParams; /**< parameters for self-shielding fitting formula */
static int Nselfshieldingtab;                   /**< length of self-shielding fitting formula table */
#endif

static PhotoCurrent pc; /**< current interpolated photo rates */

static int NheattabUVB; /**< length of UVB photo table */

static DoCoolData DoCool; /**< cooling data */

#ifdef GFM_COOLING_METAL
static double MetallicityFloor, LogMetallicityFloorInSolar, LogMinMetalTemp;
#endif

#ifdef RADCOOL

/** \brief Initializes the units at the start of the simulation.
 */
void init_radcool_units(void)
{
  fac_rad_mass        = All.UnitMass_in_g / SOLARMASS_in_g / All.HubbleParam;
  fac_rad_length_init = All.UnitLength_in_cm / KPC_in_cm / All.HubbleParam;
  fac_rad_timeon      = TIMEON_NEWSTARS * GYR_to_YR;
  if(All.ComovingIntegrationOn == 0)
    {
      fac_rad_length = fac_rad_length_init;
      fac_rad        = fac_rad_mass / (fac_rad_length * fac_rad_length);
      fac_rad_ns     = fac_rad / fac_rad_timeon;
#ifdef RADCOOL_HOTHALO
      fac_rad_dens             = fac_rad / fac_rad_length;
      fac_rad_massdens_mp2_cm5 = fac_rad * fac_rad_dens * FACTOR_Mp2_cm5 /
                                 (4 * 3.1415); /* rho*mass/r2 [Msun^2/kpc^5] to rho*mass / 4 pi r^2 [protonmass^2/m^5] */
#endif
    }
  mpi_printf("RADCOOL:INIT ONCE RADCOOL : mass factor = %le\tlength factor = %le\ttime on factor = %le\n", fac_rad_mass,
             fac_rad_length_init, fac_rad_timeon);
}

/** \brief Initializes the units as a fucntion of present time  - updates every time step.
 */
void set_radcool_units_for_current_time(void)
{
  double expfactor_hh;
  expfactor_hh = 1.0 / (All.cf_redshift + 1.0);

  if(All.ComovingIntegrationOn)
    {
      fac_rad_length = fac_rad_length_init * expfactor_hh;
      fac_rad        = fac_rad_mass / (fac_rad_length * fac_rad_length);
      fac_rad_ns     = fac_rad / fac_rad_timeon;
#ifdef RADCOOL_HOTHALO
      fac_rad_dens             = fac_rad / fac_rad_length;
      fac_rad_massdens_mp2_cm5 = fac_rad * fac_rad_dens * FACTOR_Mp2_cm5 /
                                 (4.0 * 3.1415); /* rho*mass/r2 [Msun^2/kpc^5] to rho*mass / 4 pi r^2 [protonmass^2/m^5] */
#endif
    }
  mpi_printf("RADCOOL:INIT EVERY TIMESTEP: time = %le\tlength factor = %le\tflux factor = %le", expfactor_hh, fac_rad_length, fac_rad);
#ifdef RADCOOL_HOTHALO
  mpi_printf("\t hothalo factor = %le", fac_rad_massdens_mp2_cm5);
#endif
  mpi_printf("\n");
}

#endif

#if(defined(GFM_AGN_RADIATION) || defined(GFM_UVB_CORRECTIONS) || defined(UVB_SELF_SHIELDING) || defined(RADCOOL)) && !defined(SFR_MCS)
/** \brief Reset the state of background radiation.
 *
 *   The function resets the fields of the struct GasState that refer
 *   to the background radiation state. If the UV background is selected,
 *   its ionization parameters are updated.
 */
void reset_radiation_state(void)
{
#ifdef GFM_AGN_RADIATION
  gs.FlagAGNBackground  = 0;
  gs.AGNBolIntensity    = GFM_MIN_AGN_BOL_INTENSITY;
  gs.LogAGNBolIntensity = log10(GFM_MIN_AGN_BOL_INTENSITY);
#endif
#ifdef RADCOOL
  gs.FlagRadcool = 0;
  gs.Phios       = 0.0;
  gs.LogPhios    = log10(PHIOS_MIN);
  gs.Phins       = 0.0;
  gs.LogPhins    = log10(PHINS_MIN);
#ifdef RADCOOL_HOTHALO
  gs.PhiT6    = 0.0;
  gs.LogPhiT6 = log10(PHIT_MIN);
  gs.PhiT7    = 0.0;
  gs.LogPhiT7 = log10(PHIT_MIN);
  gs.PhiT8    = 0.0;
  gs.LogPhiT8 = log10(PHIT_MIN);
#endif
#endif
  IonizeParamsUVB();  // FIXME: make faster through saved values
}

/** \brief Update the state of background radiation.
 *
 *   The function updates the fields of the struct GasState that refer
 *   to the background radiation state. In particular these fields are
 *   first reset by calling the function reset_radiation_state() and
 *   then updated according to the type of ionizing backgruond selected.
 *
 *   \param rho   the proper density of the gas cell
 *   \param xh    hydrogen mass fraction of the gas
 *   \param bolIntensity AGN bolometric intensity
 */
#ifdef RADCOOL
void update_radiation_state(MyFloat rho, MyFloat xh, MyFloat Phios, MyFloat Phins
#ifdef RADCOOL_HOTHALO
                            ,
                            MyFloat PhiT6, MyFloat PhiT7, MyFloat PhiT8
#endif
)
#else
void update_radiation_state(MyFloat rho, MyFloat xh, MyFloat bolIntensity)
#endif
{
  reset_radiation_state();

#if defined(GFM_AGN_RADIATION) || defined(UVB_SELF_SHIELDING) || defined(RADCOOL)
  double nH = xh * rho / PROTONMASS * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
#endif

#ifdef GFM_AGN_RADIATION
  double ionization_parameter = bolIntensity / nH;

  if(ionization_parameter > GFM_MAX_IONIZATION_PARAMETER && nH < All.SelfShieldingDensity)
    warn("GFM_AGN_RADIATION: ionization_parameter=%g > GFM_MAX_IONIZATION_PARAMETER=%g; bolIntensity=%g nH=%g", ionization_parameter,
         GFM_MAX_IONIZATION_PARAMETER, bolIntensity, nH);

  if(bolIntensity > 0 && ionization_parameter >= GFM_MIN_IONIZATION_PARAMETER && nH < All.SelfShieldingDensity)
    {
      /* cap bolIntensity such that ionization_parameter is at its maximum reliable value: GFM_MAX_IONIZATION_PARAMETER */
      if(ionization_parameter > GFM_MAX_IONIZATION_PARAMETER)
        bolIntensity *= (GFM_MAX_IONIZATION_PARAMETER / ionization_parameter);

      gs.FlagAGNBackground = 1;
      CellsWithAGNRadiation += 1;
      gs.AGNBolIntensity    = bolIntensity;
      gs.LogAGNBolIntensity = log10(bolIntensity);
      IonizeParamsAGN();
    }
#endif

#ifdef RADCOOL

  gs.FlagRadcool = 1;
  if((All.SelfShieldingOn == 1) && (nH >= All.SelfShieldingDensity))
    {
      gs.Phios = 0.0;
      gs.Phins = 0.0;

      gs.LogPhios = log10(PHIOS_MIN);
      gs.LogPhins = log10(PHINS_MIN);

#ifdef RADCOOL_HOTHALO
      gs.PhiT6 = 0.0;
      gs.PhiT7 = 0.0;
      gs.PhiT8 = 0.0;

      gs.LogPhiT6 = log10(PHIT_MIN);
      gs.LogPhiT7 = log10(PHIT_MIN);
      gs.LogPhiT8 = log10(PHIT_MIN);

#endif
    }
  else
    {
#ifndef TEST_COOLING_METAL

      gs.Phios = Phios * fac_rad;    /*Msun/kpc^2 */
      gs.Phins = Phins * fac_rad_ns; /*Msun/kpc^2/yr = SFR/kpc^2 */

#ifdef RADCOOL_HOTHALO
      gs.PhiT6 = PhiT6 * fac_rad_massdens_mp2_cm5;
      gs.PhiT7 = PhiT7 * fac_rad_massdens_mp2_cm5;
      gs.PhiT8 = PhiT8 * fac_rad_massdens_mp2_cm5;
#endif

#else
      gs.Phios = Phios;
      gs.Phins = Phins;
#ifdef RADCOOL_HOTHALO
      gs.PhiT6 = PhiT6;
      gs.PhiT7 = PhiT7;
      gs.PhiT8 = PhiT8;
#endif
#endif
      if(All.OldStarsOn)
        {
          if(gs.Phios > PHIOS_MAX)
            {
              gs.Phios    = PHIOS_MAX * (1.0 - EPS);
              gs.LogPhios = log10(gs.Phios);
            }
          else if(gs.Phios < PHIOS_MINVAL)
            {
              gs.Phios    = 0.0;
              gs.LogPhios = log10(PHIOS_MIN);
            }
          else
            gs.LogPhios = log10(gs.Phios);
        }
      else
        {
          gs.Phios    = 0.0;
          gs.LogPhios = log10(PHIOS_MIN);
        }

      if(All.NewStarsOn)
        {
          if(gs.Phins > PHINS_MAX)
            {
              gs.Phins    = PHINS_MAX * (1.0 - EPS);
              gs.LogPhins = log10(gs.Phins);
            }
          else if(gs.Phins < PHINS_MINVAL)
            {
              gs.Phins    = 0.0;
              gs.LogPhins = log10(PHINS_MIN);
            }
          else
            gs.LogPhins = log10(gs.Phins);
        }
      else
        {
          gs.Phins    = 0.0;
          gs.LogPhins = log10(PHINS_MIN);
        }

#ifdef RADCOOL_HOTHALO
      if((All.HotHaloOn == 1) && (All.cf_redshift < HOTHALO_ON_zFACTOR))
        {
          if(gs.PhiT6 > PHIT_MAX)
            {
              gs.PhiT6    = PHIT_MAX * (1.0 - EPS);
              gs.LogPhios = log10(gs.PhiT6);
            }
          else if(gs.PhiT6 < PHIT_MINVAL)
            {
              gs.PhiT6    = 0.0;
              gs.LogPhiT6 = log10(PHIT_MIN);
            }
          else
            gs.LogPhiT6 = log10(gs.PhiT6);

          if(gs.PhiT7 > PHIT_MAX)
            {
              gs.PhiT7    = PHIT_MAX * (1.0 - EPS);
              gs.LogPhiT7 = log10(gs.PhiT7);
            }
          else if(gs.PhiT7 < PHIT_MINVAL)
            {
              gs.PhiT7    = 0.0;
              gs.LogPhiT7 = log10(PHIT_MIN);
            }
          else
            gs.LogPhiT7 = log10(gs.PhiT7);

          if(gs.PhiT8 > PHIT_MAX)
            {
              gs.PhiT8    = PHIT_MAX * (1.0 - EPS);
              gs.LogPhiT8 = log10(gs.PhiT8);
            }
          else if(gs.PhiT8 < PHIT_MINVAL)
            {
              gs.PhiT8    = 0.0;
              gs.LogPhiT8 = log10(PHIT_MIN);
            }
          else
            gs.LogPhiT8 = log10(gs.PhiT8);

          if(gs.Phios < PHIOS_MINVAL_LZ)
            {
              gs.Phios    = 0.0;
              gs.LogPhios = log10(PHIOS_MIN);
            }

          if(gs.Phins < PHINS_MINVAL_LZ)
            {
              gs.Phins    = 0.0;
              gs.LogPhins = log10(PHINS_MIN);
            }
        }
      else
        {
          gs.PhiT6 = 0.0;
          gs.PhiT7 = 0.0;
          gs.PhiT8 = 0.0;

          gs.LogPhiT6 = log10(PHIT_MIN);
          gs.LogPhiT7 = log10(PHIT_MIN);
          gs.LogPhiT8 = log10(PHIT_MIN);
        }

#endif
    }
  IonizeParamsRADCOOL();
#endif

#ifdef GFM_UVB_CORRECTIONS
  if(rho > 0)
    {
      double delta = rho / (gs.RhoGasMean * pow(All.Time, 3));
      if(delta > All.UV_HeII_threshold)
        delta = All.UV_HeII_threshold;
      gs.UV_HeII_Factor = All.UV_HeII_alpha * pow(delta, All.UV_HeII_beta);
    }
#endif

#ifdef UVB_SELF_SHIELDING
  double selfshielding_factor = calc_self_shielding_factor(nH);
  pc.gJH0 *= selfshielding_factor;
  pc.gJHe0 *= selfshielding_factor;
  pc.gJHep *= selfshielding_factor;
  pc.epsH0 *= selfshielding_factor;
  pc.epsHe0 *= selfshielding_factor;
  pc.epsHep *= selfshielding_factor;
#endif
}
#endif

#ifdef GFM_LAMBDA
void cooling_set_partindex(int i) { gs.partindex = i; }
#endif

#if defined(GFM_COOLING_METAL) && !defined(MRT_METAL_COOLING) && !defined(SFR_MCS)
/** \brief Update the gas state for metal line cooling.
 *
 *   The function updates several properties of a given gas cell (the H mass fraction,
 *   the metallicity, the hydrogen and helium number density) needed
 *   to apply metal line cooling. Arguments are passed in code units.
 *
 *   \param rho   the proper density of the gas cell
 *   \param xh    hydrogen mass fraction of the gas
 *   \param metallicity gas metallicity
 */
void update_gas_state(MyFloat rho, MyFloat xh, MyFloat metallicity)
{
  gs.XH      = xh;
  gs.yhelium = (1 - gs.XH - metallicity) / (4. * gs.XH);

  if(!gsl_finite(gs.XH) || !gsl_finite(gs.yhelium) || !gsl_finite(metallicity))
    terminate("xh=%g, yhelium=%g, metallicity=%g", xh, gs.yhelium, metallicity);

  if(metallicity > MetallicityFloor)
    /* log_MetallicityInSolar = log10(Z/Z_solar) */
    gs.log_MetallicityInSolar = log10(metallicity / GFM_SOLAR_METALLICITY);
  else
    /* log_MetallicityInSolar = log10(Z), does not matter */
    gs.log_MetallicityInSolar = LogMetallicityFloorInSolar;

  /* convert to physical cgs units, is already physical at this point */
  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
  gs.log_HydrogenNumberDensity = log10(gs.XH * rho / PROTONMASS);
#ifdef GFM_LAMBDA
  gs.partindex = -1;
#endif
}
#endif

#ifdef EXPLICIT_COOLING
/** \brief Compute the new internal energy per unit mass with the Heun scheme.
 *
 *   The function uses the Heun scheme to compute the new internal energy per unit mass of the gas.
 *   Arguments are passed in code units.
 *
 *   \param u_old the initial (before cooling is applied) internal energy per unit mass of the gas cell
 *   \param rho   the proper density of the gas cell
 *   \param dt    the duration of the time step
 *   \param ne_guess electron number density relative to hydrogen number density (for molecular weight computation)
 *   \return the new internal energy per unit mass of the gas cell
 */
inline double DoCoolingHeun(double u_old, double rho, double dt, double *ne_guess)
{
  double u, u_prime, nHcgs;
  double ratefact, LambdaNet, LambdaNetPrime;

  nHcgs    = gs.XH * rho / PROTONMASS; /* hydrogen number dens in cgs units */
  ratefact = nHcgs * nHcgs / rho;

  LambdaNet = CoolingRateFromU(u_old, rho, ne_guess);
  u_prime   = u_old + ratefact * LambdaNet * dt;

  if(u_prime < All.MinEgySpec)
    return 0.0; /* in this way implicit solver is activated */

  LambdaNetPrime = CoolingRateFromU(u_prime, rho, ne_guess);
  u              = u_old + 0.5 * ratefact * (LambdaNet + LambdaNetPrime) * dt;

  return u;
}
#endif

/*! \brief Computes the new internal energy per unit mass.
 *
 *  The function solves for the new internal energy per unit mass of the gas
 *  by integrating the equation for the internal energy with an implicit
 *  Euler scheme. The root of resulting non linear equation,
 *  which gives tnew internal energy, is found with the bisection method.
 *  Arguments are passed in code units.
 *
 *  \param[in] u_old the initial (before cooling is applied) internal energy
 *             per unit mass of the gas cell.
 *  \param[in] rho   the proper density of the gas cell.
 *  \param[in] dt    the duration of the time step.
 *  \param[in] ne_guess electron number density relative to hydrogen number
 *             density (for molecular weight computation).
 *
 *  \return The new internal energy per unit mass of the gas cell.
 */
double DoCooling(double u_old, double rho, double dt, double *ne_guess)
{
  double u, du;
  double u_lower, u_upper;
  double ratefact;
  double LambdaNet;

#ifdef EXPLICIT_COOLING
  double ne_guess_expl = *ne_guess;
#endif

  int iter = 0;

  DoCool.u_old_input    = u_old;
  DoCool.rho_input      = rho;
  DoCool.dt_input       = dt;
  DoCool.ne_guess_input = *ne_guess;

  if(!gsl_finite(u_old))
    terminate("invalid input: u_old=%g rho=%g dt=%g ne =%g All.MinEgySpec=%g\n", u_old, rho, dt, *ne_guess, All.MinEgySpec);

  if(u_old < 0 || rho < 0)
    terminate("invalid input: task=%d u_old=%g  rho=%g  dt=%g  All.MinEgySpec=%g\n", ThisTask, u_old, rho, dt, All.MinEgySpec);

#ifdef SMUGGLE_RADIATION_FEEDBACK
  if(dt <= 0.0 && All.Time > All.TimeBegin)
    return u_old;
#endif

  /* convert to physical cgs units */
  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
  u_old *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
  dt *= All.UnitTime_in_s / All.HubbleParam;

#ifdef EXPLICIT_COOLING
  u = DoCoolingHeun(u_old, rho, dt, &ne_guess_expl);

  if(fabs(1.0 - u / u_old) < COOLING_TOLERANCE) /* use u given by the explicit scheme */
    {
      *ne_guess = ne_guess_expl;
      u *= All.UnitDensity_in_cgs / All.UnitPressure_in_cgs; /* to internal units */
      return u;
    }
#endif

  gs.nHcgs = gs.XH * rho / PROTONMASS; /* hydrogen number dens in cgs units */
  ratefact = gs.nHcgs * gs.nHcgs / rho;

  u       = u_old;
  u_lower = u;
  u_upper = u;

  LambdaNet = CoolingRateFromU(u, rho, ne_guess);

  /* bracketing */

  if(u - u_old - ratefact * LambdaNet * dt < 0) /* heating */
    {
      u_upper *= sqrt(1.1);
      u_lower /= sqrt(1.1);
      while(u_upper - u_old - ratefact * CoolingRateFromU(u_upper, rho, ne_guess) * dt < 0)
        {
          u_upper *= 1.1;
          u_lower *= 1.1;
        }
    }

  if(u - u_old - ratefact * LambdaNet * dt > 0)
    {
      u_lower /= sqrt(1.1);
      u_upper *= sqrt(1.1);
      while(u_lower - u_old - ratefact * CoolingRateFromU(u_lower, rho, ne_guess) * dt > 0)
        {
          u_upper /= 1.1;
          u_lower /= 1.1;
        }
    }

  do
    {
      u = 0.5 * (u_lower + u_upper);

      LambdaNet = CoolingRateFromU(u, rho, ne_guess);

      if(u - u_old - ratefact * LambdaNet * dt > 0)
        {
          u_upper = u;
        }
      else
        {
          u_lower = u;
        }

      du = u_upper - u_lower;

      iter++;

      if(iter >= (MAXITER - 10))
        printf("u= %g\n", u);
    }
  while(fabs(du / u) > 1.0e-6 && iter < MAXITER);

  if(iter >= MAXITER)
    terminate(
        "failed to converge in DoCooling(): DoCool.u_old_input=%g\nDoCool.rho_input= %g\nDoCool.dt_input= %g\nDoCool.ne_guess_input= "
        "%g\n",
        DoCool.u_old_input, DoCool.rho_input, DoCool.dt_input, DoCool.ne_guess_input);

  u *= All.UnitDensity_in_cgs / All.UnitPressure_in_cgs; /* to internal units */

  return u;
}

/*! \brief Returns the cooling time.
 *
 *  If we actually have heating, a cooling time of 0 is returned.
 *
 *  \param[in] u_old The initial (before cooling is applied) internal energy
 *             per unit mass of the gas cell.
 *  \param[in] rho The proper density of the gas cell.
 *  \param[in] ne_guess Electron number density relative to hydrogen number
 *             density (for molecular weight computation).
 *
 *  \return Cooling time; 0 if heating.
 */
double GetCoolingTime(double u_old, double rho, double *ne_guess)
{
  double u;
  double ratefact;
  double LambdaNet, coolingtime;

  DoCool.u_old_input    = u_old;
  DoCool.rho_input      = rho;
  DoCool.ne_guess_input = *ne_guess;

  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam; /* convert to physical cgs units */
  u_old *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;

  gs.nHcgs = gs.XH * rho / PROTONMASS; /* hydrogen number dens in cgs units */
  ratefact = gs.nHcgs * gs.nHcgs / rho;

  u = u_old;

  LambdaNet = CoolingRateFromU(u, rho, ne_guess);
  /* bracketing */

  if(LambdaNet >= 0) /* ups, we have actually heating due to UV background */
    return 0;

  coolingtime = u_old / (-ratefact * LambdaNet);

  coolingtime *= All.HubbleParam / All.UnitTime_in_s;

  return coolingtime;
}

/*! \brief Compute gas temperature from internal energy per unit mass.
 *
 *   This function determines the electron fraction, and hence the mean
 *   molecular weight. With it arrives at a self-consistent temperature.
 *   Element abundances and the rates for the emission are also computed.
 *
 *  \param[in] u   internal energy per unit mass.
 *  \param[in] rho gas density.
 *  \param[in, out] ne_guess electron number density relative to hydrogen
 *                  number density
 *
 *  \return The gas temperature.
 */
double convert_u_to_temp(double u, double rho, double *ne_guess)
{
  double temp, temp_old, temp_new, max = 0, ne_old;
  double mu;
  int iter = 0;

  double u_input, rho_input, ne_input;

  u_input   = u;
  rho_input = rho;
  ne_input  = *ne_guess;

  mu   = (1 + 4 * gs.yhelium) / (1 + gs.yhelium + *ne_guess); /* NOTE: slight inconsistency for GFM_COOLING_METAL */
  temp = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;

  do
    {
      ne_old = *ne_guess;

      if(!gsl_finite(log10(temp)))
        terminate("task=%d, iter=%d, u=%g, mu=%g, helium=%g, temp=%g, rho=%g, ne_guess=%g, temp_old=%g, temp_new=%g, max=%g", ThisTask,
                  iter, u, mu, gs.yhelium, temp, rho, *ne_guess, temp_old, temp_new, max);
      find_abundances_and_rates(log10(temp), rho, ne_guess);
      temp_old = temp;

      mu = (1 + 4 * gs.yhelium) / (1 + gs.yhelium + *ne_guess); /* NOTE: slight inconsistency for GFM_COOLING_METAL */

      temp_new = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;

      max = fmax(max,
                 temp_new / (1 + gs.yhelium + *ne_guess) *
                     fabs((*ne_guess - ne_old) / (temp_new - temp_old + 1.0))); /* NOTE: slight inconsistency for GFM_COOLING_METAL */

      temp = temp_old + (temp_new - temp_old) / (1 + max);
      iter++;

      if(iter > (MAXITER - 10))
        printf("-> temp= %g ne=%g\n", temp, *ne_guess);
    }
  while(fabs(temp - temp_old) > 1.0e-3 * temp && iter < MAXITER);

  if(iter >= MAXITER)
    {
      printf("failed to converge in convert_u_to_temp()\n");
      printf("u_input= %g\nrho_input=%g\n ne_input=%g\n", u_input, rho_input, ne_input);
      printf("DoCool.u_old_input=%g\nDoCool.rho_input= %g\nDoCool.dt_input= %g\nDoCool.ne_guess_input= %g\n", DoCool.u_old_input,
             DoCool.rho_input, DoCool.dt_input, DoCool.ne_guess_input);
#ifdef GFM_COOLING_METAL
      printf("log Z=%g\nlog nH=%g\nH mass frac=%g\n", gs.log_MetallicityInSolar * GFM_SOLAR_METALLICITY, gs.log_HydrogenNumberDensity,
             gs.XH);
#endif
      terminate("convergence failure");
    }
  gs.mu = mu;

  return temp;
}

/*! \brief Computes the actual abundance ratios.
 *
 *  The chemical composition of the gas is primordial (no metals are present).
 *
 *  \param[in] logT log10 of gas temperature.
 *  \param[in] rho Gas density.
 *  \param[in, out] ne_guess Electron number density relative to hydrogen
 *                  number density.
 *
 *  \return void
 */
void find_abundances_and_rates(double logT, double rho, double *ne_guess)
{
  double neold, nenew;
  int j, niter;
  double flow, fhi, t;

  double logT_input, rho_input, ne_input;

  logT_input = logT;
  rho_input  = rho;
  ne_input   = *ne_guess;

  if(!gsl_finite(logT))
    terminate("logT=%g, rho=%g, ne_guess=%g\n", logT, rho, *ne_guess);

  if(logT <= Tmin) /* everything neutral */
    {
      gs.nH0    = 1.0;
      gs.nHe0   = gs.yhelium;
      gs.nHp    = 0;
      gs.nHep   = 0;
      gs.nHepp  = 0;
      gs.ne     = 0;
      *ne_guess = 0;
      return;
    }

  if(logT >= Tmax) /* everything is ionized */
    {
      gs.nH0    = 0;
      gs.nHe0   = 0;
      gs.nHp    = 1.0;
      gs.nHep   = 0;
      gs.nHepp  = gs.yhelium;
      gs.ne     = gs.nHp + 2.0 * gs.nHepp;
      *ne_guess = gs.ne; /* note: in units of the hydrogen number density */
      return;
    }

  t    = (logT - Tmin) / deltaT;
  j    = (int)t;
  fhi  = t - j;
  flow = 1 - fhi;

  if(*ne_guess == 0)
    *ne_guess = 1.0;

  gs.nHcgs = gs.XH * rho / PROTONMASS; /* hydrogen number dens in cgs units */

  gs.ne    = *ne_guess;
  neold    = gs.ne;
  niter    = 0;
  gs.necgs = gs.ne * gs.nHcgs;

  /* evaluate number densities iteratively (cf KWH eqns 33-38) in units of nH */
  do
    {
      niter++;

      gs.aHp   = flow * RateT[j].AlphaHp + fhi * RateT[j + 1].AlphaHp;
      gs.aHep  = flow * RateT[j].AlphaHep + fhi * RateT[j + 1].AlphaHep;
      gs.aHepp = flow * RateT[j].AlphaHepp + fhi * RateT[j + 1].AlphaHepp;
      gs.ad    = flow * RateT[j].Alphad + fhi * RateT[j + 1].Alphad;
      gs.geH0  = flow * RateT[j].GammaeH0 + fhi * RateT[j + 1].GammaeH0;
      gs.geHe0 = flow * RateT[j].GammaeHe0 + fhi * RateT[j + 1].GammaeHe0;
      gs.geHep = flow * RateT[j].GammaeHep + fhi * RateT[j + 1].GammaeHep;

      if(gs.necgs <= 1.e-25 || pc.J_UV == 0)
        {
          gs.gJH0ne = gs.gJHe0ne = gs.gJHepne = 0;
        }
      else
        {
          gs.gJH0ne  = pc.gJH0 / gs.necgs;
          gs.gJHe0ne = pc.gJHe0 / gs.necgs;
          gs.gJHepne = pc.gJHep / gs.necgs;
        }

      gs.nH0 = gs.aHp / (gs.aHp + gs.geH0 + gs.gJH0ne); /* eqn (33) */
      gs.nHp = 1.0 - gs.nH0;                            /* eqn (34) */

      if((gs.gJHe0ne + gs.geHe0) <= SMALLNUM) /* no ionization at all */
        {
          gs.nHep  = 0.0;
          gs.nHepp = 0.0;
          gs.nHe0  = gs.yhelium;
        }
      else
        {
          gs.nHep =
              gs.yhelium / (1.0 + (gs.aHep + gs.ad) / (gs.geHe0 + gs.gJHe0ne) + (gs.geHep + gs.gJHepne) / gs.aHepp); /* eqn (35) */
          gs.nHe0  = gs.nHep * (gs.aHep + gs.ad) / (gs.geHe0 + gs.gJHe0ne);                                          /* eqn (36) */
          gs.nHepp = gs.nHep * (gs.geHep + gs.gJHepne) / gs.aHepp;                                                   /* eqn (37) */
        }

      neold = gs.ne;

      gs.ne    = gs.nHp + gs.nHep + 2 * gs.nHepp; /* eqn (38) */
      gs.necgs = gs.ne * gs.nHcgs;

      if(pc.J_UV == 0)
        break;

      nenew    = 0.5 * (gs.ne + neold);
      gs.ne    = nenew;
      gs.necgs = gs.ne * gs.nHcgs;

      if(fabs(gs.ne - neold) < 1.0e-4)
        break;

      if(niter > (MAXITER - 10))
        printf("ne= %g  niter=%d\n", gs.ne, niter);
    }
  while(niter < MAXITER);

  if(niter >= MAXITER)
    {
      printf("gs.aHp = %le\n", gs.aHp);
      char buff[MAXLEN_PATH];
      file_path_sprintf(buff, "%s/cooling_task%d.dat", All.OutputDir, ThisTask);
      FILE *fp = fopen(buff, "w");
      fwrite(&All.Time, sizeof(double), 1, fp);
      fwrite(&logT_input, sizeof(double), 1, fp);
      fwrite(&rho_input, sizeof(double), 1, fp);
      fwrite(&ne_input, sizeof(double), 1, fp);
      fclose(fp);
      terminate(
          "no convergence reached in find_abundances_and_rates(): logT_input= %g  rho_input= %g  ne_input= %g "
          "DoCool.u_old_input=%g\nDoCool.rho_input= %g\nDoCool.dt_input= %g\nDoCool.ne_guess_input= %g\n",
          logT_input, rho_input, ne_input, DoCool.u_old_input, DoCool.rho_input, DoCool.dt_input, DoCool.ne_guess_input);
    }
  gs.bH0  = flow * RateT[j].BetaH0 + fhi * RateT[j + 1].BetaH0;
  gs.bHep = flow * RateT[j].BetaHep + fhi * RateT[j + 1].BetaHep;
  gs.bff  = flow * RateT[j].Betaff + fhi * RateT[j + 1].Betaff;

  *ne_guess = gs.ne;
}

/*! \brief Get cooling rate from gas internal energy.
 *
 *  This function first computes the self-consistent temperature
 *  and abundance ratios, and then it calculates
 *  (heating rate-cooling rate)/n_h^2 in cgs units.
 *
 *  \param[in] u Gas internal energy per unit mass.
 *  \param[in] rho Gas density.
 *  \param[in, out] ne_guess Electron number density relative to hydrogen
 *                  number density.
 *
 *  \return Cooling rate.
 */
double CoolingRateFromU(double u, double rho, double *ne_guess)
{
  double temp;

  temp = convert_u_to_temp(u, rho, ne_guess);

  return CoolingRate(log10(temp), rho, ne_guess);
}

/*! \brief  This function computes the self-consistent temperature and
 *          abundance ratios.
 *
 *  Used only in io_fields.c for calculating output fields.
 *
 *  \param[in] i index into SphP for gas cell to consider.
 *  \param[in, out] ne_guess pointer to electron number density relative to
 *                  hydrogen number density (modified).
 *  \param[out] nH0 Pointer to the neutral hydrogen fraction (set to value in
 *              the GasState struct).
 *  \param[out] coolrate Pointer to cooling rate (set to value from
 *              CoolingRateFromU).
 *
 *  \return void
 */
void SetOutputGasState(int i, double *ne_guess, double *nH0, double *coolrate)
{
  double rho = SphP[i].Density * All.cf_a3inv;
  double u   = fmax(All.MinEgySpec, SphP[i].Utherm);
#if !(defined(SFR_MCS) && defined(GFM_COOLING_METAL))
  double sfr = 0;
  /* update GasState as appropriate given compile-time options and cell properties */
#if defined(USE_SFR) && (!defined(LOCAL_FEEDBACK) && !defined(SMUGGLE_STAR_FEEDBACK) && !defined(SFR_MCS))
  sfr = get_starformation_rate(i);
#endif

  if(sfr == 0)
    {
#if defined(GFM_AGN_RADIATION) || defined(GFM_UVB_CORRECTIONS) || defined(UVB_SELF_SHIELDING)
#ifdef GFM_AGN_RADIATION
      update_radiation_state(rho, get_hydrogen_abundances_of_local_cell(i), SphP[i].AGNBolIntensity);
#else

#ifdef GFM_STELLAR_EVOLUTION
      update_radiation_state(rho, get_hydrogen_abundances_of_local_cell(i), 0);
#else
      update_radiation_state(rho, HYDROGEN_MASSFRAC, 0);
#endif

#endif
#endif

#if defined(GFM_COOLING_METAL) && !defined(MRT_METAL_COOLING)
      update_gas_state(rho, get_hydrogen_abundances_of_local_cell(i), SphP[i].Metallicity);
#endif
    }
  else
    {
#if defined(GFM_AGN_RADIATION) || defined(GFM_UVB_CORRECTIONS) || defined(UVB_SELF_SHIELDING)
      update_radiation_state(rho, HYDROGEN_MASSFRAC, 0);
#endif
#if defined(GFM_COOLING_METAL) && !defined(MRT_METAL_COOLING)
      update_gas_state(rho, GFM_INITIAL_ABUNDANCE_HYDROGEN, 0.0);
#endif
    }
#else
  /*SFR_MCS and GFM_COOLING_METAL are defined in this block */
#ifdef UVB_SELF_SHIELDING
  update_radiation_state(rho, SphP[i].Metallicity, &pc);
#endif
  update_gas_state(rho, SphP[i].Metallicity, &gs);
#endif /* !(defined(SFR_MCS) && defined(GFM_COOLING_METAL))*/
  /* update DoCool */
  DoCool.u_old_input    = u;
  DoCool.rho_input      = rho;
  DoCool.ne_guess_input = *ne_guess;

  /* convert to physical cgs units */
  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
  u *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;

  /* calculate cooling rate (and so ne_guess and all of gs including nH0, nHeII) */
  *coolrate = CoolingRateFromU(u, rho, ne_guess);

  *nH0 = gs.nH0;
}

void getHeState(double *nHe0, double *nHep, double *nHepp)
{
  *nHe0  = gs.nHe0;
  *nHep  = gs.nHep;
  *nHepp = gs.nHepp;
}

/*! \brief  Calculate (heating rate-cooling rate)/n_h^2 in cgs units.
 *
 *  \param[in] logT log10 of gas temperature.
 *  \param[in] rho Gas density.
 *  \param[in, out] nelec Electron number density relative to hydrogen number
 *                  density.
 *
 *  \return (heating rate-cooling rate)/n_h^2.
 */
double CoolingRate(double logT, double rho, double *nelec)
{
  double Lambda, Heat;
  double LambdaExc, LambdaIon, LambdaRec, LambdaFF, LambdaCmptn = 0.0;
  double LambdaExcH0, LambdaExcHep, LambdaIonH0, LambdaIonHe0, LambdaIonHep;
  double LambdaRecHp, LambdaRecHep, LambdaRecHepp, LambdaRecHepd;
  double redshift;
  double T;
  double LambdaPrim = 0.0, LambdaMet = 0.0, LambdaDust = 0.0, LambdaMol = 0.0;

  if(logT <= Tmin)
    logT = Tmin + 0.5 * deltaT; /* floor at Tmin */

  gs.nHcgs = gs.XH * rho / PROTONMASS; /* hydrogen number dens in cgs units */

  if(logT < Tmax)
    {
      find_abundances_and_rates(logT, rho, nelec);

      /* Compute cooling and heating rate (cf KWH Table 1) in units of nH**2 */
      T = pow(10.0, logT);

      LambdaExcH0   = gs.bH0 * gs.ne * gs.nH0;
      LambdaExcHep  = gs.bHep * gs.ne * gs.nHep;
      LambdaExc     = LambdaExcH0 + LambdaExcHep; /* excitation */
      LambdaIonH0   = 2.18e-11 * gs.geH0 * gs.ne * gs.nH0;
      LambdaIonHe0  = 3.94e-11 * gs.geHe0 * gs.ne * gs.nHe0;
      LambdaIonHep  = 8.72e-11 * gs.geHep * gs.ne * gs.nHep;
      LambdaIon     = LambdaIonH0 + LambdaIonHe0 + LambdaIonHep; /* ionization */
      LambdaRecHp   = 1.036e-16 * T * gs.ne * (gs.aHp * gs.nHp);
      LambdaRecHep  = 1.036e-16 * T * gs.ne * (gs.aHep * gs.nHep);
      LambdaRecHepp = 1.036e-16 * T * gs.ne * (gs.aHepp * gs.nHepp);
      LambdaRecHepd = 6.526e-11 * gs.ad * gs.ne * gs.nHep;
      LambdaRec     = LambdaRecHp + LambdaRecHep + LambdaRecHepp + LambdaRecHepd;
      LambdaFF      = gs.bff * (gs.nHp + gs.nHep + 4 * gs.nHepp) * gs.ne;
      LambdaPrim    = LambdaExc + LambdaIon + LambdaRec + LambdaFF;

#if defined(GFM_COOLING_METAL) && !defined(MRT_METAL_COOLING)
      if(logT > LogMinMetalTemp)
        {
#ifdef GFM_AGN_RADIATION
          LambdaMet = get_CoolingMetalRate(gs.LogAGNBolIntensity, gs.log_MetallicityInSolar, gs.log_HydrogenNumberDensity, logT);
#else
#ifdef RADCOOL
          LambdaMet = get_CoolingMetalRate(gs.LogPhios, gs.LogPhins,
#ifdef RADCOOL_HOTHALO
                                           gs.LogPhiT6, gs.LogPhiT7, gs.LogPhiT8,
#endif
                                           gs.log_MetallicityInSolar, gs.log_HydrogenNumberDensity, logT);
#else
          LambdaMet = get_CoolingMetalRate(gs.log_MetallicityInSolar, gs.log_HydrogenNumberDensity, logT);
#endif
#endif
        }
#endif

#ifdef GFM_DUST_COOLING
      LambdaDust = get_CoolingDustRate(logT, gs.dust_to_gas_ratio, gs.ne, gs.dust_cool);
#endif

#ifdef SMUGGLE_DUST_HEATING_COOLING
      LambdaDust = get_DustHeatingCoolingRate(logT, gs.XH, gs.log_MetallicityInSolar);
#endif

      if(All.ComovingIntegrationOn)
        {
          redshift    = 1 / All.Time - 1;
          LambdaCmptn = 5.65e-36 * gs.ne * (T - 2.73 * (1. + redshift)) * pow(1. + redshift, 4.) / gs.nHcgs;
        }
      else
        LambdaCmptn = 0;

#ifdef SMUGGLE_MOLEC_COOLING
      LambdaMol = get_MolecularCoolingRate(logT, gs.nHcgs, gs.log_MetallicityInSolar, pc.J_UV);
#endif

#ifdef GFM_LAMBDA
      if(gs.partindex >= 0)
        {
          SphP[gs.partindex].Lambda[0] = LambdaPrim;
          SphP[gs.partindex].Lambda[1] = LambdaMet;
          SphP[gs.partindex].Lambda[2] = LambdaDust;
          SphP[gs.partindex].Lambda[3] = LambdaCmptn;
          SphP[gs.partindex].Lambda[4] = LambdaMol;
          SphP[gs.partindex].Lambda[5] = 0.0;
        }
#endif
      Lambda = LambdaPrim + LambdaMet + LambdaDust + LambdaCmptn + LambdaMol;

      Heat = 0;
      if(pc.J_UV != 0)
#ifdef GFM_UVB_CORRECTIONS
        Heat += All.UV_HeatBoost * (gs.nH0 * pc.epsH0 + gs.nHe0 * pc.epsHe0 + gs.UV_HeII_Factor * gs.nHep * pc.epsHep) / gs.nHcgs;
#else
        Heat += (gs.nH0 * pc.epsH0 + gs.nHe0 * pc.epsHe0 + gs.nHep * pc.epsHep) / gs.nHcgs;
#endif

/* see Hopkins+ 2018 Appendix B */
#ifdef SMUGGLE_COSMIC_RAY_HEATING
      double fac         = 0.0;
      double MeanBarDens = gs.RhoGasMean * All.cf_a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;

      if(All.ComovingIntegrationOn)
        {
          if(rho > 1000. * MeanBarDens)
            fac = gs.nHcgs / (0.01 + gs.nHcgs);
        }
      else
        fac = fmin(1., 1e2 * gs.nHcgs);

      Heat += 9.e-28 * 2.63 * gs.ne * fac / gs.nHcgs;
#endif

#ifdef SMUGGLE_PHOTOELECTRIC_HEATING
      double epe = 1.7;
      double xpe = 2 * sqrt(T) * epe / (gs.nHcgs * gs.ne);
      Heat += 1.3e-24 * epe * pow(10., gs.log_MetallicityInSolar) *
              (0.049 / (1. + pow(xpe / 1925, 0.73)) + 0.037 * pow(1e-4 * T, 0.7) / (1. + (xpe / 5000.))) / gs.nHcgs;
#endif
    }
  else /* here we're outside of tabulated rates, T>Tmax K */
    {
      /* at high T (fully ionized); only free-free and Compton cooling are present. Assumes no heating. */

      Heat = 0;

      LambdaExcH0 = LambdaExcHep = LambdaIonH0 = LambdaIonHe0 = LambdaIonHep = LambdaRecHp = LambdaRecHep = LambdaRecHepp =
          LambdaRecHepd                                                                                   = 0;

      /* very hot: H and He both fully ionized */
      gs.nHp   = 1.0;
      gs.nHep  = 0;
      gs.nHepp = gs.yhelium;
      gs.ne    = gs.nHp + 2.0 * gs.nHepp;
      *nelec   = gs.ne; /* note: in units of the hydrogen number density */

      T        = pow(10.0, logT);
      LambdaFF = 1.42e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - logT) * (5.5 - logT) / 3)) * (gs.nHp + 4 * gs.nHepp) * gs.ne;

      if(All.ComovingIntegrationOn)
        {
          redshift = 1 / All.Time - 1;
          /* add inverse Compton cooling off the microwave background */
          LambdaCmptn = 5.65e-36 * gs.ne * (T - 2.73 * (1. + redshift)) * pow(1. + redshift, 4.) / gs.nHcgs;
        }
      else
        LambdaCmptn = 0;

#ifdef GFM_DUST_COOLING
      LambdaDust = get_CoolingDustRate(logT, gs.dust_to_gas_ratio, gs.ne, gs.dust_cool);
#endif

#ifdef SMUGGLE_DUST_HEATING_COOLING
      LambdaDust = get_DustHeatingCoolingRate(logT, gs.XH, gs.log_MetallicityInSolar);
#endif

#ifdef GFM_LAMBDA
      if(gs.partindex >= 0)
        {
          SphP[gs.partindex].Lambda[0] = 0.0;
          SphP[gs.partindex].Lambda[1] = 0.0;
          SphP[gs.partindex].Lambda[2] = LambdaDust;
          SphP[gs.partindex].Lambda[3] = LambdaCmptn;
          SphP[gs.partindex].Lambda[4] = 0.0;
          SphP[gs.partindex].Lambda[5] = LambdaFF;
        }
#endif

      Lambda = LambdaFF + LambdaCmptn;
    }

  return (Heat - Lambda);
}

/** \brief  Compute the log10 of the gas temperature.
 *
 *  the function is unused
 *
 *  \param u     gas internal energy per unit mass
 *  \param ne    electron number density relative to hydrogen number density
 *  \return      log10 of gas temperature
 */
double LogTemp(double u, double ne) /* ne= electron density in terms of hydrogen density */
{
  double T;

  if(u < gs.ethmin)
    u = gs.ethmin;

  T = log10(GAMMA_MINUS1 * u * gs.mhboltz * (1 + 4 * gs.yhelium) /
            (1 + ne + gs.yhelium)); /* NOTE: slight inconsistency for GFM_COOLING_METAL */

  return T;
}

/*! \brief Make cooling rates interpolation table.
 *
 *  Set up interpolation tables in T for cooling rates given in
 *  KWH, ApJS, 105, 19.
 *
 *  \return void
 */
void MakeRateTable(void)
{
  int i;
  double T;
  double Tfact;
#ifdef GFM_PRIMORDIAL_RATES
  double VF96a, VF96b, VF96T0, VF96T1, VF96T;
  double Vor97dE, Vor97P, Vor97A, Vor97X, Vor97K, Vor97U, Vor97T_eV;
#endif

  gs.yhelium = (1 - gs.XH) / (4 * gs.XH);
  gs.mhboltz = PROTONMASS / BOLTZMANN;
  if(All.MinGasTemp > 0.0)
    Tmin = log10(0.1 * All.MinGasTemp);
  else
    Tmin = 1.0;
  deltaT    = (Tmax - Tmin) / NCOOLTAB;
  gs.ethmin = pow(10.0, Tmin) * (1. + gs.yhelium) / ((1. + 4. * gs.yhelium) * gs.mhboltz * GAMMA_MINUS1);
  /* minimum internal energy for neutral gas */

  for(i = 0; i <= NCOOLTAB; i++)
    {
      RateT[i].BetaH0 = RateT[i].BetaHep = RateT[i].Betaff = RateT[i].AlphaHp = RateT[i].AlphaHep = RateT[i].AlphaHepp =
          RateT[i].Alphad = RateT[i].GammaeH0 = RateT[i].GammaeHe0 = RateT[i].GammaeHep = 0;

      T     = pow(10.0, Tmin + deltaT * i);
      Tfact = 1.0 / (1 + sqrt(T / 1.0e5));

      /* collisional excitation */
      /* Cen 1992 */
      if(118348 / T < 70)
        RateT[i].BetaH0 = 7.5e-19 * exp(-118348 / T) * Tfact;
      if(473638 / T < 70)
        RateT[i].BetaHep = 5.54e-17 * pow(T, -0.397) * exp(-473638 / T) * Tfact;

      /* free-free */
      RateT[i].Betaff = 1.43e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - log10(T)) * (5.5 - log10(T)) / 3));

      /* recombination */
#ifdef GFM_PRIMORDIAL_RATES
      /* Verner & Ferland 1996, eq.4 */
      /* Hydrogen II */
      VF96a            = 7.982e-11;
      VF96b            = 7.480e-01;
      VF96T0           = 3.148e+00;
      VF96T1           = 7.036e+05;
      VF96T            = sqrt(T / VF96T0);
      RateT[i].AlphaHp = VF96a / (VF96T * pow(1.0 + VF96T, 1.0 - VF96b) * pow(1.0 + sqrt(T / VF96T1), 1.0 + VF96b));
      /* Helium II */
      VF96a             = 9.356e-10;
      VF96b             = 7.892e-01;
      VF96T0            = 4.266e-02;
      VF96T1            = 4.677e+06;
      VF96T             = sqrt(T / VF96T0);
      RateT[i].AlphaHep = VF96a / (VF96T * pow(1.0 + VF96T, 1.0 - VF96b) * pow(1.0 + sqrt(T / VF96T1), 1.0 + VF96b));
      /* Helium III */
      VF96a              = 1.891e-10;
      VF96b              = 7.524e-01;
      VF96T0             = 9.370e+00;
      VF96T1             = 2.774e+06;
      VF96T              = sqrt(T / VF96T0);
      RateT[i].AlphaHepp = VF96a / (VF96T * pow(1.0 + VF96T, 1.0 - VF96b) * pow(1.0 + sqrt(T / VF96T1), 1.0 + VF96b));
#else
      /* Cen 1992 */
      /* Hydrogen II */
      RateT[i].AlphaHp = 8.4e-11 * pow(T / 1000, -0.2) / (1. + pow(T / 1.0e6, 0.7)) / sqrt(T);
      /* Helium II */
      RateT[i].AlphaHep = 1.5e-10 * pow(T, -0.6353);
      /* Helium III */
      RateT[i].AlphaHepp = 4. * RateT[i].AlphaHp;
#endif
      /* Cen 1992 */
      /* dielectric recombination */
      if(470000 / T < 70)
        RateT[i].Alphad = 1.9e-3 * pow(T, -1.5) * exp(-470000 / T) * (1. + 0.3 * exp(-94000 / T));

        /* collisional ionization */
#ifdef GFM_PRIMORDIAL_RATES
      /* Voronov 1997 */
      Vor97T_eV = T / eV_to_K;
      /* Hydrogen */
      Vor97dE           = 13.6;
      Vor97P            = 0.0;
      Vor97A            = 0.291e-7;
      Vor97X            = 0.232;
      Vor97K            = 0.39;
      Vor97U            = Vor97dE / Vor97T_eV;
      RateT[i].GammaeH0 = Vor97A * (1.0 + Vor97P * sqrt(Vor97U)) * pow(Vor97U, Vor97K) * exp(-Vor97U) / (Vor97X + Vor97U);
      /* Helium */
      Vor97dE            = 24.6;
      Vor97P             = 0.0;
      Vor97A             = 0.175e-7;
      Vor97X             = 0.18;
      Vor97K             = 0.35;
      Vor97U             = Vor97dE / Vor97T_eV;
      RateT[i].GammaeHe0 = Vor97A * (1.0 + Vor97P * sqrt(Vor97U)) * pow(Vor97U, Vor97K) * exp(-Vor97U) / (Vor97X + Vor97U);
      /* Hellium II */
      Vor97dE            = 54.4;
      Vor97P             = 1.0;
      Vor97A             = 0.205e-8;
      Vor97X             = 0.265;
      Vor97K             = 0.25;
      Vor97U             = Vor97dE / Vor97T_eV;
      RateT[i].GammaeHep = Vor97A * (1.0 + Vor97P * sqrt(Vor97U)) * pow(Vor97U, Vor97K) * exp(-Vor97U) / (Vor97X + Vor97U);
#else
      /* Cen 1992 */
      /* Hydrogen */
      if(157809.1 / T < 70)
        RateT[i].GammaeH0 = 5.85e-11 * sqrt(T) * exp(-157809.1 / T) * Tfact;
      /* Helium */
      if(285335.4 / T < 70)
        RateT[i].GammaeHe0 = 2.38e-11 * sqrt(T) * exp(-285335.4 / T) * Tfact;
      /* Hellium II */
      if(631515.0 / T < 70)
        RateT[i].GammaeHep = 5.68e-12 * sqrt(T) * exp(-631515.0 / T) * Tfact;
#endif
    }
}

/*! \brief Read table input for ionizing parameters.
 *
 *  \param[in] fname Name of file that contains the tabulated parameters.
 *  \param[in] which Flag used to identify the type of the ionizing background
 *                   (0 = UV background, 1 = AGN background, 2=RADCOOL).
 *
 *  \return void
 */
void ReadIonizeParams(char *fname, int which)
{
  int iter, i;
  FILE *fdcool;
  float dummy;

  if(which == 0)
    {
      NheattabUVB = 0;

      for(iter = 0, i = 0; iter < 2; iter++)
        {
          if(!(fdcool = fopen(fname, "r")))
            terminate("COOLING: cannot read ionization table in file `%s'\n", fname);
          if(iter == 0)
            while(fscanf(fdcool, "%g %g %g %g %g %g %g", &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy) != EOF)
              NheattabUVB++;
          if(iter == 1)
            while(fscanf(fdcool, "%g %g %g %g %g %g %g", &PhotoTUVB[i].variable, &PhotoTUVB[i].gH0, &PhotoTUVB[i].gHe,
                         &PhotoTUVB[i].gHep, &PhotoTUVB[i].eH0, &PhotoTUVB[i].eHe, &PhotoTUVB[i].eHep) != EOF)
              i++;
          fclose(fdcool);

          if(iter == 0)
            {
              PhotoTUVB = (PhotoTable *)mymalloc("PhotoT", NheattabUVB * sizeof(PhotoTable));
              mpi_printf("COOLING: read ionization table with %d entries in file `%s'.\n", NheattabUVB, fname);
            }
        }
      /* ignore zeros at end of treecool file */
      for(i = 0; i < NheattabUVB; ++i)
        if(PhotoTUVB[i].gH0 == 0.0)
          break;

      NheattabUVB = i;
      mpi_printf("COOLING: using %d ionization table entries from file `%s'.\n", NheattabUVB, fname);
    }
#ifdef GFM_AGN_RADIATION
  if(which == 1)
    {
      PhotoTAGN = (PhotoTable *)mymalloc("PhotoTAGN", 1 * sizeof(PhotoTable));
      if(!(fdcool = fopen(fname, "r")))
        terminate("COOLING: cannot read ionization table in file `%s'\n", fname);
      fscanf(fdcool, "%g %g %g %g %g %g %g", &PhotoTAGN[0].variable, &PhotoTAGN[0].gH0, &PhotoTAGN[0].gHe, &PhotoTAGN[0].gHep,
             &PhotoTAGN[0].eH0, &PhotoTAGN[0].eHe, &PhotoTAGN[0].eHep);
      fclose(fdcool);
      mpi_printf("COOLING/GFM_AGN_RADIATION read `%s': %g %g %g %g %g %g %g\n", fname, PhotoTAGN[0].variable, PhotoTAGN[0].gH0,
                 PhotoTAGN[0].gHe, PhotoTAGN[0].gHep, PhotoTAGN[0].eH0, PhotoTAGN[0].eHe, PhotoTAGN[0].eHep);
    }
#endif
#ifdef RADCOOL
  if(which == 2)
    {
      PhotoTRAD_OS = (PhotoTable *)mymalloc("PhotoTRAD_OS", 1 * sizeof(PhotoTable));
      PhotoTRAD_NS = (PhotoTable *)mymalloc("PhotoTRAD_NS", 1 * sizeof(PhotoTable));
#ifdef RADCOOL_HOTHALO
      PhotoTRAD_T6 = (PhotoTable *)mymalloc("PhotoTRAD_T6", 1 * sizeof(PhotoTable));
      PhotoTRAD_T7 = (PhotoTable *)mymalloc("PhotoTRAD_T7", 1 * sizeof(PhotoTable));
      PhotoTRAD_T8 = (PhotoTable *)mymalloc("PhotoTRAD_T8", 1 * sizeof(PhotoTable));
#endif
      if(!(fdcool = fopen(fname, "r")))
        terminate("COOLING: cannot read ionization table in file `%s'\n", fname);
      fscanf(fdcool, "%g %g %g %g %g %g %g \n", &PhotoTRAD_OS[0].variable, &PhotoTRAD_OS[0].gH0, &PhotoTRAD_OS[0].gHe,
             &PhotoTRAD_OS[0].gHep, &PhotoTRAD_OS[0].eH0, &PhotoTRAD_OS[0].eHe, &PhotoTRAD_OS[0].eHep);
      fscanf(fdcool, "%g %g %g %g %g %g %g \n", &PhotoTRAD_NS[0].variable, &PhotoTRAD_NS[0].gH0, &PhotoTRAD_NS[0].gHe,
             &PhotoTRAD_NS[0].gHep, &PhotoTRAD_NS[0].eH0, &PhotoTRAD_NS[0].eHe, &PhotoTRAD_NS[0].eHep);
#ifdef RADCOOL_HOTHALO
      fscanf(fdcool, "%g %g %g %g %g %g %g \n", &PhotoTRAD_T6[0].variable, &PhotoTRAD_T6[0].gH0, &PhotoTRAD_T6[0].gHe,
             &PhotoTRAD_T6[0].gHep, &PhotoTRAD_T6[0].eH0, &PhotoTRAD_T6[0].eHe, &PhotoTRAD_T6[0].eHep);
      fscanf(fdcool, "%g %g %g %g %g %g %g \n", &PhotoTRAD_T7[0].variable, &PhotoTRAD_T7[0].gH0, &PhotoTRAD_T7[0].gHe,
             &PhotoTRAD_T7[0].gHep, &PhotoTRAD_T7[0].eH0, &PhotoTRAD_T7[0].eHe, &PhotoTRAD_T7[0].eHep);
      fscanf(fdcool, "%g %g %g %g %g %g %g \n", &PhotoTRAD_T8[0].variable, &PhotoTRAD_T8[0].gH0, &PhotoTRAD_T8[0].gHe,
             &PhotoTRAD_T8[0].gHep, &PhotoTRAD_T8[0].eH0, &PhotoTRAD_T8[0].eHe, &PhotoTRAD_T8[0].eHep);
#endif
      fclose(fdcool);
      mpi_printf("COOLING/RADCOOL_OLD_STARS read `%s': %g %g %g %g %g %g %g\n", fname, PhotoTRAD_OS[0].variable, PhotoTRAD_OS[0].gH0,
                 PhotoTRAD_OS[0].gHe, PhotoTRAD_OS[0].gHep, PhotoTRAD_OS[0].eH0, PhotoTRAD_OS[0].eHe, PhotoTRAD_OS[0].eHep);
      mpi_printf("COOLING/RADCOOL_NEW_STARS read `%s': %g %g %g %g %g %g %g\n", fname, PhotoTRAD_NS[0].variable, PhotoTRAD_NS[0].gH0,
                 PhotoTRAD_NS[0].gHe, PhotoTRAD_NS[0].gHep, PhotoTRAD_NS[0].eH0, PhotoTRAD_NS[0].eHe, PhotoTRAD_NS[0].eHep);
#ifdef RADCOOL_HOTHALO
      mpi_printf("COOLING/RADCOOL_HOTHALO_T6 read `%s': %g %g %g %g %g %g %g\n", fname, PhotoTRAD_T6[0].variable, PhotoTRAD_T6[0].gH0,
                 PhotoTRAD_T6[0].gHe, PhotoTRAD_T6[0].gHep, PhotoTRAD_T6[0].eH0, PhotoTRAD_T6[0].eHe, PhotoTRAD_T6[0].eHep);
      mpi_printf("COOLING/RADCOOL_HOTHALO_T7 read `%s': %g %g %g %g %g %g %g\n", fname, PhotoTRAD_T7[0].variable, PhotoTRAD_T7[0].gH0,
                 PhotoTRAD_T7[0].gHe, PhotoTRAD_T7[0].gHep, PhotoTRAD_T7[0].eH0, PhotoTRAD_T7[0].eHe, PhotoTRAD_T7[0].eHep);
      mpi_printf("COOLING/RADCOOL_HOTHALO_T8 read `%s': %g %g %g %g %g %g %g\n", fname, PhotoTRAD_T8[0].variable, PhotoTRAD_T8[0].gH0,
                 PhotoTRAD_T8[0].gHe, PhotoTRAD_T8[0].gHep, PhotoTRAD_T8[0].eH0, PhotoTRAD_T8[0].eHe, PhotoTRAD_T8[0].eHep);
#endif
    }
#endif
}

#ifdef UVB_SELF_SHIELDING
/** \brief Read table input for self shielding fitting formula parameters.
 *
 *  \param file that contains the tabulated parameters
 */
void ReadSelfShieldingParams(char *fname)
{
  int iter, i;
  FILE *fdselfshielding;
  float dummy;

  Nselfshieldingtab = 0;

  for(iter = 0, i = 0; iter < 2; iter++)
    {
      if(!(fdselfshielding = fopen(fname, "r")))
        terminate("UVB_SELF_SHIELDING: cannot read self shielding table in file `%s'\n", fname);
      if(iter == 0)
        while(fscanf(fdselfshielding, "%g %g %g %g %g %g %g", &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy) != EOF)
          Nselfshieldingtab++;
      if(iter == 1)
        while(fscanf(fdselfshielding, "%g %g %g %g %g %g %g", &SelfShieldingParams[i].redshift, &SelfShieldingParams[i].alpha1,
                     &SelfShieldingParams[i].alpha2, &SelfShieldingParams[i].beta, &SelfShieldingParams[i].xi,
                     &SelfShieldingParams[i].n0, &SelfShieldingParams[i].f) != EOF)
          i++;
      fclose(fdselfshielding);

      if(iter == 0)
        {
          SelfShieldingParams = (SelfShieldingTable *)mymalloc("SelfShieldingT", Nselfshieldingtab * sizeof(SelfShieldingTable));
          mpi_printf("UVB_SELF_SHIELDING: read self shielding table with %d entries in file `%s'.\n", Nselfshieldingtab, fname);
        }
    }
  /* ignore zeros at end of self shielding parameters file */
  for(i = 0; i < Nselfshieldingtab; ++i)
    if(SelfShieldingParams[i].alpha1 == 0.0)
      break;
  Nselfshieldingtab = i;
  mpi_printf("UVB_SELF_SHIELDING: using %d self shielding table entries from file `%s'.\n", Nselfshieldingtab, fname);
}
#endif

#ifdef GFM_AGN_RADIATION
/** \brief Set the ionization parameters for the AGN background.
 */
void IonizeParamsAGN(void)
{
  double fac = gs.AGNBolIntensity / PhotoTAGN[0].variable;
  /* make sure photo heating/ionization is turned on */
  pc.J_UV = 1;
  /* add scaled AGN photo rates */
  pc.gJH0 += fac * PhotoTAGN[0].gH0;
  pc.gJHe0 += fac * PhotoTAGN[0].gHe;
  pc.gJHep += fac * PhotoTAGN[0].gHep;
  pc.epsH0 += fac * PhotoTAGN[0].eH0;
  pc.epsHe0 += fac * PhotoTAGN[0].eHe;
  pc.epsHep += fac * PhotoTAGN[0].eHep;
}
#endif

#ifdef RADCOOL
/** \brief Set the ionization parameters for the Local Sources.
 */
void IonizeParamsRADCOOL(void)
{
  // double fac = P / PhotoTAGN[0].variable;
  /* make sure photo heating/ionization is turned on */
  pc.J_UV = 1;

#ifdef RADCOOL_HOTHALO
  double fac = 1.0 / PhotoTRAD_T6[0].variable;
#endif

  /* add scaled stellar  photo rates */
  pc.gJH0 = pc.gJH0 + gs.Phios * PhotoTRAD_OS[0].gH0 + gs.Phins * PhotoTRAD_NS[0].gH0
#ifdef RADCOOL_HOTHALO
            + ((gs.PhiT6 * PhotoTRAD_T6[0].gH0 + gs.PhiT7 * PhotoTRAD_T7[0].gH0 + gs.PhiT8 * PhotoTRAD_T8[0].gH0) * fac)
#endif
      ;
  pc.gJHe0 = pc.gJHe0 + gs.Phios * PhotoTRAD_OS[0].gHe + gs.Phins * PhotoTRAD_NS[0].gHe
#ifdef RADCOOL_HOTHALO
             + ((gs.PhiT6 * PhotoTRAD_T6[0].gHe + gs.PhiT7 * PhotoTRAD_T7[0].gHe + gs.PhiT8 * PhotoTRAD_T8[0].gHe) * fac)
#endif
      ;
  pc.gJHep = pc.gJHep + gs.Phios * PhotoTRAD_OS[0].gHep + gs.Phins * PhotoTRAD_NS[0].gHep
#ifdef RADCOOL_HOTHALO
             + ((gs.PhiT6 * PhotoTRAD_T6[0].gHep + gs.PhiT7 * PhotoTRAD_T7[0].gHep + gs.PhiT8 * PhotoTRAD_T8[0].gHep) * fac)
#endif
      ;
  pc.epsH0 = pc.epsH0 + gs.Phios * PhotoTRAD_OS[0].eH0 + gs.Phins * PhotoTRAD_NS[0].eH0
#ifdef RADCOOL_HOTHALO
             + ((gs.PhiT6 * PhotoTRAD_T6[0].eH0 + gs.PhiT7 * PhotoTRAD_T7[0].eH0 + gs.PhiT8 * PhotoTRAD_T8[0].eH0) * fac)
#endif
      ;
  pc.epsHe0 = pc.epsHe0 + gs.Phios * PhotoTRAD_OS[0].eHe + gs.Phins * PhotoTRAD_NS[0].eHe
#ifdef RADCOOL_HOTHALO
              + ((gs.PhiT6 * PhotoTRAD_T6[0].eHe + gs.PhiT7 * PhotoTRAD_T7[0].eHe + gs.PhiT8 * PhotoTRAD_T8[0].eHe) * fac)
#endif
      ;
  pc.epsHep = pc.epsHep + gs.Phios * PhotoTRAD_OS[0].eHep + gs.Phins * PhotoTRAD_NS[0].eHep
#ifdef RADCOOL_HOTHALO
              + ((gs.PhiT6 * PhotoTRAD_T6[0].eHep + gs.PhiT7 * PhotoTRAD_T7[0].eHep + gs.PhiT8 * PhotoTRAD_T8[0].eHep) * fac)
#endif
      ;
}
#endif

/*! \brief Set the ionization parameters for the UV background.
 *
 *  \return void
 */
void IonizeParamsUVB(void)
{
#ifdef UVB_OFF
  SetZeroIonization();
  return;
#else
  int i, ilow;
  double logz, dzlow, dzhi;
  double redshift;

  if(All.ComovingIntegrationOn)
    {
      redshift = 1 / All.Time - 1;
#ifdef UVB_START
      if(redshift > All.UVBStartRedshift)
        {
          SetZeroIonization();
          return;
        }
#endif
    }
  else
    {
      redshift = 0.0;
      // SetZeroIonization();
      // return;
    }

  logz = log10(redshift + 1.0);
  ilow = 0;
  for(i = 0; i < NheattabUVB; i++)
    {
      if(PhotoTUVB[i].variable < logz)
        ilow = i;
      else
        break;
    }

  dzlow = logz - PhotoTUVB[ilow].variable;
  dzhi  = PhotoTUVB[ilow + 1].variable - logz;

  if(NheattabUVB == 0 || logz > PhotoTUVB[NheattabUVB - 1].variable || PhotoTUVB[ilow].gH0 == 0 || PhotoTUVB[ilow + 1].gH0 == 0)
    {
      SetZeroIonization();
      return;
    }
  else
    pc.J_UV = 1;

  pc.gJH0   = pow(10., (dzhi * log10(PhotoTUVB[ilow].gH0) + dzlow * log10(PhotoTUVB[ilow + 1].gH0)) / (dzlow + dzhi));
  pc.gJHe0  = pow(10., (dzhi * log10(PhotoTUVB[ilow].gHe) + dzlow * log10(PhotoTUVB[ilow + 1].gHe)) / (dzlow + dzhi));
  pc.gJHep  = pow(10., (dzhi * log10(PhotoTUVB[ilow].gHep) + dzlow * log10(PhotoTUVB[ilow + 1].gHep)) / (dzlow + dzhi));
  pc.epsH0  = pow(10., (dzhi * log10(PhotoTUVB[ilow].eH0) + dzlow * log10(PhotoTUVB[ilow + 1].eH0)) / (dzlow + dzhi));
  pc.epsHe0 = pow(10., (dzhi * log10(PhotoTUVB[ilow].eHe) + dzlow * log10(PhotoTUVB[ilow + 1].eHe)) / (dzlow + dzhi));
  pc.epsHep = pow(10., (dzhi * log10(PhotoTUVB[ilow].eHep) + dzlow * log10(PhotoTUVB[ilow + 1].eHep)) / (dzlow + dzhi));

  return;
#endif
}

#ifdef UVB_SELF_SHIELDING
/** \brief Calculate the self-shielding factor of the radiation.
 *
 *  \param number density of hydrogen atoms in cm^-3
 */
double calc_self_shielding_factor(double nH)
{
  int i, ilow;
  double dzlow, dzhi;
  double redshift;

  float alpha1, alpha2, beta, xi, n0, f;

  if(All.ComovingIntegrationOn)
    redshift = 1 / All.Time - 1;
  else
    {
#if defined(SMUGGLE_MOLEC_COOLING) || defined(WINDTUNNEL_FIXVARIABLESININJECTIONREGION)
      redshift = 0.0;
#else
      return 1.0;
#endif
    }

  ilow = 0;
  for(i = 0; i < Nselfshieldingtab; i++)
    {
      if(SelfShieldingParams[i].redshift < redshift)
        ilow = i;
      else
        break;
    }

  dzlow = redshift - SelfShieldingParams[ilow].redshift;
  dzhi  = SelfShieldingParams[ilow + 1].redshift - redshift;

  if(redshift > SelfShieldingParams[Nselfshieldingtab - 1].redshift || Nselfshieldingtab == 0)
    return 1.0;

  alpha1 = (dzhi * SelfShieldingParams[ilow].alpha1 + dzlow * SelfShieldingParams[ilow + 1].alpha1) / (dzlow + dzhi);
  alpha2 = (dzhi * SelfShieldingParams[ilow].alpha2 + dzlow * SelfShieldingParams[ilow + 1].alpha2) / (dzlow + dzhi);
  beta   = (dzhi * SelfShieldingParams[ilow].beta + dzlow * SelfShieldingParams[ilow + 1].beta) / (dzlow + dzhi);
  xi     = (dzhi * SelfShieldingParams[ilow].xi + dzlow * SelfShieldingParams[ilow + 1].xi) / (dzlow + dzhi);
  n0     = (dzhi * SelfShieldingParams[ilow].n0 + dzlow * SelfShieldingParams[ilow + 1].n0) / (dzlow + dzhi);
  f      = (dzhi * SelfShieldingParams[ilow].f + dzlow * SelfShieldingParams[ilow + 1].f) / (dzlow + dzhi);

  double self_shielding_factor = (1 - f) * pow(1 + pow(nH / pow(10, n0), beta), alpha1) + f * pow(1 + nH / pow(10, n0), alpha2);

  if(self_shielding_factor > 1)
    {
      warn("UVB_SELF_SHIELDING: self_shielding_factor=%g > 1. alpha1=%g alpha2=%g beta=%g xi=%g n0=%g f=%g, and nH=%g",
           self_shielding_factor, alpha1, alpha2, beta, xi, n0, f, nH);
      self_shielding_factor = 1;
    }

  self_shielding_factor = 1 - xi * (1 - self_shielding_factor);

  return self_shielding_factor;
}
#endif

/*! \brief Reset the ionization parameters.
 *
 *  \return void
 */
void SetZeroIonization(void) { memset(&pc, 0, sizeof(PhotoCurrent)); }

/*! \brief Wrapper function to set the ionizing background.
 *
 *  \return void
 */
void IonizeParams(void) { IonizeParamsUVB(); }

/*! \brief Initialize the cooling module.
 *
 *  This function initializes the cooling module. In particular,
 *  it allocates the memory for the cooling rate and ionization tables
 *  and initializes them.
 *
 *  \return void
 */
void InitCool(void)
{
  /* set default hydrogen mass fraction */
  gs.XH = HYDROGEN_MASSFRAC;

#if defined(GFM_COOLING_METAL) && !defined(MRT_METAL_COOLING)
  /* set metallicity floor */
  MetallicityFloor           = pow(10.0, GFM_MIN_METAL);
  LogMetallicityFloorInSolar = log10(MetallicityFloor / GFM_SOLAR_METALLICITY);

  /* convert lower metal line cooling limit to log */
  if(All.MinMetalTemp > 0)
    LogMinMetalTemp = log10(All.MinMetalTemp);
  else
    LogMinMetalTemp = -MAX_REAL_NUMBER;

  mpi_printf("GFM_COOLING_METAL: All.MinMetalTemp=%g LogMinMetalTemp=%g MetallicityFloor=%g LogMetallicityFloorInSolar=%g\n",
             All.MinMetalTemp, LogMinMetalTemp, MetallicityFloor, LogMetallicityFloorInSolar);

#ifdef SFR_MCS
  metal_cool_mcs_init(MetallicityFloor, LogMetallicityFloorInSolar);
#endif
#endif

  /* zero photo-ionization/heating rates */
  SetZeroIonization();

  /* allocate and construct rate table */
  RateT = (RateTable *)mymalloc("RateT", (NCOOLTAB + 1) * sizeof(RateTable));
  ;
  MakeRateTable();

  /* read photo tables */
  ReadIonizeParams(All.TreecoolFile, 0);

#ifdef GFM_AGN_RADIATION
  ReadIonizeParams(All.TreecoolFileAGN, 1);
#endif

#ifdef RADCOOL
  ReadIonizeParams(All.TreecoolFileRAD, 2);
#endif

#ifdef UVB_SELF_SHIELDING
  ReadSelfShieldingParams(All.SelfShieldingFile);
#endif

  mpi_printf("GFM_COOLING: time, time begin = %le\t%le\n", All.Time, All.TimeBegin);
  All.Time = All.TimeBegin;
  set_cosmo_factors_for_current_time();

  IonizeParams();

#if defined(GFM_UVB_CORRECTIONS) || defined(SMUGGLE_COSMIC_RAY_HEATING)
  gs.RhoGasMean = All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
#endif

#if(defined(GFM_AGN_RADIATION) || defined(GFM_UVB_CORRECTIONS) || defined(UVB_SELF_SHIELDING) || defined(RADCOOL)) && !defined(SFR_MCS)
  reset_radiation_state();
#endif
}

/** \brief Apply the isochoric cooling to all the active gas cells.
 *
 */
void cooling_only(void) /* normal cooling routine when star formation is disabled */
{
  int idx, i;

#ifndef MRT
  CPU_Step[CPU_MISC] += measure_time();
#endif

#ifdef MRT
#ifdef MRT_EQUIL_CHEM_COOL
  mpi_printf_rt(0, "RT: MRT_EQUIL_CHEM_COOL doing chemistry and cooling...\n");
#else
  mpi_printf_rt(0, "RT: Caclulating the isochoric cooling rates...\n");
#endif
#endif

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i >= 0)
        {
          if(P[i].Mass == 0 && P[i].ID == 0)
            continue; /* skip cells that have been swallowed or eliminated */

#ifdef MRT
          // #ifdef MRT_EQUIL_CHEM_COOL
          if(do_cooling_mrt(i))
// #endif
#endif
            cool_cell(i);
        }
    }
#ifndef MRT
  CPU_Step[CPU_COOLINGSFR] += measure_time();
#else
  CPU_Step[CPU_RT_COOLING] += measure_time();
#endif
}

/*! \brief Apply the isochoric cooling to a given gas cell.
 *
 *  This function applies the normal isochoric cooling to a single gas cell.
 *  Once the cooling has been applied according to one of the cooling models
 *  implemented, the internal energy per unit mass, the total energy and the
 *  pressure of the cell are updated.
 *
 *  \param[in] i Index of the gas cell to which cooling is applied.
 *
 *  \return void
 */
void cool_cell(int i)
{
  double dt, dtime, ne = 1;
  double unew, dens, dtcool;

  dens = SphP[i].Density;

  dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

  dtime = All.cf_atime * dt / All.cf_time_hubble_a;

#ifdef MRT

#ifdef MRT_SUBCYCLE
  dtime /= ((double)(All.RTNumSubCycles));
#else
  dtime *= 0.5;
#endif

#endif

  dtcool = dtime;

#if defined(GFM_AGN_RADIATION)
  update_radiation_state(dens * All.cf_a3inv, get_hydrogen_abundances_of_local_cell(i), SphP[i].AGNBolIntensity);
#elif defined(RADCOOL) && defined(RADCOOL_HOTHALO)
  update_radiation_state(dens * All.cf_a3inv, get_hydrogen_abundances_of_local_cell(i), SphP[i].Phios, SphP[i].Phins, SphP[i].PhiT6,
                         SphP[i].PhiT7, SphP[i].PhiT8);
#elif defined(RADCOOL)
  update_radiation_state(dens * All.cf_a3inv, get_hydrogen_abundances_of_local_cell(i), SphP[i].Phios, SphP[i].Phins);
#elif(defined(GFM_UVB_CORRECTIONS) || defined(UVB_SELF_SHIELDING)) && defined(GFM_STELLAR_EVOLUTION)
  update_radiation_state(dens * All.cf_a3inv, get_hydrogen_abundances_of_local_cell(i), 0);
#elif defined(UVB_SELF_SHIELDING) && defined(SFR_MCS) && defined(GFM_COOLING_METAL)
  update_radiation_state(dens * All.cf_a3inv, SphP[i].Metallicity, &pc);
#elif defined(GFM_UVB_CORRECTIONS) || defined(UVB_SELF_SHIELDING)
  update_radiation_state(dens * All.cf_a3inv, HYDROGEN_MASSFRAC, 0);
#endif

#ifdef MRT
#ifdef MRT_EQUIL_CHEM_COOL
  update_radiation_state_mrt(i, &pc);
#endif
#endif

#if defined(GFM_COOLING_METAL) && !defined(MRT_METAL_COOLING)
#ifdef SFR_MCS
  update_gas_state(dens * All.cf_a3inv, SphP[i].Metallicity, &gs);
#else
  update_gas_state(dens * All.cf_a3inv, get_hydrogen_abundances_of_local_cell(i), SphP[i].Metallicity);
#endif
#endif

#ifdef GFM_LAMBDA
  gs.partindex = i;
#endif

#ifdef GFM_DUST_COOLING
  get_CoolingDustState(i, &gs.dust_to_gas_ratio, &gs.dust_cool);
#endif

#if defined(SMUGGLE_RADIATION_FEEDBACK)
  if(SphP[i].GasRadCoolShutoffTime > 0)
    dtcool = update_radiation_cooling_shutoff_time(i, dtime);
#endif

#ifndef SIMPLE_COOLING
#if(!defined(OTVET_COOLING_HEATING) && !defined(MRT_COOLING_HEATING)) || defined(MRT_EQUIL_CHEM_COOL)
#ifndef RT_COOLING_PHOTOHEATING /* if RT_COOLING is set, do not do the cooling here */
  ne = SphP[i].Ne;              /* electron abundance (gives ionization state and mean molecular weight) */
  if(ne < 0)
    terminate("SphP.Ne has gone below 0!!!\n");

  unew = DoCooling(fmax(All.MinEgySpec, SphP[i].Utherm), dens * All.cf_a3inv, dtcool, &ne);
  /* Xiaohan: keep the two lines below for the Stromgren sphere test */
  // convert_u_to_temp(fmax(All.MinEgySpec, SphP[i].Utherm) * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs, dens * All.cf_a3inv *
  // All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam, &ne); unew = fmax(All.MinEgySpec, SphP[i].Utherm);

  SphP[i].Ne = ne;

  if(unew < 0)
    terminate("invalid temperature: Thistask=%d i=%d unew=%g\n", ThisTask, i, unew);

#ifdef MRT_EQUIL_CHEM_COOL /* update abundances and photon density */
  update_chem_mrt(i, dtcool, &gs);
#endif

  double du = unew - SphP[i].Utherm;
#else
  SphP[i].Utherm += rt_DoHeating(i, dtcool);
  double du = rt_DoCooling(i, dtcool);
#endif
#else /*OTVET or MRT */

#ifdef MRT_COOLING_HEATING

#ifndef MRT_NO_UV
  double du_heat = mrt_DoHeating(i, dtcool);
#else
  double du_heat = 0.0;
#endif

#else
  double du_heat = otvet_DoHeating(i, dtcool);
#endif

#ifdef MRT_COOLING_HEATING

#ifndef MRT_NO_UV
  double du      = mrt_DoCooling(i, dtcool);
#else
  double du      = 0.0;
#endif

#else
  double du = otvet_DoCooling(i, dtcool);
#endif

#if defined(MRT_IR_LTE) && !defined(MRT_COUPLED_THERMOCHEMISTRY)
  double du_IR   = mrt_update_IR_cooling(i, dtcool);
  du += du_IR;
#endif

  if(isnan(du))
    {
      print_particle_info(i);
#ifdef OTVET
      printf("OTVET error DoCooling %g %g \n", du, SphP[i].n_gamma[0]);
#else
      printf("MRT error DoCooling %g %g \n", du, SphP[i].DensPhot[0]);
#endif

      terminate("991");
    }
  du   = du_heat + du;
  unew = SphP[i].Utherm + du;
  // if(P[i].Pos[0]<All.BoxSize/30.0)
  // printf("%g\t%g\t%g\t%g\n", du_heat, du, unew, SphP[i].Utherm) ;

  if(unew < 0.0)
    terminate("OTVET/MRT invalid temperature: Thistask=%d i=%d unew=%g %g %g\n", ThisTask, i, unew, du_heat, du);
#endif
#else /* SIMPLE COOLING */
  unew = SphP[i].Utherm + DoSimpleCoolingHeating(i, dtcool);
  if(unew < 0)
    terminate("invalid temperature: Thistask=%d i=%d unew=%g\n", ThisTask, i, unew);

  double du = unew - SphP[i].Utherm;
#endif

#if defined(SMUGGLE_RADIATION_FEEDBACK)
  /* LVS: keep T = 10^4 */
  if(dtcool <= 0.0 && All.Time > All.TimeBegin)
    {
      if((SphP[i].Utherm + du) < All.PhotoionizationEgySpec)
        {
          du   = All.PhotoionizationEgySpec - SphP[i].Utherm;
          unew = SphP[i].Utherm + du;
          ne   = SphP[i].Ne;
          convert_u_to_temp(unew, dens * All.cf_a3inv, &ne);
          SphP[i].Ne = ne;
        }
    }
#endif

  if(unew < All.MinEgySpec)
    du = All.MinEgySpec - SphP[i].Utherm;

#ifdef MRT
#ifndef MRT_UPDATE_AT_END_OF_STEP
  SphP[i].Utherm += du;
  SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;
#else
  SphP[i].RT_dutherm += du;
  SphP[i].RT_utherm += du;
#endif
#else
  SphP[i].Utherm += du;
  SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;
#endif

#ifdef OUTPUT_COOLHEAT
  if(dtime > 0)
    SphP[i].CoolHeat = du * P[i].Mass / dtime;
#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  SphP[i].A       = (GAMMA - 1.0) * SphP[i].Utherm / pow(dens * All.cf_a3inv, GAMMA - 1);
  SphP[i].Entropy = log(SphP[i].A) * P[i].Mass;
#endif

#ifdef MRT
#ifndef MRT_UPDATE_AT_END_OF_STEP
  set_pressure_of_cell(i);
#endif
#else
  set_pressure_of_cell(i);
#endif
}

#ifdef PREHEATING
void impose_preheating(void)
{
  if((All.HighestActiveTimeBin == All.HighestOccupiedTimeBin) && (All.Time > All.TimePreheating) && (!All.FlagPreheating))
    {
      double ulimit_save = All.MinEgySpec;
      double heliumfraction;
      double meanweight;

      heliumfraction = (1.0 - HYDROGEN_MASSFRAC) / (4.0 * HYDROGEN_MASSFRAC);

      meanweight = (1.0 + 4.0) * heliumfraction / (1 + 1 + 3 * heliumfraction); /* note: assuming fully ionized GAS */

      All.MinEgySpec = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempPreheating;
      All.MinEgySpec *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

      update_primitive_variables(); /* these effectively closes off the hydro step */

      All.MinEgySpec     = ulimit_save;
      All.FlagPreheating = 1;
    }
}
#endif

#ifdef MRT_INIT_IONIZATION
void SetInitGasState(int i)
{
  double sfr = 0;
  double rho = SphP[i].Density * All.cf_a3inv;
  double u   = fmax(All.MinEgySpec, SphP[i].Utherm);
  double ne_guess;
  ne_guess = 1.0;

  gs.XH      = HYDROGEN_MASSFRAC;
  gs.yhelium = (1 - gs.XH - 0.0) / (4. * gs.XH);

  /* update GasState as appropriate given compile-time options and cell properties */
#if defined(USE_SFR) && (!defined(LOCAL_FEEDBACK) && !defined(SMUGGLE_STAR_FEEDBACK))
  sfr = get_starformation_rate(i);
#endif

  if(sfr == 0)
    {
#if defined(GFM_AGN_RADIATION) || defined(GFM_UVB_CORRECTIONS) || defined(UVB_SELF_SHIELDING)
#ifdef GFM_AGN_RADIATION
      update_radiation_state(rho, get_hydrogen_abundances_of_local_cell(i), SphP[i].AGNBolIntensity);
#else

#ifdef GFM_STELLAR_EVOLUTION
      update_radiation_state(rho, get_hydrogen_abundances_of_local_cell(i), 0);
#else
      update_radiation_state(rho, HYDROGEN_MASSFRAC, 0);
#endif

#endif
#endif

#if defined(GFM_COOLING_METAL) && !defined(MRT_METAL_COOLING)
      update_gas_state(rho, get_hydrogen_abundances_of_local_cell(i), SphP[i].Metallicity);
#endif
    }
  else
    {
#if defined(GFM_AGN_RADIATION) || defined(GFM_UVB_CORRECTIONS) || defined(UVB_SELF_SHIELDING)
      update_radiation_state(rho, HYDROGEN_MASSFRAC, 0);
#endif
#if defined(GFM_COOLING_METAL) && !defined(MRT_METAL_COOLING)
      update_gas_state(rho, GFM_INITIAL_ABUNDANCE_HYDROGEN, 0.0);
#endif
    }

  /* update DoCool */
  DoCool.u_old_input    = u;
  DoCool.rho_input      = rho;
  DoCool.ne_guess_input = ne_guess;

  /* convert to physical cgs units */
  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
  u *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;

  /* calculate cooling rate (and so ne_guess and all of gs including nH0, nHeII) */
  CoolingRateFromU(u, rho, &ne_guess);

  double nH_times_volume = P[i].Mass;
  double y_fac           = (1.0 - HYDROGEN_MASSFRAC) / 4.0 / HYDROGEN_MASSFRAC;
  double nh_fac          = 1.0 / (1.0 + y_fac);
  SphP[i].HI             = gs.nH0;
  SphP[i].HII            = gs.nHp;
  SphP[i].nHI            = SphP[i].HI * nH_times_volume;
  SphP[i].nHII           = SphP[i].HII * nH_times_volume;
#ifdef MRT_INCLUDE_HE
  SphP[i].HeI    = gs.nHe0;
  SphP[i].HeII   = gs.nHep;
  SphP[i].HeIII  = gs.nHepp;
  SphP[i].Ne     = SphP[i].HII + SphP[i].HeII + 2.0 * SphP[i].HeIII;
  SphP[i].nHeI   = SphP[i].HeI * nH_times_volume;
  SphP[i].nHeII  = SphP[i].HeII * nH_times_volume;
  SphP[i].nHeIII = SphP[i].HeIII * nH_times_volume;
#else
  SphP[i].Ne = SphP[i].HII;
#endif
  SphP[i].ne = SphP[i].Ne * nH_times_volume;
}
#endif

#endif
