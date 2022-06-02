#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../allvars.h"
#include "../proto.h"
#include "f2c.h"
#include "sgchem_proto.h"

/* See README file for a description, authors list and usage policy */

void init_chemistry(void)
{
  /* Initialize parameters stored in coolr and cooli common blocks */
  INIT_CHEMISTRY_PARAMETERS(&All.DeutAbund,
#ifndef SGCHEM_VARIABLE_Z
                            &All.CarbAbund, &All.OxyAbund, &All.MAbund, &All.ZAtom,
#endif
                            &All.InitDustTemp,
#ifndef SGCHEM_VARIABLE_ISRF
                            &All.UVFieldStrength,
#endif
                            &All.LWBGType, &All.LWBGStartRedsh,
#ifndef SGCHEM_VARIABLE_Z
                            &All.DustToGasRatio,
#endif
#ifndef SGCHEM_VARIABLE_CRION
                            &All.CosmicRayIonRate,
#endif
#ifdef SGCHEM_TEMPERATURE_FLOOR
                            &All.SGChemTemperatureFloor,
#endif
                            &All.InitRedshift, &All.ExternalDustExtinction, &All.H2FormEx, &All.H2FormKin, &All.PhotoApprox,
                            &All.ISRFOption, &All.AtomicCoolOption, &All.H2OpacityOption);

  /* Initialize chemical rates */
  COOLINMO();
  CHEMINMO();

  /* ODE integrator tolerances */
  INIT_TOLERANCES();

  /* MCMA setup, if using NL99 network */
#ifdef MCMA
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  CMA_INIT();
#endif
#endif

  /* Gamma and temperature lookup, if using variable gamma */
#ifdef VARIABLE_GAMMA
  INIT_TEMPERATURE_LOOKUP();
#endif

  /* H2 cooling rate opacity correction, if running with primordial chemistry */
#if CHEMISTRYNETWORK == 1
  LOAD_H2_TABLE();
#endif

  return;
}

void evolve_chemistry(void)
{
  double dt;
  int i;

  mpi_printf("SGCHEM: Looping over the chemistry for active cells \n");
  /* Loop over active particles, and do chemical evolution for each one */
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

#ifdef SINK_PHOTOION_FEEDBACK
      if(SphP[i].CoolingFlag == 0)
        {
          SphP[i].CoolingFlag = 1;
          continue;
        }
#endif

      dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

#ifdef MRT_CHEM_SG
#ifdef MRT_SUBCYCLE
      dt /= ((double)(All.RTNumSubCycles));
#else
      dt *= 0.5;
#endif
#endif

      evolve_chemistry_for_single_cell(dt, i);
    }

  return;
}

/* Chemical evolution, radiative heating & cooling */
void evolve_chemistry_for_single_cell(double dt, int index)
{
  double time, rho, energy, energy_cgs, divv, dl, yn;
  double rho_cgs, redshift, dust_temp;
  double non_eq_abundances[SGCHEM_NUM_SPECIES];
  double column_density_projection[NPIX];
  double column_density_projection_H2[NPIX];
  double column_density_projection_CO[NPIX];
  double column_density_projection_C[NPIX];
  double photorates_direct[NDIRECT];
  double a, a3inv, hubble_a, hubble_param, hubble_param2;
  double column_correction_cosmic, column_correction_cgs;
  double column_correction_NH, column_correction_NH2, column_correction_NCO, column_correction_NC;
  double Utherm_diff, Utherm_new, gamma, gamma_minus1;
  double abh2, abe, ekn, ntot;
  double dust_to_gas, carb_abund, oxy_abund, m_abund, Z_atom;
  double breakdown_of_thermal_rates[SGCHEM_NUM_THERMAL_RATES];
  double d2sink, msink, mdotsink, rMS, p1, p2, rsink, sink_flux;
  double sink_accretion_luminosity;
  double solar_radius = 6.957e10; /* in cgs units */
  double fsuppress;
#ifdef SGCHEM_SUPPRESS_COOLING_ABOVE_RHOMAX
  double rho_max, density_ratio;
#endif
  int i, j, id;

#if NO_CHEMISTRY
  return;
#else

  /* Check size of timestep */
  if(dt < 0)
    {
      terminate("Error: negative timestep in evolve_chemistry");
    }
  else if(dt == 0)
    {
      return;
    }

  /* Set current index and ID in cool.h */
  id = P[index].ID;
  SET_INDEX_ID_FOR_CHEMISTRY(&index, &id);

  /* Starting values */
  time = dt;
  rho  = SphP[index].Density;
  divv = SphP[index].DivVel;
#ifdef STATIC_CHEMISTRY_TEST
  /* Allow direct comparison with one-zone model */
  divv = 0;
#endif

  /* Estimate of typical distance from centre to edge of local grid-cell, used for computing
   * local contribution to shielding.
   *
   * XXX: should we make it possible to switch this on or off?
   */
  dl = 0.5 * pow(SphP[index].Volume, 1. / 3.);

  for(i = 0; i < SGCHEM_NUM_SPECIES; i++)
    {
      non_eq_abundances[i] = SphP[index].TracAbund[i];
    }

    /* Get local element abundances and dust-to-gas ratio, if they're spatially variable */
#ifdef SGCHEM_VARIABLE_Z
  carb_abund = SphP[index].CarbAbund;
  oxy_abund  = SphP[index].OxyAbund;
  m_abund    = SphP[index].MAbund;
  Z_atom     = SphP[index].ZAtom;
  SET_LOCAL_ELEMENT_ABUNDANCES(&carb_abund, &oxy_abund, &m_abund, &Z_atom);

  dust_to_gas = SphP[index].DustToGasRatio;
  SET_LOCAL_DUST_ABUNDANCE(&dust_to_gas);
#else
  /* Setting this here avoids lots of messy conditional code later. We could also set the other
     values (oxy_abund, m_abund, etc.) in the same way, but right now they're not used later on */
  carb_abund = All.CarbAbund;
  oxy_abund  = All.OxyAbund;
  m_abund    = All.MAbund;
#endif

#ifdef SGCHEM_VARIABLE_ISRF
  find_and_set_local_ISRF(index);
#endif
#ifdef SGCHEM_VARIABLE_CRION
  find_and_set_local_CRION(index);
#endif

  /* If TREECOLV2 and friends are not in use, then these are all single-element arrays */
  column_density_projection[0] = column_density_projection_H2[0] = column_density_projection_CO[0] = column_density_projection_C[0] =
      0.0;

#ifdef TREECOLV2
  for(i = 0; i < NPIX; i++)
    {
      column_density_projection[i] = SphP[index].Projection[i];

      /*    if (isnan(column_density_projection[i])){
         printf("SGChem: NaN coldensity i %d index %d id %d treecol projection%g\n",i,index,id,SphP[index].Projection[i]);
         } */
#ifdef TREECOLV2_H2
      column_density_projection_H2[i] = SphP[index].ProjectionH2[i];
      if(isnan(column_density_projection[i]))
        {
          printf("SGChem: NaN H2 coldensity i %d index %d id %d treecol projection%g\n", i, index, id, SphP[index].Projection[i]);
        }
#else
      column_density_projection_H2[i] = 0.0;
#endif
#ifdef TREECOLV2_CO
      column_density_projection_CO[i] = SphP[index].ProjectionCO[i];
#else
      column_density_projection_CO[i] = 0.0;
#endif
#ifdef TREECOLV2_C
      column_density_projection_C[i]  = SphP[index].ProjectionC[i];
#else
      column_density_projection_C[i]  = 0.0;
#endif
    }
#endif

  /* Convert from comoving to physical units */
  if(All.ComovingIntegrationOn)
    {
      a             = All.Time;
      a3inv         = 1 / (a * a * a);
      hubble_a      = hubble_function(a);
      hubble_param  = All.HubbleParam;
      hubble_param2 = hubble_param * hubble_param;
      redshift      = (1 / a) - 1;
    }
  else
    {
      a = a3inv     = 1;
      hubble_a      = All.cf_hubble_a;
      hubble_param  = All.HubbleParam;
      hubble_param2 = hubble_param * hubble_param;
      redshift      = All.InitRedshift;
    }

  time *= 1 / hubble_a / hubble_param;
  rho *= a3inv * hubble_param2;

#ifdef MRT_CHEM_SG
  double vel2 = SphP[index].Momentum[0] * SphP[index].Momentum[0] + SphP[index].Momentum[1] * SphP[index].Momentum[1] +
                SphP[index].Momentum[2] * SphP[index].Momentum[2];
  double utherm_particle = (SphP[index].Energy / (All.cf_atime * All.cf_atime) - 0.5 * vel2 / P[index].Mass) / P[index].Mass;
  energy                 = utherm_particle * rho; /* Energy here is energy density, AREPO evolves specific energy */
#else
  energy     = SphP[index].Utherm * rho; /* Energy here is energy density, AREPO evolves specific energy */
#endif

  // TODO: This should print out more info if code crashes on "Calc H2 shielding", but it is not a fix!!
  if(energy <= 0)
    {
      printf("sgchem.c: particle with ID %d has zero energy\n", P[index].ID);
      printf("          Energy %.03e Utherm %.03e Volume %.03e \n", SphP[index].Energy, SphP[index].Utherm, SphP[index].Volume);
      printf("          Density %.03e Pressure %.03e \n", SphP[index].Density, SphP[index].Pressure);
      printf("          Momentum [ %.03e %.03e %.03e ] \n", SphP[index].Momentum[0], SphP[index].Momentum[1], SphP[index].Momentum[2]);
      printf("          Pos          [%.03e %.03e %.03e]\n", P[index].Pos[0], P[index].Pos[1], P[index].Pos[2]);
      printf("          TracAbund    [ ");
      for(i = 0; i < SGCHEM_NUM_SPECIES; i++)
        printf("%.03e ", SphP[index].TracAbund[i]);
      printf("]\n");
#ifdef SIMPLEX
      printf("          PhotonRates  [ ");
      for(i = 0; i < NDIRECT; i++)
        printf("%.03e ", SphP[index].sxPhotonRates[i]);
      printf("]\n");
#endif
    }

  divv *= 1.0 / a;
  if(All.ComovingIntegrationOn)
    divv += 3 * hubble_a;
  dl *= a;

  column_correction_cosmic = hubble_param / (a * a);
  for(i = 0; i < NPIX; i++)
    {
      column_density_projection[i] *= column_correction_cosmic;
      column_density_projection_H2[i] *= column_correction_cosmic;
      column_density_projection_CO[i] *= column_correction_cosmic;
      column_density_projection_C[i] *= column_correction_cosmic;
    }

  /* Convert from AREPO code units to cgs */
  time *= All.UnitTime_in_s;
  rho_cgs    = rho * All.UnitDensity_in_cgs; /* We keep rho in code units as we need it later */
  energy_cgs = energy * All.UnitEnergy_in_cgs / pow(All.UnitLength_in_cm, 3);
  divv *= 1.0 / All.UnitTime_in_s;
  dl *= All.UnitLength_in_cm;

  /* Calculate cooling suppression factor, if requested. This allows us to artificially suppress
   * cooling above some maximum density, which is set in the parameter file. */
  fsuppress     = 1.0;
#ifdef SGCHEM_SUPPRESS_COOLING_ABOVE_RHOMAX
  rho_max       = All.SGChemMaxDensityForCooling * a3inv * hubble_param2 * All.UnitDensity_in_cgs;
  density_ratio = rho_cgs / rho_max;
  if(density_ratio < 1)
    {
      fsuppress = 1.0;
    }
  else if(density_ratio > 10)
    {
      fsuppress = 0.0;
    }
  else
    {
      fsuppress = exp(1.0 - (density_ratio * density_ratio));
    }
#endif

  column_correction_cgs = All.UnitDensity_in_cgs * All.UnitLength_in_cm;
  for(i = 0; i < NPIX; i++)
    {
      column_density_projection[i] *= column_correction_cgs;
      column_density_projection_H2[i] *= column_correction_cgs;
      column_density_projection_CO[i] *= column_correction_cgs;
      column_density_projection_C[i] *= column_correction_cgs;
    }

  /* Convert from mass surface density to particle column
   * density
   */
  column_correction_NH  = 1.0 / ((1.0 + 4.0 * ABHE) * PROTONMASS);
  column_correction_NH2 = 1.0 / (2.0 * PROTONMASS);
  column_correction_NCO = 1.0 / (28.0 * PROTONMASS);
  column_correction_NC  = 1.0 / (12.0 * PROTONMASS);
  for(i = 0; i < NPIX; i++)
    {
      column_density_projection[i] *= column_correction_NH;
      column_density_projection_H2[i] *= column_correction_NH2;
      column_density_projection_CO[i] *= column_correction_NCO;
      column_density_projection_C[i] *= column_correction_NC;
    }

  /* Compute derived units */
  yn = rho_cgs / ((1.0 + 4.0 * ABHE) * PROTONMASS);

  sink_flux = 0.0;
#if defined(SGCHEM_ACCRETION_LUMINOSITY) && defined(SINK_PARTICLES)
  if(All.SGChemAccretionLuminosityOn && NSinksAllTasks > 0)
    {
      for(i = 0; i < NSinksAllTasks; i++)
        {
          if(SinkP[i].ID == 0 && SinkP[i].Mass == 0)
            continue;

          d2sink = 0.0;
          for(j = 0; j < 3; j++)
            {
              d2sink += pow((SinkP[i].Pos[j] - P[index].Pos[j]), 2);
            }
          msink    = SinkP[i].Mass;
          mdotsink = SinkP[i].AccretionRate;
          /* Now convert to cgs units */

          msink *= All.UnitMass_in_g / SOLAR_MASS;
          mdotsink *= (All.UnitMass_in_g / All.UnitTime_in_Megayears) * (1e-6 / SOLAR_MASS); /* Accretion rate in Msun/yr */
          if(mdotsink >= 0.01)
            {
              mdotsink = 0.01; /* Limit maximum value used for stability */
            }

          if(All.ComovingIntegrationOn)
            {
              d2sink *= a * a;
            }
          d2sink *= pow(All.UnitLength_in_cm, 2);

          /* If accretion rate is zero, then accretion luminosity will be zero */
          if(mdotsink == 0.0)
            {
              sink_accretion_luminosity = 0.0;
            }
          else
            {
              /* Calculate current radius of protostar, using same approach as in Smith et al (2011).
               * This code was ported from the FLASH implementation by K. Wollenberg */

              rMS = 0.28 * pow(msink, 0.61); /* Radius (in Rsolar) once on main sequence */

              p1 = 5.0 * pow((mdotsink * 1e3), 0.27);
              p2 = 7.0 * pow((mdotsink * 1e3), 0.27);

              if(msink <= p1)
                {
                  rsink = 26.0 * pow(msink, 0.27) * pow(mdotsink * 1e3, 0.41);
                }
              else if(msink < p2)
                {
                  rsink = 26.0 * pow(p1, -2.73) * pow(mdotsink * 1e3, 0.41) * pow(msink, 3);
                }
              else
                {
                  rsink = 26.0 * pow(p1, -2.73) * pow(mdotsink * 1e3, 0.41) * pow(p2, 5) * pow(msink, -2);
                }

              /*  If accretion rate is very small, this can yield a sink radius that is less than the MS
               *  radius of the star, which is unphysical. Therefore, use rMS as a lower limit on the radius
               */
              if(rsink < rMS)
                rsink = rMS;

              msink *= SOLAR_MASS;
              mdotsink *= SOLAR_MASS / SEC_PER_YEAR;
              rsink *= solar_radius;

              sink_accretion_luminosity = GRAVITY * msink * mdotsink / rsink;
            }

          /* Given this accretion luminosity, contribution to total flux then follows from L / 4 pi r**2 */
          sink_flux += sink_accretion_luminosity / (4.0 * M_PI * d2sink);
        }
    }
#endif

  /* Check that the chemistry has sane values after the advection. */
  for(i = 0; i < SGCHEM_NUM_SPECIES; i++)
    {
      if(non_eq_abundances[i] < -1e-12)
        {
          printf("sgchem.c negative abundance from advection, species = %d abundance = %g\n", i, non_eq_abundances[i]);
          printf("sgchem.c: Setting abundance to +1e-20 (this might not help!)\n");
          non_eq_abundances[i] = 1e-20;
        }
      if(non_eq_abundances[IH2] > 0.5)
        non_eq_abundances[IH2] = 0.5;
      if(non_eq_abundances[IHP] > 1)
        non_eq_abundances[IHP] = 1;
#if CHEMISTRYNETWORK != 1
      if(non_eq_abundances[ICO] > carb_abund)
        non_eq_abundances[ICO] = carb_abund;
#endif
    }

#if defined(SGCHEM_RT) || defined(MRT_CHEM_SG)

  for(i = 0; i < NDIRECT; i++)
    {
      photorates_direct[i] = SphP[index].sxPhotonRates[i];
    }

#if CHEMISTRYNETWORK == 1
  /* No need to track photoelectric heating when running with primordial chemistry, as we know the rate must be zero */
  photorates_direct[5] = 0.0;
#endif
#else
  /* If we're not using Simplex, we set all of the rates to zero */
  for(i = 0; i < NDIRECT; i++)
    {
      photorates_direct[i] = 0.0;
    }
#endif

#if defined(SGCHEM_DEBUG_PARTICLE)
  if(P[index].ID == SGCHEM_DEBUG_PARTICLE)
    {
      printf("SGCHEM - DEGBUG PARTICLE: id %d mass %g dl %g n %g divv %g energy_cgs %g dust_temp %g gas temp estimate %g \n", id,
             P[index].Mass, dl, yn, divv, energy_cgs, SphP[index].DustTemp,
             SphP[index].Utherm * All.UnitEnergy_in_cgs / All.UnitMass_in_g / 1.5 / 8.25e7 * 2.);
      printf("SGCHEM - DEGBUG PARTICLE: Abundances %g %g %g %g %g %g %g %g %g \n", SphP[index].TracAbund[0], SphP[index].TracAbund[1],
             SphP[index].TracAbund[2], SphP[index].TracAbund[3], SphP[index].TracAbund[4], SphP[index].TracAbund[5],
             SphP[index].TracAbund[6], SphP[index].TracAbund[7], SphP[index].TracAbund[8]);
      for(i = 0; i < NPIX; i++)
        printf("SGCHEM - DEGBUG PARTICLE: TreeCol pix %d total column %g H2 column %g CO column %g \n", i,
               column_density_projection[i], column_density_projection_H2[i], column_density_projection_CO[i]);
#ifdef SINK_PARTICLES
      double dist_to_closest_sink = 1e33;
      for(i = 0; i < NSinksAllTasks; i++)
        {
          double rad_to_sink = sqrt(pow(SinkP[i].Pos[0] - P[index].Pos[0], 2) + pow(SinkP[i].Pos[1] - P[index].Pos[1], 2) +
                                    pow(SinkP[i].Pos[2] - P[index].Pos[2], 2));
          if(rad_to_sink < dist_to_closest_sink)
            dist_to_closest_sink = rad_to_sink;
        }
      printf("SGCHEM - DEGBUG PARTICLE: sinks present. Distance to closest sink: %g \n", dist_to_closest_sink);
#endif
    }
#endif

  EVOLVE_ABUNDANCES(&time, &dl, &yn, &divv, &energy_cgs, &redshift, non_eq_abundances, breakdown_of_thermal_rates,
                    column_density_projection, column_density_projection_H2, column_density_projection_CO, column_density_projection_C,
                    photorates_direct, &dust_temp, &sink_flux, &id);

#ifdef SGCHEM_TEMPERATURE_FLOOR
  abh2 = non_eq_abundances[IH2];
#if CHEMISTRYNETWORK == 1 || CHEMISTRYNETWORK == 2
  abe  = non_eq_abundances[IHP] + non_eq_abundances[IHEP] + 2.0 * non_eq_abundances[IHEPP];
#else
  abe = non_eq_abundances[IHP];
#endif

  ntot        = yn * (1.0 + ABHE + abe - abh2);
  double temp = GAMMA_MINUS1 * energy_cgs / (ntot * BOLTZMANN);

  if(temp < All.SGChemTemperatureFloor)
    energy_cgs = All.SGChemTemperatureFloor * (ntot * BOLTZMANN) / GAMMA_MINUS1;
#endif /*SGCHEM_TEMPERATURE_FLOOR*/

  /* Convert from cgs to AREPO code units */
  energy = energy_cgs * pow(All.UnitLength_in_cm, 3) / All.UnitEnergy_in_cgs; /* Energy density in Arepo code units */

  /* Set new dust temperature */
  SphP[index].DustTemp = dust_temp;

  /* Update evolved values */
  Utherm_new  = energy / rho;

#ifdef MRT_CHEM_SG
  Utherm_diff = Utherm_new - utherm_particle;
  SphP[index].Energy += a * a * Utherm_diff * P[index].Mass;
#else
  Utherm_diff = Utherm_new - SphP[index].Utherm;
  SphP[index].Utherm += Utherm_diff;
#ifdef MHD_THERMAL_ENERGY_SWITCH
  SphP[index].Etherm = SphP[index].Utherm * P[index].Mass;
#endif
  SphP[index].Energy += a * a * Utherm_diff * P[index].Mass;
#endif

#ifdef VARIABLE_GAMMA
  /* Update gamma for this cell */
  abh2 = non_eq_abundances[IH2];
#if CHEMISTRYNETWORK == 1 || CHEMISTRYNETWORK == 2
  abe  = non_eq_abundances[IHP] + non_eq_abundances[IHEP] + 2.0 * non_eq_abundances[IHEPP];
#elif CHEMISTRYNETWORK == 7
  abe = non_eq_abundances[IHP] + non_eq_abundances[IHEP];
#else
  abe = non_eq_abundances[IHP];
#endif
  ntot = yn * (1.0 + ABHE + abe - abh2);
  ekn  = energy_cgs / (ntot * BOLTZMANN);
  CALC_GAMMA(&abh2, &ekn, &gamma);

  SphP[index].GammaE = gamma;
  SphP[index].GammaC = gamma;
  gamma_minus1       = gamma - 1.0;
#else
  gamma_minus1 = GAMMA_MINUS1;
#endif

#ifdef TREECOLV2_VEL
  /* Compute square of thermal velocity of H atom. Note: used for velocity-dependent molecular shielding.
   * Note: value here is in code units */
  SphP[index].Vth2 = 2.0 * SphP[index].Utherm * gamma_minus1 * (1.0 + 4.0 * ABHE) /
                     (1.0 + ABHE - SphP[index].TracAbund[IH2] + SphP[index].TracAbund[IHP]);
#endif

  /* Update cooling time */
#ifdef SGCHEM_OUTPUT_COOLTIME

#ifdef MRT_CHEM_SG
  SphP[index].CoolTime = utherm_particle * dt / Utherm_diff;
#else
  SphP[index].CoolTime = SphP[index].Utherm * dt / Utherm_diff;
#endif

#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  SphP[index].A        = gamma_minus1 * SphP[index].Utherm / pow(SphP[index].Density * a3inv, GAMMA_MINUS1);
  SphP[index].Entropy  = log(SphP[index].A) * P[index].Mass;
#endif

  for(i = 0; i < SGCHEM_NUM_SPECIES; i++)
    {
      SphP[index].TracAbund[i]     = non_eq_abundances[i];
      SphP[index].MassTracAbund[i] = P[index].Mass * SphP[index].TracAbund[i];
    }
#if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5 || CHEMISTRYNETWORK == 7
  SphP[index].TracAbund[IHATOM] = 1.0 - 2.0 * SphP[index].TracAbund[IH2] - SphP[index].TracAbund[IHP];
  if(SphP[index].TracAbund[IHATOM] < 0.0)
    {
      SphP[index].TracAbund[IHATOM] = 0.0;
    }
  SphP[index].MassTracAbund[IHATOM] = SphP[index].TracAbund[IHATOM] * P[index].Mass;
#endif
#if CHEMISTRYNETWORK == 5 || CHEMISTRYNETWORK == 7
  SphP[index].TracAbund[ICP]        = carb_abund - SphP[index].TracAbund[ICO];
  if(SphP[index].TracAbund[ICP] < 0.0)
    {
      SphP[index].TracAbund[ICP] = 0.0;
    }
  SphP[index].MassTracAbund[ICP] = SphP[index].TracAbund[ICP] * P[index].Mass;

  /*if(index == 10000)
    printf("sgchem:Final single cell %d abundances %g %g %g\n", index, SphP[index].TracAbund[0], SphP[index].TracAbund[1],
    SphP[index].TracAbund[2]);*/
#endif
#if CHEMISTRYNETWORK == 7
  SphP[index].TracAbund[IHEATOM] = ABHE - SphP[index].TracAbund[IHEP];
  if(SphP[index].TracAbund[IHEATOM] < 0.0)
    {
      SphP[index].TracAbund[IHEATOM] = 0.0;
    }
  SphP[index].MassTracAbund[IHEATOM] = SphP[index].TracAbund[IHEATOM] * P[index].Mass;
#endif
#if CHEMISTRYNETWORK == 1
  SphP[index].TracAbund[IHATOM] = 1.0 - 2.0 * SphP[index].TracAbund[IH2] - SphP[index].TracAbund[IHP] - SphP[index].TracAbund[IHD];
  if(SphP[index].TracAbund[IHATOM] < 0.0)
    {
      SphP[index].TracAbund[IHATOM] = 0.0;
    }
  SphP[index].MassTracAbund[IHATOM] = SphP[index].TracAbund[IHATOM] * P[index].Mass;

  SphP[index].TracAbund[IHEATOM] = ABHE - SphP[index].TracAbund[IHEP] - SphP[index].TracAbund[IHEPP];
  if(SphP[index].TracAbund[IHEATOM] < 0.0)
    {
      SphP[index].TracAbund[IHEATOM] = 0.0;
    }
  SphP[index].MassTracAbund[IHEATOM] = SphP[index].TracAbund[IHEATOM] * P[index].Mass;

  SphP[index].TracAbund[IDATOM] = All.DeutAbund - SphP[index].TracAbund[IDP] - SphP[index].TracAbund[IHD];
  if(SphP[index].TracAbund[IDATOM] < 0.0)
    {
      SphP[index].TracAbund[IDATOM] = 0.0;
    }
  SphP[index].MassTracAbund[IDATOM] = SphP[index].TracAbund[IDATOM] * P[index].Mass;
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  SphP[index].TracAbund[IHATOM] = 1.0 - 2.0 * SphP[index].TracAbund[IH2] - SphP[index].TracAbund[IHP] - SphP[index].TracAbund[ICHX] -
                                  SphP[index].TracAbund[IOHX] - SphP[index].TracAbund[IHCOP];

  SphP[index].TracAbund[ICATOM] = carb_abund - SphP[index].TracAbund[ICP] - SphP[index].TracAbund[ICHX] - SphP[index].TracAbund[ICO] -
                                  SphP[index].TracAbund[IHCOP];

  SphP[index].TracAbund[IOATOM] = oxy_abund - SphP[index].TracAbund[IOHX] - SphP[index].TracAbund[ICO] - SphP[index].TracAbund[IHCOP];

  SphP[index].TracAbund[IHEATOM] = ABHE - SphP[index].TracAbund[IHEP];

  SphP[index].TracAbund[IMATOM] = m_abund - SphP[index].TracAbund[IMP]; /* Actually Si for network 16 */

  for(i = 0; i < SGCHEM_NUM_ADVECTED_SPECIES; i++)
    {
      if(SphP[index].TracAbund[i] < 0)
        {
          SphP[index].TracAbund[i] = 0.0;
        }
      SphP[index].MassTracAbund[i] = SphP[index].TracAbund[i] * P[index].Mass;
    }
#endif

#ifdef SGCHEM_DUMP_THERMAL_RATES
  /* Pass out the heating / cooling rate block to the main code for writing to snapshots */
  for(i = 0; i < SGCHEM_NUM_THERMAL_RATES; i++)
    SphP[index].HeatCoolRates[i] = breakdown_of_thermal_rates[i];
#endif

    /* Update pressure */
#ifndef MRT_CHEM_SG
  set_pressure_of_cell(index);
#endif

  return;
#endif /* NO_CHEMISTRY */
}

#ifdef SGCHEM_VARIABLE_ISRF
/* Compute local strength of ISRF, using a simple Galactic model from Wolfire et al 2003 */
void find_and_set_local_ISRF(int index)
{
  double local_ISRF;

  double xi = P[index].Pos[0] - boxHalf_X;
  double yi = P[index].Pos[1] - boxHalf_Y;
  double zi = P[index].Pos[2] - boxHalf_Z;

  double R = sqrt(xi * xi + yi * yi) * All.UnitLength_in_cm / KILOPARSEC;

  /* Radial dependence of the ISRF is taken from Wolfire+ 2003 */
  double R1  = 4.;
  double HRJ = 4.1;
  double R0  = 8.;

  if(R < R1)
    local_ISRF = All.UVFieldStrength * exp(-(R1 - R0) / HRJ);
  else
    local_ISRF = All.UVFieldStrength * exp(-(R - R0) / HRJ);

  SET_LOCAL_ISRF(&local_ISRF);
  return;
}
#endif

#ifdef SGCHEM_VARIABLE_CRION
/* Compute local cosmic ray ionization rate. If we're using the COSMIC_RAYS option, we scale this according to the
 * cosmic ray energy density; otherwise, we use a simple Galactic scaling model from Wolfire et al 2003
 */
void find_and_set_local_CRION(int index)
{
  double local_CRIR, cr_energy_density;

#ifdef COSMIC_RAYS
  cr_energy_density = SphP[index].CR_SpecificEnergy * SphP[index].Density;
  cr_energy_density *=
      All.cf_a2inv * All.cf_a2inv * All.HubbleParam * All.HubbleParam * All.UnitEnergy_in_cgs / pow(All.UnitLength_in_cm, 3);
  local_CRIR = 0;
#else
  double xi = P[index].Pos[0] - boxHalf_X;
  double yi = P[index].Pos[1] - boxHalf_Y;
  double zi = P[index].Pos[2] - boxHalf_Z;

  double R = sqrt(xi * xi + yi * yi) * All.UnitLength_in_cm / KILOPARSEC;

  /*radial dependence of the CRIR is taken from Wolfire+ 2003*/
  double R0 = 8.;
  double HROB = 3.5;

  local_CRIR = All.CosmicRayIonRate * exp(-(R - R0) / HROB);
  cr_energy_density = 0;
#endif
  SET_LOCAL_CR_ION_RATE(&cr_energy_density, &local_CRIR);
}
#endif
