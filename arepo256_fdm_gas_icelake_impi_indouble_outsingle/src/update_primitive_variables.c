/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/update_primitive_variables.c
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

#include <gsl/gsl_linalg.h>

#include "allvars.h"
#include "helm_eos.h"
#include "opal_eos.h"
#include "proto.h"

/*! \brief Main routine to update the primitive hydrodynamics variables from
 *         the conserved ones.
 *
 *  Note that the primitive variables are inconsistent with the (new)
 *  conserved variables after the hydro integration up to the point this
 *  function is called.
 *
 *  \return void
 */
void update_primitive_variables(void)
{
  TIMER_START(CPU_CELL_UPDATES);

#ifdef COSMIC_RAYS_EXTRA_DIAGNOSTICS
  double InitialCREnergy = 0;
  for(int i = 0; i < NumGas; i++)
    if(P[i].Mass != 0 && P[i].ID != 0 && P[i].Type == 0)
      InitialCREnergy += SphP[i].CR_Energy;
#endif

  struct pv_update_data pvd;
  int idx, i;

  if(All.ComovingIntegrationOn)
    {
      pvd.atime    = All.Time;
      pvd.hubble_a = hubble_function(All.Time);
      pvd.a3inv    = 1 / (All.Time * All.Time * All.Time);
    }
  else
    pvd.atime = pvd.hubble_a = pvd.a3inv = 1.0;

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  pvd.count_keep_entropy = pvd.count_update_entropy = 0;
#endif

#pragma omp parallel for private(idx, i)
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      do_validity_checks(P, SphP, i, &pvd);

      do_special_boundaries(P, SphP, i, &pvd);

#ifdef SPECIAL_RELATIVITY
      update_primitive_variables_special_relativity(P, SphP, i, &pvd);
#else
#ifdef GENERAL_RELATIVITY
      update_primitive_variables_general_relativity(P, SphP, i, &pvd);
#else
      update_primitive_variables_single(P, SphP, i, &pvd);

      update_internal_energy(P, SphP, i, &pvd);

      set_pressure_of_cell_internal(P, SphP, i); /* calculate the pressure from Density and Utherm (and composition) */
#endif
#endif

#ifdef SGCHEM
      do_chemical_abund_checks(P, SphP, i, &pvd);
#endif

      SphP[i].OldMass = P[i].Mass;

#ifdef MRT_LSF_GRADIENTS
      //      for(int num1=0;num1<MRT_BINS;num1++)
      // SphP[i].OldCons_DensPhot[num1] = SphP[i].Cons_DensPhot[num1] ;
#ifdef MRT_COMOVING
      SphP[i].Old_Vel[0] = P[i].Vel[0];
      SphP[i].Old_Vel[1] = P[i].Vel[1];
      SphP[i].Old_Vel[2] = P[i].Vel[2];
#endif
#endif

      SphP[i].TimeLastPrimUpdate = All.Time;
    }

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  long long tot_count_keep_entropy, tot_count_update_entropy;

  sumup_large_ints(1, &pvd.count_keep_entropy, &tot_count_keep_entropy);
  sumup_large_ints(1, &pvd.count_update_entropy, &tot_count_update_entropy);

  mpi_printf("Keep-Entropy=%llu  Update-Entropy=%llu\n", tot_count_keep_entropy, tot_count_update_entropy);
#endif

#ifdef COSMIC_RAYS_EXTRA_DIAGNOSTICS
  double FinalCREnergy = 0;
  for(int i = 0; i < NumGas; i++)
    if(P[i].Mass != 0 && P[i].ID != 0 && P[i].Type == 0)
      FinalCREnergy += SphP[i].CR_Energy;

  double dCREnergy = FinalCREnergy - InitialCREnergy;
  All.TotalCREnergyUpdatePrims += dCREnergy;
#endif

  TIMER_STOP(CPU_CELL_UPDATES);
}

/*! \brief Wrapper function to calculate pressure of a cell from its internal
 *         energy.
 *
 *  \param[in] i Index of cell in P and SphP arrays.
 *
 *  \return void
 */
void set_pressure_of_cell(int i) { set_pressure_of_cell_internal(P, SphP, i); }

/*! \brief Function to calculate pressure from other hydrodynamics quantities.
 *
 *  How this is done depends on the adiabatic index and potentially on sub-
 *  resolution physics. Note that this is just the thermal pressure (i.e. not
 *  including magnetic fields).
 *
 *  \param[in] localP Pointer to particle data array.
 *  \param[in,out] localSphP Pointer to cell data array.
 *  \param[in] i Index in localP and localSphP arrays.
 *
 *  \return void
 */
void set_pressure_of_cell_internal(struct particle_data *localP, struct sph_particle_data *localSphP, int i)
{
#ifdef DG
  localSphP[i].Pressure = GAMMA_MINUS1 * localSphP[i].Density * localSphP[i].Utherm;

#ifndef FIX_MEAN_VALUES
  if(localSphP[i].Pressure < Epsilon_p)
    {
      printf("cell: %d, pressure: %g\n", i, localSphP[i].Pressure);
      terminate("Negative pressure detected!\n");
    }
#endif

  return;
#endif

#ifdef ISOTHERM_EQS
  localSphP[i].Pressure = localSphP[i].Density * All.IsoSoundSpeed * All.IsoSoundSpeed;
#else

#ifdef LOCALLY_ISOTHERM_DISK
  double csnd;
  csnd = get_isotherm_disk_sound_speed(i);
  localSphP[i].Pressure = localSphP[i].Density * csnd * csnd;
#else

  if(localSphP[i].Utherm >= 0)
    {
#if !defined(USE_SFR) || defined(LOCAL_FEEDBACK) || defined(SFR_MCS)

#ifdef VARIABLE_GAMMA
      localSphP[i].GammaE = GAMMA;
      localSphP[i].GammaC = GAMMA;
#endif

#ifdef TGCHEM
      localSphP[i].Pressure = (localSphP[i].Gamma - 1) * localSphP[i].Density * localSphP[i].Utherm;
#else
#ifdef EOS_DEGENERATE
      struct eos_result res;

      if(eos_calc_egiven(localSphP[i].Density, localSphP[i].Composition, localSphP[i].Utherm, &localSphP[i].EOSTemperature, &res) ==
         -2)
        {
          localSphP[i].EOSTemperature = -1.;
          if(eos_calc_egiven(localSphP[i].Density, localSphP[i].Composition, localSphP[i].Utherm, &localSphP[i].EOSTemperature,
                             &res) != 0)
            {
              if(res.temp < 1000.)
                localSphP[i].EOSTemperature = 1000;
              else if(res.temp > 1e10)
                localSphP[i].EOSTemperature = 1e10;

              eos_calc_tgiven(localSphP[i].Density, localSphP[i].Composition, localSphP[i].EOSTemperature, &res);
            }
        }

      if(localSphP[i].EOSTemperature <= 1.001 * 1000. || localSphP[i].EOSTemperature >= 1e10 * 0.999)
        {
          /* looks like the equation of state corrected the internal energy, as we run out of the table */
          EgyInjection += res.e.v - localSphP[i].Utherm;
          localSphP[i].Energy += (res.e.v - localSphP[i].Utherm) * localP[i].Mass;
          localSphP[i].Utherm = res.e.v;
        }

#ifdef OUTPUT_ENTROPY
      localSphP[i].Entropy = res.s.v;
#endif

      localSphP[i].Pressure = res.p.v;
      localSphP[i].GammaE   = localSphP[i].Pressure / localSphP[i].Utherm / localSphP[i].Density + 1.0;
      localSphP[i].GammaC   = (res.p.drho + res.temp * gsl_pow_2(res.p.dtemp / localSphP[i].Density) / res.e.dtemp) *
                            localSphP[i].Density / localSphP[i].Pressure;
      localSphP[i].cv = res.cv;

      if(localSphP[i].GammaC <= 0)
        {
          printf("Pressure=%g, Temperature=%g, Density=%g, dpdr=%g, dedt=%g.\n", res.p.v, res.temp, localSphP[i].Density, res.p.drho,
                 res.e.dtemp);
          print_particle_info(i);
          terminate("bla");
        }
#else
#ifdef EOS_OPAL
      struct opal_eos_result res;

      /* limit and renormalize composition */
      double xsum = 0.0;
      int j;
      /* X can be 0.8 maximum */
      localSphP[i].Composition[0] = fmin(0.8, localSphP[i].Composition[0]);
      for(j = 0; j < EOS_NSPECIES; j++)
        {
          localSphP[i].Composition[j] = fmin(1.0, fmax(0.0, localSphP[i].Composition[j]));
          xsum += localSphP[i].Composition[j];
        }
      for(j = 0; j < EOS_NSPECIES; j++)
        localSphP[i].Composition[j] /= xsum;

      /* call EOS */
      if(opaleos_egiven(localSphP[i].Composition[0], localSphP[i].Density, localSphP[i].Utherm, &localSphP[i].EOSTemperature,
                        EOS_PRESSURE | EOS_ENERGY | EOS_CHIT | EOS_CHIR | EOS_CV | EOS_GAMMA1 | EOS_RADIATION, &res) < 0)
        {
          /* printf("Correcting internal energy from %e to %e (total energy: %e;\
            rel. change: %e).\n", localSphP[i].Utherm, res.e, localSphP[i].Energy, res.e / localSphP[i].Energy); */
          localSphP[i].Utherm = res.e;
        }

      localSphP[i].Pressure = res.p;

      localSphP[i].GammaE = localSphP[i].Pressure / localSphP[i].Utherm / localSphP[i].Density + 1.0;
      localSphP[i].GammaC = res.gamma1;

      /* check results */
      if(!gsl_finite(localSphP[i].GammaC) || !gsl_finite(localSphP[i].GammaE))
        {
          printf("Error in EOS for particle %d.\n", i);
          print_particle_info(i);
          terminate("Infinity encountered in EOS.");
        }

#else
#if defined(SGCHEM) && defined(VARIABLE_GAMMA)
      double abh2, abe, energy, ntot, ekn, gamma, rho;
      double a, a3inv, hubble_param2;

      abh2 = localSphP[i].TracAbund[IH2];
      abe  = localSphP[i].TracAbund[IHP];
#if CHEMISTRYNETWORK == 1 || CHEMISTRYNETWORK == 2
      abe += localSphP[i].TracAbund[IHEP] + 2.0 * localSphP[i].TracAbund[IHEPP];
#endif
#if CHEMISTRYNETWORK == 7
      abe += localSphP[i].TracAbund[IHEP];
#endif
      if(All.ComovingIntegrationOn)
        {
          a             = All.Time;
          a3inv         = 1.0 / (a * a * a);
          hubble_param2 = All.HubbleParam * All.HubbleParam;
        }
      else
        {
          a = a3inv = hubble_param2 = 1;
        }

      /* Convert internal energy to cgs units */
      rho    = localSphP[i].Density * a3inv * hubble_param2;
      energy = SphP[i].Utherm * rho * All.UnitEnergy_in_cgs / pow(All.UnitLength_in_cm, 3);

      /* Get total particle number density */
      ntot = rho * All.UnitDensity_in_cgs / ((1.0 + 4.0 * HE_ABUND) * PROTONMASS);
      ntot *= (1.0 + HE_ABUND + abe - abh2);

      ekn = energy / (ntot * BOLTZMANN);
      if(abh2 < 0 || abh2 > 0.5)
        {
          printf("WARNING: ID %d abh2=%.3e ekn=%.3e ", P[i].ID, abh2, ekn);
          CALC_GAMMA(&abh2, &ekn, &gamma);
          printf(" gamma=%.3e\n", gamma);
        }
      else
        {
          CALC_GAMMA(&abh2, &ekn, &gamma);
        }

      localSphP[i].GammaE   = gamma;
      localSphP[i].GammaC   = gamma;
      localSphP[i].Pressure = (gamma - 1.0) * localSphP[i].Density * localSphP[i].Utherm;
#else
      /* this block applies to ordinary gas dynamics, without SFR */
      localSphP[i].Pressure = GAMMA_MINUS1 * localSphP[i].Density * localSphP[i].Utherm;
#endif /* SGCHEM && VARIABLE_GAMMA */

#if defined(MRT_IR_PHOTON_TRAPPING) && defined(MRT_RADIATION_PRESSURE)
      for(int num1 = UV_BINS; num1 < (UV_BINS + IR_BINS); num1++)
        localSphP[i].Pressure +=
            localSphP[i].Trapped_DensPhot[num1] * c_internal_units / (2.99792458e10 / All.UnitVelocity_in_cm_per_s) / 3.0;
#endif

#endif  /* EOS_OPAL */
#endif  /* EOS_DEGENERATE */
#endif  // TGCHEM

#ifdef JEANS_PRESSURE_LIMIT
      double jeans_pressure, celldim, ncells;
      ncells = JEANS_PRESSURE_LIMIT;
      celldim = 2.0 * pow(All.MinVolume * 3.0 / 4.0 / M_PI / 2.0, 1.0 / 3.0);
#ifdef VARIABLE_GAMMA
      jeans_pressure =
          ncells * ncells * All.G * celldim * celldim * localSphP[i].Density * localSphP[i].Density / M_PI / localSphP[i].GammaE;
#else
      jeans_pressure = ncells * ncells * All.G * celldim * celldim * localSphP[i].Density * localSphP[i].Density / M_PI / GAMMA;
#endif
      localSphP[i].Pressure = fmax(localSphP[i].Pressure, jeans_pressure);

#ifdef MAKE_PRES_UTHERM_CONSISTENT
#ifdef VARIABLE_GAMMA
      double utherm_from_pres = localSphP[i].Pressure / (localSphP[i].GammaE - 1.0) / localSphP[i].Density;
#else
      double utherm_from_pres = localSphP[i].Pressure / GAMMA_MINUS1 / localSphP[i].Density;
#endif
      double utherm_diff = utherm_from_pres - localSphP[i].Utherm;
      localSphP[i].Utherm += utherm_diff;
      localSphP[i].Energy += utherm_diff * localP[i].Mass;
#endif

#endif

#ifdef JEANS_TOTPRESSURE_LIMIT

      double jeans_pressure, celldim, ncells;
      ncells = JEANS_PRESSURE_LIMIT;
      // celldim = 2.0 * get_cell_radius(i);
      celldim = 2.0 * pow(All.MinVolume * 3.0 / 4.0 / M_PI / 2.0, 1.0 / 3.0);
#ifdef VARIABLE_GAMMA
      jeans_pressure =
          ncells * ncells * All.G * celldim * celldim * localSphP[i].Density * localSphP[i].Density / M_PI / localSphP[i].GammaE;

      double tot_pressure = (localSphP[i].GammaE - 1.0) * localSphP[i].Energy / localSphP[i].Volume;
#else
      jeans_pressure = ncells * ncells * All.G * celldim * celldim * localSphP[i].Density * localSphP[i].Density / M_PI / GAMMA;

      double tot_pressure   = GAMMA_MINUS1 * localSphP[i].Energy / localSphP[i].Volume;
#endif
      tot_pressure = fmax(jeans_pressure, tot_pressure);
#endif

#ifdef JEANS_PRESSURE_LIMIT_MCS
      double jeans_pressure;
      double ncells = JEANS_PRESSURE_LIMIT_MCS;
#ifndef JEANS_MASS_PRESSURE_LIMIT_MCS
      double celldim = 2.0 * get_cell_radius(i);
      jeans_pressure =
          ncells * ncells * All.G * celldim * celldim * localSphP[i].Density * localSphP[i].Density * All.cf_atime / M_PI / GAMMA;
#else
      double mass_rho_term  = ncells * localP[i].Mass * localSphP[i].Density * localSphP[i].Density;
      jeans_pressure        = (All.G / GAMMA) * cbrt(6.0 * mass_rho_term * mass_rho_term / pow(M_PI, 5.0));
#endif
#if defined(GRACKLE) && !defined(GRACKLE_TAB)
      localSphP[i].Pressure = fmax(get_grackle_pressure(i), jeans_pressure);
#else
      localSphP[i].Pressure = fmax(GAMMA_MINUS1 * localSphP[i].Density * localSphP[i].Utherm, jeans_pressure);
#endif
#endif  // JEANS_PRESSURE_LIMIT_MCS

#ifdef SMAUG_PRESSURE_FLOOR
      /* Standard thermal pressure */
      localSphP[i].Pressure = GAMMA_MINUS1 * localSphP[i].Density * localSphP[i].Utherm;
      /* Now add polytrope contribution, where T_poly / mu = Tstar * (n / nstar)^(gstar - 1). PolytropeFactor
      contains all constants, calculated at start of run. Note that these pressures are added (rather
      than taking max), this is equivalent to Ramses prescription */
      localSphP[i].Pressure += All.PolytropeFactor * pow(localSphP[i].Density, All.Polytrope_gstar);
#endif

#else /* now comes the SFR treatment */

      localSphP[i].Pressure = GAMMA_MINUS1 * localSphP[i].Density * localSphP[i].Utherm;

#endif  // USE_SFR
    }  // end utherm >= 0
  else
    localSphP[i].Pressure = 0;
#endif  // LOCAL_ISOTHERM_DISK
#endif  // ISOTHERM_EQS

#ifdef COSMIC_RAYS
  localSphP[i].CR_Pressure = (All.GammaCR - 1.0) * localSphP[i].CR_SpecificEnergy * localSphP[i].Density / All.cf_atime;
#endif

#if defined(GRACKLE) && !defined(GRACKLE_TAB) && !defined(JEANS_PRESSURE_LIMIT_MCS)
  localSphP[i].Pressure = get_grackle_pressure(i);
#endif

#ifdef ENFORCE_JEANS_STABILITY_OF_CELLS
#ifndef ENFORCE_JEANS_STABILITY_OF_CELLS_EEOS
#if defined(USE_SFR) && !defined(SMUGGLE_SFR) && !defined(LOCAL_FEEDBACK)
  if(get_starformation_rate(i) == 0)
#endif
    {
#endif

#ifdef ADAPTIVE_HYDRO_SOFTENING
      double cell_soft = All.ForceSoftening[localP[i].SofteningType];
#else
  double cell_soft = All.GasSoftFactor * get_cell_radius(i);
#endif

#if defined(REFINEMENT_HIGH_RES_GAS) && !defined(TGSET)
      if(SphP[i].HighResMass > HIGHRESMASSFAC * P[i].Mass)
#endif
        localSphP[i].Pressure =
            fmax(localSphP[i].Pressure, GAMMA_MINUS1 * localSphP[i].Density * 2 * All.G * localP[i].Mass / (All.cf_atime * cell_soft));

#ifndef ENFORCE_JEANS_STABILITY_OF_CELLS_EEOS
    }
#endif
#endif

#ifdef SGS_TURBULENCE
  localSphP[i].SgsTData.Pressure = (All.SgsTConst.GammaSgsT - 1.0) * localSphP[i].SgsTData.SpecificEnergy * localSphP[i].Density;
#endif
}

/*! \brief Validity checks for a gas cell.
 *
 *  So far, only a positive mass constraint implemented. Terminates if not
 *  successful.
 *
 *  \param[in] localP Pointer to particle data array
 *  \param[in,out] localSphP Pointer to cell data array
 *  \param[in] i Index in localP and localSphP arrays
 *  \param[in] pvd (unused)
 *
 *  \return void
 */
void do_validity_checks(struct particle_data *localP, struct sph_particle_data *localSphP, int i, struct pv_update_data *pvd)
{
  if(localP[i].Mass < 0)
    {
      printf("very bad...i=%d ID=%d mass=%g oldMass=%g utherm=%g pos=%g|%g|%g\n", i, (int)localP[i].ID, localP[i].Mass,
             localSphP[i].OldMass, localSphP[i].Utherm, localP[i].Pos[0], localP[i].Pos[1], localP[i].Pos[2]);

      terminate("stop");
    }
}

#ifdef SGCHEM
void do_chemical_abund_checks(struct particle_data *localP, struct sph_particle_data *localSphP, int index, struct pv_update_data *pvd)
{
  double non_eq_abundances_i, carb_abund, oxy_abund;
  int i;

#ifdef SGCHEM_VARIABLE_Z
  carb_abund = localSphP[index].CarbAbund;
  oxy_abund  = localSphP[index].OxyAbund;
#else
  carb_abund       = All.CarbAbund;
  oxy_abund        = All.OxyAbund;
#endif

  // Check that the chemistry has sane values after the advection.
  for(i = 0; i < SGCHEM_NUM_SPECIES; i++)
    {
      non_eq_abundances_i = SphP[index].TracAbund[i];

      if(non_eq_abundances_i < 0.0)
        {
#ifdef DEBUG_SGCHEM
          printf("update_primitive_variables.c negative abundance from advection, species = %d abundance = %g\n", i,
                 non_eq_abundances_i);
          printf("update_primitive_variables.c: Setting abundance to +1e-20 (this might not help!)\n");
#endif
          non_eq_abundances_i = 1e-20;
        }

      if(i == IH2 && non_eq_abundances_i > 0.5)
        {
#ifdef DEBUG_SGCHEM
          printf("update_primitive_variables.c H2 abundance greater than 0.5; abundance = %g\n", non_eq_abundances_i);
          non_eq_abundances_i = 0.5;
#endif
        }

      if(i == IHP && non_eq_abundances_i > 1.0)
        {
#ifdef DEBUG_SGCHEM
          printf("update_primitive_variables.c HP abundance greater than 1.0; abundance = %g\n", non_eq_abundances_i);
#endif
          non_eq_abundances_i = 1.0;
        }

      if(i == ICO && non_eq_abundances_i > fmin(carb_abund, oxy_abund))
        {
#ifdef DEBUG_SGCHEM
          printf("update_primitive_variables.c CO abundance greater than Carbon & Oxygen abundances; abundance = %g\n",
                 non_eq_abundances_i);
#endif
          non_eq_abundances_i = fmin(carb_abund, oxy_abund);
        }

#if CHEMISTRYNETWORK == 1
      if(i == IHEPP && non_eq_abundances_i + SphP[index].TracAbund[IHEP] > ABHE)
        {
#ifdef DEBUG_SGCHEM
          printf("update_primitive_variables.c Too much He+/He++; abundances = %g / %g \n", SphP[index].TracAbund[IHEP],
                 non_eq_abundances_i);
#endif
          non_eq_abundances_i = ABHE - SphP[index].TracAbund[IHEP];
        }
#endif
      SphP[index].TracAbund[i]     = non_eq_abundances_i;
      SphP[index].MassTracAbund[i] = P[index].Mass * SphP[index].TracAbund[i];
    }
}
#endif

/*! \brief Update the gas boundary cells if the boundary is special.
 *
 *  For now unused.
 *
 *  \param[in] localP Pointer to particle data array
 *  \param[in,out] localSphP Pointer to cell data array
 *  \param[in] i Index in localP and localSphP arrays
 *  \param[in] pvd additional data that is needed for update (e.g. cosmological
 *             factors)
 *
 *  \return void
 */
void do_special_boundaries(struct particle_data *localP, struct sph_particle_data *localSphP, int i, struct pv_update_data *pvd)
{
#ifdef BOUNDARY_FLAG
  if(localSphP[i].BoundaryFlag & 2)
    {
      localP[i].Mass = localSphP[i].Volume * All.AGBWindDensity;

      double d[3];
      d[0]        = localSphP[i].Center[0] - All.AGBWindCenterX;
      d[1]        = localSphP[i].Center[1] - All.AGBWindCenterY;
      d[2]        = localSphP[i].Center[2] - All.AGBWindCenterZ;
      double dist = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);

      for(int k = 0; k < 3; k++)
        localSphP[i].Momentum[k] = localP[i].Mass * All.AGBWindVelocity * d[k] / dist;

      localSphP[i].Energy = localP[i].Mass * All.AGBWindSpecificEnergy + 0.5 / localP[i].Mass *
                                                                             (localSphP[i].Momentum[0] * localSphP[i].Momentum[0] +
                                                                              localSphP[i].Momentum[1] * localSphP[i].Momentum[1] +
                                                                              localSphP[i].Momentum[2] * localSphP[i].Momentum[2]);
    }
#endif

#ifdef SPECIAL_BOUNDARY
  double dt;

  // ID <-3 cells are buffer cells and are never updated
  if((localP[i].ID <= -3) && (All.Time > All.TimeBegin))
    return;

  if((localP[i].ID == -2) && (All.SpecialBoundaryType == 3 || All.SpecialBoundaryType == 4))
    {
      dt = (localP[i].TimeBinHydro ? (((integertime)1) << localP[i].TimeBinHydro) : 0) * All.Timebase_interval;
      boundary_get_velocity(localP[i].Pos[0], localP[i].Pos[1], localP[i].Pos[2], &localP[i].Vel[0], &localP[i].Vel[1],
                            &localP[i].Vel[2], dt);
      localSphP[i].Momentum[0] = localP[i].Vel[0] * localP[i].Mass;
      localSphP[i].Momentum[1] = localP[i].Vel[1] * localP[i].Mass;
      localSphP[i].Momentum[2] = localP[i].Vel[2] * localP[i].Mass;
      /*update vectorial quantities only, not the scalars */
    }
#endif

#ifdef WINDTUNNEL_READ_IN_BFIELD
  int nx, ny, nz;
  float x_d, y_d, z_d;
  double Bx, By, Bz;
  double P_total_cgs, P_hydro_cgs, Bmag;
  double tmp;
#endif

#ifdef WINDTUNNEL_FIXVARIABLESININJECTIONREGION

  if(localP[i].Pos[WINDTUNNEL_COORD] < All.InjectionRegion)
    {
      localSphP[i].Density            = All.InjectionDensity;
      localP[i].Vel[0]                = 0;
      localP[i].Vel[1]                = 0;
      localP[i].Vel[2]                = 0;
      localP[i].Vel[WINDTUNNEL_COORD] = All.InjectionVelocity;
      localSphP[i].Utherm             = All.InjectionUtherm;

      localP[i].Mass           = localSphP[i].Density * localSphP[i].Volume;
      localSphP[i].Momentum[0] = localP[i].Vel[0] * localP[i].Mass;
      localSphP[i].Momentum[1] = localP[i].Vel[1] * localP[i].Mass;
      localSphP[i].Momentum[2] = localP[i].Vel[2] * localP[i].Mass;
      localSphP[i].Energy =
          localP[i].Mass * pvd->atime * pvd->atime * localSphP[i].Utherm +
          0.5 * localP[i].Mass *
              (localP[i].Vel[0] * localP[i].Vel[0] + localP[i].Vel[1] * localP[i].Vel[1] + localP[i].Vel[2] * localP[i].Vel[2]);

      localSphP[i].Pressure = GAMMA_MINUS1 * localSphP[i].Density * localSphP[i].Utherm;

#ifdef MHD
#ifdef WINDTUNNEL_READ_IN_BFIELD

      WindtunnelReadIn_CalculateIndex(localP[i].Pos[0], localP[i].Pos[1], localP[i].Pos[2], &nx, &ny, &nz, &x_d, &y_d, &z_d);
      Bx = WindtunnelReadIn_TrilinearInterpolation(All.WindtunnelReadIn_Bx, nx, ny, nz, x_d, y_d, z_d);
      By = WindtunnelReadIn_TrilinearInterpolation(All.WindtunnelReadIn_By, nx, ny, nz, x_d, y_d, z_d);
      Bz = WindtunnelReadIn_TrilinearInterpolation(All.WindtunnelReadIn_Bz, nx, ny, nz, x_d, y_d, z_d);

      Bx *= sqrt(All.InjectionBx_InGauss * All.InjectionBx_InGauss + All.InjectionBy_InGauss * All.InjectionBy_InGauss +
                 All.InjectionBz_InGauss * All.InjectionBz_InGauss);
      By *= sqrt(All.InjectionBx_InGauss * All.InjectionBx_InGauss + All.InjectionBy_InGauss * All.InjectionBy_InGauss +
                 All.InjectionBz_InGauss * All.InjectionBz_InGauss);
      Bz *= sqrt(All.InjectionBx_InGauss * All.InjectionBx_InGauss + All.InjectionBy_InGauss * All.InjectionBy_InGauss +
                 All.InjectionBz_InGauss * All.InjectionBz_InGauss);

      // Enforce pressure equilibrium such that magnetic pressure plus hydro pressure is constant (constant is determined based on
      // injection region values). We determine the "target total pressure" based on the parameters for the injection region
      P_total_cgs = 1.0 / (8.0 * M_PI) *
                        (All.InjectionBx_InGauss * All.InjectionBx_InGauss + All.InjectionBy_InGauss * All.InjectionBy_InGauss +
                         All.InjectionBz_InGauss * All.InjectionBz_InGauss) +
                    GAMMA_MINUS1 * (All.InjectionDensity * All.UnitDensity_in_cgs) *
                        (All.InjectionUtherm * All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s);  // P_Magnetic + P_Hydro

      // We determine the actual hydro pressure in cgs units
      P_hydro_cgs = P_total_cgs - 1.0 / (8.0 * M_PI) * (Bx * Bx + By * By + Bz * Bz);

      // we do not allow magnetic fields stronger than PlasmaBeta=1 for any cell. If PlasmaBeta<1 we se P(magnetic)=P(hydro)
      if(1.0 / (8.0 * M_PI) * (Bx * Bx + By * By + Bz * Bz) > 0.5 * P_total_cgs)
        {
          Bmag = sqrt(Bx * Bx + By * By + Bz * Bz);
          Bx *= sqrt(0.5 * P_total_cgs * 8.0 * M_PI) / Bmag;
          By *= sqrt(0.5 * P_total_cgs * 8.0 * M_PI) / Bmag;
          Bz *= sqrt(0.5 * P_total_cgs * 8.0 * M_PI) / Bmag;
          P_hydro_cgs = 0.5 * P_total_cgs;

          // printf("P123 %f %d %20.20lf %20.20lf %20.20lf\n",All.Time,i,P_hydro_cgs,1.0/(8.0*M_PI)*(Bx*Bx+By*By+Bz*Bz),P_total_cgs);
        }
      // else
      // printf("P234 %d %20.20lf %20.20lf %20.20lf\n",i,P_hydro_cgs,1.0/(8.0*M_PI)*(Bx*Bx+By*By+Bz*Bz),P_total_cgs);

      localSphP[i].Pressure = P_hydro_cgs / All.UnitDensity_in_cgs / (All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s);

      // We maximum allow a factor of 2 change in Utherm:
      if(localSphP[i].Pressure / (GAMMA_MINUS1 * localSphP[i].Density) < 0.5 * localSphP[i].Utherm)
        {
          localSphP[i].Utherm = 0.5 * localSphP[i].Utherm;
          mpi_printf("P123 %d %20.20lf \n", i, localSphP[i].Utherm);
          terminate("Utherm is too small. Should not happen");
        }
      else if(localSphP[i].Pressure / (GAMMA_MINUS1 * localSphP[i].Density) > 2.0 * localSphP[i].Utherm)
        {
          localSphP[i].Utherm = 2.0 * localSphP[i].Utherm;
          mpi_printf("P234 %d %20.20lf \n", i, localSphP[i].Utherm);
          terminate("Utherm is too large. Should not happen");
        }
      else
        {
          localSphP[i].Utherm = localSphP[i].Pressure / (GAMMA_MINUS1 * localSphP[i].Density);
        }

      // Overwrite energy with new value of Utherm:
      localSphP[i].Energy =
          localP[i].Mass * pvd->atime * pvd->atime * localSphP[i].Utherm +
          0.5 * localP[i].Mass *
              (localP[i].Vel[0] * localP[i].Vel[0] + localP[i].Vel[1] * localP[i].Vel[1] + localP[i].Vel[2] * localP[i].Vel[2]);

      // in the following bfac is written out with the same formula as in init.c, where it is defined with "double bfac = 1. /
      // (sqrt(All.UnitMass_in_g / All.UnitLength_in_cm) / (All.UnitTime_in_s / All.HubbleParam))";
      localSphP[i].B[0] =
          Bx / sqrt(4. * M_PI) / (sqrt(All.UnitMass_in_g / All.UnitLength_in_cm) / (All.UnitTime_in_s / All.HubbleParam));
      localSphP[i].B[1] =
          By / sqrt(4. * M_PI) / (sqrt(All.UnitMass_in_g / All.UnitLength_in_cm) / (All.UnitTime_in_s / All.HubbleParam));
      localSphP[i].B[2] =
          Bz / sqrt(4. * M_PI) / (sqrt(All.UnitMass_in_g / All.UnitLength_in_cm) / (All.UnitTime_in_s / All.HubbleParam));
      localSphP[i].BConserved[0] = Bx / sqrt(4. * M_PI) * localSphP[i].Volume /
                                   (sqrt(All.UnitMass_in_g / All.UnitLength_in_cm) / (All.UnitTime_in_s / All.HubbleParam));
      localSphP[i].BConserved[1] = By / sqrt(4. * M_PI) * localSphP[i].Volume /
                                   (sqrt(All.UnitMass_in_g / All.UnitLength_in_cm) / (All.UnitTime_in_s / All.HubbleParam));
      localSphP[i].BConserved[2] = Bz / sqrt(4. * M_PI) * localSphP[i].Volume /
                                   (sqrt(All.UnitMass_in_g / All.UnitLength_in_cm) / (All.UnitTime_in_s / All.HubbleParam));
#else

      localSphP[i].B[0] = All.InjectionBx_InGauss / sqrt(4. * M_PI) /
                          (sqrt(All.UnitMass_in_g / All.UnitLength_in_cm) / (All.UnitTime_in_s / All.HubbleParam));
      localSphP[i].B[1] = All.InjectionBy_InGauss / sqrt(4. * M_PI) /
                          (sqrt(All.UnitMass_in_g / All.UnitLength_in_cm) / (All.UnitTime_in_s / All.HubbleParam));
      localSphP[i].B[2] = All.InjectionBz_InGauss / sqrt(4. * M_PI) /
                          (sqrt(All.UnitMass_in_g / All.UnitLength_in_cm) / (All.UnitTime_in_s / All.HubbleParam));
      localSphP[i].BConserved[0] = All.InjectionBx_InGauss / sqrt(4. * M_PI) * localSphP[i].Volume /
                                   (sqrt(All.UnitMass_in_g / All.UnitLength_in_cm) / (All.UnitTime_in_s / All.HubbleParam));
      localSphP[i].BConserved[1] = All.InjectionBy_InGauss / sqrt(4. * M_PI) * localSphP[i].Volume /
                                   (sqrt(All.UnitMass_in_g / All.UnitLength_in_cm) / (All.UnitTime_in_s / All.HubbleParam));
      localSphP[i].BConserved[2] = All.InjectionBz_InGauss / sqrt(4. * M_PI) * localSphP[i].Volume /
                                   (sqrt(All.UnitMass_in_g / All.UnitLength_in_cm) / (All.UnitTime_in_s / All.HubbleParam));
#endif /* WINDTUNNEL_READ_IN_BFIELD */

      localSphP[i].Energy +=
          0.5 *
          (localSphP[i].B[0] * localSphP[i].B[0] + localSphP[i].B[1] * localSphP[i].B[1] + localSphP[i].B[2] * localSphP[i].B[2]) *
          localSphP[i].Volume * All.cf_atime;

#endif /* MHD */

#if defined(GFM_NORMALIZED_METAL_ADVECTION) && (defined(GFM_SET_METALLICITY) || defined(GFM_PREENRICH))

#ifndef GFM_SET_METALLICITY
      localSphP[i].Metallicity = All.metallicity;
#else
      localSphP[i].Metallicity = All.GasMetallicityInSolar * GFM_SOLAR_METALLICITY;
#endif
      localSphP[i].MassMetallicity = SphP[i].Metallicity * localP[i].Mass;
      for(int j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
        {
          localSphP[i].MetalsFraction[j] = All.mass_fractions[j];
          localSphP[i].MassMetals[j]     = localSphP[i].MetalsFraction[j] * localP[i].Mass;
        }
#endif
    }
#endif /* WINDTUNNEL_FIXVARIABLESININJECTIONREGION */

    /*Check if we have buffer (inflow-outflow)cells with prescribed values
     */
#if defined(BOUNDARY_INFLOWOUTFLOW_MINID) && defined(BOUNDARY_INFLOWOUTFLOW_MAXID) && \
    !defined(WINDTUNNEL_FIXVARIABLESININJECTIONREGION)
  if(localP[i].ID >= BOUNDARY_INFLOWOUTFLOW_MINID && localP[i].ID < BOUNDARY_INFLOWOUTFLOW_MAXID)
    {
#ifdef WINDTUNNEL
      localSphP[i].Density            = All.InjectionDensity;
      localP[i].Vel[0]                = 0;
      localP[i].Vel[1]                = 0;
      localP[i].Vel[2]                = 0;
      localP[i].Vel[WINDTUNNEL_COORD] = All.InjectionVelocity;
      localSphP[i].Utherm             = All.InjectionUtherm;

#ifdef GFM_NORMALIZED_METAL_ADVECTION
      localSphP[i].Metallicity     = All.metallicity;
      localSphP[i].MassMetallicity = SphP[i].Metallicity * localP[i].Mass;
      for(int j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
        {
          localSphP[i].MetalsFraction[j] = All.mass_fractions[j];
          localSphP[i].MassMetals[j]     = localSphP[i].MetalsFraction[j] * localP[i].Mass;
        }
#endif

#endif
    }
#endif

#ifdef SPECIAL_BOUNDARY
  /*If the special boundaries are at the same time buffer cells */
  if((localP[i].ID == -1) || (localP[i].ID == -2))
    {
      dt = (localP[i].TimeBinHydro ? (((integertime)1) << localP[i].TimeBinHydro) : 0) * All.Timebase_interval;
      boundary_get_velocity(localP[i].Pos[0], localP[i].Pos[1], localP[i].Pos[2], &localP[i].Vel[0], &localP[i].Vel[1],
                            &localP[i].Vel[2]);
    }
#endif

#if(defined(WINDTUNNEL) || defined(SPECIAL_BOUNDARY)) && !defined(WINDTUNNEL_FIXVARIABLESININJECTIONREGION)
  if(localP[i].ID >= BOUNDARY_INFLOWOUTFLOW_MINID && localP[i].ID < BOUNDARY_INFLOWOUTFLOW_MAXID)
    {
      localP[i].Mass           = localSphP[i].Density * localSphP[i].Volume;
      localSphP[i].Momentum[0] = localP[i].Vel[0] * localP[i].Mass;
      localSphP[i].Momentum[1] = localP[i].Vel[1] * localP[i].Mass;
      localSphP[i].Momentum[2] = localP[i].Vel[2] * localP[i].Mass;
#ifndef ISOTHERM_EQS
      localSphP[i].Energy =
          pvd->atime * pvd->atime * localP[i].Mass * localSphP[i].Utherm +
          0.5 * localP[i].Mass *
              (localP[i].Vel[0] * localP[i].Vel[0] + localP[i].Vel[1] * localP[i].Vel[1] + localP[i].Vel[2] * localP[i].Vel[2]);
      localSphP[i].Pressure = GAMMA_MINUS1 * localSphP[i].Density * localSphP[i].Utherm;

#endif
    }
#endif
}

/*! \brief Updates primitive variables in a specified cell.
 *
 *  \param[in] localP Pointer to particle data array.
 *  \param[in,out] localSphP Pointer to cell data array.
 *  \param[in] i Index of cell in localP and localSphP arrays.
 *  \param[in] pvd additional data that is needed for update (e.g. cosmological
 *             factors).
 *
 *  \return void
 */
void update_primitive_variables_single(struct particle_data *localP, struct sph_particle_data *localSphP, int i,
                                       struct pv_update_data *pvd)
{
#ifdef DG

  localSphP[i].Density = localP[i].Mass / localSphP[i].Volume;

#ifndef FIX_MEAN_VALUES
  if(localSphP[i].Density < Epsilon_rho)
    {
      printf("cell: %d, density: %g\n", i, localSphP[i].Density);
      terminate("Negative density detected!\n");
    }
#endif

  localP[i].Vel[0] = localSphP[i].Momentum[0] / localP[i].Mass;
  localP[i].Vel[1] = localSphP[i].Momentum[1] / localP[i].Mass;
  localP[i].Vel[2] = localSphP[i].Momentum[2] / localP[i].Mass;

  return;
#endif

#ifdef MHD_POWELL_ENERGYLIMITER
  if(localSphP[i].Powell_Energy / localP[i].Mass / (pvd->atime * pvd->atime) > -0.3 * localSphP[i].Utherm)
    localSphP[i].Energy += localSphP[i].Powell_Energy;

  localSphP[i].Powell_Energy = 0;
#endif

  localSphP[i].Density = localP[i].Mass / localSphP[i].Volume;

#ifdef RT_ADVECT
  int j;
  for(j = 0; j < RT_N_DIR; j++)
    localSphP[i].DensPhot[j] = localSphP[i].Photons[j] / localSphP[i].Volume;
#endif

#ifdef MRT
#ifdef MRT_COMOVING
  cell_do_lorentz_boost(i, localSphP, P[i].Vel[0] - localSphP[i].Old_Vel[0], P[i].Vel[1] - localSphP[i].Old_Vel[1],
                        P[i].Vel[2] - localSphP[i].Old_Vel[2]);
#endif

#endif

  if(localP[i].Mass > 0)
    {
#ifdef MRT_UPDATE_AT_END_OF_STEP
      double KE_old = 0.5 *
                      (localSphP[i].Momentum[0] * localSphP[i].Momentum[0] + localSphP[i].Momentum[1] * localSphP[i].Momentum[1] +
                       localSphP[i].Momentum[2] * localSphP[i].Momentum[2]) /
                      localP[i].Mass;

#if defined(MRT_COOLING_HEATING) || defined(MRT_CHEM_SG)

      double unow = (localSphP[i].Energy - KE_old) / localP[i].Mass;

      unow += localSphP[i].RT_dutherm;  //* localSphP[i].RT_mass/localP[i].Mass ;
      if(unow < All.MinEgySpec)
        unow = All.MinEgySpec;

      localSphP[i].Energy = unow * localP[i].Mass + KE_old;

#endif

#ifdef MRT_RADIATION_PRESSURE
      localSphP[i].Momentum[0] += localSphP[i].RT_mominj[0];
      localSphP[i].Momentum[1] += localSphP[i].RT_mominj[1];
      localSphP[i].Momentum[2] += localSphP[i].RT_mominj[2];

      double KE_new = 0.5 *
                      (localSphP[i].Momentum[0] * localSphP[i].Momentum[0] + localSphP[i].Momentum[1] * localSphP[i].Momentum[1] +
                       localSphP[i].Momentum[2] * localSphP[i].Momentum[2]) /
                      localP[i].Mass;

      localSphP[i].Energy += KE_new - KE_old;
#endif

#endif

      localP[i].Vel[0] = localSphP[i].Momentum[0] / localP[i].Mass;
      localP[i].Vel[1] = localSphP[i].Momentum[1] / localP[i].Mass;
      localP[i].Vel[2] = localSphP[i].Momentum[2] / localP[i].Mass;

#ifdef COSMIC_RAYS
      if(SphP[i].CR_Energy / SphP[i].Volume < All.MinimumCREnergyDensity)
        {
          SphP[i].CR_Energy = All.MinimumCREnergyDensity * SphP[i].Volume;
        }
#endif

#ifdef SGS_TURBULENCE
      if(localSphP[i].SgsTData.Energy < All.SgsTConst.MinimumSgsTSpecificEnergy * localP[i].Mass)
        {
          /*printf("SGS_TURBULENCE: artificially injected sgs turbulence energy\n");*/
          localSphP[i].SgsTData.Energy = All.SgsTConst.MinimumSgsTSpecificEnergy * localP[i].Mass;
        }
#endif

#ifdef MAXSCALARS
      for(int k = 0; k < N_Scalar; k++)
        {
#ifdef MHD_THERMAL_ENERGY_SWITCH
          if(k != ScalarIndex.Etherm) /* don't always update thermal energy */
#endif
            *(MyFloat *)(((char *)(&localSphP[i])) + scalar_elements[k].offset) =
                *(MyFloat *)(((char *)(&localSphP[i])) + scalar_elements[k].offset_mass) / localP[i].Mass;
        }
#endif

#ifdef GFM_NORMALIZED_METAL_ADVECTION
      double mass_metals = 0;

      for(int j = 2; j < GFM_N_CHEM_ELEMENTS; j++)
        mass_metals += localSphP[i].MassMetals[j];

      localSphP[i].MassMetallicity = mass_metals;
      localSphP[i].Metallicity     = mass_metals / localP[i].Mass;
#endif

#ifdef GFM_CHEMTAGS
      for(int k = 0; k < GFM_N_CHEM_TAGS; k++)
        localSphP[i].MassMetalsChemTagsFraction[k] = localSphP[i].MassMetalsChemTags[k] / localP[i].Mass;
#endif  // en

#ifdef ACTIVE_CELL_SPIN
      {
        /* Spin = I * Omega -> solve for Omega */
        int j, k;
        double momInertia[9];

        /* we have to move the tensor from the center of mass to the center of rotation of the cell */
        {
          double dx = NEAREST_X(localSphP[i].CenterOffset[0]);
          double dy = NEAREST_Y(localSphP[i].CenterOffset[1]);
          double dz = NEAREST_Z(localSphP[i].CenterOffset[2]);

          momInertia[0] = localSphP[i].MomentOfInertia[0][0] + localSphP[i].Volume * (dy * dy + dz * dz);
          momInertia[4] = localSphP[i].MomentOfInertia[1][1] + localSphP[i].Volume * (dx * dx + dz * dz);
          momInertia[8] = localSphP[i].MomentOfInertia[2][2] + localSphP[i].Volume * (dx * dx + dy * dy);
          momInertia[1] = localSphP[i].MomentOfInertia[0][1] - localSphP[i].Volume * dx * dy;
          momInertia[3] = localSphP[i].MomentOfInertia[1][0] - localSphP[i].Volume * dx * dy;
          momInertia[2] = localSphP[i].MomentOfInertia[0][2] - localSphP[i].Volume * dx * dz;
          momInertia[6] = localSphP[i].MomentOfInertia[2][0] - localSphP[i].Volume * dx * dz;
          momInertia[5] = localSphP[i].MomentOfInertia[1][2] - localSphP[i].Volume * dy * dz;
          momInertia[7] = localSphP[i].MomentOfInertia[2][1] - localSphP[i].Volume * dy * dz;
        }

        for(j = 0; j < 3; j++)
          for(k = 0; k < 3; k++)
            momInertia[j * 3 + k] *= localSphP[i].Density;

        gsl_matrix_view A = gsl_matrix_view_array(momInertia, 3, 3);
        gsl_vector_view b = gsl_vector_view_array(SphP[i].Spin, 3);
        gsl_vector *x     = gsl_vector_alloc(3);

        int s;
        gsl_permutation *p = gsl_permutation_alloc(3);
        gsl_linalg_LU_decomp(&A.matrix, p, &s);
        gsl_linalg_LU_solve(&A.matrix, p, &b.vector, x);

        for(j = 0; j < 3; j++)
          localSphP[i].Omega[j] = gsl_vector_get(x, j);

        gsl_permutation_free(p);
        gsl_vector_free(x);
      }
#endif

#if defined(GFM_STELLAR_EVOLUTION) && (GFM_STELLAR_EVOLUTION == 1)
      /* re-normalizing the primitive metallicity variables, because in this mode
         the gas MassMetallicity & MassMetals are updated (from stellar evolution) but not gas Mass */
      double mass_by_metals = 0;
      for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
        {
          mass_by_metals += localSphP[i].MassMetals[k];
        }
      localSphP[i].Metallicity *= (localP[i].Mass / mass_by_metals);
      for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
        {
          localSphP[i].MetalsFraction[k] *= (localP[i].Mass / mass_by_metals);
        }
#ifdef GFM_DUST
      /* TODO? */
      for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
        {
          localSphP[i].MetalsDustFraction[k] *= (localP[i].Mass / mass_by_metals);
        }
#endif

#endif

#ifdef MHD
      localSphP[i].B[0] = localSphP[i].BConserved[0] / localSphP[i].Volume;
      localSphP[i].B[1] = localSphP[i].BConserved[1] / localSphP[i].Volume;
      localSphP[i].B[2] = localSphP[i].BConserved[2] / localSphP[i].Volume;
#ifdef MHD_DEDNER
#ifdef MHD_DEDNER_VARIABLE_SPEED
      localSphP[i].Psi = localSphP[i].PsiConserved / localP[i].Mass;
#else
      localSphP[i].Psi         = localSphP[i].PsiConserved / localSphP[i].Volume;
#endif
#endif
#endif

#ifdef MHD_CT
      localSphP[i].A[0] = localSphP[i].AConserved[0] / localSphP[i].Volume;
      localSphP[i].A[1] = localSphP[i].AConserved[1] / localSphP[i].Volume;
      localSphP[i].A[2] = localSphP[i].AConserved[2] / localSphP[i].Volume;
#endif

#ifdef TRACER_FIELD
      localSphP[i].Tracer = localSphP[i].ConservedTracer / localP[i].Mass;
#endif
    }
  else /* P[i].Mass <= 0 */
    {
      localP[i].Vel[0] = 0;
      localP[i].Vel[1] = 0;
      localP[i].Vel[2] = 0;

#ifdef MAXSCALARS
      for(int k = 0; k < N_Scalar; k++)
        *(MyFloat *)(((char *)(&localSphP[i])) + scalar_elements[k].offset) = 0;
#endif
#ifdef TRACER_FIELD
      localSphP[i].Tracer = 0;
#endif
    }
}

/*! \brief Updates the internal energy field in a specified cell
 *
 *  \param[in] localP Pointer to particle data array
 *  \param[in,out] localSphP Pointer to cell data array
 *  \param[in] i Index of cell in localP and localSphP arrays
 *  \param[in] pvd additional data that is needed for update (e.g. cosmological
 *             factors)
 *
 *  \return void
 */
void update_internal_energy(struct particle_data *localP, struct sph_particle_data *localSphP, int i, struct pv_update_data *pvd)
{
#ifdef DG
  // rho>=Epsilon_rho
  localSphP[i].Utherm =
      (localSphP[i].Energy / localP[i].Mass -
       0.5 * (localP[i].Vel[0] * localP[i].Vel[0] + localP[i].Vel[1] * localP[i].Vel[1] + localP[i].Vel[2] * localP[i].Vel[2])) /
      (pvd->atime * pvd->atime);
  return;
#endif

#ifndef ISOTHERM_EQS
  double ulimit;
#if(defined(MEASURE_DISSIPATION_RATE) && defined(USE_ENTROPY_FOR_COLD_FLOWS)) || defined(TGCHEM)
  double dt = (localP[i].TimeBinHydro ? (((integertime)1) << localP[i].TimeBinHydro) : 0) * All.Timebase_interval / pvd->hubble_a;
#endif

#ifdef TGCHEM
  double uold = localSphP[i].Utherm;
#endif

  if(localP[i].Mass > 0)
    {
#ifdef MESHRELAX
      localSphP[i].Utherm = localSphP[i].Energy / localP[i].Mass;
#else
      localSphP[i].Utherm =
          (localSphP[i].Energy / localP[i].Mass -
           0.5 * (localP[i].Vel[0] * localP[i].Vel[0] + localP[i].Vel[1] * localP[i].Vel[1] + localP[i].Vel[2] * localP[i].Vel[2])) /
          (pvd->atime * pvd->atime);

#ifdef MRT_UPDATE_AT_END_OF_STEP
#if defined(MRT_COOLING_HEATING) || defined(MRT_CHEM_SG)
      localSphP[i].RT_utherm    = localSphP[i].Utherm;
      localSphP[i].RT_dutherm   = 0.0;
      localSphP[i].RT_mass      = localP[i].Mass;
#endif
#ifdef MRT_RADIATION_PRESSURE
      localSphP[i].RT_mominj[0] = localSphP[i].RT_mominj[1] = localSphP[i].RT_mominj[2] = 0.0;
#endif
#endif
#endif

#ifdef FLD
      MyFloat mu               = 2.33;
      MyFloat u                = SphP[i].Utherm * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
      MyFloat temp             = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;
      localSphP[i].Temperature = temp;
#endif

#if defined(VS_TURB) || defined(AB_TURB)
      localSphP[i].Utherm -= get_turb_pot(localP[i].Pos[0], localP[i].Pos[1], localP[i].Pos[2]);
#endif

#ifdef MHD
      localSphP[i].Utherm -=
          0.5 *
          (localSphP[i].B[0] * localSphP[i].B[0] + localSphP[i].B[1] * localSphP[i].B[1] + localSphP[i].B[2] * localSphP[i].B[2]) /
          localSphP[i].Density / pvd->atime;
#endif

#ifdef MHD_THERMAL_ENERGY_SWITCH
      if(localSphP[i].Utherm < 0.)
        {
          EgyInjection -= localSphP[i].Energy;
          localSphP[i].Energy -= localSphP[i].Utherm * localP[i].Mass * pvd->atime * pvd->atime;
          localSphP[i].Utherm = localSphP[i].Etherm / localP[i].Mass / (pvd->atime * pvd->atime);
          localSphP[i].Energy += localSphP[i].Utherm * localP[i].Mass * pvd->atime * pvd->atime;
          EgyInjection += localSphP[i].Energy;
        }
#endif

#ifdef ACTIVE_CELL_SPIN
      {
        int j, k;
        double Erot = 0;
        for(j = 0; j < 3; j++)
          for(k = 0; k < 3; k++)
            Erot += localSphP[i].Omega[j] * localSphP[i].Omega[k] * localSphP[i].MomentOfInertia[j][k] / localSphP[i].Volume;
        Erot *= 0.5;

        /* still undecided on including this term -> more testing required */
        // localSphP[i].Utherm -= Erot;
      }
#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
      double Anew;
      int flag = 1;

#ifdef ENTROPY_MACH_THRESHOLD
      if(localSphP[i].MaxMach > ENTROPY_MACH_THRESHOLD)
        flag = 0;
#else
      terminate("not implemented any more");
#endif

      localSphP[i].A = exp(localSphP[i].Entropy / localP[i].Mass);

#ifdef TGCHEM
      Anew = (localSphP[i].Gamma - 1) * localSphP[i].Utherm / pow(localSphP[i].Density * pvd->a3inv, GAMMA_MINUS1);
#else
#if defined(SGCHEM) && defined(VARIABLE_GAMMA)
      Anew = (localSphP[i].GammaE - 1) * localSphP[i].Utherm / pow(localSphP[i].Density * pvd->a3inv, GAMMA_MINUS1);
#else
      Anew = GAMMA_MINUS1 * localSphP[i].Utherm / pow(localSphP[i].Density * pvd->a3inv, GAMMA_MINUS1);
#endif
#endif

#ifndef MEASURE_DISSIPATION_RATE
      if(flag != 0 || Anew < 0.25 * localSphP[i].A) /* here we keep the entropy, and reinitialize the thermal energy */
        {
#ifdef TGCHEM
          localSphP[i].Utherm = localSphP[i].A * pow(localSphP[i].Density * pvd->a3inv, GAMMA_MINUS1) / (localSphP[i].Gamma - 1);
#else
#if defined(SGCHEM) && defined(VARIABLE_GAMMA)
          localSphP[i].Utherm = localSphP[i].A * pow(localSphP[i].Density * pvd->a3inv, GAMMA_MINUS1) / (localSphP[i].GammaE - 1);
#else
          localSphP[i].Utherm = localSphP[i].A * pow(localSphP[i].Density * pvd->a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;
#endif
#endif
          localSphP[i].Energy =
              localP[i].Mass * pvd->atime * pvd->atime * localSphP[i].Utherm +
              0.5 * localP[i].Mass *
                  (localP[i].Vel[0] * localP[i].Vel[0] + localP[i].Vel[1] * localP[i].Vel[1] + localP[i].Vel[2] * localP[i].Vel[2]);
#ifdef MHD
          localSphP[i].Energy +=
              0.5 *
              (localSphP[i].B[0] * localSphP[i].B[0] + localSphP[i].B[1] * localSphP[i].B[1] + localSphP[i].B[2] * localSphP[i].B[2]) *
              localSphP[i].Volume * pvd->atime;
#endif

          pvd->count_keep_entropy++;
        }
      else /* here the entropy is reinitialized, the thermal energy is kept */
#endif     /* end of not MEASURE_DISSIPATION_RATE */
        {
#ifdef MEASURE_DISSIPATION_RATE
          if(dt)
            localSphP[i].DuDt =
                (localSphP[i].Utherm - localSphP[i].A * pow(localSphP[i].Density * pvd->a3inv, GAMMA_MINUS1) / GAMMA_MINUS1) / dt;
          else
            localSphP[i].DuDt = 0;
#endif
          localSphP[i].A = Anew;

          localSphP[i].Entropy = log(localSphP[i].A) * localP[i].Mass;

          pvd->count_update_entropy++;

          if(!gsl_finite(localSphP[i].A))
            {
              print_particle_info(i);
              terminate("stop");
            }
        }
#endif /* end of USE_ENTROPY_FOR_COLD_FLOWS */

      ulimit = All.MinEgySpec;
#ifdef JEANS_UTHERM_LIMIT
      // double celldim = 2.0 * get_cell_radius(i);
      // double jeans_length = JEANS_UTHERM_LIMIT*celldim;
      // double ulimit_jeans = All.G * localSphP[i].Density*jeans_length*jeans_length/M_PI/GAMMA/GAMMA_MINUS1;
      // ulimit = fmax(ulimit,ulimit_jeans);
      double max_rho = 4.0 * All.ReferenceGasPartMass / All.MinVolume;
      ulimit *= pow(localSphP[i].Density / max_rho, 2.0 / 3.0);
#endif

      if(localSphP[i].Utherm < ulimit)
        {
          EgyInjection -= localSphP[i].Energy;

          localSphP[i].Utherm = ulimit;

#ifdef MESHRELAX
          localSphP[i].Energy = localP[i].Mass * localSphP[i].Utherm;
#else
          localSphP[i].Energy =
              pvd->atime * pvd->atime * localP[i].Mass * localSphP[i].Utherm +
              0.5 * localP[i].Mass *
                  (localP[i].Vel[0] * localP[i].Vel[0] + localP[i].Vel[1] * localP[i].Vel[1] + localP[i].Vel[2] * localP[i].Vel[2]);
#endif

#ifdef MHD
          localSphP[i].Energy +=
              0.5 *
              (localSphP[i].B[0] * localSphP[i].B[0] + localSphP[i].B[1] * localSphP[i].B[1] + localSphP[i].B[2] * localSphP[i].B[2]) *
              localSphP[i].Volume * pvd->atime;
#endif

#ifdef ACTIVE_CELL_SPIN
          {
            int j, k;
            double Erot = 0;
            for(j = 0; j < 3; j++)
              for(k = 0; k < 3; k++)
                Erot += localSphP[i].Omega[j] * localSphP[i].Omega[k] * localSphP[i].MomentOfInertia[j][k] * localSphP[i].Density;
            Erot *= 0.5;

            /* still undecided on including this term -> more testing required */
            // localSphP[i].Energy += Erot;
          }
#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS

#ifdef TGCHEM
          localSphP[i].A = (localSphP[i].Gamma - 1) * localSphP[i].Utherm / pow(localSphP[i].Density * pvd->a3inv, GAMMA_MINUS1);
#else
#if defined(SGCHEM) && defined(VARIABLE_GAMMA)
          localSphP[i].A = (localSphP[i].GammaE - 1) * localSphP[i].Utherm / pow(localSphP[i].Density * pvd->a3inv, GAMMA_MINUS1);
#else
          localSphP[i].A = GAMMA_MINUS1 * localSphP[i].Utherm / pow(localSphP[i].Density * pvd->a3inv, GAMMA_MINUS1);
#endif
#endif
          localSphP[i].Entropy = log(localSphP[i].A) * localP[i].Mass;

          if(!gsl_finite(localSphP[i].A))
            {
              print_particle_info(i);
              terminate("stop");
            }
#endif
          EgyInjection += localSphP[i].Energy;
        }
    }
  else
    localSphP[i].Utherm = 0;

  if(localSphP[i].Density < All.LimitUBelowThisDensity && localSphP[i].Utherm > All.LimitUBelowCertainDensityToThisValue)
    {
      localSphP[i].Utherm = All.LimitUBelowCertainDensityToThisValue;
      localSphP[i].Energy =
          pvd->atime * pvd->atime * localP[i].Mass * localSphP[i].Utherm +
          0.5 * localP[i].Mass *
              (localP[i].Vel[0] * localP[i].Vel[0] + localP[i].Vel[1] * localP[i].Vel[1] + localP[i].Vel[2] * localP[i].Vel[2]);
#ifdef MHD
      localSphP[i].Energy +=
          0.5 *
          (localSphP[i].B[0] * localSphP[i].B[0] + localSphP[i].B[1] * localSphP[i].B[1] + localSphP[i].B[2] * localSphP[i].B[2]) *
          localSphP[i].Volume * pvd->atime;
#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
#ifdef TGCHEM
      localSphP[i].A = (localSphP[i].Gamma - 1.0) * localSphP[i].Utherm / pow(localSphP[i].Density * pvd->a3inv, GAMMA_MINUS1);
#else
#if defined(SGCHEM) && defined(VARIABLE_GAMMA)
      localSphP[i].A = (localSphP[i].GammaE - 1.0) * localSphP[i].Utherm / pow(localSphP[i].Density * pvd->a3inv, GAMMA_MINUS1);
#else
      localSphP[i].A = GAMMA_MINUS1 * localSphP[i].Utherm / pow(localSphP[i].Density * pvd->a3inv, GAMMA_MINUS1);
#endif
#endif
      localSphP[i].Entropy = log(localSphP[i].A) * localP[i].Mass;
#endif
    }

  if(localSphP[i].Utherm < 0)
    {
      printf("negative utherm %g\n", localSphP[i].Utherm);
      terminate("stop");
    }

#ifdef ENTROPY_MACH_THRESHOLD
#ifdef OUTPUT_MACHNUM
  localSphP[i].MaxMachNumber = localSphP[i].MaxMach;
#endif
  localSphP[i].MaxMach = 0; /* reset */
#endif

#ifdef TGCHEM
  if(dt)
    localSphP[i].HydroHeatRate = (localSphP[i].Utherm - uold) / dt;
  else
    localSphP[i].HydroHeatRate = 0.;
#endif

#endif /* end of not ISOTHERM_EQS */

    // temperature calculation

#ifdef RADCOOL_HOTHALO
#ifdef COOLING
  double meanweight = 4. / (3 * HYDROGEN_MASSFRAC + 1 + 4 * HYDROGEN_MASSFRAC * SphP[i].Ne);
#else
  double meanweight = 0.5882352941176471;  // fully ionized
#endif
  localSphP[i].Temperature =
      localSphP[i].Utherm * All.UnitEnergy_in_cgs / All.UnitMass_in_g * GAMMA_MINUS1 * meanweight * PROTONMASS / BOLTZMANN;
#endif
}

/*! \brief Calculates the sound speed of a specified cell
 *
 *  Depends on equation of state and potential sub-resoluiton physics.
 *
 *  \param[in] p Index of gas cell in P and SphP arrays
 *
 *  \return Sound speed
 */
double get_sound_speed(int p)
{
  double csnd;

#ifdef ISOTHERM_EQS
  csnd = All.IsoSoundSpeed;
#else
#ifdef LOCALLY_ISOTHERM_DISK
  csnd = get_isotherm_disk_sound_speed(p);
#else

  double gamma;
#ifdef TGCHEM
  gamma = SphP[p].Gamma;
#else
#ifdef VARIABLE_GAMMA
  gamma = SphP[p].GammaC;
#else  /* This corresponds to ordinary gas dynamics, not ISOTHERM   */
  gamma = GAMMA;
#endif /* VARIABLE_GAMMA */
#endif /* TGCHEM */

  if(SphP[p].Density > 0)
    csnd = sqrt(gamma * SphP[p].Pressure / SphP[p].Density);
  else
    csnd = 0;

#endif /* LOCALLY_ISOTHERM_DISK */
#endif /* ISOTHERM_EQS */

#ifdef MHD
  /* for MHD, this is an upper bound to the signal velocity
     to do it more precisely, the magnet field in normal direction to the interfaces
     has to be taken into account */
  double Bsqr = SphP[p].B[0] * SphP[p].B[0] + SphP[p].B[1] * SphP[p].B[1] + SphP[p].B[2] * SphP[p].B[2];
  if(All.ComovingIntegrationOn)
    Bsqr /= All.Time;
  csnd = sqrt(csnd * csnd + Bsqr / SphP[p].Density);
#endif

#ifdef COSMIC_RAYS
  csnd = sqrt(csnd * csnd + All.GammaCR * SphP[p].CR_Pressure / SphP[p].Density);
#endif

#ifdef SGS_TURBULENCE_RIEMANN_PRESSURE
  csnd = sqrt(csnd * csnd + All.SgsTConst.GammaSgsT * SphP[p].SgsTData.Pressure / SphP[p].Density);
#endif

#ifdef FLD
  csnd =
      sqrt(csnd * csnd + (4. / 9.) * SphP[p].n_gamma * (1. - exp(-SphP[p].Kappa_R * amr_length[Mesh.DP[p].level])) / SphP[p].Density);
#endif

  return csnd;
}

#ifdef LOCALLY_ISOTHERM_DISK
double get_isotherm_disk_sound_speed(int p)
{
  double dx, dy, r;
  double csnd = 0;
  dx          = P[p].Pos[0] - boxHalf_X;
  dy          = P[p].Pos[1] - boxHalf_Y;

  r = sqrt(dx * dx + dy * dy);

  if(r <= All.InnerRadius)
    csnd = All.AspectRatio * sqrt(All.G * All.CentralMass / All.InnerRadius);
  if(r > All.InnerRadius && r < All.OuterRadius)
    csnd = All.AspectRatio * sqrt(All.G * All.CentralMass / r);
  if(r >= All.InnerRadius)
    csnd = All.AspectRatio * sqrt(All.G * All.CentralMass / All.OuterRadius);

  return csnd;
}

int get_isotherm_disk_flag(int p)
{
  double dx, dy, r;
  dx = P[p].Pos[0] - boxHalf_X;
  dy = P[p].Pos[1] - boxHalf_Y;

  r = sqrt(dx * dx + dy * dy);

  if(r > All.InnerRadius && r < All.OuterRadius)
    return 1;
  else
    return 0;
}
#endif
