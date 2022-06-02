/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/grackle/grackle.c
 * \date        12/2017
 * \author      BW Keller
 * \brief       Cooling using the Grackle library
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - 11.12.2017 From scratch rewrite of this module (orig. M. Smith) for grackle 3.0
 * - 12.03.2019 Added photoelectric heating (M. Smith)
 * - 17.01.2022 Updated implementation for compatibility with 2022 Arepo codebase.
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef COOLING
#ifdef GRACKLE

#define ABUNDANCE_CON_TOLERANCE 0.01
#define ABUNDANCE_CON_EPS 1e-6
#define MAX_ABUNDANCE_CON_ITER 1e4
#define ABUNDANCE_EPS 1e-20

/** \brief Initialize the grackle units and parameters structures.
 *
 * This method initializes the my_grackle_units and gracke_data structs
 * provided by grackle.h.  my_grackle_units is defined in allvars.h, and is a
 * code_units struct that contains the code units and the cosmological scale
 * factor a.  grackle_data is a chemistry_data struct that contains the
 * parameter flags for running grackle.
 */
void initialise_grackle(void)
{
#ifdef GRACKLE_VERBOSE
  grackle_verbose = 1;
#endif
  // Initialize units struct my_grackle_units
  my_grackle_units.comoving_coordinates = All.ComovingIntegrationOn;
  my_grackle_units.density_units        = All.UnitDensity_in_cgs * All.cf_a3inv * All.HubbleParam * All.HubbleParam;
  my_grackle_units.length_units         = All.UnitLength_in_cm * All.cf_atime / All.HubbleParam;
  my_grackle_units.time_units           = All.UnitTime_in_s / All.HubbleParam;
  my_grackle_units.velocity_units       = All.UnitVelocity_in_cm_per_s;
  my_grackle_units.a_units              = 1.0;
  my_grackle_units.a_value              = All.cf_atime;

  // Initialize chemistry grackle_data struct.  Contains solver parameters
  chemistry_data *my_grackle_data = (chemistry_data *)malloc(sizeof(chemistry_data));
  if(set_default_chemistry_parameters(my_grackle_data) == 0)
    terminate("Grackle: Error in set_default_chemistry_parameters.\n");
  // Set parameter values for chemistry.
  grackle_data->use_grackle            = All.GrackleOn;                // chemistry on
  grackle_data->with_radiative_cooling = All.GrackleRadiativeCooling;  // cooling on
#ifdef GRACKLE_TAB
  grackle_data->primordial_chemistry = 0;  // Tabulated mode
#ifdef HYDROGEN_ONLY
  grackle_data->HydrogenFractionByMass = 1.0;
#endif
#elif defined(GRACKLE_D)
  grackle_data->primordial_chemistry = 3;  // H, He, H2, D
#elif defined(GRACKLE_H2)
  grackle_data->primordial_chemistry = 2;  // H, He, H2
#else
  grackle_data->primordial_chemistry = 1;  // H, He
#endif
  grackle_data->metal_cooling         = All.GrackleMetalCooling;  // metal cooling on
  grackle_data->UVbackground          = All.GrackleUVB;           // UV background on
  grackle_data->self_shielding_method = All.GrackleSelfShieldingMethod;
  grackle_data->grackle_data_file     = All.GrackleDataFile;  // data file for the tabulated cooling rates (ie, for metals, UVB)

#ifdef GRACKLE_PHOTOELECTRIC
  grackle_data->photoelectric_heating      = 1;                                    // Apply global PE heating rate
  grackle_data->photoelectric_heating_rate = All.GracklePhotoelectricHeatingRate;  // erg cm^-3 s^-1 (n_H/cm^{-3})^-1
#endif

#ifdef PE_MCS
  grackle_data->use_volumetric_heating_rate = 1;
#endif

#ifdef HII_MCS_LR
  grackle_data->UVbackground_boost_on = 1;
#endif

  if(grackle_data->Gamma != GAMMA)
    terminate("The value of Gamma in AREPO and GRACKLE must be the same (almost certainly 5/3)\n");

  if(initialize_chemistry_data(&my_grackle_units) == 0)
    terminate("Grackle: Error in initialize_chemistry_data.\n");
}

/** \brief Update the grackle units as the scale factor evolves
 *
 * When running in non-comoving units, we need to change the values of the
 * units as a goes to 1.
 */
static inline void update_grackle_comoving_units(void)
{
  set_cosmo_factors_for_current_time();
  my_grackle_units.density_units = All.UnitDensity_in_cgs * All.cf_a3inv * All.HubbleParam * All.HubbleParam;
  my_grackle_units.length_units  = All.UnitLength_in_cm * All.cf_atime / All.HubbleParam;
  my_grackle_units.time_units    = All.UnitTime_in_s / All.HubbleParam;
  my_grackle_units.a_units       = 1.0;
  my_grackle_units.a_value       = All.cf_atime;
}

/** \brief Construct the grackle_field_data that will be fed to the grackle integrator
 *
 * This function does the memory allocation required for the construction of a
 * new grackle_field_data struct, which will hold the species, density, and
 * energy information used by the grackle integrator. The teardown method that
 * will free() all of these allocations is destroy_fields(). The
 * grackle_field_data struct is mostly a big blob of pointers to gr_float
 * arrays.
 *
 * \param field_size the number of particles to pass to the integrator.
 * \return pointer to the grackle_field_data
 */
static grackle_field_data *build_fields(int field_size)
{
  grackle_field_data *my_fields = (grackle_field_data *)malloc(sizeof(grackle_field_data));
  // Set grid dimension and size.
  // grid_start and grid_end are used to ignore ghost zones.
  my_fields->grid_rank      = 3;
  my_fields->grid_dimension = (int *)malloc(3 * sizeof(int));
  my_fields->grid_start     = (int *)malloc(3 * sizeof(int));
  my_fields->grid_end       = (int *)malloc(3 * sizeof(int));
  for(int i = 0; i < 3; i++)
    {
      my_fields->grid_dimension[i] = 1;
      my_fields->grid_start[i]     = 0;
      my_fields->grid_end[i]       = 0;
    }
  my_fields->grid_dimension[0] = field_size;
  my_fields->grid_end[0]       = field_size - 1;
  // Set field arrays.
  my_fields->density         = (gr_float *)mymalloc("grackle_density", field_size * sizeof(gr_float));
  my_fields->internal_energy = (gr_float *)mymalloc("grackle_energy", field_size * sizeof(gr_float));
  my_fields->metal_density   = (gr_float *)mymalloc("grackle_metal_density", field_size * sizeof(gr_float));
  my_fields->x_velocity      = (gr_float *)mymalloc("grackle_dummy_velocity", field_size * sizeof(gr_float));
  my_fields->y_velocity      = my_fields->x_velocity;
  my_fields->z_velocity      = my_fields->x_velocity;
#ifdef GRACKLE_TAB
  my_fields->HI_density    = NULL;
  my_fields->HII_density   = NULL;
  my_fields->HeI_density   = NULL;
  my_fields->HeII_density  = NULL;
  my_fields->HeIII_density = NULL;
  my_fields->e_density     = NULL;
#else
  my_fields->HI_density              = (gr_float *)mymalloc("grackle_HI_density", field_size * sizeof(gr_float));
  my_fields->HII_density             = (gr_float *)mymalloc("grackle_HII_density", field_size * sizeof(gr_float));
  my_fields->HeI_density             = (gr_float *)mymalloc("grackle_HeI_density", field_size * sizeof(gr_float));
  my_fields->HeII_density            = (gr_float *)mymalloc("grackle_HeII_density", field_size * sizeof(gr_float));
  my_fields->HeIII_density           = (gr_float *)mymalloc("grackle_HeIII_density", field_size * sizeof(gr_float));
  my_fields->e_density               = (gr_float *)mymalloc("grackle_e_density", field_size * sizeof(gr_float));
#endif
#ifdef GRACKLE_H2
  my_fields->HM_density   = (gr_float *)mymalloc("grackle_HM_density", field_size * sizeof(gr_float));
  my_fields->H2I_density  = (gr_float *)mymalloc("grackle_H2I_density", field_size * sizeof(gr_float));
  my_fields->H2II_density = (gr_float *)mymalloc("grackle_H2II_density", field_size * sizeof(gr_float));
#else
  my_fields->HM_density              = NULL;
  my_fields->H2I_density             = NULL;
  my_fields->H2II_density            = NULL;
#endif
#ifdef GRACKLE_D
  my_fields->DI_density  = (gr_float *)mymalloc("grackle_DI_density", field_size * sizeof(gr_float));
  my_fields->DII_density = (gr_float *)mymalloc("grackle_DII_density", field_size * sizeof(gr_float));
  my_fields->HDI_density = (gr_float *)mymalloc("grackle_HDI_density", field_size * sizeof(gr_float));
#else
  my_fields->DI_density              = NULL;
  my_fields->DII_density             = NULL;
  my_fields->HDI_density             = NULL;
#endif
#ifdef PE_MCS
  my_fields->volumetric_heating_rate = (gr_float *)mymalloc("volumetric_heating_rate", field_size * sizeof(gr_float));
#endif
  return my_fields;
}

/** \brief Free all of the memory used in a grackle_field_data struct
 *
 * Since the grackle_field_data is filled with numerous allocations, we want to
 * make sure each member of it gets free'd, and since myfree must be ordered,
 * this helper function can be used.
 *
 * \param fields The struct pointer we wish to free the content of.
 */
static inline void destroy_fields(grackle_field_data *fields)
{
  free(fields->grid_dimension);
  free(fields->grid_start);
  free(fields->grid_end);
#ifdef PE_MCS
  myfree(fields->volumetric_heating_rate);
#endif
#ifdef GRACKLE_D
  myfree(fields->HDI_density);
  myfree(fields->DII_density);
  myfree(fields->DI_density);
#endif
#ifdef GRACKLE_H2
  myfree(fields->H2II_density);
  myfree(fields->H2I_density);
  myfree(fields->HM_density);
#endif
#ifndef GRACKLE_TAB
  myfree(fields->e_density);
  myfree(fields->HeIII_density);
  myfree(fields->HeII_density);
  myfree(fields->HeI_density);
  myfree(fields->HII_density);
  myfree(fields->HI_density);
#endif
  myfree(fields->x_velocity);
  myfree(fields->metal_density);
  myfree(fields->internal_energy);
  myfree(fields->density);
  free(fields);
}

/** \brief Set the fields in grackle_field_data using the SphP array.
 *
 * This function uses a pair of indices (one for the field arrays, one for the
 * SphP array) to set the field values in a grackle_field_data struct.  This
 * function may need to be optimized in the future to use a memcpy rather than
 * assigning values one-by-one.
 *
 * \param fields The grackle_field_data pointer that will get values set for it.
 * \param field_idx
 * \param sph_idx  The index in the
 */
static inline void set_field_values(grackle_field_data *fields, int field_idx, int sph_idx)
{
  fields->density[field_idx]         = SphP[sph_idx].Density;
  fields->x_velocity[field_idx]      = 0;
  fields->internal_energy[field_idx] = SphP[sph_idx].Utherm;
#ifdef METALS
  fields->metal_density[field_idx] = SphP[sph_idx].Metallicity * SphP[sph_idx].Density;
#else
  fields->metal_density[field_idx]   = All.GrackleInitialMetallicity * SphP[sph_idx].Density;
#endif
#ifndef GRACKLE_TAB
  fields->HI_density[field_idx]    = SphP[sph_idx].GrackleSpeciesFraction[0] * SphP[sph_idx].Density;
  fields->HII_density[field_idx]   = SphP[sph_idx].GrackleSpeciesFraction[1] * SphP[sph_idx].Density;
  fields->HeI_density[field_idx]   = SphP[sph_idx].GrackleSpeciesFraction[2] * SphP[sph_idx].Density;
  fields->HeII_density[field_idx]  = SphP[sph_idx].GrackleSpeciesFraction[3] * SphP[sph_idx].Density;
  fields->HeIII_density[field_idx] = SphP[sph_idx].GrackleSpeciesFraction[4] * SphP[sph_idx].Density;
  fields->e_density[field_idx]     = SphP[sph_idx].e_frac * SphP[sph_idx].Density * (PROTONMASS / ELECTRONMASS);
#endif
#ifdef GRACKLE_H2
  fields->HM_density[field_idx]   = SphP[sph_idx].GrackleSpeciesFraction[5] * SphP[sph_idx].Density;
  fields->H2I_density[field_idx]  = SphP[sph_idx].GrackleSpeciesFraction[6] * SphP[sph_idx].Density;
  fields->H2II_density[field_idx] = SphP[sph_idx].GrackleSpeciesFraction[7] * SphP[sph_idx].Density;
#endif
#ifdef GRACKLE_D
  fields->DI_density[field_idx]  = SphP[sph_idx].GrackleSpeciesFraction[8] * SphP[sph_idx].Density;
  fields->DII_density[field_idx] = SphP[sph_idx].GrackleSpeciesFraction[9] * SphP[sph_idx].Density;
  fields->HDI_density[field_idx] = SphP[sph_idx].GrackleSpeciesFraction[10] * SphP[sph_idx].Density;
#endif
#ifdef PE_MCS
  gr_float pe_rate                           = (gr_float)calculate_pe_heating_rate(sph_idx);
  fields->volumetric_heating_rate[field_idx] = pe_rate;
#endif
#ifdef HII_MCS_LR
  gr_float uvb_boost                 = (gr_float)calculate_uv_background_boost_factor(sph_idx);
  fields->UV_background_boost_factor = uvb_boost;
#endif
}

/** \brief Set values in an array using fields from a grackle_field_data struct
 *
 * This function is essentially the inverse of set_field_values.  It takes the
 * values out of the grackle field arrays and stores them in a gr_float array
 * passed to this function.  It too may be improved with a memcpy refactoring.
 *
 * \param fields The grackle_field_data pointer we will get values from
 * \param field_idx The index in the grackle_field_data structs for the abundance arrays
 * \param abundance_fields The gr_float array for the species abundance
 * \param electron_field The gr_float pointer for the electron abundance
 */
static inline void get_field_values(grackle_field_data *fields, int field_idx, gr_float *abundance_fields, gr_float *electron_field)
{
  abundance_fields[0] = fields->HI_density[field_idx] / fields->density[field_idx];
  abundance_fields[1] = fields->HII_density[field_idx] / fields->density[field_idx];
  abundance_fields[2] = fields->HeI_density[field_idx] / fields->density[field_idx];
  abundance_fields[3] = fields->HeII_density[field_idx] / fields->density[field_idx];
  abundance_fields[4] = fields->HeIII_density[field_idx] / fields->density[field_idx];
  *electron_field     = fields->e_density[field_idx] / fields->density[field_idx] / (PROTONMASS / ELECTRONMASS);
#ifdef GRACKLE_H2
  abundance_fields[5] = fields->HM_density[field_idx] / fields->density[field_idx];
  abundance_fields[6] = fields->H2I_density[field_idx] / fields->density[field_idx];
  abundance_fields[7] = fields->H2II_density[field_idx] / fields->density[field_idx];
#endif
#ifdef GRACKLE_D
  abundance_fields[8]  = fields->DI_density[field_idx] / fields->density[field_idx];
  abundance_fields[9]  = fields->DII_density[field_idx] / fields->density[field_idx];
  abundance_fields[10] = fields->HDI_density[field_idx] / fields->density[field_idx];
#endif
}

/** \brief Call grackle to cool without star formation.
 *
 * This function calls the cool_active_cells() function to actually apply the
 * cooling, as well as handle the timing for diagnostic purposes.
 */
void cooling_only(void)
{
  CPU_Step[CPU_MISC] += measure_time();
  mpi_printf("GRACKLE: Cooling active cells\n");
  cool_active_cells();
  CPU_Step[CPU_COOLINGSFR] += measure_time();
}

/** \brief Use grackle to self-consistently calculate the pressure of a given cell
 *
 * This function is used by the set_pressure_of_cell_internal() function in
 * update_primitive_variables.c to self-consistently set the pressure of a gas
 * cell after Grackle has been run.
 */
MyFloat get_grackle_pressure(int i)
{
  update_grackle_comoving_units();
  gr_float pressure;
  grackle_field_data *fields = build_fields(1);
  fields->grid_dimension[0]  = 1;
  fields->grid_end[0]        = 0;
  set_field_values(fields, 0, i);
  if(calculate_pressure(&my_grackle_units, fields, &pressure) == 0)
    terminate("GRACKLE: Error in calculate_pressure\n");
  destroy_fields(fields);
  return pressure;
}

/** \brief Use grackle to cool all active gas cells.
 *
 * This function does everything required to set up and call grackle to cool
 * all of the active gas cells, as well as update the abundances of all
 * relevant species.
 */
void cool_active_cells(void)
{
  if(All.GrackleOn == 0)
    return;
  update_grackle_comoving_units();
  double dt, dtime, du;
  grackle_field_data *fields = build_fields(1);
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
      if(P[i].Mass == 0 && P[i].ID == 0)
        continue;  // skip cells that have been swallowed or eliminated
      if(SphP[i].Density < All.LimitUBelowThisDensity)
        continue;
      fields->grid_dimension[0] = 1;
      fields->grid_end[0]       = 0;
      set_field_values(fields, 0, i);
      dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
      if(dt == 0)
        continue;
      dtime = All.cf_atime * dt / All.cf_time_hubble_a;
      if(solve_chemistry(&my_grackle_units, fields, dtime) == 0)
        terminate("GRACKLE: Error in solve_chemistry\n");
      du = fields->internal_energy[0] - SphP[i].Utherm;
      SphP[i].Utherm += du;
      SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;
      assert(SphP[i].Utherm > 0);
      assert(SphP[i].Energy > 0);
#ifndef GRACKLE_TAB
      get_field_values(fields, 0, SphP[i].GrackleSpeciesFraction, &SphP[i].e_frac);

      for(int j = 0; j < GRACKLE_SPECIES_NUMBER; j++)
        SphP[i].GrackleSpeciesMass[j] = P[i].Mass * SphP[i].GrackleSpeciesFraction[j];

      SphP[i].e_mass = P[i].Mass * SphP[i].e_frac;
#endif

#if defined(GRACKLE_TEMPERATURE_FLOOR) || defined(HII_MCS)
#ifdef HII_MCS
      if((SphP[i].StromgrenSourceID > 0) && (SphP[i].StromgrenSourceID != HII_MCS_IGNORE_FLAG))
        {
          du = All.PhotoionizationEgySpec - SphP[i].Utherm;
#ifndef GRACKLE_TAB
          /* Fully ionise the Hydrogen and increase the free electron abundance appropriately */
          SphP[i].e_frac += SphP[i].GrackleSpeciesFraction[0] * ELECTRONMASS / PROTONMASS;
          SphP[i].e_mass = P[i].Mass * SphP[i].e_frac;
          SphP[i].GrackleSpeciesFraction[1] += SphP[i].GrackleSpeciesFraction[0] * (1.0 - ELECTRONMASS / PROTONMASS);
          SphP[i].GrackleSpeciesMass[1]     = P[i].Mass * SphP[i].GrackleSpeciesFraction[1];
          SphP[i].GrackleSpeciesFraction[0] = ABUNDANCE_EPS;  // Set to very small non-zero number to avoid advection errors
          SphP[i].GrackleSpeciesMass[0]     = ABUNDANCE_EPS;
#endif
        }
      else
        du = All.MinEgySpec - SphP[i].Utherm;
#else
      du = All.MinEgySpec - SphP[i].Utherm;
#endif

      if(du > 0)
        {
          SphP[i].Utherm += du;
          SphP[i].Energy += P[i].Mass * du * All.cf_atime * All.cf_atime;
        }

#endif  // defined(GRACKLE_TEMPERATURE_FLOOR) || defined(HII_MCS)
      set_pressure_of_cell(i);
    }
  destroy_fields(fields);
}

#ifndef GRACKLE_TAB
/** \brief Converge the initial species abundances.
 *
 * We need reasonable values for the abundances if we want our simulation to
 * start in equilibrium.  This function iterates the grackle chemistry solver
 * until all the species have converged on equilibrium abundances.
 */
void grackle_converge_abundances(void)
{
  int ngas_con, ngas_uncon, iter, con_local;
  int *con_list;
  gr_float new_species_frac[GRACKLE_SPECIES_NUMBER];
  gr_float new_e_frac;
  update_grackle_comoving_units();
  double dtime = 1e-2 * SEC_PER_MEGAYEAR / All.UnitTime_in_s;
  mpi_printf("GRACKLE: Converging abundances\n");
  grackle_field_data *fields = build_fields(NumGas);

  con_list = (int *)mymalloc("grackle_converged_list", NumGas * sizeof(int));

  ngas_con = 0;
  iter     = 0;
  while(ngas_con < NumGas)
    {
      ngas_uncon = 0;
      for(int i = 0; i < NumGas; i++)
        {
          if(P[i].Type != 0)
            continue;
          if(iter == 0)
            {
              con_list[i] = 0;                                  // On first pass assume all abundances are unconverged
              if(SphP[i].Density < All.LimitUBelowThisDensity)  // Except if we should be ignoring this cell
                {
                  con_list[i] = 1;  // We don't care about this cell
                  continue;
                }
            }
          else if(con_list[i] == 1)
            continue;  // This cell already has converged abundances. Don't bother passing to GRACKLE again
          set_field_values(fields, ngas_uncon, i);
          ngas_uncon++;
        }
      fields->grid_dimension[0] = ngas_uncon;
      fields->grid_end[0]       = ngas_uncon - 1;
      if(solve_chemistry(&my_grackle_units, fields, dtime) == 0)
        terminate("GRACKLE: Error in solve_chemistry\n");
      ngas_uncon = 0;
      for(int i = 0; i < NumGas; i++)
        {
          if(con_list[i] == 1)
            continue;
          get_field_values(fields, ngas_uncon, new_species_frac, &new_e_frac);
          con_local = 1;
          if(fabs(new_e_frac - SphP[i].e_frac) > ABUNDANCE_CON_TOLERANCE * SphP[i].e_frac && (new_e_frac > ABUNDANCE_CON_EPS) &&
             (SphP[i].e_frac > ABUNDANCE_CON_EPS))
            {
              con_local = 0;
            }
          else
            for(int j = 0; j < GRACKLE_SPECIES_NUMBER; j++)
              {
                if(fabs(new_species_frac[j] - SphP[i].GrackleSpeciesFraction[j]) >
                       ABUNDANCE_CON_TOLERANCE * SphP[i].GrackleSpeciesFraction[j] &&
                   (new_species_frac[j] > ABUNDANCE_CON_EPS) && (SphP[i].GrackleSpeciesFraction[j] > ABUNDANCE_CON_EPS))
                  {
                    con_local = 0;
                    break;
                  }
              }
          if(con_local == 1)
            {
              ngas_con++;
              con_list[i] = 1;
            }
          for(int j = 0; j < GRACKLE_SPECIES_NUMBER; j++)
            {
              SphP[i].GrackleSpeciesFraction[j] = new_species_frac[j];
              SphP[i].GrackleSpeciesMass[j]     = P[i].Mass * new_species_frac[j];
            }
          SphP[i].e_frac = new_e_frac;
          SphP[i].e_mass = P[i].Mass * new_e_frac;
          ngas_uncon++;
        }
      if(iter++ > MAX_ABUNDANCE_CON_ITER)
        {
#ifdef GRACKLE_UNCONVERGED_IGNORE
          printf("GRACKLE: Task %d failed to converge all cells. %d remain unconverged. Continuing anyway\n", ThisTask, ngas_uncon);
          fflush(stdout);
          break;
#else
          terminate("GRACKLE: unable to converge abundances\n");
#endif
        }
    }
  myfree(con_list);
  destroy_fields(fields);
}
#endif

#if !defined(GRACKLE_TAB) && !defined(GRACKLE_ABUNDANCE_IN_ICS)
/** \brief Set the values for the grackle species for each gas cell
 *
 * If we are not using the tabulated grackle, and our intial conditions don't
 * contain species, we need to start out with some reasonable guesses for the
 * abundances for each species.
 */
void grackle_initialise_abundances(void)
{
  MyFloat metallicity;
  mpi_printf("GRACKLE: Initialising abundances\n");
#ifndef METALS
  metallicity = All.GrackleInitialMetallicity;
#endif
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
#ifdef METALS
      metallicity = SphP[i].Metallicity;
#endif
      SphP[i].GrackleSpeciesFraction[0] = HYDROGEN_MASSFRAC * (1.0 - metallicity);          // HI
      SphP[i].GrackleSpeciesFraction[1] = ABUNDANCE_EPS;                                    // HII
      SphP[i].GrackleSpeciesFraction[2] = (1.0 - HYDROGEN_MASSFRAC) * (1.0 - metallicity);  // HeI
      SphP[i].GrackleSpeciesFraction[3] = ABUNDANCE_EPS;                                    // HeII
      SphP[i].GrackleSpeciesFraction[4] = ABUNDANCE_EPS;                                    // HeIII
#ifdef GRACKLE_H2
      SphP[i].GrackleSpeciesFraction[5] = ABUNDANCE_EPS;  // HM
      SphP[i].GrackleSpeciesFraction[6] = ABUNDANCE_EPS;  // H2I
      SphP[i].GrackleSpeciesFraction[7] = ABUNDANCE_EPS;  // H2II
#endif
#ifdef GRACKLE_D
      SphP[i].GrackleSpeciesFraction[8]  = ABUNDANCE_EPS;  // DI
      SphP[i].GrackleSpeciesFraction[9]  = ABUNDANCE_EPS;  // DII
      SphP[i].GrackleSpeciesFraction[10] = ABUNDANCE_EPS;  // HDI
#endif
      SphP[i].e_frac = ABUNDANCE_EPS;
    }
}
#endif

/** \brief Get the temperature of a gas cell.
 *
 * Use grackle to calculate the self-consistent temperature of a given cell.
 *
 * \param i The index of the cell
 * \return The temperature of the cell
 */
double get_temp_individual_cell_grackle(int i)
{
  update_grackle_comoving_units();
  gr_float temperature;
  grackle_field_data *fields = build_fields(1);
  fields->grid_dimension[0]  = 1;
  fields->grid_end[0]        = 0;
  set_field_values(fields, 0, i);
  if(calculate_temperature(&my_grackle_units, fields, &temperature) == 0)
    terminate("GRACKLE: Error in calculate_temperature.\n");
  destroy_fields(fields);
  return temperature;
}

/** \brief Get the cooling time of a gas cell.
 *
 * Use grackle to calculate the self-consistent cooling time of a given cell.
 *
 * \param i The index of the cell
 * \return The cooling time of the cell
 */
double get_cooling_time_individual_cell_grackle(int i)
{
  update_grackle_comoving_units();
  gr_float cooling_time;
  grackle_field_data *fields = build_fields(1);
  fields->grid_dimension[0]  = 1;
  fields->grid_end[0]        = 0;
  set_field_values(fields, 0, i);
  if(calculate_cooling_time(&my_grackle_units, fields, &cooling_time) == 0)
    terminate("GRACKLE: Error in calculate_cooling_time.\n");
  destroy_fields(fields);
  return cooling_time;
}

/** \brief Get the appropriate timestep, taking in to account the cooling rate
 *
 * Use grackle to calculate the self-consistent cooling time of a given cell.
 *
 * \param i The index of the cell
 * \return The cooling time of the cell
 */
double grackle_get_timestep(int i) { return fabs(get_cooling_time_individual_cell_grackle(i)); }

#endif  // COOLING
#endif  // GRACKLE
