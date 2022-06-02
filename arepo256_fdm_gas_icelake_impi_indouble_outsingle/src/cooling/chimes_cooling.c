#define _POSIX_C_SOURCE 199309L

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef CHIMES
struct globalVariables ChimesGlobalVars;
double *ChimesDustGArr;
double *ChimesH2DissocJArr;

double *ChimesAbundances;
double *ChimesPhotonDensity;
double *ChimesDustG;
double *ChimesH2dissocJ;
#ifdef CHIMES_ADVECT_ABUNDANCES
double *ChimesIonAdvect;
#endif

#ifdef CHIMES_PTHREADS
#include <time.h>
int ThisTask_node;
int NTask_node;
MPI_Comm node_comm;
int buf_index;
int N_active_tot;
pthread_mutex_t mutexcool, mutex_chimes_malloc;
struct gasVariables *ChimesGasVars_buf;
double *ChimesAbundances_buf;
double *ChimesPhotonDensity_buf;
double *ChimesH2dissocJ_buf;
double *ChimesDustG_buf;
double *ChimesTemperature_buf;
struct All_rate_variables_structure **ChimesAllRatesTh;
struct Reactions_Structure **ChimesAllReactionsRootTh;
struct Reactions_Structure **ChimesNonMolecularReactionsRootTh;
#else
struct All_rate_variables_structure *ChimesAllRates;
struct Reactions_Structure *ChimesAllReactionsRoot;
struct Reactions_Structure *ChimesNonMolecularReactionsRoot;
#endif
#ifdef CHIMES_REDSHIFT_DEPENDENT_UVB
double *Chimes_Redshifts;
int Chimes_N_Redshifts;
struct PhotoIonTables_UVB *ChimesPhotoIonTable;
int Chimes_N_Elements_in_Bens_tables;
#endif

#include "../chimes/proto.h"

/* This function initialises the CHIMES cooling module.
 * It replaces the InitCool() function that is used
 * with the standard cooling module.
 */
void ChimesInitCool(void)
{
  mpi_printf("CHIMES: initialising chemistry and cooling module.\n");
  ChimesGlobalVars.updatePhotonFluxOn     = 0;
  ChimesGlobalVars.InitIonState           = 1;
  ChimesGlobalVars.print_debug_statements = 0;
  sprintf(ChimesGlobalVars.BenTablesPath, "%s/bens_tables/", All.ChimesDataPath);
  sprintf(ChimesGlobalVars.AdditionalRatesTablesPath, "%s/additional_rates.hdf5", All.ChimesDataPath);
  sprintf(ChimesGlobalVars.MolecularTablePath, "%s/molecular_cooling_table.hdf5", All.ChimesDataPath);

  /* Currently, we only support a single UV spectrum.
   * We will add further options later. */
  ChimesGlobalVars.N_spectra = 1;

  /* The following arrays will store the dust_G and H2_dissocJ
   * parameters from the spectrum data files. */
  ChimesDustGArr     = (double *)mymalloc("Chimes_dustGarr", ChimesGlobalVars.N_spectra * sizeof(double));
  ChimesH2DissocJArr = (double *)mymalloc("Chimes_H2dissocArr", ChimesGlobalVars.N_spectra * sizeof(double));

#ifdef CHIMES_PTHREADS
  int i;
  ChimesAllRatesTh = (struct All_rate_variables_structure **)mymalloc("Chimes_allRates",
                                                                      CHIMES_PTHREADS * sizeof(struct All_rate_variables_structure *));
  ChimesAllReactionsRootTh =
      (struct Reactions_Structure **)mymalloc("Chimes_allReactions", CHIMES_PTHREADS * sizeof(struct Reactions_Structure *));
  ChimesNonMolecularReactionsRootTh =
      (struct Reactions_Structure **)mymalloc("Chimes_nonMol", CHIMES_PTHREADS * sizeof(struct Reactions_Structure *));

  if(NTask < 11)
    init_chimes(&ChimesGlobalVars, &ChimesAllRatesTh[0], &ChimesAllReactionsRootTh[0], &ChimesNonMolecularReactionsRootTh[0],
                ChimesDustGArr, ChimesH2DissocJArr);
  else
    init_chimes_parallel(&ChimesGlobalVars, &ChimesAllRatesTh[0], &ChimesAllReactionsRootTh[0], &ChimesNonMolecularReactionsRootTh[0],
                         ChimesDustGArr, ChimesH2DissocJArr);
  for(i = 1; i < CHIMES_PTHREADS; i++)
    init_chimes_pthreads(&ChimesGlobalVars, &ChimesAllRatesTh[i], &ChimesAllReactionsRootTh[i], &ChimesNonMolecularReactionsRootTh[i]);
#else
  /* If we have more that 11 tasks, then we can use the parallel
   * version of the init_chimes routine, where the rates tables
   * for each element are read in by a different task. */
  if(NTask < 11)
    init_chimes(&ChimesGlobalVars, &ChimesAllRates, &ChimesAllReactionsRoot, &ChimesNonMolecularReactionsRoot, ChimesDustGArr,
                ChimesH2DissocJArr);
  else
    init_chimes_parallel(&ChimesGlobalVars, &ChimesAllRates, &ChimesAllReactionsRoot, &ChimesNonMolecularReactionsRoot, ChimesDustGArr,
                         ChimesH2DissocJArr);
#endif /* CHIMES_PTHREADS */

#ifdef CHIMES_REDSHIFT_DEPENDENT_UVB
  char RedshiftFilename[500], set_name[500];
  hid_t file_id, dataset;
  herr_t status;
  sprintf(RedshiftFilename, "%s/HM12_cross_sections/redshifts.hdf5", All.ChimesDataPath);
  if(ThisTask == 0)
    {
      file_id = H5Fopen(RedshiftFilename, H5F_ACC_RDONLY, H5P_DEFAULT);
      sprintf(set_name, "/N_Redshifts");
      dataset = H5Dopen(file_id, set_name);
      status  = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Chimes_N_Redshifts);
      status  = H5Dclose(dataset);
    }

  MPI_Bcast(&Chimes_N_Redshifts, 1, MPI_INT, 0, MPI_COMM_WORLD);

  Chimes_Redshifts = (double *)mymalloc("Chimes_redshifts", Chimes_N_Redshifts * sizeof(double));

  if(ThisTask == 0)
    {
      sprintf(set_name, "/Redshifts");
      dataset = H5Dopen(file_id, set_name);
      status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Chimes_Redshifts);
      status  = H5Fclose(dataset);
      status  = H5Fclose(file_id);
    }

  MPI_Bcast(Chimes_Redshifts, Chimes_N_Redshifts, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  chimes_allocate_memory_to_photoion_tables(&ChimesGlobalVars, &ChimesPhotoIonTable, Chimes_N_Elements_in_Bens_tables);

  ChimesLoadPhotoIonTables();
#endif /* CHIMES_REDSHIFT_DEPENDENT_UVB */
}

#ifdef CHIMES_PTHREADS
/* This function loops over active particles and does cooling
 * and chemistry only, using the CHIMES cooling module. This
 * version uses the PTHREADS multithreading option. */
void chimes_cooling_only_pthreads(void)
{
  int idx, i, j;
  int N_active = 0;
  int *N_active_node, *buf_offset;

  mpi_printf("CHIMES: doing chemistry and cooling.\n");

  CPU_Step[CPU_MISC] += measure_time();

#ifdef CHIMES_REDSHIFT_DEPENDENT_UVB
  ChimesLoadPhotoIonTables();
#endif

  chimes_update_all_pointers();

  /* Determine number of active particles and update gas vars. */
  N_active_tot = 0;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i >= 0)
        {
          if(P[i].Mass == 0 && P[i].ID == 0)
            continue; /* skip cells that have been swallowed or eliminated */

          chimes_update_gas_vars(i);
          N_active++;
        }
    }

  N_active_node = (int *)mymalloc_movable(&N_active_node, "N_active", NTask_node * sizeof(int));
  buf_offset    = (int *)mymalloc_movable(&buf_offset, "buf_offset", NTask_node * sizeof(int));

  MPI_Gather(&N_active, 1, MPI_INT, N_active_node, 1, MPI_INT, 0, node_comm);

  if(ThisTask_node == 0)
    {
      for(j = 0; j < NTask_node; j++)
        {
          buf_offset[j] = N_active_tot;
          N_active_tot += N_active_node[j];
        }
    }

  chimes_create_pthreads_buffers(N_active);

  chimes_send_gasvars_to_node_root(N_active, N_active_node, buf_offset);

  chimes_run_threads();

  chimes_send_abundances_to_original_task(N_active, N_active_node, buf_offset);

  chimes_free_pthreads_buffers(N_active);

  myfree_movable(N_active_node);
  myfree_movable(buf_offset);

  CPU_Step[CPU_COOLINGSFR] += measure_time();
}
#else /* CHIMES_PTHREADS */
/* This function loops over active particles and does cooling
 * and chemistry only, using the CHIMES cooling module. */
void chimes_cooling_only(void)
{
  int idx, i;

  mpi_printf("CHIMES: doing chemistry and cooling.\n");

  CPU_Step[CPU_MISC] += measure_time();

#ifdef CHIMES_REDSHIFT_DEPENDENT_UVB
  ChimesLoadPhotoIonTables();
#endif

  chimes_update_all_pointers();

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i >= 0)
        {
          if(P[i].Mass == 0 && P[i].ID == 0)
            continue; /* skip cells that have been swallowed or eliminated */

          chimes_update_gas_vars(i);

          if(SphP[i].ChimesGasVars.hydro_timestep > 0.0)
            chimes_network(&(SphP[i].ChimesGasVars), &ChimesGlobalVars, ChimesAllRates, ChimesAllReactionsRoot,
                           ChimesNonMolecularReactionsRoot);

          chimes_update_particle_energies(i);
        }
    }

  CPU_Step[CPU_COOLINGSFR] += measure_time();
}
#endif /* CHIMES_PTHREADS */

/* This function updates the ChimesGasVars
 * structure for particle i. */
void chimes_update_gas_vars(int i)
{
  double dt      = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
  double dtime   = All.cf_atime * dt / All.cf_time_hubble_a;
  double u_cgs   = SphP[i].Utherm * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
  double rho_cgs = SphP[i].Density * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv;

#ifdef GFM_STELLAR_EVOLUTION
  double H_mass_fraction = SphP[i].MetalsFraction[element_index("Hydrogen")];
#else
  double H_mass_fraction = HYDROGEN_MASSFRAC;
#endif

  SphP[i].ChimesGasVars.temperature = chimes_convert_ucgs_to_temp(u_cgs, i);
  SphP[i].ChimesGasVars.nH_tot      = H_mass_fraction * rho_cgs / PROTONMASS;

  SphP[i].ChimesGasVars.TempFloor   = All.MinGasTemp;
  SphP[i].ChimesGasVars.ThermEvolOn = All.ChimesThermEvolOn;
  SphP[i].ChimesGasVars.ForceEqOn   = All.ChimesForceEqOn;

  /* Extragalactic UV background */
  SphP[i].ChimesGasVars.isotropic_photon_density[0] = All.ChimesIsotropicPhotonDensity;
  SphP[i].ChimesGasVars.dust_G_parameter[0]         = ChimesDustGArr[0];
  SphP[i].ChimesGasVars.H2_dissocJ[0]               = ChimesH2DissocJArr[0];

  SphP[i].ChimesGasVars.cr_rate               = All.ChimesCrRate;
  SphP[i].ChimesGasVars.hydro_timestep        = dtime * All.UnitTime_in_s / All.HubbleParam;
  SphP[i].ChimesGasVars.divVel                = fabs((All.HubbleParam / All.UnitTime_in_s) * SphP[i].DivVel);
  SphP[i].ChimesGasVars.constant_heating_rate = 0.0;

#if defined(CHIMES_JEANS_SHIELDING)
  /* Use the Jeans length to calculate the shielding column
   * density. The normalisation can be varied with the
   * ChimesShieldingLengthFactor parameter. */
  double mol_weight = calculate_mean_molecular_weight(&(SphP[i].ChimesGasVars), &ChimesGlobalVars);
  SphP[i].ChimesGasVars.cell_size =
      min(All.ChimesShieldingLengthFactor *
              sqrt(PI * GAMMA * BOLTZMANNCGS * SphP[i].ChimesGasVars.temperature / (mol_weight * 6.674e-8 * rho_cgs * PROTONMASS)),
          All.ChimesMaxShieldingLength_kpc * 3.086e21); /* cgs */
#elif defined(CHIMES_SOBOLEV_SHIELDING)
  /* Impose a floor to prevent division by zero */
  double grad_rho_mag =
      max(sqrt(pow(SphP[i].Grad.drho[0], 2.0) + pow(SphP[i].Grad.drho[1], 2.0) + pow(SphP[i].Grad.drho[2], 2.0)), 1.0e-40);

  SphP[i].ChimesGasVars.cell_size =
      min(All.ChimesShieldingLengthFactor * (SphP[i].Density / grad_rho_mag) * All.UnitLength_in_cm / All.HubbleParam,
          All.ChimesMaxShieldingLength_kpc * 3.086e21);
#else
  SphP[i].ChimesGasVars.cell_size = 1.0;
#endif

  SphP[i].ChimesGasVars.doppler_broad =
      7.1; /* km/s. For now, just set this constant. Thermal broadening is also added within CHIMES. */

#ifdef GFM_STELLAR_EVOLUTION
  chimes_update_element_abundances(i);
#endif
}

/* Once a particle's chemistry and cooling have been
 * evolved for this time-step, its new temperature
 * in the ChimesGasVars structure is used to update
 * its energies and pressure in SphP. */
void chimes_update_particle_energies(int i)
{
  double unew = chimes_convert_temp_to_ucgs(SphP[i].ChimesGasVars.temperature, i);

  /* Convert from cgs to code units. */
  unew /= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;

  double du = unew - SphP[i].Utherm;

  /* Note that CHIMES already imposes a temperature
   * floor, but we also impose the minimum thermal
   * energy, MinEgySpec, for consistency with the
   * rest of the code. */
  if(unew < All.MinEgySpec)
    du = All.MinEgySpec - SphP[i].Utherm;

  SphP[i].Utherm += du;
  SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;

#ifdef OUTPUT_COOLHEAT
  if(dtime > 0)
    SphP[i].CoolHeat = du * P[i].Mass / dtime;
#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  double dens     = SphP[i].Density;
  SphP[i].A       = (GAMMA - 1.0) * SphP[i].Utherm / pow(dens * All.cf_a3inv, GAMMA - 1);
  SphP[i].Entropy = log(SphP[i].A) * P[i].Mass;
#endif

  set_pressure_of_cell(i);
}

double chimes_convert_ucgs_to_temp(double ucgs, int i)
{
  double mu;
#if defined(USE_SFR) && !defined(SMUGGLE_SFR) && !defined(LOCAL_FEEDBACK) && !defined(SFR_MCS)
  /* Uses the Springel & Hernquist (2003) effective EoS. For this, we use
   * a mean molecular weight assuming fully ionised above the density threshold */

  double eos_dens_threshold = All.PhysDensThresh;
#ifdef MODIFIED_EOS
  eos_dens_threshold *= All.FactorDensThresh;
#endif

  double dens = SphP[i].Density;
  if(dens * All.cf_a3inv >= eos_dens_threshold)
    mu = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC)); /* Assumes fully ionised */
  else
    mu = calculate_mean_molecular_weight(&(SphP[i].ChimesGasVars), &ChimesGlobalVars);
#else
  mu = calculate_mean_molecular_weight(&(SphP[i].ChimesGasVars), &ChimesGlobalVars);
#endif /* USE_SFR && !(SMUGGLE_SFR) && !(LOCAL_FEEDBACK) && !(SFR_MCS) */

  chimes_set_pointers(i);
  return ucgs * GAMMA_MINUS1 * PROTONMASS * mu / BOLTZMANNCGS;
}

double chimes_convert_temp_to_ucgs(double temp, int i)
{
  double mu;
#if defined(USE_SFR) && !defined(SMUGGLE_SFR) && !defined(LOCAL_FEEDBACK) && !defined(SFR_MCS)
  /* Uses the Springel & Hernquist (2003) effective EoS. For this, we use
   * a mean molecular weight assuming fully ionised above the density threshold */

  double eos_dens_threshold = All.PhysDensThresh;
#ifdef MODIFIED_EOS
  eos_dens_threshold *= All.FactorDensThresh;
#endif

  double dens = SphP[i].Density;
  if(dens * All.cf_a3inv >= eos_dens_threshold)
    mu = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC)); /* Assumes fully ionised */
  else
    mu = calculate_mean_molecular_weight(&(SphP[i].ChimesGasVars), &ChimesGlobalVars);
#else
  mu = calculate_mean_molecular_weight(&(SphP[i].ChimesGasVars), &ChimesGlobalVars);
#endif /* USE_SFR && !(SMUGGLE_SFR) && !(LOCAL_FEEDBACK) && !(SFR_MCS) */

  chimes_set_pointers(i);
  return BOLTZMANNCGS * temp / (GAMMA_MINUS1 * mu * PROTONMASS);
}

/* This routine replaces GetCoolingTime() when CHIMES
 * is switched on. If we have heating, we instead
 * return 0.
 */
double ChimesGetCoolingTime(int i, struct All_rate_variables_structure *this_all_rates)
{
  double ucgs, rho_cgs, lambda_cool, tcool;
  double NH_tot      = 0.0;
  double HI_column   = 0.0;
  double H2_column   = 0.0;
  double HeI_column  = 0.0;
  double HeII_column = 0.0;
  double CO_column   = 0.0;
  double H2O_column  = 0.0;
  double OH_column   = 0.0;
  double extinction  = 0.0;

  chimes_set_pointers(i);

  /* If cell self-shielding is on, update
   * column densities */
  if(ChimesGlobalVars.cellSelfShieldingOn == 1)
    {
      NH_tot      = SphP[i].ChimesGasVars.cell_size * SphP[i].ChimesGasVars.nH_tot;
      HI_column   = SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[HI]] * NH_tot;
      H2_column   = SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[H2]] * NH_tot;
      HeI_column  = SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[HeI]] * NH_tot;
      HeII_column = SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[HeII]] * NH_tot;

      if(ChimesGlobalVars.speciesIndices[CO] > -1)
        CO_column = max(SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[CO]], 0.0) * NH_tot;

      if(ChimesGlobalVars.speciesIndices[H2O] > -1)
        H2O_column = max(SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[H2O]], 0.0) * NH_tot;

      if(ChimesGlobalVars.speciesIndices[OH] > -1)
        OH_column = max(SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[OH]], 0.0) * NH_tot;

      extinction = DUSTEFFSIZE * NH_tot * SphP[i].ChimesGasVars.metallicity;
    }

  /* Update the rate coefficients */
  set_constant_rates(&(SphP[i].ChimesGasVars), &ChimesGlobalVars, this_all_rates);
  update_rates(&(SphP[i].ChimesGasVars), &ChimesGlobalVars, HI_column, H2_column, HeI_column, HeII_column, CO_column, extinction,
               this_all_rates);
  update_T_dependent_rates(&(SphP[i].ChimesGasVars), &ChimesGlobalVars, this_all_rates);

  /* Compute cooling rate in erg/s/cm^3 */
  lambda_cool = calculate_total_cooling_rate(&(SphP[i].ChimesGasVars), &ChimesGlobalVars, HI_column, HeI_column, HeII_column,
                                             H2_column, CO_column, H2O_column, OH_column, extinction, this_all_rates);

  /* Compute cooling time */
  ucgs    = chimes_convert_temp_to_ucgs(SphP[i].ChimesGasVars.temperature, i); /* thermal energy per unit mass */
  rho_cgs = SphP[i].Density * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv;
  tcool   = max(rho_cgs * ucgs / lambda_cool, 0.0);

  /* Convert tcool to code units */
  tcool *= All.HubbleParam / All.UnitTime_in_s;

  return tcool;
}

#ifdef GFM_STELLAR_EVOLUTION
/* This routine looks at the change in element abundances
 * compared to the previous timestep and updates the
 * individual ion abundances accordingly. Note that we
 * only include metals in CHIMES when GFM_STELLAR_EVOLUTION
 * is switched on. */
void chimes_update_element_abundances(int i)
{
  double old_abundances[10];
  double delta_zi;
  int j, ion_index;

  double H_mass_fraction = SphP[i].MetalsFraction[element_index("Hydrogen")];

  /* Record the element abundances from the previous time step. */
  for(j = 0; j < 10; j++)
    old_abundances[j] = SphP[i].ChimesGasVars.element_abundances[j];

  /* Update the element abundances in ChimesGasVars. */
  SphP[i].ChimesGasVars.element_abundances[0] = SphP[i].MetalsFraction[element_index("Helium")] / (4.0 * H_mass_fraction);     /* He */
  SphP[i].ChimesGasVars.element_abundances[1] = SphP[i].MetalsFraction[element_index("Carbon")] / (12.0 * H_mass_fraction);    /* C */
  SphP[i].ChimesGasVars.element_abundances[2] = SphP[i].MetalsFraction[element_index("Nitrogen")] / (14.0 * H_mass_fraction);  /* N */
  SphP[i].ChimesGasVars.element_abundances[3] = SphP[i].MetalsFraction[element_index("Oxygen")] / (16.0 * H_mass_fraction);    /* O */
  SphP[i].ChimesGasVars.element_abundances[4] = SphP[i].MetalsFraction[element_index("Neon")] / (20.0 * H_mass_fraction);      /* Ne */
  SphP[i].ChimesGasVars.element_abundances[5] = SphP[i].MetalsFraction[element_index("Magnesium")] / (24.0 * H_mass_fraction); /* Mg */
  SphP[i].ChimesGasVars.element_abundances[6] = SphP[i].MetalsFraction[element_index("Silicon")] / (28.0 * H_mass_fraction);   /* Si */
  SphP[i].ChimesGasVars.element_abundances[9] = SphP[i].MetalsFraction[element_index("Iron")] / (56.0 * H_mass_fraction);      /* Fe */

  /* S and Ca are not explicitly tracked, so assume constant ratios w.r.t. Si */
  SphP[i].ChimesGasVars.element_abundances[7] =
      0.6054160 * SphP[i].MetalsFraction[element_index("Silicon")] / (32.0 * H_mass_fraction); /* S */
  SphP[i].ChimesGasVars.element_abundances[8] =
      0.0941736 * SphP[i].MetalsFraction[element_index("Silicon")] / (40.0 * H_mass_fraction); /* Ca */

  SphP[i].ChimesGasVars.metallicity = SphP[i].Metallicity / 0.0129; /* In Zsol. CHIMES uses Zsol = 0.0129. */

  /* Finally, calculate the change in each element
   * abundance and assume that it is injected into
   * all atomic/ionic/molecular states, preserving the
   * ion and molecule fractions. If this particle
   * previously had no metals, we put the new metals
   * into the lowest ionisation state (i.e. neutral).
   * NOTE: Since the CHIMES abundances are in terms of
   * n_i / n_Htot, we don't have a hydrogen
   * abundance (which would be 1 by definition). */

  /* Helium */
  delta_zi = SphP[i].ChimesGasVars.element_abundances[0] - old_abundances[0];
  for(ion_index = ChimesGlobalVars.speciesIndices[HeI]; ion_index <= ChimesGlobalVars.speciesIndices[HeIII]; ion_index++)
    SphP[i].ChimesGasVars.abundances[ion_index] += (SphP[i].ChimesGasVars.abundances[ion_index] / old_abundances[0]) * delta_zi;

  /* Carbon */
  if(ChimesGlobalVars.element_included[0] == 1)
    {
      delta_zi = SphP[i].ChimesGasVars.element_abundances[1] - old_abundances[1];
      if(old_abundances[1] > 0.0)
        {
          for(ion_index = ChimesGlobalVars.speciesIndices[CI]; ion_index <= ChimesGlobalVars.speciesIndices[CVII]; ion_index++)
            SphP[i].ChimesGasVars.abundances[ion_index] +=
                (SphP[i].ChimesGasVars.abundances[ion_index] / old_abundances[1]) * delta_zi;

          if(SphP[i].ChimesGasVars.temperature < ChimesGlobalVars.T_mol)
            {
              SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[C2]] +=
                  (SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[C2]] / old_abundances[1]) * delta_zi;
              SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[HCOp]] +=
                  (SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[HCOp]] / old_abundances[1]) * delta_zi;
              SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[CH]] +=
                  (SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[CH]] / old_abundances[1]) * delta_zi;
              SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[CH2]] +=
                  (SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[CH2]] / old_abundances[1]) * delta_zi;
              SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[CH3p]] +=
                  (SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[CH3p]] / old_abundances[1]) * delta_zi;
              SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[CO]] +=
                  (SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[CO]] / old_abundances[1]) * delta_zi;
              SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[CHp]] +=
                  (SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[CHp]] / old_abundances[1]) * delta_zi;
              SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[CH2p]] +=
                  (SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[CH2p]] / old_abundances[1]) * delta_zi;
              SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[COp]] +=
                  (SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[COp]] / old_abundances[1]) * delta_zi;
              SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[HOCp]] +=
                  (SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[HOCp]] / old_abundances[1]) * delta_zi;
            }
        }
      else
        SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[CI]] += delta_zi;
    }

  /* Nitrogen */
  if(ChimesGlobalVars.element_included[1] == 1)
    {
      delta_zi = SphP[i].ChimesGasVars.element_abundances[2] - old_abundances[2];
      if(old_abundances[2] > 0.0)
        {
          for(ion_index = ChimesGlobalVars.speciesIndices[NI]; ion_index <= ChimesGlobalVars.speciesIndices[NVIII]; ion_index++)
            SphP[i].ChimesGasVars.abundances[ion_index] +=
                (SphP[i].ChimesGasVars.abundances[ion_index] / old_abundances[2]) * delta_zi;
        }
      else
        SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[NI]] += delta_zi;
    }

  /* Oxygen */
  if(ChimesGlobalVars.element_included[2] == 1)
    {
      delta_zi = SphP[i].ChimesGasVars.element_abundances[3] - old_abundances[3];
      if(old_abundances[3] > 0.0)
        {
          for(ion_index = ChimesGlobalVars.speciesIndices[OI]; ion_index <= ChimesGlobalVars.speciesIndices[OIX]; ion_index++)
            SphP[i].ChimesGasVars.abundances[ion_index] +=
                (SphP[i].ChimesGasVars.abundances[ion_index] / old_abundances[3]) * delta_zi;

          if(SphP[i].ChimesGasVars.temperature < ChimesGlobalVars.T_mol)
            {
              SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[OH]] +=
                  (SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[OH]] / old_abundances[3]) * delta_zi;
              SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[H2O]] +=
                  (SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[H2O]] / old_abundances[3]) * delta_zi;
              SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[O2]] +=
                  (SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[O2]] / old_abundances[3]) * delta_zi;
              SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[HCOp]] +=
                  (SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[HCOp]] / old_abundances[3]) * delta_zi;
              SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[CO]] +=
                  (SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[CO]] / old_abundances[3]) * delta_zi;
              SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[OHp]] +=
                  (SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[OHp]] / old_abundances[3]) * delta_zi;
              SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[H2Op]] +=
                  (SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[H2Op]] / old_abundances[3]) * delta_zi;
              SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[H3Op]] +=
                  (SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[H3Op]] / old_abundances[3]) * delta_zi;
              SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[COp]] +=
                  (SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[COp]] / old_abundances[3]) * delta_zi;
              SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[HOCp]] +=
                  (SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[HOCp]] / old_abundances[3]) * delta_zi;
              SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[O2p]] +=
                  (SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[O2p]] / old_abundances[3]) * delta_zi;
            }
        }
      else
        SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[OI]] += delta_zi;
    }

  /* Neon */
  if(ChimesGlobalVars.element_included[3] == 1)
    {
      delta_zi = SphP[i].ChimesGasVars.element_abundances[4] - old_abundances[4];
      if(old_abundances[4] > 0.0)
        {
          for(ion_index = ChimesGlobalVars.speciesIndices[NeI]; ion_index <= ChimesGlobalVars.speciesIndices[NeXI]; ion_index++)
            SphP[i].ChimesGasVars.abundances[ion_index] +=
                (SphP[i].ChimesGasVars.abundances[ion_index] / old_abundances[4]) * delta_zi;
        }
      else
        SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[NeI]] += delta_zi;
    }

  /* Magnesium */
  if(ChimesGlobalVars.element_included[4] == 1)
    {
      delta_zi = SphP[i].ChimesGasVars.element_abundances[5] - old_abundances[5];
      if(old_abundances[5] > 0.0)
        {
          for(ion_index = ChimesGlobalVars.speciesIndices[MgI]; ion_index <= ChimesGlobalVars.speciesIndices[MgXIII]; ion_index++)
            SphP[i].ChimesGasVars.abundances[ion_index] +=
                (SphP[i].ChimesGasVars.abundances[ion_index] / old_abundances[5]) * delta_zi;
        }
      else
        SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[MgI]] += delta_zi;
    }

  /* Silicon */
  if(ChimesGlobalVars.element_included[5] == 1)
    {
      delta_zi = SphP[i].ChimesGasVars.element_abundances[6] - old_abundances[6];
      if(old_abundances[6] > 0.0)
        {
          for(ion_index = ChimesGlobalVars.speciesIndices[SiI]; ion_index <= ChimesGlobalVars.speciesIndices[SiXV]; ion_index++)
            SphP[i].ChimesGasVars.abundances[ion_index] +=
                (SphP[i].ChimesGasVars.abundances[ion_index] / old_abundances[6]) * delta_zi;
        }
      else
        SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[SiI]] += delta_zi;
    }

  /* Sulphur */
  if(ChimesGlobalVars.element_included[6] == 1)
    {
      delta_zi = SphP[i].ChimesGasVars.element_abundances[7] - old_abundances[7];
      if(old_abundances[7] > 0.0)
        {
          for(ion_index = ChimesGlobalVars.speciesIndices[SI]; ion_index <= ChimesGlobalVars.speciesIndices[SXVII]; ion_index++)
            SphP[i].ChimesGasVars.abundances[ion_index] +=
                (SphP[i].ChimesGasVars.abundances[ion_index] / old_abundances[7]) * delta_zi;
        }
      else
        SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[SI]] += delta_zi;
    }

  /* Calcium */
  if(ChimesGlobalVars.element_included[7] == 1)
    {
      delta_zi = SphP[i].ChimesGasVars.element_abundances[8] - old_abundances[8];
      if(old_abundances[8] > 0.0)
        {
          for(ion_index = ChimesGlobalVars.speciesIndices[CaI]; ion_index <= ChimesGlobalVars.speciesIndices[CaXXI]; ion_index++)
            SphP[i].ChimesGasVars.abundances[ion_index] +=
                (SphP[i].ChimesGasVars.abundances[ion_index] / old_abundances[8]) * delta_zi;
        }
      else
        SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[CaI]] += delta_zi;
    }

  /* Iron */
  if(ChimesGlobalVars.element_included[8] == 1)
    {
      delta_zi = SphP[i].ChimesGasVars.element_abundances[9] - old_abundances[9];
      if(old_abundances[9] > 0.0)
        {
          for(ion_index = ChimesGlobalVars.speciesIndices[FeI]; ion_index <= ChimesGlobalVars.speciesIndices[FeXXVII]; ion_index++)
            SphP[i].ChimesGasVars.abundances[ion_index] +=
                (SphP[i].ChimesGasVars.abundances[ion_index] / old_abundances[9]) * delta_zi;
        }
      else
        SphP[i].ChimesGasVars.abundances[ChimesGlobalVars.speciesIndices[FeI]] += delta_zi;
    }
}
#endif /* GFM_STELLAR_EVOLUTION */

#ifdef USE_SFR

#if !defined(SMUGGLE_SFR) && !defined(LOCAL_FEEDBACK)
#if !defined(QUICK_LYALPHA)

#ifdef CHIMES_PTHREADS
/* This routine loops over active gas particles and integrates
 * their cooling and chemistry, unless they are above the EoS
 * density and colder than the EoS. It is used within the
 * cooling_and_starformation() routine in sfr_eEOS.c. This
 * version uses the PTHREADS multithreading option.
 */
void chimes_cooling_for_starformation_pthreads(void)
{
  int idx, i, j, bin;
  double dt, dtime;
  double unew, du;
  double factorEVP, dens;
  double tsfr;
  double egyeff, egyhot, x, u_save;
  int N_active = 0;
  int *N_active_node, *buf_offset;

#ifdef BH_THERMALFEEDBACK
  double temp;
#endif

  double eos_dens_threshold = All.PhysDensThresh;
#ifdef MODIFIED_EOS
  eos_dens_threshold *= All.FactorDensThresh;
#endif

  mpi_printf("CHIMES: doing chemistry and cooling.\n");

#ifdef CHIMES_REDSHIFT_DEPENDENT_UVB
  ChimesLoadPhotoIonTables();
#endif

  chimes_update_all_pointers();

  /* clear the SFR stored in the active timebins */
  for(bin = 0; bin < TIMEBINS; bin++)
    if(TimeBinSynchronized[bin])
      TimeBinSfr[bin] = 0;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Mass == 0 && P[i].ID == 0)
        continue; /* skip cells that have been swallowed or eliminated */

      dens = SphP[i].Density;

      dt    = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
      dtime = All.cf_atime * dt / All.cf_time_hubble_a;

      /* apply the temperature floor (and, if relevant, BH thermal feedback) */
      unew = fmax(All.MinEgySpec, SphP[i].Utherm);

#ifdef BH_THERMALFEEDBACK
      if(SphP[i].Injected_BH_Energy)
        {
          if(SphP[i].Injected_BH_Energy < 0)
            terminate("strange feedback energy: Thistask=%d i=%d SphP[i].Injected_BH_Energy=%g\n", ThisTask, i,
                      SphP[i].Injected_BH_Energy);

          unew += SphP[i].Injected_BH_Energy / P[i].Mass;

          temp = chimes_convert_ucgs_to_temp(unew * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs, i);

          if(temp > 5.0e9)
            unew = chimes_convert_temp_to_ucgs(5.0e9, i) * All.UnitDensity_in_cgs / All.UnitPressure_in_cgs;

          AGNEnergyT_Is += SphP[i].Injected_BH_Energy;
          SphP[i].Injected_BH_Energy = 0;
        }
#endif

      if(unew < 0)
        terminate("Invalid Temperature: Task=%d i=%d unew=%g\n", ThisTask, i, unew);

      du = unew - SphP[i].Utherm;
      SphP[i].Utherm += du;
      SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;

      egyeff = 0.;
      /* calculate the effective equation of state for gas above the density threshold */
      if(dens * All.cf_a3inv >= eos_dens_threshold)
        egyeff = chimes_calc_egyeff(i, dens * All.cf_a3inv, &x, &tsfr, &factorEVP);

      if(dens * All.cf_a3inv < eos_dens_threshold || (dens * All.cf_a3inv >= eos_dens_threshold && SphP[i].Utherm > egyeff))
        {
          chimes_update_gas_vars(i);

          if(dens * All.cf_a3inv >= eos_dens_threshold)
            {
              SphP[i].ChimesGasVars.ThermEvolOn = All.ChimesThermEvolOn;
              SphP[i].ChimesGasVars.ForceEqOn   = 1;

              SphP[i].ChimesGasVars.TempFloor =
                  chimes_convert_ucgs_to_temp(egyeff * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs, i);

              if(SphP[i].ChimesGasVars.hydro_timestep > 0.0)
                {
                  chimes_network(&(SphP[i].ChimesGasVars), &ChimesGlobalVars, ChimesAllRatesTh[0], ChimesAllReactionsRootTh[0],
                                 ChimesNonMolecularReactionsRootTh[0]);
                  chimes_update_particle_energies(i);
                }
            }
          else
            {
              SphP[i].ChimesGasVars.ThermEvolOn = All.ChimesThermEvolOn;
              SphP[i].ChimesGasVars.ForceEqOn   = All.ChimesForceEqOn;
              SphP[i].ChimesGasVars.TempFloor   = All.MinGasTemp;
            }
        }
      else
        {
          /* Temporarily set particle i's thermal energy to egyeff while we set the chemistry. */
          u_save         = SphP[i].Utherm;
          SphP[i].Utherm = egyeff;

          chimes_update_gas_vars(i);
          SphP[i].ChimesGasVars.ThermEvolOn = 0; /* Fixed temperature */
          SphP[i].ChimesGasVars.ForceEqOn   = 1; /* Eqm abundances */
          SphP[i].ChimesGasVars.TempFloor   = All.MinGasTemp;

          chimes_network(&(SphP[i].ChimesGasVars), &ChimesGlobalVars, ChimesAllRatesTh[0], ChimesAllReactionsRootTh[0],
                         ChimesNonMolecularReactionsRootTh[0]);

          /* Now reset Utherm */
          SphP[i].Utherm = u_save;
        }

      /* Note that we only record how many active
       * particles are to be evolved with non-eq
       * chemistry here. */
      if(SphP[i].ChimesGasVars.ForceEqOn == 0)
        N_active++;
    } /* end of main loop over active particles */

  N_active_node = (int *)mymalloc_movable(&N_active_node, "N_active", NTask_node * sizeof(int));
  buf_offset    = (int *)mymalloc_movable(&buf_offset, "buf_offset", NTask_node * sizeof(int));

  MPI_Gather(&N_active, 1, MPI_INT, N_active_node, 1, MPI_INT, 0, node_comm);

  MPI_Barrier(node_comm);

  N_active_tot = 0;

  if(ThisTask_node == 0)
    {
      for(j = 0; j < NTask_node; j++)
        {
          buf_offset[j] = N_active_tot;
          N_active_tot += N_active_node[j];
        }
    }

  chimes_create_pthreads_buffers(N_active);

  chimes_send_gasvars_to_node_root(N_active, N_active_node, buf_offset);

  chimes_run_threads();

  chimes_send_abundances_to_original_task(N_active, N_active_node, buf_offset);

  chimes_free_pthreads_buffers(N_active);

  myfree_movable(N_active_node);
  myfree_movable(buf_offset);
}
#else /* CHIMES_PTHREADS */
/* This routine loops over active gas particles and integrates
 * their cooling and chemistry, unless they are above the EoS
 * density and colder than the EoS. It is used within the
 * cooling_and_starformation() routine in sfr_eEOS.c.
 */
void chimes_cooling_for_starformation(void)
{
  int idx, i, bin;
  double dt, dtime;
  double unew, du;
  double factorEVP, dens;
  double tsfr;
  double egyeff, egyhot, x, u_save;

#ifdef BH_THERMALFEEDBACK
  double temp;
#endif

  double eos_dens_threshold = All.PhysDensThresh;
#ifdef MODIFIED_EOS
  eos_dens_threshold *= All.FactorDensThresh;
#endif

  mpi_printf("CHIMES: doing chemistry and cooling.\n");

#ifdef CHIMES_REDSHIFT_DEPENDENT_UVB
  ChimesLoadPhotoIonTables();
#endif

  chimes_update_all_pointers();

  /* clear the SFR stored in the active timebins */
  for(bin = 0; bin < TIMEBINS; bin++)
    if(TimeBinSynchronized[bin])
      TimeBinSfr[bin] = 0;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Mass == 0 && P[i].ID == 0)
        continue; /* skip cells that have been swallowed or eliminated */

      dens = SphP[i].Density;

      dt    = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
      dtime = All.cf_atime * dt / All.cf_time_hubble_a;

      /* apply the temperature floor (and, if relevant, BH thermal feedback) */
      unew = fmax(All.MinEgySpec, SphP[i].Utherm);

#ifdef BH_THERMALFEEDBACK
      if(SphP[i].Injected_BH_Energy)
        {
          if(SphP[i].Injected_BH_Energy < 0)
            terminate("strange feedback energy: Thistask=%d i=%d SphP[i].Injected_BH_Energy=%g\n", ThisTask, i,
                      SphP[i].Injected_BH_Energy);

          unew += SphP[i].Injected_BH_Energy / P[i].Mass;

          temp = chimes_convert_ucgs_to_temp(unew * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs, i);

          if(temp > 5.0e9)
            unew = chimes_convert_temp_to_ucgs(5.0e9, i) * All.UnitDensity_in_cgs / All.UnitPressure_in_cgs;

          AGNEnergyT_Is += SphP[i].Injected_BH_Energy;
          SphP[i].Injected_BH_Energy = 0;
        }
#endif

      if(unew < 0)
        terminate("Invalid Temperature: Task=%d i=%d unew=%g\n", ThisTask, i, unew);

      du = unew - SphP[i].Utherm;
      SphP[i].Utherm += du;
      SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;

      egyeff = 0.;
      /* calculate the effective equation of state for gas above the density threshold */
      if(dens * All.cf_a3inv >= eos_dens_threshold)
        egyeff = chimes_calc_egyeff(i, dens * All.cf_a3inv, &x, &tsfr, &factorEVP);

      /* do non-eqm chemistry and cooling, except for gas above the EOS density threshold that is colder than the eEOS */
      if(dens * All.cf_a3inv < eos_dens_threshold || (dens * All.cf_a3inv >= eos_dens_threshold && SphP[i].Utherm > egyeff))
        {
          chimes_update_gas_vars(i);

          if(dens * All.cf_a3inv >= eos_dens_threshold)
            {
              /* For particles above the density threshold, evolve with equilibrium
               * cooling, and set the temperature floor to the EoS effective temperature. */
              SphP[i].ChimesGasVars.ThermEvolOn = All.ChimesThermEvolOn;
              SphP[i].ChimesGasVars.ForceEqOn   = 1;
              SphP[i].ChimesGasVars.TempFloor =
                  chimes_convert_ucgs_to_temp(egyeff * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs, i);
            }

          if(SphP[i].ChimesGasVars.hydro_timestep > 0.0)
            {
              chimes_network(&(SphP[i].ChimesGasVars), &ChimesGlobalVars, ChimesAllRates, ChimesAllReactionsRoot,
                             ChimesNonMolecularReactionsRoot);
              chimes_update_particle_energies(i);
            }

          if(dens * All.cf_a3inv >= eos_dens_threshold)
            {
              SphP[i].ChimesGasVars.ForceEqOn = All.ChimesForceEqOn;
              SphP[i].ChimesGasVars.TempFloor = All.MinGasTemp;
            }
        }
      else
        {
          /* Set abundances to eqm at the effective temperature of the eEOS
           * We will temporarily set particle i's thermal energy
           * to egyeff while we set the chemistry. */
          u_save         = SphP[i].Utherm;
          SphP[i].Utherm = egyeff;

          chimes_update_gas_vars(i);
          SphP[i].ChimesGasVars.ForceEqOn   = 1; /* Eqm abundances */
          SphP[i].ChimesGasVars.ThermEvolOn = 0; /* Fixed temperature */

          chimes_network(&(SphP[i].ChimesGasVars), &ChimesGlobalVars, ChimesAllRates, ChimesAllReactionsRoot,
                         ChimesNonMolecularReactionsRoot);

          /* Now reset Utherm and CHIMES flags. */
          SphP[i].Utherm                    = u_save;
          SphP[i].ChimesGasVars.ForceEqOn   = All.ChimesForceEqOn;
          SphP[i].ChimesGasVars.ThermEvolOn = All.ChimesThermEvolOn;
        }
    } /* end of main loop over active particles */
}
#endif /* CHIMES_PTHREADS */

/* Calculate the effective energy of the multi-phase model, using
 * CHIMES cooling and equilibrium abundances.
 */
double chimes_calc_egyeff(int i, double gasdens, double *x, double *tsfr, double *factorEVP)
{
  double egyhot, egyeff, tcool, y, u_save;
  double *abundances_save;
  int j;
  double rho = gasdens;

  chimes_set_pointers(i);

#ifndef MODIFIED_EOS
  if(rho < All.PhysDensThresh)
    terminate("Error in gas cell %d, Effective equation of state can be computed only for gas above SF density threshold", i);
#else
  rho = fmax(rho, All.PhysDensThresh);
#endif

  *tsfr = sqrt(All.PhysDensThresh / rho) * All.MaxSfrTimescale;

  *factorEVP = pow(rho / All.PhysDensThresh, -0.8) * All.FactorEVP;

  egyhot = All.EgySpecSN / (1 + *factorEVP) + All.EgySpecCold;

  /* Compute equilibrium abundances of the hot phase */

  /* We will temporarily set particle i's thermal energy to egyhot while we set the chemistry. */
  u_save         = SphP[i].Utherm;
  SphP[i].Utherm = egyhot; /* egyhot already in code units. */

  /* We also need to save the particle's abundances. We only
   * set them to eqm here to compute the egyeff, we don't
   * want to actually set the particle to chemical eqm here. */

  abundances_save = (double *)mymalloc_movable(&abundances_save, "abun_save", ChimesGlobalVars.totalNumberOfSpecies * sizeof(double));
  for(j = 0; j < ChimesGlobalVars.totalNumberOfSpecies; j++)
    abundances_save[j] = SphP[i].ChimesGasVars.abundances[j];

  chimes_update_gas_vars(i);
  SphP[i].ChimesGasVars.ForceEqOn   = 1; /* Eqm abundances */
  SphP[i].ChimesGasVars.ThermEvolOn = 0; /* Fixed temperature */

#ifdef CHIMES_DISABLE_SHIELDING_ON_EOS
  SphP[i].ChimesGasVars.cell_size = 1.0;
#endif

#ifdef CHIMES_DISABLE_ZCOOL_ON_EOS
  double save_element_abundances[9];
  for(j = 1; j < 10; j++)
    {
      save_element_abundances[j - 1]              = SphP[i].ChimesGasVars.element_abundances[j];
      SphP[i].ChimesGasVars.element_abundances[j] = 0.0;
    }
#endif

#ifdef CHIMES_PTHREADS
  chimes_network(&(SphP[i].ChimesGasVars), &ChimesGlobalVars, ChimesAllRatesTh[0], ChimesAllReactionsRootTh[0],
                 ChimesNonMolecularReactionsRootTh[0]);
#else
  chimes_network(&(SphP[i].ChimesGasVars), &ChimesGlobalVars, ChimesAllRates, ChimesAllReactionsRoot, ChimesNonMolecularReactionsRoot);
#endif

  /* Compute the cooling time of the hot phase. */
#ifdef CHIMES_PTHREADS
  tcool = ChimesGetCoolingTime(i, ChimesAllRatesTh[0]);
#else
  tcool = ChimesGetCoolingTime(i, ChimesAllRates);
#endif

  /* Now reset Utherm and CHIMES flags, and abundances. */
  SphP[i].Utherm                    = u_save;
  SphP[i].ChimesGasVars.ForceEqOn   = All.ChimesForceEqOn;
  SphP[i].ChimesGasVars.ThermEvolOn = All.ChimesThermEvolOn;

  for(j = 0; j < ChimesGlobalVars.totalNumberOfSpecies; j++)
    SphP[i].ChimesGasVars.abundances[j] = abundances_save[j];
  myfree_movable(abundances_save);

#ifdef CHIMES_DISABLE_ZCOOL_ON_EOS
  for(j = 1; j < 10; j++)
    SphP[i].ChimesGasVars.element_abundances[j] = save_element_abundances[j - 1];
#endif

  y = *tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);

  *x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

  egyeff = egyhot * (1 - *x) + All.EgySpecCold * (*x);

#ifdef SOFTEREQS
  /* use an intermediate EQS, between isothermal and the full multiphase model */
  egyeff = All.FactorForSofterEQS * egyeff + (1 - All.FactorForSofterEQS) * All.UForSofterEQS;
#endif

#ifdef STEEPER_SFR_FOR_STARBURST
  /* modify sfr at high densities, but don't change the energy calculated above */
  if(rho > All.PhysDensThreshStarburst)
    *tsfr = pow(All.PhysDensThreshStarburst / rho, All.StarburstPowerLawIndex) *
            sqrt(All.PhysDensThresh / All.PhysDensThreshStarburst) * All.MaxSfrTimescale;
#endif

  return egyeff;
}

#endif /* !(QUICK_LYALPHA) */
#endif /* !(SMUGGLE_SFR) && !(LOCAL_FEEDBACK) */
#endif /* USE_SFR */

#ifdef CHIMES_PTHREADS
/* Allocate memory to the CHIMES buffers, and read gas vars
 * of active particles into these buffers. */
void chimes_create_pthreads_buffers(int N_active)
{
  int idx, i, k;

  if(ThisTask_node == 0)
    {
      if(N_active_tot > 0)
        {
          ChimesGasVars_buf =
              (struct gasVariables *)mymalloc_movable(&ChimesGasVars_buf, "gasVars_buf", N_active_tot * sizeof(struct gasVariables));
          ChimesAbundances_buf    = (double *)mymalloc_movable(&ChimesAbundances_buf, "abun_buf",
                                                               N_active_tot * ChimesGlobalVars.totalNumberOfSpecies * sizeof(double));
          ChimesPhotonDensity_buf = (double *)mymalloc_movable(&ChimesPhotonDensity_buf, "photon_buf",
                                                               N_active_tot * ChimesGlobalVars.N_spectra * sizeof(double));
          ChimesH2dissocJ_buf     = (double *)mymalloc_movable(&ChimesH2dissocJ_buf, "H2dissocJ_buf",
                                                               N_active_tot * ChimesGlobalVars.N_spectra * sizeof(double));
          ChimesDustG_buf =
              (double *)mymalloc_movable(&ChimesDustG_buf, "dustG_buf", N_active_tot * ChimesGlobalVars.N_spectra * sizeof(double));
          ChimesTemperature_buf = (double *)mymalloc_movable(&ChimesTemperature_buf, "T_buf", N_active_tot * sizeof(double));
        }
    }
  else
    {
      if(N_active > 0)
        {
          ChimesGasVars_buf =
              (struct gasVariables *)mymalloc_movable(&ChimesGasVars_buf, "gasVars_buf", N_active * sizeof(struct gasVariables));
          ChimesAbundances_buf    = (double *)mymalloc_movable(&ChimesAbundances_buf, "abun_buf",
                                                               N_active * ChimesGlobalVars.totalNumberOfSpecies * sizeof(double));
          ChimesPhotonDensity_buf = (double *)mymalloc_movable(&ChimesPhotonDensity_buf, "photon_buf",
                                                               N_active * ChimesGlobalVars.N_spectra * sizeof(double));
          ChimesH2dissocJ_buf     = (double *)mymalloc_movable(&ChimesH2dissocJ_buf, "H2dissocJ_buf",
                                                               N_active * ChimesGlobalVars.N_spectra * sizeof(double));
          ChimesDustG_buf =
              (double *)mymalloc_movable(&ChimesDustG_buf, "dustG_buf", N_active * ChimesGlobalVars.N_spectra * sizeof(double));
          ChimesTemperature_buf = (double *)mymalloc_movable(&ChimesTemperature_buf, "T_buf", N_active * sizeof(double));
        }
    }

  if(N_active > 0)
    {
      /* Read this task's active particles into the buffers */
      buf_index = 0;
      for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i >= 0)
            {
              if(P[i].Mass == 0 && P[i].ID == 0)
                continue; /* skip cells that have been swallowed or eliminated */

#ifdef USE_SFR
              if(SphP[i].ChimesGasVars.ForceEqOn == 1)
                continue;
#endif

              ChimesGasVars_buf[buf_index] = SphP[i].ChimesGasVars;
              for(k = 0; k < ChimesGlobalVars.totalNumberOfSpecies; k++)
                ChimesAbundances_buf[(ChimesGlobalVars.totalNumberOfSpecies * buf_index) + k] = SphP[i].ChimesGasVars.abundances[k];
              for(k = 0; k < ChimesGlobalVars.N_spectra; k++)
                {
                  ChimesPhotonDensity_buf[(ChimesGlobalVars.N_spectra * buf_index) + k] =
                      SphP[i].ChimesGasVars.isotropic_photon_density[k];
                  ChimesH2dissocJ_buf[(ChimesGlobalVars.N_spectra * buf_index) + k] = SphP[i].ChimesGasVars.H2_dissocJ[k];
                  ChimesDustG_buf[(ChimesGlobalVars.N_spectra * buf_index) + k]     = SphP[i].ChimesGasVars.dust_G_parameter[k];
                }
              buf_index++;
            }
        }
    }
}

/* Send the CHIMES buffers to the root MPI task
 * on this node. */
void chimes_send_gasvars_to_node_root(int N_active, int *N_active_node, int *buf_offset)
{
  int j, k;
  int n_requests = 0;
  MPI_Request *requests;

  if(ThisTask_node == 0)
    {
      if(N_active_tot > 0)
        {
          /* Recv gas vars from other tasks */
          requests = (MPI_Request *)mymalloc_movable(&requests, "chimes_requests", NTask_node * 5 * sizeof(MPI_Request));
          for(j = 1; j < NTask_node; j++)
            {
              if(N_active_node[j] > 0)
                {
                  MPI_Irecv(&ChimesGasVars_buf[buf_offset[j]], N_active_node[j] * sizeof(struct gasVariables), MPI_BYTE, j,
                            TAG_CHIMESDATA, node_comm, &requests[n_requests++]);
                  MPI_Irecv(&ChimesAbundances_buf[ChimesGlobalVars.totalNumberOfSpecies * buf_offset[j]],
                            ChimesGlobalVars.totalNumberOfSpecies * N_active_node[j], MPI_DOUBLE, j, TAG_ABUNDATA, node_comm,
                            &requests[n_requests++]);
                  MPI_Irecv(&ChimesPhotonDensity_buf[ChimesGlobalVars.N_spectra * buf_offset[j]],
                            ChimesGlobalVars.N_spectra * N_active_node[j], MPI_DOUBLE, j, TAG_PHOTDATA, node_comm,
                            &requests[n_requests++]);
                  MPI_Irecv(&ChimesH2dissocJ_buf[ChimesGlobalVars.N_spectra * buf_offset[j]],
                            ChimesGlobalVars.N_spectra * N_active_node[j], MPI_DOUBLE, j, TAG_PHOTDATA, node_comm,
                            &requests[n_requests++]);
                  MPI_Irecv(&ChimesDustG_buf[ChimesGlobalVars.N_spectra * buf_offset[j]],
                            ChimesGlobalVars.N_spectra * N_active_node[j], MPI_DOUBLE, j, TAG_PHOTDATA, node_comm,
                            &requests[n_requests++]);
                }
            }
        }
    }
  else
    {
      if(N_active > 0)
        {
          /* Send gas vars to root */
          requests = (MPI_Request *)mymalloc_movable(&requests, "chimes_requests", 5 * sizeof(MPI_Request));
          MPI_Isend(ChimesGasVars_buf, N_active * sizeof(struct gasVariables), MPI_BYTE, 0, TAG_CHIMESDATA, node_comm,
                    &requests[n_requests++]);
          MPI_Isend(ChimesAbundances_buf, ChimesGlobalVars.totalNumberOfSpecies * N_active, MPI_DOUBLE, 0, TAG_ABUNDATA, node_comm,
                    &requests[n_requests++]);
          MPI_Isend(ChimesPhotonDensity_buf, ChimesGlobalVars.N_spectra * N_active, MPI_DOUBLE, 0, TAG_PHOTDATA, node_comm,
                    &requests[n_requests++]);
          MPI_Isend(ChimesH2dissocJ_buf, ChimesGlobalVars.N_spectra * N_active, MPI_DOUBLE, 0, TAG_PHOTDATA, node_comm,
                    &requests[n_requests++]);
          MPI_Isend(ChimesDustG_buf, ChimesGlobalVars.N_spectra * N_active, MPI_DOUBLE, 0, TAG_PHOTDATA, node_comm,
                    &requests[n_requests++]);
        }
    }

  /* Wait for all tasks to finish sending to root. */
  if(n_requests > 0)
    MPI_Waitall(n_requests, requests, MPI_STATUSES_IGNORE);

  if(ThisTask_node == 0)
    {
      if(N_active_tot > 0)
        {
          /* Allocate arrays within GasVars to corresponding section of the buffers */
          for(j = 0; j < N_active_tot; j++)
            {
              ChimesGasVars_buf[j].abundances               = &(ChimesAbundances_buf[ChimesGlobalVars.totalNumberOfSpecies * j]);
              ChimesGasVars_buf[j].isotropic_photon_density = &(ChimesPhotonDensity_buf[ChimesGlobalVars.N_spectra * j]);
              ChimesGasVars_buf[j].H2_dissocJ               = &(ChimesH2dissocJ_buf[ChimesGlobalVars.N_spectra * j]);
              ChimesGasVars_buf[j].dust_G_parameter         = &(ChimesDustG_buf[ChimesGlobalVars.N_spectra * j]);
            }
        }
    }

  if((ThisTask_node == 0 && N_active_tot > 0) || (ThisTask_node != 0 && N_active > 0))
    myfree_movable(requests);
}

void chimes_run_threads(void)
{
  int N_thread, j;
  int thread_id[CHIMES_PTHREADS];
  pthread_t callThread[CHIMES_PTHREADS];
  pthread_attr_t thread_attr;
  struct chimes_thread_data my_data[CHIMES_PTHREADS];
  int node_buffer = 0;
  int node_flag   = 0;
  struct timespec time_pause;
  MPI_Request *node_req;

  time_pause.tv_sec  = 0;
  time_pause.tv_nsec = 10000000L; /* 10 ms */
  node_req           = (MPI_Request *)mymalloc_movable(&node_req, "node_req", NTask_node * sizeof(MPI_Request));

  if(ThisTask_node == 0)
    {
      if(N_active_tot > 0)
        {
          /* Create threads and do the chemistry and cooling for all active particles in the buffer. */
          pthread_mutex_init(&mutexcool, NULL);
          pthread_mutex_init(&mutex_chimes_malloc, NULL);
          pthread_attr_init(&thread_attr);
          pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);

          N_thread  = min(N_active_tot, CHIMES_PTHREADS);
          buf_index = 0;
          for(j = 0; j < N_thread; j++)
            {
              my_data[j].thread_id = j;

              if(pthread_create(&callThread[j], &thread_attr, chimes_worker_thread_routine, (void *)&my_data[j]) != 0)
                terminate("Error: pthread_create failed on task %d, thread %d\n", ThisTask, j);
            }

          for(j = 0; j < N_thread; j++)
            pthread_join(callThread[j], NULL);

          pthread_attr_destroy(&thread_attr);

          pthread_mutex_destroy(&mutexcool);
          pthread_mutex_destroy(&mutex_chimes_malloc);
        } /* N_active_tot > 0 */

      /* We need to ensure that the other tasks on the node wait while
       * ThisTask_node == 0 loops through the particles. However, a simple
       * MPI_Barrier will usually still use up 100% of the CPU. Below, we
       * periodically use MPI_Iprobe to determine when the chemistry
       * has finished. */
      for(j = 1; j < NTask_node; j++)
        MPI_Isend(&node_buffer, 1, MPI_INT, j, TAG_CHIMES_PTHREAD, node_comm, &node_req[j]);
    } /* ThisTask_node == 0 */
  else
    {
      node_flag = 0;
      while(!node_flag)
        {
          MPI_Iprobe(0, TAG_CHIMES_PTHREAD, node_comm, &node_flag, MPI_STATUS_IGNORE);
          nanosleep(&time_pause, NULL);
        }

      MPI_Irecv(&node_buffer, 1, MPI_INT, 0, TAG_CHIMES_PTHREAD, node_comm, &node_req[ThisTask_node]);
    }

  MPI_Barrier(node_comm);

  myfree_movable(node_req);
}

void *chimes_worker_thread_routine(void *arg)
{
  struct chimes_thread_data *my_data = (struct chimes_thread_data *)arg;
  int thread_id                      = my_data->thread_id;
  int this_thread_index;

  while(1)
    {
      pthread_mutex_lock(&mutexcool);
      if(buf_index >= N_active_tot)
        {
          pthread_mutex_unlock(&mutexcool);
          break;
        }
      else
        {
          this_thread_index = buf_index;
          buf_index++;
        }

      pthread_mutex_unlock(&mutexcool);

      if(ChimesGasVars_buf[this_thread_index].hydro_timestep > 0.0)
        chimes_network(&(ChimesGasVars_buf[this_thread_index]), &(ChimesGlobalVars), ChimesAllRatesTh[thread_id],
                       ChimesAllReactionsRootTh[thread_id], ChimesNonMolecularReactionsRootTh[thread_id]);
    }

  pthread_exit((void *)0);
}

/* Send CHIMES abundances from the root task on the node
 * back to their original tasks. */
void chimes_send_abundances_to_original_task(int N_active, int *N_active_node, int *buf_offset)
{
  int n_requests, idx, i, j, k;
  MPI_Request *requests;

  n_requests = 0;
  if(ThisTask_node == 0)
    {
      if(N_active_tot > 0)
        {
          requests = (MPI_Request *)mymalloc_movable(&requests, "chimes_requests", NTask_node * 2 * sizeof(MPI_Request));

          /* Read temperatures into buffer */
          /* Abundances were already updated within their buffer */
          for(j = 0; j < N_active_tot; j++)
            ChimesTemperature_buf[j] = ChimesGasVars_buf[j].temperature;

          /* Send back to original tasks */
          n_requests = 0;
          for(j = 1; j < NTask_node; j++)
            {
              if(N_active_node[j] > 0)
                {
                  MPI_Isend(&ChimesAbundances_buf[ChimesGlobalVars.totalNumberOfSpecies * buf_offset[j]],
                            ChimesGlobalVars.totalNumberOfSpecies * N_active_node[j], MPI_DOUBLE, j, TAG_ABUNDATA, node_comm,
                            &requests[n_requests++]);
                  MPI_Isend(&ChimesTemperature_buf[buf_offset[j]], N_active_node[j], MPI_DOUBLE, j, TAG_TEMPDATA, node_comm,
                            &requests[n_requests++]);
                }
            }
        } /* N_active_tot > 0 */
    }     /* ThisTask_node == 0 */
  else
    {
      if(N_active > 0)
        {
          requests = (MPI_Request *)mymalloc_movable(&requests, "chimes_requests", NTask_node * 2 * sizeof(MPI_Request));

          MPI_Irecv(ChimesAbundances_buf, ChimesGlobalVars.totalNumberOfSpecies * N_active, MPI_DOUBLE, 0, TAG_ABUNDATA, node_comm,
                    &requests[n_requests++]);
          MPI_Irecv(ChimesTemperature_buf, N_active, MPI_DOUBLE, 0, TAG_TEMPDATA, node_comm, &requests[n_requests++]);
        }
    }

  if(n_requests > 0)
    MPI_Waitall(n_requests, requests, MPI_STATUSES_IGNORE);

  /* Now write abundances and temperatures back into the full gas vars array,
   * and update the particle energies */
  buf_index = 0;
  if(N_active > 0)
    {
      for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i >= 0)
            {
              if(P[i].Mass == 0 && P[i].ID == 0)
                continue; /* skip cells that have been swallowed or eliminated */

#ifdef USE_SFR
              if(SphP[i].ChimesGasVars.ForceEqOn == 1)
                continue;
#endif

              SphP[i].ChimesGasVars.temperature = ChimesTemperature_buf[buf_index];
              for(k = 0; k < ChimesGlobalVars.totalNumberOfSpecies; k++)
                SphP[i].ChimesGasVars.abundances[k] = ChimesAbundances_buf[(ChimesGlobalVars.totalNumberOfSpecies * buf_index) + k];
              buf_index++;

              chimes_update_particle_energies(i);
            }
        }
    }

  if((ThisTask_node == 0 && N_active_tot > 0) || (ThisTask_node != 0 && N_active > 0))
    myfree_movable(requests);
}

void chimes_free_pthreads_buffers(int N_active)
{
  int j;

  if(ThisTask_node == 0)
    {
      if(N_active_tot > 0)
        {
          myfree_movable(ChimesGasVars_buf);
          myfree_movable(ChimesAbundances_buf);
          myfree_movable(ChimesPhotonDensity_buf);
          myfree_movable(ChimesH2dissocJ_buf);
          myfree_movable(ChimesDustG_buf);
          myfree_movable(ChimesTemperature_buf);
        }
    }
  else
    {
      if(N_active > 0)
        {
          myfree_movable(ChimesGasVars_buf);
          myfree_movable(ChimesAbundances_buf);
          myfree_movable(ChimesPhotonDensity_buf);
          myfree_movable(ChimesH2dissocJ_buf);
          myfree_movable(ChimesDustG_buf);
          myfree_movable(ChimesTemperature_buf);
        }
    }
}
#endif /* CHIMES_PTHREADS */

#ifdef CHIMES_REDSHIFT_DEPENDENT_UVB
void chimes_redshift_index(double *table, int ntable, double x, int *i, double *dx)
{
  if(x <= table[0])
    {
      *i  = 0;
      *dx = 0;
    }
  else if(x >= table[ntable - 1])
    {
      *i  = ntable - 1;
      *dx = 1;
    }
  else
    {
      *i = 0;
      while(table[*i] < x)
        *i += 1;

      *i -= 1;
      *dx = (x - table[*i]) / (table[(*i) + 1] - table[*i]);
    }
}

/* This routine determines whether we need to load
 * new PhotoIon tables for the current redshift. */
void ChimesLoadPhotoIonTables(void)
{
  static int ChimesTablesIndex_lowz  = -1;
  static int ChimesTablesIndex_highz = -1;
  char table_path[500];
  int z_index;
  double dz;

  if(All.ComovingIntegrationOn == 0)
    terminate("CHIMES_REDSHIFT_DEPENDENT_UVB is switched on, but ComovingIntegrationOn == 0. Aborting.\n");

  chimes_redshift_index(Chimes_Redshifts, Chimes_N_Redshifts, All.cf_redshift, &z_index, &dz);

  if(ChimesTablesIndex_lowz == z_index)
    {
      ChimesInterpolatePhotoIonTables(&ChimesGlobalVars, ChimesPhotoIonTable, z_index, dz);
      return;
    }

  if(z_index >= Chimes_N_Redshifts - 1)
    {
      /* Current redshift is before highest redshift
       * table. Load highest redshift table into both
       * lowz and hiz tables. */
      sprintf(table_path, "%s/HM12_cross_sections/z%1.3f_cross_sections_v8.hdf5", All.ChimesDataPath,
              Chimes_Redshifts[Chimes_N_Redshifts - 1]);

      if(NTask < 11)
        ChimesReadPhotoIonTables_UVB(&ChimesGlobalVars, table_path, ChimesPhotoIonTable, Chimes_N_Elements_in_Bens_tables, 0);
      else
        ChimesReadPhotoIonTables_UVB_parallel(&ChimesGlobalVars, table_path, ChimesPhotoIonTable, Chimes_N_Elements_in_Bens_tables, 0);

      ChimesCopyPhotoIonTables(ChimesPhotoIonTable, 0, 1);

      ChimesTablesIndex_highz = z_index;
      ChimesTablesIndex_lowz  = z_index;
    }
  else
    {
      if(ChimesTablesIndex_lowz > -1)
        {
          /* Copy lowz table to hiz, then
           * read in new lowz table. */
          ChimesCopyPhotoIonTables(ChimesPhotoIonTable, 0, 1);
          sprintf(table_path, "%s/HM12_cross_sections/z%1.3f_cross_sections_v8.hdf5", All.ChimesDataPath, Chimes_Redshifts[z_index]);

          if(NTask < 11)
            ChimesReadPhotoIonTables_UVB(&ChimesGlobalVars, table_path, ChimesPhotoIonTable, Chimes_N_Elements_in_Bens_tables, 0);
          else
            ChimesReadPhotoIonTables_UVB_parallel(&ChimesGlobalVars, table_path, ChimesPhotoIonTable, Chimes_N_Elements_in_Bens_tables,
                                                  0);
        }
      else
        {
          /* This is the first time any PhotoIon tables have been
           * read in. Read in both lowz and hiz tables. */
          sprintf(table_path, "%s/HM12_cross_sections/z%1.3f_cross_sections_v8.hdf5", All.ChimesDataPath, Chimes_Redshifts[z_index]);

          if(NTask < 11)
            ChimesReadPhotoIonTables_UVB(&ChimesGlobalVars, table_path, ChimesPhotoIonTable, Chimes_N_Elements_in_Bens_tables, 0);
          else
            ChimesReadPhotoIonTables_UVB_parallel(&ChimesGlobalVars, table_path, ChimesPhotoIonTable, Chimes_N_Elements_in_Bens_tables,
                                                  0);

          sprintf(table_path, "%s/HM12_cross_sections/z%1.3f_cross_sections_v8.hdf5", All.ChimesDataPath,
                  Chimes_Redshifts[z_index + 1]);

          if(NTask < 11)
            ChimesReadPhotoIonTables_UVB(&ChimesGlobalVars, table_path, ChimesPhotoIonTable, Chimes_N_Elements_in_Bens_tables, 1);
          else
            ChimesReadPhotoIonTables_UVB_parallel(&ChimesGlobalVars, table_path, ChimesPhotoIonTable, Chimes_N_Elements_in_Bens_tables,
                                                  1);
        }

      ChimesTablesIndex_lowz  = z_index;
      ChimesTablesIndex_highz = z_index + 1;
    }

  ChimesInterpolatePhotoIonTables(&ChimesGlobalVars, ChimesPhotoIonTable, z_index, dz);
}

void chimes_allocate_memory_to_photoion_tables(struct globalVariables *myGlobalVars, struct PhotoIonTables_UVB **myPhotoIonTable,
                                               int N_Elements_in_Bens_tables)
{
  /* The myPhotoIonTable structure will contain two sets of tables, one for
   * the lower redshift and one for the higher redshift, that bracket the
   * current redshift (hence N_spectra = 2). */
  int ns, l, j, k, m, r;
  int N_arrayCells;
  int N_spectra = 2;

  *myPhotoIonTable = (struct PhotoIonTables_UVB *)mymalloc("Chimes_photoIonTable", sizeof(struct PhotoIonTables_UVB));
  (*myPhotoIonTable)->tables_per_element = (struct PhotoIonTables_per_element *)mymalloc(
      "Chimes_N_elem", N_Elements_in_Bens_tables * sizeof(struct PhotoIonTables_per_element));

  ns = 0;
  for(l = 0; l < N_Elements_in_Bens_tables; l++)
    {
      if(l < 2 || myGlobalVars->element_included[max(l - 2, 0)] == 1)
        {
          if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
            N_arrayCells = 3; /* Hydrogen photoionisation also includes H- & H2*/
          else
            N_arrayCells = chimesRateTables.NonEqIon->N_Ions[ns] - 1;

          (*myPhotoIonTable)->tables_per_element[ns].sigmaphot =
              (double ***)mymalloc("Chimes_sigmaphot_dim0", N_spectra * sizeof(double **));
          (*myPhotoIonTable)->tables_per_element[ns].sigmaphot[0] =
              (double **)mymalloc("Chimes_sigmaphot_dim1", N_spectra * N_arrayCells * sizeof(double *));
          (*myPhotoIonTable)->tables_per_element[ns].sigmaphot[0][0] = (double *)mymalloc(
              "Chimes_sigmaphot_dim2", N_spectra * N_arrayCells * chimesRateTables.NonEqIon->N_Auger[ns] * sizeof(double));
          for(j = 0; j < N_spectra; j++)
            {
              (*myPhotoIonTable)->tables_per_element[ns].sigmaphot[j] =
                  &((*myPhotoIonTable)->tables_per_element[ns].sigmaphot[0][j * N_arrayCells]);
              for(k = 0; k < N_arrayCells; k++)
                (*myPhotoIonTable)->tables_per_element[ns].sigmaphot[j][k] =
                    &((*myPhotoIonTable)
                          ->tables_per_element[ns]
                          .sigmaphot[0][0][((j * N_arrayCells) + k) * chimesRateTables.NonEqIon->N_Auger[ns]]);
            }

          (*myPhotoIonTable)->tables_per_element[ns].shieldFactor1D =
              (float ****)mymalloc("Chimes_shield1D_dim0", N_spectra * sizeof(float ***));
          (*myPhotoIonTable)->tables_per_element[ns].shieldFactor1D[0] =
              (float ***)mymalloc("Chimes_shield1D_dim1", N_spectra * 3 * sizeof(float **));
          (*myPhotoIonTable)->tables_per_element[ns].shieldFactor1D[0][0] =
              (float **)mymalloc("Chimes_shield1D_dim2", N_spectra * 3 * N_arrayCells * sizeof(float *));
          (*myPhotoIonTable)->tables_per_element[ns].shieldFactor1D[0][0][0] =
              (float *)mymalloc("Chimes_shield1D_dim3", N_spectra * 3 * N_arrayCells *
                                                            chimesRateTables.NonEqIon->shieldingColumnDimensions[0] * sizeof(float));
          for(j = 0; j < N_spectra; j++)
            {
              (*myPhotoIonTable)->tables_per_element[ns].shieldFactor1D[j] =
                  &((*myPhotoIonTable)->tables_per_element[ns].shieldFactor1D[0][j * 3]);
              for(k = 0; k < 3; k++)
                {
                  (*myPhotoIonTable)->tables_per_element[ns].shieldFactor1D[j][k] =
                      &((*myPhotoIonTable)->tables_per_element[ns].shieldFactor1D[0][0][((j * 3) + k) * N_arrayCells]);
                  for(m = 0; m < N_arrayCells; m++)
                    (*myPhotoIonTable)->tables_per_element[ns].shieldFactor1D[j][k][m] =
                        &((*myPhotoIonTable)
                              ->tables_per_element[ns]
                              .shieldFactor1D[0][0][0][((((j * 3) + k) * N_arrayCells) + m) *
                                                       chimesRateTables.NonEqIon->shieldingColumnDimensions[0]]);
                }
            }

          (*myPhotoIonTable)->tables_per_element[ns].shieldFactor2D =
              (float *****)mymalloc("Chimes_shield2D_dim0", N_spectra * sizeof(float ****));
          (*myPhotoIonTable)->tables_per_element[ns].shieldFactor2D[0] =
              (float ****)mymalloc("Chimes_shield2D_dim1", N_spectra * 6 * sizeof(float ***));
          (*myPhotoIonTable)->tables_per_element[ns].shieldFactor2D[0][0] =
              (float ***)mymalloc("Chimes_shield2D_dim2", N_spectra * 6 * N_arrayCells * sizeof(float **));
          (*myPhotoIonTable)->tables_per_element[ns].shieldFactor2D[0][0][0] = (float **)mymalloc(
              "Chimes_shield2D_dim3",
              N_spectra * 6 * N_arrayCells * chimesRateTables.NonEqIon->shieldingColumnDimensions[0] * sizeof(float *));
          (*myPhotoIonTable)->tables_per_element[ns].shieldFactor2D[0][0][0][0] = (float *)mymalloc(
              "Chimes_shield2D_dim4", N_spectra * 6 * N_arrayCells * chimesRateTables.NonEqIon->shieldingColumnDimensions[0] *
                                          chimesRateTables.NonEqIon->shieldingColumnDimensions[0] * sizeof(float));
          for(m = 0; m < N_spectra; m++)
            {
              (*myPhotoIonTable)->tables_per_element[ns].shieldFactor2D[m] =
                  &((*myPhotoIonTable)->tables_per_element[ns].shieldFactor2D[0][m * 6]);
              for(j = 0; j < 6; j++)
                {
                  (*myPhotoIonTable)->tables_per_element[ns].shieldFactor2D[m][j] =
                      &((*myPhotoIonTable)->tables_per_element[ns].shieldFactor2D[0][0][((m * 6) + j) * N_arrayCells]);
                  for(k = 0; k < N_arrayCells; k++)
                    {
                      (*myPhotoIonTable)->tables_per_element[ns].shieldFactor2D[m][j][k] =
                          &((*myPhotoIonTable)
                                ->tables_per_element[ns]
                                .shieldFactor2D[0][0][0][((((m * 6) + j) * N_arrayCells) + k) *
                                                         chimesRateTables.NonEqIon->shieldingColumnDimensions[0]]);
                      for(r = 0; r < chimesRateTables.NonEqIon->shieldingColumnDimensions[0]; r++)
                        (*myPhotoIonTable)->tables_per_element[ns].shieldFactor2D[m][j][k][r] =
                            &((*myPhotoIonTable)
                                  ->tables_per_element[ns]
                                  .shieldFactor2D[0][0][0][0][((((((m * 6) + j) * N_arrayCells) + k) *
                                                                chimesRateTables.NonEqIon->shieldingColumnDimensions[0]) +
                                                               r) *
                                                              chimesRateTables.NonEqIon->shieldingColumnDimensions[0]]);
                    }
                }
            }

          (*myPhotoIonTable)->tables_per_element[ns].epsilon = (double **)mymalloc("Chimes_eps_dim0", N_spectra * sizeof(double *));
          (*myPhotoIonTable)->tables_per_element[ns].epsilon[0] =
              (double *)mymalloc("Chimes_eps_dim1", N_spectra * N_arrayCells * sizeof(double));
          for(j = 0; j < N_spectra; j++)
            (*myPhotoIonTable)->tables_per_element[ns].epsilon[j] =
                &((*myPhotoIonTable)->tables_per_element[ns].epsilon[0][j * N_arrayCells]);

          (*myPhotoIonTable)->dust_G_parameter         = (double *)mymalloc("Chimes_dustG", N_spectra * sizeof(double));
          (*myPhotoIonTable)->H2_dissocJ               = (double *)mymalloc("Chimes_H2_dissocJ", N_spectra * sizeof(double));
          (*myPhotoIonTable)->isotropic_photon_density = (double *)mymalloc("Chimes_iso_phot", N_spectra * sizeof(double));

          ns += 1;
        }
    }
}

void ChimesReadPhotoIonTables_UVB(struct globalVariables *myGlobalVars, char *photoIonTablePath,
                                  struct PhotoIonTables_UVB *myPhotoIonTable, int N_Elements_in_Bens_tables, int current_spectrum)
{
  /* NOTE: These tables are stored as logs. This makes it quicker
   * when we interpolate the tables, because we don't have to
   * take the log each time. However, the shield factors were
   * already stored as logs. */
  hid_t file_id, dataset, dataspace_id, memspace_id;
  herr_t status;
  hsize_t dims[2];
  hsize_t dims2D[2], count2D[2], offset2D[2];
  hsize_t dims3D[3], count3D[3], offset3D[3];
  hsize_t dims4D[4], count4D[4], offset4D[4];
  int rank, i, j, k, l, m, ns, r, N_arrayCells;
  char set_name[500];
  double *sigmaphot, *epsilon;
  double n_phot, H2_dissocJ, dust_G_parameter;
  float *shieldFac_1D, *shieldFac_2D;
  dims[1] = 1;

  if(ThisTask == 0)
    printf("CHIMES: Reading PhotoIon table: %s \n", photoIonTablePath);

  file_id = H5Fopen(photoIonTablePath, H5F_ACC_RDONLY, H5P_DEFAULT);

  ns = 0;
  for(l = 0; l < N_Elements_in_Bens_tables; l++)
    {
      if(l < 2 || myGlobalVars->element_included[max(l - 2, 0)] == 1)
        {
          if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
            N_arrayCells = 3; /* Hydrogen photoionisation also includes H- & H2*/
          else
            N_arrayCells = chimesRateTables.NonEqIon->N_Ions[ns] - 1;

          sigmaphot    = (double *)mymalloc("Chimes_sigma", N_arrayCells * sizeof(double));
          epsilon      = (double *)mymalloc("Chimes_eps", N_arrayCells * sizeof(double));
          shieldFac_1D = (float *)mymalloc("Chimes_S1d", 3 * sizeof(float));
          shieldFac_2D = (float *)mymalloc("Chimes_S2d", 6 * sizeof(float));

          sprintf(set_name, "/%s/sigmaPhot", chimesRateTables.NonEqIon->ElementName[ns]);
          dataset = H5Dopen(file_id, set_name);

          if(ns == 0 || ns == 1)
            {
              status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, sigmaphot);
              for(i = 0; i < N_arrayCells; i++)
                myPhotoIonTable->tables_per_element[ns].sigmaphot[current_spectrum][i][0] = log10(sigmaphot[i]);
              status = H5Dclose(dataset);
            }
          else
            {
              dataspace_id = H5Dget_space(dataset);
              status       = H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
              dims[0]      = N_arrayCells;
              rank         = 1;
              memspace_id  = H5Screate_simple(rank, dims, NULL);

              for(j = 0; j < chimesRateTables.NonEqIon->N_Auger[ns]; j++)
                {
                  offset2D[0] = 0;
                  offset2D[1] = j;
                  count2D[0]  = N_arrayCells;
                  count2D[1]  = 1;
                  status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D, NULL);

                  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, sigmaphot);

                  for(i = 0; i < N_arrayCells; i++)
                    myPhotoIonTable->tables_per_element[ns].sigmaphot[current_spectrum][i][j] = log10(sigmaphot[i]);
                }
              H5Sclose(memspace_id);
              H5Sclose(dataspace_id);
              status = H5Dclose(dataset);
            }

          sprintf(set_name, "/%s/epsilonPhot", chimesRateTables.NonEqIon->ElementName[ns]);
          dataset = H5Dopen(file_id, set_name);
          status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, epsilon);
          status  = H5Dclose(dataset);

          sprintf(set_name, "/%s/shieldFactor_1D", chimesRateTables.NonEqIon->ElementName[ns]);
          dataset      = H5Dopen(file_id, set_name);
          dataspace_id = H5Dget_space(dataset);
          status       = H5Sget_simple_extent_dims(dataspace_id, dims3D, NULL);
          dims[0]      = 3;
          rank         = 1;
          memspace_id  = H5Screate_simple(rank, dims, NULL);

          for(i = 0; i < N_arrayCells; i++)
            {
              for(j = 0; j < chimesRateTables.NonEqIon->shieldingColumnDimensions[0]; j++)
                {
                  offset3D[0] = 0;
                  offset3D[1] = i;
                  offset3D[2] = j;
                  count3D[0]  = 3;
                  count3D[1]  = 1;
                  count3D[2]  = 1;
                  status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset3D, NULL, count3D, NULL);
                  status      = H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT, shieldFac_1D);

                  for(k = 0; k < 3; k++)
                    myPhotoIonTable->tables_per_element[ns].shieldFactor1D[current_spectrum][k][i][j] =
                        (float)max(shieldFac_1D[k], -300.0);
                }
            }
          H5Sclose(memspace_id);
          H5Sclose(dataspace_id);
          H5Dclose(dataset);

          sprintf(set_name, "/%s/shieldFactor_2D", chimesRateTables.NonEqIon->ElementName[ns]);
          dataset      = H5Dopen(file_id, set_name);
          dataspace_id = H5Dget_space(dataset);
          status       = H5Sget_simple_extent_dims(dataspace_id, dims4D, NULL);
          dims[0]      = 6;
          rank         = 1;
          memspace_id  = H5Screate_simple(rank, dims, NULL);

          for(m = 0; m < N_arrayCells; m++)
            {
              for(i = 0; i < chimesRateTables.NonEqIon->shieldingColumnDimensions[0]; i++)
                {
                  for(j = 0; j < chimesRateTables.NonEqIon->shieldingColumnDimensions[0]; j++)
                    {
                      offset4D[0] = 0;
                      offset4D[1] = m;
                      offset4D[2] = i;
                      offset4D[3] = j;
                      count4D[0]  = 6;
                      count4D[1]  = 1;
                      count4D[2]  = 1;
                      count4D[3]  = 1;
                      status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset4D, NULL, count4D, NULL);
                      status      = H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT, shieldFac_2D);

                      for(k = 0; k < 6; k++)
                        myPhotoIonTable->tables_per_element[ns].shieldFactor2D[current_spectrum][k][m][i][j] =
                            (float)max(shieldFac_2D[k], -300.0);
                    }
                }
            }
          H5Sclose(memspace_id);
          H5Sclose(dataspace_id);
          H5Dclose(dataset);

          for(i = 0; i < N_arrayCells; i++)
            {
              myPhotoIonTable->tables_per_element[ns].epsilon[current_spectrum][i] = log10(epsilon[i]);
            }

          myfree(shieldFac_2D);
          myfree(shieldFac_1D);
          myfree(epsilon);
          myfree(sigmaphot);

          ns += 1;
        }
    }

  /* Read in isotropic_photon_density, H2_dissocJ and dust_G_parameter. */
  sprintf(set_name, "/isotropic_photon_density");
  dataset = H5Dopen(file_id, set_name);
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &n_phot);
  status  = H5Dclose(dataset);

  sprintf(set_name, "/H2_dissocJ");
  dataset = H5Dopen(file_id, set_name);
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &H2_dissocJ);
  status  = H5Dclose(dataset);

  sprintf(set_name, "/dust_G_parameter");
  dataset = H5Dopen(file_id, set_name);
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dust_G_parameter);
  status  = H5Dclose(dataset);

  myPhotoIonTable->isotropic_photon_density[current_spectrum] = log10(n_phot);
  myPhotoIonTable->H2_dissocJ[current_spectrum]               = log10(H2_dissocJ);
  myPhotoIonTable->dust_G_parameter[current_spectrum]         = log10(dust_G_parameter);

  status = H5Fclose(file_id);
}

void ChimesReadPhotoIonTables_UVB_parallel(struct globalVariables *myGlobalVars, char *photoIonTablePath,
                                           struct PhotoIonTables_UVB *myPhotoIonTable, int N_Elements_in_Bens_tables,
                                           int current_spectrum)
{
  /* NOTE: These tables are stored as logs. This makes it quicker
   * when we interpolate the tables, because we don't have to
   * take the log each time. However, the shield factors were
   * already stored as logs. */
  hid_t file_id, dataset, dataspace_id, memspace_id;
  herr_t status;
  hsize_t dims[2];
  hsize_t dims2D[2], count2D[2], offset2D[2];
  hsize_t dims3D[3], count3D[3], offset3D[3];
  hsize_t dims4D[4], count4D[4], offset4D[4];
  int rank, i, j, k, l, m, ns, r, N_arrayCells;
  char set_name[500];
  double *sigmaphot, *epsilon;
  double n_phot, H2_dissocJ, dust_G_parameter;
  float *shieldFac_1D, *shieldFac_2D;
  dims[1] = 1;

  if(ThisTask == 0)
    printf("CHIMES: Reading PhotoIon table: %s \n", photoIonTablePath);

  /* Each element only needs to be read in by 1 task. */
  if(ThisTask < chimesRateTables.NonEqIon->N_Elements)
    file_id = H5Fopen(photoIonTablePath, H5F_ACC_RDONLY, H5P_DEFAULT);

  ns = 0;
  for(l = 0; l < N_Elements_in_Bens_tables; l++)
    {
      if(l < 2 || myGlobalVars->element_included[max(l - 2, 0)] == 1)
        {
          if(ThisTask == ns)
            {
              if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
                N_arrayCells = 3; /* Hydrogen photoionisation also includes H- & H2*/
              else
                N_arrayCells = chimesRateTables.NonEqIon->N_Ions[ns] - 1;

              sigmaphot    = (double *)mymalloc("Chimes_sigma", N_arrayCells * sizeof(double));
              epsilon      = (double *)mymalloc("Chimes_eps", N_arrayCells * sizeof(double));
              shieldFac_1D = (float *)mymalloc("Chimes_S1d", 3 * sizeof(float));
              shieldFac_2D = (float *)mymalloc("Chimes_S2d", 6 * sizeof(float));

              sprintf(set_name, "/%s/sigmaPhot", chimesRateTables.NonEqIon->ElementName[ns]);
              dataset = H5Dopen(file_id, set_name);

              if(ns == 0 || ns == 1)
                {
                  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, sigmaphot);
                  for(i = 0; i < N_arrayCells; i++)
                    myPhotoIonTable->tables_per_element[ns].sigmaphot[current_spectrum][i][0] = log10(sigmaphot[i]);
                  status = H5Dclose(dataset);
                }
              else
                {
                  dataspace_id = H5Dget_space(dataset);
                  status       = H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
                  dims[0]      = N_arrayCells;
                  rank         = 1;
                  memspace_id  = H5Screate_simple(rank, dims, NULL);

                  for(j = 0; j < chimesRateTables.NonEqIon->N_Auger[ns]; j++)
                    {
                      offset2D[0] = 0;
                      offset2D[1] = j;
                      count2D[0]  = N_arrayCells;
                      count2D[1]  = 1;
                      status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D, NULL);

                      status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, sigmaphot);

                      for(i = 0; i < N_arrayCells; i++)
                        myPhotoIonTable->tables_per_element[ns].sigmaphot[current_spectrum][i][j] = log10(sigmaphot[i]);
                    }
                  H5Sclose(memspace_id);
                  H5Sclose(dataspace_id);
                  status = H5Dclose(dataset);
                }

              sprintf(set_name, "/%s/epsilonPhot", chimesRateTables.NonEqIon->ElementName[ns]);
              dataset = H5Dopen(file_id, set_name);
              status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, epsilon);
              status  = H5Dclose(dataset);

              sprintf(set_name, "/%s/shieldFactor_1D", chimesRateTables.NonEqIon->ElementName[ns]);
              dataset      = H5Dopen(file_id, set_name);
              dataspace_id = H5Dget_space(dataset);
              status       = H5Sget_simple_extent_dims(dataspace_id, dims3D, NULL);
              dims[0]      = 3;
              rank         = 1;
              memspace_id  = H5Screate_simple(rank, dims, NULL);

              for(i = 0; i < N_arrayCells; i++)
                {
                  for(j = 0; j < chimesRateTables.NonEqIon->shieldingColumnDimensions[0]; j++)
                    {
                      offset3D[0] = 0;
                      offset3D[1] = i;
                      offset3D[2] = j;
                      count3D[0]  = 3;
                      count3D[1]  = 1;
                      count3D[2]  = 1;
                      status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset3D, NULL, count3D, NULL);
                      status      = H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT, shieldFac_1D);

                      for(k = 0; k < 3; k++)
                        myPhotoIonTable->tables_per_element[ns].shieldFactor1D[current_spectrum][k][i][j] =
                            (float)max(shieldFac_1D[k], -300.0);
                    }
                }
              H5Sclose(memspace_id);
              H5Sclose(dataspace_id);
              H5Dclose(dataset);

              sprintf(set_name, "/%s/shieldFactor_2D", chimesRateTables.NonEqIon->ElementName[ns]);
              dataset      = H5Dopen(file_id, set_name);
              dataspace_id = H5Dget_space(dataset);
              status       = H5Sget_simple_extent_dims(dataspace_id, dims4D, NULL);
              dims[0]      = 6;
              rank         = 1;
              memspace_id  = H5Screate_simple(rank, dims, NULL);

              for(m = 0; m < N_arrayCells; m++)
                {
                  for(i = 0; i < chimesRateTables.NonEqIon->shieldingColumnDimensions[0]; i++)
                    {
                      for(j = 0; j < chimesRateTables.NonEqIon->shieldingColumnDimensions[0]; j++)
                        {
                          offset4D[0] = 0;
                          offset4D[1] = m;
                          offset4D[2] = i;
                          offset4D[3] = j;
                          count4D[0]  = 6;
                          count4D[1]  = 1;
                          count4D[2]  = 1;
                          count4D[3]  = 1;
                          status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset4D, NULL, count4D, NULL);
                          status      = H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT, shieldFac_2D);

                          for(k = 0; k < 6; k++)
                            myPhotoIonTable->tables_per_element[ns].shieldFactor2D[current_spectrum][k][m][i][j] =
                                (float)max(shieldFac_2D[k], -300.0);
                        }
                    }
                }
              H5Sclose(memspace_id);
              H5Sclose(dataspace_id);
              H5Dclose(dataset);

              for(i = 0; i < N_arrayCells; i++)
                {
                  myPhotoIonTable->tables_per_element[ns].epsilon[current_spectrum][i] = log10(epsilon[i]);
                }

              myfree(shieldFac_2D);
              myfree(shieldFac_1D);
              myfree(epsilon);
              myfree(sigmaphot);

              ns += 1;
            }
          else
            ns += 1;
        }
    }

  /* Now broadcast each element to all tasks. */
  ns = 0;
  for(l = 0; l < N_Elements_in_Bens_tables; l++)
    {
      if(l < 2 || myGlobalVars->element_included[max(l - 2, 0)] == 1)
        {
          if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
            N_arrayCells = 3; /* Hydrogen photoionisation also includes H- & H2*/
          else
            N_arrayCells = chimesRateTables.NonEqIon->N_Ions[ns] - 1;

          /* For multi-D arrays, use a 1D buffer. */
          epsilon = (double *)mymalloc("Chimes_eps", N_arrayCells * sizeof(double));
          if(ThisTask == ns)
            {
              for(i = 0; i < N_arrayCells; i++)
                epsilon[i] = myPhotoIonTable->tables_per_element[ns].epsilon[current_spectrum][i];
            }

          MPI_Bcast(epsilon, N_arrayCells, MPI_DOUBLE, ns, MPI_COMM_WORLD);

          if(ThisTask != ns)
            {
              for(i = 0; i < N_arrayCells; i++)
                myPhotoIonTable->tables_per_element[ns].epsilon[current_spectrum][i] = epsilon[i];
            }

          if(ns == 0 || ns == 1)
            sigmaphot = (double *)mymalloc("Chimes_sigma", N_arrayCells * sizeof(double));
          else
            sigmaphot = (double *)mymalloc("Chimes_sigma", N_arrayCells * chimesRateTables.NonEqIon->N_Auger[ns] * sizeof(double));

          shieldFac_1D = (float *)mymalloc("Chimes_S1d",
                                           3 * N_arrayCells * chimesRateTables.NonEqIon->shieldingColumnDimensions[0] * sizeof(float));
          shieldFac_2D = (float *)mymalloc("Chimes_S2d", 6 * N_arrayCells * chimesRateTables.NonEqIon->shieldingColumnDimensions[0] *
                                                             chimesRateTables.NonEqIon->shieldingColumnDimensions[0] * sizeof(float));

          if(ThisTask == ns)
            {
              if(ns == 0 || ns == 1)
                {
                  for(i = 0; i < N_arrayCells; i++)
                    sigmaphot[i] = myPhotoIonTable->tables_per_element[ns].sigmaphot[current_spectrum][i][0];
                }
              else
                {
                  for(i = 0; i < N_arrayCells; i++)
                    {
                      for(j = 0; j < chimesRateTables.NonEqIon->N_Auger[ns]; j++)
                        sigmaphot[(i * chimesRateTables.NonEqIon->N_Auger[ns]) + j] =
                            myPhotoIonTable->tables_per_element[ns].sigmaphot[current_spectrum][i][j];
                    }
                }

              for(i = 0; i < 3; i++)
                {
                  for(j = 0; j < N_arrayCells; j++)
                    {
                      for(k = 0; k < chimesRateTables.NonEqIon->shieldingColumnDimensions[0]; k++)
                        shieldFac_1D[(i * chimesRateTables.NonEqIon->shieldingColumnDimensions[0] * N_arrayCells) +
                                     (j * chimesRateTables.NonEqIon->shieldingColumnDimensions[0]) + k] =
                            myPhotoIonTable->tables_per_element[ns].shieldFactor1D[current_spectrum][i][j][k];
                    }
                }

              for(i = 0; i < 6; i++)
                {
                  for(j = 0; j < N_arrayCells; j++)
                    {
                      for(k = 0; k < chimesRateTables.NonEqIon->shieldingColumnDimensions[0]; k++)
                        {
                          for(m = 0; m < chimesRateTables.NonEqIon->shieldingColumnDimensions[0]; m++)
                            shieldFac_2D[(i * chimesRateTables.NonEqIon->shieldingColumnDimensions[0] *
                                          chimesRateTables.NonEqIon->shieldingColumnDimensions[0] * N_arrayCells) +
                                         (j * chimesRateTables.NonEqIon->shieldingColumnDimensions[0] *
                                          chimesRateTables.NonEqIon->shieldingColumnDimensions[0]) +
                                         (k * chimesRateTables.NonEqIon->shieldingColumnDimensions[0]) + m] =
                                myPhotoIonTable->tables_per_element[ns].shieldFactor2D[current_spectrum][i][j][k][m];
                        }
                    }
                }
            }

          if(ns == 0 || ns == 1)
            MPI_Bcast(sigmaphot, N_arrayCells, MPI_DOUBLE, ns, MPI_COMM_WORLD);
          else
            MPI_Bcast(sigmaphot, N_arrayCells * chimesRateTables.NonEqIon->N_Auger[ns], MPI_DOUBLE, ns, MPI_COMM_WORLD);

          MPI_Bcast(shieldFac_1D, 3 * N_arrayCells * chimesRateTables.NonEqIon->shieldingColumnDimensions[0], MPI_FLOAT, ns,
                    MPI_COMM_WORLD);
          MPI_Bcast(shieldFac_2D,
                    6 * N_arrayCells * chimesRateTables.NonEqIon->shieldingColumnDimensions[0] *
                        chimesRateTables.NonEqIon->shieldingColumnDimensions[0],
                    MPI_FLOAT, ns, MPI_COMM_WORLD);

          MPI_Barrier(MPI_COMM_WORLD);

          if(ThisTask != ns)
            {
              if(ns == 0 || ns == 1)
                {
                  for(i = 0; i < N_arrayCells; i++)
                    myPhotoIonTable->tables_per_element[ns].sigmaphot[current_spectrum][i][0] = sigmaphot[i];
                }
              else
                {
                  for(i = 0; i < N_arrayCells; i++)
                    {
                      for(j = 0; j < chimesRateTables.NonEqIon->N_Auger[ns]; j++)
                        myPhotoIonTable->tables_per_element[ns].sigmaphot[current_spectrum][i][j] =
                            sigmaphot[(i * chimesRateTables.NonEqIon->N_Auger[ns]) + j];
                    }
                }

              for(i = 0; i < 3; i++)
                {
                  for(j = 0; j < N_arrayCells; j++)
                    {
                      for(k = 0; k < chimesRateTables.NonEqIon->shieldingColumnDimensions[0]; k++)
                        myPhotoIonTable->tables_per_element[ns].shieldFactor1D[current_spectrum][i][j][k] =
                            shieldFac_1D[(i * chimesRateTables.NonEqIon->shieldingColumnDimensions[0] * N_arrayCells) +
                                         (j * chimesRateTables.NonEqIon->shieldingColumnDimensions[0]) + k];
                    }
                }

              for(i = 0; i < 6; i++)
                {
                  for(j = 0; j < N_arrayCells; j++)
                    {
                      for(k = 0; k < chimesRateTables.NonEqIon->shieldingColumnDimensions[0]; k++)
                        {
                          for(m = 0; m < chimesRateTables.NonEqIon->shieldingColumnDimensions[0]; m++)
                            myPhotoIonTable->tables_per_element[ns].shieldFactor2D[current_spectrum][i][j][k][m] =
                                shieldFac_2D[(i * chimesRateTables.NonEqIon->shieldingColumnDimensions[0] *
                                              chimesRateTables.NonEqIon->shieldingColumnDimensions[0] * N_arrayCells) +
                                             (j * chimesRateTables.NonEqIon->shieldingColumnDimensions[0] *
                                              chimesRateTables.NonEqIon->shieldingColumnDimensions[0]) +
                                             (k * chimesRateTables.NonEqIon->shieldingColumnDimensions[0]) + m];
                        }
                    }
                }
            }

          MPI_Barrier(MPI_COMM_WORLD);

          myfree(shieldFac_2D);
          myfree(shieldFac_1D);
          myfree(sigmaphot);
          myfree(epsilon);

          ns += 1;
        }
    }

  if(ThisTask == 0)
    {
      /* Read in isotropic_photon_density, H2_dissocJ and dust_G_parameter. */
      sprintf(set_name, "/isotropic_photon_density");
      dataset = H5Dopen(file_id, set_name);
      status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &n_phot);
      status  = H5Dclose(dataset);

      sprintf(set_name, "/H2_dissocJ");
      dataset = H5Dopen(file_id, set_name);
      status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &H2_dissocJ);
      status  = H5Dclose(dataset);

      sprintf(set_name, "/dust_G_parameter");
      dataset = H5Dopen(file_id, set_name);
      status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dust_G_parameter);
      status  = H5Dclose(dataset);
    }

  MPI_Bcast(&n_phot, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&H2_dissocJ, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&dust_G_parameter, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  myPhotoIonTable->isotropic_photon_density[current_spectrum] = log10(n_phot);
  myPhotoIonTable->H2_dissocJ[current_spectrum]               = log10(H2_dissocJ);
  myPhotoIonTable->dust_G_parameter[current_spectrum]         = log10(dust_G_parameter);

  if(ThisTask < chimesRateTables.NonEqIon->N_Elements)
    status = H5Fclose(file_id);
}

void ChimesCopyPhotoIonTables(struct PhotoIonTables_UVB *myPhotoIonTable, int i, int j)
{
  /* Copy table i to table j. */
  int ns, k, l, m, r, N_species;

  for(ns = 0; ns < chimesRateTables.NonEqIon->N_Elements; ns++)
    {
      if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
        N_species = 3;
      else
        N_species = chimesRateTables.NonEqIon->N_Ions[ns] - 1;

      for(k = 0; k < N_species; k++)
        {
          for(l = 0; l < chimesRateTables.NonEqIon->N_Auger[ns]; l++)
            myPhotoIonTable->tables_per_element[ns].sigmaphot[j][k][l] = myPhotoIonTable->tables_per_element[ns].sigmaphot[i][k][l];
        }
    }

  for(ns = 0; ns < chimesRateTables.NonEqIon->N_Elements; ns++)
    {
      if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
        N_species = 3;
      else
        N_species = chimesRateTables.NonEqIon->N_Ions[ns] - 1;

      for(k = 0; k < 3; k++)
        {
          for(l = 0; l < N_species; l++)
            {
              for(m = 0; m < chimesRateTables.NonEqIon->shieldingColumnDimensions[0]; m++)
                myPhotoIonTable->tables_per_element[ns].shieldFactor1D[j][k][l][m] =
                    myPhotoIonTable->tables_per_element[ns].shieldFactor1D[i][k][l][m];
            }
        }
    }

  for(ns = 0; ns < chimesRateTables.NonEqIon->N_Elements; ns++)
    {
      if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
        N_species = 3;
      else
        N_species = chimesRateTables.NonEqIon->N_Ions[ns] - 1;

      for(k = 0; k < 6; k++)
        {
          for(l = 0; l < N_species; l++)
            {
              for(m = 0; m < chimesRateTables.NonEqIon->shieldingColumnDimensions[0]; m++)
                {
                  for(r = 0; r < chimesRateTables.NonEqIon->shieldingColumnDimensions[0]; r++)
                    myPhotoIonTable->tables_per_element[ns].shieldFactor2D[j][k][l][m][r] =
                        myPhotoIonTable->tables_per_element[ns].shieldFactor2D[i][k][l][m][r];
                }
            }
        }
    }

  for(ns = 0; ns < chimesRateTables.NonEqIon->N_Elements; ns++)
    {
      if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
        N_species = 3;
      else
        N_species = chimesRateTables.NonEqIon->N_Ions[ns] - 1;

      for(k = 0; k < N_species; k++)
        myPhotoIonTable->tables_per_element[ns].epsilon[j][k] = myPhotoIonTable->tables_per_element[ns].epsilon[i][k];
    }

  myPhotoIonTable->dust_G_parameter[j]         = myPhotoIonTable->dust_G_parameter[i];
  myPhotoIonTable->H2_dissocJ[j]               = myPhotoIonTable->H2_dissocJ[i];
  myPhotoIonTable->isotropic_photon_density[j] = myPhotoIonTable->isotropic_photon_density[i];
}

void ChimesInterpolatePhotoIonTables(struct globalVariables *myGlobalVars, struct PhotoIonTables_UVB *myPhotoIonTable, int z_index,
                                     double dz)
{
  /* Recall that the PhotoIon tables were stored as logs.
   * However, the shield factor arrays can be kept as
   * logs, as that is how CHIMES uses them. */
  int ns, k, l, m, r, N_species;
  double log_x;
  static int first_call = 0;

  if(z_index >= Chimes_N_Redshifts - 1)
    {
      /* Before highest redshift table. Take this to be
       * before reionisation, and set UVB to zero. Then set
       * cross sections etc. to one of the high-z table.
       * Only needs to be done once. */
      if(first_call == 0)
        {
          first_call = 1;

          for(ns = 0; ns < chimesRateTables.NonEqIon->N_Elements; ns++)
            {
              if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
                N_species = 3;
              else
                N_species = chimesRateTables.NonEqIon->N_Ions[ns] - 1;

              for(k = 0; k < N_species; k++)
                {
                  for(l = 0; l < chimesRateTables.NonEqIon->N_Auger[ns]; l++)
                    chimesRateTables.NonEqIon->NonEqRates[ns].sigmaphot[0][k][l] =
                        pow(10.0, myPhotoIonTable->tables_per_element[ns].sigmaphot[1][k][l]);
                }
            }

          for(ns = 0; ns < chimesRateTables.NonEqIon->N_Elements; ns++)
            {
              if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
                N_species = 3;
              else
                N_species = chimesRateTables.NonEqIon->N_Ions[ns] - 1;

              for(k = 0; k < 3; k++)
                {
                  for(l = 0; l < N_species; l++)
                    {
                      for(m = 0; m < chimesRateTables.NonEqIon->shieldingColumnDimensions[0]; m++)
                        chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor1D[0][k][l][m] =
                            myPhotoIonTable->tables_per_element[ns].shieldFactor1D[1][k][l][m];
                    }
                }
            }

          for(ns = 0; ns < chimesRateTables.NonEqIon->N_Elements; ns++)
            {
              if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
                N_species = 3;
              else
                N_species = chimesRateTables.NonEqIon->N_Ions[ns] - 1;

              for(k = 0; k < 6; k++)
                {
                  for(l = 0; l < N_species; l++)
                    {
                      for(m = 0; m < chimesRateTables.NonEqIon->shieldingColumnDimensions[0]; m++)
                        {
                          for(r = 0; r < chimesRateTables.NonEqIon->shieldingColumnDimensions[0]; r++)
                            chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[0][k][l][m][r] =
                                myPhotoIonTable->tables_per_element[ns].shieldFactor2D[1][k][l][m][r];
                        }
                    }
                }
            }

          for(ns = 0; ns < chimesRateTables.NonEqIon->N_Elements; ns++)
            {
              if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
                N_species = 3;
              else
                N_species = chimesRateTables.NonEqIon->N_Ions[ns] - 1;

              for(k = 0; k < N_species; k++)
                chimesRateTables.NonEqIon->NonEqRates[ns].epsilon[0][k] =
                    pow(10.0, myPhotoIonTable->tables_per_element[ns].epsilon[1][k]);
            }

          ChimesDustGArr[0]                = pow(10.0, myPhotoIonTable->dust_G_parameter[1]);
          ChimesH2DissocJArr[0]            = pow(10.0, myPhotoIonTable->H2_dissocJ[1]);
          All.ChimesIsotropicPhotonDensity = 0.0; /* Before reionisation. */
        }
    }
  else
    {
      for(ns = 0; ns < chimesRateTables.NonEqIon->N_Elements; ns++)
        {
          if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
            N_species = 3;
          else
            N_species = chimesRateTables.NonEqIon->N_Ions[ns] - 1;

          for(k = 0; k < N_species; k++)
            {
              for(l = 0; l < chimesRateTables.NonEqIon->N_Auger[ns]; l++)
                {
                  log_x = (myPhotoIonTable->tables_per_element[ns].sigmaphot[1][k][l] * dz) +
                          (myPhotoIonTable->tables_per_element[ns].sigmaphot[0][k][l] * (1.0 - dz));
                  chimesRateTables.NonEqIon->NonEqRates[ns].sigmaphot[0][k][l] = pow(10.0, log_x);
                }
            }
        }

      for(ns = 0; ns < chimesRateTables.NonEqIon->N_Elements; ns++)
        {
          if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
            N_species = 3;
          else
            N_species = chimesRateTables.NonEqIon->N_Ions[ns] - 1;

          for(k = 0; k < 3; k++)
            {
              for(l = 0; l < N_species; l++)
                {
                  for(m = 0; m < chimesRateTables.NonEqIon->shieldingColumnDimensions[0]; m++)
                    {
                      log_x = (myPhotoIonTable->tables_per_element[ns].shieldFactor1D[1][k][l][m] * dz) +
                              (myPhotoIonTable->tables_per_element[ns].shieldFactor1D[0][k][l][m] * (1.0 - dz));
                      chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor1D[0][k][l][m] = log_x;
                    }
                }
            }
        }

      for(ns = 0; ns < chimesRateTables.NonEqIon->N_Elements; ns++)
        {
          if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
            N_species = 3;
          else
            N_species = chimesRateTables.NonEqIon->N_Ions[ns] - 1;

          for(k = 0; k < 6; k++)
            {
              for(l = 0; l < N_species; l++)
                {
                  for(m = 0; m < chimesRateTables.NonEqIon->shieldingColumnDimensions[0]; m++)
                    {
                      for(r = 0; r < chimesRateTables.NonEqIon->shieldingColumnDimensions[0]; r++)
                        {
                          log_x = (myPhotoIonTable->tables_per_element[ns].shieldFactor2D[1][k][l][m][r] * dz) +
                                  (myPhotoIonTable->tables_per_element[ns].shieldFactor2D[0][k][l][m][r] * (1.0 - dz));
                          chimesRateTables.NonEqIon->NonEqRates[ns].shieldFactor2D[0][k][l][m][r] = log_x;
                        }
                    }
                }
            }
        }

      for(ns = 0; ns < chimesRateTables.NonEqIon->N_Elements; ns++)
        {
          if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
            N_species = 3;
          else
            N_species = chimesRateTables.NonEqIon->N_Ions[ns] - 1;

          for(k = 0; k < N_species; k++)
            {
              log_x = (myPhotoIonTable->tables_per_element[ns].epsilon[1][k] * dz) +
                      (myPhotoIonTable->tables_per_element[ns].epsilon[0][k] * (1.0 - dz));
              chimesRateTables.NonEqIon->NonEqRates[ns].epsilon[0][k] = pow(10.0, log_x);
            }
        }

      log_x             = (myPhotoIonTable->dust_G_parameter[1] * dz) + (myPhotoIonTable->dust_G_parameter[0] * (1.0 - dz));
      ChimesDustGArr[0] = pow(10.0, log_x);

      log_x                 = (myPhotoIonTable->H2_dissocJ[1] * dz) + (myPhotoIonTable->H2_dissocJ[0] * (1.0 - dz));
      ChimesH2DissocJArr[0] = pow(10.0, log_x);

      log_x = (myPhotoIonTable->isotropic_photon_density[1] * dz) + (myPhotoIonTable->isotropic_photon_density[0] * (1.0 - dz));
      All.ChimesIsotropicPhotonDensity = pow(10.0, log_x);
    }
}
#endif /* CHIMES_REDSHIFT_DEPENDENT_UVB */

/* The following routine sets the pointers to the dynamcally
 * allocated arrays within ChimesGasVars to the correct
 * position within the corresponding global arrays. */
void chimes_set_pointers(int i)
{
  /* Check that the index stored within ChimesGasVars
   * matches the particle index. If not, the particle
   * may have been moved within SphP[] without correctly
   * moving its corresponding abundance array within
   * the global ChimesAbundances[] array. */
  if(SphP[i].ChimesGasVars.index != i)
    terminate("i=%d, SphP[i].ChimesGasVars.index = %d \n", i, SphP[i].ChimesGasVars.index);

  SphP[i].ChimesGasVars.abundances               = &(ChimesAbundances[ChimesGlobalVars.totalNumberOfSpecies * i]);
  SphP[i].ChimesGasVars.isotropic_photon_density = &(ChimesPhotonDensity[ChimesGlobalVars.N_spectra * i]);
  SphP[i].ChimesGasVars.dust_G_parameter         = &(ChimesDustG[ChimesGlobalVars.N_spectra * i]);
  SphP[i].ChimesGasVars.H2_dissocJ               = &(ChimesH2dissocJ[ChimesGlobalVars.N_spectra * i]);
#ifdef CHIMES_ADVECT_ABUNDANCES
  SphP[i].ChimesGasVars.IonAdvect = &(ChimesIonAdvect[ChimesGlobalVars.totalNumberOfSpecies * i]);
#endif
}

/* The following routine loops through all gas particles
 * and updates their various pointers within ChimesGasVars
 * to the correct positions within the corresponding global
 * arrays. This is needed in case the addresses of the
 * global arrays have changed (e.g. if memory blocks
 * have been moved around etc. */
void chimes_update_all_pointers(void)
{
  int i;

  for(i = 0; i < NumGas; i++)
    chimes_set_pointers(i);
}

/* The following routine moves the abundances of
 * particle j to position i, and updates the
 * various pointers in the ChimesGasVars of i
 * accordingly. */
void chimes_move_abundances(int i, int j)
{
  int k;
  for(k = 0; k < ChimesGlobalVars.totalNumberOfSpecies; k++)
    ChimesAbundances[(ChimesGlobalVars.totalNumberOfSpecies * i) + k] =
        ChimesAbundances[(ChimesGlobalVars.totalNumberOfSpecies * j) + k];

  chimes_set_pointers(i);
}

/* The following routine sets the pointers to the dynamcally
 * allocated arrays within ChimesGasVars to NULL. This is used
 * when a gas particle is deleted, and ensures that we don't
 * try to accidentally access a part of the global abundances
 * array that is no longer in use. */
void chimes_set_pointers_to_null(int i)
{
  SphP[i].ChimesGasVars.abundances               = NULL;
  SphP[i].ChimesGasVars.isotropic_photon_density = NULL;
  SphP[i].ChimesGasVars.dust_G_parameter         = NULL;
  SphP[i].ChimesGasVars.H2_dissocJ               = NULL;
#ifdef CHIMES_ADVECT_ABUNDANCES
  SphP[i].ChimesGasVars.IonAdvect = NULL;
#endif
}

#endif /* CHIMES */
