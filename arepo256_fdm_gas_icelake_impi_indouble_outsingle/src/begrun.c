/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/begrun.c
 * \date        05/2018
 * \author      Volker Springel
 * \brief       Initial set-up of a simulation run
 * \details     This file contains various functions to initialize a simulation
 *              run. In particular, the parameter file is read in and parsed
 *              and global variables are initialized to their proper values.
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 03.05.2018 Function documentation and indent -- Rainer Weinberger
 */

#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "MRT/RT.h"
#include "allvars.h"
#include "compiler-command-line-args.h"
#include "domain.h"
#include "gitversion/version.h"
#include "proto.h"
#include "voronoi.h"

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

static void delete_end_file(void);

#ifdef HIGH_FREQUENCY_OUTPUT_STARS
static void init_high_freq_star_outputs(void);
#endif

/*! \brief Prints a welcome message.
 *
 *  \return void
 */
void hello(void)
{
#ifndef DG
  mpi_printf(
      "\n   __    ____  ____  ____  _____\n  /__\\  (  _ \\( ___)(  _ \\(  _  )\n /(__)\\  )   / )__)  )___/ "
      ")(_)(\n(__)(__)(_)\\_)(____)(__)  (_____)\n\n");
#else  /* #ifndef DG */
  mpi_printf(
      "  _____   U _____ u   _   _     U _____ u   _____      \n |_ \" _|  \\| ___\"|/  | \\ |\"|    \\| ___\"|/  |_ \" _|     \n"
      "   | |     |  _|\"   <|  \\| |>    |  _|\"      | |       \n  /| |\\    | |___   U| |\\  |u    | |___     /| |\\      \n"
      " u |_|U    |_____|   |_| \\_|     |_____|   u |_|U      \n _//\\\\_    <<   >>   ||   \\\\,-.  <<   >>   _// \\\\_     \n"
      "(__) (__) (__) (__)  (_\")  (_/  (__) (__) (__) (__)    \n");
#endif /* #ifndef DG #else */
}

/*! \brief Prints used compile options.
 *
 *  \return void
 */
void begrun0(void)
{
#ifdef DG
  const char *program_name    = "Tenet";
  const char *program_version = TENET_VERSION;
#else
  const char *program_name    = "Arepo";
  const char *program_version = AREPO_VERSION;
#endif

  mpi_printf(
      "\nThis is %s, version %s (git: %s).\n\nRunning with %d MPI tasks.\n\nApparently we're using %d compute nodes (we have a "
      "minimum of %d MPI tasks per node, and a maximum of %d)\n\n",
      program_name, program_version, GIT_COMMIT, NTask, NumNodes, MinTasksPerNode, MaxTasksPerNode);

  mpi_printf("Code was compiled with the following compiler and flags:\n%s\n\n\n", compiler_flags);

  mpi_printf("Code was compiled with settings:\n\n");

  if(ThisTask == 0)
    output_compile_time_options();
}

/*! \brief Initial setup of the simulation.
 *
 *  First, the parameter file is read by read_parameter_file(),
 *  then routines for setting units, etc are called. This function only does
 *  the setup necessary to load the IC file. After the IC file has been loaded
 *  and prepared by init(), setup continues with begrun2(). This splitting is
 *  done so that we can return cleanly from operations that don't actually
 *  start the simulation (converting snapshots, making projected images, etc.)
 *
 * \return void
 */
void begrun1(void)
{
#if defined(PMGRID) && !defined(ONLY_PM) && defined(ASMTH)
  if(ASMTH == 0)
    mpi_terminate("ASMTH == 0 cannot be used with TreePM");
#endif

  /* ... read in parameters for this run */
  read_parameter_file(ParameterFile);

  /* consistency check of parameters */
  check_parameters();

#ifndef TGSET
  healthtest();
#endif

#ifdef HAVE_HDF5
  H5Eset_auto(my_hdf5_error_handler, NULL);
#endif

  gsl_set_error_handler(my_gsl_error_handler);

#ifdef CUDA
  cuda_init();
#endif

#if defined(DM_WINDTUNNEL) && defined(DM_WINDTUNNEL_EXTERNAL_SOURCE)
  read_dmwindtunnel_file();
#endif

#if defined(WINDTUNNEL) && defined(WINDTUNNEL_EXTERNAL_SOURCE)
  read_windtunnel_file();
#endif

#if defined(WINDTUNNEL) && defined(WINDTUNNEL_READ_IN_BFIELD)
  WindtunnelReadIn_InitialiseGlobals();
#endif

#ifdef DEBUG
  enable_core_dumps_and_fpu_exceptions();
#endif

  mpi_printf("BEGRUN: Size of particle structure       %3d  [bytes]\n", (int)sizeof(struct particle_data));
  mpi_printf("BEGRUN: Size of SPH particle structure   %3d  [bytes]\n", (int)sizeof(struct sph_particle_data));
  mpi_printf("BEGRUN: Size of gravity tree node        %3d  [bytes]\n", (int)sizeof(struct NODE));
#ifdef MULTIPLE_NODE_SOFTENING
  mpi_printf("BEGRUN: Size of auxiliary gravity node   %3d  [bytes]\n", (int)sizeof(struct ExtNODE));
#endif
#if defined(GFM) || defined(SFR_MCS)
  mpi_printf("BEGRUN: Size of star particle structure  %3d  [bytes]\n", (int)sizeof(struct star_particle_data));
#endif
#ifdef BLACK_HOLES
  mpi_printf("BEGRUN: Size of BH particle structure    %3d  [bytes]\n\n", (int)sizeof(struct bh_particle_data));
#endif
#ifdef SINKS
  mpi_printf("BEGRUN: Size of SinkP particle structure %3d  [bytes]\n\n", (int)sizeof(struct sink_particle_data));
#endif
#ifdef DUST_LIVE
  mpi_printf("BEGRUN: Size of dust particle structure  %3d  [bytes]\n\n", (int)sizeof(struct dust_particle_data));
#endif

#if defined(DARKENERGY) && defined(TIMEDEPDE)
  /* set up table needed for hubble functions with time dependent w */
  fwa_init();
#endif

  set_units();

#ifdef MRT
  init_RT();

#if(defined(MRT_SOURCES) && defined(MRT_STARS)) || (defined(MRT_LOCAL_FEEDBACK))
  mrt_load_spectrum_table();
#endif /* #if defined(MRT_SOURCES) && defined(MRT_STARS) */

#ifndef MRT_NO_UV
  mrt_get_sigma();
  mrt_get_luminosities();
#endif /* #ifndef MRT_NO_UV */

#if(defined(MRT_SOURCES) && defined(MRT_STARS)) || (defined(MRT_LOCAL_FEEDBACK))
  mrt_free_spectrum_table();
#endif /* #if defined(MRT_SOURCES) && defined(MRT_STARS) */
#endif /* #ifdef MRT */

  /* this is needed here to allow domain decomposition right after restart */
  if(RestartFlag == RESTART_RESTART)
    if(All.ComovingIntegrationOn)
      init_drift_table();
#ifndef CHIMES
  /* When we include CHIMES, we will need to run
   * init_io_fields() after ChimesInitCool(). */
  init_io_fields();
#endif

  force_short_range_init();

#if defined(FORCETEST) && !defined(FORCETEST_TESTFORCELAW) && !defined(GRAVITY_TALLBOX)
  forcetest_ewald_init();
#endif

  /* set up random number generators */
  random_generator     = gsl_rng_alloc(gsl_rng_ranlxd1);
  random_generator_aux = gsl_rng_alloc(gsl_rng_ranlxd1);

  /* individual start-up seed */
#ifdef RANDOM_REALIZATION
  gsl_rng_set(random_generator, 42 + ThisTask + RANDOM_REALIZATION);
  gsl_rng_set(random_generator_aux, 31452 + ThisTask + RANDOM_REALIZATION);
#else
  gsl_rng_set(random_generator, 42 + ThisTask);
  gsl_rng_set(random_generator_aux, 31452 + ThisTask);
#endif

#ifdef SNE_FEEDBACK
  sne_rng = gsl_rng_alloc(gsl_rng_ranlxd2);
#endif

#ifdef PERFORMANCE_TEST_SPARSE_MPI_ALLTOALL
  test_mpi_alltoall_performance();
#endif

  timebins_init(&TimeBinsHydro, "Hydro", &All.MaxPartSph);
  timebins_init(&TimeBinsGravity, "Gravity", &All.MaxPart);

#ifdef TRACER_PARTICLE
  /* assume that there are less tracers than cells */
  timebins_init(&TimeBinsTracer, "Tracer", &All.MaxPart);
#endif

#ifdef BLACK_HOLES
  timebins_init(&TimeBinsBHAccretion, "BHAccretion", &All.MaxPart);
#endif

#ifdef SINKS
  timebins_init(&TimeBinsSinksAccretion, "SinksAccretion", &All.MaxPart);
#endif

#ifdef DUST_LIVE
  timebins_init(&TimeBinsDust, "Dust", &All.MaxPart);
#endif

#ifdef SUBBOX_SNAPSHOTS
  read_subbox_coordinates(All.SubboxCoordinatesPath);
#endif

#ifdef MRT_METAL_COOLING
  init_cooling_metal();
  mpi_printf("GFM_COOLING_METAL: Metal line cooling rates initialized.\n");
#endif

#ifdef GFM_STELLAR_EVOLUTION
#ifndef GFM_STELLAR_EVOLUTION_NO_ELEMENTS
  int mc = 0;
  /* copy element names to ElementNames */
  strcpy(ElementNames[mc++], "Hydrogen");
  strcpy(ElementNames[mc++], "Helium");
  strcpy(ElementNames[mc++], "Carbon");
  strcpy(ElementNames[mc++], "Nitrogen");
  strcpy(ElementNames[mc++], "Oxygen");
  strcpy(ElementNames[mc++], "Neon");
  strcpy(ElementNames[mc++], "Magnesium");
  strcpy(ElementNames[mc++], "Silicon");
  strcpy(ElementNames[mc++], "Iron");

#ifdef GFM_SPROCESS
  strcpy(ElementNames[mc++], "Yttrium");
  strcpy(ElementNames[mc++], "Strontium");
  strcpy(ElementNames[mc++], "Zirconium");
  strcpy(ElementNames[mc++], "Barium");
  strcpy(ElementNames[mc++], "Lead");
#endif /* #ifdef GFM_SPROCESS */

#ifdef GFM_NORMALIZED_METAL_ADVECTION
  strcpy(ElementNames[mc++], "OtherMetals");
#endif /* #ifdef GFM_NORMALIZED_METAL_ADVECTION */
#endif

  init_imf();
  mpi_printf("GFM_STELLAR_EVOLUTION: IMF initialized.\n");

  init_yields();
  mpi_printf("GFM_STELLAR_EVOLUTION: Yields initialized.\n");

#ifdef RADCOOL
  init_radcool_units();
  mpi_printf("RADCOOL, RADCOOL_HOTHALO : Mass and length units converted to Msun and kpc respectively.\n");
#endif /* #ifdef RADCOOL */

#ifdef GFM_COOLING_METAL
  init_cooling_metal();
  mpi_printf("GFM_COOLING_METAL: Metal line cooling rates initialized.\n");
#endif /* #ifdef GFM_COOLING_METAL */

#ifdef GFM_WINDS
  init_winds();
#endif /* #ifdef GFM_WINDS */

#ifdef GFM_WINDS_VARIABLE
  if(ThisTask == 0)
    init_variable_winds();
#endif /* #ifdef GFM_WINDS_VARIABLE */

#ifdef GFM_PREENRICH
  gfm_read_preenrich_table(All.PreEnrichAbundanceFile);
#endif /* #ifdef GFM_PREENRICH */

#if defined(WINDTUNNEL_FIXVARIABLESININJECTIONREGION) && defined(GFM_SET_METALLICITY)
  get_initial_mass_fractions(&All.mass_fractions[0], All.GasMetallicityInSolar);
#endif

#endif /* #ifdef GFM_STELLAR_EVOLUTION */

#if defined(GFM_COOLING_METAL) && defined(SFR_MCS)
  /* Enables use of GFM_COOLING_METAL without rest of GFM */
  init_cooling_metal();
  mpi_printf("GFM_COOLING_METAL: Metal line cooling rates initialized.\n");
#endif

#ifdef GFM_STELLAR_PHOTOMETRICS
  init_stellar_photometrics();
  mpi_printf("GFM_STELLAR_PHOTOMETRICS: Stellar photometrics initialized.\n");
#endif

#if defined(BLACK_HOLES) && defined(GFM_AGN_RADIATION)
  init_agn_radiation();
  mpi_printf("GFM_AGN_RADIATION: initialized.\n");
#endif

#if defined(TRACER_PART_NUM_FLUID_QUANTITIES) || defined(TRACER_MC_NUM_FLUID_QUANTITIES)
  set_tracer_part_indices();
#endif

#ifdef GROWING_DISK_POTENTIAL
  growing_disk_init();
#endif

#if defined(COOLING) & !defined(SIMPLE_COOLING) & !defined(GRACKLE)
  All.Time = All.TimeBegin;
  set_cosmo_factors_for_current_time();
#ifdef CHIMES
  ChimesInitCool();
  init_io_fields();
#else
  InitCool();
#endif /* CHIMES */
#endif /* #if defined(COOLING) & !defined(SIMPLE_COOLING) & !defined(GRACKLE) */

#if defined(COOLING) && defined(GRACKLE)
  All.Time = All.TimeBegin;
  set_cosmo_factors_for_current_time();
  initialise_grackle();
#endif

#ifdef ATOMIC_DM
  All.Time = All.TimeBegin;
  set_cosmo_factors_for_current_time();
  ADM_InitCool();
#endif

#ifdef TGSET
  tgset_begrun();
#endif

#ifdef HEALRAY
  /* Needs to be before TGCHEM! */
  healray_begrun();
#endif

#ifdef TGCHEM
  tgchem_begrun();
#endif

#ifdef SGCHEM
  init_chemistry();
#endif

#ifdef TREECOLV2
  int success_read_treecolv2_tables;
  success_read_treecolv2_tables = treecolv2_read_and_setup_lookup_table();
  if(success_read_treecolv2_tables == 0)
    terminate("Problem reading in TreeCol lookup tables");
#endif

#ifdef SINKS
  sinks_begrun();
#endif

#ifdef SINK_PARTICLES
  init_sink_particles(0);
#endif

#if !defined(PMGRID) && defined(SELFGRAVITY) && !defined(GRAVITY_NOT_PERIODIC) && !defined(ONEDIMS_SPHERICAL)
  ewald_init();
#endif

#ifdef TILE_ICS
  All.BoxSize *= All.TileICsFactor;
#endif

  for(size_t i = 0; i < sizeof(All.BoxSizes) / sizeof(All.BoxSizes[0]); i++)
    All.BoxSizes[i] = All.BoxSize;

#ifdef LONG_X
  All.BoxSizes[0] = All.BoxSize * LONG_X;
#endif
#ifdef LONG_Y
  All.BoxSizes[1] = All.BoxSize * LONG_Y;
#endif
#ifdef LONG_Z
  All.BoxSizes[2] = All.BoxSize * LONG_Z;
#endif

  EgyInjection = 0;

#ifdef PMGRID
  if(RestartFlag != RESTART_FOF_SUBFIND && RestartFlag != RESTART_SLICE && RestartFlag != RESTART_SNAP_CONVERSION &&
     RestartFlag != RESTART_GAS_VELOCITY_POWER_SPECTRUM && RestartFlag != RESTART_TRACER_POWER_SPECTRA)
    long_range_init();
#endif

  if(RestartFlag <= RESTART_SNAPSHOT)
    open_logfiles();

  All.TimeLastRestartFile = CPUThisRun;

#ifdef REDUCE_FLUSH
  All.FlushLast = CPUThisRun;
#endif

#ifdef EOS_DEGENERATE
#ifndef VARIABLE_GAMMA
#error "EOS_DEGENERATE requires VARIABLE_GAMMA"
#endif /* #ifndef VARIABLE_GAMMA */
  eos_init(All.EosTable, All.EosSpecies);
#endif /* #ifdef EOS_DEGENERATE */

#ifdef NUCLEAR_NETWORK
  network_init(All.EosSpecies, All.NetworkRates, All.NetworkPartFunc, All.NetworkMasses, All.NetworkWeakrates, &All.nd);

  {
    for(int k = 0; k < NUM_THREADS; k++)
      network_workspace_init(&All.nd, &All.nw[k]);
  }

#ifdef NETWORK_NSE
  network_nse_init(All.EosSpecies, All.NetworkRates, All.NetworkPartFunc, All.NetworkMasses, All.NetworkWeakrates, &All.nd_nse,
                   All.nw_nse);
#endif /* #ifdef NETWORK_NSE */

#endif /* #ifdef NUCLEAR_NETWORK */

#ifdef EOS_OPAL
#ifndef VARIABLE_GAMMA
#error "EOS_OPAL requires VARIABLE_GAMMA"
#endif /* #ifndef VARIABLE_GAMMA */
  if(opaleos_init(All.EosOpalTable) < 0)
    terminate("Error in OPAL EOS initialization!");
  /* get rho boundaries */
  opaleos_get_rho_limits(&opal_rhomin, &opal_rhomax);
#endif /* #ifdef EOS_OPAL */

#ifdef GENERAL_RELATIVITY
#if METRIC_TYPE == 3 || METRIC_TYPE == 4
  if((read_fixed_numerical_1d_metric()) != 0)
    {
      printf("problems with reading fixed metric");
      terminate("bled");
    }
#endif /* #if METRIC_TYPE==3 || METRIC_TYPE==4 */
#endif /* #ifdef GENERAL_RELATIVITY */

#ifdef COSMIC_RAYS
  init_cosmic_rays();
#endif

#ifdef DUST_LIVE
  init_dust();
#endif

  init_scalars();

  init_gradients();

#ifdef MRT
  init_gradients_RT();
#endif

#ifdef RT_ADVECT
  rt_init_gradients();
#endif

#ifdef SECOND_DERIVATIVES
  init_hessians();
#endif

#ifdef AMR
  amr_init();
#endif

#if defined(VS_TURB) || defined(AB_TURB)
  init_turb();
#ifdef POWERSPEC_GRID
  powersepc_turb_init();
#endif /* #ifdef POWERSPEC_GRID */
#endif /* #if defined (VS_TURB) || defined (AB_TURB) */

#ifdef SIDM
  sidm_Init_CrossSection();
#endif

#ifdef DG
  dg_initialize();
#endif

#ifdef GRAVITY_TABLE
  grav_table_init();
#endif

#ifdef RELAXOBJECT_COOLING2
  load_temperature_profil();
#endif

#ifdef HCOUTPUT
  hcoutput_init();
#endif

#ifdef AURIGA_MOVIE
  auriga_movie_init();
#endif

#ifdef SFR_MCS_LOG
  setup_sf_log();
#endif

#ifdef SGS_TURBULENCE
  init_sgs_turbulence();
#endif

#ifdef SN_MCS_LOG
  setup_sn_log();
#endif

#if((defined(SN_MCS) && !defined(SN_MCS_SINGLE_INJECTION)) || (defined(HII_MCS) && !defined(HII_MCS_TEST)) || defined(PE_MCS)) && \
    !(defined(SN_MCS_INITIAL_DRIVING) || defined(IMF_SAMPLING_MCS))
  read_sb99_tables();
#endif

#ifdef IMF_SAMPLING_MCS
  init_star_properties_tables();
#endif

#ifdef HIGH_FREQUENCY_OUTPUT_STARS
  init_high_freq_star_outputs();
#endif
}

/*! \brief Late setup, after the IC file has been loaded but before run() is
 *  called.
 *
 *  The output files are opened and various modules are initialized. The next
 *  output time is determined by find_next_outputtime() and various timers are
 *  set.
 *
 *  \return void
 */
void begrun2(void)
{
  char contfname[MAXLEN_PATH];
  file_path_sprintf(contfname, "%s/cont", All.OutputDir);
  unlink(contfname);

  delete_end_file();

  if(RestartFlag > RESTART_SNAPSHOT)
    open_logfiles();

#if defined(DG_SET_IC_FROM_AVERAGES) && defined(DG)
  load_weights_from_averages();
#endif /* #if defined(DG_SET_IC_FROM_AVERAGES) && defined(DG) */

#if defined(USE_SFR) && !defined(LOCAL_FEEDBACK)
  sfr_init();
#endif /* #if defined(USE_SFR) && !defined(LOCAL_FEEDBACK) */

#ifdef GFM_STELLAR_EVOLUTION
  init_SNIa_rates();
  mpi_printf("GFM_STELLAR_EVOLUTION: Type Ia rates initialized.\n");
#endif /* #ifdef GFM_STELLAR_EVOLUTION */

#ifdef PMGRID
  long_range_init_regionsize();
#endif /* #ifdef PMGRID */

#ifdef RT_ADVECT
  if(RestartFlag == RESTART_IC)
    {
      rt_set_simple_inits();
      rt_init_sourceid();
#ifdef RT_STELLAR_SOURCES
      rt_create_source_list();
#endif /* #ifdef RT_STELLAR_SOURCES */
    }
#ifdef RT_HEALPIX_NSIDE
  rt_get_vectors();
#endif /* #ifdef RT_HEALPIX_NSIDE */
#endif /* #ifdef RT_ADVECT */

#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
  special_particle_create_list();
#endif /* #ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE */

#ifdef REFINEMENT_AROUND_DM
  dm_particle_create_list();
#endif /* #ifdef REFINEMENT_AROUND_DM */

#if(defined(CIRCUMSTELLAR_IRRADIATION) || defined(ALPHA_VISCOSITY) || defined(CIRCUMSTELLAR_REFINEMENTS)) && !defined(EXTERNALGRAVITY)
  source_particle_create_list();
#endif /* #if (defined(CIRCUMSTELLAR_IRRADIATION) || defined(ALPHA_VISCOSITY) || defined(CIRCUMSTELLAR_REFINEMENTS)) && !defined \
          (EXTERNALGRAVITY) */

  /* this needs to be (re)done here because All.TimeBegin may have changed its value when the restart file was read in (if earlier a
   * restart from a snaphot with option 2 was done) */
  if(All.ComovingIntegrationOn)
    init_drift_table();

  if(RestartFlag == RESTART_SNAPSHOT)
    All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 100);
  else
    All.Ti_nextoutput = find_next_outputtime(All.Ti_Current);

#ifdef OTVET
  if(RestartFlag == RESTART_IC)
    otvet_set_simple_inits();

  ot_get_sigma();

#ifdef OTVET_MULTI_FREQUENCY
#if defined(EDDINGTON_TENSOR_STARS)
  ot_get_lum_stars();
#endif /* #if defined(EDDINGTON_TENSOR_STARS) */
#endif /* #ifdef OTVET_MULTI_FREQUENCY */
#endif /* #ifdef OTVET */

#ifdef TRACER_TRAJECTORY
  tracer_init_output_configuration();
#endif /* #ifdef TRACER_TRAJECTORY */

  All.TimeLastRestartFile = CPUThisRun;

#ifdef REDUCE_FLUSH
  All.FlushLast = CPUThisRun;
#endif /* #ifdef REDUCE_FLUSH */

#if defined(FORCETEST) && defined(FORCETEST_TESTFORCELAW)
  gravity_forcetest_testforcelaw();
#endif /* #if defined(FORCETEST) && defined(FORCETEST_TESTFORCELAW) */

#if defined(SHOCK_FINDER_POST_PROCESSING) || defined(SHOCK_FINDER_BEFORE_OUTPUT) || defined(SHOCK_FINDER_ON_THE_FLY)
  SfVars.limit_gradients = 1;
#endif

#ifdef SHOCK_FINDER_ON_THE_FLY /* create full mesh and run shock finder */
  if(RestartFlag == RESTART_RESTART || RestartFlag == RESTART_SNAPSHOT)
    {
      int k;
      short int *buTimeBin = (short int *)mymalloc_movable(&buTimeBin, "buTimeBin", NumGas * sizeof(short int));
      static int buTimeBinActive[TIMEBINS];

      for(k = 0; k < NumGas; k++)
        {
          buTimeBin[k]      = P[k].TimeBinHydro;
          P[k].TimeBinHydro = 0;
        }

      for(k = 0; k < TIMEBINS; k++)
        {
          buTimeBinActive[k] = TimeBinSynchronized[k];

          TimeBinSynchronized[k] = 1;
        }

      reconstruct_timebins();

      create_mesh();
      mesh_setup_exchange();
      shock_finder_on_the_fly();
      free_mesh();

      for(k = 0; k < TIMEBINS; k++)
        TimeBinSynchronized[k] = buTimeBinActive[k];

      for(k = 0; k < NumGas; k++)
        P[k].TimeBinHydro = buTimeBin[k];

      reconstruct_timebins();

      myfree_movable(buTimeBin);
    }
#endif /* #ifdef SHOCK_FINDER_ON_THE_FLY */

#ifdef SIMPLEX
  if(RestartFlag == RESTART_IC || RestartFlag == RESTART_SNAPSHOT)
    sx_initialize();
#endif /* #ifdef SIMPLEX */

#ifdef RELAXOBJECT_COOLING2
  load_temperature_profil();
#endif /* #ifdef RELAXOBJECT_COOLING2 */

#ifdef BH_BASED_CGM_ZOOM
  // set initial state, which will be kept if no BH has yet been seeded
  All.BlackHolePosition[0] = -1.0;
  All.BlackHolePosition[1] = -1.0;
  All.BlackHolePosition[2] = -1.0;
  All.BlackHoleTask        = -1;
  bh_based_cgm_zoom_update();
#endif /* #ifdef BH_BASED_CGM_ZOOM */
}

/*! \brief Computes conversion factors between internal code units and the
 *  cgs-system.
 *
 *  In addition constants like the gravitation constant are set.
 *
 *  \return void
 */
void set_units(void)
{
  double meanweight;

  All.UnitTime_in_s         = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
  All.UnitTime_in_Megayears = All.UnitTime_in_s / SEC_PER_MEGAYEAR;

  if(All.GravityConstantInternal == 0)
    All.G = GRAVITY / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g * pow(All.UnitTime_in_s, 2);
  else
    All.G = All.GravityConstantInternal;

#ifdef BECDM
  All.hbar   = HBAR / pow(All.UnitLength_in_cm, 2) / All.UnitMass_in_g * All.UnitTime_in_s;
  All.mAxion = All.AxionMassEv * ELECTRONVOLT_IN_ERGS / pow(CLIGHT_REAL, 2) / All.UnitMass_in_g;
  if(All.ComovingIntegrationOn)
    {
      /* need to convert to internal code units (with “little h”) */
      All.hbar *= pow(All.HubbleParam, 2);
      All.mAxion *= All.HubbleParam;
    }
#endif

  All.UnitDensity_in_cgs     = All.UnitMass_in_g / pow(All.UnitLength_in_cm, 3);
  All.UnitPressure_in_cgs    = All.UnitMass_in_g / All.UnitLength_in_cm / pow(All.UnitTime_in_s, 2);
  All.UnitCoolingRate_in_cgs = All.UnitPressure_in_cgs / All.UnitTime_in_s;
  All.UnitEnergy_in_cgs      = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);

  /* convert some physical input parameters to internal units */

  All.Hubble = HUBBLE * All.UnitTime_in_s;

  mpi_printf("BEGRUN: Hubble (internal units)   = %g\n", All.Hubble);
  mpi_printf("BEGRUN: G (internal units)        = %g\n", All.G);
  mpi_printf("BEGRUN: UnitMass_in_g             = %g\n", All.UnitMass_in_g);
  mpi_printf("BEGRUN: UnitLength_in_cm          = %g\n", All.UnitLength_in_cm);
  mpi_printf("BEGRUN: UnitTime_in_s             = %g\n", All.UnitTime_in_s);
  mpi_printf("BEGRUN: UnitVelocity_in_cm_per_s  = %g\n", All.UnitVelocity_in_cm_per_s);
  mpi_printf("BEGRUN: UnitDensity_in_cgs        = %g\n", All.UnitDensity_in_cgs);
  mpi_printf("BEGRUN: UnitEnergy_in_cgs         = %g\n", All.UnitEnergy_in_cgs);
  mpi_printf("\n");

  /* note: assuming NEUTRAL GAS */
  meanweight = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC);

  if(All.MinEgySpec == 0)
    {
      All.MinEgySpec = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.MinGasTemp;
      All.MinEgySpec *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

      mpi_printf("BEGRUN: MinEgySpec set to %g based on MinGasTemp=%g\n", All.MinEgySpec, All.MinGasTemp);
    }

#ifdef SMUGGLE_RADIATION_FEEDBACK
  All.PhotoionizationEgySpec = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.PhotoionizationGasTemp;
  All.PhotoionizationEgySpec *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
#endif

#if defined(USE_SFR) && !defined(SMUGGLE_SFR) && !defined(LOCAL_FEEDBACK) && !defined(SFR_MCS)
  set_units_sfr();
#endif /* #if defined(USE_SFR) && !defined(SMUGGLE_SFR) && !defined(LOCAL_FEEDBACK) */

#ifdef CONDUCTION  // VITALI
  init_spitzer_conductivity();
#endif

#ifdef MONOTONE_CONDUCTION
  init_conductivity();
#endif /* #ifdef MONOTONE_CONDUCTION */

#ifdef IMPLICIT_OHMIC_DIFFUSION
  init_ohm_conductivity();
#endif /* #ifdef IMPLICIT_OHMIC_DIFFUSION */

#ifdef BRAGINSKII_VISCOSITY
  init_braginskii_viscosity();
#endif /* #ifdef BRAGINSKII_VISCOSITY */

#ifdef STATICNFW
#ifdef NFW_h
  All.Hubble *= NFW_h;
#endif /* #ifdef NFW_h */
  R200    = pow(NFW_M200 * All.G / (100 * All.Hubble * All.Hubble), 1.0 / 3);
  Rs      = R200 / NFW_C;
  Dc      = 200.0 / 3 * NFW_C * NFW_C * NFW_C / (log(1 + NFW_C) - NFW_C / (1 + NFW_C));
  RhoCrit = 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
  V200    = 10 * All.Hubble * R200;
  mpi_printf("V200= %g\n", V200);

  fac         = 1.0;
  double Mtot = enclosed_mass(R200);
  mpi_printf("M200= %g\n", Mtot);
  fac  = V200 * V200 * V200 / (10 * All.G * All.Hubble) / Mtot;
  Mtot = enclosed_mass(R200);
  mpi_printf("M200= %g\n", Mtot);
#endif /* #ifdef STATICNFW */
}

/*! \brief deletes the end file if it exists.
 *
 *  This is needed in case a already completed simulation is extended or
 *  overwritten. Note that the end file is completely passive.
 *
 *  \return void
 */
static void delete_end_file(void)
{
  if(RestartFlag > RESTART_SNAPSHOT) /* no simulation happening */
    return;

  char endfname[MAXLEN_PATH];
  file_path_sprintf(endfname, "%s/end", All.OutputDir);
  unlink(endfname);
}

#ifdef HIGH_FREQUENCY_OUTPUT_STARS
void init_high_freq_star_outputs(void)
{
  FILE *fd;
  char buf[512];

  if(!(fd = fopen(All.HighFreqStarsPath, "r")))
    {
      terminate("can't read output list in file '%s'", All.HighFreqStarsPath);
    }

  All.HighFreqStarsSnapshotCount = 0;
  All.HighFreqStarsSnapshotNum   = 0;

  while(1)
    {
      if(fgets(buf, 500, fd) != buf)
        break;

      if(All.HighFreqStarsSnapshotNum >= 10000)
        terminate("too many entries in %s", All.HighFreqStarsPath);

      sscanf(buf, " %lg", &All.HighFreqStarsOutputTimes[All.HighFreqStarsSnapshotNum]);
      All.HighFreqStarsSnapshotNum++;
    }

  fclose(fd);

  mpi_printf("\nBEGRUN: found %d times in %s.", All.HighFreqStarsSnapshotNum, All.HighFreqStarsPath);
  mpi_printf("\nBEGRUN: next star output (%d) will be written at t=%g.\n\n", All.HighFreqStarsSnapshotCount,
             All.HighFreqStarsOutputTimes[All.HighFreqStarsSnapshotNum]);
}
#endif
