/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \date        22/10/2018
 * \author      Ondrej Jaura
 * \brief       Main functionality of the SPRAI (SimpleX Photon Radiation Arepo Implementation) 
 * \details     
 * 
 * Workflow of the code:
 * -> Initialization - at the beginning of the simulation
 * -> Evolve HD/PP - main wrappers for the hydro-dynamical and post-processing runs, at full hydro-steps
 *    Evolve - general RT routine that depends on the parameters "a", "t", and "dt" 
 *    -> Run - one "dt" long radiational transfer starting at "a"/"t"
 *       -> Rotation - one realization of the radiational field
 *                     different realization are integrated at the end to form a single field 
 *          -> Step - one movement of the photons between neighbouring cells
 * -> SGChem - chemistry module that reads the radiational field and calculates abundances
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

//======================================================================================================//
//                                           Initialization                                             //
//======================================================================================================//

/*!
 * \brief Initialization of SimpleX after restart
 */
void sx_restart(void)
{
  mpi_printf("SX: Restarting Simplex \n");

  // we need to initialize random number generator, otherwise we will get a segmentation fault
  sxRand = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(sxRand, All.sxSeed * (ThisTask + 1));
}

/*!
 * \brief Initialization of SimpleX computation
 */
void sx_initialize(void)
{
  int f;

  All.sxPrevTime = 0.0;                  // physical time of the last SimpleX RT [s]
  All.sxCurrentRun = 0;                  // initialize current run variable  
  All.sxSeed = 3;                        // random seed for the random number generator

  // check SPRAI requirements
#ifndef DO_NOT_RANDOMIZE_DOMAINCENTER
  // DEBUG: 
  // For now we need to ensure, that the domain center is not randomized!!!
  // Some preliminary fixes are made in sx_photon_sources.c, however we also need 
  // to properly implement transfer of photons at the edges of periodical box.
  // Without this, the photons are dumped at random positions within the box.
#error "SX-ERROR: SPRAI will not work without a compiler flag DO_NOT_RANDOMIZE_DOMAINCENTER."
#endif
#if SX_RUNNING_MODE != 2
#error "SX-ERROR: SX_RUNNING_MODE != 2 is currently not available"
#endif

  // initialize random number generator
  sxRand = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(sxRand, All.sxSeed * (ThisTask + 1));

  // DEBUG: this should be somehow generalized or set by a particular star
  All.sxTeff = pow(10., 5.);           // we set the effective temperature of the source to 10^5 K

  All.sxMaxStraightAngleCos = cos(SX_MAX_STRAIGHT_ANGLE * M_PI / 180.0);  // cosine of the Maximal Straight Angle

  All.sxNdirInv = 1. / (double)SX_NDIR;  // inverse of the number of directional bins

  mpi_printf("SX: initializing SPRAI radiative transfer\n");

  mpi_fprintf(FdSimplex,"SPRAI - SimpleX Photon Radiation in the Arepo Implementation \n");

  mpi_fprintf(FdSimplex,"\n>--- Initialization\n");
  mpi_fprintf(FdSimplex,"T_eff %.3e \nMinPhot %.3e \nMaxAngle %.1f \nFactorSA %.1f \n", 
	      All.sxTeff, All.sxMinNumPhotons, SX_MAX_STRAIGHT_ANGLE, SX_FACTOR_SA);
  // calculate the memory consumption of the SphP 
  int mSphP = NumGas * sizeof(struct sph_particle_data), mSphPall;
  MPI_Reduce(&mSphP, &mSphPall, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  mpi_fprintf(FdSimplex,"SphP structure size:     particle=%.1f KB, task=%.1f MB, total=%.1f MB \n",
	      (float)sizeof(struct sph_particle_data) * 9.765625e-4, (float)mSphPall * 9.536743164e-7 / (float)NTask,
	      (float)mSphPall * 9.536743164e-7 ); 
  // calculate an approximate number of protons in each cell
  double totalMass = 0, localMass = 0;
  int totalNumGas = 0;
  for(int i = 0; i < NumGas; i++)
      localMass += P[i].Mass;
  if (All.ComovingIntegrationOn)
    localMass *= All.UnitMass_in_g / (PROTONMASS * All.HubbleParam);
  else
    localMass *= All.UnitMass_in_g / (PROTONMASS);
  MPI_Reduce(&localMass, &totalMass, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&NumGas, &totalNumGas, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  mpi_fprintf(FdSimplex,"Mean num. of protons per gas cell: %.04e \n",totalMass/(double)totalNumGas);
  mpi_fprintf(FdSimplex,"Memory allocation ABPP [ %.2f %.2f ] QPP [ %.2f %.2f ] \n",
	      SX_ABPP_INIT,SX_ABPP_REALLOC,SX_QPP_INIT,SX_QPP_REALLOC);

  // initialize ionisation physics
  sx_chem_initialize();

#if SX_CHEMISTRY == 3
  mpi_fprintf(FdSimplex,"Frequenci bins\n Hz    [ ");
  for( f=0; f<SX_NFREQ+1; f++)
      mpi_fprintf(FdSimplex,"%.2e ", All.sxFreqBins[f]);
  mpi_fprintf(FdSimplex,"] \n eV    [ ");
  for( f=0; f<SX_NFREQ+1; f++)
      mpi_fprintf(FdSimplex,"%.1f ", All.sxFreqBins[f]/ELECTRONVOLT_IN_HZ);
  mpi_fprintf(FdSimplex,"] \n");
  mpi_fprintf(FdSimplex,"Cross-sections\n cm^2  [ ");
  for( f=0; f<SX_NSIGMA; f++)
      mpi_fprintf(FdSimplex,"%.2e ", All.sxSigma[f]);
  mpi_fprintf(FdSimplex,"] \n");
  mpi_fprintf(FdSimplex,"Photon energies\n erg   [ ");
  for( f=0; f<SX_NENERGY; f++)
      mpi_fprintf(FdSimplex,"%.2e ", All.sxEnergy[f]);
  mpi_fprintf(FdSimplex,"]\n eV    [ ");
  for( f=0; f<SX_NENERGY; f++)
      mpi_fprintf(FdSimplex,"%.2e ", All.sxEnergy[f]/ELECTRONVOLT_IN_ERGS);
  mpi_fprintf(FdSimplex,"]\n");
#endif

  // initialize photon trasnfer
  sx_phot_initialize();

  // DEBUG: Abort the simulation
  // endrun();
}

//======================================================================================================//
//                                         Run flow functions                                           //
//======================================================================================================//

/*!
 * \brief Do the initialization before the run
 * \params a - expansion factor at the beginning of the RT run
 *         t - physical time at the beginning of the RT run [s]
 *         dt - time duration of the RT run [s]
 */
void sx_start_run(double a, double t, double dt)
{
  int i, idx;

  // increase current run number
  All.sxCurrentRun++;

  // allocate run data structure and initialize variables
  sxRunData = (struct sxRunData_struct *)mymalloc("sxRunData", sizeof(struct sxRunData_struct));

  // initialize debug and start the main timer
  sx_debug_start_run();
  sx_debug_timer_start(TRUN);
#ifdef SX_DISPLAY_TIMERS
  sx_debug_timer_start(TRUNs);
#endif
  
  // initialize current step and rotation counter
  sxRunData->currentRot  = 0;
  sxRunData->currentStep = 0;
  sxRunData->nSteps = 0;

  // calculate unit conversion factors and cosmological parameters
  if(All.ComovingIntegrationOn)
    {
      sxRunData->mass2cgs = All.UnitMass_in_g  / All.HubbleParam;                                // [g]
      sxRunData->density2cgs = All.UnitDensity_in_cgs * pow( All.HubbleParam, 2 ) / pow( a, 3 ); // [g/cm^3]
      sxRunData->length2cgs = All.UnitLength_in_cm * a / All.HubbleParam;                        // [cm] 
      sxRunData->velocity2cgs = All.UnitVelocity_in_cm_per_s / a ;                               // [cm/s]
    }
  else
    {
      sxRunData->mass2cgs = All.UnitMass_in_g;                        // [g]  
      sxRunData->density2cgs = All.UnitDensity_in_cgs;                // [g/cm^3]
      sxRunData->length2cgs = All.UnitLength_in_cm;                   // [cm] 
      sxRunData->velocity2cgs = All.UnitVelocity_in_cm_per_s;         // [cm/s]
    }
  sxRunData->volume2cgs = pow(sxRunData->length2cgs, 3);  // [cm^3]

  sxRunData->time            = t;            // current time [s]
  sxRunData->expansionFactor = a;            // current expansion factor
  sxRunData->redshift        = 1.0 / a - 1;  // current redshift from expansion factor
#ifdef SGCHEM
  if(!All.ComovingIntegrationOn)
    sxRunData->redshift = All.InitRedshift;  // we use initial value from SGChem module
#endif
#ifdef SX_CMB_TEMP_LIMIT
  sxRunData->tempCMB = 2.725 * (1 + sxRunData->redshift);   // temperature of the CMB [K]
#endif
  sxRunData->runTime = dt;                                  // duration time of the RT run [s]
  sxRunData->runTimeInv = 1.0/sxRunData->runTime;           // inverse of the run time
#if (SX_NUM_STEPS>0)
  sxRunData->timeStep = sxRunData->runTime / SX_NUM_STEPS;  // time of a RT step [s]
#else
  sxRunData->timeStep = sxRunData->runTime;
#endif

#ifdef SX_RADIATION_PRESSURE
  sxRunData->cInv = 1.0/CLIGHT;
#endif
  
  mpi_printf("SX: RUN %d starting SimpleX radiation transfer\n", All.sxCurrentRun);

  mpi_fprintf(FdSimplex,"\n>--- RUN %d --- a=%f z=%f t=%f --- dt=%.03e s (%.03e yr)\n",
	     All.sxCurrentRun, sxRunData->expansionFactor, sxRunData->redshift,
	     sxRunData->time/All.UnitTime_in_s,
	     sxRunData->runTime, sxRunData->runTime/SEC_PER_YEAR );

  // allocate array for gas properties used only during this run
  // DEBUG: This should be in the future somehow allocated only for the active particles
  sxCell = (struct sxCell_struct *)mymalloc("sxCell", NumGas * sizeof(struct sxCell_struct) );

  // initialize chemistry
  sx_chem_start_run();

  // initialize radiative transfer
  sx_phot_start_run();

#ifdef SX_DISPLAY_TIMERS
  sx_debug_timer_switch(TRUNs,TRUNi);
#endif
}

/*!
 * \brief Initialize the beginning of the RT step
 */
void sx_start_step(void)
{
#ifdef SX_DISPLAY_TIMERS
  sx_debug_timer_start(TSTEP);
  sx_debug_timer_start(TSTEPs);
#endif
  
  // register a new step
  sxRunData->nSteps++;
  
#ifdef SX_DISPLAY_STATS
  mpi_fprintf(FdSimplex,"\n>--- RUN %d",All.sxCurrentRun);
#if (SX_NUM_ROT>0)
  mpi_fprintf(FdSimplex," ROT %d/%d", sxRunData->currentRot+1, SX_NUM_ROT);
#endif
#if (SX_NUM_STEPS>0)
  mpi_fprintf(FdSimplex," STEP %d/%d",sxRunData->currentStep+1, SX_NUM_STEPS);
#else
  mpi_fprintf(FdSimplex," STEP %d",sxRunData->currentStep+1);
#endif
  mpi_fprintf(FdSimplex," ---");
#if (SX_RUNNING_MODE==2) || (SX_RUNNING_MODE==1)
#if (SX_NUM_STEPS>0)
  mpi_fprintf(FdSimplex," runTime=%.2e s (%.02e yr) timeStep=%.2e s (%.02e yr)", 
	     sxRunData->runTime, sxRunData->runTime/SEC_PER_YEAR,
	     sxRunData->timeStep, sxRunData->timeStep/SEC_PER_YEAR);
#else
  mpi_fprintf(FdSimplex," runTime=%.2e s (%.02e yr)", 
	     sxRunData->runTime, sxRunData->runTime/SEC_PER_YEAR);
#endif
#endif
  mpi_fprintf(FdSimplex,"\n");
#endif

  // initialize radiative transport before each step
  sx_phot_start_step();

#ifdef SX_DISPLAY_TIMERS
  sx_debug_timer_switch(TSTEPs,TSTEPi);
#endif

}

/*!
 * \brief Solve the site for incoming photons and photons from the source
 */
void sx_evolve_sites(void)
{
#ifdef SX_DISPLAY_TIMERS
  sx_debug_timer_switch(TSTEPi,TEVOL);
#endif

#ifdef SX_DISPLAY_STATS
  sxRunData->nDifTransp = 0;
  sxRunData->nBalTransp = 0;
#endif

  // Loop over all active gas particles
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      struct sxSite site;  
#ifdef SX_DISPLAY_TIMERS_SITE
      sx_debug_timer_start(TEVOLsi);
      _Bool doTransport = sx_phot_initialize_site( &site , i );
      sx_debug_timer_end(TEVOLsi);
#else
      _Bool doTransport = sx_phot_initialize_site( &site , i );
#endif 
      if (doTransport)
	{
	  //solve the rate equation to determine ionisations and recombinations
#ifdef SX_DISPLAY_TIMERS_SITE
          sx_debug_timer_start(TEVOLss);
	  sx_chem_solve_site( &site );
          sx_debug_timer_end(TEVOLss);
#else
	  sx_chem_solve_site( &site );
#endif	  
	  //check if there are enough photons for the transport
	  _Bool doTransport = sx_phot_check_n_photons( &site );
	  if (doTransport)
	    {	  
	      if(site.isSource)
		{
	          // redistribute the photons to all neighbours
	          sx_phot_diffuse_transport( &site );
#ifdef SX_DISPLAY_STATS
                  sxRunData->nDifTransp++;
#endif
	        }
	      //redistribute the photons to most straight-forwards neighbours
#ifdef SX_DISPLAY_TIMERS_SITE
	      sx_debug_timer_start(TEVOLtb);
	      sx_phot_ballistic_transport( &site );
	      sx_debug_timer_end(TEVOLtb);
#else
	      sx_phot_ballistic_transport( &site );
#endif
#ifdef SX_DISPLAY_STATS
              sxRunData->nBalTransp++;
#endif
	    }
	  // make final calculations and update values in the cell
	  sx_chem_finalize_site( &site );
	}
    }

#ifdef SX_DISPLAY_TIMERS
  sx_debug_timer_switch(TEVOL,TSTEPi);
#endif

}

/*!
 * \brief Exchange photons after each step and convert them to current
 */
_Bool sx_end_step(void)
{
  _Bool nextStep = 1;
  int f;
  double *nPhot = NULL;

#ifdef SX_DISPLAY_TIMERS
  sx_debug_timer_switch(TSTEPi,TSTEPe);
#endif

  // finish radiative step
  sx_phot_end_step();

#ifdef SX_DISPLAY_STATS
  sx_debug_display_stats();
#endif

#if(SX_NUM_STEPS > 0)
  if(sxRunData->currentStep + 1 >= SX_NUM_STEPS)
    nextStep = 0;
#else

  // count if there are some photon packets, if not, terminate the run
  int numTotAPP = 0;
  MPI_Allreduce(&sxNumAPP, &numTotAPP, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  nextStep = (numTotAPP>0) ? 1 : 0;
#endif

  sxRunData->currentStep++;

#ifdef SX_DISPLAY_TIMERS
  sx_debug_timer_end(TSTEPe);
  sx_debug_timer_end(TSTEP);
#endif

  return nextStep;
}

/*!
 * \brief Finalize after the run
 */
void sx_end_run(void)
{
  int i, idx;

#ifdef SX_DISPLAY_TIMERS  
  sx_debug_timer_switch(TRUNi,TRUNe);
#endif

  // Finalize chemistry
  sx_chem_end_run();

  // Finalize radiation trasfer
  sx_phot_end_run();

  // Free sxCell data
  myfree( sxCell );

#ifdef SX_RADIATION_PRESSURE
  /*
  // DEBUG: I thought update of pressure primitive variables below is necessary, but I have suspicion that 
  //        it it is not, because RT does not change primitive variables of SphP directly
  //        Only radiation pressure changes SphP.Momentum and SphP.Energy
  //        I think that these updates set here and there zero energies to cells and crash the chemistry
  */
  // Loop over active particles
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      // update pressure of the cell
      set_pressure_of_cell(i);
    }

  // Update primitive variables
  update_primitive_variables();
  exchange_primitive_variables();
  calculate_gradients();
  exchange_primitive_variables_and_gradients();
#endif

  // end the main timer and close debugging
#ifdef SX_DISPLAY_TIMERS
  sx_debug_timer_end(TRUNe);
#endif
  double trun = sx_debug_timer_end(TRUN);
  mpi_fprintf(FdSimplex,"\n>--- RUN %d END --- %d steps took %e s\n",
	      All.sxCurrentRun, sxRunData->nSteps, trun);
  mpi_printf("SX: RUN %d took %e s \n",All.sxCurrentRun, trun);
  sx_debug_end_run();     

  // free memory
  myfree(sxRunData);
  
}

//======================================================================================================//
//                                  Main Radiation Transfer Routines                                    //
//======================================================================================================//

/*! 
 * Main loop of the Simplex Radiative Transfer (SRT)
 *
 * \brief Performs the radiative transfer calculation with SimpleX
 *
 * \params a  - expansion factor at the begining of the run (a=1 for non-cosmological simulations)
 *         t  - physical time at the begining of the run [s]
 *         dt - time duration of the RT run [s]
 */
void sx_evolve(double a, double t, double dt)
{
  // initialize the run
  sx_start_run(a, t, dt);

#if(SX_NUM_ROT > 1)
  // we will divide RT into SX_NUM_ROT rotation steps
  while(sxRunData->currentRot < SX_NUM_ROT)
    {
      sxRunData->currentStep = 0;
      sx_phot_start_rot();
#endif

      MPI_Barrier(MPI_COMM_WORLD);

      do
        {
          // initialize step
          sx_start_step();

          // synchronize all processors
          MPI_Barrier(MPI_COMM_WORLD);

          // evolve photons in gas particles
          sx_evolve_sites();

	  MPI_Barrier(MPI_COMM_WORLD);
	  
        } 
      while ( sx_end_step() ); // finish step
	
      MPI_Barrier(MPI_COMM_WORLD);

#if (SX_NUM_ROT > 1)
      sx_phot_end_rot();
      sxRunData->currentRot++;
    }
#endif

  // finalize the run
  sx_end_run();
  
  MPI_Barrier(MPI_COMM_WORLD);

}

/*
 * \brief wrapper function that evolves RT within the hydrodynamical run
 */
void sx_evolve_hd(void)
{
  double a, t, dt;

  // skipping the first synchronization step
  if(All.NumCurrentTiStep == 0)
    return;

  // calculating the expansion factor and the current time
  if(All.ComovingIntegrationOn)
    {
      a  = All.Time;
      t  = 0.0;
      dt = (All.Ti_Current - All.sxPrevTime) * All.Timebase_interval * All.UnitTime_in_s / hubble_function(All.Time) /
           All.HubbleParam;  // time of the last full-hydro-step
    }
  else
    {
      a  = 1.0;
      t  = All.Time * All.UnitTime_in_s;
      dt = (All.Ti_Current - All.sxPrevTime) * All.Timebase_interval * All.UnitTime_in_s;
    }

  // starting RT
  sx_evolve(a, t, dt);

  // save the last run time
  All.sxPrevTime = All.Ti_Current;
}

/*
 * \brief wrapper function that evolves RT as post processing on a snapshot
 */
void sx_evolve_pp(void)
{
  int NumRuns, Run, SnapNum, nSub;
  double TimeBegin, TimeMax, Time, TimeInterval;
  double a, time_s, t, dt;
  double nSubInv;

  // we have to open the log file 'simplex.txt' here because Arepo doesn't at this point
  // NOTE: only one task can open the file at the time
  if (ThisTask==0)
    {
      char buf[MAXLEN_PATH];
      file_path_sprintf(buf, "%s/simplex.txt", All.OutputDir);
      if(!(FdSimplex = fopen(buf, "w")))
	terminate("error in opening file '%s'\n", buf);
    }

  sx_initialize();

  mpi_fprintf(FdSimplex,"\n>--- Start SimpleX post-processing \n");
  
  if (Argv[4]==0||Argv[5]==0)
    mpi_terminate("SX-ERROR: Missing argument <TimeMax> or <numRuns>!");

  Time = All.TimeBegin;
  TimeBegin = All.TimeBegin;
  TimeMax = atof(Argv[4]);
  Run = 0;

  if (RestartSnapNum<1e6)
    {
      NumRuns = atoi(Argv[5]);
      SnapNum = 1000000+RestartSnapNum*1000;
    }
  else
    {
      NumRuns = atoi(Argv[5]) - RestartSnapNum%1000;
      SnapNum = RestartSnapNum;
    }

  if(All.ComovingIntegrationOn)
    {
      a            = Time;
      TimeInterval = (log(TimeMax) - log(TimeBegin)) / (double)NumRuns;
      nSub         = ((double)NumRuns < 1000.0) ? (int)(1000.0 / (double)NumRuns) : (int)NumRuns;
      nSubInv      = 1.0 / nSub;
    }
  else
    {
      a            = 1.0;
      TimeInterval = (TimeMax - TimeBegin) / (double)NumRuns;
    }

  mpi_fprintf(FdSimplex,"Snap %d NumRuns %d TimeBegin %.03e TimeMax %.03e TimeInterval %.03e \n",
	     RestartSnapNum, NumRuns, TimeBegin, TimeMax, TimeInterval); 
  mpi_printf("SX: PP %d/%d | a=%f z=%f dt=%.03e s (%.01e yr) | t=%.03e s (%.01e yr)\n",
	     Run, NumRuns, a, 1.0/a-1.0, 0.0, 0.0, All.sxPrevTime, 0.0);
  
#ifdef SX_PP_ZERO_SNAPSHOT
  mpi_printf("SX: Saving the first snapshot \n");
  savepositions(SnapNum, 0);
#if defined(VORONOI_IMAGES_FOREACHSNAPSHOT)
  make_voronoi_image(SnapNum);
  make_voronoi_image_slice(SnapNum);
#endif
#endif

  // simulate the main Arepo loop
  for(Run = 1, time_s = 0.0; Run <= NumRuns; Run++)
    {
      // calculating the expansion factor and the current time
      if(All.ComovingIntegrationOn)
        {
          // we do substepping to get converging values from the hubble_function
          dt = 0.0;
          for(int i = 1; i <= nSub; i++)
            {
              Time = TimeBegin * exp((Run - 1.0 + i * nSubInv) * TimeInterval);
              dt += TimeInterval * nSubInv * All.UnitTime_in_s / hubble_function(Time) / All.HubbleParam;
            }
          a = Time;
          t = 0.0;
        }
      else
        {
          Time = TimeBegin + Run * TimeInterval;
          a    = 1.0;
          t    = Time * All.UnitTime_in_s;
          dt   = TimeInterval * All.UnitTime_in_s;
        }
      All.Time = Time;
      time_s += dt;

      mpi_printf("\nSX: PP %d/%d | a=%f z=%f dt=%.03e s (%.01e yr) | t=%.03e s (%.01e yr) \n",
                 Run, NumRuns, a, 1.0/a-1.0, dt, dt/SEC_PER_YEAR, time_s, time_s/SEC_PER_YEAR);

      sx_evolve(a, t, dt);

#if(SX_CHEMISTRY == 3)
      // evolve chemistry for each active cell
      mpi_printf("SX: evolving SGChem chemistry \n");
      for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          int i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
            continue;
          evolve_chemistry_for_single_cell(TimeInterval, i);
        }
#endif

      update_primitive_variables();
      exchange_primitive_variables();
      calculate_gradients();
      exchange_primitive_variables_and_gradients();

      savepositions(SnapNum + Run, 0);
#if defined(VORONOI_IMAGES_FOREACHSNAPSHOT)
      make_voronoi_image(SnapNum + Run);
      make_voronoi_image_slice(SnapNum + Run);
#endif

    }

  mpi_fprintf(FdSimplex,"\n>--- End SimpleX post-processing \n");  
  if (ThisTask==0)
    {      
      // close the log file "simplex.txt"
      fclose(FdSimplex);
    }
}
