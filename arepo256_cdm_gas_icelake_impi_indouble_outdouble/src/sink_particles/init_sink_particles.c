#include "../allvars.h"
#include "../proto.h"

/* PCC - 02.01.2014
        Takes care of the initialization of the ITA-style sink
        particles. Called early in AREPO, with different modes
        depending on where it's called (before/after ICs/snapshot)
        is read
*/
void init_sink_particles(int mode)
{
  int LastTaskNewSink;
  char msg[100];

  mpi_printf("SINK_PARTICLES: Initializing the global variables. In mode %d \n", mode);

  if(mode == 0)
    {
      /* This branch is called from begrun.c */

      NSinksAllTasksOld          = 0;
      NSinksAllTasks             = 0;
      NSinkBufferSize            = SINKLISTLENGTH;
      LastTaskNewSink            = 0;
      SinksFormedSinceLastDomain = 0;
#ifdef SINK_PARTICLES_OUTPUT_EVERY_NEW_SINK
      NSinksLastSnapshot = 0;
#endif

      /* Broadcast the variables set in the parameter file to the other CPUs */
      MPI_Bcast(&SinkCreationDensityCodeUnits, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&SinkFormationRadius, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&SinkEvolutionDumpRateYears, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      set_sink_particle_parameters();

      mpi_printf("SINK_PARTICLES: Evolution dump rate in years %g \n", SinkEvolutionDumpRateYears);
      mpi_printf("SINK_PARTICLES: ... variables set!\n");
#ifdef DEBUG_SINK_PARTICLES
      printf("SINK PARTICLES: debug formation radius %g acc rad %g \n", SinkFormationRadius, SinkAccretionRadius);
#endif
    }
  else
    {
      /* This branch is called from init.c, once the particle structures have been read, etc. */
      SinksLastEvolutionDumpTime = All.Time; /* this ensures that we get a dump on startup */
      /* Read the histories of the sink particles... init() only gets called if RestartFlag != 1
       * so it's safe to call it here (it won't overwrite the correct histories in the restart files) */

      MPI_Bcast(&SinkCreationDensityCodeUnits, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&SinkFormationRadius, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&SinkEvolutionDumpRateYears, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      set_sink_particle_parameters();

      mpi_printf("SINK_PARTICLES: Evolution dump rate in years %g \n", SinkEvolutionDumpRateYears);
      mpi_printf("SINK_PARTICLES: ... variables set!\n");

      int success;
      if(RestartFlag == 2)
        success = load_sink_particle_snapshot(All.SnapshotFileCount - 1);
      else if(RestartFlag == 5)
        success = load_sink_particle_snapshot(RestartSnapNum);
      else if(RestartFlag == 0)
        {
          mpi_printf("Starting from IC\n");
          return;
        }

      if(!success)
        {
          warn(
              "SINK_PARTICLES: we could not locate the file holding the sink particles properties,"
              "we will therefore reset all of their properties to some sensible value\n");
          int totalsinks = get_all_sink_particle_info(2);
        }
#ifdef SINK_PARTICLES_OUTPUT_EVERY_NEW_SINK
      NSinksLastSnapshot = get_all_sink_particle_info(0);
#endif
    }
}

/*! \brief Setting sink particle parameters such as creation density in correct units
 *
 *  This is needed only for initialization in case of a fixed unit-system
 *  but every time step for comoving units.
 *
 */
void set_sink_particle_parameters(void)
{
  if(All.ComovingIntegrationOn)
    {
      SinkCreationDensityCurrent = SinkCreationDensityCodeUnits / All.cf_a3inv / (All.HubbleParam * All.HubbleParam);
      SinkFormationRadiusCurrent = SinkFormationRadius / All.Time * All.HubbleParam;
    }
  else
    {
      SinkCreationDensityCurrent = SinkCreationDensityCodeUnits / (All.HubbleParam * All.HubbleParam);
      SinkFormationRadiusCurrent = SinkFormationRadius * All.HubbleParam;
    }

  SinkFormationRadiusSquared    = SinkFormationRadiusCurrent * SinkFormationRadiusCurrent;
  SinkAccretionRadius           = SinkFormationRadiusCurrent;
  SinkSofteningLength           = SinkAccretionRadius * 0.2;
  SinkTestsGiveupDensityCurrent = SinkCreationDensityCurrent * 100;

  mpi_printf("SINK_PARTICLES: Creation density %g  and radius %g \n", SinkCreationDensityCurrent, SinkFormationRadiusCurrent);
}
