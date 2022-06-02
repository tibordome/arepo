#include "../allvars.h"
#include "../proto.h"

/* PCC - 19.12.2013
        This function loads the data for all sinks in the simulation into
        external structure SinkP. Returns the number of sinks.
        The current version tries out a different model of the MPI communication
        than is common for our previous sink particle algorithms. Instead of
        having the usual send/rec setup, where all tasks communicate the
        information in their buffers, we've adopted a Bcast setup, where only those
        tasks that have sinks, Bcast their information to the others. The rational
        is that for many applications where these sinks will be used, only a handful
        of tasks will hold active sinks. In the future, we should allow this function
        to implement the old communication model if it's more appropriate.
*/
int get_all_sink_particle_info(int mode)
{
  int *NSinksEachTask;
  int NSinksThisTask;
  int SinkList[SINKLISTLENGTH];
#ifdef SINK_PARTICLES_VARIABLE_ACC_RADIUS
  double radius;
#endif

  struct unique_to_sinks
  {
    double FormationMass;
    double FormationTime;
    int FormationOrder;
    long long ID;
#ifdef SGCHEM_ACCRETION_LUMINOSITY
    double AccretionRate;
    double TimeOld;
#endif
#ifdef SINK_PARTICLES_FEEDBACK
    int N_sne;
    double StellarMass;
    double explosion_time[MAXSNE];
    double MassStillToConvert[MAXACCRETIONEVENTS];
    double AccretionTime[MAXACCRETIONEVENTS];
#ifdef POPIII_SNE
    double stellar_mass[MAXSNE];
#endif
#endif
#ifdef SINK_PARTICLES_VARIABLE_ACC_RADIUS
    double AccretionRadius;
#endif
#ifdef SINK_SIMPLEX
    double MassOld;
#endif
#ifdef SINK_PARTICLES_VARIABLE_CREATION
    double CreationDensity;
#ifdef SINK_PARTICLES_FEEDBACK
    double SFeff;
#endif
#endif
#ifdef STORE_SINK_PARTICLE_SPIN
    double AngularMomentum[3];
#endif
  } * store_SinkP;

  /* Find the number of sinks on this Task and their location in
     the particle list
  */
  for(int i = NSinksThisTask = 0; i < NumPart; i++)
    {
      if(P[i].Type == PTYPE_BNDRY)
        {
          if(P[i].ID == 0 && P[i].Mass == 0)
            continue;

          SinkList[NSinksThisTask] = i;
          NSinksThisTask++;
          /*          printf("SINK_PARTICLES: in get_sinks, found sink particle with ID %d \n", P[i].ID);*/
          /* XXX-pcc: need to be careful here that we don't run out of space in the
             the SinkList array! Fix this later.
          */
        }
    }

  int eqSinkList[SINKLISTLENGTH];
  for(int i = 0; i < NSinksThisTask; i++)
    {
      int ind = SinkList[i];

      for(int j = 0; j < NSinksAllTasks; j++)
        if(SinkP[j].ID == P[ind].ID)
          {
            eqSinkList[i] = j;
            break;
          }
    }

  /* Now create an array that holds the number of sinks on all Tasks */
  NSinksEachTask = (int *)mymalloc("NSinksEachTask", NTask * sizeof(int));
  MPI_Allgather(&NSinksThisTask, 1, MPI_INT, NSinksEachTask, 1, MPI_INT, MPI_COMM_WORLD);

  /* Get the total number of sinks in simulation, and work out the offsets
     for the communication of the sink data.
  */
  int nSinksAllTasks_old = NSinksAllTasks;

  for(int i = NSinksAllTasks = 0; i < NTask; i++)
    NSinksAllTasks += NSinksEachTask[i];

  if(mode == 0)
    {
      myfree(NSinksEachTask);
      return (NSinksAllTasks);
    }

  if(NSinksAllTasks == 0)
    {
      mpi_printf("SINK_PARTICLES: No sinks in domain yet!\n");
      myfree(NSinksEachTask);
      return NSinksAllTasks;
    }
  mpi_printf("SINK_PARTICLES: %d sinks present in domain \n", NSinksAllTasks);

  /* Check that we still have space in the SinkP array. In the future we could be smart
   * and allow the code to dynamically reallocate. Would require the function to be called again. */
  if(NSinksAllTasks > NSinkBufferSize)
    terminate("Ran out of sink buffer, currently set at %d", NSinkBufferSize);

  /* Now make a copy of the SinkP array since we're going to overwriting it in the loop below */
  store_SinkP = (struct unique_to_sinks *)mymalloc("store_SinkP", nSinksAllTasks_old * sizeof(struct unique_to_sinks));
  for(int i = 0; i < nSinksAllTasks_old; i++)
    {
      store_SinkP[i].FormationMass  = SinkP[i].FormationMass;
      store_SinkP[i].FormationTime  = SinkP[i].FormationTime;
      store_SinkP[i].FormationOrder = SinkP[i].FormationOrder;
      store_SinkP[i].ID             = SinkP[i].ID;

#ifdef SGCHEM_ACCRETION_LUMINOSITY
      store_SinkP[i].AccretionRate = SinkP[i].AccretionRate;
      store_SinkP[i].TimeOld       = SinkP[i].TimeOld;
#endif
#ifdef SINK_PARTICLES_FEEDBACK
      if(mode == 2)
        {
          store_SinkP[i].N_sne = 0;
          memset(SinkP[i].MassStillToConvert, 0, MAXACCRETIONEVENTS * sizeof(double));
          store_SinkP[i].StellarMass = 0.;
        }
      else
        {
          store_SinkP[i].N_sne = SinkP[i].N_sne;
#ifdef POPIII_SNE
          for(int k = 0; k < SinkP[i].N_sne; k++)
            store_SinkP[i].stellar_mass[k] = SinkP[i].stellar_mass[k];
#endif
          for(int k = 0; k < SinkP[i].N_sne; k++)
            store_SinkP[i].explosion_time[k] = SinkP[i].explosion_time[k];
          for(int k = 0; k < MAXACCRETIONEVENTS; k++)
            {
              store_SinkP[i].MassStillToConvert[k] = SinkP[i].MassStillToConvert[k];
              store_SinkP[i].AccretionTime[k]      = SinkP[i].AccretionTime[k];
            }
          store_SinkP[i].StellarMass = SinkP[i].StellarMass;
        }
#endif
#ifdef SINK_PARTICLES_VARIABLE_ACC_RADIUS
      store_SinkP[i].AccretionRadius = SinkP[i].AccretionRadius;
#endif
#ifdef SINK_SIMPLEX
      store_SinkP[i].MassOld = SinkP[i].MassOld;
#endif
#ifdef SINK_PARTICLES_VARIABLE_CREATION
      store_SinkP[i].CreationDensity = SinkP[i].CreationDensity;
#ifdef SINK_PARTICLES_FEEDBACK
      store_SinkP[i].SFeff = SinkP[i].SFeff;
#endif
#endif
#ifdef STORE_SINK_PARTICLE_SPIN
      store_SinkP[i].AngularMomentum[0] = SinkP[i].AngularMomentum[0];
      store_SinkP[i].AngularMomentum[1] = SinkP[i].AngularMomentum[1];
      store_SinkP[i].AngularMomentum[2] = SinkP[i].AngularMomentum[2];
#endif
    }
#ifdef DEBUG_SINK_PARTICLES
  printf("SINK_PARTICLES: GET -- just set the store_Sink array. TASK %d \n", ThisTask);
  mpi_printf("SINK_PARTICLES: GET -- Sink formass %g formtime %g \n", store_SinkP[0].FormationMass, store_SinkP[0].FormationTime);
#endif

  /* In the loop below, we do the sink communication between the Tasks.
     Looping around each Task, we allow those containing sinks to Bcast
     their data to the other tasks. */

  export_SinkP =
      (struct global_sink_particle_data *)mymalloc("export_SinkP", NSinksThisTask * sizeof(struct global_sink_particle_data));

  /* Ready the communication buffer for this set of sinks */
  for(int i = 0; i < NSinksThisTask; i++)
    {
      int j  = SinkList[i];
      int ii = eqSinkList[i];
      if(mode == 2)
        {
          export_SinkP[i].Pos[0]         = P[j].Pos[0];
          export_SinkP[i].Pos[1]         = P[j].Pos[1];
          export_SinkP[i].Pos[2]         = P[j].Pos[2];
          export_SinkP[i].Vel[0]         = P[j].Vel[0];
          export_SinkP[i].Vel[1]         = P[j].Vel[1];
          export_SinkP[i].Vel[2]         = P[j].Vel[2];
          export_SinkP[i].Accel[0]       = 0.;
          export_SinkP[i].Accel[1]       = 0.;
          export_SinkP[i].Accel[2]       = 0.;
          export_SinkP[i].Mass           = P[j].Mass;
          export_SinkP[i].ID             = (long long)P[j].ID;
          export_SinkP[i].HomeTask       = ThisTask;
          export_SinkP[i].Index          = j;
          export_SinkP[i].FormationMass  = P[j].Mass;
          export_SinkP[i].FormationTime  = All.Time;
          export_SinkP[i].FormationOrder = j;
#ifdef SGCHEM_ACCRETION_LUMINOSITY
          export_SinkP[i].AccretionRate = 0.0;
          export_SinkP[i].TimeOld       = All.Time;
#endif
#ifdef SINK_PARTICLES_FEEDBACK
          export_SinkP[i].N_sne = 0;
          memset(export_SinkP[i].MassStillToConvert, 0, MAXACCRETIONEVENTS * sizeof(double));
          export_SinkP[i].StellarMass = 0.;
#endif
#ifdef SINK_PARTICLES_VARIABLE_ACC_RADIUS
          export_SinkP[i].AccretionRadius = SinkAccretionRadius;
#endif
#ifdef SINK_SIMPLEX
          export_SinkP[i].MassOld = 0;
#endif
#ifdef SINK_PARTICLES_VARIABLE_CREATION
          export_SinkP[i].CreationDensity = SinkCreationDensityCodeUnits;
#ifdef SINK_PARTICLES_FEEDBACK
          export_SinkP[i].SFeff = All.SINKStarFormationEfficiency;
#endif
#endif
        }
      else
        {
          export_SinkP[i].Pos[0]   = P[j].Pos[0];
          export_SinkP[i].Pos[1]   = P[j].Pos[1];
          export_SinkP[i].Pos[2]   = P[j].Pos[2];
          export_SinkP[i].Vel[0]   = P[j].Vel[0];
          export_SinkP[i].Vel[1]   = P[j].Vel[1];
          export_SinkP[i].Vel[2]   = P[j].Vel[2];
          export_SinkP[i].Accel[0] = P[j].GravAccel[0];
          export_SinkP[i].Accel[1] = P[j].GravAccel[1];
          export_SinkP[i].Accel[2] = P[j].GravAccel[2];
          export_SinkP[i].Mass     = P[j].Mass;
          export_SinkP[i].ID       = (long long)P[j].ID;
          export_SinkP[i].HomeTask = ThisTask;
          export_SinkP[i].Index    = j;
#ifdef SINK_PARTICLES_VARIABLE_ACC_RADIUS
          radius = pow(3.0 * P[j].Mass / 4. / M_PI / SinkCreationDensityCurrent, 0.3333333);
          if(SinkAccretionRadius > radius)
            radius = SinkAccretionRadius;
          export_SinkP[i].AccretionRadius = radius;
#endif
          /* Keep track of particle properties that are ONLY recorded on the
             SinkP array, such as formation mass and the order in which the
             sinks formed. Could also keep age and stellar properties here too.
           */
          if(P[j].ID != store_SinkP[ii].ID)
            terminate("SINK ID mismatch %" MYIDTYPE_PRI " %lld\n", P[j].ID, store_SinkP[ii].ID);
          export_SinkP[i].FormationMass  = store_SinkP[ii].FormationMass;
          export_SinkP[i].FormationTime  = store_SinkP[ii].FormationTime;
          export_SinkP[i].FormationOrder = store_SinkP[ii].FormationOrder;
#ifdef STORE_SINK_PARTICLE_SPIN
          export_SinkP[i].AngularMomentum[0] = store_SinkP[ii].AngularMomentum[0];
          export_SinkP[i].AngularMomentum[1] = store_SinkP[ii].AngularMomentum[1];
          export_SinkP[i].AngularMomentum[2] = store_SinkP[ii].AngularMomentum[2];
#endif
#ifdef SGCHEM_ACCRETION_LUMINOSITY
          export_SinkP[i].AccretionRate = store_SinkP[ii].AccretionRate;
          export_SinkP[i].TimeOld       = store_SinkP[ii].TimeOld;
#endif
#ifdef SINK_PARTICLES_FEEDBACK

          export_SinkP[i].N_sne = store_SinkP[ii].N_sne;
          for(int k = 0; k < export_SinkP[i].N_sne; k++)
            export_SinkP[i].explosion_time[k] = store_SinkP[ii].explosion_time[k];
#ifdef POPIII_SNE
          for(int k = 0; k < export_SinkP[i].N_sne; k++)
            export_SinkP[i].stellar_mass[k] = store_SinkP[ii].stellar_mass[k];
#endif
          for(int k = 0; k < MAXACCRETIONEVENTS; k++)
            {
              export_SinkP[i].MassStillToConvert[k] = store_SinkP[ii].MassStillToConvert[k];
              export_SinkP[i].AccretionTime[k]      = store_SinkP[ii].AccretionTime[k];
            }
          export_SinkP[i].StellarMass = store_SinkP[ii].StellarMass;
#endif
#ifdef SINK_SIMPLEX
          export_SinkP[i].MassOld = store_SinkP[ii].MassOld;
#endif
#ifdef SINK_PARTICLES_VARIABLE_CREATION
          export_SinkP[i].CreationDensity = store_SinkP[ii].CreationDensity;
#ifdef SINK_PARTICLES_FEEDBACK
          export_SinkP[i].SFeff = store_SinkP[ii].SFeff;
#endif
#ifdef SINK_PARTICLES_VARIABLE_ACC_RADIUS
          radius = pow(3.0 * P[j].Mass / 4. / M_PI / store_SinkP[ii].CreationDensity, 0.3333333);

          if(SinkAccretionRadius > radius)
            radius = SinkAccretionRadius;
          export_SinkP[i].AccretionRadius = radius;
#endif
#endif
        }
    }

  int *bytecounts = (int *)mymalloc("bytecounts", sizeof(int) * NTask);
  int *byteoffset = (int *)mymalloc("byteoffset", sizeof(int) * NTask);

  int itask;
  for(itask = 0; itask < NTask; itask++)
    bytecounts[itask] = NSinksEachTask[itask] * sizeof(struct global_sink_particle_data);

  for(itask = 1, byteoffset[0] = 0; itask < NTask; itask++)
    byteoffset[itask] = byteoffset[itask - 1] + bytecounts[itask - 1];

  MPI_Allgatherv(&export_SinkP[0], NSinksThisTask * sizeof(struct global_sink_particle_data), MPI_BYTE, &SinkP[0], bytecounts,
                 byteoffset, MPI_BYTE, MPI_COMM_WORLD);

  myfree(byteoffset);
  myfree(bytecounts);

  myfree(export_SinkP);

  /* De-allocate anything used above (needs to be done in the reverse order
     in which they were allocated.) */
  myfree(store_SinkP);
  myfree(NSinksEachTask);

  /* Should be finished! Update the 'old' number of sinks in domain, and
     return the number of sinks to the calling function */
  NSinksAllTasksOld = NSinksAllTasks;

  return (NSinksAllTasks);
}
