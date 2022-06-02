
#include "../allvars.h"
#include "../proto.h"

void sink_particles(void)
{
  CPU_Step[CPU_MISC] += measure_time();

  if(All.Time == All.TimeBegin || SinkAccretionRadius == 0)
    return;

  mpi_printf("SINK_PARTICLES: entering sink particle routine \n");
  double t0 = second();

  /* for comoving units e.g. creation density needs to updated every step
   */
  if(All.ComovingIntegrationOn)
    set_sink_particle_parameters();
  /* get the information about all the sinks that are present
     in the simulation
  */
  int num_sinks = get_all_sink_particle_info(1);
  mpi_printf("SINK_PARTICLES: number of sinks on all tasks %d\n", num_sinks);

  double t1 = second();
  mpi_printf("SINK_PARTICLES: collected all sinks, took %g sec.\n", timediff(t0, t1));

#ifdef ALLOW_MULTIPLE_SINK_CREATION_PER_TIMESTEP
  /* Allow new sinks to form
   */
  mpi_printf("SINK_PARTICLES: creating new sink particles \n");
  int success = 1;
  while(success >= 0)
    success = create_sink_particles();

  double t2 = second();
  mpi_printf("SINK_PARTICLES: creation done, took %g sec.\n", timediff(t1, t2));

  /* Allow the current sinks to accrete
   */
  mpi_printf("SINK_PARTICLES: accreting sink particles \n");
  double mass_accreted = accrete_onto_sink_particles();
  mpi_printf("SINK_PARTICLES: mass accreted by all sinks this timestep %g \n", mass_accreted);

  double t3 = second();
  mpi_printf("SINK_PARTICLES: accretion done, took %g sec.\n", timediff(t2, t3));

#else

  mpi_printf("SINK_PARTICLES: Entering accretion function... NSinksAllTasks %d \n", NSinksAllTasks);
  double mass_accreted1 = 0;
  double mass_accreted2 = 0;
  if(NSinksAllTasks > 0)
    mass_accreted1 = accrete_onto_sink_particles();
  mpi_printf("SINK_PARTICLES: mass accreted by all sinks this timestep %g \n", mass_accreted1);

  /* Allow new sinks to form
   */
  int success;
  success = create_sink_particles();
  if(success <= 0 && ThisTask == 0)
    printf("SINK_PARTICLES: success flag %d \n", success);
  if(success > 0)
    {
      /* New sink present, so we need to updated the SinkP lists on all tasks
         before letting the new sink accrete mass
      */
      num_sinks = get_all_sink_particle_info(1);
      mpi_printf("SINK_PARTICLES: number of sinks on all tasks after creation %d \n", num_sinks);
      mass_accreted2 = accrete_onto_sink_particles();
      mpi_printf("SINK_PARTICLES: mass accreted after sink creation %g \n", mass_accreted2);
    }

  if(mass_accreted1 + mass_accreted2 != 0)
    num_sinks = get_all_sink_particle_info(1);

#endif /*ALLOW_MULTIPLE_SINK_CREATION_PER_TIMESTEP*/

#ifdef SINK_MERGERS
  perform_mergers();
  num_sinks = get_all_sink_particle_info(1);
  mpi_printf("SINK MERGERS new Nsink = %d \n", num_sinks);
#endif

#ifdef DUMP_SINK_PARTICLE_INFO
  /* Dump the sink particle data to a file
   */
  if(NSinksAllTasks > 0)
    dump_sink_particle_info(success);
#endif

#ifdef SINK_PARTICLES_FEEDBACK
  eliminate_old_sinks();
#endif

  MPI_Barrier(MPI_COMM_WORLD);

  double t4 = second();

  mpi_printf("SINK_PARTICLES: done, took %g sec.\n", timediff(t0, t4));

  CPU_Step[CPU_SINKS] += measure_time();
}

#ifdef SINK_PARTICLES_FEEDBACK
void eliminate_old_sinks(void)
{
  double t0 = second();

  double toleranceTime = 10. * SEC_PER_MEGAYEAR / All.UnitTime_in_s;

  for(int isink = 0; isink < NSinksAllTasks; isink++)
    {
      if(SinkP[isink].N_sne == 0)
        if((All.Time - SinkP[isink].FormationTime) > toleranceTime)
          {
            reinject_gas_mass(isink);
            remove_sink_particle(isink);
          }
    }

  double t1 = second();

  mpi_printf("SINK_PARTICLES: removed old inactive sinks, took %g sec.\n", timediff(t0, t1));
}

void reinject_gas_mass(int isink)
{
  double rReinjection = 100. * PARSEC / All.UnitLength_in_cm;
  int local_n, total_n, *indices;

  /*collect all particles within RsMax from the sink*/
  while(1)
    {
      find_particles_within_a_sphere(SinkP[isink].Pos, rReinjection, &local_n, &total_n, &indices);
      if(total_n != 0)
        break;

      rReinjection *= 2.;
      myfree(indices);
    }

  double additional_mass = (SinkP[isink].Mass - SinkP[isink].StellarMass) / total_n;

  for(int i = 0; i < local_n; i++)
    {
      int idx = indices[i];

      double volume = P[idx].Mass / SphP[idx].Density;
      P[idx].Mass += additional_mass;

      SphP[idx].Density = P[idx].Mass / volume;

      SphP[idx].Momentum[0] = P[idx].Mass * P[idx].Vel[0];
      SphP[idx].Momentum[1] = P[idx].Mass * P[idx].Vel[1];
      SphP[idx].Momentum[2] = P[idx].Mass * P[idx].Vel[2];

      double abs_Mom2  = pow(SphP[idx].Momentum[0], 2) + pow(SphP[idx].Momentum[1], 2) + pow(SphP[idx].Momentum[2], 2);
      SphP[idx].Energy = SphP[idx].Utherm * P[idx].Mass + 0.5 * abs_Mom2 / P[idx].Mass;
#ifdef MHD
      SphP[idx].Energy +=
          0.5 * (SphP[idx].B[0] * SphP[idx].B[0] + SphP[idx].B[1] * SphP[idx].B[1] + SphP[idx].B[2] * SphP[idx].B[2]) * volume;
#endif

#ifdef MAXSCALARS
      for(int k = 0; k < N_Scalar; k++)
        *(MyFloat *)(((char *)(&SphP[idx])) + scalar_elements[k].offset_mass) =
            *(MyFloat *)(((char *)(&SphP[idx])) + scalar_elements[k].offset) * P[idx].Mass;
#endif
    }
  SinkP[isink].Mass -= additional_mass * total_n;

  if(isnan(SinkP[isink].Mass))
    printf("SinkP[isink].Mass=%g, additional_mass=%g, total_n=%d\n", SinkP[isink].Mass, additional_mass, total_n);

  int candidate_task = SinkP[isink].HomeTask;
  if(ThisTask == candidate_task)
    {
      int particle_index = SinkP[isink].Index;
      if(P[particle_index].ID != SinkP[isink].ID)
        terminate("Something is wrong, the sink is indexing the wrong particle\n");
      P[particle_index].Mass = SinkP[isink].Mass;
    }

  myfree(indices);
}
#endif /*SINK_PARTICLES_FEEDBACK*/
