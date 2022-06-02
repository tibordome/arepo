#include "../allvars.h"
#include "../proto.h"

/* PCC - 08.01.2014
         This function deals with the creation of new sinks.
*/
int create_sink_particles(void)
{
  int i, j;
  int potential_peak, still_a_peak;
  int candidate_task;
  double dx, dy, dz, dist;
  double dvx, dvy, dvz;
  double dax, day, daz;
  double egrav, esupport;
  double density_potential;
  double dist_to_sink_boundary, tenc, tff_candidate, vrad;
  double density_min_in_region;
  double rad_sink_closest;

  struct search_buffer
  {
    double density;
    int rank;
  } search_in, search_out;

  struct candidate_info
  {
    double Pos[3];
    double Vel[3];
    double Accel[3];
    long long ID;
    double Mass;
    double Utherm;
    double Potential;
    double Density;
#ifdef SINK_PARTICLES_VARIABLE_CREATION
    double FormR2;
    double creationdensity;
#ifdef SINK_PARTICLES_FEEDBACK
    double TargetMass;
#endif
#endif

  } candidate;

  struct sink_environment
  {
    double mass;
    double ekin;
    double etherm;
    double divv;
    double diva;
  } environment_in, environment_out;

  /* Initialize */
  int sink_creation_success = -1; /*-1 = no sink candidates
                                   * 0 = sink candidate failed
                                   * 1 = found valid new sink */
  double SinkAccretionRadiusSquared = SinkAccretionRadius * SinkAccretionRadius;

  /* Find our candidate sink for this timestep on THIS task
   * Candidate needs to be above the sink creation threshold
   * density, and be the most dense potential peak on this task */
  search_in.density   = -1;
  search_in.rank      = ThisTask;
  int candidate_index = -1;
  int candidate_idx   = -1;
  density_potential   = 1;
  int nfound          = 0;
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
      if(P[i].Type > 0 || P[i].ID == 0 || P[i].Mass == 0)
        continue;

#ifdef SINK_PARTICLES_VARIABLE_CREATION
      double SinkVariableCreationDensity;
      SinkVariableCreationDensity = variable_getcreationdensity(SphP[i].RefineTarget);
      /*      mpi_printf(" SinkVariableCreationDensity %g /n",SinkVariableCreationDensity );*/

      if(SphP[i].Density > SinkVariableCreationDensity && SphP[i].PotentialPeak == 1 && SphP[i].DivVel < 0.)

#else
      if(SphP[i].Density > SinkCreationDensityCurrent && SphP[i].PotentialPeak == 1 && SphP[i].DivVel < 0.)
#endif
        {
          nfound++;
          if(SphP[i].Density > search_in.density)
            {
              search_in.density = SphP[i].Density;
              candidate_index   = i;
              candidate_idx     = idx;
              density_potential = P[i].Potential;
            }
        }
    }
  int nfound_allTasks;
  MPI_Allreduce(&nfound, &nfound_allTasks, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(nfound_allTasks == 0)
    {
      mpi_printf("SINK_PARTICLES: No more candidates to check\n");
      sink_creation_success = -1;
      return sink_creation_success;
    }

  mpi_printf("SINK_PARTICLES: found %d sink candidates\n", nfound_allTasks);

  /* Check that the candidate isn't inside the accretion radius of another sink.
   * This can be done BEFORE Bcast of candidate, since we have the sink positions
   * from the SinkP struct. */
  if(NSinksAllTasks > 0)
    {
      if(nfound > 0)
        {
          for(i = 0; i < NSinksAllTasks; i++)
            {
              if(SinkP[i].ID == 0 && SinkP[i].Mass == 0)
                continue;
#ifdef SINK_PARTICLES_FEEDBACK
              if(SinkP[i].StellarMass > All.MaxStellarMassPerSink)
                continue;
#endif
              dx   = GRAVITY_NEAREST_X(SinkP[i].Pos[0] - P[candidate_index].Pos[0]);
              dy   = GRAVITY_NEAREST_Y(SinkP[i].Pos[1] - P[candidate_index].Pos[1]);
              dz   = GRAVITY_NEAREST_Z(SinkP[i].Pos[2] - P[candidate_index].Pos[2]);
              dist = dx * dx + dy * dy + dz * dz;
              if(i == 0)
                rad_sink_closest = dist;
              else
                {
                  if(dist < rad_sink_closest)
                    rad_sink_closest = dist;
                }
#ifdef SINK_PARTICLES_VARIABLE_ACC_RADIUS
              SinkAccretionRadiusSquared = SinkP[i].AccretionRadius * SinkP[i].AccretionRadius;
#endif
              if(dist < SinkAccretionRadiusSquared)
                {
                  SphP[candidate_index].PotentialPeak = 0;
                  search_in.density                   = -1;
                  /* this is never actually used in this case, as we break immedidately */
                  candidate_index = -1;
                  break;
                }
            }
        }
    }

  /* Our candidate now needs to be the DENSEST of the possible potential minima */
  MPI_Allreduce(&search_in, &search_out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

  /* If there are no candidates (no particle above rho_crit and in isolation), then we return. */
  if(search_out.density < 0)
    {
      mpi_printf("SINK_PARTICLES: near existing sinks \n");
      sink_creation_success = 0;
      return sink_creation_success;
    }

  /* BCast the information required for the energy checks */
  candidate_task = search_out.rank;
  if(ThisTask == candidate_task)
    {
      SphP[candidate_index].PotentialPeak = 0;

      candidate.Pos[0]    = P[candidate_index].Pos[0];
      candidate.Pos[1]    = P[candidate_index].Pos[1];
      candidate.Pos[2]    = P[candidate_index].Pos[2];
      candidate.Vel[0]    = P[candidate_index].Vel[0];
      candidate.Vel[1]    = P[candidate_index].Vel[1];
      candidate.Vel[2]    = P[candidate_index].Vel[2];
      candidate.Accel[0]  = P[candidate_index].GravAccel[0];
      candidate.Accel[1]  = P[candidate_index].GravAccel[1];
      candidate.Accel[2]  = P[candidate_index].GravAccel[2];
      candidate.ID        = (long long)P[candidate_index].ID;
      candidate.Mass      = P[candidate_index].Mass;
      candidate.Utherm    = SphP[candidate_index].Utherm;
      candidate.Potential = P[candidate_index].Potential;
      candidate.Density   = SphP[candidate_index].Density;
#ifdef SINK_PARTICLES_VARIABLE_CREATION
      candidate.creationdensity = variable_getcreationdensity(SphP[candidate_index].RefineTarget);
      candidate.FormR2          = variable_getformationr2(SphP[candidate_index].RefineTarget, candidate.creationdensity);
#ifdef SINK_PARTICLES_FEEDBACK
      candidate.TargetMass = SphP[candidate_index].RefineTarget;
#endif
#endif
    }

  MPI_Bcast(&candidate, sizeof(struct candidate_info), MPI_BYTE, candidate_task, MPI_COMM_WORLD);

  /* Now check that the region inside the candidate's interaction radius is
   * indeed both bound (energy check) and collapsing (divv & diva checks).
   * Once again, each CPU will run in parallel. Each task bins its gas particles
   * radially, centred on the candidate sink. We then do a MPI_All sum to get
   * the full radial profile. This is then used to make an energy profile, which
   * will be used to determine whether this gas cell can become a sink. If so,
   * the sink details are then added to the SinkP array immediately, and obviously
   * a new particle is created on the host task with P.Type = 5.
   * The newly-formed sink is then given a chance to accrete. */

  memset(&environment_in, 0, sizeof(struct sink_environment));

  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Mass == 0 || P[i].ID == 0)
        continue;
      dx   = GRAVITY_NEAREST_X(P[i].Pos[0] - candidate.Pos[0]);
      dy   = GRAVITY_NEAREST_Y(P[i].Pos[1] - candidate.Pos[1]);
      dz   = GRAVITY_NEAREST_Z(P[i].Pos[2] - candidate.Pos[2]);
      dist = dx * dx + dy * dy + dz * dz;
#ifdef SINK_PARTICLES_VARIABLE_CREATION
      if(dist > candidate.FormR2)
#else
      if(dist > SinkFormationRadiusSquared)
#endif
        continue;
      if(ThisTask == candidate_task && i == candidate_index)
        continue;
      dist = sqrt(dist);
      dvx  = P[i].Vel[0] - candidate.Vel[0];
      dvy  = P[i].Vel[1] - candidate.Vel[1];
      dvz  = P[i].Vel[2] - candidate.Vel[2];
      dax  = P[i].GravAccel[0] - candidate.Accel[0];
      day  = P[i].GravAccel[1] - candidate.Accel[1];
      daz  = P[i].GravAccel[2] - candidate.Accel[2];
      /* Now work out the partial sums (mass weighting!) for ekin, etherm
       * divv, and diva */
      environment_in.mass += P[i].Mass;
      environment_in.ekin += 0.5 * (dvx * dvx + dvy * dvy + dvz * dvz) * P[i].Mass;
      environment_in.etherm += SphP[i].Utherm * P[i].Mass;
      environment_in.divv += P[i].Mass * (dvx * dx + dvy * dy + dvz * dz) / dist;
      environment_in.diva += P[i].Mass * (dax * dx + day * dy + daz * dz) / dist;
    }

  /* get the full sum of the above values in the region */
  memset(&environment_out, 0, sizeof(struct sink_environment));
  MPI_Allreduce(&environment_in, &environment_out, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  /* Now get what we need from the sums. Need to add the candiate to the mass and
   * the thermal energy calculation. We don't need to add it to the kinetic energy though */
  environment_out.mass += candidate.Mass;
  environment_out.etherm += (candidate.Mass * candidate.Utherm);
  environment_out.divv /= environment_out.mass;
  environment_out.diva /= environment_out.mass;
#ifdef SINK_PARTICLES_VARIABLE_CREATION
  egrav = All.G * environment_out.mass * environment_out.mass / All.HubbleParam / candidate.FormR2;
#else
  egrav = All.G * environment_out.mass * environment_out.mass / All.HubbleParam / SinkFormationRadius;
#endif
  /* we use SinkFormationRadius and not SinkFormationRadiusCurrent here because egrav should be in proper units
   * the /h factor is to remove the h dependence from the second mass, so all energies scale with 1/h */
  if(All.ComovingIntegrationOn)
    {
      /* conversion of kinetc energy to proper units */
      esupport = environment_out.ekin / (All.Time * All.Time) + environment_out.etherm;
    }
  else
    {
      esupport = environment_out.ekin + environment_out.etherm;
    }

    /* Is this candidate the centre of a collapsing structure? If so, we have a sink,
     * if not, we return with the appropriate value of the success flag */

#ifdef SINK_PARTICLES_VARIABLE_CREATION
  if(candidate.Density > candidate.creationdensity * 100. && environment_out.divv <= 0 && environment_out.diva <= 0)
#else
  if(candidate.Density > SinkTestsGiveupDensityCurrent && environment_out.divv <= 0 && environment_out.diva <= 0)
#endif
    {
      mpi_printf(
          "SINK_PARTICLES: Candidate density of %g is now very dense, but converging... Forcing sink formation regardless of energy "
          "checks! \n",
          candidate.Density);
    }
  else
    {
      if(egrav <= 2. * esupport || environment_out.divv >= 0 || environment_out.diva >= 0)
        {
          mpi_printf("SINK_PARTICLES: Candidate failed the energy checks: egrav %g esupport %g divv %g diva %g \n", egrav, esupport,
                     environment_out.divv, environment_out.diva);
          mpi_printf(
              "SINK_PARTICLES: Number density cgs %g Temperature %g average Temperature %g \n",
              candidate.Density * All.UnitDensity_in_cgs * All.cf_a3inv * (All.HubbleParam * All.HubbleParam) / PROTONMASS / 1.22,
              candidate.Utherm * All.UnitEnergy_in_cgs / All.UnitMass_in_g * 2.0 / 3.0 * 1.22 * PROTONMASS / BOLTZMANN,
              environment_out.etherm / environment_out.mass * All.UnitEnergy_in_cgs / All.UnitMass_in_g * 2.0 / 3.0 * 1.22 *
                  PROTONMASS / BOLTZMANN);
#ifndef SINK_PARTICLES_FORCE_FORMATION
          sink_creation_success = 0;
          return sink_creation_success;
#else
#ifdef SINK_PARTICLES_VARIABLE_CREATION
          if(candidate.Density > candidate.creationdensity * 10000.0)
#else
          if(candidate.Density > SinkTestsGiveupDensityCurrent * 100.0)
#endif
            {
              warn(
                  "SINK_PARTICLES: Candidate density of %g is now super dense but not converging....... Forcing sink formation "
                  "regardless of all checks! \n",
                  candidate.Density);
            }
          else
            {
              sink_creation_success = 0;
              return sink_creation_success;
            }
#endif
        }
    }

#ifdef SINK_PARTICLE_FREE_FALL_TEST
  /* Final test! Check that the candidate isn't going to fall into another sink's accretion radius
   * before its is able to collapse locally: i.e. tff(candidate) < interaction time
   * This test is useful if you have discs forming around existing sinks, and there's nothing to
   * stabilise the inner disc, or if you have very dense clusters. XXX Untested! */
  if(NSinksAllTasks > 0 && candidate.Density < SinkTestsGiveupDensityCurrent)
    {
#ifdef DEBUG_SINK_PARTICLES
      mpi_printf("SINK_PARTICLES--FORM: checking the free fall time criterion \n");
#endif
      for(i = 0; i < NSinksAllTasks; i++)
        {
          if(SinkP[i].ID == 0 && SinkP[i].Mass == 0)
            continue;

          dx   = GRAVITY_NEAREST_X(candidate.Pos[0] - SinkP[i].Pos[0]);
          dy   = GRAVITY_NEAREST_Y(candidate.Pos[1] - SinkP[i].Pos[1]);
          dz   = GRAVITY_NEAREST_Z(candidate.Pos[2] - SinkP[i].Pos[2]);
          dist = sqrt(dx * dx + dy * dy + dz * dz);
#ifdef SINK_PARTICLES_VARIABLE_ACC_RADIUS
          SinkAccretionRadius = SinkP[i].AccretionRadius;
#endif
          dist_to_sink_boundary = dist - SinkAccretionRadius;
          if(dist_to_sink_boundary > 0)
            {
              dvx  = candidate.Vel[0] - SinkP[i].Vel[0];
              dvy  = candidate.Vel[1] - SinkP[i].Vel[1];
              dvz  = candidate.Vel[2] - SinkP[i].Vel[2];
              vrad = (dvx * dx + dvy * dy + dvz * dz) / dist;
#ifdef DEBUG_SINK_PARTICLES
              mpi_printf("SINK_PARTICLES--FORM: free-fall check dist %g dist to sink %g vrad %g \n", dist, dist_to_sink_boundary,
                         vrad);
#endif
              if(vrad <= 0)
                {
                  tenc          = dist_to_sink_boundary / (-vrad);
                  tff_candidate = sqrt(3 * M_PI / 32.0 / All.G / candidate.Density);
                  if(All.ComovingIntegrationOn)
                    tenc *= sqrt(All.Time);
#ifdef DEBUG_SINK_PARTICLES
                  mpi_printf("SINK_PARTICLES--FORM: free fall check tff %g tenc %g \n", tff_candidate, tenc);
#endif
                  if(tenc < tff_candidate)
                    {
#ifdef DEBUG_SINK_PARTICLES
                      mpi_printf("SINK_PARTICLES--FORM: candidate will fall into sink before formation on Task %d \n", ThisTask);
                      mpi_printf("SINK_PARTICLES--FORM: dist_to_sink_boundary %g tenc %g tff_candidate %g \n", dist_to_sink_boundary,
                                 tenc, tff_candidate);
#endif
                      sink_creation_success = -4;
                      return sink_creation_success;
                    }
                }
            }
          else
            {
              terminate("SINK_PARTICLES--FORM: candidate is inside another sink on Task %d! This shouldn't happen! Aborting",
                        ThisTask);
            }
        }
    }
#endif /* SINK_PARTICLE_FREE_FALL_TEST */

  /* Joy! Bang a drum... */
  mpi_printf("SINK_PARTICLES: SINK CREATION SUCCESSFUL!\n");
  mpi_printf("SINK_PARTICLES: CREATION egrav %g esupport %g divv %g diva %g \n", egrav, esupport, environment_out.divv,
             environment_out.diva);
  mpi_printf("SINK_PARTICLES: CREATION formation mass %g cell mass %g \n", environment_out.mass, candidate.Mass);
  mpi_printf("SINK_PARTICLES: CREATION number density cgs %g Temperature %g \n",
             candidate.Density * All.UnitDensity_in_cgs * All.cf_a3inv * (All.HubbleParam * All.HubbleParam) / PROTONMASS / 1.22,
             candidate.Utherm * All.UnitEnergy_in_cgs / All.UnitMass_in_g * 2.0 / 3.0 * 1.22 * PROTONMASS / BOLTZMANN);
#ifdef SINK_PARTICLES_VARIABLE_CREATION
  mpi_printf("SINK_PARTICLES: CREATION density %g radius %g \n", candidate.creationdensity, pow(candidate.FormR2, 0.5));
#endif

  sink_creation_success = 1;

  /* Update the SinkP array. The ID and Index are only needed on the host Task, so they
   * don't need to be Bcast to the other Tasks. */
  for(i = 0; i < 3; i++)
    SinkP[NSinksAllTasks].Pos[i] = candidate.Pos[i];
  for(i = 0; i < 3; i++)
    SinkP[NSinksAllTasks].Vel[i] = candidate.Vel[i];
  for(i = 0; i < 3; i++)
    SinkP[NSinksAllTasks].Accel[i] = candidate.Accel[i];
  SinkP[NSinksAllTasks].Mass          = candidate.Mass;
  SinkP[NSinksAllTasks].FormationMass = environment_out.mass; /* Used in early accretion */
  SinkP[NSinksAllTasks].FormationTime = All.Time;
  SinkP[NSinksAllTasks].HomeTask      = candidate_task;
  SinkP[NSinksAllTasks].ID            = candidate.ID;
  if(ThisTask == candidate_task)
    SinkP[NSinksAllTasks].Index = candidate_index;
  else
    SinkP[NSinksAllTasks].Index = -1;
  SinkP[NSinksAllTasks].FormationOrder = NSinksAllTasks + 1;

#ifdef SINK_PARTICLES_VARIABLE_CREATION
  SinkP[NSinksAllTasks].CreationDensity = candidate.creationdensity;
#ifdef SINK_PARTICLES_FEEDBACK
  double feedback_density     = environment_out.mass / (4. / 3. * M_PI * pow(candidate.FormR2, 1.5));
  double SFeff_estimate       = variable_geteff(candidate.TargetMass, feedback_density, candidate.creationdensity);
  SinkP[NSinksAllTasks].SFeff = SFeff_estimate;
  mpi_printf("SinkVAR: filled structure with SFeff %g, env_mass %g\n", SinkP[NSinksAllTasks].SFeff, environment_out.mass);
#endif
#endif /*SINK_PARTICLES_VARIABLE_CREATION*/

#ifdef SINK_PARTICLES_FEEDBACK
  SinkP[NSinksAllTasks].N_sne = 0; /*initialize to zero*/
  memset(SinkP[NSinksAllTasks].MassStillToConvert, 0, MAXACCRETIONEVENTS * sizeof(double));

#ifdef SINK_PARTICLES_VARIABLE_CREATION
  SinkP[NSinksAllTasks].StellarMass = candidate.Mass * SinkP[NSinksAllTasks].SFeff;
#else
  SinkP[NSinksAllTasks].StellarMass = candidate.Mass * All.SINKStarFormationEfficiency;
#endif

  assign_sne(NSinksAllTasks, SinkP[NSinksAllTasks].StellarMass, candidate_task);
#endif /*SINK_PARTICLES_FEEDBACK*/

  NSinksAllTasks++;
  NSinksThisTask++;

  /* Convert cell to particle */
  add_sink_to_particle_structure(candidate_idx, candidate_index, candidate_task);

  /* Keep track of which task was responsible for forming this new sink... Might be useful! */
  LastTaskNewSink = candidate_task;

  return sink_creation_success;
}

/* This function turns the gas cell ('particle' Type 0) into a sink particle (Type 5).
 * The code is based on the technique in star_formation.c (functions make_star and
 * convert_cell_into_star). Unlike the old version, where we made a new Type 5 particle,
 * and then ripped the Type 0 particle out of the mesh, here we simply change its type!
 * We also free up the memory that the particle is taking in the mesh, and remove it from
 * the list of particle in the hydro time bin.
 */
void add_sink_to_particle_structure(int idx, int candidate_index, int candidate_task)
{
  mpi_printf("SINK_PARTICLES: Adding particle to P array and cleaning up! \n");

  if(ThisTask == candidate_task)
    {
      P[candidate_index].Type          = 5;
      P[candidate_index].SofteningType = All.SofteningTypeOfPartType[P[candidate_index].Type];
#ifdef INDIVIDUAL_GRAVITY_SOFTENING
      if(((1 << P[candidate_index].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
        {
          P[candidate_index].SofteningType = get_softening_type_from_mass(P[candidate_index].Mass);
          printf("SINK_PARTICLES: Doing the weird mass softening thing... \n");
        }
#endif

      printf("SINK_PARTICLES: Task %d converting particle with ID %d to sink. the idx is %d \n", ThisTask, P[candidate_index].ID, idx);

      /* Can remove the candidate gas cell from the mesh now */
      timebin_remove_particle(&TimeBinsHydro, idx, P[candidate_index].TimeBinHydro);
#ifdef VORONOI_DYNAMIC_UPDATE
      voronoi_remove_connection(candidate_index);
#endif
      printf("SINK_PARTICLES: Finshed the update of the P-array \n");
    }
  SinksFormedSinceLastDomain++;
}

#ifdef SINK_PARTICLES_VARIABLE_CREATION

/* A simple step function that increases the sink creation density in highly refined regions */

double variable_getcreationdensity(double TargetMass)
{
  double creation_density;

  creation_density = 1.0e+8; /* ~1000 */

  if(TargetMass < 1000. && TargetMass >= 100.)
    creation_density = 1.0e+9; /* ~1e4 */
  else if(TargetMass < 100. && TargetMass >= 10.)
    creation_density = 1.0e+10; /* ~1.e5 */
  else if(TargetMass < 10.)
    creation_density = 1.0e+11; /* ~1.e6 */

  return creation_density;
}

/* For the unrefined cells we use the value from the parameter file - otherwise use the cell radius for a volume containing 8 cells are
 * our creation density */

double variable_getformationr2(double TargetMass, double CreationDensity)
{
  double var_r2, radius;

  var_r2 = SinkFormationRadius * SinkFormationRadius;
  if(TargetMass < 1000.)
    {
      radius = pow(3.0 * TargetMass * 8. / 4. / M_PI / CreationDensity, 0.3333333);
      var_r2 = radius * radius;
    }

  return var_r2;
}

#ifdef SINK_PARTICLES_FEEDBACK
double variable_geteff(double TargetMass, double Density, double CreationDensity)
{
  double eff;

  if(TargetMass <= 10.)
    {
      eff = All.SINKStarFormationEfficiency;
      mpi_printf("SINK_PARTICLES: TARGET REACHED Assigning eff of %g, %g \n", eff, TargetMass);
    }
  else
    {
      double rho_eff, tlife, tff;

      if(Density > CreationDensity)
        {
          rho_eff = CreationDensity;
        }
      else
        {
          rho_eff = Density;
        }

      tlife = 1.e7 * SEC_PER_YEAR / All.UnitTime_in_s; /* assumed cloud lifetime */
      tff   = sqrt(3 * M_PI / 32.0 / All.G / rho_eff);
      eff   = 0.01 * tlife / tff;

      mpi_printf("SINK_PARTICLES: Assigning eff of %g, %g, %g, %g \n", eff, TargetMass, Density, tff);
    }

  return eff;
}
#endif

#endif
