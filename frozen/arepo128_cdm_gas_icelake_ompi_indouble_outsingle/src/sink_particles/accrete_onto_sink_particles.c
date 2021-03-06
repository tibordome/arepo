#include "../allvars.h"
#include "../proto.h"

/* PCC - 06.01.2014
	Algorithm to handle the accertion onto the sink particles that are
	present in the domain. Since all Tasks now know the properties of
	all the sinks in the simulation, this process can run in isolation
	on each Task. 
*/
double accrete_onto_sink_particles(void)
{
  int idx, icell, isink, jsink;
  int iloc;
  int i_active, i_pass_density;
  int num_accreted_this_timestep;
  int NumberSinksAccretedThisStep;
  int isink_min; 
  double cell_radius, dist;
  double dx, dy, dz;
  double dvx, dvy, dvz;
  double dax, day, daz;
  double vrad, arad, dv;
  double sink_mass;
  double energy_total, energy_min, energy;
  double egrav, ekin, etherm;
  double dt;
  double length_this_timestep;
  double mass_sum;
  double MassAccretedViaCellMunch;
  double MassAccretedViaBondiHoyle;
  double TotalMassAccreted;
  double SinkAccretionRadiusSquared;
  double mass_accreted;
#ifdef SGCHEM_ACCRETION_LUMINOSITY
  double Msmooth;
  double sinkn, tff;
#endif
  struct accretion_buffer
    {
      double com[3];
      double mom[3];
      double mass;
    } *acc_in, *acc_out;

#if (defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)) || defined(GRAVITY_TALLBOX)
  double xtmp, ytmp, ztmp;
#endif


  /* Initialize variables
  */
  num_accreted_this_timestep = 0;
  TotalMassAccreted = 0;
  MassAccretedViaCellMunch = 0;
  MassAccretedViaBondiHoyle = 0;
  SinkAccretionRadiusSquared = SinkAccretionRadius*SinkAccretionRadius;

  acc_in = mymalloc("acc_in", NSinksAllTasks * sizeof(struct accretion_buffer));
  memset(acc_in, 0, NSinksAllTasks * sizeof(struct accretion_buffer));

  acc_out = mymalloc("acc_out", NSinksAllTasks * sizeof(struct accretion_buffer));
  memset(acc_out, 0, NSinksAllTasks * sizeof(struct accretion_buffer));  

  /* The main loop around all gas cells on this Task.
     IMPORTANT: Both the gas cell AND the sink particle need to be on
     the current timestep if we are to accrete from the cell. Since the
     sinks are forced to be on the global minimumm timestep, this amounts
     to checking whether the gas is active.
     Note that we have a problem that the sinks might MISS cells, if they 
     are on very different timesteps. Especially true if the cells are larger than the
     sink have a density lower then the creation density -- the sinks could pass
     straight through the cell, without seeing them. Even Bondi-Hoyle accretion would be missed 
     Easiest way around this would be to simply enforce refinement and lower timesteps
     near the sink. This would capture any B-H type accretion with the standard algorithm.
  */
  i_active = 0;
  i_pass_density=0;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
#ifndef SINK_PARTICLES_SKIM_CELL_MASS
      if(num_accreted_this_timestep > 50)
        break; // limit the destruction to the mesh in one timestep! 
#endif
      icell = TimeBinsHydro.ActiveParticleList[idx];
      if(icell < 0)
        continue;
      i_active++;
      if(!(TimeBinSynchronized[P[icell].TimeBinHydro])) 
        continue;
      if(SphP[icell].Density < SinkCreationDensityCodeUnits)
        continue;
      i_pass_density++;
      cell_radius = pow(SphP[icell].Volume*3./(4.*M_PI), 1./3.);
      /* Loop around all sinks in domain 
      */
      energy_min = 1e33;
      isink_min = -1;
      for(isink = 0; isink < NSinksAllTasks; isink++)
        {
          dx = GRAVITY_NEAREST_X(P[icell].Pos[0] - SinkP[isink].Pos[0]);
          dy = GRAVITY_NEAREST_Y(P[icell].Pos[1] - SinkP[isink].Pos[1]);
          dz = GRAVITY_NEAREST_Z(P[icell].Pos[2] - SinkP[isink].Pos[2]);
          dist = dx * dx  +  dy * dy  +  dz * dz;
#ifdef SINK_PARTICLES_VARIABLE_ACC_RADIUS
          SinkAccretionRadiusSquared =  SinkP[isink].AccretionRadius * SinkP[isink].AccretionRadius;
#endif
          if(dist < SinkAccretionRadiusSquared) 
            {
              /* Minimum requirement for accretion passed.
                 Now do the more sophisticated checks!
              */
              dist = sqrt(dist);

              dvx = P[icell].Vel[0] - SinkP[isink].Vel[0];
              dvy = P[icell].Vel[1] - SinkP[isink].Vel[1];
              dvz = P[icell].Vel[2] - SinkP[isink].Vel[2];
              vrad = (dvx * dx  +  dvy * dy  +  dvz * dz) / dist;

              dax = P[icell].GravAccel[0] - SinkP[isink].Accel[0];
              day = P[icell].GravAccel[1] - SinkP[isink].Accel[1];
              daz = P[icell].GravAccel[2] - SinkP[isink].Accel[2];
              arad = (dax * dx  +  day * dy  +  daz * dz) / dist;

              dv = dvx * dvx  +  dvy * dvy  +  dvz * dvz;
              ekin = 0.5 * P[icell].Mass * dv;
              /* The following is to ensure that newly created sinks are able to
                 to accrete the material inside their formation radius upon formation
              */
              if(SinkP[isink].Mass > SinkP[isink].FormationMass)
                sink_mass = SinkP[isink].Mass;
              else
                sink_mass = SinkP[isink].FormationMass;

              egrav = P[icell].Mass * sink_mass / dist;
              etherm = P[icell].Mass * SphP[icell].Utherm;
              energy_total = ekin + etherm - egrav;

              if(energy_total < energy_min)
                {
                  if( ((energy_total < 0)  &&  (vrad < 0)  &&  (arad < 0)) || (SphP[icell].Density > SinkTestsGiveupDensityCodeUnits) )
                    {
                      energy_min = energy_total;
                      isink_min = isink;
                    }
                }
            }
        }

      /* Now check if we can accrete this cell
      */
      if(isink_min > -1) 
        { 
          isink = isink_min;
#ifdef SINK_PARTICLES_SKIM_CELL_MASS
          /* We are skimming mass from the cell. Try to keep the cell density 
          ~ sink creation density. But limit the skimmed mass!
          */
          mass_accreted = P[icell].Mass - (SinkCreationDensityCodeUnits * P[icell].Mass/SphP[icell].Density);
          if(mass_accreted > 0.9 * P[icell].Mass)
            mass_accreted = 0.9 * P[icell].Mass;
          if(mass_accreted < 0)
            {
              printf("SINK WARNING: Negative accreted mass. cell mass %g required mass %g SphP.vol %g vol from m/rho %g \n", P[icell].Mass, SinkCreationDensityCodeUnits * SphP[icell].Volume, SphP[icell].Volume, P[icell].Mass/SphP[icell].Density);
              mass_accreted = 0;
            }
#else
          /* We are just eating the entire cell, so take all the mass */
          mass_accreted = P[icell].Mass;
#endif 
                  
          /* Load the momentum information into the communication buffer. Although the
             position looks a little odd below, it accounts for the periodic box. This isn't
             required for the velocities.
          */
          dx = GRAVITY_NEAREST_X(P[icell].Pos[0] - SinkP[isink].Pos[0]);
          dy = GRAVITY_NEAREST_Y(P[icell].Pos[1] - SinkP[isink].Pos[1]);
          dz = GRAVITY_NEAREST_Z(P[icell].Pos[2] - SinkP[isink].Pos[2]);
          acc_in[isink].com[0] += (SinkP[isink].Pos[0] + dx) * mass_accreted;
          acc_in[isink].com[1] += (SinkP[isink].Pos[1] + dy) * mass_accreted;
          acc_in[isink].com[2] += (SinkP[isink].Pos[2] + dz) * mass_accreted;
          acc_in[isink].mom[0] += P[icell].Vel[0] * mass_accreted;
          acc_in[isink].mom[1] += P[icell].Vel[1] * mass_accreted;
          acc_in[isink].mom[2] += P[icell].Vel[2] * mass_accreted;
          acc_in[isink].mass += mass_accreted;


          /* Now accrete from the cell
          */
          num_accreted_this_timestep++;
#ifdef SINK_PARTICLES_SKIM_CELL_MASS
          /* Skim the accreted mass from the cell and update any variables that depend on mass 
          */
          double fac = (P[icell].Mass - mass_accreted)/P[icell].Mass;
          if(fac > 1)
            fac = 1;
          P[icell].Mass *= fac;
          SphP[icell].Density *= fac;
          SphP[icell].Momentum[0] *= fac;
          SphP[icell].Momentum[1] *= fac;
          SphP[icell].Momentum[2] *= fac;
          SphP[icell].Energy *= fac;
          SphP[icell].Pressure *= fac;
          SphP[icell].OldMass *= fac;
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
          //SphP[icell].A = GAMMA_MINUS1 * SphP[icell].Utherm / pow(SphP[icell].Density, GAMMA_MINUS1);
          //SphP[icell].Entropy = log(SphP[icell].A) * P[icell].Mass;
          // XXX the method below is what the blackhole code does, but I don't think it's correct...
          SphP[icell].Entropy *= fac;
#endif
#ifdef DEBUG_SINK_PARTICLES
          printf("SINK ACCRETE:  %g mass %g skimmass %g fac %g \n", SphP[icell].Density, P[icell].Mass, mass_accreted, fac);
#endif          
#ifdef SGCHEM
          for(int k = 0; k < SGCHEM_NUM_ADVECTED_SPECIES; k++)
            SphP[icell].MassTracAbund[k] *= fac;
#endif
#else
          /* Eat the whole cell... nom nom...
          */
          MassAccretedViaCellMunch += P[icell].Mass;
          P[icell].Mass = 0;
          P[icell].ID = 0;
          P[icell].Vel[0] = 0;
          P[icell].Vel[1] = 0;
          P[icell].Vel[2] = 0;
          timebin_remove_particle(&TimeBinsHydro, idx, P[icell].TimeBinHydro);
#ifdef VORONOI_DYNAMIC_UPDATE
          voronoi_remove_connection(icell);
#endif
#endif
        }
    }

  /* Now sum up the contributions to the sinks from all Tasks.
  */
  MPI_Allreduce(acc_in, acc_out, 7 * NSinksAllTasks, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


  /* Loop around the sinks and find those on this Task. Update their
     properties in the real (P.) particle structure.
  */
  for(isink = 0; isink < NSinksAllTasks; isink++)
    {
      /* Update the SinkP array
      */
      mass_sum = SinkP[isink].Mass + acc_out[isink].mass;
#ifdef SGCHEM_ACCRETION_LUMINOSITY
      if (All.SinkAccretionRateSmoothingMass > 0) {
      /* Smooth accretion rate over time taken to accrete number of solar masses specified by
       * SinkAccretionRateSmoothingMass parameter. At early times, when less than this mass
       * has been accreted, average over time since sink formation. If SinkAccretionRateSmoothingMass == 0,
       * don't smooth.
       *
       * Note: to recover behaviour of Smith et al (2012) approach, set SinkAccretionRateSmoothingMass = 3.0
       */
        Msmooth = All.SinkAccretionRateSmoothingMass * SOLAR_MASS / All.UnitMass_in_g;
        if (mass_sum < Msmooth) {
	  SinkP[isink].TimeOld = SinkP[isink].FormationTime;
	  double sink_age = All.Time - SinkP[isink].FormationTime;
	  /* At very early times, treat the sink as if its age were the free-fall time at the sink creation density.
	   * Otherwise, this calculation of the accretion rate diverges as we approach SinkP[isink].FormationTime
	   */
	  sinkn = SinkCreationDensityCodeUnits * All.UnitDensity_in_cgs / ((1.0 + 4.0 * ABHE) * PROTONMASS);
	  tff = 1.5e15 / sqrt(sinkn) / All.UnitTime_in_s;
	  if (sink_age < tff) {
	    sink_age = tff;
	  }
	  SinkP[isink].AccretionRate = mass_sum / sink_age;
        }
        else {
	  dt = mass_sum / SinkP[isink].AccretionRate;
	  SinkP[isink].TimeOld += dt;
	  if (SinkP[isink].TimeOld > All.Time) {
	    printf("Warning: problem with smoothed accretion rate! mass_sum = %e, dt = %e, TimeOld = %e, isink = %d\n",
                   mass_sum, dt, SinkP[isink].TimeOld, isink);
	  }
	  SinkP[isink].AccretionRate = Msmooth / (All.Time - SinkP[isink].TimeOld);
        }
      }
      else {
	if (All.Time == SinkP[isink].FormationTime) {
          sinkn = SinkCreationDensityCodeUnits * All.UnitDensity_in_cgs /((1.0 + 4.0 * ABHE) *PROTONMASS);
          dt = 1.5e15 / sqrt(sinkn) / All.UnitTime_in_s;
	}
	else {
          dt = All.Time - SinkP[isink].TimeOld;
	}
	if (dt > 0.0) {
          SinkP[isink].AccretionRate = acc_out[isink].mass / dt;
	}
	SinkP[isink].TimeOld = All.Time;
      }
#endif
      SinkP[isink].Pos[0] = (SinkP[isink].Mass*SinkP[isink].Pos[0] + acc_out[isink].com[0])/mass_sum;
      SinkP[isink].Pos[1] = (SinkP[isink].Mass*SinkP[isink].Pos[1] + acc_out[isink].com[1])/mass_sum;
      SinkP[isink].Pos[2] = (SinkP[isink].Mass*SinkP[isink].Pos[2] + acc_out[isink].com[2])/mass_sum;
      SinkP[isink].Vel[0] = (SinkP[isink].Mass*SinkP[isink].Vel[0] + acc_out[isink].mom[0])/mass_sum;
      SinkP[isink].Vel[1] = (SinkP[isink].Mass*SinkP[isink].Vel[1] + acc_out[isink].mom[1])/mass_sum;
      SinkP[isink].Vel[2] = (SinkP[isink].Mass*SinkP[isink].Vel[2] + acc_out[isink].mom[2])/mass_sum;
      SinkP[isink].Mass += acc_out[isink].mass;
#ifdef SINK_PARTICLES_FEEDBACK
      if (acc_out[isink].mass != 0.)
          assign_sne(isink, acc_out[isink].mass, SinkP[isink].HomeTask);
#endif

      /* Update the real particle array
      */
      if(SinkP[isink].HomeTask == ThisTask)
        {
          iloc = SinkP[isink].Index;
          P[iloc].Pos[0] = SinkP[isink].Pos[0];
          P[iloc].Pos[1] = SinkP[isink].Pos[1];
          P[iloc].Pos[2] = SinkP[isink].Pos[2];
          P[iloc].Vel[0] = SinkP[isink].Vel[0];
          P[iloc].Vel[1] = SinkP[isink].Vel[1];
          P[iloc].Vel[2] = SinkP[isink].Vel[2];
          P[iloc].Mass = SinkP[isink].Mass;
        } 
    }


  /* Get the statistics on what has been accreted this timestep.
  */
  TotalMassAccreted = 0;
  NumberSinksAccretedThisStep = 0;
  for(isink = 0; isink < NSinksAllTasks; isink++)
    {
      TotalMassAccreted += acc_out[isink].mass;
      if(acc_out[isink].mass > 0)
        NumberSinksAccretedThisStep++;
    }


  /* purge accreted cells from list of active gravity cells. They've
     already been removed from the Hydro list above.
     Don't need this if SINK_PARTICLES_SKIM_CELL_MASS is on. 
  */
#ifndef SINK_PARTICLES_SKIM_CELL_MASS
  if(num_accreted_this_timestep > 0)
    {
      printf("SINK_PARTICLES: Accrete on Task %d -- cleaning the gravity timebins for %d particles \n", ThisTask, num_accreted_this_timestep);
      timebin_cleanup_list_of_active_particles(&TimeBinsGravity);
    }
#endif

  /* Clean the memory
  */
  myfree(acc_out);
  myfree(acc_in);


  /* Should be finished! Return the mass accreted during this step
  */
  mpi_printf("SINK_PARTICLES: ACC -- %d sink(s) accreted mass this timestep \n", NumberSinksAccretedThisStep++);
  return(TotalMassAccreted);
}

