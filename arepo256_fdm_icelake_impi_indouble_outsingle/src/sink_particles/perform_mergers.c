#include "../allvars.h"
#include "../proto.h"

void perform_mergers(void)
{
  mpi_printf("SINK_MERGERS NSinks = %d \n", NSinksAllTasks);

#ifdef SINK_MERGERS_DEBUG
  mpi_printf("SINK_MERGERS Masses");
  for(int j = 0; j < NSinksAllTasks; j++)
    {
      mpi_printf(" %g ", SinkP[j].Mass);
    }
  mpi_printf("\n");
#endif

  MPI_Barrier(MPI_COMM_WORLD);
  if(NSinksAllTasks < 2)
    mpi_printf("SINK_MERGERS not enough sinks for merger routine \n");
  else
    {
      mpi_printf("SINK_MERGERS initialising \n");
      int merged[NSinksAllTasks]; /* keep track of deleted sinks, -1 means not merged */
      int i;
      int j;
      int z;
      double dx;
      double dy;
      double dz;
      double dist;
      double rad_sink_closest;
      double dvx;
      double dvy;
      double dvz;
      double vrad;
      double dax;
      double day;
      double daz;
      double arad;
      double dv;
      double sink_mass_i;
      double sink_mass_j;
      double ekin;
      double egrav;
      double SinkAccretionRadiusSquared = SinkAccretionRadius * SinkAccretionRadius;
      int index_merge;
      int index_survive;

      /*initialise merged array*/
      for(j = 0; j < NSinksAllTasks; j++)
        {
          merged[j] = -1;
        }

      /*for each sink, loop through all other sinks to check:*/
      /* if they lie within thier accretion radius, are moving towards eachother, and are bound to eachother*/
      /*if all 3 criteria are met then the smaller sink will be marked to merge  with the more massive sink in the merged array*/
      mpi_printf("SINK_MERGERS searching for candidates \n");
      for(j = 0; j < NSinksAllTasks; j++)
        {
          mpi_printf("SINK_MERGERS j = %d \n", j);
          if(merged[j] == -1)
            /* check it still exists */
            {
              for(i = 0; i < NSinksAllTasks; i++)
                {
                  mpi_printf("SINK_MERGERS i = %d \n", i);
                  if(merged[i] == -1)
                    /* can’t merge with a sink that doesn’t exist anymore */
                    {
                      if(SinkP[i].ID != SinkP[j].ID)
                        /* can’t merge with itself */
                        {
                          if(merged[j] == -1)
                            /*the current j sink could have been merged in the last i loop if mass[i]>mass[j]*/
                            {
#ifdef SINK_MERGERS_DEBUG
                              mpi_printf("SINK_MERGERS j = %d, i=%d \n", j, i);
                              mpi_printf("SINK_MERGERS merger array ");
                              for(z = 0; z < NSinksAllTasks; z++)
                                {
                                  mpi_printf(" %d ", merged[z]);
                                }
                              mpi_printf("\n");
#endif
                              /* Check separation*/
                              dx   = GRAVITY_NEAREST_X(SinkP[i].Pos[0] - SinkP[j].Pos[0]);
                              dy   = GRAVITY_NEAREST_Y(SinkP[i].Pos[1] - SinkP[j].Pos[1]);
                              dz   = GRAVITY_NEAREST_Z(SinkP[i].Pos[2] - SinkP[j].Pos[2]);
                              dist = dx * dx + dy * dy + dz * dz;
#ifdef SINK_MERGERS_DEBUG
                              mpi_printf("SINK_MERGERS distance from sink %d to %d = %g, r_accrete = %g \n", j, i, dist,
                                         SinkAccretionRadiusSquared);
#endif
                              if((i == 1) && (j == 0))
                                rad_sink_closest = dist;
                              else
                                {
                                  if(dist < rad_sink_closest)
                                    rad_sink_closest = dist;
                                }

#ifdef SINK_PARTICLES_VARIABLE_ACC_RADIUS
                              SinkAccretionRadiusSquared = SinkP[i].AccretionRadius * SinkP[i].AccretionRadius;
#endif

                              if(dist > SinkAccretionRadiusSquared) /*close enough to merge*/

#ifdef SINK_MERGERS_DEBUG
                                mpi_printf("SINK_MERGERS not close enough to merge \n");
#else
                                ;
#endif
                              else
                                {
#ifdef SINK_MERGERS_DEBUG
                                  mpi_printf("SINK_MERGERS  pos (%g,%g,%g) and (%g,%g,%g) \n", SinkP[j].Pos[0], SinkP[j].Pos[1],
                                             SinkP[j].Pos[2], SinkP[i].Pos[0], SinkP[i].Pos[1], SinkP[i].Pos[2]);
#endif
                                  /* Check divergence*/
                                  dist = sqrt(dist);

                                  dvx  = SinkP[i].Vel[0] - SinkP[j].Vel[0];
                                  dvy  = SinkP[i].Vel[1] - SinkP[j].Vel[1];
                                  dvz  = SinkP[i].Vel[2] - SinkP[j].Vel[2];
                                  vrad = (dvx * dx + dvy * dy + dvz * dz) / dist;

                                  dax  = SinkP[i].Accel[0] - SinkP[j].Accel[0];
                                  day  = SinkP[i].Accel[1] - SinkP[j].Accel[1];
                                  daz  = SinkP[i].Accel[2] - SinkP[j].Accel[2];
                                  arad = (dax * dx + day * dy + daz * dz) / dist;

                                  if((vrad > 0) || (arad > 0)) /*moving towards each other*/
#ifdef SINK_MERGERS_DEBUG
                                    mpi_printf("SINK_MERGERS not moving towards eachother \n");
#else
                                    ;
#endif
                                  else
                                    {
                                      /* Check if bound*/
                                      dv = dvx * dvx + dvy * dvy + dvz * dvz;
                                      if(SinkP[i].FormationTime == All.Time)
                                        sink_mass_i = SinkP[i].FormationMass;
                                      else
                                        sink_mass_i = SinkP[i].Mass;

                                      if(SinkP[j].FormationTime == All.Time)
                                        sink_mass_j = SinkP[j].FormationMass;
                                      else
                                        sink_mass_j = SinkP[j].Mass;

                                      ekin  = 0.5 * sink_mass_i * dv; /*energies*/
                                      egrav = All.G * sink_mass_i * sink_mass_j / dist;
                                      if(All.ComovingIntegrationOn) /*Some cosmo stuff idk*/
                                        {
                                          /*converting energies to physical units*/
                                          ekin /= (All.Time * All.Time);
                                          egrav /= All.Time;
                                        }

                                      int e_total = ekin - egrav;
                                      if(e_total > 0) /*bound*/
#ifdef SINK_MERGERS_DEBUG
                                        mpi_printf("SINK_MERGERS not bound \n");
#else
                                        ;
#endif
                                      else
                                        {
#ifdef SINK_MERGERS_DEBUG
                                          mpi_printf("SINK_MERGERS a sink has passed the test \n");
                                          mpi_printf("SINK_MERGERS closest sinks at d = %g while r_accrete = %g \n", rad_sink_closest,
                                                     SinkAccretionRadiusSquared);
#endif
                                          /*passed all merger tests*/
                                          /*going to merge massive sink with smaller sink*/
                                          if(sink_mass_i >= sink_mass_j)
                                            {
#ifdef SINK_MERGERS_DEBUG
                                              mpi_printf("SINK_MERGERS sink larger than current sink - destroying current sink \n");
#endif
                                              int keep        = i;
                                              int destroy     = j;
                                              merged[destroy] = keep;
#ifdef SINK_MERGERS_DEBUG
                                              mpi_printf("SINK_MERGERS merger array ");
                                              for(z = 0; z < NSinksAllTasks; z++)
                                                {
                                                  mpi_printf(" %d ", merged[z]);
                                                }
                                              mpi_printf("\n");
#endif
                                              /*merged array is -1 for surviving sink, or the argument of the sink it will be eaten
                                               * by*/
                                              for(z = 0; z < NSinksAllTasks; z++)
                                                /*i inherits all of the existing mergers to j  (otherwise they'll merge with nothing)*/
                                                {
                                                  if(merged[z] == j)
                                                    {
#ifdef SINK_MERGERS_DEBUG
                                                      mpi_printf("SINK_MERGERS z = %d, j = %d, merged[z] = %d \n", z, j, merged[z]);
                                                      mpi_printf("SINK_MERGERS sink %d mergers transferred to sink  %d \n", j, i);
#endif
                                                      merged[z] = i;
                                                    }
                                                }
                                            }
                                          else
                                            {
#ifdef SINK_MERGERS_DEBUG
                                              mpi_printf("SINK_MERGERS sink smaller than current sink - keeping current sink \n");
#endif
                                              int keep        = j;
                                              int destroy     = i;
                                              merged[destroy] = keep;
                                              /*merged array is -1 for surviving sink, or the argument of the sink it will be eaten
                                               * by*/
                                            }
#ifdef SINK_MERGERS_DEBUG
                                          mpi_printf("SINK_MERGERS merger array ");

                                          for(z = 0; z < NSinksAllTasks; z++)
                                            {
                                              mpi_printf(" %d ", merged[z]);
                                            }
                                          mpi_printf("\n");
#endif
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

      /*create arrays containing only the sinks involved in a merger*/
      /* 1) array with P array index for pray sinks*/
      /* 2) array with SinkP index for pray sinks*/
      /* 3) array with the HomeTasks of the pray sinks*/
      /* 4) array with P array index of the sink they are being eaten by */
      /* 5) array with SinkP index of the sink they are being eaten by*/

      MPI_Barrier(MPI_COMM_WORLD);
      int Nmerge = 0;
      for(int i = 0; i < NSinksAllTasks; i++)
        {
          if(merged[i] > -1)
            Nmerge += 1;
        }
      mpi_printf("SINK_MERGERS performing %d mergers \n", Nmerge);

      int toMerge[Nmerge];
      int toMerge_sinkNo[Nmerge];
      int mergeTasks[Nmerge];
      int eatenBy[Nmerge];
      int eatenBy_sinkNo[Nmerge];
      int merge_no = 0;

      if(Nmerge > 0)
        {
          for(z = 0; z < NSinksAllTasks; z++)
            {
              if(merged[z] > -1)
                {
                  toMerge[merge_no]        = SinkP[z].Index;
                  toMerge_sinkNo[merge_no] = z;
                  mergeTasks[merge_no]     = SinkP[z].HomeTask;
                  eatenBy[merge_no]        = SinkP[merged[z]].Index;
                  eatenBy_sinkNo[merge_no] = merged[z];
                  merge_no += 1;
                }
            }
        }

      /*going to transfer mass, momentum of pray sinks to predators & move predator to the center of mass*/
      /*go to home task of predator sink to find it in the P array, transfer properties using SinkP array*/
      /*then go to home task of the pray sink to delete it (change it to type 3)*/
      if(Nmerge > 0)
        {
          for(i = 0; i < Nmerge; i++)
            {
              if(SinkP[eatenBy_sinkNo[i]].HomeTask == ThisTask) /*predator home task*/
                {
#ifdef SINK_MERGERS_DEBUG
                  printf("SINK_MERGERS found home task of predator sink %d \n", eatenBy_sinkNo[j]);
                  printf("SINK_MERGERS check that IDs %d & %d match  \n", P[eatenBy[i]].ID, SinkP[eatenBy_sinkNo[i]].ID);
#endif
                  int index_survive_P     = eatenBy[i];
                  int index_survive_SinkP = eatenBy_sinkNo[i];
                  int index_merge_SinkP   = toMerge_sinkNo[i];

                  /*move to centre of mass */
                  P[index_survive_P].Pos[0] = (P[index_survive_P].Mass * P[index_survive_P].Pos[0] +
                                               SinkP[index_merge_SinkP].Mass * SinkP[index_merge_SinkP].Pos[0]) /
                                              (P[index_survive_P].Mass + SinkP[index_merge_SinkP].Mass);

                  P[index_survive_P].Pos[1] = (P[index_survive_P].Mass * P[index_survive_P].Pos[1] +
                                               SinkP[index_merge_SinkP].Mass * SinkP[index_merge_SinkP].Pos[1]) /
                                              (P[index_survive_P].Mass + SinkP[index_merge_SinkP].Mass);

                  P[index_survive_P].Pos[2] = (P[index_survive_P].Mass * P[index_survive_P].Pos[2] +
                                               SinkP[index_merge_SinkP].Mass * SinkP[index_merge_SinkP].Pos[2]) /
                                              (P[index_survive_P].Mass + SinkP[index_merge_SinkP].Mass);

                  /*transfer linear momentum*/
                  P[index_survive_P].Vel[0] = (P[index_survive_P].Mass * P[index_survive_P].Vel[0] +
                                               SinkP[index_merge_SinkP].Mass * SinkP[index_merge_SinkP].Vel[0]) /
                                              (P[index_survive_P].Mass + SinkP[index_merge_SinkP].Mass);

                  P[index_survive_P].Vel[1] = (P[index_survive_P].Mass * P[index_survive_P].Vel[1] +
                                               SinkP[index_merge_SinkP].Mass * SinkP[index_merge_SinkP].Vel[1]) /
                                              (P[index_survive_P].Mass + SinkP[index_merge_SinkP].Mass);

                  P[index_survive_P].Vel[2] = (P[index_survive_P].Mass * P[index_survive_P].Vel[2] +
                                               SinkP[index_merge_SinkP].Mass * SinkP[index_merge_SinkP].Vel[2]) /
                                              (P[index_survive_P].Mass + SinkP[index_merge_SinkP].Mass);

                  /*transfer mass*/
                  P[index_survive_P].Mass += SinkP[index_merge_SinkP].Mass;

#ifdef STORE_SINK_PARTICLE_SPIN
                  SinkP[index_survive_SinkP].AngularMomentum[0] += SinkP[index_merge_SinkP].AngularMomentum[0];
                  SinkP[index_survive_SinkP].AngularMomentum[1] += SinkP[index_merge_SinkP].AngularMomentum[1];
                  SinkP[index_survive_SinkP].AngularMomentum[2] += SinkP[index_merge_SinkP].AngularMomentum[2];
#endif

#ifdef SINK_MERGERS_DEBUG
                  printf("SINK_MERGERS sink %d has eaten sink %d \n", eatenBy_sinkNo[i], toMerge_sinkNo[i]);
                  printf("SINK_MERGERS old pos (%g,%g,%g) & (%g,%g,%g) \n", SinkP[index_survive_SinkP].Pos[0],
                         SinkP[index_survive_SinkP].Pos[1], SinkP[index_survive_SinkP].Pos[2], SinkP[index_merge_SinkP].Pos[0],
                         SinkP[index_merge_SinkP].Pos[1], SinkP[index_merge_SinkP].Pos[2]);
                  printf("SINK_MERGERS new pos (%g,%g,%g) \n", P[index_survive_P].Pos[0], P[index_survive_P].Pos[1],
                         P[index_survive_P].Pos[2]);
                  printf("SINK_MERGERS old vel (%g,%g,%g) & (%g,%g,%g) \n", SinkP[index_survive_SinkP].Vel[0],
                         SinkP[index_survive_SinkP].Vel[1], SinkP[index_survive_SinkP].Vel[2], SinkP[index_merge_SinkP].Vel[0],
                         SinkP[index_merge_SinkP].Vel[1], SinkP[index_merge_SinkP].Vel[2]);
                  printf("SINK_MERGERS new vels (%g,%g,%g) \n", P[index_survive_P].Vel[0], P[index_survive_P].Vel[1],
                         P[index_survive_P].Vel[0]);
                  printf("SINK_MERGERS old mass  %g & %g \n", SinkP[index_survive_SinkP].Mass, SinkP[index_merge_SinkP].Mass);
                  printf("SINK_MERGERS new mass  %g  \n", P[index_survive_P].Mass);
#endif

                  mpi_printf("SINK_MERGERS checking it made it 1 \n");
                  /*update the SinkP array for the next round*/
                }
            }
        }

      mpi_printf("SINK_MERGERS about to enter MPI Barrier 1\n");
      MPI_Barrier(MPI_COMM_WORLD);
      mpi_printf("SINK_MERGERS passed MPI Barrier 1\n");

      /*remove the eaten sinks in the P array by changing type to 3 and set mass & vel to 0 */
      /*first go to the task of the eaten sink to remove it*/
      if(Nmerge > 0)
        {
          mpi_printf("SINK_MERGERS testing 1 \n");
          for(i = 0; i < Nmerge; i++)
            {
              mpi_printf("i = %d  \n", i);
              mpi_printf("SINK_MERGERS testing 2 \n");
              //            printf("Task %d checking in, looking for Task %d  \n",ThisTask,mergeTasks[i]);
              if(mergeTasks[i] == ThisTask)
                {
#ifdef SINK_MERGERS_DEBUG
                  printf("SINK_MERGERS testing 3 \n");
                  printf("SINK_MERGERS found sink %d to keep  \n", toMerge_sinkNo[i]);
                  printf("SINK_MERGERS check that IDs %d & %d match \n", P[toMerge[i]].ID, SinkP[toMerge_sinkNo[i]].ID);
#endif

                  int index_merge_P       = toMerge[i];
                  P[index_merge_P].Type   = 3;
                  P[index_merge_P].Mass   = 0;
                  P[index_merge_P].Vel[0] = 0;
                  P[index_merge_P].Vel[1] = 0;
                  P[index_merge_P].Vel[2] = 0;
                  int index_merge_SinkP   = toMerge_sinkNo[i];
#ifdef STORE_SINK_PARTICLE_SPIN
                  SinkP[index_merge_SinkP].AngularMomentum[0] = 0;
                  SinkP[index_merge_SinkP].AngularMomentum[1] = 0;
                  SinkP[index_merge_SinkP].AngularMomentum[2] = 0;
#endif

#ifdef SINK_MERGERS_DEBUG
                  printf("delted sink mass = %g, Type = %d  \n", P[index_merge_P].Mass, P[index_merge_P].Type);
#endif
                }
            }
          mpi_printf("SINK_MERGERS about to enter MPI Barrier 2 \n");
          MPI_Barrier(MPI_COMM_WORLD);
          mpi_printf("SINK_MERGERS passed MPI Barrier 2 \n");
          int num_sinks = get_all_sink_particle_info(1);
        }
    }

#ifdef SINK_MERGERS_DEBUG
  mpi_printf("SINK_MERGERS Masses");
  for(int j = 0; j < NSinksAllTasks; j++)
    {
      mpi_printf(" %g ", SinkP[j].Mass);
    }
  mpi_printf("\n");
#endif
  return;
}
