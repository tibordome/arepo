/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sne/sne_injection_criteria.c
 * \date        MM/YYYY
 * \author      Robin Tress, Rowan Smith, Andre Bubel
 * \brief       This file holds the functions that control where and when a new Supernova has
 *              to be injected based on the injection criterion used
 * \details     Please contact the authors at robin.tress@uni-heidelberg.de
 *              and rowan.smith@manchester.ac.uk before using this to avoid overlapping
 *              of projects. And please report any issue you may encounter by using this
 *              routine.
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"

/*! \brief based on the injection scheme used this function decides if it is
 *         time for a new Supernova
 *
 *   returns 1 if it is time for a new SN
 *
 */
enum SNE_type is_it_time_for_a_new_sn(void)
{
  enum SNE_type kaboom;
  switch(All.SNEInjectionCriterion)
    {
      case 1:
        kaboom = sne_only_once_at_center_TIMING();
        break;
      case 2:
        kaboom = sne_random_TIMING();
        break;
      case 3:
        kaboom = sne_at_cluster_position_TIMING();  // sne at the position of a single cluster particle representing a stellar cluster
        break;
      case 4:
        kaboom = sne_random_in_thin_disc_TIMING();
        break;
      case 5:
        kaboom = sne_following_a_given_general_distribution_TIMING();
        break;
      case 6:
        kaboom = sne_at_cluster_position_and_random_TIMING();
        break;
      case 7:
        kaboom = sne_sink_feedback_TIMING();
        break;
      case 8:
        kaboom = sne_sinks_and_SNIa_TIMING();
        break;
      case 9:
        kaboom = sne_sinks_and_SNIa_McMillan_TIMING();
        break;
      case 10:
        kaboom = sne_poisson_sampling_SNII_AND_SNIa_TIMING();
        break;

      default:
        mpi_printf("Wrong injection criterion chosen, no supernova will be created\n");
        kaboom = NO_SNE;
        break;
    }

  return kaboom;
}

/*! \brief based on the injection scheme used this function decides where the SN has to be positioned.
 *         It first finds a candidate position based on the injection scheme used. Then, with find_injection_cells(), finds the
 *         cells around this position where the energy has later to be deposited (the indices of these are stored into indices).
 *         If nobody complains about the position (complaints can be raised by is_this_position_viable()) the function returns,
 *         otherwise a new candidate is chosen.
 *
 *  \param sne_pos[3] variable where the position of the SN will be stored
 *  \param *radius pointer to the variable where the radius of the injection region will be stored
 *  \param *local_n will store the numer of injection particles found on the local processor
 *  \param *total_n will store the total number on any processor of injectin particles
 *  \param **indices will hold the indices of the local injection cells found
 *
 */
void determine_position_and_injection_region_of_sn(double sne_pos[3], enum SNE_type sne_type, double *radius, int *local_n,
                                                   int *total_n, int **indices, int sn_size)
{
  int viablePosition;

  while(1)
    {
      switch(sne_type)
        {
          case SNE_ONLY_ONCE:
            sne_only_once_at_center_POSITIONING(sne_pos);
            break;
          case SNE_RANDOM:
            sne_random_POSITIONING(sne_pos);
            break;
          case SNE_CLUSTER:
            sne_at_cluster_position_POSITIONING(sne_pos);
            break;
          case SNE_IN_DISC:
            sne_random_in_thin_disc_POSITIONING(sne_pos);
            break;
          case SNE_GENERAL_DISTRIBUTION:
            sne_following_a_given_general_distribution_POSITIONING(sne_pos);
            break;
          case SNE_SINK:
            sne_sink_feedback_POSITIONING(sne_pos);
            break;
          case SNE_POISSON:
            sne_poisson_sampling_SNII_AND_SNIa_POSITIONING(sne_pos);
            break;
          case SNE_DISC_PARTICLE:
            sne_at_disc_particle_position_POSITIONING(sne_pos);
            break;

          default:
            terminate("Wrong injection criterion chosen\n");
        }
      /*find injection region*/
      find_injection_cells(sne_pos, radius, local_n, total_n, indices, sn_size);
      /*is the position and the injection region fine? or should I choose a different position*/
      viablePosition = is_this_position_viable(sne_pos, *radius);
      if(viablePosition)
        {
          mpi_printf("\t found a viable position \n");
          return;
        }
      else
        myfree(*indices);
    }
}

/*! \brief based on the injeciton scheme used this function will ask if
 *         the chosen position for the SN is fine.
 *
 *  returns 1 if the position is ok, 0 otherwise.
 *
 *  \param sne_pos[3] position of the supernova
 *  \param radius radius of the injection region
 *
 */
int is_this_position_viable(double sne_pos[3], double radius)
{
  int viable_position = 1;
  if(All.SNEInjectionCriterion == 3)
    viable_position = sne_random_in_thin_disc_VALIDATE_POSITION(sne_pos, radius);
  return viable_position;
}

/*! \brief decides when it is time for a new supernova
 *         given the injection criterion All.SNEInjectionCriterion = 0 (Single SN at t = 0 at the center of the box)
 *         i.e. at t=0
 *
 *  returns 1 if it is time for a new SN
 */
enum SNE_type sne_only_once_at_center_TIMING(void)
{
  // Flag: Have we run yet?
  static int run_once = 0;
  // Only run once and only if time == 0.0
  if(run_once == 0)
    run_once = 1;
  else
    return NO_SNE;

  if(All.Time > 0.0)
    return NO_SNE;

  return SNE_ONLY_ONCE;
}

/*! \brief decides where to inject the new supernova
 *         given the injection criterion All.SNEInjectionCriterion = 0 (Single SN at t = 0 at the center of the box)
 *         i.e. at the center of the box.
 *
 */
void sne_only_once_at_center_POSITIONING(double sne_pos[3])
{
  sne_pos[0] = boxHalf_X;
  sne_pos[1] = boxHalf_Y;
  sne_pos[2] = boxHalf_Z;
}

/*! \brief decides when it is time for a new supernova
 *         given the injection criterion All.SNEInjectionCriterion = 1 (SNe randomly distributed in space, with a given rate in time).
 *         i.e. a new SN is due if we exeed a given time wich is updated based on a distribution peaked around a given rate.
 *
 */
enum SNE_type sne_random_TIMING(void)
{
  /* Not compatible with comoving coordinates */
  myassert(!All.ComovingIntegrationOn);

  double rate = 1. / (All.SNEPeriodInYears * SEC_PER_YEAR / All.UnitTime_in_s);

  /* if we are restarting from a snapshot, skip the SNe at time 0.0*/
  if(RestartFlag == RESTART_SNAPSHOT && SNERandomNextTime == 0.0)
    SNERandomNextTime = All.Time + random_exponential_distributed(rate);

  if(All.Time >= SNERandomNextTime)
    {
      SNERandomNextTime += random_exponential_distributed(rate);
      mpi_printf("SNe: Next supernova around t = %g \n", SNERandomNextTime);
      return SNE_RANDOM;
    }

  return NO_SNE;
}

/*! \brief decides where to inject the new supernova
 *         given the injection criterion All.SNEInjectionCriterion = 1 (SNe randomly distributed in space, with a given rate in time).
 *         I.e. randomly within the domain.
 *
 */
void sne_random_POSITIONING(double sne_pos[3])
{
  sne_pos[0] = boxSize_X * gsl_rng_uniform(sne_rng);
  sne_pos[1] = boxSize_Y * gsl_rng_uniform(sne_rng);
  sne_pos[2] = boxSize_Z * gsl_rng_uniform(sne_rng);
}

/*! \brief
 *
 */
enum SNE_type sne_at_cluster_position_TIMING(void)
{
#ifdef CLUSTERED_SNE
  /*supernovae are equally spaced in time*/
  double deltat = (All.SNEClusterTfin - All.SNEClusterTinit) / (All.SNENumber - 1);
  /*number of supernovae that already exploded*/
  int nt = ceil((All.Time - All.SNEClusterTinit) / deltat);

  /*if not already set, set the next time for supernova injection*/
  if(SNEClusterNextTime == 0.0)
    SNEClusterNextTime = All.SNEClusterTinit + ((nt < 0) ? 0. : nt) * deltat;

  if(All.Time >= SNEClusterNextTime && All.Time < (All.SNEClusterTfin + deltat))
    {
      SNEClusterNextTime += deltat;
      return SNE_CLUSTER;
    }

#else
  terminate(
      "If you want SN feedback from a single stellar cluster recompile the code with CLUSTERED_SNE, otherwise change injection "
      "criterion.\n");
#endif

  return NO_SNE;
}

/*! \brief
 *
 */
void sne_at_cluster_position_POSITIONING(double sne_pos[3])
{
  double sink_pos[3];

  // CAREFUL: will work only if you have a single cluster and if that is only on a single task, how it should be
  if(NumPart != NumGas)  // i.e. the cluster is on this task
    {
      sink_pos[0] = P[NumPart - 1].Pos[0];  // the sink particle index is the last in the structure
      sink_pos[1] = P[NumPart - 1].Pos[1];
      sink_pos[2] = P[NumPart - 1].Pos[2];

      for(int task = 0; task < NTask; task++)
        if(task != ThisTask)
          MPI_Send(&sink_pos, 3, MPI_DOUBLE, task, 0, MPI_COMM_WORLD);
    }
  else
    MPI_Recv(&sink_pos, 3, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  sne_pos[0] = sink_pos[0];
  sne_pos[1] = sink_pos[1];
  sne_pos[2] = sink_pos[2];
}

/*! \brief
 *
 */
enum SNE_type sne_random_in_thin_disc_TIMING(void)
{
  if(sne_random_TIMING())
    return SNE_IN_DISC;

  return NO_SNE;
}

/*! \brief
 *
 */
void sne_random_in_thin_disc_POSITIONING(double sne_pos[3])
{
  /* Not compatible with comoving coords. */
  myassert(!All.ComovingIntegrationOn);

  sne_pos[0] = boxSize_X * gsl_rng_uniform(sne_rng);
  sne_pos[1] = boxSize_Y * gsl_rng_uniform(sne_rng);

  double sne_height = log((100. * PARSEC / All.UnitLength_in_cm) / gsl_rng_uniform(sne_rng));

  double sign = 1.0;
  if(gsl_rng_uniform(sne_rng) > 0.5)
    sign = -1.0;
  sne_pos[2] = boxHalf_Z + sign * sne_height;
}

/*! \brief
 *
 */
int sne_random_in_thin_disc_VALIDATE_POSITION(double sne_pos[3], double radius)
{
  /* Not compatible with comoving coords. */
  myassert(!All.ComovingIntegrationOn);

  double rdisc_inner = 4.0 * KILOPARSEC / All.UnitLength_in_cm;
  double rdisc_outer = 10.0 * KILOPARSEC / All.UnitLength_in_cm;

  /*Add the bubble size to the supernova central position*/
  double xminus          = sne_pos[0] - radius;
  double xplus           = sne_pos[0] + radius;
  double yminus          = sne_pos[1] - radius;
  double yplus           = sne_pos[1] + radius;
  double zminus          = sne_pos[2] - radius;
  double zplus           = sne_pos[2] + radius;
  double galactic_radius = sqrt(pow(sne_pos[0] - boxHalf_X, 2) + pow(sne_pos[1] - boxHalf_Y, 2));

  /* Ensure SNe only go off within central 500 pc in the z-direction for computational reasons- this should be rare anyway due to the
   * exponential function.*/

  double snzmin = boxHalf_Z - 250. * PARSEC / All.UnitLength_in_cm;
  double snzmax = boxHalf_Z + 250. * PARSEC / All.UnitLength_in_cm;

  if(xminus < 0.0 || yminus < 0.0 || zminus < snzmin || xplus > boxSize_X || yplus > boxSize_Y || zplus > snzmax ||
     galactic_radius < rdisc_inner * 1.05 || galactic_radius > rdisc_outer * 0.95)
    return 0; /*Rejected*/
  else
    return 1; /*Accepted*/
}

/*! \brief
 *
 */
enum SNE_type sne_following_a_given_general_distribution_TIMING(void)
{
  if(sne_random_TIMING())
    return SNE_GENERAL_DISTRIBUTION;

  return NO_SNE;
}

/*! \brief positions of the supernovae follow a general distribution.
 *         The CDF of the radial and z distribution are defined elsewhere
 *         (see inverse_radial_CDF() and inverse_z_CDF() in sne_utility.c)
 */
void sne_following_a_given_general_distribution_POSITIONING(double sne_pos[3])
{
  // is this the first call of this function? then initialize the inverse_radial_CDF
  static int first_call = 1;
  if(first_call)
    {
      first_call = 0;
      sne_initialize_radial_CDF();
    }

  double rr    = inverse_radial_CDF(gsl_rng_uniform(sne_rng));
  double theta = gsl_rng_uniform(sne_rng) * 2. * M_PI;
  double z     = inverse_z_CDF(gsl_rng_uniform(sne_rng));

  sne_pos[0] = rr * cos(theta) + boxHalf_X;
  sne_pos[1] = rr * sin(theta) + boxHalf_Y;
  sne_pos[2] = z + boxHalf_Z;
}

enum SNE_type sne_at_cluster_position_and_random_TIMING(void)
{
  if(sne_at_cluster_position_TIMING())
    return SNE_CLUSTER;

  if(sne_random_TIMING())
    return SNE_GENERAL_DISTRIBUTION;

  return NO_SNE;
}

int check_for_additional_sne(int sinkidx, int idx, int n_thistimeAlreadyDone)
{
#if defined(SINK_PARTICLES) && defined(SINK_PARTICLES_FEEDBACK)
  mpi_printf("Checking sink %lld for additional SNe\n", SinkP[sinkidx].ID);

  int n_sne, n_sneLeft = 0;
  double massarr[MAXSNE];
  double explosion_time[MAXSNE];
  int idx_sneLeft[MAXSNE];
  int foundAnotherSN = 0;

  find_number_of_sne_drawing_from_IMF(SinkP[sinkidx].MassStillToConvert[idx], &n_sne, &massarr[0]);
  // count how many are left
  for(int i = 0; i < n_sne; i++)
    {
      explosion_time[i] = SinkP[sinkidx].AccretionTime[idx] + stellar_lifetime(massarr[i]);
      if(explosion_time[i] > All.LastSNTime)
        {
          idx_sneLeft[n_sneLeft] = i;
          n_sneLeft += 1;
        }
    }
  mpi_printf("Mass to convert = %g, n_sne = %d, n_sneLeft = %d\n", SinkP[sinkidx].MassStillToConvert[idx], n_sne, n_sneLeft);

  // if the remaining SNe fit into the standard SN vector copy them there
  if(SinkP[sinkidx].N_sne + n_sneLeft <= MAXSNE)
    {
      mpi_printf("Sink ID = %lld, N_sne = %d, n_sneLeft = %d\n", SinkP[sinkidx].ID, SinkP[sinkidx].N_sne, n_sneLeft);
      mpi_printf("The additional SNe fit into the standard SN vector now, we copy them there\n");
      for(int i = 0; i < n_sneLeft; i++)
        SinkP[sinkidx].explosion_time[SinkP[sinkidx].N_sne + i] = explosion_time[idx_sneLeft[i]];

      SinkP[sinkidx].N_sne += n_sneLeft;
      SinkP[sinkidx].MassStillToConvert[idx] = 0;
      return 0;  // check
    }

  // check if this buffer has some additional SNe
  int n_thistime = 0;

  for(int i = 0; i < n_sneLeft; i++)
    if(All.Time > explosion_time[idx_sneLeft[i]])
      n_thistime += 1;

  if((n_thistime - n_thistimeAlreadyDone == 1) && (n_sneLeft - n_thistimeAlreadyDone == 1)) /*the last SN in this chunk*/
    {
      mpi_printf("last additional SN for this accretion event\n");
      SinkP[sinkidx].MassStillToConvert[idx] = 0;
    }

  return n_thistime;
#else
  return 0;
#endif
}

enum SNE_type sne_sink_feedback_TIMING(void)
{
#if defined(SINK_PARTICLES) && defined(SINK_PARTICLES_FEEDBACK)
  /*Loop over all the sinks*/
  static int n_thistimeAlreadyDone = 0;

  for(int i = Sink_about_to_explode; i < NSinksAllTasks; i++)
    {
      if(SinkP[i].ID == 0 && SinkP[i].Mass == 0)
        continue;

      /*check within the massive stars of the sink wether there is one that exceeded its life time*/
      for(int j = ((i == Sink_about_to_explode) ? Star_about_to_explode : 0); j < SinkP[i].N_sne; j++)
        {
          /*if the time has come for the star it goes supernova */
          if(All.Time > SinkP[i].explosion_time[j])
            {
              Sink_about_to_explode = i;
              Star_about_to_explode = j;
              return SNE_SINK;
            }
        }

      /*check if there are additional sne in the buffer vector*/
      Wherewasi = (Wherewasi == -1) ? 0 : Wherewasi;
      for(int j = Wherewasi; j < MAXACCRETIONEVENTS; j++)
        {
          if(SinkP[i].MassStillToConvert[j] != 0)
            {
              int foundAnotherSN = check_for_additional_sne(i, j, n_thistimeAlreadyDone);
              if(foundAnotherSN)
                {
                  mpi_printf("Found %d additional SNe this turn, %d of them already exploded\n", foundAnotherSN,
                             n_thistimeAlreadyDone);

                  n_thistimeAlreadyDone += 1;

                  if(foundAnotherSN == n_thistimeAlreadyDone)
                    {
                      mpi_printf("this is the last one for this chunk.\n");
                      Wherewasi             = j + 1;
                      n_thistimeAlreadyDone = 0;
                    }

                  Sink_about_to_explode = i;

                  return SNE_SINK;
                }
            }
        }
      Wherewasi = -1;
    }
  /*if no SNe were found reset flags*/
  Sink_about_to_explode = 0;
  Star_about_to_explode = 0;

#else
  terminate(
      "If you want SN feedback coupled with sink particles recompile the code with SINK_PARTICLES and SINK_PARTICLES_FEEDBACK, "
      "otherwise change injection criterion.");
#endif
  return NO_SNE;
}

void sne_sink_feedback_POSITIONING(double sne_pos[3])
{
#if defined(SINK_PARTICLES) && defined(SINK_PARTICLES_FEEDBACK)
  /*find position of the SN associated to the flagged sink*/
  sne_pos[0] = SinkP[Sink_about_to_explode].Pos[0] + gsl_ran_gaussian(sne_rng, All.SNEScatterAroundSink);
  sne_pos[1] = SinkP[Sink_about_to_explode].Pos[1] + gsl_ran_gaussian(sne_rng, All.SNEScatterAroundSink);
  sne_pos[2] = SinkP[Sink_about_to_explode].Pos[2] + gsl_ran_gaussian(sne_rng, All.SNEScatterAroundSink);

  if(Wherewasi < 0)
    {
      /*Safety check*/
      if(All.Time < SinkP[Sink_about_to_explode].explosion_time[Star_about_to_explode])
        terminate("Something went wrong, the wrong sink/star is flagged as about to explode");

      /*remove the dead star from the star list associated to the sink*/
      SinkP[Sink_about_to_explode].explosion_time[Star_about_to_explode] =
          SinkP[Sink_about_to_explode].explosion_time[SinkP[Sink_about_to_explode].N_sne - 1];
#ifdef POPIII_SNE
      SinkP[Sink_about_to_explode].stellar_mass[Star_about_to_explode] =
          SinkP[Sink_about_to_explode].stellar_mass[SinkP[Sink_about_to_explode].N_sne - 1];
#endif

      SinkP[Sink_about_to_explode].N_sne -= 1;

      mpi_printf("SINK_PARTICLES_FEEDBACK EXPLOSION: we just exploded sink %lld which now has still %d sne to go\n",
                 SinkP[Sink_about_to_explode].ID, SinkP[Sink_about_to_explode].N_sne);
    }
  else
    mpi_printf("SINK_PARTICLES_FEEDBACK EXPLOSION: we just exploded an additional SN of sink %lld\n", SinkP[Sink_about_to_explode].ID);

#endif
}

void sne_at_disc_particle_position_POSITIONING(double sne_pos[3])
{
  int nDiscStarsThisTask = 0;
  int selectedTask;

  // find total number of disc particles on each processor
  for(int i = NumGas; i < NumPart; i++)
    {
      if(P[i].Type == 2)
        nDiscStarsThisTask++;
    }

  int *nDiscStarsEachTask = (int *)mymalloc("nDiscStarsEachTask", NTask * sizeof(int));
  MPI_Allgather(&nDiscStarsThisTask, 1, MPI_INT, nDiscStarsEachTask, 1, MPI_INT, MPI_COMM_WORLD);

  int nTotalDiscStars = header.npart[2];

  // get random number between 0 and total number of disc particles
  int selectedStarIndex = floor(gsl_rng_uniform(sne_rng) * nTotalDiscStars);

  // find on which processor the selected particle lies, if it is this one then get the position of this particle
  int cumSumStars = 0;
  for(int i = 0; i < NTask; i++)
    {
      if(selectedStarIndex < (cumSumStars + nDiscStarsEachTask[i]))
        {
          selectedTask = i;
          break;
        }
      cumSumStars += nDiscStarsEachTask[i];
    }

  if(ThisTask == selectedTask)
    {
      int localSelectedStarIndex = selectedStarIndex - cumSumStars;

      int nDiscStarsSelectedTask = 0;
      for(int i = NumGas; i < NumPart; i++)
        {
          if(P[i].Type == 2)
            nDiscStarsSelectedTask++;

          if(localSelectedStarIndex == nDiscStarsSelectedTask)
            {
              sne_pos[0] = P[i].Pos[0];
              sne_pos[1] = P[i].Pos[1];
              sne_pos[2] = P[i].Pos[2];
            }
        }
    }

  // comunicate the position just found to the other processors
  MPI_Bcast(sne_pos, 3, MPI_DOUBLE, selectedTask, MPI_COMM_WORLD);

  myfree(nDiscStarsEachTask);

  return;
}

enum SNE_type sne_sinks_and_SNIa_TIMING(void)
{
  if(sne_random_TIMING())
    return SNE_DISC_PARTICLE;

  if(sne_sink_feedback_TIMING())
    return SNE_SINK;

  return NO_SNE;
}

enum SNE_type sne_poisson_sampling_SNII_AND_SNIa_TIMING(void)
{
#ifdef GALPOT
  if(NSNetot > 0)
    return SNE_POISSON;
  if(NSNetot == 0)
    {
      NSNetot = -1;
      return NO_SNE;
    }

  double rateSNII = .8 / (All.SNEPeriodInYears * SEC_PER_YEAR / All.UnitTime_in_s);
  double rateSNIa = .2 / (All.SNEPeriodInYears * SEC_PER_YEAR / All.UnitTime_in_s);

  int nSNeEstimate = 10 * ceil(rateSNIa * All.TimeStep + rateSNII * All.TimeStep);
  struct snPositions *localSNePos =
      (struct snPositions *)mymalloc_movable(&localSNePos, "localSNePos", nSNeEstimate * sizeof(struct snPositions));

  int nSNeloc = 0;
  for(int i = 0; i < NumGas; i++)
    {
      double muSNII = rateSNII * All.TimeStep * P[i].Mass / GlobalGasMass;

      double xi = P[i].Pos[0] - boxHalf_X;
      double yi = P[i].Pos[1] - boxHalf_Y;
      double zi = P[i].Pos[2] - boxHalf_Z;

      double stellarMass = galpot_find_stellar_mass_of_particle(xi, yi, zi, All.Time, SphP[i].Volume);
      double muSNIa      = rateSNIa * All.TimeStep * stellarMass / GlobalStellarMass;

      int nSNII = gsl_ran_poisson(random_generator, muSNII);
      int nSNIa = gsl_ran_poisson(random_generator, muSNIa);

      int nSNeHere = nSNII + nSNIa;

      for(int j = 0; j < nSNeHere; j++)
        {
          double hi                       = get_cell_radius(i) * 0.00001;
          localSNePos[nSNeloc + j].pos[0] = P[i].Pos[0] + hi;
          localSNePos[nSNeloc + j].pos[1] = P[i].Pos[1] + hi;
          localSNePos[nSNeloc + j].pos[2] = P[i].Pos[2] + hi;
        }

      nSNeloc += nSNeHere;
    }

  int *nSNeEachTask;

  nSNeEachTask = mymalloc_movable(&nSNeEachTask, "nSNeEachTask", NTask * sizeof(int));
  MPI_Allgather(&nSNeloc, 1, MPI_INT, nSNeEachTask, 1, MPI_INT, MPI_COMM_WORLD);

  NSNetot = 0;
  for(int i = 0; i < NTask; i++)
    NSNetot += nSNeEachTask[i];

  mpi_printf("FOUND %d total SNe this timestep\n", NSNetot);
  int *nSNeOffset = mymalloc_movable(&nSNeOffset, "nSNeOffset", NTask * sizeof(int));
  nSNeOffset[0]   = 0;
  for(int i = 1; i < NTask; i++)
    nSNeOffset[i] = nSNeOffset[i - 1] + nSNeEachTask[i - 1];

  SNePositions = mymalloc_movable(&SNePositions, "SNePositions", NSNetot * sizeof(struct snPositions));

  for(int itask = 0; itask < NTask; itask++)
    {
      if(nSNeEachTask[itask] > 0)
        {
          struct snPositions *export_snPos =
              (struct snPositions *)mymalloc("export_SinkP", nSNeEachTask[itask] * sizeof(struct snPositions));

          if(itask == ThisTask)
            {
              for(int i = 0; i < nSNeloc; i++)
                {
                  export_snPos[i].pos[0] = localSNePos[i].pos[0];
                  export_snPos[i].pos[1] = localSNePos[i].pos[1];
                  export_snPos[i].pos[2] = localSNePos[i].pos[2];
                }
            }

          MPI_Bcast(export_snPos, nSNeEachTask[itask] * sizeof(struct snPositions), MPI_BYTE, itask, MPI_COMM_WORLD);

          for(int i = 0; i < nSNeEachTask[itask]; i++)
            {
              int ioffset                  = nSNeOffset[itask] + i;
              SNePositions[ioffset].pos[0] = export_snPos[i].pos[0];
              SNePositions[ioffset].pos[1] = export_snPos[i].pos[1];
              SNePositions[ioffset].pos[2] = export_snPos[i].pos[2];
            }
          myfree(export_snPos);
        }
    }

  myfree_movable(nSNeOffset);
  myfree_movable(nSNeEachTask);
  myfree_movable(localSNePos);

  if(NSNetot > 0)
    return SNE_POISSON;

  NSNetot = -1;
  myfree_movable(SNePositions);
#else
  terminate("This SN injection criterium requires GALPOT to be active");
#endif
  return NO_SNE;
}

void sne_poisson_sampling_SNII_AND_SNIa_POSITIONING(double sne_pos[3])
{
#ifdef GALPOT
  sne_pos[0] = SNePositions[NSNetot - 1].pos[0];
  sne_pos[1] = SNePositions[NSNetot - 1].pos[1];
  sne_pos[2] = SNePositions[NSNetot - 1].pos[2];

  NSNetot--;
  if(NSNetot == 0)
    myfree_movable(SNePositions);
#else
  terminate("This SN injection criterium requires GALPOT to be active\n");
#endif
}

enum SNE_type sne_sinks_and_SNIa_McMillan_TIMING(void)
{
  if(sne_random_TIMING())
    return SNE_GENERAL_DISTRIBUTION;

  if(sne_sink_feedback_TIMING())
    return SNE_SINK;

  return NO_SNE;
}
