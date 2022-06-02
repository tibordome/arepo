#include "../allvars.h"
#include "../proto.h"

/*! \brief Master function to perform simplified UV ionization feedback around the sink particles
 *
 *         Based on the number of massive stars within the sink particles, we get an ionising
 *         photon flux. We then look for which cells surrounding the sink particle this photon
 *         flux is able to keep ionised (for details on how this is done see function find_stromgren_radius_around_sink()).
 *         These cells are then heated up to All.StromgrenTemperature and ionised.
 *         Also these cells are prevented from cooling during the next timestep as the photon flux is supposed to be
 *         continuous and thus keeping the surrounding gas hot and ionised.
 *
 */
void photoionisation_feedback_from_sinks()
{
  mpi_printf("SINK_FEEDBACK: entering photo ionization routine\n");

  double t0 = second();

  int local_n, total_n, *indices;

  for(int isink = 0; isink < NSinksAllTasks; isink++)
    {
      if(SinkP[isink].N_sne > 0)
        {
          double strom_radius = find_stromgren_radius_around_sink(isink, &local_n, &total_n, &indices);

          heat_and_ionise_stromgren_cells(local_n, indices);

          myfree(indices);
        }
    }

  double t1 = second();

  mpi_printf("SINK_FEEDBACK: photo ionisation around sinks done, took %g sec.\n", timediff(t0, t1));
}

/*! \brief finds the region around the sink that the massive stars within the
 *         sink are able to ionise.
 *
 *         We find the cells within R_max to the sink. We then walk these cells
 *         in increasing radius order and flag it if the photon budget is still
 *         enough to fully ionise the cell.
 *
 *  \param isink index of the sink in the SinkP structure
 *  \param local_n number of particles on the local process to be ionised
 *  \param total_n number of total particles on any processor to be ionised
 *  \param **indices indices of the SphP structure of the local particles
 *                   to be ionised.
 *
 *  returns the maximum distance the photons manage to travel before they are all
 *  consumed by ionisations.
 */
double find_stromgren_radius_around_sink(int isink, int *local_n, int *total_n, int **indices)
{
  double RsMax = 50. * PARSEC / All.UnitLength_in_cm;

  double beta = 2.56e-13;

  /*collect all particles within RsMax from the sink*/
  find_particles_within_a_sphere(SinkP[isink].Pos, RsMax, local_n, total_n, indices);

  /*gather these cells on all processors on the structure storeGlobalData*/
  struct collect_data
  {
    double r2_ParticleSink;
    double localRecombinations;

    int homeTask;
    int localIndex;

  } * storeLocalData, *storeGlobalData;

  storeGlobalData = (struct collect_data *)mymalloc("storeGlobalData", *total_n * sizeof(struct collect_data));
  storeLocalData  = (struct collect_data *)mymalloc("storeLocalData", *local_n * sizeof(struct collect_data));

  int *nfoundEachTask = (int *)mymalloc("nfoundEachTask", NTask * sizeof(int));
  MPI_Allgather(local_n, 1, MPI_INT, nfoundEachTask, 1, MPI_INT, MPI_COMM_WORLD);

  double sinkx = SinkP[isink].Pos[0];
  double sinky = SinkP[isink].Pos[1];
  double sinkz = SinkP[isink].Pos[2];

  for(int i = 0; i < *local_n; i++)
    {
      int idx = (*indices)[i];

      storeLocalData[i].r2_ParticleSink =
          pow((P[idx].Pos[0] - sinkx), 2) + pow((P[idx].Pos[1] - sinky), 2) + pow((P[idx].Pos[2] - sinkz), 2);

      double yn = (SphP[idx].Density * All.UnitDensity_in_cgs) / ((1.0 + 4.0 * ABHE) * PROTONMASS);
      double Yn = P[idx].Mass * All.UnitMass_in_g / ((1.0 + 4.0 * ABHE) * PROTONMASS);

      double n_HHe = (1. + ABHE) * yn;
      double N_HHe = (1. + ABHE) * Yn;

      storeLocalData[i].localRecombinations = beta * n_HHe * N_HHe;

      storeLocalData[i].homeTask = ThisTask;

      storeLocalData[i].localIndex = idx;
    }

  int *bytecounts = (int *)mymalloc("bytecounts", sizeof(int) * NTask);
  int *byteoffset = (int *)mymalloc("byteoffset", sizeof(int) * NTask);

  int itask;
  for(itask = 0; itask < NTask; itask++)
    bytecounts[itask] = nfoundEachTask[itask] * sizeof(struct collect_data);

  for(itask = 1, byteoffset[0] = 0; itask < NTask; itask++)
    byteoffset[itask] = byteoffset[itask - 1] + bytecounts[itask - 1];

  MPI_Allgatherv(&storeLocalData[0], *local_n * sizeof(struct collect_data), MPI_BYTE, &storeGlobalData[0], bytecounts, byteoffset,
                 MPI_BYTE, MPI_COMM_WORLD);

  myfree(byteoffset);
  myfree(bytecounts);
  myfree(nfoundEachTask);
  myfree(storeLocalData);

  /*Walk the particle list until photons are exhausted*/
  *local_n  = 0;
  int tot_n = 0;

  double Rs                  = 0.;
  double totalRecombinations = 0.;

  double ionPhotonRate_singleStar = 1.e49;
  double ionPhotonRate            = SinkP[isink].N_sne * ionPhotonRate_singleStar;

  /*sort storeGlobalData with increasing distance to the sink*/
  for(int i = 0; i < *total_n - 1; i++)
    {
      for(int j = i + 1; j < *total_n; j++)
        {
          if(storeGlobalData[j].r2_ParticleSink < storeGlobalData[i].r2_ParticleSink)
            {
              struct collect_data tmp = storeGlobalData[i];
              storeGlobalData[i]      = storeGlobalData[j];
              storeGlobalData[j]      = tmp;
            }
        }

      /*walk until photons are exhausted*/
      totalRecombinations += storeGlobalData[i].localRecombinations;

      if(totalRecombinations > ionPhotonRate)
        break;

      Rs = sqrt(storeGlobalData[i].r2_ParticleSink);
      if(ThisTask == storeGlobalData[i].homeTask)
        (*indices)[(*local_n)++] = storeGlobalData[i].localIndex;

      tot_n++;
    }

  *total_n = tot_n;

  myfree(storeGlobalData);
  return Rs;
}

/*! \brief heats the selected particles to a temperature of All.StromgrenTemperature,
 *         ionises the gas and flaggs the cell to be skipped during the next sgchem call
 *
 *  \param local_n numer of particles on the local processor flagged to be ionised
 *  \param indices[] indices of these particles
 *
 */
void heat_and_ionise_stromgren_cells(int local_n, int indices[])
{
  for(int i = 0; i < local_n; i++)
    {
      int idx = indices[i];

      /*ionise everything but He and metals*/
#if CHEMISTRYNETWORK == 1
      SphP[idx].TracAbund[IH2]    = 0.;
      SphP[idx].TracAbund[IHP]    = 1.;
      SphP[idx].TracAbund[IDP]    = All.DeutAbund;
      SphP[idx].TracAbund[IHD]    = 0.;
      SphP[idx].TracAbund[IHATOM] = 0.;
      SphP[idx].TracAbund[IDATOM] = 0.;
#endif
#if CHEMISTRYNETWORK == 4
      SphP[idx].TracAbund[IH2]    = 0.;
      SphP[idx].TracAbund[IHP]    = 1.;
      SphP[idx].TracAbund[IHATOM] = 0.;
#endif
#if CHEMISTRYNETWORK == 5
      SphP[idx].TracAbund[IH2]    = 0.;
      SphP[idx].TracAbund[IHP]    = 1.;
      SphP[idx].TracAbund[IHATOM] = 0.;
#ifdef SGCHEM_VARIABLE_Z
      SphP[idx].TracAbund[ICP] = SphP[i].CarbAbund;
      SphP[idx].TracAbund[ICO] = 0.;
#else
      SphP[idx].TracAbund[ICP]    = All.CarbAbund;
      SphP[idx].TracAbund[ICO]    = 0.;
#endif
#endif /* CHEMISTRYNETWORK == 5 */
#if CHEMISTRYNETWORK == 15
      SphP[idx].TracAbund[IH2]    = 0.;
      SphP[idx].TracAbund[IHP]    = 1.;
      SphP[idx].TracAbund[ICHX]   = 0.;
      SphP[idx].TracAbund[IOHX]   = 0.;
      SphP[idx].TracAbund[ICO]    = 0.;
      SphP[idx].TracAbund[IHCOP]  = 0.;
      SphP[idx].TracAbund[IHATOM] = 0.;
#ifdef SGCHEM_VARIABLE_Z
      SphP[idx].TracAbund[ICP]    = SphP[idx].CarbAbund;
      SphP[idx].TracAbund[ICATOM] = 0.;
      SphP[idx].TracAbund[IOATOM] = SphP[idx].OxyAbund;
#else
      SphP[idx].TracAbund[ICP]    = All.CarbAbund;
      SphP[idx].TracAbund[ICATOM] = 0.;
      SphP[idx].TracAbund[IOATOM] = All.OxyAbund;
#endif
#endif /* CHEMISTRYNETWORK == 15 */
      for(int j = 0; j < SGCHEM_NUM_SPECIES; j++)
        SphP[idx].MassTracAbund[j] = P[idx].Mass * SphP[idx].TracAbund[j];

      /*Heat to T = All.StromgrenTemperature K*/
      double gas_temp = find_temperature_of_cell(idx);

      if(gas_temp < All.StromgrenTemperature)
        {
          double yn = (SphP[idx].Density * All.UnitDensity_in_cgs) / ((1.0 + 4.0 * ABHE) * PROTONMASS);

#if CHEMISTRYNETWORK == 1
          double yntot = (1.0 + ABHE - SphP[idx].TracAbund[IH2] + SphP[idx].TracAbund[IHP] + SphP[idx].TracAbund[IHEP] +
                          2. * SphP[idx].TracAbund[IHEPP]) *
                         yn;
#elif CHEMISTRYNETWORK == 15
          double yntot = (1.0 + ABHE - SphP[idx].TracAbund[IH2] + SphP[idx].TracAbund[IHP] + SphP[idx].TracAbund[IHEP]) * yn;
#else
          double yntot = (1.0 + ABHE - SphP[idx].TracAbund[IH2] + SphP[idx].TracAbund[IHP]) * yn;
#endif
          double gas_energy_density = 1.5 * yntot * BOLTZMANN * All.StromgrenTemperature;
          double unew               = gas_energy_density * pow(All.UnitLength_in_cm, 3) / (SphP[idx].Density * All.UnitEnergy_in_cgs);
          SphP[idx].Utherm          = unew;
#ifdef MHD_THERMAL_ENERGY_SWITCH
          SphP[idx].Etherm = SphP[idx].Utherm * P[idx].Mass;
#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
          // Update Entropy, Taken from SGChem
          SphP[idx].A       = GAMMA_MINUS1 * SphP[idx].Utherm / pow(SphP[idx].Density * All.cf_a3inv, GAMMA_MINUS1);
          SphP[idx].Entropy = log(SphP[idx].A) * P[idx].Mass;
#endif
          double mom2 = SphP[idx].Momentum[0] * SphP[idx].Momentum[0] + SphP[idx].Momentum[1] * SphP[idx].Momentum[1] +
                        SphP[idx].Momentum[2] * SphP[idx].Momentum[2];

          SphP[idx].Energy = SphP[idx].Utherm * P[idx].Mass + 0.5 * mom2 / P[idx].Mass;
#ifdef MHD
          SphP[idx].Energy +=
              0.5 * (SphP[idx].B[0] * SphP[idx].B[0] + SphP[idx].B[1] * SphP[idx].B[1] + SphP[idx].B[2] * SphP[idx].B[2]) * volume;
#endif
          /*flag them as cells that shall not to be cooled*/
          SphP[idx].CoolingFlag = 0;
        }
    }
}
