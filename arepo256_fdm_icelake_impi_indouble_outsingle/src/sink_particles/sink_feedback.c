#include <gsl/gsl_randist.h>
#include "../allvars.h"
#include "../proto.h"

#ifdef POPIII_SNE
#define N_IMF_BINS 8
#else
#define N_IMF_BINS 10 /*we use 10 mass bins for the range between 8 and 120 Msun*/
#endif
/*! \brief
 */

void assign_sne(int sink_index, double stellar_mass, int candidate_task)
{
  int i;

  struct sne_buf
  {
    int n_sne;
    double massarr[MAXSNE];
  } sne;
  /* Divide by little h here to avoid having to keep track of it in the sink feedback routines */

  stellar_mass = stellar_mass / All.HubbleParam;

  if(ThisTask == candidate_task)
    find_number_of_sne_drawing_from_IMF(stellar_mass, &sne.n_sne, &sne.massarr[0]);

  MPI_Bcast(&sne, sizeof(struct sne_buf), MPI_BYTE, candidate_task, MPI_COMM_WORLD);

  if((SinkP[sink_index].N_sne + sne.n_sne) > MAXSNE)
    {
      for(i = 0; i < MAXACCRETIONEVENTS; i++)
        {
          if(SinkP[sink_index].MassStillToConvert[i] == 0.)
            {
              SinkP[sink_index].MassStillToConvert[i] = stellar_mass;
              SinkP[sink_index].AccretionTime[i]      = All.Time;
              mpi_printf("SINK_PARTICLES_FEEDBACK: sink %lld has N_sne = %d and just created another %d sne\n", SinkP[sink_index].ID,
                         SinkP[sink_index].N_sne, sne.n_sne);
              mpi_printf("SINK_PARTICLES_FEEDBACK: N_sne > MAXSNE, %g mass still to convert\n", stellar_mass);
              break;
            }
          if(i == (MAXACCRETIONEVENTS - 1))
            warn(
                "WARNING SINK_PARTICLES_FEEDBACK: way too many SNe for the sink particle, try to change MAXSNE and/or "
                "MAXACCRETIONEVENTS, we are throwing away n = %d SNe\n",
                sne.n_sne);
        }
    }
  else
    {
      for(i = 0; i < sne.n_sne; i++)
        {
          SinkP[sink_index].explosion_time[SinkP[sink_index].N_sne + i] = All.Time + stellar_lifetime(sne.massarr[i]);
#ifdef POPIII_SNE
          SinkP[sink_index].stellar_mass[SinkP[sink_index].N_sne + i] = sne.massarr[i];
#endif
        }

      SinkP[sink_index].N_sne += sne.n_sne;

      mpi_printf("SINK_PARTICLE_FEEDBACK: given an accreted mass of m = %g onto the sink ID = %lld, I just created %d massive stars\n",
                 stellar_mass, SinkP[sink_index].ID, sne.n_sne);
    }
}

/*! \brief given the mass of the star this function returns its life-time
 *
 *   TODO: introduce some exceptions
 */
double stellar_lifetime(double mass_of_star)
{
  /*See table 25.6 of Maeder 2009*/

#ifdef POPIII_SNE
  double m[N_IMF_BINS]    = {9., 15., 25., 40., 60., 80., 120., 200.};
  double t_ms[N_IMF_BINS] = {20220000., 10400000., 6459000., 3864000., 3464000., 3012000., 2521000., 2204000.};
#else
  double m[N_IMF_BINS]    = {7., 9., 12., 15., 20., 25., 40., 60., 85., 120.};
  double t_ms[N_IMF_BINS] = {4.319e7, 2.639e7, 1.6e7, 1.158e7, 8.141e6, 6.408e6, 4.303e6, 3.447e6, 2.823e6, 2.561e6};
#endif
  mass_of_star *= All.UnitMass_in_g / SOLAR_MASS; /*stellar mass in solar masses*/

  if(mass_of_star < m[0])
    mass_of_star = m[0];
  else if(mass_of_star > m[N_IMF_BINS - 1])
    mass_of_star = m[N_IMF_BINS - 1];

  int i = 0;
  while(mass_of_star > m[i + 1])
    i++;

  /*Linear interpolation between the values*/
  double lifetime = t_ms[i] + ((t_ms[i + 1] - t_ms[i]) / (m[i + 1] - m[i])) * (mass_of_star - m[i]);

  if(All.ComovingIntegrationOn)
    {
      /* comoving case: compute atime difference */
      double t_Gyr = get_time_difference_in_Gyr(0.0, All.Time);
      t_Gyr += lifetime / 1.0e9;
      lifetime = get_atime_of_time(t_Gyr) - All.Time;
    }
  else
    lifetime *= SEC_PER_YEAR / All.UnitTime_in_s; /*lifetime in internal units*/
  return lifetime;
}

/*  \brief given the mass of formed stars this function returns
 *   the number of massive stars.
 *
 *   Given a Kroupa IMF we know the fraction of stars for each
 *   mass bin. We then get the number of stars within each mass bin
 *   by drawing from a Poisson distribution with mean related to this
 *   fraction (see Sormani+ 2016).
 *
 */
void find_number_of_sne_drawing_from_IMF(double mass_to_convert, int *n_sne, double *massarr)
{
  *n_sne = 0;

#ifdef POPIII_SNE /* log flat IMF*/
  double m_bins[N_IMF_BINS + 1] = {8., 12., 20., 32.5, 50., 70., 100., 160., 200.};
  double m[N_IMF_BINS]          = {9., 15., 25., 40., 60., 80., 120., 200.};
  /*f is the fraction of the total mass associated to stars within the given bin, calculated given a Kroupa IMF*/

  double f[N_IMF_BINS] = {0.12596482, 0.1586969, 0.15083148, 0.13383024, 0.10453098, 0.1108073, 0.14601484, 0.06932344};
#else
  /*we use equispaced in log10 scale mass bins and the characteristic mass for
    each bin is the mean mass of stars within that bin given a Kroupa IMF*/
  double m_bins[N_IMF_BINS + 1] = {8.,          10.48815538, 13.75017542, 18.02674705, 23.63341551, 30.98386677,
                                   40.62045114, 53.25420041, 69.8172911,  91.53182469, 120.};
  double m[N_IMF_BINS]          = {9.11535456,  11.95040688, 15.66721553, 20.54002387, 26.92837024,
                          35.30361642, 46.28372683, 60.67886485, 79.55117038, 104.2931295};
  /*f is the fraction of the total mass associated to stars within the given bin, calculated given a Kroupa IMF*/

  double f[N_IMF_BINS] = {0.02914699, 0.02687268, 0.02477583, 0.0228426, 0.02106021,
                          0.0194169,  0.01790182, 0.01650496, 0.0152171, 0.01402972};
#endif

  int i, j;
  /* Converting mass arrays to code units */
  for(i = 0; i < N_IMF_BINS; i++)
    m[i] = m[i] * SOLAR_MASS / All.UnitMass_in_g;
  for(i = 0; i < N_IMF_BINS + 1; i++)
    m_bins[i] = m_bins[i] * SOLAR_MASS / All.UnitMass_in_g;

    /*for each mass bin draw from a Poisson distribution with mean mu = fi * Mtot / mi */
#ifndef SINK_FEEDBACK_SINGLE_STAR
  unsigned long int seed;
  memcpy(&seed, &mass_to_convert, sizeof(double));
  gsl_rng_set(random_generator, seed);

  double mu; /*mean of the poisson distribution*/
  int ni;    /*number of stars in given bin*/

  for(int i = 0; i < N_IMF_BINS; i++)
    {
      mu = f[i] * mass_to_convert / m[i];
      ni = gsl_ran_poisson(random_generator, mu);
      if(*n_sne + ni > MAXSNE)
        {
          int lost_sne = *n_sne + ni - MAXSNE;
          warn("Ups, too massive accretion event, no way this is going to fit in the SN vector, %d SNe will be lost\n", lost_sne);

          continue;
        }

      for(int j = *n_sne; j < (*n_sne + ni); j++)
        {
#ifdef POPIII_SNE
          massarr[j] = m_bins[i] * exp(gsl_ran_flat(random_generator, 0.0, 1.0) * log(m_bins[i + 1] / m_bins[i])); /*LOG flat IMF*/
#else
          massarr[j] = gsl_ran_flat(random_generator, m_bins[i],
                                    m_bins[i + 1]); /*the mass of the star varies uniformly within the given mass bin*/
#endif
        }
      *n_sne += ni;
      mpi_printf(
          "SINK_PARTICLE_FEEDBACK: drawing from IMF in bin %d mass in bin %g, no. of supernova %d, last entry mass array %g, mass to "
          "convert %g %g \n",
          i, mu, ni, massarr[*n_sne - 1], mass_to_convert, f[i]);
    }
#else
  massarr[0]           = SINK_FEEDBACK_SINGLE_STAR * SOLAR_MASS / All.UnitMass_in_g;
  // massarr[1] = 160.0 * SOLAR_MASS / All.UnitMass_in_g;
  *n_sne += 1;
#endif

  return;
}

double get_atime_of_time(double t_Gyr)
{
  double a_15, time, factor1;
  time = t_Gyr * (HUBBLE * All.HubbleParam);
  time *= SEC_PER_MEGAYEAR * 1000;
  factor1 = exp(3.0 * sqrt(All.OmegaLambda) / 2.0 * time);
  a_15    = (factor1 * factor1 - 1) / (2 * factor1 * sqrt(All.OmegaLambda / All.Omega0));
  return pow(a_15, 2.0 / 3.0);
}
