#include "../allvars.h"
#include "../proto.h"
#include <gsl/gsl_randist.h>

#define N_IMF_BINS 10 /*we use 10 mass bins for the range between 8 and 120 Msun*/

/*! \brief 
 */ 
void assign_sne(int sink_ID, double mass_to_convert, int candidate_task)
{  
  int i;

  struct sne_buf
    {
      int n_sne;
      double massarr[MAXSNE];
    } sne;
  
  double effective_stellar_mass = All.SINKStarFormationEfficiency * mass_to_convert;

  if(ThisTask == candidate_task)
    find_number_of_sne_drawing_from_IMF(effective_stellar_mass, &sne.n_sne, &sne.massarr[0]);

  MPI_Bcast(&sne, sizeof(struct sne_buf), MPI_BYTE, candidate_task, MPI_COMM_WORLD);

  if ((SinkP[sink_ID].N_sne + sne.n_sne) > MAXSNE)
    terminate("SINK_PARTICLES_FEEDBACK ERROR: N_sne > MAXSNE, increase MAXSNE\n");

  for(i=0; i<sne.n_sne; i++)
      SinkP[sink_ID].explosion_time[SinkP[sink_ID].N_sne + i] = All.Time + stellar_lifetime(sne.massarr[i]);

  SinkP[sink_ID].N_sne += sne.n_sne;

  mpi_printf("SINK_PARTICLE_FEEDBACK: given an accreted mass of m = %g onto the sink ID = %g, I just created %d massive stars\n", mass_to_convert, sink_ID, sne.n_sne);
}

/*! \brief given the mass of the star this function returns its life-time
 *   
 *   TODO: introduce some exceptions
 */ 
double stellar_lifetime(double mass_of_star)
{
  /*See table 25.6 of Maeder 2009*/
  double m[N_IMF_BINS]    = {7.,      9.,      12.,   15.,     20.,     25.,     40.,     60.,     85.,     120.};
  double t_ms[N_IMF_BINS] = {4.319e7, 2.639e7, 1.6e7, 1.158e7, 8.141e6, 6.408e6, 4.303e6, 3.447e6, 2.823e6, 2.561e6};


  mass_of_star *= All.UnitMass_in_g / SOLAR_MASS; /*stellar mass in solar masses*/

  int i=0;
  while (mass_of_star > m[i+1])
    i++; 

  /*Linear interpolation between the values*/
  double lifetime = t_ms[i] + ((t_ms[i+1]-t_ms[i]) / (m[i+1]-m[i])) * (mass_of_star - m[i]);
  
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

 /*we use equispaced in log10 scale mass bins and the characteristic mass for
   each bin is the mean mass of stars within that bin given a Kroupa IMF*/

  double m_bins[N_IMF_BINS+1] = {8., 10.48815538, 13.75017542, 18.02674705, 23.63341551, 30.98386677, 40.62045114, 53.25420041, 69.8172911 , 91.53182469, 120.}; 

  double m[N_IMF_BINS] = {9.11535456, 11.95040688, 15.66721553, 20.54002387, 26.92837024, 35.30361642, 46.28372683, 60.67886485, 79.55117038, 104.2931295};

 /*f is the fraction of the total mass associated to stars within the given bin, calculated given a Kroupa IMF*/

  double f[N_IMF_BINS] = {0.02914699, 0.02687268, 0.02477583, 0.0228426, 0.02106021, 0.0194169, 0.01790182, 0.01650496, 0.0152171, 0.01402972};


 /*for each mass bin draw from a Poisson distribution with mean mu = fi * Mtot / mi */
  srand(time(NULL));
  unsigned long randSeed = rand();
  gsl_rng_set(random_generator, randSeed);

  int i,j;
  double mu; /*mean of the poisson distribution*/
  int ni; /*number of stars in given bin*/

  for(i=0; i<N_IMF_BINS; i++)
    {
      mu = f[i] * mass_to_convert/m[i];
      ni = gsl_ran_poisson(random_generator, mu);
      for (j=*n_sne; j<(*n_sne+ni); j++)
        massarr[j] = gsl_ran_flat(random_generator, m_bins[i], m_bins[i+1]); /*the mass of the star varies uniformly within the given mass bin*/
      *n_sne += ni;
    }

  return;
}
