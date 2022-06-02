/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \date        22/10/2018
 * \author      Ondrej Jaura
 * \brief       Main photon transport functionality
 * \details     
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"
#include "sx_direction_bins.h"

//======================================================================================================//
//                                     Other common functions                                           //
//======================================================================================================//

/*!
 * \biref Initialize the site structure
 */
_Bool sx_phot_initialize_site(struct sxSite *site, int index)
{
  _Bool hasPhotons = 0;
  int sid = sxCell[index].Source;
  int spdOffset = index*SX_NDIR; 
  
  // Initialize the site
  site->index = index;
  site->isSource = (sid<0) ? 0: 1;

  // DEBUG: this timing scheme should be simplified somehow
#if(SX_RUNNING_MODE == 1) || (SX_RUNNING_MODE == 2)
  site->dt = sxRunData->timeStep;
#else
  site->dt = (P[index].TimeBinHydro ? (((integertime)1) << P[index].TimeBinHydro) : 0) * All.Timebase_interval * All.UnitTime_in_s;
#endif
  
  for(int f=0; f<SX_NFREQ; f++)
    {
      site->nPhotons[f]       = 0.;
      site->nPhotonsSource[f] = 0.;

      // if source, add photons from ionisation rate
#ifdef SX_SINGLE_SHOT
      if (site->isSource && sxRunData->currentStep==0) // all at the first step only
	{	  
	  double nPhotons = sxSources[sid].ionRate[f] * sxRunData->runTime;
#else
      if (site->isSource)                              // continualy 
	{	  
	  double nPhotons = sxSources[sid].ionRate[f] * sxRunData->timeStep;
#endif
#if(SX_NUM_ROT > 1)
          nPhotons /= SX_NUM_ROT;
#endif
	  site->nPhotons[f] += nPhotons;
	  site->nPhotonsSource[f] += nPhotons;
	}
#ifdef SX_SINGLE_SHOT
      else
        {
#endif
          // add current photons to the photon bundle
          for(int d=0; d < SX_NDIR; d++ )
	    {
	      int app = sxSPDtoAPP[ spdOffset + d ];
	      if (app<0) continue;
	      site->nPhotons[f] += (double)sxAPP[app].nPhotons[f];
	    }
#ifdef SX_SINGLE_SHOT
        }
#endif

#ifdef SX_OUTPUT_FLUX
// Add up the flux values in code units
#ifdef SX_SKIP_RADIUS
      if (sxCell[site->index].skip==0)
	SphP[index].sxPhotonFlux[f] += site->nPhotons[f] * sxRunData->runTimeInv;
#else
      SphP[index].sxPhotonFlux[f] += site->nPhotons[f] * sxRunData->runTimeInv;
#endif
#endif

      // convert to physical units
      site->nPhotons[f] *= All.UnitPhotons_per_s;
      site->nPhotonsSource[f] *= All.UnitPhotons_per_s;

      // save the original number of photons
      site->nPhotonsOriginal[f] = site->nPhotons[f];

#ifdef SX_DISPLAY_STATS
      // Update statistics
      sxRunData->nEmitPhotons[f] += site->nPhotonsSource[f];
#endif

      // check if there are some photons waiting for trasfer
      if(site->nPhotons[f] > 0.)
        hasPhotons = 1;
    }

    // initialize chemistry specific variables
#if(SX_CHEMISTRY == 4)
  // For SX_CHEMISTRY==4 we need to treat the recombination and recombination cooling for every cell
  sx_chem_initialize_site(site);
#else
  if(hasPhotons)
    sx_chem_initialize_site(site);
#endif

  return hasPhotons;
}

/*!
 * \brief compute the mean delaunay neighbour distance of the particle
 */
void sx_phot_compute_neigh_distance(int index)
{
  double dist = 0.0;
  int nDist = 0;

  int q = SphP[index].first_connection;
  while(q >= 0)
    {
      // get the neighbour index
      int dp = DC[q].dp_index;
      int j = Mesh.DP[dp].index;

      /* Check if we have a valid particle, otherwise go to the next */
      if(j < 0)
        {
          q = DC[q].next;
          continue;
        }

#ifndef SX_OUTPUT_IMAGE_FLAGS
      // DEBUG: for now filter out periodic border particles (discuss with Volker Springel)
      if (DC[q].image_flags == SX_DC_IMAGE_FLAGS) 
      {
#endif

          // calculate the distance of the two points
          // NOTE: we have to use Mesh.DP.x/y/z values because some neighbours are on different processors
          // DEBUG: This works only with SX_NDIM == 3
          myassert(SX_NDIM == 3);
          dist += (double)sqrt((P[index].Pos[0] - Mesh.DP[dp].x) * (P[index].Pos[0] - Mesh.DP[dp].x) +
                               (P[index].Pos[1] - Mesh.DP[dp].y) * (P[index].Pos[1] - Mesh.DP[dp].y) +
                               (P[index].Pos[2] - Mesh.DP[dp].z) * (P[index].Pos[2] - Mesh.DP[dp].z));
          nDist++;

#ifndef SX_OUTPUT_IMAGE_FLAGS      
      }
#endif

      // end with the last neighbour
      if(q == SphP[index].last_connection)
        break;

      q = DC[q].next;
    
    }/* while q >= 0 */
  
  //set the mean neighbour distance to the particle [code comoving units]
  sxCell[index].MeanNeighDistance = dist / (double)nDist; 
  
}


/*!
 * \brief Check if the total number of photons is bigger than allowed minimum
 */
_Bool sx_phot_check_n_photons(struct sxSite *site)
{
  // bool determines whether photons need to be transported
  _Bool doTransport = 0;
  
  // do not perform ballistic transport if there is not enough photons
  for(int f=0; f<SX_NFREQ; f++ )
    {
      if(isnan(site->nPhotons[f]))
        {
          printf("SX-ERROR: number of photons is NaN for the following site:\n");
          sx_debug_particle(site->index, P[site->index].ID);
          sx_debug_site(site, P[site->index].ID);
          terminate("SX-ERROR: number of photons is NaN");
        }
      if(site->nPhotons[f] > 0.0)
        {
          doTransport = 1;
          break;
        }
    }

  return doTransport;
}

//======================================================================================================//
//                                     Directional bins                                                 //
//======================================================================================================//

/*!
 * \brief Initialize global set of directional bins
 *
 * There are two possible variants for setting the global set of directional bins.
 * Either we will use a fixed set of direction for the whole simulation (SX_NUM_ROT=0),
 * or we will randomly rotate the set before each run (SX_NUM_ROT>0)
 * NOTE: In a case when we want to store actual number of photons in the cell during the hydrostep
 *       and subsequently use them in the following RT run we need to keep global set fixed at all times!!!
 */
void sx_phot_initialize_global_directions(void)
{
#if (SX_NUM_ROT>0)
  double M[3][3];    // random rotation matrix 3x3
  double base[3][3]; // list of global base vectors (M * v = rotated vector) of a form: [base vector x,y,z][coordinate]
  sx_math_compute_rotational_matrix(M, M_PI);
#elif defined(SX_ROT_EVERY_STEP)
  double M[3][3];
  double base[3][3];
  double angle=0;
  for (int c=0; c<SX_NDIR; c++)
    angle += sxSolidAngles_[c];
  angle = acos(1.0 - angle / SX_NDIR / (2.0 * M_PI));
  sx_math_compute_rotational_matrix(M, angle);  // rotating only within the solid angle of the directional bin
                                                // mpi_printf("SX: rotation angle %f M [[%f %f %f],[%f %f %f],[%f %f %f]] \n",angle,
                                                //	 M[0][0],M[0][1],M[0][2],
                                                //	 M[1][0],M[1][1],M[1][2],
                                                //	 M[2][0],M[2][1],M[2][2]);
#endif

  for(int d=0; d<SX_NDIR; d++ ) // loop through all directions
    {
#if(SX_NUM_ROT > 0) || defined(SX_ROT_EVERY_STEP)
      // set the orthogonal base
      for(int c=0; c<SX_NDIM; c++) // loop through coordinates
	{
	  base[0][c] = sxOrthoBase_[d][0][c]; // x-vector
	  base[1][c] = sxOrthoBase_[d][1][c]; // y-vector
          base[2][c] = sxDirectionBins_[d][c]; // z-vector = direction
	}
      // rotate the base using the random rotation matrix
      for(int v=0; v<SX_NDIM; v++) // loop through vectors
	{
          sxRunData->orthoBase[d][v][0] = M[0][0] * base[v][0] + M[0][1] * base[v][1] + M[0][2] * base[v][2];
          sxRunData->orthoBase[d][v][1] = M[1][0] * base[v][0] + M[1][1] * base[v][1] + M[1][2] * base[v][2];
          sxRunData->orthoBase[d][v][2] = M[2][0] * base[v][0] + M[2][1] * base[v][1] + M[2][2] * base[v][2];
        }
#else
      // set the orthogonal base
      for(int c=0; c<SX_NDIM; c++) // loop through coordinates
	{
	  sxRunData->orthoBase[d][0][c] = sxOrthoBase_[d][0][c]; // x-vector 
	  sxRunData->orthoBase[d][1][c] = sxOrthoBase_[d][1][c]; // y-vector
          sxRunData->orthoBase[d][2][c] = sxDirectionBins_[d][c]; // z-vector = direction
	}
#endif      
    }
}

//======================================================================================================//
//                                         Run flow functions                                           //
//======================================================================================================//

void sx_phot_initialize(void)
{
}

void sx_phot_start_run(void)
{
  // set solid angles
  for(int d=0; d<SX_NDIR; d++ ) 
    {
      sxRunData->solidAngles[d] = sxSolidAngles_[d];
    }
#ifdef SX_DISPLAY_STATS
  mpi_fprintf(FdSimplex,"SA theta is %.02f * %.1f = %.02f (%.02f deg) \n", 
	     sxRunData->solidAngles[0], SX_FACTOR_SA, 
	     sxRunData->solidAngles[0] * SX_FACTOR_SA,
	     sxRunData->solidAngles[0] * SX_FACTOR_SA * 180 / M_PI);
#endif
  
  // initialize a global set of directional bins
  sx_phot_initialize_global_directions();

  // loop over active particles
  for(int d = 0; d < TimeBinsHydro.NActiveParticles; d++)
    {
      int i = TimeBinsHydro.ActiveParticleList[d];
      if(i < 0)
        continue;

      /* Calculate the mean neighbour distance of every gas particle */
      sx_phot_compute_neigh_distance(i);

#ifdef SX_OUTPUT_FLUX
      memset(SphP[i].sxPhotonFlux, 0, SX_NFREQ * sizeof(double));
#endif
    }
  // endrun();

  // initialize photon transport
  sx_phot_transport_start_run();

  // calculate ionistion rates of sources and assign them to the corresponding gas cells
  sx_phot_sources_start_run();

  // initialize variables for photon exchange between processors
  sx_phot_exchange_start_run();

  // Allocate all movable memory at once
  // NOTE: we need to do it at this point, because movable arrays have to be allocated at the end
  // NOTE: there are some non-movable arrays allocated below (during steps), however, they are freed before any reallocation
  sxRealAPP = (struct sxABPP_struct *)mymalloc_movable(&sxRealAPP,"sxAPP", sxMaxNumABPP * sizeof(struct sxABPP_struct));
  sxRealBPP = (struct sxABPP_struct *)mymalloc_movable(&sxRealBPP,"sxBPP", sxMaxNumABPP * sizeof(struct sxABPP_struct));
  sxAPP = sxRealAPP; 
  sxBPP = sxRealBPP;
  sxQPP = (struct sxQPP_struct *)mymalloc_movable(&sxQPP,"sxQPP", sxMaxNumQPP * sizeof(struct sxQPP_struct));

#ifdef SX_DISPLAY_STATS
  // clear the photon statistics
  memset(sxRunData->nInitPhotons, 0, SX_NFREQ * sizeof(double));
  memset(sxRunData->nEmitPhotons, 0, SX_NFREQ * sizeof(double));
  memset(sxRunData->nAbsPhotons, 0, SX_NFREQ * sizeof(double));
  memset(sxRunData->nLostPhotons, 0, SX_NFREQ * sizeof(double));

  // count initial number of photons in the system at the beginning of the run
  for(int spd=0; spd<sxNumSPD; spd++)
    {
      int app = sxSPDtoAPP[spd];
      if (app<0) continue;
      for(int f=0; f<SX_NFREQ; f++ )
	  sxRunData->nInitPhotons[f] += (double)sxAPP[app].nPhotons[f];
    }
  for(int f=0; f<SX_NFREQ; f++ )
    sxRunData->nInitPhotons[f] *= All.UnitPhotons_per_s;
#endif
  
}

#if (SX_NUM_ROT>1)
void sx_phot_start_rot(void)
{
  sx_debug_timer_start(TROT);
#ifdef SX_DISPLAY_TIMERS 
  sx_debug_timer_start(TROTs);
#endif

#ifdef SX_DISPLAY_STATS
  mpi_fprintf(FdSimplex,"\n>--- RUN %d ROT %d/%d \n", 
	      All.sxCurrentRun, sxRunData->currentRot+1, SX_NUM_ROT); 
#endif

  // initialize a global set of directional bins
  sx_phot_initialize_global_directions();

  // Reset the APP and BPP arrays
  memset(sxSPDtoAPP, -1, sxNumSPD * sizeof(int));
  sxNumAPP = 0;
  memset(sxSPDtoBPP, -1, sxNumSPD * sizeof(int));
  sxNumBPP = 0;

#ifdef SX_DISPLAY_TIMERS
  sx_debug_timer_end(TROTs);
#endif
}
#endif

void sx_phot_start_step(void)
{
#ifdef SX_ROT_EVERY_STEP
  sx_phot_initialize_global_directions();
#endif
}

void sx_phot_end_step(void)
{
#ifdef SX_DISPLAY_TIMERS
  sx_debug_timer_start(TEXCHc);
#endif

   // transport photon packages
  sx_phot_transport_copy();

#ifdef SX_DISPLAY_TIMERS
  sx_debug_timer_switch(TEXCHc,TEXCHe);
#endif

  // exchange photon packages between processors
  sx_phot_exchange_photons();

#ifdef SX_DISPLAY_TIMERS
  sx_debug_timer_end(TEXCHe);
#endif

#ifdef SX_DISPLAY_STATS
  // count remaining number of photons in the system
  memset(sxRunData->nFinPhotons, 0, SX_NFREQ * sizeof(double));
  for(int spd=0; spd<sxNumSPD; spd++)
    {
      int app = sxSPDtoAPP[spd];
      if (app<0) continue;
      for(int f=0; f<SX_NFREQ; f++ )
	sxRunData->nFinPhotons[f] += (double)sxAPP[app].nPhotons[f];
    }
  for(int f=0; f<SX_NFREQ; f++ )
    sxRunData->nFinPhotons[f] *= All.UnitPhotons_per_s;
#endif
}

#if (SX_NUM_ROT>1)
void sx_phot_end_rot(void)
{
#ifdef SX_DISPLAY_TIMERS
  sx_debug_timer_start(TROTe);
#endif

  // end the timer and print the information
#ifdef SX_DISPLAY_TIMERS
  sx_debug_timer_end(TROTe);
#endif
  double trot = sx_debug_timer_end(TROT);
  mpi_fprintf(FdSimplex,"\n>--- RUN %d ROT %d/%d --- %d steps took %e s \n", 
	      All.sxCurrentRun, sxRunData->currentRot+1, SX_NUM_ROT, sxRunData->currentStep, trot); 
}
#endif

void sx_phot_end_run(void)
{

  // free movable memory
  myfree_movable(sxQPP);
  myfree_movable(sxRealBPP);
  myfree_movable(sxRealAPP);

  // free non-movable memory
  sx_phot_exchange_end_run();
  sx_phot_sources_end_run();
  sx_phot_transport_end_run();

}
