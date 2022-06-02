/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \date        22/10/2018
 * \author      Ondrej Jaura
 * \brief       Transport of the photons on the Arepo Voronoi mesh
 * \details     
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

//======================================================================================================//
//                                     Photon Package functions                                         //
//======================================================================================================//

/*!
 * \brief Reallocate ABPP memory
 */
void sx_phot_reallocate_ABPP(void)
{
  int newNum = (int)( sxMaxNumABPP * SX_ABPP_REALLOC );
  if (newNum>sxNumSPD) newNum = sxNumSPD;
  //printf("SX: (%d) ABPP realloc max %d old %d new %d sxSwapABPP %d \n", ThisTask, 
  // 	 sxNumSPD, sxMaxNumABPP, newNum, sxSwapABPP);
  sxMaxNumABPP = newNum;
  sxRealAPP = (struct sxABPP_struct *)myrealloc_movable(sxRealAPP, sxMaxNumABPP * sizeof(struct sxABPP_struct));
  sxRealBPP = (struct sxABPP_struct *)myrealloc_movable(sxRealBPP, sxMaxNumABPP * sizeof(struct sxABPP_struct));
  sxAPP = (sxSwapABPP) ? sxRealBPP : sxRealAPP; 
  sxBPP = (sxSwapABPP) ? sxRealAPP : sxRealBPP;
#ifdef SX_DISPLAY_MEMORY
  sxNumReallocABPP++;
#endif
}

/*!
 * \brief Copy PP value to the APP array
 */
void sx_phot_copy_PP_to_APP( int i, int d, MySxDouble nPhotons[SX_NFREQ] )
{
  int f, spd, app;
  spd = i * SX_NDIR + d;
  if (sxSPDtoAPP[spd]<0) 
    {
      app = sxNumAPP++;	  
      if (app==sxMaxNumABPP) sx_phot_reallocate_ABPP();
      sxSPDtoAPP[spd] = app;	  	  
      for (f=0; f<SX_NFREQ; f++)
	  sxAPP[app].nPhotons[f] = nPhotons[f]; // assigning a new value !!
    }
  else
    {
      app = sxSPDtoAPP[spd];	  	  
      for (f=0; f<SX_NFREQ; f++)
	  sxAPP[app].nPhotons[f] += nPhotons[f]; // adding photons !!
    }
  /*
    printf("SX: (%d) IN  i %d d %d app %d nPhotons [ ",ThisTask,i,d,app);
    for (f=0; f<SX_NFREQ; f++)
      printf("%.03e ",nPhotons[f]);
    printf("] \n");
  */
}

/*!
 * \breif Copy PP value to the BPP array
 */
void sx_phot_copy_PP_to_BPP( int i, int d, MySxDouble nPhotons[SX_NFREQ] )
{
  int f, spd, bpp;
  spd = i * SX_NDIR + d;
  if (sxSPDtoBPP[spd]<0) 
    {
      bpp = sxNumBPP++;	  
      if (bpp==sxMaxNumABPP) sx_phot_reallocate_ABPP();
      sxSPDtoBPP[spd] = bpp;	  	  
      for (f=0; f<SX_NFREQ; f++)
	  sxBPP[bpp].nPhotons[f] = nPhotons[f]; // assigning a new value !!
    }
  else
    {
      bpp = sxSPDtoBPP[spd];	  	  
      for (f=0; f<SX_NFREQ; f++)
	  sxBPP[bpp].nPhotons[f] += nPhotons[f]; // adding photons !!
    }
  /*
    printf("SX: (%d) IN  i %d d %d bpp %d nPhotons [ ",ThisTask,i,d,bpp);
    for (f=0; f<SX_NFREQ; f++)
      printf("%.03e ",nPhotons[f]);
    printf("] \n");
  */
}

/*!
 * \brief Add photons to the local/external photon packet list
 */
void sx_phot_add_photons_next( int dp, int d, MySxDouble nPhotons[SX_NFREQ] )
{
  int i = Mesh.DP[dp].index;
  int task = Mesh.DP[dp].task;

  if( task == ThisTask )
    sx_phot_copy_PP_to_BPP( i, d, nPhotons );
  else
    sx_phot_copy_PP_to_QPP( dp, d, nPhotons );

}

/*!
 * \brief Copy/swap photons from sxBPP to sxAPP
 */
void sx_phot_transport_copy(void)
{
  int i, d, f, bpp;

  // Swap the sxAPP/sxBPP pointers
  struct sxABPP_struct *tmp;
  tmp = sxAPP; 
  sxAPP = sxBPP;
  sxBPP = tmp;

  int *tmp2;
  tmp2 = sxSPDtoAPP;
  sxSPDtoAPP = sxSPDtoBPP;
  sxSPDtoBPP = tmp2;

  sxNumAPP = sxNumBPP;

  sxSwapABPP = (sxSwapABPP) ? 0 : 1;

  // Reset SPD and BPP arrays
  memset(sxSPDtoBPP, -1, sxNumSPD * sizeof(int));
  sxNumBPP = 0;

  //MPI_Barrier(MPI_COMM_WORLD);
  //MPI_Abort(MPI_COMM_WORLD, 1);

}

//======================================================================================================//
//                                     Common functions                                                 //
//======================================================================================================//

/*! 
 * \brief Computes a random direction in the opening angle of the direction bin
 *
 * Computes a random number on a circle centred on the direction bin, with radius 
 * equal to f * opening_angle_dir_bin. Here, f is a fraction that can be set by 
 * the user to ensure the opening angle is properly filled. 
 * Recommended value is f ~ 5.
 *
 * \param[in] binIndex  Index of the direction bin
 *
 * \param[out] directionBin  Coordinates of a random direction 
 * within the opening angle
 */
void sx_phot_compute_random_direction_bin( int binIndex, double directionBin[SX_NDIM] )
{
  
  int dim;
  double phi, cosTheta, randTheta, x, y, z, xy;
  double binVector[SX_NDIM];
        
  for( dim=0; dim<SX_NDIM; dim++ )
    {
      binVector[dim] = sxRunData->orthoBase[binIndex][2][dim]; //sxDirectionBins_[binIndex][dim];
    }

  // NOTE: "solidAngles" actually contains theta angles in radians, 
  //        however, the values for small number of bins are not accurate
  phi = 2 * M_PI * gsl_rng_uniform( sxRand );                          // random from [ 0, 2*pi ]
  //cosTheta = cos( All.sxDirConeAngle * 0.017453292519943295 );       // cos( theta_deg * pi/180 )
  cosTheta = cos( SX_FACTOR_SA * sxRunData->solidAngles[binIndex] );   // cos( fSA * theta_rad )
  z = cosTheta + (1.0 - cosTheta) * gsl_rng_uniform( sxRand );         // random from [ cos(theta), 1 ]
  xy = sqrt( 1.0 - pow(z,2) );                                         //  <=  1^2 = z^2 + xy^2
  y = xy*sin(phi);                                                     // sqrt(1-z^2)*sin(phi)
  x = xy*cos(phi);                                                     // sqrt(1-z^2)*cos(phi)

  //use orthogonal base to get vector randomly rotated in cone
  for( dim=0; dim<SX_NDIM; dim++ )
    {
      directionBin[dim] = x*sxRunData->orthoBase[binIndex][0][dim] + //sxOrthoBase_[binIndex][0][dim] + 
	y*sxRunData->orthoBase[binIndex][1][dim] + //sxOrthoBase_[binIndex][1][dim] + 
	z*sxRunData->orthoBase[binIndex][2][dim]; //binVector[dim];
    }

}

//======================================================================================================//
//                                     Diffuse transport                                                //
//======================================================================================================//

/*!
 * \brief Compute the neighbour index associated with direction bin
 */
int sx_phot_compute_neighbours_index( struct sxSite *site , double directionBin[SX_NDIM] )
{
  //printf("SX: (%d:%d:%d) pos [ %.03e %.03e %.03e ] \n",
  //	 ThisTask, site->index, P[site->index].ID,
  //	 P[site->index].Pos[0], P[site->index].Pos[1], P[site->index].Pos[2]);

  //maximum inner product between this direction bin and neighbours
  double maxInnerProduct = 0.0;
  
  //local neighbour index belonging to maximum inner product
  int neighDPIndex = -1;

  int dim, q, dp, j;
  
  q = SphP[site->index].first_connection;
  while(q >= 0)
    {
      dp = DC[q].dp_index;
      j = Mesh.DP[dp].index;
  
      /* Check if we have a valid particle, otherwise go to the next */
      if(j < 0)
        {
	  q = DC[q].next;
	  continue;
        }

      // create a neighbour vector
      // DEBUG: This works only if SX_NDIM == 3
      double neighVector[SX_NDIM];
      myassert(SX_NDIM == 3);  // If SX_NDIM != 3, better to fail here than compute garbage
      neighVector[0] = Mesh.DP[dp].x - P[site->index].Pos[0];
      neighVector[1] = Mesh.DP[dp].y - P[site->index].Pos[1];
      neighVector[2] = Mesh.DP[dp].z - P[site->index].Pos[2];
      
      // calculate inner product between the refVector and the current neighbour
      double innerProduct = sx_math_compute_inner_product( neighVector, directionBin );

      // store the highest inner product and its delaunay point index
      if ( innerProduct > maxInnerProduct )
	{
          maxInnerProduct = innerProduct;
#ifndef SX_OUTPUT_IMAGE_FLAGS
	  // DEBUG: for now we filter out any ghost particles (outside of the domain, etc.)....
	  // TODO: discuss flags with Volker Springel
	  if (DC[q].image_flags > SX_DC_IMAGE_FLAGS)
	    neighDPIndex = -1;
	  else
	    neighDPIndex = dp;
#else
	  neighDPIndex = dp;
#endif      
        }
      
      // end with the last neighbour
      if(q == SphP[site->index].last_connection)
	break;

      q = DC[q].next;
    
    }/* while q >= 0 */

  /*
  if( neighDPIndex < 0 )
    {
      printf("SX: (%d:%d:%d) Neighbour not found!!! Pos [ %.03e %.03e %.03e ] DirBin [ %.03e %.03e %.03e ]\n",
	     ThisTask,site->index,P[site->index].ID,
	     P[site->index].Pos[0], P[site->index].Pos[1], P[site->index].Pos[2],
	     directionBin[0], directionBin[1], directionBin[2]);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  */
 
  return neighDPIndex;
        
}

/*! 
 * \brief Diffuse transport of photons
 *
 * Transports photons to all neighbours in the Delaunay triangulation
 *
 * \param[in] particleIndex  place in the SphP array of the particle 
 * containing the photons
 * \param[in] nPhotons       number of photons to distribute
 *
 */
void sx_phot_diffuse_transport( struct sxSite *site )
{	  
  int d, f;
  MySxDouble NIon[SX_NFREQ];
  double splitFactor = 1. / ( (double)SX_NDIR * All.UnitPhotons_per_s );

  // Compute the number of photons to be sent in every direction bin 
  for( f=0; f<SX_NFREQ; f++ )
    {     
      if (site->nPhotonsOriginal[f]==0)
	NIon[f] = 0.0;
      else
	NIon[f] = (MySxDouble)(site->nPhotonsSource[f] * site->nPhotonsRatio[f] * splitFactor);
    }

  // Distribute the photons in every direction bin to the closest 
  //   Delaunay connection. Randomize the direction bins in their
  //   opening angle.  
  for( d=0; d<SX_NDIR; d++ )
    {
      //vector of coordinates of direction bin
      double directionBin[SX_NDIM];

      //get a direction in the solid angle of the initial bin
      sx_phot_compute_random_direction_bin( d, directionBin );
    
      //associate this direction with neighbours
      int neighDPIndex = sx_phot_compute_neighbours_index( site, directionBin );

      if (neighDPIndex>-1)
	{
	  //add the radiation to the site
	  sx_phot_add_photons_next( neighDPIndex, d, NIon );
	}
#ifdef SX_DISPLAY_STATS
      else
	{
	  // all photons are lost
	  for( f=0; f<SX_NFREQ; f++ )
	    {
	      sxRunData->nLostPhotons[f] += (double)NIon[f];
	    }
	}
#endif	
    }

}

//======================================================================================================//
//                                     Ballistic transport                                              //
//======================================================================================================//

/*!
 * \brief compute the most STRAIGHTFORWARD neighbours with respect to the bin
 */
int sx_phot_compute_straight_neighbours( struct sxSite *site , double directionBin[SX_NDIM], int neighDPIndex[SX_NDIM] )
{

  int dim, n=0;
  double max[SX_NDIM]; // array to hold the neighbour vector and the maximum inner product
  int DPIndex[SX_NDIM];  // array to hold the indices of the most straight forward neighbours

  // initialise values
  for( dim=0; dim<SX_NDIM; dim++ )
    {
      max[dim] = -FLT_MAX;
      DPIndex[dim] = -1;
      neighDPIndex[dim] = -1;
    }

  // loop over all neighbours to calculate the most strait-forward ones with respect to the direction bin
  int q, dp, j;
  
  q = SphP[site->index].first_connection;
  while(q >= 0)
    {
      dp = DC[q].dp_index;
      j = Mesh.DP[dp].index;
  
      /* Check if we have a valid particle, otherwise go to the next */
      if(j < 0)
        {
	  q = DC[q].next;
	  continue;
        }

      // create a neighbour vector
      // DEBUG: This works only if SX_NDIM == 3
      double neighVector[SX_NDIM];
      myassert(SX_NDIM == 3);  // If SX_NDIM != 3, better to fail here than compute garbage
      neighVector[0] = Mesh.DP[dp].x - P[site->index].Pos[0];
      neighVector[1] = Mesh.DP[dp].y - P[site->index].Pos[1];
      neighVector[2] = Mesh.DP[dp].z - P[site->index].Pos[2];
      
      // calculate inner product between both vectors
      double innerProduct = sx_math_compute_inner_product( neighVector, directionBin );

#ifndef SX_OUTPUT_IMAGE_FLAGS
      // DEBUG: for now we filter out any ghost particles (outside of the domain, etc.)....
      // TODO: discuss flags with Volker Springel
      if (DC[q].image_flags > SX_DC_IMAGE_FLAGS)
	{
	  dp = -1;
	}
#endif
      
      // store d largest inner products
      if ( innerProduct > max[0] )
	{
	  max[2] = max[1];
	  max[1] = max[0];
	  max[0] = innerProduct;
	  DPIndex[2] = DPIndex[1];
	  DPIndex[1] = DPIndex[0];
	  DPIndex[0] = dp;
	}
      else if( innerProduct > max[1] )
	{
	  max[2] = max[1];
	  max[1] = innerProduct;
	  DPIndex[2] = DPIndex[1];
	  DPIndex[1] = dp;
	}
      else if( innerProduct > max[2] )
	{
	  max[2] = innerProduct;
	  DPIndex[2] = dp;
	}


      // end with the last neighbour
      if(q == SphP[site->index].last_connection)
	break;

      q = DC[q].next;
    
    }/* while q >= 0 */

  //store most straightforward neighbour
  //DEBUG: restricting now only on one most straightforward neighbor!! to see the how ray works
  //for( dim=0; dim<SX_NDIM; dim++ )
  for( dim=0; dim<SX_NDIM; dim++ )
    {
      //add the neighbour if the maximum inner product is larger than the cosine of the largest angle
      if( max[dim] > All.sxMaxStraightAngleCos )
	{
	  neighDPIndex[dim] = DPIndex[dim];
	  n++;
	}
      else
	{
	  break;
	}
    }

  return n;
  
}

/*! 
 * \brief Ballistic transport of photons
 *
 * Transports photons to the NDIM neighbours that form the most straightforward
 * direction relative to the incoming direction. Incoming directions are sampled
 * randomly in the opening angle of the tessellation of the unit sphere. 
 *
 * \param[in] particleIndex          place in the SphP array of the particle containing the photons
 * \param[in] nPhotons       number of photons to distribute
 *
 */
void sx_phot_ballistic_transport( struct sxSite *site )
{
  for(int d=0; d<SX_NDIR; d++ )
    {
      // Check if there is a APP in a particular direction
      int app = sxSPDtoAPP[ site->index * SX_NDIR + d ];
      if ( app >= 0 )
	{
	  //get a direction in the solid angle of the initial bin
	  double directionBin[SX_NDIM]; //vector of coordinates of direction bin
          sx_phot_compute_random_direction_bin( d, directionBin );

	  // compute most straightforward neighbours
	  int neighDPIndex[SX_NDIM]; //indexes of the most straight neighbours (negative index -> no neighbour)
	  int numNeigh = sx_phot_compute_straight_neighbours( site, directionBin, neighDPIndex );
	  double invNumNeigh = 1.0/(double)numNeigh;

	  //DEBUG: if there are no straight neighbours (i.e. at the box border), we do not continue
          if (numNeigh > 0)
	    {
	      MySxDouble NIon[SX_NFREQ];
  
	      //compute the intensity to send out
	      for(int f=0; f<SX_NFREQ; f++ )
		{
		  //get correct part from all added intensities and divide by the number of straight neighbours
		  NIon[f] = (MySxDouble)( (double)sxAPP[app].nPhotons[f] * invNumNeigh * site->nPhotonsRatio[f] );
		}

	      //add the raditation to straightforward neighbours
	      for(int n=0; n<numNeigh; n++ )
		{
#ifndef SX_OUTPUT_IMAGE_FLAGS
		  // filtering out particles on the border
		  if (neighDPIndex[n]<0)
		    {
#ifdef SX_DISPLAY_STATS
		      for(int f=0; f<SX_NFREQ; f++ )
			{
			  sxRunData->nLostPhotons[f] += (double)(NIon[f] * All.UnitPhotons_per_s);
			}
#endif
		      continue;
		    }
#endif
		  //add the radiation to the site
		  sx_phot_add_photons_next( neighDPIndex[n], d, NIon );
		}
	    }
#ifdef SX_DISPLAY_STATS
	  else
	    {
	      // all photons are lost
	      for(int f=0; f<SX_NFREQ; f++ )
		{
		  sxRunData->nLostPhotons[f] += (double)sxAPP[app].nPhotons[f] * site->nPhotonsRatio[f] * All.UnitPhotons_per_s;
		}
	    }
#endif	
	}
    }

}

//======================================================================================================//
//                                     Runtime functions                                                //
//======================================================================================================//

void sx_phot_transport_start_run(void)
{
  // Number of all SphP particle Directions
  sxNumSPD = NumGas * SX_NDIR;

  // Initialize counters
  sxMaxNumABPP = (int)( sxNumSPD * SX_ABPP_INIT );
  if (All.sxCurrentRun>1 && All.sxLastMaxNumABPP > sxMaxNumABPP && All.sxLastMaxNumABPP < sxNumSPD)
      sxMaxNumABPP = All.sxLastMaxNumABPP;
  sxNumBPP = 0;  
  sxNumAPP = 0;  
  sxSwapABPP = 0;
#ifdef SX_DISPLAY_MEMORY
  sxNumReallocABPP = 0;
#endif

  // Initialize an array that will keep indexes of the transported photon packages
  sxRealSPDtoAPP = (int *)mymalloc("sxRealSPDtoAPP", sxNumSPD * sizeof(int));
  sxRealSPDtoBPP = (int *)mymalloc("sxRealSPDtoBPP", sxNumSPD * sizeof(int));
  sxSPDtoAPP = sxRealSPDtoAPP;
  sxSPDtoBPP = sxRealSPDtoBPP;
  memset(sxSPDtoAPP, -1, sxNumSPD * sizeof(int));  
  memset(sxSPDtoBPP, -1, sxNumSPD * sizeof(int));

}

void sx_phot_transport_end_run(void)
{
  // Save the current sxMaxNumABPP for the next run
  All.sxLastMaxNumABPP = sxMaxNumABPP;

  // Free all memory
  myfree(sxRealSPDtoBPP);
  myfree(sxRealSPDtoAPP);
}
