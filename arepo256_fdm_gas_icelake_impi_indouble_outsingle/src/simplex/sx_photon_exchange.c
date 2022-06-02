/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \date        22/10/2018
 * \author      Ondrej Jaura
 * \brief       Photon exchange between processors
 * \details     
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

//======================================================================================================//
//                                     Photon exchange                                                  //
//======================================================================================================//

/*!
 * \brief Add photon packet to exchange queue
 */
void sx_phot_copy_PP_to_QPP( int dp, int d, MySxDouble nPhotons[SX_NFREQ] )
{
  int f, app, newNum;
  int task = Mesh.DP[dp].task;
  int edpd = ( sxOffEDP[task] + sxDPtoEDP[dp] ) * SX_NDIR + d;

  //if (sxNumExportPP>5) return;

  if ( sxEDPDtoQPP[edpd] < 0 ) // create a new packet
    {
      app = sxNumExportPP++;

      if (app==sxMaxNumQPP) // increase the memory if needed
	{
	  newNum = (int)( sxMaxNumQPP * SX_QPP_REALLOC );
	  if (newNum>sxNumTotEDPD)  newNum = sxNumTotEDPD; // constrain the maximum
	  //printf("SX: (%d) QPP realloc max %d old %d new %d \n",ThisTask,
	  // 	 sxNumTotEDPD, sxMaxNumQPP, newNum);
	  sxMaxNumQPP = newNum;
	  sxQPP = (struct sxQPP_struct *)myrealloc_movable(sxQPP, sxMaxNumQPP * sizeof(struct sxQPP_struct));
#ifdef SX_DISPLAY_MEMORY
	  sxNumReallocQPP++;
#endif
	  //MPI_Abort(MPI_COMM_WORLD, 1);
	}
      sxEDPDtoQPP[edpd] = app;
      for (f=0; f<SX_NFREQ; f++ )
	{
	  sxQPP[app].nPhotons[f] = nPhotons[f];  // assigning a new value !!
	}
    }
  else // add to the existing packet
    {
      app = sxEDPDtoQPP[edpd];
      for (f=0; f<SX_NFREQ; f++ )
	{
	  sxQPP[app].nPhotons[f] += nPhotons[f]; // adding new photons !!
	}
    }

  sxQPP[app].index = Mesh.DP[dp].originalindex;

  /*
  if (sxRunData->currentStep==0)
    {
      printf("SX: (%d) dp %d d %d task %d sxOffEDP %d sxDPtoEDP %d SX_NDIR %d edpd %d app %d sxNumExportPP %d origIndex %d \n",
	     ThisTask,
	     dp, d, task, sxOffEDP[task], sxDPtoEDP[dp], SX_NDIR, edpd, app, sxNumExportPP, Mesh.DP[dp].originalindex);
    }
  */
  /*
  if (sxNumExportPP>8)
    MPI_Abort(MPI_COMM_WORLD, 1);
  */
}

/*!
 * \brief Exchange photons between processors
 */
void sx_phot_exchange_photons(void)
{
  int edpd, app, i, t, d, p, f;
  int ngrp;
  // Reset counter
  for(t = 0; t < NTask; t++){
    Send_count[t] = 0;
  }

  // Create an export photon packet list
  sxExportPP = (struct sxEPP_struct *)mymalloc("sxExportPP", sxNumExportPP * sizeof(struct sxEPP_struct));
  for( t=0, i=0; t<NTask; t++)
    {
      for( d=0; d<SX_NDIR; d++)
	{
	  for( p=0; p<sxNumEDP[t]; p++)
	    {
	      edpd = ( sxOffEDP[t] + p ) * SX_NDIR + d;
	      app = sxEDPDtoQPP[edpd];

	      if (app<0) continue;

	      sxExportPP[i].dir = d;
	      sxExportPP[i].index = sxQPP[app].index;
	      for ( f=0; f<SX_NFREQ; f++ )
		  sxExportPP[i].nPhotons[f] = sxQPP[app].nPhotons[f];
	      
	      Send_count[t]++;
	      i++;
	      /*
	      if (sxRunData->currentStep==0)
		{
		  printf("SX: (%d) t %d d %d p %d edpd %d app %d sxNumEDP[t] %d \n",ThisTask,
			 t, d, p, edpd, app, sxNumEDP[t]);
		}
	      */
	    }
	}
    }

  // Count items that we will receive
  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  // Calculate exchange variables
  for(t = 0, sxNumImportPP = 0, Recv_offset[0] = 0, Send_offset[0] = 0; t < NTask; t++)
    {
      sxNumImportPP += Recv_count[t];
      if(t > 0)
        {
          Send_offset[t] = Send_offset[t - 1] + Send_count[t - 1];
          Recv_offset[t] = Recv_offset[t - 1] + Recv_count[t - 1];
        }
    }

  /*
  if (sxRunData->currentStep==0 && (ThisTask==10||ThisTask==4||ThisTask==8||ThisTask==15))
    {
      printf("SX: (%d) sxNumImportPP %d sxNumExportPP %d \n      Send_count [ ",ThisTask,sxNumImportPP,sxNumExportPP);
      for(t = 0; t < NTask; t++){
	printf("%d ",Send_count[t]);
      }
      printf("] \n      Recv_count [ ");
      for(t = 0; t < NTask; t++){
	printf("%d ",Recv_count[t]);
      }
      printf("] \n      Send_offset [ ");
      for(t = 0; t < NTask; t++){
	printf("%d ",Send_offset[t]);
      }
      printf("] \n      Recv_offset [ ");
      for(t = 0; t < NTask; t++){
	printf("%d ",Recv_offset[t]);
      }
      printf("] \n");
    }
  */

  // Create an array that will receive packages from another processors
  struct sxEPP_struct sxImportPP[sxNumImportPP];

  // Exchange photon packets
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              // get the particles
              MPI_Sendrecv(&sxExportPP[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct sxEPP_struct),
                           MPI_BYTE, recvTask, TAG_DENS_A,
                           &sxImportPP[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct sxEPP_struct),
                           MPI_BYTE, recvTask, TAG_DENS_A,
                           MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  // Free memory
  myfree(sxExportPP);

  // Add photons to the corresponding particles
  for( p=0; p<sxNumImportPP; p++)
    {
      i = sxImportPP[p].index;
      d = sxImportPP[p].dir;
      // add up the photons
      sx_phot_copy_PP_to_APP( i, d, sxImportPP[p].nPhotons );
    }

#ifdef SX_DISPLAY_STATS
  // save statistics about photon exchange
  sxRunData->nRecPackets = sxNumImportPP;
  sxRunData->nSendPackets = sxNumExportPP;
#endif

  // Reset the the packet list
  memset(sxEDPDtoQPP, -1, sxNumTotEDPD * sizeof(int));
  sxNumExportPP = 0;
  sxNumImportPP = 0;

}

//======================================================================================================//
//                                     Runtime functions                                                //
//======================================================================================================//

/*!
 * \brief Initialize variables used for photon exchange between processors
 */
void sx_phot_exchange_start_run(void)
{
  // Find all External DP and sort them by the Task
  sxNumEDP = (int *)mymalloc("sxNumEDP", NTask * sizeof(int));
  memset(sxNumEDP, 0, NTask * sizeof(int));
  sxDPtoEDP = (int *)mymalloc("sxDPtoEDP", Mesh.Ndp * sizeof(int));
  memset(sxDPtoEDP, -1, NTask * sizeof(int));

  for( int dp=0; dp<Mesh.Ndp; dp++ )
    {
#ifndef SX_OUTPUT_IMAGE_FLAGS
      // DEBUG: exclude all ghost particles
      if (Mesh.DP[dp].image_flags > SX_DP_IMAGE_FLAGS)
        continue;
#endif

      int task = Mesh.DP[dp].task;
      int index = Mesh.DP[dp].index;

      if (index < 0 || task == ThisTask)
        continue;

      sxDPtoEDP[dp] = sxNumEDP[task]++;  
    }

  // Calculate total numbers of EDP/EDPD and their offsets with respect to the Task
  sxOffEDP = (int *)mymalloc("sxOffEDP", NTask * sizeof(int));
  memset(sxOffEDP, 0, NTask * sizeof(int));
  sxNumTotEDP = 0;
  for( int i=0; i<NTask; i++ )
    {
      sxNumTotEDP += sxNumEDP[i];
      for( int j=0; j<i; j++ )
	{
	  sxOffEDP[i] += sxNumEDP[j];
	}
    }
  sxNumTotEDPD = sxNumTotEDP * SX_NDIR;

  // Create an array that will keep track of the connection between EDPD and QPP
  sxEDPDtoQPP = (int *)mymalloc("sxEDPDtoQPP", sxNumTotEDPD * sizeof(int));
  memset(sxEDPDtoQPP, -1, sxNumTotEDPD * sizeof(int));
  sxNumExportPP = 0;
  sxNumImportPP = 0;

  // Initialize variables for the QPP
  sxMaxNumQPP = (int)( sxNumTotEDPD * SX_QPP_INIT ); // setting the initial value
  if (All.sxCurrentRun>1 && All.sxLastMaxNumQPP > sxMaxNumQPP && All.sxLastMaxNumQPP < sxNumTotEDPD) 
    sxMaxNumQPP = All.sxLastMaxNumQPP;
#ifdef SX_DISPLAY_MEMORY
  sxNumReallocQPP = 0;
#endif
  
  //printf("SX: (%d) sxMaxNumQPP %d \n",ThisTask,sxMaxNumQPP);
  /*
  if (ThisTask==0)
    {
      printf("SX: (%d) sxNumTotEDP %d sxNumTotEDPD %d int %d \n",ThisTask,sxNumTotEDP,sxNumTotEDPD,sizeof(int));
      
      printf("SX: (%d) sxNumEDP [ ",ThisTask);
      for( i=0; i<NTask; i++ )
	{
	  printf("%d ",sxNumEDP[i]);
	}
      printf("] \n");
      
      printf("SX: (%d) sxOffEDP [ ",ThisTask);
      for( i=0; i<NTask; i++ )
	{
	  printf("%d ",sxOffEDP[i]);
	}
      printf("] \n");
    }
  */
  //MPI_Abort(MPI_COMM_WORLD, 1);
  //terminate(0);
}

/*!
 * \brief Finalize photon exchange and clean the memory
 */
void sx_phot_exchange_end_run(void)
{

  //printf("SX: (%d) sxMaxNumQPP %d sxNumExportPP %d \n",ThisTask,sxMaxNumQPP,sxNumExportPP);

  // Update number of Added PP
  All.sxLastMaxNumQPP = sxMaxNumQPP;

  // Free the memory
  myfree(sxEDPDtoQPP);
  myfree(sxOffEDP);
  myfree(sxDPtoEDP);
  myfree(sxNumEDP);

  //MPI_Abort(MPI_COMM_WORLD, 1);

}
