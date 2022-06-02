/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \date        22/10/2018
 * \author      Ondrej Jaura
 * \brief       Loading of photon sources and injecting of the photons to the cells
 * \details     
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

// Structures that store source properties
struct sxSourceIon_struct
{
  double Rate[SX_NFREQ];
} *sxSourceIon;
mesh_search_data * sxSourceCoord;

//======================================================================================================//
//                                     Ionisation from sources                                          //
//======================================================================================================//

/*
 * \brief Load source data
 * SX_SOURCES == 4 or 5  - positions are those of the stars/sink particles
 * SX_SOURCES == 10      - sources are taken from the binary file "test_sources.bin"
 * SX_SOURCES > 10       - sources are taken from the header file "sx_test_sources.h"
 */
#if (SX_SOURCES==10)

int sx_phot_sources_load(void)
{  
  // Open a binary file with sources
  //mpi_printf("SX: Reading test sources from the binary file \n");
  FILE *fd;
  if(!(fd = fopen(All.sxTestSrcFile, "r")))
    terminate("SX-ERROR: Could not find the binary file %s", All.sxTestSrcFile);

  // Reading offset
  int binOffset = 0;

  // Read the header
  int tsNumSrc=0, tsNumFreq=0, tsNumSigma=0, tsNumEnergy=0;
  my_fread(&tsNumSigma, sizeof(int), 1, fd);  // number of crossections (can be zero)
  my_fread(&tsNumEnergy, sizeof(int), 1, fd); // number of energies (can be zero)
  my_fread(&tsNumSrc, sizeof(int), 1, fd);    // number of sources (non-zero)
  my_fread(&tsNumFreq, sizeof(int), 1, fd);   // number of frequencies (non-zero)
  binOffset = 4*sizeof(int); // update offset
  if (tsNumFreq!=SX_NFREQ)
    terminate("SX-ERROR: Number of frequency bins of a test source %d has to be equal to %d", tsNumFreq, SX_NFREQ);
  if (tsNumSigma>0 && tsNumSigma!=SX_NSIGMA)
    terminate("SX-ERROR: Number of cross sections in the test source %d has to be equal to %d", tsNumSigma, SX_NSIGMA);
  if (tsNumEnergy>0 && tsNumEnergy!=SX_NENERGY)
    terminate("SX-ERROR: Number of energies in the test source %d has to be equal to %d", tsNumEnergy, SX_NENERGY);

#if SX_CHEMISTRY == 3
  // Read and change cross sections if set in the source file
  if (tsNumSigma>0)
    {
      fseek(fd, binOffset, SEEK_SET); 
      my_fread(All.sxSigma, sizeof(double), SX_NSIGMA, fd);  
      binOffset += sizeof(double)*SX_NSIGMA;  // update offset

      mpi_fprintf(FdSimplex,"Cross-sections\n cm^2  [ ");
      for(int f=0; f<SX_NSIGMA; f++)
	mpi_fprintf(FdSimplex,"%.2e ", All.sxSigma[f]);
      mpi_fprintf(FdSimplex,"] \n");
    }

  // Read and change energies if set in the source file
  if (tsNumEnergy>0)
    {
      fseek(fd, binOffset, SEEK_SET); 
      my_fread(All.sxEnergy, sizeof(double), SX_NENERGY, fd);  
      binOffset += sizeof(double)*SX_NENERGY;  // update offset

      mpi_fprintf(FdSimplex,"Photon energies\n erg   [ ");
      for(int f=0; f<SX_NENERGY; f++)
	mpi_fprintf(FdSimplex,"%.2e ", All.sxEnergy[f]);
      mpi_fprintf(FdSimplex,"]\n eV    [ ");
      for(int f=0; f<SX_NENERGY; f++)
	mpi_fprintf(FdSimplex,"%.2e ", All.sxEnergy[f]/ELECTRONVOLT_IN_ERGS);
      mpi_fprintf(FdSimplex,"]\n");
    }
#endif

  // Redistribute the sources to all tasks
  int tsOffset=0, nexport=0;
  for(int i=0; i < tsNumSrc; i++)
    if (i%NTask<ThisTask)
      tsOffset++;  // count memory offset of the local task
    else if (i%NTask == ThisTask)
      nexport++;   // count number of sources on the local task

  // Read coordinates
  double tsCoord[nexport][SX_NDIM];
  fseek(fd, binOffset+sizeof(double)*SX_NDIM*tsOffset, SEEK_SET); 
  my_fread(&tsCoord, sizeof(double)*SX_NDIM, nexport, fd);  
  binOffset += sizeof(double)*SX_NDIM*tsNumSrc;  // update offset

  // Read rates
  double tsRates[nexport][SX_NFREQ];
  fseek(fd, binOffset+sizeof(double)*SX_NFREQ*tsOffset, SEEK_SET); 
  my_fread(&tsRates, sizeof(double)*SX_NFREQ, nexport, fd);
  binOffset += sizeof(double)*SX_NFREQ*tsNumSrc;  // update offset

  // Allocate the memory and copy the data
  sxSourceCoord = (mesh_search_data *)mymalloc("sxSourceCoord", nexport * sizeof(mesh_search_data) );
  sxSourceIon = (struct sxSourceIon_struct *)mymalloc("sxSourceIon", nexport * sizeof(struct sxSourceIon_struct) );
  for(int i=0; i<nexport; i++)
    {
      sxSourceCoord[i].Pos[0] = tsCoord[i][0];
      sxSourceCoord[i].Pos[1] = tsCoord[i][1];
      sxSourceCoord[i].Pos[2] = tsCoord[i][2];
#ifndef DO_NOT_RANDOMIZE_DOMAINCENTER
      // we have to correct the source coordinates according to the current domain displacement
      // DEBUG: this still does not solve all problems, because we need to first implement
      //        periodic photon radiation
      domain_displacePosition( sxSourceCoord[i].Pos, DISPLACE_POSITION_FORWARD );
#endif
      for(int f=0; f<SX_NFREQ; f++)
	{
	  sxSourceIon[i].Rate[f] = tsRates[i][f];  // ionisation rate (ph/s)
	}
    }  

  // Close the file
  fclose(fd);

  return nexport;
}

#elif (SX_SOURCES>10) 

int sx_phot_sources_load(void)
{
  int nexport=0;

  // Include source settings
#include "sx_test_sources.h"

  // Count number of sources for export
  for(int i = 0; i < sxTS_num; i++)
    if (i%NTask == ThisTask) // distribute evenly on all tasks
      nexport++;

  // Allocate the memory and copy the data
  sxSourceCoord = (mesh_search_data *)mymalloc("sxSourceCoord", nexport * sizeof(mesh_search_data) );
  sxSourceIon = (struct sxSourceIon_struct *)mymalloc("sxSourceIon", nexport * sizeof(struct sxSourceIon_struct) );
  for(int s = 0; s < nexport ;s++)
    {
      int i = s*NTask+ThisTask; // calculate the correct test source index
      sxSourceCoord[s].Pos[0] = sxTS_positions[i][0];
      sxSourceCoord[s].Pos[1] = sxTS_positions[i][1];
      sxSourceCoord[s].Pos[2] = sxTS_positions[i][2];
#ifndef DO_NOT_RANDOMIZE_DOMAINCENTER
      // we have to correct the source coordinates according to the current domain displacement
      // DEBUG: this still does not solve all problems, because we need to first implement
      //        periodic photon radiation
      domain_displacePosition( sxSourceCoord[i].Pos, DISPLACE_POSITION_FORWARD );
#endif
      for(int f=0; f<SX_NFREQ; f++)
	  sxSourceIon[s].Rate[f] = sxTS_ionizationRates[i][f];
    }

  return nexport;
}

#elif (SX_SOURCES==4)

int sx_phot_sources_load(void)
{
  // Count number of sources for export
  int nexport = 0;
  for(int i = 0; i < NumPart; i++)
    if(P[i].Type == SX_SOURCES)
      nexport++;

  // Allocate the memory and copy the data
  sxSourceCoord = (mesh_search_data *)mymalloc("sxSourceCoord", nexport * sizeof(mesh_search_data) );
  sxSourceIon = (struct sxSourceIon_struct *)mymalloc("sxSourceIon", nexport * sizeof(struct sxSourceIon_struct) );
  int s = 0;
  for(int i = 0; i < NumPart; i++)
    {                                              
      if(P[i].Type == SX_SOURCES)
	{
	  sxSourceCoord[s].Pos[0] = P[i].Pos[0];                                                            
	  sxSourceCoord[s].Pos[1] = P[i].Pos[1];
	  sxSourceCoord[s].Pos[2] = P[i].Pos[2];                                                                                
	  sx_chem_compute_ionisation_rates( sxSourceIon[s].Rate, i );
	  s++;                                                                   
	}                                                                                                     
    }                

  return nexport;
}

#elif (SX_SOURCES==5) 

#if defined(SX_SKIP_RADIUS)
/*! \brief Finds the particles within a sphere
 *
 *  remember to deallocate indices after this function is called
 */
void find_particles_within_a_sphere(double center[3], double radius, int *local_n, int *total_n, int **indices)
{
  // only search on local threads
  // kind of naive use of the function ngb_treefind_variable_threads(), but it seems to work
  int n, j;
  int thread_id = get_thread_num();


  Thread[thread_id].R2list  = (double *)mymalloc("R2list", NumPart * sizeof(double));
  Thread[thread_id].Ngblist = (int *)mymalloc("Ngblist", NumPart * sizeof(int));

  int nfound = ngb_treefind_variable_threads(center, radius, -1, MODE_LOCAL_PARTICLES, thread_id, 1, NULL);

  int index[nfound];
  double r2;
  *local_n = 0;

  for(int n = 0; n < nfound; n++)
    {
      j = Thread[thread_id].Ngblist[n];

      // r2 = Thread[thread_id].R2list[n];
      r2 = pow(center[0] - P[j].Pos[0], 2) + pow(center[1] - P[j].Pos[1], 2) + pow(center[2] - P[j].Pos[2], 2);
      if((r2 < radius * radius) && (P[j].ID != 0) && (P[j].Mass > 0) && (P[j].Type == 0))
        {
          index[(*local_n)++] = j;
        }
    }

  myfree(Thread[thread_id].Ngblist);
  myfree(Thread[thread_id].R2list);
  *indices = (int *)mymalloc("indices", *local_n * sizeof(int));

  memcpy(*indices, &index, *local_n * sizeof(int));
  MPI_Allreduce(local_n, total_n, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}
void skip_sink_radius(int index)
{
  int local_n, ntot, *indices;
  find_particles_within_a_sphere(SinkP[index].Pos, (double)SinkFormationRadius, &local_n, &ntot, &indices);
  if (local_n>1)
    {
      for (int c=0; c<local_n; c++)
	{
	  sxCell[indices[c]].skip = 1;
	  //printf("SX: (%d) i %d local_n %d ntot %d c %d \n",ThisTask,i,local_n,ntot,indices[c]);
	}
    }
  myfree(indices);
}
#endif

int sx_phot_sources_load(void)
{
  int nexport=0;

  // Count numbger of sources for export
  for(int i = 0, s=0; i < NSinksAllTasks; i++)
      if(SinkP[i].HomeTask == ThisTask)
	nexport++;

#if defined(SX_SKIP_RADIUS) // "un-skip" all cells
  for(int i = 0; i < NumGas; i++)
      sxCell[i].skip = 0;
  int maxIndex = 0;
#if SX_SKIP_RADIUS>0
  double solMass = SX_SKIP_RADIUS*SOLAR_MASS/All.UnitMass_in_g;  // solar mass in code units
#else
  double maxMass = 0;
#endif
#endif

  // Allocate the memory and copy the data
  sxSourceCoord = (mesh_search_data *)mymalloc("sxSourceCoord", nexport * sizeof(mesh_search_data) );
  sxSourceIon = (struct sxSourceIon_struct *)mymalloc("sxSourceIon", nexport * sizeof(struct sxSourceIon_struct) );
  for(int i = 0, s=0; i < NSinksAllTasks; i++)
    { 
      if(SinkP[i].HomeTask == ThisTask)
	{
	  sxSourceCoord[s].Pos[0] = SinkP[i].Pos[0];                                                            
	  sxSourceCoord[s].Pos[1] = SinkP[i].Pos[1];
	  sxSourceCoord[s].Pos[2] = SinkP[i].Pos[2];                                                                                
#ifdef POPIII_SNE
          sx_chem_compute_popIII_ionisation_rates(sxSourceIon[s].Rate, i );
#else
	  sx_chem_compute_ionisation_rates( sxSourceIon[s].Rate, i );
#endif
	  s++;                                                                   
	}
#if defined(SX_SKIP_RADIUS)  // find the most massive sink
#if SX_SKIP_RADIUS>0
      // skipping radius of all sinks with smaller than 10 M_sol
      if (SinkP[i].Mass>solMass)  
	{
	  skip_sink_radius(i);
	}
#else
      // finding an index of the most massive sink particle
      if (SinkP[i].Mass>maxMass)
	{
	  maxIndex = i;
	  maxMass = SinkP[i].Mass;
	}
#endif
#endif      
    }

#if defined(SX_SKIP_RADIUS)&&(SX_SKIP_RADIUS==0)  
  // skip all cells around the sources
  if (maxMass>0)
    {
      skip_sink_radius(maxIndex);
    }
#endif

  return nexport;
}
#endif

/*
 * \brief Calculate an ionization rate in every gas cell with a source
 */
void sx_phot_sources_start_run(void)
{
  int i, s, p, t, f;
  int nimport=0, nexport=0, ngrp;
  int dummy[NTask]; // dummy source counter
  double conversionFactor = 1.0 / All.UnitPhotons_per_s;

  // Load sources
  nexport = sx_phot_sources_load();

#ifdef SX_DISPLAY_SOURCES
  // Print out a list of loaded sources
  mpi_printf("SX: Source positions\n");
  MPI_Barrier(MPI_COMM_WORLD);
  for(int i=0; i<nexport; i++)
    {
      char buff[50];
      for(int f=0; f<SX_NFREQ; f++)
	{
	  sprintf(&buff[f*9],"%.02e ",sxSourceIon[i].Rate[f]);
	}
      printf("SX: (%d) coord [ %.03e %.03e %.03e ] rates [ %s]\n",ThisTask,
	     sxSourceCoord[i].Pos[0],sxSourceCoord[i].Pos[1],sxSourceCoord[i].Pos[2],buff);
    }  
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // Search for the nearest gas particles
  find_nearest_meshpoint_global(sxSourceCoord, nexport, 0, 0);

  // Count items that we want to send
  for(t = 0; t < NTask; t++){
    Send_count[t] = 0;
  }
  for(s = 0; s<nexport; s++)
    {
      //p = sxSourceCoord[s].u.Index;
      t = sxSourceCoord[s].Task;
      Send_count[t]++;
      /*
      if (t==ThisTask){
	printf("SX: (%d) Task %d Index %d ID %d SrcPos [ %.03e %.03e %.03e ] PartPos [ %.03e %.03e %.03e ] \n", ThisTask,
	       t, p, P[p].ID, sxSourceCoord[s].Pos[0], sxSourceCoord[s].Pos[1], sxSourceCoord[s].Pos[2],
	       P[p].Pos[0], P[p].Pos[1], P[p].Pos[2]);
      }else{
	printf("SX: (%d) Task %d Index %d SrcPos [ %.03e %.03e %.03e ] \n", ThisTask,
	       t, p, sxSourceCoord[s].Pos[0], sxSourceCoord[s].Pos[1], sxSourceCoord[s].Pos[2] );
	}
      */
      
    }

  // Count items that we will receive
  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  // Calculate exchange variables
  for(t = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; t < NTask; t++)
    {
      nimport += Recv_count[t];
      if(t > 0)
        {
          Send_offset[t] = Send_offset[t - 1] + Send_count[t - 1];
          Recv_offset[t] = Recv_offset[t - 1] + Recv_count[t - 1];
        }
    }

  /*
  printf("SX: (%d) nimport %d nexport %d Send_count [ ",ThisTask,nimport,nexport);
  for(t = 0; t < NTask; t++){
    printf("%d ",Send_count[t]);
  }
  printf("] Recv_count [ ");
  for(t = 0; t < NTask; t++){
    printf("%d ",Recv_count[t]);
  }
  printf("] Send_offset [ ");
  for(t = 0; t < NTask; t++){
    printf("%d ",Send_offset[t]);
  }
  printf("] Recv_offset [ ");
  for(t = 0; t < NTask; t++){
    printf("%d ",Recv_offset[t]);
  }
  printf("] \n");
  */

  // Prepare array with sources for the sending
  struct sxSource_struct SourceSend[nexport];
  memset(dummy, 0, NTask * sizeof(int));
  for(s = 0; s<nexport; s++)
    {
      t = sxSourceCoord[s].Task;          // task of the gas particle
      p = sxSourceCoord[s].u.Index;       // index of the gas particle on the task "t"
      i = Send_offset[t]+dummy[t]; // index in the SourceSend
      //printf("SX (%d) %d = %d + %d \n",ThisTask, i, Send_offset[t], dummy[t]);
      SourceSend[i].index = p;
      //printf("SX: (%d) source %d index %d ionRate [ ",ThisTask,s,p);
      for(f=0; f<SX_NFREQ; f++)
	{
	  SourceSend[i].ionRate[f] = sxSourceIon[s].Rate[f] * conversionFactor;
	  //printf("%.03e ", ionRates[s][f]);
	}
      //printf("] \n");
      dummy[t]++;
    }
  /*
  for(s = 0; s<nexport; s++)
    {
      printf("SX: (%d) SEND: source %d index %d ionRate [ ",
	     ThisTask, s, sxSourceSend[s].index);
      for(f=0; f<SX_NFREQ; f++)
	{
	  printf("%.03e ", sxSourceSend[s].ionRate[f]);
	}      
      printf("] \n");
    }
  */

  // Free the search data
  myfree(sxSourceIon);
  myfree(sxSourceCoord);
  
  // Exchange sources
  struct sxSource_struct SourceRecv[nimport];
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              // get the particles
              MPI_Sendrecv(&SourceSend[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct sxSource_struct),
			   MPI_BYTE, recvTask, TAG_DENS_A,
                           &SourceRecv[Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(struct sxSource_struct),
			   MPI_BYTE, recvTask, TAG_DENS_A,
			   MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  // Count number of source cells and set indexes to ionisation rates
  int nSources = 0;
  for(i = 0; i < NumGas; i++)
    {
      sxCell[i].Source = -1;
    }
  for (i = 0; i < nimport; i++)
    {
      p = SourceRecv[i].index;
      if (sxCell[p].Source<0)
	{
          sxCell[p].Source = nSources++;
        }
    }
  
  // Create source array and fill in the ionisation rates
  sxSources = (struct sxSource_struct *)mymalloc("sxSources", nSources * sizeof(struct sxSource_struct));
  for (i = 0; i < nSources; i++)
    memset(sxSources[i].ionRate, 0, SX_NFREQ * sizeof(double));
  for (i = 0; i < nimport; i++)
    {
      p = SourceRecv[i].index;
      s = sxCell[p].Source;
      sxSources[s].index = p;
      for(f=0; f<SX_NFREQ; f++)
	  sxSources[s].ionRate[f] += SourceRecv[i].ionRate[f];
    }

#ifdef SX_DISPLAY_SOURCES
  // Print out a list of source cells
  mpi_printf("SX: Source cells\n");
  MPI_Barrier(MPI_COMM_WORLD);
  for(int i=0; i<NumGas; i++)
    {
      s = sxCell[i].Source;
      if (s>=0)
	{
	  char buff[60];
	  for(int f=0; f<SX_NFREQ; f++)
	      sprintf(&buff[f*9],"%.02e ",sxSources[s].ionRate[f] / conversionFactor);
	  printf("SX: (%d) coord [ %.03e %.03e %.03e ] rates [ %s]\n",ThisTask,
		  P[i].Pos[0],P[i].Pos[1],P[i].Pos[2],buff);
	}
    }  
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // update number of sources
  sxRunData->nSourcesLocal = nSources;
  sxRunData->nSourcesExport = nexport;
  MPI_Reduce(&nSources, &sxRunData->nSources, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  // Display information about the sources
  sx_debug_sources();

}

/*
 * \brief free the source memory
 */
void sx_phot_sources_end_run(void)
{
  
  // free the source list
  myfree(sxSources);  

}
