/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \date        22/10/2018
 * \author      Ondrej Jaura
 * \brief       Debugging functionality of the SPRAI
 * \details     
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

static int nLoc = 2;
//static double loc[][3] = { {0.5,0.5,0.5}, {0.55,0.5,0.5}, {0.6,0.5,0.5}, {0.65,0.5,0.5},
//                         {0.7,0.5,0.5}, {0.75,0.5,0.5}, {0.8,0.5,0.5}, {0.85,0.5,0.5} };
static double locCoord[][3] = { {0.5,0.5,0.5}, {0.8,0.5,0.5} };
mesh_search_data * locSearch;

// Structure that keeps record of different imers
struct sxTimer_struct
{
  clock_t TimeStart[SX_NTIMERS];
  double TimeSum[SX_NTIMERS];
} * sxTimer;

//======================================================================================================//
//                                     Particle location finder                                         //
//======================================================================================================//

/*!
 * \brief find debug particles closest to the given locations
 */
void sx_debug_find_locations(void)
{
  int i, t, p;
  locSearch = (mesh_search_data *)mymalloc("locSearch", nLoc * sizeof(mesh_search_data) );
  for(i = 0; i < nLoc; i++)
    {
      locSearch[i].Pos[0] = locCoord[i][0];
      locSearch[i].Pos[1] = locCoord[i][1];
      locSearch[i].Pos[2] = locCoord[i][2];
    }
  find_nearest_meshpoint_global(locSearch, nLoc, 0, 0);
  mpi_fprintf(FdSimplex,"SX: Found debug particles: \n");
  for(i = 0; i < nLoc; i++)
    {
      t = locSearch[i].Task;
      if (t!=ThisTask) continue;
      p = locSearch[i].u.Index;
      mpi_fprintf(FdSimplex,"SX:    (%d:%d:%d) xyz [ %.03e %.03e %.03e ] \n",
             t,p,P[p].ID, locSearch[i].Pos[0], locSearch[i].Pos[1], locSearch[i].Pos[2]);
    }
}

/*!
 * \brief This function tells you if the particle with the given ID is one of the given locations
 */
_Bool sx_debug_is_location( int pid )
{
  int i,p,t;
  for(i = 0; i < nLoc; i++)
    {
      t = locSearch[i].Task;
      if (t!=ThisTask) continue;
      p = locSearch[i].u.Index;
      if (P[p].ID==pid) return 1;
    }
  return 0;
}

/*!
 * \brief initialize variables at the beginning of the RT run
 */
/*
void sx_debug_start_run(void)
{
  int i;

  // DEBUG: we want to find index of a particle with ID==SXPID
  SXPID = 32769;  // 19382;
  SXID  = -1;

  // look for a tracked particle with SXPID
  for(i = 0; i < NumGas; i++)
    {
      if(P[i].ID == SXPID)
        SXID = i;
    }

  // DEBUG: check if tracked particle was found at some processor
  int GlobalSXID;
  MPI_Reduce(&SXID, &GlobalSXID, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  if(GlobalSXID < 0)
    mpi_printf("SX: tracking particle with ID %d was not found \n", SXPID);
  if(SXID >= 0)
    printf("SX: tracking particle with ID %d has index %d on task %d \n", SXPID, SXID, ThisTask);

  sx_debug_find_locations();

}
*/

/*!
 * \brief finalize variables at the end of the RT run
 */
/*
void sx_debug_end_run(void)
{
  myfree(locSearch);
}
*/

//======================================================================================================//
//                                     Site debugging functions                                         //
//======================================================================================================//

void sx_debug_sources(void)
{
  int t, f;
  int dummy[NTask]; // dummy source counter

  // print source statistics
  mpi_fprintf(FdSimplex,"Number of source gas particles is %d \n", sxRunData->nSources);

#ifdef SX_DISPLAY_STATS
  // Number of source particles on each task
  MPI_Gather(&sxRunData->nSourcesExport, 1, MPI_INT, &dummy, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (ThisTask==0)
    {
      mpi_fprintf(FdSimplex,"Source locations per task:     [ ");
      for(t = 0; t < NTask; t++){
	mpi_fprintf(FdSimplex,"%d ",dummy[t]);
      }
      mpi_fprintf(FdSimplex,"] \n");
    }
  // Number of source sites (gas) on each task
  MPI_Gather(&sxRunData->nSourcesLocal, 1, MPI_INT, &dummy, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (ThisTask==0)
    {
      mpi_fprintf(FdSimplex,"Source gas particles per task: [ ");
      for(t = 0; t < NTask; t++){
	mpi_fprintf(FdSimplex,"%d ",dummy[t]);
      }
      mpi_fprintf(FdSimplex,"] \n");
    }
  // Number of memory used to store sources on each task
  int mSources = sxRunData->nSourcesLocal * sizeof(struct sxSource_struct);
  MPI_Gather(&mSources, 1, MPI_INT, &dummy, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (ThisTask==0)
    {
      mpi_fprintf(FdSimplex,"Source memory per task: [ ");
      for(t = 0; t < NTask; t++){
	mpi_fprintf(FdSimplex,"%d ",dummy[t]);
      }
      mpi_fprintf(FdSimplex,"] byte \n");
    }
#endif

#ifdef SX_DISPLAY_SOURCES
  // List sources on each task (to the output.log file)
  for(t = 0; t < sxRunData->nSourcesLocal; t++)
    {
      printf("SX: (%d) source %d ion [ ",ThisTask,t);	 
      for(f = 0; f < SX_NFREQ; f++)
	{
	  printf("%.03e ",sxSources[t].ionRate[f]*All.UnitPhotons_per_s);
	}
      f = sxSources[t].index;
      printf("] pos [ %.03e %.03e %.03e ] \n",P[f].Pos[0],P[f].Pos[1],P[f].Pos[2]);
    }
#endif
}

/*!
 * \brief show information about SXPID particle
 */
void sx_debug_particle(int index, int ID)
{
  int f;
  double radius = 0.0;
  if(P[index].ID == ID)
    {
      printf("SX: (%d:%d:%d) Pos [%.03e %.03e %.03e] Energy %.03e Density %.03e Pressure %.03e Utherm %.03e \n",
	     ThisTask,index,P[index].ID,
	     P[index].Pos[0], P[index].Pos[0], P[index].Pos[0],
	     SphP[index].Energy, SphP[index].Density, SphP[index].Pressure, SphP[index].Utherm);
      printf("SX: (%d:%d:%d) Ionisation  [ ", ThisTask,index,P[index].ID);
      for(f=0; f<SX_NFREQ; f++)
	{
	  printf("%.3e ", sxSources[sxCell[index].Source].ionRate[f]);
	}
      printf("] \n");
      /*
      printf("SX:   Fluxes      [ ");
      for(f=0; f<SX_NFREQ; f++)
	{
	  printf("%.3e ", pow(10,SphP[index].sxPhotonFlux[f]) );
	}
      printf("] \n");
      */
#if (SX_CHEMISTRY==3)
      printf("SX: (%d:%d:%d) Rates       [ ", ThisTask,index,P[index].ID);
      for(f=0; f<SX_NRATES; f++)
	{
	  printf("%.3e ", SphP[index].sxPhotonRates[f]);
	}
      printf("] \n");
      printf("SX: (%d:%d:%d) TracedAbund [ ", ThisTask, index, P[index].ID);
      for(f = 0; f < SGCHEM_NUM_ADVECTED_SPECIES; f++)
        {
          printf("%.3e ", SphP[index].TracAbund[f]);
        }
      printf("] \n");
#elif(SX_CHEMISTRY == 4)
      printf("SX: (%d:%d:%d) MassFractions [ ", ThisTask, index, P[index].ID);
      for(f = 0; f < SX_NMASSFRACT; f++)
        {
          printf("%.3e ", SphP[index].MassFract[f]);
        }
      printf("] \n");
#endif
    }
}
void sx_debug_site(struct sxSite *site, int ID)
{
  int f;
  if(P[site->index].ID == ID)
    {
      printf("SX: >--------------  Site: Task %d Index %d ID %d source %d ------------------ \n", ThisTask, site->index,
             P[site->index].ID, site->isSource);
      printf("SX:   PhotonsOrig [ ");
      for(f = 0; f < SX_NFREQ; f++)
        {
          printf("%.3e ", site->nPhotonsOriginal[f]);
        }
      printf("] \n");
      printf("SX:   Photons     [ ");
      for(f = 0; f < SX_NFREQ; f++)
        {
          printf("%.3e ", site->nPhotons[f]);
        }
      printf("] \n");
      /*
      printf("SX:   Fluxes      [ ");
      for(f=0; f<SX_NFREQ; f++)
        {
          printf("%.3e ", pow(10,site->fluxes[f]) );
        }
      printf("] \n");
      */
#if(SX_CHEMISTRY == 3)
      printf("SX:   Rates       [ ");
      for(f = 0; f < SX_NRATES; f++)
        {
          printf("%.3e ", site->rates[f]);
        }
      printf("] \n");
      printf("SX:   xH %.3e xH2 %.3e xHp %.3e volume %.3e numdens %.3e N %.3e density %.3e \n", site->xH, site->xH2, site->xHp,
             site->volume, site->numdens, site->volume * site->numdens, site->density);
#elif(SX_CHEMISTRY == 4)
      printf("SX:   nFract [ ");
      for(f = 0; f < SX_NMASSFRACT; f++)
        {
          printf("%.3e ", site->nFract[f]);
        }
      printf("] \n");
      printf("SX:   nNucleons %.3e density %.3e volume %.3e numdens %.3e numdr %.3e temp %.3e heating %.3e cooling %.3e \n",
             site->nNucleons, site->density, site->volume, site->numdens, site->numdr, site->temp, site->heating, site->cooling);
#endif
      printf("SX: >-------------endsite-------------< \n");
    }
}

//======================================================================================================//
//                                     Debuging statistics                                              //
//======================================================================================================//

#ifdef SX_DISPLAY_MEMORY
/*!
 * \brief Display memory statistics
 */
void sx_debug_display_memory(void)
{
  // Memory allocation of the ABPP and QPP arrays
  float maxAllocABPP, locAllocABPP = 100.0*(float)sxMaxNumABPP/(float)sxNumSPD;
  float maxAllocQPP, locAllocQPP = 100.0*(float)sxMaxNumQPP/(float)sxNumTotEDPD;
  MPI_Reduce(&locAllocABPP, &maxAllocABPP, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&locAllocQPP, &maxAllocQPP, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
  mpi_fprintf(FdSimplex,"Maximum allocation: ABPP %.2f %% QPP %.2f %% \n", maxAllocABPP, maxAllocQPP);

  int maxNumABPP, maxNumQPP;
  MPI_Reduce(&sxMaxNumABPP, &maxNumABPP, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sxMaxNumQPP, &maxNumQPP, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  mpi_fprintf(FdSimplex,"Maximum memory: ABPP %.04e b QPP %.04e b \n", 
	      (double)(maxNumABPP*sizeof(struct sxABPP_struct)*2),
	      (double)(maxNumQPP*sizeof(struct sxQPP_struct)) );

  int maxReallocABPP, maxReallocQPP;
  MPI_Reduce(&sxNumReallocABPP, &maxReallocABPP, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sxNumReallocQPP, &maxReallocQPP, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  mpi_fprintf(FdSimplex,"Maximum re-allocation: ABPP %d QPP %d \n", maxReallocABPP, maxReallocQPP);
}
#endif

#ifdef SX_DISPLAY_STATS
/*!
 * \brief Display detail information about current RT step
 */
void sx_debug_display_stats(void)
{
  int t, *nBalTr=NULL, *nDifTr=NULL, *nRecPack=NULL, *nSendPack=NULL;
  double *nInitPh=NULL, *nEmitPh=NULL, *nAbsPh=NULL, *nLostPh=NULL, *nFinPh=NULL;

  nInitPh = (double *)mymalloc("nInitPh", SX_NFREQ * sizeof(double));
  nEmitPh = (double *)mymalloc("nEmitPh", SX_NFREQ * sizeof(double));
  nAbsPh  = (double *)mymalloc("nAbsPh", SX_NFREQ * sizeof(double));
  nLostPh = (double *)mymalloc("nLostPh", SX_NFREQ * sizeof(double));
  nFinPh = (double *)mymalloc("nPhot", SX_NFREQ * sizeof(double));

  if(ThisTask == 0)
    {
      nBalTr    = (int *)mymalloc("nBalTr", NTask * sizeof(int));
      nDifTr    = (int *)mymalloc("nDifTr", NTask * sizeof(int));
      nRecPack  = (int *)mymalloc("nRecPack", NTask * sizeof(int));
      nSendPack = (int *)mymalloc("nSendPack", NTask * sizeof(int));
    }

  MPI_Gather(&sxRunData->nBalTransp, 1, MPI_INT, nBalTr, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&sxRunData->nDifTransp, 1, MPI_INT, nDifTr, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&sxRunData->nRecPackets, 1, MPI_INT, nRecPack, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&sxRunData->nSendPackets, 1, MPI_INT, nSendPack, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Allreduce(&sxRunData->nInitPhotons, nInitPh, SX_NFREQ, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&sxRunData->nEmitPhotons, nEmitPh, SX_NFREQ, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&sxRunData->nAbsPhotons, nAbsPh, SX_NFREQ, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&sxRunData->nLostPhotons, nLostPh, SX_NFREQ, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&sxRunData->nFinPhotons,  nFinPh,  SX_NFREQ, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      mpi_fprintf(FdSimplex,"DiffuseTransfer   [ ");
      for( t=0; t<NTask; t++ )
	{
	  mpi_fprintf(FdSimplex,"%d ",nDifTr[t]);
	}
      mpi_fprintf(FdSimplex,"]\n");
      mpi_fprintf(FdSimplex,"BallisticTransfer [ ");
      for( t=0; t<NTask; t++ )
	{
	  mpi_fprintf(FdSimplex,"%d ",nBalTr[t]);
	}
      mpi_fprintf(FdSimplex,"]\n");
      mpi_fprintf(FdSimplex,"SentPackets       [ ");
      for( t=0; t<NTask; t++ )
	{
	  mpi_fprintf(FdSimplex,"%d ",nSendPack[t]);
	}
      mpi_fprintf(FdSimplex,"]\n");
      mpi_fprintf(FdSimplex,"ReceivedPackets   [ ");
      for( t=0; t<NTask; t++ )
	{
	  mpi_fprintf(FdSimplex,"%d ",nRecPack[t]);
	}
      mpi_fprintf(FdSimplex,"] \n");
      mpi_fprintf(FdSimplex,"InitialPhotons    [ ");
      for( t=0; t<SX_NFREQ; t++ )
	{
	  mpi_fprintf(FdSimplex,"%.4e ",nInitPh[t]);
	}
      mpi_fprintf(FdSimplex,"] \n");
      mpi_fprintf(FdSimplex,"EmitedPhotons     [ ");
      for( t=0; t<SX_NFREQ; t++ )
	{
	  mpi_fprintf(FdSimplex,"%.4e ",nEmitPh[t]);
	}
      mpi_fprintf(FdSimplex,"] \n");
      mpi_fprintf(FdSimplex,"AbsorbedPhotons   [ ");
      for( t=0; t<SX_NFREQ; t++ )
	{
	  mpi_fprintf(FdSimplex,"%.4e ",nAbsPh[t]);
	}
      mpi_fprintf(FdSimplex,"] \n");
      mpi_fprintf(FdSimplex,"LostPhotons       [ ");
      for( t=0; t<SX_NFREQ; t++ )
	{
	  mpi_fprintf(FdSimplex,"%.4e ",nLostPh[t]);
	}
      mpi_fprintf(FdSimplex,"] \n");
      mpi_fprintf(FdSimplex,"RemainingPhotons  [ ");
      for( t=0; t<SX_NFREQ; t++ )
	{
	  mpi_fprintf(FdSimplex,"%.4e ",nFinPh[t]);
	}
      mpi_fprintf(FdSimplex,"] \n");
      
      myfree(nSendPack);
      myfree(nRecPack);
      myfree(nDifTr);
      myfree(nBalTr);
    }

  myfree(nFinPh);
  myfree(nLostPh);
  myfree(nAbsPh);
  myfree(nEmitPh);
  myfree(nInitPh);
}
#endif

//======================================================================================================//
//                                      Timers                                                          //
//======================================================================================================//

void sx_debug_timer_start( int tid )
{
  sxTimer->TimeStart[tid] = clock();
}

double sx_debug_timer_end( int tid )
{
  double dt = (double)(clock()-sxTimer->TimeStart[tid])/CLOCKS_PER_SEC;
  sxTimer->TimeSum[tid] += dt;
  return dt;
}

double sx_debug_timer_switch( int tid_old, int tid_new )
{
  sx_debug_timer_start( tid_new );
  return sx_debug_timer_end( tid_old );
}

//======================================================================================================//
//                                     Runtime functions                                                //
//======================================================================================================//

void sx_debug_start_run(void)
{
  int i;

  // allocate all timers
  sxTimer = (struct sxTimer_struct *)mymalloc("sxTimer", sizeof(struct sxTimer_struct) );
  for (i=0;i<SX_NTIMERS;i++)
    {
      sxTimer->TimeSum[i] = 0.;
    }
}

void sx_debug_end_run(void)
{
#ifdef SX_DISPLAY_TIMERS
  int i;
  double times[SX_NTIMERS];

  MPI_Reduce(&sxTimer->TimeSum, &times, SX_NTIMERS, MPI_DOUBLE, MPI_SUM, 0,  MPI_COMM_WORLD);
  if (ThisTask==0)
    {
      mpi_fprintf(FdSimplex, "Timers s [ ");
      for (i=0;i<SX_NTIMERS;i++)
	{
	  mpi_fprintf(FdSimplex,"%.03e ",times[i]/NTask);
	}
      mpi_fprintf(FdSimplex, "] \n");
      mpi_fprintf(FdSimplex, "Timers %% [ ");
      for (i=0;i<SX_NTIMERS;i++)
	{
	  mpi_fprintf(FdSimplex,"%.1f ",100*times[i]/times[TRUN]);
	}
      mpi_fprintf(FdSimplex, "] \n");
    }
#endif

#ifdef SX_DISPLAY_MEMORY
  // print out memory statistics
  sx_debug_display_memory();
#endif

  // free the memory
  myfree(sxTimer);
}
