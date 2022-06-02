/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/GFM/stellar_evolution_logs.c
 * \date        MM/YYYY
 * \author
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef GFM_STELLAR_EVOLUTION

void output_stellar_evolution_statistics(void)
{
  mpi_printf("GFM_STELLAR_EVOLUTION: Writing stellar evolution statistics.\n");

  int i, j;

  MyFloat metal_mass_gas[GFM_N_CHEM_ELEMENTS + 1], metal_mass_stars[GFM_N_CHEM_ELEMENTS + 1],
      metal_mass_total[GFM_N_CHEM_ELEMENTS + 1], global_metal_mass[GFM_N_CHEM_ELEMENTS + 1];
  MyFloat SNIaRate = 0.0, global_SNIaRate = 0.0;
  MyFloat SNIIRate = 0.0, global_SNIIRate = 0.0;

  /* metal data */
  for(i = 0; i < GFM_N_CHEM_ELEMENTS + 1; i++)
    {
      metal_mass_gas[i]   = 0.0;
      metal_mass_stars[i] = 0.0;
    }

  for(j = 0; j < NumPart; j++)
    {
      if(P[j].Type == 0) /* gas particle */
        {
          for(i = 0; i < GFM_N_CHEM_ELEMENTS; i++)
            metal_mass_gas[i] += SphP[j].MassMetals[i];

          metal_mass_gas[GFM_N_CHEM_ELEMENTS] += P[j].Mass;
        }
      if(P[j].Type == 4 && P[j].Mass > 0 && STP(j).BirthTime > 0) /* star particle */
        {
          for(i = 0; i < GFM_N_CHEM_ELEMENTS; i++)
            metal_mass_stars[i] += STP(j).MassMetals[i];

          metal_mass_stars[GFM_N_CHEM_ELEMENTS] += P[j].Mass;
        }
    }

  for(i = 0; i < GFM_N_CHEM_ELEMENTS; i++)
    metal_mass_total[i] = metal_mass_gas[i] + metal_mass_stars[i];

  metal_mass_total[GFM_N_CHEM_ELEMENTS] = metal_mass_gas[GFM_N_CHEM_ELEMENTS] + metal_mass_stars[GFM_N_CHEM_ELEMENTS];

  /* output total metal mass in gas particles */
  MPI_Reduce(metal_mass_gas, global_metal_mass, GFM_N_CHEM_ELEMENTS + 1, MPI_MYFLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fprintf(FdMetalsGas, "%e", All.Time);
      for(i = 0; i < GFM_N_CHEM_ELEMENTS + 1; i++)
        fprintf(FdMetalsGas, " %e", global_metal_mass[i]);
      fprintf(FdMetalsGas, "\n");
      myflush(FdMetalsGas);
    }

  /* output total metal mass in star particles */
  MPI_Reduce(metal_mass_stars, global_metal_mass, GFM_N_CHEM_ELEMENTS + 1, MPI_MYFLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fprintf(FdMetalsStars, "%e", All.Time);
      for(i = 0; i < GFM_N_CHEM_ELEMENTS + 1; i++)
        fprintf(FdMetalsStars, " %e", global_metal_mass[i]);
      fprintf(FdMetalsStars, "\n");
      myflush(FdMetalsStars);
    }

  /* output total metal mass in all particles */
  MPI_Reduce(metal_mass_total, global_metal_mass, GFM_N_CHEM_ELEMENTS + 1, MPI_MYFLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fprintf(FdMetalsTot, "%e", All.Time);
      for(i = 0; i < GFM_N_CHEM_ELEMENTS + 1; i++)
        fprintf(FdMetalsTot, " %e", global_metal_mass[i]);
      fprintf(FdMetalsTot, "\n");
      myflush(FdMetalsTot);
    }

  /* SN data */
  for(i = 0; i < N_star; i++)
    {
      SNIaRate += StarP[i].SNIaRate;
      SNIIRate += StarP[i].SNIIRate;
    }

  MPI_Reduce(&SNIaRate, &global_SNIaRate, 1, MPI_MYFLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&SNIIRate, &global_SNIIRate, 1, MPI_MYFLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fprintf(FdSN, "%e %e %e\n", All.Time, global_SNIaRate, global_SNIIRate);
      myflush(FdSN);
    }
}

#endif
