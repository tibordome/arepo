/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/calculate_quantities_in_postprocess/calculate_quantities.c
 * \date        MM/YYYY
 * \author      Rahul Kannan
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../domain.h"
#include "../proto.h"
#include "../voronoi.h"

#include "calculate_quantities.h"

#ifdef CALCULATE_QUANTITIES_IN_POSTPROCESS

void calculate_quantities_in_postprocess(void)
{
  mpi_printf("CALCULATE_QUANTITIES_IN_POSTPROCESS: Calculating Vorticities\n");
  calculate_cell_centered_vorticity();
}

void calculate_cell_centered_vorticity(void)
{
  int i, sph_idx;

  for(i = 0; i < NumGas; i++)
    {
      sph_idx = i;
      //      SphP[i].Vorticity[0] = SphP[sph_idx].Grad.dvel[2][1] - SphP[sph_idx].Grad.dvel[1][2];
      // SphP[i].Vorticity[1] = SphP[sph_idx].Grad.dvel[0][2] - SphP[sph_idx].Grad.dvel[2][0];
      // SphP[i].Vorticity[2] = SphP[sph_idx].Grad.dvel[1][0] - SphP[sph_idx].Grad.dvel[0][1];
      //      mpi_printf("%le\t%le\t%le\t", P[i].Vel[0], P[i].Vel[1], P[i].Vel[2]) ;
      // mpi_printf("%le\t%le\t%le\n", SphP[i].Vorticity[0], SphP[i].Vorticity[1], SphP[i].Vorticity[2]) ;
      mpi_printf("%le || \t", SphP[i].Utherm);
      mpi_printf("%le\t%le\t%le\n", SphP[i].Grad.dutherm[0], SphP[i].Grad.dutherm[1], SphP[i].Grad.dutherm[2]);
    }

  return;
}

void save_additional_quantities_data(int RestartSnapNum)
{
  char buf[MAXLEN_PATH];

  file_path_sprintf(buf, "%s_%03d_addquant_%03d.dat", All.SnapshotFileBase, RestartSnapNum, ThisTask);

  FILE *fp = fopen(buf, "w");
  fwrite(&NumGas, 4, 1, fp);
  int i;
  for(i = 0; i < NumGas; i++)
    fwrite(&P[i].ID, sizeof(MyIDType), 1, fp);

  for(i = 0; i < NumGas; i++)
    fwrite(&SphP[i].Vorticity[0], sizeof(MyFloat), 1, fp);

  for(i = 0; i < NumGas; i++)
    fwrite(&SphP[i].Vorticity[1], sizeof(MyFloat), 1, fp);

  for(i = 0; i < NumGas; i++)
    fwrite(&SphP[i].Vorticity[2], sizeof(MyFloat), 1, fp);

  for(i = 0; i < NumGas; i++)
    fwrite(&SphP[i].Grad.dutherm[0], sizeof(MyFloat), 1, fp);

  for(i = 0; i < NumGas; i++)
    fwrite(&SphP[i].Grad.dutherm[1], sizeof(MyFloat), 1, fp);

  for(i = 0; i < NumGas; i++)
    fwrite(&SphP[i].Grad.dutherm[2], sizeof(MyFloat), 1, fp);

  fclose(fp);
}

#endif
