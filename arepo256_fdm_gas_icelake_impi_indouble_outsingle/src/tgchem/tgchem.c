/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/tgchem/tgchem.c
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       Primordial chemistry and cooling network
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"

void tgchem(void)
{
  CPU_Step[CPU_MISC] += measure_time();

  TGCD.DebugFlag = -1;
  // TGCD.DebugFlag = 1225861435;

  if(TGCD.ChemMode >= 0)
    {
      double t0 = second();

      mpi_printf("TGCHEM: Doing chemistry step...\n");

      int i, idx;

      step_data step;

      tgchem_step_data(&step);

      for(i = 0; i < TimeBinsHydro.NActiveParticles; i++)
        {
          idx = TimeBinsHydro.ActiveParticleList[i];

          if(idx < 0)
            continue;

          // continue;

          cell_data cell;

          tgchem_cell_data(&step, &cell, idx);

          var_data var;

          tgchem_var_data(&cell, &var);

          int eq_flag_h2   = cell.eq_flag[1];
          int eq_flag_hii  = cell.eq_flag[2];
          double nh_bef    = cell.nh;
          double temp_bef  = var.temp;
          double abhii_bef = var.abhii;

          // if(cell.cell_id == TGCD.DebugFlag)
          //{
          // printf("nh = %g, temp = %g, abhii = %g, id = %d\n", nh_bef, temp_bef, abhii_bef, cell.cell_id);
          // terminate("");
          //}
          // if(nh_bef > 3e14 && abhii > 1e-4)

          tgchem_step(&step, &cell);

          tgchem_var_data(&cell, &var);

          double nh_af    = cell.nh;
          double temp_af  = var.temp;
          double abhii_af = var.abhii;

          // if(nh_bef > 3e14 && abhii_bef < 1e-4 && abhii_af > 1e-4)
          //{
          // if(cell.cell_id == TGCD.DebugFlag)
          // printf("eq_flag_h2 = %d, eq_flag_hii = %d, nh = %g, temp_bef = %g, temp_af = %g, abhii_bef = %g, abhii_af = %g, id =
          // %d\n", eq_flag_h2, eq_flag_hii, nh_bef, temp_bef, temp_af, abhii_bef, abhii_af, cell.cell_id); terminate("");
          //}

          // if(cell.cell_id == TGCD.DebugFlag)
          // if(nh_bef > 8e14 && log10(temp_af) < 3.6)
          // if(temp_af < 0.5 * temp_bef)
          // printf("MOEP! nh = %g, temp_bef = %g, temp_af = %g, id = %d\n", nh_bef, temp_bef, temp_af, cell.cell_id);
        }

      double t1 = second();

      double tdiff = timediff(t0, t1);

      long tot_num_rate_calls, tot_num_cells;

      MPI_Allreduce(&step.num_rate_calls, &tot_num_rate_calls, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&step.num_cells, &tot_num_cells, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

      long tot_num_neq, tot_num_eq_h2, tot_num_eq_hii;

      MPI_Allreduce(&step.num_neq, &tot_num_neq, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&step.num_eq_h2, &tot_num_eq_h2, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&step.num_eq_hii, &tot_num_eq_hii, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

      long tot_num_sub_steps = tot_num_neq + tot_num_eq_h2 + tot_num_eq_hii;

      double frac_cells_neq    = 0.;
      double frac_cells_eq_h2  = 0.;
      double frac_cells_eq_hii = 0.;

      if(tot_num_sub_steps)
        {
          frac_cells_neq    = (double)tot_num_neq / tot_num_sub_steps;
          frac_cells_eq_h2  = (double)tot_num_eq_h2 / tot_num_sub_steps;
          frac_cells_eq_hii = (double)tot_num_eq_hii / tot_num_sub_steps;
        }

      long avg_num_sub_steps  = 0;
      long avg_num_rate_calls = 0;

      if(tot_num_cells)
        {
          avg_num_sub_steps  = tot_num_sub_steps / tot_num_cells;
          avg_num_rate_calls = tot_num_rate_calls / tot_num_cells;
        }

      double tdiff_neq, tdiff_eq_h2, tdiff_eq_hii, tdiff_avg, tdiff_max;

      MPI_Allreduce(&step.dt_neq, &tdiff_neq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&step.dt_eq_h2, &tdiff_eq_h2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&step.dt_eq_hii, &tdiff_eq_hii, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&tdiff, &tdiff_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&tdiff, &tdiff_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      double frac_tdiff_neq    = tdiff_neq / tdiff_avg;
      double frac_tdiff_eq_h2  = tdiff_eq_h2 / tdiff_avg;
      double frac_tdiff_eq_hii = tdiff_eq_hii / tdiff_avg;

      double cells_per_s = tot_num_cells / tdiff_avg;

      tdiff_avg /= NTask;

      double imba = (tdiff_max - tdiff_avg) / tdiff_avg;

      mpi_printf(
          "TGCHEM: Done! Calls: %d|%d, Eq fractions: %g|%g|%g (%g|%g|%g), Imbalance: %g, Cells/task/second: %g, Took %g seconds\n",
          avg_num_sub_steps, avg_num_rate_calls, frac_cells_neq, frac_cells_eq_h2, frac_cells_eq_hii, frac_tdiff_neq, frac_tdiff_eq_h2,
          frac_tdiff_eq_hii, imba, cells_per_s, tdiff_max);

      // if(All.NumCurrentTiStep > 3)
      // endrun();

#ifdef HEALRAY
      if(HRD.SourceFlag >= 0 && !HRD.RayTimeFac)
        {
          produce_dump();

          endrun();
        }
#endif
    }

  CPU_Step[CPU_TGCHEM] += measure_time();
}
