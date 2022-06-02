/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/dust_timestep.c
 * \date        MM/YYYY
 * \author      Ryan McKinnon
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
#include "../proto.h"

#ifdef DUST_LIVE

/* Update dust time bins of active dust particles using minimum hydro time bin
 * of nearby gas cells.
 */
void update_timebins_dust(void)
{
  start_dust();
  find_drag_cells(Ndust);
  drag_kernel();
  for(int i = 0; i < Ndust; i++)
    {
      const int p               = DustParticle[i].index;
      const int bin_old         = P[p].TimeBinHydro;
      const integertime ti_step = get_timestep_dust(p);

      /* Since we have a time bin directly, and not dt, we could just handle
       * synchronization manually, but we'll use the helper for its validity
       * checks and error handling. */
      const int bin_new = timebins_get_bin_and_do_validity_checks(ti_step, bin_old);

      timebin_move_particle(&TimeBinsDust, p, bin_old, bin_new);
      P[p].TimeBinHydro = bin_new;
    }
  end_dust();
}

integertime get_timestep_dust(int p)
{
  double dt;
  integertime ti_step;

  /* Start with timestep corresponding to minimum of nearby gas cells. */
  double dt_gas_min = (DTP(p).MinGasTimeBin ? (((integertime)1) << DTP(p).MinGasTimeBin) : 0) * All.Timebase_interval;
  dt_gas_min /= All.cf_hubble_a;

  dt = dt_gas_min;

  double rel_vel2 = 0.0;
  for(int j = 0; j < 3; j++)
    {
      double dvel = (DTP(p).LocalGasVelocity[j] - P[p].Vel[j]) / All.cf_atime;
      rel_vel2 += (dvel * dvel);
    }

  double csnd_eff = sqrt(pow(DTP(p).LocalSoundSpeed, 2) + rel_vel2);
  if(csnd_eff <= 0.0)
    csnd_eff = 1.0e-30;

  double dt_courant = All.CourantFac * (DTP(p).Hsml * All.cf_atime) / csnd_eff;

  if(dt_courant < dt)
    dt = dt_courant;

#if !defined(DL_DRAG_SEMI_IMPLICIT) && !defined(DL_NODRAG)
  /* If we require explicit drag timesteps, make sure to resolve the stopping timescale. */
  if((DTP(p).StoppingTime > 0.0) && (dt > DTP(p).StoppingTime * All.StoppingTimeFrac))
    {
      dt = DTP(p).StoppingTime * All.StoppingTimeFrac;
    }
#endif

#ifdef DL_GRAIN_BINS
#if defined(DL_GROWTH) || defined(DL_SPUTTERING) || defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION)
  /* If we are evolving a grain size distribution, choose a timestep limit
   * based on the maximum allowable fractional change in mass in the grain size
   * bins. */
  double dt_gsd = get_dust_dt_gsd(p);
#ifdef DL_SUBCYCLE
  /* If we are subcycling the grain size evolution updates, allow particle
   * timesteps to be greater than the grain size evolution timescale. */
  dt_gsd *= All.DustSubcycleFac;
#endif
  if(dt > dt_gsd)
    dt = dt_gsd;

  DTP(p).OldBinMassChgTau = DTP(p).BinMassChgTau;

  /* Reset in preparation for next timestep. */
  DTP(p).BinMassChgTau = MAX_REAL_NUMBER;

  /* Store current mass for use in timescale calculations taking place during
   * grain size evolution. */
  DTP(p).OrigMass = P[p].Mass;
#endif
#endif

  /* Convert physical timestep to dloga if needed. */
  dt *= All.cf_hubble_a;

  if(dt >= All.MaxSizeTimestep)
    dt = All.MaxSizeTimestep;

#ifdef PMGRID
  if(dt >= All.DtDisplacement)
    dt = All.DtDisplacement;
#endif

  ti_step = (integertime)(dt / All.Timebase_interval);
  validate_timestep(dt, ti_step, p, TSTP_DUST_LIVE);

  return ti_step;
}

#ifdef DL_GRAIN_BINS
#if defined(DL_GROWTH) || defined(DL_SPUTTERING) || defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION)
/* Return dust particle's allowable timestep for grain size evolution in
 * internal time units. */
double get_dust_dt_gsd(int p) { return All.MaxBinFracMassChg * DTP(p).BinMassChgTau; }
#endif
#endif

#endif
