/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/criterion_derefinement.c
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

#if defined(REFINEMENT_MERGE_CELLS) && !defined(ONEDIMS)

static int derefine_criterion_windtunnel(int i);
static int derefine_criterion_jeans_ref(int i);
static int derefine_criterion_default(int i);
static int derefine_criterion_special_boundary(int i);
static int derefine_criterion_gmcs(int i);
static int derefine_criterion_anna(int i);
static int derefine_criterion_shellresol(int i);
static int jeans_derefinement_criteria(int i);
static int anna_derefinement_criteria(int i);
static int derefinement_criteria_ion(int i);
static int derefinement_criteria_simplex(int i);

#ifdef ROTATING_HIGHRES_REGION
static int derefine_rotating_highres_region(int i);
#endif

#if defined(SPIRAL) && defined(RAMP_REFINE)
static int derefine_disczoom(int i);
#endif

#if defined(REFINE_MCTR) && defined(TRACER_MC)
static int derefine_criterion_mctr(int i);
#endif

#ifdef BH_BASED_CGM_ZOOM
static int derefine_criterion_cgm(int i);
#endif

/*! \brief Should this cell be dissolved?
 *
 *  This function signals whether a cell should be dissolved. This needs to be
 *  adjusted according to the needs of the simulation in question. One may also
 *  set the SphP[].Flag variable beforehand, these cells will also be
 *  dissolved.
 *
 *  \param[in] i Index of cell in P and SphP arrays.
 *  \param[in] flag If this is nonzero, flag is returned.
 *
 *  \return Flag if this cell should be dissolved.
 */
int derefine_should_this_cell_be_merged(int i, int flag)
{
#ifdef REFINEMENT_HIGH_RES_GAS
  if(SphP[i].AllowRefinement == 0)
    return 0;
#endif

#ifdef NODEREFINE_BACKGROUND_GRID
  if(SphP[i].Volume > 0.1 * All.MeanVolume)
    return 0;
#endif

#ifdef BOUNDARY_FLAG
  if(SphP[i].BoundaryFlag & 1)
    return 0;
#endif

#if defined(BOUNDARY_INFLOWOUTFLOW_MINID) && defined(BOUNDARY_INFLOWOUTFLOW_MAXID)
  if(P[i].ID >= BOUNDARY_INFLOWOUTFLOW_MINID && P[i].ID < BOUNDARY_INFLOWOUTFLOW_MAXID)
    return 0;
#endif

#if defined(BOUNDARY_STICKY_MINID) && defined(BOUNDARY_STICKY_MAXID)
  if(P[i].ID >= BOUNDARY_STICKY_MINID && P[i].ID < BOUNDARY_STICKY_MAXID)
    return 0;
#endif

#ifdef STICKYFLAGS
  if(SphP[i].StickyFlag > 0)
    return 0;

  double boundary_dist;
  boundary_dist = fmax(boxSize_X, fmax(boxSize_Y, boxSize_Z));

#if(REFLECTIVE_X == 2)
  boundary_dist = fmin(boundary_dist, fmin(boxSize_X - P[i].Pos[0], P[i].Pos[0]));
#endif
#if(REFLECTIVE_Y == 2)
  boundary_dist = fmin(boundary_dist, fmin(boxSize_Y - P[i].Pos[1], P[i].Pos[1]));
#endif
#if(REFLECTIVE_Z == 2)
  boundary_dist = fmin(boundary_dist, fmin(boxSize_Z - P[i].Pos[2], P[i].Pos[2]));
#endif

  if(boundary_dist < All.StickyLayerMaxDist)
    return 1;

#endif

#if defined(BOUNDARY_REFL_FLUIDSIDE_MINID) && defined(BOUNDARY_REFL_FLUIDSIDE_MAXID)
  if(P[i].ID >= BOUNDARY_REFL_FLUIDSIDE_MINID && P[i].ID < BOUNDARY_REFL_FLUIDSIDE_MAXID)
    return 0;
#endif

#if defined(BOUNDARY_REFL_SOLIDSIDE_MINID) && defined(BOUNDARY_REFL_SOLIDSIDE_MAXID)
  if(P[i].ID >= BOUNDARY_REFL_SOLIDSIDE_MINID && P[i].ID < BOUNDARY_REFL_SOLIDSIDE_MAXID)
    return 0;
#endif

#ifdef SPECIAL_BOUNDARY
  if(P[i].ID == -1 || P[i].ID == -2 || P[i].ID <= -3)
    return 0;
#endif

#ifdef DEREFINE_GENTLY
  if(SphP[i].DoNotDerefFlag)
    return 0;
#endif

#if defined(REFINEMENT_KEEP_INITIAL_VOLUME)
  if(SphP[i].Volume > 0.5 * SphP[i].InitialVolume)
    return 0;

  if(SphP[i].Volume < 0.5 * SphP[i].InitialVolume)
    return 1;
#endif

#if defined(REFINEMENT_VOLUME_LIMIT) && !defined(WINDTUNNEL_REFINEMENT_VOLUME_LIMIT)
  double maxvolume = All.MaxVolume;
  double minvolume = All.MinVolume;
#ifdef REFINE_ONLY_WITH_TRACER
  if(P[i].TracerHead != -1)
    {
      maxvolume = All.MaxTracerVolume;
      minvolume = All.MinTracerVolume;
    }
#endif

  if(SphP[i].Volume > 0.5 * maxvolume)
    return 0;

  if(SphP[i].Volume < 0.5 * minvolume)
    return 1;

#ifdef REFINEMENT_VOLUME_LIMIT_MASS_LIMIT
  if(P[i].Mass < 0.1 * All.TargetGasMass)
    return 1;
#endif

  if(All.MaxVolumeDiff > 0 && SphP[i].Volume > 0.3 * All.MaxVolumeDiff * SphP[i].MinNgbVolume)
    return 0;
#endif

#if defined(CIRCUMSTELLAR) && defined(CIRCUMSTELLAR_REFINEMENTS)
  double StellarDistance;
  StellarDistance = get_circumstellar_distance(i);
  if((StellarDistance < All.CircumstellarDerefinementDistance) &&
     (SphP[i].Volume < 0.5 * 0.75 / M_PI / pow(All.CircumstellarDerefinementDistance, 3)))
    return 1;
#endif

#ifdef REFINEMENT_AROUND_BH
  if(SphP[i].RefBHFlag)
    {
      if(SphP[i].RefBHMaxRad / All.RefBHLowerFactorC > get_cell_radius(i))
        return 1;
      return 0; /* We are in refining region, but the correct size */
    }
#endif

#ifdef REFINEMENT_AROUND_DM
  int j;
  double dx, dy, dz, dist;
  double cuberoot2 = 1.2599210498948732;
  for(j = 0; j < All.TotPartDM; j++)
    {
      dx   = P[i].Pos[0] - DMPartListGlobal[j].pos[0];
      dy   = P[i].Pos[1] - DMPartListGlobal[j].pos[1];
      dz   = P[i].Pos[2] - DMPartListGlobal[j].pos[2];
      dist = sqrt(dx * dx + dy * dy + dz * dz);
      if(dist < 1.0 * DMPartListGlobal[j].softening)
        {
          if(get_cell_radius(i) > DMPartListGlobal[j].softening / (All.RefinementCellsPerSoftening * cuberoot2) * 0.8)
            return 0;
          else
            return 1;
        }
    }
#endif

  if(flag)
    return flag;

  switch(All.DerefinementCriterion)
    {
      case 0:
        return 0;
        break;

      case 1:
        return derefine_criterion_default(i);
        break;

      case 2:
        return derefine_criterion_jeans_ref(i);
        break;

      case 3:
        return derefine_criterion_windtunnel(i);
        break;

      case 4:
        return derefine_criterion_special_boundary(i);
        break;

      case 5:
        return derefine_criterion_gmcs(i);
        break;

#ifdef RAMP_REFINE
      case 6:
        return derefine_disczoom(i);
        break;
#endif

#ifdef PASSIVE_SCALARS
      case 8:
        return derefine_criterion_shellresol(i);
        break;
#endif

#ifdef SGCHEM
      case 7:
        return derefine_criterion_anna(i);
        break;
#endif

#ifdef DISC_REFINE_ONLY
      case 8:
        return derefine_disconly(i);
        break;
#endif

#ifdef ROTATING_HIGHRES_REGION
      case 9:
        return derefine_rotating_highres_region(i);
        break;
#endif

#if defined(REFINE_MCTR) && defined(TRACER_MC)
      case 10:
        return derefine_criterion_mctr(i);
#endif
#ifdef BH_BASED_CGM_ZOOM
      case 10:
        return derefine_criterion_cgm(i);
        break;
#endif

      case 11:
        return derefinement_criteria_ion(i);
        break;

      case 12:
        return derefinement_criteria_simplex(i);
        break;

      default:
        terminate("invalid derefinement criterion specified");
        break;
    }

  return 0;
}

static int derefine_criterion_special_boundary(int i)
{
#ifdef SPECIAL_BOUNDARY
  if(SphP[i].MinDistBoundaryCell < All.BoundaryLayerScaleFactor && P[i].ID > 0)
    return 1;
#endif

  return 0;
}

static int derefine_criterion_default(int i)
{
#if defined(REFINEMENT_SPLIT_CELLS) && defined(REFINEMENT_MERGE_CELLS)

#ifdef DEREFINE_ONLY_DENSE_GAS

  if(P[i].Mass < 0.5 * All.TargetGasMass && SphP[i].Density * All.cf_a3inv >= 1.0e-6 * All.PhysDensThresh)
#else

  double TargetGasMass      = All.TargetGasMass;
#ifdef REFINEMENT_LIMIT_STARFORMING_GAS
  double eos_dens_threshold = All.PhysDensThresh;
#ifdef MODIFIED_EOS
  eos_dens_threshold *= All.FactorDensThresh;
#endif

  double dens = SphP[i].Density * All.cf_a3inv;
  if(dens >= 4. * eos_dens_threshold)
    {
      TargetGasMass = All.TargetGasMass * dens / (4. * eos_dens_threshold);
      if(TargetGasMass > All.HighDensityMaxGasDerefinementFactor * All.TargetGasMass)
        TargetGasMass = All.HighDensityMaxGasDerefinementFactor * All.TargetGasMass;
    }
#endif

#ifdef REFINEMENT_CGM
  if(SphP[i].HighResMassCGM > HIGHRESMASSFAC * P[i].Mass)
    {
      if(P[i].Mass < 0.5 * TargetGasMass && SphP[i].Volume < 0.5 * All.TargetGasVolume * All.cf_a3inv)
        return 1;
    }
  else
    {
      if(P[i].Mass < 0.5 * TargetGasMass)
        return 1;
    }
#else
  if(P[i].Mass < 0.5 * TargetGasMass)
    return 1;
#endif

#endif

#endif

    return 0;
}

#ifdef BH_BASED_CGM_ZOOM
/*! \brief Derefinement criterion based on distance from a single central blackhole.
 *
 * \param[in] i Gas cell index into SphP.
 *
 * \return void
 */
static int derefine_criterion_cgm(int i)
{
#if defined(REFINEMENT_SPLIT_CELLS) && defined(REFINEMENT_MERGE_CELLS)
  double tgm;

  // if there isn't yet a flagged BH, then All.BlackHolePosition is set to [-1,-1,-1]
  if(All.BlackHolePosition[0] < 0.0)
    {
      tgm = All.TargetGasMass;
    }
  else
    {
      // calculate distance to All.BlackHolePosition
      double dx = (P[i].Pos[0] - All.BlackHolePosition[0]);
      double dy = (P[i].Pos[1] - All.BlackHolePosition[1]);
      double dz = (P[i].Pos[2] - All.BlackHolePosition[2]);
      double rr = sqrt(dx * dx + dy * dy + dz * dz);

      // get distance dependent target gas mass
      tgm = distance_dependent_target_mass(rr);
    }

#ifdef DEREFINE_ONLY_DENSE_GAS
  if(P[i].Mass < 0.5 * tgm && SphP[i].Density * All.cf_a3inv >= 1.0e-6 * All.PhysDensThresh)
#else
  if(P[i].Mass < 0.5 * tgm)
#endif
    return 1;
#endif  // REFINEMENT_SPLIT_CELLS && REFINEMENT_MERGE_CELLS

  return 0;
}
#endif /* #ifdef BH_BASED_CGM_ZOOM */

static int derefine_criterion_jeans_ref(int i)
{
#ifdef TGSET
  if(tgset_jeans_ref(1, i))
    return 1;
/*
#ifdef SINKS
    if(SphP[i].InAccrRadius)
      if(P[i].Mass < MIN_TARGET_MASS_FACTOR_FOR_ACC * All.TargetGasMass)
        return 1;
#endif
*/
#else
#ifdef JEANS_REFINEMENT
  return jeans_derefinement_criteria(i);
#endif
#endif

  return 0;
}

static int derefine_criterion_windtunnel(int i)
{
#ifdef WINDTUNNEL
#ifdef WINDTUNNEL_REFINEMENT_VOLUME_LIMIT

  if(SphP[i].Volume > 0.5 * All.MaxVolume)
    return 0;

  if(SphP[i].Volume < 0.5 * All.MinVolume)
    return 1;

  if(P[i].Pos[WINDTUNNEL_COORD] < All.InjectionRegion)  // injection region
    {
      // if(SphP[i].Volume < 0.5 * All.InjectionVolume)
      //    return 1;
    }
  else if(P[i].Pos[WINDTUNNEL_COORD] > All.InjectionRegion &&
          All.BoxSizes[WINDTUNNEL_COORD] - P[i].Pos[WINDTUNNEL_COORD] > All.InjectionRegion)  // the "normal region"
    {
      if(All.MaxVolumeDiff > 0 && SphP[i].Volume > 0.3 * All.MaxVolumeDiff * SphP[i].MinNgbVolume)
        return 0;
    }
  else if(All.BoxSizes[WINDTUNNEL_COORD] - P[i].Pos[WINDTUNNEL_COORD] < All.InjectionRegion)  // outflow region
    {
      // if(SphP[i].Volume < 0.8 * All.InjectionVolume)
      //    return 1;
    }

#endif

  if(P[i].Mass < 0.5 * All.TargetGasMass && P[i].Pos[WINDTUNNEL_COORD] > All.InjectionRegion)
    return 1;

#endif

  return 0;
}

static int derefine_criterion_gmcs(int i)
{
#ifdef GMC_REFINEMENT
  return gmc_derefinement_criteria(i);
#else
  return 0;
#endif
}

#if defined(REFINE_MCTR) && defined(TRACER_MC)
int derefine_criterion_mctr(int i)
{
  double dx    = P[i].Pos[0] - boxHalf_X;
  double dy    = P[i].Pos[1] - boxHalf_Y;
  double dz    = P[i].Pos[2] - boxHalf_Z;
  double erre2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);

  double zheight = 100.;
  double zoomtarget;
  double MCtarget = 1.;

  if(P[i].NumberOfTracers > 0 && fabs(dz) < zheight)
    zoomtarget = MCtarget;
  else if(erre2 < pow(20. * KILOPARSEC / All.UnitLength_in_cm, 2))
    zoomtarget = All.TargetGasMass;
  else if(erre2 < pow(30. * KILOPARSEC / All.UnitLength_in_cm, 2))
    zoomtarget = 10. * All.TargetGasMass;
  else if(erre2 < pow(40. * KILOPARSEC / All.UnitLength_in_cm, 2))
    zoomtarget = 100. * All.TargetGasMass;
  else
    zoomtarget = 1000. * All.TargetGasMass;

  if(P[i].Mass < 0.5 * zoomtarget)
    {
      return 1;
    }

  return 0;
}
#endif

static int derefine_criterion_anna(int i)
{
  if(P[i].Mass < 0.5 * All.TargetGasMass)
    return 1;
#ifdef SGCHEM
  return anna_derefinement_criteria(i); /* jeans criterion here for high densities*/
#else
  return 0;
#endif
}

static int derefine_criterion_shellresol(int i)
{
#if defined(REFINEMENT_SPLIT_CELLS) && defined(REFINEMENT_MERGE_CELLS)
#ifdef PASSIVE_SCALARS
  if(SphP[i].PScalars[0] > 0.01)
    {
      if(can_this_cell_be_split(i) && P[i].Mass < 0.5 * All.TargetGasMass / 10)
        return 1;
    }
  else
    {
      if(can_this_cell_be_split(i) && P[i].Mass < 0.5 * All.TargetGasMass)
        return 1;
    }
#endif
  return 0;
#endif
}

int jeans_derefinement_criteria(int i)
{
#ifndef NO_TARGET_MASS_CONDITION
  if(P[i].Mass < 0.5 * All.TargetGasMass)
    return 1;
#endif

#ifdef JEANS_DEREFINEMENT_DENSITY_THRESHOLD
  if(SphP[i].Density < JEANS_DEREFINEMENT_DENSITY_THRESHOLD)
    return 0;
#endif

#ifdef JEANS_REFINEMENT
  double jeans_number, jeans_length, sound_speed, dx;
#ifdef VARIABLE_GAMMA
  sound_speed = sqrt(SphP[i].GammaE * SphP[i].Pressure / SphP[i].Density);
#else
  sound_speed = sqrt(GAMMA * SphP[i].Pressure / SphP[i].Density);
#endif
  jeans_length = sqrt(M_PI / All.G / SphP[i].Density) * sound_speed;
  dx           = 2.0 * get_cell_radius(i);
  jeans_number = jeans_length / dx;
  /*For comoving integration, the a-dependencies do not vanish! */
  if(All.ComovingIntegrationOn)
    {
      jeans_number = jeans_number * sqrt(All.Time);
    }

  if(jeans_number > 1.5 * JEANS_REFINEMENT && P[i].Mass < 0.5 * All.TargetGasMass)
    return 1;
#endif
  return 0;
}

int anna_derefinement_criteria(int i)
{
  if(P[i].Mass < 0.5 * All.TargetGasMass)
    return 1;

#ifdef SGCHEM
  double n_dens_anna;
  double n_dens_anna_min;

  n_dens_anna     = SphP[i].Density * All.UnitDensity_in_cgs * All.cf_a3inv / (All.HubbleParam * All.HubbleParam);
  n_dens_anna     = n_dens_anna * 0.81967 / PROTONMASS;
  n_dens_anna_min = 1.0e2; /* for number densities larger than 100, use jeans refinement! */
  if(n_dens_anna > n_dens_anna_min)
    {
#ifdef JEANS_REFINEMENT
      return jeans_derefinement_criteria(i);
#endif
    }
  return derefine_criterion_default(i);
#endif
  return 0;
}

int derefinement_criteria_ion(int i)
{
  if(P[i].Mass > 0.5 * All.TargetGasMass)
    return 0;

  double phys_dens = SphP[i].Density * All.UnitDensity_in_cgs * All.cf_a3inv * (All.HubbleParam * All.HubbleParam);
  double M_strom   = 6.2025e+10 / phys_dens / All.UnitMass_in_g * All.HubbleParam;
  if(All.ComovingIntegrationOn && All.Time < 0.025)
    M_strom = M_strom * 1.0e20;  // deactive stromgren refinement at high-redshift
  if(P[i].Mass > 0.5 * M_strom)
    return 0;
#ifdef SINK_PARTICLES_FEEDBACK
  for(int j = 0; j < NSinksAllTasks; j++)
    {
      double dist = (P[i].Pos[0] - SinkP[j].Pos[0]) * (P[i].Pos[0] - SinkP[j].Pos[0]) +
                    (P[i].Pos[1] - SinkP[j].Pos[1]) * (P[i].Pos[1] - SinkP[j].Pos[1]) +
                    (P[i].Pos[2] - SinkP[j].Pos[2]) * (P[i].Pos[2] - SinkP[j].Pos[2]);
      if(dist < pow((12 * PARSEC / All.UnitLength_in_cm * All.HubbleParam / All.cf_atime), 2))
        {
          if(get_cell_radius(i) > 0.075 * PARSEC / All.UnitLength_in_cm * All.HubbleParam / All.cf_atime)
            return 0;
        }
    }
#endif
#ifdef JEANS_REFINEMENT
  return jeans_derefinement_criteria(i);
#else
  return 1;
#endif
}

int derefinement_criteria_simplex(int i)
{
  if(P[i].Mass > 0.5 * All.TargetGasMass)
    return 0;
#ifdef SINK_PARTICLES
  for(int j = 0; j < NSinksAllTasks; j++)
    {
      double dist = (P[i].Pos[0] - SinkP[j].Pos[0]) * (P[i].Pos[0] - SinkP[j].Pos[0]) +
                    (P[i].Pos[1] - SinkP[j].Pos[1]) * (P[i].Pos[1] - SinkP[j].Pos[1]) +
                    (P[i].Pos[2] - SinkP[j].Pos[2]) * (P[i].Pos[2] - SinkP[j].Pos[2]);
      if(dist < 1.2 * SinkFormationRadius * SinkFormationRadius)
        {
          double phys_dens = SphP[i].Density * All.UnitDensity_in_cgs * All.cf_a3inv * (All.HubbleParam * All.HubbleParam);
          double M_strom   = 6.2025e+10 / phys_dens / All.UnitMass_in_g * All.HubbleParam;
          if(P[i].Mass > 0.5 * M_strom)
            return 0;
        }
    }
#endif
#ifdef JEANS_REFINEMENT
  return jeans_derefinement_criteria(i);
#else
  return 1;
#endif
}

#ifdef ROTATING_HIGHRES_REGION
static int derefine_rotating_highres_region(int i)
{
  double dx    = P[i].Pos[0] - boxHalf_X;
  double dy    = P[i].Pos[1] - boxHalf_Y;
  double dz    = P[i].Pos[2] - boxHalf_Z;
  double erre2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);

  double zoomtarget;

  if(erre2 < pow(20. * KILOPARSEC / All.UnitLength_in_cm, 2))
    zoomtarget = All.TargetGasMass;
  else if(erre2 < pow(30. * KILOPARSEC / All.UnitLength_in_cm, 2))
    zoomtarget = 10. * All.TargetGasMass;
  else if(erre2 < pow(40. * KILOPARSEC / All.UnitLength_in_cm, 2))
    zoomtarget = 100. * All.TargetGasMass;
  else
    zoomtarget = 1000. * All.TargetGasMass;

  if(fabs(dz) < All.Highres_deltaz)
    {
      double Rref = sqrt(pow((All.Highres_y0 - boxHalf_Y), 2) +
                         pow((All.Highres_x0 - boxHalf_X), 2));  // radius at which the high res region is centered

      double dR = sqrt(dx * dx + dy * dy);

      if(fabs(Rref - dR) < All.Highres_deltaR)
        {
          double theta_0 = atan2((All.Highres_y0 - boxHalf_Y), (All.Highres_x0 - boxHalf_X));  // initial angle of the high res region
          double wr      = All.Highres_vrot / Rref;  // angular velocity of the high res region

          double thref = theta_0 + wr * (All.Time - All.Highres_t0);  // current angle of the high res region
          thref -= (thref < 0 ? -1. : 1.) * floor((fabs(thref) + M_PI) / (2. * M_PI)) * 2. *
                   M_PI;  // so that thref is always between -pi and pi
          double dtheta = atan2(dy, dx);

          double delta_th = dtheta - thref;

          if(delta_th > M_PI)
            delta_th = 2. * M_PI - delta_th;
          else if(delta_th < -M_PI)
            delta_th = 2. * M_PI + delta_th;

          delta_th = fabs(delta_th);

          if((delta_th * dR) < All.Highres_deltaThetaxR)
            zoomtarget = All.Highres_targetmass;
        }
    }

  if(P[i].Mass < 0.5 * zoomtarget)
    return 1;

  return 0;
}
#endif

#ifdef SPIRAL
static int derefine_disczoom(int i)
{
  /* This routine refines a segment of a galactic disc that moves with the
     mean gas velocity at 7.5 kpc. By Rowan Smith 2013. */
  double zoomtarget;
  double wr, thref, dtheta, theta_range, theta_range2, zrange, rinner, Rref, normalise;
  double rfloat, sp_ang;
  int revolve;

  Rref    = 7.5 * KILOPARSEC;
  wr      = -220.0 * 1.e5 / Rref; /* vr of 220 kms at 7.5 kpc */
  thref   = wr * All.Time * All.UnitTime_in_s;
  revolve = (int)(fabs(thref) / M_PI);
  rfloat  = (float)revolve;
  thref   = thref + rfloat * 2.0 * M_PI; /*Note this assumes the rotation is clockwise */

  theta_range  = M_PI / 4.;
  theta_range2 = theta_range / 2.;
  zrange       = KILOPARSEC / All.UnitLength_in_cm;      /*Only refine 1kpc from disc plane*/
  rinner       = 4. * KILOPARSEC / All.UnitLength_in_cm; /*Inner radius set to 4 kpc*/
  normalise    = 2.3e-13;                                /*2.3e-13 gives 5 Msun mass resolution*/

  double dx = P[i].Pos[0] - boxHalf_X;
  double dy = P[i].Pos[1] - boxHalf_Y;
  double dz = P[i].Pos[2] - boxHalf_Z;

  double rr = sqrt(dx * dx + dy * dy);

  zoomtarget = All.TargetGasMass;

  if(dz < boxHalf_Z + zrange && dz > boxHalf_Z - zrange && rr > rinner)
    {
      dtheta = atan2(dy, dx) - thref;

      if(dtheta > M_PI)
        dtheta = 2. * M_PI - dtheta;
      if(dtheta < -1.0 * M_PI)
        dtheta = 2.0 * M_PI + dtheta;

      sp_ang = fabs(dtheta);

#ifdef RAMP_REFINE
      if(sp_ang < theta_range && sp_ang > theta_range2)
        {
          zoomtarget = All.TargetGasMass / 1000.0 * (2490.0 * (sp_ang - M_PI / 8.0) + 5.0);
#ifdef SNE_RAMP_REFINE
          zoomtarget = All.TargetGasMass / 1000.0 * 2.0 * (1245.0 * (sp_ang - M_PI / 8.0) + 5.0);
#endif
        }
      if(zoomtarget < 100. * SOLAR_MASS / All.UnitMass_in_g) /* to avoid sharp contrasts use exponential */
        zoomtarget = normalise * exp(78.2 * sp_ang) * SOLAR_MASS / All.UnitMass_in_g;

      /*Now assign the inner region where resolution is high */
      if(sp_ang <= theta_range2)
        {
          zoomtarget = normalise * exp(78.2 * theta_range2) * SOLAR_MASS / All.UnitMass_in_g;
#ifdef SNE_RAMP_REFINE
          zoomtarget = 10.0 * SOLAR_MASS / All.UnitMass_in_g;
#endif
        }
#else
      if(sp_ang < theta_range)
        zoomtarget = All.TargetGasMass / 10.0;
#endif
    }

  if(P[i].Mass < 0.5 * zoomtarget)
    return 1;

  return 0;
}
#endif

#ifdef GMC_REFINEMENT
int gmc_derefinement_criteria(int i)
{
  /* PCC - 08.05.2013
     This file contains the different refinement criteria for
     modelling molecular clouds. At the moment we have a simple Jeans
     criterion, but we should eventually add a shock criterion. Note
     that the Jeans length definition here is different (smaller) than
     that used in the standard Jeans refinement. This code also ignores
     the cell mass, unlike the other version.
  */

  double jeans_temperature, jeans_length;
  double cells_per_jeans_length;
  double density;
  double three_over_four_pi, temp_prefactor;
  double cellrad;
  double numdens_ref_min, numdens_ref_max;
  double yn, energy, yntot;

  three_over_four_pi = 0.238732414637843;
  /* 5.*kb/2./gg/2.33/mp * T */
  temp_prefactor = 1.32740615361024e+15;

  /* make sure density is within range */
  density = SphP[i].Density * All.UnitDensity_in_cgs;
  if((density < All.GMCDerefMinDensity) || (density > All.GMCRefMaxDensity))
    return 0;

    /* Get the temperature of the gas  */
#ifdef SGCHEM
  yn    = density / ((1.0 + 4.0 * ABHE) * PROTONMASS);
  yntot = (1 + ABHE - SphP[i].TracAbund[IH2] + SphP[i].TracAbund[IHP]) * yn;
#else
  yntot             = density / 2.33 / PROTONMASS;
#endif

  energy = (SphP[i].Utherm * All.UnitEnergy_in_cgs / All.UnitMass_in_g) * density;
#ifdef VARIABLE_GAMMA
  jeans_temperature = (SphP[i].GammaE - 1.0) * energy / (yntot * BOLTZMANN);
#else
  jeans_temperature = 2 * energy / (3 * yntot * BOLTZMANN);
#endif

  /* Calculate the number of cells per Jeans length */
  jeans_length           = 2.0 * sqrt(three_over_four_pi / density) * sqrt(temp_prefactor * jeans_temperature);
  cellrad                = 2.0 * get_cell_radius(i) * All.UnitLength_in_cm;
  cells_per_jeans_length = jeans_length / cellrad;

  /* Fix cells per jeans length here... */
  if(cells_per_jeans_length > 1.5 * All.GMCRefCellsPerJeansLength && P[i].Mass < 0.5 * All.TargetGasMass)
    return (1);

  return (0);
}
#endif

#ifdef DISC_REFINE_ONLY
int derefine_disconly(int i)
{
  double zoomtarget;
  double dx = P[i].Pos[0] - boxHalf_X;
  double dy = P[i].Pos[1] - boxHalf_Y;
  double dz = P[i].Pos[2] - boxHalf_Z;
  double rr = sqrt(dx * dx + dy * dy);

  double zrange = KILOPARSEC / All.UnitLength_in_cm; /*Only refine 1kpc from plane*/
  double rinner = 3. * KILOPARSEC / All.UnitLength_in_cm;
  double router = 11. * KILOPARSEC / All.UnitLength_in_cm;

  zoomtarget = All.TargetGasMass;
  if(fabs(dz) < zrange && rr > rinner && rr < router)
    {
      zoomtarget = All.TargetGasMass / 100.;
    }

  if(P[i].Mass < 0.5 * zoomtarget)
    return 1;
  return 0;
}
#endif

#endif
