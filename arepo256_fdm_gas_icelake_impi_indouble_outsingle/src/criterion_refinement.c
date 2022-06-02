/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/criterion_refinement.c
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

#if defined(REFINEMENT_SPLIT_CELLS) && !defined(ONEDIMS)

static int refine_criterion_anna(int i);
static int refine_criterion_gmcs(int i);
static int refine_criterion_shellresol(int i);
static int anna_refinement_criteria(int i);
static int jeans_refinement_criteria(int i);
static int refinement_criteria_ion(int i);
static int refinement_criteria_simplex(int i);

#ifdef ROTATING_HIGHRES_REGION
static int refine_rotating_highres_region(int i);
#endif

#if defined(SPIRAL) && defined(RAMP_REFINE)
static int refine_disczoom(int i);
#endif

#if defined(REFINE_MCTR) && defined(TRACER_MC)
static int refine_criterion_mctr(int i);
#endif

/*! \brief Should this cell be refined?
 *
 *  This function signals whether a cell needs further refinement. This needs
 *  to be adjusted according to the needs of the simulation in question.
 *
 *  \param[in] i Index of cell in P and SphP arrays.
 *
 *  \return Flag if this cell should be split.
 */
int should_this_cell_be_split(int i)
{
#ifdef REFINEMENT_MERGE_CELLS
  if(FlagDoNotRefine[i])
    return 0;
#endif

  if(P[i].Mass == 0 && P[i].ID == 0) /* skip cells that have been swallowed or dissolved */
    return 0;

#ifdef BOUNDARY_FLAG
  if(SphP[i].BoundaryFlag & 1)
    return 0;
#endif

#ifdef BOUNDARY_INFLOWOUTFLOW_MINID
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
  if(P[i].ID == -1 || P[i].ID == -2)
    return 0;
#endif

#if defined(REFINEMENT_KEEP_INITIAL_VOLUME)
  if(SphP[i].Volume > 2.0 * SphP[i].InitialVolume)
    if(can_this_cell_be_split(i))
      return 1;

  if(SphP[i].Volume < 2.0 * SphP[i].InitialVolume)
    return 0;
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

  if(SphP[i].Volume > 2. * maxvolume)
    if(can_this_cell_be_split(i))
      return 1;

  if(SphP[i].Volume < 2. * minvolume)
    return 0;

#ifdef REFINEMENT_VOLUME_LIMIT_MASS_LIMIT
  if(P[i].Mass < 0.1 * All.TargetGasMass)
    return 0;
#endif

  if(refine_criterion_volume(i))
    if(can_this_cell_be_split(i))
      return 1;
#endif

#ifdef REFINEMENT_AROUND_BH
  if(SphP[i].RefBHFlag)
    {
      if(SphP[i].RefBHMaxRad < get_cell_radius(i))
        {
          if(P[i].Mass < All.RefBHMinCellMass)
            {
              return 0;
            }
#if(REFINEMENT_AROUND_BH == 0)
          if(can_this_cell_be_split(i))
            return 1;
#endif
#if(REFINEMENT_AROUND_BH == 1)
          return 1;
#endif
        }
    }
#endif

#ifdef REFINEMENT_AROUND_DM
  int j;
  for(j = 0; j < All.TotPartDM; j++)
    {
      double dx   = P[i].Pos[0] - DMPartListGlobal[j].pos[0];
      double dy   = P[i].Pos[1] - DMPartListGlobal[j].pos[1];
      double dz   = P[i].Pos[2] - DMPartListGlobal[j].pos[2];
      double dist = sqrt(dx * dx + dy * dy + dz * dz);
      if(dist < 1.0 * DMPartListGlobal[j].softening)
        if(get_cell_radius(i) > DMPartListGlobal[j].softening / All.RefinementCellsPerSoftening && can_this_cell_be_split(i))
          {
#ifdef DEBUG_REFINE
            printf("Refining particle ID %d with radius %e around DM particle %d with softening length %e\n", P[i].ID,
                   get_cell_radius(i), DMPartListGlobal[j].ID, DMPartListGlobal[j].softening);
#endif
            return 1;
          }
    }
#endif

#if defined(CIRCUMSTELLAR) && defined(CIRCUMSTELLAR_REFINEMENTS)
  double StellarDistance;
  StellarDistance = get_circumstellar_distance(i);
  if((StellarDistance < All.CircumstellarDerefinementDistance) &&
     (SphP[i].Volume < 0.75 / M_PI / pow(All.CircumstellarDerefinementDistance, 3)))
    return 0;
#endif

#if defined SINK_PARTICLES_REFINEMENT_LIMIT
  int isink;
#ifdef JEANS_REFINEMENT
  double ncells_per_sink = JEANS_REFINEMENT;
#else
  double ncells_per_sink = 4;
#endif
  double dx, dy, dz, dist, dist_between_sink_cell_surface;
  double cell_rad = get_cell_radius(i);
  int number_of_sink_radii;
  if(can_this_cell_be_split(i))
    {
      if(NSinksAllTasks > 0)
        {
          for(isink = 0; isink < NSinksAllTasks; isink++)
            {
              /* Ensure that the regions around the sink particles are well resolved, and that the
                 cells refine smoothly around the sinks.
              */
              dx                             = P[i].Pos[0] - SinkP[isink].Pos[0];
              dy                             = P[i].Pos[1] - SinkP[isink].Pos[1];
              dz                             = P[i].Pos[2] - SinkP[isink].Pos[2];
              dist                           = sqrt(dx * dx + dy * dy + dz * dz);
              dist_between_sink_cell_surface = dist - cell_rad - SinkAccretionRadius;
              if(dist_between_sink_cell_surface < 0)
                dist_between_sink_cell_surface = 0;
              number_of_sink_radii = dist_between_sink_cell_surface / SinkAccretionRadius;
              if((number_of_sink_radii + 1) * dist_between_sink_cell_surface / ncells_per_sink > cell_rad)
                return 1;
            }
        }
    }
#endif

#if defined(MHD) && defined(MHD_REFINE_ON_DIVB_FACTOR)
  if(can_this_cell_be_split(i))
    {
      double Bmag                    = sqrt(SphP[i].B[0] * SphP[i].B[0] + SphP[i].B[1] * SphP[i].B[1] + SphP[i].B[2] * SphP[i].B[2]);
      double cell_rad                = get_cell_radius(i);
      double spurious_to_real_factor = SphP[i].DivB * 2.0 * cell_rad / Bmag;
      if(spurious_to_real_factor > MHD_REFINE_ON_DIVB_FACTOR)
        return 1;
    }
#endif

  switch(All.RefinementCriterion) /* select the function that evaluates the refinement criterion */
    {
      case 0:
        return 0;
        break;

      case 1:
        return refine_criterion_default(i);
        break;

      case 2:
        return refine_criterion_jeans_ref(i);
        break;

      case 3:
        return refine_criterion_windtunnel(i);
        break;

      case 4:
        return refine_criterion_special_boundary(i);
        break;

      case 5:
        return refine_criterion_gmcs(i);
        break;

#ifdef RAMP_REFINE
      case 6:
        return refine_disczoom(i);
        break;
#endif

#ifdef PASSIVE_SCALARS
      case 8:
        return refine_criterion_shellresol(i);
        break;
#endif

#ifdef SGCHEM
      case 7:
        return refine_criterion_anna(i);
        break;
#endif

#ifdef DISC_REFINE_ONLY
      case 8:
        return refine_disconly(i);
        break;
#endif

#ifdef ROTATING_HIGHRES_REGION
      case 9:
        return refine_rotating_highres_region(i);
#endif

#if defined(REFINE_MCTR) && defined(TRACER_MC)
      case 10:
        return refine_criterion_mctr(i);
#endif

#ifdef BH_BASED_CGM_ZOOM
      case 10:
        return refine_criterion_cgm(i);
        break;
#endif

      case 11:
        return refinement_criteria_ion(i);

      case 12:
        return refinement_criteria_simplex(i);

      default:
        terminate("invalid refinement criterion specified");
        break;
    }

  return 0;
}

int can_this_cell_be_split(int i)
{
#ifdef REFINE_ABOVE_WNM_DENSITY
  /* ~ 1 cm^-3 */
  if(SphP[i].Density * All.UnitDensity_in_cgs / PROTONMASS / 1.4 < 2.0)
    return 0;
#endif
#ifdef REGULARIZE_MESH_FACE_ANGLE
  if(SphP[i].MaxFaceAngle < 1.5 * All.CellMaxAngleFactor)
#else
  double dx      = nearest_x(P[i].Pos[0] - SphP[i].Center[0]);
  double dy      = nearest_y(P[i].Pos[1] - SphP[i].Center[1]);
  double dz      = nearest_z(P[i].Pos[2] - SphP[i].Center[2]);
  double d       = sqrt(dx * dx + dy * dy + dz * dz);
  double cellrad = get_cell_radius(i);

  if(d < 2.0 * All.CellShapingFactor * cellrad) /* only refine cells which are reasonably 'round' */
#endif
    return 1;

  return 0;
}

int refine_criterion_default(int i)
{
#ifdef REFINEMENT_HIGH_RES_GAS
  if(SphP[i].AllowRefinement != 0)
#ifndef TGSET
    if(SphP[i].HighResMass > HIGHRESMASSFAC * P[i].Mass)
#endif
#endif
      {
        double TargetGasMass = All.TargetGasMass;

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

        if(can_this_cell_be_split(i))
          {
            if(P[i].Mass > 2.0 * TargetGasMass)
              return 1;
#ifdef GRADIENTREFINEMENT
            double deltaBx_totsquare =
                pow(SphP[i].Grad.dB[0][0], 2.0) + pow(SphP[i].Grad.dB[0][1], 2.0) + pow(SphP[i].Grad.dB[0][2], 2.0);
            double deltaBy_totsquare =
                pow(SphP[i].Grad.dB[1][0], 2.0) + pow(SphP[i].Grad.dB[1][1], 2.0) + pow(SphP[i].Grad.dB[1][2], 2.0);
            double deltaBz_totsquare =
                pow(SphP[i].Grad.dB[2][0], 2.0) + pow(SphP[i].Grad.dB[2][1], 2.0) + pow(SphP[i].Grad.dB[2][2], 2.0);

            double deltaB_tot = sqrt(deltaBx_totsquare + deltaBy_totsquare + deltaBz_totsquare);
            double B_total    = sqrt(pow(SphP[i].B[0], 2.0) + pow(SphP[i].B[1], 2.0) + pow(SphP[i].B[2], 2.0));

            double deltavelx_totsquare =
                pow(SphP[i].Grad.dvel[0][0], 2.0) + pow(SphP[i].Grad.dvel[0][1], 2.0) + pow(SphP[i].Grad.dvel[0][2], 2.0);
            double deltavely_totsquare =
                pow(SphP[i].Grad.dvel[1][0], 2.0) + pow(SphP[i].Grad.dvel[1][1], 2.0) + pow(SphP[i].Grad.dvel[1][2], 2.0);
            double deltavelz_totsquare =
                pow(SphP[i].Grad.dvel[2][0], 2.0) + pow(SphP[i].Grad.dvel[2][1], 2.0) + pow(SphP[i].Grad.dvel[2][2], 2.0);

            double delta_v_tot = sqrt(deltavelx_totsquare + deltavely_totsquare + deltavelz_totsquare);
            double v_tot       = sqrt(pow(P[i].Vel[0], 2.0) + pow(P[i].Vel[1], 2.0) + pow(P[i].Vel[2], 2.0));

            double delta_rho_tot =
                sqrt(pow(SphP[i].Grad.drho[0], 2.0) + pow(SphP[i].Grad.drho[1], 2.0) + pow(SphP[i].Grad.drho[2], 2.0));
            double rho_tot = P[i].Mass / SphP[i].Volume;
            double critfac = 1.2;  // corresponds to an increase of factor 10 over 10 cells

            if((delta_rho_tot / rho_tot * get_cell_radius(i)) > (critfac))
              {
                return 1;
              }
            else if(((delta_v_tot / v_tot * get_cell_radius(i)) > (critfac)) && (All.GradVelRefinement == 1))
              {
                return 1;
              }
            else if(((deltaB_tot / B_total * get_cell_radius(i)) > (critfac)) && (All.GradBfldRefinement == 1))
              {
                return 1;
              }
#endif

#ifdef REFINEMENT_CGM
            if(SphP[i].HighResMassCGM > HIGHRESMASSFAC * P[i].Mass)
              if(SphP[i].Volume > 2.0 * All.TargetGasVolume * All.cf_a3inv)
                return 1;
#endif
          }
      }
  return 0; /* default is not to refine */
}

#ifdef BH_BASED_CGM_ZOOM
/*! \brief Define the distance dependent target gas mass (interpolant) function.
 *
 * \param[in] r Distance, code units.
 *
 * \return Target gas mass for this distance, code units.
 */
double distance_dependent_target_mass(double r)
{
  double CGM_MinRadius, CGM_MaxRadius, IGM_Radius;
  double redshift, m1, b1, m2, b2;

  // convert parameter file values from physical lengths at the specified redshift to (comoving) code
  CGM_MinRadius = All.CGM_MinRadius * (1.0 + All.CGM_RadiiRedshift) * All.HubbleParam;
  CGM_MaxRadius = All.CGM_MaxRadius * (1.0 + All.CGM_RadiiRedshift) * All.HubbleParam;
  IGM_Radius    = All.IGM_Radius * (1.0 + All.CGM_RadiiRedshift) * All.HubbleParam;

  double massratio = All.TargetGasMass / All.CGM_RefinementFactor;

  // interpolate
  m1 = (massratio - All.TargetGasMass) / CGM_MinRadius;
  b1 = All.TargetGasMass;
  m2 = (All.TargetGasMass - massratio) / (IGM_Radius - CGM_MaxRadius);
  b2 = All.TargetGasMass + IGM_Radius * (massratio - All.TargetGasMass) / (IGM_Radius - CGM_MaxRadius);

  if(r < CGM_MinRadius)
    {
      return m1 * r + b1;
    }
  else if(r >= CGM_MinRadius && r <= CGM_MaxRadius)
    {
      return massratio;
    }
  else if(r >= CGM_MaxRadius && r <= IGM_Radius)
    {
      return m2 * r + b2;
    }
  else if(r > IGM_Radius)
    {
      return All.TargetGasMass;
    }
}

/*! \brief Refinement criterion based on distance from a single central blackhole.
 *
 * \param[in] i Gas cell index into SphP.
 *
 * \return void
 */
int refine_criterion_cgm(int i)
{
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

#ifdef REFINEMENT_HIGH_RES_GAS
  if(SphP[i].AllowRefinement != 0)
#endif
    if(can_this_cell_be_split(i) && P[i].Mass > 2.0 * tgm)
      return 1;

  return 0;
}
#endif /* # ifdef BH_BASED_CGM_ZOOM */

int refine_criterion_jeans_ref(int i)
{
#ifdef REFINEMENT_HIGH_RES_GAS
  if(SphP[i].AllowRefinement != 0)
#ifndef TGSET
    if(SphP[i].HighResMass > HIGHRESMASSFAC * P[i].Mass)
#endif
#endif
      if(can_this_cell_be_split(i))
        {
          if(P[i].Mass > 2.0 * All.TargetGasMass)
            return 1;
#ifdef TGSET
          return tgset_jeans_ref(0, i);
#else
#ifdef JEANS_REFINEMENT
      return jeans_refinement_criteria(i);
#else
      return 0;
#endif
#endif
        }

  return 0;
}

int refine_criterion_gmcs(int i)
{
  if(can_this_cell_be_split(i))
    {
#ifdef GMC_REFINEMENT
      return gmc_refinement_criteria(i);  // currently only jeans criterion here */
#else
      return 0;
#endif
    }
  return 0;
}

int refine_criterion_anna(int i)
{
  if(can_this_cell_be_split(i))
    {
#ifdef SGCHEM
      return anna_refinement_criteria(i);  // jeans criterion here, applied to high densities */
#else
      return 0;
#endif
    }
  return 0;
}

int refine_criterion_shellresol(int i)
{
#ifdef PASSIVE_SCALARS
  if(SphP[i].PScalars[0] > 0.01)
    {
      if(can_this_cell_be_split(i) && P[i].Mass > 2.0 * All.TargetGasMass / 10)
        return 1;
    }
  else
    {
      if(can_this_cell_be_split(i) && P[i].Mass > 2.0 * All.TargetGasMass)
        return 1;
    }
#endif
  return 0;
}

int jeans_refinement_criteria(int i)
{
#ifdef SINK_PARTICLES
  if(SphP[i].Density > 0.1 * SinkCreationDensityCurrent)
    return 0;
#endif
#ifdef JEANS_REFINEMENT
  if(can_this_cell_be_split(i))
    {
      double jeans_number, jeans_length, sound_speed, dx;
#ifdef VARIABLE_GAMMA
      sound_speed = sqrt(SphP[i].GammaE * SphP[i].Pressure / SphP[i].Density);
#else
      sound_speed = sqrt(GAMMA * SphP[i].Pressure / SphP[i].Density);
#endif
      jeans_length = sqrt(M_PI / All.G / SphP[i].Density) * sound_speed;
      dx           = 2.0 * get_cell_radius(i);
      jeans_number = jeans_length / dx;
#ifdef EXTERNALSHEARBOX
      double bscale0    = 61.0 * (All.ShearBoxFg / 0.1 / All.ShearBoxMu) / (All.ShearBoxSigma0 / 10.0);
      double rho_thresh = All.ShearBoxSigma0 / 2.0 / bscale0;
      if(SphP[i].Density < rho_thresh)
        return 0;
#endif
      /*For comoving integration, the a-dependencies do not vanish!*/
      if(All.ComovingIntegrationOn)
        {
          jeans_number = jeans_number * sqrt(All.Time);
        }

      if(jeans_number < JEANS_REFINEMENT)
        {
          return 1;
        }
    }
#endif
  return 0;
}

#if defined(REFINE_MCTR) && defined(TRACER_MC)

int refine_criterion_mctr(int i)
{
  if(P[i].Type == 0)
    {
      double dx    = P[i].Pos[0] - boxHalf_X;
      double dy    = P[i].Pos[1] - boxHalf_Y;
      double dz    = P[i].Pos[2] - boxHalf_Z;
      double erre2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);

      double zheight = 100.;
      double zoomtarget;
      double MCtarget = 0.25;

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

#ifdef SINK_PARTICLES_VARIABLE_CREATION
      SphP[i].RefineTarget = zoomtarget;
#endif
      if(can_this_cell_be_split(i) && P[i].Mass > 2.0 * zoomtarget)
        return 1;
    }

  return 0;
}
#endif

int anna_refinement_criteria(int i)
{
  if(can_this_cell_be_split(i))
    {
      double n_dens_anna;
      double n_dens_anna_min;

      n_dens_anna     = SphP[i].Density * All.UnitDensity_in_cgs * All.cf_a3inv / (All.HubbleParam * All.HubbleParam);
      n_dens_anna     = n_dens_anna * 0.81967 / PROTONMASS;
      n_dens_anna_min = 1.0e2; /* for number densities larger than 100, use jeans refinement! */
#ifdef JEANS_REFINEMENT
      if(n_dens_anna > n_dens_anna_min)
        {
          return jeans_refinement_criteria(i);
        }

#endif
      return refine_criterion_default(i);
    }
  return 0;
}

int refinement_criteria_ion(int i)
{
  if(can_this_cell_be_split(i))
    {
#ifdef SINK_PARTICLES
      if(SphP[i].Density > SinkCreationDensityCurrent)
        // do not apply advanced refinement criteria above sink creation density threshold
        return 0;
#endif
#ifdef JEANS_REFINEMENT
      if(jeans_refinement_criteria(i) == 1)
        return 1;
#endif
      if(refine_criterion_default(i) == 1)
        return 1;
      double phys_dens = SphP[i].Density * All.UnitDensity_in_cgs * All.cf_a3inv * (All.HubbleParam * All.HubbleParam);
      double M_strom   = 6.2025e+10 / phys_dens / All.UnitMass_in_g * All.HubbleParam;
      if(All.ComovingIntegrationOn && All.Time < 0.025)
        M_strom = M_strom * 1.0e20;  // deactive stromgren refinement at high-redshift
      if(P[i].Mass > 2.0 * M_strom)
        return 1;
#ifdef SINK_PARTICLES_FEEDBACK
      for(int j = 0; j < NSinksAllTasks; j++)
        {
          double dist = (P[i].Pos[0] - SinkP[j].Pos[0]) * (P[i].Pos[0] - SinkP[j].Pos[0]) +
                        (P[i].Pos[1] - SinkP[j].Pos[1]) * (P[i].Pos[1] - SinkP[j].Pos[1]) +
                        (P[i].Pos[2] - SinkP[j].Pos[2]) * (P[i].Pos[2] - SinkP[j].Pos[2]);
          if(dist < pow((10 * PARSEC / All.UnitLength_in_cm * All.HubbleParam / All.cf_atime), 2))
            {
              if(get_cell_radius(i) > 0.3 * PARSEC / All.UnitLength_in_cm * All.HubbleParam / All.cf_atime)
                return 1;
            }
        }
#endif
    }
  return 0;
}

int refinement_criteria_simplex(int i)
{
  if(can_this_cell_be_split(i))
    {
#if defined(JEANS_REFINEMENT) && defined(SINK_PARTICLES)
      if(SphP[i].Density < SinkCreationDensityCurrent)
        if(jeans_refinement_criteria(i) == 1)
          return 1;
#endif
      if(refine_criterion_default(i) == 1)
        return 1;
#ifdef SINK_PARTICLES
      for(int j = 0; j < NSinksAllTasks; j++)
        {
          double dist = (P[i].Pos[0] - SinkP[j].Pos[0]) * (P[i].Pos[0] - SinkP[j].Pos[0]) +
                        (P[i].Pos[1] - SinkP[j].Pos[1]) * (P[i].Pos[1] - SinkP[j].Pos[1]) +
                        (P[i].Pos[2] - SinkP[j].Pos[2]) * (P[i].Pos[2] - SinkP[j].Pos[2]);
          if(dist < SinkFormationRadius * SinkFormationRadius)
            {
              double phys_dens = SphP[i].Density * All.UnitDensity_in_cgs * All.cf_a3inv * (All.HubbleParam * All.HubbleParam);
              double M_strom   = 6.2025e+10 / phys_dens / All.UnitMass_in_g * All.HubbleParam;
              if(P[i].Mass > 2.0 * M_strom)
                return 1;
            }
        }
#endif
    }
  return 0;
}

int refine_criterion_windtunnel(int i)
{
#ifdef WINDTUNNEL
  if(P[i].Type == PTYPE_GAS)
    {
#ifdef WINDTUNNEL_REFINEMENT_VOLUME_LIMIT

      if(SphP[i].Volume > 2. * All.MaxVolume)
        if(can_this_cell_be_split(i))
          return 1;

      if(SphP[i].Volume < 2. * All.MinVolume)
        return 0;

      if(P[i].Pos[WINDTUNNEL_COORD] < All.InjectionRegion)  // injection region
        {
          // if(SphP[i].Volume > 1.5 * All.InjectionVolume)//|| P[i].Mass > 2.0 * All.TargetGasMass
          //    if(can_this_cell_be_split(i))
          //        return 1;
        }
      else if(P[i].Pos[WINDTUNNEL_COORD] > All.InjectionRegion &&
              All.BoxSizes[WINDTUNNEL_COORD] - P[i].Pos[WINDTUNNEL_COORD] > All.InjectionRegion)  // the "normal region"
        {
          if(refine_criterion_volume(i))
            if(can_this_cell_be_split(i))
              return 1;
        }
      else if(All.BoxSizes[WINDTUNNEL_COORD] - P[i].Pos[WINDTUNNEL_COORD] < All.InjectionRegion)  // outflow region
        {
          // if (SphP[i].Volume > 2.0 * All.InjectionVolume)
          //    if(can_this_cell_be_split(i))
          //        return 1;
        }
#endif

      if(can_this_cell_be_split(i))
        {
          if((P[i].Pos[WINDTUNNEL_COORD] < All.InjectionRegion && SphP[i].Volume > 1.5 * All.InjectionVolume) ||
             (SphP[i].Volume > 2.0 * All.InjectionVolume))
            return 1;

          if(P[i].Mass > 2.0 * All.TargetGasMass)
            return 1;
        }
    }

#endif

  return 0; /* default is not to refine */
}

int refine_criterion_special_boundary(int i)
{
#ifdef SPECIAL_BOUNDARY
  if(SphP[i].MinDistBoundaryCell < 4 * All.BoundaryLayerScaleFactor && P[i].ID != -1 && P[i].ID != -2)
    return 0;
    // else
    //  return refine_criterion_default(i);
#endif
  return 0;
}

#ifdef REFINEMENT_VOLUME_LIMIT
int refine_criterion_volume(int i)
{
  if(All.MaxVolumeDiff > 0 && SphP[i].Volume > All.MaxVolumeDiff * SphP[i].MinNgbVolume)
    {
#ifdef REGULARIZE_MESH_FACE_ANGLE
      if(SphP[i].MaxFaceAngle < 1.5 * All.CellMaxAngleFactor)
        return 1;
#else
      double dx   = nearest_x(P[i].Pos[0] - SphP[i].Center[0]);
      double dy   = nearest_y(P[i].Pos[1] - SphP[i].Center[1]);
      double dz   = nearest_z(P[i].Pos[2] - SphP[i].Center[2]);

      double d = sqrt(dx * dx + dy * dy + dz * dz);

      double cellrad = get_cell_radius(i);

      if(d < 2.0 * All.CellShapingFactor * cellrad) /* only refine cells which are reasonably 'round' */
        return 1;
#endif
    }

  return 0;
}
#endif

#ifdef ROTATING_HIGHRES_REGION
int refine_rotating_highres_region(int i)
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

#ifdef SINK_PARTICLES_VARIABLE_CREATION
  SphP[i].RefineTarget = zoomtarget;
#endif

  if(can_this_cell_be_split(i) && P[i].Mass > 2.0 * zoomtarget)
    return 1;

  return 0; /* default is not to refine */
}
#endif

#ifdef SPIRAL
int refine_disczoom(int i)
{
  /* This routine refines a segment of a galactic disc that moves with the
     mean gas. By Rowan Smith 2013. */
  double zoomtarget;
  double wr, thref, dtheta, theta_range, theta_range2, zrange, rinner, Rref, normalise;
  float rfloat, sp_ang;
  int revolve;

  Rref    = 7.5 * KILOPARSEC;
  wr      = -220.0 * 1.e5 / Rref; /* vr of 220 kms at 7.5 kpc */
  thref   = wr * All.Time * All.UnitTime_in_s;
  revolve = (int)(fabs(thref) / (2.0 * M_PI));
  rfloat  = (float)revolve;
  thref   = thref + rfloat * 2.0 * M_PI; /*Note this assumes the rotation is clockwise */

  theta_range  = M_PI / 4.;
  theta_range2 = theta_range / 2.;
  zrange       = 1. * KILOPARSEC / All.UnitLength_in_cm; /*Only refine 1kpc from disc plane*/
  rinner       = 4. * KILOPARSEC / All.UnitLength_in_cm; /*Inner radius set to 4 kpc*/
  normalise    = 2.3e-13;                                /*2.3e-13 gives 5 solar mass resolution*/

  double dx = P[i].Pos[0] - boxHalf_X;
  double dy = P[i].Pos[1] - boxHalf_Y;
  double dz = P[i].Pos[2] - boxHalf_Z;
  double rr = sqrt(dx * dx + dy * dy);

  zoomtarget = All.TargetGasMass;

  if(dz < zrange && dz > -1.0 * zrange && rr > rinner)
    {
      /*  mpi_printf("Test values boxhalf x %g & z %g zrange %g, rinner %g, rr %g, dx %g, dz
       * %g\n",boxHalf_X,boxHalf_Z,zrange,rinner,rr,dx,dz);*/

      dtheta = atan2(dy, dx) - thref;

      if(dtheta > 2.0 * M_PI)
        dtheta = dtheta - 2.0 * M_PI;

      if(dtheta > M_PI)
        dtheta = 2. * M_PI - dtheta;

      /*      if(dtheta < -1.0 * M_PI)
                {
                  dtheta = 2.0 * M_PI + dtheta;
                }*/
      sp_ang = fabs(dtheta);
#ifdef RAMP_REFINE
      if(sp_ang < theta_range && sp_ang > theta_range2)
        {
          zoomtarget = All.TargetGasMass / 1000.0 * (2490.0 * (sp_ang - M_PI / 8.0) + 5.0);
#ifdef SNE_RAMP_REFINE
          zoomtarget = All.TargetGasMass / 1000.0 * 2.0 * (1245.0 * (sp_ang - M_PI / 8.0) + 5.0);
#endif
        }
      if(zoomtarget < 100.0 * SOLAR_MASS / All.UnitMass_in_g) /* to avoid sharp contrasts use exponential */
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

#ifdef SINK_PARTICLES_VARIABLE_CREATION
  SphP[i].RefineTarget = zoomtarget;
#endif

  if(can_this_cell_be_split(i) && P[i].Mass > 2.0 * zoomtarget)
    return 1;

  return 0; /* default is not to refine */
}

#endif

#ifdef DISC_REFINE_ONLY
int refine_disconly(int i)
{
  double zoomtarget;
  double dx = P[i].Pos[0] - boxHalf_X;
  double dy = P[i].Pos[1] - boxHalf_Y;
  double dz = P[i].Pos[2] - boxHalf_Z;
  double rr = sqrt(dx * dx + dy * dy);

  double zrange = KILOPARSEC / All.UnitLength_in_cm; /*Only refine 1kpc from plane*/
  double rinner = 1. * KILOPARSEC / All.UnitLength_in_cm;
  double router = 20. * KILOPARSEC / All.UnitLength_in_cm;

  zoomtarget = All.TargetGasMass;
  if(fabs(dz) < zrange && rr > rinner && rr < router)
    {
      zoomtarget = All.TargetGasMass / 10.;
    }

  if(can_this_cell_be_split(i) && P[i].Mass > 2.0 * zoomtarget)
    return 1;
  return 0; /* default is not to refine */
}
#endif

#ifdef REFINEMENT_AROUND_DM
void dm_particle_create_list(void)
{
  struct refine_dm_data *DMPartList;
  DMPartList = (struct refine_dm_data *)mymalloc("DMPartList", All.TotPartDM * sizeof(struct refine_dm_data));

  int i, j, nsrc, nimport, ngrp;
  for(i = 0, nsrc = 0; i < NumPart; i++)
    {
      if(P[i].Type == 1)
        {
          DMPartList[nsrc].ID = P[i].ID;

          DMPartList[nsrc].pos[0] = P[i].Pos[0];
          DMPartList[nsrc].pos[1] = P[i].Pos[1];
          DMPartList[nsrc].pos[2] = P[i].Pos[2];

          // DMPartList[nsrc++].softening = All.SofteningTable[P[i].SofteningType];
          DMPartList[nsrc++].softening = All.ForceSoftening[P[i].SofteningType];
        }
    }

  for(j = 0; j < NTask; j++)
    Send_count[j] = nsrc;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = 0;
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  /* exchange particle data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&DMPartList[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct refine_dm_data), MPI_BYTE,
                           recvTask, TAG_DENS_A, &DMPartListGlobal[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct refine_dm_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(DMPartList);

#ifdef REFINEMENT_VOLUME_LIMIT
  /* now set minimum volume to refine to */
  /* get minimum softening length of all DM particles */
  double minrad = MAX_REAL_NUMBER;
  for(j = 0; j < All.TotPartDM; j++)
    if(minrad > DMPartListGlobal[j].softening / All.RefinementCellsPerSoftening)
      minrad = DMPartListGlobal[j].softening / All.RefinementCellsPerSoftening;
  /* set volume to 1% of corresponding sphere volume */
  All.MinVolume = 0.01 * 4. * M_PI / 3. * minrad * minrad * minrad;
  mpi_printf("Setting minimum volume to %e\n", All.MinVolume);
#endif
}

void dm_particle_update_list(void)
{
  struct refine_dm_data *DMPartList;
  DMPartList = (struct refine_dm_data *)mymalloc("DMPartList", All.TotPartDM * sizeof(struct refine_dm_data));

  int i, j, nsrc, nimport, ngrp;
  for(i = 0, nsrc = 0; i < NumPart; i++)
    {
      if(P[i].Type == 1)
        {
          DMPartList[nsrc].ID = P[i].ID;

          DMPartList[nsrc].pos[0] = P[i].Pos[0];
          DMPartList[nsrc].pos[1] = P[i].Pos[1];
          DMPartList[nsrc].pos[2] = P[i].Pos[2];

          // DMPartList[nsrc++].softening = All.SofteningTable[P[i].SofteningType];
          DMPartList[nsrc++].softening = All.ForceSoftening[P[i].SofteningType];
        }
    }

  for(j = 0; j < NTask; j++)
    Send_count[j] = nsrc;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = 0;
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  /* exchange particle data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&DMPartList[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct refine_dm_data), MPI_BYTE,
                           recvTask, TAG_DENS_A, &DMPartListGlobal[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct refine_dm_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(DMPartList);
}
#endif

#ifdef GMC_REFINEMENT
int gmc_refinement_criteria(int i)
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
  if((density < All.GMCRefMinDensity) || (density > All.GMCRefMaxDensity))
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
  if(cells_per_jeans_length < All.GMCRefCellsPerJeansLength)
    return (1);

  return (0);
}
#endif

#endif
