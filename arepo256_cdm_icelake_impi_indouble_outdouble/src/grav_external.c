/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/grav_external.c
 * \date        MM/YYYY
 * \author
 * \brief       special gravity routines for external forces
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "allvars.h"
#include "domain.h"
#include "proto.h"

#ifdef EXTERNALGRAVITY
static void gravity_external_get_force(double pos[3], int type, MyIDType ID, double acc[3], double *pot, int *flag_set);

/*! \brief Main routine to add contribution of external gravitational potential
 *  to accelerations.
 *
 *  Function is called in gravity() (in accel.c). Function also evaluates
 *  the gradient of the accelerations which is needed for the timestep
 *  criterion due to the external potential.
 *
 *  \return void
 */
void gravity_external(void)
{
  mpi_printf("EXTERNALGRAVITY: execute\n");

  TIMER_START(CPU_TREE);

  for(int idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      int i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double *pos;

#ifdef CELL_CENTER_GRAVITY
      if(P[i].Type == 0)
        pos = SphP[i].Center;
      else
#endif
        pos = P[i].Pos;

      double acc[3], pot;
      int flag_set = 0;
      gravity_external_get_force(pos, P[i].Type, P[i].ID, acc, &pot, &flag_set);

      if(flag_set)
        {
          for(int k = 0; k < NUMDIMS; k++)
            P[i].GravAccel[k] = acc[k];
          for(int k = NUMDIMS; k < 3; k++)
            P[i].GravAccel[k] = 0;
#ifndef SPIRAL
          P[i].ExtPotential = pot;
#endif
        }
      else
        {
          for(int k = 0; k < NUMDIMS; k++)
            P[i].GravAccel[k] += acc[k];
#ifdef EVALPOTENTIAL
          P[i].Potential += pot;
#endif
#ifndef SPIRAL
          P[i].ExtPotential += pot;
#endif
        }

      double dx;
      if(P[i].Type == 0)
        dx = 0.1 * get_cell_radius(i);
      else
        dx = 0.1 * All.ForceSoftening[P[i].SofteningType];

      P[i].dGravAccel = 0;
      for(int dim = 0; dim < NUMDIMS; dim++)
        {
          double accL[3], posL[3];
          for(int k = 0; k < 3; k++)
            posL[k] = pos[k];
          posL[dim] -= dx;
          gravity_external_get_force(posL, P[i].Type, P[i].ID, accL, &pot, &flag_set);

          double accR[3], posR[3];
          for(int k = 0; k < 3; k++)
            posR[k] = pos[k];
          posR[dim] += dx;
          gravity_external_get_force(posR, P[i].Type, P[i].ID, accR, &pot, &flag_set);

          for(int k = 0; k < NUMDIMS; k++)
            {
              double dGrav = accR[k] - accL[k];
              P[i].dGravAccel += dGrav * dGrav;
            }
        }
      P[i].dGravAccel = sqrt(P[i].dGravAccel) / (2. * dx);
    }

  TIMER_STOP(CPU_TREE);
}

static void gravity_external_get_force(double pos[3], int type, MyIDType ID, double acc[3], double *pot, int *flag_set)
{
  for(int k = 0; k < 3; k++)
    acc[k] = 0;

#ifdef THERMAL_INSTABILITY_TEST
  acc[1]    = 1.0;
  *flag_set = 1;
#endif

#ifdef MRT
#ifdef MRT_LEVITATION_TEST
  acc[1]    = -1.46e-6 * pow(All.UnitTime_in_s, 2.0) / All.UnitLength_in_cm;
  *flag_set = 1;
#endif
#endif

#ifdef DG_TEST_PROBLEM
  dg_acceleration(pos[0], pos[1], pos[2], acc);
  *flag_set = 1;
#endif

  *pot = 0;

#ifdef EXTERNALGY
  acc[1] += EXTERNALGY;
  *pot = -(EXTERNALGY)*pos[1];
#endif

#ifdef EXTERNALDISKPOTENTIAL
  double dx, dy, r;

  dx = pos[0] - boxHalf_X;
  dy = pos[1] - boxHalf_Y;
  r  = sqrt(dx * dx + dy * dy);

  double Sigma0 = 1.0 / (2 * M_PI);
  double y      = r / (2);

  double dphidR = externaldisk_dphidR(r, NULL);
  double pot    = externaldisk_potential(r);

  acc[0] -= dphidR * dx / r;
  acc[1] -= dphidR * dy / r;
  *pot = pot;
#endif

#ifdef EXTERNALSHEARBOX
  double dz, b;

  dz = pos[2] - boxHalf_Z;

  double Sigma0 = All.ShearBoxSigma0;
  double fg     = All.ShearBoxFg;
  double mu     = All.ShearBoxMu;

  b = 61.0 * (fg / 0.1 / mu) / (Sigma0 / 10.0);
  double f;

#ifdef SELFGRAVITY
  f = 1.0 / (1.0 / fg - 1.0);
#else
  f = fg;
#endif

  acc[2] -= 2.0 * M_PI * All.G * Sigma0 * tanh(dz / b) / f;
  *pot = 2.0 * M_PI * All.G * b * Sigma0 * log(cosh(dz / b)) / f;
#endif

#ifdef GRAVITY_TABLE
  double xi = pos[0] - boxHalf_X;
  double yi = pos[1] - boxHalf_Y;
  double zi = pos[2] - boxHalf_Z;

  sp_acc[0] = 0.0;
  sp_acc[1] = 0.0;
  sp_acc[2] = 0.0;

  /* contributions from a multi component Galaxy model defined in an external file */
  grav_table_find_grav_acceleration(xi, yi, zi, sp_acc);
#ifdef SPIRAL
  /* Now the contributions from the time-dependent spiral
     given by Cox & Gomez (2002) */
  double hi;
  if(NumPart != NumGas && i >= NumPart - 1)
    hi = 1. * PARSEC / All.UnitLength_in_cm;  // collisionless particles have no volume field in their structure so just set it to a
                                              // small enough number
  else
    hi = get_cell_radius(i);

  double dhi = hi / 1000.;

  double potent1 = spiral_potential_calc(xi, yi, zi, All.Time);

  double potent2 = spiral_potential_calc(xi + dhi, yi, zi, All.Time);
  sp_acc[0] += -(potent2 - potent1) / dhi;

  potent2 = spiral_potential_calc(xi, yi + dhi, zi, All.Time);
  sp_acc[1] += -(potent2 - potent1) / dhi;

  potent2 = spiral_potential_calc(xi, yi, zi + dhi, All.Time);
  sp_acc[2] += -(potent2 - potent1) / dhi;
#endif
  P[i].GravAccel[0] = sp_acc[0];
  P[i].GravAccel[1] = sp_acc[1];
  P[i].GravAccel[2] = sp_acc[2];
#endif

#ifdef GALPOT
  double xi = pos[0] - boxHalf_X;
  double yi = pos[1] - boxHalf_Y;
  double zi = pos[2] - boxHalf_Z;

  double dPhidx[3];
  galpot_dPhidx(xi, yi, zi, All.Time, &dPhidx[0]);

  acc[0] -= dPhidx[0];
  acc[1] -= dPhidx[1];
  acc[2] -= dPhidx[2];

  double R = sqrt(xi * xi + yi * yi + zi * zi);

  if(R == 0)
    {
      acc[0] = 0.0;
      acc[1] = 0.0;
      acc[2] = 0.0;
    }
#endif

#ifdef EXTERNALSHEETY
  if(type == 0)
    acc[1] = (-2 * M_PI * tanh((pos[1] - 1.0) / 0.1));
#endif

#ifdef STATICISO
#if !defined(ISO_Eps) || !defined(ISO_M200) || !defined(ISO_R200)
#error "STATICISO requires NFW_Eps, ISO_M200 and ISO_R200 to be defined"
#endif
  {
    double r, m;
    double dx, dy, dz;

#ifdef AXISYMMETRY
    dx = pos[0];
    dy = pos[1] - boxHalf_Y;
    dz = 0;
#else
    dx = pos[0] - boxHalf_X;
    dy = pos[1] - boxHalf_Y;
    dz = pos[2] - boxHalf_Z;
#endif
    r = sqrt(dx * dx + dy * dy + dz * dz);

    if(r > ISO_R200)
      m = ISO_M200;
    else
      m = ISO_M200 * r / ISO_R200;

#ifdef ISO_FRACTION
    m *= ISO_FRACTION;
#endif

    if(r > 0)
      {
        acc[0] += -All.G * m * dx / r / (r * r + ISO_Eps * ISO_Eps);
        acc[1] += -All.G * m * dy / r / (r * r + ISO_Eps * ISO_Eps);
        acc[2] += -All.G * m * dz / r / (r * r + ISO_Eps * ISO_Eps);
      }
  }
#endif

#ifdef GROWING_DISK_POTENTIAL
  {
    double mdisk, dx, dy, dz, r, z, aR, az;

    growing_disk_init();

    mdisk = get_disk_mass(All.Time);

    dx = pos[0] - boxHalf_X;
    dy = pos[1] - boxHalf_Y;
    dz = pos[2] - boxHalf_Z;

    r = sqrt(dx * dx + dy * dy);
    z = fabs(dz);

    get_disk_forces(r, z, &aR, &az);

    aR *= mdisk;
    az *= mdisk;

    if(r > 0)
      {
        acc[0] += -dx / r * aR;
        acc[1] += -dy / r * aR;
        acc[2] += -dz / z * az;
      }
  }
#endif

#ifdef STATICNFW
#if !defined(NFW_C) || !defined(NFW_Eps) || !defined(NFW_M200)
#error "STATICNFW requires NFW_C, NFW_Eps and NFW_M200 to be defined"
#endif
  {
    double r, m;
    double dx, dy, dz;

    dx = pos[0] - boxHalf_X;
    dy = pos[1] - boxHalf_Y;
    dz = pos[2] - boxHalf_Z;

    r = sqrt(dx * dx + dy * dy + dz * dz);
    m = enclosed_mass(r);
#ifdef NFW_DARKFRACTION
    m *= NFW_DARKFRACTION;
#endif
    if(r > 0)
      {
        acc[0] += -All.G * m * dx / (r * r * r);
        acc[1] += -All.G * m * dy / (r * r * r);
        acc[2] += -All.G * m * dz / (r * r * r);
      }
#if defined(VS_TURB) || defined(AB_TURB)
    *pot += get_turb_pot(pos[0], pos[1], pos[2]);
#endif
  }
#endif

#ifdef STATIC_ISOTHERMAL_CLUSTER
  {
    /*
    A simple potential for an isothermal intracluster medium in a galaxy cluster.
    See details in the example rising_bubble_2d where this profile is used.
    */
    double r;
    double dx, dy, dz;
    double r0 = 1.0, cs = 1.0;

    dx = pos[0] - boxHalf_X;
    dy = pos[1] - boxHalf_Y;
    dz = pos[2] - boxHalf_Z;

#ifdef TWODIMS
    dz = 0.0;
#endif

    double g_acc;
    r     = sqrt(dx * dx + dy * dy + dz * dz);
    g_acc = cs * cs * 3.0 / 2.0 * r / (r * r + r0 * r0);
    *pot  = cs * cs * 3.0 / 4.0 * log(r * r + r0 * r0);
    acc[0] += -g_acc * dx / r;
    acc[1] += -g_acc * dy / r;
    acc[2] += -g_acc * dz / r;
  }
#endif

#ifdef STATICHQ
#if !defined(HQ_C) || !defined(HQ_M200)
#error "STATICHQ requires HQ_C and HQ_M200 to be defined"
#endif
  {
    double r, m, a;
    double dx, dy, dz;

    dx = pos[0] - boxHalf_X;
    dy = pos[1] - boxHalf_Y;
    dz = pos[2] - boxHalf_Z;

    r = sqrt(dx * dx + dy * dy + dz * dz);

#ifdef HQ_A
    a = HQ_A;
#else
    a  = pow(All.G * HQ_M200 / (100 * All.Hubble * All.Hubble), 1.0 / 3) / HQ_C * sqrt(2 * (log(1 + HQ_C) - HQ_C / (1 + HQ_C)));
#endif

    m = HQ_M200 * pow(r / (r + a), 2);
#ifdef HQ_DARKFRACTION
    m *= HQ_DARKFRACTION;
#endif
    if(r > 0)
      {
        acc[0] += -All.G * m * dx / (r * r * r);
        acc[1] += -All.G * m * dy / (r * r * r);
        acc[2] += -All.G * m * dz / (r * r * r);
      }
  }
#endif

#ifdef CONSTANT_GRAVITY
  acc[1] += -1.0;
#endif

#ifdef CENTRAL_MASS_POTENTIAL
  double fac, wp;
  double dx, dy, r, r2;
  double dz = 0;
  double h, h_inv, h3_inv, u;

  h = All.SofteningCentral * 2.8;

  dx = pos[0] - boxHalf_X;
  dy = pos[1] - boxHalf_Y;
#ifndef TWODIMS
  dz = pos[2] - boxHalf_Z;
#endif

  r2 = dx * dx + dy * dy + dz * dz;
  r  = sqrt(r2);

  // using spline softening
  if(r >= h)
    {
      fac = 1 / (r2 * r);
      wp  = -1 / r;
    }
  else
    {
      h_inv  = 1.0 / h;
      h3_inv = h_inv * h_inv * h_inv;
      u      = r * h_inv;

      if(u < 0.5)
        {
          fac = h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
          wp  = h_inv * (-2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6)));
        }
      else
        {
          fac = h3_inv * (21.333333333333 - 48.0 * u + 38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
          wp  = h_inv * (-3.2 + 0.066666666667 / u + u * u * (10.666666666667 + u * (-16.0 + u * (9.6 - 2.133333333333 * u))));
        }
    }

  acc[0] -= All.G * All.CentralMass * fac * dx;
  acc[1] -= All.G * All.CentralMass * fac * dy;
#ifndef TWODIMS
  acc[2] -= All.G * All.CentralMass * fac * dz;
#endif

  *pot = All.G * All.CentralMass * wp;

#ifdef SPECIAL_BOUNDARY
  if(ID <= -3)
    {
      acc[0]    = 0.0;
      acc[1]    = 0.0;
      acc[2]    = 0.0;
      *flag_set = 1;
    }
#endif
#endif

#ifdef STAR_PLANET_POTENTIAL
  double dx, dy, r;
  double first_mass, second_mass, dx_2, dy_2, r_2, soft;
  double indirect;

  soft        = All.SofteningPlanet * 2.8;
  first_mass  = (1 - All.MassRatio);
  second_mass = All.MassRatio;

  if(All.Time <= All.PlanetGrowthTime)
    second_mass *= sin(0.5 * M_PI * All.Time / All.PlanetGrowthTime) * sin(0.5 * M_PI * All.Time / All.PlanetGrowthTime);

  dx       = pos[0] - boxHalf_X;
  dy       = pos[1] - boxHalf_Y;
  dx_2     = pos[0] - boxHalf_X - cos(All.Time);
  dy_2     = pos[1] - boxHalf_Y - sin(All.Time);
  indirect = All.G * second_mass * (dx * cos(All.Time) + dy * sin(All.Time));

  r   = sqrt(dx * dx + dy * dy);
  r_2 = sqrt(dx_2 * dx_2 + dy_2 * dy_2 + soft * soft);

  // or using spline softening
  r_2 = sqrt(dx_2 * dx_2 + dy_2 * dy_2);
  double u_2, fac_2, wp_2, soft_inv, soft3_inv;
  if(r_2 >= soft)
    {
      fac_2 = 1 / (r_2 * r_2 * r_2);
      wp_2  = -1 / r_2;
    }
  else
    {
      soft_inv  = 1.0 / soft;
      soft3_inv = soft_inv * soft_inv * soft_inv;
      u_2       = r_2 * soft_inv;

      if(u_2 < 0.5)
        {
          fac_2 = soft3_inv * (10.666666666667 + u_2 * u_2 * (32.0 * u_2 - 38.4));
          wp_2  = soft_inv * (-2.8 + u_2 * u_2 * (5.333333333333 + u_2 * u_2 * (6.4 * u_2 - 9.6)));
        }
      else
        {
          fac_2 = soft3_inv * (21.333333333333 - 48.0 * u_2 + 38.4 * u_2 * u_2 - 10.666666666667 * u_2 * u_2 * u_2 -
                               0.066666666667 / (u_2 * u_2 * u_2));
          wp_2  = soft_inv *
                 (-3.2 + 0.066666666667 / u_2 + u_2 * u_2 * (10.666666666667 + u_2 * (-16.0 + u_2 * (9.6 - 2.133333333333 * u_2))));
        }
    }

  *pot = -All.G * first_mass / r + All.G * second_mass * wp_2 + indirect;

  acc[0] -= All.G * first_mass / r / r / r * dx + All.G * second_mass * fac_2 * dx_2 + All.G * second_mass * cos(All.Time);
  acc[1] -= All.G * first_mass / r / r / r * dy + All.G * second_mass * fac_2 * dy_2 + All.G * second_mass * sin(All.Time);

#ifdef SPECIAL_BOUNDARY
  if(ID <= -3)
    {
      acc[0]    = 0.0;
      acc[1]    = 0.0;
      acc[2]    = 0.0;
      *flag_set = 1;
    }
#endif
#endif

#ifdef BINARY_POTENTIAL
  double dx_1, dy_1, dr_1, dx_2, dy_2, dr_2;
  double first_mass, second_mass;

  double soft;
  double soft1 = All.BinarySoftening * 2.8;
  double soft2 = soft1;

  double indirect = 0.0;

  double mu   = All.BinaryMassRatio / (1.0 + All.BinaryMassRatio);
  first_mass  = (1.0 - mu);
  second_mass = mu;
  if(All.Time < All.BinaryGrowthTime)
    second_mass *= sin(0.5 * M_PI * All.Time / All.BinaryGrowthTime) * sin(0.5 * M_PI * All.Time / All.BinaryGrowthTime);

  if(ID == 5)
    second_mass = 0;

  double x, y, vx, vy;
  circumstellar_solve_kepler(All.Time + M_PI, All.BinaryEccentricity, &x, &y, &vx, &vy);

  double x1, x2, y1, y2;
  if(All.BinaryBarycentricCoord)
    {
      x1 = boxHalf_X + mu * x;
      y1 = boxHalf_Y + mu * y;

      x2 = boxHalf_X - (1.0 - mu) * x;
      y2 = boxHalf_Y - (1.0 - mu) * y;
    }
  else
    {
      x1 = boxHalf_X;
      y1 = boxHalf_Y;

      x2 = boxHalf_X - x;
      y2 = boxHalf_Y - y;
    }
#ifdef CIRCUMSTELLAR_WBOUNDARIES
  if(All.BinaryBarycentricCoord == 0)
    soft1 = 0.0;
#endif

  dx_1       = pos[0] - x1;
  dy_1       = pos[1] - y1;
  dx_2       = pos[0] - x2;
  dy_2       = pos[1] - y2;
  double r_2 = sqrt(x2 * x2 + y2 * y2);

  if(All.BinaryBarycentricCoord)
    indirect = 0;
  else
    indirect = All.G * second_mass * (dx_1 * x2 + dy_1 * y2) / r_2 / r_2 / r_2;

  dr_1 = sqrt(dx_1 * dx_1 + dy_1 * dy_1);
  dr_2 = sqrt(dx_2 * dx_2 + dy_2 * dy_2);

  int k;
  double r, u, fac, wp, soft_inv, soft3_inv;
  double wp_1, wp_2, fac_1, fac_2;
  for(k = 0; k < 2; k++)
    {
      if(k == 0)
        {
          r    = dr_1;
          soft = soft1;
        }
      else
        {
          r    = dr_2;
          soft = soft2;
        }

      if(r >= soft)
        {
          fac = 1 / (r * r * r);
          wp  = -1 / r;
        }
      else
        {
          soft_inv  = 1.0 / soft;
          soft3_inv = soft_inv * soft_inv * soft_inv;
          u         = r * soft_inv;

          if(u < 0.5)
            {
              fac = soft3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
              wp  = soft_inv * (-2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6)));
            }
          else
            {
              fac =
                  soft3_inv * (21.333333333333 - 48.0 * u + 38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
              wp = soft_inv * (-3.2 + 0.066666666667 / u + u * u * (10.666666666667 + u * (-16.0 + u * (9.6 - 2.133333333333 * u))));
            }
        }

      if(k == 0)
        {
          wp_1  = wp;
          fac_1 = fac;
        }
      else
        {
          wp_2  = wp;
          fac_2 = fac;
        }
    }

  *pot = -All.G * first_mass * wp_1 - All.G * second_mass * wp_2 + indirect;

  acc[0] -= All.G * first_mass * fac_1 * dx_1 + All.G * second_mass * fac_2 * dx_2;
  acc[1] -= All.G * first_mass * fac_1 * dy_1 + All.G * second_mass * fac_2 * dy_2;

  if(All.BinaryBarycentricCoord == 0)
    {
      acc[0] -= All.G * second_mass * x2 / r_2 / r_2 / r_2;
      acc[1] -= All.G * second_mass * y2 / r_2 / r_2 / r_2;
    }

#ifdef SPECIAL_BOUNDARY
  if(ID <= -3)
    {
      acc[0]    = 0.0;
      acc[1]    = 0.0;
      acc[2]    = 0.0;
      *flag_set = 1;
    }
#endif
#endif

#if defined(CIRCUMSTELLAR) && defined(GRAVITY_FROM_STARS_PLANETS_ONLY)
  circumstellar_calc_gravity_from_stars_planets_only();
#endif
}
#endif

#ifdef ONEDIMS_SPHERICAL
/*! \brief One-dimensional gravity in the spherically symmetric case.
 *
 *  \return void
 */
void gravity_monopole_1d_spherical(void)
{
  printf("Doing 1D gravity...\n");

  int i;
  double msum = All.CoreMass;

  for(i = 0; i < NumGas; i++)
    {
      double dm  = 0.5 * P[i].Mass;
      double rad = SphP[i].Center[0];

      P[i].GravAccel[0] = -(msum + dm) * All.G / (rad * rad);

#ifdef EVALPOTENTIAL
      P[i].Potential = -(msum + dm) * All.G / rad;
#endif

      msum += P[i].Mass;

      P[i].GravAccel[1] = 0;
      P[i].GravAccel[2] = 0;
    }

  printf("... 1D gravity done.\n");
}
#endif

#ifdef STATICNFW
/*! \brief Auxiliary function for static NFW potential.
 *
 *  \param[in] R Radius from center of potential.
 *
 *  \return Enclosed mass (which causes the external potential).
 */
double enclosed_mass(double R)
{
  /* Eps is in units of Rs !!!! */

  if(R > Rs * NFW_C)
    R = Rs * NFW_C;

  return fac * 4 * M_PI * RhoCrit * Dc *
         (-(Rs * Rs * Rs * (1 - NFW_Eps + log(Rs) - 2 * NFW_Eps * log(Rs) + NFW_Eps * NFW_Eps * log(NFW_Eps * Rs))) /
              ((NFW_Eps - 1) * (NFW_Eps - 1)) +
          (Rs * Rs * Rs *
           (Rs - NFW_Eps * Rs - (2 * NFW_Eps - 1) * (R + Rs) * log(R + Rs) + NFW_Eps * NFW_Eps * (R + Rs) * log(R + NFW_Eps * Rs))) /
              ((NFW_Eps - 1) * (NFW_Eps - 1) * (R + Rs)));
}

#ifdef DG_EXTERNAL_ACCELERATION
void dg_acceleration(double x, double y, double z, double *acc)
{
  double r, M, a, dx, dy, dz;
  dx = x - boxHalf_X;
  dy = y - boxHalf_Y;
  dz = z - boxHalf_Z;
  r  = sqrt(dx * dx + dy * dy + dz * dz);
  M  = enclosed_mass(r);
  a  = -All.G * M / r / r / r;

  acc[0] = a * dx;
  acc[1] = a * dy;
  acc[2] = a * dz;
}
#endif
#endif

#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
/*! \brief Routine that computes gravitational force by direct summation.
 *
 *  Called by gravity() (in accel.c).
 *
 *  \return void
 */
void calc_exact_gravity_for_particle_type(void)
{
  int i, idx;
#ifdef EXACT_GRAVITY_REACTION
  double *accx, *accy, *accz;
  accx = (double *)mymalloc("accx", All.MaxPartSpecial * sizeof(double));
  accy = (double *)mymalloc("accy", All.MaxPartSpecial * sizeof(double));
  accz = (double *)mymalloc("accz", All.MaxPartSpecial * sizeof(double));
#ifdef EVALPOTENTIAL
  double *pot;
  pot = (double *)mymalloc("pot", All.MaxPartSpecial * sizeof(double));
#endif
  int n;
  for(n = 0; n < All.MaxPartSpecial; n++)
    {
      accx[n] = accy[n] = accz[n] = 0.0;
#ifdef EVALPOTENTIAL
      pot[n] = 0.0;
#endif
    }
#endif

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double fac, wp;
      double dx, dy, dz, r, r2;
      double h, h_inv, h3_inv, u;
      int k;

      /* set softening to corresponding particle's softening length */
      h = All.ForceSoftening[All.SofteningTypeOfPartType[EXACT_GRAVITY_FOR_PARTICLE_TYPE]];

      for(k = 0; k < All.MaxPartSpecial; k++)
        {
          if(PartSpecialListGlobal[k].ID == P[i].ID)
            continue;
#ifdef REFINEMENT_AROUND_DM
          /* reset softening length to the length for this particle */
          h = DMPartListGlobal[k].softening;
#endif

          dx = P[i].Pos[0] - PartSpecialListGlobal[k].pos[0];
          dy = P[i].Pos[1] - PartSpecialListGlobal[k].pos[1];
          dz = P[i].Pos[2] - PartSpecialListGlobal[k].pos[2];

          r2 = dx * dx + dy * dy + dz * dz;
          r  = sqrt(r2);

          // using spline softening
          if(r >= h)
            {
              fac = 1 / (r2 * r);
              wp  = -1 / r;
            }
          else
            {
              h_inv  = 1.0 / h;
              h3_inv = h_inv * h_inv * h_inv;
              u      = r * h_inv;

              if(u < 0.5)
                {
                  fac = h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
                  wp  = h_inv * (-2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6)));
                }
              else
                {
                  fac = h3_inv *
                        (21.333333333333 - 48.0 * u + 38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
                  wp = h_inv * (-3.2 + 0.066666666667 / u + u * u * (10.666666666667 + u * (-16.0 + u * (9.6 - 2.133333333333 * u))));
                }
            }

          P[i].GravAccel[0] -= All.G * PartSpecialListGlobal[k].mass * fac * dx;
          P[i].GravAccel[1] -= All.G * PartSpecialListGlobal[k].mass * fac * dy;
          P[i].GravAccel[2] -= All.G * PartSpecialListGlobal[k].mass * fac * dz;

#ifdef EVALPOTENTIAL
          P[i].Potential += All.G * PartSpecialListGlobal[k].mass * wp;
#endif
#ifdef EXACT_GRAVITY_REACTION
          /* avoid double counting */
          if(P[i].Type != EXACT_GRAVITY_FOR_PARTICLE_TYPE)
            {
              accx[k] += All.G * P[i].Mass * fac * dx;
              accy[k] += All.G * P[i].Mass * fac * dy;
              accz[k] += All.G * P[i].Mass * fac * dz;
#ifdef EVALPOTENTIAL
              pot[k] += All.G * P[i].Mass * wp;
#endif
            }
#endif
        }
    }
#ifdef EXACT_GRAVITY_REACTION
  double *buf = (double *)mymalloc("buf", All.MaxPartSpecial * sizeof(double));

  MPI_Allreduce(accx, buf, All.MaxPartSpecial, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for(n = 0; n < All.MaxPartSpecial; n++)
    accx[n] = buf[n];
  MPI_Allreduce(accy, buf, All.MaxPartSpecial, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for(n = 0; n < All.MaxPartSpecial; n++)
    accy[n] = buf[n];
  MPI_Allreduce(accz, buf, All.MaxPartSpecial, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for(n = 0; n < All.MaxPartSpecial; n++)
    accz[n] = buf[n];
#ifdef EVALPOTENTIAL
  MPI_Allreduce(pot, buf, All.MaxPartSpecial, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for(n = 0; n < All.MaxPartSpecial; n++)
    pot[n] = buf[n];
#endif
  myfree(buf);

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;
      for(n = 0; n < All.MaxPartSpecial; n++)
        {
          if(PartSpecialListGlobal[n].ID == P[i].ID)
            {
              P[i].GravAccel[0] += accx[n];
              P[i].GravAccel[1] += accy[n];
              P[i].GravAccel[2] += accz[n];
#ifdef EVALPOTENTIAL
              P[i].Potential += pot[n];
#endif
            }
        }
    }

#ifdef EVALPOTENTIAL
  myfree(pot);
#endif
  myfree(accz);
  myfree(accy);
  myfree(accx);
#endif
}
#endif

#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
/*! \brief Creates list of special particles, i.e. particles for which gravity
 *  is calculated by direct summation.
 *
 *  Called in begrund2() (begrun.c), i.e. only at startup of the simulation.
 *
 *  \return void
 */
void special_particle_create_list(void)
{
  struct special_particle_data *SpecialPartList;
  SpecialPartList =
      (struct special_particle_data *)mymalloc("SpecialPartList", All.MaxPartSpecial * sizeof(struct special_particle_data));

  int i, j, nsrc, nimport, ngrp;
  for(i = 0, nsrc = 0; i < NumPart; i++)
    {
      if(P[i].Type == EXACT_GRAVITY_FOR_PARTICLE_TYPE)
        {
          SpecialPartList[nsrc].ID = P[i].ID;

          SpecialPartList[nsrc].pos[0] = P[i].Pos[0];
          SpecialPartList[nsrc].pos[1] = P[i].Pos[1];
          SpecialPartList[nsrc].pos[2] = P[i].Pos[2];

          SpecialPartList[nsrc++].mass = P[i].Mass;
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
              MPI_Sendrecv(&SpecialPartList[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct special_particle_data),
                           MPI_BYTE, recvTask, TAG_DENS_A, &PartSpecialListGlobal[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct special_particle_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(SpecialPartList);
}

/*! \brief Updates list of special particles, i.e. particles for which gravity
 *  is calculated by direct summation.
 *
 *  Called in run() (run.c).
 *
 *  \return void
 */
void special_particle_update_list(void)
{
  struct special_particle_data *SpecialPartList;
  SpecialPartList =
      (struct special_particle_data *)mymalloc("SpecialPartList", All.MaxPartSpecial * sizeof(struct special_particle_data));

  int i, j, nsrc, nimport, ngrp;
  for(i = 0, nsrc = 0; i < NumPart; i++)
    {
      if(P[i].Type == EXACT_GRAVITY_FOR_PARTICLE_TYPE)
        {
          SpecialPartList[nsrc].ID = P[i].ID;

          SpecialPartList[nsrc].pos[0] = P[i].Pos[0];
          SpecialPartList[nsrc].pos[1] = P[i].Pos[1];
          SpecialPartList[nsrc].pos[2] = P[i].Pos[2];

          SpecialPartList[nsrc++].mass = P[i].Mass;
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
              MPI_Sendrecv(&SpecialPartList[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct special_particle_data),
                           MPI_BYTE, recvTask, TAG_DENS_A, &PartSpecialListGlobal[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct special_particle_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(SpecialPartList);
}
#endif

#if defined(ACCRETE_ONTO_CENTRAL_POTENTIAL) && defined(CENTRAL_MASS_POTENTIAL)
/*! \brief Routine that allows gas accretion onto a central potential.
 *
 *  \return void
 */
void accrete_onto_central_potential(void)
{
  int idx, i;
  double dx, dy, dz, r2, r, racc2, vr, new_racc, nma;
  double racc_global;
  double dm_local = 0.0;
  double dm_global;
  int count_local = 0;
  int count_global;

  mpi_printf("CENTRALPOTENTIAL: Accreting all cells within = %g\n", All.CentralAccretionRadius);
  racc2    = All.CentralAccretionRadius * All.CentralAccretionRadius;
  new_racc = 2 * sqrt(racc2);

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      dx = P[i].Pos[0] - boxHalf_X;
      dy = P[i].Pos[1] - boxHalf_Y;
      dz = P[i].Pos[2] - boxHalf_Z;

      r2 = dx * dx + dy * dy + dz * dz;

      // Accretion criterion based on density and position
      if(r2 < racc2 && P[i].Mass > 0.1 * All.TargetGasMass)
        {
          count_local += 1;

          P[i].Mass *= 0.5;
          dm_local += P[i].Mass;

#ifdef MHD
          double Emag = 0.5 * (SphP[i].B[0] * SphP[i].B[0] + SphP[i].B[1] * SphP[i].B[1] + SphP[i].B[2] * SphP[i].B[2]) *
                        SphP[i].Volume * All.cf_atime;
          SphP[i].Energy -= Emag;
#endif
          SphP[i].Energy *= 0.5;
#ifdef MHD
          SphP[i].Energy += Emag;
#endif

          SphP[i].Momentum[0] *= 0.5;
          SphP[i].Momentum[1] *= 0.5;
          SphP[i].Momentum[2] *= 0.5;

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
          SphP[i].Entropy *= 0.5;
#endif
#ifdef MAXSCALARS
          for(int s = 0; s < N_Scalar; s++)
            *(MyFloat *)(((char *)(&SphP[i])) + scalar_elements[s].offset_mass) *= 0.5;
#endif
        }
      else if(All.HighestOccupiedTimeBin == All.HighestActiveTimeBin && r2 < 4.0 * racc2)
        {
          r   = sqrt(r2);
          vr  = (dx * P[i].Vel[0] + dy * P[i].Vel[1] + dz * P[i].Vel[2]) / r;
          nma = vr / SphP[i].Csnd;
          if(nma > -2.0 && r < new_racc)
            new_racc = r;
        }
    }

  if(All.HighestOccupiedTimeBin == All.HighestActiveTimeBin)
    {
      MPI_Allreduce(&new_racc, &racc_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      if(racc_global > All.CentralAccretionRadius)
        {
          All.CentralAccretionRadius = racc_global;
          mpi_printf("CENTRALPOTENTIAL: Adjusting new accretion radius to %g\n", All.CentralAccretionRadius);
        }
    }

  // Sum accreted mass
  MPI_Allreduce(&count_local, &count_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&dm_local, &dm_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  All.CentralMass += dm_global;

  mpi_printf("CENTRALPOTENTIAL: Accreted %g from %i cells onto central mass potential, total = %g\n", dm_global, count_global,
             All.CentralMass);
}
#endif

#ifdef SPIRAL

/*
  The function below adds the spiral contribution. Called several times
  to get gradients.
*/
double spiral_potential_calc(double x, double y, double z, double ti)
{
  double Szq    = 0.7;
  double Salpha = 0.261799;
  double Src2   = 2.019e21 / All.UnitLength_in_cm;
  double Sp0    = 2.12889e-24 * pow(All.UnitLength_in_cm, 3) / All.UnitMass_in_g;
  double Sp1    = 6.19e-31 * pow(All.UnitLength_in_cm, 3) / All.UnitMass_in_g;
  double SCo    = 2.31e14 * pow(All.UnitTime_in_s, 2) / pow(All.UnitLength_in_cm, 2);
  double SRc    = 3.09398e20 / All.UnitLength_in_cm;
  double Srs    = 2.16578e22 / All.UnitLength_in_cm;
  double Sphir  = 6.3e-16 * All.UnitTime_in_s;
  double SHz    = 5.56916e20 / All.UnitLength_in_cm;
  double SNarms = 4.0;

  double Cz[3];
  double phi, r, gamma, sum, Sr0, St0;
  double Kn, Bn, Dn;
  double np1float;
  double spiral;
  int n;

  Cz[0] = 8. / (3. * M_PI);
  Cz[1] = 0.5;
  Cz[2] = 8. / (15. * M_PI);
  Sr0   = 2.47518e22 / All.UnitLength_in_cm;
  St0   = 3.153e+15 / All.UnitTime_in_s;

  /* Convert to polar co-ords */

  if(x != 0 && y != 0)
    {
      phi = atan(y / x);
      if(x < 0)
        phi += M_PI;
    }
  else
    {
      if(x == 0 && y > 0)
        phi = M_PI / 2.;
      if(x == 0 && y < 0)
        phi = -M_PI / 2.;
      if(y == 0 && x > 0)
        phi = 0;
      if(y == 0 && x < 0)
        phi = M_PI;
    }

  r = sqrt(x * x + y * y);

  /* Calculate the potential */

  gamma = SNarms * (phi + Sphir * (St0 + ti) - log(r / Sr0) / tan(Salpha));

  sum = 0;
  for(n = 0; n < 3; n++)
    {
      np1float = 1 + (float)n;
      Kn       = np1float * SNarms / (r * sin(Salpha));
      Bn       = Kn * SHz * (1. + 0.4 * Kn * SHz);
      Dn       = (1. + Kn * SHz + 0.3 * pow(Kn * SHz, 2)) / (1. + 0.3 * Kn * SHz);

      sum += (Cz[n] / (Dn * Kn)) * cos(np1float * gamma) * pow(1. / cosh((Kn * z) / Bn), Bn);
    }

  spiral = -4. * M_PI * All.G * SHz * Sp0 * exp(-(r - Sr0) / Srs) * sum;

  return spiral;
}

#endif
