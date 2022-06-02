/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/coord_util.h
 * \date        02/2021
 * \author      Simon May
 * \brief       Inline definitions routines for coordinate handling, e.g.
 *              periodic wrapping
 */

#ifndef COORD_UTIL_H
#define COORD_UTIL_H

#include "allvars.h"
#include "dtypes.h"
#include "forcetree.h"

/* functions related to periodic wrapping */
static inline double coord_wrap(const double box_size, double x)
{
  while(x < 0)
    x += box_size;
  while(x >= box_size)
    x -= box_size;
  return x;
}

static inline double coord_nearest(const double box_size, const double box_half, const double x)
{
  return x > box_half ? x - box_size : (x < -box_half ? x + box_size : x);
}

static inline double coord_nearest_abs(const double box_size, const double box_half, const double x_abs)
{
  return x_abs > box_half ? box_size - x_abs : x_abs;
}

static inline double coord_nearest_gravity(const double box_size, const double box_half, const double x)
{
#if !defined(GRAVITY_NOT_PERIODIC) || defined(GRAVITY_TALLBOX)
  return coord_nearest(box_size, box_half, x);
#else
  return x;
#endif
}

static inline double coord_ngb_periodic_long(const double box_size, const double box_half, const double x)
{
  return coord_nearest_abs(box_size, box_half, fabs(x));
}

static inline double coord_fof_nearest_long(const double box_size, const double box_half, const double x)
{
  const double xtmp = fabs(x);
#ifdef GRAVITY_NOT_PERIODIC
  return xtmp;
#else
  return coord_nearest_abs(box_size, box_half, xtmp);
#endif
}

/*! \brief Nearest distance in x direction, accounting for periodicity.
 *
 *  \param[in] d Distance to be checked.
 *
 *  \return Nearest distance.
 */
static inline double nearest_x(const double d)
{
#ifdef REFLECTIVE_X
  return d;
#else
  return coord_nearest(boxSize_X, boxHalf_X, d);
#endif
}

/*! \brief Nearest distance in y direction, accounting for periodicity.
 *
 *  \param[in] d Distance to be checked.
 *
 *  \return Nearest distance.
 */
static inline double nearest_y(const double d)
{
#ifdef REFLECTIVE_Y
  return d;
#else
  return coord_nearest(boxSize_Y, boxHalf_Y, d);
#endif
}

/* \brief Nearest distance in z direction, accounting for periodicity.
 *
 * \param[in] d Distance to be checked.
 *
 * \return Nearest distance.
 */
static inline double nearest_z(const double d)
{
#ifdef REFLECTIVE_Z
  return d;
#else
  return coord_nearest(boxSize_Z, boxHalf_Z, d);
#endif
}

static inline double WRAP(const int dim, const double x) { return coord_wrap(All.BoxSizes[dim], x); }
static inline double WRAP_X(const double x) { return WRAP(0, x); }
static inline double WRAP_Y(const double y) { return WRAP(1, y); }
static inline double WRAP_Z(const double z) { return WRAP(2, z); }

static inline double NEAREST(const int dim, const double x) { return coord_nearest(All.BoxSizes[dim], 0.5 * All.BoxSizes[dim], x); }
static inline double NEAREST_X(const double x) { return NEAREST(0, x); }
static inline double NEAREST_Y(const double y) { return NEAREST(1, y); }
static inline double NEAREST_Z(const double z) { return NEAREST(2, z); }

static inline double GRAVITY_NEAREST(const int dim, const double x)
{
#ifdef GRAVITY_TALLBOX
  if(GRAVITY_TALLBOX == dim)
    return x;
#else
  return coord_nearest_gravity(All.BoxSizes[dim], 0.5 * All.BoxSizes[dim], x);
#endif
}
static inline double GRAVITY_NEAREST_X(const double x) { return GRAVITY_NEAREST(0, x); }
static inline double GRAVITY_NEAREST_Y(const double y) { return GRAVITY_NEAREST(1, y); }
static inline double GRAVITY_NEAREST_Z(const double z) { return GRAVITY_NEAREST(2, z); }

static inline double NGB_PERIODIC_LONG(const int dim, const double x)
{
  return coord_ngb_periodic_long(All.BoxSizes[dim], 0.5 * All.BoxSizes[dim], x);
}
static inline double NGB_PERIODIC_LONG_X(const double x) { return NGB_PERIODIC_LONG(0, x); }
static inline double NGB_PERIODIC_LONG_Y(const double y) { return NGB_PERIODIC_LONG(1, y); }
static inline double NGB_PERIODIC_LONG_Z(const double z) { return NGB_PERIODIC_LONG(2, z); }

static inline double FOF_NEAREST_LONG(const int dim, const double x)
{
  return coord_fof_nearest_long(All.BoxSizes[dim], 0.5 * All.BoxSizes[dim], x);
}
static inline double FOF_NEAREST_LONG_X(const double x) { return FOF_NEAREST_LONG(0, x); }
static inline double FOF_NEAREST_LONG_Y(const double y) { return FOF_NEAREST_LONG(1, y); }
static inline double FOF_NEAREST_LONG_Z(const double z) { return FOF_NEAREST_LONG(2, z); }

static inline double domain_shiftCoordinate(const double pos, const int dim, const enum domain_displace_mode mode)
{
  switch(mode)
    {
      case DISPLACE_POSITION_FORWARD:
        return pos + All.GlobalDisplacementVector[dim];
      case DISPLACE_POSITION_BACKWARD:
        return pos - All.GlobalDisplacementVector[dim];
      default:
        terminate("Unkown mode %d.", (int)mode);
    }
}

static inline double domain_displaceCoordinate(const double pos, const int dim, const enum domain_displace_mode mode)
{
  return WRAP(dim, domain_shiftCoordinate(pos, dim, mode));
}

static inline double wrap_position_shift(const double pos, const int dim)
{
#ifdef REFLECTIVE_X
  if(dim == 0)
    return pos;
#endif
#ifdef REFLECTIVE_Y
  if(dim == 1)
    return pos;
#endif
#ifdef REFLECTIVE_Z
  if(dim == 2)
    return pos;
#endif

  return domain_displaceCoordinate(pos, dim, DISPLACE_POSITION_BACKWARD);
}

/* \brief Periodic wrapping of point
 *
 *  If pos is more than a half-box away from ref, wrap it so that they
 *  are in the same octant.
 *
 *  \param[in] pos Position of point of interest
 *  \param[in] ref Position of reference point
 *  \param[in] dim Coordinate axis (0: x, 1: y, 2: z)
 *
 *  \return (potentially wrapped) point of interest
 */
static inline double periodic_wrap_coordinate(const double pos, const double ref, const int dim)
{
  const double box_size = All.BoxSizes[dim];
  const double box_half = 0.5 * All.BoxSizes[dim];
  const double d        = pos - ref;
  return d > box_half ? pos - box_size : (d < -box_half ? pos + box_size : pos);
}

/* \brief Periodic wrapping of point
 *
 *  If pos is more than a half-box away from ref, wrap it so that they
 *  are in the same octant.
 *
 *  \param[in, out] pos Position of point of interest
 *  \param[in] ref Position of reference point
 *
 *  \return void
 */
static inline void periodic_wrap_point(double *const pos, const double *const ref)
{
#ifndef REFLECTIVE_X
  pos[0] = periodic_wrap_coordinate(pos[0], ref[0], 0);
#endif
#ifndef REFLECTIVE_Y
  pos[1] = periodic_wrap_coordinate(pos[1], ref[1], 1);
#endif
#ifndef REFLECTIVE_Z
  pos[2] = periodic_wrap_coordinate(pos[2], ref[2], 2);
#endif
}

/* \brief Periodic wraping of point
 *
 *  If pos is more than a half-box away from ref, wrap it so that they
 *  are in the same octant.
 *
 *  \param[in, out] pos Position of point of interest
 *  \param[in] ref Position of reference point
 *
 *  \return void
 */
static inline void periodic_wrap_point_MyDouble(MyDouble *const pos, const MyDouble *const ref)
{
#ifndef REFLECTIVE_X
  pos[0] = periodic_wrap_coordinate(pos[0], ref[0], 0);
#endif
#ifndef REFLECTIVE_Y
  pos[1] = periodic_wrap_coordinate(pos[1], ref[1], 1);
#endif
#ifndef REFLECTIVE_Z
  pos[2] = periodic_wrap_coordinate(pos[2], ref[2], 2);
#endif
}

/* definitions functions related to PM computations */
#ifdef PMGRID

#define GRIDX (PMGRID * STRETCHX * DBX + DBX_EXTRA)
#define GRIDY (PMGRID * STRETCHY * DBY + DBY_EXTRA)
#define GRIDZ (PMGRID * STRETCHZ * DBZ + DBZ_EXTRA)

#define GRIDz (GRIDZ / 2 + 1)
#define GRID2 (2 * GRIDz)

#if defined(PLACEHIGHRESREGION) || defined(GRAVITY_NOT_PERIODIC)
#ifndef GRIDBOOST
#define GRIDBOOST 2
#endif

#define GRID_NP (GRIDBOOST * PMGRID)
#define GRIDz_NP (GRID_NP / 2 + 1)
#define GRID2_NP (2 * GRIDz_NP)
#endif

#if GRIDX > 1024 || GRIDY > 1024 || GRIDZ > 1024 || (defined(GRID_NP) && GRID_NP > 1024)
/* use a larger data type in this case so that we can always address all cells
 * of the 3D grid with a single index */
typedef long long large_array_offset;
#else
typedef unsigned int large_array_offset;
#endif

#ifdef NUMPART_PER_TASK_LARGE
/* if there is a risk that the local particle number times 8 overflows a 32-bit
 * integer, this data type should be used */
typedef long long large_numpart_type;
#else
typedef int large_numpart_type;
#endif

/* short-cut functions for accessing different 3D arrays */
static inline large_array_offset PM_FI(const int gridy, const int grid2, const int x, const int y, const int z)
{
  return (large_array_offset)grid2 * (gridy * x + y) + z;
}

static inline large_array_offset PM_FC(const int base_firstcol, const int grid2, const int c, const int z)
{
  return (large_array_offset)grid2 * (c - base_firstcol) + z;
}

static inline large_array_offset PM_NI(const int nslab_y, const int gridz, const int x, const int y, const int z)
{
  return (large_array_offset)gridz * (y + x * nslab_y) + z;
}

static inline large_array_offset PM_TI(const int nslab_x, const int gridz, const int x, const int y, const int z)
{
  return (large_array_offset)gridz * (x + y * nslab_x) + z;
}

#endif

#endif
