/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/domain_box.c
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
#include <strings.h>

#include "allvars.h"
#include "domain.h"
#include "proto.h"
#include "voronoi.h"

void domain_shiftPosition(double *const pos, const enum domain_displace_mode mode)
{
  for(int i = 0; i < 3; i++)
    pos[i] = domain_shiftCoordinate(pos[i], i, mode);
}

void domain_displacePosition(double *const pos, const enum domain_displace_mode mode)
{
  for(int i = 0; i < 3; i++)
    pos[i] = domain_displaceCoordinate(pos[i], i, mode);
}

void domain_displacePositions(const enum domain_displace_mode mode)
{
  for(int i = 0; i < NumPart; i++)
    {
      if(P[i].ID == 0 && P[i].Mass == 0) /* derefined */
        continue;

      domain_displacePosition(P[i].Pos, mode);

      if(i < NumGas)
        {
          domain_displacePosition(SphP[i].Center, mode);
          domain_displacePosition(SphP[i].Grad.Center, mode);
        }

#if defined(BLACK_HOLES) && defined(BH_FRICTION)
      if(P[i].Type == PTYPE_BNDRY)
        {
          domain_displacePosition(BPP(i).BH_MinPotPos_Extended, mode);
          domain_displacePosition(BPP(i).BH_MinPotPos_Previous, mode);
        }
#endif

#if defined(BLACK_HOLES) && defined(MEASURE_POTMIN_AROUND_BH)
      if(P[i].Type == PTYPE_BNDRY)
        {
          domain_displacePosition(BPP(i).BH_MinPotPos, mode);
        }
#endif
    }

#ifdef BH_BASED_CGM_ZOOM
  domain_displacePosition(All.BlackHolePosition, mode);
#endif

#ifdef PLACEHIGHRESREGION
  domain_displacePosition(All.Xmintot[1], mode);
  domain_displacePosition(All.Xmaxtot[1], mode);
  domain_displacePosition(All.Corner[1], mode);
  domain_displacePosition(All.UpperCorner[1], mode);
#endif
}

/*! \brief Finds the extent of the global domain grid.
 *
 *  The minimum extent is the box size.
 *
 *  \return void
 */
void domain_findExtent(void)
{
  /* determine local extension */
  /* preset to simulation box, taking care of stretched box */
  double xmin[] = {0, 0, 0};
  double xmax[] = {boxSize_X, boxSize_Y, boxSize_Z};

  for(int i = 0; i < NumPart; i++)
    {
#ifdef ADDBACKGROUNDGRID
      if(P[i].Type != 0)
        continue;
#endif
      for(int j = 0; j < 3; j++)
        {
          if(xmin[j] > P[i].Pos[j])
            xmin[j] = P[i].Pos[j];

          if(xmax[j] < P[i].Pos[j])
            xmax[j] = P[i].Pos[j];
        }
    }

  double xmin_glob[3], xmax_glob[3];
  MPI_Allreduce(xmin, xmin_glob, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(xmax, xmax_glob, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

#ifdef ADDBACKGROUNDGRID
  for(int j = 0; j < 3; j++)
    if(xmax_glob[j] < All.BoxSize)
      xmax_glob[j] = All.BoxSize;

  for(int j = 0; j < 3; j++)
    if(xmin_glob[j] > 0)
      xmin_glob[j] = 0;
#endif

  double len = 0;
  for(int j = 0; j < 3; j++)
    if(xmax_glob[j] - xmin_glob[j] > len)
      len = xmax_glob[j] - xmin_glob[j];

#if defined(GRAVITY_NOT_PERIODIC) && !defined(ADDBACKGROUNDGRID)
  len *= 1.2; /* enlarge box a bit to avoid triggering of an out of box recovery */
#else
  len *= 1.00001;
#endif

#ifndef AMR

#if defined(DO_NOT_RANDOMIZE_DOMAINCENTER) || !defined(GRAVITY_NOT_PERIODIC) || defined(ONEDIMS) || defined(TWODIMS)
  for(int j = 0; j < 3; j++)
    {
      DomainCenter[j] = 0.5 * (xmin_glob[j] + xmax_glob[j]);
      DomainCorner[j] = 0.5 * (xmin_glob[j] + xmax_glob[j]) - 0.5 * len;
    }
#else
  for(int j = 0; j < 3; j++)
    {
      DomainCenter[j] = 0.5 * (xmin_glob[j] + xmax_glob[j]);
      DomainCenter[j] += (2. * get_random_number() - 1.) * 0.5 * len;
    }

  MPI_Bcast(DomainCenter, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  len *= 2;

  for(int j = 0; j < 3; j++)
    DomainCorner[j] = DomainCenter[j] - 0.5 * len;
#endif

  DomainLen = len;

#else

  len             = fmax(boxSize_X, fmax(boxSize_Y, boxSize_Z));
  DomainLen       = len;
  DomainCorner[0] = 0.;
  DomainCorner[1] = 0.;
  DomainCorner[2] = 0.;

  DomainCenter[0] = DomainLen / 2.;
  DomainCenter[1] = DomainLen / 2.;
  DomainCenter[2] = DomainLen / 2.;
#endif

  DomainInverseLen = 1.0 / DomainLen;
  DomainFac        = 1.0 / len * (((peanokey)1) << (BITS_PER_DIMENSION));
  DomainBigFac     = (DomainLen / (((long long)1) << 52));
}

/*! \brief Makes sure all particles are within the box.
 *
 *  This function makes sure that all particle coordinates (Pos) are
 *  periodically mapped onto the interval [0, BoxSize). After this function
 *  has been called, a new domain decomposition should be done, which will
 *  also force a new tree construction.
 *
 *  \return void
 */
void do_box_wrapping(void)
{
#ifdef ADDBACKGROUNDGRID
  return;
#endif

#if !defined(GRAVITY_NOT_PERIODIC) && !defined(DO_NOT_RANDOMIZE_DOMAINCENTER) && defined(SELFGRAVITY) && (NUMDIMS > 2)
  domain_displacePositions(DISPLACE_POSITION_BACKWARD);

  if(ThisTask == 0)
    {
      double prefac = 1.;
#ifdef PLACEHIGHRESREGION
      prefac = 0.5;
#endif
      for(int j = 0; j < 3; j++)
        All.GlobalDisplacementVector[j] = (get_random_number() - 0.5) * All.BoxSizes[j] * prefac;
    }

  mpi_printf("DOMAIN: New global displacement vector: %g, %g, %g, box=%g, %g, %g\n", All.GlobalDisplacementVector[0],
             All.GlobalDisplacementVector[1], All.GlobalDisplacementVector[2], All.BoxSizes[0], All.BoxSizes[1], All.BoxSizes[2]);
  MPI_Bcast(All.GlobalDisplacementVector, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  domain_displacePositions(DISPLACE_POSITION_FORWARD);
#endif

  for(int i = 0; i < NumPart; i++)
    {
#if defined(VORONOI_DYNAMIC_UPDATE) || defined(AMR_CONNECTIONS)
      if(i < NumGas)
        trans_table[i].wrapped = 0;
#endif

#ifdef GRAVITY_NOT_PERIODIC
      if(P[i].Type != 0)
        continue;
#endif

#ifdef VORONOI_DYNAMIC_UPDATE
      if(i < NumGas)
        {
          int flag = 1;
          for(int j = 0; j < 3; j++)
            {
              if(P[i].Pos[j] < 0)
                trans_table[i].wrapped |= flag;
              flag <<= 1;
              if(P[i].Pos[j] >= All.BoxSizes[j])
                trans_table[i].wrapped |= flag;
              flag <<= 1;
            }
        }
#endif

#ifdef REFLECTIVE_X
      if(P[i].Pos[0] < 0 || P[i].Pos[0] >= All.BoxSizes[0])
        terminate("i=%d ID=%" MYIDTYPE_PRI " type=%d moved out of box. x=%g", i, P[i].ID, P[i].Type, P[i].Pos[0]);
#else
      P[i].Pos[0] = WRAP_X(P[i].Pos[0]);
#endif
#ifdef REFLECTIVE_Y
      if(P[i].Pos[1] < 0 || P[i].Pos[1] >= All.BoxSizes[1])
        terminate("i=%d ID=%" MYIDTYPE_PRI " type=%d moved out of box. y=%g", i, P[i].ID, P[i].Type, P[i].Pos[1]);
#else
      P[i].Pos[1] = WRAP_Y(P[i].Pos[1]);
#endif
#ifdef REFLECTIVE_Z
      if(P[i].Pos[2] < 0 || P[i].Pos[2] >= All.BoxSizes[2])
        terminate("i=%d ID=%" MYIDTYPE_PRI " type=%d moved out of box. z=%g", i, P[i].ID, P[i].Type, P[i].Pos[2]);
#else
      P[i].Pos[2] = WRAP_Z(P[i].Pos[2]);
#endif
    }
}
