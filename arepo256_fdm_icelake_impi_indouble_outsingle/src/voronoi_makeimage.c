/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/voronoi_makeimage.c
 * \date        MM/YYYY
 * \author
 * \brief       Routines to create images such as projection through a Voronoi-mesh
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
#include <sys/stat.h>
#include <sys/types.h>

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"

#ifndef GFM_N_CHEM_ELEMENTS
#define GFM_N_CHEM_ELEMENTS 0
#endif

/* \brief Limits of image
 *
 *  Extracts the limits and resolutions for the image slice and
 *  projection routines. These are largely the same for both
 *  projection and slicing, so this means less duplication, and that
 *  the slicing routine gets the same centering options as the
 *  projection ones.
 *
 *  \param argc[in] Number of parameters handed over to Arepo
 *  \param argv[in] Pointer to parameters handed over to Arepo
 *  \param restart_flag[in] Mode in which Arepo is called
 *  \param pixels_x[out] Number of pixels in x direction
 *  \param pixels_y[out] Number of pixels in y direction
 *  \param xaxis[out] Projection x axis in simulation coordinates
 *  \param yaxis[out] Projeciton y axis in simulation coordinates
 *  \param zaxis[out] Projeciton z axis in simulation coordinates
 *  \param xmin[out] Minimum x coordinate in projection system
 *  \param xmax[out] Maximum x coordinate in projection system
 *  \param ymin[out] Minimum y coordinate in projection system
 *  \param ymax[out] Maximum y coordinate in projection system
 *  \param zmin[out] Minimum z coordinate in projection system
 *  \param zmax[out] Maximum z coordinate in projection system
 *  \param weight_flag[out] Weight flag: (1: dens^2, else: volume)
 *
 *  \return 0: wrong number of parameters; -1 normal exit
 */
int get_image_limits(int argc, char **argv, int restart_flag, int *pixels_x, int *pixels_y, int *xaxis, int *yaxis, int *zaxis,
                     double *xmin, double *xmax, double *ymin, double *ymax, double *zmin, double *zmax, int *weight_flag)
{
#ifdef TGSET
  if(Argc != 13 && Argc != 16)
    {
      mpi_printf(
          "\nwrong parameters found: Expecting:\n Call with: <ParameterFile> %d <SnapNum> <pixelsX> <pixelsY> <axisX> <axisY> <axisZ> "
          "<center_flag> <weight_flag>\n <width_flag> <width> <center_x> <center_y> <center_z> (The last three only needed for "
          "center_flag = 1)\n\n",
          restart_flag);
      return 0;
    }
#else
  if(!((restart_flag == RESTART_SLICE && argc == 14) || (restart_flag == RESTART_PROJECTION && (argc == 15 || argc == 16))
#ifdef SUBBOX_SNAPSHOTS
       || (restart_flag == RESTART_PROJECTION && (argc == 16 || argc == 17))
#endif
       || (restart_flag == RESTART_PROJECTION_GRID_RAYTRACING && argc == 15)))
    {
      mpi_printf(
          "\nwrong parameters found: Expecting:\n Call with: <ParameterFile> %d <SnapNum> <pixelsX> <pixelsY> <axisX> <axisY> <axisZ> "
          " <xmin> <xmax> <ymin> <ymax> %s\n",
          restart_flag, restart_flag == RESTART_SLICE ? "<zval>" : "<zmin> <zmax>");
      if(restart_flag == RESTART_PROJECTION)
        mpi_printf("Optionally also with <weight_flag> as last argument (default is 0)\n\n");
      return 0;
    }
#endif

  *pixels_x = atoi(argv[4]);
  *pixels_y = atoi(argv[5]);

  *xaxis = atoi(argv[6]);
  *yaxis = atoi(argv[7]);
  *zaxis = atoi(argv[8]);

#ifdef TGSET
  tgset_get_image_limits(argv, xaxis, yaxis, zaxis, xmin, xmax, ymin, ymax, zmin, zmax, weight_flag);
#else
  *xmin = atof(argv[9]);
  *xmax = atof(argv[10]);
  *ymin = atof(argv[11]);
  *ymax = atof(argv[12]);
  *zmin = atof(argv[13]);

  if(restart_flag == RESTART_PROJECTION || restart_flag == RESTART_PROJECTION_GRID_RAYTRACING)
    {
#ifndef SUBBOX_SNAPSHOTS
      *zmax = atof(argv[14]);
      if(Argc == 16)
        {
          *weight_flag = atoi(Argv[15]);
          mpi_printf("VORONOI_MAKEIMAGE: weight_flag = %d\n", *weight_flag);
        }
#else
      *zmax = atof(argv[14]);
      if(Argc == 17)
        {
          *weight_flag = atoi(Argv[16]);
          mpi_printf("VORONOI_MAKEIMAGE: weight_flag = %d\n", *weight_flag);
        }
#endif
    }
#endif

  return -1;
}

#ifdef VORONOI

#if defined(VORONOI_IMAGES_FOREACHSNAPSHOT) || defined(VORONOI_FREQUENT_IMAGES)
/*! \brief Creates a projected image.
 *
 *  \param[in] num Index in filename.
 *
 *  \return void
 */
void make_voronoi_image(int num)
{
#if !defined(TWODIMS) && !defined(ONEDIMS)

  make_3d_voronoi_projected_image(num, 1, All.PicXpixels, All.PicYpixels, All.PicXaxis, All.PicYaxis, All.PicZaxis, All.PicXmin,
                                  All.PicXmax, All.PicYmin, All.PicYmax, All.PicZmin, All.PicZmax, 0);
#ifdef VORONOI_MULTIPLE_PROJECTIONS
  make_3d_voronoi_projected_image(num, 1, All.PicXpixels, All.PicYpixels, All.PicXaxis, All.PicZaxis, All.PicYaxis, All.PicXmin,
                                  All.PicXmax, All.PicZmin, All.PicZmax, All.PicYmin, All.PicYmax, 0);
#endif

#ifdef VORONOI_NOGRADS
  make_3d_voronoi_projected_image(num, 0, All.PicXpixels, All.PicYpixels, All.PicXaxis, All.PicYaxis, All.PicZaxis, All.PicXmin,
                                  All.PicXmax, All.PicYmin, All.PicYmax, All.PicZmin, All.PicZmax, 0);
#ifdef VORONOI_MULTIPLE_PROJECTIONS
  make_3d_voronoi_projected_image(num, 0, All.PicXpixels, All.PicYpixels, All.PicXaxis, All.PicXaxis, All.PicYaxis, All.PicXmin,
                                  All.PicXmax, All.PicZmin, All.PicZmax, All.PicYmin, All.PicYmax, 0);
#endif
#endif /* #ifdef VORONOI_NOGRADS */
#endif /* #if !defined(TWODIMS) && !defined (ONEDIMS) */

#ifdef TWODIMS
  make_2d_voronoi_image(num, All.PicXpixels, All.PicYpixels);
#endif /* #ifdef TWODIMS */
}

/*! \brief Creates a slice image.
 *
 *  \param[in] num Index in filename.
 *
 *  \return void
 */
void make_voronoi_image_slice(int num)
{
#if !defined(TWODIMS) && !defined(ONEDIMS)
  make_3d_voronoi_slice_image(num, 1, All.PicXpixels, All.PicYpixels, All.PicXaxis, All.PicYaxis, All.PicZaxis, All.PicXmin,
                              All.PicXmax, All.PicYmin, All.PicYmax, All.PicZmin);
#ifdef VORONOI_NOGRADS
  make_3d_voronoi_slice_image(num, 0, All.PicXpixels, All.PicYpixels, All.PicXaxis, All.PicYaxis, All.PicZaxis, All.PicXmin,
                              All.PicXmax, All.PicYmin, All.PicYmax, All.PicZmin);
#endif /* #ifdef VORONOI_NOGRADS */
#endif /* #if !defined(TWODIMS) && !defined (ONEDIMS) */
}
#endif /* #if defined(VORONOI_IMAGES_FOREACHSNAPSHOT) || defined(VORONOI_FREQUENT_IMAGES) */

#if !defined(TWODIMS) && !defined(ONEDIMS) /* will only be compiled in 3D case */

#define INSIDE_EPS 1.0e-6
#define MAX_COUNT_MOVES 500000

/*! \brief Probes if there is an intersection and returns intersecting elements
 *
 *  Tests the line segment between ppstart and ppend for intersections
 *  with the tetrahedron tt. Count is set to the number of
 *  intersections found. face, edge, and corner is set to the
 *  intersecting elements. orientations is a int[6] array with the
 *  orientations of the points. (If the line segment intersects more
 *  than one face, only one of them will be returned. What do we do in
 *  that case?)
 *
 *  \param[in] T Pointer to tessellation
 *  \param[in] tt index in DT array
 *  \param[in] pstart Pointer to starting point
 *  \param[in] pend Pointer to end point
 *  \param[out] count Number of intersections
 *  \param[out] f Intersecting face
 *  \param[out] edge Intersecting edge
 *  \param[out] corner Intersecting corner
 *  \param[out] orientation Int[6] array with the orientations of the points
 *
 *  \return void
 */
void voronoi_probe_intersections(tessellation *T, int tt, point *pstart, point *pend, int *count, int *f, int *edge, int *corner,
                                 int *orientations)
{
  tetra *DT = T->DT;
  point *DP = T->DP;

  int o1, o2, o3, o4, o5, o6;

  tetra *t = &DT[tt];

  point *p0 = &DP[t->p[0]];
  point *p1 = &DP[t->p[1]];
  point *p2 = &DP[t->p[2]];
  point *p3 = &DP[t->p[3]];

  o1 = Orient3d(p0, p1, pend, pstart);
  o2 = Orient3d(p1, p2, pend, pstart);
  o3 = Orient3d(p2, p0, pend, pstart);
  o4 = Orient3d(p1, p3, pend, pstart);
  o5 = Orient3d(p3, p2, pend, pstart);
  o6 = Orient3d(p3, p0, pend, pstart);

  /*
     o1 = Orient3d_Exact(p0, p1, pend, pstart);
     o2 = Orient3d_Exact(p1, p2, pend, pstart);
     o3 = Orient3d_Exact(p2, p0, pend, pstart);
     o4 = Orient3d_Exact(p1, p3, pend, pstart);
     o5 = Orient3d_Exact(p3, p2, pend, pstart);
     o6 = Orient3d_Exact(p3, p0, pend, pstart);
   */

  *f      = -1;
  *edge   = -1;
  *corner = -1;
  *count  = 0;

  if(o4 == 1 && o5 == 1 && o2 == -1)
    {
      *f     = 0;
      *count = *count + 1;
    }
  if(o3 == -1 && o5 == -1 && o6 == 1)
    {
      *f     = 1;
      *count = *count + 1;
    }
  if(o4 == -1 && o1 == -1 && o6 == -1)
    {
      *f     = 2;
      *count = *count + 1;
    }
  if(o1 == 1 && o2 == 1 && o3 == 1)
    {
      *f     = 3;
      *count = *count + 1;
    }

  /* edge between face 0 and face 1 (is edge between corners 2 & 3) */
  if((o4 == 1 && o5 == 0 && o2 == -1 && o3 == -1 && o5 == 0 && o6 == 1) ||
     (o4 == 1 && o5 == 0 && o2 == -1 && o3 == 0 && o5 == 0 && o6 == 0) ||
     (o4 == 0 && o5 == 0 && o2 == 0 && o3 == -1 && o5 == 0 && o6 == 1))
    {
      *edge  = 5;
      *count = *count + 1;
    }

  /* edge between face 0 and face 2 (is edge between corners 1 & 3) */
  if((o4 == 0 && o5 == 1 && o2 == -1 && o4 == 0 && o1 == -1 && o6 == -1) ||
     (o4 == 0 && o5 == 1 && o2 == -1 && o4 == 0 && o1 == 0 && o6 == 0) ||
     (o4 == 0 && o5 == 0 && o2 == 0 && o4 == 0 && o1 == -1 && o6 == -1))
    {
      *edge  = 4;
      *count = *count + 1;
    }

  /* edge between face 0 and face 3 (is edge between corners 1 & 2) */
  if((o4 == 1 && o5 == 1 && o2 == 0 && o1 == 1 && o2 == 0 && o3 == 1) ||
     (o4 == 1 && o5 == 1 && o2 == 0 && o1 == 0 && o2 == 0 && o3 == 0) ||
     (o4 == 0 && o5 == 0 && o2 == 0 && o1 == 1 && o2 == 0 && o3 == 1))
    {
      *edge  = 3;
      *count = *count + 1;
    }

  /* edge between face 1 and face 2 (is edge between corners 0 & 3) */
  if((o3 == -1 && o5 == -1 && o6 == 0 && o4 == -1 && o1 == -1 && o6 == 0) ||
     (o3 == -1 && o5 == -1 && o6 == 0 && o4 == 0 && o1 == 0 && o6 == 0) ||
     (o3 == 0 && o5 == 0 && o6 == 0 && o4 == -1 && o1 == -1 && o6 == 0))
    {
      *edge  = 2;
      *count = *count + 1;
    }

  /* edge between face 1 and face 3 (is edge between corners 0 & 2) */
  if((o3 == 0 && o5 == -1 && o6 == 1 && o1 == 1 && o2 == 1 && o3 == 0) ||
     (o3 == 0 && o5 == -1 && o6 == 1 && o1 == 0 && o2 == 0 && o3 == 0) ||
     (o3 == 0 && o5 == 0 && o6 == 0 && o1 == 1 && o2 == 1 && o3 == 0))
    {
      *edge  = 1;
      *count = *count + 1;
    }

  /* edge between face 2 and face 3 (is edge between corners 0 & 1) */
  if((o4 == -1 && o1 == 0 && o6 == -1 && o1 == 0 && o2 == 1 && o3 == 1) ||
     (o4 == -1 && o1 == 0 && o6 == -1 && o1 == 0 && o2 == 0 && o3 == 0) ||
     (o4 == 0 && o1 == 0 && o6 == 0 && o1 == 0 && o2 == 1 && o3 == 1))
    {
      *edge  = 0;
      *count = *count + 1;
    }

  /* corner 0 (true if all adjacent edges are cut */

  if(o3 == 0 && o5 == -1 && o6 == 0 && o4 == -1 && o1 == 0 && o2 == 1)
    {
      *corner = 0;
      *count  = *count + 1;
    }

  /* corner 1 (true if all adjacent edges are cut */
  if(o4 == 0 && o1 == 0 && o6 == -1 && o2 == 0 && o3 == 1 && o5 == 1)
    {
      *corner = 1;
      *count  = *count + 1;
    }

  /* corner 2 (true if all adjacent edges are cut */
  if(o3 == 0 && o5 == 0 && o6 == 1 && o1 == 1 && o2 == 0 && o4 == 1)
    {
      *corner = 2;
      *count  = *count + 1;
    }

  /* corner 3 (true if all adjacent edges are cut */
  if(o5 == 0 && o6 == 0 && o4 == 0 && o1 == -1 && o2 == -1 && o3 == -1)
    {
      *corner = 3;
      *count  = *count + 1;
    }

  orientations[0] = o1;
  orientations[1] = o2;
  orientations[2] = o3;
  orientations[3] = o4;
  orientations[4] = o5;
  orientations[5] = o6;
}

/*! \brief Gets next tetrahedron
 *
 *  Tests whether the line defined by [pstart, pend] has an
 *  intersection with the tetrahedron tt. IF yes, 1 is returned and the
 *  coordinates of the point ppexit are updated with the coordinates of
 *  the exit point.  If no (or if it intersects a corner?) 0 is returned
 *  and ppexit is updated to the index of the corner point. If the tetra
 *  contains infinity points and the segment did not cross a valid
 *  (finity) face, 0 is returned and ppexit is unset. previous_tetra
 *  indicates a previously visited tetrahedron in case a ray is followed
 *  from pstart to pend.
 *
 *  \param[in] T Pointer to tessellation
 *  \param[in] tt The DT index of the tetrahedron to test
 *  \param[in] pstart The DP index of the starting point
 *  \param[in] pend The DP index of the end point
 *  \param[out] nexttetra The DT index of the next tetrahedron
 *  \param[out] pexit DP index of the point in which the exit coordinates will
 *  be placed
 *  \param[in] previous_tetra previous tetrahedron
 *
 *  \return (0,1)
 */
int image_get_next_tetra(tessellation *T,
                         /// The DT index of the tetrahedron to test
                         int tt,
                         /// The DP index of the starting point
                         point *pstart,
                         /// The DP index of the end point
                         point *pend,
                         /// [out] The DT index of the next tetrahedron
                         int *nexttetra,
                         /// [out] DP index of the point in which the exit coordinates will be placed
                         point *pexit, int *previous_tetra)
{
  tetra *DT = T->DT;
  point *DP = T->DP;
  tetra *t  = &DT[tt];

  double w, u, v, s = 0;
  double denom;
  int status;

  // int pp0, pp1, pp2, pp3;
  point *p0, *p1, *p2, *p3;

  /*
     pp0 = t->p[0];
     pp1 = t->p[1];
     pp2 = t->p[2];
     pp3 = t->p[3];
   */

  p0 = &DP[t->p[0]];
  p1 = &DP[t->p[1]];
  p2 = &DP[t->p[2]];
  p3 = &DP[t->p[3]];

  // This is not a valid end condition! We could go through the valid face!
  // XXX not true anymore, is it?
#ifndef SUNRISE
  if(isInfinity(p0) || isInfinity(p1) || isInfinity(p2) || isInfinity(p3))
    {
      printf("we are in a tetraeder with in infinity point. tetra=%d\n", tt);
      terminate("tetraeder with in infinity point");
    }
#endif

  int face, edge, corner, count, orientations[6];

  voronoi_probe_intersections(T, tt, pstart, pend, &count, &face, &edge, &corner, orientations);

  if(count != 1)
    {
      if(count == 0)
        // no intersection found, that means (I think) we're in a
        // tetra with infinity points and did not cross a valid face
        return 0;

      printf("count=%d\n", count);

      printf("orientations: %d %d %d %d %d %d\n", orientations[0], orientations[1], orientations[2], orientations[3], orientations[4],
             orientations[5]);
      printf("tetra: %d   adjacent ones %d %d %d %d\n", tt, t->t[0], t->t[1], t->t[2], t->t[3]);
      printf("points: %ld %ld %ld %ld\n", p0 - DP, p1 - DP, p2 - DP, p3 - DP);
      printf("face-cuts:  %d\n", face);
      printf("edge-cuts:  %d\n", edge);
      printf("corner-cuts: %d\n", corner);

      printf("\n");

      printf("p0:     %g %g %g\n", p0->x - pstart->x, p0->y - pstart->y, p0->z - pstart->z);
      printf("p1:     %g %g %g\n", p1->x - pstart->x, p1->y - pstart->y, p1->z - pstart->z);
      printf("p2:     %g %g %g\n", p2->x - pstart->x, p2->y - pstart->y, p2->z - pstart->z);
      printf("p3:     %g %g %g\n", p3->x - pstart->x, p3->y - pstart->y, p3->z - pstart->z);
      printf("\n");
      printf("pstart: %g %g %g\n", pstart->x, pstart->y, pstart->z);
      printf("pend:   %g %g %g\n", pend->x, pend->y, pend->z);

      terminate("we have not found one exit intersection with the tetrahedron");
    }

#ifndef OPTIMIZE_MEMORY_USAGE
  double ax = p1->xx - p0->xx;
  double ay = p1->yy - p0->yy;
  double az = p1->zz - p0->zz;

  double bx = p2->xx - p0->xx;
  double by = p2->yy - p0->yy;
  double bz = p2->zz - p0->zz;

  double cx = p3->xx - p0->xx;
  double cy = p3->yy - p0->yy;
  double cz = p3->zz - p0->zz;

  double qx = pend->xx - p0->xx;
  double qy = pend->yy - p0->yy;
  double qz = pend->zz - p0->zz;
#else
  double ax, ay, az, bx, by, bz, cx, cy, cz, qx, qy, qz;
  IntegerMapType pA_ixyz[3], pB_ixyz[3];
  double pA_xyz[3], pB_xyz[3];

  get_integers_for_point(p0, pA_ixyz, pA_xyz);

  get_integers_for_point(p1, pB_ixyz, pB_xyz);
  ax = pB_xyz[0] - pA_xyz[0];
  ay = pB_xyz[1] - pA_xyz[1];
  az = pB_xyz[2] - pA_xyz[2];

  get_integers_for_point(p2, pB_ixyz, pB_xyz);
  bx = pB_xyz[0] - pA_xyz[0];
  by = pB_xyz[1] - pA_xyz[1];
  bz = pB_xyz[2] - pA_xyz[2];

  get_integers_for_point(p3, pB_ixyz, pB_xyz);
  cx = pB_xyz[0] - pA_xyz[0];
  cy = pB_xyz[1] - pA_xyz[1];
  cz = pB_xyz[2] - pA_xyz[2];

  get_integers_for_point(pend, pB_ixyz, pB_xyz);
  qx = pB_xyz[0] - pA_xyz[0];
  qy = pB_xyz[1] - pA_xyz[1];
  qz = pB_xyz[2] - pA_xyz[2];
#endif

  double mv_data[] = {ax, bx, cx, qx, ay, by, cy, qy, az, bz, cz, qz};
  double xend[3];

  status = solve_linear_equations(mv_data, xend);

  if(status < 0)
    {
      /*
         printf("warning. Potentially inaccurate solution of linear system\n");
         printf("a:     %g %g %g\n", ax, ay, az);
         printf("b:     %g %g %g\n", bx, by, bz);
         printf("b:     %g %g %g\n", cx, cy, cz);
         printf("q:     %g %g %g\n", qx, qy, qz);
         printf("vol=%g test_orient=%d exact_orient=%d\n",
         calculate_tetra_volume(p0, p1, p2, p3),
         test_tetra_orientation(p0, p1, p2, p3),
         Orient3d_Exact(p0, p1, p2, p3));
       */
      //      terminate("status < 0");
    }

#ifndef OPTIMIZE_MEMORY_USAGE
  double qstartx = pstart->xx - p0->xx;
  double qstarty = pstart->yy - p0->yy;
  double qstartz = pstart->zz - p0->zz;
#else
  get_integers_for_point(pstart, pB_ixyz, pB_xyz);

  double qstartx = pB_xyz[0] - pA_xyz[0];
  double qstarty = pB_xyz[1] - pA_xyz[1];
  double qstartz = pB_xyz[2] - pA_xyz[2];
#endif

  double mv_start[] = {ax, bx, cx, qstartx, ay, by, cy, qstarty, az, bz, cz, qstartz};
  double xstart[3];

  status = solve_linear_equations(mv_start, xstart);

  if(status < 0)
    {
      /*
         printf("warning. Potentially inaccurate solution of linear system\n");
         printf("a:     %g %g %g\n", ax, ay, az);
         printf("b:     %g %g %g\n", bx, by, bz);
         printf("b:     %g %g %g\n", cx, cy, cz);
         printf("q:     %g %g %g\n", qx, qy, qz);
         printf("vol=%g test_orient=%d exact_orient=%d\n",
         calculate_tetra_volume(p0, p1, p2, p3),
         test_tetra_orientation(p0, p1, p2, p3),
         Orient3d_Exact(p0, p1, p2, p3));
       */
      //      terminate("status < 0");
    }

  *nexttetra = *previous_tetra;

  if(face >= 0 && face <= 3)
    {
      if(face == 0) /* cut with face 0 */
        {
          denom = (xend[0] + xend[1] + xend[2] - (xstart[0] + xstart[1] + xstart[2]));
          if(denom != 0)
            s = (1 - (xstart[0] + xstart[1] + xstart[2])) / denom;
          else
            terminate("denom==0");
        }

      if(face == 1) /* cut with face 1 */
        {
          denom = (xstart[0] - xend[0]);

          if(denom != 0)
            s = xstart[0] / denom;
          else
            terminate("denom==0");
        }

      if(face == 2) /* cut with face 2 */
        {
          denom = (xstart[1] - xend[1]);

          if(denom != 0)
            s = xstart[1] / denom;
          else
            terminate("denom == 0");
        }

      if(face == 3) /* cut with face 3 */
        {
          denom = (xstart[2] - xend[2]);

          if(denom != 0)
            s = xstart[2] / denom;
          else
            terminate("denom == 0");
        }
      *nexttetra = t->t[face];
    }

  if(edge >= 0 && face <= 5)
    {
      int o1, o2, o3, o4, o5, o6;

      o1 = orientations[0];
      o2 = orientations[1];
      o3 = orientations[2];
      o4 = orientations[3];
      o5 = orientations[4];
      o6 = orientations[5];

      if(edge == 0)
        {
          if(o4 != 0 || o1 != 0 || o6 != 0)
            {
              denom = (xstart[1] - xend[1]);
              if(denom != 0)
                s = xstart[1] / denom;
              else
                terminate("denom == 0");
            }
          else
            {
              denom = (xstart[2] - xend[2]);
              if(denom != 0)
                s = xstart[2] / denom;
              else
                terminate("denom == 0");
            }
        }
      else if(edge == 1)
        {
          if(o3 != 0 || o5 != 0 || o6 != 0)
            {
              denom = (xstart[0] - xend[0]);
              if(denom != 0)
                s = xstart[0] / denom;
              else
                terminate("denom == 0");
            }
          else
            {
              denom = (xstart[2] - xend[2]);
              if(denom != 0)
                s = xstart[2] / denom;
              else
                {
                  printf("xstart[0]-xend[0]=%g\n", xstart[0] - xend[0]);
                  printf("xstart[1]-xend[1]=%g\n", xstart[1] - xend[1]);
                  printf("xstart[2]-xend[2]=%g\n", xstart[2] - xend[2]);
                  terminate("denom == 0");
                }
            }
        }
      else if(edge == 2)
        {
          if(o3 != 0 || o5 != 0 || o6 != 0)
            {
              denom = (xstart[0] - xend[0]);
              if(denom != 0)
                s = xstart[0] / denom;
              else
                terminate("denom == 0");
            }
          else
            {
              denom = (xstart[1] - xend[1]);
              if(denom != 0)
                s = xstart[1] / denom;
              else
                terminate("denom == 0");
            }
        }
      else if(edge == 3)
        {
          if(o4 != 0 || o5 != 0 || o2 != 0)
            {
              denom = (xend[0] + xend[1] + xend[2] - (xstart[0] + xstart[1] + xstart[2]));
              if(denom != 0)
                s = (1 - (xstart[0] + xstart[1] + xstart[2])) / denom;
              else
                terminate("denom == 0");
            }
          else
            {
              denom = (xstart[2] - xend[2]);
              if(denom != 0)
                s = xstart[2] / denom;
              else
                terminate("denom == 0");
            }
        }
      else if(edge == 4)
        {
          if(o4 != 0 || o5 != 0 || o2 != 0)
            {
              denom = (xend[0] + xend[1] + xend[2] - (xstart[0] + xstart[1] + xstart[2]));
              if(denom != 0)
                s = (1 - (xstart[0] + xstart[1] + xstart[2])) / denom;
              else
                terminate("denom == 0");
            }
          else
            {
              denom = (xstart[1] - xend[1]);
              if(denom != 0)
                s = xstart[1] / denom;
              else
                terminate("denom == 0");
            }
        }
      else if(edge == 5)
        {
          if(o4 != 0 || o5 != 0 || o2 != 0)
            {
              denom = (xend[0] + xend[1] + xend[2] - (xstart[0] + xstart[1] + xstart[2]));
              if(denom != 0)
                s = (1 - (xstart[0] + xstart[1] + xstart[2])) / denom;
              else
                terminate("denom == 0");
            }
          else
            {
              denom = (xstart[0] - xend[0]);
              if(denom != 0)
                s = xstart[0] / denom;
              else
                terminate("denom==0");
            }
        }

      /* now we circle around the edge to find the correct next tetrahedron */

      int i, j, k, l, m, ii, jj, kk, ll, pp, nn, iter;
      tetra *prev, *next;

      i = edge_start[edge];
      j = edge_end[edge];
      k = edge_opposite[edge];
      l = edge_nexttetra[edge];

      iter = 0;

      pp   = tt;
      prev = &DT[pp];
      do
        {
          nn   = prev->t[l];
          next = &DT[nn];

          //      printf("nn=%d tt=%d *previous_tetra=%d\n", nn, tt, *previous_tetra);

          if(nn != tt && nn != (*previous_tetra))
            {
              int face2, edge2, corner2, count2, orientations2[6];

              voronoi_probe_intersections(T, nn, pstart, pend, &count2, &face2, &edge2, &corner2, orientations2);
              /*
                 printf("count2=%d\n", count2);
                 printf("orientations2: %d %d %d %d %d %d\n", orientations2[0], orientations2[1], orientations2[2], orientations2[3],
                 orientations2[4], orientations2[5]);
               */
              if(count2 > 1)
                terminate("count2 > 1");

              if(count2 == 1)
                {
                  *nexttetra = nn;
                  break;
                }
            }

          for(m = 0, ll = ii = jj = -1; m < 4; m++)
            {
              if(next->p[m] == prev->p[k])
                ll = m;
              if(next->p[m] == prev->p[i])
                ii = m;
              if(next->p[m] == prev->p[j])
                jj = m;
            }

          if(ll < 0 || ii < 0 || jj < 0)
            terminate("inconsistency");

          kk = 6 - (ll + ii + jj);

          prev = next;
          pp   = nn;
          i    = ii;
          l    = ll;
          j    = jj;
          k    = kk;

          iter++;

          if(iter > 1000)
            terminate("iter is too large");
        }
      while(nn != tt);
    }

  if(corner >= 0 && corner <= 3)
    {
      *pexit = DP[t->p[corner]];

      if(pow(pexit->x - pstart->x, 2) + pow(pexit->y - pstart->y, 2) + pow(pexit->z - pstart->z, 2) >
             pow(pend->x - pstart->x, 2) + pow(pend->y - pstart->y, 2) + pow(pend->z - pstart->z, 2)
#ifdef SUNRISE
         || t->p[corner] == DPinfinity
#endif
      )
        {
          *pexit = *pend;
          return 0;
        }

      /* now we look among all tetras that share the point for a suitable next tetrahedron */

      int i;
      int pp = t->p[corner];

      for(i = 0; i < T->Ndt; i++)
        {
          if(DT[i].p[0] == pp || DT[i].p[1] == pp || DT[i].p[2] == pp || DT[i].p[3] == pp)
            if(i != tt && i != (*previous_tetra))
              {
                if(DT[i].t[0] < 0) /* don't consider deleted ones */
                  continue;

                int face2, edge2, corner2, count2, orientations2[6];

                voronoi_probe_intersections(T, i, pstart, pend, &count2, &face2, &edge2, &corner2, orientations2);
                /*
                   printf("count2=%d\n", count2);
                   printf("orientations2: %d %d %d %d %d %d\n", orientations2[0], orientations2[1], orientations2[2], orientations2[3],
                   orientations2[4], orientations2[5]);
                 */
                if(count2 > 1)
                  terminate("count2 > 1");

                if(count2 == 1)
                  {
                    *nexttetra = i;
                    break;
                  }
              }
        }
    }

  if(face >= 0 || edge >= 0)
    {
      if(s >= 1)
        {
          *pexit = *pend;
          return 0;
        }

      if(status >= 0)
        {
          u = xstart[0] + s * (xend[0] - xstart[0]);
          v = xstart[1] + s * (xend[1] - xstart[1]);
          w = xstart[2] + s * (xend[2] - xstart[2]);

#ifndef OPTIMIZE_MEMORY_USAGE
          pexit->xx = u * ax + v * bx + w * cx + p0->xx;
          pexit->yy = u * ay + v * by + w * cy + p0->yy;
          pexit->zz = u * az + v * bz + w * cz + p0->zz;
          pexit->x  = (pexit->xx - 1.0) / ConversionFac + CentralOffsetX;
          pexit->y  = (pexit->yy - 1.0) / ConversionFac + CentralOffsetY;
          pexit->z  = (pexit->zz - 1.0) / ConversionFac + CentralOffsetZ;
#else
          double pexit_xyz[3];
          pexit_xyz[0] = u * ax + v * bx + w * cx + pA_xyz[0];
          pexit_xyz[1] = u * ay + v * by + w * cy + pA_xyz[1];
          pexit_xyz[2] = u * az + v * bz + w * cz + pA_xyz[2];

          pexit->x = (pexit_xyz[0] - 1.0) / ConversionFac + CentralOffsetX;
          pexit->y = (pexit_xyz[1] - 1.0) / ConversionFac + CentralOffsetY;
          pexit->z = (pexit_xyz[2] - 1.0) / ConversionFac + CentralOffsetZ;
#endif
        }
    }

  if(*nexttetra == *previous_tetra)
    {
      printf("orientations: %d %d %d %d %d %d\n", orientations[0], orientations[1], orientations[2], orientations[3], orientations[4],
             orientations[5]);

      printf("face-cuts:  %d\n", face);
      printf("edge-cuts:  %d\n", edge);
      printf("corner-cuts: %d\n", corner);

      printf("\n");
      printf("p0:     %g %g %g\n", p0->x - pstart->x, p0->y - pstart->y, p0->z - pstart->z);
      printf("p1:     %g %g %g\n", p1->x - pstart->x, p1->y - pstart->y, p1->z - pstart->z);
      printf("p2:     %g %g %g\n", p2->x - pstart->x, p2->y - pstart->y, p2->z - pstart->z);
      printf("p3:     %g %g %g\n", p3->x - pstart->x, p3->y - pstart->y, p3->z - pstart->z);
      printf("\n");
      printf("pstart: %g %g %g\n", pstart->x, pstart->y, pstart->z);
      printf("pend:   %g %g %g\n", pend->x, pend->y, pend->z);

      terminate("nexttetra == previous_tetra");
    }

  return 1;
}

/* if this is defined, the definition is in voronoi_makeimage_new.c */
#ifndef VORONOI_NEW_IMAGE
/*! \brief Create a 3d projected image
 *
 *  \param[in] num Index in output file
 *  \param[in] gradients_flag Flag if gradients should be used in projection
 *  \param[in] pixels_x Number of pixels in x direction
 *  \param[in] pixels_y Number of pixels in y direction
 *  \param[in] xaxis Projection x axis in simulation coordinate system
 *  \param[in] yaxis Projection y axis in simulation coordinate system
 *  \param[in] zaxis Projection z axis in simulation coordinate system
 *  \param[in] xmin Minimum x coordinate in projection coordinate system
 *  \param[in] xmax Maximum x coordinate in projection coordinate system
 *  \param[in] ymin Minimum y coordinate in projection coordinate system
 *  \param[in] ymax Maximum y coordinate in projection coordinate system
 *  \param[in] zmin Minimum z coordinate in projection coordinate system
 *  \param[in] zmax Maximum z coordinate in projection coordinate system
 *  \param[in] weight_flag Flag for weighting; density^weight_flag
 */
void make_3d_voronoi_projected_image(int num, int gradients_flag, int pixels_x, int pixels_y, int xaxis, int yaxis, int zaxis,
                                     double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int weight_flag)
{
  CPU_Step[CPU_MISC] += measure_time();
  tessellation *T = &Mesh;

  char buf[MAXLEN_PATH];
  float *dens, *denssum, *dp, *temp, *tempsum, *tp, *weight, *weightsum, *wp;
#if defined(COOLING) && !defined(GRACKLE)
  float *szy, *szysum, *sp;
#endif
#ifdef GFM_STELLAR_EVOLUTION
  float *metal, *metalsum, *mp;
#endif
#ifdef GFM_AGN_RADIATION
  float *agnbol, *agnbolsum, *ap;
#endif
#ifdef TRACER_MC
  float *tracernum, *tracernumsum, *trp, *trweight, *trweightsum, *trwp;
#endif
#ifdef CHEM_IMAGE
  float *dust, *xH2, *xHP, *xCO, *dustsum, *h2sum, *hpsum, *cosum, *dst, *dh2, *dhp, *dco;
#endif
#ifdef MRT
  float *photon, *photonsum, *pp;
  float *h1, *h1sum, *h1p;
#endif

#ifdef SX_OUTPUT_IMAGE
  float *rih, *hrih;
  float *rihsum, *hrihsum;
  float *drih, *dhrih;
#endif

  FILE *fd = 0;
  point p, pold, pnew, pend, pstart;
  int tt0, tt, ttstart, ttrow, next_tetra, previous_tetra;
  double sigma, sigmatemp, sigmaweight;
#if defined(COOLING) && !defined(GRACKLE)
  double sigmaszy;
#endif
#ifdef GFM_STELLAR_EVOLUTION
  double sigmametal;
#endif
#ifdef GFM_AGN_RADIATION
  double sigmaagnbol;
#endif
#ifdef TRACER_MC
  double sigmatracernum, sigmatrweight;
#endif
#ifdef CHEM_IMAGE
  double sigmadust, sigmah2, sigmahp, sigmaco;
#endif
#ifdef MRT
  double sigmaphoton;
  double sigmah1;
#endif

#ifdef SX_OUTPUT_IMAGE
  double sigmarih, sigmahrih;
#endif
  int i, j, moves, ret, count_moves;

  mpi_printf(
      "Warning: Generating projected image without VORONOI_NEW_IMAGE.\nThe resulting column densities are known to be incorrect.\n");
  fflush(stdout);

  mpi_printf("we start to project the image... gradients_flag=%d\n", gradients_flag);

  T->DTF    = (char *)mymalloc_movable(&T->DTF, "DTF", T->MaxNdt * sizeof(char));
  tetra *DT = T->DT;

  for(i = 0; i < T->Ndt; i++)
    T->DTF[i] = 0;

  compute_circumcircles(T);

  /* check that we have enough space in DP for the points we
   * need. Note that we do NOT increase Ndp -- this means that this
   * code is not reentrant, though we can pass separate DP arrays to
   * different threads */
  while((T->Ndp >= (T->MaxNdp - 5)))
    {
      T->Indi.AllocFacNdp *= ALLOC_INCREASE_FACTOR;
      T->MaxNdp = T->Indi.AllocFacNdp;
#ifdef VERBOSE
      printf("Task=%d: increase memory allocation, MaxNdp=%d Indi.AllocFacNdp=%g\n", ThisTask, T->MaxNdp, T->Indi.AllocFacNdp);
#endif
      T->DP -= 5;
      T->DP = (point *)myrealloc_movable(T->DP, (T->MaxNdp + 5) * sizeof(point));
      T->DP += 5;

      if(T->Ndp >= (T->MaxNdp - 5) && NumGas == 0)
        terminate("(Ndp >= (MaxNdp - 5)");
    }

#ifndef VORONOI_MULTIPLE_PROJECTIONS
  if(gradients_flag == 1)
    file_path_sprintf(buf, "%s/proj_density_field_%03d", All.OutputDir, num);
  else if(gradients_flag == 0)
    file_path_sprintf(buf, "%s/proj_density_field_nograds_%03d", All.OutputDir, num);
#else
  if(gradients_flag == 1)
    file_path_sprintf(buf, "%s/proj_density_field_x%d_y%d_%03d", All.OutputDir, num, xaxis, yaxis);
  else if(gradients_flag == 0)
    file_path_sprintf(buf, "%s/proj_density_field_nograds_x%d_y%d_%03d", All.OutputDir, num, xaxis, yaxis);
#endif
  else
    terminate("gradients_flag != 0 && gradients_flag != 1");

  if(ThisTask == 0)
    {
      if(!(fd = fopen(buf, "w")))
        terminate("can't open file `%s' for writing snapshot.\n", buf);

      my_fwrite(&pixels_x, sizeof(int), 1, fd);
      my_fwrite(&pixels_y, sizeof(int), 1, fd);
    }

  /* allocate images and zero them out */
  dens    = (float *)mymalloc("dens", pixels_x * pixels_y * sizeof(float));
  denssum = (float *)mymalloc("denssum", pixels_x * pixels_y * sizeof(float));

  temp    = (float *)mymalloc("temp", pixels_x * pixels_y * sizeof(float));
  tempsum = (float *)mymalloc("tempsum", pixels_x * pixels_y * sizeof(float));

  weight    = (float *)mymalloc("weight", pixels_x * pixels_y * sizeof(float));
  weightsum = (float *)mymalloc("weightsum", pixels_x * pixels_y * sizeof(float));

#if defined(COOLING) && !defined(GRACKLE)
  szy    = (float *)mymalloc("szy", pixels_x * pixels_y * sizeof(float));
  szysum = (float *)mymalloc("szysum", pixels_x * pixels_y * sizeof(float));
#endif

#ifdef GFM_STELLAR_EVOLUTION
  metal    = (float *)mymalloc("metal", pixels_x * pixels_y * sizeof(float));
  metalsum = (float *)mymalloc("metalsum", pixels_x * pixels_y * sizeof(float));
#endif

#ifdef GFM_AGN_RADIATION
  agnbol    = (float *)mymalloc("agnbol", pixels_x * pixels_y * sizeof(float));
  agnbolsum = (float *)mymalloc("agnbolsum", pixels_x * pixels_y * sizeof(float));
#endif

#ifdef TRACER_MC
  tracernum    = (float *)mymalloc("tracernum", pixels_x * pixels_y * sizeof(float));
  tracernumsum = (float *)mymalloc("tracernumsum", pixels_x * pixels_y * sizeof(float));

  trweight    = (float *)mymalloc("trweight", pixels_x * pixels_y * sizeof(float));
  trweightsum = (float *)mymalloc("trweightsum", pixels_x * pixels_y * sizeof(float));
#endif

#ifdef CHEM_IMAGE
  dust    = (float *)mymalloc("dust", pixels_x * pixels_y * sizeof(float));
  dustsum = (float *)mymalloc("dustsum", pixels_x * pixels_y * sizeof(float));

  xH2   = (float *)mymalloc("xH2", pixels_x * pixels_y * sizeof(float));
  h2sum = (float *)mymalloc("h2sum", pixels_x * pixels_y * sizeof(float));

  xHP   = (float *)mymalloc("xHP", pixels_x * pixels_y * sizeof(float));
  hpsum = (float *)mymalloc("hpsum", pixels_x * pixels_y * sizeof(float));

  xCO   = (float *)mymalloc("xCO", pixels_x * pixels_y * sizeof(float));
  cosum = (float *)mymalloc("cosum", pixels_x * pixels_y * sizeof(float));
#endif

#ifdef SX_OUTPUT_IMAGE
  rih     = (float *)mymalloc("rih", pixels_x * pixels_y * sizeof(float));
  rihsum  = (float *)mymalloc("rihsum", pixels_x * pixels_y * sizeof(float));
  hrih    = (float *)mymalloc("hrih", pixels_x * pixels_y * sizeof(float));
  hrihsum = (float *)mymalloc("hrihsum", pixels_x * pixels_y * sizeof(float));
#endif

#ifdef MRT
  photon    = (float *)mymalloc("photon", pixels_x * pixels_y * sizeof(float));
  photonsum = (float *)mymalloc("photonsum", pixels_x * pixels_y * sizeof(float));
  h1        = (float *)mymalloc("h1", pixels_x * pixels_y * sizeof(float));
  h1sum     = (float *)mymalloc("h1sum", pixels_x * pixels_y * sizeof(float));
#endif

  for(i = 0, dp = dens; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *dp++ = 0;

  for(i = 0, tp = temp; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *tp++ = 0;

  for(i = 0, wp = weight; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *wp++ = 0;

#if defined(COOLING) && !defined(GRACKLE)
  for(i = 0, sp = szy; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *sp++ = 0;
#endif

#ifdef GFM_STELLAR_EVOLUTION
  for(i = 0, mp = metal; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *mp++ = 0;
#endif

#ifdef GFM_AGN_RADIATION
  for(i = 0, ap = agnbol; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *ap++ = 0;
#endif

#ifdef TRACER_MC
  for(i = 0, trp = tracernum; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *trp++ = 0;

  for(i = 0, trwp = trweight; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *trwp++ = 0;
#endif

#ifdef CHEM_IMAGE
  for(i = 0, dst = dust; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *dst++ = 0;

  for(i = 0, dh2 = xH2; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *dh2++ = 0;

  for(i = 0, dhp = xHP; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *dhp++ = 0;

  for(i = 0, dco = xCO; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *dco++ = 0;
#endif

#ifdef SX_OUTPUT_IMAGE
  for(i = 0, drih = rih; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *drih++ = 0;

  for(i = 0, dhrih = hrih; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *dhrih++ = 0;
#endif

#ifdef MRT
  for(i = 0, pp = photon; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *pp++ = 0;

  for(i = 0, h1p = h1; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *h1p++ = 0;
#endif

  // Get a suitable start tetrahedron?
  ttrow = 0;
  while(DT[ttrow].t[0] < 0 || DT[ttrow].p[0] == DPinfinity || DT[ttrow].p[1] == DPinfinity || DT[ttrow].p[2] == DPinfinity ||
        DT[ttrow].p[3] == DPinfinity)
    ttrow++;

  // Loop over x (columns)
  for(i = 0; i < pixels_x; i++)
    {
      ttstart = ttrow;

      // Loop over y (rows)
      for(j = 0; j < pixels_y; j++)
        {
          // The logic below sets the points p and pend correctly
          // depending on which axes we have chosen.
          if(xaxis == 0 && yaxis == 1)
            {
              p.x    = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
              p.y    = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
              p.z    = zmin;
              pend.x = p.x;
              pend.y = p.y;
              pend.z = zmax;
            }
          else if(xaxis == 1 && yaxis == 0)
            {
              p.x    = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
              p.y    = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
              p.z    = zmin;
              pend.x = p.x;
              pend.y = p.y;
              pend.z = zmax;
            }
          else if(xaxis == 0 && yaxis == 2)
            {
              p.x    = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
              p.y    = zmin;
              p.z    = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
              pend.x = p.x;
              pend.y = zmax;
              pend.z = p.z;
            }
          else if(xaxis == 2 && yaxis == 0)
            {
              p.x    = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
              p.y    = zmin;
              p.z    = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
              pend.x = p.x;
              pend.y = zmax;
              pend.z = p.z;
            }
          else if(xaxis == 1 && yaxis == 2)
            {
              p.x    = zmin;
              p.y    = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
              p.z    = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
              pend.x = zmax;
              pend.y = p.y;
              pend.z = p.z;
            }
          else if(xaxis == 2 && yaxis == 1)
            {
              p.x    = zmin;
              p.y    = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
              p.z    = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
              pend.x = zmax;
              pend.y = p.y;
              pend.z = p.z;
            }
          else
            terminate("invalid combination of axes");
#ifndef OPTIMIZE_MEMORY_USAGE
          set_integers_for_pointer(&p);
          set_integers_for_pointer(&pend);
#endif
          if(ttstart < 0)
            terminate("ttstart < 0");

          if(DT[ttstart].t[0] < 0)
            {
              printf("this start tetra (ttstart=%d) was deleted\n", ttstart);
              terminate("deleted start tetrahedron");
            }

          // Find first tetra that contains point pp
          tt0 = get_tetra(T, &p, &moves, ttstart, &ret, &ret);

          // Update starting tetras in both row and column, if applicable
          if(j == 0)
            ttrow = tt0;

          ttstart = tt0;

          tt = tt0;

          pstart         = p; /* this is the starting point */
          pold           = p;
          previous_tetra = -1;

          count_moves = 0;

          // Now walk through the tetras along the ray. ppnew is set
          // to the exit point from this/entry point to next tetra
          while((ret = image_get_next_tetra(T, tt, &pstart, &pend, &next_tetra, &pnew, &previous_tetra)))
            {
#ifndef OPTIMIZE_MEMORY_USAGE
              set_integers_for_pointer(&pold);
              set_integers_for_pointer(&pnew);
#endif
              if(count_moves > MAX_COUNT_MOVES - 10)
                {
                  printf("i/j=(%d|%d)  tetra=%d next=%d   x=%g y=%g z=%g\n", i, j, (int)(tt), (int)(next_tetra), pnew.x, pnew.y,
                         pnew.z);
                }

              count_moves++;

              if(count_moves > MAX_COUNT_MOVES)
                terminate("count_moves > MAX_COUNT_MOVES");

              // calculate the contribution to this pixel for this tetrahedron
              calc_picture_contribution(T, tt, &pold, &pnew, &sigma, &sigmatemp, &sigmaweight, weight_flag, gradients_flag
#if defined(COOLING) && !defined(GRACKLE)
                                        ,
                                        &sigmaszy
#endif
#ifdef GFM_STELLAR_EVOLUTION
                                        ,
                                        &sigmametal
#endif
#ifdef GFM_AGN_RADIATION
                                        ,
                                        &sigmaagnbol
#endif
#ifdef TRACER_MC
                                        ,
                                        &sigmatracernum, &sigmatrweight
#endif
#ifdef CHEM_IMAGE
                                        ,
                                        &sigmadust, &sigmah2, &sigmahp, &sigmaco
#endif
#ifdef MRT
                                        ,
                                        &sigmaphoton, &sigmah1
#endif

#ifdef SX_OUTPUT_IMAGE
                                        ,
                                        &sigmarih, &sigmahrih
#endif
              );

              dens[i * pixels_y + j] += sigma;
              temp[i * pixels_y + j] += sigmatemp;
              weight[i * pixels_y + j] += sigmaweight;
#if defined(COOLING) && !defined(GRACKLE)
              szy[i * pixels_y + j] += sigmaszy * BOLTZMANN * THOMPSON / (ELECTRONMASS * CLIGHT * CLIGHT);
#endif
#ifdef GFM_STELLAR_EVOLUTION
              metal[i * pixels_y + j] += sigmametal;
#endif
#ifdef GFM_AGN_RADIATION
              agnbol[i * pixels_y + j] += sigmaagnbol;
#endif
#ifdef TRACER_MC
              tracernum[i * pixels_y + j] += sigmatracernum;
              trweight[i * pixels_y + j] += sigmatrweight;
#endif
#ifdef CHEM_IMAGE
              dust[i * pixels_y + j] += sigmadust;
              xH2[i * pixels_y + j] += sigmah2;
              xHP[i * pixels_y + j] += sigmahp;
              xCO[i * pixels_y + j] += sigmaco;
#endif
#ifdef MRT
              photon[i * pixels_y + j] += sigmaphoton;
              h1[i * pixels_y + j] += sigmah1;
#endif

#ifdef SX_OUTPUT_IMAGE
              rih[i * pixels_y + j] += sigmarih;
              hrih[i * pixels_y + j] += sigmahrih;
#endif
              /* update "old" values */
              previous_tetra = tt;
              tt             = next_tetra;

              pold = pnew;
            }

          calc_picture_contribution(T, tt, &pold, &pnew, &sigma, &sigmatemp, &sigmaweight, weight_flag, gradients_flag
#if defined(COOLING) && !defined(GRACKLE)
                                    ,
                                    &sigmaszy
#endif
#ifdef GFM_STELLAR_EVOLUTION
                                    ,
                                    &sigmametal
#endif
#ifdef GFM_AGN_RADIATION
                                    ,
                                    &sigmaagnbol
#endif
#ifdef TRACER_MC
                                    ,
                                    &sigmatracernum, &sigmatrweight
#endif
#ifdef CHEM_IMAGE
                                    ,
                                    &sigmadust, &sigmah2, &sigmahp, &sigmaco
#endif
#ifdef MRT
                                    ,
                                    &sigmaphoton, &sigmah1
#endif

#ifdef SX_OUTPUT_IMAGE
                                    ,
                                    &sigmarih, &sigmahrih
#endif
          );

          dens[i * pixels_y + j] += sigma;
          temp[i * pixels_y + j] += sigmatemp;
          weight[i * pixels_y + j] += sigmaweight;
#if defined(COOLING) && !defined(GRACKLE)
          szy[i * pixels_y + j] += sigmaszy * BOLTZMANN * THOMPSON / (ELECTRONMASS * CLIGHT * CLIGHT);
#endif
#ifdef GFM_STELLAR_EVOLUTION
          metal[i * pixels_y + j] += sigmametal;
#endif
#ifdef GFM_AGN_RADIATION
          agnbol[i * pixels_y + j] += sigmaagnbol;
#endif
#ifdef TRACER_MC
          tracernum[i * pixels_y + j] += sigmatracernum;
          trweight[i * pixels_y + j] += sigmatrweight;
#endif
#ifdef CHEM_IMAGE
          dust[i * pixels_y + j] += sigmadust;
          xH2[i * pixels_y + j] += sigmah2;
          xHP[i * pixels_y + j] += sigmahp;
          xCO[i * pixels_y + j] += sigmaco;
#endif
#ifdef MRT
          photon[i * pixels_y + j] += sigmaphoton;
          h1[i * pixels_y + j] += sigmah1;
#endif

#ifdef SX_OUTPUT_IMAGE
          rih[i * pixels_y + j] += sigmarih;
          hrih[i * pixels_y + j] += sigmahrih;
#endif
        }
    }

  MPI_Reduce(dens, denssum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(temp, tempsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(weight, weightsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#if defined(COOLING) && !defined(GRACKLE)
  MPI_Reduce(szy, szysum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
#ifdef GFM_STELLAR_EVOLUTION
  MPI_Reduce(metal, metalsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
#ifdef GFM_AGN_RADIATION
  MPI_Reduce(agnbol, agnbolsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
#ifdef TRACER_MC
  MPI_Reduce(tracernum, tracernumsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(trweight, trweightsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
#ifdef CHEM_IMAGE
  MPI_Reduce(dust, dustsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(xH2, h2sum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(xHP, hpsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(xCO, cosum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

#ifdef SX_OUTPUT_IMAGE
  MPI_Reduce(rih, rihsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(hrih, hrihsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

#ifdef MRT
  MPI_Reduce(photon, photonsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(h1, h1sum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

  if(ThisTask == 0)
    {
      for(i = 0, wp = weightsum, tp = tempsum, dp = denssum
#if defined(COOLING) && !defined(GRACKLE)
          ,
      sp = szysum
#endif
#ifdef GFM_STELLAR_EVOLUTION
          ,
      mp = metalsum
#endif
#ifdef GFM_AGN_RADIATION
          ,
      ap = agnbolsum
#endif
#ifdef TRACER_MC
          ,
      trp = tracernumsum, trwp = trweightsum
#endif
#ifdef CHEM_IMAGE
          ,
      dst = dustsum, dh2 = h2sum, dhp = hpsum, dco = cosum
#endif
#ifdef MRT
          ,
      pp = photonsum, h1p = h1sum
#endif

#ifdef SX_OUTPUT_IMAGE
          ,
      drih = rihsum, dhrih = hrihsum
#endif
          ;
          i < pixels_x; i++)

        for(j = 0; j < pixels_y; j++, dp++, tp++, wp++
#if defined(COOLING) && !defined(GRACKLE)
                                 ,
        sp++
#endif
#ifdef GFM_STELLAR_EVOLUTION
                                 ,
        mp++
#endif
#ifdef GFM_AGN_RADIATION
                                 ,
        ap++
#endif
#ifdef TRACER_MC
                                 ,
        trp++, trwp++
#endif
#ifdef CHEM_IMAGE
                                 ,
        dst++, dh2++, dhp++, dco++
#endif
#ifdef MRT
                                 ,
        pp++, h1p++
#endif

#ifdef SX_OUTPUT_IMAGE
                                 ,
        drih++, dhrih++
#endif
        )
          {
            if(*dp > 0)
              {
                *tp /= *dp;
              }

#ifdef CHEM_IMAGE
            if(*dp > 0)
              {
                *dst /= *dp;
                *dh2 /= *dp;
                *dhp /= *dp;
                *dco /= *dp;
              }
#endif
#ifdef MRT
            if(*dp > 0)
              {
                *pp /= *dp;
                *h1p /= *dp;
              }
#endif

#ifdef SX_OUTPUT_IMAGE
            if(*dp > 0)
              {
                *drih /= *dp;
                *dhrih /= *dp;
              }
#endif

            if(*wp > 0)
              *dp /= *wp;

#ifdef GFM_STELLAR_EVOLUTION
            if(*dp > 0)
              *mp /= *dp;
#endif
#ifdef GFM_AGN_RADIATION
            if(*dp > 0)
              *ap /= *dp;
#endif
#ifdef TRACER_MC
            if(*trwp > 0)
              *trp /= *trwp;
#endif

            if(weight_flag > 0)
              {
                *dp *= (zmax - zmin);
#ifdef TRACER_MC
                *trp *= (zmax - zmin);
#endif
              }
          }

      my_fwrite(denssum, sizeof(float), pixels_x * pixels_y, fd);
      my_fwrite(tempsum, sizeof(float), pixels_x * pixels_y, fd);
#if defined(COOLING) && !defined(GRACKLE)
      my_fwrite(szysum, sizeof(float), pixels_x * pixels_y, fd);
#endif
#ifdef GFM_STELLAR_EVOLUTION
      my_fwrite(metalsum, sizeof(float), pixels_x * pixels_y, fd);
#endif
#ifdef GFM_AGN_RADIATION
      my_fwrite(agnbolsum, sizeof(float), pixels_x * pixels_y, fd);
#endif
#ifdef TRACER_MC
      my_fwrite(tracernumsum, sizeof(float), pixels_x * pixels_y, fd);
#endif
#ifdef CHEM_IMAGE
      my_fwrite(dustsum, sizeof(float), pixels_x * pixels_y, fd);
      my_fwrite(h2sum, sizeof(float), pixels_x * pixels_y, fd);
      my_fwrite(hpsum, sizeof(float), pixels_x * pixels_y, fd);
      my_fwrite(cosum, sizeof(float), pixels_x * pixels_y, fd);
#endif
#ifdef MRT
      my_fwrite(photonsum, sizeof(float), pixels_x * pixels_y, fd);
      my_fwrite(h1sum, sizeof(float), pixels_x * pixels_y, fd);
#endif

#ifdef SX_OUTPUT_IMAGE
      my_fwrite(rihsum, sizeof(float), pixels_x * pixels_y, fd);
      my_fwrite(hrihsum, sizeof(float), pixels_x * pixels_y, fd);
#endif
      fclose(fd);
    }

#ifdef MRT
  myfree(h1sum);
  myfree(h1);
  myfree(photonsum);
  myfree(photon);
#endif

#ifdef SX_OUTPUT_IMAGE
  myfree(hrihsum);
  myfree(hrih);
  myfree(rihsum);
  myfree(rih);
#endif

#ifdef CHEM_IMAGE
  myfree(cosum);
  myfree(xCO);
  myfree(hpsum);
  myfree(xHP);
  myfree(h2sum);
  myfree(xH2);
  myfree(dustsum);
  myfree(dust);
#endif
#ifdef TRACER_MC
  myfree(trweightsum);
  myfree(trweight);
  myfree(tracernumsum);
  myfree(tracernum);
#endif
#ifdef GFM_AGN_RADIATION
  myfree(agnbolsum);
  myfree(agnbol);
#endif
#ifdef GFM_STELLAR_EVOLUTION
  myfree(metalsum);
  myfree(metal);
#endif
#if defined(COOLING) && !defined(GRACKLE)
  myfree(szysum);
  myfree(szy);
#endif
  myfree(weightsum);
  myfree(weight);
  myfree(tempsum);
  myfree(temp);
  myfree(denssum);
  myfree(dens);

  myfree(T->DTF);

  CPU_Step[CPU_MAKEIMAGES] += measure_time();
}
#endif

/* this is the code for simple tetrahedra interpolation */
#ifdef SIMPLE_TETRAHEDRA_PROJECTION
/*! \brief Simple approximation to the contribution of a cell to a projection
 *
 *  \param[in] T Pointer to tessellation
 *  \param[in] tt Index in DT
 *  \param[in] p0 Array with coordinates of start point
 *  \param[in] p1 Array with coordinates of end point
 *  \param[out] sigma Output array for surface density
 *  \param[out] sigmatemp Output array for temperature
 *  \param[out] sigmaweight Output array for weight
 *  \param[in] weight_flag Flag for projection weighting: rho^weight_flag
 *  \param[in] gradients_flag Flag whether gradients are used in projection
 *  \param[out] sigmaszy Output array for He abunance (not functional)
 *  \param[out] sigmatracernum Output array for number of tracers
 *  \param[out] sigmatrweight Output arrray for weighted tracer sum
 *  \param[out] sigmarih Output array
 *  \param[out] sigmahrih Output array
 *
 *  \return void
 */
void calc_picture_contribution(tessellation *T, tetra *t, point *p0, point *p1, double *sigma, double *sigmatemp, double *sigmaweight,
                               int weight_flag, int gradients_flag
#if defined(COOLING) && !defined(GRACKLE)
                               ,
                               double *sigmaszy
#endif
)
{
  double rho = 0;
  int k, li;

  for(k = 0; k < 4; k++)
    {
      if(t->p[k]->task == ThisTask)
        {
          li = t->p[k]->index;

          if(li >= NumGas)
            li -= NumGas;

          rho += SphP[li].Density;
        }
    }

  rho /= 4;

  double dx = p0->x - p1->x;
  double dy = p0->y - p1->y;
  double dz = p0->z - p1->z;

  double r = sqrt(dx * dx + dy * dy + dz * dz);

  if(weight_flag > 0)
    {
      *sigma       = pow(rho, weight_flag) * r;
      *sigmaweight = pow(rho, weight_flag - 1) * r;
    }
  else
    {
      *sigma       = rho * r;
      *sigmaweight = 0;
    }

  *sigmatemp = 0;
#if defined(COOLING) && !defined(GRACKLE)
  *sigmaszy = 0;
#endif
}
#else

/** These 6 permutations of 0,1,2,3 define the intersection tests
    below. The first two numbers are the 6 edges. The final two are
    just the two remaining points that must be tested against. */
const char pairs[6][4] = {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2}, {1, 2, 0, 3}, {1, 3, 0, 2}, {2, 3, 0, 1}};

/*! \brief Checks that the specified position is on the inside of all the
 *  face planes that make up the specified cell.
 *
 *  \param[in] T Pointer to tessellation
 *  \param[in] cell Index of cell in P and SphP array
 *  \param[in] p0 Position to be checked
 *
 *  \return minimum distance to a face
 */
double assert_contains(tessellation *T, int cell, MyDouble p0[3])
{
#ifdef VORONOI_DYNAMIC_UPDATE
  point *DP = T->DP;

  MyDouble cell_p[3];
  cell_p[0] = P[cell].Pos[0];
  cell_p[1] = P[cell].Pos[1];
  cell_p[2] = P[cell].Pos[2];

  // if cell center is across the boundary, wrap it
  periodic_wrap_point_MyDouble(cell_p, p0);

  MyDouble nb_p[3];
  double m[3];
  double c[3];
  double q[3];

  int edge = SphP[cell].first_connection;
  int last_edge = SphP[cell].last_connection;

  double mindist = MAX_DOUBLE_NUMBER;

  while(1)
    {
      const int neighbor = DC[edge].dp_index;
      // myassert((DC[edge].task != ThisTask) || (DC[edge].index != cell));   // this assert is not fulfilled for reflective boundary
      // conditions

      nb_p[0] = DP[neighbor].x;
      nb_p[1] = DP[neighbor].y;
      nb_p[2] = DP[neighbor].z;
      // if neighbor is across the boundary, wrap it
      periodic_wrap_point_MyDouble(nb_p, p0);

      int i;
      for(i = 0; i < 3; ++i)
        {
          // m is the edge midpoint, which is a point on the face
          m[i] = 0.5 * (nb_p[i] + cell_p[i]);
          // c is the vector from entry point p0 to the point on the face
          c[i] = m[i] - p0[i];
          /* q is the edge vector to the neighboring cell, which is a
             normal vector of the face plane, pointing outward. */
          q[i] = nb_p[i] - cell_p[i];
        }

      double cdotq = c[0] * q[0] + c[1] * q[1] + c[2] * q[2];
      // distance from point to face, counting POSITIVE INSIDE. (If
      // q and c point in the same direction, point is inside and
      // cdotq is positive.)
      double dist_to_face = cdotq / sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
      if(dist_to_face < mindist)
        mindist = dist_to_face;

      if(edge == last_edge)
        break;

      edge = DC[edge].next;
    }
  return mindist;
#else
  myassert(0);
  return 0;
#endif
}

/*! \brief Calculates the intersections between a ray and the internal Voronoi
 *  faces in a Delaunay tetrahedron.
 *
 *  The intersections are returned in
 *  list (which must be an array of 8 elements), with nlist indicating
 *  how many intersections there are. (There are always at least 2, the
 *  start and end points.)
 *
 *  \param[in] T Pointer to tessellation
 *  \param[in] tt Index in DT
 *  \param[in] p0 Starting point
 *  \param[in] p1 End point
 *  \param[out] list Intersection list
 *  \param[out] nlist Number of elements in Intersection list
 *
 *  \return void
 */
void calc_delaunay_intersections(tessellation *T, int tt, point *p0, point *p1, intersection_list *list, int *nlist)
{
  tetra *DT = T->DT;
  point *DP = T->DP;
  tetra_center *DTC = T->DTC;
  tetra *t = &DT[tt];
  tetra_center *tc = &DTC[tt];

  int i, j, k, l;
  double qx, qy, qz;
  double dx, dy, dz;
  double cx, cy, cz;
  double px, py, pz;
  double rx, ry, rz;
  double s, dist2, dist2max;

  *nlist = 2;

  // initialize list[0] to be the starting point p0. indA and indB are
  // both set to the index (0-4) of the closest point. s=0.
  list[0].s = 0;
  list[0].p.x = p0->x;
  list[0].p.y = p0->y;
  list[0].p.z = p0->z;

  for(i = 0, dist2max = 1.0e30; i < 4; i++)
    {
      dx = p0->x - DP[t->p[i]].x;
      dy = p0->y - DP[t->p[i]].y;
      dz = p0->z - DP[t->p[i]].z;

      dist2 = dx * dx + dy * dy + dz * dz;
      if(dist2 < dist2max)
        {
          dist2max = dist2;
          list[0].indA = list[0].indB = i;
        }
    }

  // initialize list[1] to be the end point p1. indA and indB are
  // both set to the index (0-4) of the closest point.
  list[1].s = 1.0;
  list[1].p.x = p1->x;
  list[1].p.y = p1->y;
  list[1].p.z = p1->z;

  for(i = 0, dist2max = 1.0e30; i < 4; i++)
    {
      dx = p1->x - DP[t->p[i]].x;
      dy = p1->y - DP[t->p[i]].y;
      dz = p1->z - DP[t->p[i]].z;

      dist2 = dx * dx + dy * dy + dz * dz;
      if(dist2 < dist2max)
        {
          dist2max = dist2;
          list[1].indA = list[1].indB = i;
        }
    }

  // c is the vector from the tetra center to tetra point 0
  // XXX CLOBBERED BELOW!
  cx = tc->cx - DP[t->p[0]].x;
  cy = tc->cy - DP[t->p[0]].y;
  cz = tc->cz - DP[t->p[0]].z;

  // d is the vector from entry to exit point (the full line segment)
  dx = p1->x - p0->x;
  dy = p1->y - p0->y;
  dz = p1->z - p0->z;

  // c is the vector from entry point p0 to tetra center
  cx = tc->cx - p0->x;
  cy = tc->cy - p0->y;
  cz = tc->cz - p0->z;

  for(i = 0; i < 6; i++)
    {
      // pairs is a [6][4] const array defined above that defines 6 different
      // permutations of the tetrahedron points 0,1,2,3.
      j = pairs[i][0];
      k = pairs[i][1];

      point *dpj = &DP[t->p[j]];
      point *dpk = &DP[t->p[k]];

      // q is the tetra edge vector from pk to pj.
      qx = dpj->x - dpk->x;
      qy = dpj->y - dpk->y;
      qz = dpj->z - dpk->z;

      // s = c.q / d.q
      // This is the standard formula for the intersection between a
      // ray and a plane. s is the point where the ray p0+s*d
      // intersects the plane which is perpendicular to q and goes
      // through c, i.e. the Voronoi face corresponding to the j-k
      // edge.
      s = (cx * qx + cy * qy + cz * qz) / (dx * qx + dy * qy + dz * qz);

      if(s > 0 && s < 1)
        {
          /* The ray intersects the j-k Voronoi face inside the
             tetrahedron, which means we'll have at least one
             split. But we need to check whether it intersects another
             face before. */

          // set p to the intersection point
          px = p0->x + s * dx;
          py = p0->y + s * dy;
          pz = p0->z + s * dz;

          // r is the vector from pj to p
          rx = px - dpj->x;
          ry = py - dpj->y;
          rz = pz - dpj->z;

          /* Now we check whether p is really within voronoi cell l,
             by comparing the distance from p to l. Since Voronoi cell
             j is defined by the points closer to j than any other
             point. By construction p is on the j-k face, so is
             equally distant from j and k. We thus only need to check
             the distance from p to the two other points in the
             tetrahedron. */
          l = pairs[i][2];
          point *dpl = &DP[t->p[l]];

          // rr is the vector from pl to p
          double rrx, rry, rrz;
          rrx = px - dpl->x;
          rry = py - dpl->y;
          rrz = pz - dpl->z;

          if((rx * rx + ry * ry + rz * rz) < (rrx * rrx + rry * rry + rrz * rrz))
            {
              /* r is shorter than rr, the intersection point p is
                 closer to pj than pl, which means p can not be in pl's
                 voronoi cell. */

              // Now, in the same way, we check whether p is really
              // within voronoi cell m.
              int m = pairs[i][3];
              point *dpm = &DP[t->p[m]];

              rrx = px - dpm->x;
              rry = py - dpm->y;
              rrz = pz - dpm->z;

              if((rx * rx + ry * ry + rz * rz) < (rrx * rrx + rry * rry + rrz * rrz))
                {
                  /* p is indeed on the j-k voronoi face. We have
                     found a new cut-point. Perform an insertion sort
                     on the list. */

                  for(l = 0; l < *nlist - 1; l++)
                    {
                      // find the place where to insert the cut point
                      if(list[l].s < s && s < list[l + 1].s)
                        {
                          // make space in the list and add the new element
                          memmove(&list[l + 2], &list[l + 1], (*nlist - (l + 1)) * sizeof(intersection_list));

                          list[l + 1].s = s;
                          list[l + 1].p.x = px;
                          list[l + 1].p.y = py;
                          list[l + 1].p.z = pz;
                          list[l + 1].indA = j;
                          list[l + 1].indB = k;
                          (*nlist)++;
                          break;
                        }
                    }
                }
            }
        }
    }
}

/*! \brief Decodes an intersection list
 *
 *  This little function returns the Voronoi cell that the line
 *  segment [i, i+1] in the intersection_list is internal to,
 *  expressed as the point in the tetra
 *
 *  \param[in] i Index in intersection list
 *  \param[in] list Array of intersection lists
 *
 *  \return Voronoi cell index that line segment is internal to
 */
int decode_intersection_list(int i, intersection_list *list)
{
  if(list[i].indA == list[i + 1].indA)
    return list[i].indA;
  else if(list[i].indA == list[i + 1].indB)
    return list[i].indA;
  else if(list[i].indB == list[i + 1].indA)
    return list[i].indB;
  else if(list[i].indB == list[i + 1].indB)
    return list[i].indB;
  else
    return -1;
}

/*! \brief get contribution of given tetarhedron for line segment
 *
 *  Calculates the surface density for tetrahedron tt for the line
 *  segment [pp0, pp1] (which is known to be inside this tetra). The
 *  strategy appears to be to split the line segment into segments
 *  that are within the 4 Voronoi cells that exist within the
 *  tetrahedron. For each of these segments, the density field is
 *  described by the density at the point (ie, Voronoi center) and the
 *  gradient. It is thus trivial to integrate the column density
 *  within that segment. Most of the logic goes into finding the
 *  splitting points.
 *
 *  \param[in] T Pointer to tessellation
 *  \param[in] tt Index in DT
 *  \param[in] p0 Array with coordinates of start point
 *  \param[in] p1 Array with coordinates of end point
 *  \param[out] sigma Output array for surface density
 *  \param[out] sigmatemp Output array for temperature
 *  \param[out] sigmaweight Output array for weight
 *  \param[in] weight_flag Flag for projection weighting: rho^weight_flag
 *  \param[in] gradients_flag Flag whether gradients are used in projection
 *  \param[out] sigmaszy Output array for He abunance (not functional)
 *  \param[out] sigmatracernum Output array for number of tracers
 *  \param[out] sigmatrweight Output arrray for weighted tracer sum
 *  \param[out] sigmarih Output array
 *  \param[out] sigmahrih Output array
 *
 *  \return void
 */
void calc_picture_contribution(tessellation *T, int tt, point *p0, point *p1, double *sigma, double *sigmatemp, double *sigmaweight,
                               int weight_flag, int gradients_flag
#if defined(COOLING) && !defined(GRACKLE)
                               ,
                               double *sigmaszy
#endif
#ifdef GFM_STELLAR_EVOLUTION
                               ,
                               double *sigmametal
#endif
#ifdef GFM_AGN_RADIATION
                               ,
                               double *sigmaagnbol
#endif
#ifdef TRACER_MC
                               ,
                               double *sigmatracernum, double *sigmatrweight
#endif
#ifdef CHEM_IMAGE
                               ,
                               double *sigmadust, double *sigmah2, double *sigmahp, double *sigmaco
#endif
#ifdef MRT
                               ,
                               double *sigmaphoton, double *sigmah1
#endif

#ifdef SX_OUTPUT_IMAGE
                               ,
                               double *sigmarih, double *sigmahrih
#endif
)
{
  tetra *DT = T->DT;
  point *DP = T->DP;

  int i, j, li;
  double dx, dy, dz;
  double mx, my, mz;

  tetra *t = &DT[tt];

  intersection_list list[8];
  int nlist;

  double rho, proj_rhotemp_sum = 0;

  calc_delaunay_intersections(T, tt, p0, p1, list, &nlist);

  double proj_sum = 0;
  double proj_weight_sum = 0;
#if defined(COOLING) && !defined(GRACKLE)
  double proj_szy_sum = 0;
#endif
#ifdef GFM_STELLAR_EVOLUTION
  double proj_rhometal_sum = 0;
#endif
#ifdef GFM_AGN_RADIATION
  double proj_rhoagnbol_sum = 0;
#endif
#ifdef TRACER_MC
  double proj_rhotracernum_sum = 0;
  double proj_trweight_sum = 0;
#endif
#ifdef CHEM_IMAGE
  double proj_dust_sum = 0;
  double proj_h2_sum = 0;
  double proj_hp_sum = 0;
  double proj_co_sum = 0;
#endif
#ifdef MRT
  double proj_photon_sum = 0;
  double proj_h1_sum = 0;
#endif

#ifdef SX_OUTPUT_IMAGE
  double proj_rih_sum = 0;
  double proj_hrih_sum = 0;
#endif

  /* now that all cuts have been established, we can integrate the
     column density for the segments */

  for(i = 0; i < nlist - 1; i++)
    {
      //      printf("(%d|%d)  s0=%g s1=%g  (%d|%d)  (%d|%d)\n", i, nlist, list[i].s, list[i+1].s, list[i].indA, list[i].indB,
      //      list[i+1].indA, list[i+1].indB);

      /* figure out which voronoi cell the segment is internal to */
      j = decode_intersection_list(i, list);

      /* not sure what could cause this to fail... */
      if(j >= 0)
        {
          /* check that this voronoi cell indeed is on this processor */
          if(DP[t->p[j]].task == ThisTask)
            {
              li = DP[t->p[j]].index;

              /* check that this point actually has a hydro quantity
               * (what are failure scenarios?) */
              if(li >= 0 && li < NumGas)
                {
                  /* now determine the interpolated density value at the
                     middle of the two points */

                  /* Is the SphP[i].Center different from DP[i] ?? Yes. */

                  mx = 0.5 * (list[i].p.x + list[i + 1].p.x) - SphP[li].Center[0];
                  my = 0.5 * (list[i].p.y + list[i + 1].p.y) - SphP[li].Center[1];
                  mz = 0.5 * (list[i].p.z + list[i + 1].p.z) - SphP[li].Center[2];

                  if(gradients_flag == 1)
                    rho = SphP[li].Density + SphP[li].Grad.drho[0] * mx + SphP[li].Grad.drho[1] * my + SphP[li].Grad.drho[2] * mz;
                  else if(gradients_flag == 0)
                    rho = SphP[li].Density;
                  else
                    terminate("gradients_flag != 0 && gradients_flag != 1");

                  /* get the column length (this seems suboptimal, we already have s) */
                  // r=(list[i+1].s-list[i].s)*l_seg;
                  dx = list[i].p.x - list[i + 1].p.x;
                  dy = list[i].p.y - list[i + 1].p.y;
                  dz = list[i].p.z - list[i + 1].p.z;

                  double r = sqrt(dx * dx + dy * dy + dz * dz);

                  if(weight_flag > 0)
                    {
                      double pow_rho = pow(rho, weight_flag), pow_rho_sum = pow(rho, weight_flag - 1);
                      proj_sum += pow_rho * r;
#ifdef VORONOI_PROJ_TEMP
                      double meanweight = 4.0 / (1. + 3. * HYDROGEN_MASSFRAC) * PROTONMASS;
#if defined(COOLING) && !defined(GRACKLE)
                      meanweight = 4.0 / (1. + 3. * HYDROGEN_MASSFRAC + 4. * HYDROGEN_MASSFRAC * SphP[li].Ne) * PROTONMASS;
#endif
#ifdef VARIABLE_GAMMA
                      double temp_in_K = (SphP[li].GammaE - 1.0) * SphP[li].Utherm / BOLTZMANN * All.UnitEnergy_in_cgs /
                                         All.UnitMass_in_g * meanweight;
#else
                      double temp_in_K =
                          GAMMA_MINUS1 * SphP[li].Utherm / BOLTZMANN * All.UnitEnergy_in_cgs / All.UnitMass_in_g * meanweight;
#endif
#ifdef SGCHEM
                      double yn, en, yntot;

                      yn = rho * All.UnitDensity_in_cgs / ((1.0 + 4.0 * ABHE) * PROTONMASS);
                      en = SphP[li].Utherm * rho * All.UnitEnergy_in_cgs / pow(All.UnitLength_in_cm, 3);
#if CHEMISTRYNETWORK == 1
                      yntot = (1.0 + ABHE - SphP[li].TracAbund[IH2] + SphP[li].TracAbund[IHP] + SphP[li].TracAbund[IHEP] +
                               SphP[li].TracAbund[IHEPP]) *
                              yn;
#else
                      yntot     = (1.0 + ABHE - SphP[li].TracAbund[IH2] + SphP[li].TracAbund[IHP]) * yn;
#endif
#ifdef VARIABLE_GAMMA
                      temp_in_K = (SphP[li].GammaE - 1.0) * en / (yntot * BOLTZMANN);
#else
                      temp_in_K = (GAMMA - 1.0) * en / (yntot * BOLTZMANN);
#endif
#endif /* SGCHEM */

                      proj_rhotemp_sum += temp_in_K * pow_rho * r;
#if defined(COOLING) && !defined(GRACKLE)
                      proj_szy_sum += temp_in_K * SphP[li].Ne * (rho * All.UnitDensity_in_cgs * HYDROGEN_MASSFRAC / PROTONMASS) *
                                      (r * All.UnitLength_in_cm);  // in pysical units
#endif
#else
                      proj_rhotemp_sum += SphP[li].Utherm * pow_rho * r;
#endif

#ifdef GFM_STELLAR_EVOLUTION
                      proj_rhometal_sum += SphP[li].Metallicity * pow_rho * r;
#endif
#ifdef GFM_AGN_RADIATION
                      proj_rhoagnbol_sum += SphP[li].AGNBolIntensity * pow_rho * r;
#endif
#ifdef TRACER_MC
                      proj_rhotracernum_sum +=
                          pow(get_number_of_tracers(li) * All.ReferenceTracerMCMass / P[li].Mass, 3.0) * pow_rho * r;
                      proj_trweight_sum +=
                          pow(get_number_of_tracers(li) * All.ReferenceTracerMCMass / P[li].Mass, 2.0) * pow_rho_sum * r;
#endif
#ifdef CHEM_IMAGE
                      proj_dust_sum += SphP[li].DustTemp * pow_rho * r;
                      proj_h2_sum += SphP[li].TracAbund[IH2] * pow_rho * r;
                      proj_hp_sum += SphP[li].TracAbund[IHP] * pow_rho * r;
#if CHEMISTRYNETWORK != 1
                      proj_co_sum += SphP[li].TracAbund[ICO] * pow_rho * r;
#else
                      proj_co_sum += 0.0;
#endif
#endif
#ifdef MRT
                      proj_photon_sum += SphP[li].DensPhot[0] * pow_rho * r;
                      proj_h1_sum += SphP[li].HI * pow_rho * r;
#endif

#ifdef SX_OUTPUT_IMAGE
                      proj_rih_sum += sx_chem_image1(li) * pow_rho * r;
                      proj_hrih_sum += sx_chem_image2(li) * pow_rho * r;
#endif

                      proj_weight_sum += pow_rho_sum * r;
                    }
                  else
                    {
                      proj_sum += rho * r;
#ifdef VORONOI_PROJ_TEMP
                      double meanweight = 4.0 / (1. + 3. * HYDROGEN_MASSFRAC) * PROTONMASS;
#if defined(COOLING) && !defined(GRACKLE)
                      meanweight = 4.0 / (1. + 3. * HYDROGEN_MASSFRAC + 4. * HYDROGEN_MASSFRAC * SphP[li].Ne) * PROTONMASS;
#endif
#ifdef VARIABLE_GAMMA
                      double temp_in_K = (SphP[li].GammaE - 1.0) * SphP[li].Utherm / BOLTZMANN * All.UnitEnergy_in_cgs /
                                         All.UnitMass_in_g * meanweight;
#else
                      double temp_in_K =
                          GAMMA_MINUS1 * SphP[li].Utherm / BOLTZMANN * All.UnitEnergy_in_cgs / All.UnitMass_in_g * meanweight;
#endif
#ifdef SGCHEM
                      double yn, en, yntot;

                      yn = rho * All.UnitDensity_in_cgs / ((1.0 + 4.0 * ABHE) * PROTONMASS);
                      en = SphP[li].Utherm * rho * All.UnitEnergy_in_cgs / pow(All.UnitLength_in_cm, 3);
#if CHEMISTRYNETWORK == 1
                      yntot = (1.0 + ABHE - SphP[li].TracAbund[IH2] + SphP[li].TracAbund[IHP] + SphP[li].TracAbund[IHEP] +
                               SphP[li].TracAbund[IHEPP]) *
                              yn;
#else
                      yntot     = (1.0 + ABHE - SphP[li].TracAbund[IH2] + SphP[li].TracAbund[IHP]) * yn;
#endif
#ifdef VARIABLE_GAMMA
                      temp_in_K = (SphP[li].GammaE - 1.0) * en / (yntot * BOLTZMANN);
#else
                      temp_in_K = (GAMMA - 1.0) * en / (yntot * BOLTZMANN);
#endif
#endif /* SGCHEM */

                      proj_rhotemp_sum += temp_in_K * rho * r;
#if defined(COOLING) && !defined(GRACKLE)
                      proj_szy_sum += temp_in_K * SphP[li].Ne * (rho * All.UnitDensity_in_cgs * HYDROGEN_MASSFRAC / PROTONMASS) *
                                      (r * All.UnitLength_in_cm);
#endif
#else
                      proj_rhotemp_sum += SphP[li].Utherm * rho * r;
#endif

#ifdef GFM_STELLAR_EVOLUTION
                      proj_rhometal_sum += SphP[li].Metallicity * rho * r;
#endif
#ifdef GFM_AGN_RADIATION
                      proj_rhoagnbol_sum += SphP[li].AGNBolIntensity * rho * r;
#endif
#ifdef TRACER_MC
                      proj_rhotracernum_sum += rho * r * get_number_of_tracers(li) * All.ReferenceTracerMCMass / P[li].Mass;
#endif
#ifdef CHEM_IMAGE
                      proj_dust_sum += SphP[li].DustTemp * rho * r;
                      proj_h2_sum += SphP[li].TracAbund[IH2] * rho * r;
                      proj_hp_sum += SphP[li].TracAbund[IHP] * rho * r;
#if CHEMISTRYNETWORK != 1
                      proj_co_sum += SphP[li].TracAbund[ICO] * rho * r;
#else
                      proj_co_sum += 0.0;
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
                      proj_rih_sum += sx_chem_image1(li) * rho * r;
                      proj_hrih_sum += sx_chem_image2(li) * rho * r;
#endif

#ifdef MRT
                      proj_photon_sum += SphP[li].DensPhot[0] * rho * r;
                      proj_h1_sum += SphP[li].HI * rho * r;
#endif
                    }
                }
            }
        }
    }

  *sigma = proj_sum;
  *sigmatemp = proj_rhotemp_sum;
  *sigmaweight = proj_weight_sum;
#if defined(COOLING) && !defined(GRACKLE)
  *sigmaszy = proj_szy_sum;
#endif
#ifdef GFM_STELLAR_EVOLUTION
  *sigmametal = proj_rhometal_sum;
#endif
#ifdef GFM_AGN_RADIATION
  *sigmaagnbol = proj_rhoagnbol_sum;
#endif
#ifdef TRACER_MC
  *sigmatracernum = proj_rhotracernum_sum;
  *sigmatrweight = proj_trweight_sum;
#endif
#ifdef CHEM_IMAGE
  *sigmadust = proj_dust_sum;
  *sigmah2 = proj_h2_sum;
  *sigmahp = proj_hp_sum;
  *sigmaco = proj_co_sum;
#endif
#ifdef MRT
  *sigmaphoton = proj_photon_sum;
  *sigmah1 = proj_h1_sum;
#endif

#ifdef SX_OUTPUT_IMAGE
  *sigmarih = proj_rih_sum;
  *sigmahrih = proj_hrih_sum;
#endif
}

#endif

/*! \brief Get position of pixel for slice imaging
 *
 *  \param[in] i Index of pixel in x direction
 *  \param[in] j Index of pixel in y direction
 *  \param[in] xaxis Slice x axis in simulation coordinate system
 *  \param[in] yaxis Slice y axis in simulation coordinate system
 *  \param[in] zaxis Slice z axis in simulaiton coordinate system
 *  \param[in] pixels_x Number of pixels in x direction
 *  \param[in] pixels_y Number of pixels in y direction
 *  \param[in] xmin Minimum x coordinate in slice coordinate system
 *  \param[in] xmax Maximum x coordinate in slice coordinate system
 *  \param[in] ymin Minimum y coordinate in slice coordinate system
 *  \param[in] ymax Maximum y coordinate in slice coordinate system
 *  \param[in] zval Z position of slice in slice coordinate system
 *  \param[out] p Array returning position of point
 *
 *  \return void
 */
void extract_position_from_axes(int i, int j, int xaxis, int yaxis, int zaxis, int pixels_x, int pixels_y, double xmin, double xmax,
                                double ymin, double ymax, double zval, double *p)
{
  /* The logic below sets the point p depending on which axes
   * we have chosen. */
  if(xaxis == 0 && yaxis == 1)
    {
      p[0] = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
      p[1] = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
      p[2] = zval;
    }
  else if(xaxis == 1 && yaxis == 0)
    {
      p[0] = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
      p[1] = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
      p[2] = zval;
    }
  else if(xaxis == 0 && yaxis == 2)
    {
      p[0] = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
      p[1] = zval;
      p[2] = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
    }
  else if(xaxis == 2 && yaxis == 0)
    {
      p[0] = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
      p[1] = zval;
      p[2] = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
    }
  else if(xaxis == 1 && yaxis == 2)
    {
      p[0] = zval;
      p[1] = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
      p[2] = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
    }
  else if(xaxis == 2 && yaxis == 1)
    {
      p[0] = zval;
      p[1] = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
      p[2] = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
    }
  else
    terminate("invalid combination of axes");
}

/*! \brief Generates an image of several fields in a slice given by the
 *  specified coordinate ranges.
 *
 *  \param[in] num Number of the image
 *  \param[in] gradients_flag Flag if gradients are used as well
 *  \param[in] pixels_x Number of pixels in x direction
 *  \param[in] pixels_y Number of pixels in y direction
 *  \param[in] xaxis Slice x axis in simulation coordinate system
 *  \param[in] yaxis Slice y axis in simulation coordinate system
 *  \param[in] zaxis Slice z axis in simulation coordinate system
 *  \param[in] xmin Minimum x coordinate in slice coordinate system
 *  \param[in] xmax Maximum x coordinate in slice coordinate system
 *  \param[in] ymin Minimum y coordinate in slice coordinate system
 *  \param[in] ymax Maximum y coordinate in slice coordinate system
 *  \param[in] zval Z value of slice in slice coordinate system
 *
 *  \return void
 */
void make_3d_voronoi_slice_image(int num, int gradients_flag, int pixels_x, int pixels_y, int xaxis, int yaxis, int zaxis, double xmin,
                                 double xmax, double ymin, double ymax, double zval)
{
  CPU_Step[CPU_MISC] += measure_time();

  float *density = 0, *temp = 0, *met = 0, *velocity = 0, *Bfield = 0, *vorticity = 0, *photon_density = 0, *chem_elements = 0,
        *density_trmc = 0;
  FILE *fd = 0, *fdtemp = 0, *fdmet = 0, *fdvel = 0, *fdmag = 0, *fdvort = 0, *fdphot = 0, *fdchem = 0, *fdtr = 0;
  char buf[MAXLEN_PATH];
#ifdef CHEM_IMAGE
  float *dust, *xH2, *xHP, *xCO;
  FILE *fddust = 0, *fdh2 = 0, *fdhp = 0, *fdco = 0;
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  float *xCHX, *xOHX, *xHCOP, *xCP, *xMP;
  FILE *fdchx = 0, *fdohx = 0, *fdhcop = 0, *fdcp = 0, *fdmp = 0;
#endif
#if CHEMISTRYNETWORK == 1 || CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  float *xHEP;
  FILE *fdhep = 0;
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
  float *rih = 0, *hrih = 0;
  FILE *fdrih = 0, *fdhrih = 0;
#endif
#ifdef SX_OUTPUT_IMAGE_ALL
  float *sxrates  = 0;
  FILE *fdsxrates = 0;
#endif

  if(gradients_flag == 1)
    file_path_sprintf(buf, "slice");
  else if(gradients_flag == 0)
    file_path_sprintf(buf, "slice_nograds");
  else
    terminate("gradients_flag != 1 && gradients_flag != 0");

  mpi_printf("we start to generate the slice image... gradients_flag=%d\n", gradients_flag);

  open_image_files(buf, num, &fd, &fdtemp, &fdmet, &fdvel,
#ifdef MHD
                   &fdmag,
#else
                   0,
#endif
#ifdef OUTPUT_VORTICITY
                   &fdvort,
#else
                   0,
#endif
#ifdef RT_ADVECT
                   &fdphot,
#else
                   0,
#endif
#ifdef GFM_STELLAR_EVOLUTION
                   &fdchem,
#else
                   0,
#endif
#ifdef TRACER_MC
                   &fdtr,
#else
                   0,
#endif

#ifdef CHEM_IMAGE
                   &fddust, &fdh2, &fdhp, &fdco
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
                   ,
                   &fdchx, &fdohx, &fdhcop, &fdcp, &fdmp, &fdhep,
#elif CHEMISTRYNETWORK == 1
                   ,
                   0, 0, 0, 0, 0, &fdhep,
#else
                   ,
                   0, 0, 0, 0, 0, 0,
#endif
#else
                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#endif

#ifdef SX_OUTPUT_IMAGE
                   &fdrih, &fdhrih,
#else
                   0, 0,
#endif
#ifdef SX_OUTPUT_IMAGE_ALL
                   &fdsxrates
#else
                   0
#endif
  );

  write_image_header(fd, pixels_x, pixels_y, 0);
  write_image_header(fdtemp, pixels_x, pixels_y, 0);
  write_image_header(fdmet, pixels_x, pixels_y, 0);
  write_image_header(fdvel, pixels_x, pixels_y, 0);
  write_image_header(fdmag, pixels_x, pixels_y, 0);
  write_image_header(fdvort, pixels_x, pixels_y, 0);
  write_image_header(fdphot, pixels_x, pixels_y, 0);
  write_image_header(fdchem, pixels_x, pixels_y, 0);
  write_image_header(fdtr, pixels_x, pixels_y, 0);
#ifdef CHEM_IMAGE
  write_image_header(fddust, pixels_x, pixels_y, 0);
  write_image_header(fdh2, pixels_x, pixels_y, 0);
  write_image_header(fdhp, pixels_x, pixels_y, 0);
  write_image_header(fdco, pixels_x, pixels_y, 0);
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  write_image_header(fdchx, pixels_x, pixels_y, 0);
  write_image_header(fdohx, pixels_x, pixels_y, 0);
  write_image_header(fdhcop, pixels_x, pixels_y, 0);
  write_image_header(fdcp, pixels_x, pixels_y, 0);
  write_image_header(fdmp, pixels_x, pixels_y, 0);
#endif
#if CHEMISTRYNETWORK == 1 || CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  write_image_header(fdhep, pixels_x, pixels_y, 0);
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
  write_image_header(fdrih, pixels_x, pixels_y, 0);
  write_image_header(fdhrih, pixels_x, pixels_y, 0);
#endif
#ifdef SX_OUTPUT_IMAGE_ALL
  write_image_header(fdsxrates, pixels_x, pixels_y, 0);
#endif

  density        = (float *)mymalloc("density", pixels_x * pixels_y * sizeof(float));
  temp           = (float *)mymalloc("temp", pixels_x * pixels_y * sizeof(float));
  met            = (float *)mymalloc("met", pixels_x * pixels_y * sizeof(float));
  velocity       = (float *)mymalloc("velocity", 3 * pixels_x * pixels_y * sizeof(float));
  Bfield         = (float *)mymalloc("Bfield", 3 * pixels_x * pixels_y * sizeof(float));
  vorticity      = (float *)mymalloc("vorticity", 3 * pixels_x * pixels_y * sizeof(float));
  photon_density = (float *)mymalloc("photon_density", RT_N_DIR * pixels_x * pixels_y * sizeof(float));
  chem_elements  = (float *)mymalloc("chem_elements", GFM_N_CHEM_ELEMENTS * pixels_x * pixels_y * sizeof(float));
  density_trmc   = (float *)mymalloc("density_trmc", pixels_x * pixels_y * sizeof(float));
#ifdef CHEM_IMAGE
  dust = (float *)mymalloc("dust", pixels_x * pixels_y * sizeof(float));
  xH2  = (float *)mymalloc("xH2", pixels_x * pixels_y * sizeof(float));
  xHP  = (float *)mymalloc("xHP", pixels_x * pixels_y * sizeof(float));
  xCO  = (float *)mymalloc("xCO", pixels_x * pixels_y * sizeof(float));
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  xCHX  = (float *)mymalloc("CHX", pixels_x * pixels_y * sizeof(float));
  xOHX  = (float *)mymalloc("OHX", pixels_x * pixels_y * sizeof(float));
  xHCOP = (float *)mymalloc("HCOP", pixels_x * pixels_y * sizeof(float));
  xCP   = (float *)mymalloc("CP", pixels_x * pixels_y * sizeof(float));
  xMP   = (float *)mymalloc("MP", pixels_x * pixels_y * sizeof(float));
#endif
#if CHEMISTRYNETWORK == 1 || CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  xHEP = (float *)mymalloc("HEP", pixels_x * pixels_y * sizeof(float));
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
  rih  = (float *)mymalloc("rih", pixels_x * pixels_y * sizeof(float));
  hrih = (float *)mymalloc("hrih", pixels_x * pixels_y * sizeof(float));
#endif
#ifdef SX_OUTPUT_IMAGE_ALL
  sxrates = (float *)mymalloc("sxrates", SX_NIMAGES * pixels_x * pixels_y * sizeof(float));
#endif

  MaxNray = pixels_x * pixels_y;
  Ray     = (ray_data *)mymalloc_movable(&Ray, "Ray", MaxNray * sizeof(ray_data));

  /* For simplicity, we use the "Ray" struct from the projection code. */

  setup_rays(pixels_x, pixels_y, xaxis, yaxis, zaxis, xmin, xmax, ymin, ymax, zval, zval);

  fill_slice(Ray, Nray, gradients_flag, pixels_x, pixels_y, density, temp, met, velocity, Bfield, vorticity, photon_density,
             chem_elements, density_trmc
#ifdef CHEM_IMAGE
             ,
             dust, xH2, xHP, xCO
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
             ,
             xCHX, xOHX, xHCOP, xCP, xMP
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16 || CHEMISTRYNETWORK == 1
             ,
             xHEP
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
             ,
             rih, hrih
#endif
#ifdef SX_OUTPUT_IMAGE_ALL
             ,
             sxrates
#endif
  );

  my_fwrite(density, sizeof(float), pixels_x * pixels_y, fd);
  my_fwrite(temp, sizeof(float), pixels_x * pixels_y, fdtemp);
  my_fwrite(met, sizeof(float), pixels_x * pixels_y, fdmet);
  my_fwrite(velocity, sizeof(float), 3 * pixels_x * pixels_y, fdvel);
  my_fwrite(Bfield, sizeof(float), 3 * pixels_x * pixels_y, fdmag);
  my_fwrite(vorticity, sizeof(float), 3 * pixels_x * pixels_y, fdvort);
  my_fwrite(photon_density, sizeof(float), RT_N_DIR * pixels_x * pixels_y, fdphot);
  my_fwrite(chem_elements, sizeof(float), GFM_N_CHEM_ELEMENTS * pixels_x * pixels_y, fdchem);
  my_fwrite(density_trmc, sizeof(float), pixels_x * pixels_y, fdtr);

#ifdef SX_OUTPUT_IMAGE
  my_fwrite(rih, sizeof(float), pixels_x * pixels_y, fdrih);
  my_fwrite(hrih, sizeof(float), pixels_x * pixels_y, fdhrih);
#endif
#ifdef SX_OUTPUT_IMAGE_ALL
  my_fwrite(sxrates, sizeof(float), SX_NIMAGES * pixels_x * pixels_y, fdsxrates);
#endif

#ifdef CHEM_IMAGE
  my_fwrite(dust, sizeof(float), pixels_x * pixels_y, fddust);
  my_fwrite(xH2, sizeof(float), pixels_x * pixels_y, fdh2);
  my_fwrite(xHP, sizeof(float), pixels_x * pixels_y, fdhp);
  my_fwrite(xCO, sizeof(float), pixels_x * pixels_y, fdco);
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  my_fwrite(xCHX, sizeof(float), pixels_x * pixels_y, fdchx);
  my_fwrite(xOHX, sizeof(float), pixels_x * pixels_y, fdohx);
  my_fwrite(xHCOP, sizeof(float), pixels_x * pixels_y, fdhcop);
  my_fwrite(xCP, sizeof(float), pixels_x * pixels_y, fdcp);
  my_fwrite(xMP, sizeof(float), pixels_x * pixels_y, fdmp);
#endif
#if CHEMISTRYNETWORK == 1 || CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  my_fwrite(xHEP, sizeof(float), pixels_x * pixels_y, fdhep);
#endif
#endif

  if(fd)
    fclose(fd);
  if(fdtemp)
    fclose(fdtemp);
  if(fdmet)
    fclose(fdmet);
  if(fdvel)
    fclose(fdvel);
  if(fdmag)
    fclose(fdmag);
  if(fdvort)
    fclose(fdvort);
  if(fdphot)
    fclose(fdphot);
  if(fdchem)
    fclose(fdchem);
  if(fdtr)
    fclose(fdtr);
#ifdef CHEM_IMAGE
  if(fddust)
    fclose(fddust);
  if(fdh2)
    fclose(fdh2);
  if(fdhp)
    fclose(fdhp);
  if(fdco)
    fclose(fdco);
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  if(fdchx)
    fclose(fdchx);
  if(fdohx)
    fclose(fdohx);
  if(fdhcop)
    fclose(fdhcop);
  if(fdcp)
    fclose(fdcp);
  if(fdmp)
    fclose(fdmp);
#endif
#if CHEMISTRYNETWORK == 1 || CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  if(fdhep)
    fclose(fdhep);
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
  if(fdrih)
    fclose(fdrih);
  if(fdhrih)
    fclose(fdhrih);
#endif
#ifdef SX_OUTPUT_IMAGE_ALL
  if(fdsxrates)
    fclose(fdsxrates);
#endif

  myfree(Ray);

#ifdef SX_OUTPUT_IMAGE_ALL
  myfree(sxrates);
#endif
#ifdef SX_OUTPUT_IMAGE
  myfree(hrih);
  myfree(rih);
#endif

#ifdef CHEM_IMAGE
#if CHEMISTRYNETWORK == 1 || CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  myfree(xHEP);
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  myfree(xMP);
  myfree(xCP);
  myfree(xHCOP);
  myfree(xOHX);
  myfree(xCHX);
#endif
  myfree(xCO);
  myfree(xHP);
  myfree(xH2);
  myfree(dust);
#endif
  myfree(density_trmc);
  myfree(chem_elements);
  myfree(photon_density);
  myfree(vorticity);
  myfree(Bfield);
  myfree(velocity);
  myfree(met);
  myfree(temp);
  myfree(density);

  CPU_Step[CPU_MAKEIMAGES] += measure_time();
}

FILE *FdFaces;

/*! \brief Writes Voronoi faces intersecting with a plane to file
 *
 *  Writes a file listing all the Voronoi faces that intersect the
 *  slice plane as a list of line segments (x0,y0,x1,y1). Note that all
 *  intersections across the entire box are written, regardless of the
 *  x/y limits made on the image slice. Also note that this only works
 *  correctly on one processor because no attempt is made to gather the
 *  tessellations from a decomposed domain.
 *
 *  \param[in] T Pointer to tessellation
 *  \param[in] num Number of the image in filename
 *  \param[in] xaxis Slice x axis
 *  \param[in] yaxis Slice y axis
 *  \param[in] zaxis Slice z axis
 *  \param[in] zval Position along Slice z axis
 *
 *  \return void
 */
void make_3d_voronoi_listfaces(tessellation *T, int num, int xaxis, int yaxis, int zaxis, double zval)
{
  CPU_Step[CPU_MISC] += measure_time();

  int i, nr, bit;
  char buf[MAXLEN_PATH];

  tetra_center *tmpDTC = T->DTC;
  T->DTC               = (tetra_center *)mymalloc_movable(&T->DTC, "DTC", T->MaxNdt * sizeof(tetra_center));
  T->DTF               = (char *)mymalloc_movable(&T->DTF, "DTF", T->MaxNdt * sizeof(char));
  for(i = 0; i < T->Ndt; i++)
    T->DTF[i] = 0;
  compute_circumcircles(T);

  Edge_visited = (unsigned char *)mymalloc("Edge_visited", T->Ndt * sizeof(unsigned char));

  for(i = 0; i < T->Ndt; i++)
    Edge_visited[i] = 0;

  file_path_sprintf(buf, "%s/faces_list_%03d.txt", All.OutputDir, num);

  if(ThisTask == 0)
    if(!(FdFaces = fopen(buf, "w")))
      terminate("can't open file `%s' for writing snapshot.\n", buf);

  for(i = 0; i < T->Ndt; i++)
    {
      if(T->DT[i].t[0] < 0) /* deleted? */
        continue;

      bit = 1;
      nr  = 0;

      /* Loop over all edges in tetra i. */
      while(Edge_visited[i] != EDGE_ALL)
        {
          if((Edge_visited[i] & bit) == 0)
            make_3d_voronoi_listfaces_check_for_cut(T, i, nr, xaxis, yaxis, zaxis, zval);

          bit <<= 1;
          nr++;
        }
    }

  if(ThisTask == 0)
    fclose(FdFaces);

  myfree(Edge_visited);

  myfree(T->DTF);
  myfree(T->DTC);
  T->DTC = tmpDTC;

  CPU_Step[CPU_MAKEIMAGES] += measure_time();
}

/*! \brief  Checks whether the triangle defined by the three points in pp
 *  intersects the slice plane, and if so writes it to the output
 *  file.
 *
 *  \param[in] pp
 *  \param[in] xaxis Projection x axis
 *  \param[in] yaxis Projection y axis
 *  \param[in] zaxis Projection z axis
 *  \param[in] Position along projeciton z axis
 *
 *  \return void
 */
void check_for_cut(double pp[3][3], int xaxis, int yaxis, int zaxis, double zval)
{
  int i, ii, count;
  double x[3], y[3], w;

  for(i = 0, count = 0; i < 3; i++)
    {
      ii = i + 1;
      if(ii > 2)
        ii = 0;

      // The edge only intersects the slice plane if the points are on
      // opposite sides.
      if((pp[i][zaxis] - zval) * (pp[ii][zaxis] - zval) < 0)
        {
          // determine the x and y of the intersection point
          w        = (zval - pp[i][zaxis]) / (pp[ii][zaxis] - pp[i][zaxis]);
          x[count] = pp[i][xaxis] + w * (pp[ii][xaxis] - pp[i][xaxis]);
          y[count] = pp[i][yaxis] + w * (pp[ii][yaxis] - pp[i][yaxis]);
          count++;
        }
    }

  if(count == 2 && ThisTask == 0)
    fprintf(FdFaces, "%g %g %g %g\n", x[0], y[0], x[1], y[1]);
}

/*! \brief Checks whether the face corresponding to edge nr in tetrahedron tt
 *  intersects the slice plane.
 *
 *  \param[in] T Pointer to tessellation
 *  \param[in] tt Index in DT array
 *  \param[in] nr Index of edge
 *  \param[in] xaxis Projection x axis
 *  \param[in] yaxis Projection y axis
 *  \param[in] zaxis Projection z axis
 *  \param[in] zval Position along projection z axis
 *
 *  \return void
 */
void make_3d_voronoi_listfaces_check_for_cut(tessellation *T, int tt, int nr, int xaxis, int yaxis, int zaxis, double zval)
{
  point *DP         = T->DP;
  tetra *DT         = T->DT;
  tetra_center *DTC = T->DTC;
  int i, j, k, l, m, ii, jj, kk, ll, count, nr_next, flag, nn;
  tetra *prev, *next;
  tetra_center *prevc, *nextc;
  double pp[3][3];

  tetra *t         = &DT[tt];
  tetra_center *tc = &DTC[tt];

  i = edge_start[nr];
  j = edge_end[nr];
  k = edge_opposite[nr];
  l = edge_nexttetra[nr];

  Edge_visited[tt] |= (1 << nr);

  if(T->Nvf + 1 >= T->MaxNvf)
    terminate("Nvf + 1 >= MaxNvf");

  pp[0][0] = tc->cx;
  pp[0][1] = tc->cy;
  pp[0][2] = tc->cz;

  count = 0;

  flag = 0;

  if(DP[t->p[i]].task == ThisTask && DP[t->p[i]].index >= 0 && DP[t->p[i]].index < NumGas)
    flag = 1;

  if(DP[t->p[j]].task == ThisTask && DP[t->p[j]].index >= 0 && DP[t->p[j]].index < NumGas)
    flag = 1;

  prev  = t;
  prevc = tc;
  do
    {
      nn    = prev->t[l];
      next  = &DT[nn];
      nextc = &DTC[nn];

      if(prev != t && next != t)
        {
          pp[1][0] = prevc->cx;
          pp[1][1] = prevc->cy;
          pp[1][2] = prevc->cz;

          pp[2][0] = nextc->cx;
          pp[2][1] = nextc->cy;
          pp[2][2] = nextc->cz;

          if(flag)
            check_for_cut(pp, xaxis, yaxis, zaxis, zval);
        }

      for(m = 0, ll = ii = jj = -1; m < 4; m++)
        {
          if(next->p[m] == prev->p[k])
            ll = m;
          if(next->p[m] == prev->p[i])
            ii = m;
          if(next->p[m] == prev->p[j])
            jj = m;
        }

      if(ll < 0 || ii < 0 || jj < 0)
        terminate("ll < 0 || ii < 0 || jj < 0");

      kk = 6 - (ll + ii + jj);

      /* need to determine the edge number to be able to flag it */

      for(nr_next = 0; nr_next < 6; nr_next++)
        if((edge_start[nr_next] == ii && edge_end[nr_next] == jj) || (edge_start[nr_next] == jj && edge_end[nr_next] == ii))
          {
            if((Edge_visited[nn] & (1 << nr_next)) && next != t)
              terminate("can't be");

            Edge_visited[nn] |= (1 << nr_next);
            break;
          }

      prevc = nextc;
      prev  = next;
      i     = ii;
      l     = ll;
      j     = jj;
      k     = kk;

      count++;

      if(count > 1000)
        terminate("count > 1000");
    }
  while(next != t);
}

#endif

#ifdef TWODIMS

/*! \brief Crate a 2d image of entire box
 *
 *  \param[in] num Number of the image in filename
 *  \param[in] pixels_x Number of pixels in x direction
 *  \param[in] pixels_y Number of pixels in y direction
 *
 *  \return void
 */
void make_2d_voronoi_image(int num, int pixels_x, int pixels_y)
{
  CPU_Step[CPU_MISC] += measure_time();

  char buf[MAXLEN_PATH];
  float *dens, *denssum, *dp;
  FILE *fd = 0;
  point *p;
  int pp;
  int tt0, ttstart, ttrow;
  tetra *t0;
  double l_dx, l_dy;
  int i, j, k, kmin, li, moves, ret, no, task;
  double r2, r2min, rho_L;
  peanokey key;

  if(Mesh.Ndp >= Mesh.MaxNdp)
    {
      terminate("Ndp >= MaxNdp");
    }

  pp = Mesh.Ndp;
  p  = &Mesh.DP[pp];

  file_path_sprintf(buf, "%s/density_field_%03d", All.OutputDir, num);

  if(ThisTask == 0)
    {
      if(!(fd = fopen(buf, "w")))
        terminate("can't open file `%s' for writing snapshot.\n", buf);

      my_fwrite(&pixels_x, sizeof(int), 1, fd);
      my_fwrite(&pixels_y, sizeof(int), 1, fd);
    }

  dens    = (float *)mymalloc("dens", pixels_x * pixels_y * sizeof(float));
  denssum = (float *)mymalloc("denssum", pixels_x * pixels_y * sizeof(float));

  for(i = 0, dp = dens; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *dp++ = 0;

  ttrow = 0;

  point *DP = Mesh.DP;

  for(i = 0; i < pixels_x; i++)
    {
      ttstart = ttrow;

      for(j = 0; j < pixels_y; j++)
        {
          p->x = (i + 0.5) / pixels_x * boxSize_X;
          p->y = (j + 0.5) / pixels_y * boxSize_Y;
          p->z = 0;

          key = peano_hilbert_key((int)((p->x - DomainCorner[0]) * DomainFac), (int)((p->y - DomainCorner[1]) * DomainFac),
                                  (int)((p->z - DomainCorner[2]) * DomainFac), BITS_PER_DIMENSION);

          no = 0;
          while(TopNodes[no].Daughter >= 0)
            no = TopNodes[no].Daughter + (key - TopNodes[no].StartKey) / (TopNodes[no].Size >> 3);

          no   = TopNodes[no].Leaf;
          task = DomainTask[no];

          if(task == ThisTask)
            {
#ifndef OPTIMIZE_MEMORY_USAGE
              set_integers_for_point(&Mesh, pp);
#endif
              tt0 = get_triangle(&Mesh, pp, &moves, &ret, ttstart);
              t0  = &Mesh.DT[tt0];

              for(k = 0, kmin = 0, r2min = 1.0e30; k < 3; k++)
                {
                  r2 = (p->x - DP[t0->p[k]].x) * (p->x - DP[t0->p[k]].x) + (p->y - DP[t0->p[k]].y) * (p->y - DP[t0->p[k]].y);
                  if(r2 < r2min)
                    {
                      r2min = r2;
                      kmin  = k;
                    }
                }

              li = DP[t0->p[kmin]].index;

              if(li >= NumGas)
                li -= NumGas;

              if(DP[t0->p[kmin]].task == ThisTask)
                {
                  l_dx = p->x - SphP[li].Center[0];
                  l_dy = p->y - SphP[li].Center[1];
                }
              else
                {
                  l_dx = p->x - PrimExch[li].Center[0];
                  l_dy = p->y - PrimExch[li].Center[1];
                }
              l_dx = nearest_x(l_dx);
              l_dy = nearest_y(l_dy);
              if(DP[t0->p[kmin]].task == ThisTask)
                {
                  rho_L = SphP[li].Density + SphP[li].Grad.drho[0] * l_dx + SphP[li].Grad.drho[1] * l_dy;
                }
              else
                {
                  rho_L = PrimExch[li].Density + GradExch[li].drho[0] * l_dx + GradExch[li].drho[1] * l_dy;
                }

              dens[i * pixels_y + j] = rho_L;

              ttstart = tt0;

              if(j == 0)
                ttrow = tt0;
            }
        }
    }

  MPI_Reduce(dens, denssum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      my_fwrite(denssum, sizeof(float), pixels_x * pixels_y, fd);
      fclose(fd);
    }

  myfree(denssum);
  myfree(dens);

  CPU_Step[CPU_MAKEIMAGES] += measure_time();
}

/*! \brief Crate a 2d zoomed image
 *
 *  \param[in] T Pointer to tessellation
 *  \param[in] num Number of the image in filename
 *  \param[in] pixels_x Number of pixels in x direction
 *  \param[in] pixels_y Number of pixels in y direction
 *  \param[in] xmin Minimum x coordinate in image box
 *  \param[in] xmax Maximum x coordinate in image box
 *  \param[in] ymin Minimum y coordinate in image box
 *  \param[in] ymax Maximum y coordinate in image box
 *
 *  \return void
 */
void make_2d_voronoi_image_zoomed(tessellation *T, int num, int pixels_x, int pixels_y, double xmin, double xmax, double ymin,
                                  double ymax)
{
  CPU_Step[CPU_MISC] += measure_time();

  float *density = 0, *temp = 0, *velocity = 0, *vorticity = 0;
  FILE *fd = 0, *fdtemp = 0, *fdvel = 0, *fdvort = 0;

  open_image_files("field", num, &fd, &fdtemp, 0, &fdvel, 0,
#ifdef OUTPUT_VORTICITY
                   &fdvort,
#else
                   0,
#endif
                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

  write_image_header(fd, pixels_x, pixels_y, 0);
  write_image_header(fdtemp, pixels_x, pixels_y, 0);
  write_image_header(fdvel, pixels_x, pixels_y, 0);
  write_image_header(fdvort, pixels_x, pixels_y, 0);

  density   = (float *)mymalloc("density", pixels_x * pixels_y * sizeof(float));
  temp      = (float *)mymalloc("temp", pixels_x * pixels_y * sizeof(float));
  velocity  = (float *)mymalloc("velocity", 3 * pixels_x * pixels_y * sizeof(float));
  vorticity = (float *)mymalloc("vorticity", 3 * pixels_x * pixels_y * sizeof(float));

  MaxNray = pixels_x * pixels_y;
  Ray     = (ray_data *)mymalloc_movable(&Ray, "Ray", MaxNray * sizeof(ray_data));

  /* For simplicity, we use the "Ray" struct from the projection
     code. */

  setup_rays(pixels_x, pixels_y, 0, 1, 2, xmin, xmax, ymin, ymax, 0.0, 0.0);

  fill_slice(Ray, Nray, 1, pixels_x, pixels_y, density, temp, 0, velocity, 0, vorticity, 0, 0, 0);

  my_fwrite(density, sizeof(float), pixels_x * pixels_y, fd);
  my_fwrite(temp, sizeof(float), pixels_x * pixels_y, fdtemp);
  my_fwrite(velocity, sizeof(float), 3 * pixels_x * pixels_y, fdvel);
  my_fwrite(vorticity, sizeof(float), 3 * pixels_x * pixels_y, fdvort);

  if(fd)
    fclose(fd);
  if(fdtemp)
    fclose(fdtemp);
  if(fdvel)
    fclose(fdvel);
  if(fdvort)
    fclose(fdvort);

  myfree(Ray);
  myfree(vorticity);
  myfree(velocity);
  myfree(temp);
  myfree(density);

  CPU_Step[CPU_MAKEIMAGES] += measure_time();
}

#endif

#if !defined(ONEDIMS) && !defined(TWODIMS)

/*! \brief Creates a 3d grid of various properties
 *
 *  \param[in] num Number of the image in filename
 *  \param[in] pixels_x Number of pixels in x direction
 *  \param[in] pixels_y Number of pixels in y direction
 *  \param[in] pixels_z Number of pixels in z direction
 *  \param[in] xmin Minimum x coordinate in grid box
 *  \param[in] xmax Maximum x coordinate in grid box
 *  \param[in] ymin Minimum y coordinate in grid box
 *  \param[in] ymax Maximum y coordinate in grid box
 *  \param[in] zmin Minimum z coordinate in grid box
 *  \param[in] zmax Maximum z coordinate in grid box
 *
 *  \return void
 */
void make_3d_voronoi_grid(int num, int pixels_x, int pixels_y, int pixels_z, double xmin, double xmax, double ymin, double ymax,
                          double zmin, double zmax)
{
  CPU_Step[CPU_MISC] += measure_time();

  float *density = 0, *temp = 0, *met = 0, *velocity = 0, *Bfield = 0, *vorticity = 0, *photon_density = 0, *chem_elements = 0,
        *density_trmc = 0;
  FILE *fd = 0, *fdvel = 0, *fdtemp = 0, *fdmet = 0, *fdmag = 0, *fdvort = 0, *fdphot = 0, *fdchem = 0, *fdtr = 0;
#ifdef CHEM_IMAGE
  float *dust = 0, *xH2 = 0, *xHP = 0, *xCO = 0;
  FILE *fddust = 0, *fdh2 = 0, *fdhp = 0, *fdco = 0;
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  float *xCHX, *xOHX, *xHCOP, *xCP, *xMP;
  FILE *fdchx = 0, *fdohx = 0, *fdhcop = 0, *fdcp = 0, *fdmp = 0;
#endif
#if CHEMISTRYNETWORK == 1 || CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  float *xHEP;
  FILE *fdhep = 0;
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
  float *rih = 0, *hrih = 0;
  FILE *fdrih = 0, *fdhrih = 0;
#endif
#ifdef SX_OUTPUT_IMAGE_ALL
  float *sxrates  = 0;
  FILE *fdsxrates = 0;
#endif

  double xstep;

#ifdef VORONOI_NOGRADS
  int gradients_flag = 0;
#else
  int gradients_flag = 1;
#endif

  open_image_files("grid", num, &fd, &fdtemp, &fdmet, &fdvel,
#ifdef MHD
                   &fdmag,
#else
                   0,
#endif
#ifdef OUTPUT_VORTICITY
                   &fdvort,
#else
                   0,
#endif
#ifdef RT_ADVECT
                   &fdphot,
#else
                   0,
#endif
#ifdef GFM_STELLAR_EVOLUTION
                   &fdchem,
#else
                   0,
#endif
#ifdef TRACER_MC
                   &fdtr,
#else
                   0,
#endif

#ifdef CHEM_IMAGE
                   &fddust, &fdh2, &fdhp, &fdco,
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
                   &fdchx, &fdohx, &fdhcop, &fdcp, &fdmp, &fdhep,
#elif CHEMISTRYNETWORK == 1
                   0, 0, 0, 0, 0, &fdhep,
#else
                   0, 0, 0, 0, 0, 0,
#endif
#else
                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#endif

#ifdef SX_OUTPUT_IMAGE
                   &fdrih, &fdhrih,
#else
                   0, 0,
#endif
#ifdef SX_OUTPUT_IMAGE_ALL
                   &fdsxrates
#else
                   0
#endif
  );

  write_image_header(fd, pixels_x, pixels_y, pixels_z);
  write_image_header(fdtemp, pixels_x, pixels_y, pixels_z);
  write_image_header(fdmet, pixels_x, pixels_y, pixels_z);
  write_image_header(fdvel, pixels_x, pixels_y, pixels_z);
  write_image_header(fdmag, pixels_x, pixels_y, pixels_z);
  write_image_header(fdvort, pixels_x, pixels_y, pixels_z);
  write_image_header(fdphot, pixels_x, pixels_y, pixels_z);
  write_image_header(fdchem, pixels_x, pixels_y, pixels_z);
  write_image_header(fdtr, pixels_x, pixels_y, pixels_z);

#ifdef CHEM_IMAGE
  write_image_header(fddust, pixels_x, pixels_y, pixels_z);
  write_image_header(fdh2, pixels_x, pixels_y, pixels_z);
  write_image_header(fdhp, pixels_x, pixels_y, pixels_z);
  write_image_header(fdco, pixels_x, pixels_y, pixels_z);
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  write_image_header(fdchx, pixels_x, pixels_y, pixels_z);
  write_image_header(fdohx, pixels_x, pixels_y, pixels_z);
  write_image_header(fdhcop, pixels_x, pixels_y, pixels_z);
  write_image_header(fdcp, pixels_x, pixels_y, pixels_z);
  write_image_header(fdmp, pixels_x, pixels_y, pixels_z);
#endif
#if CHEMISTRYNETWORK == 1 || CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  write_image_header(fdhep, pixels_x, pixels_y, pixels_z);
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
  write_image_header(fdrih, pixels_x, pixels_y, pixels_z);
  write_image_header(fdhrih, pixels_x, pixels_y, pixels_z);
#endif
#ifdef SX_OUTPUT_IMAGE_ALL
  write_image_header(fdsxrates, pixels_x, pixels_y, pixels_z);
#endif

  MaxNray = pixels_y * pixels_z;
  Ray     = (ray_data *)mymalloc_movable(&Ray, "Ray", MaxNray * sizeof(ray_data));

  xstep = (xmax - xmin) / pixels_x;

  /* First create rays at the starting positions in the xmin plane,
     with length appropriate to take them to the next x-slice */
  setup_rays(pixels_y, pixels_z, 1, 2, 0, ymin, ymax, zmin, zmax, xmin, xmin + xstep);

  density        = (float *)mymalloc("density", pixels_y * pixels_z * sizeof(float));
  temp           = (float *)mymalloc("temp", pixels_y * pixels_z * sizeof(float));
  met            = (float *)mymalloc("met", pixels_y * pixels_z * sizeof(float));
  velocity       = (float *)mymalloc("velocity", 3 * pixels_y * pixels_z * sizeof(float));
  Bfield         = (float *)mymalloc("Bfield", 3 * pixels_y * pixels_z * sizeof(float));
  vorticity      = (float *)mymalloc("vorticity", 3 * pixels_y * pixels_z * sizeof(float));
  photon_density = (float *)mymalloc("photon_density", RT_N_DIR * pixels_y * pixels_z * sizeof(float));
  chem_elements  = (float *)mymalloc("chem_elements", GFM_N_CHEM_ELEMENTS * pixels_y * pixels_z * sizeof(float));
  density_trmc   = (float *)mymalloc("density_trmc", pixels_y * pixels_z * sizeof(float));

#ifdef CHEM_IMAGE
  dust = (float *)mymalloc("dust", pixels_y * pixels_z * sizeof(float));
  xH2  = (float *)mymalloc("H2", pixels_y * pixels_z * sizeof(float));
  xHP  = (float *)mymalloc("HP", pixels_y * pixels_z * sizeof(float));
  xCO  = (float *)mymalloc("CO", pixels_y * pixels_z * sizeof(float));
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  xCHX  = (float *)mymalloc("CHX", pixels_y * pixels_z * sizeof(float));
  xOHX  = (float *)mymalloc("OHX", pixels_y * pixels_z * sizeof(float));
  xHCOP = (float *)mymalloc("HCOP", pixels_y * pixels_z * sizeof(float));
  xCP   = (float *)mymalloc("CP", pixels_y * pixels_z * sizeof(float));
  xMP   = (float *)mymalloc("MP", pixels_y * pixels_z * sizeof(float));
#endif
#if CHEMISTRYNETWORK == 1 || CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  xHEP = (float *)mymalloc("HEP", pixels_y * pixels_z * sizeof(float));
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
  rih  = (float *)mymalloc("rih", pixels_y * pixels_z * sizeof(float));
  hrih = (float *)mymalloc("hrih", pixels_y * pixels_z * sizeof(float));
#endif
#ifdef SX_OUTPUT_IMAGE_ALL
  sxrates = (float *)mymalloc("sxrates", SX_NIMAGES * pixels_y * pixels_z * sizeof(float));
#endif

  mpi_printf("Extracting grid info...\n");

  /* loop over x-slices */
  for(int i = 0; i < pixels_x; i++)
    {
      /* reset ray propagation lengths */
      for(int j = 0; j < Nray; ++j)
        {
          Ray[j].len        = 0;
          Ray[j].target_len = xstep;
        }

      /* now we simply use fill_slice for the rays */
      fill_slice(Ray, Nray, gradients_flag, pixels_y, pixels_z, density, temp, met, velocity, Bfield, vorticity, photon_density,
                 chem_elements, density_trmc
#ifdef CHEM_IMAGE
                 ,
                 dust, xH2, xHP, xCO
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
                 ,
                 xCHX, xOHX, xHCOP, xCP, xMP
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16 || CHEMISTRYNETWORK == 1
                 ,
                 xHEP
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
                 ,
                 rih, hrih
#endif
#ifdef SX_OUTPUT_IMAGE_ALL
                 ,
                 sxrates
#endif
      );

#ifdef VORONOI_DYNAMIC_UPDATE
      /* step rays to the next x-position */
      int left_this_task, rays_left;
      do
        {
          left_this_task = advance_rays_for_one_cell(0, 0, 1);
          exchange_rays();

          MPI_Allreduce(&left_this_task, &rays_left, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        }
      while(rays_left);
#else
      /* without being able to trace the rays, we must redo the mesh
         search from scratch */
      setup_rays(pixels_y, pixels_z, 1, 2, 0, ymin, ymax, zmin, zmax, xmin + (i + 1) * xstep, xmax);
#endif

      /* write the slice */
      if(ThisTask == 0)
        {
          mpi_printf("Writing slice %d/%d      \r", i, pixels_x);
          int b, j, k;
          for(j = 0; j < pixels_y; j++)
            for(k = 0; k < pixels_z; k++)
              for(b = 0; b < 3; b++)
                {
                  myassert(gsl_finite(velocity[3 * (j * pixels_z + k) + b]));
                }

          /* these don't write if the fd is null */
          my_fwrite(density, sizeof(float), pixels_y * pixels_z, fd);
          my_fwrite(temp, sizeof(float), pixels_y * pixels_z, fdtemp);
          my_fwrite(met, sizeof(float), pixels_y * pixels_z, fdmet);
          my_fwrite(velocity, sizeof(float), 3 * pixels_y * pixels_z, fdvel);
          my_fwrite(Bfield, sizeof(float), 3 * pixels_y * pixels_z, fdmag);
          my_fwrite(vorticity, sizeof(float), 3 * pixels_y * pixels_z, fdvort);
          my_fwrite(photon_density, sizeof(float), RT_N_DIR * pixels_y * pixels_z, fdvort);
          my_fwrite(chem_elements, sizeof(float), GFM_N_CHEM_ELEMENTS * pixels_y * pixels_z, fdchem);
          my_fwrite(density_trmc, sizeof(float), pixels_y * pixels_z, fdtr);
#ifdef CHEM_IMAGE
          my_fwrite(dust, sizeof(float), pixels_y * pixels_z, fddust);
          my_fwrite(xH2, sizeof(float), pixels_y * pixels_z, fdh2);
          my_fwrite(xHP, sizeof(float), pixels_y * pixels_z, fdhp);
          my_fwrite(xCO, sizeof(float), pixels_y * pixels_z, fdco);
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
          my_fwrite(xCHX, sizeof(float), pixels_y * pixels_z, fdchx);
          my_fwrite(xOHX, sizeof(float), pixels_y * pixels_z, fdohx);
          my_fwrite(xHCOP, sizeof(float), pixels_y * pixels_z, fdhcop);
          my_fwrite(xCP, sizeof(float), pixels_y * pixels_z, fdcp);
          my_fwrite(xMP, sizeof(float), pixels_y * pixels_z, fdmp);
#endif
#if CHEMISTRYNETWORK == 1 || CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
          my_fwrite(xHEP, sizeof(float), pixels_y * pixels_z, fdhep);
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
          my_fwrite(rih, sizeof(float), pixels_y * pixels_z, fdrih);
          my_fwrite(hrih, sizeof(float), pixels_y * pixels_z, fdhrih);
#endif
#ifdef SX_OUTPUT_IMAGE_ALL
          my_fwrite(sxrates, sizeof(float), SX_NIMAGES * pixels_y * pixels_z, fdsxrates);
#endif
        }
    }

#ifdef IMAGE_FOOTERS
  write_image_footer(fd, xmin, xmax, ymin, ymax, zmin, zmax);
  write_image_footer(fdtemp, xmin, xmax, ymin, ymax, zmin, zmax);
  write_image_footer(fdmet, xmin, xmax, ymin, ymax, zmin, zmax);
  write_image_footer(fdvel, xmin, xmax, ymin, ymax, zmin, zmax);
#ifdef CHEM_IMAGE
  write_image_footer(fddust, xmin, xmax, ymin, ymax, zmin, zmax);
  write_image_footer(fdh2, xmin, xmax, ymin, ymax, zmin, zmax);
  write_image_footer(fdhp, xmin, xmax, ymin, ymax, zmin, zmax);
  write_image_footer(fdco, xmin, xmax, ymin, ymax, zmin, zmax);
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  write_image_footer(fdchx, xmin, xmax, ymin, ymax, zmin, zmax);
  write_image_footer(fdohx, xmin, xmax, ymin, ymax, zmin, zmax);
  write_image_footer(fdhcop, xmin, xmax, ymin, ymax, zmin, zmax);
  write_image_footer(fdcp, xmin, xmax, ymin, ymax, zmin, zmax);
  write_image_footer(fdmp, xmin, xmax, ymin, ymax, zmin, zmax);
#endif
#if CHEMISTRYNETWORK == 1 || CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  write_image_footer(fdhep, xmin, xmax, ymin, ymax, zmin, zmax);
#endif
#endif
#endif

  if(fd)
    fclose(fd);
  if(fdtemp)
    fclose(fdtemp);
  if(fdmet)
    fclose(fdmet);
  if(fdvel)
    fclose(fdvel);
  if(fdmag)
    fclose(fdmag);
  if(fdvort)
    fclose(fdvort);
  if(fdphot)
    fclose(fdphot);
  if(fdchem)
    fclose(fdchem);
  if(fdtr)
    fclose(fdtr);
#ifdef CHEM_IMAGE
  if(fddust)
    fclose(fddust);
  if(fdh2)
    fclose(fdh2);
  if(fdhp)
    fclose(fdhp);
  if(fdco)
    fclose(fdco);
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  if(fdchx)
    fclose(fdchx);
  if(fdohx)
    fclose(fdohx);
  if(fdhcop)
    fclose(fdhcop);
  if(fdcp)
    fclose(fdcp);
  if(fdmp)
    fclose(fdmp);
#endif
#if CHEMISTRYNETWORK == 1 || CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  if(fdhep)
    fclose(fdhep);
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
  if(fdrih)
    fclose(fdrih);
  if(fdhrih)
    fclose(fdhrih);
#endif
#ifdef SX_OUTPUT_IMAGE_ALL
  if(fdsxrates)
    fclose(fdsxrates);
#endif

#ifdef SX_OUTPUT_IMAGE_ALL
  myfree(sxrates);
#endif
#ifdef SX_OUTPUT_IMAGE
  myfree(hrih);
  myfree(rih);
#endif

#ifdef CHEM_IMAGE
#if CHEMISTRYNETWORK == 1 || CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  myfree(xHEP);
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  myfree(xMP);
  myfree(xCP);
  myfree(xHCOP);
  myfree(xOHX);
  myfree(xCHX);
#endif
  myfree(xCO);
  myfree(xHP);
  myfree(xH2);
  myfree(dust);
#endif
  myfree(density_trmc);
  myfree(chem_elements);
  myfree(photon_density);
  myfree(vorticity);
  myfree(Bfield);
  myfree(velocity);
  myfree(met);
  myfree(temp);
  myfree(density);
  myfree(Ray);

  CPU_Step[CPU_MAKEIMAGES] += measure_time();
}
#endif /* #if !defined(ONEDIMS) && !defined(TWODIMS) */

/*! \brief Fills an image slice with information for the rays
 *
 *   Is used by both the slice image and "grid of slices" code. The data is
 *   returned in the output arrays, which should be allocated with
 *   pixels_x*pixels_y entries (or 3x that in the case of the vector
 *   quantities).
 *
 *   \param[in] ray Array with rays of interest
 *   \param[in] Nray Number of rays; length of ray array
 *   \param[in] gradients_flag Flag whether to take into account gradients
 *   \param[in] pixels_x Number of pixels in x direction
 *   \param[in] pixels_y Number of pixels in y direction
 *   \param[out] density Array for density data
 *   \param[out] temperature Array for temperature data
 *   \param[out] metallicity Array for metallicity data
 *   \param[out] velocity Array for velocity data
 *   \param[out] Bfield Array for magnetic field data
 *   \param[out] vorticity Array for vorticity data
 *   \param[out] photon_density Array for photon density data
 *   \param[out] chem_elements Array for chemical elements data
 *   \param[out] density_trmc Array for MC tracer data
 *   \param[out] rih (required only with SX_OUTPUT_IMAGE) Array for further image data
 *   \param[out] hrih (required only with SX_OUTPUT_IMAGE) Array for further image data
 *   \param[out] sxrates (required only with SX_OUTPUT_IMAGE_ALL) Array for further
 *   image data
 *
 *   \return void
 */
void fill_slice(ray_data *ray, int n_ray, int gradients_flag, int pixels_x, int pixels_y, float *density, float *temperature,
                float *metallicity, float *velocity, float *Bfield, float *vorticity, float *photon_density, float *chem_elements,
                float *density_trmc
#ifdef CHEM_IMAGE
                ,
                float *dust, float *xH2, float *xHP, float *xCO
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
                ,
                float *xCHX, float *xOHX, float *xHCOP, float *xCP, float *xMP
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16 || CHEMISTRYNETWORK == 1
                ,
                float *xHEP
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
                ,
                float *rih, float *hrih
#endif
#ifdef SX_OUTPUT_IMAGE_ALL
                ,
                float *sxrates
#endif
)
{
  float *dens, *temp, *vel, *mag, *vort, *met, *phot, *chem_elem, *denstrmc;
  int i, sph_idx, pix;
#if defined(RT_ADVECT) || defined(GFM_STELLAR_EVOLUTION) || defined(SX_OUTPUT_IMAGE_ALL)
  int j;
#endif
#ifdef CHEM_IMAGE
  float *dst, *H2, *HP, *CO;
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  float *CHX, *OHX, *HCOP, *CP, *MP;
#endif
#if CHEMISTRYNETWORK == 1 || CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  float *HEP;
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
  float *rih2, *hrih2;
#endif
#ifdef SX_OUTPUT_IMAGE_ALL
  float *sxrates2;
#endif
  double l_dx, l_dy, l_dz;

  dens      = (float *)mymalloc("dens", pixels_x * pixels_y * sizeof(float));
  temp      = (float *)mymalloc("temp", pixels_x * pixels_y * sizeof(float));
  met       = (float *)mymalloc("met", pixels_x * pixels_y * sizeof(float));
  vel       = (float *)mymalloc("vel", 3 * pixels_x * pixels_y * sizeof(float));
  mag       = (float *)mymalloc("mag", 3 * pixels_x * pixels_y * sizeof(float));
  vort      = (float *)mymalloc("vort", 3 * pixels_x * pixels_y * sizeof(float));
  phot      = (float *)mymalloc("phot", RT_N_DIR * pixels_x * pixels_y * sizeof(float));
  chem_elem = (float *)mymalloc("chem_elem", GFM_N_CHEM_ELEMENTS * pixels_x * pixels_y * sizeof(float));
  denstrmc  = (float *)mymalloc("denstrmc", pixels_x * pixels_y * sizeof(float));
#ifdef CHEM_IMAGE
  dst = (float *)mymalloc("dust", pixels_x * pixels_y * sizeof(float));
  H2  = (float *)mymalloc("xH2", pixels_x * pixels_y * sizeof(float));
  HP  = (float *)mymalloc("xHP", pixels_x * pixels_y * sizeof(float));
  CO  = (float *)mymalloc("xCO", pixels_x * pixels_y * sizeof(float));
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  CHX  = (float *)mymalloc("CHX", pixels_x * pixels_y * sizeof(float));
  OHX  = (float *)mymalloc("OHX", pixels_x * pixels_y * sizeof(float));
  HCOP = (float *)mymalloc("HCOP", pixels_x * pixels_y * sizeof(float));
  CP   = (float *)mymalloc("CP", pixels_x * pixels_y * sizeof(float));
  MP   = (float *)mymalloc("MP", pixels_x * pixels_y * sizeof(float));
#endif
#if CHEMISTRYNETWORK == 1 || CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  HEP = (float *)mymalloc("HEP", pixels_x * pixels_y * sizeof(float));
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
  rih2  = (float *)mymalloc("rih2", pixels_x * pixels_y * sizeof(float));
  hrih2 = (float *)mymalloc("hrih2", pixels_x * pixels_y * sizeof(float));
#endif
#ifdef SX_OUTPUT_IMAGE_ALL
  sxrates2 = (float *)mymalloc("sxrates2", SX_NIMAGES * pixels_x * pixels_y * sizeof(float));
#endif

  memset(dens, 0, pixels_x * pixels_y * sizeof(float));
  memset(temp, 0, pixels_x * pixels_y * sizeof(float));
  memset(met, 0, pixels_x * pixels_y * sizeof(float));
  memset(vel, 0, 3 * pixels_x * pixels_y * sizeof(float));
  memset(mag, 0, 3 * pixels_x * pixels_y * sizeof(float));
  memset(vort, 0, 3 * pixels_x * pixels_y * sizeof(float));
  memset(phot, 0, RT_N_DIR * pixels_x * pixels_y * sizeof(float));
  memset(chem_elem, 0, GFM_N_CHEM_ELEMENTS * pixels_x * pixels_y * sizeof(float));
  memset(denstrmc, 0, pixels_x * pixels_y * sizeof(float));
#ifdef CHEM_IMAGE
  memset(dst, 0, pixels_x * pixels_y * sizeof(float));
  memset(H2, 0, pixels_x * pixels_y * sizeof(float));
  memset(HP, 0, pixels_x * pixels_y * sizeof(float));
  memset(CO, 0, pixels_x * pixels_y * sizeof(float));
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  memset(CHX, 0, pixels_x * pixels_y * sizeof(float));
  memset(OHX, 0, pixels_x * pixels_y * sizeof(float));
  memset(HCOP, 0, pixels_x * pixels_y * sizeof(float));
  memset(CP, 0, pixels_x * pixels_y * sizeof(float));
  memset(MP, 0, pixels_x * pixels_y * sizeof(float));
#endif
#if CHEMISTRYNETWORK == 1 || CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  memset(HEP, 0, pixels_x * pixels_y * sizeof(float));
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
  memset(rih2, 0, pixels_x * pixels_y * sizeof(float));
  memset(hrih2, 0, pixels_x * pixels_y * sizeof(float));
#endif
#ifdef SX_OUTPUT_IMAGE_ALL
  memset(sxrates2, 0, SX_NIMAGES * pixels_x * pixels_y * sizeof(float));
#endif

  /* Now we simply loop over the pixels *we own* */
  for(i = 0; i < n_ray; ++i)
    {
      sph_idx = ray[i].index;
      pix     = ray[i].pixel;

      l_dx = ray[i].pos[0] - SphP[sph_idx].Center[0];
      l_dy = ray[i].pos[1] - SphP[sph_idx].Center[1];
      l_dz = ray[i].pos[2] - SphP[sph_idx].Center[2];

      l_dx = nearest_x(l_dx);
      l_dy = nearest_y(l_dy);
#if !defined(REFLECTIVE_X) || !defined(REFLECTIVE_Y)
      if(l_dz > boxHalf_Z)
        l_dz -= boxSize_Z;
#endif

      if(gradients_flag == 1)
        dens[pix] = SphP[sph_idx].Density + SphP[sph_idx].Grad.drho[0] * l_dx + SphP[sph_idx].Grad.drho[1] * l_dy +
                    SphP[sph_idx].Grad.drho[2] * l_dz;
      else if(gradients_flag == 0)
        dens[pix] = SphP[sph_idx].Density;
      else
        terminate("gradients_flag != 1 && gradients_flag != 0");

      if(temperature)
        {
#ifdef VORONOI_PROJ_TEMP
          double meanweight = 4.0 / (1. + 3. * HYDROGEN_MASSFRAC) * PROTONMASS;

#if defined(COOLING) && !defined(GRACKLE)
          meanweight = 4.0 / (1. + 3. * HYDROGEN_MASSFRAC + 4. * HYDROGEN_MASSFRAC * SphP[pix].Ne) * PROTONMASS;
#endif
#ifdef VARIABLE_GAMMA
          double temp_in_K =
              (SphP[sph_idx].GammaE - 1.0) * SphP[sph_idx].Utherm / BOLTZMANN * All.UnitEnergy_in_cgs / All.UnitMass_in_g * meanweight;
#else
          double temp_in_K = GAMMA_MINUS1 * SphP[sph_idx].Utherm / BOLTZMANN * All.UnitEnergy_in_cgs / All.UnitMass_in_g * meanweight;
#endif
#ifdef SGCHEM
          double yn, en, yntot;

          yn = SphP[sph_idx].Density * All.UnitDensity_in_cgs / ((1.0 + 4.0 * ABHE) * PROTONMASS);
          en = SphP[sph_idx].Utherm * SphP[sph_idx].Density * All.UnitEnergy_in_cgs / pow(All.UnitLength_in_cm, 3);
#if CHEMISTRYNETWORK == 1
          yntot = (1.0 + ABHE - SphP[sph_idx].TracAbund[IH2] + SphP[sph_idx].TracAbund[IHP] + SphP[sph_idx].TracAbund[IHEP] +
                   SphP[sph_idx].TracAbund[IHEPP]) *
                  yn;
#else
          yntot     = (1.0 + ABHE - SphP[sph_idx].TracAbund[IH2] + SphP[sph_idx].TracAbund[IHP]) * yn;
#endif
#ifdef VARIABLE_GAMMA
          temp_in_K = (SphP[sph_idx].GammaE - 1.0) * en / (yntot * BOLTZMANN);
#else
          temp_in_K = (GAMMA - 1.0) * en / (yntot * BOLTZMANN);
#endif
#endif /* SGCHEM */

          temp[pix] = temp_in_K;
#else
          temp[pix] = SphP[sph_idx].Utherm;
#endif
        }

#ifdef METALS
      if(metallicity)
        met[pix] = SphP[sph_idx].Metallicity;
#endif

      if(velocity)
        {
          if(gradients_flag == 1)
            {
              vel[3 * pix + 0] = P[sph_idx].Vel[0] + SphP[sph_idx].Grad.dvel[0][0] * l_dx + SphP[sph_idx].Grad.dvel[0][1] * l_dy +
                                 SphP[sph_idx].Grad.dvel[0][2] * l_dz;
              vel[3 * pix + 1] = P[sph_idx].Vel[1] + SphP[sph_idx].Grad.dvel[1][0] * l_dx + SphP[sph_idx].Grad.dvel[1][1] * l_dy +
                                 SphP[sph_idx].Grad.dvel[1][2] * l_dz;
              vel[3 * pix + 2] = P[sph_idx].Vel[2] + SphP[sph_idx].Grad.dvel[2][0] * l_dx + SphP[sph_idx].Grad.dvel[2][1] * l_dy +
                                 SphP[sph_idx].Grad.dvel[2][2] * l_dz;
            }
          else
            {
              vel[3 * pix + 0] = P[sph_idx].Vel[0];
              vel[3 * pix + 1] = P[sph_idx].Vel[1];
              vel[3 * pix + 2] = P[sph_idx].Vel[2];
            }
        }

#ifdef MHD
      if(Bfield)
        {
          double mag_x, mag_y, mag_z;
          if(gradients_flag == 1)
            {
              mag_x = SphP[sph_idx].B[0] + SphP[sph_idx].Grad.dB[0][0] * l_dx + SphP[sph_idx].Grad.dB[0][1] * l_dy +
                      SphP[sph_idx].Grad.dB[0][2] * l_dz;
              mag_y = SphP[sph_idx].B[1] + SphP[sph_idx].Grad.dB[1][0] * l_dx + SphP[sph_idx].Grad.dB[1][1] * l_dy +
                      SphP[sph_idx].Grad.dB[1][2] * l_dz;
              mag_z = SphP[sph_idx].B[2] + SphP[sph_idx].Grad.dB[2][0] * l_dx + SphP[sph_idx].Grad.dB[2][1] * l_dy +
                      SphP[sph_idx].Grad.dB[2][2] * l_dz;
            }
          else
            {
              mag_x = SphP[sph_idx].B[0];
              mag_y = SphP[sph_idx].B[1];
              mag_z = SphP[sph_idx].B[2];
            }
          mag[3 * pix + 0] = mag_x;
          mag[3 * pix + 1] = mag_y;
          mag[3 * pix + 2] = mag_z;
        }
#endif

      vort[3 * pix + 0] = SphP[sph_idx].Grad.dvel[2][1] - SphP[sph_idx].Grad.dvel[1][2];
      vort[3 * pix + 1] = SphP[sph_idx].Grad.dvel[0][2] - SphP[sph_idx].Grad.dvel[2][0];
      vort[3 * pix + 2] = SphP[sph_idx].Grad.dvel[1][0] - SphP[sph_idx].Grad.dvel[0][1];

#ifdef RT_ADVECT
      if(photon_density)
        {
          for(j = 0; j < RT_N_DIR; ++j)
            phot[RT_N_DIR * pix + j] = SphP[sph_idx].DensPhot[j];
        }
#endif

#ifdef GFM_STELLAR_EVOLUTION
      if(chem_elements)
        {
          for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
            chem_elem[GFM_N_CHEM_ELEMENTS * pix + j] = SphP[sph_idx].MetalsFraction[j];
        }
#endif

#ifdef TRACER_MC
      if(denstrmc)
        denstrmc[pix] = dens[pix] * get_number_of_tracers(sph_idx) * All.ReferenceTracerMCMass / P[sph_idx].Mass;
#endif

#ifdef CHEM_IMAGE
      dst[pix] += SphP[sph_idx].DustTemp;
      H2[pix] += SphP[sph_idx].TracAbund[IH2];
      HP[pix] += SphP[sph_idx].TracAbund[IHP];
#if CHEMISTRYNETWORK != 1
      CO[pix] += SphP[sph_idx].TracAbund[ICO];
#else
      CO[pix] += 0.0;
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
      CHX[pix] += SphP[sph_idx].TracAbund[ICHX];
      OHX[pix] += SphP[sph_idx].TracAbund[IOHX];
      HCOP[pix] += SphP[sph_idx].TracAbund[IHCOP];
      CP[pix] += SphP[sph_idx].TracAbund[ICP];
      MP[pix] += SphP[sph_idx].TracAbund[IMP];
#endif
#if CHEMISTRYNETWORK == 1 || CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
      HEP[pix] += SphP[sph_idx].TracAbund[IHEP];
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
      rih2[pix] += sx_chem_image1(sph_idx);
      hrih2[pix] += sx_chem_image2(sph_idx);
#endif
#ifdef SX_OUTPUT_IMAGE_ALL
      if(sxrates)
        {
          for(j = 0; j < SX_NIMAGES; j++)
            sxrates2[SX_NIMAGES * pix + j] += sx_chem_image_all(sph_idx, j);
        }
#endif
    }

  /* Now assemble the densities. Since all points except the ones that
   * were set by the requisite tasks are zero, we just sum the contribution */
  MPI_Reduce(dens, density, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(temperature)
    MPI_Reduce(temp, temperature, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(metallicity)
    MPI_Reduce(met, metallicity, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(velocity)
    MPI_Reduce(vel, velocity, 3 * pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(Bfield)
    MPI_Reduce(mag, Bfield, 3 * pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(vorticity)
    MPI_Reduce(vort, vorticity, 3 * pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(photon_density)
    MPI_Reduce(phot, photon_density, RT_N_DIR * pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(chem_elements)
    MPI_Reduce(chem_elem, chem_elements, GFM_N_CHEM_ELEMENTS * pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(density_trmc)
    MPI_Reduce(denstrmc, density_trmc, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#ifdef CHEM_IMAGE
  if(dust)
    MPI_Reduce(dst, dust, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(xH2)
    MPI_Reduce(H2, xH2, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(xHP)
    MPI_Reduce(HP, xHP, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(CO)
    MPI_Reduce(CO, xCO, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  if(xCHX)
    MPI_Reduce(CHX, xCHX, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(xOHX)
    MPI_Reduce(OHX, xOHX, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(xHCOP)
    MPI_Reduce(HCOP, xHCOP, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(xCP)
    MPI_Reduce(CP, xCP, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(xMP)
    MPI_Reduce(MP, xMP, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
#if CHEMISTRYNETWORK == 1 || CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  if(xHEP)
    MPI_Reduce(HEP, xHEP, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
  if(rih)
    MPI_Reduce(rih2, rih, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(hrih)
    MPI_Reduce(hrih2, hrih, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
#ifdef SX_OUTPUT_IMAGE_ALL
  if(sxrates)
    MPI_Reduce(sxrates2, sxrates, SX_NIMAGES * pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

#ifdef SX_OUTPUT_IMAGE_ALL
  myfree(sxrates2);
#endif
#ifdef SX_OUTPUT_IMAGE
  myfree(hrih2);
  myfree(rih2);
#endif

#ifdef CHEM_IMAGE
#if CHEMISTRYNETWORK == 1 || CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  myfree(HEP);
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  myfree(MP);
  myfree(CP);
  myfree(HCOP);
  myfree(OHX);
  myfree(CHX);
#endif
  myfree(CO);
  myfree(HP);
  myfree(H2);
  myfree(dst);
#endif
  myfree(denstrmc);
  myfree(chem_elem);
  myfree(phot);
  myfree(vort);
  myfree(mag);
  myfree(vel);
  myfree(met);
  myfree(temp);
  myfree(dens);
}

/*! \brief Finds next Voronoi cell from a given position in a given direction
 *
 *  Calculates the intersections between a ray at position p0 and
 *  direction dir and the Voronoi faces across all the connections to
 *  the cell and returns the element in the DC array representing the
 *  face the cell will cross. Cell is the index of the primary mesh
 *  cell the ray is in. Previous is the cell that we came from so, if
 *  set, it is ignored in the intersection test (making sure we don't
 *  detect a face we have just entered through due to numerical
 *  truncation.) The distance to the face is returned in length.
 *
 *  Note that unlike the function in Sunrise that this was back-ported
 *  from, it DOES use periodic boundary conditions. If an edge to a
 *  point on the opposite side of the box is encountered, it will be
 *  wrapped to get the correct intersections.
 *
 *  \param[in] T Pointer to tessellation
 *  \param[in] cell Index of cell that point is currently in
 *  \param[in] p0 Starting point
 *  \param[in] dir Direction of next cell search
 *  \param[in] previous The index of the cell on the original task, where the cell is not imported.
 *  \param[out] length Distance traveled of ray in this cell
 *
 *  \return Index of next Voronoi cell
 */
int find_next_voronoi_cell(tessellation *T, int cell, MyDouble p0[3], double dir[3], int previous, double *length)
{
#ifdef VORONOI_DYNAMIC_UPDATE
  point *DP = T->DP;
  // myassert(DP[cell].index >= 0);
  // myassert(DP[cell].index < T->Ndp);

  MyDouble cell_p[3];
  cell_p[0] = P[cell].Pos[0];
  cell_p[1] = P[cell].Pos[1];
  cell_p[2] = P[cell].Pos[2];

  /* if mesh point is across the boundary, wrap it */
  periodic_wrap_point_MyDouble(cell_p, p0);

  MyDouble nb_p[3];
  double m[3];
  double c[3];
  double q[3];
  double s;

  /* init next to -1, which means it didn't cross the face */
  int next = -1;
  /* initialize length to huge */
  *length = HUGE_VAL;

  int edge      = SphP[cell].first_connection;
  int last_edge = SphP[cell].last_connection;
  int iter      = 0;

  while(1)
    {
      ++iter;
      const int neighbor = DC[edge].dp_index;

#if(defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z))
      if(DC[edge].image_flags == 1)  // only ignore the edge we entered through if the edge does not cross the box
        {
          // ignore the edge we entered through
          if((DC[edge].index == previous) && (DC[edge].task == ThisTask))
            {
              if(edge == last_edge)
                break;
              edge = DC[edge].next;
              continue;
            }
        }
#else
      myassert((DC[edge].task != ThisTask) || (DC[edge].index != cell));

      /* ignore the edge we entered through */
      if((DC[edge].index == previous) && (DC[edge].task == ThisTask))
        {
          if(edge == last_edge)
            break;
          edge = DC[edge].next;
          continue;
        }
#endif

      nb_p[0] = DP[neighbor].x;
      nb_p[1] = DP[neighbor].y;
      nb_p[2] = DP[neighbor].z;

      /* if neighbor is across the boundary, wrap it */
      periodic_wrap_point_MyDouble(nb_p, p0);

      int i;
      for(i = 0; i < 3; ++i)
        {
          /* m is the edge midpoint, which is a point on the face plane */
          m[i] = 0.5 * (nb_p[i] + cell_p[i]);
          /* c is the vector from point p0 to m. */
          c[i] = m[i] - p0[i];
          /* q is the edge vector to the neighboring cell, which is a
           * normal vector of the plane */
          q[i] = nb_p[i] - cell_p[i];
        }

      /* sanity check: by construction we know that the point is inside
       * the cell, which means that c.q > 0. If this is not true, it's
       * because some numerical error has put the point on the other
       * side of the face. In this case we short-circuit the process and
       * say it is on the face, the propagation distance is zero (if
       * it's heading towards the face, ie d.q > 0). We must also guard
       * against the case where the ray is perfectly on the face, in
       * which case we will get NaN. In that case we simply ignore the
       * face. */
      double cdotq = c[0] * q[0] + c[1] * q[1] + c[2] * q[2];
      double ddotq = dir[0] * q[0] + dir[1] * q[1] + dir[2] * q[2];

      if(cdotq > 0)
        {
          /* s = c.nq / d.q
           * This is the standard formula for the intersection between a
           * ray and a plane. s is the point where the ray p0 + s * dir
           * intersects the plane which is perpendicular to q and goes
           * through c, i.e. the Voronoi face corresponding to the j-k
           * edge. */
          s = cdotq / ddotq;
        }
      else
        {
          /* point on wrong (outside) side of face. Something's up.
           * If this is due to numerical truncation, distance to face
           * should be small. If it's large, it's likely we aren't in
           * the cell we think we are in. */

          // double dist_to_face = cdotq / sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
          // myassert(dist_to_face > -1e-6);

          if(ddotq > 0)
            {
              /* it's heading away from the cell, so it must have been
               * supposed to intersect this face
               * set s = 0. */
              s = 0;
            }
          else
            /* it's heading into the cell, so it must have come in on
             * this face. ignore the face */
            s = HUGE_VAL;
        }

      if(s >= 0 && s < *length)
        {
          /* The ray intersects the Voronoi face closer to the
           * starting point than any previous points. This is our
           * current candidate exit face. */
          next    = edge;
          *length = s;
        }

      if(edge == last_edge)
        break;

      myassert(edge != DC[edge].next);
      edge = DC[edge].next;
    }

  /* set length to the physical length instead of the fractional length */
  *length *= sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
  return next;
#else
  myassert(0);
  return 0;
#endif
}

/*! \brief Finds next Voronoi cell from a given position in a given direction
 *
 *  This function works exactly as find_next_voronoi_cell, but in addition it
 *  also skips the previous cell when the ray changes the domain.
 *
 *  \param[in] T Pointer to tessellation
 *  \param[in] cell Index of cell that point is currently in
 *  \param[in] p0 Starting point
 *  \param[in] dir Direction of next cell search
 *  \param[in] previous The index of the cell on the original task, where the cell is not imported.
 *  \param[in] task_of_previous The task of the cell 'previous' / the original task.
 *  \param[out] length Distance traveled of ray in this cell
 *
 *  \return Index of next Voronoi cell
 */
int find_next_voronoi_cell2(tessellation *T, int cell, MyDouble p0[3], double dir[3], int previous, int task_of_previous,
                            double *length)
{
#ifdef VORONOI_DYNAMIC_UPDATE
  point *DP = T->DP;
  // myassert(DP[cell].index >= 0);
  // myassert(DP[cell].index < T->Ndp);

  MyDouble cell_p[3];
  cell_p[0] = P[cell].Pos[0];
  cell_p[1] = P[cell].Pos[1];
  cell_p[2] = P[cell].Pos[2];

  /* if mesh point is across the boundary, wrap it */
  periodic_wrap_point_MyDouble(cell_p, p0);

  MyDouble nb_p[3];
  double m[3];
  double c[3];
  double q[3];
  double s;

  /* init next to -1, which means it didn't cross the face */
  int next = -1;
  /* initialize length to huge */
  *length = HUGE_VAL;

  int edge      = SphP[cell].first_connection;
  int last_edge = SphP[cell].last_connection;
  int iter      = 0;

  int ignored = 0;

  while(1)
    {
      ++iter;
      const int neighbor = DC[edge].dp_index;

#if(defined(REFLECTIVE_X) && defined(REFLECTIVE_Y) && defined(REFLECTIVE_Z))

      if(DC[edge].image_flags == 1)  // only ignore the edge we entered through if the edge does not cross the box
        {
          // ignore the edge we entered through
          if((DC[edge].index == previous) && (DC[edge].task == task_of_previous))
            {
              ignored = 1;

              if(edge == last_edge)
                break;
              edge = DC[edge].next;
              continue;
            }
        }
#else

      myassert((DC[edge].task != ThisTask) || (DC[edge].index != cell));

      /* ignore the edge we entered through */
      if((DC[edge].index == previous) && (DC[edge].task == task_of_previous))
        {
          ignored = 1;

          if(edge == last_edge)
            break;
          edge = DC[edge].next;
          continue;
        }

#endif

      nb_p[0] = DP[neighbor].x;
      nb_p[1] = DP[neighbor].y;
      nb_p[2] = DP[neighbor].z;

#if !(defined(REFLECTIVE_X) && defined(REFLECTIVE_Y) && defined(REFLECTIVE_Z))
      /* if neighbor is across the boundary, wrap it */
      periodic_wrap_point_MyDouble(nb_p, p0);
#endif

      int i;
      for(i = 0; i < 3; ++i)
        {
          /* m is the edge midpoint, which is a point on the face plane */
          m[i] = 0.5 * (nb_p[i] + cell_p[i]);
          //*c is the vector from point p0 to m. */
          c[i] = m[i] - p0[i];
          /* q is the edge vector to the neighboring cell, which is a
             normal vector of the plane */
          q[i] = nb_p[i] - cell_p[i];
        }

      /* sanity check: by construction we know that the point is inside
       * the cell, which means that c.q > 0. If this is not true, it's
       * because some numerical error has put the point on the other
       * side of the face. In this case we short-circuit the process and
       * say it is on the face, the propagation distance is zero (if
       * it's heading towards the face, ie d.q > 0). We must also guard
       * against the case where the ray is perfectly on the face, in
       * which case we will get NaN. In that case we simply ignore the
       * face. */
      double cdotq = c[0] * q[0] + c[1] * q[1] + c[2] * q[2];
      double ddotq = dir[0] * q[0] + dir[1] * q[1] + dir[2] * q[2];

      if(cdotq > 0)
        {
          /* s = c.nq / d.q
           * This is the standard formula for the intersection between a
           * ray and a plane. s is the point where the ray p0 + s * dir
           * intersects the plane which is perpendicular to q and goes
           * through c, i.e. the Voronoi face corresponding to the j-k
           * edge. */
          s = cdotq / ddotq;
        }
      else
        {
          /* point on wrong (outside) side of face. Something's up.
           * If this is due to numerical truncation, distance to face
           * should be small. If it's large, it's likely we aren't in
           * the cell we think we are in. */

          // double dist_to_face = cdotq / sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
          // myassert(dist_to_face > -1e-6);

          if(ddotq > 0)
            {
              /* it's heading away from the cell, so it must have been
               * supposed to intersect this face
               * set s=0. */
              s = 0;
            }
          else
            /* it's heading into the cell, so it must have come in on
             * this face. ignore the face */
            s = HUGE_VAL;
        }

      if(s >= 0 && s < *length)
        {
          /* The ray intersects the Voronoi face closer to the
           * starting point than any previous points. This is our
           * current candidate exit face. */
          next    = edge;
          *length = s;
        }

      if(edge == last_edge)
        break;

      myassert(edge != DC[edge].next);
      edge = DC[edge].next;
    }

#if defined(DEBUG) && !(defined(REFLECTIVE_X) && defined(REFLECTIVE_Y) && defined(REFLECTIVE_Z))
  /* assert that the previous cell gets ignored unless the ray just started */
  if(!ignored && (!(SphP[cell].Center[0] == p0[0] && SphP[cell].Center[1] == p0[1] && SphP[cell].Center[2] == p0[2])))
    {
      myassert(0);
    }
#endif

  /* set length to the physical length instead of the fractional length */
  *length *= sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
  return next;
#else
  myassert(0);
  return 0;
#endif
}

#endif /* VORONOI */
