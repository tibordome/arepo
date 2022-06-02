/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/voronoi_makeimage_new.c
 * \date        MM/YYYY
 * \author
 * \brief       More routines for image generation
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <math.h>
#include <mpi.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"

#ifdef VORONOI_NEW_IMAGE
#ifndef VORONOI_DYNAMIC_UPDATE
#error VORONOI_NEW_IMAGE requires VORONOI_DYNAMIC_UPDATE
#endif
#endif

/* arrays for image data */
float *Dens, *Temp, *Met, *Vel, *Weight;
#ifdef MHD
float *Bfield;
#endif
#ifdef TRACER_MC
float *DensMCtr, *WeightMCtr;
#endif
#ifdef CHEM_IMAGE
float *Dust, *xH2, *xHP, *xCO;
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
float *xCHX, *xOHX, *xHCOP, *xCP, *xMP;
#endif
#if CHEMISTRYNETWORK == 1 || CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
float *xHEP;
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
float *rih, *hrih;
#endif

static void get_ray_start_and_end(int i, int j, ray_data *ray, int pixels_x, int pixels_y, int xaxis, int yaxis, int zaxis,
                                  double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);

/*! \brief Initializes projection rays
 *
 *  Creates the rays from the specified position grid, finds the
 *  cells/tasks containing them, and sends them to the correct task.
 *
 *  \param[in] pixels_x Number of pixels in x direction
 *  \param[in] pixels_y Number of pixels in y direction
 *  \param[in] xaxis X-axis of projection
 *  \param[in] yaxis Y-axis of projection
 *  \param[in] zaxis Z-axis of projection
 *  \param[in] xmin Minimum (projection system) x coordinate for projection
 *  \param[in] xmax Maximum (projection system) x coordinate for projection
 *  \param[in] ymin Minimum (projection system) y coordinate for projection
 *  \param[in] ymax Maximum (projection system) y coordinate for projection
 *  \param[in] zmin Minimum (projection system) z coordinate for projection
 *  \param[in] zmax Maximum (projection system) z coordinate for projection
 *
 *  \return void
 */
void setup_rays(int pixels_x, int pixels_y, int xaxis, int yaxis, int zaxis, double xmin, double xmax, double ymin, double ymax,
                double zmin, double zmax)
{
  int i, j;
  mesh_search_data *searchdata;

  /* first all tasks create the rays that lie in their domain */
  Nray = 0;
  for(i = 0; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      {
        if(Nray >= MaxNray)
          terminate("need to realloc Ray-Field");

        /* set position and direction */
        get_ray_start_and_end(i, j, &Ray[Nray], pixels_x, pixels_y, xaxis, yaxis, zaxis, xmin, xmax, ymin, ymax, zmin, zmax);

        Ray[Nray].len   = 0;
        Ray[Nray].pixel = i * pixels_y + j;
        Ray[Nray].index = -1;
        Ray[Nray].prev  = -1;

        /* now check domain */
        peanokey key   = position_to_peanokey(Ray[Nray].pos);
        int no         = peanokey_to_topnode(key);
        Ray[Nray].task = DomainTask[no];
        if(Ray[Nray].task == ThisTask)
          {
            // if ray is on our domain, keep it. otherwise it will be
            // overwritten by the next
            Nray++;
          }
      }

  /* determine the index and task of the mesh cells containing the rays */
  searchdata = (mesh_search_data *)mymalloc("searchdata", Nray * sizeof(mesh_search_data));

  /* fill search array */
  for(i = 0; i < Nray; i++)
    {
      searchdata[i].Pos[0] = Ray[i].pos[0];
      searchdata[i].Pos[1] = Ray[i].pos[1];
      searchdata[i].Pos[2] = Ray[i].pos[2];
    }

  find_nearest_meshpoint_global(searchdata, Nray, 0, 1);

  /* set ray index/task and send rays to correct tasks */
  mpi_distribute_items_from_search(searchdata, Ray, &Nray, &MaxNray, sizeof(*Ray), TAG_DENS_A, offsetof(ray_data, task),
                                   offsetof(ray_data, index));

  myfree(searchdata);
}

#ifdef VORONOI_NEW_IMAGE

// if this is not defined, the definition in voronoi_makeimage.c
void make_3d_voronoi_projected_image(int num, int gradients_flag, int pixels_x, int pixels_y, int xaxis, int yaxis, int zaxis,
                                     double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int weight_flag)
{
  CPU_Step[CPU_MISC] += measure_time();

  float *denssum, *tempsum, *metsum, *velsum, *bmagsum, *weightsum;
  FILE *fd = 0, *fdtemp = 0, *fdmet = 0, *fdvel = 0, *fdmag = 0;

#ifdef TRACER_MC
  float *tracersum, *trweightsum;
  FILE *fdtr = 0;
#endif

#ifdef CHEM_IMAGE
  float *dustsum, *h2sum, *hpsum, *cosum;
  FILE *fddust = 0, *fdh2 = 0, *fdhp = 0, *fdco = 0;
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  float *chxsum, *ohxsum, *hcopsum, *cpsum, *mpsum;
  FILE *fdchx = 0, *fdohx = 0, *fdhcop = 0, *fdcp = 0, *fdmp = 0;
#endif
#if CHEMISTRYNETWORK == 1 || CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  float *hepsum;
  FILE *fdhep = 0;
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
  float *rihsum, *hrihsum;
  FILE *fdrih = 0, *fdhrih = 0;
#endif

  int i;
  char buf[1000];

#ifndef VORONOI_MULTIPLE_PROJECTIONS
  if(gradients_flag == 1)
    sprintf(buf, "proj");
  else if(gradients_flag == 0)
    sprintf(buf, "proj_nograds");
#else
  if(gradients_flag == 1)
    sprintf(buf, "proj_x%d_y%d", xaxis, yaxis);
  else if(gradients_flag == 0)
    sprintf(buf, "proj_nograds_x%d_y%d", xaxis, yaxis);
#endif
  else
    terminate("gradients_flag != 1 && gradients_flag != 0");

  mpi_printf("Generating projected images, %dx%d pixels... gradients_flag=%d\n", pixels_x, pixels_y, gradients_flag);

  open_image_files(buf, num, &fd, &fdtemp, &fdmet, &fdvel,
#ifdef MHD
                   &fdmag,
#else
                   0,
#endif
                   0, 0, 0,
#ifdef TRACER_MC
                   &fdtr,
#else
                   0,
#endif

#ifdef CHEM_IMAGE
                   &fddust, &fdh2, &fdhp, &fdco
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
                   ,
                   &fdchx, &fdohx, &fdhcop, &fdcp, &fdmp, &fdhep
#elif CHEMISTRYNETWORK == 1
                   ,
                   0, 0, 0, 0, 0, &fdhep
#else
                   ,
                   0, 0, 0, 0, 0, 0
#endif
#else
                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#endif

#ifdef SX_OUTPUT_IMAGE
                   ,
                   &fdrih, &fdhrih, 0
#else
                   ,
                   0, 0, 0
#endif
  );

  write_image_header(fd, pixels_x, pixels_y, 0);
  write_image_header(fdtemp, pixels_x, pixels_y, 0);
  write_image_header(fdmet, pixels_x, pixels_y, 0);
  write_image_header(fdvel, pixels_x, pixels_y, 0);
#ifdef MHD
  write_image_header(fdmag, pixels_x, pixels_y, 0);
#endif

#ifdef TRACER_MC
  write_image_header(fdtr, pixels_x, pixels_y, 0);
#endif

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

  /* prepare the image fields that we want to construct */

  Dens = (float *)mymalloc("dens", pixels_x * pixels_y * sizeof(float));
  Temp = (float *)mymalloc("temp", pixels_x * pixels_y * sizeof(float));
  Met  = (float *)mymalloc("met", pixels_x * pixels_y * sizeof(float));
  Vel  = (float *)mymalloc("vel", 3 * pixels_x * pixels_y * sizeof(float));
#ifdef MHD
  Bfield = (float *)mymalloc("bfield", 3 * pixels_x * pixels_y * sizeof(float));
#endif
  Weight = (float *)mymalloc("weight", pixels_x * pixels_y * sizeof(float));

#ifdef TRACER_MC
  DensMCtr   = (float *)mymalloc("densMCtr", pixels_x * pixels_y * sizeof(float));
  WeightMCtr = (float *)mymalloc("weightMCtr", pixels_x * pixels_y * sizeof(float));
#endif

#ifdef CHEM_IMAGE
  Dust = (float *)mymalloc("dust", pixels_x * pixels_y * sizeof(float));
  xH2  = (float *)mymalloc("H2", pixels_x * pixels_y * sizeof(float));
  xHP  = (float *)mymalloc("HP", pixels_x * pixels_y * sizeof(float));
  xCO  = (float *)mymalloc("CO", pixels_x * pixels_y * sizeof(float));
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  xCHX  = (float *)mymalloc("CHX", pixels_x * pixels_y * sizeof(float));
  xOHX  = (float *)mymalloc("OHX", pixels_x * pixels_y * sizeof(float));
  xHCOP = (float *)mymalloc("HCOP", pixels_x * pixels_y * sizeof(float));
  xCP   = (float *)mymalloc("CP", pixels_x * pixels_y * sizeof(float));
  xMP   = (float *)mymalloc("MP", pixels_x * pixels_y * sizeof(float));
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16 || CHEMISTRYNETWORK == 1
  xHEP = (float *)mymalloc("HEP", pixels_x * pixels_y * sizeof(float));
#endif
#endif
#ifdef SX_OUTPUT_IMAGE
  rih  = (float *)mymalloc("rih", pixels_x * pixels_y * sizeof(float));
  hrih = (float *)mymalloc("hrih", pixels_x * pixels_y * sizeof(float));
#endif

  memset(Dens, 0, pixels_x * pixels_y * sizeof(float));
  memset(Temp, 0, pixels_x * pixels_y * sizeof(float));
  memset(Met, 0, pixels_x * pixels_y * sizeof(float));
  memset(Vel, 0, 3 * pixels_x * pixels_y * sizeof(float));
#ifdef MHD
  memset(Bfield, 0, 3 * pixels_x * pixels_y * sizeof(float));
#endif
  memset(Weight, 0, pixels_x * pixels_y * sizeof(float));

#ifdef TRACER_MC
  memset(DensMCtr, 0, pixels_x * pixels_y * sizeof(float));
  memset(WeightMCtr, 0, pixels_x * pixels_y * sizeof(float));
#endif

#ifdef CHEM_IMAGE
  memset(Dust, 0, pixels_x * pixels_y * sizeof(float));
  memset(xH2, 0, pixels_x * pixels_y * sizeof(float));
  memset(xHP, 0, pixels_x * pixels_y * sizeof(float));
  memset(xCO, 0, pixels_x * pixels_y * sizeof(float));
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  memset(xCHX, 0, pixels_x * pixels_y * sizeof(float));
  memset(xOHX, 0, pixels_x * pixels_y * sizeof(float));
  memset(xHCOP, 0, pixels_x * pixels_y * sizeof(float));
  memset(xCP, 0, pixels_x * pixels_y * sizeof(float));
  memset(xMP, 0, pixels_x * pixels_y * sizeof(float));
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16 || CHEMISTRYNETWORK == 1
  memset(xHEP, 0, pixels_x * pixels_y * sizeof(float));
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
  memset(rih, 0, pixels_x * pixels_y * sizeof(float));
  memset(hrih, 0, pixels_x * pixels_y * sizeof(float));
#endif

  MaxNray = pixels_x * pixels_y;
  Ray     = (ray_data *)mymalloc_movable(&Ray, "Ray", MaxNray * sizeof(ray_data));

  measure_time();

  setup_rays(pixels_x, pixels_y, xaxis, yaxis, zaxis, xmin, xmax, ymin, ymax, zmin, zmax);

  MPI_Barrier(MPI_COMM_WORLD);

  printf("task %d has %d rays\n", ThisTask, Nray);

  {
    double t = measure_time();
    mpi_printf("Finding initial mesh cells took %f s.\n", t);
  }

  // double check the positions are correct after exchanging
  // for(i = 0; i < Nray; i++)
  //  {
  //    assert_contains(&Mesh, Ray[i].index, Ray[i].pos);
  //  }

  int left_this_task, rays_left;
  i = 0;
  do
    {
      left_this_task = advance_rays_for_one_cell(1, weight_flag, gradients_flag);
      exchange_rays();
      ++i;

      MPI_Allreduce(&left_this_task, &rays_left, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      /* mpi_printf("iteration %d rays left: %d        \n", i, rays_left); */
    }
  while(rays_left);

  {
    double t = measure_time();
    mpi_printf("integration done after %d iterations and %f s\n", i, t);
  }

  /* let's add up the image contributions and write them to the file we opened */

  denssum = (float *)mymalloc("denssum", pixels_x * pixels_y * sizeof(float));
  tempsum = (float *)mymalloc("tempsum", pixels_x * pixels_y * sizeof(float));
  metsum  = (float *)mymalloc("metsum", pixels_x * pixels_y * sizeof(float));
  velsum  = (float *)mymalloc("velsum", 3 * pixels_x * pixels_y * sizeof(float));
#ifdef MHD
  bmagsum = (float *)mymalloc("bmagsum", 3 * pixels_x * pixels_y * sizeof(float));
#endif
  weightsum = (float *)mymalloc("weightsum", pixels_x * pixels_y * sizeof(float));

#ifdef TRACER_MC
  tracersum   = (float *)mymalloc("tracersum", pixels_x * pixels_y * sizeof(float));
  trweightsum = (float *)mymalloc("trweightsum", pixels_x * pixels_y * sizeof(float));
#endif

#ifdef CHEM_IMAGE
  dustsum = (float *)mymalloc("dustsum", pixels_x * pixels_y * sizeof(float));
  h2sum   = (float *)mymalloc("h2sum", pixels_x * pixels_y * sizeof(float));
  hpsum   = (float *)mymalloc("hpsum", pixels_x * pixels_y * sizeof(float));
  cosum   = (float *)mymalloc("cosum", pixels_x * pixels_y * sizeof(float));
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  chxsum  = (float *)mymalloc("chxsum", pixels_x * pixels_y * sizeof(float));
  ohxsum  = (float *)mymalloc("ohxsum", pixels_x * pixels_y * sizeof(float));
  hcopsum = (float *)mymalloc("hcopsum", pixels_x * pixels_y * sizeof(float));
  cpsum   = (float *)mymalloc("cpsum", pixels_x * pixels_y * sizeof(float));
  mpsum   = (float *)mymalloc("mpsum", pixels_x * pixels_y * sizeof(float));
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16 || CHEMISTRYNETWORK == 1
  hepsum = (float *)mymalloc("hepsum", pixels_x * pixels_y * sizeof(float));
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
  rihsum  = (float *)mymalloc("rihsum", pixels_x * pixels_y * sizeof(float));
  hrihsum = (float *)mymalloc("hrihsum", pixels_x * pixels_y * sizeof(float));
#endif

  MPI_Reduce(Dens, denssum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Temp, tempsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Met, metsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(Vel, velsum, 3 * pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#ifdef MHD
  MPI_Reduce(Bfield, bmagsum, 3 * pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
  MPI_Reduce(Weight, weightsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

#ifdef TRACER_MC
  MPI_Reduce(DensMCtr, tracersum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(WeightMCtr, trweightsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

#ifdef CHEM_IMAGE
  MPI_Reduce(Dust, dustsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(xH2, h2sum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(xHP, hpsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(xCO, cosum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  MPI_Reduce(xCHX, chxsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(xOHX, ohxsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(xHCOP, hcopsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(xCP, cpsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(xMP, mpsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16 || CHEMISTRYNETWORK == 1
  MPI_Reduce(xHEP, hepsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
  MPI_Reduce(rih, rihsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(hrih, hrihsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

  if(ThisTask == 0)
    {
      // divide out the column density weighting

      for(i = 0; i < pixels_x * pixels_y; i++)
        {
          if(denssum[i] > 0)
            {
              tempsum[i] /= denssum[i];
              metsum[i] /= denssum[i];
              velsum[3 * i + 0] /= denssum[i];
              velsum[3 * i + 1] /= denssum[i];
              velsum[3 * i + 2] /= denssum[i];
#ifdef MHD
              bmagsum[3 * i + 0] /= denssum[i];
              bmagsum[3 * i + 1] /= denssum[i];
              bmagsum[3 * i + 2] /= denssum[i];
#ifndef MHD_DONT_PRINT_BMAGSUM
              printf("bmagsum: pixel %d x %g y %g z %g \n", i, bmagsum[3 * i + 0], bmagsum[3 * i + 1], bmagsum[3 * i + 2]);
#endif
#endif

#ifdef CHEM_IMAGE
              dustsum[i] /= denssum[i];
              h2sum[i] /= denssum[i];
              hpsum[i] /= denssum[i];
              cosum[i] /= denssum[i];
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
              chxsum[i] /= denssum[i];
              ohxsum[i] /= denssum[i];
              hcopsum[i] /= denssum[i];
              cpsum[i] /= denssum[i];
              mpsum[i] /= denssum[i];
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16 || CHEMISTRYNETWORK == 1
              hepsum[i] /= denssum[i];
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
              rihsum[i] /= denssum[i];
              hrihsum[i] /= denssum[i];
#endif
            }

          if(weightsum[i] > 0)
            denssum[i] /= weightsum[i];

#ifdef TRACER_MC
          if(trweightsum[i] > 0)
            tracersum[i] /= trweightsum[i];
#endif
        }

      my_fwrite(denssum, sizeof(float), pixels_x * pixels_y, fd);
      my_fwrite(tempsum, sizeof(float), pixels_x * pixels_y, fdtemp);
      my_fwrite(metsum, sizeof(float), pixels_x * pixels_y, fdmet);
      my_fwrite(velsum, 3 * sizeof(float), pixels_x * pixels_y, fdvel);
#ifdef MHD
      my_fwrite(bmagsum, 3 * sizeof(float), pixels_x * pixels_y, fdmag);
#endif
#ifdef TRACER_MC
      my_fwrite(tracersum, sizeof(float), pixels_x * pixels_y, fdtr);
#endif

#ifdef CHEM_IMAGE
      my_fwrite(dustsum, sizeof(float), pixels_x * pixels_y, fddust);
      my_fwrite(h2sum, sizeof(float), pixels_x * pixels_y, fdh2);
      my_fwrite(hpsum, sizeof(float), pixels_x * pixels_y, fdhp);
      my_fwrite(cosum, sizeof(float), pixels_x * pixels_y, fdco);
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
      my_fwrite(chxsum, sizeof(float), pixels_x * pixels_y, fdchx);
      my_fwrite(ohxsum, sizeof(float), pixels_x * pixels_y, fdohx);
      my_fwrite(hcopsum, sizeof(float), pixels_x * pixels_y, fdhcop);
      my_fwrite(cpsum, sizeof(float), pixels_x * pixels_y, fdcp);
      my_fwrite(mpsum, sizeof(float), pixels_x * pixels_y, fdmp);
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16 || CHEMISTRYNETWORK == 1
      my_fwrite(hepsum, sizeof(float), pixels_x * pixels_y, fdhep);
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
      my_fwrite(rihsum, sizeof(float), pixels_x * pixels_y, fdrih);
      my_fwrite(hrihsum, sizeof(float), pixels_x * pixels_y, fdhrih);
#endif

#ifdef IMAGE_FOOTERS
      write_image_footer(fd, xmin, xmax, ymin, ymax, zmin, zmax);
      write_image_footer(fdtemp, xmin, xmax, ymin, ymax, zmin, zmax);
      write_image_footer(fdmet, xmin, xmax, ymin, ymax, zmin, zmax);
      write_image_footer(fdvel, xmin, xmax, ymin, ymax, zmin, zmax);
#ifdef MHD
      write_image_footer(fdmag, xmin, xmax, ymin, ymax, zmin, zmax);
#endif

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
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16 || CHEMISTRYNETWORK == 1
      write_image_footer(fdhep, xmin, xmax, ymin, ymax, zmin, zmax);
#endif
#endif
#endif

      fclose(fd);
      fclose(fdtemp);
      fclose(fdmet);
      fclose(fdvel);
#ifdef MHD
      fclose(fdmag);
#endif
#ifdef TRACER_MC
      fclose(fdtr);
#endif

#ifdef CHEM_IMAGE
      fclose(fddust);
      fclose(fdh2);
      fclose(fdhp);
      fclose(fdco);
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
      fclose(fdchx);
      fclose(fdohx);
      fclose(fdhcop);
      fclose(fdcp);
      fclose(fdmp);
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16 || CHEMISTRYNETWORK == 1
      fclose(fdhep);
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
      fclose(fdrih);
      fclose(fdhrih);
#endif
    }

#ifdef SX_OUTPUT_IMAGE
  myfree(hrihsum);
  myfree(rihsum);
#endif

#ifdef CHEM_IMAGE
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16 || CHEMISTRYNETWORK == 1
  myfree(hepsum);
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  myfree(mpsum);
  myfree(cpsum);
  myfree(hcopsum);
  myfree(ohxsum);
  myfree(chxsum);
#endif
  myfree(cosum);
  myfree(hpsum);
  myfree(h2sum);
  myfree(dustsum);
#endif

#ifdef TRACER_MC
  myfree(trweightsum);
  myfree(tracersum);
#endif

  myfree(weightsum);
#ifdef MHD
  myfree(bmagsum);
#endif
  myfree(velsum);
  myfree(metsum);
  myfree(tempsum);
  myfree(denssum);

  myfree(Ray);

#ifdef SX_OUTPUT_IMAGE
  myfree(hrih);
  myfree(rih);
#endif

#ifdef CHEM_IMAGE
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16 || CHEMISTRYNETWORK == 1
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
  myfree(Dust);
#endif

#ifdef TRACER_MC
  myfree(WeightMCtr);
  myfree(DensMCtr);
#endif

  myfree(Weight);
#ifdef MHD
  myfree(Bfield);
#endif
  myfree(Vel);
  myfree(Met);
  myfree(Temp);
  myfree(Dens);

  CPU_Step[CPU_MAKEIMAGES] += measure_time();
}
#endif

/*! \brief Send rays to their tasks if they aren't there already.
 *
 *  \return void
 */
void exchange_rays(void)
{
  // reset prev field for rays that will be moved to another task
  int i;
  for(i = 0; i < Nray; i++)
    {
      if(Ray[i].task != ThisTask)
        Ray[i].prev = -1;
    }

  mpi_distribute_items_to_tasks(Ray, offsetof(ray_data, task), &Nray, &MaxNray, sizeof(*Ray), TAG_DENS_A);
}

/*! \brief Updates the column values for the specified pixel.
 *
 *  \param[in] cell Index of cell in P and SphP array
 *  \param[in] rho Density
 *  \param[in] len Length of ray in this cell
 *  \param[in] pix Index of pixel of interest
 *  \param[in] weight_flag Density^2 weighted (==1) or volume weighted (==else)
 *
 *  \return void
 */
void update_column_integral(int cell, MyFloat rho, double len, int pix, int weight_flag)
{
  myassert(Dens[pix] == Dens[pix]);
  myassert(rho == rho);
  myassert(len >= 0);

  if(weight_flag == 1)
    {
      Dens[pix] += rho * rho * rho * len;
#ifdef VORONOI_PROJ_TEMP
#ifdef GRACKLE
      double temp_in_K = get_temp_individual_cell_grackle(cell);
#else
      double meanweight = 4. / (1. + 3. * HYDROGEN_MASSFRAC) * PROTONMASS;
#ifdef COOLING
      meanweight        = 4.0 / (1. + 3. * HYDROGEN_MASSFRAC + 4. * HYDROGEN_MASSFRAC * SphP[cell].Ne) * PROTONMASS;
#endif
#ifdef VARIABLE_GAMMA
      double temp_in_K =
          (SphP[cell].GammaE - 1.0) * SphP[cell].Utherm / BOLTZMANN * All.UnitEnergy_in_cgs / All.UnitMass_in_g * meanweight;
#else
      double temp_in_K = GAMMA_MINUS1 * SphP[cell].Utherm / BOLTZMANN * All.UnitEnergy_in_cgs / All.UnitMass_in_g * meanweight;
#endif
#endif

#ifdef SGCHEM
      double yn, en, yntot;

      yn    = rho * All.UnitDensity_in_cgs / ((1.0 + 4.0 * ABHE) * PROTONMASS);
      en    = SphP[cell].Utherm * rho * All.UnitEnergy_in_cgs / pow(All.UnitLength_in_cm, 3);
      yntot = (1.0 + ABHE - SphP[cell].TracAbund[0] + SphP[cell].TracAbund[1]) * yn;
#ifdef VARIABLE_GAMMA
      temp_in_K = (SphP[cell].GammaE - 1.0) * en / (yntot * BOLTZMANN);
#else
      temp_in_K = 2.0 * en / (3.0 * yntot * BOLTZMANN);
#endif
#endif

      Temp[pix] += temp_in_K * rho * rho * rho * len;
#else
      Temp[pix] += SphP[cell].Utherm * rho * rho * rho * len;
#endif
      Weight[pix] += rho * rho * len;
#ifdef TRACER_MC
      DensMCtr[pix] += pow(get_number_of_tracers(cell) * All.ReferenceTracerMCMass / P[cell].Mass, 3.0) * rho * rho * rho * len;
      WeightMCtr[pix] += pow(get_number_of_tracers(cell) * All.ReferenceTracerMCMass / P[cell].Mass, 2.0) * rho * rho * len;
#endif
    }
  else
    {
      Dens[pix] += rho * len;
#ifdef VORONOI_PROJ_TEMP
#ifdef GRACKLE
      double temp_in_K = get_temp_individual_cell_grackle(cell);
#else
      double meanweight = 4 / (1 + 3 * HYDROGEN_MASSFRAC) * PROTONMASS;
#ifdef COOLING
      meanweight        = 4.0 / (1. + 3. * HYDROGEN_MASSFRAC + 4. * HYDROGEN_MASSFRAC * SphP[cell].Ne) * PROTONMASS;
#endif
#ifdef VARIABLE_GAMMA
      double temp_in_K =
          (SphP[cell].GammaE - 1.0) * SphP[cell].Utherm / BOLTZMANN * All.UnitEnergy_in_cgs / All.UnitMass_in_g * meanweight;
#else
      double temp_in_K = GAMMA_MINUS1 * SphP[cell].Utherm / BOLTZMANN * All.UnitEnergy_in_cgs / All.UnitMass_in_g * meanweight;
#endif
#endif

#ifdef SGCHEM
      double yn, en, yntot;
      yn    = rho * All.UnitDensity_in_cgs / ((1.0 + 4.0 * ABHE) * PROTONMASS);
      en    = SphP[cell].Utherm * rho * All.UnitEnergy_in_cgs / pow(All.UnitLength_in_cm, 3);
      yntot = (1.0 + ABHE - SphP[cell].TracAbund[IH2] + SphP[cell].TracAbund[IHP]) * yn;
#ifdef VARIABLE_GAMMA
      temp_in_K = (SphP[cell].GammaE - 1.0) * en / (yntot * BOLTZMANN);
#else
      temp_in_K = 2.0 * en / (3.0 * yntot * BOLTZMANN);
#endif
#endif

      Temp[pix] += temp_in_K * rho * len;
#else
      Temp[pix] += SphP[cell].Utherm * rho * len;
#endif

#ifdef TRACER_MC
      DensMCtr[pix] += get_number_of_tracers(cell) * All.ReferenceTracerMCMass / P[cell].Mass * rho * len;
#endif
    }

#if defined(METALS) || defined(GFM_STELLAR_EVOLUTION)
  Met[pix] += SphP[cell].Metallicity * rho * len;
#endif

#ifdef SGCHEM
#ifdef CHEM_IMAGE
  Dust[pix] += SphP[cell].DustTemp * rho * len;
  xH2[pix] += SphP[cell].TracAbund[IH2] * rho * len;
  xHP[pix] += SphP[cell].TracAbund[IHP] * rho * len;
#if CHEMISTRYNETWORK != 1
  xCO[pix] += SphP[cell].TracAbund[ICO] * rho * len;
#else
  xCO[pix] += 0.0;
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  xCHX[pix] += SphP[cell].TracAbund[ICHX] * rho * len;
  xOHX[pix] += SphP[cell].TracAbund[IOHX] * rho * len;
  xHCOP[pix] += SphP[cell].TracAbund[IHCOP] * rho * len;
  xCP[pix] += SphP[cell].TracAbund[ICP] * rho * len;
  xMP[pix] += SphP[cell].TracAbund[IMP] * rho * len;
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16 || CHEMISTRYNETWORK == 1
  xHEP[pix] += SphP[cell].TracAbund[IHEP] * rho * len;
#endif
#endif
#endif

#ifdef SX_OUTPUT_IMAGE
  rih[pix] += sx_chem_image1(cell) * rho * len;
  hrih[pix] += sx_chem_image2(cell) * rho * len;
#endif

  Vel[3 * pix + 0] += P[cell].Vel[0] * rho * len;
  Vel[3 * pix + 1] += P[cell].Vel[1] * rho * len;
  Vel[3 * pix + 2] += P[cell].Vel[2] * rho * len;

#ifdef MHD
  Bfield[3 * pix + 0] += SphP[cell].B[0] * rho * len;
  Bfield[3 * pix + 1] += SphP[cell].B[1] * rho * len;
  Bfield[3 * pix + 2] += SphP[cell].B[2] * rho * len;
#endif
  myassert(Dens[pix] == Dens[pix]);
}

/*! \brief Advance the rays
 *
 *  Step the rays through the current cell and find next cell they
 *  will enter. If integrate is true, it also updates the column
 *  integral values. Returns the number of rays that have not yet
 *  reached their target_len.
 *
 *  \param[in] integrate Flag if this is a column integral
 *  \param[in] weight_flag Density^2 weighted (==1) or volume weighted (==else)
 *  \param[in] gradients_flag Use gradients for projection?
 *
 *  \return Number of rays that are left for next iteration (= did not yet hit
 *  a boundary)
 */
int advance_rays_for_one_cell(int integrate, int weight_flag, int gradients_flag)
{
#if !defined(VORONOI) || !defined(VORONOI_DYNAMIC_UPDATE)
  terminate("Ray tracing only works with VORONOI and DYNAMIC_UPDATE enabled.");
#else
  int i, j, nleft = 0;

  for(i = 0; i < Nray; i++)
    {
      myassert(Ray[i].task == ThisTask);

      if(Ray[i].len >= Ray[i].target_len)
        /* ray done */
        continue;

      /* this is the index of the cell in which we currently are */
      int sph_idx = Ray[i].index;

      myassert((sph_idx >= 0) && (sph_idx < NumGas));

      double l;

      int next_edge = find_next_voronoi_cell(&Mesh, sph_idx, Ray[i].pos, Ray[i].dir, Ray[i].prev, &l);

      if(Ray[i].len + l > Ray[i].target_len)
        {
          /* The next cell boundary is past the target length. truncate
             propagation and mark ray as done. do not change the
             cell/task. set prev to -1 since we are inside the cell */
          l = Ray[i].target_len - Ray[i].len;
          Ray[i].len = Ray[i].target_len;
          Ray[i].prev = -1;
        }
      else
        {
          /* we will enter next cell. update cell info */
          Ray[i].task = DC[next_edge].task;
          Ray[i].prev = Ray[i].index;
          Ray[i].index = DC[next_edge].index;
          Ray[i].len += l;
          ++nleft;
        }

      myassert(next_edge >= 0);

#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
      double out_of_box_flag = 0;
#endif

#ifdef REFLECTIVE_X
      double posx_old = Ray[i].pos[0];
#endif
#ifdef REFLECTIVE_Y
      double posy_old = Ray[i].pos[1];
#endif
#ifdef REFLECTIVE_Z
      double posz_old = Ray[i].pos[2];
#endif

      /* now determine the interpolated value at midpoint of segment */
      double delta_center[3];

      for(j = 0; j < 3; j++)
        {
          delta_center[j] = NEAREST_X(Ray[i].pos[j] + 0.5 * l * Ray[i].dir[j] - SphP[sph_idx].Center[j]);
          Ray[i].pos[j] += l * Ray[i].dir[j];
        }

#ifdef REFLECTIVE_X
      if(Ray[i].pos[0] < 0)
        Ray[i].pos[0] = 0;
      if(Ray[i].pos[0] > boxSize_X)
        Ray[i].pos[0] = boxSize_X;

      delta_center[0] = 0.5 * (posx_old + Ray[i].pos[0]);
      out_of_box_flag = 1;
#endif
#ifdef REFLECTIVE_Y
      if(Ray[i].pos[1] < 0)
        Ray[i].pos[1] = 0;
      if(Ray[i].pos[1] > boxSize_Y)
        Ray[i].pos[1] = boxSize_Y;

      delta_center[1] = 0.5 * (posy_old + Ray[i].pos[1]);
      out_of_box_flag = 1;
#endif
#ifdef REFLECTIVE_Z
      if(Ray[i].pos[2] < 0)
        Ray[i].pos[2] = 0;
      if(Ray[i].pos[2] > boxSize_Z)
        Ray[i].pos[2] = boxSize_Z;

      delta_center[2] = 0.5 * (posz_old + Ray[i].pos[2]);
      out_of_box_flag = 1;
#endif

#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
      if(out_of_box_flag)
        {
          Ray[i].target_len = Ray[i].len;
          nleft--;
        }
#endif

      /* check that ray is still in cell */
      // assert_contains(&Mesh, Ray[i].index, Ray[i].pos);

      MyFloat rho;
      if(gradients_flag == 1)
        rho = SphP[sph_idx].Density + SphP[sph_idx].Grad.drho[0] * delta_center[0] + SphP[sph_idx].Grad.drho[1] * delta_center[1] +
              SphP[sph_idx].Grad.drho[2] * delta_center[2];
      else if(gradients_flag == 0)
        rho = SphP[sph_idx].Density;
      else
        terminate("gradients_flag != 0 && gradients_flag != 1");

      if(integrate)
        update_column_integral(sph_idx, rho, l, Ray[i].pixel, weight_flag);

      /* see if we need to wrap ray around boundary. (This is hoaky,
       * because if the point struct is changed it will break.) */
      MyDouble ref[3];
      for(j = 0; j < 3; j++)
        ref[j] = (&Mesh.DP[DC[next_edge].dp_index].x)[j];
      periodic_wrap_point_MyDouble(Ray[i].pos, ref);
    }
  return nleft;
#endif
}

/*! \brief Set start and end points for projection
 *
 *  So far this only works along simulation box coordinates
 *
 *  \param[in] i Pixel index in X direction
 *  \param[in] j Pixel index in Y direction
 *  \param[out] ray Ray data this is initialized in
 *  \param[in] pixels_x Number of pixels in x direction
 *  \param[in] pixels_y Number of pixels in y direction
 *  \param[in] xaxis Index of the projection x axis in simulation coordinates
 *  \param[in] yaxis Index of the projection y axis in simulation coordinates
 *  \param[in] zaxis Index of the projection z axis in simulation coordinates
 *  \param[in] xmin Minimum x coordinate in projection coordinates
 *  \param[in] xmax Maximum x coordinate in projection coordinates
 *  \param[in] ymin Minimum y coordinate in projection coordinates
 *  \param[in] ymax Maximum y coordinate in projection coordinates
 *  \param[in] zmin Minumum z coordinate in projection coordinates
 *  \param[in] zmax Maximum z coordinate in projection coordinates
 *
 *  \return void
 */
void get_ray_start_and_end(int i, int j, ray_data *ray, int pixels_x, int pixels_y, int xaxis, int yaxis, int zaxis, double xmin,
                           double xmax, double ymin, double ymax, double zmin, double zmax)
{
  double dir      = zmax > zmin ? 1.0 : -1.0;
  ray->target_len = fabs(zmax - zmin);

  if(xaxis == 0 && yaxis == 1)
    {
      ray->pos[0] = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
      ray->pos[1] = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
      ray->pos[2] = zmin;

      ray->dir[0] = 0;
      ray->dir[1] = 0;
      ray->dir[2] = dir;
    }
  else if(xaxis == 1 && yaxis == 0)
    {
      ray->pos[0] = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
      ray->pos[1] = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
      ray->pos[2] = zmin;

      ray->dir[0] = 0;
      ray->dir[1] = 0;
      ray->dir[2] = dir;
    }
  else if(xaxis == 0 && yaxis == 2)
    {
      ray->pos[0] = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
      ray->pos[1] = zmin;
      ray->pos[2] = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;

      ray->dir[0] = 0;
      ray->dir[1] = dir;
      ray->dir[2] = 0;
    }
  else if(xaxis == 2 && yaxis == 0)
    {
      ray->pos[0] = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
      ray->pos[1] = zmin;
      ray->pos[2] = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;

      ray->dir[0] = 0;
      ray->dir[1] = dir;
      ray->dir[2] = 0;
    }
  else if(xaxis == 1 && yaxis == 2)
    {
      ray->pos[0] = zmin;
      ray->pos[1] = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
      ray->pos[2] = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;

      ray->dir[0] = dir;
      ray->dir[1] = 0;
      ray->dir[2] = 0;
    }
  else
    {
      ray->pos[0] = zmin;
      ray->pos[1] = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
      ray->pos[2] = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;

      ray->dir[0] = dir;
      ray->dir[1] = 0;
      ray->dir[2] = 0;
    }
}
