/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/movie_auriga/movie_projection.c
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

#include <stdio.h>
#include <sys/stat.h>

#include "../allvars.h"
#include "../proto.h"
#include "movie.h"

#ifdef AURIGA_MOVIE
static void auriga_movie_projection_gas(struct auriga_movie_projection *self);
static void auriga_movie_projection_gas_fast(struct auriga_movie_projection *self);
static void auriga_movie_projection_dark_matter(struct auriga_movie_projection *self);
static void auriga_movie_projection_stars(struct auriga_movie_projection *self);

static void auriga_movie_projection_vector_shift(double v_new[3], double v_orig[3], double v_shift[3])
{
  v_new[0] = NEAREST_X(v_orig[0] - v_shift[0]);
  v_new[1] = NEAREST_Y(v_orig[1] - v_shift[1]);
  v_new[2] = NEAREST_Z(v_orig[2] - v_shift[2]);
}

static void auriga_movie_projection_vector_shift_non_periodic(double v_new[3], double v_orig[3], double v_shift[3])
{
  v_new[0] = v_orig[0] - v_shift[0];
  v_new[1] = v_orig[1] - v_shift[1];
  v_new[2] = v_orig[2] - v_shift[2];
}

static void auriga_movie_projection_vector_rotate(double vector[3], double rotation[3][3])
{
  double tx = vector[0];
  double ty = vector[1];
  double tz = vector[2];

  vector[0] = tx * rotation[0][0] + ty * rotation[0][1] + tz * rotation[0][2];
  vector[1] = tx * rotation[1][0] + ty * rotation[1][1] + tz * rotation[1][2];
  vector[2] = tx * rotation[2][0] + ty * rotation[2][1] + tz * rotation[2][2];
}

static double auriga_movie_projection_vector_radius(double vector[3])
{
  return sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
}

static void auriga_movie_projection_vector_normalize(double vector[3])
{
  double length = auriga_movie_projection_vector_radius(vector);
  int k;
  if(length > 0)
    for(k = 0; k < 3; k++)
      vector[k] /= length;
}

static double auriga_movie_projection_dot_product(double vector1[3], double vector2[3])
{
  return vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2];
}

inline double _getkernel(double h, double r2)
{
  double coeff1, coeff2, coeff5;
  double hinv, hinv3, u;
  coeff1 = 8.0 / M_PI;
  coeff2 = coeff1 * 6.0;
  coeff5 = coeff1 * 2.0;

  hinv  = 1.0 / h;
  hinv3 = hinv * hinv * hinv;
  u     = sqrt(r2) * hinv;

  if(u > 1.0)
    return 0.;

  if(u < 0.5)
    {
      return hinv3 * (coeff1 + coeff2 * (u - 1.0) * u * u);
    }
  else
    {
      return hinv3 * coeff5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
    }
}

void auriga_movie_projection_init(struct auriga_movie_projection *self, double center[3], double rotation[3][3], double cmvel[3],
                                  double width, double depth, int nx, int ny)
{
  int k, ix, iy;

  for(k = 0; k < 3; k++)
    self->center[k] = center[k];

  if(rotation)
    {
      for(ix = 0; ix < 3; ix++)
        for(iy = 0; iy < 3; iy++)
          self->rotation[ix][iy] = rotation[ix][iy];
    }
  else
    {
      for(ix = 0; ix < 3; ix++)
        for(iy = 0; iy < 3; iy++)
          {
            if(ix == iy)
              self->rotation[ix][iy] = 1;
            else
              self->rotation[ix][iy] = 0;
          }
    }

  if(cmvel)
    {
      for(k = 0; k < 3; k++)
        self->cmvel[k] = cmvel[k];
    }
  else
    {
      for(k = 0; k < 3; k++)
        self->cmvel[k] = 0;
    }

  self->nx = nx;
  self->ny = ny;

  self->width  = width;
  self->height = width / nx * ny;
  self->depth  = depth;

  if(width == 0)
    {
      self->allsky = 1;
      self->ny     = self->nx / 2;
      self->width  = 2. * M_PI;
      self->height = M_PI;
    }
  else
    self->allsky = 0;

  for(k = 0; k < aum_count; k++)
    {
      self->projection[k] = (double *)mymalloc("projection", self->nx * self->ny * sizeof(double));
      memset(self->projection[k], 0, self->nx * self->ny * sizeof(double));
    }
}

void auriga_movie_projection_project_particles(struct auriga_movie_projection *self)
{
  double dx = self->width / self->nx;
  double dy = self->height / self->ny;

  mpi_printf(
      "AURIGA MOVIE: Doing projections, center=%g,%g,%g, vel=%g,%g,%g, width=%g, height=%g, depth=%g, nx=%d, ny=%d, dx=%g, dy=%g.\n",
      self->center[0], self->center[1], self->center[2], self->cmvel[0], self->cmvel[1], self->cmvel[2], self->width, self->height,
      self->depth, self->nx, self->ny, dx, dy);

  double t0 = second();

  /* gas first */
  mpi_printf("AURIGA MOVIE: Doing gas projection.\n");
  auriga_movie_projection_gas_fast(self);
  double t1 = second();
  mpi_printf("AURIGA MOVIE: Gas projection done in %gs.\n", timediff(t0, t1));

  /* then dark matter */
  mpi_printf("AURIGA MOVIE: Doing DM projection.\n");
  auriga_movie_projection_dark_matter(self);
  double t2 = second();
  mpi_printf("AURIGA MOVIE: DM projection done in %gs.\n", timediff(t1, t2));

  /* then stars */
  mpi_printf("AURIGA MOVIE: Doing stellar projection.\n");
  auriga_movie_projection_stars(self);
  double t3 = second();
  mpi_printf("AURIGA MOVIE: Stellar projection done in %gs.\n", timediff(t2, t3));

  mpi_printf("AURIGA MOVIE: Collecting projections on task 0.\n");
  int k;
  for(k = 0; k < aum_count; k++)
    {
      if(ThisTask == 0)
        MPI_Reduce(MPI_IN_PLACE, self->projection[k], self->nx * self->ny, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      else
        MPI_Reduce(self->projection[k], NULL, self->nx * self->ny, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

  double t4 = second();
  mpi_printf("AURIGA MOVIE: Collected projections on task 0 in %gs.\n", timediff(t3, t4));

  mpi_printf("AURIGA MOVIE: Doing final operations on projections.\n");
  if(ThisTask == 0)
    {
      double area = dx * dy;

      int k;
      for(k = 0; k < self->nx * self->ny; k++)
        {
          if(self->projection[aum_gas_rho][k] > 0)
            {
              self->projection[aum_gas_vel][k] = sqrt(self->projection[aum_gas_vel][k] / self->projection[aum_gas_rho][k]);
              self->projection[aum_gas_velx][k] /= self->projection[aum_gas_rho][k];
              self->projection[aum_gas_vely][k] /= self->projection[aum_gas_rho][k];
              self->projection[aum_gas_velz][k] /= self->projection[aum_gas_rho][k];

              self->projection[aum_gas_temp][k] /= self->projection[aum_gas_rho][k];
              self->projection[aum_gas_metallicity][k] /= self->projection[aum_gas_rho][k];
            }

          self->projection[aum_gas_pres][k] = sqrt(self->projection[aum_gas_pres][k] / self->depth);

#ifdef MHD
          self->projection[aum_gas_bfld][k] = sqrt(self->projection[aum_gas_bfld][k] / self->depth);
          self->projection[aum_gas_bx][k] /= self->depth;
          self->projection[aum_gas_by][k] /= self->depth;
          self->projection[aum_gas_bz][k] /= self->depth;
#endif
          self->projection[aum_dm_rho][k] /= area;

          if(self->projection[aum_stars_rho][k] > 0)
            self->projection[aum_stars_metallicity][k] /= self->projection[aum_stars_rho][k];

          self->projection[aum_stars_rho][k] /= area;
        }
    }
  double t5 = second();
  mpi_printf("AURIGA MOVIE: Final operations on projections done in %gs.\n", timediff(t4, t5));
}

void auriga_movie_projection_gas_fast(struct auriga_movie_projection *self)
{
  double timetotal        = 0;
  double timetonextupdate = 0;
  double previous         = second();

  double dx = self->width / self->nx;
  double dy = self->height / self->ny;

  double *len = (double *)mymalloc("len", self->nx * self->ny * sizeof(double));
  memset(len, 0, self->nx * self->ny * sizeof(double));

  double xoffset = 0;
  if(self->allsky)
    xoffset = -8e-3 / All.Time * All.HubbleParam;

  int igas;
  for(igas = 0; igas < NumGas; igas++)
    {
      if(P[igas].Type != 0)
        continue;

      double pos[3];
      auriga_movie_projection_vector_shift(pos, P[igas].Pos, self->center);
      auriga_movie_projection_vector_rotate(pos, self->rotation);

      pos[0] -= xoffset;

      double maxrad;
      if(self->allsky)
        {
          double radius = auriga_movie_projection_vector_radius(pos);
          maxrad        = fmin(asin(SphP[igas].MaxDelaunayRadius / radius), 2. * M_PI);

          if(radius - SphP[igas].MaxDelaunayRadius > self->depth)
            continue;

          if(SphP[igas].MaxDelaunayRadius > radius * 0.5)
            continue;

          double theta = acos(pos[2] / radius);
          double phi   = atan2(pos[1], pos[0]);

          pos[0] = phi;
          pos[1] = theta - 0.5 * M_PI;
        }
      else
        {
          maxrad = SphP[igas].MaxDelaunayRadius;

          if(pos[0] + maxrad < -0.5 * self->width || pos[0] - maxrad > +0.5 * self->width || pos[1] + maxrad < -0.5 * self->height ||
             pos[1] - maxrad > +0.5 * self->height || pos[2] + maxrad < -0.5 * self->depth || pos[2] - maxrad > +0.5 * self->depth)
            continue;
        }

      int xmin = (int)floor((pos[0] - maxrad + 0.5 * self->width) / dx);
      int xmax = (int)ceil((pos[0] + maxrad + 0.5 * self->width) / dx);
      int ymin = (int)floor((pos[1] - maxrad + 0.5 * self->height) / dy);
      int ymax = (int)ceil((pos[1] + maxrad + 0.5 * self->height) / dy);

      xmin = imax(xmin, 0);
      ymin = imax(ymin, 0);
      xmax = imin(xmax, self->nx - 1);
      ymax = imin(ymax, self->ny - 1);

      int nx = xmax - xmin + 1;
      int ny = ymax - ymin + 1;

      if(nx <= 0 || ny <= 0)
        continue;

      double *lmin = (double *)mymalloc("lmin", nx * ny * sizeof(double));
      double *lmax = (double *)mymalloc("lmax", nx * ny * sizeof(double));

      int ix, iy;
      if(self->allsky)
        {
          for(ix = 0; ix < nx; ix++)
            for(iy = 0; iy < ny; iy++)
              {
                int idx   = iy * nx + ix;
                lmin[idx] = 0;
                lmax[idx] = +self->depth;
              }
        }
      else
        {
          for(ix = 0; ix < nx; ix++)
            for(iy = 0; iy < ny; iy++)
              {
                int idx   = iy * nx + ix;
                lmin[idx] = -0.5 * self->depth;
                lmax[idx] = +0.5 * self->depth;
              }
        }

      double line_dir[3];
      double line_base[3];

      int k;
      for(k = 0; k < 3; k++)
        {
          line_dir[k]  = 0;
          line_base[k] = 0;
        }

      if(self->allsky)
        {
          line_base[0] = xoffset;
        }
      else
        {
          line_dir[0] = 0;
          line_dir[1] = 0;
          line_dir[2] = 1;
        }

      double *pcell = P[igas].Pos;
      int q         = SphP[igas].first_connection;
      while(q >= 0)
        {
          int dp       = DC[q].dp_index;
          int particle = Mesh.DP[dp].index;
          if(particle < 0)
            {
              q = DC[q].next;
              continue;
            }

          double pother[3];
          pother[0] = Mesh.DP[dp].x;
          pother[1] = Mesh.DP[dp].y;
          pother[2] = Mesh.DP[dp].z;

          int vf = DC[q].vf_index;

          double face_center[3];
          face_center[0] = Mesh.VF[vf].cx;
          face_center[1] = Mesh.VF[vf].cy;
          face_center[2] = Mesh.VF[vf].cz;

          auriga_movie_projection_vector_shift(face_center, face_center, self->center);
          auriga_movie_projection_vector_rotate(face_center, self->rotation);

          double normal[3];
          auriga_movie_projection_vector_shift(normal, pother, pcell);
          auriga_movie_projection_vector_rotate(normal, self->rotation);
          auriga_movie_projection_vector_normalize(normal);

          double ln = 0;
          if(!self->allsky)
            ln = auriga_movie_projection_dot_product(normal, line_dir);

          for(ix = xmin; ix <= xmax; ix++)
            {
              if(ix < 0 || ix >= self->nx)
                continue;

              double cx = -0.5 * self->width + (ix + 0.5) * dx;
              for(iy = ymin; iy <= ymax; iy++)
                {
                  if(iy < 0 || iy >= self->ny)
                    continue;

                  int idx = (iy - ymin) * nx + (ix - xmin);
                  if(lmax[idx] <= lmin[idx])
                    continue;

                  double cy = -0.5 * self->height + (iy + 0.5) * dy;

                  if(self->allsky)
                    {
                      // cx += 0.5*self->width;
                      cy += 0.5 * self->height;

                      /* phi = cx, theta = cy */
                      double nx = sin(cy) * cos(cx);
                      double ny = sin(cy) * sin(cx);
                      double nz = cos(cy);

                      line_dir[0] = nx;
                      line_dir[1] = ny;
                      line_dir[2] = nz;

                      ln = auriga_movie_projection_dot_product(normal, line_dir);
                    }
                  else
                    {
                      /* the line is (cx,cy,0)b + a * (0,0,1) */
                      line_base[0] = cx;
                      line_base[1] = cy;
                      line_base[2] = 0;
                    }

                  double diff[3];
                  auriga_movie_projection_vector_shift(diff, face_center, line_base);
                  double l = auriga_movie_projection_dot_product(diff, normal);

                  if(ln != 0)
                    {
                      l /= ln;

                      if(ln > 0 && l < lmax[idx])
                        lmax[idx] = l;
                      if(ln < 0 && l > lmin[idx])
                        lmin[idx] = l;
                    }
                  else
                    {
                      /* the whole line is outside the cell */
                      if(l < 0)
                        lmax[idx] = lmin[idx];
                    }
                }
            }

          if(q == SphP[igas].last_connection)
            break;

          q = DC[q].next;
        }

      for(ix = xmin; ix <= xmax; ix++)
        {
          if(ix < 0 || ix >= self->nx)
            continue;
          for(iy = ymin; iy <= ymax; iy++)
            {
              if(iy < 0 || iy >= self->ny)
                continue;

              int idx = (iy - ymin) * nx + (ix - xmin);
              if(lmax[idx] > lmin[idx])
                {
                  double l  = lmax[idx] - lmin[idx];
                  int pixel = iy * self->nx + ix;

                  double rho = SphP[igas].Density;
                  self->projection[aum_gas_rho][pixel] += l * rho;

                  double v[3];
                  auriga_movie_projection_vector_shift_non_periodic(v, P[igas].Vel, self->cmvel);
                  auriga_movie_projection_vector_rotate(v, self->rotation);
                  self->projection[aum_gas_vel][pixel] += l * rho * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
                  self->projection[aum_gas_velx][pixel] += l * rho * v[0];
                  self->projection[aum_gas_vely][pixel] += l * rho * v[1];
                  self->projection[aum_gas_velz][pixel] += l * rho * v[2];

#ifdef MHD
                  double B[3];
                  for(k = 0; k < 3; k++)
                    B[k] = SphP[igas].B[k];
                  auriga_movie_projection_vector_rotate(B, self->rotation);
                  self->projection[aum_gas_bfld][pixel] += l * (B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
                  self->projection[aum_gas_bx][pixel] += l * B[0];
                  self->projection[aum_gas_by][pixel] += l * B[1];
                  self->projection[aum_gas_bz][pixel] += l * B[2];
#endif

                  double meanweight = 4.0 / (1. + 3. * HYDROGEN_MASSFRAC + 4. * HYDROGEN_MASSFRAC * SphP[igas].Ne) * PROTONMASS;
                  double temp_in_K =
                      GAMMA_MINUS1 * SphP[igas].Utherm / BOLTZMANN * All.UnitEnergy_in_cgs / All.UnitMass_in_g * meanweight;
                  self->projection[aum_gas_temp][pixel] += l * rho * temp_in_K;

                  self->projection[aum_gas_pres][pixel] += l * SphP[igas].Pressure;

                  self->projection[aum_gas_metallicity][pixel] += l * rho * SphP[igas].Metallicity;

                  len[pixel] += l;
                }
            }
        }

      myfree(lmax);
      myfree(lmin);

      if(ThisTask == 0)
        {
          double current = second();
          double seconds = timediff(previous, current);
          previous       = current;

          timetotal += seconds;
          timetonextupdate += seconds;

          if(timetonextupdate > 60.)
            {
              printf("AURIGA MOVIE: Done with %d of %d cells after %gs on Task 0.\n", igas + 1, NumGas, timetotal);
              timetonextupdate -= 60.;
            }
        }
    }

  if(ThisTask == 0)
    MPI_Reduce(MPI_IN_PLACE, len, self->nx * self->ny, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  else
    MPI_Reduce(len, NULL, self->nx * self->ny, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      double lmin = +MAX_REAL_NUMBER;
      double lmax = 0;
      double lav  = 0;
      int k;
      for(k = 0; k < self->nx * self->ny; k++)
        {
          if(len[k] > lmax)
            lmax = len[k];
          if(len[k] < lmin)
            lmin = len[k];
          lav += len[k];
        }
      printf("AURIGA MOVIE: lmin=%g, lmax=%g, lav=%g, depth=%g\n", lmin, lmax, lav / (self->nx * self->ny), self->depth);
    }

  myfree(len);
}

void auriga_movie_projection_gas(struct auriga_movie_projection *self)
{
  double dx = self->width / self->nx;
  double dy = self->height / self->ny;

  int MaxNRay   = self->nx * self->ny;
  ray_data *Ray = (ray_data *)mymalloc("Rays", MaxNRay * sizeof(ray_data));

  int Nray = 0;
  int k, i;

  if(self->allsky)
    {
      double xoffset = -8e-3 / All.Time * All.HubbleParam;

      double pos[3];
      for(k = 0; k < 3; k++)
        pos[k] = self->center[k] + xoffset * self->rotation[0][k];

      mpi_printf("AURIGA MOVIE: center for allsky projection: %g,%g,%g\n", pos[0], pos[1], pos[2]);

      peanokey key = position_to_peanokey(pos);
      int no       = peanokey_to_topnode(key);
      if(DomainTask[no] == ThisTask)
        {
          int itheta, iphi;
          for(iphi = 0; iphi < self->nx; iphi++)
            {
              double phi = (iphi + 0.5) / self->nx * self->width - M_PI;
              for(itheta = 0; itheta < self->ny; itheta++)
                {
                  double theta = (itheta + 0.5) / self->ny * self->height;
                  double nx    = sin(theta) * cos(phi);
                  double ny    = sin(theta) * sin(phi);
                  double nz    = cos(theta);

                  for(k = 0; k < 3; k++)
                    {
                      Ray[Nray].pos[k] = pos[k];
                      Ray[Nray].dir[k] = nx * self->rotation[0][k] + ny * self->rotation[1][k] + nz * self->rotation[2][k];
                    }

                  Ray[Nray].len        = 0;
                  Ray[Nray].target_len = self->depth;
                  Ray[Nray].pixel      = itheta * self->nx + iphi;
                  Ray[Nray].index      = -1;
                  Ray[Nray].prev       = -1;
                  Ray[Nray].task       = ThisTask;

                  Nray++;
                }
            }
        }
    }
  else
    {
      int ix, iy;
      double pz = -0.5 * self->depth;
      for(ix = 0; ix < self->nx; ix++)
        {
          double px = (ix + 0.5 - self->nx / 2) * dx;
          for(iy = 0; iy < self->ny; iy++)
            {
              double py = (iy + 0.5 - self->ny / 2) * dy;
              for(k = 0; k < 3; k++)
                {
                  Ray[Nray].pos[k] =
                      self->center[k] + px * self->rotation[0][k] + py * self->rotation[1][k] + pz * self->rotation[2][k];
                  Ray[Nray].dir[k] = self->rotation[2][k];
                }

              peanokey key   = position_to_peanokey(Ray[Nray].pos);
              int no         = peanokey_to_topnode(key);
              Ray[Nray].task = DomainTask[no];
              if(Ray[Nray].task != ThisTask)
                continue;

              Ray[Nray].len        = 0;
              Ray[Nray].target_len = self->depth;
              Ray[Nray].pixel      = iy * self->nx + ix;
              Ray[Nray].index      = -1;
              Ray[Nray].prev       = -1;

              Nray++;
            }
        }
    }

  int NrayAll;
  MPI_Reduce(&Nray, &NrayAll, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  mpi_printf("AURIGA MOVIE: Created %d of %d rays.\n", NrayAll, self->nx * self->ny);

  mesh_search_data *searchdata = (mesh_search_data *)mymalloc("searchdata", Nray * sizeof(mesh_search_data));
  for(i = 0; i < Nray; i++)
    for(k = 0; k < 3; k++)
      searchdata[i].Pos[k] = Ray[i].pos[k];

  find_nearest_meshpoint_global(searchdata, Nray, 0, 1);

  mpi_distribute_items_from_search(searchdata, Ray, &Nray, &MaxNRay, sizeof(ray_data), TAG_DENS_B, offsetof(ray_data, task),
                                   offsetof(ray_data, index));
  myfree(searchdata);

  mpi_printf("AURIGA MOVIE: Distributed rays do tasks.\n");

  /* advance all rays until they have reached their final point */
  int rays_left, iter = 0;
  do
    {
      int i, j, nleft = 0;

      for(i = 0; i < Nray; i++)
        {
          if(Ray[i].len >= Ray[i].target_len)
            continue;

          /* this is the index of the cell in which we currently are */
          int sph_idx = Ray[i].index;
          double l;
          int next_edge = find_next_voronoi_cell(&Mesh, sph_idx, Ray[i].pos, Ray[i].dir, Ray[i].prev, &l);

          if(Ray[i].len + l > Ray[i].target_len)
            {
              l           = Ray[i].target_len - Ray[i].len;
              Ray[i].len  = Ray[i].target_len;
              Ray[i].prev = -1;
            }
          else
            {
              Ray[i].task  = DC[next_edge].task;
              Ray[i].prev  = Ray[i].index;
              Ray[i].index = DC[next_edge].index;
              Ray[i].len += l;
              ++nleft;
            }

          double rho = SphP[sph_idx].Density;
          self->projection[aum_gas_rho][Ray[i].pixel] += l * rho;

          double vx = P[sph_idx].Vel[0] - self->cmvel[0];
          double vy = P[sph_idx].Vel[1] - self->cmvel[1];
          double vz = P[sph_idx].Vel[2] - self->cmvel[2];
          self->projection[aum_gas_vel][Ray[i].pixel] += l * rho * (vx * vx + vy * vy + vz * vz);
          self->projection[aum_gas_velx][Ray[i].pixel] += l * rho * vx;
          self->projection[aum_gas_vely][Ray[i].pixel] += l * rho * vy;
          self->projection[aum_gas_velz][Ray[i].pixel] += l * rho * vz;

#ifdef MHD
          double bx = SphP[sph_idx].B[0];
          double by = SphP[sph_idx].B[1];
          double bz = SphP[sph_idx].B[2];
          self->projection[aum_gas_bfld][Ray[i].pixel] += l * (bx * bx + by * by + bz * bz);
          self->projection[aum_gas_bx][Ray[i].pixel] += l * bx;
          self->projection[aum_gas_by][Ray[i].pixel] += l * by;
          self->projection[aum_gas_bz][Ray[i].pixel] += l * bz;
#endif

          double meanweight = 4.0 / (1. + 3. * HYDROGEN_MASSFRAC + 4. * HYDROGEN_MASSFRAC * SphP[sph_idx].Ne) * PROTONMASS;
          double temp_in_K  = GAMMA_MINUS1 * SphP[sph_idx].Utherm / BOLTZMANN * All.UnitEnergy_in_cgs / All.UnitMass_in_g * meanweight;
          self->projection[aum_gas_temp][Ray[i].pixel] += l * rho * temp_in_K;

          self->projection[aum_gas_pres][Ray[i].pixel] += l * SphP[sph_idx].Pressure;

          self->projection[aum_gas_metallicity][Ray[i].pixel] += l * rho * SphP[sph_idx].Metallicity;

          for(j = 0; j < 3; j++)
            Ray[i].pos[j] += l * Ray[i].dir[j];

          MyDouble ref[3];
          for(j = 0; j < 3; j++)
            ref[j] = (&Mesh.DP[DC[next_edge].dp_index].x)[j];
          periodic_wrap_point_MyDouble(Ray[i].pos, ref);
        }

      for(i = 0; i < Nray; i++)
        {
          if(Ray[i].task != ThisTask)
            Ray[i].prev = -1;
        }
      mpi_distribute_items_to_tasks(Ray, offsetof(ray_data, task), &Nray, &MaxNRay, sizeof(*Ray), TAG_DENS_A);

      MPI_Allreduce(&nleft, &rays_left, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      if(iter % 50 == 0)
        mpi_printf("AURIGA_MOVIE: iteration %d rays left: %5d\n", iter, rays_left);
      iter++;
    }
  while(rays_left);

  myfree(Ray);
}

void auriga_movie_projection_dark_matter(struct auriga_movie_projection *self)
{
  double dx = self->width / self->nx;
  double dy = self->height / self->ny;

  double xoffset = 0;
  if(self->allsky)
    xoffset = -8e-3 / All.Time * All.HubbleParam;

  int idm;
  for(idm = 0; idm < NumPart; idm++)
    {
      if(P[idm].Type == 0 || P[idm].Type > 3)
        continue;

      double px = NEAREST_X(P[idm].Pos[0] - self->center[0]);
      double py = NEAREST_Y(P[idm].Pos[1] - self->center[1]);
      double pz = NEAREST_Z(P[idm].Pos[2] - self->center[2]);

      double tx = px;
      double ty = py;
      double tz = pz;
      px        = tx * self->rotation[0][0] + ty * self->rotation[0][1] + tz * self->rotation[0][2];
      py        = tx * self->rotation[1][0] + ty * self->rotation[1][1] + tz * self->rotation[1][2];
      pz        = tx * self->rotation[2][0] + ty * self->rotation[2][1] + tz * self->rotation[2][2];

      px -= xoffset;

      double hsml;
      if(self->allsky)
        {
          double rad = sqrt(px * px + py * py + pz * pz);
          hsml       = asin(P[idm].Auriga_Movie_Hsml / rad);
          hsml       = fmax(hsml, 1.001 * fmax(dx, dy) * 0.5);

          double theta = acos(pz / rad);
          double phi   = atan2(py, px);

          px = phi;
          py = theta - 0.5 * M_PI;

          if(rad > self->depth)
            continue;
        }
      else
        {
          hsml = fmax(P[idm].Auriga_Movie_Hsml, 1.001 * fmax(dx, dy) * 0.5);
          if(px + hsml < -0.5 * self->width || px - hsml > +0.5 * self->width || py + hsml < -0.5 * self->height ||
             py - hsml > +0.5 * self->height || pz < -0.5 * self->depth || pz > +0.5 * self->depth)
            continue;
        }

      int xmin = (int)floor((px - hsml + 0.5 * self->width) / dx);
      int xmax = (int)ceil((px + hsml + 0.5 * self->width) / dx);
      int ymin = (int)floor((py - hsml + 0.5 * self->height) / dy);
      int ymax = (int)ceil((py + hsml + 0.5 * self->height) / dy);

      double sum = 0;

      int ix, iy;
      for(ix = xmin; ix <= xmax; ix++)
        {
          double cx = -0.5 * self->width + (ix + 0.5) * dx;
          for(iy = ymin; iy <= ymax; iy++)
            {
              double cy = -0.5 * self->height + (iy + 0.5) * dy;
              double r2 = (cx - px) * (cx - px) + (cy - py) * (cy - py);
              sum += hsml * _getkernel(hsml, r2);
            }
        }

      for(ix = xmin; ix <= xmax; ix++)
        {
          double cx = -0.5 * self->width + (ix + 0.5) * dx;
          for(iy = ymin; iy <= ymax; iy++)
            {
              if(ix < 0 || iy < 0 || ix >= self->nx || iy >= self->ny)
                continue;
              double cy = -0.5 * self->height + (iy + 0.5) * dy;
              double r2 = (cx - px) * (cx - px) + (cy - py) * (cy - py);
              double wk = hsml * _getkernel(hsml, r2) / sum;

              if(wk > 0.)
                {
                  int idx = iy * self->nx + ix;
                  self->projection[aum_dm_rho][idx] += wk * P[idm].Mass;
                  self->projection[aum_dm_rho2][idx] += wk * P[idm].Mass * P[idm].Auriga_Movie_Density;
                }
            }
        }
    }
}

void auriga_movie_projection_stars(struct auriga_movie_projection *self)
{
  double dx = self->width / self->nx;
  double dy = self->height / self->ny;

  double xoffset = 0;
  if(self->allsky)
    xoffset = -8e-3 / All.Time * All.HubbleParam;

  int istar;
  for(istar = 0; istar < NumPart; istar++)
    {
      if(P[istar].Type != 4 || StarP[P[istar].AuxDataID].BirthTime < 0 || P[istar].Mass > 2. * All.TargetGasMass)
        continue;

      double px = NEAREST_X(P[istar].Pos[0] - self->center[0]);
      double py = NEAREST_Y(P[istar].Pos[1] - self->center[1]);
      double pz = NEAREST_Z(P[istar].Pos[2] - self->center[2]);

      double tx = px;
      double ty = py;
      double tz = pz;
      px        = tx * self->rotation[0][0] + ty * self->rotation[0][1] + tz * self->rotation[0][2];
      py        = tx * self->rotation[1][0] + ty * self->rotation[1][1] + tz * self->rotation[1][2];
      pz        = tx * self->rotation[2][0] + ty * self->rotation[2][1] + tz * self->rotation[2][2];

      px -= xoffset;

      double hsml;
      if(self->allsky)
        {
          double rad = sqrt(px * px + py * py + pz * pz);
          hsml       = asin(P[istar].Auriga_Movie_Hsml / rad);
          hsml       = fmin(hsml, fmax(4. * fmax(dx, dy), 2. * M_PI / 720.));
          hsml       = fmax(hsml, 1.001 * fmax(dx, dy) * 0.5);

          double theta = acos(pz / rad);
          double phi   = atan2(py, px);

          px = phi;
          py = theta - 0.5 * M_PI;

          if(rad > self->depth)
            continue;
        }
      else
        {
          hsml = fmin(P[istar].Auriga_Movie_Hsml, fmax(4. * fmax(dx, dy), 2e-5 / All.Time * All.HubbleParam));
          hsml = fmax(hsml, 1.001 * fmax(dx, dy) * 0.5);
          if(px + hsml < -0.5 * self->width || px - hsml > +0.5 * self->width || py + hsml < -0.5 * self->height ||
             py - hsml > +0.5 * self->height || pz < -0.5 * self->depth || pz > +0.5 * self->depth)
            continue;
        }

      int xmin = floor((px - hsml + 0.5 * self->width) / dx);
      int xmax = ceil((px + hsml + 0.5 * self->width) / dx);
      int ymin = floor((py - hsml + 0.5 * self->height) / dy);
      int ymax = ceil((py + hsml + 0.5 * self->height) / dy);

      double sum = 0;

      int ix, iy;
      for(ix = xmin; ix <= xmax; ix++)
        {
          double cx = -0.5 * self->width + (ix + 0.5) * dx;
          for(iy = ymin; iy <= ymax; iy++)
            {
              double cy = -0.5 * self->height + (iy + 0.5) * dy;
              double r2 = (cx - px) * (cx - px) + (cy - py) * (cy - py);
              sum += hsml * _getkernel(hsml, r2);
            }
        }

      stellar_photometrics st_photo;
      assign_stellar_photometrics(istar, &st_photo);
      double M_U = pow(10., -0.4 * st_photo.Magnitude_U);
      double M_B = pow(10., -0.4 * st_photo.Magnitude_B);
      double M_K = pow(10., -0.4 * st_photo.Magnitude_K);

      xmin = imax(xmin, 0);
      xmax = imin(xmax, self->nx - 1);
      ymin = imax(ymin, 0);
      ymax = imin(ymax, self->ny - 1);

      for(ix = xmin; ix <= xmax; ix++)
        {
          double cx = -0.5 * self->width + (ix + 0.5) * dx;
          for(iy = ymin; iy <= ymax; iy++)
            {
              // if(ix < 0 || iy < 0 || ix >= self->nx || iy >= self->ny)
              //  continue;

              double cy = -0.5 * self->height + (iy + 0.5) * dy;
              double r2 = (cx - px) * (cx - px) + (cy - py) * (cy - py);
              double wk = hsml * _getkernel(hsml, r2) / sum;

              if(wk > 0.)
                {
                  self->projection[aum_stars_rho][iy * self->nx + ix] += wk * P[istar].Mass;
                  self->projection[aum_stars_metallicity][iy * self->nx + ix] +=
                      wk * P[istar].Mass * StarP[P[istar].AuxDataID].Metallicity;

                  self->projection[aum_stars_u][iy * self->nx + ix] += wk * M_U;
                  self->projection[aum_stars_g][iy * self->nx + ix] += wk * M_B;
                  self->projection[aum_stars_r][iy * self->nx + ix] += wk * M_K;
                }
            }
        }
    }
}

void auriga_movie_projection_save(struct auriga_movie_projection *self, int outNum, char *movieDir, char *prefix)
{
  if(ThisTask != 0)
    return;

  mkdir(movieDir, MKDIR_MODE);

  char buf[1000];
  sprintf(buf, "%s/%s", movieDir, prefix);
  mkdir(buf, MKDIR_MODE);

  /* write general information */
  sprintf(buf, "%s/%s/projection_%05d.txt", movieDir, prefix, outNum);

  FILE *fp = fopen(buf, "w");
  fprintf(fp, "time: %g\n", All.Time);
  fprintf(fp, "redshift: %g\n", 1. / All.Time - 1.);
  fprintf(fp, "width: %g\n", self->width);
  fprintf(fp, "height: %g\n", self->height);
  fprintf(fp, "nx: %d\n", self->nx);
  fprintf(fp, "ny: %d\n", self->ny);
  fprintf(fp, "center: %g %g %g\n", self->center[0], self->center[1], self->center[2]);
  fprintf(fp, "rotation: %g %g %g %g %g %g %g %g %g\n", self->rotation[0][0], self->rotation[0][1], self->rotation[0][2],
          self->rotation[1][0], self->rotation[1][1], self->rotation[1][2], self->rotation[2][0], self->rotation[2][1],
          self->rotation[2][2]);
  fclose(fp);

  /* write projections */
  int k;
  for(k = 0; k < aum_count; k++)
    {
      printf("Writing projection %s of type %s.\n", aum_projection_names[k], prefix);
      sprintf(buf, "%s/%s/%s_%05d.dat", movieDir, prefix, aum_projection_names[k], outNum);

      fp = fopen(buf, "w");
      fwrite(self->projection[k], sizeof(double), self->nx * self->ny, fp);
      fclose(fp);
    }
}

void auriga_movie_projection_free(struct auriga_movie_projection *self)
{
  int k;
  for(k = aum_count - 1; k >= 0; k--)
    myfree(self->projection[k]);
}

#endif
