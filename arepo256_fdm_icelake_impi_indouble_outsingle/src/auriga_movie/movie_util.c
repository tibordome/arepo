/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/movie_auriga/movie_util.c
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

#include <hdf5.h>

#include "../allvars.h"
#include "../proto.h"
#include "movie.h"

#ifdef AURIGA_MOVIE

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>

int auriga_movie_get_center_guess_groupcatalogue(double center[3])
{
  int halo = 0;

  /* we want to know the center of halo 0 */
  if(ThisTask == 0)
    {
      char fname[MAXLEN_PATH];

      if(RestartSnapNum < 1000)
        file_path_sprintf(fname, "%s/groups_%03d/fof_subhalo_tab_%03d.0.hdf5", All.OutputDir, RestartSnapNum, RestartSnapNum);
      else
        file_path_sprintf(fname, "%s/groups_%03d/fof_subhalo_tab_%d.0.hdf5", All.OutputDir, RestartSnapNum, RestartSnapNum);

      printf("AURIGA MOVIE: Trying to read group catalogue file `%s`.\n", fname);
      hid_t hdf5_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
      if(hdf5_file >= 0)
        {
          hid_t hdf5_grp = my_H5Gopen(hdf5_file, "/Subhalo");

          hid_t hdf5_dataset = my_H5Dopen_if_existing(hdf5_grp, "SubhaloPos");
          if(hdf5_dataset >= 0)
            {
              hid_t hdf5_space = H5Dget_space(hdf5_dataset);

              hsize_t data_offset[2];
              data_offset[0] = 0;
              data_offset[1] = 0;
              hsize_t data_count[2];
              data_count[0] = 1;
              data_count[1] = 3;
              H5Sselect_hyperslab(hdf5_space, H5S_SELECT_SET, data_offset, NULL, data_count, NULL);

              hsize_t memory_dim[2];
              memory_dim[0]  = 1;
              memory_dim[1]  = 3;
              hid_t memspace = H5Screate_simple(2, memory_dim, NULL);

              my_H5Dread(hdf5_dataset, H5T_NATIVE_DOUBLE, memspace, hdf5_space, H5P_DEFAULT, center, "SubhaloPos");
              printf("Found subhalo 0 at coordinates %g,%g,%g.\n", center[0], center[1], center[2]);

              domain_displacePosition(center, DISPLACE_POSITION_FORWARD);

              my_H5Dclose(hdf5_dataset, "SubhaloPos");
              my_H5Gclose(hdf5_grp, "/Subhalo");
              my_H5Fclose(hdf5_file, fname);
              halo = 1;
            }
          else
            {
              printf("Did not find dataset SubhaloPos.\n");
            }
        }
      else
        printf("Could not open file %s.\n", fname);
    }

  MPI_Bcast(center, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&halo, 1, MPI_INT, 0, MPI_COMM_WORLD);
  return halo;
}

void auriga_movie_get_center_guess_onthefly(double center[3])
{
  All.Auriga_Movie_Halo_Initialized = 1;
  /*
  if(!All.Auriga_Movie_Halo_Initialized)
    {
      int k;
      for(k = 0; k < 3; k++)
        center[k] = All.Auriga_Movie_Halo_Center[k];
      All.Auriga_Movie_Halo_Initialized = 1;
    }
  else
  */
  {
    double cm[3];
    double mass;

    int k;
    for(k = 0; k < 3; k++)
      cm[k] = 0;
    mass = 0;

    int i;
    for(i = 0; i < NumPart; i++)
      {
        if(P[i].Auriga_Movie_Center_Marker)
          {
            double dx, dy, dz;

            if(P[i].Type == 0)
              {
                dx = NEAREST_X(SphP[i].Center[0] - 0.5 * All.BoxSize);
                dy = NEAREST_Y(SphP[i].Center[1] - 0.5 * All.BoxSize);
                dz = NEAREST_Z(SphP[i].Center[2] - 0.5 * All.BoxSize);
              }
            else
              {
                dx = NEAREST_X(P[i].Pos[0] - 0.5 * All.BoxSize);
                dy = NEAREST_Y(P[i].Pos[1] - 0.5 * All.BoxSize);
                dz = NEAREST_Z(P[i].Pos[2] - 0.5 * All.BoxSize);
              }

            cm[0] += dx * P[i].Mass;
            cm[1] += dy * P[i].Mass;
            cm[2] += dz * P[i].Mass;
            mass += P[i].Mass;
          }
      }

    double allcm[3], allmass;
    MPI_Allreduce(&mass, &allmass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(cm, allcm, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if(allmass == 0)
      terminate("Fail.");

    for(k = 0; k < 3; k++)
      {
        center[k] = 0.5 * All.BoxSize + allcm[k] / allmass;
      }
  }
}

int auriga_movie_compute_center_and_rotation(double center[3], double rotation[3][3], double cmvel[3], int niter)
{
  double radius  = All.Auriga_Movie_CenterRadius / All.HubbleParam * All.Time;  // convert from physical to internal code units
  double radius2 = radius * radius;

  /* recenter, radius 20kpc physical */
  int iter;
  for(iter = 0; iter < niter; iter++)
    {
      double cm[3], mass;
      int j;

      mass = 0;
      for(j = 0; j < 3; j++)
        cm[j] = 0;

      int i;
      for(i = 0; i < NumPart; i++)
        {
          double dx, dy, dz;

          if(P[i].Type == 0)
            {
              dx = NEAREST_X(SphP[i].Center[0] - center[0]);
              dy = NEAREST_Y(SphP[i].Center[1] - center[1]);
              dz = NEAREST_Z(SphP[i].Center[2] - center[2]);
            }
          else
            {
              dx = NEAREST_X(P[i].Pos[0] - center[0]);
              dy = NEAREST_Y(P[i].Pos[1] - center[1]);
              dz = NEAREST_Z(P[i].Pos[2] - center[2]);
            }
          double r2 = dx * dx + dy * dy + dz * dz;

          if(r2 < radius2)
            {
              mass += P[i].Mass;
              cm[0] += dx * P[i].Mass;
              cm[1] += dy * P[i].Mass;
              cm[2] += dz * P[i].Mass;
            }
        }

      double allcm[3], allmass;
      MPI_Allreduce(&mass, &allmass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      if(allmass == 0)
        {
          /* radius is too small, stop with current center */
          mpi_printf("AURIGA MOVIE: iter=%d, no mass in sphere of radius %g, using current center estimate of %g,%g,%g.\n", iter,
                     sqrt(radius2), center[0], center[1], center[2]);
          MPI_Bcast(center, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
          break;
        }

      MPI_Reduce(cm, allcm, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

      if(ThisTask == 0)
        {
          for(j = 0; j < 3; j++)
            center[j] += allcm[j] / allmass;

          mpi_printf("AURIGA MOVIE: iter=%d, Mass in CenterRadius: %g, delta=%g,%g,%g, new center: %g,%g,%g\n", iter, allmass,
                     allcm[0] / allmass, allcm[1] / allmass, allcm[2] / allmass, center[0], center[1], center[2]);
        }

      MPI_Bcast(center, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      radius2 *= 0.7 * 0.7;
    }

  radius2 = radius * radius;

  double mass = 0;
  double massvel[3];
  double L[3], tensor[3][3];
  int k, l;
  for(k = 0; k < 3; k++)
    {
      massvel[k] = 0;
      L[k]       = 0;
      for(l = 0; l < 3; l++)
        tensor[k][l] = 0;
    }

  /* compute angular momentum and moment tensor */
  int i;
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type != 4)
        continue;

      double dx, dy, dz;

      if(P[i].Type == 0)
        {
          dx = NEAREST_X(SphP[i].Center[0] - center[0]);
          dy = NEAREST_Y(SphP[i].Center[1] - center[1]);
          dz = NEAREST_Z(SphP[i].Center[2] - center[2]);
        }
      else
        {
          dx = NEAREST_X(P[i].Pos[0] - center[0]);
          dy = NEAREST_Y(P[i].Pos[1] - center[1]);
          dz = NEAREST_Z(P[i].Pos[2] - center[2]);
        }
      double r2 = dx * dx + dy * dy + dz * dz;

      if(r2 < radius2)
        {
          L[0] += P[i].Mass * (dy * P[i].Vel[2] - dz * P[i].Vel[1]);
          L[1] += P[i].Mass * (dz * P[i].Vel[0] - dx * P[i].Vel[2]);
          L[2] += P[i].Mass * (dx * P[i].Vel[1] - dy * P[i].Vel[0]);

          for(k = 0; k < 3; k++)
            massvel[k] += P[i].Mass * P[i].Vel[k];
          mass += P[i].Mass;

          tensor[0][0] += P[i].Mass * (dy * dy + dz * dz);
          tensor[1][1] += P[i].Mass * (dx * dx + dz * dz);
          tensor[2][2] += P[i].Mass * (dx * dx + dy * dy);

          tensor[0][1] += -P[i].Mass * dx * dy;
          tensor[1][0] += -P[i].Mass * dx * dy;

          tensor[0][2] += -P[i].Mass * dx * dz;
          tensor[2][0] += -P[i].Mass * dx * dz;

          tensor[1][2] += -P[i].Mass * dy * dz;
          tensor[2][1] += -P[i].Mass * dy * dz;
        }

      P[i].Auriga_Movie_Center_Marker = 0;
      if(r2 < 0.01 * radius2)
        {
          P[i].Auriga_Movie_Center_Marker = 1;
        }
    }

  double massall;
  MPI_Allreduce(&mass, &massall, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  /* check if there are any stars for meaningful estimates */
  if(massall == 0)
    {
      int ix, iy;
      for(ix = 0; ix < 3; ix++)
        for(iy = 0; iy < 3; iy++)
          {
            if(ix == iy)
              rotation[ix][iy] = 1.;
            else
              rotation[ix][iy] = 0.;
          }

      int k;
      for(k = 0; k < 3; k++)
        cmvel[k] = 0;

      return 2;
    }

  double Ltot[3];
  MPI_Reduce(L, Ltot, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  double tensorAll[3][3];
  MPI_Reduce(tensor, tensorAll, 9, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      double Lsum = sqrt(Ltot[0] * Ltot[0] + Ltot[1] * Ltot[1] + Ltot[2] * Ltot[2]);
      for(k = 0; k < 3; k++)
        Ltot[k] /= Lsum;

      printf("Ltot: %g, %g, %g, Mass: %g\n", Ltot[0], Ltot[1], Ltot[2], massall);

      gsl_matrix *m = gsl_matrix_alloc(3, 3);
      for(k = 0; k < 3; k++)
        for(l = 0; l < 3; l++)
          gsl_matrix_set(m, k, l, tensorAll[k][l]);

      gsl_vector *eval             = gsl_vector_alloc(3);
      gsl_matrix *evec             = gsl_matrix_alloc(3, 3);
      gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(3);

      gsl_eigen_symmv(m, eval, evec, w);
      gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_DESC);

      int idx        = 0;
      double dot_old = 0;
      for(l = 0; l < 3; l++)
        dot_old += gsl_matrix_get(evec, l, idx) * Ltot[l];
      dot_old = fabs(dot_old);

      for(k = 1; k < 3; k++)
        {
          double dot_new = 0;
          for(l = 0; l < 3; l++)
            dot_new += gsl_matrix_get(evec, l, k) * Ltot[l];
          dot_new = fabs(dot_new);

          if(dot_new > dot_old)
            {
              idx     = k;
              dot_old = dot_new;
            }
        }

      double prefac = 1.;
      if(dot_old < 0.)
        prefac = -1.;

      for(l = 0; l < 3; l++)
        Ltot[l] = prefac * gsl_matrix_get(evec, l, idx);

      gsl_eigen_symmv_free(w);
      gsl_matrix_free(evec);
      gsl_vector_free(eval);
      gsl_matrix_free(m);

      for(k = 0; k < 3; k++)
        rotation[2][k] = Ltot[k];

      printf("Eigenvector along L (chose %d): %g, %g, %g\n", idx, rotation[2][0], rotation[2][1], rotation[2][2]);

      if(rotation[2][0] != 0 || rotation[2][1] != 0)
        {
          rotation[0][0] = -rotation[2][1];
          rotation[0][1] = +rotation[2][0];
          rotation[0][2] = 0;
        }
      else
        {
          rotation[0][0] = 1.;
          rotation[0][1] = 0.;
          rotation[0][2] = 0.;
        }

      double vsum = sqrt(rotation[0][0] * rotation[0][0] + rotation[0][1] * rotation[0][1] + rotation[0][2] * rotation[0][2]);
      for(k = 0; k < 3; k++)
        rotation[0][k] /= vsum;

      rotation[1][0] = rotation[2][1] * rotation[0][2] - rotation[2][2] * rotation[0][1];
      rotation[1][1] = rotation[2][2] * rotation[0][0] - rotation[2][0] * rotation[0][2];
      rotation[1][2] = rotation[2][0] * rotation[0][1] - rotation[2][1] * rotation[0][0];

      printf("AURIGA MOVIE Rotation matrix: %8g %8g %8g\n", rotation[0][0], rotation[1][0], rotation[2][0]);
      printf("AURIGA MOVIE Rotation matrix: %8g %8g %8g\n", rotation[0][1], rotation[1][1], rotation[2][1]);
      printf("AURIGA MOVIE Rotation matrix: %8g %8g %8g\n", rotation[0][2], rotation[1][2], rotation[2][2]);
    }

  MPI_Bcast(rotation, 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Reduce(massvel, cmvel, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(ThisTask == 0)
    for(k = 0; k < 3; k++)
      cmvel[k] /= massall;
  MPI_Bcast(cmvel, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  return 1;
}

void auriga_movie_get_center_highres(double center[3], double cmvel[3])
{
  int i;
  if(ThisTask == 0)
    {
      for(i = 0; i < NumPart; i++)
        {
          if(P[i].Type == 1)
            {
              int k;
              for(k = 0; k < 3; k++)
                center[k] = P[i].Pos[k];
              break;
            }
        }
    }
  MPI_Bcast(center, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  double dist[3];
  double mass;
  double massvel[3];

  int k;
  for(k = 0; k < 3; k++)
    {
      dist[k]    = 0;
      massvel[k] = 0;
    }
  mass = 0;

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 1)
        {
          dist[0] += NEAREST_X(P[i].Pos[0] - center[0]) * P[i].Mass;
          dist[1] += NEAREST_Y(P[i].Pos[1] - center[1]) * P[i].Mass;
          dist[2] += NEAREST_Z(P[i].Pos[2] - center[2]) * P[i].Mass;
          mass += P[i].Mass;
          for(k = 0; k < 3; k++)
            massvel[k] += P[i].Mass * P[i].Vel[k];
        }
    }

  double alldist[3], allmass, allmassvel[3];
  MPI_Reduce(&mass, &allmass, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(dist, alldist, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(massvel, allmassvel, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      for(k = 0; k < 3; k++)
        {
          center[k] += alldist[k] / allmass;
          cmvel[k] = allmassvel[k] / allmass;
        }
    }
  MPI_Bcast(center, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(cmvel, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  mpi_printf("AURIGA MOVIE: Center of the highres region: %g,%g,%g vel: %g,%g,%g\n", center[0], center[1], center[2], cmvel[0],
             cmvel[1], cmvel[2]);
}

#endif
