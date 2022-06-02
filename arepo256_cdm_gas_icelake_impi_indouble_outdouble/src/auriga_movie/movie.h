/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/movie_auriga/movie.h
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

#ifndef AURIGA_MOVIE_H
#define AURIGA_MOVIE_H

#ifdef AURIGA_MOVIE

#define MAX_MOVIE_OUTPUT_TIMES 8192

void auriga_movie_init(void);
void auriga_movie_check_output(void);
void auriga_movie_make(void);

/* utils */
int auriga_movie_get_center_guess_groupcatalogue(double center[3]);
void auriga_movie_get_center_guess_onthefly(double center[3]);
int auriga_movie_compute_center_and_rotation(double center[3], double rotation[3][3], double cmvel[3], int niter);
void auriga_movie_get_center_highres(double center[3], double cmvel[3]);

/* density */
void auriga_movie_calculate_smoothing_lenghts(void);

/* projection */
enum
{
  aum_gas_rho = 0,
  aum_gas_vel,
  aum_gas_velx,
  aum_gas_vely,
  aum_gas_velz,
  aum_gas_bfld,
  aum_gas_bx,
  aum_gas_by,
  aum_gas_bz,
  aum_gas_temp,
  aum_gas_pres,
  aum_gas_metallicity,
  aum_dm_rho,
  aum_dm_rho2,
  aum_stars_rho,
  aum_stars_u,
  aum_stars_g,
  aum_stars_r,
  aum_stars_metallicity,
  aum_count
};

extern char *aum_projection_names[aum_count];

struct auriga_movie_projection
{
  double center[3];
  double rotation[3][3];
  double cmvel[3];

  int nx, ny;
  double width, height, depth;
  int allsky;

  double *projection[aum_count];
};

void auriga_movie_projection_init(struct auriga_movie_projection *self, double center[3], double rotation[3][3], double cmvel[3],
                                  double width, double depth, int nx, int ny);
void auriga_movie_projection_project_particles(struct auriga_movie_projection *self);
void auriga_movie_projection_save(struct auriga_movie_projection *self, int outNum, char *movieDir, char *prefix);
void auriga_movie_projection_free(struct auriga_movie_projection *self);

#endif

#endif
