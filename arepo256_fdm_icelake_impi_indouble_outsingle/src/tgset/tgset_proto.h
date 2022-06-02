/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/tgset/tgset_proto.h
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       Special settings for primordial runs
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

void tgset_begrun(void);
void tgset_dens_init(int i);
integertime tgset_snap_free_fall(integertime ti_curr);
void tgset_get_min_timestep(void);
void tgset_limit_timestep(integertime *pti_step, int *pbin, int binold);
void tgset_get_nh_max(void);
int tgset_jeans_ref(int mode, int i);
void tgset_get_image_limits(char **argv, int *xaxis, int *yaxis, int *zaxis, double *xmin, double *xmax, double *ymin, double *ymax,
                            double *zmin, double *zmax, int *weight_flag);
