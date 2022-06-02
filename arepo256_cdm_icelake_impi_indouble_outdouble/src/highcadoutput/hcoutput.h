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

#ifndef HCOUTPUT_H
#define HCOUTPUT_H

#ifdef HCOUTPUT

#define MAX_HCSNIPS_OUTPUT_TIMES 8192

void hcoutput_init(void);
void hcoutput_check_output(void);
void hcoutput_make(void);
void hcoutput_dump(void);

void hcoutput_get_center_guess_onthefly(double center[3]);
void hcoutput_compute_center(double center[3], double cmvel[3], int niter);

#endif

#endif
