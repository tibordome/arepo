/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/SMUGGLE/sfr_proto.h
 * \date        03/2020
 * \author
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#ifndef FM_STAR_FORM_PROTO_H
#define FM_STAR_FORM_PROTO_H

void init_star_formation(void);

#ifdef SMUGGLE_OUTPUT_SF_PROBABILITY
double compute_sf_probability(int index);
#endif

#endif
