/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/dust_vars.c
 * \date        MM/YYYY
 * \author      Ryan McKinnon
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"

#ifdef DUST_LIVE
int Ndust = 0;

struct dust_particle *DustParticle;

#ifdef DL_GRAIN_BINS
struct grain_size_distribution GSD;
#endif

#if(defined(DL_GRAIN_BINS) && (defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION))) || \
    defined(DL_DRAG_BACKREACTION)
int NforcesDust = 0;

struct shatter_data_p *PShatter;
#endif

#ifdef DL_RADIATION
struct dust_radiation_pressure DustRP_gra, DustRP_sil;
struct dust_rp_helpers DRT;
#endif

#endif
