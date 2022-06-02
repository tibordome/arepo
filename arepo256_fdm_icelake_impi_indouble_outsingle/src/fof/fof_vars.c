/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/fof/fof_vars.c
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

#include <gsl/gsl_math.h>
#include <inttypes.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "../allvars.h"
#include "../domain.h"
#include "../proto.h"
#include "../subfind/subfind.h"
#include "fof.h"

/*! \file fof_vars.c
 *  \brief In this file we create instances for the global variables used by FOF, which are declared in fof.h
 */

#ifdef FOF

int Ngroups, NgroupsExt, MaxNgroups, TotNgroups, Nsubgroups, TotNsubgroups;
int Nids;
long long TotNids;

double LinkL = 0;

int fof_OldMaxPart;
int fof_OldMaxPartSph;
#if defined(GFM) || defined(SFR_MCS)
int fof_OldMaxPartStar;
#endif
#ifdef BLACK_HOLES
int fof_OldMaxPartBHs;
#endif
#ifdef DUST_LIVE
int fof_OldMaxPartDust;
#endif

unsigned char *flag_node_inside_linkinglength;

struct group_properties *Group;

#ifdef ADD_GROUP_PROPERTIES
int Ngroups_eff;
struct group_properties *GroupAll;
struct group_catalogue *GroupCat, *GroupCatLocal, *GroupCatSend;
#endif

struct fofdata_in *FoFDataIn, *FoFDataGet;

struct fofdata_out *FoFDataResult, *FoFDataOut;

struct fof_particle_list *FOF_PList;

struct fof_group_list *FOF_GList;

struct id_list *ID_list;

struct bit_flags *Flags;

struct fof_subfind_header catalogue_header;

#endif /* of FOF */
