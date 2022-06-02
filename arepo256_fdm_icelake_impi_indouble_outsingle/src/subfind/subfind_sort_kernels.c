/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/subfind/subfind_sort_kernels.c
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

#include <gsl/gsl_rng.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "../allvars.h"
#include "../domain.h"
#include "../fof/fof.h"
#include "../proto.h"
#include "subfind.h"

#ifdef SUBFIND

/*! \brief Comparison function for proc_assign_data objects.
 *
 *  Sorting kernel comparing element GrNr.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int subfind_compare_procassign_GrNr(const void *a, const void *b)
{
  if(((struct proc_assign_data *)a)->GrNr < ((struct proc_assign_data *)b)->GrNr)
    return -1;

  if(((struct proc_assign_data *)a)->GrNr > ((struct proc_assign_data *)b)->GrNr)
    return +1;

  return 0;
}

/*! \brief Comparison function for submp_data objects.
 *
 *  Sorting kernel comparing element (most important first):
 *  GrNr, DM_Density.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b, except for DM density, where -1 if a > b
 */
int subfind_compare_submp_GrNr_DM_Density(const void *a, const void *b)
{
  if(((struct submp_data *)a)->GrNr < ((struct submp_data *)b)->GrNr)
    return -1;

  if(((struct submp_data *)a)->GrNr > ((struct submp_data *)b)->GrNr)
    return +1;

#ifdef ADD_GROUP_PROPERTIES
  if(((struct submp_data *)a)->OriginalSubNr < ((struct submp_data *)b)->OriginalSubNr)
    return -1;

  if(((struct submp_data *)a)->OriginalSubNr > ((struct submp_data *)b)->OriginalSubNr)
    return +1;
#endif

  if(((struct submp_data *)a)->DM_Density > ((struct submp_data *)b)->DM_Density)
    return -1;

  if(((struct submp_data *)a)->DM_Density < ((struct submp_data *)b)->DM_Density)
    return +1;

  return 0;
}

/*! \brief Comparison function for submp_data objects.
 *
 *  Sorting kernel comparing element OldIndex.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int subfind_compare_submp_OldIndex(const void *a, const void *b)
{
  if(((struct submp_data *)a)->OldIndex < ((struct submp_data *)b)->OldIndex)
    return -1;

  if(((struct submp_data *)a)->OldIndex > ((struct submp_data *)b)->OldIndex)
    return +1;

  return 0;
}

/*! \brief Comparison function for id_list objects.
 *
 *  Sorting kernel comparing elements (most important first):
 *  GrNr, SubNr, Type, BindingEgy.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int subfind_compare_ID_list(const void *a, const void *b)
{
  if(((struct id_list *)a)->GrNr < ((struct id_list *)b)->GrNr)
    return -1;

  if(((struct id_list *)a)->GrNr > ((struct id_list *)b)->GrNr)
    return +1;

  if(((struct id_list *)a)->SubNr < ((struct id_list *)b)->SubNr)
    return -1;

  if(((struct id_list *)a)->SubNr > ((struct id_list *)b)->SubNr)
    return +1;

  if(((struct id_list *)a)->Type < ((struct id_list *)b)->Type)
    return -1;

  if(((struct id_list *)a)->Type > ((struct id_list *)b)->Type)
    return +1;

  if(((struct id_list *)a)->BindingEgy < ((struct id_list *)b)->BindingEgy)
    return -1;

  if(((struct id_list *)a)->BindingEgy > ((struct id_list *)b)->BindingEgy)
    return +1;

  return 0;
}

/*! \brief Comparison function for subgroup_properties objects.
 *
 *  Sorting kernel comparing elements (most important first):
 *  GrNr and SubNr.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int subfind_compare_SubGroup_GrNr_SubNr(const void *a, const void *b)
{
  if(((struct subgroup_properties *)a)->GrNr < ((struct subgroup_properties *)b)->GrNr)
    return -1;

  if(((struct subgroup_properties *)a)->GrNr > ((struct subgroup_properties *)b)->GrNr)
    return +1;

  if(((struct subgroup_properties *)a)->SubNr < ((struct subgroup_properties *)b)->SubNr)
    return -1;

  if(((struct subgroup_properties *)a)->SubNr > ((struct subgroup_properties *)b)->SubNr)
    return +1;

  return 0;
}

/*! \brief Comparison function for sort_r2list objects.
 *
 *  Sorting kernel comparing element r.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int subfind_compare_dist_rotcurve(const void *a, const void *b)
{
  if(((sort_r2list *)a)->r < ((sort_r2list *)b)->r)
    return -1;

  if(((sort_r2list *)a)->r > ((sort_r2list *)b)->r)
    return +1;

  return 0;
}

#if defined(MHD) && defined(ADD_MAGNETIC_GROUP_PROPERTIES)
int subfind_compare_rlist_mhd(const void *a, const void *b)
{
  if(((struct rlist_mhd *)a)->r < ((struct rlist_mhd *)b)->r)
    return -1;

  if(((struct rlist_mhd *)a)->r > ((struct rlist_mhd *)b)->r)
    return +1;

  return 0;
}
#endif

/*! \brief Comparison function for variables of type double.
 *
 *  Sorting kernel.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int subfind_compare_binding_energy(const void *a, const void *b)
{
  if(*((double *)a) > *((double *)b))
    return -1;

  if(*((double *)a) < *((double *)b))
    return +1;

  return 0;
}

/*! \brief Comparison function for cand_dat objects.
 *
 *  Sorting kernel comparing elements (most important first):
 *  bound_length and rank.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b, excpet bound length, where -1 if a > b.
 */
int subfind_compare_serial_candidates_boundlength(const void *a, const void *b)
{
  if(((struct cand_dat *)a)->bound_length > ((struct cand_dat *)b)->bound_length)
    return -1;

  if(((struct cand_dat *)a)->bound_length < ((struct cand_dat *)b)->bound_length)
    return +1;

  if(((struct cand_dat *)a)->rank < ((struct cand_dat *)b)->rank)
    return -1;

  if(((struct cand_dat *)a)->rank > ((struct cand_dat *)b)->rank)
    return +1;

  return 0;
}

/*! \brief Comparison function for cand_dat objects.
 *
 *  Sorting kernel comparing elements (most important first):
 *  rank and len.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b, except for len where -1 if a>b.
 */
int subfind_compare_serial_candidates_rank(const void *a, const void *b)
{
  if(((struct cand_dat *)a)->rank < ((struct cand_dat *)b)->rank)
    return -1;

  if(((struct cand_dat *)a)->rank > ((struct cand_dat *)b)->rank)
    return +1;

  if(((struct cand_dat *)a)->len > ((struct cand_dat *)b)->len)
    return -1;

  if(((struct cand_dat *)a)->len < ((struct cand_dat *)b)->len)
    return +1;

  return 0;
}

/*! \brief Comparison function for cand_dat objects.
 *
 *  Sorting kernel comparing element subnr.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int subfind_compare_serial_candidates_subnr(const void *a, const void *b)
{
  if(((struct cand_dat *)a)->subnr < ((struct cand_dat *)b)->subnr)
    return -1;

  if(((struct cand_dat *)a)->subnr > ((struct cand_dat *)b)->subnr)
    return +1;

  return 0;
}

/*! \brief Comparison function for coll_cand_dat objects.
 *
 *  Sorting kernel comparing element subnr.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int subfind_compare_coll_candidates_subnr(const void *a, const void *b)
{
  if(((struct coll_cand_dat *)a)->subnr < ((struct coll_cand_dat *)b)->subnr)
    return -1;

  if(((struct coll_cand_dat *)a)->subnr > ((struct coll_cand_dat *)b)->subnr)
    return +1;

  return 0;
}

/*! \brief Comparison function for coll_cand_dat objects.
 *
 *  Sorting kernel comparing element nsub.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int subfind_compare_coll_candidates_nsubs(const void *a, const void *b)
{
  if(((struct coll_cand_dat *)a)->nsub < ((struct coll_cand_dat *)b)->nsub)
    return -1;

  if(((struct coll_cand_dat *)a)->nsub > ((struct coll_cand_dat *)b)->nsub)
    return +1;

  return 0;
}

/*! \brief Comparison function for coll_cand_dat objects.
 *
 *  Sorting kernel comparing elements (most important first):
 *  bound_length, rank.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b, except for bound length where -1 if a > b.
 */
int subfind_compare_coll_candidates_boundlength(const void *a, const void *b)
{
  if(((struct coll_cand_dat *)a)->bound_length > ((struct coll_cand_dat *)b)->bound_length)
    return -1;

  if(((struct coll_cand_dat *)a)->bound_length < ((struct coll_cand_dat *)b)->bound_length)
    return +1;

  if(((struct coll_cand_dat *)a)->rank < ((struct coll_cand_dat *)b)->rank)
    return -1;

  if(((struct coll_cand_dat *)a)->rank > ((struct coll_cand_dat *)b)->rank)
    return +1;

  return 0;
}

/*! \brief Comparison function for coll_cand_dat objects.
 *
 *  Sorting kernel comparing elements (most important first):
 *  rank and len.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b, except for len, where -1 if a > b
 */
int subfind_compare_coll_candidates_rank(const void *a, const void *b)
{
  if(((struct coll_cand_dat *)a)->rank < ((struct coll_cand_dat *)b)->rank)
    return -1;

  if(((struct coll_cand_dat *)a)->rank > ((struct coll_cand_dat *)b)->rank)
    return +1;

  if(((struct coll_cand_dat *)a)->len > ((struct coll_cand_dat *)b)->len)
    return -1;

  if(((struct coll_cand_dat *)a)->len < ((struct coll_cand_dat *)b)->len)
    return +1;

  return 0;
}

#endif
