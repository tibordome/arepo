/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sne/sne.h
 * \date        MM/YYYY
 * \author      Robin Tress, Rowan Smith, Andre Bubel
 * \brief
 * \details     Please contact the authors at robin.tress@uni-heidelberg.de
 *              and rowan.smith@manchester.ac.uk before using this to avoid overlapping
 *              of projects. And please report any issue you may encounter by using this
 *              routine.
 *
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <stdio.h>

extern gsl_rng *sne_rng; /**< a random number generator used for sne concerns, seeded the same on all tasks */
extern FILE *FdSNe;      /** Supernova logfile **/

extern double SNERandomNextTime;
#ifdef CLUSTERED_SNE
extern double SNEClusterNextTime;
#endif

#ifdef GALPOT
double GlobalStellarMass;
double GlobalGasMass;
int NSNetot;
struct snPositions
{
  double pos[3];
} * SNePositions;
#endif

#if defined(SINK_PARTICLES) && defined(SINK_PARTICLES_FEEDBACK)
extern int Sink_about_to_explode;
extern int Star_about_to_explode;
extern int Wherewasi;
#endif

/*variables and parameters needed for SNEInjectionCriterion = 4
 * i.e. to define the cumulative distribution function of the radial supernova distribution*/
#define SNEPDF_NBINS 1000
#define SNEPDF_RMAX 20. * KILOPARSEC / All.UnitLength_in_cm

enum SNE_type
{
  NO_SNE = 0,

  SNE_ONLY_ONCE            = 1,
  SNE_RANDOM               = 2,
  SNE_CLUSTER              = 4,
  SNE_IN_DISC              = 8,
  SNE_GENERAL_DISTRIBUTION = 16,
  SNE_SINK                 = 32,
  SNE_POISSON              = 64,
  SNE_DISC_PARTICLE        = 128
};
