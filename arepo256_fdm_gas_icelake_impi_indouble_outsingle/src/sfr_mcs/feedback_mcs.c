/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sfr_mcs/feedback_mcs.c
 * \date        02/2019
 * \author     	Matthew C Smith
 * \brief
 * \details     Wrapper for stellar feedback modules
 *
 *
 * \par Major modifications and contributions:
 *
 * - 17.01.2022 Ported into current codebase
 */

#include "../allvars.h"
#include "../proto.h"

#if defined(SN_MCS) || defined(HII_MCS) || defined(PE_MCS)

void init_stellar_feedback(void)
{
  All.MaxFBStarMass = All.TargetGasMass * All.MaxFBStarMassFac;

  All.MaxFBStarEvalTimestep *= (SEC_PER_MEGAYEAR / All.UnitTime_in_s * All.HubbleParam);
#ifndef IMF_SAMPLING_MCS
  All.MaxFBStarEvalTimestepCut *= 1e6;  // Now in years
#endif

#ifdef HII_MCS
  init_hii();
#endif

#ifdef PE_MCS
  init_pe();
#endif

#ifdef SN_MCS
  init_sne();
#endif
}

void do_stellar_feedback(void)
{
#ifdef IMF_SAMPLING_MCS
  check_for_dead_stars();
#endif

#ifdef HII_MCS
  place_hii_regions();
#endif

#ifdef SN_MCS
  do_supernovae();
#endif

#ifdef PE_MCS
  update_FUV_luminosities();
#endif
}
#endif
