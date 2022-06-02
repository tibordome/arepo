/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/GFM/helper.c
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

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"
#include "../voronoi.h"

void check_AuxDataID_references(void)
{
  mpi_printf("GFM: GFM Checks...\n");

  int i;
  for(i = 0; i < NumPart; i++)
    {
#ifdef GFM
      if(P[i].Type == 4)
        if(StarP[P[i].AuxDataID].PID != i)
          {
            printf("StarP broken: %llu %d %d %d %d %d\n", (long long)P[i].AuxDataID, i, StarP[P[i].AuxDataID].PID, NumGas, N_star,
                   NumPart);
            terminate("StarP[P[i].AuxDataID].PID!=i\n");
          }
#endif

#ifdef BLACK_HOLES
      if(P[i].Type == 5)
        if(BHP[P[i].AuxDataID].PID != i)
          {
            printf("BHP broken: %llu %d %d %d %d %d\n", (long long)P[i].AuxDataID, i, BHP[P[i].AuxDataID].PID, NumGas, NumBHs,
                   NumPart);
            terminate("BHP[P[i].AuxDataID].PID!=i\n");
          }
#endif
    }

  mpi_printf("GFM: done.\n");
}

#ifdef GFM
void gfm_inject_into_cell(int j, double inj_mass, double inj_thermalenergy, double inj_mom[3])
{
  if(inj_mass <= 0)
    terminate("GFM: inj_mass<=0");

  /* remove current kinetic energy from total energy */
  SphP[j].Energy -= 0.5 *
                    (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                     SphP[j].Momentum[2] * SphP[j].Momentum[2]) /
                    P[j].Mass;

  /* add injected momentum to cell */
  SphP[j].Momentum[0] += inj_mom[0];
  SphP[j].Momentum[1] += inj_mom[1];
  SphP[j].Momentum[2] += inj_mom[2];

  /* add injected mass to cell */
  P[j].Mass += inj_mass;

  /* add injected thermal energy to cell */
  SphP[j].Energy += inj_thermalenergy;

  /* add new kinetic energy to total energy */
  SphP[j].Energy += 0.5 *
                    (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                     SphP[j].Momentum[2] * SphP[j].Momentum[2]) /
                    P[j].Mass;

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  SphP[j].Utherm = (SphP[j].Energy - 0.5 *
                                         (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                                          SphP[j].Momentum[2] * SphP[j].Momentum[2]) /
                                         P[j].Mass) /
                   P[j].Mass / (All.cf_atime * All.cf_atime);
  SphP[j].A       = (GAMMA - 1.0) * SphP[j].Utherm / pow(SphP[j].Density * All.cf_a3inv, GAMMA - 1);
  SphP[j].Entropy = log(SphP[j].A) * P[j].Mass;
#endif
}
#endif

#ifdef GFM
/* add a star particle to StarP */
void gfm_add_star(int i, int j, MyDouble mass_of_star, MyFloat birthtime, MyFloat hsml_guess)
{
  if(N_star >= All.MaxPartStar)
    terminate("There is no space left to creat new stars. N_star = %d, MaxPartStar = %d", N_star, All.MaxPartStar);

#ifdef GFM_STELLAR_EVOLUTION
  int k_elem;
#endif
  P[i].AuxDataID = N_star;

  /* zero StarP[] entries */
  memset(&StarP[N_star], 0, sizeof(struct star_particle_data));

  /* set values */
  StarP[N_star].PID       = i;
  StarP[N_star].BirthTime = birthtime;

  StarP[N_star].BirthPos[0] = P[i].Pos[0];
  StarP[N_star].BirthPos[1] = P[i].Pos[1];
  StarP[N_star].BirthPos[2] = P[i].Pos[2];
  domain_displacePosition(StarP[N_star].BirthPos, DISPLACE_POSITION_BACKWARD);

  StarP[N_star].BirthVel[0]  = P[i].Vel[0];
  StarP[N_star].BirthVel[1]  = P[i].Vel[1];
  StarP[N_star].BirthVel[2]  = P[i].Vel[2];
  StarP[N_star].BirthDensity = SphP[j].Density;

#ifdef GFM_RPROCESS_CHANNELS_NS_KICKS
  StarP[N_star].NSNS_channel = -1;
#endif

#ifdef GFM_STELLAR_EVOLUTION
  StarP[N_star].InitialMass = mass_of_star;
  StarP[N_star].Hsml        = hsml_guess;

  for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
    StarP[N_star].MassMetals[k_elem] = (mass_of_star / P[j].Mass) * SphP[j].MassMetals[k_elem];

#ifdef GFM_STELLAR_EVOLUTION_NO_ELEMENTS
  StarP[N_star].Metallicity = SphP[j].Metallicity;
#else
  double metmass = 0;
  for(k_elem = 2; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
    metmass += StarP[N_star].MassMetals[k_elem];

  StarP[N_star].Metallicity = metmass / mass_of_star;
#endif

#ifdef GFM_CHEMTAGS
  for(k_elem = 0; k_elem < GFM_N_CHEM_TAGS; k_elem++)
    StarP[N_star].MassMetalsChemTags[k_elem] = (mass_of_star / P[j].Mass) * SphP[j].MassMetalsChemTags[k_elem];
#endif

#ifdef GFM_RPROCESS_CHANNELS
  for(k_elem = 0; k_elem < GFM_RPROCESS_CHANNELS; k_elem++)
    {
      StarP[N_star].MassRProcess[k_elem]        = mass_of_star * SphP[j].MassRProcess[k_elem] / P[j].Mass;
      StarP[N_star].NRProcessInjections[k_elem] = 0;
    }
#endif
#ifdef GFM_SNIA_ENERGY_INJECTION
  StarP[N_star].NumSNIa = 0;
#endif

#ifdef GFM_DUST
  /* Assume that a star particle also starts with initial metals */
  /* due to dust in the ISM. */
  int chan;
  for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
    {
      for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
        {
          StarP[N_star].MassMetals[k_elem] += mass_of_star * SphP[j].MetalsDustFraction[chan][k_elem];
          StarP[N_star].Metallicity += SphP[j].MetalsDustFraction[chan][k_elem];
        }
    }

  for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
    {
      double tot_mass = 0.0;
      tot_mass += SphP[j].MassMetals[k_elem];
      for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
        {
          tot_mass += SphP[j].MassMetalsDust[chan][k_elem];
        }

      /* Normalize the star's initial metal and dust fractions.  If not */
      /* normalizable, just set to zero: the  ISM had no gas-phase or dust */
      /* metals. */
      if(tot_mass > 0.0)
        {
          StarP[N_star].InitialMetalFractions[k_elem] = SphP[j].MassMetals[k_elem] / tot_mass;
          for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
            {
              StarP[N_star].InitialDustFractions[chan][k_elem] = SphP[j].MassMetalsDust[chan][k_elem] / tot_mass;
            }
        }
      else
        {
          StarP[N_star].InitialMetalFractions[k_elem] = 0.0;
          for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
            {
              StarP[N_star].InitialDustFractions[chan][k_elem] = 0.0;
            }
        }
    }
#endif

  StarP[N_star].lastEnrichTime = birthtime;

#ifdef GFM_DISCRETE_ENRICHMENT
#ifdef SMUGGLE_DISCRETE_SN
  StarP[N_star].SNTime = birthtime;
#endif
#endif

#ifdef SMUGGLE_AGB_WINDS
  StarP[N_star].OBTime = birthtime;
#endif

#if defined(GFM_VARIABLE_IMF) && (GFM_VARIABLE_IMF == 0)
  StarP[N_star].DMVelDisp = SphP[j].w.DMVelDisp;
#endif
#endif

#ifdef SMUGGLE_STAR_FEEDBACK
#ifdef SMUGGLE_OUTPUT_STELLAR_FEEDBACK
  StarP[N_star].SNII_Num                   = 0;
  StarP[N_star].SNIa_Num                   = 0;
  StarP[N_star].FeedbackEnergy             = 0;
  StarP[N_star].FeedbackKinEnergy          = 0;
  StarP[N_star].FeedbackThEnergy           = 0;
  StarP[N_star].FeedbackMomentum           = 0;
  StarP[N_star].FeedbackMomentumAGB        = 0;
  StarP[N_star].TotalMassReleased          = 0;
  StarP[N_star].TotalMassToEnrich          = 0;
  StarP[N_star].Cum_SNII_Num               = 0;
  StarP[N_star].Cum_SNIa_Num               = 0;
  StarP[N_star].Cum_FeedbackEnergy         = 0;
  StarP[N_star].Cum_FeedbackMomentum       = 0;
  StarP[N_star].Cum_InjFeedbackMomentum    = 0;
  StarP[N_star].Cum_InjFeedbackMomentumAGB = 0;
  StarP[N_star].MaxFeedRadius              = -1;
#endif

#ifdef SMUGGLE_RADIATION_FEEDBACK
  StarP[N_star].RadFeed_Flag     = 0; /* it changes to -1 once the momentum due to radiation is input */
  StarP[N_star].StromgrenRadius  = 0;
  StarP[N_star].GasColumnDensity = 0;
#ifdef SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION
  StarP[N_star].StromgrenMass = 0;
#endif
  StarP[N_star].RadFeed_NormSphRad = 0;
  StarP[N_star].RadFeedTau         = 0;
  StarP[N_star].RadFeed_NumNgb     = 0;
#ifdef SMUGGLE_RADIATION_FEEDBACK_DEBUG /* #LVS: TODO --> once tested, keep variables in "debug" for only StarParticle */
  StarP[N_star].RadiationMomentumReleased     = 0;
  StarP[N_star].NormSphRadFeedback            = 0;
  StarP[N_star].Cum_RadiationMomentumReleased = 0;
  StarP[N_star].Cum_RadMomentumRealInjected   = 0;
#ifdef SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION
  StarP[N_star].NormSphRadFeedback_cold = 0;
#endif
  StarP[N_star].RadCoolShutoffTime = 0;
#endif
#endif
#endif

  N_star++;
}

#ifdef GFM_DUST_CAP
void gfm_cap_dust(void)
{
  for(int i = 0; i < NumGas; i++)
    {
      double sum = 0.0;
      for(int chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
        {
          for(int k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
            {
              if(SphP[i].MetalsDustFraction[chan][k_elem] < 0)
                {
                  SphP[i].MetalsDustFraction[chan][k_elem] = 0.0;
                }
              if(SphP[i].MassMetalsDust[chan][k_elem] < 0)
                {
                  sum += -SphP[i].MassMetalsDust[chan][k_elem];
                  SphP[i].MassMetalsDust[chan][k_elem] = 0.0;
                }
            }
        }
      SphP[i].DustMassCap += sum;
    }
}

void gfm_check_dust(int checkid)
{
  for(int i = 0; i < NumGas; i++)
    {
      for(int chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
        {
          for(int k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
            {
              if(SphP[i].MetalsDustFraction[chan][k_elem] < 0)
                terminate("GFM_DUST: %d negative dust %g\n", checkid, SphP[i].MetalsDustFraction[chan][k_elem]);
              if(SphP[i].MassMetalsDust[chan][k_elem] < 0)
                terminate("GFM_DUST: %d negative dust %g\n", checkid, SphP[i].MassMetalsDust[chan][k_elem]);
            }
        }
    }
}

#endif

#endif  // GFM
