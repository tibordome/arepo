/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/GFM/stellar_evolution_main.c
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
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef GFM_STELLAR_EVOLUTION

#ifdef TRACER_MC
#define MAXLEN_TRACERTARGET_LIST 6
typedef struct
{
  int attach_to_task;
  int attach_to_index;
#ifdef TRACER_MC_CHECKS
  MyIDType attach_to_ID;
#endif
} tracertarget_data;
#endif

static int enrich_evaluate(int target, int mode, int threadid);

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Vel[3];
  MyFloat Hsml;
  MyFloat NormSph;
  MyFloat MetalsReleased[GFM_N_CHEM_ELEMENTS];
#ifdef GFM_DUST
  MyFloat DustReleased[GFM_DUST_N_CHANNELS][GFM_N_CHEM_ELEMENTS];
#endif
#if defined(GFM_DUST) || defined(DUST_LIVE)
  MyFloat NumSNII;
#endif
#ifdef GFM_CHEMTAGS
  MyFloat MetalReleasedChemTags[GFM_N_CHEM_TAGS];
#endif
#ifdef GFM_RPROCESS_CHANNELS
  MyFloat MassReleasedRProcess[GFM_RPROCESS_CHANNELS];
#endif
#ifdef GFM_SNIA_ENERGY_INJECTION
  MyIDType ID;
  MyFloat SNIaEnergy;
#endif
  MyFloat TotalMetalMassReleased;
  MyFloat TotalMassReleased;
#ifdef TRACER_MC
  MyFloat Mass;
  int NumberOfTracers;
#endif
#ifdef GFM_INJECT_B_FROM_SN
  MyFloat m[3];
#endif

  int Firstnode;
} data_in;

static data_in *DataGet;

/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in *in, int i, int firstnode)
{
  int iel;

  in->Pos[0]                 = P[StarParticle[i].index].Pos[0];
  in->Pos[1]                 = P[StarParticle[i].index].Pos[1];
  in->Pos[2]                 = P[StarParticle[i].index].Pos[2];
  in->Vel[0]                 = P[StarParticle[i].index].Vel[0];
  in->Vel[1]                 = P[StarParticle[i].index].Vel[1];
  in->Vel[2]                 = P[StarParticle[i].index].Vel[2];
  in->Hsml                   = STP(StarParticle[i].index).Hsml;
  in->NormSph                = StarParticle[i].NormSph;
  in->TotalMassReleased      = StarParticle[i].TotalMassReleased;
  in->TotalMetalMassReleased = StarParticle[i].TotalMetalMassReleased;
  for(iel = 0; iel < GFM_N_CHEM_ELEMENTS; iel++)
    {
      in->MetalsReleased[iel] = StarParticle[i].MetalMassReleased[iel];
    }
#ifdef GFM_DUST
  for(int chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
    {
      for(iel = 0; iel < GFM_N_CHEM_ELEMENTS; iel++)
        {
          in->DustReleased[chan][iel] = StarParticle[i].DustMassReleased[chan][iel];
        }
    }
#endif
#ifdef GFM_SNIA_ENERGY_INJECTION
  in->ID         = P[StarParticle[i].index].ID;
  in->SNIaEnergy = 1e51 / All.UnitEnergy_in_cgs * StarParticle[i].NumSNIa;
#endif
#if defined(GFM_DUST) || defined(DUST_LIVE)
  in->NumSNII = StarParticle[i].NumSNII;
#endif
#ifdef GFM_CHEMTAGS
  for(iel = 0; iel < GFM_N_CHEM_TAGS; iel++)
    {
      in->MetalReleasedChemTags[iel] = StarParticle[i].MetalMassReleasedChemTags[iel];
    }
#endif
#ifdef GFM_RPROCESS_CHANNELS
  for(iel = 0; iel < GFM_RPROCESS_CHANNELS; iel++)
    {
      in->MassReleasedRProcess[iel] = StarParticle[i].MassReleasedRProcess[iel];
    }
#endif
#ifdef TRACER_MC
  in->Mass            = P[StarParticle[i].index].Mass;
  in->NumberOfTracers = P[StarParticle[i].index].NumberOfTracers;
#endif
#ifdef GFM_INJECT_B_FROM_SN
  in->m[0] = StarParticle[i].m[0];
  in->m[1] = StarParticle[i].m[1];
  in->m[2] = StarParticle[i].m[2];
#endif

  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
#ifdef GFM_SNIA_ENERGY_INJECTION
  int CountSNIa;
#endif
#ifdef TRACER_MC
  int ntracertargets;
  tracertarget_data tracertargets[MAXLEN_TRACERTARGET_LIST];
#else
  char dummy;
#endif
} data_out;

static data_out *DataResult;

/* routine to store or combine result data */
static void out2particle(data_out *out, int i, int mode)
{
#ifdef GFM_SNIA_ENERGY_INJECTION
  STP(StarParticle[i].index).NumSNIa += out->CountSNIa;
#endif

  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
#ifdef TRACER_MC
      if(out->ntracertargets)
        {
          for(int k = 0; k < out->ntracertargets; k++)
            {
#ifdef TRACER_MC_CHECKS
              move_one_tracer(StarParticle[i].index, out->tracertargets[k].attach_to_task, out->tracertargets[k].attach_to_index,
                              out->tracertargets[k].attach_to_ID);
#else
              move_one_tracer(StarParticle[i].index, out->tracertargets[k].attach_to_task, out->tracertargets[k].attach_to_index, 0);
#endif
            }
        }
#endif
    }
  else /* merge */
    {
#ifdef TRACER_MC
      if(out->ntracertargets)
        {
          for(int k = 0; k < out->ntracertargets; k++)
            {
#ifdef TRACER_MC_CHECKS
              move_one_tracer(StarParticle[i].index, out->tracertargets[k].attach_to_task, out->tracertargets[k].attach_to_index,
                              out->tracertargets[k].attach_to_ID);
#else
              move_one_tracer(StarParticle[i].index, out->tracertargets[k].attach_to_task, out->tracertargets[k].attach_to_index, 0);
#endif
            }
        }
#endif
    }
}

#include "../generic_comm_helpers2.h"

static int Ncount;

static void kernel_local(void)
{
  int i;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif

#pragma omp parallel private(i)
  {
    int j, threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              if(generic_polling_primary(count, Ncount))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= Ncount)
          break;

        enrich_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}

static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, cnt = 0;
#pragma omp parallel private(i)
  {
    int threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    while(1)
      {
#pragma omp atomic capture
        i = cnt++;

        if(i >= Nimport)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              generic_polling_secondary();
          }

        count++;
#endif

        enrich_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void do_chemical_enrichment(void)
{
  Ncount = Nstar;

  long long ntot;
  sumup_large_ints(1, &Nstar, &ntot);
  if(ntot == 0)
    return;

  generic_set_MaxNexport();

  double t0 = second();

#ifdef TRACER_MC
  start_MC_tracer(N_tracer); /* allocate buffer for tracer exchange */
#endif

  generic_comm_pattern(Nstar, kernel_local, kernel_imported);

#ifdef TRACER_MC
  finish_MC_tracer();
#endif

  double t1 = second();
  mpi_printf("GFM_STELLAR_EVOLUTION: enriching from stars took %g sec\n", timediff(t0, t1));
}

#ifdef GFM_WINDS_STRIPPING
void strip_active_winds(void)
{
  int i, iel;
  double dm;

  /* do local star particles */
  for(i = 0; i < Nwinds_to_strip; i++)
    {
      /* do not strip the wind particle if it has nagative metallicity */
      if(STP(StarParticle[i].index).Metallicity < 0.0)
        {
          StarParticle[i].TotalMassReleased      = 0;
          StarParticle[i].TotalMetalMassReleased = 0;
          continue;
        }

#ifndef GFM_DUST
      dm            = 0;
      double summet = 0;

#ifndef GFM_STELLAR_EVOLUTION_NO_ELEMENTS
      /* do individual heavy elements */
      for(iel = 2; iel < GFM_N_CHEM_ELEMENTS; iel++)
        {
          StarParticle[i].MetalMassReleased[iel] = All.WindDumpFactor * STP(StarParticle[i].index).MassMetals[iel];
          dm += All.WindDumpFactor * STP(StarParticle[i].index).MassMetals[iel];
          STP(StarParticle[i].index).MassMetals[iel] *= (1 - All.WindDumpFactor);
          summet += STP(StarParticle[i].index).MassMetals[iel];
        }

      P[StarParticle[i].index].Mass -= dm;
      STP(StarParticle[i].index).Metallicity = summet / P[StarParticle[i].index].Mass;

      StarParticle[i].TotalMassReleased      = dm;
      StarParticle[i].TotalMetalMassReleased = dm;
#else
      dm = All.WindDumpFactor * STP(StarParticle[i].index).Metallicity * (MyFloat)P[StarParticle[i].index].Mass;
      P[StarParticle[i].index].Mass -= dm;
      STP(StarParticle[i].index).Metallicity *=
          (1.0 - All.WindDumpFactor) / (1.0 - All.WindDumpFactor * STP(StarParticle[i].index).Metallicity);

      StarParticle[i].TotalMassReleased      = dm;
      StarParticle[i].TotalMetalMassReleased = dm;
#endif
#else /* ifdef GFM_DUST */

      dm = All.WindDumpFactor * STP(StarParticle[i].index).Metallicity * (MyFloat)P[StarParticle[i].index].Mass;
      P[StarParticle[i].index].Mass -= dm;
      STP(StarParticle[i].index).Metallicity *=
          (1.0 - All.WindDumpFactor) / (1.0 - All.WindDumpFactor * STP(StarParticle[i].index).Metallicity);

      StarParticle[i].TotalMassReleased = dm;
      StarParticle[i].TotalMetalMassReleased = dm;

      for(iel = 2; iel < GFM_N_CHEM_ELEMENTS; iel++)
        {
          double released_metal = All.WindDumpFactor * STP(StarParticle[i].index).MassMetals[iel];
          STP(StarParticle[i].index).MassMetals[iel] *= (1 - All.WindDumpFactor);

          StarParticle[i].MetalMassReleased[iel] = released_metal * STP(StarParticle[i].index).InitialMetalFractions[iel];
          for(int chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
            {
              double released_dust = released_metal * STP(StarParticle[i].index).InitialDustFractions[chan][iel];
              StarParticle[i].DustMassReleased[chan][iel] = released_dust;
              /* Total metal mass released should only refer to gas-phase. */
              StarParticle[i].TotalMetalMassReleased -= released_dust;
            }
        }
#endif
#ifdef GFM_CHEMTAGS
      for(iel = 0; iel < GFM_N_CHEM_TAGS; iel++)
        {
          StarParticle[i].MetalMassReleasedChemTags[iel] = All.WindDumpFactor * STP(StarParticle[i].index).MassMetalsChemTags[iel];
          STP(StarParticle[i].index).MassMetalsChemTags[iel] *= (1 - All.WindDumpFactor);
        }
#endif
#ifdef GFM_RPROCESS_CHANNELS
      for(iel = 0; iel < GFM_RPROCESS_CHANNELS; iel++)
        {
          StarParticle[i].MassReleasedRProcess[iel] = All.WindDumpFactor * STP(StarParticle[i].index).MassRProcess[iel];
          STP(StarParticle[i].index).MassRProcess[iel] *= (1 - All.WindDumpFactor);
        }
#endif
#ifdef GFM_SNIA_ENERGY_INJECTION
      StarParticle[i].NumSNIa = 0;
#endif
    }
}

void do_chemical_stripping(void)
{
  Ncount = Nwinds_to_strip;

  generic_set_MaxNexport();

  double t0 = second();

#ifdef TRACER_MC
  start_MC_tracer(N_tracer); /* allocate buffer for tracer exchange */
#endif

  generic_comm_pattern(Nwinds_to_strip, kernel_local, kernel_imported);

#ifdef TRACER_MC
  finish_MC_tracer();
#endif

  double t1 = second();

  mpi_printf("GFM_STELLAR_EVOLUTION: stripping from winds took %g sec\n", timediff(t0, t1));
}
#endif

#ifdef TRACER_MC
static int sort_probs_kernel(const void *a, const void *b)
{
  if(*((double *)a) < *((double *)b))
    return -1;

  if(*((double *)a) > *((double *)b))
    return +1;

  return 0;
}
#endif

static int enrich_evaluate(int target, int mode, int threadid)
{
  int j, n, iel;
  int numnodes, *firstnode;
  double h;
#if !defined(GFM_TOPHAT_KERNEL) && !defined(SMUGGLE_OMEGA_WEIGHT_SN)
  double wk, u, r, hinv, hinv3;
#endif
  double weight_fac;
  MyDouble *pos;
  MyFloat *vel;
  MyFloat dm_total, dm_metals, dm_metal[GFM_N_CHEM_ELEMENTS];
  MyFloat TotalMassReleased, TotalMetalMassReleased, MetalMassReleased[GFM_N_CHEM_ELEMENTS];
#ifdef GFM_DUST
  MyFloat DustMassReleased[GFM_DUST_N_CHANNELS][GFM_N_CHEM_ELEMENTS], dm_metal_dust[GFM_DUST_N_CHANNELS][GFM_N_CHEM_ELEMENTS];
#endif
#if defined(GFM_DUST) || defined(DUST_LIVE)
  MyFloat NumSNII, dNumSNII;
#endif
#ifdef GFM_CHEMTAGS
  MyFloat MetalMassReleasedChemTags[GFM_N_CHEM_TAGS], dm_metal_tags[GFM_N_CHEM_TAGS];
#endif
#ifdef GFM_RPROCESS_CHANNELS
  MyFloat MassReleasedRProcess[GFM_RPROCESS_CHANNELS], dm_metal_rprocess[GFM_RPROCESS_CHANNELS];
#endif
#ifdef GFM_SNIA_ENERGY_INJECTION
  MyFloat SNIaEnergy;
#endif
#ifdef GFM_INJECT_B_FROM_SN
  MyDouble *m;
#endif
  MyFloat normsph;

  data_in local, *in;
  data_out out;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      in = &local;

      numnodes  = 1;
      firstnode = NULL;
    }
  else
    {
      in = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

#ifdef GFM_DISCRETE_ENRICHMENT
  /* if this particle has no mass to distribute, quit before we actually do the neighbor search */
  if(in->TotalMassReleased == 0.0)
    return 0;
#endif

  pos                    = in->Pos;
  vel                    = in->Vel;
  h                      = in->Hsml;
  normsph                = in->NormSph;
  TotalMassReleased      = in->TotalMassReleased;
  TotalMetalMassReleased = in->TotalMetalMassReleased;
  for(iel = 0; iel < GFM_N_CHEM_ELEMENTS; iel++)
    {
      MetalMassReleased[iel] = in->MetalsReleased[iel];
    }
#ifdef GFM_DUST
  for(int chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
    {
      for(iel = 0; iel < GFM_N_CHEM_ELEMENTS; iel++)
        {
          DustMassReleased[chan][iel] = in->DustReleased[chan][iel];
        }
    }
#endif
#if defined(GFM_DUST) || defined(DUST_LIVE)
  NumSNII = in->NumSNII;
#endif
#ifdef GFM_CHEMTAGS
  for(iel = 0; iel < GFM_N_CHEM_TAGS; iel++)
    {
      MetalMassReleasedChemTags[iel] = in->MetalReleasedChemTags[iel];
    }
#endif
#ifdef GFM_RPROCESS_CHANNELS
  for(iel = 0; iel < GFM_RPROCESS_CHANNELS; iel++)
    {
      MassReleasedRProcess[iel] = in->MassReleasedRProcess[iel];
    }
#endif
#ifdef GFM_SNIA_ENERGY_INJECTION
  SNIaEnergy    = in->SNIaEnergy;
  out.CountSNIa = 0;
#endif
#ifdef GFM_INJECT_B_FROM_SN
  m = in->m;
#endif

#if !defined(GFM_TOPHAT_KERNEL) && !defined(SMUGGLE_OMEGA_WEIGHT_SN)
  hinv = 1.0 / h;
#ifndef TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif
#endif

#ifdef TRACER_MC
  out.ntracertargets = 0;
#else
  out.dummy = 0;
#endif

#ifdef TRACER_MC
  double *prob_tracers = (double *)mymalloc("prob_tracers", in->NumberOfTracers * sizeof(double));
  for(int k = 0; k < in->NumberOfTracers; k++)
    prob_tracers[k] = get_random_number();
  double pcum1 = 0, pcum2 = 0;

  mysort(prob_tracers, in->NumberOfTracers, sizeof(double), sort_probs_kernel);

  int kpoint = 0;
#endif

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0) /* skip cells that have been swallowed or dissolved */
        {
#ifdef SMUGGLE_STAR_FEEDBACK
#ifdef BH_BASED_CGM_ZOOM
          if(P[j].Mass < 0.3 * All.TargetGasMass / All.CGM_RefinementFactor)
            continue;
#else
          if(P[j].Mass < 0.3 * All.TargetGasMass)
            continue;
#endif
#endif

#ifdef SMUGGLE_OMEGA_WEIGHT_SN
          double r2          = Thread[threadid].R2list[n];
          double cell_radius = get_cell_radius(j);
          double cell_area   = cell_radius * cell_radius;
          double omega       = 0.5 * (1.0 - sqrt(r2 / (r2 + cell_area)));
          weight_fac         = omega / normsph;
#elif !defined(GFM_TOPHAT_KERNEL)
          double r2 = Thread[threadid].R2list[n];

          r = sqrt(r2);
          u = r * hinv;

          if(u < 0.5)
            wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
          else
            wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

          weight_fac = SphP[j].Volume * wk / normsph;
#else
          weight_fac = SphP[j].Volume / normsph;
#endif

          dm_total  = weight_fac * TotalMassReleased;
          dm_metals = weight_fac * TotalMetalMassReleased;
          for(iel = 0; iel < GFM_N_CHEM_ELEMENTS; iel++)
            {
              dm_metal[iel] = weight_fac * MetalMassReleased[iel];
            }
#ifdef GFM_DUST
          for(int chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
            {
              for(iel = 0; iel < GFM_N_CHEM_ELEMENTS; iel++)
                {
                  dm_metal_dust[chan][iel] = weight_fac * DustMassReleased[chan][iel];
                }
            }
#endif
#if defined(GFM_DUST) || defined(DUST_LIVE)
          dNumSNII = weight_fac * NumSNII;
#endif
#ifdef GFM_CHEMTAGS
          for(iel = 0; iel < GFM_N_CHEM_TAGS; iel++)
            {
              dm_metal_tags[iel] = weight_fac * MetalMassReleasedChemTags[iel];
            }
#endif
#ifdef GFM_RPROCESS_CHANNELS
          for(iel = 0; iel < GFM_RPROCESS_CHANNELS; iel++)
            {
              dm_metal_rprocess[iel] = weight_fac * MassReleasedRProcess[iel];
            }
#endif

          /* add metals (primitive variables updated below) */
          SphP[j].MassMetallicity += dm_metals;
          for(iel = 0; iel < GFM_N_CHEM_ELEMENTS; iel++)
            {
              SphP[j].MassMetals[iel] += dm_metal[iel];
            }
#ifdef GFM_DUST
          for(int chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
            {
              for(iel = 0; iel < GFM_N_CHEM_ELEMENTS; iel++)
                {
                  SphP[j].MassMetalsDust[chan][iel] += dm_metal_dust[chan][iel];
                }
            }
#endif
#if defined(GFM_DUST) || defined(DUST_LIVE)
          SphP[j].NumSNII += dNumSNII;
#endif
#ifdef GFM_CHEMTAGS
          for(iel = 0; iel < GFM_N_CHEM_TAGS; iel++)
            {
              SphP[j].MassMetalsChemTags[iel] += dm_metal_tags[iel];
            }
#endif
#ifdef GFM_RPROCESS_CHANNELS
          for(iel = 0; iel < GFM_RPROCESS_CHANNELS; iel++)
            {
              SphP[j].MassRProcess[iel] += dm_metal_rprocess[iel];
            }
#endif

#if defined(REFINEMENT_HIGH_RES_GAS) || defined(REFINEMENT_CGM)
          double mass_fac = (P[j].Mass + dm_total) / P[j].Mass;
#endif

#ifndef SMUGGLE_SN_COOLING_RADIUS_BOOST /* mass enrichment is done in the feedback loop in this case */
          if(dm_total > 0)
            {
              double inj_mass, inj_thermalenergy, vrelsq, x, inj_mom[3], vrel[3];
              inj_mass   = dm_total;
              inj_mom[0] = dm_total * vel[0];
              inj_mom[1] = dm_total * vel[1];
              inj_mom[2] = dm_total * vel[2];

              /* compute thermal energy variation due to momentum injection */
              vrel[0]           = SphP[j].Momentum[0] / P[j].Mass - vel[0];
              vrel[1]           = SphP[j].Momentum[1] / P[j].Mass - vel[1];
              vrel[2]           = SphP[j].Momentum[2] / P[j].Mass - vel[2];
              vrelsq            = vrel[0] * vrel[0] + vrel[1] * vrel[1] + vrel[2] * vrel[2];
              x                 = inj_mass / P[j].Mass;
              inj_thermalenergy = 0.5 * P[j].Mass * x * vrelsq / (1. + x);

              /* compute thermal energy variation due to thermal energy of ejecta  */
              /* note: assuming FULL ionization of ejecta */
              double u_to_temp_fac = (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 *
                                     All.UnitEnergy_in_cgs / All.UnitMass_in_g;
              inj_thermalenergy += All.cf_atime * All.cf_atime * dm_total * GFM_STELLAR_EJECTA_TEMPERATURE / u_to_temp_fac;

              gfm_inject_into_cell(j, inj_mass, inj_thermalenergy, inj_mom);
            }
#endif

#ifdef GFM_SNIA_ENERGY_INJECTION
          // don't inject into starforming gas
          if(SphP[j].Sfr == 0 && weight_fac > 0 && SNIaEnergy > 0)
            {
              SphP[j].Energy += All.cf_atime * All.cf_atime * weight_fac * SNIaEnergy;

              double pos[3];
              for(int k = 0; k < 3; k++)
                pos[k] = P[j].Pos[k];
              domain_displacePosition(pos, DISPLACE_POSITION_BACKWARD);
              printf("SNIA: Task %d, t=%g, star ID=%lld adding %g erg to cell ID=%lld, Mass=%g, Rho=%g, Utherm=%g, Pos=%g,%g,%g\n",
                     ThisTask, All.Time, (long long)in->ID, SNIaEnergy * All.UnitEnergy_in_cgs, (long long)P[j].ID, P[j].Mass,
                     SphP[j].Density, SphP[j].Utherm, pos[0], pos[1], pos[2]);
              out.CountSNIa++;
            }
#endif

#ifdef REFINEMENT_HIGH_RES_GAS
          /* mass scale factor to new total mass */
          if(dm_total > 0)
            SphP[j].HighResMass *= mass_fac;
#endif
#ifdef REFINEMENT_CGM
          if(dm_total > 0)
            SphP[j].HighResMassCGM *= mass_fac;
#endif

#ifdef GFM_INJECT_B_FROM_SN
          double mstrength = sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
          if(mstrength > 0.0)
            {
              double dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
              double dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
              double dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);
              double r2 = dx * dx + dy * dy + dz * dz;

              double soft = get_default_softening_of_particletype(0) * All.cf_atime;

              double r = sqrt(r2);
              if(r > 0.0)
                {
                  double dir[3];
                  dir[0] = -dx / r;
                  dir[1] = -dy / r;
                  dir[2] = -dz / r;

                  double mdotn = m[0] * dir[0] + m[1] * dir[1] + m[2] * dir[2];
                  SphP[j].BConserved[0] += (3.0 * mdotn * dir[0] - m[0]) / (r * r2 + soft * soft * soft) * SphP[j].Volume;
                  SphP[j].BConserved[1] += (3.0 * mdotn * dir[1] - m[1]) / (r * r2 + soft * soft * soft) * SphP[j].Volume;
                  SphP[j].BConserved[2] += (3.0 * mdotn * dir[2] - m[2]) / (r * r2 + soft * soft * soft) * SphP[j].Volume;
                }
              else
                {
                  double dir[3];
                  /* B field is aligned with the magnetic moment vector at the very center of the dipole */
                  dir[0] = m[0] / mstrength;
                  dir[1] = m[1] / mstrength;
                  dir[2] = m[2] / mstrength;

                  double mdotn = m[0] * dir[0] + m[1] * dir[1] + m[2] * dir[2];
                  SphP[j].BConserved[0] += (3.0 * mdotn * dir[0] - m[0]) / (soft * soft * soft) * SphP[j].Volume;
                  SphP[j].BConserved[1] += (3.0 * mdotn * dir[1] - m[1]) / (soft * soft * soft) * SphP[j].Volume;
                  SphP[j].BConserved[2] += (3.0 * mdotn * dir[2] - m[2]) / (soft * soft * soft) * SphP[j].Volume;
                }
            }
#endif

#ifdef TRACER_MC
          if(in->Mass > 0)
            {
              double prob = dm_total / in->Mass;
              pcum2 += prob;

              while(kpoint < in->NumberOfTracers && (prob_tracers[kpoint] < pcum1))
                kpoint++;

              while(kpoint < in->NumberOfTracers && (pcum1 < prob_tracers[kpoint] && prob_tracers[kpoint] < pcum2))
                {
                  if(out.ntracertargets < MAXLEN_TRACERTARGET_LIST)
                    {
                      int n                                = out.ntracertargets++;
                      out.tracertargets[n].attach_to_task  = ThisTask;
                      out.tracertargets[n].attach_to_index = j;
#ifdef TRACER_MC_CHECKS
                      out.tracertargets[n].attach_to_ID = P[j].ID;
#endif
                    }
                  else
                    warn("reached MAXLEN_TRACERTARGET_LIST");

                  kpoint++;
                }

              pcum1 = pcum2;
            }
#endif
        }
    }

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

#ifdef TRACER_MC
  myfree(prob_tracers);
#endif

  return 0;
}

#endif
