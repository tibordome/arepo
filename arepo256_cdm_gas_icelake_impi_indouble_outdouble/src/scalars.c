/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/scalars.c
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"

#ifdef MAXSCALARS
int N_Scalar = 0;
struct scalar_elements scalar_elements[MAXSCALARS];
struct scalar_index ScalarIndex;
#endif

/*! \brief Main routine to initialize passive scalar quantities.
 *
 *  \return void
 */
void init_scalars(void)
{
#ifdef MAXSCALARS

#if defined(EOS_DEGENERATE) || defined(EOS_OPAL) || defined(EOS_PASSIVE)
  for(int i = 0; i < EOS_NSPECIES; i++)
    {
      scalar_init(&SphP[0].Composition[i], &SphP[0].MassComposition[i], SCALAR_TYPE_SPECIES);
    }
#endif

#if defined(REFINEMENT_HIGH_RES_GAS) && !defined(TGSET)
  ScalarIndex.HighResMass = scalar_init(&SphP[0].HighResDensity, &SphP[0].HighResMass, SCALAR_TYPE_PASSIVE);
  if(ScalarIndex.HighResMass == -1)
    terminate("ScalarIndex.HighResMass initialized incorrectly\n");
#endif

#if defined(REFINEMENT_CGM)
  ScalarIndex.HighResMassCGM = scalar_init(&SphP[0].HighResDensityCGM, &SphP[0].HighResMassCGM, SCALAR_TYPE_PASSIVE);
#endif

#ifdef METALS
  scalar_init(&SphP[0].Metallicity, &SphP[0].MassMetallicity, SCALAR_TYPE_PASSIVE);
#endif

#if defined(MRT) && !defined(MRT_NO_UV)
#ifndef MRT_EQUIL_CHEM_COOL
  scalar_init(&SphP[0].HI, &SphP[0].nHI, SCALAR_TYPE_HYDROGEN);
  scalar_init(&SphP[0].HII, &SphP[0].nHII, SCALAR_TYPE_HYDROGEN);
  //  scalar_init(&SphP[0].Ne, &SphP[0].ne, SCALAR_TYPE_PASSIVE);
#ifdef MRT_INCLUDE_HE
  scalar_init(&SphP[0].HeI, &SphP[0].nHeI, SCALAR_TYPE_HELIUM);
  scalar_init(&SphP[0].HeII, &SphP[0].nHeII, SCALAR_TYPE_HELIUM);
  scalar_init(&SphP[0].HeIII, &SphP[0].nHeIII, SCALAR_TYPE_HELIUM);
#endif
#endif
#endif

#ifdef SGCHEM
#if CHEMISTRYNETWORK == 1
  scalar_init(&SphP[0].TracAbund[IH2], &SphP[0].MassTracAbund[IH2], SCALAR_TYPE_H2);
  scalar_init(&SphP[0].TracAbund[IHP], &SphP[0].MassTracAbund[IHP], SCALAR_TYPE_HYDROGEN);
  scalar_init(&SphP[0].TracAbund[IDP], &SphP[0].MassTracAbund[IDP], SCALAR_TYPE_DEUTERIUM);
  scalar_init(&SphP[0].TracAbund[IHD], &SphP[0].MassTracAbund[IHD], SCALAR_TYPE_HD);
  scalar_init(&SphP[0].TracAbund[IHEP], &SphP[0].MassTracAbund[IHEP], SCALAR_TYPE_HELIUM);
  scalar_init(&SphP[0].TracAbund[IHEPP], &SphP[0].MassTracAbund[IHEPP], SCALAR_TYPE_HELIUM);
  scalar_init(&SphP[0].TracAbund[IHATOM], &SphP[0].MassTracAbund[IHATOM], SCALAR_TYPE_HYDROGEN);
  scalar_init(&SphP[0].TracAbund[IDATOM], &SphP[0].MassTracAbund[IDATOM], SCALAR_TYPE_DEUTERIUM);
  scalar_init(&SphP[0].TracAbund[IHEATOM], &SphP[0].MassTracAbund[IHEATOM], SCALAR_TYPE_HELIUM);
#endif
#if CHEMISTRYNETWORK == 4
  scalar_init(&SphP[0].TracAbund[IH2], &SphP[0].MassTracAbund[IH2], SCALAR_TYPE_H2);
  scalar_init(&SphP[0].TracAbund[IHP], &SphP[0].MassTracAbund[IHP], SCALAR_TYPE_HYDROGEN);
  scalar_init(&SphP[0].TracAbund[IHATOM], &SphP[0].MassTracAbund[IHATOM], SCALAR_TYPE_HYDROGEN);
#endif
#if CHEMISTRYNETWORK == 5 || CHEMISTRYNETWORK == 7
  scalar_init(&SphP[0].TracAbund[IH2], &SphP[0].MassTracAbund[IH2], SCALAR_TYPE_H2);
  scalar_init(&SphP[0].TracAbund[IHP], &SphP[0].MassTracAbund[IHP], SCALAR_TYPE_HYDROGEN);
  scalar_init(&SphP[0].TracAbund[IHATOM], &SphP[0].MassTracAbund[IHATOM], SCALAR_TYPE_HYDROGEN);
  scalar_init(&SphP[0].TracAbund[ICP], &SphP[0].MassTracAbund[ICP], SCALAR_TYPE_CARBON);
  scalar_init(&SphP[0].TracAbund[ICO], &SphP[0].MassTracAbund[ICO], SCALAR_TYPE_CARBON);
#if CHEMISTRYNETWORK == 7
  scalar_init(&SphP[0].TracAbund[IHEP], &SphP[0].MassTracAbund[IHEP], SCALAR_TYPE_HELIUM);
  scalar_init(&SphP[0].TracAbund[IHEATOM], &SphP[0].MassTracAbund[IHEATOM], SCALAR_TYPE_HELIUM);
#endif
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
  scalar_init(&SphP[0].TracAbund[IH2], &SphP[0].MassTracAbund[IH2], SCALAR_TYPE_MCMA);
  scalar_init(&SphP[0].TracAbund[IHP], &SphP[0].MassTracAbund[IHP], SCALAR_TYPE_MCMA);
  scalar_init(&SphP[0].TracAbund[ICP], &SphP[0].MassTracAbund[ICP], SCALAR_TYPE_MCMA);
  scalar_init(&SphP[0].TracAbund[ICHX], &SphP[0].MassTracAbund[ICHX], SCALAR_TYPE_MCMA);
  scalar_init(&SphP[0].TracAbund[IOHX], &SphP[0].MassTracAbund[IOHX], SCALAR_TYPE_MCMA);
  scalar_init(&SphP[0].TracAbund[ICO], &SphP[0].MassTracAbund[ICO], SCALAR_TYPE_MCMA);
  scalar_init(&SphP[0].TracAbund[IHCOP], &SphP[0].MassTracAbund[IHCOP], SCALAR_TYPE_MCMA);
  scalar_init(&SphP[0].TracAbund[IHEP], &SphP[0].MassTracAbund[IHEP], SCALAR_TYPE_HELIUM);
  scalar_init(&SphP[0].TracAbund[IMP], &SphP[0].MassTracAbund[IMP], SCALAR_TYPE_METAL);
  scalar_init(&SphP[0].TracAbund[IHATOM], &SphP[0].MassTracAbund[IHATOM], SCALAR_TYPE_MCMA);
  scalar_init(&SphP[0].TracAbund[IHEATOM], &SphP[0].MassTracAbund[IHEATOM], SCALAR_TYPE_HELIUM);
  scalar_init(&SphP[0].TracAbund[ICATOM], &SphP[0].MassTracAbund[ICATOM], SCALAR_TYPE_MCMA);
  scalar_init(&SphP[0].TracAbund[IOATOM], &SphP[0].MassTracAbund[IOATOM], SCALAR_TYPE_MCMA);
  scalar_init(&SphP[0].TracAbund[IMATOM], &SphP[0].MassTracAbund[IMATOM], SCALAR_TYPE_METAL);
#endif

#ifdef SGCHEM_VARIABLE_Z
  scalar_init(&SphP[0].CarbAbund, &SphP[0].CarbMass, SCALAR_TYPE_CELEM);
  scalar_init(&SphP[0].OxyAbund, &SphP[0].OxyMass, SCALAR_TYPE_OELEM);
  scalar_init(&SphP[0].MAbund, &SphP[0].MMass, SCALAR_TYPE_MELEM);
  scalar_init(&SphP[0].ZAtom, &SphP[0].ZMass, SCALAR_TYPE_PASSIVE);
  scalar_init(&SphP[0].DustToGasRatio, &SphP[0].ScaledDustMass, SCALAR_TYPE_PASSIVE);
#endif
#endif

#ifdef REFINEMENT_RPS
  scalar_init(&SphP[0].RPSGalaxyDensity, &SphP[0].RPSGalaxyMass, SCALAR_TYPE_PASSIVE);
#endif

#ifdef GFM_STELLAR_EVOLUTION
  ScalarIndex.Metallicity = scalar_init(&SphP[0].Metallicity, &SphP[0].MassMetallicity, SCALAR_TYPE_PASSIVE);
  if(ScalarIndex.Metallicity == -1)
    terminate("ScalarIndex.Metallicity initialized incorrectly\n");
#ifdef GFM_NORMALIZED_METAL_ADVECTION
  for(int i = 0; i < GFM_N_CHEM_ELEMENTS; i++)
    scalar_init(&SphP[0].MetalsFraction[i], &SphP[0].MassMetals[i], SCALAR_TYPE_SPECIES);
#else
  for(int i = 0; i < GFM_N_CHEM_ELEMENTS; i++)
    scalar_init(&SphP[0].MetalsFraction[i], &SphP[0].MassMetals[i], SCALAR_TYPE_PASSIVE);
#endif
#ifdef GFM_DUST
  /* Keep dust as SCALAR_TYPE_PASSIVE since the normalization constraint from
   * SCALAR_TYPE_SPECIES is not applicable like for metals. */
  for(int i = 0; i < GFM_DUST_N_CHANNELS; i++)
    {
      for(int j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
        {
          scalar_init(&SphP[0].MetalsDustFraction[i][j], &SphP[0].MassMetalsDust[i][j], SCALAR_TYPE_PASSIVE);
        }
    }
#endif
#endif

#ifdef GFM_RPROCESS_CHANNELS
  for(int i = 0; i < GFM_RPROCESS_CHANNELS; i++)
    {
      scalar_init(&SphP[0].MassFractionRProcess[i], &SphP[0].MassRProcess[i], SCALAR_TYPE_PASSIVE);
    }
#endif

#if defined(GRACKLE) && !defined(GRACKLE_TAB)
  for(int i = 0; i < GRACKLE_SPECIES_NUMBER; i++)
    scalar_init(&SphP[0].GrackleSpeciesFraction[i], &SphP[0].GrackleSpeciesMass[i], SCALAR_TYPE_PASSIVE);

  scalar_init(&SphP[0].e_frac, &SphP[0].e_mass, SCALAR_TYPE_PASSIVE);
#endif

#ifdef PASSIVE_SCALARS
  for(int i = 0; i < PASSIVE_SCALARS; i++)
    {
      scalar_init(&SphP[0].PScalars[i], &SphP[0].PConservedScalars[i], SCALAR_TYPE_PASSIVE);
    }
#endif

#ifdef COSMIC_RAYS
  ScalarIndex.CR_Energy = scalar_init(&SphP[0].CR_SpecificEnergy, &SphP[0].CR_Energy, SCALAR_TYPE_PASSIVE);
  if(ScalarIndex.CR_Energy == -1)
    terminate("ScalarIndex.CR_Energy initialized incorrectly\n");
#endif

#ifdef MHD_THERMAL_ENERGY_SWITCH
  ScalarIndex.Etherm = scalar_init(&SphP[0].Utherm, &SphP[0].Etherm, SCALAR_TYPE_PASSIVE);
#endif

#ifdef SGS_TURBULENCE
  ScalarIndex.SGS_T_Energy = scalar_init(&SphP[0].SgsTData.SpecificEnergy, &SphP[0].SgsTData.Energy, SCALAR_TYPE_PASSIVE);
  if(ScalarIndex.SGS_T_Energy == -1)
    terminate("ScalarIndex.SGS_T_Energy initialized incorrectly\n");
#endif

#ifdef TURB_APPROX_MCS
  scalar_init(&SphP[0].TurbSpecEnergy, &SphP[0].TurbEnergy, SCALAR_TYPE_PASSIVE);
#endif

  mpi_printf("INIT: %d/%d Scalars used.\n", N_Scalar, MAXSCALARS);
#endif /* MAXSCALARS */
}

/*! \brief Initialize a specific scalar property.
 *
 *  \param[in] addr Pointer to (primitive) scalar in SphP[0] struct.
 *  \param[in] addr_mass Pointer to conserved scalar quantity in SphP[0].
 *  \param[in] type Type of scalar (e.g. SCALAR_TYPE_PASSIVE for passive
 *             scalar)
 *
 *  \return Number of scalars - 1
 */
int scalar_init(MyFloat *addr, MyFloat *addr_mass, int type)
{
#ifdef MAXSCALARS
  if(N_Scalar == MAXSCALARS)
    {
      mpi_printf("Failed to register scalar, maximum of %d already reached\n", MAXSCALARS);
      terminate("MAXSCALARS reached");
    }

  /* save type and relative address */
  scalar_elements[N_Scalar].type        = type;
  scalar_elements[N_Scalar].offset      = ((char *)addr) - ((char *)&SphP[0]);
  scalar_elements[N_Scalar].offset_mass = ((char *)addr_mass) - ((char *)&SphP[0]);

  N_Scalar++;

  return N_Scalar - 1;
  /* note: gradients are initialized in init_gradients */
#else
  return -1;
#endif
}
