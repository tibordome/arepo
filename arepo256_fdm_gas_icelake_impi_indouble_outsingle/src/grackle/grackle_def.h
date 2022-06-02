/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/grackle/grackle_proto.h
 * \date        MM/YYYY
 * \author     	Matthew C Smith
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#ifndef GRACKLE_H
#define GRACKLE_H

#if defined(COOLING) && defined(GRACKLE)

#define CONFIG_BFLOAT_8

#ifdef __cplusplus
extern "C"
{
#include <grackle.h>
}
#else
#include <grackle.h>
#endif

#ifdef GRACKLE_D
#define GRACKLE_SPECIES_NUMBER 11
#define GRACKLE_H2
#elif defined(GRACKLE_H2)
#define GRACKLE_SPECIES_NUMBER 8
#elif !defined(GRACKLE_TAB)
#define GRACKLE_SPECIES_NUMBER 5
#else
#define GRACKLE_SPECIES_NUMBER 0
#endif

MyFloat get_grackle_pressure(int i);
void initialise_grackle(void);
void cooling_only(void);
void cool_active_cells(void);
double get_temp_individual_cell_grackle(int i);
double get_cooling_time_individual_cell_grackle(int i);
double grackle_get_timestep(int i);
void grackle_initialise_abundances(void);
void grackle_converge_abundances(void);

#endif
#endif
