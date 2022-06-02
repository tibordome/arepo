/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dmwindtunnel/dmwindtunnel.h
 * \date        12/2020
 * \author      Philip Mocz, Martin Sparre, Hayden Foote
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - 08.02.2021 HF added functions for reading wind parameters from a file
 * - 01.11.2022 HF added star particles to the windtunnel
 */

#ifndef DM_WINDTUNNEL_H
#define DM_WINDTUNNEL_H

#if(defined(DM_WINDTUNNEL) && defined(DM_WINDTUNNEL_EXTERNAL_SOURCE))
/* DM windtunnel parameters from an external file*/
void interpolate_from_dm_wind_table(double t, double *rho, double *vel, double *sigma);
void read_dmwindtunnel_file(void);
#endif

#ifdef DM_WINDTUNNEL

/* DM Windtunnel Main Methods */
void apply_windtunnel_bcs(void);

#endif

#endif /* DM_WINDTUNNEL_H */
