/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sinks/sinks_proto.h
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       function prototypes for sink particles
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#ifdef SINKS

void sinks(void);

void sinks_dmass(void);
void sinks_accrete(void);

void sinks_set_constants(void);
void write_sink_data(int xaxis, int yaxis, int zaxis, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
double get_accretion_radius(double mass, double temp);

void sinks_create(void);

void sinks_begrun(void);

void sinks_get_num_sinks(void);
void sinks_get_active_sinks(void);
void sinks_begin(void);
void sinks_end(void);

void sinks_check_AuxDataID_references(void);

#ifdef SINKS_MERGERS
void sinks_mergers(void);
void sinks_find_neighboring_sinks(void);
void sinks_do_mergers(void);
#endif

#endif
