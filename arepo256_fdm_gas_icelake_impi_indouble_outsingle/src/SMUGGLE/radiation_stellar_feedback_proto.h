/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/SMUGGLE/radiation_stellar_feedback_proto.h
 * \date        03/2020
 * \author
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#if defined(SMUGGLE_RADIATION_FEEDBACK)
MyFloat compute_stromgren_radius(double lum, double avg_gas_dens, double minGasDist, double hsml);
#ifdef SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION
MyFloat compute_stromgren_mass(double lum, double avg_gas_dens, double hsml);
#endif
MyFloat update_radiation_cooling_shutoff_time(int i, MyFloat dt);
#endif

#ifdef SMUGGLE_RADIATION_FEEDBACK
void do_radiation_stellar_feedback(void);
int radiation_feedback_evaluate(int target, int mode, int thread_id);
MyFloat compute_radiation_feedback_velocity(int n, MyFloat *stromgren_mass, MyFloat stromgren_radius, MyFloat radiatmom, MyFloat slenght);
void find_radiation_feedback_cells(void);
int find_radiation_feedback_cells_evaluate(int target, int mode, int thread_id);
#endif
